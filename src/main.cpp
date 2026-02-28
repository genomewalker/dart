#include "cli/args.hpp"
#include "cli/subcommand.hpp"
#include "dart/damage_model.hpp"
#include "dart/adaptive_damage.hpp"
#include "dart/sequence_io.hpp"
#include "dart/frame_selector.hpp"
#include "dart/hexamer_tables.hpp"
#include "dart/codon_tables.hpp"
#include "dart/unified_codon_scorer.hpp"
#include "dart/damage_index_writer.hpp"
#include "dart/log_utils.hpp"
#include "dart/version.h"
#include <iostream>
#include <fstream>
#include <chrono>
#include <iomanip>
#include <mutex>
#include <atomic>
#include <vector>
#include <cstring>
#include <unistd.h>
#include <sstream>

#ifdef _OPENMP
#include <omp.h>
#endif

using Options = dart::cli::Options;

// Convert CLI library type to SampleDamageProfile library type
inline dart::SampleDamageProfile::LibraryType to_sample_library_type(dart::cli::LibraryType lt) {
    switch (lt) {
        case dart::cli::LibraryType::DOUBLE_STRANDED:
            return dart::SampleDamageProfile::LibraryType::DOUBLE_STRANDED;
        case dart::cli::LibraryType::SINGLE_STRANDED:
            return dart::SampleDamageProfile::LibraryType::SINGLE_STRANDED;
        default:
            return dart::SampleDamageProfile::LibraryType::UNKNOWN;
    }
}

namespace dart {
namespace cli {

// Thread-local output buffer for parallel gene prediction
struct ThreadOutputBuffer {
    std::ostringstream gff;
    std::ostringstream fasta_nt;
    std::ostringstream fasta_aa;
    std::ostringstream fasta_aa_masked;
    std::vector<std::tuple<std::string, Gene, std::string>> damage_records;
    size_t gene_count = 0;

    void clear() {
        gff.str("");
        gff.clear();
        fasta_nt.str("");
        fasta_nt.clear();
        fasta_aa.str("");
        fasta_aa.clear();
        fasta_aa_masked.str("");
        fasta_aa_masked.clear();
        damage_records.clear();
        gene_count = 0;
    }

    // Format GFF line for a gene
    void write_gff(const std::string& seq_id, const Gene& gene, int gene_num,
                   const std::string& source = "AncientGenePredictor") {
        gff << seq_id << "\t"
            << source << "\t"
            << "CDS\t"
            << (gene.start + 1) << "\t"
            << gene.end << "\t"
            << std::fixed << std::setprecision(3) << gene.score << "\t"
            << (gene.is_forward ? "+" : "-") << "\t"
            << gene.frame << "\t"
            << "ID=" << seq_id << "_gene" << gene_num
            << ";damage_pct=" << std::fixed << std::setprecision(1) << gene.damage_score
            << ";p_damaged=" << std::fixed << std::setprecision(3) << gene.p_read_damaged
            << ";aa_corr=" << gene.aa_corrections
            << ";stop_fix=" << gene.stop_restorations
            << "\n";
    }

    // Format FASTA entry for nucleotide
    void write_fasta_nt(const std::string& seq_id, const Gene& gene, int gene_num) {
        fasta_nt << ">" << seq_id << "_" << (gene.is_forward ? "+" : "-") << "_" << gene.frame
                 << " " << (gene.start + 1) << ".." << gene.end
                 << " length=" << gene.sequence.length() << "nt"
                 << "\n" << gene.sequence << "\n";
    }

    // Format FASTA entry for protein
    void write_fasta_aa(const std::string& seq_id, const Gene& gene, int gene_num) {
        fasta_aa << ">" << seq_id << "_" << (gene.is_forward ? "+" : "-") << "_" << gene.frame
                 << " " << (gene.start + 1) << ".." << gene.end
                 << " length=" << gene.protein.length() << "aa"
                 << " damage_pct=" << static_cast<int>(gene.damage_score)
                 << " p_damaged=" << static_cast<int>(gene.p_read_damaged * 100) / 100.0
                 << " log_lr=" << static_cast<int>(gene.damage_signal * 100) / 100.0
                 << " dmg_info=" << gene.dna_corrections
                 << " aa_corr=" << gene.aa_corrections
                 << " stop_fix=" << gene.stop_restorations
                 << " corr_score=" << gene.orf_rank
                 << "\n" << gene.protein << "\n";
    }

    // Format FASTA entry for search protein (X-masked)
    void write_fasta_aa_masked(const std::string& seq_id, const Gene& gene, int gene_num) {
        const std::string& prot = gene.search_protein.empty() ? gene.protein : gene.search_protein;
        fasta_aa_masked << ">" << seq_id << "_" << (gene.is_forward ? "+" : "-") << "_" << gene.frame
                        << " " << (gene.start + 1) << ".." << gene.end
                        << " length=" << prot.length() << "aa"
                        << " x_count=" << gene.x_count
                        << "\n" << prot << "\n";
    }
};

int cmd_predict(int argc, char* argv[]) {
    try {
        auto run_start = std::chrono::steady_clock::now();
        Options opts = dart::cli::parse_args(argc, argv);

        // Set up threading
        int num_threads = opts.num_threads;
#ifdef _OPENMP
        if (num_threads == 0) {
            num_threads = omp_get_max_threads();
        }
        omp_set_num_threads(num_threads);
#else
        num_threads = 1;
#endif

        bool is_tty = isatty(fileno(stderr));

        // Set active domain for hexamer scoring
        dart::Domain active_domain = dart::parse_domain(opts.domain_name);
        dart::set_active_domain(active_domain);

        // Header
        std::cerr << "DART v" << DART_VERSION << "\n";
        std::cerr << "Input: " << opts.input_file << "\n";
        std::cerr << "Output: " << opts.output_file << "\n";
        std::cerr << "Domain: " << dart::domain_name(active_domain) << "\n";
        std::cerr << "Threads: " << num_threads << "\n";

        if (opts.verbose) {
            std::cerr << "Min ORF length: " << opts.orf_min_aa << " aa\n";
            std::cerr << "ORF mode: " << (opts.adaptive_orf ? "adaptive" : "6-frame") << "\n";
            std::cerr << "Damage detection: " << (opts.use_damage ? "ON" : "OFF") << "\n";
            if (opts.aggregate_damage) std::cerr << "Aggregate damage: ON (two-pass)\n";
        }
        std::cerr << "\n";

        // Initialize damage model
        dart::DamageModel damage_model;

        // Open output files - using direct streams for parallel-friendly writes
        std::ofstream gff_file(opts.output_file);
        if (!gff_file) {
            throw std::runtime_error("Cannot open output file: " + opts.output_file);
        }
        // Write GFF3 header
        gff_file << "##gff-version 3\n";

        // Direct file streams for parallel output buffering
        std::ofstream fasta_nt_file;
        std::ofstream fasta_aa_file;
        std::ofstream fasta_aa_masked_file;

        // Flags to check if outputs are enabled
        bool write_fasta_nt = !opts.fasta_nt.empty();
        bool write_fasta_aa = !opts.fasta_aa.empty();
        bool write_fasta_aa_masked = !opts.fasta_aa_masked.empty();

        // Keep corrected writer for backward compatibility (not optimized)
        std::unique_ptr<dart::FastaWriter> fasta_nt_corr_writer;

        if (write_fasta_nt) {
            fasta_nt_file.open(opts.fasta_nt);
            if (!fasta_nt_file) {
                throw std::runtime_error("Cannot open FASTA NT file: " + opts.fasta_nt);
            }
        }
        if (!opts.fasta_nt_corrected.empty()) {
            fasta_nt_corr_writer = std::make_unique<dart::FastaWriter>(opts.fasta_nt_corrected);
        }
        if (write_fasta_aa) {
            fasta_aa_file.open(opts.fasta_aa);
            if (!fasta_aa_file) {
                throw std::runtime_error("Cannot open FASTA AA file: " + opts.fasta_aa);
            }
        }
        if (write_fasta_aa_masked) {
            fasta_aa_masked_file.open(opts.fasta_aa_masked);
            if (!fasta_aa_masked_file) {
                throw std::runtime_error("Cannot open FASTA AA masked file: " + opts.fasta_aa_masked);
            }
        }

        // Sample-level damage profile
        dart::SampleDamageProfile sample_profile;
        sample_profile.forced_library_type = to_sample_library_type(opts.forced_library_type);

        // Adaptive damage calibrator
        dart::AdaptiveDamageCalibrator calibrator;

        // Pass 1: Damage detection
        if (opts.aggregate_damage && opts.use_damage) {
            auto pass1_start = std::chrono::high_resolution_clock::now();
            std::cerr << "[Pass 1] Scanning for damage patterns...\n";

            dart::SequenceReader first_pass(opts.input_file);
            const size_t BATCH_SIZE = 50000;
            std::vector<dart::SequenceRecord> batch;
            batch.reserve(BATCH_SIZE);
            size_t count = 0;

            // Thread-local profiles
            std::vector<dart::SampleDamageProfile> thread_profiles(num_threads);
            for (auto& tp : thread_profiles) {
                tp.forced_library_type = to_sample_library_type(opts.forced_library_type);
            }

            // EM training sequences
            const size_t EM_MAX_SEQS = 50000;
            std::vector<std::string> em_training_seqs(EM_MAX_SEQS);
            std::atomic<size_t> em_training_count{0};

            while (true) {
                batch.clear();
                dart::SequenceRecord record;
                while (batch.size() < BATCH_SIZE && first_pass.read_next(record)) {
                    batch.push_back(std::move(record));
                }
                if (batch.empty()) break;

                #pragma omp parallel for schedule(static)
                for (size_t i = 0; i < batch.size(); ++i) {
                    int tid = 0;
                    #ifdef _OPENMP
                    tid = omp_get_thread_num();
                    #endif
                    std::string seq = dart::SequenceUtils::clean(batch[i].sequence);
                    if (seq.length() >= opts.min_length) {
                        dart::FrameSelector::update_sample_profile(thread_profiles[tid], seq);

                        if (seq.length() >= 60 && seq.length() <= 200) {
                            size_t slot = em_training_count.fetch_add(1, std::memory_order_relaxed);
                            if (slot < EM_MAX_SEQS) {
                                em_training_seqs[slot] = seq;
                            }
                        }
                    }
                }

                count += batch.size();
                if (count % 1000000 < BATCH_SIZE) {
                    std::cerr << "  Scanned " << count / 1000000 << "M sequences...\n";
                }
            }

            // Merge thread profiles
            for (int t = 0; t < num_threads; ++t) {
                dart::FrameSelector::merge_sample_profiles(sample_profile, thread_profiles[t]);
            }

            size_t kept_training = std::min(em_training_count.load(std::memory_order_relaxed), EM_MAX_SEQS);
            em_training_seqs.resize(kept_training);

            dart::FrameSelector::finalize_sample_profile(sample_profile);

            // =========================================================================
            // Pass 1.5: Compute d_metamatch for high-damage samples
            // This re-weights reads by P(damaged) to better estimate the damage rate
            // in the ancient DNA fraction, matching metaDMG's selection bias.
            // =========================================================================
            if (sample_profile.damage_validated && sample_profile.d_max_combined > 0.15f) {
                bool needs_metamatch = false;
                if (sample_profile.inverted_pattern_5prime || sample_profile.inverted_pattern_3prime) {
                    needs_metamatch = true;
                } else if (sample_profile.position_0_artifact_5prime) {
                    needs_metamatch = true;
                } else {
                    float gap = std::abs(sample_profile.d_max_from_channel_b - sample_profile.d_max_combined);
                    if (gap > 0.05f) {
                        needs_metamatch = true;
                    }
                }

                if (needs_metamatch) {
                    std::cerr << "  Computing d_metamatch (weighted profile)...\n";

                    dart::SequenceReader reader2(opts.input_file);
                    std::vector<dart::SampleDamageProfile> thread_profiles2(num_threads);
                    std::vector<dart::SequenceRecord> batch2;
                    batch2.reserve(50000);

                    while (true) {
                        batch2.clear();
                        dart::SequenceRecord rec;
                        while (batch2.size() < 50000 && reader2.read_next(rec)) {
                            batch2.push_back(std::move(rec));
                        }
                        if (batch2.empty()) break;

                        #pragma omp parallel for schedule(static)
                        for (size_t i = 0; i < batch2.size(); ++i) {
                            int tid = 0;
                            #ifdef _OPENMP
                            tid = omp_get_thread_num();
                            #endif
                            std::string seq = dart::SequenceUtils::clean(batch2[i].sequence);
                            if (seq.length() >= 30) {
                                float p_damaged = dart::FrameSelector::compute_per_read_damage_prior(seq, sample_profile);
                                if (p_damaged > 0.01f) {
                                    dart::FrameSelector::update_sample_profile_weighted(
                                        thread_profiles2[tid], seq, p_damaged);
                                }
                            }
                        }
                    }

                    dart::SampleDamageProfile weighted_profile;
                    for (int t = 0; t < num_threads; ++t) {
                        dart::FrameSelector::merge_sample_profiles(weighted_profile, thread_profiles2[t]);
                    }
                    dart::FrameSelector::finalize_sample_profile(weighted_profile);

                    // Update d_max_combined with d_metamatch for prediction
                    float d_metamatch = weighted_profile.d_max_combined;
                    std::cerr << "  d_metamatch: " << std::fixed << std::setprecision(1)
                              << (d_metamatch * 100.0f) << "% (was "
                              << (sample_profile.d_max_combined * 100.0f) << "%)\n";
                    sample_profile.d_metamatch = d_metamatch;
                    sample_profile.d_max_combined = d_metamatch;  // Use for prediction
                }
            }

            // Train periodic model
            if (!em_training_seqs.empty()) {
                dart::codon::DamageParams em_dmg;
                em_dmg.d_max_5p = sample_profile.max_damage_5prime;
                em_dmg.d_max_3p = sample_profile.max_damage_3prime;
                em_dmg.lambda_5p = sample_profile.lambda_5prime;
                em_dmg.lambda_3p = sample_profile.lambda_3prime;
                dart::codon::train_periodic_model(em_training_seqs, em_dmg, 3);
            }

            auto pass1_end = std::chrono::high_resolution_clock::now();
            auto pass1_duration = std::chrono::duration_cast<std::chrono::milliseconds>(pass1_end - pass1_start);

            std::cerr << "  Reads: " << sample_profile.n_reads
                     << " | 5' damage: " << std::fixed << std::setprecision(1)
                     << sample_profile.max_damage_5prime * 100.0f << "%"
                     << " | 3' damage: " << sample_profile.max_damage_3prime * 100.0f << "%"
                     << " | " << dart::log_utils::format_duration_ms(pass1_duration.count()) << "\n\n";

            damage_model.update_from_sample_profile(sample_profile);
            calibrator.initialize(sample_profile);
        }

        // Exit early if damage-only mode
        if (opts.damage_only) {
            std::cerr << "Done (damage-only mode).\n";

            // Write summary if requested
            if (!opts.summary_file.empty()) {
                std::ofstream summary(opts.summary_file);
                summary << "{\n";
                summary << "  \"d_max_5prime\": " << sample_profile.max_damage_5prime << ",\n";
                summary << "  \"d_max_3prime\": " << sample_profile.max_damage_3prime << ",\n";
                summary << "  \"lambda_5prime\": " << sample_profile.lambda_5prime << ",\n";
                summary << "  \"lambda_3prime\": " << sample_profile.lambda_3prime << ",\n";
                summary << "  \"n_reads\": " << sample_profile.n_reads << "\n";
                summary << "}\n";
            }
            return 0;
        }

        // Pass 2: Gene prediction via ORF enumeration
        std::cerr << "[Pass 2] Gene prediction...\n";
        auto pass2_start = std::chrono::high_resolution_clock::now();

        dart::SequenceReader reader(opts.input_file);
        const size_t BATCH_SIZE = 10000;
        std::vector<dart::SequenceRecord> batch;
        batch.reserve(BATCH_SIZE);

        std::atomic<size_t> total_seqs{0};
        std::atomic<size_t> total_genes{0};

        // Optional damage index writer for post-mapping annotation
        std::unique_ptr<dart::DamageIndexWriter> damage_index_writer;
        if (!opts.damage_index.empty()) {
            damage_index_writer = std::make_unique<dart::DamageIndexWriter>(
                opts.damage_index, sample_profile);
        }

        // Thread-local output buffers for parallel formatting
        std::vector<ThreadOutputBuffer> thread_buffers(num_threads);

        while (true) {
            batch.clear();
            dart::SequenceRecord record;
            while (batch.size() < BATCH_SIZE && reader.read_next(record)) {
                batch.push_back(std::move(record));
            }
            if (batch.empty()) break;

            // Clear thread buffers for this batch
            for (auto& buf : thread_buffers) {
                buf.clear();
            }

            // Process batch in parallel - format output to thread-local buffers
            #pragma omp parallel for schedule(dynamic, 100)
            for (size_t i = 0; i < batch.size(); ++i) {
                std::string seq = dart::SequenceUtils::clean(batch[i].sequence);
                if (seq.length() < opts.min_length) continue;

                std::string id = batch[i].id;

                // Compute per-read damage prior BEFORE ORF enumeration
                // This allows damage-aware scoring during ORF selection
                float per_read_damage_prior = opts.use_damage
                    ? dart::FrameSelector::compute_per_read_damage_prior(seq, sample_profile)
                    : 0.0f;

                // ORF enumeration with damage-aware scoring
                // Uses per_read_damage_prior to adjust stop penalties and X-mask thresholds
                auto orfs = dart::FrameSelector::enumerate_orf_fragments(
                    seq, sample_profile,
                    opts.orf_min_aa,
                    opts.adaptive_orf,
                    per_read_damage_prior);

                float damage_pct = opts.use_damage
                    ? dart::FrameSelector::compute_damage_percentage(seq, sample_profile)
                    : 0.0f;

                // Refine per-read damage using ORF-specific evidence (stops, pre-stops)
                float per_read_damage = dart::FrameSelector::infer_per_read_aa_damage(seq, orfs, sample_profile);

                // Get thread-local buffer
                int tid = 0;
                #ifdef _OPENMP
                tid = omp_get_thread_num();
                #endif
                auto& buf = thread_buffers[tid];

                // Pre-compute reverse complement once (not per-ORF)
                std::string rc;
                int gene_num = 0;
                for (size_t orf_idx = 0; orf_idx < orfs.size(); ++orf_idx) {
                    const auto& orf = orfs[orf_idx];
                    // Use effective length (max of protein and search_protein) for filtering
                    // This allows short proteins that were extended via search_protein to pass
                    size_t effective_length = std::max(orf.protein.length(), orf.search_protein.length());
                    if (effective_length < opts.orf_min_aa) continue;

                    dart::Gene gene;
                    size_t L = seq.length();

                    if (orf.is_forward) {
                        gene.start = orf.nt_start;
                        gene.end = orf.nt_end;
                        gene.sequence = seq.substr(orf.nt_start, orf.nt_end - orf.nt_start);
                    } else {
                        gene.start = L - orf.nt_end;
                        gene.end = L - orf.nt_start;
                        if (rc.empty()) rc = dart::reverse_complement(seq);
                        gene.sequence = rc.substr(orf.nt_start, orf.nt_end - orf.nt_start);
                    }

                    gene.is_forward = orf.is_forward;
                    gene.is_fragment = true;
                    gene.frame = orf.frame;

                    float normalized_score = orf.length > 0
                        ? orf.score / static_cast<float>(orf.length)
                        : 0.0f;
                    gene.frame_score = orf.score;
                    gene.coding_prob = std::max(0.0f, std::min(1.0f, normalized_score / 1.5f + 0.5f));
                    gene.score = gene.coding_prob;

                    gene.protein = orf.protein;
                    gene.search_protein = orf.search_protein;
                    gene.corrected_protein = orf.observed_protein;
                    gene.corrected_sequence = orf.corrected_nt;
                    gene.internal_stops = 0;
                    gene.damage_score = damage_pct;
                    gene.damage_signal = per_read_damage;
                    gene.p_read_damaged = per_read_damage;

                    gene.dna_corrections = orf.rescued_stops;
                    gene.aa_corrections = orf.aa_corrections;
                    gene.stop_restorations = orf.rescued_stops;

                    gene.orf_rank = orf_idx;
                    gene.passed_stops = orf.passed_stops;
                    gene.x_count = orf.x_count;
                    gene.orf_length_score = static_cast<float>(orf.length);

                    if (orfs.size() > 1 && orf_idx == 0) {
                        gene.frame_margin = static_cast<float>(orfs[0].length - orfs[1].length);
                    } else {
                        gene.frame_margin = 0.0f;
                    }
                    gene.frame_pmax = gene.coding_prob;

                    // Write to thread-local buffers (no locking needed)
                    buf.write_gff(id, gene, gene_num);
                    if (write_fasta_nt) buf.write_fasta_nt(id, gene, gene_num);
                    if (write_fasta_aa) buf.write_fasta_aa(id, gene, gene_num);
                    if (write_fasta_aa_masked) buf.write_fasta_aa_masked(id, gene, gene_num);
                    if (damage_index_writer) {
                        buf.damage_records.emplace_back(id, gene, gene.sequence);
                    }
                    buf.gene_count++;
                    gene_num++;
                }
            }

            // Merge and write thread buffers (sequential I/O, but no formatting overhead)
            for (const auto& buf : thread_buffers) {
                if (buf.gene_count == 0) continue;

                // Write pre-formatted strings directly to files
                gff_file << buf.gff.str();
                if (write_fasta_nt) {
                    fasta_nt_file << buf.fasta_nt.str();
                }
                if (write_fasta_aa) {
                    fasta_aa_file << buf.fasta_aa.str();
                }
                if (write_fasta_aa_masked) {
                    fasta_aa_masked_file << buf.fasta_aa_masked.str();
                }
                if (damage_index_writer) {
                    for (const auto& [id, gene, seq] : buf.damage_records) {
                        damage_index_writer->add_record(id, gene, seq);
                    }
                }
                total_genes += buf.gene_count;
            }

            total_seqs += batch.size();
            if (total_seqs % 1000000 < BATCH_SIZE) {
                std::cerr << "  Processed " << total_seqs / 1000000 << "M sequences...\n";
            }
        }

        // Finalize damage index if enabled
        if (damage_index_writer) {
            damage_index_writer->finalize();
            std::cerr << "  Damage index: " << damage_index_writer->record_count() << " records\n";
        }

        auto pass2_end = std::chrono::high_resolution_clock::now();
        auto pass2_duration = std::chrono::duration_cast<std::chrono::milliseconds>(pass2_end - pass2_start);

        std::cerr << "  Sequences: " << total_seqs.load()
                  << " | Genes: " << total_genes.load()
                  << " | " << dart::log_utils::format_duration_ms(pass2_duration.count())
                  << " (" << (total_seqs.load() * 1000 / std::max(1LL, (long long)pass2_duration.count())) << " seq/s)\n";

        // Write summary if requested
        if (!opts.summary_file.empty()) {
            std::ofstream summary(opts.summary_file);
            summary << "{\n";
            summary << "  \"version\": \"" << DART_VERSION << "\",\n";
            summary << "  \"input\": \"" << opts.input_file << "\",\n";
            summary << "  \"sequences\": " << total_seqs.load() << ",\n";
            summary << "  \"genes\": " << total_genes.load() << ",\n";
            summary << "  \"damage\": {\n";
            summary << "    \"d_max_5prime\": " << sample_profile.max_damage_5prime << ",\n";
            summary << "    \"d_max_3prime\": " << sample_profile.max_damage_3prime << ",\n";
            summary << "    \"lambda_5prime\": " << sample_profile.lambda_5prime << ",\n";
            summary << "    \"lambda_3prime\": " << sample_profile.lambda_3prime << "\n";
            summary << "  }\n";
            summary << "}\n";
        }

        auto run_end = std::chrono::steady_clock::now();
        std::cerr << "\nDone. Total runtime: "
                  << dart::log_utils::format_elapsed(run_start, run_end) << "\n";
        return 0;

    } catch (const dart::cli::ParseArgsExit& e) {
        if (e.exit_code() != 0) {
            if (std::strlen(e.what()) > 0) {
                std::cerr << e.what() << "\n\n";
            }
            dart::cli::print_usage(argv[0]);
        }
        return e.exit_code();
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
}

}  // namespace cli
}  // namespace dart

// Register the predict subcommand
namespace {
    struct PredictRegistrar {
        PredictRegistrar() {
            dart::cli::SubcommandRegistry::instance().register_command(
                "predict", "Predict genes from sequencing reads",
                dart::cli::cmd_predict, 20);
        }
    } predict_registrar;
}
