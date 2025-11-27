#include "agp/damage_model.hpp"
#include "agp/sequence_io.hpp"
#include "agp/frame_selector.hpp"
#include "agp/multi_domain_hexamer.hpp"
#include "agp/version.h"
#include <iostream>
#include <fstream>
#include <chrono>
#include <iomanip>
#include <mutex>
#include <vector>
#include <array>
#include <future>
#include <unistd.h>

#ifdef _OPENMP
#include <omp.h>
#endif

void print_version() {
    std::cout << "agp " << AGP_VERSION << "\n";
}

void print_usage(const char* program_name) {
    std::cout << "Ancient Gene Predictor v" << AGP_VERSION << "\n\n";
    std::cout << "Usage: " << program_name << " -i <input> -o <output> [options]\n";
    std::cout << "       " << program_name << " [options] <input.fasta|fastq>\n\n";
    std::cout << "Options:\n";
    std::cout << "  -i, --input <file>       Input FASTA/FASTQ file (or .gz)\n";
    std::cout << "  -o, --output <file>      Output GFF3 file (default: predictions.gff)\n";
    std::cout << "  --fasta-nt <file>        Output nucleotide FASTA\n";
    std::cout << "  --fasta-aa <file>        Output amino acid FASTA\n";
    std::cout << "  --min-length <int>       Minimum sequence length (default: 30)\n";
    std::cout << "  --min-coding-prob <f>    Minimum coding probability (default: 0.3)\n";
    std::cout << "  --threads <int>          Number of threads (default: auto)\n";
    std::cout << "  --no-damage              Disable damage detection\n";
    std::cout << "  --single-strand          Score only forward strand\n";
    std::cout << "  --both-strands           Output predictions for both strands\n";
    std::cout << "  --no-aggregate           Disable two-pass damage aggregation\n";
    std::cout << "  --stop-only              Use stop-codon-only frame selection (97% accuracy)\n";
    std::cout << "  --domain <name>          Taxonomic domain for hexamer scoring (default: gtdb)\n";
    std::cout << "                           Options: gtdb, fungi, plant, protozoa, invertebrate,\n";
    std::cout << "                                    viral, mammal, vertebrate, meta\n";
    std::cout << "                           Use 'meta' for weighted ensemble scoring across all domains\n";
    std::cout << "  -v, --verbose            Verbose output\n";
    std::cout << "  -V, --version            Show version and exit\n";
    std::cout << "  -h, --help               Show this help message\n";
    std::cout << "\n";
    std::cout << "Defaults:\n";
    std::cout << "  - Best strand output (1 prediction per read, use --both-strands for both)\n";
    std::cout << "  - Aggregate damage detection (two-pass for sample-level statistics)\n";
    std::cout << "  - damage_prob threshold 0.5 recommended (F1=0.986, 97% recall)\n";
    std::cout << "\n";
    std::cout << "Examples:\n";
    std::cout << "  # Standard usage (best strand + aggregate damage)\n";
    std::cout << "  " << program_name << " reads.fastq -o predictions.gff\n\n";
    std::cout << "  # Output both strands per read\n";
    std::cout << "  " << program_name << " --both-strands reads.fastq\n\n";
    std::cout << "  # Undamaged DNA (skip damage detection)\n";
    std::cout << "  " << program_name << " --no-damage reads.fastq\n\n";
}

struct Options {
    std::string input_file;
    std::string output_file = "predictions.gff";
    std::string fasta_nt;
    std::string fasta_aa;
    std::string domain_name = "gtdb";  // Default to bacteria/archaea for ancient DNA
    bool use_damage = true;
    size_t min_length = 30;
    float min_coding_prob = 0.3f;
    int num_threads = 0;
    bool verbose = false;
    bool dual_strand = true;       // Default ON for short aDNA
    bool best_strand_only = true;  // Default ON: output only best-scoring strand (cleaner output)
    bool aggregate_damage = true;  // Default ON for better damage detection
    bool stop_only = false;        // Use stop-codon-only frame selection (experimental)
    bool metagenome_mode = false;   // Use weighted ensemble scoring across all domains
};

Options parse_args(int argc, char* argv[]) {
    Options opts;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];

        if (arg == "-h" || arg == "--help") {
            print_usage(argv[0]);
            exit(0);
        } else if (arg == "-V" || arg == "--version") {
            print_version();
            exit(0);
        } else if (arg == "-i" || arg == "--input") {
            if (i + 1 < argc) opts.input_file = argv[++i];
        } else if (arg == "-o" || arg == "--output") {
            if (i + 1 < argc) opts.output_file = argv[++i];
        } else if (arg == "--fasta-nt") {
            if (i + 1 < argc) opts.fasta_nt = argv[++i];
        } else if (arg == "--fasta-aa") {
            if (i + 1 < argc) opts.fasta_aa = argv[++i];
        } else if (arg == "--no-damage") {
            opts.use_damage = false;
        } else if (arg == "--min-coding-prob") {
            if (i + 1 < argc) opts.min_coding_prob = std::stof(argv[++i]);
        } else if (arg == "--min-length") {
            if (i + 1 < argc) opts.min_length = std::stoul(argv[++i]);
        } else if (arg == "--threads") {
            if (i + 1 < argc) opts.num_threads = std::stoi(argv[++i]);
        } else if (arg == "-v" || arg == "--verbose") {
            opts.verbose = true;
        } else if (arg == "--single-strand") {
            opts.dual_strand = false;
        } else if (arg == "--dual-strand") {
            opts.dual_strand = true;
        } else if (arg == "--best-strand") {
            opts.best_strand_only = true;
        } else if (arg == "--both-strands") {
            opts.best_strand_only = false;
        } else if (arg == "--no-aggregate") {
            opts.aggregate_damage = false;
        } else if (arg == "--aggregate-damage") {
            opts.aggregate_damage = true;
        } else if (arg == "--stop-only") {
            opts.stop_only = true;
        } else if (arg == "--domain") {
            if (i + 1 < argc) {
                opts.domain_name = argv[++i];
                // "meta" enables weighted ensemble scoring across all domains
                if (opts.domain_name == "meta" || opts.domain_name == "metagenome" || opts.domain_name == "all") {
                    opts.metagenome_mode = true;
                }
            }
        } else if (arg[0] != '-') {
            opts.input_file = arg;
        } else {
            std::cerr << "Unknown option: " << arg << "\n";
            exit(1);
        }
    }

    if (opts.input_file.empty()) {
        std::cerr << "Error: No input file specified\n\n";
        print_usage(argv[0]);
        exit(1);
    }

    return opts;
}

int main(int argc, char* argv[]) {
    try {
        Options opts = parse_args(argc, argv);

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

        // Check if stderr is a TTY for inline progress updates
        bool is_tty = isatty(fileno(stderr));

        // Set active domain for hexamer scoring
        agp::Domain active_domain = agp::parse_domain(opts.domain_name);
        agp::set_active_domain(active_domain);

        // Always show basic info
        std::cerr << "Ancient Gene Predictor v" << AGP_VERSION << "\n";
        std::cerr << "Input: " << opts.input_file << "\n";
        std::cerr << "Output: " << opts.output_file << "\n";
        if (opts.metagenome_mode) {
            std::cerr << "Domain: metagenome (weighted ensemble scoring)\n";
        } else {
            std::cerr << "Domain: " << agp::domain_name(active_domain) << "\n";
        }
        std::cerr << "Threads: " << num_threads << "\n";

        if (opts.verbose) {
            std::cerr << "Min coding prob: " << opts.min_coding_prob << "\n";
            std::cerr << "Damage detection: " << (opts.use_damage ? "ON" : "OFF") << "\n";
            if (opts.aggregate_damage) std::cerr << "Aggregate damage: ON (two-pass)\n";
            if (opts.stop_only) std::cerr << "Frame selection: STOP-ONLY (97% accuracy mode)\n";
            if (opts.dual_strand) {
                std::cerr << "Dual-strand: ON";
                if (opts.best_strand_only) std::cerr << " (best only)";
                std::cerr << "\n";
            }
        }
        std::cerr << "\n";

        // Initialize damage model
        agp::DamageModel damage_model;

        // Open output files
        agp::GeneWriter writer(opts.output_file);
        std::unique_ptr<agp::FastaWriter> fasta_nt_writer;
        std::unique_ptr<agp::FastaWriter> fasta_aa_writer;

        if (!opts.fasta_nt.empty()) {
            fasta_nt_writer = std::make_unique<agp::FastaWriter>(opts.fasta_nt);
        }
        if (!opts.fasta_aa.empty()) {
            fasta_aa_writer = std::make_unique<agp::FastaWriter>(opts.fasta_aa);
        }

        // Sample-level damage profile
        agp::SampleDamageProfile sample_profile;

        // First pass: collect damage statistics (if aggregate mode)
        if (opts.aggregate_damage && opts.use_damage) {
            auto pass1_start = std::chrono::high_resolution_clock::now();
            std::cerr << "[Pass 1] Scanning for damage patterns...\n";

            agp::SequenceReader first_pass(opts.input_file);
            const size_t PASS1_BATCH_SIZE = 50000;  // Larger batches for Pass 1
            std::vector<agp::SequenceRecord> batch;
            batch.reserve(PASS1_BATCH_SIZE);
            size_t count = 0;

            // Thread-local profiles for parallel aggregation
            std::vector<agp::SampleDamageProfile> thread_profiles(num_threads);

            while (true) {
                // Read a batch
                batch.clear();
                agp::SequenceRecord record;
                while (batch.size() < PASS1_BATCH_SIZE && first_pass.read_next(record)) {
                    batch.push_back(std::move(record));
                }
                if (batch.empty()) break;

                // Process batch in parallel with thread-local profiles
                #pragma omp parallel for schedule(static)
                for (size_t i = 0; i < batch.size(); ++i) {
                    int tid = 0;
                    #ifdef _OPENMP
                    tid = omp_get_thread_num();
                    #endif
                    std::string seq = agp::SequenceUtils::clean(batch[i].sequence);
                    if (seq.length() >= opts.min_length) {
                        agp::FrameSelector::update_sample_profile(thread_profiles[tid], seq);
                    }
                }

                count += batch.size();
                if (count % 1000000 < PASS1_BATCH_SIZE) {
                    if (is_tty) {
                        std::cerr << "\r  Scanned " << count / 1000000 << "M sequences..." << std::flush;
                    } else {
                        std::cerr << "  Scanned " << count / 1000000 << "M sequences..." << std::endl;
                    }
                }
            }

            // Merge all thread-local profiles
            for (int t = 0; t < num_threads; ++t) {
                agp::FrameSelector::merge_sample_profiles(sample_profile, thread_profiles[t]);
            }

            // Store raw counts before finalization for debug
            double raw_t_5prime = sample_profile.t_freq_5prime[0];
            double raw_c_5prime = sample_profile.c_freq_5prime[0];
            double raw_a_3prime = sample_profile.a_freq_3prime[0];
            double raw_g_3prime = sample_profile.g_freq_3prime[0];

            agp::FrameSelector::finalize_sample_profile(sample_profile);

            auto pass1_end = std::chrono::high_resolution_clock::now();
            auto pass1_duration = std::chrono::duration_cast<std::chrono::milliseconds>(pass1_end - pass1_start);

            if (is_tty) std::cerr << "\r                                        \r";

            // Always show summary
            std::cerr << "  Reads: " << sample_profile.n_reads
                     << " | 5' damage: " << std::fixed << std::setprecision(1)
                     << sample_profile.max_damage_5prime * 100.0f << "%"
                     << " | 3' damage: " << sample_profile.max_damage_3prime * 100.0f << "%"
                     << " | ";
            // Classification based on average damage level
            float avg_damage = (sample_profile.max_damage_5prime + sample_profile.max_damage_3prime) / 2.0f;
            if (avg_damage > 0.10f) {
                std::cerr << "HIGH DAMAGE";      // >10% = clearly ancient
            } else if (avg_damage > 0.05f) {
                std::cerr << "MODERATE DAMAGE";  // 5-10% = likely ancient
            } else if (avg_damage > 0.02f) {
                std::cerr << "LOW DAMAGE";       // 2-5% = possibly ancient
            } else {
                std::cerr << "MINIMAL DAMAGE";   // <2% = likely modern or very degraded
            }
            std::cerr << " | " << std::setprecision(1) << pass1_duration.count() / 1000.0 << "s\n";

            // Verbose: show debug details
            if (opts.verbose) {
                std::cerr << "  [DEBUG] Raw counts at pos 0: T=" << std::fixed << std::setprecision(0)
                         << raw_t_5prime << " C=" << raw_c_5prime
                         << " A=" << raw_a_3prime << " G=" << raw_g_3prime << "\n";
                std::cerr << "  [DEBUG] 5' T/(T+C): " << std::setprecision(1)
                         << sample_profile.t_freq_5prime[0] * 100.0 << "%"
                         << " | 3' A/(A+G): " << sample_profile.a_freq_3prime[0] * 100.0 << "%"
                         << " | Baseline T/(T+C): "
                         << (sample_profile.baseline_t_freq / (sample_profile.baseline_t_freq + sample_profile.baseline_c_freq + 0.001)) * 100.0 << "%\n";
                std::cerr << "  [DEBUG] Decay λ: 5'=" << std::setprecision(2)
                         << sample_profile.lambda_5prime
                         << " 3'=" << sample_profile.lambda_3prime
                         << " (half-life: " << std::setprecision(1)
                         << 0.693f / sample_profile.lambda_5prime << " / "
                         << 0.693f / sample_profile.lambda_3prime << " bp)\n";
                std::cerr << "  [DEBUG] Library type: " << sample_profile.library_type_str() << "\n";
            }

            // Update damage model from Pass 1 profile
            damage_model.update_from_sample_profile(sample_profile);

            std::cerr << "\n";
        }

        // Main processing with double-buffering
        if (opts.aggregate_damage && opts.use_damage) {
            std::cerr << "[Pass 2] Gene prediction...\n";
        } else {
            std::cerr << "[Processing] Gene prediction...\n";
        }
        auto start_time = std::chrono::high_resolution_clock::now();
        agp::SequenceReader reader(opts.input_file);

        const size_t BATCH_SIZE = 50000;  // Larger batches for better thread utilization

        // Double-buffer: read next batch while processing current
        std::vector<agp::SequenceRecord> batch_a, batch_b;
        batch_a.reserve(BATCH_SIZE);
        batch_b.reserve(BATCH_SIZE);

        // Lambda to read a batch
        auto read_batch = [&reader](std::vector<agp::SequenceRecord>& batch) {
            batch.clear();
            agp::SequenceRecord record;
            while (batch.size() < BATCH_SIZE && reader.read_next(record)) {
                batch.push_back(std::move(record));
            }
            return batch.size();
        };

        size_t seq_count = 0;
        size_t total_genes = 0;

        // Damage probability distribution bins: [0-0.1), [0.1-0.2), ..., [0.9-1.0]
        std::array<size_t, 10> damage_prob_hist = {};
        // Damage percentage distribution bins: [0-10), [10-20), ..., [90-100]
        std::array<size_t, 10> damage_pct_hist = {};
        double total_damage_pct = 0.0;  // For computing mean

        // Correction statistics
        size_t total_genes_corrected = 0;
        size_t total_dna_corrections = 0;
        size_t total_aa_corrections = 0;

        std::mutex hist_mutex;

        // Read first batch synchronously
        read_batch(batch_a);
        if (batch_a.empty()) {
            std::cerr << "No sequences found in input file\n";
            return 0;
        }

        // Process with prefetching
        std::vector<agp::SequenceRecord>* current_batch = &batch_a;
        std::vector<agp::SequenceRecord>* next_batch = &batch_b;
        std::future<size_t> prefetch_future;
        bool prefetching = false;

        while (true) {
            // Start prefetching next batch asynchronously
            if (!prefetching) {
                prefetch_future = std::async(std::launch::async, read_batch, std::ref(*next_batch));
                prefetching = true;
            }

            auto& batch = *current_batch;
            if (batch.empty()) break;

            // Results storage
            std::vector<std::pair<std::string, std::vector<agp::Gene>>> results(batch.size());

            // Process batch in parallel
            #pragma omp parallel for schedule(static)
            for (size_t i = 0; i < batch.size(); ++i) {
                const auto& rec = batch[i];
                std::string seq = agp::SequenceUtils::clean(rec.sequence);

                if (seq.length() < opts.min_length) continue;

                // Use quality if available and lengths match (no cleaning removed bases)
                const std::string& quality = (rec.quality.length() == seq.length()) ? rec.quality : "";

                // Enable ensemble mode for metagenome scoring (weighted across all domains)
                if (opts.metagenome_mode && seq.length() >= 30) {
                    agp::MultiDomainResult domain_result = agp::score_all_domains(seq, 0);
                    agp::set_ensemble_mode(true);
                    agp::set_domain_probs(domain_result);
                } else {
                    agp::set_ensemble_mode(false);
                    agp::set_active_domain(active_domain);
                }

                // Create damage profile
                const agp::DamageProfile* profile_ptr = nullptr;
                agp::DamageProfile local_profile;
                if (opts.use_damage) {
                    local_profile = damage_model.create_profile(seq.length());
                    profile_ptr = &local_profile;
                }

                std::vector<agp::Gene> genes;

                if (opts.dual_strand) {
                    // Output best frame from BOTH strands
                    // Use appropriate frame selection method based on options
                    std::pair<agp::FrameScore, agp::FrameScore> strand_results;
                    if (opts.stop_only) {
                        // Stop-only frame selection: 97% accuracy, uses ONLY stop codon count
                        strand_results = agp::FrameSelector::select_best_per_strand_stop_only(seq);
                    } else {
                        // Standard scoring: strict stop penalties
                        // Use this regardless of damage detection settings
                        // (Damage-aware frame selection was found to hurt accuracy)
                        // Hexamer corrections are applied after frame selection
                        strand_results = agp::FrameSelector::select_best_per_strand(seq, profile_ptr);
                    }
                    auto [best_fwd, best_rev] = strand_results;

                    // Calculate ancient_prob using codon-aware detection for each strand
                    // This uses the reading frame information for better accuracy
                    float ancient_prob_fwd = 0.5f;
                    float ancient_prob_rev = 0.5f;

                    if (best_fwd.total_score >= opts.min_coding_prob) {
                        // Use quality-aware detection if quality available, else codon-aware
                        if (!quality.empty()) {
                            ancient_prob_fwd = agp::FrameSelector::estimate_ancient_prob_with_quality(
                                seq, quality, best_fwd.frame, sample_profile);
                        } else {
                            ancient_prob_fwd = agp::FrameSelector::estimate_ancient_prob_codon_aware(
                                seq, best_fwd.frame, sample_profile);
                        }
                    }

                    if (best_rev.total_score >= opts.min_coding_prob) {
                        // For reverse strand: damage is on ORIGINAL read, not RC
                        // Use original seq for damage scoring, but adjust frame for RC
                        // The frame in RC corresponds to (len - frame) % 3 in original
                        int orig_frame = (seq.length() - best_rev.frame) % 3;
                        if (!quality.empty()) {
                            ancient_prob_rev = agp::FrameSelector::estimate_ancient_prob_with_quality(
                                seq, quality, orig_frame, sample_profile);
                        } else {
                            ancient_prob_rev = agp::FrameSelector::estimate_ancient_prob_codon_aware(
                                seq, orig_frame, sample_profile);
                        }
                    }

                    // Determine which strands to output
                    bool output_fwd = (best_fwd.total_score >= opts.min_coding_prob);
                    bool output_rev = (best_rev.total_score >= opts.min_coding_prob);

                    // If --best-strand, only output the higher-scoring strand
                    if (opts.best_strand_only && output_fwd && output_rev) {
                        // Use combined score: coding_prob * (1 + ancient_prob) for ranking
                        float score_fwd = best_fwd.total_score * (1.0f + ancient_prob_fwd);
                        float score_rev = best_rev.total_score * (1.0f + ancient_prob_rev);
                        if (score_fwd >= score_rev) {
                            output_rev = false;
                        } else {
                            output_fwd = false;
                        }
                    }

                    // Compute damage percentage using sample-level distribution
                    // (damage is always on original read, not reverse complement)
                    float damage_pct = agp::FrameSelector::compute_damage_percentage(seq, sample_profile);

                    if (output_fwd) {
                        agp::Gene gene;
                        gene.start = 0;
                        gene.end = seq.length();
                        gene.is_forward = true;
                        gene.is_fragment = true;
                        gene.coding_prob = best_fwd.total_score;
                        gene.score = best_fwd.total_score;
                        gene.frame = best_fwd.frame;
                        gene.sequence = seq;
                        gene.protein = std::move(best_fwd.protein);
                        gene.ancient_prob = ancient_prob_fwd;
                        gene.damage_score = damage_pct;
                        genes.push_back(std::move(gene));
                    }

                    if (output_rev) {
                        agp::Gene gene;
                        gene.start = 0;
                        gene.end = seq.length();
                        gene.is_forward = false;
                        gene.is_fragment = true;
                        gene.coding_prob = best_rev.total_score;
                        gene.score = best_rev.total_score;
                        gene.frame = best_rev.frame;
                        // Store reverse complement as the coding sequence
                        gene.sequence = agp::FrameSelector::reverse_complement(seq);
                        gene.protein = std::move(best_rev.protein);
                        gene.ancient_prob = ancient_prob_rev;
                        gene.damage_score = damage_pct;
                        genes.push_back(std::move(gene));
                    }
                } else {
                    // Single best frame
                    agp::FrameScore best = agp::FrameSelector::select_best_frame(seq, profile_ptr);

                    if (best.total_score >= opts.min_coding_prob) {
                        agp::Gene gene;
                        gene.start = 0;
                        gene.end = seq.length();
                        gene.is_forward = best.forward;
                        gene.is_fragment = true;
                        gene.coding_prob = best.total_score;
                        gene.score = best.total_score;
                        gene.frame = best.frame;

                        // Use quality-aware damage detection if available, else codon-aware
                        // Store coding sequence (reverse complement for reverse strand)
                        std::string working_seq = best.forward ? seq :
                            agp::FrameSelector::reverse_complement(seq);
                        if (!quality.empty()) {
                            gene.ancient_prob = agp::FrameSelector::estimate_ancient_prob_with_quality(
                                working_seq, quality, best.frame, sample_profile);
                        } else {
                            gene.ancient_prob = agp::FrameSelector::estimate_ancient_prob_codon_aware(
                                working_seq, best.frame, sample_profile);
                        }

                        // Compute damage percentage using sample-level distribution
                        gene.damage_score = agp::FrameSelector::compute_damage_percentage(seq, sample_profile);

                        gene.sequence = std::move(working_seq);
                        gene.protein = std::move(best.protein);
                        genes.push_back(std::move(gene));
                    }
                }

                // Apply damage correction if confident
                // Use sample-profile-aware correction that scales by ancient_prob
                size_t batch_dna_corrections = 0;
                size_t batch_aa_corrections = 0;
                size_t batch_genes_corrected = 0;

                if (!genes.empty() && profile_ptr && opts.use_damage) {
                    for (auto& gene : genes) {
                        // Apply correction to reads with moderate-to-high damage probability
                        if (gene.ancient_prob > 0.6f) {
                            // Use new sample-profile-aware correction
                            // This uses position-specific damage rates from Phase 1
                            // and scales threshold by ancient_prob
                            auto [corrected_dna, dna_corr] = damage_model.correct_damage_with_profile(
                                gene.sequence, sample_profile, gene.ancient_prob, gene.frame);
                            gene.corrected_sequence = corrected_dna;
                            gene.dna_corrections = dna_corr;

                            // Re-translate corrected DNA to get corrected protein
                            if (dna_corr > 0) {
                                batch_dna_corrections += dna_corr;
                                batch_genes_corrected++;

                                gene.corrected_protein.clear();
                                for (size_t i = gene.frame; i + 2 < corrected_dna.length(); i += 3) {
                                    char c1 = corrected_dna[i];
                                    char c2 = corrected_dna[i + 1];
                                    char c3 = corrected_dna[i + 2];
                                    char aa = agp::CodonTable::translate_codon(c1, c2, c3);
                                    gene.corrected_protein += aa;
                                }

                                // Count amino acid changes
                                gene.aa_corrections = 0;
                                for (size_t i = 0; i < std::min(gene.protein.length(), gene.corrected_protein.length()); ++i) {
                                    if (gene.protein[i] != gene.corrected_protein[i]) {
                                        gene.aa_corrections++;
                                    }
                                }
                                batch_aa_corrections += gene.aa_corrections;
                            }
                        }
                    }
                }

                results[i] = {rec.id, std::move(genes)};
            }

            // Write results and track damage probability distribution
            size_t batch_total_dna_corr = 0;
            size_t batch_total_aa_corr = 0;
            size_t batch_total_genes_corr = 0;

            for (const auto& [id, genes] : results) {
                if (!genes.empty()) {
                    writer.write_genes(id, genes);
                    if (fasta_nt_writer) fasta_nt_writer->write_genes_nucleotide(id, genes);
                    if (fasta_aa_writer) fasta_aa_writer->write_genes_protein(id, genes);

                    total_genes += genes.size();

                    // Update damage probability and percentage histograms
                    for (const auto& gene : genes) {
                        int bin = std::min(9, static_cast<int>(gene.ancient_prob * 10));
                        damage_prob_hist[bin]++;

                        // Damage percentage: 0-100, bin into 10% increments
                        int pct_bin = std::min(9, static_cast<int>(gene.damage_score / 10.0));
                        damage_pct_hist[pct_bin]++;
                        total_damage_pct += gene.damage_score;

                        // Track corrections
                        if (gene.dna_corrections > 0) {
                            batch_total_dna_corr += gene.dna_corrections;
                            batch_total_aa_corr += gene.aa_corrections;
                            batch_total_genes_corr++;
                        }
                    }
                }
            }

            // Update global correction stats
            total_dna_corrections += batch_total_dna_corr;
            total_aa_corrections += batch_total_aa_corr;
            total_genes_corrected += batch_total_genes_corr;

            seq_count += batch.size();
            if (seq_count % 1000000 < batch.size()) {
                if (is_tty) {
                    std::cerr << "\r  Processed " << seq_count / 1000000 << "M sequences..." << std::flush;
                } else {
                    std::cerr << "  Processed " << seq_count / 1000000 << "M sequences..." << std::endl;
                }
            }

            // Wait for prefetch to complete and swap buffers
            if (prefetching) {
                prefetch_future.wait();
                prefetching = false;
            }
            std::swap(current_batch, next_batch);
        }

        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);

        if (is_tty) std::cerr << "\r                                                  \r";  // Clear progress line

        // Always show summary
        double speed = (duration.count() > 0) ? (seq_count * 1000.0 / duration.count()) : 0;
        std::cerr << "  Sequences: " << seq_count
                 << " | Genes: " << total_genes
                 << " | " << std::fixed << std::setprecision(1) << duration.count() / 1000.0 << "s"
                 << " (" << std::setprecision(0) << speed << " seq/s)\n";

        // Show damage distributions with ASCII bar charts
        if (opts.use_damage && total_genes > 0) {
            // Helper lambda to draw histogram bar
            auto draw_bar = [](double pct, int max_width = 20) -> std::string {
                int filled = static_cast<int>(pct / 100.0 * max_width + 0.5);
                filled = std::min(filled, max_width);
                return std::string(filled, '#') + std::string(max_width - filled, ' ');
            };

            // Ancient probability distribution
            std::cerr << "\n  Ancient Probability Distribution:\n";
            const char* prob_labels[] = {"0.0-0.1", "0.1-0.2", "0.2-0.3", "0.3-0.4", "0.4-0.5",
                                         "0.5-0.6", "0.6-0.7", "0.7-0.8", "0.8-0.9", "0.9-1.0"};
            for (int i = 0; i < 10; ++i) {
                double pct = 100.0 * damage_prob_hist[i] / total_genes;
                std::cerr << "    " << prob_labels[i] << " |" << draw_bar(pct)
                         << "| " << std::fixed << std::setprecision(1) << std::setw(5) << pct << "%\n";
            }

            // Damage percentage distribution
            double mean_dmg_pct = total_damage_pct / total_genes;
            std::cerr << "\n  Per-Read Damage Score (mean: " << std::setprecision(1) << mean_dmg_pct << "%):\n";
            const char* pct_labels[] = {"  0-10%", " 10-20%", " 20-30%", " 30-40%", " 40-50%",
                                        " 50-60%", " 60-70%", " 70-80%", " 80-90%", "90-100%"};
            for (int i = 0; i < 10; ++i) {
                double pct = 100.0 * damage_pct_hist[i] / total_genes;
                std::cerr << "    " << pct_labels[i] << " |" << draw_bar(pct)
                         << "| " << std::fixed << std::setprecision(1) << std::setw(5) << pct << "%\n";
            }
        }

        if (opts.verbose) {
            std::cerr << "\n  [DEBUG] Gene/seq ratio: " << std::setprecision(3)
                     << (seq_count > 0 ? (double)total_genes / seq_count : 0.0) << "\n";

            // Correction statistics
            if (total_genes_corrected > 0) {
                double corr_pct = 100.0 * total_genes_corrected / total_genes;
                double avg_dna_corr = (double)total_dna_corrections / total_genes_corrected;
                double avg_aa_corr = (double)total_aa_corrections / total_genes_corrected;
                std::cerr << "\n  [CORRECTION] Genes corrected: " << total_genes_corrected
                         << " / " << total_genes << " (" << std::setprecision(1) << corr_pct << "%)\n";
                std::cerr << "  [CORRECTION] Total DNA corrections: " << total_dna_corrections
                         << " (avg " << std::setprecision(2) << avg_dna_corr << " per gene)\n";
                std::cerr << "  [CORRECTION] Total AA corrections: " << total_aa_corrections
                         << " (avg " << std::setprecision(2) << avg_aa_corr << " per gene)\n";
            } else {
                std::cerr << "\n  [CORRECTION] No genes corrected (0 / " << total_genes << ")\n";
            }
        }

        std::cerr << "\nDone.\n";

        writer.close();
        return 0;

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
}
