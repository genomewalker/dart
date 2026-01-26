#include "cli/args.hpp"
#include "agp/damage_model.hpp"
#include "agp/adaptive_damage.hpp"
#include "agp/sequence_io.hpp"
#include "agp/frame_selector.hpp"
#include "agp/hexamer_tables.hpp"
#include "agp/hexamer_fast.hpp"
#include "agp/codon_tables.hpp"
#include "agp/bdamage_reader.hpp"
#include "agp/score_calibration_2d.hpp"
#include "agp/score_calibration_4d.hpp"
#include "agp/decoy_generator.hpp"
#include "agp/fdr_estimator.hpp"
#include "agp/bayesian_frame_v2.hpp"
#include "agp/version.h"
#include <iostream>
#include <fstream>
#include <chrono>
#include <iomanip>
#include <mutex>
#include <atomic>
#include <vector>
#include <array>
#include <future>
#include <unistd.h>
#include <unordered_map>
#include <algorithm>
#include <cstring>

#ifdef _OPENMP
#include <omp.h>
#endif

// Use CLI types from the cli module
using Options = agp::cli::Options;

// Convert CLI library type to SampleDamageProfile library type
inline agp::SampleDamageProfile::LibraryType to_sample_library_type(agp::cli::LibraryType lt) {
    switch (lt) {
        case agp::cli::LibraryType::DOUBLE_STRANDED:
            return agp::SampleDamageProfile::LibraryType::DOUBLE_STRANDED;
        case agp::cli::LibraryType::SINGLE_STRANDED:
            return agp::SampleDamageProfile::LibraryType::SINGLE_STRANDED;
        default:
            return agp::SampleDamageProfile::LibraryType::UNKNOWN;
    }
}

int main(int argc, char* argv[]) {
    try {
        Options opts = agp::cli::parse_args(argc, argv);

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
        agp::Domain base_domain = active_domain;

        // Per-domain damage models for meta blending
        std::array<agp::DamageModel, 8> domain_damage_models;

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
            if (opts.iterative_damage) {
                std::cerr << "Iterative damage: ON (4-pass refinement)\n";
                std::cerr << "Iterative threshold: " << opts.iterative_threshold << "\n";
            }
            if (opts.all_frames) std::cerr << "Output mode: ALL-FRAMES (6 frames per read)\n";
            if (opts.use_bayesian) std::cerr << "Frame selection: BAYESIAN (hexamer posteriors)\n";
            if (opts.dual_strand) {
                std::cerr << "Dual-strand: ON";
                if (opts.orf_confidence_threshold >= 0.0f) {
                    std::cerr << " (adaptive: threshold=" << opts.orf_confidence_threshold << ")";
                } else if (opts.best_strand_only) {
                    std::cerr << " (best only)";
                }
                std::cerr << "\n";
            }
            if (opts.forced_library_type != agp::cli::LibraryType::UNKNOWN) {
                std::cerr << "Library type: FORCED "
                         << agp::cli::library_type_to_string(opts.forced_library_type) << "\n";
            }
            if (opts.compute_fdr) {
                std::cerr << "FDR estimation: ON (decoys=" << opts.fdr_decoys
                          << ", strategy=" << opts.fdr_strategy
                          << ", seed=" << opts.fdr_seed << ")\n";
                if (opts.fdr_threshold >= 0.0f) {
                    std::cerr << "FDR threshold: " << (opts.fdr_threshold * 100.0f) << "%\n";
                }
            }
        }
        std::cerr << "\n";

        // Initialize decoy generator for FDR estimation
        std::unique_ptr<agp::DecoyGenerator> decoy_gen;
        std::vector<float> fdr_target_scores;
        std::vector<float> fdr_decoy_scores;
        std::mutex fdr_mutex;
        agp::DecoyStrategy fdr_decoy_strategy = agp::DecoyStrategy::DINUCLEOTIDE_SHUFFLE;
        if (opts.compute_fdr) {
            decoy_gen = std::make_unique<agp::DecoyGenerator>(opts.fdr_seed);
            fdr_target_scores.reserve(100000);  // Pre-allocate for efficiency
            fdr_decoy_scores.reserve(100000 * opts.fdr_decoys);

            // Map strategy string to enum
            if (opts.fdr_strategy == "synonymous") {
                fdr_decoy_strategy = agp::DecoyStrategy::SYNONYMOUS_SHUFFLE;
            } else if (opts.fdr_strategy == "dinucleotide") {
                fdr_decoy_strategy = agp::DecoyStrategy::DINUCLEOTIDE_SHUFFLE;
            } else if (opts.fdr_strategy == "random") {
                fdr_decoy_strategy = agp::DecoyStrategy::RANDOM_FRAME;
            }
        }

        // Initialize damage model
        agp::DamageModel damage_model;

        // Open output files
        agp::GeneWriter writer(opts.output_file);
        std::unique_ptr<agp::FastaWriter> fasta_nt_writer;
        std::unique_ptr<agp::FastaWriter> fasta_nt_corr_writer;
        std::unique_ptr<agp::FastaWriter> fasta_aa_writer;
        std::unique_ptr<agp::FastaWriter> fasta_corrected_writer;

        if (!opts.fasta_nt.empty()) {
            fasta_nt_writer = std::make_unique<agp::FastaWriter>(opts.fasta_nt);
        }
        if (!opts.fasta_nt_corrected.empty()) {
            fasta_nt_corr_writer = std::make_unique<agp::FastaWriter>(opts.fasta_nt_corrected);
        }
        if (!opts.fasta_aa.empty()) {
            fasta_aa_writer = std::make_unique<agp::FastaWriter>(opts.fasta_aa);
        }
        if (!opts.fasta_corrected.empty()) {
            fasta_corrected_writer = std::make_unique<agp::FastaWriter>(opts.fasta_corrected);
        }

        // Sample-level damage profile
        agp::SampleDamageProfile sample_profile;
        sample_profile.forced_library_type = to_sample_library_type(opts.forced_library_type);
        // Per-domain damage profiles for meta mode blending
        std::array<agp::SampleDamageProfile, 8> domain_profiles;
        for (auto& dp : domain_profiles) {
            dp.forced_library_type = to_sample_library_type(opts.forced_library_type);
        }

        // Adaptive damage calibrator (initialized after Pass 1/2)
        agp::AdaptiveDamageCalibrator calibrator;

        // Flag to track if damage profile was loaded externally
        bool bdamage_loaded = false;

        // Load external bdamage file if provided (skips Pass 1)
        if (!opts.bdamage_file.empty() && opts.use_damage) {
            std::cerr << "[External] Loading damage profile from metaDMG bdamage file...\n";
            std::cerr << "  File: " << opts.bdamage_file << "\n";

            auto bdamage_result = agp::load_bdamage(opts.bdamage_file);
            if (!bdamage_result.success) {
                std::cerr << "  [ERROR] Failed to load bdamage file: " << bdamage_result.error_message << "\n";
                return 1;
            }

            sample_profile = bdamage_result.profile;
            sample_profile.forced_library_type = to_sample_library_type(opts.forced_library_type);
            bdamage_loaded = true;

            std::cerr << "  Loaded damage profile from " << bdamage_result.num_references << " reference(s)\n";
            std::cerr << "  5' damage (C→T): " << std::fixed << std::setprecision(1)
                      << (sample_profile.max_damage_5prime * 100.0f) << "%\n";
            std::cerr << "  3' damage (G→A): " << std::fixed << std::setprecision(1)
                      << (sample_profile.max_damage_3prime * 100.0f) << "%\n";

            // Determine damage level
            float d_max = std::max(sample_profile.max_damage_5prime, sample_profile.max_damage_3prime);
            std::string level;
            if (d_max >= 0.10f) level = "high";
            else if (d_max >= 0.05f) level = "moderate";
            else if (d_max >= 0.02f) level = "low";
            else level = "undetectable";
            std::cerr << "  Damage level: " << level << " (D_max=" << std::setprecision(1)
                      << (d_max * 100.0f) << "%)\n";

            // Initialize calibrator with external profile
            calibrator.initialize(sample_profile);

            // Update damage model with external profile
            damage_model.update_from_sample_profile(sample_profile);
        }

                // First pass: collect damage statistics (always runs; blends with external if provided)
        if (opts.aggregate_damage && opts.use_damage) {
            auto pass1_start = std::chrono::high_resolution_clock::now();
            std::cerr << "[Pass 1] Scanning for damage patterns...\n";

            // Use a temporary reference-free profile when external damage is present
            agp::SampleDamageProfile ref_profile;
            ref_profile.forced_library_type = to_sample_library_type(opts.forced_library_type);
            agp::SampleDamageProfile& pass_profile = bdamage_loaded ? ref_profile : sample_profile;

            agp::SequenceReader first_pass(opts.input_file);
            const size_t PASS1_BATCH_SIZE = 50000;  // Larger batches for Pass 1
            std::vector<agp::SequenceRecord> batch;
            batch.reserve(PASS1_BATCH_SIZE);
            size_t count = 0;

            // Thread-local profiles for parallel aggregation
            std::vector<agp::SampleDamageProfile> thread_profiles(num_threads);
            for (auto& tp : thread_profiles) {
                tp.forced_library_type = to_sample_library_type(opts.forced_library_type);
            }
            // Per-domain profiles for meta blending (only if not using external)
            std::vector<std::array<agp::SampleDamageProfile, 8>> thread_domain_profiles;
            if (opts.metagenome_mode && !bdamage_loaded) {
                thread_domain_profiles.resize(num_threads);
                for (int t = 0; t < num_threads; ++t) {
                    for (auto& dp : thread_domain_profiles[t]) {
                        dp.forced_library_type = to_sample_library_type(opts.forced_library_type);
                    }
                }
            }

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

                        if (!bdamage_loaded && opts.metagenome_mode && seq.length() >= 30) {
                            agp::MultiDomainResult dom = agp::score_all_domains_fast(seq, 0);
                            const float probs[8] = {
                                dom.gtdb_prob, dom.fungi_prob, dom.protozoa_prob, dom.invertebrate_prob,
                                dom.plant_prob, dom.vertebrate_mammalian_prob, dom.vertebrate_other_prob,
                                dom.viral_prob
                            };
                            for (int d = 0; d < 8; ++d) {
                                if (probs[d] > 0.01f) {
                                    agp::FrameSelector::update_sample_profile_weighted(
                                        thread_domain_profiles[tid][d], seq, probs[d]);
                                }
                            }
                        }
                    }
                }

                count += batch.size();
                if (count % 1000000 < PASS1_BATCH_SIZE) {
                    if (is_tty) {
                        std::cerr << "  Scanned " << count / 1000000 << "M sequences..." << std::flush;
                    } else {
                        std::cerr << "  Scanned " << count / 1000000 << "M sequences..." << std::endl;
                    }
                }
            }

            // Merge all thread-local profiles
            for (int t = 0; t < num_threads; ++t) {
                agp::FrameSelector::merge_sample_profiles(pass_profile, thread_profiles[t]);
                if (!bdamage_loaded && opts.metagenome_mode) {
                    for (int d = 0; d < 8; ++d) {
                        agp::FrameSelector::merge_sample_profiles(domain_profiles[d], thread_domain_profiles[t][d]);
                    }
                }
            }

            // Store raw counts before finalization for debug
            double raw_t_5prime = pass_profile.t_freq_5prime[0];
            double raw_c_5prime = pass_profile.c_freq_5prime[0];
            double raw_a_3prime = pass_profile.a_freq_3prime[0];
            double raw_g_3prime = pass_profile.g_freq_3prime[0];

            agp::FrameSelector::finalize_sample_profile(pass_profile);

            auto pass1_end = std::chrono::high_resolution_clock::now();
            auto pass1_duration = std::chrono::duration_cast<std::chrono::milliseconds>(pass1_end - pass1_start);

            if (is_tty) std::cerr << "                                        ";

            if (bdamage_loaded) {
                // Blend external (sample_profile) with reference-free (pass_profile)
                // When reference-free shows inverted pattern or very low damage but external shows high,
                // trust external more heavily
                std::array<float, 15> ext5{}, ext3{};
                for (int i = 0; i < 15; ++i) {
                    ext5[i] = sample_profile.damage_rate_5prime[i];
                    ext3[i] = sample_profile.damage_rate_3prime[i];
                }
                float ext_lambda5 = sample_profile.lambda_5prime;
                float ext_lambda3 = sample_profile.lambda_3prime;
                float ext_max5 = sample_profile.max_damage_5prime;
                float ext_max3 = sample_profile.max_damage_3prime;

                // Calculate blending weights
                // Base weight for reference-free: up to 0.3 based on read count
                float base_w_ref = std::min(0.3f, static_cast<float>(pass_profile.n_reads) / 1000000.0f);

                // Reduce reference-free weight when:
                // 1. Reference-free shows inverted pattern (unreliable)
                // 2. Reference-free shows much lower damage than external
                bool ref_free_unreliable = pass_profile.terminal_inversion ||
                    pass_profile.inverted_pattern_5prime || pass_profile.inverted_pattern_3prime;
                float ref_max_damage = std::max(pass_profile.max_damage_5prime, pass_profile.max_damage_3prime);
                float ext_max_damage = std::max(ext_max5, ext_max3);
                bool damage_mismatch = (ext_max_damage > 0.10f) && (ref_max_damage < ext_max_damage * 0.3f);

                float w_ref = base_w_ref;
                if (ref_free_unreliable || damage_mismatch) {
                    w_ref = 0.0f;  // Trust external fully when reference-free is unreliable
                }
                float w_ext = 1.0f - w_ref;

                for (int i = 0; i < 15; ++i) {
                    sample_profile.damage_rate_5prime[i] = w_ext * ext5[i] + w_ref * pass_profile.damage_rate_5prime[i];
                    sample_profile.damage_rate_3prime[i] = w_ext * ext3[i] + w_ref * pass_profile.damage_rate_3prime[i];
                    sample_profile.t_freq_5prime[i] = pass_profile.t_freq_5prime[i];
                    sample_profile.c_freq_5prime[i] = pass_profile.c_freq_5prime[i];
                    sample_profile.a_freq_3prime[i] = pass_profile.a_freq_3prime[i];
                    sample_profile.g_freq_3prime[i] = pass_profile.g_freq_3prime[i];
                }

                sample_profile.max_damage_5prime = w_ext * ext_max5 + w_ref * pass_profile.max_damage_5prime;
                sample_profile.max_damage_3prime = w_ext * ext_max3 + w_ref * pass_profile.max_damage_3prime;
                sample_profile.lambda_5prime = w_ext * ext_lambda5 + w_ref * pass_profile.lambda_5prime;
                sample_profile.lambda_3prime = w_ext * ext_lambda3 + w_ref * pass_profile.lambda_3prime;

                if (pass_profile.n_reads > 0) {
                    sample_profile.baseline_t_freq = pass_profile.baseline_t_freq;
                    sample_profile.baseline_c_freq = pass_profile.baseline_c_freq;
                    sample_profile.baseline_a_freq = pass_profile.baseline_a_freq;
                    sample_profile.baseline_g_freq = pass_profile.baseline_g_freq;
                }

                sample_profile.n_reads = 10000 + pass_profile.n_reads;
                sample_profile.terminal_shift_5prime = pass_profile.terminal_shift_5prime;
                sample_profile.terminal_shift_3prime = pass_profile.terminal_shift_3prime;
                sample_profile.terminal_z_5prime = pass_profile.terminal_z_5prime;
                sample_profile.terminal_z_3prime = pass_profile.terminal_z_3prime;
                // Don't copy terminal_inversion from pass_profile - trust external for this

                damage_model.update_from_sample_profile(sample_profile);

                if (opts.verbose) {
                    std::cerr << "  [Blend] External bdamage + ref-free: w_ext=" << std::fixed << std::setprecision(2)
                             << w_ext << " w_ref=" << w_ref;
                    if (ref_free_unreliable) std::cerr << " (ref-free unreliable)";
                    if (damage_mismatch) std::cerr << " (damage mismatch)";
                    std::cerr << "\n";
                }
            } else {
                // Always show summary
                std::cerr << "  Reads: " << sample_profile.n_reads
                         << " | 5' damage: " << std::fixed << std::setprecision(1)
                         << sample_profile.max_damage_5prime * 100.0f << "%"
                         << " | 3' damage: " << sample_profile.max_damage_3prime * 100.0f << "%"
                         << " | ";
                // Classification based on average damage level
                float avg_damage = (sample_profile.max_damage_5prime + sample_profile.max_damage_3prime) / 2.0f;
                if (avg_damage > 0.10f) {
                    std::cerr << "HIGH DAMAGE";
                } else if (avg_damage > 0.05f) {
                    std::cerr << "MODERATE DAMAGE";
                } else if (avg_damage > 0.02f) {
                    std::cerr << "LOW DAMAGE";
                } else {
                    std::cerr << "MINIMAL DAMAGE";
                }
                std::cerr << " | " << std::setprecision(1) << pass1_duration.count() / 1000.0 << "s\n";

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

                damage_model.update_from_sample_profile(sample_profile);

                std::cerr << "\n";
            }

            calibrator.initialize(sample_profile);

            if (bdamage_loaded && opts.verbose) {
                std::cerr << "  [INFO] Blended external bdamage with reference-free estimate.\n\n";
            }
        }
// ========================================================================
        // ITERATIVE DAMAGE REFINEMENT (Pass 2 + Pass 3)
        // When enabled, performs:
        //   Pass 2: Quick gene scoring (coding probability only)
        //   Pass 3: Re-estimate damage using only high-scoring sequences
        // This helps extract damage signal from metagenomes with mixed content
        // ========================================================================
        agp::SampleDamageProfile refined_profile;  // Will hold refined damage profile
        refined_profile.forced_library_type = to_sample_library_type(opts.forced_library_type);
        if (opts.iterative_damage && opts.aggregate_damage && opts.use_damage) {
            // Store initial profile for comparison
            float initial_5prime = sample_profile.max_damage_5prime;
            float initial_3prime = sample_profile.max_damage_3prime;

            std::cerr << "[Pass 2] Codon-damage scoring for damage enrichment...\n";
            auto pass2_start = std::chrono::high_resolution_clock::now();

            // Re-read file for scoring pass
            agp::SequenceReader scoring_pass(opts.input_file);
            const size_t PASS2_BATCH_SIZE = 50000;
            std::vector<agp::SequenceRecord> batch;
            batch.reserve(PASS2_BATCH_SIZE);

            // Thread-local profiles for refined damage estimation
            std::vector<agp::SampleDamageProfile> thread_profiles(num_threads);
            size_t high_score_count = 0;
            size_t pass2_count = 0;

            while (true) {
                batch.clear();
                agp::SequenceRecord record;
                while (batch.size() < PASS2_BATCH_SIZE && scoring_pass.read_next(record)) {
                    batch.push_back(std::move(record));
                }
                if (batch.empty()) break;

                // Thread-local counters for high-scoring sequences
                std::vector<size_t> thread_high_count(num_threads, 0);

                #pragma omp parallel for schedule(static)
                for (size_t i = 0; i < batch.size(); ++i) {
                    int tid = 0;
                    #ifdef _OPENMP
                    tid = omp_get_thread_num();
                    #endif

                    std::string seq = agp::SequenceUtils::clean(batch[i].sequence);
                    if (seq.length() < opts.min_length) continue;

                    // Use damage codon score to identify reads with damage patterns
                    // This scores reads by damage-diagnostic codon patterns:
                    // - Terminal T enrichment (C→T damage)
                    // - Damage-induced stop codons (TAA/TAG/TGA from CAA/CAG/CGA)
                    // - Within-read T/(T+C) gradient
                    float damage_codon_score = agp::FrameSelector::compute_damage_codon_score(seq);

                    // Include sequences with damage signal above threshold
                    // Higher threshold = more stringent filtering for damage-indicative reads
                    if (damage_codon_score >= opts.iterative_threshold) {
                        // Weight by damage codon score: higher scores = more damage signal
                        float weight = damage_codon_score;
                        agp::FrameSelector::update_sample_profile_weighted(
                            thread_profiles[tid], seq, weight);
                        thread_high_count[tid]++;
                    }
                }

                // Sum up high-score counts
                for (int t = 0; t < num_threads; ++t) {
                    high_score_count += thread_high_count[t];
                    thread_high_count[t] = 0;  // Reset for next batch
                }

                pass2_count += batch.size();
                if (pass2_count % 1000000 < PASS2_BATCH_SIZE) {
                    if (is_tty) {
                        std::cerr << "\r  Scanned " << pass2_count / 1000000 << "M sequences ("
                                 << high_score_count / 1000 << "K damage-indicative)..." << std::flush;
                    } else {
                        std::cerr << "  Scanned " << pass2_count / 1000000 << "M sequences ("
                                 << high_score_count / 1000 << "K damage-indicative)..." << std::endl;
                    }
                }
            }

            // Merge thread-local profiles into refined profile
            for (int t = 0; t < num_threads; ++t) {
                agp::FrameSelector::merge_sample_profiles(refined_profile, thread_profiles[t]);
            }

            // Finalize the refined profile
            agp::FrameSelector::finalize_sample_profile(refined_profile);

            auto pass2_end = std::chrono::high_resolution_clock::now();
            auto pass2_duration = std::chrono::duration_cast<std::chrono::milliseconds>(pass2_end - pass2_start);

            if (is_tty) std::cerr << "\r                                                            \r";

            // Report refinement results
            float refined_avg = (refined_profile.max_damage_5prime + refined_profile.max_damage_3prime) / 2.0f;
            float initial_avg = (initial_5prime + initial_3prime) / 2.0f;

            std::cerr << "  Damage-indicative reads: " << high_score_count << " / " << pass2_count
                     << " (" << std::fixed << std::setprecision(1)
                     << (100.0 * high_score_count / (pass2_count + 1)) << "%)\n";
            std::cerr << "  Refined damage: 5'=" << std::setprecision(1)
                     << refined_profile.max_damage_5prime * 100.0f << "% 3'="
                     << refined_profile.max_damage_3prime * 100.0f << "%";

            // Show improvement or note if refinement helped
            if (refined_avg > initial_avg + 0.01f) {
                std::cerr << " (+" << std::setprecision(1)
                         << (refined_avg - initial_avg) * 100.0f << "% vs all-reads)\n";
            } else if (refined_avg < initial_avg - 0.01f) {
                std::cerr << " (" << std::setprecision(1)
                         << (refined_avg - initial_avg) * 100.0f << "% vs all-reads)\n";
            } else {
                std::cerr << " (similar to all-reads)\n";
            }

            std::cerr << "  | " << std::setprecision(1) << pass2_duration.count() / 1000.0 << "s\n";

            // Determine if we should trust the refined profile
            // The codon-damage enrichment can produce false positives on modern DNA
            // by selecting naturally AT-biased sequences. We detect this by checking
            // whether the enrichment produced a physically plausible damage signal.

            float initial_t_ratio = 0.5f;
            float baseline_t_ratio = 0.5f;
            if (sample_profile.t_freq_5prime[0] + sample_profile.c_freq_5prime[0] > 0) {
                initial_t_ratio = sample_profile.t_freq_5prime[0] /
                    (sample_profile.t_freq_5prime[0] + sample_profile.c_freq_5prime[0]);
            }
            if (sample_profile.baseline_t_freq + sample_profile.baseline_c_freq > 0) {
                baseline_t_ratio = sample_profile.baseline_t_freq /
                    (sample_profile.baseline_t_freq + sample_profile.baseline_c_freq);
            }

            // Check if refined profile shows DECAY pattern (hallmark of real damage)
            // Real damage: high at position 0, decreasing toward interior
            // False positive: high at ALL positions (just AT-biased reads)
            float refined_pos0_t = 0.5f;
            float refined_pos4_t = 0.5f;  // Position 4-5 for decay check
            bool decay_data_available = false;
            if (refined_profile.t_freq_5prime[0] + refined_profile.c_freq_5prime[0] > 0) {
                refined_pos0_t = refined_profile.t_freq_5prime[0] /
                    (refined_profile.t_freq_5prime[0] + refined_profile.c_freq_5prime[0]);
            }
            if (refined_profile.t_freq_5prime[4] + refined_profile.c_freq_5prime[4] > 0) {
                refined_pos4_t = refined_profile.t_freq_5prime[4] /
                    (refined_profile.t_freq_5prime[4] + refined_profile.c_freq_5prime[4]);
                decay_data_available = true;
            }

            // Decay signal: position 0 should be notably higher than position 4
            // If no decay, we're just selecting AT-biased sequences
            float decay_signal = refined_pos0_t - refined_pos4_t;
            bool has_decay_pattern = decay_signal > 0.03f;  // At least 3% decay over 4 positions

            // Check decay at position 10 (interior) for stronger validation
            float refined_pos10_t = 0.5f;
            bool pos10_valid = false;
            if (refined_profile.t_freq_5prime[10] + refined_profile.c_freq_5prime[10] > 0) {
                refined_pos10_t = refined_profile.t_freq_5prime[10] /
                    (refined_profile.t_freq_5prime[10] + refined_profile.c_freq_5prime[10]);
                pos10_valid = true;
            }
            float decay_to_interior = refined_pos0_t - refined_pos10_t;
            bool has_strong_decay = pos10_valid && decay_to_interior > 0.05f;  // >5% decay to interior

            // Check R² of refined profile - good fit indicates real exponential decay
            bool has_good_fit = refined_profile.r_squared_5prime > 0.7f ||
                                refined_profile.r_squared_3prime > 0.7f;

            // Check if 3' also shows decay (for ds libraries, both ends should decay)
            float refined_pos0_a = 0.5f;
            float refined_pos4_a = 0.5f;
            if (refined_profile.a_freq_3prime[0] + refined_profile.g_freq_3prime[0] > 0) {
                refined_pos0_a = refined_profile.a_freq_3prime[0] /
                    (refined_profile.a_freq_3prime[0] + refined_profile.g_freq_3prime[0]);
            }
            if (refined_profile.a_freq_3prime[4] + refined_profile.g_freq_3prime[4] > 0) {
                refined_pos4_a = refined_profile.a_freq_3prime[4] /
                    (refined_profile.a_freq_3prime[4] + refined_profile.g_freq_3prime[4]);
            }
            float decay_3prime = refined_pos0_a - refined_pos4_a;
            bool has_3prime_decay = decay_3prime > 0.02f;

            // Initial profile must show SOME damage signal for refinement to be trusted
            // Use computed d_max values (which use correct C->T formula) rather than raw T/(T+C)
            bool initial_has_damage_signal = (initial_t_ratio > baseline_t_ratio + 0.01f) ||
                                              (initial_5prime > 0.005f) || (initial_3prime > 0.005f);

            // CRITICAL FIX: Check if initial d_max is essentially zero
            // If raw estimate shows 0% damage, iterative mode is likely selecting AT-rich reads
            // that happen to have T at terminal positions (compositional bias, not damage)
            bool initial_is_zero = (initial_avg < 0.005f);  // <0.5% raw damage = essentially zero

            // Calculate enrichment ratio: how much did iterative boost the signal?
            float enrichment_ratio = (initial_avg > 0.001f) ? (refined_avg / initial_avg) : 999.0f;

            // For truly damaged samples, iterative should boost by 2-10x (diluted damage)
            // For modern samples, iterative can boost by 30x+ (selecting compositional bias)
            bool suspicious_enrichment = (enrichment_ratio > 20.0f) && (initial_avg < 0.02f);

            // Use refined profile if:
            // 1. Sufficient reads and signal, AND
            // 2. Either initial profile had signal, OR refined shows decay pattern
            // 3. NOT a suspicious case (zero initial + massive enrichment)
            bool use_refined = refined_profile.n_reads >= 1000 && refined_avg > 0.01f;

            // Smarter decision logic for the middle case:
            // - dfb2272499: initial_avg=54.7%, enrichment ~1.2x → TRUST RAW (high R², strong signal)
            // - 521588e724: initial_avg=1.97%, enrichment ~15x → TRUST ITERATIVE (real diluted damage)
            // - 37e4738093: initial_avg=0.0%, enrichment=∞ → TRUST RAW (no signal = modern)
            // - 5fa278026b: initial_avg=0.0%, but refined shows decay + R² → TRUST ITERATIVE

            if (initial_is_zero) {
                // Initial profile showed essentially NO damage signal
                // BUT: if refined profile shows clear decay pattern AND good R² fit,
                // this indicates the damage was masked by metagenomic composition
                // and iterative refinement successfully recovered the signal.

                // Criteria for trusting refined when initial is zero:
                // 1. Clear decay pattern at 5' (pos0 - pos4 > 3%)
                // 2. Strong decay to interior (pos0 - pos10 > 5%)
                // 3. Good exponential fit (R² > 0.7) OR symmetric 3' decay (for ds)
                // 4. Refined d_max in plausible range (5-70%)
                bool refined_is_plausible = (refined_avg > 0.05f && refined_avg < 0.70f);
                bool has_validated_decay = has_decay_pattern &&
                                           (has_strong_decay || !pos10_valid) &&
                                           (has_good_fit || has_3prime_decay);

                if (refined_is_plausible && has_validated_decay) {
                    // Refined profile passes validation - trust it
                    use_refined = true;
                    if (opts.verbose) {
                        std::cerr << "  [INFO] Initial d_max is near-zero (" << std::fixed << std::setprecision(2)
                                  << initial_avg * 100.0f << "%), but refined profile passes decay validation\n";
                        std::cerr << "  [INFO] Refined decay: pos0=" << std::setprecision(1) << refined_pos0_t * 100.0f
                                  << "% pos4=" << refined_pos4_t * 100.0f
                                  << "% pos10=" << refined_pos10_t * 100.0f << "%\n";
                        std::cerr << "  [INFO] R²: 5'=" << std::setprecision(3) << refined_profile.r_squared_5prime
                                  << " 3'=" << refined_profile.r_squared_3prime << "\n";
                        std::cerr << "  [INFO] Accepting refined profile (damage likely masked by composition)\n";
                    }
                } else {
                    // Refined profile fails validation - reject
                    use_refined = false;
                    if (opts.verbose) {
                        std::cerr << "  [INFO] Initial d_max is near-zero (" << std::fixed << std::setprecision(2)
                                  << initial_avg * 100.0f << "%), likely modern sample\n";
                        if (!refined_is_plausible) {
                            std::cerr << "  [INFO] Refined d_max=" << refined_avg * 100.0f
                                      << "% outside plausible range (5-70%)\n";
                        }
                        if (!has_validated_decay) {
                            std::cerr << "  [INFO] Refined profile lacks validated decay pattern\n";
                            std::cerr << "  [INFO] decay_pattern=" << (has_decay_pattern ? "yes" : "no")
                                      << " strong_decay=" << (has_strong_decay ? "yes" : "no")
                                      << " good_fit=" << (has_good_fit ? "yes" : "no")
                                      << " 3prime_decay=" << (has_3prime_decay ? "yes" : "no") << "\n";
                        }
                        std::cerr << "  [INFO] Rejecting refined profile (enrichment=" << std::setprecision(0)
                                  << enrichment_ratio << "x likely from compositional bias)\n";
                    }
                }
            } else if (!initial_has_damage_signal && !has_decay_pattern) {
                // Sample showed weak initial damage AND refined has no decay pattern
                // This indicates we're just selecting AT-biased reads, not damaged ones
                if (refined_avg > 0.05f) {
                    use_refined = false;
                    if (opts.verbose) {
                        std::cerr << "  [WARNING] Initial profile shows no damage signal (T/(T+C)="
                                  << std::fixed << std::setprecision(1) << initial_t_ratio * 100.0f
                                  << "% vs baseline=" << baseline_t_ratio * 100.0f << "%)\n";
                        std::cerr << "  [WARNING] Refined profile has no decay pattern (pos0="
                                  << refined_pos0_t * 100.0f << "% pos4=" << refined_pos4_t * 100.0f
                                  << "% decay=" << decay_signal * 100.0f << "%)\n";
                        std::cerr << "  [WARNING] Refined profile likely false positive from AT-biased reads\n";
                    }
                }
            } else if (suspicious_enrichment) {
                // Weak initial signal but massive enrichment - be cautious
                // Only trust if there's a clear decay pattern
                if (!has_decay_pattern) {
                    use_refined = false;
                    if (opts.verbose) {
                        std::cerr << "  [WARNING] Suspicious enrichment (" << std::setprecision(0)
                                  << enrichment_ratio << "x) with no decay pattern\n";
                        std::cerr << "  [WARNING] Initial d_max=" << std::fixed << std::setprecision(2)
                                  << initial_avg * 100.0f << "%, refined=" << refined_avg * 100.0f << "%\n";
                    }
                } else if (opts.verbose) {
                    std::cerr << "  [INFO] High enrichment (" << std::setprecision(0) << enrichment_ratio
                              << "x) but decay pattern present - trusting refined\n";
                }
            } else if (!initial_has_damage_signal && has_decay_pattern && opts.verbose) {
                std::cerr << "  [INFO] Initial profile has weak signal but refined shows decay pattern\n";
                std::cerr << "  [INFO] (pos0=" << std::fixed << std::setprecision(1)
                          << refined_pos0_t * 100.0f << "% pos4=" << refined_pos4_t * 100.0f
                          << "% decay=" << decay_signal * 100.0f << "%) - using refined\n";
            }

            if (use_refined) {
                // CRITICAL: Preserve hexamer fields from initial profile before overwriting
                // Hexamer analysis from Pass 1 characterizes the sample composition
                // and is needed for inverted pattern detection and d_max calculation
                auto saved_hexamer_terminal_tc = sample_profile.hexamer_terminal_tc;
                auto saved_hexamer_interior_tc = sample_profile.hexamer_interior_tc;
                auto saved_hexamer_excess_tc = sample_profile.hexamer_excess_tc;
                auto saved_hexamer_damage_llr = sample_profile.hexamer_damage_llr;
                auto saved_inverted_5prime = sample_profile.inverted_pattern_5prime;
                auto saved_inverted_3prime = sample_profile.inverted_pattern_3prime;

                // CRITICAL: Also save Pass 1's d_max values
                // When inverted pattern is detected, the refined profile's d_max is unreliable
                // because it's based on a biased subset (high-confidence coding genes)
                auto saved_d_max_5prime = sample_profile.d_max_5prime;
                auto saved_d_max_3prime = sample_profile.d_max_3prime;
                auto saved_d_max_combined = sample_profile.d_max_combined;
                auto saved_d_max_source = sample_profile.d_max_source;

                // CRITICAL: Save Pass 1's Channel B decision
                // Channel B (stop conversion) is our "smoking gun" for real damage
                // If Pass 1 identified the signal as artifact, preserve that decision
                auto saved_damage_artifact = sample_profile.damage_artifact;
                auto saved_damage_validated = sample_profile.damage_validated;
                auto saved_channel_b_valid = sample_profile.channel_b_valid;
                auto saved_stop_decay_llr = sample_profile.stop_decay_llr_5prime;
                auto saved_stop_baseline = sample_profile.stop_conversion_rate_baseline;
                auto saved_stop_amplitude = sample_profile.stop_amplitude_5prime;

                // Update the main sample_profile with refined estimates
                sample_profile = refined_profile;

                // Restore hexamer fields from initial analysis
                sample_profile.hexamer_terminal_tc = saved_hexamer_terminal_tc;
                sample_profile.hexamer_interior_tc = saved_hexamer_interior_tc;
                sample_profile.hexamer_excess_tc = saved_hexamer_excess_tc;
                sample_profile.hexamer_damage_llr = saved_hexamer_damage_llr;
                sample_profile.inverted_pattern_5prime = saved_inverted_5prime;
                sample_profile.inverted_pattern_3prime = saved_inverted_3prime;

                // CRITICAL: Restore Channel B fields from Pass 1
                // Pass 2+ doesn't track convertible codons, so we must preserve Pass 1's decision
                sample_profile.damage_artifact = saved_damage_artifact;
                sample_profile.damage_validated = saved_damage_validated;
                sample_profile.channel_b_valid = saved_channel_b_valid;
                sample_profile.stop_decay_llr_5prime = saved_stop_decay_llr;
                sample_profile.stop_conversion_rate_baseline = saved_stop_baseline;
                sample_profile.stop_amplitude_5prime = saved_stop_amplitude;

                // CRITICAL: If Pass 1 identified the signal as ARTIFACT via Channel B,
                // override the refined d_max values back to 0
                // This is the key fix for compositional artifacts like sample 0267130b40
                if (saved_damage_artifact) {
                    sample_profile.d_max_5prime = 0.0f;
                    sample_profile.d_max_3prime = 0.0f;
                    sample_profile.d_max_combined = 0.0f;
                    sample_profile.d_max_source = agp::SampleDamageProfile::DmaxSource::NONE;
                    if (opts.verbose) {
                        std::cerr << "  [ARTIFACT] Pass 1 Channel B identified compositional artifact - d_max set to 0\n";
                    }
                }
                // ALSO: If Pass 1 showed NO damage (both channels low/negative) AND
                // Channel B strongly contradicts damage (very negative stop_llr),
                // override to 0. Iterative refinement can re-inflate false positives.
                else if (!saved_damage_validated && saved_channel_b_valid &&
                         saved_stop_decay_llr < -100.0f && saved_d_max_combined < 0.01f) {
                    sample_profile.d_max_5prime = 0.0f;
                    sample_profile.d_max_3prime = 0.0f;
                    sample_profile.d_max_combined = 0.0f;
                    sample_profile.d_max_source = agp::SampleDamageProfile::DmaxSource::NONE;
                    if (opts.verbose) {
                        std::cerr << "  [NO-DAMAGE] Pass 1 showed no damage (stop_llr=" << saved_stop_decay_llr
                                  << ") - rejecting refined d_max\n";
                    }
                }
                // If both ends showed inverted pattern, use Pass 1's d_max
                else if (saved_inverted_5prime && saved_inverted_3prime) {
                    sample_profile.d_max_5prime = saved_d_max_5prime;
                    sample_profile.d_max_3prime = saved_d_max_3prime;
                    sample_profile.d_max_combined = saved_d_max_combined;
                    sample_profile.d_max_source = saved_d_max_source;
                }
                // ALWAYS use Pass 1's d_max when damage was validated
                // Pass 2+ selects damage-indicative reads, which biases d_max upward
                // Pass 1's estimate on ALL reads is unbiased and matches metaDMG methodology
                else if (saved_damage_validated && saved_d_max_combined > 0.01f) {
                    sample_profile.d_max_5prime = saved_d_max_5prime;
                    sample_profile.d_max_3prime = saved_d_max_3prime;
                    sample_profile.d_max_combined = saved_d_max_combined;
                    sample_profile.d_max_source = saved_d_max_source;
                    if (opts.verbose) {
                        std::cerr << "  [D_MAX] Using Pass 1 estimate (unbiased): "
                                  << saved_d_max_combined * 100.0f << "%\n";
                    }
                }
                // Otherwise keep d_max from refined profile (for low-damage or unvalidated samples)

                // Update damage model with refined profile
                damage_model.update_from_sample_profile(sample_profile);

                // Update per-domain profiles/models for meta mode
                if (opts.metagenome_mode) {
                    for (int d = 0; d < 8; ++d) {
                        agp::FrameSelector::finalize_sample_profile(domain_profiles[d]);
                        domain_damage_models[d].update_from_sample_profile(domain_profiles[d]);
                    }
                }

                // Re-initialize adaptive calibrator with refined profile
                calibrator.initialize(sample_profile);

                if (opts.verbose) {
                    std::cerr << "  [DEBUG] Using refined damage profile for gene prediction\n";
                    std::cerr << "  [DEBUG] Library type: " << refined_profile.library_type_str() << "\n";
                }
            } else if (opts.verbose) {
                if (!initial_has_damage_signal && refined_avg > 0.05f) {
                    std::cerr << "  [DEBUG] Rejecting refined profile (false positive safeguard)\n";
                } else {
                    std::cerr << "  [DEBUG] Refined profile has insufficient signal, keeping original\n";
                }
            }

            std::cerr << "\n";
        }

        // Main processing with double-buffering
        if (opts.aggregate_damage && opts.use_damage) {
            if (opts.iterative_damage) {
                std::cerr << "[Pass 3] Gene prediction with refined damage model...\n";
            } else {
                std::cerr << "[Pass 2] Gene prediction...\n";
            }
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
        size_t total_stop_restorations_all = 0;
        size_t total_corrections_reverted = 0;
        size_t total_low_confidence_skipped = 0;
        std::atomic<size_t> total_preframe_corrections{0};  // Pre-frame DNA corrections
        std::atomic<size_t> total_reads_precorrected{0};    // Reads with pre-frame corrections

        // Frame-aware wobble position damage tracking
        // For correct wobble_ratio, we need to use predicted reading frame
        std::array<std::atomic<size_t>, 3> wobble_t_count = {};  // T counts at codon positions 0,1,2
        std::array<std::atomic<size_t>, 3> wobble_c_count = {};  // C counts at codon positions 0,1,2

        // Per-domain terminal damage stats (meta mode only)
        std::array<size_t, 8> global_domain_terminal_t = {};  // T counts per domain
        std::array<size_t, 8> global_domain_terminal_c = {};  // C counts per domain
        std::array<size_t, 8> global_domain_gene_counts = {}; // Gene counts per domain

        std::mutex hist_mutex;

        // Thread-local terminal stats for adaptive calibration
        std::vector<agp::TerminalKmerStats> thread_terminal_stats(num_threads);

        // Thread-local domain hexamer stats for ensemble weight calibration
        struct ThreadDomainStats {
            std::array<double, 8> hexamer_sums = {};
            std::array<size_t, 8> counts = {};
            double coding_score_before = 0.0;
            double coding_score_after = 0.0;
            double hexamer_llr_before = 0.0;
            double hexamer_llr_after = 0.0;
            size_t genes_processed = 0;
            size_t stop_restorations = 0;  // Stops corrected back to sense codons
            size_t corrections_reverted = 0;  // Corrections reverted due to frame instability
            size_t low_confidence_skipped = 0;  // Corrections skipped due to low frame confidence
            // Per-domain terminal damage tracking (5' T/(T+C) ratio proxy)
            std::array<size_t, 8> domain_terminal_t = {};  // T counts at 5' terminal (first 3bp)
            std::array<size_t, 8> domain_terminal_c = {};  // C counts at 5' terminal (first 3bp)
        };
        std::vector<ThreadDomainStats> thread_domain_stats(num_threads);

        // Thread-local FDR score storage to avoid mutex contention
        struct ThreadFDRScores {
            std::vector<float> target_scores;
            std::vector<float> decoy_scores;
        };
        std::vector<ThreadFDRScores> thread_fdr_scores(num_threads);

        // Buffer for deferred output (FDR filtering or percentile computation)
        // Stores (id, genes) pairs for deferred writing
        std::vector<std::pair<std::string, std::vector<agp::Gene>>> fdr_gene_buffer;
        bool defer_output_for_fdr = (opts.fdr_threshold >= 0.0f && opts.compute_fdr) ||
                                    (opts.min_damage_pctile > 0.0f);

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

                // Update thread-local terminal stats for adaptive calibration
                if (opts.use_damage && seq.length() >= 30) {
                    int tid = 0;
                    #ifdef _OPENMP
                    tid = omp_get_thread_num();
                    #endif
                    auto& tstats = thread_terminal_stats[tid];

                    // 5' terminal (positions 0-5)
                    for (size_t j = 0; j < 6 && j < seq.length(); ++j) {
                        char b = agp::fast_upper(seq[j]);
                        if (b == 'T') tstats.terminal_5prime_t[j]++;
                        else if (b == 'C') tstats.terminal_5prime_c[j]++;
                    }

                    // 3' terminal (positions 0-5 from end)
                    for (size_t j = 0; j < 6 && j < seq.length(); ++j) {
                        char b = agp::fast_upper(seq[seq.length() - 1 - j]);
                        if (b == 'A') tstats.terminal_3prime_a[j]++;
                        else if (b == 'G') tstats.terminal_3prime_g[j]++;
                    }

                    // Interior (positions 15-25)
                    for (size_t j = 15; j < 25 && j < seq.length(); ++j) {
                        char b = agp::fast_upper(seq[j]);
                        if (b == 'T') tstats.interior_t++;
                        else if (b == 'C') tstats.interior_c++;
                        else if (b == 'A') tstats.interior_a++;
                        else if (b == 'G') tstats.interior_g++;
                    }
                    tstats.n_reads++;
                }

                // Enable ensemble mode for metagenome scoring (weighted across all domains)
                agp::MultiDomainResult domain_result;
                if (opts.metagenome_mode && seq.length() >= 30) {
                    domain_result = agp::score_all_domains_fast(seq, 0);
                    agp::set_ensemble_mode(true);
                    agp::set_domain_probs(domain_result);
                    // Use best domain for position-dependent damage profile
                    agp::set_active_domain(domain_result.best_domain);
                } else {
                    agp::set_ensemble_mode(false);
                    agp::set_active_domain(base_domain);
                }

                // Create damage profile
                const agp::DamageProfile* profile_ptr = nullptr;
                agp::DamageProfile local_profile;
                agp::DamageProfile blended_profile;
                if (opts.use_damage) {
                    if (opts.metagenome_mode && seq.length() >= 30) {
                        // Blend per-domain profiles using ensemble weights
                        // Optimized: array-based cache indexed by length bucket (avoids hash map)
                        static thread_local struct {
                            std::array<agp::DamageProfile, 64> profiles;  // Buckets for lengths 0-31, 32-63, ...
                            std::array<size_t, 64> cached_len;
                        } profile_cache = {};

                        const size_t len = seq.length();
                        const size_t bucket = std::min(size_t(63), len / 8);  // 8bp per bucket
                        agp::DamageProfile& bp = profile_cache.profiles[bucket];

                        // Resize vectors only if needed (reuse allocation)
                        if (bp.ct_prob_5prime.size() != len) {
                            bp.ct_prob_5prime.resize(len);
                            bp.ga_prob_3prime.resize(len);
                            bp.ga_prob_5prime.resize(len);
                            bp.ct_prob_3prime.resize(len);
                        }

                        // Zero scalars and vectors
                        bp.lambda_5prime = 0.0f;
                        bp.lambda_3prime = 0.0f;
                        bp.delta_max = 0.0f;
                        bp.delta_background = 0.0f;
                        std::memset(bp.ct_prob_5prime.data(), 0, len * sizeof(float));
                        std::memset(bp.ga_prob_3prime.data(), 0, len * sizeof(float));

                        const auto& probs = agp::get_domain_probs();
                        const float weights[8] = {
                            probs.gtdb_prob, probs.fungi_prob, probs.protozoa_prob, probs.invertebrate_prob,
                            probs.plant_prob, probs.vertebrate_mammalian_prob, probs.vertebrate_other_prob,
                            probs.viral_prob
                        };

                        // Find best domain (skip blending if one dominates)
                        int best_d = 0;
                        float best_w = weights[0];
                        for (int d = 1; d < 8; ++d) {
                            if (weights[d] > best_w) { best_w = weights[d]; best_d = d; }
                        }

                        // If one domain dominates (>90%), just use that profile directly
                        if (best_w > 0.90f) {
                            blended_profile = domain_damage_models[best_d].create_profile_cached(len);
                            profile_ptr = &blended_profile;
                        } else {
                            // Blend only significant domains (>1%)
                            float weight_sum = 0.0f;
                            for (int d = 0; d < 8; ++d) {
                                if (weights[d] <= 0.01f) continue;
                                weight_sum += weights[d];
                                const auto& dom_profile = domain_damage_models[d].create_profile_cached(len);
                                bp.lambda_5prime += weights[d] * dom_profile.lambda_5prime;
                                bp.lambda_3prime += weights[d] * dom_profile.lambda_3prime;
                                bp.delta_max += weights[d] * dom_profile.delta_max;
                                bp.delta_background += weights[d] * dom_profile.delta_background;
                                const float* ct_src = dom_profile.ct_prob_5prime.data();
                                const float* ga_src = dom_profile.ga_prob_3prime.data();
                                float* ct_dst = bp.ct_prob_5prime.data();
                                float* ga_dst = bp.ga_prob_3prime.data();
                                const float w = weights[d];
                                for (size_t p = 0; p < len; ++p) {
                                    ct_dst[p] += w * ct_src[p];
                                    ga_dst[p] += w * ga_src[p];
                                }
                            }
                            if (weight_sum > 0.0f) {
                                float inv = 1.0f / weight_sum;
                                bp.lambda_5prime *= inv;
                                bp.lambda_3prime *= inv;
                                bp.delta_max *= inv;
                                bp.delta_background *= inv;
                                for (size_t p = 0; p < len; ++p) {
                                    bp.ct_prob_5prime[p] *= inv;
                                    bp.ga_prob_3prime[p] *= inv;
                                }
                            }
                            blended_profile = bp;
                            profile_ptr = &blended_profile;
                        }
                    } else {
                        local_profile = damage_model.create_profile(seq.length());
                        profile_ptr = &local_profile;
                    }
                }

                std::vector<agp::Gene> genes;

                // Pre-frame-selection damage correction
                // Apply conservative corrections to restore hexamer signal before frame selection
                // This helps frame selection by undoing damage that corrupts the hexamer patterns
                std::string seq_for_frame = seq;  // Default: use original
                size_t preframe_corrections = 0;
                if (opts.use_damage && sample_profile.max_damage_5prime >= 0.10f) {
                    auto [precorrected, n_precorr] = agp::DamageModel::correct_damage_preframe(seq, sample_profile);
                    if (n_precorr > 0) {
                        seq_for_frame = std::move(precorrected);
                        preframe_corrections = n_precorr;
                        total_preframe_corrections += n_precorr;
                        total_reads_precorrected++;
                    }
                }

                // Handle --all-frames mode: output all 6 reading frames
                if (opts.all_frames) {
                    std::vector<agp::FrameScore> all_frames = agp::FrameSelector::score_all_frames_full(seq_for_frame, profile_ptr);
                    float damage_pct = opts.use_damage
                        ? agp::FrameSelector::compute_damage_percentage(seq, sample_profile)
                        : 0.0f;

                    for (const auto& frame_score : all_frames) {
                        // Apply same min_coding_prob threshold as other modes for consistency
                        if (frame_score.total_score < opts.min_coding_prob) continue;

                        agp::Gene gene;
                        gene.start = 0;
                        gene.end = seq.length();
                        gene.is_forward = frame_score.forward;
                        gene.is_fragment = true;
                        gene.coding_prob = frame_score.total_score;
                        gene.score = agp::normalize_coding_score(frame_score.total_score);
                        gene.frame = frame_score.frame;
                        gene.damage_score = damage_pct;

                        if (frame_score.forward) {
                            gene.sequence = seq;
                            gene.damage_signal = agp::FrameSelector::estimate_damage_signal_codon_aware(
                                seq, frame_score.frame, sample_profile);
                        } else {
                            gene.sequence = agp::FrameSelector::reverse_complement(seq);
                            int orig_frame = (seq.length() - frame_score.frame) % 3;
                            gene.damage_signal = agp::FrameSelector::estimate_damage_signal_codon_aware(
                                seq, orig_frame, sample_profile);
                        }
                        gene.protein = frame_score.protein;
                        gene.p_correct = agp::calibrate_score_4d(frame_score.total_score, gene.sequence.size(),
                                                                  gene.damage_score, gene.damage_signal);
                        genes.push_back(std::move(gene));
                    }
                } else if (opts.use_bayesian) {
                    // Bayesian V2 frame selection using dicodon phase hexamers
                    // Uses phase-specific frequencies for better frame discrimination
                    // Includes damage-aware stop codon scoring and correction
                    agp::BayesianFrameSelectorV2 bayesian(active_domain);

                    // Set damage parameters from sample profile
                    if (sample_profile.is_valid()) {
                        bayesian.set_damage_params(
                            sample_profile.max_damage_5prime,
                            sample_profile.max_damage_3prime,
                            sample_profile.lambda_5prime,
                            sample_profile.lambda_3prime);
                    }

                    auto output = bayesian.select_frames(seq_for_frame);

                    // Compute damage percentage (terminal pattern-based)
                    float damage_pct = opts.use_damage
                        ? agp::FrameSelector::compute_damage_percentage(seq, sample_profile)
                        : 0.0f;

                    // Create genes from selected frames
                    for (const auto& frame_result : output.selected_frames) {
                        agp::Gene gene;
                        gene.start = 0;
                        gene.end = seq.length();
                        gene.is_forward = !frame_result.is_reverse;
                        gene.is_fragment = true;
                        gene.coding_prob = frame_result.posterior;
                        gene.score = frame_result.posterior;
                        gene.frame = frame_result.frame;
                        gene.damage_score = damage_pct;

                        // Compute codon-aware per-read damage probability
                        // Output log-LR and informativeness for downstream aggregation
                        if (opts.use_damage && sample_profile.max_damage_5prime > 0) {
                            auto dmg_result = bayesian.compute_damage_probability(
                                frame_result.sequence, frame_result.frame, 0.5f, frame_result.is_reverse);
                            // Store all three: p_damaged (human-readable), log_lr (for aggregation), info (weight)
                            gene.p_read_damaged = 1.0f / (1.0f + std::exp(-dmg_result.log_lr));
                            gene.damage_log_lr = dmg_result.log_lr;
                            gene.damage_info = dmg_result.informativeness;
                        }

                        // Use the oriented sequence from V2
                        gene.sequence = frame_result.sequence;
                        gene.protein = frame_result.protein;
                        gene.corrected_sequence = frame_result.corrected_sequence;
                        gene.corrected_protein = frame_result.corrected_protein;
                        gene.dna_corrections = frame_result.nt_corrections;
                        gene.stop_restorations = frame_result.stops_corrected;

                        // Estimate ancient probability
                        if (gene.is_forward) {
                            gene.damage_signal = agp::FrameSelector::estimate_damage_signal_codon_aware(
                                seq, frame_result.frame, sample_profile);
                        } else {
                            int orig_frame = (seq.length() - frame_result.frame) % 3;
                            gene.damage_signal = agp::FrameSelector::estimate_damage_signal_codon_aware(
                                seq, orig_frame, sample_profile);
                        }

                        // Use posterior as p_correct (already calibrated by softmax)
                        gene.p_correct = frame_result.posterior;

                        genes.push_back(std::move(gene));
                    }
                } else if (opts.dual_strand) {
                    // Output best frame from BOTH strands
                    auto strand_results = agp::FrameSelector::select_best_per_strand(seq_for_frame, profile_ptr, &sample_profile);
                    auto [best_fwd, best_rev] = strand_results;

                    // Calculate damage_signal using codon-aware detection for each strand
                    // This uses the reading frame information for better accuracy
                    float damage_signal_fwd = 0.5f;
                    float damage_signal_rev = 0.5f;

                    if (best_fwd.total_score >= opts.min_coding_prob) {
                        // Use quality-aware detection if quality available, else codon-aware
                        if (!quality.empty()) {
                            damage_signal_fwd = agp::FrameSelector::estimate_damage_signal_with_quality(
                                seq, quality, best_fwd.frame, sample_profile);
                        } else {
                            damage_signal_fwd = agp::FrameSelector::estimate_damage_signal_codon_aware(
                                seq, best_fwd.frame, sample_profile);
                        }
                    }

                    if (best_rev.total_score >= opts.min_coding_prob) {
                        // For reverse strand: damage is on ORIGINAL read, not RC
                        // Use original seq for damage scoring, but adjust frame for RC
                        // The frame in RC corresponds to (len - frame) % 3 in original
                        int orig_frame = (seq.length() - best_rev.frame) % 3;
                        if (!quality.empty()) {
                            damage_signal_rev = agp::FrameSelector::estimate_damage_signal_with_quality(
                                seq, quality, orig_frame, sample_profile);
                        } else {
                            damage_signal_rev = agp::FrameSelector::estimate_damage_signal_codon_aware(
                                seq, orig_frame, sample_profile);
                        }
                    }

                    // Determine which strands to output
                    bool output_fwd = (best_fwd.total_score >= opts.min_coding_prob);
                    bool output_rev = (best_rev.total_score >= opts.min_coding_prob);

                    // Adaptive confidence-based ORF output
                    // When enabled, uses score gap to decide between 1 or 2 ORFs
                    if (opts.orf_confidence_threshold >= 0.0f && output_fwd && output_rev) {
                        // Use combined score with configurable damage_signal weight
                        float score_fwd = best_fwd.total_score * (1.0f + opts.strand_weight * damage_signal_fwd);
                        float score_rev = best_rev.total_score * (1.0f + opts.strand_weight * damage_signal_rev);
                        
                        // Calculate score gap between top-1 and top-2
                        float score_gap = std::abs(score_fwd - score_rev);
                        
                        // If gap >= threshold, we're confident enough to output only 1 ORF
                        if (score_gap >= opts.orf_confidence_threshold) {
                            if (score_fwd >= score_rev) {
                                output_rev = false;
                            } else {
                                output_fwd = false;
                            }
                        }
                        // Otherwise, output both ORFs (low confidence)
                    }
                    // If --best-strand (legacy mode), only output the higher-scoring strand
                    else if (opts.best_strand_only && output_fwd && output_rev) {
                        // Use combined score with configurable damage_signal weight
                        // strand_weight=0: pure coding, strand_weight=1: full damage bias
                        float score_fwd = best_fwd.total_score * (1.0f + opts.strand_weight * damage_signal_fwd);
                        float score_rev = best_rev.total_score * (1.0f + opts.strand_weight * damage_signal_rev);
                        if (score_fwd >= score_rev) {
                            output_rev = false;
                        } else {
                            output_fwd = false;
                        }
                    }

                    // Compute damage percentage using sample-level distribution
                    // (damage is always on original read, not reverse complement)
                    float damage_pct = opts.use_damage
                        ? agp::FrameSelector::compute_damage_percentage(seq, sample_profile)
                        : 0.0f;

                    if (output_fwd) {
                        agp::Gene gene;
                        gene.start = 0;
                        gene.end = seq.length();
                        gene.is_forward = true;
                        gene.is_fragment = true;
                        gene.coding_prob = best_fwd.total_score;
                        gene.score = agp::normalize_coding_score(best_fwd.total_score);
                        gene.frame = best_fwd.frame;
                        gene.sequence = seq;
                        gene.protein = std::move(best_fwd.protein);
                        gene.damage_signal = damage_signal_fwd;
                        gene.damage_score = damage_pct;
                        gene.p_correct = agp::calibrate_score_4d(best_fwd.total_score, gene.sequence.size(),
                                                                  damage_pct, damage_signal_fwd);
                        genes.push_back(std::move(gene));
                    }

                    if (output_rev) {
                        agp::Gene gene;
                        gene.start = 0;
                        gene.end = seq.length();
                        gene.is_forward = false;
                        gene.is_fragment = true;
                        gene.coding_prob = best_rev.total_score;
                        gene.score = agp::normalize_coding_score(best_rev.total_score);
                        gene.frame = best_rev.frame;
                        // Store reverse complement as the coding sequence
                        gene.sequence = agp::FrameSelector::reverse_complement(seq);
                        gene.protein = std::move(best_rev.protein);
                        gene.damage_signal = damage_signal_rev;
                        gene.damage_score = damage_pct;
                        gene.p_correct = agp::calibrate_score_4d(best_rev.total_score, gene.sequence.size(),
                                                                  damage_pct, damage_signal_rev);
                        genes.push_back(std::move(gene));
                    }
                } else {
                    // Single best frame
                    // Use pre-corrected sequence for frame selection to restore hexamer signal
                    agp::FrameScore best = agp::FrameSelector::select_best_frame(seq_for_frame, profile_ptr);

                    if (best.total_score >= opts.min_coding_prob) {
                        agp::Gene gene;
                        gene.start = 0;
                        gene.end = seq.length();
                        gene.is_forward = best.forward;
                        gene.is_fragment = true;
                        gene.coding_prob = best.total_score;
                        gene.score = agp::normalize_coding_score(best.total_score);
                        gene.frame = best.frame;

                        // Use quality-aware damage detection if available, else codon-aware
                        // Store coding sequence (reverse complement for reverse strand)
                        // Use cached RC to avoid allocation in hot loop
                        std::string working_seq;
                        if (best.forward) {
                            working_seq = seq;
                        } else {
                            const std::string& rc = agp::FrameSelector::reverse_complement_cached(seq);
                            working_seq = rc;  // Copy from cached buffer
                        }
                        if (!quality.empty()) {
                            gene.damage_signal = agp::FrameSelector::estimate_damage_signal_with_quality(
                                working_seq, quality, best.frame, sample_profile);
                        } else {
                            gene.damage_signal = agp::FrameSelector::estimate_damage_signal_codon_aware(
                                working_seq, best.frame, sample_profile);
                        }

                        // Compute damage percentage using sample-level distribution
                        gene.damage_score = opts.use_damage
                            ? agp::FrameSelector::compute_damage_percentage(seq, sample_profile)
                            : 0.0f;

                        // Calibrate score to probability of correctness (4D: score, length, damage, ancient)
                        gene.p_correct = agp::calibrate_score_4d(best.total_score, seq.length(),
                                                                  gene.damage_score, gene.damage_signal);

                        gene.sequence = std::move(working_seq);
                        gene.protein = std::move(best.protein);
                        genes.push_back(std::move(gene));
                    }
                }

                // Apply damage correction if confident
                // Use adaptive damage correction with calibrator for smarter decisions
                size_t batch_dna_corrections = 0;
                size_t batch_aa_corrections = 0;
                size_t batch_genes_corrected = 0;

                // Get thread ID for domain stats tracking
                int tid = 0;
                #ifdef _OPENMP
                tid = omp_get_thread_num();
                #endif
                auto& dstats = thread_domain_stats[tid];

                if (!genes.empty() && profile_ptr && opts.use_damage) {
                    for (auto& gene : genes) {
                        // Track coding score before correction using frame total_score
                        // This gives consistent scaling with after-correction scores
                        float score_before = gene.coding_prob;  // Default fallback
                        float frame_confidence = 1.0f;  // Margin between best and 2nd best frame score
                        if (gene.sequence.length() >= 30) {
                            auto before_scores = agp::FrameSelector::score_all_frames_full(
                                gene.sequence, profile_ptr);

                            // Find best and 2nd best scores to compute confidence margin
                            float best_score = -1e9f;
                            float second_best = -1e9f;
                            for (const auto& fs : before_scores) {
                                if (fs.frame == gene.frame && fs.forward == gene.is_forward) {
                                    score_before = fs.total_score;
                                }
                                if (fs.total_score > best_score) {
                                    second_best = best_score;
                                    best_score = fs.total_score;
                                } else if (fs.total_score > second_best) {
                                    second_best = fs.total_score;
                                }
                            }
                            // Frame confidence = margin between best and 2nd best
                            // High margin = confident prediction, safe to correct
                            // Low margin = marginal prediction, risky to correct
                            frame_confidence = best_score - second_best;
                        }
                        dstats.coding_score_before += score_before;

                        // Track per-domain hexamer LLRs for ensemble weight calibration
                        if (gene.sequence.length() >= 30) {
                            static const std::array<agp::Domain, 8> domains = {
                                agp::Domain::GTDB, agp::Domain::FUNGI, agp::Domain::PROTOZOA,
                                agp::Domain::INVERTEBRATE, agp::Domain::PLANT,
                                agp::Domain::VERTEBRATE_MAMMALIAN, agp::Domain::VERTEBRATE_OTHER,
                                agp::Domain::VIRAL
                            };

                            // Compute hexamer LLRs using same formula as score_all_domains()
                            // log2(freq * 4096 + 1e-10) for consistent scaling
                            // Optimized: encode hexamer once, look up freq across domains
                            std::array<float, 8> domain_llrs = {};
                            size_t n_hex = 0;
                            size_t start = (gene.frame >= 0) ? static_cast<size_t>(gene.frame) : 0;
                            for (size_t k = start; k + 6 <= gene.sequence.length(); k += 3) {
                                // Encode hexamer once
                                uint32_t code = agp::encode_hexamer(gene.sequence.data() + k);
                                if (code >= 4096) continue;
                                // Look up freq for each domain
                                for (size_t d = 0; d < 8; ++d) {
                                    float freq = agp::get_hexamer_freq(code, domains[d]);
                                    domain_llrs[d] += std::log2(freq * 4096.0f + 1e-10f);
                                }
                                n_hex++;
                            }
                            if (n_hex > 0) {
                                float inv_n = 1.0f / static_cast<float>(n_hex);
                                size_t best_domain = 0;
                                float best_llr = domain_llrs[0];
                                for (size_t d = 0; d < 8; ++d) {
                                    dstats.hexamer_sums[d] += domain_llrs[d] * inv_n;
                                    dstats.counts[d]++;
                                    if (domain_llrs[d] > best_llr) {
                                        best_llr = domain_llrs[d];
                                        best_domain = d;
                                    }
                                }
                                // Track terminal damage for best-matching domain
                                // Count T and C in first 3 positions for 5' T/(T+C) ratio
                                const std::string& seq = gene.sequence;
                                for (size_t p = 0; p < std::min<size_t>(3, seq.length()); ++p) {
                                    char nt = agp::fast_upper(seq[p]);
                                    if (nt == 'T') dstats.domain_terminal_t[best_domain]++;
                                    else if (nt == 'C') dstats.domain_terminal_c[best_domain]++;
                                }
                            }
                            dstats.hexamer_llr_before += calibrator.get_ensemble_hexamer_llr(
                                gene.sequence, gene.frame);
                        }

                        // Frame confidence gating: only attempt corrections if frame prediction is confident
                        // Low confidence (small margin between best and 2nd best) = risky to correct
                        // Threshold: 0.05 = 5% score margin required for correction
                        constexpr float min_frame_confidence = 0.05f;
                        bool skip_correction = (frame_confidence < min_frame_confidence);

                        std::string corrected_dna;
                        size_t dna_corr = 0;

                        if (!skip_correction) {
                            // Use adaptive correction that validates each correction
                            // Calibrator handles: per-read LLR gating, max corrections,
                            // and correction validation (stop codon check, hexamer quality)
                            // IMPORTANT: Pass strand orientation for correct damage profile application
                            // Forward strand: C→T at 5' end, G→A at 3' end
                            // Reverse strand: G→A at 5' end, C→T at 3' end (inverted due to revcomp)
                            auto result = damage_model.correct_damage_adaptive(
                                gene.sequence, sample_profile, gene.damage_signal, gene.frame,
                                calibrator, gene.is_forward);
                            corrected_dna = std::move(result.first);
                            dna_corr = result.second;
                        } else {
                            corrected_dna = gene.sequence;
                            dstats.low_confidence_skipped++;
                        }
                        gene.corrected_sequence = corrected_dna;
                        gene.dna_corrections = dna_corr;

                        // Re-translate corrected DNA to get corrected protein
                        if (dna_corr > 0) {
                            batch_dna_corrections += dna_corr;
                            batch_genes_corrected++;

                            gene.corrected_protein.clear();
                            for (size_t j = gene.frame; j + 2 < corrected_dna.length(); j += 3) {
                                char c1 = corrected_dna[j];
                                char c2 = corrected_dna[j + 1];
                                char c3 = corrected_dna[j + 2];
                                char aa = agp::CodonTable::translate_codon(c1, c2, c3);
                                gene.corrected_protein += aa;
                            }

                            // Count amino acid changes and stop restorations
                            gene.aa_corrections = 0;
                            gene.stop_restorations = 0;
                            for (size_t j = 0; j < std::min(gene.protein.length(), gene.corrected_protein.length()); ++j) {
                                if (gene.protein[j] != gene.corrected_protein[j]) {
                                    gene.aa_corrections++;
                                    // Track stop restorations (original was stop, corrected is not)
                                    if (gene.protein[j] == '*' && gene.corrected_protein[j] != '*') {
                                        gene.stop_restorations++;
                                        dstats.stop_restorations++;
                                    }
                                }
                            }
                            batch_aa_corrections += gene.aa_corrections;

                            // Track hexamer LLR after correction
                            if (corrected_dna.length() >= 30) {
                                dstats.hexamer_llr_after += calibrator.get_ensemble_hexamer_llr(
                                    corrected_dna, gene.frame);
                            }

                            // Re-score corrected sequence and verify frame stability
                            if (corrected_dna.length() >= 30) {
                                auto corrected_scores = agp::FrameSelector::score_all_frames_full(
                                    corrected_dna, profile_ptr);

                                // Find the best frame/strand from corrected scores
                                int best_frame = gene.frame;
                                bool best_forward = gene.is_forward;
                                float best_score = -1e9f;
                                float original_frame_score = 0.0f;

                                for (const auto& fs : corrected_scores) {
                                    if (fs.total_score > best_score) {
                                        best_score = fs.total_score;
                                        best_frame = fs.frame;
                                        best_forward = fs.forward;
                                    }
                                    if (fs.frame == gene.frame && fs.forward == gene.is_forward) {
                                        original_frame_score = fs.total_score;
                                    }
                                }

                                // Frame stability check: revert if corrections change best frame/strand
                                if (best_frame != gene.frame || best_forward != gene.is_forward) {
                                    // Corrections destabilized frame selection - revert
                                    gene.corrected_sequence = gene.sequence;
                                    gene.dna_corrections = 0;
                                    gene.corrected_protein = gene.protein;
                                    gene.aa_corrections = 0;
                                    gene.corrected_coding_prob = score_before;
                                    dstats.coding_score_after += score_before;
                                    dstats.corrections_reverted++;
                                    // Undo batch counters (already incremented)
                                    batch_dna_corrections -= dna_corr;
                                    batch_genes_corrected--;
                                    batch_aa_corrections -= gene.aa_corrections;
                                } else {
                                    // Frame is stable - keep corrections
                                    gene.corrected_coding_prob = original_frame_score;
                                    dstats.coding_score_after += original_frame_score;
                                }
                            } else {
                                gene.corrected_coding_prob = score_before;
                                dstats.coding_score_after += score_before;  // Short seq: use same as before
                            }
                        } else {
                            // No correction - after equals before
                            gene.corrected_coding_prob = score_before;  // Same as original
                            if (gene.sequence.length() >= 30) {
                                dstats.hexamer_llr_after += calibrator.get_ensemble_hexamer_llr(
                                    gene.sequence, gene.frame);
                            }
                            dstats.coding_score_after += score_before;  // No change: after == before
                        }
                        dstats.genes_processed++;
                    }
                }

                // FDR score collection and decoy scoring
                if (opts.compute_fdr && !genes.empty()) {
                    auto& fdr_scores = thread_fdr_scores[tid];

                    // Collect target scores
                    for (const auto& gene : genes) {
                        fdr_scores.target_scores.push_back(gene.score);
                    }

                    // Generate and score decoys using selected strategy
                    // Use the original sequence (not pre-corrected) for fair comparison
                    for (size_t d = 0; d < opts.fdr_decoys; ++d) {
                        std::string decoy_seq = decoy_gen->generate(seq, fdr_decoy_strategy);

                        // Score the decoy using the same method as targets
                        // Apply same normalization as targets for fair comparison
                        if (opts.dual_strand) {
                            auto [best_fwd, best_rev] = agp::FrameSelector::select_best_per_strand(
                                decoy_seq, profile_ptr);
                            float raw_score = std::max(best_fwd.total_score, best_rev.total_score);
                            fdr_scores.decoy_scores.push_back(agp::normalize_coding_score(raw_score));
                        } else {
                            agp::FrameScore best = agp::FrameSelector::select_best_frame(decoy_seq, profile_ptr);
                            fdr_scores.decoy_scores.push_back(agp::normalize_coding_score(best.total_score));
                        }
                    }
                }

                results[i] = {rec.id, std::move(genes)};
            }

            // Write results and track damage probability distribution
            size_t batch_total_dna_corr = 0;
            size_t batch_total_aa_corr = 0;
            size_t batch_total_genes_corr = 0;
            size_t batch_total_stop_restor = 0;

            for (auto& [id, genes] : results) {
                if (!genes.empty()) {
                    // Apply confidence filtering if min_confidence > 0
                    if (opts.min_confidence > 0.0f) {
                        genes.erase(
                            std::remove_if(genes.begin(), genes.end(),
                                [&opts](const agp::Gene& g) {
                                    return g.p_correct < opts.min_confidence;
                                }),
                            genes.end());
                    }

                    if (genes.empty()) continue;

                    // Either buffer for FDR filtering or write immediately
                    if (defer_output_for_fdr) {
                        fdr_gene_buffer.push_back({id, genes});
                    } else {
                        writer.write_genes(id, genes);
                        if (fasta_nt_writer) fasta_nt_writer->write_genes_nucleotide(id, genes);
                        if (fasta_nt_corr_writer) fasta_nt_corr_writer->write_genes_nucleotide_corrected(id, genes);
                        if (fasta_aa_writer) fasta_aa_writer->write_genes_protein(id, genes);
                        if (fasta_corrected_writer) fasta_corrected_writer->write_genes_protein_corrected(id, genes);
                    }

                    total_genes += genes.size();

                    // Update damage probability and percentage histograms
                    for (const auto& gene : genes) {
                        int bin = std::min(9, static_cast<int>(gene.damage_signal * 10));
                        damage_prob_hist[bin]++;

                        if (opts.use_damage) {
                            // Damage percentage: 0-100, bin into 10% increments
                            int pct_bin = std::min(9, static_cast<int>(gene.damage_score / 10.0));
                            damage_pct_hist[pct_bin]++;
                            total_damage_pct += gene.damage_score;
                        }

                        // Track corrections
                        if (gene.dna_corrections > 0) {
                            batch_total_dna_corr += gene.dna_corrections;
                            batch_total_aa_corr += gene.aa_corrections;
                            batch_total_stop_restor += gene.stop_restorations;
                            batch_total_genes_corr++;
                        }

                        // Frame-aware wobble counting (first 15 bases from 5' end)
                        // Use original read orientation for codon positions even for reverse hits
                        int wobble_frame = gene.frame;
                        const std::string* seq_ptr = &gene.sequence;
                        if (!gene.is_forward) {
                            // Use cached RC to avoid allocation (buffer valid until next RC call)
                            const std::string& rc_seq = agp::FrameSelector::reverse_complement_cached(gene.sequence);
                            // Map frame back to original coordinates
                            wobble_frame = (static_cast<int>(rc_seq.length()) - gene.frame) % 3;
                            seq_ptr = &rc_seq;
                        }

                        const std::string& seq = *seq_ptr;
                        int frame = wobble_frame;
                        for (size_t i = 0; i < std::min<size_t>(15, seq.length()); ++i) {
                            // Frame-adjusted codon position: (i + frame) % 3
                            // frame=0: pos 0,1,2,3,4,5... -> codon_pos 0,1,2,0,1,2...
                            // frame=1: pos 0,1,2,3,4,5... -> codon_pos 1,2,0,1,2,0...
                            // frame=2: pos 0,1,2,3,4,5... -> codon_pos 2,0,1,2,0,1...
                            int codon_pos = (i + frame) % 3;
                            char base = agp::fast_upper(seq[i]);
                            if (base == 'T') {
                                wobble_t_count[codon_pos].fetch_add(1, std::memory_order_relaxed);
                            } else if (base == 'C') {
                                wobble_c_count[codon_pos].fetch_add(1, std::memory_order_relaxed);
                            }
                        }
                    }
                }
            }

            // Update global correction stats
            total_dna_corrections += batch_total_dna_corr;
            total_aa_corrections += batch_total_aa_corr;
            total_stop_restorations_all += batch_total_stop_restor;
            total_genes_corrected += batch_total_genes_corr;

            // Update adaptive calibrator with batch statistics
            // This allows the calibrator to adjust thresholds based on observed correction quality
            if (opts.use_damage && batch.size() > 0) {
                // Merge thread-local terminal stats into calibrator
                agp::TerminalKmerStats merged_stats;
                for (int t = 0; t < num_threads; ++t) {
                    auto& tstats = thread_terminal_stats[t];
                    for (size_t j = 0; j < 6; ++j) {
                        merged_stats.terminal_5prime_t[j] += tstats.terminal_5prime_t[j];
                        merged_stats.terminal_5prime_c[j] += tstats.terminal_5prime_c[j];
                        merged_stats.terminal_3prime_a[j] += tstats.terminal_3prime_a[j];
                        merged_stats.terminal_3prime_g[j] += tstats.terminal_3prime_g[j];
                    }
                    merged_stats.interior_t += tstats.interior_t;
                    merged_stats.interior_c += tstats.interior_c;
                    merged_stats.interior_a += tstats.interior_a;
                    merged_stats.interior_g += tstats.interior_g;
                    merged_stats.n_reads += tstats.n_reads;
                }

                // Merge batch terminal stats into calibrator's cumulative stats
                calibrator.merge_terminal_stats(merged_stats);

                // Compute decay slopes for decay monitoring
                float decay_5 = merged_stats.compute_decay_slope_5prime();
                float decay_3 = merged_stats.compute_decay_slope_3prime();
                calibrator.add_decay_sample(decay_5, decay_3);

                // Merge thread-local domain stats for ensemble weight calibration
                std::array<double, 8> merged_domain_sums = {};
                std::array<size_t, 8> merged_domain_counts = {};
                double total_coding_before = 0.0, total_coding_after = 0.0;
                double total_hexamer_before = 0.0, total_hexamer_after = 0.0;
                size_t total_genes_in_batch = 0;
                size_t total_stop_restorations = 0;
                size_t batch_corrections_reverted = 0;
                size_t batch_low_confidence_skipped = 0;

                for (int t = 0; t < num_threads; ++t) {
                    auto& ds = thread_domain_stats[t];
                    for (size_t d = 0; d < 8; ++d) {
                        merged_domain_sums[d] += ds.hexamer_sums[d];
                        merged_domain_counts[d] += ds.counts[d];
                        // Accumulate per-domain terminal damage stats
                        global_domain_terminal_t[d] += ds.domain_terminal_t[d];
                        global_domain_terminal_c[d] += ds.domain_terminal_c[d];
                        global_domain_gene_counts[d] += ds.counts[d];
                    }
                    total_coding_before += ds.coding_score_before;
                    total_coding_after += ds.coding_score_after;
                    total_hexamer_before += ds.hexamer_llr_before;
                    total_hexamer_after += ds.hexamer_llr_after;
                    total_genes_in_batch += ds.genes_processed;
                    total_stop_restorations += ds.stop_restorations;
                    batch_corrections_reverted += ds.corrections_reverted;
                    batch_low_confidence_skipped += ds.low_confidence_skipped;
                }

                // Feed batch stats to calibrator with real metrics
                size_t batch_reads = batch.size();
                size_t batch_corrected = batch_total_genes_corr;
                size_t batch_corrections = batch_total_dna_corr;

                // Create batch stats with real per-domain hexamer data
                agp::BatchStats batch_stats;
                batch_stats.n_reads = batch_reads;
                batch_stats.n_corrected = batch_corrected;
                batch_stats.total_corrections = batch_corrections;
                batch_stats.stop_restorations = total_stop_restorations;
                batch_stats.total_coding_score_before = total_coding_before;
                batch_stats.total_coding_score_after = total_coding_after;
                batch_stats.total_hexamer_llr_before = total_hexamer_before;
                batch_stats.total_hexamer_llr_after = total_hexamer_after;
                batch_stats.domain_hexamer_sums = merged_domain_sums;
                batch_stats.domain_counts = merged_domain_counts;

                // Update calibrator with real batch stats
                calibrator.update_batch_stats(
                    batch_reads, batch_corrected, batch_corrections, total_stop_restorations,
                    total_coding_before, total_coding_after,
                    total_hexamer_before, total_hexamer_after);

                // Update per-domain hexamer stats for ensemble weight learning
                calibrator.update_domain_stats(merged_domain_sums, merged_domain_counts);

                // Calibrate thresholds at end of each batch
                calibrator.calibrate_after_batch();

                // Update global counters for reverted and skipped corrections
                total_corrections_reverted += batch_corrections_reverted;
                total_low_confidence_skipped += batch_low_confidence_skipped;

                // Reset thread-local stats for next batch
                for (int t = 0; t < num_threads; ++t) {
                    thread_terminal_stats[t] = agp::TerminalKmerStats();
                    thread_domain_stats[t] = ThreadDomainStats();
                }
            }

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

        // Merge thread-local FDR scores after all processing
        agp::FDREstimator::FDRReport fdr_report;
        float fdr_score_threshold = 0.0f;
        if (opts.compute_fdr) {
            for (int t = 0; t < num_threads; ++t) {
                auto& fdr_scores = thread_fdr_scores[t];
                fdr_target_scores.insert(fdr_target_scores.end(),
                    fdr_scores.target_scores.begin(), fdr_scores.target_scores.end());
                fdr_decoy_scores.insert(fdr_decoy_scores.end(),
                    fdr_scores.decoy_scores.begin(), fdr_scores.decoy_scores.end());
            }

            // Compute FDR report
            if (!fdr_target_scores.empty() && !fdr_decoy_scores.empty()) {
                fdr_report = agp::FDREstimator::generate_report(fdr_target_scores, fdr_decoy_scores);

                // Determine threshold for filtering
                if (opts.fdr_threshold >= 0.0f) {
                    fdr_score_threshold = agp::FDREstimator::find_threshold_for_fdr(
                        opts.fdr_threshold, fdr_target_scores, fdr_decoy_scores);
                }

                std::cerr << "\n[FDR Estimation]\n";
                std::cerr << "  Targets: " << fdr_report.n_targets
                          << " | Decoys: " << fdr_report.n_decoys << "\n";
                std::cerr << "  Null distribution: mean=" << std::fixed << std::setprecision(3)
                          << fdr_report.empirical_null_mean << " sd=" << fdr_report.empirical_null_sd << "\n";
                std::cerr << "  At 1% FDR: " << fdr_report.genes_01 << " genes"
                          << " (threshold=" << std::setprecision(3) << fdr_report.threshold_01 << ")\n";
                std::cerr << "  At 5% FDR: " << fdr_report.genes_05 << " genes"
                          << " (threshold=" << fdr_report.threshold_05 << ")\n";
                std::cerr << "  At 10% FDR: " << fdr_report.genes_10 << " genes"
                          << " (threshold=" << fdr_report.threshold_10 << ")\n";

                if (opts.fdr_threshold >= 0.0f) {
                    size_t genes_passing = 0;
                    for (float s : fdr_target_scores) {
                        if (s >= fdr_score_threshold) genes_passing++;
                    }
                    std::cerr << "  At " << std::setprecision(1) << (opts.fdr_threshold * 100.0f) << "% FDR: "
                              << genes_passing << " genes (threshold=" << std::setprecision(3)
                              << fdr_score_threshold << ")\n";
                }
            }
        }

        // =====================================================================
        // PERCENTILE RANK COMPUTATION
        // Compute damage_signal percentile ranks for all buffered genes
        // This enables contaminant filtering based on relative damage levels
        // =====================================================================
        if (defer_output_for_fdr && !fdr_gene_buffer.empty()) {
            // Collect all damage_signal values
            std::vector<float> all_damage_signals;
            all_damage_signals.reserve(total_genes);
            for (const auto& [id, genes] : fdr_gene_buffer) {
                for (const auto& gene : genes) {
                    all_damage_signals.push_back(gene.damage_signal);
                }
            }

            // Sort for percentile computation
            std::vector<float> sorted_signals = all_damage_signals;
            std::sort(sorted_signals.begin(), sorted_signals.end());

            // Create lookup function for percentile (0-100 scale)
            auto compute_percentile = [&sorted_signals](float value) -> float {
                if (sorted_signals.empty()) return 50.0f;
                auto it = std::lower_bound(sorted_signals.begin(), sorted_signals.end(), value);
                size_t rank = std::distance(sorted_signals.begin(), it);
                return 100.0f * static_cast<float>(rank) / sorted_signals.size();
            };

            // Assign percentile ranks to all genes
            for (auto& [id, genes] : fdr_gene_buffer) {
                for (auto& gene : genes) {
                    gene.damage_pctile = compute_percentile(gene.damage_signal);
                }
            }

            if (opts.verbose) {
                std::cerr << "[Percentile] Computed damage_signal percentiles for "
                          << all_damage_signals.size() << " genes\n";
                if (!sorted_signals.empty()) {
                    std::cerr << "  Range: " << std::setprecision(3)
                              << sorted_signals.front() << " - " << sorted_signals.back()
                              << " | Median: " << sorted_signals[sorted_signals.size()/2] << "\n";
                }
            }
        }

        // Write buffered genes with FDR filtering applied
        size_t fdr_filtered_genes = 0;
        if (defer_output_for_fdr && fdr_score_threshold > 0.0f) {
            std::cerr << "\n[FDR Filtering] Writing genes passing " << std::setprecision(1)
                      << (opts.fdr_threshold * 100.0f) << "% FDR threshold (score >= "
                      << std::setprecision(3) << fdr_score_threshold << ")...\n";

            // Reset histograms for filtered data
            std::fill(damage_prob_hist.begin(), damage_prob_hist.end(), 0);
            std::fill(damage_pct_hist.begin(), damage_pct_hist.end(), 0);
            total_damage_pct = 0.0;

            for (auto& [id, genes] : fdr_gene_buffer) {
                // Filter genes by FDR score threshold
                std::vector<agp::Gene> filtered;
                filtered.reserve(genes.size());
                for (auto& gene : genes) {
                    if (gene.score >= fdr_score_threshold) {
                        filtered.push_back(std::move(gene));
                        // Update histograms for filtered genes
                        int prob_bucket = static_cast<int>(gene.damage_signal * 10);
                        if (prob_bucket >= 10) prob_bucket = 9;
                        if (prob_bucket >= 0) damage_prob_hist[prob_bucket]++;
                        int dmg_bucket = std::min(9, static_cast<int>(gene.damage_score / 10.0));
                        damage_pct_hist[dmg_bucket]++;
                        total_damage_pct += gene.damage_score;
                    }
                }

                if (!filtered.empty()) {
                    writer.write_genes(id, filtered);
                    if (fasta_nt_writer) fasta_nt_writer->write_genes_nucleotide(id, filtered);
                    if (fasta_nt_corr_writer) fasta_nt_corr_writer->write_genes_nucleotide_corrected(id, filtered);
                    if (fasta_aa_writer) fasta_aa_writer->write_genes_protein(id, filtered);
                    if (fasta_corrected_writer) fasta_corrected_writer->write_genes_protein_corrected(id, filtered);
                    fdr_filtered_genes += filtered.size();
                }
            }

            std::cerr << "  Wrote " << fdr_filtered_genes << " genes passing FDR threshold"
                      << " (filtered " << (total_genes - fdr_filtered_genes) << " genes)\n";
            total_genes = fdr_filtered_genes;  // Update for summary statistics
        }
        // Write buffered genes with percentile filtering (when FDR not enabled)
        else if (defer_output_for_fdr && opts.min_damage_pctile > 0.0f) {
            std::cerr << "\n[Percentile Filtering] Writing genes with damage_pctile >= "
                      << std::setprecision(1) << opts.min_damage_pctile << "%...\n";

            // Reset histograms for filtered data
            std::fill(damage_prob_hist.begin(), damage_prob_hist.end(), 0);
            std::fill(damage_pct_hist.begin(), damage_pct_hist.end(), 0);
            total_damage_pct = 0.0;
            size_t pctile_filtered_genes = 0;

            for (auto& [id, genes] : fdr_gene_buffer) {
                std::vector<agp::Gene> filtered;
                filtered.reserve(genes.size());
                for (auto& gene : genes) {
                    if (gene.damage_pctile >= opts.min_damage_pctile) {
                        filtered.push_back(std::move(gene));
                        // Update histograms
                        int prob_bucket = static_cast<int>(gene.damage_signal * 10);
                        if (prob_bucket >= 10) prob_bucket = 9;
                        if (prob_bucket >= 0) damage_prob_hist[prob_bucket]++;
                        int dmg_bucket = std::min(9, static_cast<int>(gene.damage_score / 10.0));
                        damage_pct_hist[dmg_bucket]++;
                        total_damage_pct += gene.damage_score;
                    }
                }

                if (!filtered.empty()) {
                    writer.write_genes(id, filtered);
                    if (fasta_nt_writer) fasta_nt_writer->write_genes_nucleotide(id, filtered);
                    if (fasta_nt_corr_writer) fasta_nt_corr_writer->write_genes_nucleotide_corrected(id, filtered);
                    if (fasta_aa_writer) fasta_aa_writer->write_genes_protein(id, filtered);
                    if (fasta_corrected_writer) fasta_corrected_writer->write_genes_protein_corrected(id, filtered);
                    pctile_filtered_genes += filtered.size();
                }
            }

            std::cerr << "  Wrote " << pctile_filtered_genes << " genes above percentile threshold"
                      << " (filtered " << (total_genes - pctile_filtered_genes) << " potential contaminants)\n";
            total_genes = pctile_filtered_genes;
        }
        // Write buffered genes without filtering (percentile computation only)
        else if (defer_output_for_fdr && !fdr_gene_buffer.empty()) {
            for (auto& [id, genes] : fdr_gene_buffer) {
                writer.write_genes(id, genes);
                if (fasta_nt_writer) fasta_nt_writer->write_genes_nucleotide(id, genes);
                if (fasta_nt_corr_writer) fasta_nt_corr_writer->write_genes_nucleotide_corrected(id, genes);
                if (fasta_aa_writer) fasta_aa_writer->write_genes_protein(id, genes);
                if (fasta_corrected_writer) fasta_corrected_writer->write_genes_protein_corrected(id, genes);
            }
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

            // Pre-frame correction statistics
            size_t preframe_total = total_preframe_corrections.load();
            size_t reads_precorrected = total_reads_precorrected.load();
            if (preframe_total > 0) {
                std::cerr << "\n  [PRE-FRAME] Reads with pre-corrections: " << reads_precorrected
                         << " / " << seq_count << " (" << std::setprecision(1)
                         << (100.0 * reads_precorrected / seq_count) << "%)\n";
                std::cerr << "  [PRE-FRAME] Total DNA corrections: " << preframe_total
                         << " (avg " << std::setprecision(2)
                         << (double)preframe_total / reads_precorrected << " per read)\n";
            }

            // Post-frame correction statistics
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
                std::cerr << "  [CORRECTION] Stop codons restored: " << total_stop_restorations_all << "\n";
                if (total_corrections_reverted > 0) {
                    std::cerr << "  [CORRECTION] Frame-destabilizing corrections reverted: " << total_corrections_reverted << "\n";
                }
                if (total_low_confidence_skipped > 0) {
                    std::cerr << "  [CORRECTION] Low-confidence reads skipped: " << total_low_confidence_skipped << "\n";
                }
            } else {
                std::cerr << "\n  [CORRECTION] No genes corrected (0 / " << total_genes << ")\n";
                if (total_low_confidence_skipped > 0) {
                    std::cerr << "  [CORRECTION] Low-confidence reads skipped: " << total_low_confidence_skipped << "\n";
                }
            }
        }

        std::cerr << "\nDone.\n";

        // Calculate frame-aware wobble ratio from accumulated counts
        // This replaces the frame-unaware values from damage_detection.cpp
        if (opts.use_damage && total_genes > 0) {
            // Calculate T/(T+C) damage rate at each frame-aware codon position
            double pos0_t = wobble_t_count[0].load();
            double pos0_c = wobble_c_count[0].load();
            double pos1_t = wobble_t_count[1].load();
            double pos1_c = wobble_c_count[1].load();
            double pos2_t = wobble_t_count[2].load();
            double pos2_c = wobble_c_count[2].load();

            // Damage = T/(T+C) for each codon position
            // Note: codon_pos 0 = 1st position, 1 = 2nd position, 2 = 3rd (wobble) position
            double pos0_damage = (pos0_t + pos0_c > 100) ? pos0_t / (pos0_t + pos0_c) : 0.0;
            double pos1_damage = (pos1_t + pos1_c > 100) ? pos1_t / (pos1_t + pos1_c) : 0.0;
            double pos2_damage = (pos2_t + pos2_c > 100) ? pos2_t / (pos2_t + pos2_c) : 0.0;

            // Update sample_profile with frame-aware values
            sample_profile.codon_pos1_damage = static_cast<float>(pos0_damage);
            sample_profile.codon_pos2_damage = static_cast<float>(pos1_damage);
            sample_profile.codon_pos3_damage = static_cast<float>(pos2_damage);

            // Calculate frame-aware wobble ratio: pos3 / avg(pos1, pos2)
            // For ancient DNA, we expect pos3 > avg(pos1, pos2) because wobble position
            // tolerates more synonymous substitutions (C->T at wobble often silent)
            double avg_pos12 = (pos0_damage + pos1_damage) / 2.0;
            if (avg_pos12 > 0.001) {
                sample_profile.wobble_ratio = static_cast<float>(pos2_damage / avg_pos12);
            } else {
                sample_profile.wobble_ratio = 1.0f;  // No damage detected
            }

            if (opts.verbose) {
                std::cerr << "\n  [WOBBLE] Frame-aware codon position damage:\n";
                std::cerr << "    Position 1: " << std::setprecision(2) << pos0_damage * 100.0 << "% T/(T+C)\n";
                std::cerr << "    Position 2: " << std::setprecision(2) << pos1_damage * 100.0 << "% T/(T+C)\n";
                std::cerr << "    Position 3 (wobble): " << std::setprecision(2) << pos2_damage * 100.0 << "% T/(T+C)\n";
                std::cerr << "    Wobble ratio: " << std::setprecision(3) << sample_profile.wobble_ratio << "\n";

                // Per-domain terminal damage (meta mode only)
                if (opts.metagenome_mode) {
                    static const char* domain_names[] = {
                        "GTDB", "Fungi", "Protozoa", "Invertebrate",
                        "Plant", "Mammal", "Vertebrate", "Viral"
                    };
                    std::cerr << "\n  [DOMAIN] Per-domain terminal T/(T+C) (5' 3bp):\n";
                    for (size_t d = 0; d < 8; ++d) {
                        size_t t = global_domain_terminal_t[d];
                        size_t c = global_domain_terminal_c[d];
                        size_t n = global_domain_gene_counts[d];
                        if (n > 10) {  // Only report domains with enough genes
                            double damage = (t + c > 0) ? 100.0 * t / (t + c) : 0.0;
                            std::cerr << "    " << std::setw(12) << std::left << domain_names[d]
                                     << " " << std::setw(6) << std::right << n << " genes  "
                                     << std::fixed << std::setprecision(1) << damage << "% T/(T+C)\n";
                        }
                    }
                }
            }
        }

        // Write JSON summary if requested
        if (!opts.summary_file.empty()) {
            std::ofstream summary(opts.summary_file);
            if (summary.is_open()) {
                double mean_dmg_pct = total_genes > 0 ? total_damage_pct / total_genes : 0.0;
                // Use calibrated d_max_combined for damage level (comparable to metaDMG)
                float calibrated_damage = sample_profile.d_max_combined;

                // Determine damage level classification using calibrated values
                std::string damage_level;
                if (sample_profile.is_detection_unreliable()) {
                    // Inverted pattern - reference-free detection failed
                    damage_level = "undetectable";
                } else if (sample_profile.high_asymmetry) {
                    // High asymmetry suggests artifact - mark as unreliable
                    damage_level = "unreliable";
                } else if (calibrated_damage > 0.10f) damage_level = "high";
                else if (calibrated_damage > 0.05f) damage_level = "moderate";
                else if (calibrated_damage > 0.02f) damage_level = "low";
                else damage_level = "minimal";

                summary << "{\n";
                summary << "  \"version\": \"" << AGP_VERSION << "\",\n";
                summary << "  \"input\": \"" << opts.input_file << "\",\n";
                summary << "  \"output\": \"" << opts.output_file << "\",\n";
                summary << "  \"domain\": \"" << (opts.metagenome_mode ? "meta" : opts.domain_name) << "\",\n";
                summary << "  \"sequences\": " << seq_count << ",\n";
                summary << "  \"genes\": " << total_genes << ",\n";
                summary << "  \"processing_time_seconds\": " << std::fixed << std::setprecision(1)
                        << duration.count() / 1000.0 << ",\n";
                summary << "  \"damage\": {\n";
                summary << "    \"detection_enabled\": " << (opts.use_damage ? "true" : "false") << ",\n";
                summary << "    \"aggregate_mode\": " << (opts.aggregate_damage ? "true" : "false") << ",\n";
                summary << "    \"level\": \"" << damage_level << "\",\n";
                summary << "    \"five_prime_percent\": " << std::setprecision(2)
                        << sample_profile.max_damage_5prime * 100.0f << ",\n";
                summary << "    \"three_prime_percent\": " << std::setprecision(2)
                        << sample_profile.max_damage_3prime * 100.0f << ",\n";
                // Calibrated D_max using exponential fit (comparable to metaDMG)
                summary << "    \"d_max\": " << std::setprecision(2)
                        << sample_profile.d_max_combined * 100.0f << ",\n";
                summary << "    \"d_max_5prime\": " << std::setprecision(2)
                        << sample_profile.d_max_5prime * 100.0f << ",\n";
                summary << "    \"d_max_3prime\": " << std::setprecision(2)
                        << sample_profile.d_max_3prime * 100.0f << ",\n";
                summary << "    \"d_max_source\": \"" << sample_profile.d_max_source_str() << "\",\n";
                summary << "    \"asymmetry\": " << std::setprecision(3)
                        << sample_profile.asymmetry << ",\n";
                summary << "    \"high_asymmetry\": "
                        << (sample_profile.high_asymmetry ? "true" : "false") << ",\n";
                // Exponential fit parameters
                summary << "    \"fit_baseline_5prime\": " << std::setprecision(4)
                        << sample_profile.fit_baseline_5prime << ",\n";
                summary << "    \"fit_baseline_3prime\": " << std::setprecision(4)
                        << sample_profile.fit_baseline_3prime << ",\n";
                summary << "    \"fit_amplitude_5prime\": " << std::setprecision(4)
                        << sample_profile.fit_amplitude_5prime << ",\n";
                summary << "    \"fit_amplitude_3prime\": " << std::setprecision(4)
                        << sample_profile.fit_amplitude_3prime << ",\n";
                summary << "    \"fit_rmse_5prime\": " << std::setprecision(4)
                        << sample_profile.fit_rmse_5prime << ",\n";
                summary << "    \"fit_rmse_3prime\": " << std::setprecision(4)
                        << sample_profile.fit_rmse_3prime << ",\n";
                summary << "    \"lambda_5prime\": " << std::setprecision(3)
                        << sample_profile.lambda_5prime << ",\n";
                summary << "    \"lambda_3prime\": " << std::setprecision(3)
                        << sample_profile.lambda_3prime << ",\n";
                // Briggs-like model parameters (single/double-stranded deamination rates)
                summary << "    \"delta_s_5prime\": " << std::setprecision(4)
                        << sample_profile.delta_s_5prime << ",\n";
                summary << "    \"delta_d_5prime\": " << std::setprecision(4)
                        << sample_profile.delta_d_5prime << ",\n";
                summary << "    \"delta_s_3prime\": " << std::setprecision(4)
                        << sample_profile.delta_s_3prime << ",\n";
                summary << "    \"delta_d_3prime\": " << std::setprecision(4)
                        << sample_profile.delta_d_3prime << ",\n";
                summary << "    \"r_squared_5prime\": " << std::setprecision(3)
                        << sample_profile.r_squared_5prime << ",\n";
                summary << "    \"r_squared_3prime\": " << std::setprecision(3)
                        << sample_profile.r_squared_3prime << ",\n";
                summary << "    \"library_type\": \"" << sample_profile.library_type_str() << "\",\n";
                // Terminal inversion detection (statistically robust: shift < -0.01 AND z < -2.0)
                summary << "    \"terminal_shift_5prime\": " << std::setprecision(4)
                        << sample_profile.terminal_shift_5prime << ",\n";
                summary << "    \"terminal_shift_3prime\": " << std::setprecision(4)
                        << sample_profile.terminal_shift_3prime << ",\n";
                summary << "    \"terminal_z_5prime\": " << std::setprecision(2)
                        << sample_profile.terminal_z_5prime << ",\n";
                summary << "    \"terminal_z_3prime\": " << std::setprecision(2)
                        << sample_profile.terminal_z_3prime << ",\n";
                summary << "    \"terminal_inversion\": "
                        << (sample_profile.terminal_inversion ? "true" : "false") << ",\n";
                // Codon-position damage analysis (unique to raw read analysis)
                summary << "    \"codon_pos1_damage\": " << std::setprecision(4)
                        << sample_profile.codon_pos1_damage << ",\n";
                summary << "    \"codon_pos2_damage\": " << std::setprecision(4)
                        << sample_profile.codon_pos2_damage << ",\n";
                summary << "    \"codon_pos3_damage\": " << std::setprecision(4)
                        << sample_profile.codon_pos3_damage << ",\n";
                summary << "    \"wobble_ratio\": " << std::setprecision(3)
                        << sample_profile.wobble_ratio << ",\n";
                summary << "    \"hexamer_damage_llr\": " << std::setprecision(4)
                        << sample_profile.hexamer_damage_llr << ",\n";
                // Likelihood-based model comparison (exponential decay vs constant)
                summary << "    \"decay_llr_5prime\": " << std::setprecision(4)
                        << sample_profile.decay_llr_5prime << ",\n";
                summary << "    \"decay_llr_3prime\": " << std::setprecision(4)
                        << sample_profile.decay_llr_3prime << ",\n";
                // Control channel decay LLR (should be ~0 for real damage)
                summary << "    \"ctrl_decay_llr_5prime\": " << std::setprecision(4)
                        << sample_profile.ctrl_decay_llr_5prime << ",\n";
                summary << "    \"ctrl_decay_llr_3prime\": " << std::setprecision(4)
                        << sample_profile.ctrl_decay_llr_3prime << ",\n";
                // Delta LLR: damage - control (positive = real damage)
                summary << "    \"delta_llr_5prime\": " << std::setprecision(4)
                        << sample_profile.delta_llr_5prime << ",\n";
                summary << "    \"delta_llr_3prime\": " << std::setprecision(4)
                        << sample_profile.delta_llr_3prime << ",\n";
                // Channel divergence (damage channel vs control channel)
                summary << "    \"channel_divergence_5prime\": " << std::setprecision(4)
                        << sample_profile.channel_divergence_5prime << ",\n";
                summary << "    \"channel_divergence_3prime\": " << std::setprecision(4)
                        << sample_profile.channel_divergence_3prime << ",\n";
                // Channel B: stop codon conversion signal (independent validator)
                summary << "    \"stop_conversion_baseline\": " << std::setprecision(4)
                        << sample_profile.stop_conversion_rate_baseline << ",\n";
                summary << "    \"stop_decay_llr_5prime\": " << std::setprecision(2)
                        << sample_profile.stop_decay_llr_5prime << ",\n";
                summary << "    \"stop_amplitude_5prime\": " << std::setprecision(4)
                        << sample_profile.stop_amplitude_5prime << ",\n";
                summary << "    \"channel_b_valid\": " << (sample_profile.channel_b_valid ? "true" : "false") << ",\n";
                summary << "    \"damage_validated\": " << (sample_profile.damage_validated ? "true" : "false") << ",\n";
                summary << "    \"damage_artifact\": " << (sample_profile.damage_artifact ? "true" : "false") << "\n";
                summary << "  },\n";
                summary << "  \"damage_signal_distribution\": [";
                for (int i = 0; i < 10; ++i) {
                    double pct = total_genes > 0 ? 100.0 * damage_prob_hist[i] / total_genes : 0.0;
                    summary << std::setprecision(2) << pct;
                    if (i < 9) summary << ", ";
                }
                summary << "],\n";
                summary << "  \"damage_pct_distribution\": [";
                for (int i = 0; i < 10; ++i) {
                    double pct = total_genes > 0 ? 100.0 * damage_pct_hist[i] / total_genes : 0.0;
                    summary << std::setprecision(2) << pct;
                    if (i < 9) summary << ", ";
                }
                summary << "],\n";
                summary << "  \"corrections\": {\n";
                summary << "    \"genes_corrected\": " << total_genes_corrected << ",\n";
                summary << "    \"total_dna_corrections\": " << total_dna_corrections << ",\n";
                summary << "    \"total_aa_corrections\": " << total_aa_corrections << ",\n";
                summary << "    \"stop_codons_restored\": " << total_stop_restorations_all << ",\n";
                summary << "    \"frame_destabilizing_reverted\": " << total_corrections_reverted << ",\n";
                summary << "    \"low_confidence_skipped\": " << total_low_confidence_skipped << ",\n";
                summary << "    \"preframe_dna_corrections\": " << total_preframe_corrections.load() << ",\n";
                summary << "    \"reads_with_preframe_corrections\": " << total_reads_precorrected.load() << "\n";
                summary << "  }";

                // Add FDR section if computed
                if (opts.compute_fdr && !fdr_target_scores.empty()) {
                    summary << ",\n";
                    summary << "  \"fdr\": {\n";
                    summary << "    \"targets\": " << fdr_report.n_targets << ",\n";
                    summary << "    \"decoys\": " << fdr_report.n_decoys << ",\n";
                    summary << "    \"decoys_per_target\": " << opts.fdr_decoys << ",\n";
                    summary << "    \"seed\": " << opts.fdr_seed << ",\n";
                    summary << "    \"null_mean\": " << std::setprecision(4) << fdr_report.empirical_null_mean << ",\n";
                    summary << "    \"null_sd\": " << std::setprecision(4) << fdr_report.empirical_null_sd << ",\n";
                    summary << "    \"threshold_1pct\": " << std::setprecision(4) << fdr_report.threshold_01 << ",\n";
                    summary << "    \"threshold_5pct\": " << std::setprecision(4) << fdr_report.threshold_05 << ",\n";
                    summary << "    \"threshold_10pct\": " << std::setprecision(4) << fdr_report.threshold_10 << ",\n";
                    summary << "    \"genes_at_1pct\": " << fdr_report.genes_01 << ",\n";
                    summary << "    \"genes_at_5pct\": " << fdr_report.genes_05 << ",\n";
                    summary << "    \"genes_at_10pct\": " << fdr_report.genes_10;
                    if (opts.fdr_threshold >= 0.0f) {
                        size_t genes_at_custom = 0;
                        for (float s : fdr_target_scores) {
                            if (s >= fdr_score_threshold) genes_at_custom++;
                        }
                        summary << ",\n";
                        summary << "    \"custom_threshold\": " << std::setprecision(4) << opts.fdr_threshold << ",\n";
                        summary << "    \"custom_score_threshold\": " << std::setprecision(4) << fdr_score_threshold << ",\n";
                        summary << "    \"genes_at_custom\": " << genes_at_custom << "\n";
                    } else {
                        summary << "\n";
                    }
                    summary << "  }\n";
                } else {
                    summary << "\n";
                }
                summary << "}\n";
                summary.close();
            }
        }

        writer.close();
        return 0;

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
}
