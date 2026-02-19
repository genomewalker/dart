// agp sample-damage: Quick sample-level damage profiling
//
// Usage: agp sample-damage <input.fq> [-o output.json] [-v]
//
// Runs Pass 1 only to compute sample-level damage profile.
// Faster than full prediction when you only need damage metrics.

#include "subcommand.hpp"
#include "args.hpp"
#include "agp/sequence_io.hpp"
#include "agp/frame_selector.hpp"
#include "agp/hexamer_tables.hpp"
#include "agp/log_utils.hpp"
#include "agp/version.h"
#include <iostream>
#include <fstream>
#include <cstring>
#include <iomanip>
#include <chrono>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace agp {
namespace cli {

int cmd_sample_damage(int argc, char* argv[]) {
    auto run_start = std::chrono::steady_clock::now();
    std::string input_file;
    std::string output_file;
    std::string domain_str = "gtdb";
    std::string library_type_str = "auto";
    bool verbose = false;
    int num_threads = 0;

    for (int i = 1; i < argc; ++i) {
        if (strcmp(argv[i], "-o") == 0 && i + 1 < argc) {
            output_file = argv[++i];
        } else if ((strcmp(argv[i], "-d") == 0 || strcmp(argv[i], "--domain") == 0) && i + 1 < argc) {
            domain_str = argv[++i];
        } else if ((strcmp(argv[i], "-t") == 0 || strcmp(argv[i], "--threads") == 0) && i + 1 < argc) {
            num_threads = std::stoi(argv[++i]);
        } else if (strcmp(argv[i], "--library-type") == 0 && i + 1 < argc) {
            library_type_str = argv[++i];
        } else if (strcmp(argv[i], "-v") == 0 || strcmp(argv[i], "--verbose") == 0) {
            verbose = true;
        } else if (strcmp(argv[i], "--help") == 0 || strcmp(argv[i], "-h") == 0) {
            std::cout << "Usage: agp sample-damage <input.fq> [options]\n\n";
            std::cout << "Quick sample-level damage profiling (Pass 1 only).\n\n";
            std::cout << "Options:\n";
            std::cout << "  -o, --output FILE       Output JSON file (default: stdout)\n";
            std::cout << "  -d, --domain DOMAIN     Domain for hexamer scoring (default: gtdb)\n";
            std::cout << "  --library-type TYPE     Library type: ds, ss, or auto (default: auto)\n";
            std::cout << "  -t, --threads N         Number of threads (default: auto)\n";
            std::cout << "  -v, --verbose           Verbose output\n";
            return 0;
        } else if (argv[i][0] != '-') {
            input_file = argv[i];
        } else {
            std::cerr << "Unknown option: " << argv[i] << "\n";
            std::cerr << "Run 'agp sample-damage --help' for usage.\n";
            return 1;
        }
    }

    if (input_file.empty()) {
        std::cerr << "Error: No input file specified.\n";
        std::cerr << "Run 'agp sample-damage --help' for usage.\n";
        return 1;
    }

#ifdef _OPENMP
    if (num_threads == 0) {
        num_threads = omp_get_max_threads();
    }
    omp_set_num_threads(num_threads);
#else
    num_threads = 1;
#endif

    // Set domain
    Domain active_domain = parse_domain(domain_str);
    set_active_domain(active_domain);

    // Parse library type
    SampleDamageProfile::LibraryType forced_library = SampleDamageProfile::LibraryType::UNKNOWN;
    if (library_type_str == "ds") {
        forced_library = SampleDamageProfile::LibraryType::DOUBLE_STRANDED;
    } else if (library_type_str == "ss") {
        forced_library = SampleDamageProfile::LibraryType::SINGLE_STRANDED;
    } else if (library_type_str != "auto") {
        std::cerr << "Error: Unknown library type '" << library_type_str << "'. Use ds, ss, or auto.\n";
        return 1;
    }

    if (verbose) {
        std::cerr << "Sample damage profiling v" << AGP_VERSION << "\n";
        std::cerr << "Input: " << input_file << "\n";
        std::cerr << "Domain: " << domain_name(active_domain) << "\n";
        std::cerr << "Threads: " << num_threads << "\n\n";
    }

    // Run Pass 1: damage detection
    auto pass1_start = std::chrono::steady_clock::now();
    SampleDamageProfile profile;
    profile.forced_library_type = forced_library;
    size_t total_reads = 0;

    try {
        SequenceReader reader(input_file);
        const size_t BATCH_SIZE = 50000;
        std::vector<SequenceRecord> batch;
        batch.reserve(BATCH_SIZE);

        // Thread-local profiles for parallel aggregation
        std::vector<SampleDamageProfile> thread_profiles(num_threads);

        while (true) {
            // Read a batch
            batch.clear();
            SequenceRecord record;
            while (batch.size() < BATCH_SIZE && reader.read_next(record)) {
                batch.push_back(std::move(record));
            }
            if (batch.empty()) break;

            // Process batch in parallel
            #pragma omp parallel for schedule(static)
            for (size_t i = 0; i < batch.size(); ++i) {
                int tid = 0;
                #ifdef _OPENMP
                tid = omp_get_thread_num();
                #endif
                std::string seq = SequenceUtils::clean(batch[i].sequence);
                if (seq.length() >= 30) {
                    FrameSelector::update_sample_profile(thread_profiles[tid], seq);
                }
            }

            total_reads += batch.size();

            if (verbose && total_reads % 100000 == 0) {
                std::cerr << "Processed " << total_reads << " sequences...\r" << std::flush;
            }
        }

        if (verbose) {
            std::cerr << "Processed " << total_reads << " sequences.    \n";
        }

        // Merge thread-local profiles
        for (int t = 0; t < num_threads; ++t) {
            FrameSelector::merge_sample_profiles(profile, thread_profiles[t]);
        }

        // Finalize profile (compute decay rates, d_max, etc.)
        FrameSelector::finalize_sample_profile(profile);
        if (verbose) {
            auto pass1_end = std::chrono::steady_clock::now();
            std::cerr << "Pass 1 runtime: "
                      << agp::log_utils::format_elapsed(pass1_start, pass1_end) << "\n";
        }

    } catch (const std::exception& e) {
        std::cerr << "Error reading input: " << e.what() << "\n";
        return 1;
    }

    // =========================================================================
    // PASS 2: Damage-weighted profile for d_metamatch
    //
    // Re-read the file and weight each read's contribution by P(damaged).
    // This mimics metaDMG's selection bias: reads showing clear damage patterns
    // are likely ancient and would align to references.
    //
    // Only run if Pass 2 would meaningfully improve the estimate:
    // 1. Damage was validated in Pass 1 (otherwise d_metamatch = 0 anyway)
    // 2. d_max > 5% (meaningful damage to filter by)
    // 3. At least one of:
    //    a) Inverted pattern detected (Channel A unreliable)
    //    b) Position-0 artifact (contamination likely)
    //    c) Large gap between d_max and Channel B (selection bias issue)
    // =========================================================================
    SampleDamageProfile weighted_profile;
    weighted_profile.forced_library_type = forced_library;
    float d_metamatch_filtered = 0.0f;
    size_t n_high_damage_reads = 0;

    // Check if Pass 2 would add value
    // Pass 2 only helps for HIGH-damage samples (>15%) where damaged reads
    // are distinguishable from undamaged. For low-damage, weighting dilutes signal.
    bool needs_pass2 = false;
    if (profile.damage_validated && profile.d_max_combined > 0.15f) {
        // For high-damage samples, run Pass 2 if there are reliability issues
        if (profile.inverted_pattern_5prime || profile.inverted_pattern_3prime) {
            needs_pass2 = true;  // Channel A unreliable
        } else if (profile.position_0_artifact_5prime) {
            needs_pass2 = true;  // Contamination likely
        } else {
            // Run if there's a meaningful gap between d_max and Channel B
            float gap = std::abs(profile.d_max_from_channel_b - profile.d_max_combined);
            if (gap > 0.05f) {  // >5% absolute difference
                needs_pass2 = true;
            }
        }
    }

    if (needs_pass2) {
        auto pass2_start = std::chrono::steady_clock::now();
        if (verbose) {
            std::cerr << "\nPass 2: Computing damage-weighted profile...\n";
            if (profile.inverted_pattern_5prime || profile.inverted_pattern_3prime) {
                std::cerr << "  Reason: inverted pattern detected\n";
            } else if (profile.position_0_artifact_5prime) {
                std::cerr << "  Reason: position-0 artifact detected\n";
            } else {
                std::cerr << "  Reason: gap between d_max and Channel B > 5%\n";
            }
        }

        try {
            SequenceReader reader2(input_file);
            const size_t BATCH_SIZE = 50000;
            std::vector<SequenceRecord> batch;
            batch.reserve(BATCH_SIZE);

            std::vector<SampleDamageProfile> thread_profiles2(num_threads);
            std::vector<size_t> thread_high_counts(num_threads, 0);
            size_t pass2_reads = 0;

            while (true) {
                batch.clear();
                SequenceRecord record;
                while (batch.size() < BATCH_SIZE && reader2.read_next(record)) {
                    batch.push_back(std::move(record));
                }
                if (batch.empty()) break;

                #pragma omp parallel for schedule(static)
                for (size_t i = 0; i < batch.size(); ++i) {
                    int tid = 0;
                    #ifdef _OPENMP
                    tid = omp_get_thread_num();
                    #endif
                    std::string seq = SequenceUtils::clean(batch[i].sequence);
                    if (seq.length() >= 30) {
                        // Compute per-read damage probability using Pass 1 profile
                        float p_damaged = FrameSelector::compute_per_read_damage_prior(seq, profile);

                        // Use soft weighting: contribution proportional to P(damaged)
                        // This gives ancient-looking reads higher influence
                        if (p_damaged > 0.01f) {
                            FrameSelector::update_sample_profile_weighted(
                                thread_profiles2[tid], seq, p_damaged);
                        }

                        // Count high-confidence damaged reads
                        if (p_damaged > 0.6f) {
                            thread_high_counts[tid]++;
                        }
                    }
                }

                pass2_reads += batch.size();
                if (verbose && pass2_reads % 500000 == 0) {
                    std::cerr << "Pass 2: " << pass2_reads << " reads...\r" << std::flush;
                }
            }

            // Merge thread-local profiles
            for (int t = 0; t < num_threads; ++t) {
                FrameSelector::merge_sample_profiles(weighted_profile, thread_profiles2[t]);
                n_high_damage_reads += thread_high_counts[t];
            }

            // Finalize weighted profile
            FrameSelector::finalize_sample_profile(weighted_profile);

            // Use weighted profile's d_max as d_metamatch
            d_metamatch_filtered = weighted_profile.d_max_combined;

            if (verbose) {
                std::cerr << "Pass 2: " << pass2_reads << " reads processed.    \n";
                std::cerr << "High-confidence damaged reads: " << n_high_damage_reads
                          << " (" << std::fixed << std::setprecision(1)
                          << (100.0f * n_high_damage_reads / pass2_reads) << "%)\n";
                std::cerr << "Weighted d_max: " << std::fixed << std::setprecision(1)
                          << (d_metamatch_filtered * 100.0f) << "%\n";
                auto pass2_end = std::chrono::steady_clock::now();
                std::cerr << "Pass 2 runtime: "
                          << agp::log_utils::format_elapsed(pass2_start, pass2_end) << "\n";
            }

        } catch (const std::exception& e) {
            if (verbose) {
                std::cerr << "Warning: Pass 2 failed: " << e.what() << "\n";
            }
            // Fall back to Pass 1 d_metamatch
            d_metamatch_filtered = profile.d_metamatch;
        }

        // Update profile's d_metamatch with the filtered estimate
        profile.d_metamatch = d_metamatch_filtered;
    } else if (verbose && profile.damage_validated && profile.d_max_combined > 0.05f) {
        // Explain why Pass 2 was skipped
        std::cerr << "\nPass 2 skipped: d_max and Channel B are within 5%\n";
        std::cerr << "Using Channel B-anchored d_metamatch: " << std::fixed
                  << std::setprecision(1) << (profile.d_metamatch * 100.0f) << "%\n";
    }

    // Determine damage level string
    auto damage_level_str = [](float d_max) -> const char* {
        if (d_max >= 0.10f) return "high";
        if (d_max >= 0.05f) return "moderate";
        if (d_max >= 0.02f) return "low";
        return "undetectable";
    };

    // Output JSON
    std::ostream* out = &std::cout;
    std::ofstream file_out;
    if (!output_file.empty()) {
        file_out.open(output_file);
        if (!file_out.is_open()) {
            std::cerr << "Error: Cannot open output file: " << output_file << "\n";
            return 1;
        }
        out = &file_out;
    }

    float d_max = profile.d_max_combined;

    *out << "{\n";
    *out << "  \"version\": \"" << AGP_VERSION << "\",\n";
    *out << "  \"input\": \"" << input_file << "\",\n";
    *out << "  \"domain\": \"" << domain_str << "\",\n";
    *out << "  \"sequences\": " << total_reads << ",\n";
    *out << "  \"damage\": {\n";
    *out << "    \"detection_enabled\": true,\n";
    *out << "    \"level\": \"" << damage_level_str(d_max) << "\",\n";
    *out << "    \"d_max\": " << std::fixed << std::setprecision(2) << (d_max * 100.0f) << ",\n";
    *out << "    \"d_max_5prime\": " << std::fixed << std::setprecision(2) << (profile.d_max_5prime * 100.0f) << ",\n";
    *out << "    \"d_max_3prime\": " << std::fixed << std::setprecision(2) << (profile.d_max_3prime * 100.0f) << ",\n";
    *out << "    \"delta_s_5prime\": " << std::fixed << std::setprecision(4) << profile.delta_s_5prime << ",\n";
    *out << "    \"delta_s_3prime\": " << std::fixed << std::setprecision(4) << profile.delta_s_3prime << ",\n";
    *out << "    \"lambda_5prime\": " << std::fixed << std::setprecision(3) << profile.lambda_5prime << ",\n";
    *out << "    \"lambda_3prime\": " << std::fixed << std::setprecision(3) << profile.lambda_3prime << ",\n";
    *out << "    \"channel_b_valid\": " << (profile.channel_b_valid ? "true" : "false") << ",\n";
    *out << "    \"channel_b_llr\": " << std::fixed << std::setprecision(2) << profile.stop_decay_llr_5prime << ",\n";
    *out << "    \"damage_validated\": " << (profile.damage_validated ? "true" : "false") << ",\n";
    *out << "    \"position_0_artifact_5prime\": " << (profile.position_0_artifact_5prime ? "true" : "false") << ",\n";
    *out << "    \"terminal_shift_5prime\": " << std::fixed << std::setprecision(4) << profile.terminal_shift_5prime << ",\n";
    *out << "    \"inverted_pattern_5prime\": " << (profile.inverted_pattern_5prime ? "true" : "false") << ",\n";
    *out << "    \"inverted_pattern_3prime\": " << (profile.inverted_pattern_3prime ? "true" : "false") << ",\n";
    *out << "    \"d_max_from_channel_b\": " << std::fixed << std::setprecision(2) << (profile.d_max_from_channel_b * 100.0f) << ",\n";
    *out << "    \"d_metamatch\": " << std::fixed << std::setprecision(2) << (profile.d_metamatch * 100.0f) << ",\n";
    *out << "    \"metamatch_high_damage_reads\": " << n_high_damage_reads << ",\n";
    *out << "    \"metamatch_high_damage_pct\": " << std::fixed << std::setprecision(2)
         << (total_reads > 0 ? 100.0f * n_high_damage_reads / total_reads : 0.0f) << ",\n";
    *out << "    \"library_type\": \"" << profile.library_type_str() << "\"\n";
    *out << "  }\n";
    *out << "}\n";

    if (verbose) {
        std::cerr << "\nDamage profile:\n";
        std::cerr << "  Level: " << damage_level_str(d_max) << "\n";
        std::cerr << "  D_max: " << std::fixed << std::setprecision(1) << (d_max * 100.0f) << "%\n";
        std::cerr << "  Channel B LLR: " << std::fixed << std::setprecision(1) << profile.stop_decay_llr_5prime << "\n";
        std::cerr << "  Validated: " << (profile.damage_validated ? "yes" : "no") << "\n";
        auto run_end = std::chrono::steady_clock::now();
        std::cerr << "  Runtime: " << agp::log_utils::format_elapsed(run_start, run_end) << "\n";
    }

    return 0;
}

// Register subcommand
namespace {
    struct SampleDamageRegistrar {
        SampleDamageRegistrar() {
            SubcommandRegistry::instance().register_command(
                "sample-damage",
                "Quick sample-level damage profiling",
                cmd_sample_damage, 10);
        }
    } sample_damage_registrar;
}

}  // namespace cli
}  // namespace agp
