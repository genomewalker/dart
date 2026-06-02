// dart profile: Reference-free ancient DNA damage profiling
//
// Usage: dart profile <input.fq> [-o output.json] [-v]
// Alias:  dart damage
//
// Mirrors fqdup's profiling pipeline:
//   Pre-scan  – first 100k reads, adapter detection via libtaph
//   Pass 1    – full scan with adapter clipping, parallel aggregation
//   Pass 2    – optional damage-weighted re-scan for d_metamatch

#include "subcommand.hpp"
#include "args.hpp"
#include "dart/sequence_io.hpp"
#include "dart/frame_selector.hpp"
#include "dart/hexamer_tables.hpp"
#include "dart/log_utils.hpp"
#include "dart/version.h"
#include <taph/damage_profiler.hpp>
#include <taph/frame_selector_decl.hpp>
#include <iostream>
#include <fstream>
#include <cstring>
#include <iomanip>
#include <chrono>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace dart {
namespace cli {

// Build a BatchFn backed by a single SequenceReader (reads once through the file).
static taph::BatchFn make_batch_fn(const std::string& path, int min_len = 30) {
    auto reader = std::make_shared<SequenceReader>(path);
    return [reader, min_len](std::vector<std::string>& out, size_t max) -> bool {
        SequenceRecord rec;
        size_t added = 0;
        while (added < max && reader->read_next(rec)) {
            std::string seq = SequenceUtils::clean(rec.sequence);
            if ((int)seq.size() >= min_len) {
                out.push_back(std::move(seq));
                ++added;
            }
        }
        return added > 0;
    };
}

int cmd_profile(int argc, char* argv[]) {
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
            std::cout << "Usage: dart profile <input.fq> [options]\n\n";
            std::cout << "Reference-free ancient DNA damage profiling.\n\n";
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
            std::cerr << "Run 'dart profile --help' for usage.\n";
            return 1;
        }
    }

    if (input_file.empty()) {
        std::cerr << "Error: No input file specified.\n";
        std::cerr << "Run 'dart profile --help' for usage.\n";
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

    Domain active_domain = parse_domain(domain_str);
    set_active_domain(active_domain);

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
        std::cerr << "dart profile v" << DART_VERSION << "\n";
        std::cerr << "Input: " << input_file << "\n";
        std::cerr << "Domain: " << domain_name(active_domain) << "\n";
        std::cerr << "Threads: " << num_threads << "\n\n";
    }

    // =========================================================================
    // PRE-SCAN: adapter detection (first 100k reads, single-threaded)
    // =========================================================================
    taph::AdapterStubs stubs;
    try {
        if (verbose) std::cerr << "Pre-scan: detecting adapters...\r" << std::flush;
        auto prescan = taph::run_prescan(
            make_batch_fn(input_file),
            forced_library, 30, 100000);
        stubs = std::move(prescan.stubs);
        if (verbose) {
            if (stubs.adapter_clipped) {
                std::cerr << "Pre-scan: adapter detected (5': " << stubs.stubs5[0]
                          << (stubs.adapter3_clipped ? ", 3': " + stubs.stubs3[0] : "") << ")\n";
            } else {
                std::cerr << "Pre-scan: no adapter detected          \n";
            }
        }
    } catch (const std::exception& e) {
        if (verbose) std::cerr << "Pre-scan warning: " << e.what() << "\n";
    }

    // Encode stub hexamer codes for fast per-read lookup
    auto stub5_codes = taph::encode_stub_codes(stubs.stubs5);
    auto stub3_codes = taph::encode_stub_codes(stubs.stubs3);

    // =========================================================================
    // PASS 1: full scan with adapter clipping
    // =========================================================================
    auto pass1_start = std::chrono::steady_clock::now();
    SampleDamageProfile profile;
    profile.forced_library_type = forced_library;
    size_t total_reads = 0;
    taph::HexEndAsymmetry hex_asym{};

    try {
        SequenceReader reader(input_file);
        const size_t BATCH_SIZE = 50000;
        std::vector<SequenceRecord> batch;
        batch.reserve(BATCH_SIZE);

        std::vector<SampleDamageProfile> thread_profiles(num_threads);

        while (true) {
            batch.clear();
            SequenceRecord record;
            while (batch.size() < BATCH_SIZE && reader.read_next(record)) {
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
                if (stubs.adapter_clipped || stubs.adapter3_clipped) {
                    seq = taph::clip_adapters(seq, stub5_codes, stub3_codes, 30);
                }
                if (seq.size() >= 30) {
                    taph::FrameSelector::update_sample_profile(thread_profiles[tid], seq);
                }
            }

            total_reads += batch.size();
            if (verbose && total_reads % 100000 == 0) {
                std::cerr << "Processed " << total_reads << " sequences...\r" << std::flush;
            }
        }

        if (verbose) std::cerr << "Processed " << total_reads << " sequences.    \n";

        for (int t = 0; t < num_threads; ++t) {
            taph::FrameSelector::merge_sample_profiles(profile, thread_profiles[t]);
        }
        taph::FrameSelector::finalize_sample_profile(profile);

        // Post-clip hex enrichment (matches fqdup's compute_hex_enriched_5prime call)
        stubs.top_enriched     = taph::compute_hex_enriched_5prime(profile);
        stubs.flag_hex_artifact = !stubs.top_enriched.empty()
                                   && stubs.top_enriched[0].log2fc > 1.5f
                                   && !stubs.top_enriched[0].damage_consistent;
        stubs.top_enriched_3prime = taph::compute_hex_enriched_3prime(profile);
        hex_asym = taph::compute_hex_end_asymmetry(
            profile, stubs.top_enriched, stubs.top_enriched_3prime);

        if (verbose) {
            auto pass1_end = std::chrono::steady_clock::now();
            std::cerr << "Pass 1 runtime: "
                      << dart::log_utils::format_elapsed(pass1_start, pass1_end) << "\n";
        }
    } catch (const std::exception& e) {
        std::cerr << "Error reading input: " << e.what() << "\n";
        return 1;
    }

    // =========================================================================
    // PASS 2: damage-weighted re-scan for d_metamatch (optional)
    // =========================================================================
    SampleDamageProfile weighted_profile;
    weighted_profile.forced_library_type = forced_library;
    float d_metamatch_filtered = 0.0f;
    size_t n_high_damage_reads = 0;

    bool needs_pass2 = false;
    if (profile.damage_validated && profile.d_max_combined > 0.15f) {
        if (profile.inverted_pattern_5prime || profile.inverted_pattern_3prime) {
            needs_pass2 = true;
        } else if (profile.position_0_artifact_5prime) {
            needs_pass2 = true;
        } else {
            float gap = std::abs(profile.d_max_from_channel_b - profile.d_max_combined);
            if (gap > 0.05f) needs_pass2 = true;
        }
    }

    if (needs_pass2) {
        auto pass2_start = std::chrono::steady_clock::now();
        if (verbose) {
            std::cerr << "\nPass 2: Computing damage-weighted profile...\n";
            if (profile.inverted_pattern_5prime || profile.inverted_pattern_3prime)
                std::cerr << "  Reason: inverted pattern detected\n";
            else if (profile.position_0_artifact_5prime)
                std::cerr << "  Reason: position-0 artifact detected\n";
            else
                std::cerr << "  Reason: gap between d_max and Channel B > 5%\n";
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
                    if (stubs.adapter_clipped || stubs.adapter3_clipped) {
                        seq = taph::clip_adapters(seq, stub5_codes, stub3_codes, 30);
                    }
                    if (seq.size() >= 30) {
                        float p_damaged = FrameSelector::compute_per_read_damage_prior(seq, profile);
                        if (p_damaged > 0.01f) {
                            taph::FrameSelector::update_sample_profile_weighted(
                                thread_profiles2[tid], seq, p_damaged);
                        }
                        if (p_damaged > 0.6f) thread_high_counts[tid]++;
                    }
                }

                pass2_reads += batch.size();
                if (verbose && pass2_reads % 500000 == 0) {
                    std::cerr << "Pass 2: " << pass2_reads << " reads...\r" << std::flush;
                }
            }

            for (int t = 0; t < num_threads; ++t) {
                taph::FrameSelector::merge_sample_profiles(weighted_profile, thread_profiles2[t]);
                n_high_damage_reads += thread_high_counts[t];
            }
            taph::FrameSelector::finalize_sample_profile(weighted_profile);
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
                          << dart::log_utils::format_elapsed(pass2_start, pass2_end) << "\n";
            }
        } catch (const std::exception& e) {
            if (verbose) std::cerr << "Warning: Pass 2 failed: " << e.what() << "\n";
            d_metamatch_filtered = profile.d_metamatch;
        }

        profile.d_metamatch = d_metamatch_filtered;
    } else if (verbose && profile.damage_validated && profile.d_max_combined > 0.05f) {
        std::cerr << "\nPass 2 skipped: d_max and Channel B are within 5%\n";
        std::cerr << "Using Channel B-anchored d_metamatch: " << std::fixed
                  << std::setprecision(1) << (profile.d_metamatch * 100.0f) << "%\n";
    }

    // =========================================================================
    // OUTPUT
    // =========================================================================
    auto damage_level_str = [](float d_max) -> const char* {
        if (d_max >= 0.10f) return "high";
        if (d_max >= 0.05f) return "moderate";
        if (d_max >= 0.02f) return "low";
        return "undetectable";
    };

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

    auto boolstr = [](bool b) { return b ? "true" : "false"; };
    *out << "{\n";
    *out << "  \"version\": \"" << DART_VERSION << "\",\n";
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
    *out << "    \"channel_b_valid\": " << boolstr(profile.channel_b_valid) << ",\n";
    *out << "    \"channel_b_llr\": " << std::fixed << std::setprecision(2) << profile.stop_decay_llr_5prime << ",\n";
    *out << "    \"channel_b3_llr\": " << std::fixed << std::setprecision(2) << profile.stop_decay_llr_3prime << ",\n";
    *out << "    \"damage_validated\": " << boolstr(profile.damage_validated) << ",\n";
    *out << "    \"position_0_artifact_5prime\": " << boolstr(profile.position_0_artifact_5prime) << ",\n";
    *out << "    \"terminal_shift_5prime\": " << std::fixed << std::setprecision(4) << profile.terminal_shift_5prime << ",\n";
    *out << "    \"inverted_pattern_5prime\": " << boolstr(profile.inverted_pattern_5prime) << ",\n";
    *out << "    \"inverted_pattern_3prime\": " << boolstr(profile.inverted_pattern_3prime) << ",\n";
    *out << "    \"d_max_from_channel_b\": " << std::fixed << std::setprecision(2) << (profile.d_max_from_channel_b * 100.0f) << ",\n";
    *out << "    \"channel_c_valid\": " << boolstr(profile.channel_c_valid) << ",\n";
    *out << "    \"ox_d_max\": " << std::fixed << std::setprecision(2) << (profile.ox_d_max * 100.0f) << ",\n";
    *out << "    \"ox_damage_detected\": " << boolstr(profile.ox_damage_detected) << ",\n";
    *out << "    \"depurination_detected\": " << boolstr(profile.depurination_detected) << ",\n";
    *out << "    \"library_type\": \"" << profile.library_type_str() << "\",\n";
    *out << "    \"d5_profile_flat\": " << boolstr(profile.d5_profile_flat) << ",\n";
    *out << "    \"d3_profile_flat\": " << boolstr(profile.d3_profile_flat) << ",\n";
    *out << "    \"d5_blunting_suspected\": " << boolstr(profile.d5_blunting_suspected) << ",\n";
    *out << std::fixed << std::setprecision(4);
    *out << "    \"d5_max_rate_pos0_4\": " << profile.d5_max_rate_pos0_4 << ",\n";
    *out << "    \"d3_max_rate_pos0_4\": " << profile.d3_max_rate_pos0_4 << ",\n";
    // hex_end_asymmetry is unreliable when position-0 hexamers are corrupted by adapter ligation
    if (!profile.position_0_artifact_5prime && !std::isnan(hex_asym.rc_excess_jsd)) {
        *out << "    \"hex_end_asymmetry\": {\n";
        *out << "      \"rc_excess_jsd\": " << hex_asym.rc_excess_jsd << ",\n";
        *out << "      \"fwd_excess_jsd\": " << hex_asym.fwd_excess_jsd << ",\n";
        *out << "      \"status\": \"" << hex_asym.status << "\"\n";
        *out << "    }\n";
    } else {
        *out << "    \"hex_end_asymmetry\": null\n";
    }
    *out << "  }\n";
    *out << "}\n";

    if (verbose) {
        std::cerr << "\nDamage profile:\n";
        std::cerr << "  Level: " << damage_level_str(d_max) << "\n";
        std::cerr << "  D_max: " << std::fixed << std::setprecision(1) << (d_max * 100.0f) << "%\n";
        std::cerr << "  Channel B LLR:  " << std::fixed << std::setprecision(1)
                  << profile.stop_decay_llr_5prime << " (5' C→T stop conversion)\n";
        std::cerr << "  Channel B3 LLR: " << std::fixed << std::setprecision(1)
                  << profile.stop_decay_llr_3prime << " (3' G→A stop conversion)\n";
        if (profile.channel_c_valid) {
            std::cerr << "  Channel C: ox_d_max=" << std::fixed << std::setprecision(1)
                      << (profile.ox_d_max * 100.0f)
                      << "% uniformity=" << std::setprecision(2) << profile.ox_uniformity_ratio
                      << (profile.ox_damage_detected ? " DETECTED" : "") << " (G→T oxidation)\n";
        }
        if (profile.depurination_detected) {
            std::cerr << "  Channel E: depurination detected (purine enrichment at 5' termini)\n";
        }
        std::cerr << "  Validated: " << (profile.damage_validated ? "yes" : "no") << "\n";
        if (profile.d5_blunting_suspected) {
            std::cerr << "  5' blunting suspected (flat 5', real 3' decay)\n";
        } else if (profile.d5_profile_flat) {
            std::cerr << "  5' profile flat (no detectable terminal signal)\n";
        }
        if (!std::isnan(hex_asym.rc_excess_jsd)) {
            std::cerr << "  Hex end asymmetry: rc_jsd=" << std::fixed << std::setprecision(3)
                      << hex_asym.rc_excess_jsd << " (" << hex_asym.status << ")\n";
        }
        auto run_end = std::chrono::steady_clock::now();
        std::cerr << "  Runtime: " << dart::log_utils::format_elapsed(run_start, run_end) << "\n";
    }

    return 0;
}

// Register subcommand
namespace {
    struct ProfileRegistrar {
        ProfileRegistrar() {
            SubcommandRegistry::instance().register_command(
                "profile",
                "Reference-free ancient DNA damage profiling",
                cmd_profile, 10);
            SubcommandRegistry::instance().register_command(
                "damage",
                "Alias for 'profile'",
                cmd_profile, 11);
        }
    } profile_registrar;
}

}  // namespace cli
}  // namespace dart
