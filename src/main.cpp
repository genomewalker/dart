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
#include "dart/libtaph_damage.hpp"
#include "taph/damage_summary.hpp"
#include "taph/damage_profiler.hpp"
#include "dart/log_utils.hpp"
#include "dart/version.h"
#include <iostream>
#include <fstream>
#include <optional>
#include <chrono>
#include <iomanip>
#include <mutex>
#include <atomic>
#include <vector>
#include <cstring>
#include <unistd.h>
#include <sstream>
#include <zlib.h>

#ifdef _OPENMP
#include <omp.h>
#endif

using Options = dart::cli::Options;

// Output file wrapper: supports both plain and gzip-compressed (.gz) output.
// Uses zlib (already linked) for in-process gz writing — no subprocess needed.
class OutputFile {
    gzFile        gz_    = nullptr;
    std::ofstream plain_;
    bool          is_gz_ = false;

public:
    OutputFile() = default;
    explicit OutputFile(const std::string& path) { open(path); }
    ~OutputFile() { close(); }

    void open(const std::string& path) {
        is_gz_ = path.size() > 3 && path.compare(path.size() - 3, 3, ".gz") == 0;
        if (is_gz_) {
            gz_ = gzopen(path.c_str(), "wb");
        } else {
            plain_.open(path);
        }
    }

    void close() {
        if (gz_) { gzclose(gz_); gz_ = nullptr; }
        if (plain_.is_open()) plain_.close();
    }

    explicit operator bool() const {
        return is_gz_ ? gz_ != nullptr : plain_.good();
    }

    OutputFile& operator<<(const std::string& s) {
        if (is_gz_) gzwrite(gz_, s.data(), static_cast<unsigned>(s.size()));
        else        plain_ << s;
        return *this;
    }
    OutputFile& operator<<(const char* s) {
        if (is_gz_) gzputs(gz_, s);
        else        plain_ << s;
        return *this;
    }
};

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
    std::ostringstream fasta_aa_corrected;
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
        fasta_aa_corrected.str("");
        fasta_aa_corrected.clear();
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

    // Format FASTA entry for Bayesian-corrected protein (ancestor AAs at damage stops)
    void write_fasta_aa_corrected(const std::string& seq_id, const Gene& gene, int gene_num) {
        const std::string& prot = gene.corrected_protein.empty() ? gene.protein : gene.corrected_protein;
        fasta_aa_corrected << ">" << seq_id << "_" << (gene.is_forward ? "+" : "-") << "_" << gene.frame
                           << " " << (gene.start + 1) << ".." << gene.end
                           << " length=" << prot.length() << "aa"
                           << " stop_fix=" << gene.stop_restorations
                           << " corr_score=" << gene.orf_rank
                           << "\n" << prot << "\n";
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

        // Library damage profile (single-engine sidecar). Fail fast: DART's self-profiling
        // diverges from fqdup on deduplicated input. On identical reads the damage FIT
        // agrees (LV7008956662: d_cooc 0.449 vs 0.447, pi 0.272 vs 0.251, d_5p 0.122 vs
        // 0.112 within noise) but the confidence STATE lands LOW_ABUNDANCE where fqdup
        // DETECTs — the divergence is in the eligibility-count / strand-symmetry inputs to
        // the shared libtaph state gate, not the fit or a collapsed denominator. The exact
        // gate input is not yet root-caused, so the canonical fqdup-fit DamageSummary is
        // required unless self-profiling is explicitly opted into for standalone use.
        std::optional<taph::DamageSummary> library_summary;
        if (!opts.library_profile.empty()) {
            library_summary = taph::read_damage_summary(opts.library_profile);
            if (!library_summary) {
                std::cerr << "[FATAL] --library-profile could not be read/parsed: "
                          << opts.library_profile << "\n";
                return 1;
            }
            std::cerr << "[library-profile] loaded " << opts.library_profile
                      << " (pi.state=" << static_cast<int>(library_summary->pi_state)
                      << " pi.point=" << library_summary->pi_point
                      << " d_max_cooc=" << library_summary->d_max_cooc_point
                      << " d_max_5p=" << library_summary->d_max_5prime << ")\n";
        } else if (opts.use_damage && !opts.allow_self_profile) {
            std::cerr << "[FATAL] predict requires --library-profile: DART self-profiling is "
                         "disabled\n        because it is unreliable on deduplicated libraries. "
                         "Pass --allow-self-profile\n        to use the tracked self-profile "
                         "fallback for standalone use.\n";
            return 1;
        }

        // Initialize damage model
        dart::DamageModel damage_model;

        // Open output files — supports plain and .gz paths transparently
        OutputFile gff_file(opts.output_file);
        if (!gff_file) {
            throw std::runtime_error("Cannot open output file: " + opts.output_file);
        }
        // Write GFF3 header
        gff_file << "##gff-version 3\n";

        // Direct file streams for parallel output buffering
        OutputFile fasta_nt_file;
        OutputFile fasta_aa_file;
        OutputFile fasta_aa_corrected_file;
        OutputFile fasta_aa_masked_file;

        // Flags to check if outputs are enabled
        bool write_fasta_nt = !opts.fasta_nt.empty();
        bool write_fasta_aa = !opts.fasta_aa.empty();
        bool write_fasta_aa_corrected = !opts.fasta_aa_corrected.empty();
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
        if (write_fasta_aa_corrected) {
            fasta_aa_corrected_file.open(opts.fasta_aa_corrected);
            if (!fasta_aa_corrected_file) {
                throw std::runtime_error("Cannot open FASTA AA corrected file: " + opts.fasta_aa_corrected);
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

        // Adapter pre-scan (single-engine libtaph): DART must clip adapter stubs before
        // profiling AND before per-read authentication, exactly as fqdup does, or terminal
        // remnants corrupt the 5'/3' deamination counts (LV7008956662: un-clipped d_5p 0.122
        // -> LOW_ABUNDANCE; clipped 0.112 -> DETECTED, matching fqdup to the digit). Hoisted
        // OUT of the aggregate block and run whenever damage is used, so the --no-aggregate
        // path (Pass 2 per-read scoring) clips read termini too. Reuses shared taph primitives.
        taph::AdapterStubs adapter_stubs;
        std::vector<uint32_t> stub5_codes, stub3_codes;
        bool clip_reads = false;
        // §8.1: damage is profiled on prof_in (the merged pre-dedup reads when
        // --profile-input is set), NOT -i (the dedup set, which lands LOW_ABUNDANCE).
        // ORF prediction always runs on -i. Adapter stubs are detected on the same
        // reads we profile (fqdup-faithful), then applied to -i in Pass 2.
        const std::string prof_in = opts.profile_input.empty() ? opts.input_file
                                                               : opts.profile_input;
        if (opts.use_damage) {
            auto scan_reader = std::make_shared<dart::SequenceReader>(prof_in);
            taph::BatchFn scan_fn =
                [scan_reader](std::vector<std::string>& out, size_t max) -> bool {
                    dart::SequenceRecord rec;
                    size_t added = 0;
                    while (added < max && scan_reader->read_next(rec)) {
                        std::string s = dart::SequenceUtils::clean(rec.sequence);
                        if ((int)s.size() >= 30) { out.push_back(std::move(s)); ++added; }
                    }
                    return added > 0;
                };
            try {
                auto prescan = taph::run_prescan(
                    scan_fn, to_sample_library_type(opts.forced_library_type), 30,
                    opts.profile_sample);   // 0 = scan all (§8.1)
                adapter_stubs = std::move(prescan.stubs);
                stub5_codes = taph::encode_stub_codes(adapter_stubs.stubs5);
                stub3_codes = taph::encode_stub_codes(adapter_stubs.stubs3);
                if (adapter_stubs.adapter_clipped || adapter_stubs.adapter3_clipped) {
                    std::cerr << "  Adapter pre-scan: clipping stubs (5': "
                              << (adapter_stubs.stubs5.empty() ? "-" : adapter_stubs.stubs5[0])
                              << ", 3': "
                              << (adapter_stubs.stubs3.empty() ? "-" : adapter_stubs.stubs3[0])
                              << ")\n";
                }
            } catch (const std::exception& e) {
                // FAIL-CLOSED: a thrown prescan means the adapter state is UNKNOWN. Silently
                // scoring raw termini would re-poison the terminal deamination counts feeding
                // authentication (the exact bug the clip fix removes), so abort rather than
                // authenticate on unclipped reads. Explicit override for non-authentication or
                // known-adapterless input: env DART_ALLOW_UNCLIPPED_AUTH=1.
                const char* allow = std::getenv("DART_ALLOW_UNCLIPPED_AUTH");
                if (allow && std::string(allow) == "1") {
                    std::cerr << "  Adapter pre-scan failed (" << e.what()
                              << "); DART_ALLOW_UNCLIPPED_AUTH=1 set — proceeding UNCLIPPED.\n";
                } else {
                    throw std::runtime_error(
                        std::string("Adapter pre-scan failed on an authentication run: ") + e.what() +
                        " — refusing to score raw termini. Set DART_ALLOW_UNCLIPPED_AUTH=1 to override.");
                }
            }
            clip_reads = adapter_stubs.adapter_clipped || adapter_stubs.adapter3_clipped;
        }

        // Pass 1: Damage detection
        if (opts.aggregate_damage && opts.use_damage) {
            auto pass1_start = std::chrono::high_resolution_clock::now();
            std::cerr << "[Pass 1] Scanning for damage patterns...\n";

            taph::ProfileConfig cfg;
            cfg.threads              = num_threads;
            cfg.forced_lib           = to_sample_library_type(opts.forced_library_type);
            cfg.min_read_len         = static_cast<int>(opts.min_length);
            cfg.short_fragment_floor = 30;                  // dart does not lower the floor
            cfg.prescan_reads        = opts.profile_sample; // 0 = all
            cfg.precomputed_stubs    = clip_reads ? &adapter_stubs : nullptr;

            // Reader factory: a fresh reader at the head of prof_in per call, cleaned
            // exactly like Pass 2 so the profiled bytes match the scored bytes. The
            // clip + length gate now live inside profile_reads (single engine).
            auto make_reader = [&]() -> taph::BatchFn {
                auto rdr = std::make_shared<dart::SequenceReader>(prof_in);
                return [rdr](std::vector<std::string>& out, size_t max) -> bool {
                    dart::SequenceRecord rec;
                    size_t added = 0;
                    while (added < max && rdr->read_next(rec)) {
                        out.push_back(dart::SequenceUtils::clean(rec.sequence));
                        ++added;
                    }
                    return added > 0;
                };
            };

            // Side channel: capture 60..200 bp clipped reads for the periodic-codon EM
            // model (feeds Pass-2 gene prediction, NOT the damage gate). Per-thread
            // buffers, flattened in tid order after the pass — no locks.
            const size_t EM_MAX_SEQS = 50000;
            std::vector<std::vector<std::string>> em_by_thread(std::max(1, num_threads));
            taph::PerReadHook em_hook =
                [&](int tid, const std::string& seq, int /*L*/) {
                    if (seq.size() >= 60 && seq.size() <= 200 &&
                        em_by_thread[tid].size() < EM_MAX_SEQS)
                        em_by_thread[tid].push_back(seq);
                };

            taph::ProfileResult pres = taph::profile_reads(make_reader, cfg, em_hook);
            sample_profile = std::move(pres.profile);

            // Single-engine override: the fqdup-fit sidecar is authoritative for the damage
            // gate (pi.state/pi.point, d_max block). Applied after finalize so it supersedes
            // DART's self-derived profile before Pass-1.5 and the per-read scoring loop.
            if (library_summary) {
                library_summary->apply_to(sample_profile);
                std::cerr << "[library-profile] applied sidecar as authoritative damage gate "
                             "(self-profile superseded)\n";
            }

            // Pass 1.5 (d_metamatch reweighting) removed in the profile_reads unification:
            // it re-derived the profile from a second manual pass and never fed the damage
            // gate when a sidecar was present (the production path). The single-engine
            // sidecar override above is authoritative.

            // Train periodic model from the hook-captured EM sequences (flattened in
            // tid order, capped at EM_MAX_SEQS). Feeds Pass-2 gene prediction.
            std::vector<std::string> em_training_seqs;
            em_training_seqs.reserve(EM_MAX_SEQS);
            for (auto& v : em_by_thread) {
                for (auto& s : v) {
                    if (em_training_seqs.size() >= EM_MAX_SEQS) break;
                    em_training_seqs.push_back(std::move(s));
                }
                if (em_training_seqs.size() >= EM_MAX_SEQS) break;
            }
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

        // If a sidecar was supplied but Pass 1 did not run (--no-aggregate or
        // --no-damage), the override above was skipped; apply it here so the AGD
        // writer and the per-read authentication loop below consume the
        // authoritative fqdup profile rather than the default/self one. Without
        // this the --library-profile flag is silently ignored under --no-aggregate.
        if (library_summary && !(opts.aggregate_damage && opts.use_damage)) {
            library_summary->apply_to(sample_profile);
            std::cerr << "[library-profile] applied sidecar (no Pass-1 profiling path)\n";
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
                // Clip adapter stubs before ORF enumeration AND per-read authentication, on
                // BOTH the aggregate and --no-aggregate paths — un-clipped termini otherwise
                // feed adapter remnants into per_read_damage_libtaph (the .agd authentication
                // LLR) and the damage-aware ORF scoring, matching fqdup's universal clip.
                // This single clipped seq is also what the emitted ORF query and .agd derive
                // from: gene.sequence is a substr() of it, and extract_terminal_codons reads
                // gene.sequence — so DO NOT split to keep an unclipped ORF path. Unclipping
                // here re-poisons the per-read .agd terminal codons with adapter bases scored
                // as terminal C->T damage. No frame desync results: six-frame ORF enumeration
                // re-frames on this clipped seq, so every downstream product (ORF query, fasta,
                // .agd terminal codons) derives from the same clipped sequence.
                if (clip_reads) seq = taph::clip_adapters(seq, stub5_codes, stub3_codes, 30);
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
                    per_read_damage_prior,
                    id);

                float damage_pct = opts.use_damage
                    ? dart::FrameSelector::compute_damage_percentage(seq, sample_profile)
                    : 0.0f;

                // Per-read damage authentication via libtaph's reference-free terminal-site
                // LLR (short-read capable, prior-free) instead of the frame_selector AA-window
                // path. Written to the .agd p_damaged_q and consumed by damage-annotate.
                float per_read_damage = opts.use_damage
                    ? dart::per_read_damage_libtaph(seq, sample_profile)
                    : 0.0f;

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

                    // Use rank_score (emit + LLR contrastive bonus) so GFF score field
                    // and FASTA corr_score=0 both reflect the same ranking criterion.
                    float normalized_score = orf.length > 0
                        ? orf.rank_score / static_cast<float>(orf.length)
                        : 0.0f;
                    gene.frame_score = orf.rank_score;
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
                    if (write_fasta_aa_corrected) buf.write_fasta_aa_corrected(id, gene, gene_num);
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
                if (write_fasta_aa_corrected) {
                    fasta_aa_corrected_file << buf.fasta_aa_corrected.str();
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
            summary << "  \"libtaph\": \"" << DART_LIBTAPH_SHA << "\",\n";
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
