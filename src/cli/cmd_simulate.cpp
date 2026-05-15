// dart simulate — generate synthetic ancient DNA reads with controlled damage
//
// Embeds random ORFs in reads, applies the Briggs (2007) mechanistic damage
// model (geometric overhangs + nick-position DS interior), and writes a
// ground-truth TSV so downstream scoring can be validated exactly.
//
// Briggs model parameters (typical ancient DNA values in parentheses):
//   nv       — nick probability in DS interior     (0.024)
//   lambda   — geometric overhang length parameter (0.36)
//   delta-s  — SS overhang deamination rate        (0.68)
//   delta-d  — DS interior deamination rate        (0.0097)

#include "subcommand.hpp"
#include "dart/codon_tables.hpp"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iostream>
#include <memory>
#include <random>
#include <string>
#include <vector>

namespace dart {
namespace cli {

namespace {

// ── Sequence primitives ───────────────────────────────────────────────────────

static bool is_stop(char c0, char c1, char c2) {
    if (c0 != 'T') return false;
    return (c1 == 'A' && (c2 == 'A' || c2 == 'G')) ||
           (c1 == 'G' &&  c2 == 'A');
}

static char rand_base(std::mt19937& rng) {
    static const char B[4] = {'A', 'C', 'G', 'T'};
    return B[std::uniform_int_distribution<int>(0, 3)(rng)];
}

static std::string generate_orf(std::mt19937& rng, size_t n_aa) {
    std::string orf;
    orf.reserve(n_aa * 3);
    orf += "ATG";
    char c[3];
    for (size_t i = 1; i < n_aa; ++i) {
        do { c[0] = rand_base(rng); c[1] = rand_base(rng); c[2] = rand_base(rng); }
        while (is_stop(c[0], c[1], c[2]));
        orf.append(c, 3);
    }
    return orf;
}

static std::string translate(const std::string& seq, size_t start = 0) {
    std::string prot;
    prot.reserve((seq.size() - start) / 3 + 1);
    for (size_t i = start; i + 2 < seq.size(); i += 3) {
        char aa = translate_codon_fast(seq[i], seq[i + 1], seq[i + 2]);
        if (aa == '*') break;
        prot += aa;
    }
    return prot;
}

// ── Briggs (2007) damage model ────────────────────────────────────────────────
//
// Mechanistic model of ancient DNA deamination:
//
//  |<-- l (SS) -->|<------- DS interior ------->|<-- r (SS) -->|
//  5' [overhang]  [left of nick] m [right of nick] [overhang]  3'
//
// SS overhang (positions 0..l-1): C→T at rate delta_s (5') / G→A (DS 3') or C→T (SS 3')
// DS interior: C→T left of nick and G→A right of nick at rate delta_d
//
// Reference: Briggs et al. 2007, PNAS; implementation adapted from NGSNGS
// (github.com/RAHenriksen/NGSNGS, Biotin_ds_454Roche).

struct DamageRecord {
    std::vector<size_t> ct_positions;
    std::vector<size_t> ga_positions;
};

static DamageRecord apply_damage_briggs(std::string& seq,
                                        float nv,
                                        float lambda,
                                        float delta_s,
                                        float delta_d,
                                        bool  single_stranded,
                                        std::mt19937& rng) {
    DamageRecord rec;
    const int L = static_cast<int>(seq.size());
    if (L < 3) return rec;

    std::uniform_real_distribution<float> U(0.0f, 1.0f);
    std::geometric_distribution<int> geom(lambda);

    // Draw overhang lengths; retry until they fit
    int l = 0, r = 0;
    do {
        l = (U(rng) > 0.5f) ? geom(rng) : 0;
        r = (U(rng) > 0.5f) ? geom(rng) : 0;
    } while (l + r > L - 2);

    // 5' single-stranded overhang: C→T at rate delta_s
    for (int i = 0; i < l; ++i) {
        if (seq[i] == 'C' && U(rng) < delta_s) {
            seq[i] = 'T';
            rec.ct_positions.push_back(static_cast<size_t>(i));
        }
    }

    // 3' single-stranded overhang
    for (int j = 0; j < r; ++j) {
        int i = L - j - 1;
        if (single_stranded) {
            // SS library: C→T at 3' end as well
            if (seq[i] == 'C' && U(rng) < delta_s) {
                seq[i] = 'T';
                rec.ct_positions.push_back(static_cast<size_t>(i));
            }
        } else {
            // DS library: G→A at 3' end
            if (seq[i] == 'G' && U(rng) < delta_s) {
                seq[i] = 'A';
                rec.ga_positions.push_back(static_cast<size_t>(i));
            }
        }
    }

    // DS interior: place nick, apply delta_d background damage
    const int ds_span = L - l - r - 1;
    if (nv > 0.0f && ds_span > 0) {
        // Nick position m drawn from approximate geometric in [l, L-r-1)
        const float p_nick = nv / (static_cast<float>(ds_span) * nv + 1.0f - nv);
        int m = l;
        float cum = p_nick;
        float u_nick = U(rng);
        while (u_nick > cum && m < L - r - 1) {
            cum += p_nick;
            ++m;
        }

        // Left of nick (incl. m): C→T at rate delta_d
        for (int i = l; i <= m; ++i) {
            if (seq[i] == 'C' && U(rng) < delta_d) {
                seq[i] = 'T';
                rec.ct_positions.push_back(static_cast<size_t>(i));
            }
        }
        // Right of nick: G→A at rate delta_d (DS only)
        if (!single_stranded) {
            for (int i = m + 1; i < L - r; ++i) {
                if (seq[i] == 'G' && U(rng) < delta_d) {
                    seq[i] = 'A';
                    rec.ga_positions.push_back(static_cast<size_t>(i));
                }
            }
        }
    }

    std::sort(rec.ct_positions.begin(), rec.ct_positions.end());
    std::sort(rec.ga_positions.begin(), rec.ga_positions.end());
    return rec;
}

// ── Ground truth helpers ──────────────────────────────────────────────────────

static std::string find_damage_stops(const std::string& orig, const std::string& dmg,
                                     size_t start) {
    std::string result;
    size_t n = std::min(orig.size(), dmg.size());
    for (size_t i = start; i + 2 < n; i += 3) {
        if (is_stop(orig[i], orig[i + 1], orig[i + 2])) continue;
        if (is_stop(dmg[i],  dmg[i + 1],  dmg[i + 2])) {
            if (!result.empty()) result += ',';
            result += std::to_string((i - start) / 3);
        }
    }
    return result.empty() ? "." : result;
}

static std::string fmt_positions(const std::vector<size_t>& v) {
    if (v.empty()) return ".";
    std::string s;
    for (size_t x : v) { if (!s.empty()) s += ','; s += std::to_string(x); }
    return s;
}

// ── Simulator ─────────────────────────────────────────────────────────────────

struct SimParams {
    size_t   n_reads         = 1000;
    size_t   read_length     = 100;
    float    damage_fraction = 1.0f;   // fraction of reads that receive damage
    float    orf_fraction    = 0.70f;  // fraction of read length used for ORF
    size_t   min_orf_aa      = 10;
    uint64_t seed            = 42;
    bool     single_stranded = false;
    std::string output_fastq;
    std::string output_truth;

    // Briggs model parameters
    float    nv      = 0.024f;   // nick probability
    float    lambda  = 0.36f;    // geometric overhang parameter
    float    delta_s = 0.68f;    // SS overhang deamination rate
    float    delta_d = 0.0097f;  // DS interior deamination rate
};

static void print_usage(const char* prog) {
    std::cerr <<
        "Usage: " << prog << " simulate [options]\n\n"
        "Generate synthetic ancient DNA reads using the Briggs (2007) mechanistic\n"
        "damage model. Embeds a random ORF in each read and writes a ground-truth\n"
        "TSV for validating damage-aware frame selection.\n\n"
        "General options:\n"
        "  -n, --n-reads N           Number of reads (default: 1000)\n"
        "  --read-length N           Read length in bp (default: 100)\n"
        "  --damage-fraction FLOAT   Fraction of reads that receive damage (default: 1.0)\n"
        "  --library-type ds|ss      ds = double-stranded, ss = single-stranded (default: ds)\n"
        "  --orf-fraction FLOAT      Fraction of read length used for ORF (default: 0.70)\n"
        "  --min-orf-aa N            Minimum ORF length in amino acids (default: 10)\n"
        "  --seed INT                Random seed (default: 42)\n"
        "  -o, --output FILE         Output FASTQ (default: stdout)\n"
        "  --truth FILE              Output ground-truth TSV (optional)\n"
        "  -h, --help                Show this message\n\n"
        "Briggs damage model parameters:\n"
        "  --nv FLOAT                Nick probability in DS interior (default: 0.024)\n"
        "  --lambda FLOAT            Geometric overhang length parameter (default: 0.36)\n"
        "  --delta-s FLOAT           SS overhang deamination rate (default: 0.68)\n"
        "  --delta-d FLOAT           DS interior background deamination rate (default: 0.0097)\n\n"
        "Typical parameter sets:\n"
        "  Low damage:   --nv 0.01 --lambda 0.5 --delta-s 0.20 --delta-d 0.001\n"
        "  Medium:       --nv 0.024 --lambda 0.36 --delta-s 0.68 --delta-d 0.0097  (default)\n"
        "  High damage:  --nv 0.05 --lambda 0.25 --delta-s 0.90 --delta-d 0.02\n\n"
        "Ground-truth TSV columns:\n"
        "  read_id  frame  orf_nt_start  orf_nt_end  is_damaged  overhang_5  overhang_3\n"
        "  original_protein  damaged_protein  damage_induced_stops\n"
        "  ct_positions  ga_positions\n";
}

static int run_simulate(const SimParams& p) {
    std::mt19937 rng(p.seed);
    std::uniform_real_distribution<float> U(0.0f, 1.0f);
    std::uniform_int_distribution<int>    frame_dist(0, 2);

    std::ostream* fq_out = &std::cout;
    std::unique_ptr<std::ofstream> fq_file;
    if (!p.output_fastq.empty()) {
        fq_file = std::make_unique<std::ofstream>(p.output_fastq);
        if (!*fq_file) {
            std::cerr << "simulate: cannot open '" << p.output_fastq << "'\n";
            return 1;
        }
        fq_out = fq_file.get();
    }

    std::unique_ptr<std::ofstream> truth_file;
    std::ostream* truth_out = nullptr;
    if (!p.output_truth.empty()) {
        truth_file = std::make_unique<std::ofstream>(p.output_truth);
        if (!*truth_file) {
            std::cerr << "simulate: cannot open '" << p.output_truth << "'\n";
            return 1;
        }
        truth_out = truth_file.get();
        *truth_out << "read_id\tframe\torf_nt_start\torf_nt_end\tis_damaged"
                      "\toriginal_protein\tdamaged_protein"
                      "\tdamage_induced_stops\tct_positions\tga_positions\n";
    }

    size_t target_aa = static_cast<size_t>(p.read_length * p.orf_fraction / 3.0f);
    if (target_aa < p.min_orf_aa) target_aa = p.min_orf_aa;

    const std::string base_qual(p.read_length, 'I');  // Q40

    size_t n_damaged = 0, n_damage_stops = 0;

    for (size_t r = 0; r < p.n_reads; ++r) {
        const std::string read_id = "sim_" + std::to_string(r + 1);
        const bool is_damaged = (U(rng) < p.damage_fraction);

        const int frame = frame_dist(rng);

        std::string seq(p.read_length, 'N');
        for (char& c : seq) c = rand_base(rng);

        // Embed ORF at frame offset
        size_t avail_aa = (p.read_length - static_cast<size_t>(frame)) / 3;
        size_t orf_aa   = std::min(target_aa, avail_aa);
        size_t orf_start = 0, orf_end = 0;

        if (orf_aa >= p.min_orf_aa) {
            std::string orf = generate_orf(rng, orf_aa);
            orf_start = static_cast<size_t>(frame);
            orf_end   = orf_start + orf.size();
            seq.replace(orf_start, orf.size(), orf);
        }

        const std::string original_seq = seq;

        DamageRecord drec;
        if (is_damaged) {
            drec = apply_damage_briggs(seq,
                                       p.nv, p.lambda, p.delta_s, p.delta_d,
                                       p.single_stranded, rng);
            ++n_damaged;
        }

        *fq_out << '@' << read_id << '\n'
                << seq << '\n'
                << "+\n"
                << base_qual << '\n';

        if (truth_out && orf_end > orf_start) {
            const std::string orig_prot = translate(original_seq, orf_start);
            const std::string dmg_prot  = translate(seq, orf_start);
            const std::string dmg_stops = find_damage_stops(original_seq, seq, orf_start);
            if (dmg_stops != ".") ++n_damage_stops;

            *truth_out << read_id    << '\t'
                       << frame      << '\t'
                       << orf_start  << '\t'
                       << orf_end    << '\t'
                       << (is_damaged ? '1' : '0') << '\t'
                       << orig_prot  << '\t'
                       << dmg_prot   << '\t'
                       << dmg_stops  << '\t'
                       << fmt_positions(drec.ct_positions) << '\t'
                       << fmt_positions(drec.ga_positions) << '\n';
        }
    }

    std::cerr << "simulate: wrote " << p.n_reads << " reads"
              << " (" << n_damaged << " damaged"
              << ", " << n_damage_stops << " with damage-induced stops)\n";
    return 0;
}

} // anonymous namespace

int cmd_simulate(int argc, char* argv[]) {
    SimParams params;

    for (int i = 1; i < argc; ++i) {
        auto need = [&](const char* flag) -> bool {
            if (i + 1 >= argc) {
                std::cerr << "simulate: " << flag << " requires an argument\n";
                return false;
            }
            return true;
        };

        if (!strcmp(argv[i], "-h") || !strcmp(argv[i], "--help")) {
            print_usage(argv[0]);
            return 0;
        } else if (!strcmp(argv[i], "-n") || !strcmp(argv[i], "--n-reads")) {
            if (!need(argv[i])) return 1;
            params.n_reads = std::stoull(argv[++i]);
        } else if (!strcmp(argv[i], "--read-length")) {
            if (!need(argv[i])) return 1;
            params.read_length = std::stoull(argv[++i]);
        } else if (!strcmp(argv[i], "--damage-fraction")) {
            if (!need(argv[i])) return 1;
            params.damage_fraction = std::stof(argv[++i]);
        } else if (!strcmp(argv[i], "--library-type")) {
            if (!need(argv[i])) return 1;
            ++i;
            if (!strcmp(argv[i], "ss") || !strcmp(argv[i], "SS"))
                params.single_stranded = true;
            else if (!strcmp(argv[i], "ds") || !strcmp(argv[i], "DS"))
                params.single_stranded = false;
            else {
                std::cerr << "simulate: --library-type must be 'ds' or 'ss'\n";
                return 1;
            }
        } else if (!strcmp(argv[i], "--orf-fraction")) {
            if (!need(argv[i])) return 1;
            params.orf_fraction = std::stof(argv[++i]);
        } else if (!strcmp(argv[i], "--min-orf-aa")) {
            if (!need(argv[i])) return 1;
            params.min_orf_aa = std::stoull(argv[++i]);
        } else if (!strcmp(argv[i], "--seed")) {
            if (!need(argv[i])) return 1;
            params.seed = std::stoull(argv[++i]);
        } else if (!strcmp(argv[i], "-o") || !strcmp(argv[i], "--output")) {
            if (!need(argv[i])) return 1;
            params.output_fastq = argv[++i];
        } else if (!strcmp(argv[i], "--truth")) {
            if (!need(argv[i])) return 1;
            params.output_truth = argv[++i];
        // Briggs parameters
        } else if (!strcmp(argv[i], "--nv")) {
            if (!need(argv[i])) return 1;
            params.nv = std::stof(argv[++i]);
        } else if (!strcmp(argv[i], "--lambda")) {
            if (!need(argv[i])) return 1;
            params.lambda = std::stof(argv[++i]);
        } else if (!strcmp(argv[i], "--delta-s")) {
            if (!need(argv[i])) return 1;
            params.delta_s = std::stof(argv[++i]);
        } else if (!strcmp(argv[i], "--delta-d")) {
            if (!need(argv[i])) return 1;
            params.delta_d = std::stof(argv[++i]);
        } else {
            std::cerr << "simulate: unknown option '" << argv[i] << "'\n";
            print_usage(argv[0]);
            return 1;
        }
    }

    return run_simulate(params);
}

namespace {
    struct SimulateRegistrar {
        SimulateRegistrar() {
            SubcommandRegistry::instance().register_command(
                "simulate",
                "Generate synthetic ancient DNA reads (Briggs damage model)",
                cmd_simulate, 5);
        }
    } simulate_registrar;
}

}  // namespace cli
}  // namespace dart
