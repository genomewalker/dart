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

// Fragment-length CDF sampler (inverse-CDF method).
// Built-in "ancient" distribution is the ngsngs Size_dist_sampling.txt empirical CDF
// from Henriksen et al. (github.com/RAHenriksen/NGSNGS), covering 34–191 bp.
struct LengthCDF {
    std::vector<int>   lengths;
    std::vector<float> cum_probs;

    bool empty() const { return lengths.empty(); }

    // Load from a two-column file: length<tab>cum_prob (ngsngs format, no header).
    static LengthCDF load(const std::string& path) {
        LengthCDF cdf;
        std::ifstream f(path);
        if (!f) { std::cerr << "simulate: cannot open length file '" << path << "'\n"; return cdf; }
        int len; float p;
        while (f >> len >> p) { cdf.lengths.push_back(len); cdf.cum_probs.push_back(p); }
        return cdf;
    }

    // Built-in ancient DNA distribution (ngsngs Size_dist_sampling.txt).
    static LengthCDF ancient() {
        static const int   L[] = {34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,186,187,188,189,190,191};
        static const float P[] = {0.0f,0.00540914f,0.01326621f,0.02248544f,0.03442894f,0.04907704f,0.06617304f,0.08449304f,0.10468844f,0.12308744f,0.14265114f,0.16255044f,0.18223264f,0.20317824f,0.22357114f,0.24354944f,0.26398174f,0.28548014f,0.30425424f,0.32255454f,0.34040074f,0.35840484f,0.37711964f,0.39389984f,0.40929814f,0.42568354f,0.44439834f,0.45944124f,0.47365504f,0.48830314f,0.50218524f,0.51625354f,0.52900044f,0.54143334f,0.55319714f,0.56450214f,0.57543434f,0.58573154f,0.59530694f,0.60476714f,0.61370274f,0.62219744f,0.63054194f,0.63852634f,0.64635984f,0.65379174f,0.66073204f,0.66759424f,0.67426534f,0.68092644f,0.68709294f,0.69319334f,0.69921284f,0.70498524f,0.71054444f,0.71588524f,0.72105194f,0.72616454f,0.73097114f,0.73575564f,0.74022214f,0.74468864f,0.74904394f,0.75323144f,0.75740094f,0.76145224f,0.76542764f,0.76918784f,0.77286594f,0.77643284f,0.77993364f,0.78332324f,0.78666674f,0.78990214f,0.79306944f,0.79614464f,0.79916574f,0.80200564f,0.80476344f,0.80741304f,0.80998054f,0.81251994f,0.81491994f,0.81729594f,0.81959984f,0.82184964f,0.82406534f,0.82619394f,0.82826044f,0.83031894f,0.83229334f,0.83422764f,0.83612784f,0.83797794f,0.83973194f,0.84148594f,0.84320584f,0.84491374f,0.84655754f,0.84819334f,0.84981514f,0.85138284f,0.85291644f,0.85440994f,0.85587534f,0.85729264f,0.85867984f,0.86003494f,0.86135994f,0.86264684f,0.86390164f,0.86511834f,0.86631094f,0.86748954f,0.86862804f,0.86976254f,0.87085894f,0.87192524f,0.87295744f,0.87398164f,0.87498174f,0.87595774f,0.87692174f,0.87783764f,0.87875354f,0.87965144f,0.88051124f,0.88135304f,0.88216874f,0.88298444f,0.88379214f,0.88457574f,0.88534134f,0.88610294f,0.88683244f,0.88756194f,0.88825534f,0.88895274f,0.88963614f,0.89030554f,0.89097494f,0.89161624f,0.89225754f,0.89289484f,0.89350004f,0.89410524f,0.89468834f,0.89526744f,0.89583254f,0.89639764f,0.89694274f,0.89747984f,0.89800894f,0.89851994f,0.89903094f,0.89953394f,0.90001484f,0.90049574f,0.90096464f,0.90143354f,0.98781947f,0.98874731f,0.98951723f,0.99030688f,0.99097809f,0.99170852f,0.99253766f,0.99312990f,0.99376162f,0.99435387f,0.99494611f,0.99565680f,0.99626878f,0.99693999f,0.99743352f,0.99792706f,0.99859826f,0.99911154f,0.99919051f,0.99924973f,0.99934844f,0.99942740f,0.99956559f,0.99962482f,0.99966430f,0.99976301f,0.99988145f,0.99992094f,0.99994068f,1.0f};
        LengthCDF cdf;
        cdf.lengths.assign(std::begin(L), std::end(L));
        cdf.cum_probs.assign(std::begin(P), std::end(P));
        return cdf;
    }

    int sample(std::mt19937& rng) const {
        float u = std::uniform_real_distribution<float>(0.0f, 1.0f)(rng);
        auto it = std::lower_bound(cum_probs.begin(), cum_probs.end(), u);
        size_t idx = static_cast<size_t>(it - cum_probs.begin());
        if (idx >= lengths.size()) idx = lengths.size() - 1;
        return lengths[idx];
    }
};

struct SimParams {
    size_t   n_reads         = 1000;
    size_t   read_length     = 100;
    std::string length_file;           // ngsngs CDF file; overrides built-in ancient dist
    bool        fixed_length = false;  // --read-length was explicit; disable ancient dist
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
        "  --length-file FILE        Two-column CDF file: length<tab>cum_prob (ngsngs format)\n"
        "                            Default: built-in ancient aDNA distribution (34–191 bp)\n"
        "  --read-length N           Use a fixed read length instead of the length distribution\n"
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

    // Build length sampler: ancient dist by default, file overrides, --read-length disables.
    LengthCDF len_cdf;
    if (!p.length_file.empty()) {
        len_cdf = LengthCDF::load(p.length_file);
        if (len_cdf.empty()) return 1;
    } else if (!p.fixed_length) {
        len_cdf = LengthCDF::ancient();
    }

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

    size_t n_damaged = 0, n_damage_stops = 0;

    for (size_t r = 0; r < p.n_reads; ++r) {
        const std::string read_id = "sim_" + std::to_string(r + 1);
        const bool is_damaged = (U(rng) < p.damage_fraction);

        const size_t read_length = len_cdf.empty()
            ? p.read_length
            : static_cast<size_t>(std::max(len_cdf.sample(rng),
                                           static_cast<int>(p.min_orf_aa * 3 + 2)));

        const int frame = frame_dist(rng);

        std::string seq(read_length, 'N');
        for (char& c : seq) c = rand_base(rng);

        size_t target_aa = static_cast<size_t>(read_length * p.orf_fraction / 3.0f);
        if (target_aa < p.min_orf_aa) target_aa = p.min_orf_aa;

        // Embed ORF at frame offset
        size_t avail_aa = (read_length - static_cast<size_t>(frame)) / 3;
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
                << std::string(read_length, 'I') << '\n';

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
            params.fixed_length = true;
        } else if (!strcmp(argv[i], "--length-file")) {
            if (!need(argv[i])) return 1;
            params.length_file = argv[++i];
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
