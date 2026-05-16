#include "cli/subcommand.hpp"
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <string>
#include <unordered_map>

// dart eval — fast evaluation of dart predict output against simulation truth.
//
// Reads the truth TSV produced by dart simulate and one or two GFF files,
// then reports top-1 accuracy (acc), stop-read accuracy (stp), and recall
// (rec/srec: is the correct frame present anywhere in the output?).
//
// Usage:
//   dart eval --truth <truth.tsv> <predictions.gff>
//   dart eval --truth <truth.tsv> --ref <ref.gff> <predictions.gff>

namespace {

struct ReadTruth {
    int8_t   strand_is_fwd;   // 1 = +, 0 = -
    int8_t   frame_offset;    // 0, 1, or 2
    bool     has_dmg_stop;
    uint64_t stop_positions;  // bitmask: bit i set if codon i was damage-induced stop
    std::string original_protein;
};

struct ReadPred {
    float   best_score;
    int8_t  best_strand_fwd;
    int8_t  best_frame;
    uint8_t seen;  // bitmask: bits 0-2 = fwd frames 0-2, bits 3-5 = rev frames 0-2
};

// Split a line on tabs, fill fields[], return count of fields found
static int split_tabs(char* line, char* fields[], int max_fields) {
    int n = 0;
    char* p = line;
    while (n < max_fields) {
        fields[n++] = p;
        p = strchr(p, '\t');
        if (!p) break;
        *p++ = '\0';
    }
    // Strip trailing newline from last field
    if (n > 0) {
        size_t len = strlen(fields[n-1]);
        while (len > 0 && (fields[n-1][len-1] == '\n' || fields[n-1][len-1] == '\r'))
            fields[n-1][--len] = '\0';
    }
    return n;
}

static std::unordered_map<std::string, ReadTruth> load_truth(const char* path) {
    std::unordered_map<std::string, ReadTruth> truth;
    truth.reserve(12000000);
    FILE* f = fopen(path, "r");
    if (!f) { fprintf(stderr, "eval: cannot open truth file: %s\n", path); return truth; }
    char buf[4096];
    char* fields[16];
    bool header = true;
    while (fgets(buf, sizeof(buf), f)) {
        if (header) { header = false; continue; }  // skip header
        int n = split_tabs(buf, fields, 10);
        if (n < 8) continue;
        int frame = atoi(fields[1]);
        bool has_stop = (fields[7][0] != '.' || fields[7][1] != '\0');
        ReadTruth t;
        t.strand_is_fwd = (frame < 3) ? 1 : 0;
        t.frame_offset  = frame % 3;
        t.has_dmg_stop  = has_stop;
        if (n >= 6) t.original_protein = fields[5];
        // Parse damage_induced_stops (field 7): comma-separated 0-based codon positions
        t.stop_positions = 0;
        if (n >= 8 && fields[7][0] != '.') {
            char* p = fields[7];
            while (*p) {
                int pos = atoi(p);
                if (pos >= 0 && pos < 64) t.stop_positions |= (1ULL << pos);
                p = strchr(p, ',');
                if (!p) break;
                p++;
            }
        }
        truth[fields[0]] = std::move(t);
    }
    fclose(f);
    return truth;
}

static std::unordered_map<std::string, ReadPred> load_gff(const char* path) {
    std::unordered_map<std::string, ReadPred> preds;
    preds.reserve(12000000);
    FILE* f = fopen(path, "r");
    if (!f) { fprintf(stderr, "eval: cannot open GFF: %s\n", path); return preds; }
    char buf[4096];
    char* fields[10];
    while (fgets(buf, sizeof(buf), f)) {
        if (buf[0] == '#') continue;
        int n = split_tabs(buf, fields, 9);
        if (n < 8) continue;
        float sc = (float)atof(fields[5]);
        int8_t strand_fwd = (fields[6][0] == '+') ? 1 : 0;
        int8_t frame = (int8_t)atoi(fields[7]);
        const char* rid = fields[0];
        auto it = preds.find(rid);
        if (it == preds.end()) {
            ReadPred p; p.best_score = sc; p.best_strand_fwd = strand_fwd;
            p.best_frame = frame; p.seen = 0;
            int bit = strand_fwd ? frame : (3 + frame);
            p.seen |= (1u << bit);
            preds[rid] = p;
        } else {
            if (sc > it->second.best_score) {
                it->second.best_score = sc;
                it->second.best_strand_fwd = strand_fwd;
                it->second.best_frame = frame;
            }
            int bit = strand_fwd ? frame : (3 + frame);
            it->second.seen |= (1u << bit);
        }
    }
    fclose(f);
    return preds;
}

// Load top-1 protein sequences from --fasta-aa output.
// Header format: >{read_id}_{strand}_{frame} coord ... corr_score={rank}
// corr_score=0 marks the top-1 prediction for each read.
// read_id is recovered by stripping the trailing _{strand}_{frame} (4 chars).
static std::unordered_map<std::string, std::string> load_proteins(const char* path) {
    std::unordered_map<std::string, std::string> proteins;
    proteins.reserve(12000000);
    FILE* f = fopen(path, "r");
    if (!f) { fprintf(stderr, "eval: cannot open proteins FASTA: %s\n", path); return proteins; }
    char buf[65536];
    std::string cur_rid;
    while (fgets(buf, sizeof(buf), f)) {
        size_t len = strlen(buf);
        while (len > 0 && (buf[len-1] == '\n' || buf[len-1] == '\r')) buf[--len] = '\0';
        if (buf[0] == '>') {
            cur_rid.clear();
            // First token is the gene ID: {read_id}_{strand}_{frame}
            const char* sp = strchr(buf + 1, ' ');
            size_t id_len = sp ? (size_t)(sp - buf - 1) : strlen(buf + 1);
            // Check corr_score=0 in rest of header (top-1 only)
            const char* rest = sp ? sp : "";
            if (!strstr(rest, "corr_score=0")) continue;
            // Strip trailing _{strand}_{frame}: always 4 chars (_[+-]_[012])
            if (id_len < 4) continue;
            cur_rid.assign(buf + 1, id_len - 4);
        } else if (!cur_rid.empty()) {
            proteins[cur_rid] += buf;
        }
    }
    fclose(f);
    return proteins;
}

struct TranslationCounts {
    long aa_total = 0, aa_correct = 0;        // all correctly-framed reads
    long stop_aa_total = 0, stop_aa_correct = 0; // at damage-induced stop positions
};

static TranslationCounts score_translation(
    const std::unordered_map<std::string, ReadTruth>& truth,
    const std::unordered_map<std::string, ReadPred>&  preds,
    const std::unordered_map<std::string, std::string>& proteins)
{
    TranslationCounts tc;
    for (const auto& [rid, t] : truth) {
        auto pit = preds.find(rid);
        if (pit == preds.end()) continue;
        const ReadPred& p = pit->second;
        if (p.best_strand_fwd != t.strand_is_fwd || p.best_frame != t.frame_offset) continue;
        auto protp = proteins.find(rid);
        if (protp == proteins.end()) continue;
        const std::string& pred = protp->second;
        const std::string& orig = t.original_protein;
        size_t cmp = std::min(pred.size(), orig.size());
        for (size_t i = 0; i < cmp; i++) {
            tc.aa_total++;
            if (pred[i] == orig[i]) tc.aa_correct++;
            if (i < 64 && (t.stop_positions >> i) & 1u) {
                tc.stop_aa_total++;
                if (pred[i] == orig[i]) tc.stop_aa_correct++;
            }
        }
    }
    return tc;
}

struct Counts {
    long total=0, correct=0, recall_total=0, recall_correct=0;
    long stop_total=0, stop_correct=0, srecall_total=0, srecall_correct=0;
};

static Counts score(
    const std::unordered_map<std::string,ReadTruth>& truth,
    const std::unordered_map<std::string,ReadPred>&  preds)
{
    Counts c;
    for (const auto& [rid, t] : truth) {
        c.total++;
        bool has_pred = false;
        bool top1_ok  = false;
        bool found    = false;
        auto it = preds.find(rid);
        if (it != preds.end()) {
            const ReadPred& p = it->second;
            has_pred = true;
            top1_ok  = (p.best_strand_fwd == t.strand_is_fwd &&
                        p.best_frame      == t.frame_offset);
            int bit  = t.strand_is_fwd ? t.frame_offset : (3 + t.frame_offset);
            found    = (p.seen >> bit) & 1u;
        }
        if (top1_ok) c.correct++;
        if (has_pred) {
            c.recall_total++;
            if (found) c.recall_correct++;
        }
        if (t.has_dmg_stop) {
            c.stop_total++;
            if (top1_ok)  c.stop_correct++;
            if (has_pred) {
                c.srecall_total++;
                if (found) c.srecall_correct++;
            }
        }
    }
    return c;
}

static void print_usage(const char* prog) {
    fprintf(stderr,
        "Usage: %s eval --truth <truth.tsv> [--ref <ref.gff>] [--proteins <faa>] <predictions.gff>\n"
        "\n"
        "Evaluate dart predict output against simulation truth (dart simulate output).\n"
        "\n"
        "Options:\n"
        "  --truth <file>     Truth TSV from dart simulate (required)\n"
        "  --ref <file>       Reference GFF to compare against (optional)\n"
        "  --proteins <file>            Protein FASTA from --fasta-aa (enables trl)\n"
        "  --proteins-corrected <file>  Corrected FASTA from --fasta-aa-corrected (enables trl_s)\n"
        "  -h, --help         Show this help\n"
        "\n"
        "Metrics:\n"
        "  acc    Top-1 accuracy (correct strand + frame offset)\n"
        "  stp    Top-1 accuracy on damage-induced-stop reads\n"
        "  rec    Recall: correct frame present anywhere in output\n"
        "  srec   Recall on damage-induced-stop reads\n"
        "  trl    AA accuracy on correctly-framed reads (requires --proteins)\n"
        "  trl_s  Ancestor AA accuracy at damage stops (requires --proteins-corrected)\n",
        prog);
}

}  // namespace

int cmd_eval(int argc, char* argv[]) {
    const char* truth_path             = nullptr;
    const char* ref_path               = nullptr;
    const char* pred_path              = nullptr;
    const char* proteins_path          = nullptr;
    const char* proteins_corrected_path = nullptr;

    for (int i = 1; i < argc; ++i) {
        if (!strcmp(argv[i], "-h") || !strcmp(argv[i], "--help")) {
            print_usage(argv[0]);
            return 0;
        } else if (!strcmp(argv[i], "--truth") && i+1 < argc) {
            truth_path = argv[++i];
        } else if (!strcmp(argv[i], "--ref") && i+1 < argc) {
            ref_path = argv[++i];
        } else if (!strcmp(argv[i], "--proteins") && i+1 < argc) {
            proteins_path = argv[++i];
        } else if (!strcmp(argv[i], "--proteins-corrected") && i+1 < argc) {
            proteins_corrected_path = argv[++i];
        } else if (argv[i][0] != '-') {
            pred_path = argv[i];
        }
    }

    if (!truth_path || !pred_path) {
        print_usage(argv[0]);
        return 1;
    }

    auto truth = load_truth(truth_path);
    auto preds = load_gff(pred_path);
    Counts c   = score(truth, preds);

    double acc  = c.total         ? 100.0*c.correct        /c.total         : 0;
    double sacc = c.stop_total    ? 100.0*c.stop_correct   /c.stop_total    : 0;
    double rec  = c.recall_total  ? 100.0*c.recall_correct /c.recall_total  : 0;
    double srec = c.srecall_total ? 100.0*c.srecall_correct/c.srecall_total : 0;

    if (ref_path) {
        auto refs   = load_gff(ref_path);
        Counts rc   = score(truth, refs);
        double racc  = rc.total      ? 100.0*rc.correct       /rc.total      : 0;
        double rsacc = rc.stop_total ? 100.0*rc.stop_correct  /rc.stop_total : 0;
        double rrec  = rc.recall_total  ? 100.0*rc.recall_correct /rc.recall_total  : 0;
        double rsrec = rc.srecall_total ? 100.0*rc.srecall_correct/rc.srecall_total : 0;
        printf("acc:  %.2f%%  (ref %.2f%%, Δ=%+.2f%%)\n",  acc,  racc,  acc-racc);
        printf("stp:  %.2f%%  (ref %.2f%%, Δ=%+.2f%%)  n_stops=%ld\n", sacc, rsacc, sacc-rsacc, c.stop_total);
        printf("rec:  %.2f%%  (ref %.2f%%, Δ=%+.2f%%)  correct frame in any output\n", rec, rrec, rec-rrec);
        printf("srec: %.2f%%  (ref %.2f%%, Δ=%+.2f%%)  recall on stop reads\n", srec, rsrec, srec-rsrec);
    } else {
        printf("acc:  %.2f%%  (n=%ld)\n",  acc,  c.total);
        printf("stp:  %.2f%%  (n_stops=%ld)\n", sacc, c.stop_total);
        printf("rec:  %.2f%%  (correct frame in any output)\n", rec);
        printf("srec: %.2f%%  (recall on stop reads)\n", srec);
    }

    if (proteins_path) {
        auto proteins = load_proteins(proteins_path);
        TranslationCounts tc = score_translation(truth, preds, proteins);
        double trl = tc.aa_total ? 100.0*tc.aa_correct/tc.aa_total : 0;
        printf("trl:  %.2f%%  AA accuracy on correctly-framed reads  (n_aa=%ld)\n",
               trl, tc.aa_total);
    }
    if (proteins_corrected_path) {
        auto proteins_corr = load_proteins(proteins_corrected_path);
        TranslationCounts tc = score_translation(truth, preds, proteins_corr);
        double trl_s = tc.stop_aa_total ? 100.0*tc.stop_aa_correct/tc.stop_aa_total : 0;
        printf("trl_s: %.2f%%  ancestor AA accuracy at damage-induced stop positions  (n_stops=%ld)\n",
               trl_s, tc.stop_aa_total);
    }
    return 0;
}

namespace {
    struct EvalRegistrar {
        EvalRegistrar() {
            dart::cli::SubcommandRegistry::instance().register_command(
                "eval", "Evaluate predictions against simulation truth",
                cmd_eval, 30);
        }
    } eval_registrar;
}
