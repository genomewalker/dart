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
    int8_t  strand_is_fwd;   // 1 = +, 0 = -
    int8_t  frame_offset;    // 0, 1, or 2
    bool    has_dmg_stop;
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
        truth[fields[0]] = t;
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
        "Usage: %s eval --truth <truth.tsv> [--ref <ref.gff>] <predictions.gff>\n"
        "\n"
        "Evaluate dart predict output against simulation truth (dart simulate output).\n"
        "\n"
        "Options:\n"
        "  --truth <file>   Truth TSV from dart simulate (required)\n"
        "  --ref <file>     Reference GFF to compare against (optional)\n"
        "  -h, --help       Show this help\n"
        "\n"
        "Metrics:\n"
        "  acc   Top-1 accuracy (correct strand + frame offset)\n"
        "  stp   Top-1 accuracy on damage-induced-stop reads\n"
        "  rec   Recall: correct frame present anywhere in output\n"
        "  srec  Recall on damage-induced-stop reads\n",
        prog);
}

}  // namespace

int cmd_eval(int argc, char* argv[]) {
    const char* truth_path = nullptr;
    const char* ref_path   = nullptr;
    const char* pred_path  = nullptr;

    for (int i = 1; i < argc; ++i) {
        if (!strcmp(argv[i], "-h") || !strcmp(argv[i], "--help")) {
            print_usage(argv[0]);
            return 0;
        } else if (!strcmp(argv[i], "--truth") && i+1 < argc) {
            truth_path = argv[++i];
        } else if (!strcmp(argv[i], "--ref") && i+1 < argc) {
            ref_path = argv[++i];
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
