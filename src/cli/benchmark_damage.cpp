// Fast benchmark for per-read damage prediction
// Computes AUC-ROC comparing AGP predictions to synthetic ground truth

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <cmath>
#include <charconv>
#include <cerrno>
#include <cctype>
#include <cstring>
#include <cstdlib>
#include <cstdio>

struct Prediction {
    float damage_signal = 0;
    float damage_pct = 0;
    float p_damaged = 0;
    float score = 0;
};

struct Truth {
    int total_damage = 0;      // All damage events (from 'damage' column)
    int terminal_damage = 0;   // Damage at positions 1-15 from either end
    int inframe_damage = 0;    // Damage in coding frame (ct + ag)
    bool has_aa_change = false; // Has amino acid change (damage_aaseq_diffs non-empty)
    bool has_stop_damage = false; // Creates premature stop codon (>*)
};

namespace {

std::string shell_escape(const std::string& arg) {
    std::string escaped = "'";
    for (char c : arg) {
        if (c == '\'') {
            escaped += "'\\''";
        } else {
            escaped += c;
        }
    }
    escaped += "'";
    return escaped;
}

bool parse_int_token(const std::string& value, int& out) {
    const char* begin = value.data();
    const char* end = begin + value.size();
    while (begin < end && std::isspace(static_cast<unsigned char>(*begin))) ++begin;
    while (begin < end && std::isspace(static_cast<unsigned char>(*(end - 1)))) --end;
    if (begin == end) return false;

    auto res = std::from_chars(begin, end, out);
    return res.ec == std::errc() && res.ptr == end;
}

bool parse_float_token(const std::string& value, float& out) {
    const char* begin = value.c_str();
    while (*begin != '\0' && std::isspace(static_cast<unsigned char>(*begin))) ++begin;
    if (*begin == '\0') return false;

    char* end = nullptr;
    errno = 0;
    float parsed = std::strtof(begin, &end);
    if (begin == end || errno == ERANGE) return false;
    while (*end != '\0' && std::isspace(static_cast<unsigned char>(*end))) ++end;
    if (*end != '\0') return false;

    out = parsed;
    return true;
}

}  // namespace

// Parse GFF attributes
float get_attr_float(const std::string& attrs, const std::string& key, float default_val = 0.0f) {
    size_t pos = attrs.find(key + "=");
    if (pos == std::string::npos) return default_val;
    pos += key.length() + 1;
    size_t end = attrs.find(';', pos);
    if (end == std::string::npos) end = attrs.length();

    float parsed = default_val;
    if (!parse_float_token(attrs.substr(pos, end - pos), parsed)) return default_val;
    return parsed;
}

// Parse AGP GFF output
std::unordered_map<std::string, Prediction> parse_gff(const std::string& path) {
    std::unordered_map<std::string, Prediction> predictions;
    predictions.reserve(10000000);

    std::ifstream f(path);
    std::string line;
    size_t count = 0;

    while (std::getline(f, line)) {
        if (line.empty() || line[0] == '#') continue;

        // Parse tab-separated fields
        size_t tab1 = line.find('\t');
        if (tab1 == std::string::npos) continue;

        std::string read_name = line.substr(0, tab1);

        // Find score (field 6, 0-indexed field 5)
        size_t pos = tab1;
        for (int i = 0; i < 4; i++) {
            pos = line.find('\t', pos + 1);
            if (pos == std::string::npos) break;
        }
        if (pos == std::string::npos) continue;

        size_t score_start = pos + 1;
        size_t score_end = line.find('\t', score_start);
        std::string score_str = line.substr(score_start, score_end - score_start);
        float score = 0.0f;
        if (!parse_float_token(score_str, score)) continue;

        // Find attributes (field 9)
        size_t attr_start = line.rfind('\t');
        if (attr_start == std::string::npos) continue;
        std::string attrs = line.substr(attr_start + 1);

        // Extract damage metrics
        Prediction pred;
        pred.damage_signal = get_attr_float(attrs, "damage_signal");
        pred.damage_pct = get_attr_float(attrs, "damage_pct");
        pred.p_damaged = get_attr_float(attrs, "p_damaged");
        pred.score = score;

        // Keep best prediction per read
        auto it = predictions.find(read_name);
        if (it == predictions.end() || score > it->second.score) {
            predictions[read_name] = pred;
        }

        if (++count % 1000000 == 0) {
            std::cerr << "  Parsed " << count / 1000000 << "M GFF lines...\r" << std::flush;
        }
    }
    std::cerr << "  Parsed " << count << " GFF lines, " << predictions.size() << " unique reads\n";

    return predictions;
}

// Parse ground truth TSV (gzipped)
std::unordered_map<std::string, Truth> parse_truth(const std::string& path) {
    std::unordered_map<std::string, Truth> truth;
    truth.reserve(10000000);

    // Use zcat/pigz for gzip
    std::string escaped_path = shell_escape(path);
    std::string cmd = "pigz -dc " + escaped_path + " 2>/dev/null || zcat " + escaped_path;
    FILE* pipe = popen(cmd.c_str(), "r");
    if (!pipe) {
        std::cerr << "Error opening: " << path << "\n";
        return truth;
    }

    // Use getline for arbitrary length lines
    char* line_buf = nullptr;
    size_t line_cap = 0;
    ssize_t line_len;

    // Read header
    line_len = getline(&line_buf, &line_cap, pipe);
    if (line_len <= 0) {
        free(line_buf);
        pclose(pipe);
        return truth;
    }

    // Find column indices from header
    std::string header(line_buf);
    std::vector<std::string> cols;
    std::istringstream hss(header);
    std::string col;
    while (std::getline(hss, col, '\t')) {
        while (!col.empty() && (col.back() == '\r' || col.back() == '\n')) col.pop_back();
        cols.push_back(col);
    }

    int read_name_idx = -1, damage_idx = -1, ct_idx = -1, ag_idx = -1;
    int aa_diffs_idx = -1, codon_diffs_idx = -1;
    for (size_t i = 0; i < cols.size(); i++) {
        if (cols[i] == "read_name") read_name_idx = i;
        else if (cols[i] == "damage") damage_idx = i;
        else if (cols[i] == "damage_inframe_ct") ct_idx = i;
        else if (cols[i] == "damage_inframe_ag") ag_idx = i;
        else if (cols[i] == "damage_aaseq_diffs") aa_diffs_idx = i;
        else if (cols[i] == "damage_codon_diffs") codon_diffs_idx = i;
    }

    if (read_name_idx < 0 || damage_idx < 0) {
        std::cerr << "Missing required columns (read_name, damage) in truth file\n";
        std::cerr << "  Found columns: ";
        for (const auto& c : cols) std::cerr << c << " ";
        std::cerr << "\n";
        free(line_buf);
        pclose(pipe);
        return truth;
    }

    int max_idx = std::max({read_name_idx, damage_idx, ct_idx, ag_idx, aa_diffs_idx, codon_diffs_idx});
    size_t count = 0;
    size_t invalid_damage_tokens = 0;

    while ((line_len = getline(&line_buf, &line_cap, pipe)) > 0) {
        // Fast tab-delimited parsing
        std::vector<std::string> parts;
        parts.reserve(35);

        char* start = line_buf;
        char* end = line_buf;
        while (*end) {
            if (*end == '\t' || *end == '\n' || *end == '\r') {
                parts.emplace_back(start, end - start);
                if (*end == '\n' || *end == '\r') break;
                start = end + 1;
            }
            end++;
        }
        if (start < end) {
            parts.emplace_back(start, end - start);
        }

        if ((int)parts.size() <= max_idx) continue;

        Truth t;

        // Parse 'damage' column: comma-separated positions, "None" if no damage
        const std::string& damage_str = parts[damage_idx];
        if (!damage_str.empty() && damage_str != "None") {
            // Count damage events and check for terminal positions
            size_t pos = 0;
            while (pos < damage_str.length()) {
                size_t comma = damage_str.find(',', pos);
                if (comma == std::string::npos) comma = damage_str.length();
                std::string num_str = damage_str.substr(pos, comma - pos);
                if (!num_str.empty()) {
                    int damage_pos = 0;
                    if (parse_int_token(num_str, damage_pos)) {
                        t.total_damage++;
                        // Terminal = position 1-15 from either end
                        if ((damage_pos > 0 && damage_pos <= 15) ||
                            (damage_pos < 0 && damage_pos >= -15)) {
                            t.terminal_damage++;
                        }
                    } else {
                        ++invalid_damage_tokens;
                    }
                }
                pos = comma + 1;
            }
        }

        // Parse inframe damage counts
        if (ct_idx >= 0 && ct_idx < (int)parts.size()) {
            const std::string& ct_str = parts[ct_idx];
            if (!ct_str.empty() && ct_str != "None") {
                // Count comma-separated values
                size_t n = 1;
                for (char c : ct_str) if (c == ',') n++;
                t.inframe_damage += n;
            }
        }
        if (ag_idx >= 0 && ag_idx < (int)parts.size()) {
            const std::string& ag_str = parts[ag_idx];
            if (!ag_str.empty() && ag_str != "None") {
                size_t n = 1;
                for (char c : ag_str) if (c == ',') n++;
                t.inframe_damage += n;
            }
        }

        // Check for AA changes
        if (aa_diffs_idx >= 0 && aa_diffs_idx < (int)parts.size()) {
            const std::string& aa_str = parts[aa_diffs_idx];
            t.has_aa_change = !aa_str.empty() && aa_str != "None";
        }

        // Check for stop-creating damage (">*" in codon_diffs)
        if (codon_diffs_idx >= 0 && codon_diffs_idx < (int)parts.size()) {
            const std::string& codon_str = parts[codon_diffs_idx];
            t.has_stop_damage = codon_str.find(">*") != std::string::npos;
        }

        truth[parts[read_name_idx]] = t;

        if (++count % 1000000 == 0) {
            std::cerr << "  Parsed " << count / 1000000 << "M truth lines...\r" << std::flush;
        }
    }
    std::cerr << "  Parsed " << count << " truth lines                    \n";
    if (invalid_damage_tokens > 0) {
        std::cerr << "  Warning: ignored " << invalid_damage_tokens
                  << " malformed damage position tokens\n";
    }

    free(line_buf);
    pclose(pipe);
    return truth;
}

// Compute AUC-ROC efficiently using sorted scores
double compute_auc(const std::vector<std::pair<float, bool>>& scored_labels) {
    // scored_labels: [(score, has_damage), ...] sorted by score descending

    long long n_pos = 0, n_neg = 0;
    for (const auto& sl : scored_labels) {
        if (sl.second) n_pos++;
        else n_neg++;
    }

    if (n_pos == 0 || n_neg == 0) return 0.5;

    // Count inversions (Mann-Whitney U statistic)
    long long sum_ranks = 0;
    long long rank = 0;

    for (const auto& sl : scored_labels) {
        rank++;
        if (sl.second) {  // positive example
            sum_ranks += rank;
        }
    }

    // AUC = (sum_ranks - n_pos*(n_pos+1)/2) / (n_pos * n_neg)
    double auc = (double)(sum_ranks - n_pos * (n_pos + 1) / 2) / ((double)n_pos * n_neg);
    return 1.0 - auc;  // We want high score = damaged, but sorted descending
}

int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <predictions.gff> <truth.tsv.gz>\n";
        return 1;
    }

    std::string gff_path = argv[1];
    std::string truth_path = argv[2];

    std::cerr << "Loading AGP predictions from: " << gff_path << "\n";
    auto predictions = parse_gff(gff_path);
    std::cerr << "  Loaded " << predictions.size() << " predictions\n";

    std::cerr << "\nLoading ground truth from: " << truth_path << "\n";
    auto truth = parse_truth(truth_path);
    std::cerr << "  Loaded " << truth.size() << " truth records\n";

    // Match and collect - separate labels for each target and metric
    // Using p_damaged
    std::vector<std::pair<float, bool>> pd_any_damage;
    std::vector<std::pair<float, bool>> pd_terminal;
    std::vector<std::pair<float, bool>> pd_aa_change;
    std::vector<std::pair<float, bool>> pd_stop_damage;
    // Using damage_signal
    std::vector<std::pair<float, bool>> ds_any_damage;
    std::vector<std::pair<float, bool>> ds_terminal;
    std::vector<std::pair<float, bool>> ds_aa_change;
    std::vector<std::pair<float, bool>> ds_stop_damage;

    std::vector<std::pair<float, int>> signal_vs_true;

    size_t matched = 0;
    size_t n_any = 0, n_terminal = 0, n_aa_change = 0, n_stop = 0;
    double sum_signal = 0, sum_p_damaged = 0;

    for (const auto& [read_name, pred] : predictions) {
        auto it = truth.find(read_name);
        if (it == truth.end()) continue;

        matched++;
        const Truth& t = it->second;

        bool has_any = t.total_damage > 0;
        bool has_terminal = t.terminal_damage > 0;

        pd_any_damage.emplace_back(pred.p_damaged, has_any);
        pd_terminal.emplace_back(pred.p_damaged, has_terminal);
        pd_aa_change.emplace_back(pred.p_damaged, t.has_aa_change);
        pd_stop_damage.emplace_back(pred.p_damaged, t.has_stop_damage);

        ds_any_damage.emplace_back(pred.damage_signal, has_any);
        ds_terminal.emplace_back(pred.damage_signal, has_terminal);
        ds_aa_change.emplace_back(pred.damage_signal, t.has_aa_change);
        ds_stop_damage.emplace_back(pred.damage_signal, t.has_stop_damage);

        signal_vs_true.emplace_back(pred.damage_signal, t.total_damage);

        sum_signal += pred.damage_signal;
        sum_p_damaged += pred.p_damaged;

        if (has_any) n_any++;
        if (has_terminal) n_terminal++;
        if (t.has_aa_change) n_aa_change++;
        if (t.has_stop_damage) n_stop++;
    }

    if (matched == 0) {
        std::cerr << "No overlapping reads between predictions and truth; cannot compute metrics.\n";
        return 1;
    }

    std::cout << "\n=== Ground Truth Distribution ===\n";
    std::cout << "Matched reads: " << matched << "\n";
    std::cout << "  Any damage:       " << n_any << " (" << 100.0 * n_any / matched << "%)\n";
    std::cout << "  Terminal damage:  " << n_terminal << " (" << 100.0 * n_terminal / matched << "%)\n";
    std::cout << "  AA change:        " << n_aa_change << " (" << 100.0 * n_aa_change / matched << "%)\n";
    std::cout << "  Stop-creating:    " << n_stop << " (" << 100.0 * n_stop / matched << "%)\n";
    std::cout << "\nMean p_damaged: " << sum_p_damaged / matched << "\n";

    // Compute correlation
    double mean_signal = sum_signal / matched;
    double mean_true = 0;
    for (const auto& st : signal_vs_true) mean_true += st.second;
    mean_true /= matched;

    double cov = 0, var_signal = 0, var_true = 0;
    for (const auto& st : signal_vs_true) {
        double ds = st.first - mean_signal;
        double dt = st.second - mean_true;
        cov += ds * dt;
        var_signal += ds * ds;
        var_true += dt * dt;
    }

    double corr = (var_signal > 0 && var_true > 0) ?
        cov / std::sqrt(var_signal * var_true) : 0;
    std::cout << "Correlation (damage_signal vs count): " << corr << "\n";

    // Compute AUC-ROC for each target
    auto sort_and_auc = [](std::vector<std::pair<float, bool>>& labels) {
        std::sort(labels.begin(), labels.end(),
                  [](const auto& a, const auto& b) { return a.first > b.first; });
        return compute_auc(labels);
    };

    // p_damaged AUC
    double pd_auc_any = sort_and_auc(pd_any_damage);
    double pd_auc_terminal = sort_and_auc(pd_terminal);
    double pd_auc_aa = sort_and_auc(pd_aa_change);
    double pd_auc_stop = sort_and_auc(pd_stop_damage);

    // damage_signal AUC
    double ds_auc_any = sort_and_auc(ds_any_damage);
    double ds_auc_terminal = sort_and_auc(ds_terminal);
    double ds_auc_aa = sort_and_auc(ds_aa_change);
    double ds_auc_stop = sort_and_auc(ds_stop_damage);

    std::cout << "\n=== AUC-ROC by metric and target ===\n";
    std::cout << "Target               | p_damaged | dmg_signal | Prevalence\n";
    std::cout << "---------------------|-----------|------------|------------\n";
    printf("Any damage           |   %7.4f |    %7.4f | %5.1f%%\n", pd_auc_any, ds_auc_any, 100.0 * n_any / matched);
    printf("Terminal damage      |   %7.4f |    %7.4f | %5.1f%%\n", pd_auc_terminal, ds_auc_terminal, 100.0 * n_terminal / matched);
    printf("AA change            |   %7.4f |    %7.4f | %5.1f%%\n", pd_auc_aa, ds_auc_aa, 100.0 * n_aa_change / matched);
    printf("Stop-creating damage |   %7.4f |    %7.4f | %5.1f%%\n", pd_auc_stop, ds_auc_stop, 100.0 * n_stop / matched);

    // Breakdown by damage level
    std::cout << "\n=== Mean p_damaged by damage characteristics ===\n";

    // By total damage count
    std::cout << "\nBy total damage count:\n";
    struct Bin { int lo, hi; const char* label; };
    Bin bins[] = {{0, 0, "None"}, {1, 2, "1-2"}, {3, 5, "3-5"}, {6, 10, "6-10"}, {11, 1000, "11+"}};

    for (const auto& bin : bins) {
        double sum_pd = 0;
        int count = 0;

        for (const auto& [read_name, pred] : predictions) {
            auto it = truth.find(read_name);
            if (it == truth.end()) continue;

            int td = it->second.total_damage;
            if (td >= bin.lo && td <= bin.hi) {
                sum_pd += pred.p_damaged;
                count++;
            }
        }

        if (count > 0) {
            printf("  %-8s: %.4f (n=%d)\n", bin.label, sum_pd / count, count);
        }
    }

    // By terminal damage
    std::cout << "\nBy terminal damage (pos 1-15):\n";
    for (int terminal = 0; terminal <= 3; terminal++) {
        double sum_pd = 0;
        int count = 0;

        for (const auto& [read_name, pred] : predictions) {
            auto it = truth.find(read_name);
            if (it == truth.end()) continue;

            int td = it->second.terminal_damage;
            if ((terminal < 3 && td == terminal) || (terminal == 3 && td >= 3)) {
                sum_pd += pred.p_damaged;
                count++;
            }
        }

        if (count > 0) {
            const char* label = (terminal == 0) ? "0" : (terminal == 1) ? "1" :
                               (terminal == 2) ? "2" : "3+";
            printf("  %-8s: %.4f (n=%d)\n", label, sum_pd / count, count);
        }
    }

    // By functional impact
    std::cout << "\nBy functional impact:\n";
    double sum_no_aa = 0, sum_aa = 0, sum_stop = 0;
    int n_no_aa = 0, n_aa_only = 0, n_stop_only = 0;

    for (const auto& [read_name, pred] : predictions) {
        auto it = truth.find(read_name);
        if (it == truth.end()) continue;

        const Truth& t = it->second;
        if (t.has_stop_damage) {
            sum_stop += pred.p_damaged;
            n_stop_only++;
        } else if (t.has_aa_change) {
            sum_aa += pred.p_damaged;
            n_aa_only++;
        } else {
            sum_no_aa += pred.p_damaged;
            n_no_aa++;
        }
    }

    if (n_no_aa > 0) printf("  No AA change:  %.4f (n=%d)\n", sum_no_aa / n_no_aa, n_no_aa);
    if (n_aa_only > 0) printf("  AA change:     %.4f (n=%d)\n", sum_aa / n_aa_only, n_aa_only);
    if (n_stop_only > 0) printf("  Stop damage:   %.4f (n=%d)\n", sum_stop / n_stop_only, n_stop_only);

    return 0;
}
