/**
 * analyze_errors.cpp - Fast error characterization for DART predictions
 *
 * Analyzes mispredicted reads to understand what makes them hard:
 * - Length distribution
 * - Score margins
 * - Damage levels
 * - Feature distributions
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <regex>
#include <algorithm>
#include <cmath>
#include <iomanip>

struct Prediction {
    float score;
    char strand;
    int frame;
    float damage_pct;
    int length;
};

struct GroundTruth {
    char strand;
    int frame;
    int start;
    int end;
    int length;
};

struct ReadAnalysis {
    int length;
    float score;
    float margin;
    float damage_pct;
    char pred_strand;
    int pred_frame;
    char true_strand;
    int true_frame;
    bool is_correct;
    bool strand_correct;
};

// Parse ground truth from read name
// Format: ...---N:ancient:STRAND:START:END:LENGTH:DAMAGE
bool parse_ground_truth(const std::string& read_name, GroundTruth& truth) {
    // Find the pattern ---N:ancient:
    size_t pos = read_name.find(":ancient:");
    if (pos == std::string::npos) return false;

    // Go back to find the ---N part
    size_t dash_pos = read_name.rfind("---", pos);
    if (dash_pos == std::string::npos) return false;

    // Parse from after "---"
    size_t start = dash_pos + 3;

    // Skip the index number and :ancient:
    pos = read_name.find(":ancient:", start);
    if (pos == std::string::npos) return false;

    // Now parse STRAND:START:END:LENGTH
    pos += 9; // skip ":ancient:"

    if (pos >= read_name.size()) return false;
    truth.strand = read_name[pos];
    pos += 2; // skip strand and ':'

    // Parse START
    size_t colon = read_name.find(':', pos);
    if (colon == std::string::npos) return false;
    truth.start = std::stoi(read_name.substr(pos, colon - pos));
    truth.frame = truth.start % 3;
    pos = colon + 1;

    // Parse END
    colon = read_name.find(':', pos);
    if (colon == std::string::npos) return false;
    truth.end = std::stoi(read_name.substr(pos, colon - pos));
    pos = colon + 1;

    // Parse LENGTH
    colon = read_name.find(':', pos);
    if (colon == std::string::npos) return false;
    truth.length = std::stoi(read_name.substr(pos, colon - pos));

    return true;
}

// Parse GFF attributes
float get_attr_float(const std::string& attrs, const std::string& key, float default_val = 0.0f) {
    std::string search = key + "=";
    size_t pos = attrs.find(search);
    if (pos == std::string::npos) return default_val;

    pos += search.size();
    size_t end = attrs.find(';', pos);
    if (end == std::string::npos) end = attrs.size();

    try {
        return std::stof(attrs.substr(pos, end - pos));
    } catch (const std::exception&) {
        return default_val;
    }
}

int get_attr_int(const std::string& attrs, const std::string& key, int default_val = 0) {
    std::string search = key + "=";
    size_t pos = attrs.find(search);
    if (pos == std::string::npos) return default_val;

    pos += search.size();
    size_t end = attrs.find(';', pos);
    if (end == std::string::npos) end = attrs.size();

    try {
        return std::stoi(attrs.substr(pos, end - pos));
    } catch (const std::exception&) {
        return default_val;
    }
}

void print_histogram(const std::string& title,
                     const std::vector<std::pair<float, float>>& bins,
                     const std::vector<ReadAnalysis>& correct,
                     const std::vector<ReadAnalysis>& mispred,
                     std::function<float(const ReadAnalysis&)> getter,
                     const std::string& unit = "") {
    std::cout << "=== " << title << " ===" << std::endl;
    std::cout << std::left << std::setw(14) << "Range"
              << std::right << std::setw(12) << "Correct"
              << std::setw(12) << "Mispred"
              << std::setw(10) << "Acc %" << std::endl;
    std::cout << std::string(50, '-') << std::endl;

    for (const auto& [lo, hi] : bins) {
        int correct_count = 0, mispred_count = 0;
        for (const auto& r : correct) {
            float val = getter(r);
            if (val >= lo && val < hi) correct_count++;
        }
        for (const auto& r : mispred) {
            float val = getter(r);
            if (val >= lo && val < hi) mispred_count++;
        }

        int total = correct_count + mispred_count;
        if (total > 0) {
            float acc = 100.0f * correct_count / total;
            std::ostringstream range;
            if (unit == "bp") {
                range << std::fixed << std::setprecision(0) << lo << "-" << hi << " " << unit;
            } else if (unit == "%") {
                range << std::fixed << std::setprecision(0) << lo << "-" << hi << unit;
            } else {
                range << std::fixed << std::setprecision(2) << lo << "-" << hi;
            }
            std::cout << std::left << std::setw(14) << range.str()
                      << std::right << std::setw(12) << correct_count
                      << std::setw(12) << mispred_count
                      << std::setw(9) << std::fixed << std::setprecision(1) << acc << "%" << std::endl;
        }
    }
    std::cout << std::endl;
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <predictions.gff>" << std::endl;
        return 1;
    }

    std::string gff_path = argv[1];
    std::ifstream file(gff_path);
    if (!file.is_open()) {
        std::cerr << "Error: Cannot open " << gff_path << std::endl;
        return 1;
    }

    // Collect predictions per read
    std::unordered_map<std::string, std::vector<Prediction>> read_predictions;

    std::string line;
    size_t line_count = 0;

    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#') continue;

        // Parse GFF line
        std::istringstream iss(line);
        std::string read_name, source, type;
        int start, end;
        float score;
        char strand;
        std::string phase, attrs;

        if (!(iss >> read_name >> source >> type >> start >> end >> score >> strand)) {
            continue;
        }
        iss >> phase;
        std::getline(iss, attrs);
        // Trim leading whitespace from attrs
        attrs.erase(0, attrs.find_first_not_of(" \t"));

        Prediction pred;
        pred.score = score;
        pred.strand = strand;
        pred.frame = get_attr_int(attrs, "frame", -1);
        pred.damage_pct = get_attr_float(attrs, "damage_pct", 0.0f);
        pred.length = end;

        read_predictions[read_name].push_back(pred);

        line_count++;
        if (line_count % 1000000 == 0) {
            std::cerr << "Processed " << line_count / 1000000 << "M lines..." << std::endl;
        }
    }

    std::cerr << "Total lines: " << line_count << ", unique reads: " << read_predictions.size() << std::endl;

    // Analyze each read
    std::vector<ReadAnalysis> correct_reads, mispred_reads;

    for (const auto& [read_name, preds] : read_predictions) {
        GroundTruth truth;
        if (!parse_ground_truth(read_name, truth)) continue;

        // Find best prediction
        const Prediction* best = &preds[0];
        for (const auto& p : preds) {
            if (p.score > best->score) best = &p;
        }

        // Calculate margin
        std::vector<float> scores;
        for (const auto& p : preds) scores.push_back(p.score);
        std::sort(scores.begin(), scores.end(), std::greater<float>());
        float margin = scores.size() > 1 ? scores[0] - scores[1] : 0.0f;

        // Create analysis record
        ReadAnalysis analysis;
        analysis.length = best->length;
        analysis.score = best->score;
        analysis.margin = margin;
        analysis.damage_pct = best->damage_pct;
        analysis.pred_strand = best->strand;
        analysis.pred_frame = best->frame;
        analysis.true_strand = truth.strand;
        analysis.true_frame = truth.frame;
        analysis.strand_correct = (best->strand == truth.strand);
        analysis.is_correct = (best->strand == truth.strand && best->frame == truth.frame);

        if (analysis.is_correct) {
            correct_reads.push_back(analysis);
        } else {
            mispred_reads.push_back(analysis);
        }
    }

    // Print summary
    size_t total = correct_reads.size() + mispred_reads.size();
    std::cout << "=== Error Analysis Summary ===" << std::endl;
    std::cout << "Total reads: " << total << std::endl;
    std::cout << "Correct: " << correct_reads.size() << " ("
              << std::fixed << std::setprecision(1) << 100.0 * correct_reads.size() / total << "%)" << std::endl;
    std::cout << "Mispredicted: " << mispred_reads.size() << " ("
              << std::fixed << std::setprecision(1) << 100.0 * mispred_reads.size() / total << "%)" << std::endl;
    std::cout << std::endl;

    // Strand vs frame errors
    int strand_errors = 0, frame_only_errors = 0;
    for (const auto& r : mispred_reads) {
        if (!r.strand_correct) strand_errors++;
        else frame_only_errors++;
    }
    std::cout << "Strand errors: " << strand_errors << " ("
              << std::fixed << std::setprecision(1) << 100.0 * strand_errors / mispred_reads.size() << "% of errors)" << std::endl;
    std::cout << "Frame-only errors: " << frame_only_errors << " ("
              << std::fixed << std::setprecision(1) << 100.0 * frame_only_errors / mispred_reads.size() << "% of errors)" << std::endl;
    std::cout << std::endl;

    // Length distribution
    std::vector<std::pair<float, float>> length_bins = {
        {30, 50}, {50, 75}, {75, 100}, {100, 150}, {150, 200}, {200, 500}
    };
    print_histogram("Length Distribution", length_bins, correct_reads, mispred_reads,
                    [](const ReadAnalysis& r) { return (float)r.length; }, "bp");

    // Score margin distribution
    std::vector<std::pair<float, float>> margin_bins = {
        {0, 0.02}, {0.02, 0.05}, {0.05, 0.10}, {0.10, 0.20}, {0.20, 0.50}, {0.50, 1.0}
    };
    print_histogram("Score Margin Distribution", margin_bins, correct_reads, mispred_reads,
                    [](const ReadAnalysis& r) { return r.margin; });

    // Damage distribution
    std::vector<std::pair<float, float>> damage_bins = {
        {0, 10}, {10, 20}, {20, 30}, {30, 40}, {40, 50}, {50, 100}
    };
    print_histogram("Damage Distribution", damage_bins, correct_reads, mispred_reads,
                    [](const ReadAnalysis& r) { return r.damage_pct; }, "%");

    // Score distribution
    std::vector<std::pair<float, float>> score_bins = {
        {0, 0.3}, {0.3, 0.5}, {0.5, 0.7}, {0.7, 0.8}, {0.8, 0.9}, {0.9, 1.0}
    };
    print_histogram("Score Distribution", score_bins, correct_reads, mispred_reads,
                    [](const ReadAnalysis& r) { return r.score; });

    // Feature comparison
    std::cout << "=== Feature Comparison (Correct vs Mispredicted) ===" << std::endl;
    std::cout << std::left << std::setw(20) << "Metric"
              << std::right << std::setw(15) << "Correct"
              << std::setw(15) << "Mispred"
              << std::setw(12) << "Diff" << std::endl;
    std::cout << std::string(65, '-') << std::endl;

    auto calc_mean = [](const std::vector<ReadAnalysis>& reads, std::function<float(const ReadAnalysis&)> getter) {
        if (reads.empty()) return 0.0;
        double sum = 0;
        for (const auto& r : reads) sum += getter(r);
        return sum / reads.size();
    };

    std::vector<std::tuple<std::string, std::function<float(const ReadAnalysis&)>>> metrics = {
        {"Length (bp)", [](const ReadAnalysis& r) { return (float)r.length; }},
        {"Score", [](const ReadAnalysis& r) { return r.score; }},
        {"Margin", [](const ReadAnalysis& r) { return r.margin; }},
        {"Damage %", [](const ReadAnalysis& r) { return r.damage_pct; }},
    };

    for (const auto& [name, getter] : metrics) {
        double correct_mean = calc_mean(correct_reads, getter);
        double mispred_mean = calc_mean(mispred_reads, getter);
        double diff = mispred_mean - correct_mean;

        std::cout << std::left << std::setw(20) << name
                  << std::right << std::setw(15) << std::fixed << std::setprecision(2) << correct_mean
                  << std::setw(15) << mispred_mean
                  << std::setw(12) << std::showpos << diff << std::noshowpos << std::endl;
    }
    std::cout << std::endl;

    // Abstention analysis
    std::cout << "=== Abstention Analysis (reject low-margin predictions) ===" << std::endl;
    std::cout << std::left << std::setw(20) << "Margin threshold"
              << std::right << std::setw(12) << "Coverage"
              << std::setw(12) << "Accuracy" << std::endl;
    std::cout << std::string(50, '-') << std::endl;

    std::vector<float> thresholds = {0.0, 0.02, 0.05, 0.10, 0.15, 0.20, 0.30, 0.40, 0.50};
    for (float thresh : thresholds) {
        int accepted_correct = 0, accepted_total = 0;
        for (const auto& r : correct_reads) {
            if (r.margin >= thresh) { accepted_correct++; accepted_total++; }
        }
        for (const auto& r : mispred_reads) {
            if (r.margin >= thresh) accepted_total++;
        }

        if (accepted_total > 0) {
            float coverage = 100.0f * accepted_total / total;
            float acc = 100.0f * accepted_correct / accepted_total;
            std::ostringstream thresh_str;
            thresh_str << ">= " << std::fixed << std::setprecision(2) << thresh;
            std::cout << std::left << std::setw(20) << thresh_str.str()
                      << std::right << std::setw(10) << std::fixed << std::setprecision(1) << coverage << "%"
                      << std::setw(11) << acc << "%" << std::endl;
        }
    }
    std::cout << std::endl;

    // Frame transition analysis
    std::cout << "=== Frame Transitions (Mispredicted) ===" << std::endl;
    std::cout << std::left << std::setw(15) << "True->Pred"
              << std::right << std::setw(12) << "Count"
              << std::setw(10) << "%" << std::endl;
    std::cout << std::string(40, '-') << std::endl;

    std::unordered_map<std::string, int> transitions;
    for (const auto& r : mispred_reads) {
        std::ostringstream key;
        key << r.true_strand << r.true_frame << "->" << r.pred_strand << r.pred_frame;
        transitions[key.str()]++;
    }

    // Sort by count
    std::vector<std::pair<std::string, int>> sorted_trans(transitions.begin(), transitions.end());
    std::sort(sorted_trans.begin(), sorted_trans.end(),
              [](const auto& a, const auto& b) { return a.second > b.second; });

    for (size_t i = 0; i < std::min(size_t(15), sorted_trans.size()); i++) {
        const auto& [key, count] = sorted_trans[i];
        float pct = 100.0f * count / mispred_reads.size();
        std::cout << std::left << std::setw(15) << key
                  << std::right << std::setw(12) << count
                  << std::setw(9) << std::fixed << std::setprecision(1) << pct << "%" << std::endl;
    }

    return 0;
}
