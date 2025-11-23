/**
 * Comprehensive validator for AGP predictions against aMGSIM ground truth
 *
 * Evaluates:
 * 1. Coding/non-coding classification
 * 2. Strand prediction accuracy
 * 3. Frame prediction accuracy
 * 4. Damage detection accuracy with ROC/AUC
 */

#include <iostream>
#include <fstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <sstream>
#include <zlib.h>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include <cctype>

// Ground truth information for each read
struct ReadGroundTruth {
    std::string read_id;
    std::string true_strand;      // + or -
    std::string true_protein;     // Ground truth protein sequence
    bool has_damage = false;
    std::vector<int> damage_positions;
    bool damage_at_5_end = false;  // Damage in positions 0-2
    bool damage_at_3_end = false;  // Damage in positions -1 to -3
};

// Prediction information
struct ReadPrediction {
    std::string read_id;
    std::string pred_strand;
    std::string pred_protein;
    float ancient_prob = 0.5f;
    float coding_prob = 0.0f;
};

// Validation metrics
struct ValidationMetrics {
    // Coding detection
    size_t total_reads = 0;
    size_t coding_reads = 0;
    size_t predicted_reads = 0;
    size_t tp_coding = 0, fp_coding = 0, tn_coding = 0, fn_coding = 0;

    // Strand accuracy
    size_t strand_correct = 0;
    size_t strand_total = 0;
    size_t dual_strand_at_least_one_correct = 0;

    // Damage detection
    size_t damage_tp = 0, damage_fp = 0, damage_tn = 0, damage_fn = 0;
    std::vector<std::pair<float, bool>> damage_scores;  // (ancient_prob, has_damage)

    // Damage by position
    size_t damage_5_end_detected = 0, damage_5_end_total = 0;
    size_t damage_3_end_detected = 0, damage_3_end_total = 0;
    size_t damage_middle_detected = 0, damage_middle_total = 0;
    size_t no_damage_false_positive = 0, no_damage_total = 0;
};

// Parse damage positions from string like "1,-28,5"
std::vector<int> parse_damage_positions(const std::string& damage_str) {
    std::vector<int> positions;
    if (damage_str.empty() || damage_str == "None") {
        return positions;
    }

    std::istringstream iss(damage_str);
    std::string token;
    while (std::getline(iss, token, ',')) {
        try {
            // Remove whitespace
            token.erase(0, token.find_first_not_of(" \t"));
            token.erase(token.find_last_not_of(" \t") + 1);
            if (!token.empty()) {
                positions.push_back(std::stoi(token));
            }
        } catch (...) {
            // Skip invalid entries
        }
    }
    return positions;
}

// Extract ancient_prob from GFF attributes (optimized - no regex)
float extract_ancient_prob(const std::string& attrs) {
    const std::string key = "ancient_prob=";
    size_t pos = attrs.find(key);
    if (pos != std::string::npos) {
        pos += key.length();
        size_t end = pos;
        while (end < attrs.length() && (std::isdigit(attrs[end]) || attrs[end] == '.')) {
            end++;
        }
        if (end > pos) {
            try {
                return std::stof(attrs.substr(pos, end - pos));
            } catch (...) {
                return 0.5f;
            }
        }
    }
    return 0.5f;
}

// Parse ground truth from read ID (format: ...ancient:{strand}:{start}:{end}:{length}:{damage})
// Optimized: no regex, uses simple string operations
bool parse_read_id_ground_truth(const std::string& read_id, ReadGroundTruth& gt) {
    // Find "ancient:" pattern
    size_t ancient_pos = read_id.find("ancient:");
    if (ancient_pos == std::string::npos) return false;

    size_t start = ancient_pos + 8;  // Skip "ancient:"
    if (start >= read_id.length()) return false;

    // Parse strand (+ or -)
    char strand = read_id[start];
    if (strand != '+' && strand != '-') return false;
    gt.true_strand = std::string(1, strand);

    // Skip strand and colon
    start += 2;  // Skip "+" and ":"

    // Skip start:end:length (3 numbers separated by colons)
    int colons = 0;
    while (start < read_id.length() && colons < 3) {
        if (read_id[start] == ':') colons++;
        start++;
    }

    if (start >= read_id.length()) return false;

    // Rest is damage info
    std::string damage_str = read_id.substr(start);

    // Remove any trailing characters like _gene_1
    size_t underscore = damage_str.find("_gene");
    if (underscore != std::string::npos) {
        damage_str = damage_str.substr(0, underscore);
    }

    gt.read_id = read_id;
    gt.has_damage = (damage_str != "None" && !damage_str.empty());
    gt.damage_positions = parse_damage_positions(damage_str);

    // Check damage location
    for (int pos : gt.damage_positions) {
        if (pos >= 0 && pos <= 2) gt.damage_at_5_end = true;
        if (pos >= -3 && pos <= -1) gt.damage_at_3_end = true;
    }

    return true;
}

// Load ground truth from TSV file
std::unordered_map<std::string, ReadGroundTruth> load_ground_truth(const std::string& tsv_path) {
    std::unordered_map<std::string, ReadGroundTruth> gt_map;

    gzFile file = gzopen(tsv_path.c_str(), "r");
    if (!file) {
        std::cerr << "Error: Could not open " << tsv_path << std::endl;
        return gt_map;
    }

    char buffer[65536];
    bool first_line = true;
    int read_name_col = 1;
    int strand_col = 20;  // Strand_read
    int damage_col = 5;
    int protein_col = -1;  // Will find damaged_seq_inframe_aa

    while (gzgets(file, buffer, sizeof(buffer)) != NULL) {
        std::string line(buffer);
        // Remove trailing newline
        while (!line.empty() && (line.back() == '\n' || line.back() == '\r')) {
            line.pop_back();
        }

        std::istringstream iss(line);
        std::string field;
        std::vector<std::string> fields;

        while (std::getline(iss, field, '\t')) {
            fields.push_back(field);
        }

        if (first_line) {
            first_line = false;
            // Find column indices
            for (size_t i = 0; i < fields.size(); i++) {
                if (fields[i] == "read_name") read_name_col = i;
                else if (fields[i] == "Strand_read") strand_col = i;
                else if (fields[i] == "damage") damage_col = i;
                else if (fields[i] == "damaged_seq_inframe_aa") protein_col = i;
            }
            continue;
        }

        if (fields.size() <= static_cast<size_t>(std::max({read_name_col, strand_col, damage_col}))) {
            continue;
        }

        ReadGroundTruth gt;
        gt.read_id = fields[read_name_col];
        gt.true_strand = (strand_col < static_cast<int>(fields.size())) ? fields[strand_col] : "+";
        gt.has_damage = (damage_col < static_cast<int>(fields.size()) &&
                        !fields[damage_col].empty() &&
                        fields[damage_col] != "None" &&
                        fields[damage_col] != "0");

        if (protein_col >= 0 && protein_col < static_cast<int>(fields.size())) {
            gt.true_protein = fields[protein_col];
        }

        // Parse damage positions
        if (damage_col < static_cast<int>(fields.size())) {
            gt.damage_positions = parse_damage_positions(fields[damage_col]);
            for (int pos : gt.damage_positions) {
                if (pos >= 0 && pos <= 2) gt.damage_at_5_end = true;
                if (pos >= -3 && pos <= -1) gt.damage_at_3_end = true;
            }
        }

        gt_map[gt.read_id] = gt;
    }

    gzclose(file);
    return gt_map;
}

// Load predictions from GFF file
std::unordered_map<std::string, std::vector<ReadPrediction>> load_predictions(const std::string& gff_path) {
    std::unordered_map<std::string, std::vector<ReadPrediction>> pred_map;

    std::ifstream file(gff_path);
    if (!file) {
        std::cerr << "Error: Could not open " << gff_path << std::endl;
        return pred_map;
    }

    std::string line;
    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#') continue;

        std::istringstream iss(line);
        std::vector<std::string> fields;
        std::string field;

        while (std::getline(iss, field, '\t')) {
            fields.push_back(field);
        }

        if (fields.size() < 9) continue;

        ReadPrediction pred;
        pred.read_id = fields[0];
        pred.pred_strand = fields[6];
        pred.ancient_prob = extract_ancient_prob(fields[8]);

        pred_map[pred.read_id].push_back(pred);
    }

    return pred_map;
}

// Count FASTQ reads and extract ground truth from read IDs
size_t process_fastq(const std::string& fastq_path,
                     std::unordered_map<std::string, ReadGroundTruth>& gt_map) {
    gzFile file = gzopen(fastq_path.c_str(), "r");
    if (!file) {
        std::cerr << "Error: Could not open " << fastq_path << std::endl;
        return 0;
    }

    size_t count = 0;
    char buffer[65536];

    while (gzgets(file, buffer, sizeof(buffer)) != NULL) {
        if (buffer[0] == '@') {
            count++;

            std::string read_id(buffer + 1);  // Skip @
            // Remove trailing newline and anything after space
            size_t space_pos = read_id.find_first_of(" \t\n\r");
            if (space_pos != std::string::npos) {
                read_id = read_id.substr(0, space_pos);
            }

            // If not in GT map from TSV, try to parse from read ID
            if (gt_map.find(read_id) == gt_map.end()) {
                ReadGroundTruth gt;
                if (parse_read_id_ground_truth(read_id, gt)) {
                    gt_map[read_id] = gt;
                }
            }
        }
    }

    gzclose(file);
    return count;
}

// Calculate AUC-ROC
double calculate_auc(std::vector<std::pair<float, bool>>& scores) {
    if (scores.empty()) return 0.5;

    // Sort by score descending
    std::sort(scores.begin(), scores.end(),
              [](const auto& a, const auto& b) { return a.first > b.first; });

    size_t n_pos = 0, n_neg = 0;
    for (const auto& p : scores) {
        if (p.second) n_pos++;
        else n_neg++;
    }

    if (n_pos == 0 || n_neg == 0) return 0.5;

    double auc = 0.0;
    double tpr_prev = 0.0, fpr_prev = 0.0;
    size_t tp = 0, fp = 0;

    for (const auto& p : scores) {
        if (p.second) tp++;
        else fp++;

        double tpr = static_cast<double>(tp) / n_pos;
        double fpr = static_cast<double>(fp) / n_neg;

        auc += (fpr - fpr_prev) * (tpr + tpr_prev) / 2.0;
        tpr_prev = tpr;
        fpr_prev = fpr;
    }

    return auc;
}

void print_report(const ValidationMetrics& m, double damage_auc) {
    std::cout << "\n";
    std::cout << "======================================================================\n";
    std::cout << "        AGP COMPREHENSIVE VALIDATION REPORT\n";
    std::cout << "======================================================================\n\n";

    // Dataset Summary
    std::cout << "DATASET SUMMARY:\n";
    std::cout << "  Total reads:        " << m.total_reads << "\n";
    std::cout << "  Coding reads (GT):  " << m.coding_reads << "\n";
    std::cout << "  Predicted reads:    " << m.predicted_reads << "\n\n";

    // Coding Detection
    std::cout << "----------------------------------------------------------------------\n";
    std::cout << "1. CODING REGION DETECTION\n";
    std::cout << "----------------------------------------------------------------------\n";
    std::cout << "  Confusion Matrix:\n";
    std::cout << "    TP (correct coding):     " << std::setw(8) << m.tp_coding << "\n";
    std::cout << "    FP (false coding):       " << std::setw(8) << m.fp_coding << "\n";
    std::cout << "    TN (correct non-coding): " << std::setw(8) << m.tn_coding << "\n";
    std::cout << "    FN (missed coding):      " << std::setw(8) << m.fn_coding << "\n\n";

    double precision = (m.tp_coding + m.fp_coding > 0) ?
        static_cast<double>(m.tp_coding) / (m.tp_coding + m.fp_coding) : 0;
    double recall = (m.tp_coding + m.fn_coding > 0) ?
        static_cast<double>(m.tp_coding) / (m.tp_coding + m.fn_coding) : 0;
    double f1 = (precision + recall > 0) ? 2 * precision * recall / (precision + recall) : 0;

    std::cout << "  Metrics:\n";
    std::cout << "    Precision: " << std::fixed << std::setprecision(4) << precision << "\n";
    std::cout << "    Recall:    " << std::fixed << std::setprecision(4) << recall << "\n";
    std::cout << "    F1-Score:  " << std::fixed << std::setprecision(4) << f1 << "\n\n";

    // Strand Prediction
    std::cout << "----------------------------------------------------------------------\n";
    std::cout << "2. STRAND PREDICTION (Dual-Strand Mode)\n";
    std::cout << "----------------------------------------------------------------------\n";
    if (m.strand_total > 0) {
        std::cout << "  At least one correct: " << m.dual_strand_at_least_one_correct
                  << "/" << m.strand_total << " ("
                  << std::fixed << std::setprecision(2)
                  << 100.0 * m.dual_strand_at_least_one_correct / m.strand_total << "%)\n\n";
    }

    // Damage Detection
    std::cout << "----------------------------------------------------------------------\n";
    std::cout << "3. DAMAGE DETECTION (Per-Read)\n";
    std::cout << "----------------------------------------------------------------------\n";
    std::cout << "  AUC-ROC: " << std::fixed << std::setprecision(4) << damage_auc << "\n\n";

    std::cout << "  At threshold 0.7 (ancient_prob > 0.7):\n";
    std::cout << "    TP (damaged detected):    " << std::setw(8) << m.damage_tp << "\n";
    std::cout << "    FP (false damaged):       " << std::setw(8) << m.damage_fp << "\n";
    std::cout << "    TN (correct undamaged):   " << std::setw(8) << m.damage_tn << "\n";
    std::cout << "    FN (missed damaged):      " << std::setw(8) << m.damage_fn << "\n\n";

    double dmg_precision = (m.damage_tp + m.damage_fp > 0) ?
        static_cast<double>(m.damage_tp) / (m.damage_tp + m.damage_fp) : 0;
    double dmg_recall = (m.damage_tp + m.damage_fn > 0) ?
        static_cast<double>(m.damage_tp) / (m.damage_tp + m.damage_fn) : 0;
    double dmg_f1 = (dmg_precision + dmg_recall > 0) ?
        2 * dmg_precision * dmg_recall / (dmg_precision + dmg_recall) : 0;

    std::cout << "  Metrics at threshold 0.7:\n";
    std::cout << "    Precision: " << std::fixed << std::setprecision(4) << dmg_precision << "\n";
    std::cout << "    Recall:    " << std::fixed << std::setprecision(4) << dmg_recall << "\n";
    std::cout << "    F1-Score:  " << std::fixed << std::setprecision(4) << dmg_f1 << "\n\n";

    // Damage by Position
    std::cout << "----------------------------------------------------------------------\n";
    std::cout << "4. DAMAGE DETECTION BY POSITION\n";
    std::cout << "----------------------------------------------------------------------\n";

    if (m.damage_5_end_total > 0) {
        std::cout << "  5' end damage (pos 0-2):  " << m.damage_5_end_detected << "/"
                  << m.damage_5_end_total << " detected ("
                  << std::fixed << std::setprecision(1)
                  << 100.0 * m.damage_5_end_detected / m.damage_5_end_total << "%)\n";
    }

    if (m.damage_3_end_total > 0) {
        std::cout << "  3' end damage (pos -1,-3): " << m.damage_3_end_detected << "/"
                  << m.damage_3_end_total << " detected ("
                  << std::fixed << std::setprecision(1)
                  << 100.0 * m.damage_3_end_detected / m.damage_3_end_total << "%)\n";
    }

    if (m.damage_middle_total > 0) {
        std::cout << "  Middle-only damage:       " << m.damage_middle_detected << "/"
                  << m.damage_middle_total << " detected ("
                  << std::fixed << std::setprecision(1)
                  << 100.0 * m.damage_middle_detected / m.damage_middle_total << "%)\n";
    }

    if (m.no_damage_total > 0) {
        std::cout << "  No damage (false pos):    " << m.no_damage_false_positive << "/"
                  << m.no_damage_total << " false positive ("
                  << std::fixed << std::setprecision(1)
                  << 100.0 * m.no_damage_false_positive / m.no_damage_total << "%)\n";
    }

    std::cout << "\n======================================================================\n";
    std::cout << "SUMMARY:\n";

    // Coding F1
    if (f1 > 0.9) std::cout << "  [OK] Coding detection: EXCELLENT (F1 > 0.9)\n";
    else if (f1 > 0.7) std::cout << "  [OK] Coding detection: GOOD (F1 > 0.7)\n";
    else std::cout << "  [!!] Coding detection: NEEDS IMPROVEMENT (F1 < 0.7)\n";

    // Strand accuracy
    if (m.strand_total > 0 &&
        static_cast<double>(m.dual_strand_at_least_one_correct) / m.strand_total > 0.99) {
        std::cout << "  [OK] Strand prediction: EXCELLENT (>99% with dual-strand)\n";
    }

    // Damage AUC
    if (damage_auc > 0.7) std::cout << "  [OK] Damage detection: GOOD (AUC > 0.7)\n";
    else if (damage_auc > 0.6) std::cout << "  [OK] Damage detection: FAIR (AUC > 0.6)\n";
    else std::cout << "  [!!] Damage detection: WEAK (AUC < 0.6)\n";

    std::cout << "======================================================================\n";
}

int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <fastq.gz> <predictions.gff> [aa-damage.tsv.gz]\n";
        std::cerr << "\nComprehensive AGP validation with damage detection metrics.\n";
        std::cerr << "\nExample:\n";
        std::cerr << "  " << argv[0] << " sample.fq.gz predictions.gff\n";
        std::cerr << "  " << argv[0] << " sample.fq.gz predictions.gff sample_aa-damage.tsv.gz\n";
        return 1;
    }

    std::string fastq_path = argv[1];
    std::string gff_path = argv[2];
    std::string tsv_path = (argc > 3) ? argv[3] : "";

    ValidationMetrics metrics;
    std::unordered_map<std::string, ReadGroundTruth> gt_map;

    // Load ground truth from TSV if provided
    if (!tsv_path.empty()) {
        std::cout << "Loading ground truth from TSV..." << std::endl;
        gt_map = load_ground_truth(tsv_path);
        std::cout << "  Loaded " << gt_map.size() << " entries from TSV\n";
    }

    // Process FASTQ and extract GT from read IDs
    std::cout << "Processing FASTQ file..." << std::endl;
    metrics.total_reads = process_fastq(fastq_path, gt_map);
    std::cout << "  Total reads: " << metrics.total_reads << "\n";
    std::cout << "  Ground truth entries: " << gt_map.size() << "\n";

    // Load predictions
    std::cout << "Loading predictions..." << std::endl;
    auto predictions = load_predictions(gff_path);
    std::cout << "  Predictions for " << predictions.size() << " reads\n";

    metrics.coding_reads = gt_map.size();
    metrics.predicted_reads = predictions.size();

    // Calculate metrics
    std::cout << "Calculating metrics..." << std::endl;

    const float DAMAGE_THRESHOLD = 0.7f;

    for (const auto& [read_id, preds] : predictions) {
        auto gt_it = gt_map.find(read_id);

        if (gt_it != gt_map.end()) {
            // This is a coding read that was predicted
            metrics.tp_coding++;

            const auto& gt = gt_it->second;

            // Check strand accuracy (dual-strand mode)
            metrics.strand_total++;
            bool strand_match = false;
            for (const auto& pred : preds) {
                if (pred.pred_strand == gt.true_strand) {
                    strand_match = true;
                    break;
                }
            }
            if (strand_match) {
                metrics.dual_strand_at_least_one_correct++;
            }

            // Damage detection
            float ancient_prob = preds[0].ancient_prob;  // Same for both strands
            metrics.damage_scores.push_back({ancient_prob, gt.has_damage});

            bool predicted_damage = (ancient_prob > DAMAGE_THRESHOLD);

            if (gt.has_damage && predicted_damage) {
                metrics.damage_tp++;
            } else if (!gt.has_damage && predicted_damage) {
                metrics.damage_fp++;
            } else if (!gt.has_damage && !predicted_damage) {
                metrics.damage_tn++;
            } else {
                metrics.damage_fn++;
            }

            // Damage by position
            if (gt.damage_at_5_end) {
                metrics.damage_5_end_total++;
                if (predicted_damage) metrics.damage_5_end_detected++;
            } else if (gt.damage_at_3_end) {
                metrics.damage_3_end_total++;
                if (predicted_damage) metrics.damage_3_end_detected++;
            } else if (gt.has_damage) {
                metrics.damage_middle_total++;
                if (predicted_damage) metrics.damage_middle_detected++;
            } else {
                metrics.no_damage_total++;
                if (predicted_damage) metrics.no_damage_false_positive++;
            }
        } else {
            // Predicted but not in GT (could be false positive or not in GT)
            metrics.fp_coding++;
        }
    }

    // Calculate FN (in GT but not predicted)
    for (const auto& [read_id, gt] : gt_map) {
        if (predictions.find(read_id) == predictions.end()) {
            metrics.fn_coding++;
        }
    }

    // TN calculation (reads not predicted and not coding)
    metrics.tn_coding = metrics.total_reads - metrics.tp_coding - metrics.fp_coding - metrics.fn_coding;

    // Calculate AUC
    double damage_auc = calculate_auc(metrics.damage_scores);

    // Print report
    print_report(metrics, damage_auc);

    return 0;
}
