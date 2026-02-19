#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <chrono>
#include <zlib.h>

#include "agp/bayesian_frame.hpp"

// Simple gzip reader
std::string read_gzip_file(const std::string& path) {
    gzFile file = gzopen(path.c_str(), "rb");
    if (!file) {
        throw std::runtime_error("Cannot open: " + path);
    }

    std::string content;
    char buffer[8192];
    int bytes_read;
    while ((bytes_read = gzread(file, buffer, sizeof(buffer))) > 0) {
        content.append(buffer, bytes_read);
    }
    gzclose(file);
    return content;
}

struct GroundTruth {
    std::string read_id;
    std::string damaged_seq;
    int true_frame;      // 0, 1, 2
    bool is_reverse;     // True if reverse strand
    std::string true_protein;
};

// Parse ground truth from TSV
std::vector<GroundTruth> load_ground_truth(const std::string& tsv_content, int max_reads = -1) {
    std::vector<GroundTruth> records;
    std::istringstream iss(tsv_content);
    std::string line;

    // Read header
    std::getline(iss, line);
    std::map<std::string, int> col_idx;
    {
        std::istringstream header_ss(line);
        std::string col;
        int idx = 0;
        while (std::getline(header_ss, col, '\t')) {
            col_idx[col] = idx++;
        }
    }

    // Required columns
    int read_name_col = col_idx["read_name"];
    int damaged_seq_col = col_idx["damaged_seq"];
    int strand_read_col = col_idx["Strand_read"];
    int strand_gene_col = col_idx["Strand_gene"];
    int protein_col = col_idx["damaged_seq_inframe_aa"];

    while (std::getline(iss, line)) {
        if (max_reads > 0 && (int)records.size() >= max_reads) break;

        std::vector<std::string> fields;
        std::istringstream line_ss(line);
        std::string field;
        while (std::getline(line_ss, field, '\t')) {
            fields.push_back(field);
        }

        if (fields.size() <= (size_t)std::max({read_name_col, damaged_seq_col, strand_read_col, strand_gene_col, protein_col})) {
            continue;
        }

        GroundTruth gt;
        gt.read_id = fields[read_name_col];
        gt.damaged_seq = fields[damaged_seq_col];
        gt.true_protein = fields[protein_col];

        // Determine if we need to RC the sequence for translation
        // is_reverse = true when Strand_read != Strand_gene (need to RC to get coding sequence)
        std::string strand_read = fields[strand_read_col];
        std::string strand_gene = fields[strand_gene_col];
        gt.is_reverse = (strand_read != strand_gene);

        // Frame is encoded in the read metadata - extract from position info
        // For simplicity, use frame 0 and rely on protein comparison
        gt.true_frame = 0;

        // Skip reads without protein or damaged sequence
        if (gt.damaged_seq.empty() || gt.true_protein.empty()) {
            continue;
        }

        // Skip very short reads
        if (gt.damaged_seq.length() < 30) {
            continue;
        }

        records.push_back(gt);
    }

    return records;
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <ground_truth.tsv.gz> [max_reads]\n";
        return 1;
    }

    std::string tsv_path = argv[1];
    int max_reads = (argc > 2) ? std::stoi(argv[2]) : 100000;

    std::cout << "=== Bayesian Frame Selector Test ===" << std::endl;
    std::cout << "Loading ground truth from: " << tsv_path << std::endl;

    auto start = std::chrono::high_resolution_clock::now();

    std::string tsv_content = read_gzip_file(tsv_path);
    std::vector<GroundTruth> ground_truth = load_ground_truth(tsv_content, max_reads);

    auto load_end = std::chrono::high_resolution_clock::now();
    auto load_ms = std::chrono::duration_cast<std::chrono::milliseconds>(load_end - start).count();

    std::cout << "Loaded " << ground_truth.size() << " reads with ground truth (" << load_ms << "ms)" << std::endl;

    // Create Bayesian frame selector with damage parameters for high-damage sample
    agp::DamageChannel damage(0.30f, 0.30f, 0.51f, 0.51f);  // 30% D_max, λ=0.51
    agp::BayesianFrameSelector selector(damage);

    // Evaluation metrics
    int total_reads = 0;
    int correct_strand = 0;
    int correct_protein_match = 0;
    int frames_emitted_0 = 0;
    int frames_emitted_1 = 0;
    int frames_emitted_2 = 0;
    int high_confidence_correct = 0;
    int high_confidence_total = 0;

    // Length bins for analysis
    std::map<std::string, std::pair<int, int>> length_bins;  // bin -> (correct, total)

    std::cout << "\nProcessing reads..." << std::endl;

    auto process_start = std::chrono::high_resolution_clock::now();

    for (const auto& gt : ground_truth) {
        auto result = selector.select_frames(gt.damaged_seq);
        total_reads++;

        // Determine length bin
        std::string bin;
        size_t len = gt.damaged_seq.length();
        if (len < 50) bin = "30-50bp";
        else if (len < 75) bin = "50-75bp";
        else if (len < 100) bin = "75-100bp";
        else if (len < 130) bin = "100-130bp";
        else if (len < 160) bin = "130-160bp";
        else bin = "160-200bp";

        if (length_bins.find(bin) == length_bins.end()) {
            length_bins[bin] = {0, 0};
        }
        length_bins[bin].second++;

        // Count frames emitted
        if (result.selected_frames.empty()) {
            frames_emitted_0++;
            continue;
        } else if (result.selected_frames.size() == 1) {
            frames_emitted_1++;
        } else {
            frames_emitted_2++;
        }

        // Check if any emitted frame matches ground truth strand
        bool strand_match = false;
        bool protein_match = false;

        for (const auto& fr : result.selected_frames) {
            if (fr.is_reverse == gt.is_reverse) {
                strand_match = true;

                // Check protein similarity (allow some mismatches due to damage)
                // Simple check: see if predicted protein contains truth or vice versa
                if (fr.protein.find(gt.true_protein) != std::string::npos ||
                    gt.true_protein.find(fr.protein) != std::string::npos ||
                    fr.protein == gt.true_protein) {
                    protein_match = true;
                }

                // Also check substring match for partial overlaps
                size_t min_len = std::min(fr.protein.length(), gt.true_protein.length());
                if (min_len >= 10) {
                    int matches = 0;
                    for (size_t i = 0; i < min_len; i++) {
                        if (fr.protein[i] == gt.true_protein[i]) matches++;
                    }
                    if (matches >= (int)(min_len * 0.8)) {
                        protein_match = true;
                    }
                }
            }
        }

        if (strand_match) {
            correct_strand++;
            length_bins[bin].first++;
        }
        if (protein_match) {
            correct_protein_match++;
        }

        // Track high-confidence predictions
        if (result.confidence > 0.3f && result.selected_frames.size() == 1) {
            high_confidence_total++;
            if (strand_match) {
                high_confidence_correct++;
            }
        }

        // Progress
        if (total_reads % 10000 == 0) {
            std::cout << "  Processed " << total_reads << " reads..." << std::endl;
        }
    }

    auto process_end = std::chrono::high_resolution_clock::now();
    auto process_ms = std::chrono::duration_cast<std::chrono::milliseconds>(process_end - process_start).count();

    // Report results
    std::cout << "\n=== RESULTS ===" << std::endl;
    std::cout << "Total reads processed: " << total_reads << std::endl;
    std::cout << "Processing time: " << process_ms << "ms (" << (1000.0 * total_reads / process_ms) << " reads/sec)" << std::endl;

    std::cout << "\n--- Frame Emission ---" << std::endl;
    std::cout << "  0 frames (noncoding/uncertain): " << frames_emitted_0 << " (" << (100.0 * frames_emitted_0 / total_reads) << "%)" << std::endl;
    std::cout << "  1 frame (confident): " << frames_emitted_1 << " (" << (100.0 * frames_emitted_1 / total_reads) << "%)" << std::endl;
    std::cout << "  2 frames (ambiguous): " << frames_emitted_2 << " (" << (100.0 * frames_emitted_2 / total_reads) << "%)" << std::endl;

    float avg_frames = (float)(frames_emitted_1 + 2 * frames_emitted_2) / total_reads;
    std::cout << "  Average frames/read: " << avg_frames << std::endl;

    int reads_with_output = frames_emitted_1 + frames_emitted_2;
    std::cout << "\n--- Accuracy (reads with output) ---" << std::endl;
    std::cout << "  Correct strand: " << correct_strand << "/" << reads_with_output
              << " (" << (100.0 * correct_strand / reads_with_output) << "%)" << std::endl;
    std::cout << "  Protein match: " << correct_protein_match << "/" << reads_with_output
              << " (" << (100.0 * correct_protein_match / reads_with_output) << "%)" << std::endl;

    if (high_confidence_total > 0) {
        std::cout << "\n--- High Confidence (margin > 0.3) ---" << std::endl;
        std::cout << "  Correct: " << high_confidence_correct << "/" << high_confidence_total
                  << " (" << (100.0 * high_confidence_correct / high_confidence_total) << "%)" << std::endl;
    }

    std::cout << "\n--- Accuracy by Read Length ---" << std::endl;
    for (const auto& kv : length_bins) {
        float acc = (kv.second.second > 0) ? (100.0f * kv.second.first / kv.second.second) : 0;
        std::cout << "  " << kv.first << ": " << kv.second.first << "/" << kv.second.second
                  << " (" << acc << "%)" << std::endl;
    }

    std::cout << "\n=== COMPARISON TARGET ===" << std::endl;
    std::cout << "  6-frame translation: ~72% recall, ~59% precision, F1=65.2%" << std::endl;
    std::cout << "  FGS: ~22% recall, ~64% precision, F1=32.5%" << std::endl;
    std::cout << "  Goal: Similar recall to 6-frame with fewer ORFs/read" << std::endl;

    return 0;
}
