/**
 * Enhanced validator for AGP predictions with sequence-level comparison
 *
 * Compares predicted amino acid sequences against ground truth sequences
 * from aMGSIM, including damage correction validation
 */

#include <iostream>
#include <fstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <sstream>
#include <zlib.h>
#include <iomanip>
#include <algorithm>
#include <vector>

struct GroundTruthRead {
    std::string read_id;
    std::string undamaged_aa;       // intersect_seq_inframe_aa - ground truth AA
    std::string damaged_aa;         // damaged_seq_inframe_aa - damaged AA
    std::string undamaged_nt;       // intersect_seq_inframe_nt - ground truth NT
    std::string damaged_nt;         // damaged_seq_inframe_nt - damaged NT
    std::string aa_diffs;           // damage_aaseq_diffs
    std::string codon_diffs;        // damage_codon_diffs
    std::string strand_gene;        // Strand_gene (+/-)
    std::string strand_read;        // Strand_read (+/-)
    int start_intersection = 0;     // Start position
    int end_intersection = 0;       // End position
    int overlap = 0;                // Overlap length
    bool is_coding = true;
};

struct SequenceMetrics {
    size_t total_predictions = 0;
    size_t total_ground_truth = 0;

    // Binary classification metrics
    size_t tp = 0;  // True positives (found coding reads)
    size_t fp = 0;  // False positives (predicted non-coding as coding)
    size_t fn = 0;  // False negatives (missed coding reads)
    size_t tn = 0;  // True negatives (correctly rejected non-coding)

    // Sequence-level metrics (original damaged predictions)
    size_t exact_matches = 0;           // Perfect sequence matches
    size_t partial_matches = 0;         // >80% identity
    size_t with_stop_codons = 0;        // Contains internal stop codons
    size_t stop_codon_corrected = 0;    // Had damage-induced stops

    double total_identity = 0.0;        // Sum of all sequence identities
    size_t sequences_compared = 0;      // Number of sequences compared

    // Corrected sequence metrics
    size_t corrected_exact_matches = 0;    // Exact matches after correction
    size_t corrected_partial_matches = 0;  // >80% identity after correction
    double corrected_total_identity = 0.0; // Sum of corrected identities
    size_t corrected_sequences = 0;        // Number with corrections

    // Frame/strand metrics
    size_t correct_strand = 0;          // Predicted same strand as ground truth
    size_t wrong_strand = 0;            // Predicted opposite strand
    size_t correct_frame = 0;           // Correct reading frame (high identity)
    size_t wrong_frame = 0;             // Wrong reading frame (low identity despite TP)
    size_t frame_shift_detected = 0;    // Likely frame shift (partial match)

    // Damage-related metrics
    size_t damage_detected = 0;         // Reads with damage signatures
    size_t nonsyn_substitutions = 0;    // Non-synonymous changes
    size_t syn_substitutions = 0;       // Synonymous changes
    size_t damage_corrections_helped = 0;   // Correction improved identity
    size_t damage_corrections_hurt = 0;     // Correction decreased identity

    size_t total_reads = 0;
};

// Calculate sequence identity (simple character-by-character comparison)
double calculate_identity(const std::string& seq1, const std::string& seq2) {
    if (seq1.empty() || seq2.empty()) return 0.0;

    size_t matches = 0;
    size_t min_len = std::min(seq1.length(), seq2.length());

    for (size_t i = 0; i < min_len; i++) {
        if (seq1[i] == seq2[i]) matches++;
    }

    // Penalize length differences
    size_t max_len = std::max(seq1.length(), seq2.length());
    return (double)matches / max_len;
}

// Count stop codons in sequence
size_t count_stop_codons(const std::string& seq) {
    return std::count(seq.begin(), seq.end(), '*');
}

// Parse damage info from TSV column (e.g., "0:H>Y,36:V>V")
void parse_damage_diffs(const std::string& diffs,
                        size_t& nonsyn, size_t& syn) {
    nonsyn = 0;
    syn = 0;

    if (diffs.empty() || diffs == "NA") return;

    std::istringstream iss(diffs);
    std::string diff;

    while (std::getline(iss, diff, ',')) {
        // Parse format: "position:original>damaged"
        size_t arrow_pos = diff.find('>');
        size_t colon_pos = diff.find(':');

        if (arrow_pos != std::string::npos && colon_pos != std::string::npos) {
            std::string original = diff.substr(colon_pos + 1, arrow_pos - colon_pos - 1);
            std::string damaged = diff.substr(arrow_pos + 1);

            if (original != damaged) {
                nonsyn++;
            } else {
                syn++;
            }
        }
    }
}

// Load ground truth sequences from TSV
std::unordered_map<std::string, GroundTruthRead> load_ground_truth(
    const std::string& tsv_path) {

    std::unordered_map<std::string, GroundTruthRead> ground_truth;

    gzFile file = gzopen(tsv_path.c_str(), "r");
    if (!file) {
        std::cerr << "Error: Could not open " << tsv_path << std::endl;
        return ground_truth;
    }

    char buffer[65536];
    bool first_line = true;

    while (gzgets(file, buffer, sizeof(buffer)) != NULL) {
        if (first_line) {
            first_line = false;
            continue;  // Skip header
        }

        std::string line(buffer);
        std::istringstream iss(line);
        std::string field;
        int col = 0;

        GroundTruthRead read;

        while (std::getline(iss, field, '\t')) {
            // NOTE: col is 0-indexed, so TSV column 2 is col==1, etc.
            switch (col) {
                case 1:  // read_name (column 2)
                    read.read_id = field;
                    break;
                case 9:  // damage_aaseq_diffs (column 10)
                    read.aa_diffs = field;
                    break;
                case 10: // damage_codon_diffs (column 11)
                    read.codon_diffs = field;
                    break;
                case 11: // Overlap (column 12)
                    read.overlap = std::atoi(field.c_str());
                    break;
                case 14: // intersect_seq_inframe_aa (column 15 - GROUND TRUTH AA)
                    read.undamaged_aa = field;
                    break;
                case 15: // intersect_seq_inframe_nt (column 16 - GROUND TRUTH NT)
                    read.undamaged_nt = field;
                    break;
                case 16: // Start (column 17)
                    read.start_intersection = std::atoi(field.c_str());
                    break;
                case 17: // End (column 18)
                    read.end_intersection = std::atoi(field.c_str());
                    break;
                case 20: // Strand_read (column 21)
                    read.strand_read = field;
                    break;
                case 24: // Strand_gene (column 25)
                    read.strand_gene = field;
                    break;
                case 26: // damaged_seq_inframe_aa (column 27)
                    read.damaged_aa = field;
                    break;
                case 27: // damaged_seq_inframe_nt (column 28)
                    read.damaged_nt = field;
                    break;
            }
            col++;
        }

        if (!read.read_id.empty() && !read.undamaged_aa.empty()) {
            ground_truth[read.read_id] = read;
        }
    }

    gzclose(file);
    return ground_truth;
}

struct PredictedSequence {
    std::string protein;
    std::string corrected_protein;
    std::string dna_sequence;
    bool is_forward = true;         // + strand
    int start = 0;
    int end = 0;
    bool has_correction = false;
};

// Load predicted sequences from FASTA (including corrected versions)
std::unordered_map<std::string, PredictedSequence> load_predictions(
    const std::string& fasta_path) {

    std::unordered_map<std::string, PredictedSequence> predictions;

    std::ifstream file(fasta_path);
    if (!file) {
        std::cerr << "Error: Could not open " << fasta_path << std::endl;
        return predictions;
    }

    std::string line, read_id, sequence;
    bool is_corrected = false;

    while (std::getline(file, line)) {
        if (line.empty()) continue;

        if (line[0] == '>') {
            // Save previous sequence
            if (!read_id.empty()) {
                if (is_corrected) {
                    predictions[read_id].corrected_protein = sequence;
                    predictions[read_id].has_correction = true;
                } else {
                    predictions[read_id].protein = sequence;
                }
            }

            // Parse new header
            // Format: >protein_N READ_ID_gene_N +/- start..end length=Xaa
            // or: >protein_N_corrected READ_ID_corrected ... damage_corrected=true
            std::string header = line.substr(1);

            // Check if this is a corrected sequence
            is_corrected = (header.find("_corrected") != std::string::npos &&
                           header.find("damage_corrected=true") != std::string::npos);

            // Extract strand info
            bool is_forward = true;
            if (header.find(" + ") != std::string::npos) {
                is_forward = true;
            } else if (header.find(" - ") != std::string::npos) {
                is_forward = false;
            }

            size_t first_space = header.find(' ');
            if (first_space != std::string::npos) {
                std::string rest = header.substr(first_space + 1);
                // Read ID ends at "_gene_" marker (AGP format)
                size_t gene_pos = rest.find("_gene_");
                if (gene_pos != std::string::npos) {
                    read_id = rest.substr(0, gene_pos);
                } else {
                    size_t minus_pos = rest.find(" - ");
                    size_t plus_pos = rest.find(" + ");
                    size_t end_pos = std::min(minus_pos, plus_pos);
                    if (end_pos != std::string::npos) {
                        read_id = rest.substr(0, end_pos);
                    } else {
                        read_id = rest;
                    }
                }

                // Store strand info for non-corrected sequences
                if (!is_corrected) {
                    predictions[read_id].is_forward = is_forward;
                }
            } else {
                // FGS format: >READ_ID_START_END_STRAND (no space)
                // Example: >...None_2_73_+ means gene from position 2 to 73 on + strand
                // Need to strip trailing _START_END_[+-] suffix
                read_id = header;

                // Check for FGS suffix pattern: _NUMBER_NUMBER_[+-] at end
                size_t len = read_id.length();
                if (len > 5) {
                    char last = read_id[len - 1];
                    if (last == '+' || last == '-') {
                        // Pattern: _START_END_STRAND where START/END are numbers
                        // Example suffix: _2_73_+
                        // We need to find: underscore before '+', underscore before '73', underscore before '2'

                        // Find underscore right before strand (the _ before +/-)
                        size_t strand_under = read_id.rfind('_');  // Last underscore
                        if (strand_under != std::string::npos && strand_under > 1) {
                            // Find underscore before END number
                            size_t end_under = read_id.rfind('_', strand_under - 1);
                            if (end_under != std::string::npos && end_under > 0) {
                                // Extract END number (between end_under and strand_under)
                                std::string end_str = read_id.substr(end_under + 1, strand_under - end_under - 1);

                                // Find underscore before START number
                                size_t start_under = read_id.rfind('_', end_under - 1);
                                if (start_under != std::string::npos) {
                                    // Extract START number (between start_under and end_under)
                                    std::string start_str = read_id.substr(start_under + 1, end_under - start_under - 1);

                                    bool start_numeric = !start_str.empty() && std::all_of(start_str.begin(), start_str.end(), ::isdigit);
                                    bool end_numeric = !end_str.empty() && std::all_of(end_str.begin(), end_str.end(), ::isdigit);

                                    if (start_numeric && end_numeric) {
                                        // Strip FGS suffix: _START_END_STRAND
                                        read_id = read_id.substr(0, start_under);
                                        // Extract strand from FGS suffix
                                        is_forward = (last == '+');
                                    }
                                }
                            }
                        }
                    }
                }
                predictions[read_id].is_forward = is_forward;
            }
            sequence.clear();
        } else {
            sequence += line;
        }
    }

    // Save last sequence
    if (!read_id.empty()) {
        if (is_corrected) {
            predictions[read_id].corrected_protein = sequence;
            predictions[read_id].has_correction = true;
        } else {
            predictions[read_id].protein = sequence;
        }
    }

    return predictions;
}

// Count total reads
size_t count_total_reads(const std::string& fastq_path) {
    gzFile file = gzopen(fastq_path.c_str(), "r");
    if (!file) {
        std::cerr << "Error: Could not open " << fastq_path << std::endl;
        return 0;
    }

    size_t count = 0;
    char buffer[65536];

    while (gzgets(file, buffer, sizeof(buffer)) != NULL) {
        if (buffer[0] == '@') count++;
    }

    gzclose(file);
    return count;
}

// Compare sequences and calculate metrics
void compare_sequences(
    const std::unordered_map<std::string, GroundTruthRead>& ground_truth,
    const std::unordered_map<std::string, PredictedSequence>& predictions,
    SequenceMetrics& metrics) {

    metrics.total_ground_truth = ground_truth.size();
    metrics.total_predictions = predictions.size();

    // For each ground truth coding sequence
    for (const auto& [read_id, gt_read] : ground_truth) {
        auto pred_it = predictions.find(read_id);

        if (pred_it != predictions.end()) {
            // True positive: Found the coding read
            metrics.tp++;

            const std::string& predicted_seq = pred_it->second.protein;
            const std::string& truth_seq = gt_read.undamaged_aa;

            // Calculate sequence identity for original (damaged) prediction
            double identity = calculate_identity(predicted_seq, truth_seq);
            metrics.total_identity += identity;
            metrics.sequences_compared++;

            if (identity > 0.999) {
                metrics.exact_matches++;
            } else if (identity > 0.80) {
                metrics.partial_matches++;
            }

            // Check for stop codons
            size_t stops = count_stop_codons(predicted_seq);
            if (stops > 0) {
                metrics.with_stop_codons++;
            }

            // Check strand correctness
            bool gt_forward = (gt_read.strand_gene == "+");
            bool pred_forward = pred_it->second.is_forward;
            if (gt_forward == pred_forward) {
                metrics.correct_strand++;
            } else {
                metrics.wrong_strand++;
            }

            // Determine frame correctness based on identity
            // High identity (>50%) suggests correct frame
            // Low identity (<20%) suggests wrong frame
            if (identity > 0.50) {
                metrics.correct_frame++;
            } else if (identity < 0.20) {
                metrics.wrong_frame++;
            } else {
                // Intermediate - possible frame shift
                metrics.frame_shift_detected++;
            }

            // Compare corrected sequence if available
            if (pred_it->second.has_correction && !pred_it->second.corrected_protein.empty()) {
                const std::string& corrected_seq = pred_it->second.corrected_protein;
                double corrected_identity = calculate_identity(corrected_seq, truth_seq);
                metrics.corrected_total_identity += corrected_identity;
                metrics.corrected_sequences++;

                // Track if correction helped or hurt
                if (corrected_identity > identity + 0.01) {
                    metrics.damage_corrections_helped++;
                } else if (corrected_identity < identity - 0.01) {
                    metrics.damage_corrections_hurt++;
                }

                if (corrected_identity > 0.999) {
                    metrics.corrected_exact_matches++;
                } else if (corrected_identity > 0.80) {
                    metrics.corrected_partial_matches++;
                }
            }

            // Check if damage was present
            if (!gt_read.aa_diffs.empty() && gt_read.aa_diffs != "NA") {
                metrics.damage_detected++;

                // Count damaged stops that might be corrected
                size_t damaged_stops = count_stop_codons(gt_read.damaged_aa);
                size_t predicted_stops = count_stop_codons(predicted_seq);

                if (damaged_stops > predicted_stops) {
                    metrics.stop_codon_corrected++;
                }
            }

            // Parse damage differences
            size_t nonsyn, syn;
            parse_damage_diffs(gt_read.aa_diffs, nonsyn, syn);
            metrics.nonsyn_substitutions += nonsyn;
            metrics.syn_substitutions += syn;

        } else {
            // False negative: Missed a coding read
            metrics.fn++;
        }
    }

    // For each prediction, check if it's a false positive
    for (const auto& [read_id, pred] : predictions) {
        if (ground_truth.find(read_id) == ground_truth.end()) {
            metrics.fp++;
        }
    }

    // Calculate true negatives
    size_t non_coding = metrics.total_reads - metrics.total_ground_truth;
    metrics.tn = non_coding - metrics.fp;
}

void print_report(const SequenceMetrics& metrics) {
    std::cout << "========================================================================" << std::endl;
    std::cout << "  AGP Sequence-Level Validation Report" << std::endl;
    std::cout << "========================================================================" << std::endl;
    std::cout << std::endl;

    std::cout << "DATASET SUMMARY:" << std::endl;
    std::cout << "  Total reads:           " << metrics.total_reads << std::endl;
    std::cout << "  Ground truth (coding): " << metrics.total_ground_truth << std::endl;
    std::cout << "  Predicted genes:       " << metrics.total_predictions << std::endl;
    std::cout << std::endl;

    std::cout << "BINARY CLASSIFICATION:" << std::endl;
    std::cout << "  True Positives (TP):   " << metrics.tp
              << "  (coding reads found)" << std::endl;
    std::cout << "  False Negatives (FN):  " << metrics.fn
              << "  (coding reads missed)" << std::endl;
    std::cout << "  False Positives (FP):  " << metrics.fp
              << "  (incorrect predictions)" << std::endl;
    std::cout << "  True Negatives (TN):   " << metrics.tn
              << "  (non-coding correctly rejected)" << std::endl;
    std::cout << std::endl;

    double precision = (double)metrics.tp / (metrics.tp + metrics.fp);
    double recall = (double)metrics.tp / (metrics.tp + metrics.fn);
    double f1 = 2.0 * precision * recall / (precision + recall);
    double accuracy = (double)(metrics.tp + metrics.tn) / metrics.total_reads;

    std::cout << "CLASSIFICATION METRICS:" << std::endl;
    std::cout << "  Precision:  " << std::fixed << std::setprecision(4) << precision << std::endl;
    std::cout << "  Recall:     " << std::fixed << std::setprecision(4) << recall << std::endl;
    std::cout << "  F1-Score:   " << std::fixed << std::setprecision(4) << f1 << std::endl;
    std::cout << "  Accuracy:   " << std::fixed << std::setprecision(4) << accuracy << std::endl;
    std::cout << std::endl;

    std::cout << "SEQUENCE-LEVEL METRICS:" << std::endl;
    std::cout << "  Sequences compared:    " << metrics.sequences_compared << std::endl;

    if (metrics.sequences_compared > 0) {
        double avg_identity = metrics.total_identity / metrics.sequences_compared;
        double exact_rate = 100.0 * metrics.exact_matches / metrics.sequences_compared;
        double partial_rate = 100.0 * metrics.partial_matches / metrics.sequences_compared;

        std::cout << "  Average identity:      " << std::fixed << std::setprecision(2)
                  << (avg_identity * 100) << "%" << std::endl;
        std::cout << "  Exact matches:         " << metrics.exact_matches
                  << " (" << std::fixed << std::setprecision(1) << exact_rate << "%)" << std::endl;
        std::cout << "  Partial matches (>80%): " << metrics.partial_matches
                  << " (" << std::fixed << std::setprecision(1) << partial_rate << "%)" << std::endl;
    }
    std::cout << std::endl;

    // Corrected sequence metrics
    if (metrics.corrected_sequences > 0) {
        std::cout << "DAMAGE-CORRECTED SEQUENCE METRICS:" << std::endl;
        double corr_avg_identity = metrics.corrected_total_identity / metrics.corrected_sequences;
        double corr_exact_rate = 100.0 * metrics.corrected_exact_matches / metrics.corrected_sequences;
        double corr_partial_rate = 100.0 * metrics.corrected_partial_matches / metrics.corrected_sequences;

        std::cout << "  Corrected sequences:   " << metrics.corrected_sequences << std::endl;
        std::cout << "  Average identity:      " << std::fixed << std::setprecision(2)
                  << (corr_avg_identity * 100) << "%" << std::endl;
        std::cout << "  Exact matches:         " << metrics.corrected_exact_matches
                  << " (" << std::fixed << std::setprecision(1) << corr_exact_rate << "%)" << std::endl;
        std::cout << "  Partial matches (>80%): " << metrics.corrected_partial_matches
                  << " (" << std::fixed << std::setprecision(1) << corr_partial_rate << "%)" << std::endl;
        std::cout << std::endl;
    }

    // Frame/strand metrics
    if (metrics.tp > 0) {
        std::cout << "FRAME/STRAND ANALYSIS:" << std::endl;
        double correct_strand_rate = 100.0 * metrics.correct_strand / metrics.tp;
        double correct_frame_rate = 100.0 * metrics.correct_frame / metrics.tp;
        double wrong_frame_rate = 100.0 * metrics.wrong_frame / metrics.tp;
        double frameshift_rate = 100.0 * metrics.frame_shift_detected / metrics.tp;

        std::cout << "  Correct strand:        " << metrics.correct_strand
                  << " (" << std::fixed << std::setprecision(1) << correct_strand_rate << "%)" << std::endl;
        std::cout << "  Wrong strand:          " << metrics.wrong_strand
                  << " (" << std::fixed << std::setprecision(1) << (100.0 - correct_strand_rate) << "%)" << std::endl;
        std::cout << "  Correct frame (>50% id): " << metrics.correct_frame
                  << " (" << std::fixed << std::setprecision(1) << correct_frame_rate << "%)" << std::endl;
        std::cout << "  Wrong frame (<20% id):   " << metrics.wrong_frame
                  << " (" << std::fixed << std::setprecision(1) << wrong_frame_rate << "%)" << std::endl;
        std::cout << "  Possible frameshift:     " << metrics.frame_shift_detected
                  << " (" << std::fixed << std::setprecision(1) << frameshift_rate << "%)" << std::endl;
        std::cout << std::endl;
    }

    std::cout << "DAMAGE-RELATED METRICS:" << std::endl;
    std::cout << "  Reads with damage:     " << metrics.damage_detected << std::endl;
    std::cout << "  Non-syn substitutions: " << metrics.nonsyn_substitutions << std::endl;
    std::cout << "  Syn substitutions:     " << metrics.syn_substitutions << std::endl;
    std::cout << "  Internal stop codons:  " << metrics.with_stop_codons << std::endl;
    std::cout << "  Stop codons corrected: " << metrics.stop_codon_corrected << std::endl;
    if (metrics.corrected_sequences > 0) {
        std::cout << "  Corrections helped:    " << metrics.damage_corrections_helped
                  << " (" << std::fixed << std::setprecision(1)
                  << (100.0 * metrics.damage_corrections_helped / metrics.corrected_sequences) << "%)" << std::endl;
        std::cout << "  Corrections hurt:      " << metrics.damage_corrections_hurt
                  << " (" << std::fixed << std::setprecision(1)
                  << (100.0 * metrics.damage_corrections_hurt / metrics.corrected_sequences) << "%)" << std::endl;
    }
    std::cout << std::endl;

    std::cout << "========================================================================" << std::endl;
    std::cout << "ASSESSMENT:" << std::endl;

    if (metrics.sequences_compared > 0) {
        double avg_identity = metrics.total_identity / metrics.sequences_compared;

        if (avg_identity > 0.95) {
            std::cout << "✓ EXCELLENT: Average sequence identity > 95%" << std::endl;
        } else if (avg_identity > 0.85) {
            std::cout << "✓ GOOD: Average sequence identity > 85%" << std::endl;
        } else if (avg_identity > 0.75) {
            std::cout << "⚠ FAIR: Average sequence identity > 75%" << std::endl;
        } else {
            std::cout << "✗ POOR: Average sequence identity < 75%" << std::endl;
        }

        double exact_rate = 100.0 * metrics.exact_matches / metrics.sequences_compared;
        if (exact_rate > 80) {
            std::cout << "✓ High exact match rate (>80%)" << std::endl;
        } else if (exact_rate > 60) {
            std::cout << "⚠ Moderate exact match rate (60-80%)" << std::endl;
        } else {
            std::cout << "✗ Low exact match rate (<60%)" << std::endl;
        }
    }

    if (f1 > 0.9) {
        std::cout << "✓ EXCELLENT: F1-score > 0.9" << std::endl;
    } else if (f1 > 0.8) {
        std::cout << "✓ GOOD: F1-score > 0.8" << std::endl;
    }

    std::cout << "========================================================================" << std::endl;
}

int main(int argc, char* argv[]) {
    if (argc != 4) {
        std::cerr << "Usage: " << argv[0] << " <fastq.gz> <aa-damage.tsv.gz> <predictions.fasta>" << std::endl;
        std::cerr << std::endl;
        std::cerr << "Example:" << std::endl;
        std::cerr << "  " << argv[0] << " reads.fq.gz ground_truth.tsv.gz genes.fna" << std::endl;
        std::cerr << std::endl;
        std::cerr << "This validator compares predicted amino acid sequences to ground truth." << std::endl;
        return 1;
    }

    std::string fastq_path = argv[1];
    std::string tsv_path = argv[2];
    std::string fasta_path = argv[3];

    SequenceMetrics metrics;

    std::cout << "Loading ground truth sequences from TSV..." << std::endl;
    auto ground_truth = load_ground_truth(tsv_path);
    std::cout << "  Loaded " << ground_truth.size() << " ground truth sequences" << std::endl;

    std::cout << "Loading predicted sequences from FASTA..." << std::endl;
    auto predictions = load_predictions(fasta_path);
    std::cout << "  Loaded " << predictions.size() << " predictions" << std::endl;

    std::cout << "Counting total reads..." << std::endl;
    metrics.total_reads = count_total_reads(fastq_path);
    std::cout << "  Total reads: " << metrics.total_reads << std::endl;

    std::cout << "Comparing sequences..." << std::endl;
    compare_sequences(ground_truth, predictions, metrics);

    std::cout << std::endl;
    print_report(metrics);

    return 0;
}
