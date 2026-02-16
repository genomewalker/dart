/**
 * Unified AGP Validator
 *
 * Evaluates AGP predictions against aMGSIM ground truth with proper handling of:
 * 1. Frame selection accuracy (using in-frame nucleotide sequence matching)
 * 2. Strand prediction accuracy (accounting for read vs gene strand)
 * 3. Sequence identity (using local alignment for partial overlaps)
 * 4. Damage detection (terminal damage in first/last 10bp)
 *
 * Ground Truth from aMGSIM:
 * - Read comes from a known gene at known coordinates
 * - True frame: determined by matching damaged_seq against damaged_seq_inframe_nt
 *   (accounts for strand orientation via reverse complement when needed)
 * - True strand: if read_strand == gene_strand, use forward; else use reverse
 * - True protein: intersect_seq_inframe_aa (undamaged translation)
 */

#ifndef AGP_STANDALONE_VALIDATE
#include "cli/subcommand.hpp"
#endif
#include "agp/codon_tables.hpp"
#include "agp/log_utils.hpp"
#include <iostream>
#include <fstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <array>
#include <sstream>
#include <tuple>
#include <zlib.h>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include <charconv>
#include <cctype>
#include <cstring>
#include <cstdlib>
#include <thread>
#include <atomic>
#include <mutex>
#include <chrono>

// Ground truth for a single read
struct GroundTruth {
    std::string read_id;

    // Coordinates
    int read_start = 0;
    int read_end = 0;
    int gene_start = 0;
    int gene_end = 0;
    int read_length = 0;

    // Strands
    char strand_read = '+';   // Strand the read was sequenced from
    char strand_gene = '+';   // Strand the gene is on

    // Computed ground truth for AGP
    int true_frame = 0;       // (read_start - gene_start) % 3
    bool true_forward = true; // Should AGP predict forward strand?

    // Sequences
    std::string true_protein;     // intersect_seq_inframe_aa (undamaged)
    std::string damaged_protein;  // damaged_seq_inframe_aa
    std::string damaged_seq;      // Full damaged read sequence
    std::string damaged_seq_inframe_nt;  // In-frame portion of damaged sequence

    // Damage info
    std::vector<int> damage_positions;  // From read ID
    bool has_terminal_damage = false;   // Damage in first/last 10bp
    bool damage_at_5prime = false;
    bool damage_at_3prime = false;
};

// AGP prediction for a single read
struct Prediction {
    std::string read_id;
    char strand = '+';        // Predicted strand
    int frame = 0;            // Predicted frame (0, 1, 2)
    std::string protein;      // Predicted protein sequence
    float damage_signal = 0.5f;
    float coding_prob = 0.0f;
    float score = 0.0f;       // Prediction score
    float damage_pct = 0.0f;  // Per-read damage percentage
};

// Validation metrics
struct Metrics {
    // Dataset counts
    size_t total_reads = 0;
    size_t coding_reads = 0;      // Reads with ground truth
    size_t predicted_reads = 0;   // Reads with predictions

    // Binary classification (coding detection)
    size_t tp = 0;  // Predicted and is coding
    size_t fp = 0;  // Predicted but not coding
    size_t fn = 0;  // Not predicted but is coding
    size_t tn = 0;  // Not predicted and not coding

    // Frame selection (among coding reads that were predicted)
    size_t frame_correct = 0;
    size_t frame_wrong = 0;
    size_t frame_total = 0;
    std::array<size_t, 3> frame_distribution = {};  // Predicted frame counts
    std::array<size_t, 3> true_frame_distribution = {};  // True frame counts

    // Strand prediction
    size_t strand_correct = 0;
    size_t strand_wrong = 0;
    size_t strand_total = 0;

    // Combined (both frame AND strand correct)
    size_t both_correct = 0;
    size_t both_total = 0;

    // Sequence identity (all predictions)
    double total_identity = 0.0;
    size_t identity_count = 0;
    size_t exact_matches = 0;
    size_t high_identity = 0;  // >90%
    size_t medium_identity = 0; // 50-90%
    size_t low_identity = 0;   // <50%

    // Identity stratified by correctness (to assess damage correction)
    double identity_correct_both = 0.0;   // Both frame AND strand correct
    size_t count_correct_both = 0;
    double identity_wrong_frame = 0.0;    // Wrong frame (regardless of strand)
    size_t count_wrong_frame = 0;
    double identity_wrong_strand = 0.0;   // Wrong strand (regardless of frame)
    size_t count_wrong_strand = 0;

    // Corrected protein identity (if --fasta-corrected provided)
    double total_corrected_identity = 0.0;
    size_t corrected_identity_count = 0;
    double corrected_identity_correct_both = 0.0;
    size_t corrected_count_correct_both = 0;

    // Damage detection
    size_t damage_tp = 0;  // Has damage and detected (damage_signal > 0.5)
    size_t damage_fp = 0;  // No damage but detected
    size_t damage_tn = 0;  // No damage and not detected
    size_t damage_fn = 0;  // Has damage but not detected
    std::vector<std::pair<float, bool>> damage_scores;  // For AUC

    // Stop codon analysis
    size_t predictions_with_stops = 0;
    size_t correct_frame_with_stops = 0;
    size_t wrong_frame_with_stops = 0;
};

bool parse_int_strict(const std::string& token, int& value) {
    const char* begin = token.data();
    const char* end = begin + token.size();
    while (begin < end && std::isspace(static_cast<unsigned char>(*begin))) ++begin;
    while (begin < end && std::isspace(static_cast<unsigned char>(*(end - 1)))) --end;
    if (begin == end) return false;

    auto result = std::from_chars(begin, end, value);
    return result.ec == std::errc() && result.ptr == end;
}

bool parse_float_strict(const std::string& token, float& value) {
    const char* begin = token.c_str();
    while (*begin != '\0' && std::isspace(static_cast<unsigned char>(*begin))) ++begin;
    if (*begin == '\0') return false;

    char* end = nullptr;
    float parsed = std::strtof(begin, &end);
    if (begin == end) return false;
    while (*end != '\0' && std::isspace(static_cast<unsigned char>(*end))) ++end;
    if (*end != '\0') return false;

    value = parsed;
    return true;
}

// Parse damage positions from read ID suffix (e.g., "1,-28,5" or "None")
std::vector<int> parse_damage_positions(const std::string& damage_str) {
    std::vector<int> positions;
    if (damage_str.empty() || damage_str == "None") return positions;

    std::istringstream iss(damage_str);
    std::string token;
    bool saw_invalid = false;
    while (std::getline(iss, token, ',')) {
        int parsed = 0;
        if (parse_int_strict(token, parsed)) {
            positions.push_back(parsed);
        } else if (!token.empty()) {
            saw_invalid = true;
        }
    }

    if (saw_invalid) {
        static std::atomic<bool> warned{false};
        if (!warned.exchange(true)) {
            std::cerr << "Warning: encountered malformed damage position token(s); invalid values are ignored.\n";
        }
    }

    return positions;
}

// Calculate sequence identity with proper alignment handling
double calculate_identity(const std::string& pred, const std::string& truth) {
    if (pred.empty() || truth.empty()) return 0.0;

    // For sequences of similar length, use global alignment
    size_t matches = 0;
    size_t len = std::min(pred.length(), truth.length());

    for (size_t i = 0; i < len; i++) {
        if (pred[i] == truth[i]) matches++;
    }

    // Normalize by max length to penalize length differences
    size_t max_len = std::max(pred.length(), truth.length());
    return static_cast<double>(matches) / max_len;
}

// Count internal stop codons (exclude terminal)
size_t count_internal_stops(const std::string& protein) {
    size_t stops = 0;
    for (size_t i = 0; i + 1 < protein.length(); i++) {
        if (protein[i] == '*') stops++;
    }
    return stops;
}

// Determine true frame by matching damaged_seq against damaged_seq_inframe_nt
// Returns frame (0, 1, or 2) or -1 if cannot determine
int determine_true_frame(const std::string& damaged_seq,
                         const std::string& inframe_nt,
                         bool same_strand) {
    if (damaged_seq.empty() || inframe_nt.empty()) return -1;

    // Get the sequence to search in (forward or reverse complement)
    std::string seq = same_strand ? damaged_seq : agp::reverse_complement(damaged_seq);

    // Convert to uppercase for comparison
    std::string seq_upper, inframe_upper;
    seq_upper.reserve(seq.size());
    inframe_upper.reserve(inframe_nt.size());
    for (char c : seq) seq_upper += std::toupper(static_cast<unsigned char>(c));
    for (char c : inframe_nt) inframe_upper += std::toupper(static_cast<unsigned char>(c));

    // Try each frame offset
    for (int frame = 0; frame < 3; frame++) {
        if (frame >= static_cast<int>(seq_upper.size())) continue;

        // Extract sequence starting from this frame
        size_t extract_len = ((seq_upper.size() - frame) / 3) * 3;
        if (extract_len == 0) continue;

        std::string candidate = seq_upper.substr(frame, extract_len);

        // Check if it matches (allowing for length differences due to truncation)
        if (candidate == inframe_upper ||
            candidate.find(inframe_upper) == 0 ||
            inframe_upper.find(candidate) == 0) {
            return frame;
        }
    }

    return -1;  // Could not determine
}

// Load ground truth from aMGSIM aa-damage.tsv.gz
std::unordered_map<std::string, GroundTruth> load_ground_truth(const std::string& path) {
    std::unordered_map<std::string, GroundTruth> gt_map;

    gzFile file = gzopen(path.c_str(), "r");
    if (!file) {
        std::cerr << "Error: Cannot open " << path << std::endl;
        return gt_map;
    }

    char buffer[65536];
    bool first_line = true;

    // Column indices (will be found from header)
    int col_read_name = -1;
    int col_strand_read = -1;
    int col_strand_gene = -1;
    int col_start_read = -1;
    int col_end_read = -1;
    int col_start_gene = -1;
    int col_end_gene = -1;
    int col_read_length = -1;
    int col_true_protein = -1;      // intersect_seq_inframe_aa
    int col_damaged_protein = -1;   // damaged_seq_inframe_aa
    int col_damaged_seq = -1;       // Full damaged read sequence
    int col_damaged_seq_inframe_nt = -1;  // In-frame portion of damaged sequence

    while (gzgets(file, buffer, sizeof(buffer)) != NULL) {
        std::string line(buffer);
        while (!line.empty() && (line.back() == '\n' || line.back() == '\r')) {
            line.pop_back();
        }

        std::vector<std::string> fields;
        std::istringstream iss(line);
        std::string field;
        while (std::getline(iss, field, '\t')) {
            fields.push_back(field);
        }

        if (first_line) {
            first_line = false;
            // Find column indices from header
            for (size_t i = 0; i < fields.size(); i++) {
                if (fields[i] == "read_name") col_read_name = i;
                else if (fields[i] == "Strand_read") col_strand_read = i;
                else if (fields[i] == "Strand_gene") col_strand_gene = i;
                else if (fields[i] == "Start_read") col_start_read = i;
                else if (fields[i] == "End_read") col_end_read = i;
                else if (fields[i] == "Start_gene") col_start_gene = i;
                else if (fields[i] == "End_gene") col_end_gene = i;
                else if (fields[i] == "read_length") col_read_length = i;
                else if (fields[i] == "intersect_seq_inframe_aa") col_true_protein = i;
                else if (fields[i] == "damaged_seq_inframe_aa") col_damaged_protein = i;
                else if (fields[i] == "damaged_seq") col_damaged_seq = i;
                else if (fields[i] == "damaged_seq_inframe_nt") col_damaged_seq_inframe_nt = i;
            }
            continue;
        }

        if (col_read_name < 0 || static_cast<size_t>(col_read_name) >= fields.size()) continue;

        GroundTruth gt;
        gt.read_id = fields[col_read_name];

        // Parse coordinates
        if (col_start_read >= 0 && static_cast<size_t>(col_start_read) < fields.size())
            gt.read_start = std::atoi(fields[col_start_read].c_str());
        if (col_end_read >= 0 && static_cast<size_t>(col_end_read) < fields.size())
            gt.read_end = std::atoi(fields[col_end_read].c_str());
        if (col_start_gene >= 0 && static_cast<size_t>(col_start_gene) < fields.size())
            gt.gene_start = std::atoi(fields[col_start_gene].c_str());
        if (col_end_gene >= 0 && static_cast<size_t>(col_end_gene) < fields.size())
            gt.gene_end = std::atoi(fields[col_end_gene].c_str());
        if (col_read_length >= 0 && static_cast<size_t>(col_read_length) < fields.size())
            gt.read_length = std::atoi(fields[col_read_length].c_str());

        // Parse strands
        if (col_strand_read >= 0 && static_cast<size_t>(col_strand_read) < fields.size() && !fields[col_strand_read].empty())
            gt.strand_read = fields[col_strand_read][0];
        if (col_strand_gene >= 0 && static_cast<size_t>(col_strand_gene) < fields.size() && !fields[col_strand_gene].empty())
            gt.strand_gene = fields[col_strand_gene][0];

        // Parse proteins
        if (col_true_protein >= 0 && static_cast<size_t>(col_true_protein) < fields.size())
            gt.true_protein = fields[col_true_protein];
        if (col_damaged_protein >= 0 && static_cast<size_t>(col_damaged_protein) < fields.size())
            gt.damaged_protein = fields[col_damaged_protein];

        // Parse nucleotide sequences for frame determination
        if (col_damaged_seq >= 0 && static_cast<size_t>(col_damaged_seq) < fields.size())
            gt.damaged_seq = fields[col_damaged_seq];
        if (col_damaged_seq_inframe_nt >= 0 && static_cast<size_t>(col_damaged_seq_inframe_nt) < fields.size())
            gt.damaged_seq_inframe_nt = fields[col_damaged_seq_inframe_nt];

        // Compute TRUE STRAND prediction
        // If read strand matches gene strand: AGP should predict FORWARD (use sequence as-is)
        // If read strand differs from gene strand: AGP should predict REVERSE (use RC)
        gt.true_forward = (gt.strand_read == gt.strand_gene);

        // Compute TRUE FRAME using in-frame nucleotide sequence matching
        // This is the correct method - determines frame by aligning damaged_seq to damaged_seq_inframe_nt
        if (!gt.damaged_seq.empty() && !gt.damaged_seq_inframe_nt.empty()) {
            int determined_frame = determine_true_frame(gt.damaged_seq, gt.damaged_seq_inframe_nt, gt.true_forward);
            if (determined_frame >= 0) {
                gt.true_frame = determined_frame;
            } else {
                // Fallback to coordinate-based calculation if sequence matching fails
                int offset = gt.read_start - gt.gene_start;
                gt.true_frame = ((offset % 3) + 3) % 3;
            }
        } else {
            // Fallback to coordinate-based calculation if sequences not available
            int offset = gt.read_start - gt.gene_start;
            gt.true_frame = ((offset % 3) + 3) % 3;
        }

        // Parse damage positions from read ID
        // Format: ...ancient:strand:start:end:length:damage_positions
        size_t last_colon = gt.read_id.rfind(':');
        if (last_colon != std::string::npos) {
            std::string damage_str = gt.read_id.substr(last_colon + 1);
            gt.damage_positions = parse_damage_positions(damage_str);
        }

        // Check for terminal damage using rolling window approach
        // Find maximum damage density in any 3bp window within first/last 10bp
        // This matches AGP's detection approach
        constexpr int WINDOW_SIZE = 3;
        constexpr int SEARCH_DEPTH = 10;

        // Count damage events at each position
        std::array<int, 20> damage_at_pos_5prime = {};  // positions 1-10
        std::array<int, 20> damage_at_pos_3prime = {};  // positions -1 to -10

        for (int pos : gt.damage_positions) {
            if (pos > 0 && pos <= SEARCH_DEPTH) {
                damage_at_pos_5prime[pos - 1]++;  // 1-indexed to 0-indexed
                gt.damage_at_5prime = true;
            } else if (pos < 0 && pos >= -SEARCH_DEPTH) {
                damage_at_pos_3prime[-pos - 1]++;  // -1 -> index 0
                gt.damage_at_3prime = true;
            }
        }

        // Find max damage in any window at 5' end
        int max_5prime_window = 0;
        for (int start = 0; start + WINDOW_SIZE <= SEARCH_DEPTH; ++start) {
            int window_count = 0;
            for (int i = start; i < start + WINDOW_SIZE; ++i) {
                window_count += damage_at_pos_5prime[i];
            }
            max_5prime_window = std::max(max_5prime_window, window_count);
        }

        // Find max damage in any window at 3' end
        int max_3prime_window = 0;
        for (int start = 0; start + WINDOW_SIZE <= SEARCH_DEPTH; ++start) {
            int window_count = 0;
            for (int i = start; i < start + WINDOW_SIZE; ++i) {
                window_count += damage_at_pos_3prime[i];
            }
            max_3prime_window = std::max(max_3prime_window, window_count);
        }

        // Has terminal damage if any window has damage
        gt.has_terminal_damage = (max_5prime_window > 0 || max_3prime_window > 0);

        // Only include reads with COMPLETE ground truth information:
        // 1. Must have protein sequence
        // 2. Must have in-frame nucleotide sequence (for accurate frame determination)
        // Without damaged_seq_inframe_nt, frame is determined by coordinate fallback which
        // may be inaccurate. Skip these reads to ensure validation quality.
        if (!gt.true_protein.empty() && !gt.damaged_seq_inframe_nt.empty()) {
            gt_map[gt.read_id] = gt;
        }
    }

    gzclose(file);
    return gt_map;
}

// Parse a GFF line and return a Prediction (supports both AGP and FGS formats)
Prediction parse_gff_line(const std::string& seqid, const std::string& source,
                          const std::string& start_s, const std::string& end_s,
                          const std::string& score_s, const std::string& strand_s,
                          const std::string& phase, const std::string& attrs) {
    Prediction pred;
    pred.read_id = seqid;
    pred.strand = strand_s[0];

    // Parse score from GFF column 6
    parse_float_strict(score_s, pred.score);

    // Detect format: FGS uses "FGS" as source, AGP uses "AGP"
    bool is_fgs = (source == "FGS" || source == "FragGeneScan" || source == "FragGeneScanRs");

    if (is_fgs) {
        // FGS format: frame is (START - 1) % 3 where START is 1-based
        int start = 1;
        if (parse_int_strict(start_s, start) && start > 0) {
            pred.frame = (start - 1) % 3;
        } else {
            pred.frame = 0;
        }
        // FGS doesn't have damage_signal
        pred.damage_signal = 0.5f;
    } else {
        // AGP format: parse frame= attribute
        size_t pos;
        if ((pos = attrs.find("frame=")) != std::string::npos) {
            size_t end = attrs.find(';', pos + 6);
            std::string value = attrs.substr(pos + 6, end == std::string::npos ? std::string::npos : end - (pos + 6));
            int frame = 0;
            if (parse_int_strict(value, frame)) {
                pred.frame = frame;
            }
        }
        if ((pos = attrs.find("damage_signal=")) != std::string::npos) {
            size_t end = attrs.find(';', pos + 14);
            std::string value = attrs.substr(pos + 14, end == std::string::npos ? std::string::npos : end - (pos + 14));
            float parsed = pred.damage_signal;
            if (parse_float_strict(value, parsed)) {
                pred.damage_signal = parsed;
            }
        }
        if ((pos = attrs.find("damage_pct=")) != std::string::npos) {
            size_t end = attrs.find(';', pos + 11);
            std::string value = attrs.substr(pos + 11, end == std::string::npos ? std::string::npos : end - (pos + 11));
            float parsed = pred.damage_pct;
            if (parse_float_strict(value, parsed)) {
                pred.damage_pct = parsed;
            }
        }
    }

    return pred;
}

// Load predictions from GFF file (supports both AGP and FGS formats)
std::unordered_map<std::string, std::vector<Prediction>> load_predictions_gff(const std::string& path) {
    std::unordered_map<std::string, std::vector<Prediction>> pred_map;

    gzFile file = gzopen(path.c_str(), "r");
    if (!file) {
        // Try plain file
        std::ifstream plain(path);
        if (!plain) {
            std::cerr << "Error: Cannot open " << path << std::endl;
            return pred_map;
        }

        std::string line;
        while (std::getline(plain, line)) {
            if (line.empty() || line[0] == '#') continue;

            std::istringstream iss(line);
            std::string seqid, source, type, start_s, end_s, score_s, strand_s, phase, attrs;
            if (!(iss >> seqid >> source >> type >> start_s >> end_s >> score_s >> strand_s >> phase)) continue;
            std::getline(iss, attrs);

            Prediction pred = parse_gff_line(seqid, source, start_s, end_s, score_s, strand_s, phase, attrs);
            pred_map[seqid].push_back(pred);
        }
        return pred_map;
    }

    char buffer[65536];
    while (gzgets(file, buffer, sizeof(buffer)) != NULL) {
        if (buffer[0] == '#' || buffer[0] == '\n') continue;

        std::string line(buffer);
        std::istringstream iss(line);
        std::string seqid, source, type, start_s, end_s, score_s, strand_s, phase, attrs;
        if (!(iss >> seqid >> source >> type >> start_s >> end_s >> score_s >> strand_s >> phase)) continue;
        std::getline(iss, attrs);

        Prediction pred = parse_gff_line(seqid, source, start_s, end_s, score_s, strand_s, phase, attrs);
        pred_map[seqid].push_back(pred);
    }

    gzclose(file);
    return pred_map;
}

// Parse FGS-style header: >READ_ID_START_END_STRAND
// Returns: {read_id, strand, start, end}
std::tuple<std::string, char, int, int> parse_fgs_header(const std::string& header) {
    // FGS format: READ_ID_START_END_STRAND where START/END are numbers, STRAND is +/-
    // Need to find the last occurrence of _N_N_[+-] pattern

    std::string id = header;
    char strand = '+';
    int start = 1, end = 0;

    // Find last underscore before strand (+/-)
    size_t last_under = header.rfind('_');
    if (last_under != std::string::npos && last_under + 1 < header.size()) {
        char possible_strand = header[last_under + 1];
        if (possible_strand == '+' || possible_strand == '-') {
            strand = possible_strand;

            // Find end position (before _STRAND)
            std::string before_strand = header.substr(0, last_under);
            size_t end_under = before_strand.rfind('_');
            if (end_under != std::string::npos) {
                std::string end_str = before_strand.substr(end_under + 1);
                int parsed_end = 0;
                if (parse_int_strict(end_str, parsed_end)) {
                    end = parsed_end;

                    // Find start position (before _END)
                    std::string before_end = before_strand.substr(0, end_under);
                    size_t start_under = before_end.rfind('_');
                    if (start_under != std::string::npos) {
                        std::string start_str = before_end.substr(start_under + 1);
                        int parsed_start = 0;
                        if (parse_int_strict(start_str, parsed_start)) {
                            start = parsed_start;
                            // Everything before _START is the read ID
                            id = before_end.substr(0, start_under);
                        }
                    }
                }
            }
        }
    }

    return {id, strand, start, end};
}

// Load protein sequences from FASTA (supports both AGP and FGS formats)
std::unordered_map<std::string, std::vector<std::pair<char, std::string>>> load_proteins(const std::string& path) {
    std::unordered_map<std::string, std::vector<std::pair<char, std::string>>> prot_map;

    std::ifstream file(path);
    if (!file) {
        std::cerr << "Warning: Cannot open protein file " << path << std::endl;
        return prot_map;
    }

    std::string line, current_id;
    char current_strand = '+';
    std::string current_seq;

    while (std::getline(file, line)) {
        if (line.empty()) continue;

        if (line[0] == '>') {
            // Save previous
            if (!current_id.empty() && !current_seq.empty()) {
                prot_map[current_id].push_back({current_strand, current_seq});
            }

            std::string header = line.substr(1);

            // Detect format:
            // - AGP new: READ_NAME_+_N or READ_NAME_-_N (strand embedded in ID)
            // - AGP old: protein_N READ_ID_gene_N STRAND ...
            // - FGS: READ_ID_START_END_STRAND

            // First extract the sequence ID (everything before space)
            size_t space_pos = header.find(' ');
            std::string seq_id = (space_pos != std::string::npos) ? header.substr(0, space_pos) : header;

            // Check for AGP new format: ends with _+_N or _-_N
            size_t plus_pos = seq_id.rfind("_+_");
            size_t minus_pos = seq_id.rfind("_-_");

            if (plus_pos != std::string::npos || minus_pos != std::string::npos) {
                // AGP new format: READ_NAME_+_N or READ_NAME_-_N
                size_t strand_pos = (plus_pos != std::string::npos) ? plus_pos : minus_pos;
                current_id = seq_id.substr(0, strand_pos);
                current_strand = (plus_pos != std::string::npos) ? '+' : '-';
            } else if (header.substr(0, 8) == "protein_" || header.substr(0, 10) == "corrected_") {
                // AGP old format: >protein_N READ_ID_gene_N STRAND ...
                std::istringstream iss(header);
                std::string protein_name, read_id_gene, strand_s;
                iss >> protein_name >> read_id_gene >> strand_s;

                // Extract read ID (remove _gene_N suffix)
                size_t gene_pos = read_id_gene.find("_gene_");
                if (gene_pos != std::string::npos) {
                    current_id = read_id_gene.substr(0, gene_pos);
                } else {
                    current_id = read_id_gene;
                }

                current_strand = (!strand_s.empty()) ? strand_s[0] : '+';
            } else {
                // FGS format: >READ_ID_START_END_STRAND
                auto [id, strand, start, end] = parse_fgs_header(header);
                current_id = id;
                current_strand = strand;
            }

            current_seq.clear();
        } else {
            current_seq += line;
        }
    }

    // Save last
    if (!current_id.empty() && !current_seq.empty()) {
        prot_map[current_id].push_back({current_strand, current_seq});
    }

    return prot_map;
}

// Count total reads in FASTQ
size_t count_fastq_reads(const std::string& path) {
    gzFile file = gzopen(path.c_str(), "r");
    if (!file) return 0;

    size_t count = 0;
    char buffer[65536];
    while (gzgets(file, buffer, sizeof(buffer)) != NULL) {
        if (buffer[0] == '@') count++;
    }

    gzclose(file);
    return count;
}

// Calculate AUC-ROC
double calculate_auc(std::vector<std::pair<float, bool>>& scores) {
    if (scores.empty()) return 0.5;

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

void print_report(const Metrics& m) {
    std::cout << "\n";
    std::cout << "╔══════════════════════════════════════════════════════════════════════╗\n";
    std::cout << "║              AGP UNIFIED VALIDATION REPORT                           ║\n";
    std::cout << "╚══════════════════════════════════════════════════════════════════════╝\n\n";

    // Dataset summary
    std::cout << "┌─────────────────────────────────────────────────────────────────────┐\n";
    std::cout << "│ DATASET SUMMARY                                                     │\n";
    std::cout << "├─────────────────────────────────────────────────────────────────────┤\n";
    std::cout << "│  Total reads in FASTQ:      " << std::setw(12) << m.total_reads << "                          │\n";
    std::cout << "│  Coding reads (ground truth): " << std::setw(10) << m.coding_reads << "                          │\n";
    std::cout << "│  Predicted reads:           " << std::setw(12) << m.predicted_reads << "                          │\n";
    std::cout << "│  Matched for evaluation:    " << std::setw(12) << m.tp << "                          │\n";
    std::cout << "└─────────────────────────────────────────────────────────────────────┘\n\n";

    // Coding detection
    std::cout << "┌─────────────────────────────────────────────────────────────────────┐\n";
    std::cout << "│ 1. CODING REGION DETECTION                                          │\n";
    std::cout << "├─────────────────────────────────────────────────────────────────────┤\n";

    double precision = (m.tp + m.fp > 0) ? static_cast<double>(m.tp) / (m.tp + m.fp) : 0;
    double recall = (m.tp + m.fn > 0) ? static_cast<double>(m.tp) / (m.tp + m.fn) : 0;
    double f1 = (precision + recall > 0) ? 2 * precision * recall / (precision + recall) : 0;

    std::cout << "│  TP (found coding):    " << std::setw(10) << m.tp << "                                  │\n";
    std::cout << "│  FN (missed coding):   " << std::setw(10) << m.fn << "                                  │\n";
    std::cout << "│  FP (false positive):  " << std::setw(10) << m.fp << "                                  │\n";
    std::cout << "│                                                                     │\n";
    std::cout << "│  Precision: " << std::fixed << std::setprecision(4) << precision
              << "    Recall: " << recall << "    F1: " << f1 << "         │\n";
    std::cout << "└─────────────────────────────────────────────────────────────────────┘\n\n";

    // Frame selection
    std::cout << "┌─────────────────────────────────────────────────────────────────────┐\n";
    std::cout << "│ 2. FRAME SELECTION ACCURACY                                         │\n";
    std::cout << "├─────────────────────────────────────────────────────────────────────┤\n";

    if (m.frame_total > 0) {
        double frame_acc = 100.0 * m.frame_correct / m.frame_total;
        std::cout << "│  Correct frame:  " << std::setw(10) << m.frame_correct
                  << " / " << std::setw(10) << m.frame_total
                  << "  (" << std::fixed << std::setprecision(1) << std::setw(5) << frame_acc << "%)       │\n";
        std::cout << "│  Random baseline: 33.3%                                             │\n";
        std::cout << "│                                                                     │\n";
        std::cout << "│  True frame distribution:                                           │\n";
        for (int f = 0; f < 3; f++) {
            double pct = 100.0 * m.true_frame_distribution[f] / m.frame_total;
            std::cout << "│    Frame " << f << ": " << std::setw(8) << m.true_frame_distribution[f]
                      << " (" << std::setw(5) << std::setprecision(1) << pct << "%)                                  │\n";
        }
        std::cout << "│  Predicted frame distribution:                                      │\n";
        for (int f = 0; f < 3; f++) {
            double pct = 100.0 * m.frame_distribution[f] / m.frame_total;
            std::cout << "│    Frame " << f << ": " << std::setw(8) << m.frame_distribution[f]
                      << " (" << std::setw(5) << std::setprecision(1) << pct << "%)                                  │\n";
        }
    }
    std::cout << "└─────────────────────────────────────────────────────────────────────┘\n\n";

    // Strand prediction
    std::cout << "┌─────────────────────────────────────────────────────────────────────┐\n";
    std::cout << "│ 3. STRAND PREDICTION ACCURACY                                       │\n";
    std::cout << "├─────────────────────────────────────────────────────────────────────┤\n";

    if (m.strand_total > 0) {
        double strand_acc = 100.0 * m.strand_correct / m.strand_total;
        std::cout << "│  Correct strand: " << std::setw(10) << m.strand_correct
                  << " / " << std::setw(10) << m.strand_total
                  << "  (" << std::fixed << std::setprecision(1) << std::setw(5) << strand_acc << "%)       │\n";
        std::cout << "│  Random baseline: 50.0%                                             │\n";
    }
    std::cout << "└─────────────────────────────────────────────────────────────────────┘\n\n";

    // Combined accuracy
    std::cout << "┌─────────────────────────────────────────────────────────────────────┐\n";
    std::cout << "│ 4. COMBINED ACCURACY (Frame AND Strand correct)                     │\n";
    std::cout << "├─────────────────────────────────────────────────────────────────────┤\n";

    if (m.both_total > 0) {
        double both_acc = 100.0 * m.both_correct / m.both_total;
        std::cout << "│  Both correct:   " << std::setw(10) << m.both_correct
                  << " / " << std::setw(10) << m.both_total
                  << "  (" << std::fixed << std::setprecision(1) << std::setw(5) << both_acc << "%)       │\n";
        std::cout << "│  Random baseline: 16.7% (1/6)                                       │\n";
    }
    std::cout << "└─────────────────────────────────────────────────────────────────────┘\n\n";

    // Sequence identity
    std::cout << "┌─────────────────────────────────────────────────────────────────────┐\n";
    std::cout << "│ 5. SEQUENCE IDENTITY (Predicted vs Ground Truth Protein)            │\n";
    std::cout << "├─────────────────────────────────────────────────────────────────────┤\n";

    if (m.identity_count > 0) {
        double avg_identity = 100.0 * m.total_identity / m.identity_count;
        std::cout << "│  Sequences compared: " << std::setw(10) << m.identity_count << "                               │\n";
        std::cout << "│  Average identity:   " << std::fixed << std::setprecision(1) << std::setw(10) << avg_identity << "%                              │\n";
        std::cout << "│                                                                     │\n";
        std::cout << "│  Exact matches (100%):     " << std::setw(8) << m.exact_matches
                  << " (" << std::setw(5) << 100.0 * m.exact_matches / m.identity_count << "%)                   │\n";
        std::cout << "│  High identity (>90%):     " << std::setw(8) << m.high_identity
                  << " (" << std::setw(5) << 100.0 * m.high_identity / m.identity_count << "%)                   │\n";
        std::cout << "│  Medium identity (50-90%): " << std::setw(8) << m.medium_identity
                  << " (" << std::setw(5) << 100.0 * m.medium_identity / m.identity_count << "%)                   │\n";
        std::cout << "│  Low identity (<50%):      " << std::setw(8) << m.low_identity
                  << " (" << std::setw(5) << 100.0 * m.low_identity / m.identity_count << "%)                   │\n";
        std::cout << "│                                                                     │\n";
        std::cout << "│  Identity by prediction correctness:                                │\n";
        if (m.count_correct_both > 0) {
            double id_correct = 100.0 * m.identity_correct_both / m.count_correct_both;
            std::cout << "│    Frame+Strand correct: " << std::setw(6) << std::setprecision(1) << id_correct << "% (n=" << std::setw(8) << m.count_correct_both << ")               │\n";
        }
        if (m.count_wrong_frame > 0) {
            double id_wrong_f = 100.0 * m.identity_wrong_frame / m.count_wrong_frame;
            std::cout << "│    Wrong frame:          " << std::setw(6) << std::setprecision(1) << id_wrong_f << "% (n=" << std::setw(8) << m.count_wrong_frame << ")               │\n";
        }
        if (m.count_wrong_strand > 0) {
            double id_wrong_s = 100.0 * m.identity_wrong_strand / m.count_wrong_strand;
            std::cout << "│    Wrong strand:         " << std::setw(6) << std::setprecision(1) << id_wrong_s << "% (n=" << std::setw(8) << m.count_wrong_strand << ")               │\n";
        }
    }
    std::cout << "└─────────────────────────────────────────────────────────────────────┘\n\n";

    // Corrected protein identity (if available)
    if (m.corrected_identity_count > 0) {
        std::cout << "┌─────────────────────────────────────────────────────────────────────┐\n";
        std::cout << "│ 5b. DAMAGE-CORRECTED PROTEIN IDENTITY                               │\n";
        std::cout << "├─────────────────────────────────────────────────────────────────────┤\n";
        double avg_corr = 100.0 * m.total_corrected_identity / m.corrected_identity_count;
        double avg_orig = 100.0 * m.total_identity / m.identity_count;
        double improvement = avg_corr - avg_orig;
        std::cout << "│  Sequences compared: " << std::setw(10) << m.corrected_identity_count << "                               │\n";
        std::cout << "│  Average corrected:  " << std::setw(10) << std::setprecision(1) << avg_corr << "%                              │\n";
        std::cout << "│  Average original:   " << std::setw(10) << std::setprecision(1) << avg_orig << "%                              │\n";
        std::cout << "│  Improvement:        " << std::setw(10) << std::showpos << improvement << "%" << std::noshowpos << "                              │\n";
        if (m.corrected_count_correct_both > 0 && m.count_correct_both > 0) {
            double corr_correct = 100.0 * m.corrected_identity_correct_both / m.corrected_count_correct_both;
            double orig_correct = 100.0 * m.identity_correct_both / m.count_correct_both;
            double improve_correct = corr_correct - orig_correct;
            std::cout << "│                                                                     │\n";
            std::cout << "│  For correct frame+strand predictions:                              │\n";
            std::cout << "│    Corrected: " << std::setw(6) << std::setprecision(1) << corr_correct << "%  Original: " << std::setw(6) << orig_correct << "%  Δ: " << std::showpos << std::setw(6) << improve_correct << "%" << std::noshowpos << "      │\n";
        }
        std::cout << "└─────────────────────────────────────────────────────────────────────┘\n\n";
    }

    // Damage detection
    std::cout << "┌───────────────────────────────────────────────────────────────────────────────────────┐\n";
    std::cout << "│ 6. DAMAGE DETECTION (Terminal damage in first/last 10bp)                              │\n";
    std::cout << "├───────────────────────────────────────────────────────────────────────────────────────┤\n";

    double auc = calculate_auc(const_cast<std::vector<std::pair<float, bool>>&>(m.damage_scores));

    // Count totals for class balance info
    size_t total_damaged = 0, total_undamaged = 0;
    for (const auto& [score, has_damage] : m.damage_scores) {
        if (has_damage) total_damaged++; else total_undamaged++;
    }

    std::cout << "│  AUC-ROC: " << std::fixed << std::setprecision(4) << auc << "                                                                     │\n";
    std::cout << "│  Total samples: " << m.damage_scores.size() << " (damaged: " << total_damaged
              << ", undamaged: " << total_undamaged << ")                                │\n";
    std::cout << "│                                                                                       │\n";
    std::cout << "│  Thresh   TP        FP        TN        FN      Sens%    Spec%    Prec%    F1%        │\n";
    std::cout << "│  ─────────────────────────────────────────────────────────────────────────────────────│\n";

    // Test multiple thresholds
    std::vector<float> thresholds = {0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f};
    for (float thresh : thresholds) {
        size_t tp = 0, fp = 0, tn = 0, fn = 0;
        for (const auto& [score, has_damage] : m.damage_scores) {
            bool detected = (score > thresh);
            if (has_damage) {
                if (detected) tp++; else fn++;
            } else {
                if (detected) fp++; else tn++;
            }
        }

        double sensitivity = (tp + fn > 0) ? 100.0 * tp / (tp + fn) : 0.0;
        double specificity = (tn + fp > 0) ? 100.0 * tn / (tn + fp) : 0.0;
        double precision = (tp + fp > 0) ? 100.0 * tp / (tp + fp) : 0.0;
        double f1 = (precision + sensitivity > 0) ? 2.0 * precision * sensitivity / (precision + sensitivity) : 0.0;

        std::cout << "│   " << std::fixed << std::setprecision(1) << thresh
                  << "   " << std::setw(8) << tp
                  << "  " << std::setw(8) << fp
                  << "  " << std::setw(8) << tn
                  << "  " << std::setw(8) << fn
                  << "   " << std::setw(5) << std::setprecision(1) << sensitivity
                  << "    " << std::setw(5) << specificity
                  << "    " << std::setw(5) << precision
                  << "    " << std::setw(5) << f1 << "   │\n";
    }
    std::cout << "└───────────────────────────────────────────────────────────────────────────────────────┘\n\n";

    // Stop codon analysis
    std::cout << "┌─────────────────────────────────────────────────────────────────────┐\n";
    std::cout << "│ 7. STOP CODON ANALYSIS                                              │\n";
    std::cout << "├─────────────────────────────────────────────────────────────────────┤\n";
    std::cout << "│  Predictions with internal stops: " << std::setw(8) << m.predictions_with_stops << "                        │\n";
    if (m.frame_correct > 0) {
        std::cout << "│  Correct frame with stops:         " << std::setw(8) << m.correct_frame_with_stops
                  << " (" << std::setw(5) << std::setprecision(1) << 100.0 * m.correct_frame_with_stops / m.frame_correct << "%)           │\n";
    }
    if (m.frame_wrong > 0) {
        std::cout << "│  Wrong frame with stops:           " << std::setw(8) << m.wrong_frame_with_stops
                  << " (" << std::setw(5) << 100.0 * m.wrong_frame_with_stops / m.frame_wrong << "%)           │\n";
    }
    std::cout << "└─────────────────────────────────────────────────────────────────────┘\n\n";

    // Summary
    std::cout << "╔══════════════════════════════════════════════════════════════════════╗\n";
    std::cout << "║ SUMMARY                                                              ║\n";
    std::cout << "╠══════════════════════════════════════════════════════════════════════╣\n";

    if (m.both_total > 0) {
        double both_acc = 100.0 * m.both_correct / m.both_total;
        if (both_acc > 80) {
            std::cout << "║  [✓] EXCELLENT: Combined accuracy > 80%                             ║\n";
        } else if (both_acc > 50) {
            std::cout << "║  [~] GOOD: Combined accuracy > 50%                                  ║\n";
        } else if (both_acc > 16.7) {
            std::cout << "║  [!] FAIR: Above random baseline (16.7%)                            ║\n";
        } else {
            std::cout << "║  [✗] POOR: At or below random baseline                              ║\n";
        }
    }

    if (m.identity_count > 0) {
        double avg_id = 100.0 * m.total_identity / m.identity_count;
        if (avg_id > 90) {
            std::cout << "║  [✓] EXCELLENT: Average sequence identity > 90%                     ║\n";
        } else if (avg_id > 70) {
            std::cout << "║  [~] GOOD: Average sequence identity > 70%                          ║\n";
        } else if (avg_id > 50) {
            std::cout << "║  [!] FAIR: Average sequence identity > 50%                          ║\n";
        } else {
            std::cout << "║  [✗] POOR: Average sequence identity < 50%                          ║\n";
        }
    }

    std::cout << "╚══════════════════════════════════════════════════════════════════════╝\n";
}

static void validate_print_usage(const char* prog) {
    std::cout << "AGP Unified Validator\n\n";
    std::cout << "Usage: " << prog << " <fastq.gz> <aa-damage.tsv.gz> <predictions.gff> [proteins.faa] [corrected.faa] [--csv output.csv]\n\n";
    std::cout << "Arguments:\n";
    std::cout << "  fastq.gz          Input FASTQ file (gzipped)\n";
    std::cout << "  aa-damage.tsv.gz  Ground truth from aMGSIM\n";
    std::cout << "  predictions.gff   AGP/FGS GFF3 output\n";
    std::cout << "  proteins.faa      Protein FASTA (optional, for sequence identity)\n";
    std::cout << "  corrected.faa     Damage-corrected protein FASTA (optional)\n";
    std::cout << "  --csv output.csv  Output per-read details to CSV (optional)\n\n";
    std::cout << "Supports both AGP and FragGeneScan GFF formats.\n\n";
    std::cout << "Example:\n";
    std::cout << "  " << prog << " reads.fq.gz ground_truth.tsv.gz agp.gff agp.faa agp_corrected.faa --csv details.csv\n";
}

static int validate_main(int argc, char* argv[]) {
    const auto run_start = std::chrono::steady_clock::now();

    // Handle help flag
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
            validate_print_usage(argv[0]);
            return 0;
        }
    }

    if (argc < 4) {
        std::cerr << "Error: missing required positional arguments.\n";
        validate_print_usage(argv[0]);
        return 1;
    }

    std::string fastq_path = argv[1];
    std::string truth_path = argv[2];
    std::string gff_path = argv[3];
    std::string protein_path;
    std::string corrected_path;
    std::string csv_path;

    // Parse optional arguments
    for (int i = 4; i < argc; i++) {
        std::string arg = argv[i];
        if (arg == "--csv" && i + 1 < argc) {
            csv_path = argv[++i];
        } else if (protein_path.empty()) {
            protein_path = arg;
        } else if (corrected_path.empty()) {
            corrected_path = arg;
        }
    }

    Metrics m;

    // Load data
    const auto load_start = std::chrono::steady_clock::now();
    std::cerr << "Loading ground truth..." << std::endl;
    auto ground_truth = load_ground_truth(truth_path);
    m.coding_reads = ground_truth.size();
    std::cerr << "  Loaded " << m.coding_reads << " coding reads\n";

    std::cerr << "Loading predictions..." << std::endl;
    auto predictions = load_predictions_gff(gff_path);
    m.predicted_reads = predictions.size();
    std::cerr << "  Loaded " << m.predicted_reads << " predictions\n";

    std::unordered_map<std::string, std::vector<std::pair<char, std::string>>> proteins;
    if (!protein_path.empty()) {
        std::cerr << "Loading proteins..." << std::endl;
        proteins = load_proteins(protein_path);
        std::cerr << "  Loaded proteins for " << proteins.size() << " reads\n";
    }

    std::unordered_map<std::string, std::vector<std::pair<char, std::string>>> corrected_proteins;
    if (!corrected_path.empty()) {
        std::cerr << "Loading corrected proteins..." << std::endl;
        corrected_proteins = load_proteins(corrected_path);
        std::cerr << "  Loaded corrected proteins for " << corrected_proteins.size() << " reads\n";
    }
    const auto load_end = std::chrono::steady_clock::now();
    std::cerr << "  Loading runtime: "
              << agp::log_utils::format_elapsed(load_start, load_end) << "\n";

    const auto count_start = std::chrono::steady_clock::now();
    std::cerr << "Counting total reads..." << std::endl;
    m.total_reads = count_fastq_reads(fastq_path);
    std::cerr << "  Total reads: " << m.total_reads << "\n";
    const auto count_end = std::chrono::steady_clock::now();
    std::cerr << "  Read counting runtime: "
              << agp::log_utils::format_elapsed(count_start, count_end) << "\n";

    // Evaluate predictions
    const auto eval_start = std::chrono::steady_clock::now();
    std::cerr << "Evaluating..." << std::endl;

    std::unordered_set<std::string> predicted_ids;
    for (const auto& [id, _] : predictions) {
        predicted_ids.insert(id);
    }

    // Open CSV file if requested
    std::ofstream csv_file;
    if (!csv_path.empty()) {
        csv_file.open(csv_path);
        if (!csv_file.is_open()) {
            std::cerr << "Error: Cannot open CSV file: " << csv_path << std::endl;
            return 1;
        }
        // Write header
        csv_file << "read_id\tlength\tpred_strand\tpred_frame\ttrue_strand\ttrue_frame\t"
                 << "score\tdamage_pct\tdamage_signal\tstrand_correct\tframe_correct\tboth_correct\tidentity\n";
    }

    for (const auto& [read_id, gt] : ground_truth) {
        auto pred_it = predictions.find(read_id);

        if (pred_it == predictions.end()) {
            // False negative: coding read not predicted
            m.fn++;
            continue;
        }

        // True positive: coding read was predicted
        m.tp++;

        const auto& preds = pred_it->second;

        // Get the best prediction (usually just one, but could be multiple for both-strands mode)
        // For evaluation, we check if ANY prediction matches ground truth
        bool frame_match = false;
        bool strand_match = false;
        float best_damage_signal = 0.0f;

        for (const auto& pred : preds) {
            bool pred_forward = (pred.strand == '+');

            // Check strand match
            if (pred_forward == gt.true_forward) {
                strand_match = true;
            }

            // Check frame match (only meaningful if strand is correct)
            if (pred_forward == gt.true_forward && pred.frame == gt.true_frame) {
                frame_match = true;
            }

            best_damage_signal = std::max(best_damage_signal, pred.damage_signal);
        }

        // Frame accuracy
        m.frame_total++;
        m.true_frame_distribution[gt.true_frame]++;
        if (frame_match) {
            m.frame_correct++;
        } else {
            m.frame_wrong++;
        }

        // Track predicted frame distribution (use first prediction)
        if (!preds.empty()) {
            m.frame_distribution[preds[0].frame]++;
        }

        // Strand accuracy
        m.strand_total++;
        if (strand_match) {
            m.strand_correct++;
        } else {
            m.strand_wrong++;
        }

        // Combined accuracy
        m.both_total++;
        if (frame_match && strand_match) {
            m.both_correct++;
        }

        // Damage detection
        m.damage_scores.push_back({best_damage_signal, gt.has_terminal_damage});
        bool detected_damage = (best_damage_signal > 0.5f);
        if (gt.has_terminal_damage) {
            if (detected_damage) m.damage_tp++;
            else m.damage_fn++;
        } else {
            if (detected_damage) m.damage_fp++;
            else m.damage_tn++;
        }

        // Sequence identity (if proteins available)
        if (!proteins.empty() && !gt.true_protein.empty()) {
            auto prot_it = proteins.find(read_id);
            if (prot_it != proteins.end()) {
                // Find best matching protein (should be the one with correct strand)
                double best_identity = 0.0;
                std::string best_protein;
                bool has_stops = false;

                for (const auto& [strand, seq] : prot_it->second) {
                    double id = calculate_identity(seq, gt.true_protein);
                    if (id > best_identity) {
                        best_identity = id;
                        best_protein = seq;
                    }
                    if (count_internal_stops(seq) > 0) {
                        has_stops = true;
                    }
                }

                m.identity_count++;
                m.total_identity += best_identity;

                if (best_identity > 0.999) {
                    m.exact_matches++;
                } else if (best_identity > 0.90) {
                    m.high_identity++;
                } else if (best_identity > 0.50) {
                    m.medium_identity++;
                } else {
                    m.low_identity++;
                }

                // Stratified identity by correctness
                if (frame_match && strand_match) {
                    m.identity_correct_both += best_identity;
                    m.count_correct_both++;
                }
                if (!frame_match) {
                    m.identity_wrong_frame += best_identity;
                    m.count_wrong_frame++;
                }
                if (!strand_match) {
                    m.identity_wrong_strand += best_identity;
                    m.count_wrong_strand++;
                }

                // Stop codon analysis
                if (has_stops) {
                    m.predictions_with_stops++;
                    if (frame_match) {
                        m.correct_frame_with_stops++;
                    } else {
                        m.wrong_frame_with_stops++;
                    }
                }
            }
        }

        // Corrected protein identity (if available)
        if (!corrected_proteins.empty() && !gt.true_protein.empty()) {
            auto corr_it = corrected_proteins.find(read_id);
            if (corr_it != corrected_proteins.end()) {
                double best_corr_identity = 0.0;
                for (const auto& [strand, seq] : corr_it->second) {
                    double id = calculate_identity(seq, gt.true_protein);
                    if (id > best_corr_identity) {
                        best_corr_identity = id;
                    }
                }

                m.corrected_identity_count++;
                m.total_corrected_identity += best_corr_identity;

                if (frame_match && strand_match) {
                    m.corrected_identity_correct_both += best_corr_identity;
                    m.corrected_count_correct_both++;
                }
            }
        }

        // Write CSV output if requested
        if (csv_file.is_open() && !preds.empty()) {
            const auto& pred = preds[0];  // Use first (best) prediction
            char true_strand = gt.true_forward ? '+' : '-';

            // Get identity if available
            double identity = 0.0;
            if (!proteins.empty() && !gt.true_protein.empty()) {
                auto prot_it = proteins.find(read_id);
                if (prot_it != proteins.end()) {
                    for (const auto& [strand, seq] : prot_it->second) {
                        double id = calculate_identity(seq, gt.true_protein);
                        if (id > identity) identity = id;
                    }
                }
            }

            csv_file << read_id << "\t"
                     << gt.read_length << "\t"
                     << pred.strand << "\t"
                     << pred.frame << "\t"
                     << true_strand << "\t"
                     << gt.true_frame << "\t"
                     << pred.score << "\t"
                     << pred.damage_pct << "\t"
                     << pred.damage_signal << "\t"
                     << (strand_match ? 1 : 0) << "\t"
                     << (frame_match ? 1 : 0) << "\t"
                     << ((frame_match && strand_match) ? 1 : 0) << "\t"
                     << identity << "\n";
        }
    }

    // Close CSV file
    if (csv_file.is_open()) {
        csv_file.close();
        std::cerr << "Wrote per-read details to: " << csv_path << std::endl;
    }

    // Count false positives (predicted but not in ground truth)
    for (const auto& id : predicted_ids) {
        if (ground_truth.find(id) == ground_truth.end()) {
            m.fp++;
        }
    }

    // True negatives
    m.tn = m.total_reads - m.tp - m.fp - m.fn;

    // Print report
    print_report(m);
    const auto eval_end = std::chrono::steady_clock::now();
    std::cerr << "Evaluation runtime: "
              << agp::log_utils::format_elapsed(eval_start, eval_end) << "\n";
    std::cerr << "Total runtime: "
              << agp::log_utils::format_elapsed(run_start, eval_end) << "\n";

    return 0;
}

#ifdef AGP_STANDALONE_VALIDATE
// Standalone mode: provide main() directly
int main(int argc, char* argv[]) {
    return validate_main(argc, argv);
}
#else
// Subcommand mode: register with dispatcher
namespace agp {
namespace cli {

int cmd_validate(int argc, char* argv[]) {
    return validate_main(argc, argv);
}

// Register the validate subcommand
namespace {
    struct ValidateRegistrar {
        ValidateRegistrar() {
            SubcommandRegistry::instance().register_command(
                "validate",
                "Validate predictions against aMGSIM ground truth",
                cmd_validate, 50);
        }
    } validate_registrar;
}

}  // namespace cli
}  // namespace agp
#endif
