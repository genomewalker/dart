/**
 * benchmark_hamming.cpp
 *
 * Fast C++ tool to benchmark AGP corrections using Hamming distance.
 * Compares corrected sequences against ground truth original sequences.
 *
 * CRITICAL: Separates frame/strand accuracy from correction quality.
 * Only computes correction Hamming distance for reads where AGP
 * selected the correct frame (matching protein length ±5%).
 *
 * Usage:
 *   ./benchmark_hamming <ground_truth.tsv.gz> <agp_output.gff> [--fasta proteins.faa]
 *
 * Metrics computed:
 *   1. Frame accuracy - did AGP select the correct reading frame?
 *   2. Strand accuracy - did AGP select the correct strand?
 *   3. Hamming(damaged, original) - baseline damage level (correct frame only)
 *   4. Hamming(corrected, original) - AGP correction quality (correct frame only)
 *   5. Improvement = baseline - corrected (positive = AGP helped)
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>
#include <algorithm>
#include <cmath>
#include <zlib.h>

// Simple gzip reader
class GzipReader {
public:
    gzFile file;
    char buffer[65536];
    std::string line_buffer;

    GzipReader(const std::string& path) {
        file = gzopen(path.c_str(), "rb");
        if (!file) {
            throw std::runtime_error("Cannot open: " + path);
        }
    }

    ~GzipReader() {
        if (file) gzclose(file);
    }

    bool getline(std::string& line) {
        line.clear();
        while (true) {
            if (gzgets(file, buffer, sizeof(buffer)) == nullptr) {
                if (line_buffer.empty()) return false;
                line = line_buffer;
                line_buffer.clear();
                return true;
            }
            line_buffer += buffer;
            size_t pos = line_buffer.find('\n');
            if (pos != std::string::npos) {
                line = line_buffer.substr(0, pos);
                line_buffer = line_buffer.substr(pos + 1);
                return true;
            }
        }
    }
};

struct GroundTruth {
    std::string read_name;
    std::string original_aa;   // intersect_seq_inframe_aa - undamaged
    std::string damaged_aa;    // damaged_seq_inframe_aa - with damage
    std::string damage_diffs;  // damage_aaseq_diffs
    char strand = '?';         // Strand_read: + or -
    int damage_count = 0;
};

struct AGPResult {
    std::string read_name;
    std::string corrected_aa;
    char strand = '?';         // From GFF column 7
    int start = 0;             // From GFF column 4
    int end = 0;               // From GFF column 5
    int aa_corrections = 0;
    float damage_signal = 0;
};

struct FASTAEntry {
    std::string sequence;
    char strand = '?';
    int start = 0;
    int end = 0;
};

// Compute Hamming distance (number of mismatches)
int hamming_distance(const std::string& a, const std::string& b) {
    int len = std::min(a.size(), b.size());
    int dist = std::abs((int)a.size() - (int)b.size()); // Length difference
    for (int i = 0; i < len; i++) {
        if (a[i] != b[i]) dist++;
    }
    return dist;
}

// Check if two proteins are likely from the same frame
// (similar length, not completely different)
bool same_frame(const std::string& a, const std::string& b, float tolerance = 0.15) {
    if (a.empty() || b.empty()) return false;

    int len_a = a.size();
    int len_b = b.size();
    int min_len = std::min(len_a, len_b);
    int max_len = std::max(len_a, len_b);

    // Length must be within tolerance
    float len_ratio = (float)min_len / max_len;
    if (len_ratio < (1.0 - tolerance)) return false;

    // Also check that they're not completely different
    // (Hamming distance should be < 50% of length if same frame)
    int dist = hamming_distance(a, b);
    float dist_ratio = (float)dist / max_len;

    return dist_ratio < 0.5; // Less than 50% different
}

// Parse TSV header to find column indices
std::unordered_map<std::string, int> parse_header(const std::string& header) {
    std::unordered_map<std::string, int> cols;
    std::istringstream ss(header);
    std::string col;
    int idx = 0;
    while (std::getline(ss, col, '\t')) {
        cols[col] = idx++;
    }
    return cols;
}

// Split string by delimiter
std::vector<std::string> split(const std::string& s, char delim) {
    std::vector<std::string> parts;
    std::istringstream ss(s);
    std::string part;
    while (std::getline(ss, part, delim)) {
        parts.push_back(part);
    }
    return parts;
}

// Load ground truth from gzipped TSV
std::unordered_map<std::string, GroundTruth> load_ground_truth(const std::string& path) {
    std::unordered_map<std::string, GroundTruth> truth;
    GzipReader reader(path);

    std::string line;
    reader.getline(line); // Header
    auto cols = parse_header(line);

    int col_name = cols.count("read_name") ? cols["read_name"] : -1;
    int col_orig = cols.count("intersect_seq_inframe_aa") ? cols["intersect_seq_inframe_aa"] : -1;
    int col_dmg = cols.count("damaged_seq_inframe_aa") ? cols["damaged_seq_inframe_aa"] : -1;
    int col_diffs = cols.count("damage_aaseq_diffs") ? cols["damage_aaseq_diffs"] : -1;
    int col_strand = cols.count("Strand_read") ? cols["Strand_read"] : -1;

    if (col_name < 0 || col_orig < 0) {
        throw std::runtime_error("Required columns not found in ground truth");
    }

    while (reader.getline(line)) {
        auto fields = split(line, '\t');
        int max_col = std::max({col_name, col_orig, col_dmg, col_diffs, col_strand});
        if (fields.size() <= (size_t)max_col) continue;

        GroundTruth gt;
        gt.read_name = fields[col_name];
        gt.original_aa = fields[col_orig];
        gt.damaged_aa = col_dmg >= 0 ? fields[col_dmg] : "";
        gt.damage_diffs = col_diffs >= 0 ? fields[col_diffs] : "";

        if (col_strand >= 0 && !fields[col_strand].empty()) {
            gt.strand = fields[col_strand][0];
        }

        // Count damage positions
        if (!gt.damage_diffs.empty() && gt.damage_diffs != "None") {
            gt.damage_count = std::count(gt.damage_diffs.begin(), gt.damage_diffs.end(), ':');
        }

        truth[gt.read_name] = gt;
    }

    return truth;
}

// Load AGP results from GFF
std::unordered_map<std::string, AGPResult> load_agp_gff(const std::string& path) {
    std::unordered_map<std::string, AGPResult> results;
    std::ifstream file(path);

    std::string line;
    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#') continue;

        auto fields = split(line, '\t');
        if (fields.size() < 9) continue;

        AGPResult r;
        r.read_name = fields[0];

        // GFF columns: seqid source type start end score strand phase attributes
        if (fields.size() > 3) r.start = std::stoi(fields[3]);
        if (fields.size() > 4) r.end = std::stoi(fields[4]);
        if (fields.size() > 6 && !fields[6].empty()) r.strand = fields[6][0];

        // Parse attributes
        for (auto& attr : split(fields[8], ';')) {
            size_t eq = attr.find('=');
            if (eq != std::string::npos) {
                std::string key = attr.substr(0, eq);
                std::string val = attr.substr(eq + 1);
                if (key == "aa_corr") r.aa_corrections = std::stoi(val);
                else if (key == "damage_signal") r.damage_signal = std::stof(val);
            }
        }

        results[r.read_name] = r;
    }

    return results;
}

// Load corrected proteins from FASTA
// Header format: >protein_N sample___readname_gene_N strand start..end ...
std::unordered_map<std::string, FASTAEntry> load_fasta(const std::string& path) {
    std::unordered_map<std::string, FASTAEntry> seqs;
    std::ifstream file(path);

    std::string line, name, seq;
    FASTAEntry entry;

    while (std::getline(file, line)) {
        if (line.empty()) continue;
        if (line[0] == '>') {
            if (!name.empty()) {
                entry.sequence = seq;
                seqs[name] = entry;
            }

            // Parse header: >protein_N readname_gene_N strand start..end ...
            std::string header = line.substr(1);
            auto parts = split(header, ' ');

            // Extract readname between first space and "_gene_"
            name.clear();
            entry = FASTAEntry();

            if (parts.size() >= 2) {
                name = parts[1];
                // Find _gene_ suffix
                size_t gene_pos = name.find("_gene_");
                if (gene_pos != std::string::npos) {
                    name = name.substr(0, gene_pos);
                }
            }

            // Extract strand (third field)
            if (parts.size() >= 3 && !parts[2].empty()) {
                entry.strand = parts[2][0];
            }

            // Extract start..end (fourth field)
            if (parts.size() >= 4) {
                size_t dots = parts[3].find("..");
                if (dots != std::string::npos) {
                    try {
                        entry.start = std::stoi(parts[3].substr(0, dots));
                        entry.end = std::stoi(parts[3].substr(dots + 2));
                    } catch (...) {}
                }
            }

            seq.clear();
        } else {
            seq += line;
        }
    }
    if (!name.empty()) {
        entry.sequence = seq;
        seqs[name] = entry;
    }

    return seqs;
}

int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <ground_truth.tsv.gz> <agp_output.gff> [--fasta proteins.faa]\n";
        return 1;
    }

    std::string truth_path = argv[1];
    std::string gff_path = argv[2];
    std::string fasta_path;

    for (int i = 3; i < argc; i++) {
        if (std::string(argv[i]) == "--fasta" && i + 1 < argc) {
            fasta_path = argv[++i];
        }
    }

    std::cerr << "Loading ground truth from " << truth_path << "...\n";
    auto truth = load_ground_truth(truth_path);
    std::cerr << "  Loaded " << truth.size() << " ground truth entries\n";

    std::cerr << "Loading AGP results from " << gff_path << "...\n";
    auto agp = load_agp_gff(gff_path);
    std::cerr << "  Loaded " << agp.size() << " AGP predictions\n";

    std::unordered_map<std::string, FASTAEntry> corrected_seqs;
    if (!fasta_path.empty()) {
        std::cerr << "Loading corrected proteins from " << fasta_path << "...\n";
        corrected_seqs = load_fasta(fasta_path);
        std::cerr << "  Loaded " << corrected_seqs.size() << " sequences\n";
    }

    // Counters
    long long total_reads = 0;
    long long reads_with_damage = 0;
    long long reads_corrected = 0;

    // Frame/strand accuracy
    long long strand_correct = 0;
    long long strand_wrong = 0;
    long long strand_unknown = 0;
    long long frame_correct = 0;   // Protein length matches within tolerance
    long long frame_wrong = 0;     // Protein length differs significantly

    // Hamming distance (ONLY for correct-frame reads)
    long long total_baseline_dist = 0;
    long long total_corrected_dist = 0;
    long long total_positions = 0;

    // Per-read improvements (correct frame only)
    long long improved_reads = 0;
    long long worsened_reads = 0;
    long long unchanged_reads = 0;
    long long positions_fixed = 0;
    long long positions_broken = 0;

    for (auto& [name, gt] : truth) {
        if (gt.original_aa.empty()) continue;
        if (agp.find(name) == agp.end()) continue;

        total_reads++;
        auto& a = agp[name];

        // Get sequences
        std::string original = gt.original_aa;
        std::string damaged = gt.damaged_aa.empty() ? original : gt.damaged_aa;
        std::string corrected;
        char agp_strand = a.strand;

        // Get corrected protein from FASTA if available
        if (corrected_seqs.count(name)) {
            corrected = corrected_seqs[name].sequence;
            if (agp_strand == '?') {
                agp_strand = corrected_seqs[name].strand;
            }
        }

        // Check strand accuracy
        if (gt.strand != '?' && agp_strand != '?') {
            if (gt.strand == agp_strand) {
                strand_correct++;
            } else {
                strand_wrong++;
            }
        } else {
            strand_unknown++;
        }

        // Check frame accuracy by comparing protein lengths
        bool is_correct_frame = false;
        if (!corrected.empty()) {
            is_correct_frame = same_frame(corrected, original);
            if (is_correct_frame) {
                frame_correct++;
            } else {
                frame_wrong++;
            }
        }

        // Count damage and corrections
        int baseline_dist = hamming_distance(damaged, original);
        if (baseline_dist > 0) {
            reads_with_damage++;
        }
        if (a.aa_corrections > 0) {
            reads_corrected++;
        }

        // Only compute correction quality for correct-frame predictions
        if (is_correct_frame) {
            int corrected_dist = hamming_distance(corrected, original);

            total_baseline_dist += baseline_dist;
            total_corrected_dist += corrected_dist;
            total_positions += original.size();

            // Track improvement
            int improvement = baseline_dist - corrected_dist;
            if (improvement > 0) {
                improved_reads++;
                positions_fixed += improvement;
            } else if (improvement < 0) {
                worsened_reads++;
                positions_broken += (-improvement);
            } else {
                unchanged_reads++;
            }
        }
    }

    // Output results
    std::cout << "============================================================\n";
    std::cout << "AGP Correction Benchmark - Hamming Distance Analysis\n";
    std::cout << "============================================================\n\n";

    std::cout << "Dataset Overview:\n";
    std::cout << "  Total reads in ground truth: " << truth.size() << "\n";
    std::cout << "  Reads matched to AGP output: " << total_reads << "\n";
    std::cout << "  Reads with damage: " << reads_with_damage;
    if (total_reads > 0) {
        std::cout << " (" << (100.0 * reads_with_damage / total_reads) << "%)";
    }
    std::cout << "\n";
    std::cout << "  Reads AGP corrected: " << reads_corrected;
    if (total_reads > 0) {
        std::cout << " (" << (100.0 * reads_corrected / total_reads) << "%)";
    }
    std::cout << "\n\n";

    std::cout << "============================================================\n";
    std::cout << "FRAME/STRAND ACCURACY (critical for valid comparison)\n";
    std::cout << "============================================================\n\n";

    std::cout << "Strand Selection:\n";
    std::cout << "  Correct strand: " << strand_correct;
    if (strand_correct + strand_wrong > 0) {
        std::cout << " (" << (100.0 * strand_correct / (strand_correct + strand_wrong)) << "%)";
    }
    std::cout << "\n";
    std::cout << "  Wrong strand: " << strand_wrong << "\n";
    std::cout << "  Unknown: " << strand_unknown << "\n\n";

    long long total_frame_checked = frame_correct + frame_wrong;
    std::cout << "Frame Selection (protein length match):\n";
    std::cout << "  Correct frame: " << frame_correct;
    if (total_frame_checked > 0) {
        std::cout << " (" << (100.0 * frame_correct / total_frame_checked) << "%)";
    }
    std::cout << "\n";
    std::cout << "  Wrong frame: " << frame_wrong;
    if (total_frame_checked > 0) {
        std::cout << " (" << (100.0 * frame_wrong / total_frame_checked) << "%)";
    }
    std::cout << "\n\n";

    std::cout << "============================================================\n";
    std::cout << "CORRECTION QUALITY (correct-frame reads only: " << frame_correct << ")\n";
    std::cout << "============================================================\n\n";

    if (frame_correct == 0) {
        std::cout << "WARNING: No correct-frame predictions found!\n";
        std::cout << "Cannot compute correction quality metrics.\n\n";
        std::cout << "This means AGP is selecting different frames than the ground truth.\n";
        std::cout << "Possible causes:\n";
        std::cout << "  1. Frame selection algorithm favors different frames\n";
        std::cout << "  2. Scoring weights need adjustment\n";
        std::cout << "  3. Ground truth uses non-standard frame assignment\n";
    } else {
        std::cout << "Hamming Distance (lower is better):\n";
        std::cout << "  Total positions: " << total_positions << "\n";
        std::cout << "  Baseline (damaged vs original): " << total_baseline_dist;
        if (total_positions > 0) {
            std::cout << " (" << (100.0 * total_baseline_dist / total_positions) << "% error)";
        }
        std::cout << "\n";
        std::cout << "  After AGP (corrected vs original): " << total_corrected_dist;
        if (total_positions > 0) {
            std::cout << " (" << (100.0 * total_corrected_dist / total_positions) << "% error)";
        }
        std::cout << "\n\n";

        double improvement_pct = 0;
        if (total_baseline_dist > 0) {
            improvement_pct = 100.0 * (total_baseline_dist - total_corrected_dist) / total_baseline_dist;
        }
        std::cout << "Improvement:\n";
        std::cout << "  Positions fixed: " << positions_fixed << "\n";
        std::cout << "  Positions broken: " << positions_broken << "\n";
        std::cout << "  Net improvement: " << (total_baseline_dist - total_corrected_dist)
                  << " positions (" << improvement_pct << "% error reduction)\n\n";

        std::cout << "Per-read breakdown (correct frame only):\n";
        std::cout << "  Improved reads: " << improved_reads << "\n";
        std::cout << "  Worsened reads: " << worsened_reads << "\n";
        std::cout << "  Unchanged reads: " << unchanged_reads << "\n";
    }

    if (!fasta_path.empty() && corrected_seqs.empty()) {
        std::cerr << "\nWARNING: FASTA file specified but no sequences loaded.\n";
        std::cerr << "Make sure the read names match between GFF and FASTA.\n";
    }

    if (corrected_seqs.empty()) {
        std::cout << "\nNOTE: No corrected protein FASTA provided.\n";
        std::cout << "Run AGP with --fasta-aa to output corrected proteins for full analysis.\n";
    }

    return 0;
}
