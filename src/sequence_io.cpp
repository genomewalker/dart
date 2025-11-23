#include "agp/sequence_io.hpp"
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <zlib.h>

namespace agp {

// Large I/O buffer for better throughput (4MB like fastq-rmdup)
constexpr size_t GZBUF_SIZE = 4 * 1024 * 1024;

// Check if a command exists in PATH
static bool command_exists(const char* cmd) {
    std::string check = "which " + std::string(cmd) + " > /dev/null 2>&1";
    return system(check.c_str()) == 0;
}

// SequenceReader implementation
class SequenceReader::Impl {
public:
    std::ifstream file_;
    gzFile gz_file_ = nullptr;
    FILE* pipe_file_ = nullptr;  // For pigz
    Format format_ = Format::UNKNOWN;
    bool is_gzipped_ = false;
    bool use_pigz_ = false;
    char buffer_[65536];  // Line buffer (separate from I/O buffer)
    std::string lookahead_line_;  // For FASTA multi-line handling
    bool has_lookahead_ = false;

    bool open(const std::string& filename) {
        // Check if gzipped
        if (filename.size() > 3 &&
            filename.substr(filename.size() - 3) == ".gz") {
            is_gzipped_ = true;

            // Try pigz first for parallel decompression (disable with AGP_NO_PIGZ=1)
            // Use AGP_PIGZ_THREADS to set threads (default: 4)
            const char* no_pigz = std::getenv("AGP_NO_PIGZ");
            if (!no_pigz && command_exists("pigz")) {
                const char* pigz_threads_env = std::getenv("AGP_PIGZ_THREADS");
                int pigz_threads = pigz_threads_env ? std::atoi(pigz_threads_env) : 4;
                if (pigz_threads < 1) pigz_threads = 4;
                std::string cmd = "pigz -dc -p " + std::to_string(pigz_threads) + " \"" + filename + "\"";
                pipe_file_ = popen(cmd.c_str(), "r");
                if (pipe_file_) {
                    use_pigz_ = true;
                    // Set large buffer for pipe reads (4MB)
                    setvbuf(pipe_file_, nullptr, _IOFBF, GZBUF_SIZE);
                    // Read first line to determine format
                    if (fgets(buffer_, sizeof(buffer_), pipe_file_)) {
                        format_ = buffer_[0] == '>' ? Format::FASTA :
                                 buffer_[0] == '@' ? Format::FASTQ : Format::UNKNOWN;
                        // Reopen to reset position
                        pclose(pipe_file_);
                        pipe_file_ = popen(cmd.c_str(), "r");
                        if (pipe_file_) {
                            setvbuf(pipe_file_, nullptr, _IOFBF, GZBUF_SIZE);
                        }
                    }
                    return pipe_file_ != nullptr;
                }
            }

            // Fallback to zlib
            use_pigz_ = false;
            gz_file_ = gzopen(filename.c_str(), "rb");
            if (!gz_file_) return false;

            // Set large buffer for zlib (4MB)
            gzbuffer(gz_file_, GZBUF_SIZE);

            // Read first line to determine format
            if (gzgets(gz_file_, buffer_, sizeof(buffer_))) {
                format_ = buffer_[0] == '>' ? Format::FASTA :
                         buffer_[0] == '@' ? Format::FASTQ : Format::UNKNOWN;
                gzrewind(gz_file_);
            }
        } else {
            file_.open(filename);
            if (!file_) return false;

            // Determine format from first character
            char c = file_.peek();
            format_ = c == '>' ? Format::FASTA :
                     c == '@' ? Format::FASTQ : Format::UNKNOWN;
        }

        return true;
    }

    // Read line into existing string (avoids allocation if string has capacity)
    bool getline(std::string& line) {
        if (is_gzipped_) {
            if (use_pigz_) {
                if (fgets(buffer_, sizeof(buffer_), pipe_file_)) {
                    size_t len = strlen(buffer_);
                    if (len > 0 && buffer_[len-1] == '\n') len--;
                    line.assign(buffer_, len);
                    return true;
                }
                return false;
            } else {
                if (gzgets(gz_file_, buffer_, sizeof(buffer_))) {
                    size_t len = strlen(buffer_);
                    if (len > 0 && buffer_[len-1] == '\n') len--;
                    line.assign(buffer_, len);
                    return true;
                }
                return false;
            }
        } else {
            return static_cast<bool>(std::getline(file_, line));
        }
    }

    // Fast path: read directly into buffer, return pointer and length (no string copy)
    const char* getline_fast(size_t& len) {
        if (is_gzipped_) {
            char* result = use_pigz_ ?
                fgets(buffer_, sizeof(buffer_), pipe_file_) :
                gzgets(gz_file_, buffer_, sizeof(buffer_));
            if (result) {
                len = strlen(buffer_);
                if (len > 0 && buffer_[len-1] == '\n') len--;
                return buffer_;
            }
            len = 0;
            return nullptr;
        }
        // For uncompressed, fall back to regular getline
        len = 0;
        return nullptr;
    }

    bool is_open() const {
        if (is_gzipped_) {
            return use_pigz_ ? (pipe_file_ != nullptr) : (gz_file_ != nullptr);
        }
        return file_.is_open();
    }

    void close() {
        if (is_gzipped_) {
            if (use_pigz_ && pipe_file_) {
                pclose(pipe_file_);
                pipe_file_ = nullptr;
            } else if (gz_file_) {
                gzclose(gz_file_);
                gz_file_ = nullptr;
            }
        } else {
            file_.close();
        }
    }

    ~Impl() {
        close();
    }
};

SequenceReader::SequenceReader(const std::string& filename)
    : impl_(std::make_unique<Impl>()) {
    if (!impl_->open(filename)) {
        throw std::runtime_error("Failed to open file: " + filename);
    }
}

SequenceReader::~SequenceReader() = default;

bool SequenceReader::read_next(SequenceRecord& record) {
    std::string line;

    if (impl_->format_ == Format::FASTA) {
        // Read FASTA record - use lookahead if available
        if (impl_->has_lookahead_) {
            line = std::move(impl_->lookahead_line_);
            impl_->has_lookahead_ = false;
        } else {
            if (!impl_->getline(line)) {
                return false;
            }
        }

        // Skip empty lines
        while (line.empty()) {
            if (!impl_->getline(line)) return false;
        }

        if (line[0] != '>') {
            return false;
        }

        // Parse header - avoid substr copies where possible
        const char* hdr = line.c_str() + 1;  // Skip '>'
        const char* space = strchr(hdr, ' ');
        if (space) {
            record.id.assign(hdr, space - hdr);
            record.description.assign(space + 1);
        } else {
            record.id.assign(hdr);
            record.description.clear();
        }

        // Read sequence lines until next header or EOF
        record.sequence.clear();
        record.sequence.reserve(256);  // Pre-allocate for typical read length
        while (impl_->getline(line)) {
            if (line.empty()) continue;
            if (line[0] == '>') {
                // Save for next call
                impl_->lookahead_line_ = std::move(line);
                impl_->has_lookahead_ = true;
                break;
            }
            record.sequence += line;
        }

        return true;

    } else if (impl_->format_ == Format::FASTQ) {
        // FASTQ: 4 lines per record
        // Try fast path for gzipped files (avoid string copies)
        size_t len;
        const char* ptr;

        // Line 1: Header (@id description)
        if (impl_->is_gzipped_) {
            ptr = impl_->getline_fast(len);
            if (!ptr || len == 0 || ptr[0] != '@') return false;

            // Parse header directly from buffer
            const char* hdr = ptr + 1;  // Skip '@'
            const char* space = static_cast<const char*>(memchr(hdr, ' ', len - 1));
            if (space) {
                record.id.assign(hdr, space - hdr);
                record.description.assign(space + 1, len - 1 - (space - ptr));
            } else {
                record.id.assign(hdr, len - 1);
                record.description.clear();
            }

            // Line 2: Sequence
            ptr = impl_->getline_fast(len);
            if (!ptr) return false;
            record.sequence.assign(ptr, len);

            // Line 3: '+' line (skip)
            ptr = impl_->getline_fast(len);
            if (!ptr) return false;

            // Line 4: Quality
            ptr = impl_->getline_fast(len);
            if (!ptr) return false;
            record.quality.assign(ptr, len);
        } else {
            // Uncompressed: use regular getline
            if (!impl_->getline(line) || line.empty() || line[0] != '@') {
                return false;
            }

            // Parse header
            const char* hdr = line.c_str() + 1;
            const char* space = strchr(hdr, ' ');
            if (space) {
                record.id.assign(hdr, space - hdr);
                record.description.assign(space + 1);
            } else {
                record.id.assign(hdr);
                record.description.clear();
            }

            // Sequence
            if (!impl_->getline(record.sequence)) return false;

            // '+' line
            if (!impl_->getline(line)) return false;

            // Quality
            if (!impl_->getline(record.quality)) return false;
        }

        return true;
    }

    return false;
}

void SequenceReader::for_each(std::function<void(const SequenceRecord&)> callback) {
    SequenceRecord record;
    while (read_next(record)) {
        callback(record);
    }
}

std::vector<SequenceRecord> SequenceReader::read_all() {
    std::vector<SequenceRecord> records;
    SequenceRecord record;
    while (read_next(record)) {
        records.push_back(record);
    }
    return records;
}

bool SequenceReader::is_open() const {
    return impl_->is_open();
}

SequenceReader::Format SequenceReader::get_format() const {
    return impl_->format_;
}

// GeneWriter implementation
class GeneWriter::Impl {
public:
    std::ofstream file_;
    FILE* pipe_file_ = nullptr;
    bool use_pigz_ = false;
    static constexpr size_t BUFFER_SIZE = 1024 * 1024;  // 1MB buffer
    char buffer_[BUFFER_SIZE];
    char fmt_buffer_[4096];  // Reusable formatting buffer

    void write_header() {
        if (use_pigz_) {
            fprintf(pipe_file_, "##gff-version 3\n");
        } else {
            file_ << "##gff-version 3\n";
        }
    }

    void write(const char* data, size_t len) {
        if (use_pigz_) {
            fwrite(data, 1, len, pipe_file_);
        } else {
            file_.write(data, len);
        }
    }

    void write(const std::string& str) {
        write(str.c_str(), str.size());
    }
};

GeneWriter::GeneWriter(const std::string& filename)
    : impl_(std::make_unique<Impl>()) {

    // Check if output should be gzipped
    bool want_gzip = (filename.size() > 3 &&
                      filename.substr(filename.size() - 3) == ".gz");

    if (want_gzip && command_exists("pigz")) {
        std::string cmd = "pigz -c > \"" + filename + "\"";
        impl_->pipe_file_ = popen(cmd.c_str(), "w");
        if (impl_->pipe_file_) {
            impl_->use_pigz_ = true;
            impl_->write_header();
            return;
        }
    }

    // Fallback to regular file (no compression if .gz requested without pigz)
    std::string out_filename = filename;
    if (want_gzip && !command_exists("pigz")) {
        // Remove .gz extension if pigz not available
        out_filename = filename.substr(0, filename.size() - 3);
        std::cerr << "Warning: pigz not found, writing uncompressed to " << out_filename << "\n";
    }

    impl_->file_.open(out_filename);
    if (!impl_->file_) {
        throw std::runtime_error("Failed to open file: " + out_filename);
    }
    // Use larger buffer for better I/O performance
    impl_->file_.rdbuf()->pubsetbuf(impl_->buffer_, GeneWriter::Impl::BUFFER_SIZE);
    impl_->write_header();
}

GeneWriter::~GeneWriter() = default;

void GeneWriter::write_genes(const std::string& sequence_id,
                             const std::vector<Gene>& genes,
                             const std::string& source) {
    int gene_num = 1;
    for (const auto& gene : genes) {
        write_gene(sequence_id, gene, gene_num++, source);
    }
}

void GeneWriter::write_gene(const std::string& sequence_id,
                            const Gene& gene,
                            int gene_num,
                            const std::string& source) {
    // Use snprintf to avoid ostringstream allocation overhead
    int len = snprintf(impl_->fmt_buffer_, sizeof(impl_->fmt_buffer_),
        "%s\t%s\tgene\t%u\t%u\t%g\t%c\t.\tID=gene%d;ancient_prob=%g;damage_pct=%g\n",
        sequence_id.c_str(),
        source.c_str(),
        gene.start + 1,  // GFF is 1-based
        gene.end,
        gene.score,
        gene.is_forward ? '+' : '-',
        gene_num,
        gene.ancient_prob,
        gene.damage_score);
    impl_->write(impl_->fmt_buffer_, len);
}

void GeneWriter::close() {
    if (impl_->use_pigz_ && impl_->pipe_file_) {
        pclose(impl_->pipe_file_);
        impl_->pipe_file_ = nullptr;
    } else {
        impl_->file_.close();
    }
}

// SequenceUtils implementation
std::string SequenceUtils::reverse_complement(const std::string& seq) {
    std::string rc = seq;
    std::reverse(rc.begin(), rc.end());
    for (char& c : rc) {
        c = complement(c);
    }
    return rc;
}

char SequenceUtils::complement(char nt) {
    switch(nt) {
        case 'A': return 'T';
        case 'T': return 'A';
        case 'C': return 'G';
        case 'G': return 'C';
        case 'a': return 't';
        case 't': return 'a';
        case 'c': return 'g';
        case 'g': return 'c';
        default: return nt;
    }
}

std::string SequenceUtils::translate(const std::string& seq, int frame) {
    std::string protein;
    for (size_t i = frame; i + 2 < seq.length(); i += 3) {
        char aa = CodonTable::translate_codon(seq[i], seq[i+1], seq[i+2]);
        if (aa == '*') break;
        protein += aa;
    }
    return protein;
}

std::vector<Gene> SequenceUtils::find_orfs(const std::string& seq,
                                          size_t min_length) {
    std::vector<Gene> orfs;

    // Find ORFs in forward strand
    for (int frame = 0; frame < 3; ++frame) {
        for (size_t i = frame; i + 2 < seq.length(); i += 3) {
            if (CodonTable::is_start_codon(seq[i], seq[i+1], seq[i+2])) {
                // Found start, look for stop
                size_t start = i;
                size_t end = start;

                for (size_t j = start + 3; j + 2 < seq.length(); j += 3) {
                    if (CodonTable::is_stop_codon(seq[j], seq[j+1], seq[j+2])) {
                        end = j + 3;
                        break;
                    }
                }

                if (end > start && (end - start) >= min_length) {
                    Gene gene;
                    gene.start = start;
                    gene.end = end;
                    gene.is_forward = true;
                    gene.sequence = seq.substr(start, end - start);
                    gene.protein = translate(gene.sequence, 0);
                    gene.score = 0.0f;
                    gene.ancient_prob = 0.5f;
                    orfs.push_back(gene);
                }
            }
        }
    }

    // Reverse strand
    std::string rc = reverse_complement(seq);
    size_t seq_len = seq.length();

    for (int frame = 0; frame < 3; ++frame) {
        for (size_t i = frame; i + 2 < rc.length(); i += 3) {
            if (CodonTable::is_start_codon(rc[i], rc[i+1], rc[i+2])) {
                size_t start = i;
                size_t end = start;

                for (size_t j = start + 3; j + 2 < rc.length(); j += 3) {
                    if (CodonTable::is_stop_codon(rc[j], rc[j+1], rc[j+2])) {
                        end = j + 3;
                        break;
                    }
                }

                if (end > start && (end - start) >= min_length) {
                    Gene gene;
                    gene.start = seq_len - end;
                    gene.end = seq_len - start;
                    gene.is_forward = false;
                    gene.sequence = seq.substr(gene.start, gene.end - gene.start);
                    gene.protein = translate(rc.substr(start, end - start), 0);
                    gene.score = 0.0f;
                    gene.ancient_prob = 0.5f;
                    orfs.push_back(gene);
                }
            }
        }
    }

    return orfs;
}

float SequenceUtils::gc_content(const std::string& seq) {
    size_t gc_count = 0;
    size_t total = 0;

    for (char c : seq) {
        if (c == 'G' || c == 'C' || c == 'g' || c == 'c') {
            gc_count++;
            total++;
        } else if (c == 'A' || c == 'T' || c == 'a' || c == 't') {
            total++;
        }
    }

    return total > 0 ? static_cast<float>(gc_count) / total : 0.0f;
}

bool SequenceUtils::is_valid_dna(const std::string& seq) {
    for (char c : seq) {
        if (c != 'A' && c != 'C' && c != 'G' && c != 'T' && c != 'N' &&
            c != 'a' && c != 'c' && c != 'g' && c != 't' && c != 'n') {
            return false;
        }
    }
    return true;
}

std::string SequenceUtils::clean(const std::string& seq) {
    // Fast-path: check if sequence is already clean (uppercase, no whitespace)
    // Most FASTQ sequences are already uppercase
    bool needs_cleaning = false;
    for (char c : seq) {
        if (c <= ' ' || (c >= 'a' && c <= 'z')) {
            needs_cleaning = true;
            break;
        }
    }

    if (!needs_cleaning) {
        return seq;  // Return copy (often elided by compiler)
    }

    // Slow path: actually clean the sequence
    std::string cleaned;
    cleaned.reserve(seq.length());

    for (char c : seq) {
        if (c <= ' ') continue;
        if (c >= 'a' && c <= 'z') {
            cleaned += (c - 32);
        } else {
            cleaned += c;
        }
    }
    return cleaned;
}

std::vector<Nucleotide> SequenceUtils::encode(const std::string& seq) {
    std::vector<Nucleotide> encoded;
    encoded.reserve(seq.size());
    for (char c : seq) {
        encoded.push_back(char_to_nt(c));
    }
    return encoded;
}

std::string SequenceUtils::decode(const std::vector<Nucleotide>& seq) {
    std::string decoded;
    decoded.reserve(seq.size());
    for (Nucleotide nt : seq) {
        decoded += nt_to_char(nt);
    }
    return decoded;
}

std::array<size_t, 4> SequenceUtils::composition(const std::string& seq) {
    std::array<size_t, 4> comp = {0, 0, 0, 0};
    for (char c : seq) {
        Nucleotide nt = char_to_nt(c);
        if (nt != Nucleotide::N) {
            comp[static_cast<size_t>(nt)]++;
        }
    }
    return comp;
}

// FastaWriter implementation
class FastaWriter::Impl {
public:
    std::ofstream file_;
    FILE* pipe_file_ = nullptr;
    bool use_pigz_ = false;
    int gene_counter_ = 0;
    static constexpr size_t BUFFER_SIZE = 1024 * 1024;  // 1MB buffer
    char buffer_[BUFFER_SIZE];
    char fmt_buffer_[4096];  // Reusable formatting buffer

    void write(const char* data, size_t len) {
        if (use_pigz_) {
            fwrite(data, 1, len, pipe_file_);
        } else {
            file_.write(data, len);
        }
    }

    void write(const std::string& str) {
        write(str.c_str(), str.size());
    }
};

FastaWriter::FastaWriter(const std::string& filename)
    : impl_(std::make_unique<Impl>()) {

    // Check if output should be gzipped
    bool want_gzip = (filename.size() > 3 &&
                      filename.substr(filename.size() - 3) == ".gz");

    if (want_gzip && command_exists("pigz")) {
        std::string cmd = "pigz -c > \"" + filename + "\"";
        impl_->pipe_file_ = popen(cmd.c_str(), "w");
        if (impl_->pipe_file_) {
            impl_->use_pigz_ = true;
            setvbuf(impl_->pipe_file_, nullptr, _IOFBF, GZBUF_SIZE);
            return;
        }
    }

    // Fallback to regular file
    std::string out_filename = filename;
    if (want_gzip && !command_exists("pigz")) {
        out_filename = filename.substr(0, filename.size() - 3);
        std::cerr << "Warning: pigz not found, writing uncompressed to " << out_filename << "\n";
    }

    impl_->file_.open(out_filename);
    if (!impl_->file_) {
        throw std::runtime_error("Failed to open FASTA file: " + out_filename);
    }
    impl_->file_.rdbuf()->pubsetbuf(impl_->buffer_, FastaWriter::Impl::BUFFER_SIZE);
}

FastaWriter::~FastaWriter() {
    if (impl_) {
        if (impl_->use_pigz_ && impl_->pipe_file_) {
            pclose(impl_->pipe_file_);
        } else if (impl_->file_.is_open()) {
            impl_->file_.close();
        }
    }
}

void FastaWriter::write_sequence(const std::string& id,
                                 const std::string& description,
                                 const std::string& sequence) {
    // Write header directly without ostringstream
    impl_->write(">", 1);
    impl_->write(id);
    if (!description.empty()) {
        impl_->write(" ", 1);
        impl_->write(description);
    }
    impl_->write("\n", 1);

    // Write sequence in lines of 80 characters - avoid substr() copies
    const char* seq_data = sequence.c_str();
    size_t seq_len = sequence.length();
    for (size_t i = 0; i < seq_len; i += 80) {
        size_t line_len = std::min(size_t(80), seq_len - i);
        impl_->write(seq_data + i, line_len);
        impl_->write("\n", 1);
    }
}

void FastaWriter::write_genes_nucleotide(const std::string& sequence_id,
                                         const std::vector<Gene>& genes) {
    for (size_t i = 0; i < genes.size(); ++i) {
        const auto& gene = genes[i];
        impl_->gene_counter_++;

        // Use snprintf to avoid ostringstream allocation
        char id_buf[64];
        snprintf(id_buf, sizeof(id_buf), "gene_%d", impl_->gene_counter_);

        int desc_len = snprintf(impl_->fmt_buffer_, sizeof(impl_->fmt_buffer_),
            "%s_gene_%zu %c %u..%u score=%g ancient_prob=%g damage_pct=%g",
            sequence_id.c_str(), i + 1,
            gene.is_forward ? '+' : '-',
            gene.start, gene.end,
            gene.score, gene.ancient_prob, gene.damage_score);

        write_sequence(id_buf, std::string(impl_->fmt_buffer_, desc_len), gene.sequence);
    }
}

void FastaWriter::write_genes_protein(const std::string& sequence_id,
                                      const std::vector<Gene>& genes) {
    for (size_t i = 0; i < genes.size(); ++i) {
        const auto& gene = genes[i];
        impl_->gene_counter_++;

        char id_buf[64];
        snprintf(id_buf, sizeof(id_buf), "protein_%d", impl_->gene_counter_);

        int desc_len;
        if (gene.dna_corrections > 0 || gene.aa_corrections > 0) {
            desc_len = snprintf(impl_->fmt_buffer_, sizeof(impl_->fmt_buffer_),
                "%s_gene_%zu %c %u..%u length=%zuaa damage_pct=%g dna_corr=%zu aa_corr=%zu",
                sequence_id.c_str(), i + 1,
                gene.is_forward ? '+' : '-',
                gene.start, gene.end,
                gene.protein.length(), gene.damage_score,
                gene.dna_corrections, gene.aa_corrections);
        } else {
            desc_len = snprintf(impl_->fmt_buffer_, sizeof(impl_->fmt_buffer_),
                "%s_gene_%zu %c %u..%u length=%zuaa damage_pct=%g",
                sequence_id.c_str(), i + 1,
                gene.is_forward ? '+' : '-',
                gene.start, gene.end,
                gene.protein.length(), gene.damage_score);
        }

        write_sequence(id_buf, std::string(impl_->fmt_buffer_, desc_len), gene.protein);

        // If corrections were made, also output the corrected protein
        if (!gene.corrected_protein.empty() && gene.aa_corrections > 0) {
            snprintf(id_buf, sizeof(id_buf), "protein_%d_corrected", impl_->gene_counter_);
            desc_len = snprintf(impl_->fmt_buffer_, sizeof(impl_->fmt_buffer_),
                "%s_gene_%zu_corrected %c %u..%u length=%zuaa damage_corrected=true",
                sequence_id.c_str(), i + 1,
                gene.is_forward ? '+' : '-',
                gene.start, gene.end,
                gene.corrected_protein.length());

            write_sequence(id_buf, std::string(impl_->fmt_buffer_, desc_len), gene.corrected_protein);
        }
    }
}

void FastaWriter::close() {
    if (impl_->use_pigz_ && impl_->pipe_file_) {
        pclose(impl_->pipe_file_);
        impl_->pipe_file_ = nullptr;
    } else if (impl_->file_.is_open()) {
        impl_->file_.close();
    }
}

// StatsWriter implementation
class StatsWriter::Impl {
public:
    std::ofstream file_;
};

StatsWriter::StatsWriter(const std::string& filename)
    : impl_(std::make_unique<Impl>()) {
    impl_->file_.open(filename);
    if (!impl_->file_) {
        throw std::runtime_error("Failed to open stats file: " + filename);
    }
    // Write header
    impl_->file_ << "# AGP Statistics Output\n";
    impl_->file_ << "# Generated by Ancient Gene Predictor\n\n";
}

StatsWriter::~StatsWriter() {
    if (impl_ && impl_->file_.is_open()) {
        impl_->file_.close();
    }
}

void StatsWriter::write_codon_frequencies(const std::map<std::string, float>& frequencies) {
    impl_->file_ << "## Codon Frequencies\n";
    impl_->file_ << "# Codon\tFrequency\n";

    for (const auto& [codon, freq] : frequencies) {
        impl_->file_ << codon << "\t" << freq << "\n";
    }
    impl_->file_ << "\n";
}

void StatsWriter::write_damage_parameters(float lambda_5p, float lambda_3p,
                                          float delta_max, float delta_bg) {
    impl_->file_ << "## Damage Model Parameters\n";
    impl_->file_ << "lambda_5prime\t" << lambda_5p << "\n";
    impl_->file_ << "lambda_3prime\t" << lambda_3p << "\n";
    impl_->file_ << "delta_max\t" << delta_max << "\n";
    impl_->file_ << "delta_background\t" << delta_bg << "\n";
    impl_->file_ << "\n";
}

void StatsWriter::write_damage_profile(const std::vector<float>& ct_profile,
                                       const std::vector<float>& ga_profile) {
    impl_->file_ << "## Position-Specific Damage Rates\n";
    impl_->file_ << "# Position\tC->T_Rate\tG->A_Rate\n";

    size_t max_len = std::max(ct_profile.size(), ga_profile.size());
    for (size_t i = 0; i < max_len; ++i) {
        impl_->file_ << i << "\t";
        if (i < ct_profile.size()) {
            impl_->file_ << ct_profile[i];
        } else {
            impl_->file_ << "NA";
        }
        impl_->file_ << "\t";
        if (i < ga_profile.size()) {
            impl_->file_ << ga_profile[i];
        } else {
            impl_->file_ << "NA";
        }
        impl_->file_ << "\n";
    }
    impl_->file_ << "\n";
}

void StatsWriter::write_organism_stats(float gc_content, float enc,
                                       const std::map<std::string, float>& rscu) {
    impl_->file_ << "## Organism Statistics\n";
    impl_->file_ << "GC_content\t" << gc_content << "\n";
    impl_->file_ << "ENC\t" << enc << "\n";
    impl_->file_ << "\n";

    impl_->file_ << "## Relative Synonymous Codon Usage (RSCU)\n";
    impl_->file_ << "# Codon\tRSCU\n";
    for (const auto& [codon, value] : rscu) {
        impl_->file_ << codon << "\t" << value << "\n";
    }
    impl_->file_ << "\n";
}

void StatsWriter::close() {
    if (impl_->file_.is_open()) {
        impl_->file_.close();
    }
}

} // namespace agp
