#include "dart/sequence_io.hpp"
#include "dart/gz_reader_base.hpp"
#include "dart/codon_tables.hpp"
#include <filesystem>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <stdexcept>
#include <vector>
#include <cstdint>
#include <zlib.h>
#include <isa-l/igzip_lib.h>

namespace dart {

// Large I/O buffer for better throughput (4MB like fastq-rmdup)
constexpr size_t GZBUF_SIZE = 4 * 1024 * 1024;

// In-process gzip sink using ISA-L igzip (isal_deflate). The single compression path
// for every DART writer -- there is no pigz subprocess anywhere. Streams: input is fed
// in chunks, compressed output flushed to the file as it is produced, gzip header and
// CRC trailer handled by IGZIP_GZIP. Matches fqdup's ISA-L backend.
class IgzipSink {
public:
    explicit IgzipSink(const std::string& path, int level = 1)
        : out_buf_(GZBUF_SIZE) {
        fp_ = std::fopen(path.c_str(), "wb");
        if (!fp_) throw std::runtime_error("Failed to open gzip output: " + path);
        isal_deflate_init(&zs_);
        zs_.gzip_flag = IGZIP_GZIP;
        static const size_t lvl_sz[] = {
            ISAL_DEF_LVL0_DEFAULT, ISAL_DEF_LVL1_DEFAULT,
            ISAL_DEF_LVL2_DEFAULT, ISAL_DEF_LVL3_DEFAULT};
        int lv = level < 0 ? 0 : (level > 3 ? 3 : level);
        zs_.level = static_cast<uint32_t>(lv);
        level_buf_.resize(lvl_sz[lv]);
        zs_.level_buf = level_buf_.data();
        zs_.level_buf_size = static_cast<uint32_t>(level_buf_.size());
    }
    IgzipSink(const IgzipSink&) = delete;
    IgzipSink& operator=(const IgzipSink&) = delete;

    void write(const char* data, size_t len) {
        const uint8_t* p = reinterpret_cast<const uint8_t*>(data);
        while (len > 0) {
            uint32_t chunk = len > 0x40000000u ? 0x40000000u : static_cast<uint32_t>(len);
            zs_.next_in = const_cast<uint8_t*>(p);
            zs_.avail_in = chunk;
            zs_.end_of_stream = 0;
            pump();
            p += chunk;
            len -= chunk;
        }
    }

    void close() {
        if (!fp_) return;
        zs_.next_in = nullptr;
        zs_.avail_in = 0;
        zs_.end_of_stream = 1;
        pump();
        std::fclose(fp_);
        fp_ = nullptr;
    }
    ~IgzipSink() { close(); }
    bool is_open() const { return fp_ != nullptr; }

private:
    void pump() {
        do {
            zs_.next_out = out_buf_.data();
            zs_.avail_out = static_cast<uint32_t>(out_buf_.size());
            isal_deflate(&zs_);
            size_t have = out_buf_.size() - zs_.avail_out;
            if (have && std::fwrite(out_buf_.data(), 1, have, fp_) != have)
                throw std::runtime_error("gzip write failed");
        } while (zs_.avail_out == 0 ||
                 (zs_.end_of_stream && zs_.internal_state.state != ZSTATE_END));
    }

    std::FILE* fp_ = nullptr;
    isal_zstream zs_{};
    std::vector<uint8_t> level_buf_;
    std::vector<uint8_t> out_buf_;
};

// SequenceReader implementation
class SequenceReader::Impl {
public:
    std::ifstream file_;
    gzFile gz_file_ = nullptr;
    std::unique_ptr<GzLineReader> rg_reader_;  // rapidgzip (ISA-L) backend — preferred
    Format format_ = Format::UNKNOWN;
    bool is_gzipped_ = false;
    char buffer_[65536];  // Line buffer (separate from I/O buffer)
    std::string lookahead_line_;  // For FASTA multi-line handling
    bool has_lookahead_ = false;

    bool open(const std::string& filename) {
        // Check if gzipped
        if (filename.size() > 3 &&
            filename.substr(filename.size() - 3) == ".gz") {
            is_gzipped_ = true;

            // Try rapidgzip first: in-process parallel decompression, no PATH dependency.
            rg_reader_ = make_gz_reader(filename);
            if (rg_reader_) {
                // Detect format from first line, then reopen to reset to BOF.
                std::string first_line;
                if (rg_reader_->readline(first_line) && !first_line.empty())
                    format_ = first_line[0] == '>' ? Format::FASTA :
                              first_line[0] == '@' ? Format::FASTQ : Format::UNKNOWN;
                rg_reader_ = make_gz_reader(filename);
                return bool(rg_reader_);
            }

            // rapidgzip (ISA-L) is the decompression path. If it could not open the
            // stream, fall back to zlib -- never a subprocess.
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
        line.clear();
        if (is_gzipped_ && rg_reader_) return rg_reader_->readline(line);
        if (is_gzipped_) {
            auto read_chunk = [&](auto getter) -> bool {
                while (true) {
                    if (!getter()) return !line.empty();
                    size_t raw_len = strlen(buffer_);
                    size_t chunk_len = strcspn(buffer_, "\n");
                    line.append(buffer_, chunk_len);
                    bool has_newline = (chunk_len < raw_len) && buffer_[chunk_len] == '\n';
                    if (has_newline || raw_len < sizeof(buffer_) - 1) return true;
                }
            };

            auto getter = [&]() -> bool { return gzgets(gz_file_, buffer_, sizeof(buffer_)) != nullptr; };
            return read_chunk(getter);
        } else {
            return static_cast<bool>(std::getline(file_, line));
        }
    }

    // Fast path: read directly into buffer, return pointer and length (no string copy)
    const char* getline_fast(size_t& len) {
        if (is_gzipped_ && rg_reader_) return rg_reader_->readline_fast(len);

        static thread_local std::string dyn_buffer;

        if (is_gzipped_) {
            auto read_into_dyn = [&](auto getter) -> const char* {
                dyn_buffer.clear();
                while (true) {
                    if (!getter()) break;
                    size_t raw_len = strlen(buffer_);
                    size_t chunk_len = strcspn(buffer_, "\n");
                    dyn_buffer.append(buffer_, chunk_len);
                    bool has_newline = (chunk_len < raw_len) && buffer_[chunk_len] == '\n';
                    if (has_newline || raw_len < sizeof(buffer_) - 1) break;
                }
                if (dyn_buffer.empty()) {
                    len = 0;
                    return nullptr;
                }
                len = dyn_buffer.size();
                return dyn_buffer.c_str();
            };

            auto getter_once = [&](auto getter) -> const char* {
                if (!getter()) {
                    len = 0;
                    return nullptr;
                }
                size_t raw_len = strlen(buffer_);
                size_t chunk_len = strcspn(buffer_, "\n");
                bool has_newline = (chunk_len < raw_len) && buffer_[chunk_len] == '\n';
                if (has_newline || raw_len < sizeof(buffer_) - 1) {
                    len = chunk_len;
                    return buffer_;
                }
                // Overflow: accumulate dynamically
                return read_into_dyn(getter);
            };

            auto getter = [&]() -> bool { return gzgets(gz_file_, buffer_, sizeof(buffer_)) != nullptr; };
            return getter_once(getter);
        }
        // For uncompressed, fall back to regular getline
        len = 0;
        return nullptr;
    }

    bool is_open() const {
        if (is_gzipped_) return rg_reader_ ? true : (gz_file_ != nullptr);
        return file_.is_open();
    }

    void close() {
        rg_reader_.reset();
        if (is_gzipped_) {
            if (gz_file_) {
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

// SequenceWriter implementation
class SequenceWriter::Impl {
public:
    std::ofstream file_;
    std::unique_ptr<IgzipSink> gz_;  // ISA-L igzip when the output is .gz
    bool is_fastq_ = false;
    bool format_locked_ = false;
    static constexpr size_t BUFFER_SIZE = 1024 * 1024;  // 1MB buffer
    char buffer_[BUFFER_SIZE];

    void write(const char* data, size_t len) {
        if (gz_) {
            gz_->write(data, len);
        } else {
            file_.write(data, len);
        }
    }

    void write(const std::string& str) {
        write(str.c_str(), str.size());
    }

    void write_fasta_wrapped(const std::string& sequence) {
        const char* seq_data = sequence.c_str();
        const size_t seq_len = sequence.length();
        for (size_t i = 0; i < seq_len; i += 80) {
            const size_t line_len = std::min(size_t(80), seq_len - i);
            write(seq_data + i, line_len);
            write("\n", 1);
        }
    }
};

SequenceWriter::SequenceWriter(const std::string& filename, bool compress)
    : impl_(std::make_unique<Impl>()) {

    const bool want_gzip = compress ||
        (filename.size() > 3 && filename.substr(filename.size() - 3) == ".gz");

    if (auto p = std::filesystem::path(filename).parent_path(); !p.empty())
        std::filesystem::create_directories(p);

    if (want_gzip) {
        impl_->gz_ = std::make_unique<IgzipSink>(filename);
        return;
    }

    impl_->file_.open(filename);
    if (!impl_->file_) {
        throw std::runtime_error("Failed to open sequence output file: " + filename);
    }
    impl_->file_.rdbuf()->pubsetbuf(impl_->buffer_, SequenceWriter::Impl::BUFFER_SIZE);
}

SequenceWriter::~SequenceWriter() {
    if (impl_) {
        close();
    }
}

void SequenceWriter::write(const SequenceRecord& record) {
    const bool record_is_fastq = !record.quality.empty();

    if (!impl_->format_locked_) {
        impl_->format_locked_ = true;
        impl_->is_fastq_ = record_is_fastq;
    } else if (impl_->is_fastq_ != record_is_fastq) {
        throw std::runtime_error("SequenceWriter cannot mix FASTA and FASTQ records");
    }

    if (impl_->is_fastq_) {
        if (record.quality.size() != record.sequence.size()) {
            throw std::runtime_error(
                "FASTQ quality length must match sequence length for record: " + record.id);
        }
        impl_->write("@", 1);
        impl_->write(record.id);
        if (!record.description.empty()) {
            impl_->write(" ", 1);
            impl_->write(record.description);
        }
        impl_->write("\n", 1);
        impl_->write(record.sequence);
        impl_->write("\n+\n", 3);
        impl_->write(record.quality);
        impl_->write("\n", 1);
        return;
    }

    impl_->write(">", 1);
    impl_->write(record.id);
    if (!record.description.empty()) {
        impl_->write(" ", 1);
        impl_->write(record.description);
    }
    impl_->write("\n", 1);
    impl_->write_fasta_wrapped(record.sequence);
}

void SequenceWriter::write_all(const std::vector<SequenceRecord>& records) {
    for (const auto& record : records) {
        write(record);
    }
}

void SequenceWriter::close() {
    if (impl_->gz_) {
        impl_->gz_->close();
        impl_->gz_.reset();
    } else if (impl_->file_.is_open()) {
        impl_->file_.close();
    }
}

// GeneWriter implementation
class GeneWriter::Impl {
public:
    std::ofstream file_;
    std::unique_ptr<IgzipSink> gz_;  // ISA-L igzip when the output is .gz
    static constexpr size_t BUFFER_SIZE = 1024 * 1024;  // 1MB buffer
    char buffer_[BUFFER_SIZE];
    char fmt_buffer_[4096];  // Reusable formatting buffer

    void write_header() {
        static const char h[] = "##gff-version 3\n";
        write(h, sizeof(h) - 1);
    }

    void write(const char* data, size_t len) {
        if (gz_) {
            gz_->write(data, len);
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

    if (auto p = std::filesystem::path(filename).parent_path(); !p.empty())
        std::filesystem::create_directories(p);

    if (want_gzip) {
        impl_->gz_ = std::make_unique<IgzipSink>(filename);
        impl_->write_header();
        return;
    }

    impl_->file_.open(filename);
    if (!impl_->file_) {
        throw std::runtime_error("Failed to open file: " + filename);
    }
    // Use larger buffer for better I/O performance
    impl_->file_.rdbuf()->pubsetbuf(impl_->buffer_, GeneWriter::Impl::BUFFER_SIZE);
    impl_->write_header();
}

GeneWriter::~GeneWriter() {
    if (impl_) {
        close();
    }
}

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
    // Include frame= attribute for validation tools
    // frame_pmax = best 6-way posterior, frame_margin = confidence gap (p_max - p_second)
    // log_lr and dmg_info for downstream aggregation across BLAST hits
    int len;

    // Check if this is an ORF mode output (has ORF-specific fields set)
    bool is_orf_mode = (gene.orf_length_score > 0.0f || gene.passed_stops > 0 || gene.x_count > 0);

    if (is_orf_mode) {
        // ORF mode: include orf_rank, passed_stops, x_count, orf_length_score, per_read_damage
        // Also include aa_corr and stop_fix for consistency with FASTA headers
        len = snprintf(impl_->fmt_buffer_, sizeof(impl_->fmt_buffer_),
            "%s\t%s\tgene\t%u\t%u\t%g\t%c\t.\tID=gene%d;frame=%d;orf_rank=%zu;orf_len=%g;aa_corr=%zu;stop_fix=%zu;passed_stops=%zu;x_count=%zu;damage_signal=%g;damage_pct=%g;damage_evidence=%g;p_damaged=%g\n",
            sequence_id.c_str(),
            source.c_str(),
            gene.start + 1,  // GFF is 1-based
            gene.end,
            gene.score,
            gene.is_forward ? '+' : '-',
            gene_num,
            gene.frame,
            gene.orf_rank,
            gene.orf_length_score,
            gene.aa_corrections,
            gene.stop_restorations,
            gene.passed_stops,
            gene.x_count,
            gene.damage_signal,
            gene.damage_score,
            gene.per_read_damage_evidence,
            gene.p_read_damaged);
    } else {
        // Standard mode
        // Calculate confidence score from protein length
        float confidence = calculate_confidence(gene.protein.length(), gene.coding_prob);

        len = snprintf(impl_->fmt_buffer_, sizeof(impl_->fmt_buffer_),
            "%s\t%s\tgene\t%u\t%u\t%g\t%c\t.\tID=gene%d;frame=%d;conf=%g;frame_pmax=%g;frame_margin=%g;damage_signal=%g;damage_pct=%g;damage_pctile=%g;p_damaged=%g;log_lr=%g;dmg_info=%g;dna_corr=%zu;aa_corr=%zu;stop_fix=%zu;int_stops=%zu;corr_score=%g\n",
            sequence_id.c_str(),
            source.c_str(),
            gene.start + 1,  // GFF is 1-based
            gene.end,
            gene.score,
            gene.is_forward ? '+' : '-',
            gene_num,
            gene.frame,
            confidence,
            gene.frame_pmax,
            gene.frame_margin,
            gene.damage_signal,
            gene.damage_score,
            gene.damage_pctile,
            gene.p_read_damaged,
            gene.damage_log_lr,
            gene.damage_info,
            gene.dna_corrections,
            gene.aa_corrections,
            gene.stop_restorations,
            gene.internal_stops,
            gene.corrected_coding_prob);
    }
    // Clamp len to buffer size to prevent read past buffer
    if (len > 0 && static_cast<size_t>(len) >= sizeof(impl_->fmt_buffer_)) {
        len = sizeof(impl_->fmt_buffer_) - 1;
    }
    impl_->write(impl_->fmt_buffer_, len);
}

void GeneWriter::close() {
    if (impl_->gz_) {
        impl_->gz_->close();
        impl_->gz_.reset();
    } else {
        impl_->file_.close();
    }
}

// SequenceUtils implementation - delegates to canonical reverse_complement
std::string SequenceUtils::reverse_complement(const std::string& seq) {
    return dart::reverse_complement(seq);
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
                    gene.damage_signal = 0.5f;
                    orfs.push_back(std::move(gene));
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
                    gene.damage_signal = 0.5f;
                    orfs.push_back(std::move(gene));
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
    std::unique_ptr<IgzipSink> gz_;  // ISA-L igzip when the output is .gz
    int gene_counter_ = 0;
    static constexpr size_t BUFFER_SIZE = 1024 * 1024;  // 1MB buffer
    char buffer_[BUFFER_SIZE];
    char fmt_buffer_[4096];  // Reusable formatting buffer

    void write(const char* data, size_t len) {
        if (gz_) {
            gz_->write(data, len);
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

    if (auto p = std::filesystem::path(filename).parent_path(); !p.empty())
        std::filesystem::create_directories(p);

    if (want_gzip) {
        impl_->gz_ = std::make_unique<IgzipSink>(filename);
        return;
    }

    impl_->file_.open(filename);
    if (!impl_->file_) {
        throw std::runtime_error("Failed to open FASTA file: " + filename);
    }
    impl_->file_.rdbuf()->pubsetbuf(impl_->buffer_, FastaWriter::Impl::BUFFER_SIZE);
}

FastaWriter::~FastaWriter() {
    if (impl_) {
        close();
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

        // Always output original (uncorrected) sequence
        // Use write_genes_nucleotide_corrected() for damage-corrected nucleotides
        const std::string& output_seq = gene.sequence;

        // Skip genes with empty sequences (no valid ORF found)
        if (output_seq.empty()) continue;

        impl_->gene_counter_++;

        // Use snprintf to avoid ostringstream allocation
        char id_buf[64];
        snprintf(id_buf, sizeof(id_buf), "gene_%d", impl_->gene_counter_);

        int desc_len;
        if (gene.dna_corrections > 0) {
            desc_len = snprintf(impl_->fmt_buffer_, sizeof(impl_->fmt_buffer_),
                "%s_gene_%zu %c %u..%u score=%g damage_signal=%g damage_pct=%g p_damaged=%g log_lr=%g dmg_info=%g dna_corr=%zu stop_fix=%zu corr_score=%g",
                sequence_id.c_str(), i + 1,
                gene.is_forward ? '+' : '-',
                gene.start, gene.end,
                gene.score, gene.damage_signal, gene.damage_score, gene.p_read_damaged,
                gene.damage_log_lr, gene.damage_info,
                gene.dna_corrections, gene.stop_restorations, gene.corrected_coding_prob);
        } else {
            desc_len = snprintf(impl_->fmt_buffer_, sizeof(impl_->fmt_buffer_),
                "%s_gene_%zu %c %u..%u score=%g damage_signal=%g damage_pct=%g p_damaged=%g log_lr=%g dmg_info=%g stop_fix=%zu",
                sequence_id.c_str(), i + 1,
                gene.is_forward ? '+' : '-',
                gene.start, gene.end,
                gene.score, gene.damage_signal, gene.damage_score, gene.p_read_damaged,
                gene.damage_log_lr, gene.damage_info, gene.stop_restorations);
        }

        write_sequence(id_buf, std::string(impl_->fmt_buffer_, desc_len), output_seq);
    }
}

// Extract the longest ORF (before first stop codon) from a protein sequence
static std::string extract_longest_orf(const std::string& protein) {
    // Find first stop codon - output only the ORF before it
    // This is standard behavior for protein FASTA output (MMseqs2, BLAST compatibility)
    size_t stop_pos = protein.find('*');
    if (stop_pos == std::string::npos) {
        return protein;  // No stops, return full sequence
    }
    return protein.substr(0, stop_pos);
}

void FastaWriter::write_genes_protein(const std::string& sequence_id,
                                      const std::vector<Gene>& genes) {
    for (size_t i = 0; i < genes.size(); ++i) {
        const auto& gene = genes[i];

        // Extract longest ORF (before first stop) for MMseqs2/BLAST compatibility
        // Internal stops indicate wrong frame or damage - output clean ORF
        const std::string output_protein = extract_longest_orf(gene.protein);

        // Skip genes with empty proteins (no valid ORF found)
        if (output_protein.empty()) continue;

        impl_->gene_counter_++;

        // Use read name + strand + gene index as ID for MMseqs2 compatibility
        // Format: READNAME_+_1 or READNAME_-_1 (allows extracting read name by splitting on _+_ or _-_)
        char id_buf[1024];
        snprintf(id_buf, sizeof(id_buf), "%s_%c_%zu",
            sequence_id.c_str(),
            gene.is_forward ? '+' : '-',
            i + 1);

        int desc_len;
        desc_len = snprintf(impl_->fmt_buffer_, sizeof(impl_->fmt_buffer_),
            "%u..%u length=%zuaa damage_pct=%g p_damaged=%g log_lr=%g dmg_info=%g aa_corr=%zu stop_fix=%zu corr_score=%g",
            gene.start, gene.end,
            output_protein.length(), gene.damage_score, gene.p_read_damaged,
            gene.damage_log_lr, gene.damage_info,
            gene.aa_corrections, gene.stop_restorations, gene.corrected_coding_prob);

        write_sequence(id_buf, std::string(impl_->fmt_buffer_, desc_len), output_protein);
    }
}

void FastaWriter::write_genes_protein_corrected(const std::string& sequence_id,
                                                const std::vector<Gene>& genes) {
    for (size_t i = 0; i < genes.size(); ++i) {
        const auto& gene = genes[i];

        // Output the corrected protein (for benchmarking damage correction)
        // Falls back to original if no correction was made
        // Extract longest ORF (before first stop) for MMseqs2/BLAST compatibility
        const std::string output_protein = extract_longest_orf(
            gene.corrected_protein.empty() ? gene.protein : gene.corrected_protein);

        // Skip genes with empty proteins (no valid ORF found)
        if (output_protein.empty()) continue;

        impl_->gene_counter_++;

        // Use read name + strand + gene index as ID for MMseqs2 compatibility
        char id_buf[1024];
        snprintf(id_buf, sizeof(id_buf), "%s_%c_%zu",
            sequence_id.c_str(),
            gene.is_forward ? '+' : '-',
            i + 1);

        // Calculate confidence score from protein length and coding probability
        float confidence = calculate_confidence(output_protein.length(), gene.coding_prob);

        int desc_len;
        desc_len = snprintf(impl_->fmt_buffer_, sizeof(impl_->fmt_buffer_),
            "%u..%u length=%zuaa conf=%.3f damage_pct=%g aa_corr=%zu stop_fix=%zu",
            gene.start, gene.end,
            output_protein.length(), confidence, gene.damage_score,
            gene.aa_corrections, gene.stop_restorations);

        write_sequence(id_buf, std::string(impl_->fmt_buffer_, desc_len), output_protein);
    }
}

void FastaWriter::write_genes_nucleotide_corrected(const std::string& sequence_id,
                                                    const std::vector<Gene>& genes) {
    for (size_t i = 0; i < genes.size(); ++i) {
        const auto& gene = genes[i];

        // Output the corrected nucleotide sequence (for benchmarking damage correction)
        // Falls back to original if no correction was made
        const std::string& output_seq = gene.corrected_sequence.empty() ?
            gene.sequence : gene.corrected_sequence;

        // Skip genes with empty sequences (no valid ORF found)
        if (output_seq.empty()) continue;

        impl_->gene_counter_++;

        char id_buf[64];
        snprintf(id_buf, sizeof(id_buf), "corrected_nt_%d", impl_->gene_counter_);

        int desc_len;
        desc_len = snprintf(impl_->fmt_buffer_, sizeof(impl_->fmt_buffer_),
            "%s_gene_%zu %c %u..%u length=%zubp damage_pct=%g dna_corr=%zu stop_fix=%zu",
            sequence_id.c_str(), i + 1,
            gene.is_forward ? '+' : '-',
            gene.start, gene.end,
            output_seq.length(), gene.damage_score,
            gene.dna_corrections, gene.stop_restorations);

        write_sequence(id_buf, std::string(impl_->fmt_buffer_, desc_len), output_seq);
    }
}

void FastaWriter::write_genes_protein_search(const std::string& sequence_id,
                                             const std::vector<Gene>& genes) {
    // Search-optimized protein output:
    // - Early damage-convertible stops are extended with 'X' instead of truncating
    // - This allows MMseqs2 to align through these positions without penalty
    // - X is interpreted as "unknown amino acid" and gets neutral scoring
    // - Better than correction (which makes false claims) or truncation (which loses signal)
    for (size_t i = 0; i < genes.size(); ++i) {
        const auto& gene = genes[i];

        // Use pre-computed search_protein from frame_selector (has extension logic)
        // Falls back to original protein if search_protein is empty
        const std::string& source_protein = gene.search_protein.empty() ?
            gene.protein : gene.search_protein;

        // Extract longest ORF (before first real stop) for output
        // X characters are preserved, only real stops (*) truncate
        std::string output_protein;
        for (char c : source_protein) {
            if (c == '*') break;  // Real stop terminates ORF
            output_protein += c;   // X characters pass through
        }

        // Skip genes with empty proteins (no valid ORF found)
        if (output_protein.empty()) continue;

        impl_->gene_counter_++;

        // Use read name + strand + gene index as ID for MMseqs2 compatibility
        char id_buf[1024];
        snprintf(id_buf, sizeof(id_buf), "%s_%c_%zu",
            sequence_id.c_str(),
            gene.is_forward ? '+' : '-',
            i + 1);

        // Count X characters for header metadata
        size_t x_count = 0;
        for (char c : output_protein) {
            if (c == 'X') x_count++;
        }

        // Calculate confidence score from protein length and coding probability
        float confidence = calculate_confidence(output_protein.length(), gene.coding_prob);

        int desc_len;
        desc_len = snprintf(impl_->fmt_buffer_, sizeof(impl_->fmt_buffer_),
            "%u..%u length=%zuaa conf=%.3f x_count=%zu damage_pct=%g p_damaged=%g",
            gene.start, gene.end,
            output_protein.length(), confidence, x_count,
            gene.damage_score, gene.p_read_damaged);

        write_sequence(id_buf, std::string(impl_->fmt_buffer_, desc_len), output_protein);
    }
}

void FastaWriter::close() {
    if (impl_->gz_) {
        impl_->gz_->close();
        impl_->gz_.reset();
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
    if (auto p = std::filesystem::path(filename).parent_path(); !p.empty())
        std::filesystem::create_directories(p);
    impl_->file_.open(filename);
    if (!impl_->file_) {
        throw std::runtime_error("Failed to open stats file: " + filename);
    }
    // Write header
    impl_->file_ << "# DART Statistics Output\n";
    impl_->file_ << "# Generated by DART\n\n";
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

} // namespace dart
