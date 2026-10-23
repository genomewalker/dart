// sequence_io_backend.cpp
// Implements make_gz_reader(): single-threaded ISA-L igzip.
//
// rapidgzip was removed 2026-10-21, as fqdup removed it on 2026-07-16: its parallel decode
// raced on large gz over NFS and silently truncated a 30 GB fqdup input (115,146,482 reads
// lost, no error), and its chunk cache added ~0.7 GB to dart profile on a 2M-read library.

#include "dart/gz_reader_base.hpp"
#include <isa-l/igzip_lib.h>
#include <memory>
#include <vector>
#include <string>
#include <cstdio>
#include <cstring>
#include <stdexcept>

namespace {

// Line splitting over a decompressor; subclasses supply decompress().
class BufferedGzLineReader : public dart::GzLineReader {
public:
    BufferedGzLineReader() : buf_pos_(0), buf_used_(0), eof_(false) {
        buf_.resize(GZBUF_SIZE);
        overflow_.reserve(256);
    }

    bool readline(std::string& line) override {
        line.clear();
        while (true) {
            for (size_t i = buf_pos_; i < buf_used_; ++i) {
                if (buf_[i] == '\n') {
                    line.append(buf_.data() + buf_pos_, i - buf_pos_);
                    buf_pos_ = i + 1;
                    return true;
                }
            }
            if (buf_pos_ < buf_used_)
                line.append(buf_.data() + buf_pos_, buf_used_ - buf_pos_);
            buf_pos_ = buf_used_ = 0;
            if (eof_) return !line.empty();
            if (!refill()) return !line.empty();
        }
    }

    size_t read_raw(char* buf, size_t n) override {
        size_t written = 0;
        // Drain internal line-buffer first.
        if (buf_pos_ < buf_used_) {
            const size_t avail = buf_used_ - buf_pos_;
            const size_t take  = std::min(avail, n);
            std::memcpy(buf, buf_.data() + buf_pos_, take);
            buf_pos_ += take;
            written  += take;
            if (written == n) return written;
        }
        // Read directly from the decompressor for the remainder.
        if (!eof_) {
            const size_t got = decompress(buf + written, n - written);
            if (got == 0) eof_ = true;
            written += got;
        }
        return written;
    }

    const char* readline_fast(size_t& len) override {
        overflow_.clear();
        // Fast path: newline found within current buffer — return direct pointer.
        for (size_t i = buf_pos_; i < buf_used_; ++i) {
            if (buf_[i] == '\n') {
                const char* ptr = buf_.data() + buf_pos_;
                len = i - buf_pos_;
                buf_pos_ = i + 1;
                return ptr;
            }
        }
        // Cross-buffer: accumulate remainder then scan into the next buffer.
        if (buf_pos_ < buf_used_)
            overflow_.append(buf_.data() + buf_pos_, buf_used_ - buf_pos_);
        buf_pos_ = buf_used_ = 0;
        while (!eof_) {
            if (!refill()) break;
            for (size_t i = 0; i < buf_used_; ++i) {
                if (buf_[i] == '\n') {
                    overflow_.append(buf_.data(), i);
                    buf_pos_ = i + 1;
                    len = overflow_.size();
                    return len > 0 ? overflow_.c_str() : nullptr;
                }
            }
            overflow_.append(buf_.data(), buf_used_);
            buf_pos_ = buf_used_ = 0;
        }
        len = overflow_.size();
        return len > 0 ? overflow_.c_str() : nullptr;
    }

protected:
    // Fill up to n decompressed bytes; 0 means end of stream.
    virtual size_t decompress(char* out, size_t n) = 0;

private:
    bool refill() {
        const size_t n = decompress(buf_.data(), buf_.size());
        if (n == 0) { eof_ = true; return false; }
        buf_used_ = n;
        buf_pos_  = 0;
        return true;
    }

    std::vector<char> buf_;
    std::string       overflow_;
    size_t            buf_pos_;
    size_t            buf_used_;
    bool              eof_;
};

// Single-threaded ISA-L inflate, as fqdup reads gz (fqdup src/fastq_io_backend.cpp). Two 4 MB
// buffers; handles multi-member files and treats a stream that stops mid-member as an error,
// not an EOF.
class GzLineReaderIgzip : public BufferedGzLineReader {
public:
    explicit GzLineReaderIgzip(const std::string& path)
        : path_(path), in_(GZBUF_SIZE) {
        fp_ = std::fopen(path.c_str(), "rb");
        if (!fp_) throw std::runtime_error("Cannot open: " + path);
        isal_inflate_init(&state_);
        state_.crc_flag = ISAL_GZIP;
    }
    ~GzLineReaderIgzip() override { std::fclose(fp_); }

protected:
    size_t decompress(char* out, size_t n) override {
        while (true) {
            if (state_.avail_in == 0 && !in_eof_) {
                const size_t got = std::fread(in_.data(), 1, in_.size(), fp_);
                if (got == 0) {
                    if (std::ferror(fp_))
                        throw std::runtime_error("Read error on '" + path_ + "'");
                    in_eof_ = true;
                } else {
                    state_.next_in  = in_.data();
                    state_.avail_in = static_cast<uint32_t>(got);
                }
            }
            if (state_.avail_in == 0 && in_eof_) {
                if (state_.block_state != ISAL_BLOCK_FINISH)
                    throw std::runtime_error("Truncated gzip input '" + path_ +
                                             "': stream ends mid-member");
                return 0;
            }
            // Previous member ended and bytes remain: start the next one, keeping the cursor.
            if (state_.block_state == ISAL_BLOCK_FINISH) {
                uint8_t*       ni = state_.next_in;
                const uint32_t ai = state_.avail_in;
                isal_inflate_reset(&state_);
                state_.crc_flag = ISAL_GZIP;
                state_.next_in  = ni;
                state_.avail_in = ai;
            }
            state_.next_out  = reinterpret_cast<uint8_t*>(out);
            state_.avail_out = static_cast<uint32_t>(n);
            if (isal_inflate(&state_) < 0)
                throw std::runtime_error("ISA-L decompression error in '" + path_ + "'");
            const size_t produced = n - state_.avail_out;
            if (produced > 0) return produced;
            // Zero bytes out at a member boundary or with input pending: keep going.
        }
    }

private:
    std::string          path_;
    std::FILE*           fp_ = nullptr;
    std::vector<uint8_t> in_;
    inflate_state        state_{};
    bool                 in_eof_ = false;
};

static bool is_gzip(const std::string& path) {
    FILE* f = fopen(path.c_str(), "rb");
    if (!f) return false;
    unsigned char magic[2] = {0, 0};
    bool gz = (fread(magic, 1, 2, f) == 2) &&
              (magic[0] == 0x1f && magic[1] == 0x8b);
    fclose(f);
    return gz;
}

}  // anonymous namespace

namespace dart {

std::unique_ptr<GzLineReader> make_gz_reader(const std::string& path) {
    if (is_gzip(path))
        return std::make_unique<GzLineReaderIgzip>(path);
    return nullptr;  // caller falls back to zlib
}

}  // namespace dart
