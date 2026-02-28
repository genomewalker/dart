// sequence_io_backend.cpp
// Implements make_gz_reader() — the only translation unit that includes
// rapidgzip headers, isolating its ODR-visible symbols from all other TUs.

#ifdef HAVE_RAPIDGZIP
#include <rapidgzip/rapidgzip.hpp>
#include <filereader/Standard.hpp>
#endif

#include "dart/gz_reader_base.hpp"
#include <memory>
#include <vector>
#include <string>
#include <cstdio>
#include <cstring>

namespace {

#ifdef HAVE_RAPIDGZIP
class GzLineReaderRapidgzip : public dart::GzLineReader {
public:
    explicit GzLineReaderRapidgzip(const std::string& path)
        : buf_pos_(0), buf_used_(0), eof_(false) {
        buf_.resize(GZBUF_SIZE);
        overflow_.reserve(256);
        reader_ = std::make_unique<rapidgzip::ParallelGzipReader<>>(
            std::make_unique<rapidgzip::StandardFileReader>(path),
            0,          // threads: 0 = auto-detect from hardware_concurrency
            GZBUF_SIZE
        );
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
            const size_t got = reader_->read(buf + written, n - written);
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

private:
    bool refill() {
        const size_t n = reader_->read(buf_.data(), buf_.size());
        if (n == 0) { eof_ = true; return false; }
        buf_used_ = n;
        buf_pos_  = 0;
        return true;
    }

    std::unique_ptr<rapidgzip::ParallelGzipReader<>> reader_;
    std::vector<char> buf_;
    std::string       overflow_;
    size_t            buf_pos_;
    size_t            buf_used_;
    bool              eof_;
};
#endif  // HAVE_RAPIDGZIP

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
#ifdef HAVE_RAPIDGZIP
    if (is_gzip(path))
        return std::make_unique<GzLineReaderRapidgzip>(path);
#else
    (void)path;
#endif
    return nullptr;  // caller falls back to zlib
}

}  // namespace dart
