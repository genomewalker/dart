#pragma once
// Abstract gzip line reader — virtual dispatch firewall for rapidgzip.
//
// GzLineReaderRapidgzip is defined in src/sequence_io_backend.cpp, which is
// the only TU that includes rapidgzip headers (ODR firewall, same pattern as
// fqdup's fastq_io_backend.cpp).  All other TUs see only this interface.

#include <memory>
#include <string>
#include <cstddef>

namespace dart {

class GzLineReader {
public:
    static constexpr size_t GZBUF_SIZE = 4 * 1024 * 1024;  // 4 MB

    virtual ~GzLineReader() = default;

    // Read one line (without trailing newline) into `line`. Returns false on EOF.
    virtual bool readline(std::string& line) = 0;

    // Zero-copy fast path: return pointer into internal buffer, set len to
    // line length (excluding newline).  Pointer is valid until the next call.
    // Returns nullptr on EOF.
    virtual const char* readline_fast(size_t& len) = 0;

    // Raw chunk read — returns bytes copied into buf (0 = EOF).
    // Drains any internal line-buffer first, then reads directly from the
    // decompressor.  Use this for bulk callers (e.g., TSV chunk parser)
    // that do NOT interleave readline() calls.
    virtual size_t read_raw(char* buf, size_t n) = 0;
};

// Factory — implemented in src/sequence_io_backend.cpp.
// Returns nullptr when HAVE_RAPIDGZIP is not defined; caller falls back to zlib.
std::unique_ptr<GzLineReader> make_gz_reader(const std::string& path);

}  // namespace dart
