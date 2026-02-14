#pragma once
/**
 * @file fast_tsv.hpp
 * @brief Parallel TSV parser using mmap + chunked processing
 *
 * Key optimizations:
 * - mmap entire file (no read syscalls)
 * - memchr for tab/newline finding (SIMD on glibc)
 * - Parallel chunk processing with OpenMP
 * - Zero-copy field extraction (string_view)
 */

#include <atomic>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <functional>
#include <string_view>
#include <vector>

#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace agp {

/**
 * @brief Memory-mapped TSV file for zero-copy parsing
 */
class MappedTSV {
public:
    explicit MappedTSV(const char* path) {
        fd_ = open(path, O_RDONLY);
        if (fd_ < 0) return;

        struct stat st;
        if (fstat(fd_, &st) < 0) {
            close(fd_);
            fd_ = -1;
            return;
        }
        size_ = static_cast<size_t>(st.st_size);

        data_ = static_cast<const char*>(
            mmap(nullptr, size_, PROT_READ, MAP_PRIVATE | MAP_POPULATE, fd_, 0));
        if (data_ == MAP_FAILED) {
            close(fd_);
            fd_ = -1;
            data_ = nullptr;
            return;
        }

        // Advise sequential access
        madvise(const_cast<char*>(data_), size_, MADV_SEQUENTIAL);
    }

    ~MappedTSV() {
        if (data_) munmap(const_cast<char*>(data_), size_);
        if (fd_ >= 0) close(fd_);
    }

    MappedTSV(const MappedTSV&) = delete;
    MappedTSV& operator=(const MappedTSV&) = delete;

    bool valid() const { return data_ != nullptr; }
    const char* data() const { return data_; }
    size_t size() const { return size_; }

private:
    const char* data_ = nullptr;
    size_t size_ = 0;
    int fd_ = -1;
};

/**
 * @brief Fast field extraction using memchr (SIMD-optimized in glibc)
 */
struct TSVLine {
    std::string_view fields[16];
    uint8_t n_fields = 0;

    // Parse a line into fields (up to 16)
    void parse(const char* start, const char* end) {
        n_fields = 0;
        const char* p = start;
        while (p < end && n_fields < 16) {
            const char* tab = static_cast<const char*>(
                memchr(p, '\t', static_cast<size_t>(end - p)));
            if (!tab) tab = end;
            fields[n_fields++] = {p, static_cast<size_t>(tab - p)};
            p = tab + 1;
        }
    }
};

/**
 * @brief Parallel TSV processing with row callback (legacy - finds line boundaries first)
 *
 * @param tsv Memory-mapped TSV file
 * @param callback Called for each row (must be thread-safe)
 * @param skip_header Skip first line
 * @param min_fields Minimum fields required (skip malformed rows)
 * @return Number of rows processed
 */
template<typename Callback>
size_t parallel_parse_tsv(
    const MappedTSV& tsv,
    Callback&& callback,
    bool skip_header = false,
    uint8_t min_fields = 1)
{
    if (!tsv.valid()) return 0;

    const char* data = tsv.data();
    const size_t size = tsv.size();

    // Find all line boundaries
    std::vector<size_t> line_starts;
    line_starts.reserve(size / 100);  // Estimate ~100 bytes per line
    line_starts.push_back(0);

    for (size_t i = 0; i < size; ++i) {
        if (data[i] == '\n' && i + 1 < size) {
            line_starts.push_back(i + 1);
        }
    }

    size_t start_idx = skip_header ? 1 : 0;
    size_t n_lines = line_starts.size();

    std::atomic<size_t> count{0};

    #pragma omp parallel
    {
        TSVLine line;

        #pragma omp for schedule(dynamic, 4096)
        for (size_t i = start_idx; i < n_lines; ++i) {
            size_t line_start = line_starts[i];
            size_t line_end = (i + 1 < n_lines) ? line_starts[i + 1] - 1 : size;

            // Skip empty lines
            if (line_end <= line_start) continue;

            // Parse fields
            line.parse(data + line_start, data + line_end);

            if (line.n_fields >= min_fields) {
                callback(line);
                count.fetch_add(1, std::memory_order_relaxed);
            }
        }
    }

    return count.load();
}

/**
 * @brief Partitioned parallel TSV parsing
 *
 * Key optimization: Each thread independently finds and processes lines
 * within its partition - no sequential line boundary scan.
 *
 * Algorithm:
 * 1. Divide file into N partitions (one per thread)
 * 2. Each thread scans to find first complete line in its partition
 * 3. Each thread processes lines until reaching the next partition
 * 4. No global synchronization during parsing
 *
 * @param tsv Memory-mapped TSV file
 * @param callback Called for each row (must be thread-safe)
 * @param skip_header Skip first line
 * @param min_fields Minimum fields required (skip malformed rows)
 * @return Number of rows processed
 */
template<typename Callback>
size_t partitioned_parse_tsv(
    const MappedTSV& tsv,
    Callback&& callback,
    bool skip_header = false,
    uint8_t min_fields = 1)
{
    if (!tsv.valid()) return 0;

    const char* data = tsv.data();
    const size_t size = tsv.size();
    if (size == 0) return 0;

    // Determine number of partitions
    #ifdef _OPENMP
    const int n_threads = omp_get_max_threads();
    #else
    const int n_threads = 1;
    #endif

    // Minimum partition size (avoid too many small partitions)
    constexpr size_t MIN_PARTITION_SIZE = 1024 * 1024;  // 1MB
    const int n_partitions = std::min(n_threads, std::max(1, static_cast<int>(size / MIN_PARTITION_SIZE)));
    const size_t partition_size = size / static_cast<size_t>(n_partitions);

    std::atomic<size_t> total_count{0};

    #pragma omp parallel num_threads(n_partitions)
    {
        #ifdef _OPENMP
        const int tid = omp_get_thread_num();
        #else
        const int tid = 0;
        #endif

        // Compute partition boundaries
        size_t part_start = static_cast<size_t>(tid) * partition_size;
        size_t part_end = (tid == n_partitions - 1) ? size : (static_cast<size_t>(tid) + 1) * partition_size;

        // Find first complete line in partition
        // Thread 0 starts at 0 (or after header), others scan for newline
        if (tid > 0) {
            // Find the first newline after part_start
            const char* nl = static_cast<const char*>(
                memchr(data + part_start, '\n', part_end - part_start));
            if (nl) {
                part_start = static_cast<size_t>(nl - data) + 1;
            } else {
                // No newline in this partition - skip it
                part_start = part_end;
            }
        } else if (skip_header) {
            // Thread 0: skip header line
            const char* nl = static_cast<const char*>(
                memchr(data, '\n', size));
            if (nl) {
                part_start = static_cast<size_t>(nl - data) + 1;
            }
        }

        // Process lines in this partition
        TSVLine line;
        size_t local_count = 0;
        const char* p = data + part_start;
        const char* end = data + part_end;

        while (p < end) {
            // Find end of line
            const char* line_end = static_cast<const char*>(
                memchr(p, '\n', static_cast<size_t>(end - p)));

            // If no newline, this is the last line (only process if we're the last partition)
            if (!line_end) {
                if (tid == n_partitions - 1 && p < data + size) {
                    line_end = data + size;
                } else {
                    break;
                }
            }

            // Parse and process line
            if (line_end > p) {
                // Handle \r\n
                const char* actual_end = (line_end > p && *(line_end - 1) == '\r')
                    ? line_end - 1 : line_end;

                line.parse(p, actual_end);

                if (line.n_fields >= min_fields) {
                    callback(line);
                    ++local_count;
                }
            }

            p = line_end + 1;
        }

        total_count.fetch_add(local_count, std::memory_order_relaxed);
    }

    return total_count.load();
}

/**
 * @brief Fast string to float (no locale, no error checking)
 */
inline float fast_stof(std::string_view s) {
    // Simple fast path for common cases
    float result = 0.0f;
    float sign = 1.0f;
    size_t i = 0;

    if (s.empty()) return 0.0f;
    if (s[0] == '-') { sign = -1.0f; i = 1; }
    else if (s[0] == '+') { i = 1; }

    // Integer part
    while (i < s.size() && s[i] >= '0' && s[i] <= '9') {
        result = result * 10.0f + static_cast<float>(s[i] - '0');
        i++;
    }

    // Fractional part
    if (i < s.size() && s[i] == '.') {
        i++;
        float frac = 0.0f;
        float div = 1.0f;
        while (i < s.size() && s[i] >= '0' && s[i] <= '9') {
            frac = frac * 10.0f + static_cast<float>(s[i] - '0');
            div *= 10.0f;
            i++;
        }
        result += frac / div;
    }

    // Exponent (E notation)
    if (i < s.size() && (s[i] == 'e' || s[i] == 'E')) {
        i++;
        float exp_sign = 1.0f;
        if (i < s.size() && s[i] == '-') { exp_sign = -1.0f; i++; }
        else if (i < s.size() && s[i] == '+') { i++; }

        int exp = 0;
        while (i < s.size() && s[i] >= '0' && s[i] <= '9') {
            exp = exp * 10 + (s[i] - '0');
            i++;
        }
        result *= std::pow(10.0f, exp_sign * static_cast<float>(exp));
    }

    return sign * result;
}

/**
 * @brief Fast string to uint32 (no locale, no error checking)
 */
inline uint32_t fast_stou(std::string_view s) {
    uint32_t result = 0;
    for (char c : s) {
        if (c >= '0' && c <= '9') {
            result = result * 10 + static_cast<uint32_t>(c - '0');
        }
    }
    return result;
}

} // namespace agp
