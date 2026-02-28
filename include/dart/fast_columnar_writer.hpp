#pragma once
/**
 * @file fast_columnar_writer.hpp
 * @brief High-throughput TSV to columnar conversion (DuckDB-style)
 *
 * Architecture (inspired by DuckDB):
 * 1. Chunked file reading (64MB chunks, no full mmap)
 * 2. Morsel-driven parallelism (10K row batches)
 * 3. Spilling dictionary (memory-bounded with disk overflow)
 * 4. Streaming row-group emission (bounded memory)
 * 5. SIMD delimiter scanning
 * 6. Single-pass: dictionary build + columnar write inline
 *
 * Memory bounds:
 * - Input buffer: 128 MB (double-buffered chunks)
 * - Dictionary: configurable cap (default 4 GB, spills to disk)
 * - Row-group buffers: ~50 MB (fixed)
 * - Total: ~4-5 GB regardless of file size
 */

#include "dart/columnar_index.hpp"
#include "dart/fast_tsv.hpp"
#include "dart/damage_index_reader.hpp"
#include "dart/damage_alignment_evidence.hpp"
#include "dart/log_utils.hpp"

#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstring>
#include <numeric>
#include <condition_variable>
#include <deque>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <exception>
#include <limits>
#include <memory>
#include <mutex>
#include <queue>
#include <string>
#include <thread>
#include <unordered_map>
#include <vector>

#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#include <zlib.h>
#include "dart/gz_reader_base.hpp"

#ifdef HAVE_ZSTD
#include <zstd.h>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef USE_AVX2
#include <immintrin.h>
#endif

namespace dart {

// ============================================================================
// Constants for streaming architecture
// ============================================================================
constexpr size_t CHUNK_SIZE_BYTES = 64 * 1024 * 1024;    // 64 MB input chunks
constexpr size_t MORSEL_SIZE = 8192;                      // 8K rows per morsel
constexpr size_t DICT_MEMORY_CAP = 4ULL * 1024 * 1024 * 1024;  // 4 GB before spill
constexpr size_t SPILL_BATCH_SIZE = 1024 * 1024;          // 1M entries per spill file

#ifdef USE_AVX2
/**
 * @brief SIMD-accelerated tab finder using AVX2
 * Finds first '\t' in buffer, returns pointer or end if not found
 */
inline const char* simd_find_tab(const char* p, const char* end) {
    const __m256i tab_vec = _mm256_set1_epi8('\t');

    // Process 32 bytes at a time
    while (p + 32 <= end) {
        __m256i chunk = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(p));
        __m256i cmp = _mm256_cmpeq_epi8(chunk, tab_vec);
        int mask = _mm256_movemask_epi8(cmp);
        if (mask != 0) {
            return p + __builtin_ctz(static_cast<unsigned>(mask));
        }
        p += 32;
    }

    // Scalar fallback for remainder
    while (p < end && *p != '\t') ++p;
    return p;
}

/**
 * @brief SIMD-accelerated newline finder using AVX2
 */
inline const char* simd_find_newline(const char* p, const char* end) {
    const __m256i nl_vec = _mm256_set1_epi8('\n');

    while (p + 32 <= end) {
        __m256i chunk = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(p));
        __m256i cmp = _mm256_cmpeq_epi8(chunk, nl_vec);
        int mask = _mm256_movemask_epi8(cmp);
        if (mask != 0) {
            return p + __builtin_ctz(static_cast<unsigned>(mask));
        }
        p += 32;
    }

    while (p < end && *p != '\n') ++p;
    return p;
}
#else
// Fallback to standard memchr
inline const char* simd_find_tab(const char* p, const char* end) {
    const char* result = static_cast<const char*>(memchr(p, '\t', static_cast<size_t>(end - p)));
    return result ? result : end;
}
inline const char* simd_find_newline(const char* p, const char* end) {
    const char* result = static_cast<const char*>(memchr(p, '\n', static_cast<size_t>(end - p)));
    return result ? result : end;
}
#endif

// ============================================================================
// Chunked File Reader - reads file in bounded chunks, handles line boundaries
// ============================================================================

/**
 * @brief Streaming file reader with bounded memory
 *
 * Reads file in 64MB chunks, handles partial lines at chunk boundaries.
 * Supports both plain and gzip files.
 */
class ChunkedFileReader {
public:
    explicit ChunkedFileReader(const std::string& path, size_t chunk_size = CHUNK_SIZE_BYTES)
        : chunk_size_(chunk_size) {
        is_gzip_ = (path.size() >= 3 && path.substr(path.size() - 3) == ".gz");

        if (is_gzip_) {
            gz_ = gzopen(path.c_str(), "rb");
            if (!gz_) throw std::runtime_error("Cannot open gzip file: " + path);
            gzbuffer(gz_, 1024 * 1024);  // 1MB internal buffer
        } else {
            fd_ = open(path.c_str(), O_RDONLY);
            if (fd_ < 0) throw std::runtime_error("Cannot open file: " + path);

            struct stat st;
            if (fstat(fd_, &st) == 0) file_size_ = static_cast<size_t>(st.st_size);

#ifdef __linux__
            posix_fadvise(fd_, 0, 0, POSIX_FADV_SEQUENTIAL);
#endif
        }

        // Allocate double buffer for overlapped I/O
        buffer_.resize(chunk_size_ * 2);
        read_pos_ = 0;
        data_end_ = 0;
    }

    ~ChunkedFileReader() {
        if (gz_) gzclose(gz_);
        if (fd_ >= 0) close(fd_);
    }

    bool valid() const { return gz_ != nullptr || fd_ >= 0; }
    size_t file_size() const { return file_size_; }
    uint64_t bytes_read() const { return total_bytes_read_; }

    /**
     * @brief Get next complete line
     * @return string_view to line (without newline), empty if EOF
     *
     * Data is valid until next call to next_line() or next_chunk()
     */
    std::string_view next_line() {
        while (true) {
            // Search for newline in current buffer
            const char* start = buffer_.data() + read_pos_;
            const char* end = buffer_.data() + data_end_;
            const char* nl = simd_find_newline(start, end);

            if (nl < end) {
                // Found newline
                std::string_view line(start, static_cast<size_t>(nl - start));
                read_pos_ = static_cast<size_t>(nl - buffer_.data()) + 1;
                return line;
            }

            // No newline found - need more data
            if (eof_) {
                // Return remaining data as final line
                if (start < end) {
                    std::string_view line(start, static_cast<size_t>(end - start));
                    read_pos_ = data_end_;
                    return line;
                }
                return {};  // True EOF
            }

            // Move partial line to beginning and read more
            size_t partial_len = static_cast<size_t>(end - start);
            if (partial_len > 0 && start != buffer_.data()) {
                std::memmove(buffer_.data(), start, partial_len);
            }
            read_pos_ = 0;
            data_end_ = partial_len;

            // Read more data
            size_t to_read = chunk_size_;
            if (data_end_ + to_read > buffer_.size()) {
                to_read = buffer_.size() - data_end_;
            }

            ssize_t bytes = read_raw(buffer_.data() + data_end_, to_read);
            if (bytes <= 0) {
                eof_ = true;
            } else {
                data_end_ += static_cast<size_t>(bytes);
                total_bytes_read_ += static_cast<size_t>(bytes);
            }
        }
    }

    /**
     * @brief Skip header line if present
     */
    void skip_header_if_present() {
        // Peek at first bytes
        if (data_end_ == 0) {
            ssize_t bytes = read_raw(buffer_.data(), std::min(chunk_size_, buffer_.size()));
            if (bytes > 0) {
                data_end_ = static_cast<size_t>(bytes);
                total_bytes_read_ = static_cast<size_t>(bytes);
            }
        }

        if (data_end_ >= 5 && std::strncmp(buffer_.data(), "query", 5) == 0) {
            next_line();  // Consume header
        }
    }

    /**
     * @brief Reset to beginning of file (for second pass)
     */
    void reset() {
        if (gz_) {
            gzrewind(gz_);
        } else if (fd_ >= 0) {
            lseek(fd_, 0, SEEK_SET);
        }
        read_pos_ = 0;
        data_end_ = 0;
        eof_ = false;
        total_bytes_read_ = 0;
    }

private:
    ssize_t read_raw(char* buf, size_t len) {
        if (gz_) {
            return gzread(gz_, buf, static_cast<unsigned>(len));
        } else {
            return ::read(fd_, buf, len);
        }
    }

    gzFile gz_ = nullptr;
    int fd_ = -1;
    bool is_gzip_ = false;
    size_t file_size_ = 0;
    size_t chunk_size_;

    std::vector<char> buffer_;
    size_t read_pos_ = 0;
    size_t data_end_ = 0;
    bool eof_ = false;
    uint64_t total_bytes_read_ = 0;
};

// ============================================================================
// Morsel - a batch of parsed rows for parallel processing
// ============================================================================

struct Morsel {
    static constexpr size_t CAPACITY = MORSEL_SIZE;

    // Raw line data (pointers into reader buffer - only valid during parsing)
    std::vector<std::string_view> lines;

    // Parsed data (stable, owned by morsel)
    std::vector<std::string> queries;    // Copied strings (stable storage)
    std::vector<std::string> targets;
    std::vector<float> bit_scores;
    std::vector<float> evalues;
    std::vector<float> fidents;
    std::vector<uint32_t> alnlens;
    std::vector<uint32_t> tstarts;
    std::vector<uint32_t> tends;
    std::vector<uint32_t> qlens;
    std::vector<uint32_t> tlens;
    std::vector<std::string> qalns;
    std::vector<std::string> talns;

    size_t size() const { return queries.size(); }

    void clear() {
        lines.clear();
        queries.clear();
        targets.clear();
        bit_scores.clear();
        evalues.clear();
        fidents.clear();
        alnlens.clear();
        tstarts.clear();
        tends.clear();
        qlens.clear();
        tlens.clear();
        qalns.clear();
        talns.clear();
    }

    void reserve(size_t n) {
        lines.reserve(n);
        queries.reserve(n);
        targets.reserve(n);
        bit_scores.reserve(n);
        evalues.reserve(n);
        fidents.reserve(n);
        alnlens.reserve(n);
        tstarts.reserve(n);
        tends.reserve(n);
        qlens.reserve(n);
        tlens.reserve(n);
        qalns.reserve(n);
        talns.reserve(n);
    }
};

// ============================================================================
// Spilling String Dictionary - memory-bounded with disk overflow
// ============================================================================

/**
 * @brief Dictionary that spills to disk when memory cap exceeded
 *
 * Uses sorted temp files for spilled entries, merged during finalization.
 */
class SpillingStringDict {
public:
    struct TransparentStringHash {
        using is_transparent = void;
        size_t operator()(std::string_view sv) const noexcept {
            return std::hash<std::string_view>{}(sv);
        }
        size_t operator()(const std::string& s) const noexcept {
            return std::hash<std::string_view>{}(s);
        }
        size_t operator()(const char* s) const noexcept {
            return std::hash<std::string_view>{}(s ? std::string_view(s) : std::string_view{});
        }
    };

    struct TransparentStringEq {
        using is_transparent = void;
        bool operator()(std::string_view a, std::string_view b) const noexcept {
            return a == b;
        }
        bool operator()(const std::string& a, const std::string& b) const noexcept {
            return a == b;
        }
    };

    using StringToIdMap = std::unordered_map<std::string, uint32_t, TransparentStringHash, TransparentStringEq>;

    SpillingStringDict(const std::string& name, size_t memory_cap = DICT_MEMORY_CAP / 2)
        : name_(name), memory_cap_(memory_cap) {
        // Create temp directory
        temp_dir_ = std::filesystem::temp_directory_path() / ("agp_dict_" + name_ + "_" + std::to_string(getpid()));
        std::filesystem::create_directories(temp_dir_);
    }

    ~SpillingStringDict() {
        // Cleanup temp files
        std::error_code ec;
        std::filesystem::remove_all(temp_dir_, ec);
    }

    /**
     * @brief Insert string, returns true if new (not seen before)
     */
    bool insert(std::string_view sv) {
        uint64_t h = hash_sv(sv);

        // Check in-memory map first
        auto it = in_memory_.find(sv);
        if (it != in_memory_.end()) return false;

        // Check if we've seen this in spilled data (bloom filter approximation)
        if (maybe_in_spilled(h)) {
            // Could be false positive - will dedupe during merge
        }

        // Add to in-memory
        std::string owned(sv);
        memory_used_ += owned.size() + 32;  // String + overhead
        in_memory_.emplace(std::move(owned), next_id_++);
        ++total_count_;

        // Check if need to spill
        if (memory_used_ > memory_cap_) {
            spill_to_disk();
        }

        return true;
    }

    /**
     * @brief Finalize dictionary - merge all spilled files + in-memory
     * @return Final dictionary with sequential IDs
     */
    std::pair<std::vector<std::string>, StringToIdMap> finalize() {
        std::vector<std::string> names;
        StringToIdMap lookup;

        if (spill_count_ == 0) {
            // No spilling - just convert in-memory map
            names.reserve(in_memory_.size());
            for (auto& [str, id] : in_memory_) {
                uint32_t new_id = static_cast<uint32_t>(names.size());
                names.push_back(std::move(str));
                lookup[names.back()] = new_id;
            }
        } else {
            // Merge spilled files with in-memory
            names = merge_all_spilled();
            for (uint32_t i = 0; i < names.size(); ++i) {
                lookup[names[i]] = i;
            }
        }

        // Clear in-memory to free memory
        in_memory_.clear();
        memory_used_ = 0;

        return {std::move(names), std::move(lookup)};
    }

    size_t size() const { return total_count_; }
    size_t memory_used() const { return memory_used_; }
    size_t spill_count() const { return spill_count_; }

private:
    static uint64_t hash_sv(std::string_view sv) {
        uint64_t h = 14695981039346656037ULL;
        for (unsigned char c : sv) {
            h ^= c;
            h *= 1099511628211ULL;
        }
        return h;
    }

    bool maybe_in_spilled(uint64_t h) {
        // Simple bloom filter check
        size_t idx = h % bloom_.size();
        return bloom_[idx];
    }

    void spill_to_disk() {
        if (in_memory_.empty()) return;

        // Sort entries
        std::vector<std::string> sorted;
        sorted.reserve(in_memory_.size());
        for (auto& [str, id] : in_memory_) {
            sorted.push_back(str);

            // Update bloom filter
            uint64_t h = hash_sv(str);
            bloom_[h % bloom_.size()] = true;
        }
        std::sort(sorted.begin(), sorted.end());

        // Write to temp file
        std::string spill_path = temp_dir_ / ("spill_" + std::to_string(spill_count_++) + ".txt");
        std::ofstream out(spill_path);
        for (const auto& s : sorted) {
            out << s << '\n';
        }
        spill_files_.push_back(spill_path);

        // Clear in-memory
        in_memory_.clear();
        memory_used_ = 0;
    }

    std::vector<std::string> merge_all_spilled() {
        // K-way merge of sorted spill files + in-memory
        std::vector<std::string> result;
        result.reserve(total_count_);

        // Open all spill files
        std::vector<std::ifstream> files;
        std::vector<std::string> current_lines;
        std::vector<bool> eof;

        for (const auto& path : spill_files_) {
            files.emplace_back(path);
            current_lines.emplace_back();
            eof.push_back(false);
            if (!std::getline(files.back(), current_lines.back())) {
                eof.back() = true;
            }
        }

        // Sort in-memory entries
        std::vector<std::string> mem_sorted;
        mem_sorted.reserve(in_memory_.size());
        for (auto& [str, id] : in_memory_) {
            mem_sorted.push_back(str);
        }
        std::sort(mem_sorted.begin(), mem_sorted.end());
        size_t mem_idx = 0;

        // K-way merge
        std::string last_added;
        while (true) {
            // Find minimum across all sources
            std::string* min_ptr = nullptr;
            size_t min_source = SIZE_MAX;

            for (size_t i = 0; i < files.size(); ++i) {
                if (!eof[i] && (!min_ptr || current_lines[i] < *min_ptr)) {
                    min_ptr = &current_lines[i];
                    min_source = i;
                }
            }

            // Check in-memory
            if (mem_idx < mem_sorted.size() && (!min_ptr || mem_sorted[mem_idx] < *min_ptr)) {
                min_ptr = &mem_sorted[mem_idx];
                min_source = SIZE_MAX;
            }

            if (!min_ptr) break;  // All sources exhausted

            // Add if not duplicate
            if (result.empty() || *min_ptr != result.back()) {
                result.push_back(*min_ptr);
            }

            // Advance source
            if (min_source == SIZE_MAX) {
                ++mem_idx;
            } else {
                if (!std::getline(files[min_source], current_lines[min_source])) {
                    eof[min_source] = true;
                }
            }
        }

        return result;
    }

    std::string name_;
    size_t memory_cap_;
    std::filesystem::path temp_dir_;

    StringToIdMap in_memory_;
    size_t memory_used_ = 0;
    uint32_t next_id_ = 0;
    uint64_t total_count_ = 0;

    std::vector<std::string> spill_files_;
    size_t spill_count_ = 0;
    std::vector<bool> bloom_ = std::vector<bool>(1024 * 1024, false);  // 1M bit bloom filter
};

// ============================================================================
// Morsel Queue - thread-safe work queue for parallel processing
// ============================================================================

class MorselQueue {
public:
    void push(Morsel&& m) {
        std::lock_guard<std::mutex> lock(mutex_);
        queue_.push(std::move(m));
        cv_.notify_one();
    }

    bool pop(Morsel& m, bool blocking = true) {
        std::unique_lock<std::mutex> lock(mutex_);
        if (blocking) {
            cv_.wait(lock, [this] { return !queue_.empty() || done_; });
        }
        if (queue_.empty()) return false;
        m = std::move(queue_.front());
        queue_.pop();
        return true;
    }

    void set_done() {
        std::lock_guard<std::mutex> lock(mutex_);
        done_ = true;
        cv_.notify_all();
    }

    bool is_done() const { return done_; }
    size_t size() const { return queue_.size(); }

private:
    std::queue<Morsel> queue_;
    std::mutex mutex_;
    std::condition_variable cv_;
    bool done_ = false;
};

/**
 * @brief FNV-1a hash for string_view (fast, no allocations)
 */
struct StringViewHash {
    size_t operator()(std::string_view sv) const noexcept {
        size_t hash = 14695981039346656037ULL;  // FNV offset basis
        for (char c : sv) {
            hash ^= static_cast<size_t>(c);
            hash *= 1099511628211ULL;  // FNV prime
        }
        return hash;
    }
};

/**
 * @brief String interning pool with stable storage (never moves memory)
 *
 * Uses chunked allocation to avoid invalidating existing string_views.
 * Each chunk is 1MB; strings never span chunks.
 */
class StringPool {
public:
    static constexpr size_t CHUNK_SIZE = 1024 * 1024;  // 1MB chunks

    StringPool() {
        chunks_.emplace_back(std::make_unique<char[]>(CHUNK_SIZE));
    }

    std::string_view intern(std::string_view s) {
        size_t needed = s.size() + 1;  // +1 for null terminator
        if (needed > CHUNK_SIZE) {
            // Oversized string gets its own allocation
            chunks_.emplace_back(std::make_unique<char[]>(needed));
            char* ptr = chunks_.back().get();
            std::memcpy(ptr, s.data(), s.size());
            ptr[s.size()] = '\0';
            total_size_ += needed;
            return {ptr, s.size()};
        }

        // Check if current chunk has room
        if (current_offset_ + needed > CHUNK_SIZE) {
            // Allocate new chunk
            chunks_.emplace_back(std::make_unique<char[]>(CHUNK_SIZE));
            current_offset_ = 0;
        }

        char* ptr = chunks_.back().get() + current_offset_;
        std::memcpy(ptr, s.data(), s.size());
        ptr[s.size()] = '\0';
        current_offset_ += needed;
        total_size_ += needed;
        return {ptr, s.size()};
    }

    size_t size() const { return total_size_; }

private:
    std::vector<std::unique_ptr<char[]>> chunks_;
    size_t current_offset_ = 0;
    size_t total_size_ = 0;
};

// ============================================================================
// Flat open-addressing hash map (string_view -> uint32_t) for parallel dictionaries
// ============================================================================

struct FlatSVU32Map {
    struct Slot {
        uint64_t hash = 0;         // 0 == empty
        std::string_view key{};
        uint32_t value = 0;
    };

    std::vector<Slot> slots_;
    size_t mask_ = 0;
    size_t size_ = 0;

    static inline uint64_t hash_sv(std::string_view sv) {
        uint64_t h = 14695981039346656037ULL; // FNV-1a 64
        for (unsigned char c : sv) {
            h ^= c;
            h *= 1099511628211ULL;
        }
        return (h == 0) ? 1 : h;
    }

    static size_t next_pow2(size_t x) {
        size_t p = 1;
        while (p < x) p <<= 1;
        return p;
    }

    void reserve(size_t expected) {
        size_t cap = next_pow2(std::max<size_t>(16, (expected * 10) / 7 + 1));
        rehash(cap);
    }

    void clear_release() {
        std::vector<Slot>().swap(slots_);
        mask_ = 0;
        size_ = 0;
    }

    bool try_emplace(std::string_view key, uint32_t value, uint64_t h) {
        if (slots_.empty()) reserve(1024);
        if ((size_ + 1) * 10 >= slots_.size() * 7) rehash(slots_.size() * 2);

        size_t i = static_cast<size_t>(h) & mask_;
        while (true) {
            Slot& s = slots_[i];
            if (s.hash == 0) {
                s.hash = h;
                s.key = key;
                s.value = value;
                ++size_;
                return true;
            }
            if (s.hash == h && s.key == key) return false;
            i = (i + 1) & mask_;
        }
    }

    uint32_t find_or_throw(std::string_view key, uint64_t h) const {
        if (slots_.empty()) throw std::runtime_error("Key missing in empty dictionary");
        size_t i = static_cast<size_t>(h) & mask_;
        while (true) {
            const Slot& s = slots_[i];
            if (s.hash == 0) break;
            if (s.hash == h && s.key == key) return s.value;
            i = (i + 1) & mask_;
        }
        throw std::runtime_error("Key missing in dictionary");
    }

private:
    void rehash(size_t new_cap) {
        new_cap = next_pow2(new_cap);
        std::vector<Slot> old;
        old.swap(slots_);
        slots_.assign(new_cap, Slot{});
        mask_ = new_cap - 1;
        size_ = 0;

        for (const Slot& s : old) {
            if (s.hash == 0) continue;
            size_t i = static_cast<size_t>(s.hash) & mask_;
            while (slots_[i].hash != 0) i = (i + 1) & mask_;
            slots_[i] = s;
            ++size_;
        }
    }
};

struct DictResult {
    std::vector<std::string_view> names;   // id -> string_view
    std::vector<FlatSVU32Map> shard_maps;  // lookup by shard
    size_t shard_mask = 0;                 // shards must be power-of-two
};

inline bool parse_first_two_fields(
    const char* p, const char* line_end,
    std::string_view& q, std::string_view& t)
{
    const char* tab1 = simd_find_tab(p, line_end);
    if (tab1 >= line_end) return false;
    const char* tab2 = simd_find_tab(tab1 + 1, line_end);
    if (tab2 > line_end) return false;
    q = std::string_view(p, static_cast<size_t>(tab1 - p));
    t = std::string_view(tab1 + 1, static_cast<size_t>(tab2 - tab1 - 1));
    return true;
}

inline DictResult build_global_dict_from_locals(
    std::vector<std::vector<std::string_view>>& local_uniques,
    size_t num_shards_pow2)
{
    DictResult out;
    out.shard_maps.resize(num_shards_pow2);
    out.shard_mask = num_shards_pow2 - 1;

    // 1) Count per-shard keys
    std::vector<size_t> shard_counts(num_shards_pow2, 0);
    for (const auto& vec : local_uniques) {
        for (std::string_view sv : vec) {
            uint64_t h = FlatSVU32Map::hash_sv(sv);
            ++shard_counts[static_cast<size_t>(h) & out.shard_mask];
        }
    }

    // 2) Pre-allocate shard vectors
    std::vector<std::vector<std::string_view>> shard_keys(num_shards_pow2);
    for (size_t s = 0; s < num_shards_pow2; ++s) {
        shard_keys[s].reserve(shard_counts[s]);
    }

    // 3) Distribute keys to shards
    for (const auto& vec : local_uniques) {
        for (std::string_view sv : vec) {
            uint64_t h = FlatSVU32Map::hash_sv(sv);
            shard_keys[static_cast<size_t>(h) & out.shard_mask].push_back(sv);
        }
    }

    // Release thread-local vectors early
    for (auto& v : local_uniques) {
        std::vector<std::string_view>().swap(v);
    }

    // 4) Per-shard sort+unique in parallel
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1)
#endif
    for (size_t s = 0; s < num_shards_pow2; ++s) {
        auto& v = shard_keys[s];
        std::sort(v.begin(), v.end());
        v.erase(std::unique(v.begin(), v.end()), v.end());
    }

    // 5) Prefix-sum shard sizes -> global ID ranges
    std::vector<uint32_t> shard_base(num_shards_pow2 + 1, 0);
    for (size_t s = 0; s < num_shards_pow2; ++s) {
        shard_base[s + 1] = shard_base[s] + static_cast<uint32_t>(shard_keys[s].size());
    }

    out.names.resize(shard_base.back());

    // 6) Build read-only per-shard lookup maps in parallel
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1)
#endif
    for (size_t s = 0; s < num_shards_pow2; ++s) {
        auto& keys = shard_keys[s];
        auto& map = out.shard_maps[s];
        map.reserve(keys.size());

        uint32_t id = shard_base[s];
        for (std::string_view sv : keys) {
            uint64_t h = FlatSVU32Map::hash_sv(sv);
            map.try_emplace(sv, id, h);
            out.names[id] = sv;
            ++id;
        }

        // Release shard vector after building map
        std::vector<std::string_view>().swap(keys);
    }

    return out;
}

inline uint32_t dict_lookup(const DictResult& d, std::string_view sv) {
    uint64_t h = FlatSVU32Map::hash_sv(sv);
    size_t shard = static_cast<size_t>(h) & d.shard_mask;
    return d.shard_maps[shard].find_or_throw(sv, h);
}

/**
 * @brief ZSTD compression utilities
 */
#ifdef HAVE_ZSTD
inline std::vector<char> compress_zstd(const void* src, size_t src_size, int level = 3) {
    size_t max_dst = ZSTD_compressBound(src_size);
    std::vector<char> dst(max_dst);
    size_t compressed = ZSTD_compress(dst.data(), max_dst, src, src_size, level);
    if (ZSTD_isError(compressed)) {
        throw std::runtime_error("ZSTD compression failed");
    }
    dst.resize(compressed);
    return dst;
}

inline std::vector<char> decompress_zstd(const void* src, size_t src_size, size_t dst_size) {
    std::vector<char> dst(dst_size);
    size_t decompressed = ZSTD_decompress(dst.data(), dst_size, src, src_size);
    if (ZSTD_isError(decompressed) || decompressed != dst_size) {
        throw std::runtime_error("ZSTD decompression failed");
    }
    return dst;
}
#endif

}  // namespace dart (compression utilities)

namespace dart {

/**
 * @brief Streaming gzip reader using ring buffer
 *
 * Prefers rapidgzip (parallel, in-process) when available; falls back to zlib.
 */
class GzipStreamReader {
public:
    static constexpr size_t CHUNK_SIZE = 16 * 1024 * 1024;  // 16MB chunks

    explicit GzipStreamReader(const std::string& path) {
        buffer_.resize(CHUNK_SIZE);
        // Try rapidgzip first (faster, parallel, no PATH dependency).
        rg_reader_ = dart::make_gz_reader(path);
        if (rg_reader_) return;

        // Fall back to zlib.
        gz_ = gzopen(path.c_str(), "rb");
        if (!gz_) return;
        gzbuffer(gz_, 1024 * 1024);  // 1MB internal buffer
    }

    ~GzipStreamReader() {
        if (gz_) gzclose(gz_);
    }

    bool valid() const { return rg_reader_ != nullptr || gz_ != nullptr; }

    // Read next chunk, returns (data, size) — data valid until next call.
    std::pair<const char*, size_t> read_chunk() {
        if (rg_reader_) {
            const size_t n = rg_reader_->read_raw(buffer_.data(), CHUNK_SIZE);
            return n > 0 ? std::pair{buffer_.data(), n} : std::pair{nullptr, size_t{0}};
        }
        if (!gz_ || gzeof(gz_)) return {nullptr, 0};
        const int bytes = gzread(gz_, buffer_.data(), static_cast<unsigned>(CHUNK_SIZE));
        if (bytes <= 0) return {nullptr, 0};
        return {buffer_.data(), static_cast<size_t>(bytes)};
    }

    // Read entire file into memory.
    std::vector<char> read_all() {
        std::vector<char> result;
        result.reserve(100 * 1024 * 1024);
        while (true) {
            auto [data, size] = read_chunk();
            if (!data || size == 0) break;
            result.insert(result.end(), data, data + size);
        }
        return result;
    }

private:
    std::unique_ptr<dart::GzLineReader> rg_reader_;
    gzFile gz_ = nullptr;
    std::vector<char> buffer_;
};

/**
 * @brief Detect if file is gzip compressed
 */
inline bool is_gzip(const std::string& path) {
    // Check extension
    if (path.size() >= 3 && path.substr(path.size() - 3) == ".gz") return true;
    if (path.size() >= 4 && path.substr(path.size() - 4) == ".bgz") return true;

    // Check magic bytes
    FILE* f = fopen(path.c_str(), "rb");
    if (!f) return false;
    unsigned char magic[2] = {0, 0};
    size_t n = fread(magic, 1, 2, f);
    fclose(f);
    return n == 2 && magic[0] == 0x1f && magic[1] == 0x8b;
}

/**
 * @brief TSV to columnar converter (DuckDB-style streaming)
 *
 * Single-pass streaming algorithm:
 * Builds dictionaries inline while writing row-groups in one file read.
 * Indices pre-assigned by main thread; row-group processing is synchronous.
 *
 * Memory bounded to ~10GB regardless of file size.
 */
class FastColumnarWriter {
public:
    struct Config {
        float d_max = 0.0f;
        float lambda = 0.3f;
        bool d_max_set_by_user = false;
        bool lambda_set_by_user = false;
        bool verbose = false;
        bool compress = true;  // Use ZSTD compression (if available)
        int compression_level = 3;  // ZSTD level (1-19, 3 = good balance)
        bool skip_alignments = false;  // Skip qaln/taln strings (faster, but damage-annotate needs them)
        size_t memory_limit = 16ULL * 1024 * 1024 * 1024;  // 16GB default memory cap
        bool force_streaming = false;  // Force streaming even for small files
        std::string damage_index_path; // Optional .agd, embed per-read p_damaged into DAMAGE_SCORE
    };

    /**
     * @brief Convert TSV to columnar index with fully streaming I/O
     *
     * DuckDB-style bounded memory processing:
     * - Chunked file reading (64MB at a time, no full mmap)
     * - Morsel-driven parallelism (8K row batches)
     * - Spilling dictionary (disk overflow when >4GB)
     * - Memory usage: ~4-5GB regardless of file size
     *
     * @param tsv_path Input TSV file (16-column MMseqs2 format, optionally .gz)
     * @param output_path Output .emi file
     * @param config Configuration options
     * @return Number of alignments written
     */
    static uint64_t convert(
        const std::string& tsv_path,
        const std::string& output_path,
        const Config& config)
    {
        // Check file size to decide mode
        struct stat st;
        if (stat(tsv_path.c_str(), &st) != 0) {
            throw std::runtime_error("Cannot stat file: " + tsv_path);
        }
        size_t file_size = static_cast<size_t>(st.st_size);

        // For gzip, estimate 5x expansion
        const bool is_gz = is_gzip(tsv_path);
        size_t estimated_size = is_gz ? file_size * 5 : file_size;

        // Prefer streaming for gzip inputs to avoid full decompression into RAM.
        const bool use_mmap = !config.force_streaming && !is_gz && estimated_size < config.memory_limit;

        if (use_mmap) {
            if (config.verbose) {
                std::cerr << "Using mmap mode (file " << (estimated_size / (1024*1024))
                          << " MB < limit " << (config.memory_limit / (1024*1024*1024)) << " GB)\n";
            }
            return convert_mmap(tsv_path, output_path, config);
        } else {
            if (config.verbose) {
                std::cerr << "Using streaming mode (file " << (estimated_size / (1024*1024))
                          << " MB, limit " << (config.memory_limit / (1024*1024*1024)) << " GB";
                if (is_gz && !config.force_streaming) std::cerr << ", gzip input";
                std::cerr << ")\n";
            }
            return convert_streaming(tsv_path, output_path, config);
        }
    }

private:
    static std::unique_ptr<DamageIndexReader> init_damage_reader(
        const Config& config,
        float& d_max_out,
        float& lambda_out)
    {
        if (config.damage_index_path.empty()) return nullptr;

        auto reader = std::make_unique<DamageIndexReader>(config.damage_index_path, false);
        if (!reader->is_valid()) {
            throw std::runtime_error("Cannot open damage index: " + config.damage_index_path);
        }

        // Force AGD into OS page cache before main loop.
        // madvise(MADV_WILLNEED) is ignored by NFS on Linux 4.18; synchronous
        // warmup reads all pages sequentially (1.5 GB/s NFS → ~21s for 31 GB).
        // Without this, every random AGD lookup causes a blocking NFS page fault,
        // keeping the process in D-state and limiting throughput to ~6 MB/s.
        if (config.verbose) {
            std::cerr << "Warming AGD page cache ("
                      << (reader->record_count() * 32 / (1024*1024)) << " MB)...\n";
        }
        reader->warmup_cache();

        if (!config.d_max_set_by_user && reader->d_max() > 0.0f) {
            d_max_out = reader->d_max();
        }
        if (!config.lambda_set_by_user && reader->lambda() > 0.0f) {
            lambda_out = reader->lambda();
        }
        return reader;
    }

    template <typename NameVec>
    static std::vector<uint8_t> build_read_damage_quantized(
        const DamageIndexReader* damage_reader,
        const NameVec& read_names)
    {
        std::vector<uint8_t> out;
        if (!damage_reader || !damage_reader->is_valid()) return out;

        out.resize(read_names.size(), 0);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
        for (int64_t i = 0; i < static_cast<int64_t>(read_names.size()); ++i) {
            const auto ridx = static_cast<size_t>(i);
            if (const AgdRecord* rec = damage_reader->find(read_names[ridx])) {
                out[ridx] = rec->p_damaged_q;
            }
        }
        return out;
    }

public:

    /**
     * @brief Fast mmap-based conversion (for files that fit in memory)
     */
    static uint64_t convert_mmap(
        const std::string& tsv_path,
        const std::string& output_path,
        const Config& config)
    {
        auto total_start = std::chrono::high_resolution_clock::now();
        float header_d_max = config.d_max;
        float header_lambda = config.lambda;
        auto damage_reader = init_damage_reader(config, header_d_max, header_lambda);

        if (is_gzip(tsv_path)) {
            throw std::runtime_error(
                "mmap mode does not support gzip input; use streaming mode");
        }

        const std::unique_ptr<MappedTSV> tsv_holder = std::make_unique<MappedTSV>(tsv_path.c_str());
        if (!tsv_holder->valid()) {
            throw std::runtime_error("Cannot open TSV: " + tsv_path);
        }
        const char* data = tsv_holder->data();
        const size_t size = tsv_holder->size();

        // Skip header if present
        size_t start_offset = 0;
        if (size > 5 && strncmp(data, "query", 5) == 0) {
            const char* nl = static_cast<const char*>(memchr(data, '\n', size));
            if (nl) start_offset = static_cast<size_t>(nl - data) + 1;
        }

        const char* file_end = data + size;

        // Count lines (fast scan)
        auto scan_start = std::chrono::high_resolution_clock::now();
        size_t n_lines = 0;
        {
            const char* p = data + start_offset;
            while (p < file_end) {
                const char* nl = simd_find_newline(p, file_end);
                ++n_lines;
                if (nl >= file_end) break;
                p = nl + 1;
            }
        }

        auto scan_end = std::chrono::high_resolution_clock::now();
        auto scan_ms = std::chrono::duration_cast<std::chrono::milliseconds>(scan_end - scan_start).count();

        if (config.verbose) {
            std::cerr << "Lines: " << n_lines << " (scan: "
                      << dart::log_utils::format_duration_ms(scan_ms)
                      << "), building dictionary...\n";
        }

        auto pass_start = std::chrono::high_resolution_clock::now();

        // Parallel dictionary building
        int n_threads = 1;
#ifdef _OPENMP
        n_threads = omp_get_max_threads();
#endif

        // Divide file into byte ranges
        std::vector<size_t> chunk_starts(n_threads + 1);
        chunk_starts[0] = start_offset;
        for (int t = 1; t < n_threads; ++t) {
            size_t target = start_offset + (size - start_offset) * t / n_threads;
            const char* p = data + target;
            while (p < file_end && *p != '\n') ++p;
            chunk_starts[t] = (p < file_end) ? static_cast<size_t>(p - data + 1) : size;
        }
        chunk_starts[n_threads] = size;

        std::vector<std::vector<std::string_view>> local_read_uniques(n_threads);
        std::vector<std::vector<std::string_view>> local_ref_uniques(n_threads);
        std::vector<FlatSVU32Map> local_read_seen(n_threads);
        std::vector<FlatSVU32Map> local_ref_seen(n_threads);

#ifdef _OPENMP
#pragma omp parallel
#endif
        {
            int tid = 0;
#ifdef _OPENMP
            tid = omp_get_thread_num();
#endif
            const char* p = data + chunk_starts[tid];
            const char* chunk_end = data + chunk_starts[tid + 1];

            auto& read_seen = local_read_seen[tid];
            auto& ref_seen = local_ref_seen[tid];
            auto& read_u = local_read_uniques[tid];
            auto& ref_u = local_ref_uniques[tid];

            size_t est_lines = (chunk_starts[tid + 1] - chunk_starts[tid]) / 200;
            read_seen.reserve(std::max<size_t>(1024, est_lines / 2));
            ref_seen.reserve(std::max<size_t>(256, est_lines / 16));
            read_u.reserve(est_lines / 3);
            ref_u.reserve(est_lines / 16);

            while (p < chunk_end) {
                const char* line_end = simd_find_newline(p, chunk_end);
                std::string_view q, t;
                if (parse_first_two_fields(p, line_end, q, t)) {
                    uint64_t hq = FlatSVU32Map::hash_sv(q);
                    if (read_seen.try_emplace(q, 0u, hq)) read_u.push_back(q);
                    uint64_t ht = FlatSVU32Map::hash_sv(t);
                    if (ref_seen.try_emplace(t, 0u, ht)) ref_u.push_back(t);
                }
                p = (line_end < chunk_end) ? line_end + 1 : chunk_end;
            }
        }

        // Release local seen maps
        for (int t = 0; t < n_threads; ++t) {
            local_read_seen[t].clear_release();
            local_ref_seen[t].clear_release();
        }

        // Build read-order list by first appearance across file chunks.
        // local_read_uniques[tid] preserves first-seen order within each contiguous chunk.
        size_t est_read_uniques = 0;
        for (const auto& v : local_read_uniques) est_read_uniques += v.size();
        std::vector<std::string_view> read_first_seen_order;
        read_first_seen_order.reserve(est_read_uniques);
        FlatSVU32Map read_first_seen_set;
        read_first_seen_set.reserve(std::max<size_t>(1024, est_read_uniques));
        for (int tid = 0; tid < n_threads; ++tid) {
            for (std::string_view sv : local_read_uniques[tid]) {
                uint64_t h = FlatSVU32Map::hash_sv(sv);
                if (read_first_seen_set.try_emplace(sv, 0u, h)) {
                    read_first_seen_order.push_back(sv);
                }
            }
        }
        read_first_seen_set.clear_release();

        // Merge into global dictionaries
        constexpr size_t READ_SHARDS = 1024;
        constexpr size_t REF_SHARDS = 256;
        auto read_dict = build_global_dict_from_locals(local_read_uniques, READ_SHARDS);
        auto ref_dict = build_global_dict_from_locals(local_ref_uniques, REF_SHARDS);

        std::vector<std::string_view> read_names = std::move(read_dict.names);
        std::vector<std::string_view> ref_names = std::move(ref_dict.names);
        std::vector<uint32_t> read_old_to_new(
            read_names.size(), std::numeric_limits<uint32_t>::max());
        std::vector<uint32_t> read_new_to_old;
        read_new_to_old.reserve(read_names.size());

        // Remap read IDs by first appearance in input so row-group local sorting
        // can preserve global (read_idx, bit_score DESC) ordering for grouped inputs.
        for (std::string_view sv : read_first_seen_order) {
            uint32_t old_idx = dict_lookup(read_dict, sv);
            uint32_t& mapped = read_old_to_new[old_idx];
            if (mapped == std::numeric_limits<uint32_t>::max()) {
                mapped = static_cast<uint32_t>(read_new_to_old.size());
                read_new_to_old.push_back(old_idx);
            }
        }
        std::vector<std::string_view>().swap(read_first_seen_order);
        if (read_new_to_old.size() != read_names.size()) {
            throw std::runtime_error("Read-ID remap incomplete while building EMI");
        }
        std::vector<std::string_view> read_names_ordered(read_names.size());
        for (uint32_t new_idx = 0; new_idx < read_new_to_old.size(); ++new_idx) {
            read_names_ordered[new_idx] = read_names[read_new_to_old[new_idx]];
        }
        read_names.swap(read_names_ordered);

        std::vector<uint8_t> read_damage_q = build_read_damage_quantized(damage_reader.get(), read_names);
        if (config.verbose && damage_reader) {
            std::cerr << "  Embedded per-read damage scores from " << config.damage_index_path
                      << " into EMI DAMAGE_SCORE\n";
            if (!config.d_max_set_by_user || !config.lambda_set_by_user) {
                std::cerr << "  Using AGD metadata for EMI header: d_max=" << header_d_max
                          << " lambda=" << header_lambda << "\n";
            }
        }

        auto dict_end = std::chrono::high_resolution_clock::now();
        auto dict_ms = std::chrono::duration_cast<std::chrono::milliseconds>(dict_end - pass_start).count();
        if (config.verbose) {
            std::cerr << "  Dictionary: " << dart::log_utils::format_duration_ms(dict_ms)
                      << " (" << read_names.size() << " reads, "
                      << ref_names.size() << " refs, " << n_threads << " threads)\n";
        }

        // Build name offset tables
        uint64_t total_rows = n_lines;
        const size_t rg_size = ROW_GROUP_SIZE;
        const size_t num_row_groups = (n_lines + rg_size - 1) / rg_size;

        std::vector<uint32_t> read_name_offsets;
        std::vector<char> read_names_data;
        read_name_offsets.reserve(read_names.size());
        for (const auto& sv : read_names) {
            read_name_offsets.push_back(static_cast<uint32_t>(read_names_data.size()));
            read_names_data.insert(read_names_data.end(), sv.begin(), sv.end());
            read_names_data.push_back('\0');
        }

        std::vector<uint32_t> ref_name_offsets;
        std::vector<char> ref_names_data;
        ref_name_offsets.reserve(ref_names.size());
        for (const auto& sv : ref_names) {
            ref_name_offsets.push_back(static_cast<uint32_t>(ref_names_data.size()));
            ref_names_data.insert(ref_names_data.end(), sv.begin(), sv.end());
            ref_names_data.push_back('\0');
        }

        std::vector<std::string_view>().swap(read_names);
        std::vector<std::string_view>().swap(ref_names);

        // Open output and write header
        std::ofstream out(output_path, std::ios::binary);
        if (!out) throw std::runtime_error("Cannot create file: " + output_path);

        EMIHeader header{};
        header.magic = EMI_MAGIC;
        header.version = EMI_VERSION;
        header.num_alignments = total_rows;
        header.num_reads = static_cast<uint32_t>(read_name_offsets.size());
        header.num_refs = static_cast<uint32_t>(ref_name_offsets.size());
        header.num_row_groups = static_cast<uint32_t>(num_row_groups);
        header.row_group_size = static_cast<uint32_t>(rg_size);
        header.d_max = header_d_max;
        header.lambda = header_lambda;
        header.flags = 0;

        out.write(reinterpret_cast<const char*>(&header), sizeof(header));

        header.row_group_dir_offset = out.tellp();
        std::vector<RowGroupMeta> row_groups(num_row_groups);
        out.write(reinterpret_cast<const char*>(row_groups.data()), row_groups.size() * sizeof(RowGroupMeta));
        std::vector<RowGroupReadIndex> row_group_read_index(num_row_groups);

        header.string_dict_offset = out.tellp();
        uint32_t read_names_count = static_cast<uint32_t>(read_name_offsets.size());
        out.write(reinterpret_cast<const char*>(&read_names_count), sizeof(read_names_count));
        out.write(reinterpret_cast<const char*>(read_name_offsets.data()), read_name_offsets.size() * sizeof(uint32_t));
        uint32_t read_names_size_val = static_cast<uint32_t>(read_names_data.size());
        out.write(reinterpret_cast<const char*>(&read_names_size_val), sizeof(read_names_size_val));
        out.write(read_names_data.data(), read_names_data.size());

        uint32_t ref_names_count = static_cast<uint32_t>(ref_name_offsets.size());
        out.write(reinterpret_cast<const char*>(&ref_names_count), sizeof(ref_names_count));
        out.write(reinterpret_cast<const char*>(ref_name_offsets.data()), ref_name_offsets.size() * sizeof(uint32_t));
        uint32_t ref_names_size_val = static_cast<uint32_t>(ref_names_data.size());
        out.write(reinterpret_cast<const char*>(&ref_names_size_val), sizeof(ref_names_size_val));
        out.write(ref_names_data.data(), ref_names_data.size());

        std::vector<uint32_t>().swap(read_name_offsets);
        std::vector<char>().swap(read_names_data);
        std::vector<uint32_t>().swap(ref_name_offsets);
        std::vector<char>().swap(ref_names_data);

        // CSR offset building: count alignments per read
        std::vector<uint32_t> read_counts(header.num_reads, 0);
        uint64_t global_row_cursor = 0;
        bool globally_sorted_by_read_then_score = true;
        bool has_prev_sort_key = false;
        uint32_t prev_read_idx = 0;
        float prev_bit_score = 0.0f;

        // Row-group buffers
        std::vector<uint32_t> rg_read_idx(rg_size);
        std::vector<uint32_t> rg_ref_idx(rg_size);
        std::vector<float> rg_bit_score(rg_size);
        std::vector<float> rg_damage_score(rg_size);
        std::vector<float> rg_evalue_log10(rg_size);
        std::vector<uint16_t> rg_dmg_k(rg_size);
        std::vector<uint16_t> rg_dmg_m(rg_size);
        std::vector<float> rg_dmg_ll_a(rg_size);
        std::vector<float> rg_dmg_ll_m(rg_size);
        std::vector<uint16_t> rg_identity_q(rg_size);
        std::vector<uint16_t> rg_aln_len(rg_size);
        std::vector<uint16_t> rg_qstart(rg_size);
        std::vector<uint16_t> rg_qend(rg_size);
        std::vector<uint16_t> rg_tstart(rg_size);
        std::vector<uint16_t> rg_tend(rg_size);
        std::vector<uint16_t> rg_tlen(rg_size);
        std::vector<uint16_t> rg_qlen(rg_size);
        std::vector<uint16_t> rg_mismatch(rg_size);
        std::vector<uint16_t> rg_gapopen(rg_size);
        std::vector<uint32_t> perm_tmp_u32(rg_size);
        std::vector<float> perm_tmp_f32(rg_size);
        std::vector<uint16_t> perm_tmp_u16(rg_size);
        std::vector<uint32_t> perm_qaln_offsets(rg_size + 1);
        std::vector<uint32_t> perm_taln_offsets(rg_size + 1);
        std::vector<uint32_t> rg_qaln_offsets(rg_size + 1);
        std::vector<char> rg_qaln_data;
        rg_qaln_data.reserve(rg_size * 50);
        std::vector<uint32_t> rg_taln_offsets(rg_size + 1);
        std::vector<char> rg_taln_data;
        rg_taln_data.reserve(rg_size * 50);

        // Process row-groups
        const char* file_pos = data + start_offset;

        for (size_t rg = 0; rg < num_row_groups; ++rg) {
            size_t rg_start = rg * rg_size;
            size_t rg_end_idx = std::min(rg_start + rg_size, n_lines);
            uint32_t rg_rows = static_cast<uint32_t>(rg_end_idx - rg_start);

            rg_qaln_data.clear();
            rg_taln_data.clear();

            // Collect line boundaries
            std::vector<const char*> rg_line_starts(rg_rows + 1);
            {
                const char* p = file_pos;
                for (uint32_t i = 0; i < rg_rows && p < file_end; ++i) {
                    rg_line_starts[i] = p;
                    const char* nl = simd_find_newline(p, file_end);
                    p = (nl < file_end) ? nl + 1 : file_end;
                }
                rg_line_starts[rg_rows] = p;
                file_pos = p;
            }

            // Parse numeric columns (parallel)
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
            for (uint32_t local_idx = 0; local_idx < rg_rows; ++local_idx) {
                const char* p = rg_line_starts[local_idx];
                const char* line_end = rg_line_starts[local_idx + 1];
                if (line_end > p && *(line_end - 1) == '\n') --line_end;

                const char* fields[16];
                size_t lens[16];
                const char* fp = p;
                for (int f = 0; f < 16; ++f) {
                    fields[f] = fp;
                    const char* tab = (f == 15) ? line_end : simd_find_tab(fp, line_end);
                    lens[f] = static_cast<size_t>(tab - fp);
                    fp = tab + 1;
                }

                std::string_view q(fields[0], lens[0]);
                std::string_view t(fields[1], lens[1]);
                uint32_t old_read_idx = dict_lookup(read_dict, q);
                uint32_t mapped_read_idx = read_old_to_new[old_read_idx];
                rg_read_idx[local_idx] =
                    (mapped_read_idx != std::numeric_limits<uint32_t>::max())
                        ? mapped_read_idx
                        : 0;
                rg_ref_idx[local_idx] = dict_lookup(ref_dict, t);

                float fident = fast_stof({fields[2], lens[2]});
                uint32_t alnlen = fast_stou({fields[3], lens[3]});
                uint32_t mismatch_v = fast_stou({fields[4], lens[4]});
                uint32_t gapopen_v = fast_stou({fields[5], lens[5]});
                uint32_t qstart_v = fast_stou({fields[6], lens[6]});
                uint32_t qend_v = fast_stou({fields[7], lens[7]});
                uint32_t tstart_v = fast_stou({fields[8], lens[8]});
                uint32_t tend_v = fast_stou({fields[9], lens[9]});
                float evalue = fast_stof({fields[10], lens[10]});
                float bits = fast_stof({fields[11], lens[11]});
                uint32_t qlen_v = fast_stou({fields[12], lens[12]});
                uint32_t tlen_v = fast_stou({fields[13], lens[13]});

                rg_bit_score[local_idx] = bits;
                rg_damage_score[local_idx] = read_damage_q.empty()
                    ? 0.0f
                    : static_cast<float>(read_damage_q[rg_read_idx[local_idx]]) * (1.0f / 255.0f);
                rg_evalue_log10[local_idx] = (evalue > 0) ? std::log10(evalue) : -300.0f;
                rg_dmg_k[local_idx] = 0;
                rg_dmg_m[local_idx] = 0;
                rg_dmg_ll_a[local_idx] = 0.0f;
                rg_dmg_ll_m[local_idx] = 0.0f;
                rg_identity_q[local_idx] = static_cast<uint16_t>(std::clamp(fident, 0.0f, 1.0f) * 65535.0f);
                rg_aln_len[local_idx] = static_cast<uint16_t>(std::min(alnlen, 65535u));
                rg_qstart[local_idx] = static_cast<uint16_t>(std::min(qstart_v, 65535u));
                rg_qend[local_idx] = static_cast<uint16_t>(std::min(qend_v, 65535u));
                rg_tstart[local_idx] = static_cast<uint16_t>(std::min(tstart_v, 65535u));
                rg_tend[local_idx] = static_cast<uint16_t>(std::min(tend_v, 65535u));
                rg_tlen[local_idx] = static_cast<uint16_t>(std::min(tlen_v, 65535u));
                rg_qlen[local_idx] = static_cast<uint16_t>(std::min(qlen_v, 65535u));
                rg_mismatch[local_idx] = static_cast<uint16_t>(std::min(mismatch_v, 65535u));
                rg_gapopen[local_idx] = static_cast<uint16_t>(std::min(gapopen_v, 65535u));
            }

            // Parse string columns (sequential)
            if (!config.skip_alignments) {
                for (uint32_t local_idx = 0; local_idx < rg_rows; ++local_idx) {
                    const char* p = rg_line_starts[local_idx];
                    const char* line_end = rg_line_starts[local_idx + 1];
                    if (line_end > p && *(line_end - 1) == '\n') --line_end;

                    const char* qid_start = p;
                    const char* qid_end = simd_find_tab(qid_start, line_end);
                    std::string_view query_id(qid_start, static_cast<size_t>(qid_end - qid_start));

                    const char* fp = (qid_end < line_end) ? (qid_end + 1) : line_end;
                    for (int f = 1; f < 14 && fp < line_end; ++f) {
                        const char* tab = simd_find_tab(fp, line_end);
                        fp = (tab < line_end) ? (tab + 1) : line_end;
                    }
                    const char* qaln_start = fp;
                    const char* qaln_end = simd_find_tab(fp, line_end);
                    fp = (qaln_end < line_end) ? (qaln_end + 1) : line_end;
                    const char* taln_start = fp;
                    const char* taln_end = simd_find_tab(fp, line_end);

                    rg_qaln_offsets[local_idx] = static_cast<uint32_t>(rg_qaln_data.size());
                    rg_qaln_data.insert(rg_qaln_data.end(), qaln_start, qaln_end);
                    rg_taln_offsets[local_idx] = static_cast<uint32_t>(rg_taln_data.size());
                    rg_taln_data.insert(rg_taln_data.end(), taln_start, taln_end);

                    const std::string_view qaln_sv(qaln_start, static_cast<size_t>(qaln_end - qaln_start));
                    const std::string_view taln_sv(taln_start, static_cast<size_t>(taln_end - taln_start));
                    const auto dmg = compute_alignment_damage_evidence(
                        query_id,
                        qaln_sv,
                        taln_sv,
                        rg_qstart[local_idx],
                        rg_qlen[local_idx],
                        header_d_max,
                        header_lambda);
                    rg_dmg_k[local_idx] = dmg.k_hits;
                    rg_dmg_m[local_idx] = dmg.m_sites;
                    rg_dmg_ll_a[local_idx] = dmg.ll_ancient;
                    rg_dmg_ll_m[local_idx] = dmg.ll_modern;
                }
            }
            rg_qaln_offsets[rg_rows] = static_cast<uint32_t>(rg_qaln_data.size());
            rg_taln_offsets[rg_rows] = static_cast<uint32_t>(rg_taln_data.size());

            // Sort row-group by (read_idx, bit_score DESC) for correct CSR offsets
            // CSR requires alignments to be grouped by read, and within each read
            // the best hit (highest bit_score) must come first for damage-annotate
            std::vector<uint32_t> perm(rg_rows);
            std::iota(perm.begin(), perm.end(), 0);
            std::sort(perm.begin(), perm.end(), [&](uint32_t a, uint32_t b) {
                if (rg_read_idx[a] != rg_read_idx[b]) {
                    return rg_read_idx[a] < rg_read_idx[b];
                }
                return rg_bit_score[a] > rg_bit_score[b];  // Descending by score
            });

            // Apply permutation with reusable scratch buffers (avoids per-column allocations)
            auto apply_perm_u32 = [&](std::vector<uint32_t>& vec) {
                for (uint32_t i = 0; i < rg_rows; ++i) perm_tmp_u32[i] = vec[perm[i]];
                std::copy_n(perm_tmp_u32.data(), rg_rows, vec.data());
            };
            auto apply_perm_f32 = [&](std::vector<float>& vec) {
                for (uint32_t i = 0; i < rg_rows; ++i) perm_tmp_f32[i] = vec[perm[i]];
                std::copy_n(perm_tmp_f32.data(), rg_rows, vec.data());
            };
            auto apply_perm_u16 = [&](std::vector<uint16_t>& vec) {
                for (uint32_t i = 0; i < rg_rows; ++i) perm_tmp_u16[i] = vec[perm[i]];
                std::copy_n(perm_tmp_u16.data(), rg_rows, vec.data());
            };
            apply_perm_u32(rg_read_idx);
            apply_perm_u32(rg_ref_idx);
            apply_perm_f32(rg_bit_score);
            apply_perm_f32(rg_damage_score);
            apply_perm_f32(rg_evalue_log10);
            apply_perm_u16(rg_dmg_k);
            apply_perm_u16(rg_dmg_m);
            apply_perm_f32(rg_dmg_ll_a);
            apply_perm_f32(rg_dmg_ll_m);
            apply_perm_u16(rg_identity_q);
            apply_perm_u16(rg_aln_len);
            apply_perm_u16(rg_qstart);
            apply_perm_u16(rg_qend);
            apply_perm_u16(rg_tstart);
            apply_perm_u16(rg_tend);
            apply_perm_u16(rg_tlen);
            apply_perm_u16(rg_qlen);
            apply_perm_u16(rg_mismatch);
            apply_perm_u16(rg_gapopen);

            // Rebuild qaln/taln with permuted order
            if (!config.skip_alignments) {
                std::vector<char> perm_qaln_data;
                std::vector<char> perm_taln_data;
                perm_qaln_data.reserve(rg_qaln_data.size());
                perm_taln_data.reserve(rg_taln_data.size());
                for (uint32_t i = 0; i < rg_rows; ++i) {
                    uint32_t old_i = perm[i];
                    uint32_t q_start = rg_qaln_offsets[old_i];
                    uint32_t q_end = rg_qaln_offsets[old_i + 1];
                    uint32_t t_start = rg_taln_offsets[old_i];
                    uint32_t t_end = rg_taln_offsets[old_i + 1];
                    perm_qaln_offsets[i] = static_cast<uint32_t>(perm_qaln_data.size());
                    perm_qaln_data.insert(
                        perm_qaln_data.end(),
                        rg_qaln_data.begin() + q_start,
                        rg_qaln_data.begin() + q_end);
                    perm_taln_offsets[i] = static_cast<uint32_t>(perm_taln_data.size());
                    perm_taln_data.insert(
                        perm_taln_data.end(),
                        rg_taln_data.begin() + t_start,
                        rg_taln_data.begin() + t_end);
                }
                perm_qaln_offsets[rg_rows] = static_cast<uint32_t>(perm_qaln_data.size());
                perm_taln_offsets[rg_rows] = static_cast<uint32_t>(perm_taln_data.size());
                std::copy_n(perm_qaln_offsets.data(), rg_rows + 1, rg_qaln_offsets.data());
                std::copy_n(perm_taln_offsets.data(), rg_rows + 1, rg_taln_offsets.data());
                rg_qaln_data.swap(perm_qaln_data);
                rg_taln_data.swap(perm_taln_data);
            }

            // Accumulate CSR counts (now sorted)
            for (uint32_t i = 0; i < rg_rows; ++i) {
                if (!has_prev_sort_key) {
                    has_prev_sort_key = true;
                } else if (globally_sorted_by_read_then_score) {
                    if (rg_read_idx[i] < prev_read_idx ||
                        (rg_read_idx[i] == prev_read_idx && rg_bit_score[i] > prev_bit_score)) {
                        globally_sorted_by_read_then_score = false;
                    }
                }
                prev_read_idx = rg_read_idx[i];
                prev_bit_score = rg_bit_score[i];
                read_counts[rg_read_idx[i]]++;
            }

            // Write row-group
            row_groups[rg].num_rows = rg_rows;
            row_groups[rg].first_read_idx = rg_read_idx[0];
            row_group_read_index[rg].first_read_idx = rg_read_idx[0];
            row_group_read_index[rg].last_read_idx = rg_read_idx[rg_rows - 1];
            row_group_read_index[rg].start_row = global_row_cursor;
            row_group_read_index[rg].end_row_exclusive = global_row_cursor + rg_rows;
            global_row_cursor += rg_rows;

            auto write_chunk = [&](ColumnID col, const void* chunk_data, size_t elem_size) {
                auto& chunk = row_groups[rg].columns[static_cast<size_t>(col)];
                chunk.offset = out.tellp();
                chunk.uncompressed_size = static_cast<uint32_t>(elem_size * rg_rows);
                chunk.compressed_size = chunk.uncompressed_size;
                chunk.codec = Codec::NONE;
                out.write(static_cast<const char*>(chunk_data), chunk.uncompressed_size);
            };
            auto write_alignment_chunk = [&](ColumnID col,
                                             const std::vector<uint32_t>& offsets,
                                             const std::vector<char>& data) {
                auto& chunk = row_groups[rg].columns[static_cast<size_t>(col)];
                chunk.offset = out.tellp();

                const size_t offsets_bytes = static_cast<size_t>(rg_rows + 1) * sizeof(uint32_t);
                const size_t raw_bytes = offsets_bytes + data.size();
                chunk.uncompressed_size = static_cast<uint32_t>(raw_bytes);

#ifdef HAVE_ZSTD
                if (config.compress && raw_bytes > 0) {
                    std::vector<char> raw(raw_bytes);
                    std::memcpy(raw.data(), offsets.data(), offsets_bytes);
                    if (!data.empty()) {
                        std::memcpy(raw.data() + offsets_bytes, data.data(), data.size());
                    }
                    auto compressed = compress_zstd(raw.data(), raw.size(), config.compression_level);
                    if (compressed.size() < raw.size()) {
                        chunk.compressed_size = static_cast<uint32_t>(compressed.size());
                        chunk.codec = Codec::ZSTD;
                        out.write(compressed.data(), static_cast<std::streamsize>(compressed.size()));
                        return;
                    }
                }
#endif
                chunk.compressed_size = chunk.uncompressed_size;
                chunk.codec = Codec::NONE;
                out.write(reinterpret_cast<const char*>(offsets.data()),
                          static_cast<std::streamsize>(offsets_bytes));
                if (!data.empty()) {
                    out.write(data.data(), static_cast<std::streamsize>(data.size()));
                }
            };

            write_chunk(ColumnID::READ_IDX, rg_read_idx.data(), sizeof(uint32_t));
            write_chunk(ColumnID::REF_IDX, rg_ref_idx.data(), sizeof(uint32_t));
            write_chunk(ColumnID::BIT_SCORE, rg_bit_score.data(), sizeof(float));
            write_chunk(ColumnID::DAMAGE_SCORE, rg_damage_score.data(), sizeof(float));
            write_chunk(ColumnID::EVALUE_LOG10, rg_evalue_log10.data(), sizeof(float));
            write_chunk(ColumnID::DMG_K, rg_dmg_k.data(), sizeof(uint16_t));
            write_chunk(ColumnID::DMG_M, rg_dmg_m.data(), sizeof(uint16_t));
            write_chunk(ColumnID::DMG_LL_A, rg_dmg_ll_a.data(), sizeof(float));
            write_chunk(ColumnID::DMG_LL_M, rg_dmg_ll_m.data(), sizeof(float));
            write_chunk(ColumnID::IDENTITY_Q, rg_identity_q.data(), sizeof(uint16_t));
            write_chunk(ColumnID::ALN_LEN, rg_aln_len.data(), sizeof(uint16_t));
            write_chunk(ColumnID::QSTART, rg_qstart.data(), sizeof(uint16_t));
            write_chunk(ColumnID::QEND, rg_qend.data(), sizeof(uint16_t));
            write_chunk(ColumnID::TSTART, rg_tstart.data(), sizeof(uint16_t));
            write_chunk(ColumnID::TEND, rg_tend.data(), sizeof(uint16_t));
            write_chunk(ColumnID::TLEN, rg_tlen.data(), sizeof(uint16_t));
            write_chunk(ColumnID::QLEN, rg_qlen.data(), sizeof(uint16_t));
            write_chunk(ColumnID::MISMATCH, rg_mismatch.data(), sizeof(uint16_t));
            write_chunk(ColumnID::GAPOPEN, rg_gapopen.data(), sizeof(uint16_t));

            if (!config.skip_alignments) {
                write_alignment_chunk(ColumnID::QALN, rg_qaln_offsets, rg_qaln_data);
                write_alignment_chunk(ColumnID::TALN, rg_taln_offsets, rg_taln_data);
            }
        }

        // Cleanup dictionaries
        for (auto& shard : read_dict.shard_maps) shard.clear_release();
        for (auto& shard : ref_dict.shard_maps) shard.clear_release();

        if (!config.skip_alignments) {
            header.flags |= EMI_FLAG_HAS_ALIGNMENT_STRINGS;
        }
        header.flags |= EMI_FLAG_HAS_ROW_GROUP_READ_INDEX;
        if (globally_sorted_by_read_then_score) {
            header.flags |= EMI_FLAG_SORTED_BY_READ_THEN_SCORE;
        } else if (config.verbose) {
            std::cerr << "Warning: EMI is not globally sorted by (read, bit_score). "
                         "Input is likely not grouped by query.\n";
        }

        header.row_group_read_index_offset = out.tellp();
        header.row_group_read_index_count = static_cast<uint32_t>(row_group_read_index.size());
        header.row_group_read_index_entry_size = sizeof(RowGroupReadIndex);
        out.write(reinterpret_cast<const char*>(row_group_read_index.data()),
                  row_group_read_index.size() * sizeof(RowGroupReadIndex));

        // Build CSR offsets from counts (prefix sum)
        header.csr_offsets_offset = out.tellp();
        uint32_t csr_count = header.num_reads + 1;
        std::vector<uint64_t> csr_offsets(csr_count);
        csr_offsets[0] = 0;
        for (uint32_t r = 0; r < header.num_reads; ++r) {
            csr_offsets[r + 1] = csr_offsets[r] + read_counts[r];
        }
        out.write(reinterpret_cast<const char*>(csr_offsets.data()), csr_count * sizeof(uint64_t));

        // Update header
        out.seekp(header.row_group_dir_offset);
        out.write(reinterpret_cast<const char*>(row_groups.data()), row_groups.size() * sizeof(RowGroupMeta));
        out.seekp(0);
        out.write(reinterpret_cast<const char*>(&header), sizeof(header));

        auto total_end = std::chrono::high_resolution_clock::now();
        auto total_ms = std::chrono::duration_cast<std::chrono::milliseconds>(total_end - total_start).count();

        if (config.verbose) {
            std::cerr << "Parsed: " << total_rows << " rows, "
                      << header.num_reads << " reads, "
                      << header.num_refs << " refs in "
                      << dart::log_utils::format_duration_ms(total_ms)
                      << " ("
                      << (total_rows * 1000 / std::max(1LL, static_cast<long long>(total_ms)))
                      << " rows/sec)\n";
        }

        return total_rows;
    }

    /**
     * @brief Streaming conversion for files larger than memory limit
     */
    static uint64_t convert_streaming(
        const std::string& tsv_path,
        const std::string& output_path,
        const Config& config)
    {
        auto total_start = std::chrono::high_resolution_clock::now();
        float header_d_max = config.d_max;
        float header_lambda = config.lambda;
        auto damage_reader = init_damage_reader(config, header_d_max, header_lambda);

        if (config.verbose) {
            std::cerr << "Streaming mode (single-pass)\n";
            if (damage_reader && (!config.d_max_set_by_user || !config.lambda_set_by_user)) {
                std::cerr << "Using AGD metadata for EMI header: d_max=" << header_d_max
                          << " lambda=" << header_lambda << "\n";
            }
        }

        // =====================================================================
        // Single-pass: build dicts + write row groups in one file read.
        // Layout: [header(provisional)] [row_group_data...] [rg_read_index]
        //         [csr_offsets] [string_dicts_v6] [row_group_dir]
        // Header overwritten at end via seekp(0) with final values.
        //
        // Previously: 2 file reads (Pass 1 dict scan + Pass 2 row-group write)
        // Now:        1 file read — dicts built inline, indices pre-assigned
        //             before each row group is processed.
        // =====================================================================

        std::ofstream out(output_path, std::ios::binary);
        if (!out) {
            throw std::runtime_error("Cannot create file: " + output_path);
        }

        // Provisional header (unknown fields are zero; overwritten at end)
        EMIHeader header{};
        header.magic = EMI_MAGIC;
        header.version = EMI_VERSION;
        header.d_max = header_d_max;
        header.lambda = header_lambda;
        out.write(reinterpret_cast<const char*>(&header), sizeof(header));

        ChunkedFileReader reader(tsv_path);
        if (!reader.valid()) {
            throw std::runtime_error("Cannot open TSV: " + tsv_path);
        }
        reader.skip_header_if_present();

        // Compact flat map for refs: key index into ref_names_vec + value id.
        struct FlatRefMap {
            struct Slot {
                uint32_t key_idx = empty_sentinel();
                uint32_t value = 0;
            };

            std::vector<Slot> slots;
            size_t mask = 0;
            size_t size = 0;

            static uint32_t empty_sentinel() { return 0xFFFFFFFFu; }

            static uint64_t hash_sv(std::string_view sv) {
                uint64_t h = 14695981039346656037ULL;
                for (unsigned char c : sv) {
                    h ^= c;
                    h *= 1099511628211ULL;
                }
                return h;
            }

            static size_t next_pow2(size_t x) {
                size_t p = 1;
                while (p < x) p <<= 1;
                return p;
            }

            void reserve(size_t expected_entries) {
                const size_t min_buckets = std::max<size_t>(16, (expected_entries * 10) / 6 + 1);
                const size_t cap = next_pow2(min_buckets);
                slots.assign(cap, Slot{});
                mask = cap - 1;
                size = 0;
            }

            const uint32_t* find(std::string_view key,
                                 const std::vector<std::string>& keys) const {
                if (slots.empty()) return nullptr;
                size_t i = static_cast<size_t>(hash_sv(key)) & mask;
                while (true) {
                    const Slot& s = slots[i];
                    if (s.key_idx == empty_sentinel()) return nullptr;
                    if (keys[s.key_idx] == key) return &s.value;
                    i = (i + 1) & mask;
                }
            }

            void emplace(uint32_t key_idx,
                         uint32_t value,
                         const std::vector<std::string>& keys) {
                if (slots.empty()) reserve(1024);
                if ((size + 1) * 10 >= slots.size() * 6) {
                    rehash(slots.size() * 2, keys);
                }

                const std::string_view key = keys[key_idx];
                size_t i = static_cast<size_t>(hash_sv(key)) & mask;
                while (true) {
                    Slot& s = slots[i];
                    if (s.key_idx == empty_sentinel()) {
                        s.key_idx = key_idx;
                        s.value = value;
                        ++size;
                        return;
                    }
                    if (keys[s.key_idx] == key) {
                        s.value = value;
                        return;
                    }
                    i = (i + 1) & mask;
                }
            }

        private:
            void rehash(size_t new_cap,
                        const std::vector<std::string>& keys) {
                const size_t cap = next_pow2(new_cap);
                std::vector<Slot> old;
                old.swap(slots);
                slots.assign(cap, Slot{});
                mask = cap - 1;
                size = 0;

                for (const Slot& old_slot : old) {
                    if (old_slot.key_idx == empty_sentinel()) continue;
                    const std::string_view key = keys[old_slot.key_idx];
                    size_t i = static_cast<size_t>(hash_sv(key)) & mask;
                    while (slots[i].key_idx != empty_sentinel()) i = (i + 1) & mask;
                    slots[i] = old_slot;
                    ++size;
                }
            }
        };

        // Inline dictionaries (first-appearance order for both reads and refs)
        std::vector<std::string> read_names_vec;
        read_names_vec.reserve(1 << 26);  // 64M slots — avoids reallocs up to ~64M reads
        FlatRefMap ref_lookup;
        ref_lookup.reserve(1 << 24);   // Pre-size for ~16M refs at ~0.6 load factor
        std::vector<std::string> ref_names_vec;
        // Quantized per-read damage score — grows in sync with read_names_vec
        std::vector<uint8_t> read_damage_q;
        read_damage_q.reserve(1 << 26);
        // Cache for consecutive identical queries (MMseqs2 output sorted by query)
        std::string last_read_name;
        uint32_t last_read_idx_cached = 0;

        const size_t rg_size = ROW_GROUP_SIZE;
        std::vector<RowGroupMeta> row_groups;
        std::vector<RowGroupReadIndex> row_group_read_index_vec;
        // Per-read alignment counts for CSR prefix sum — grows with read_names_vec
        std::vector<uint32_t> read_counts;
        read_counts.reserve(1 << 26);

        uint64_t n_lines = 0;
        uint64_t global_row_cursor = 0;
        bool globally_sorted = true;
        bool has_prev = false;
        uint32_t prev_read_idx_s = 0;
        float prev_bit_score_s = 0.0f;

        // RowGroupTask stores pre-parsed columns built in the main loop so
        // process_row_group only needs to sort + finalize write buffers.
        struct RowGroupTask {
            size_t rg_idx = 0;
            std::vector<uint32_t> read_idx, ref_idx;
            std::vector<float> bit_score, damage_score, evalue_log10;
            std::vector<uint16_t> dmg_k, dmg_m;
            std::vector<float> dmg_ll_a, dmg_ll_m;
            std::vector<uint16_t> identity_q, aln_len, qstart, qend, tstart, tend, tlen, qlen;
            std::vector<uint16_t> mismatch, gapopen;
            std::vector<uint32_t> qaln_offsets, taln_offsets;
            std::vector<char> qaln_data, taln_data;
        };
        struct RowGroupResult {
            size_t rg_idx = 0;
            uint32_t rg_rows = 0;
            std::vector<uint32_t> read_idx, ref_idx;
            std::vector<float> bit_score, damage_score, evalue_log10;
            std::vector<uint16_t> dmg_k, dmg_m;
            std::vector<float> dmg_ll_a, dmg_ll_m;
            std::vector<uint16_t> identity_q, aln_len, qstart, qend, tstart, tend, tlen, qlen;
            std::vector<uint16_t> mismatch, gapopen;
            std::vector<uint32_t> qaln_offsets, taln_offsets;
            std::vector<char> qaln_data, taln_data;
            std::vector<std::pair<uint32_t, uint32_t>> read_runs;
            uint32_t first_read_idx = 0;
            uint32_t last_read_idx = 0;
        };

        auto make_row_group_task = [&](size_t rg_idx) -> RowGroupTask {
            RowGroupTask task;
            task.rg_idx = rg_idx;
            task.read_idx.reserve(rg_size);
            task.ref_idx.reserve(rg_size);
            task.bit_score.reserve(rg_size);
            task.damage_score.reserve(rg_size);
            task.evalue_log10.reserve(rg_size);
            task.dmg_k.reserve(rg_size);
            task.dmg_m.reserve(rg_size);
            task.dmg_ll_a.reserve(rg_size);
            task.dmg_ll_m.reserve(rg_size);
            task.identity_q.reserve(rg_size);
            task.aln_len.reserve(rg_size);
            task.qstart.reserve(rg_size);
            task.qend.reserve(rg_size);
            task.tstart.reserve(rg_size);
            task.tend.reserve(rg_size);
            task.tlen.reserve(rg_size);
            task.qlen.reserve(rg_size);
            task.mismatch.reserve(rg_size);
            task.gapopen.reserve(rg_size);
            if (!config.skip_alignments) {
                task.qaln_offsets.reserve(rg_size + 1);
                task.taln_offsets.reserve(rg_size + 1);
                task.qaln_data.reserve(rg_size * 50);
                task.taln_data.reserve(rg_size * 50);
            }
            return task;
        };

        auto process_row_group = [&](RowGroupTask&& task) -> RowGroupResult {
            RowGroupResult res;
            res.rg_idx = task.rg_idx;
            res.rg_rows = static_cast<uint32_t>(task.read_idx.size());
            const uint32_t rg_rows = res.rg_rows;

            res.read_idx = std::move(task.read_idx);
            res.ref_idx = std::move(task.ref_idx);
            res.bit_score = std::move(task.bit_score);
            res.damage_score = std::move(task.damage_score);
            res.evalue_log10 = std::move(task.evalue_log10);
            res.dmg_k = std::move(task.dmg_k);
            res.dmg_m = std::move(task.dmg_m);
            res.dmg_ll_a = std::move(task.dmg_ll_a);
            res.dmg_ll_m = std::move(task.dmg_ll_m);
            res.identity_q = std::move(task.identity_q);
            res.aln_len = std::move(task.aln_len);
            res.qstart = std::move(task.qstart);
            res.qend = std::move(task.qend);
            res.tstart = std::move(task.tstart);
            res.tend = std::move(task.tend);
            res.tlen = std::move(task.tlen);
            res.qlen = std::move(task.qlen);
            res.mismatch = std::move(task.mismatch);
            res.gapopen = std::move(task.gapopen);
            res.qaln_offsets = std::move(task.qaln_offsets);
            res.taln_offsets = std::move(task.taln_offsets);
            res.qaln_data = std::move(task.qaln_data);
            res.taln_data = std::move(task.taln_data);

            if (!config.skip_alignments) {
                if (res.qaln_offsets.size() == rg_rows) {
                    res.qaln_offsets.push_back(static_cast<uint32_t>(res.qaln_data.size()));
                }
                if (res.taln_offsets.size() == rg_rows) {
                    res.taln_offsets.push_back(static_cast<uint32_t>(res.taln_data.size()));
                }
                if (res.qaln_offsets.size() != rg_rows + 1 ||
                    res.taln_offsets.size() != rg_rows + 1) {
                    throw std::runtime_error("Invalid alignment offset vector size in row-group task");
                }
            }

            std::vector<uint32_t> perm(rg_rows);
            std::iota(perm.begin(), perm.end(), 0);
            std::sort(perm.begin(), perm.end(), [&](uint32_t a, uint32_t b) {
                if (res.read_idx[a] != res.read_idx[b]) return res.read_idx[a] < res.read_idx[b];
                return res.bit_score[a] > res.bit_score[b];
            });

            std::vector<uint32_t> tmp_u32(rg_rows);
            std::vector<float>    tmp_f32(rg_rows);
            std::vector<uint16_t> tmp_u16(rg_rows);
            auto apply_perm_u32 = [&](std::vector<uint32_t>& vec) {
                for (uint32_t i = 0; i < rg_rows; ++i) tmp_u32[i] = vec[perm[i]];
                std::copy_n(tmp_u32.data(), rg_rows, vec.data());
            };
            auto apply_perm_f32 = [&](std::vector<float>& vec) {
                for (uint32_t i = 0; i < rg_rows; ++i) tmp_f32[i] = vec[perm[i]];
                std::copy_n(tmp_f32.data(), rg_rows, vec.data());
            };
            auto apply_perm_u16 = [&](std::vector<uint16_t>& vec) {
                for (uint32_t i = 0; i < rg_rows; ++i) tmp_u16[i] = vec[perm[i]];
                std::copy_n(tmp_u16.data(), rg_rows, vec.data());
            };
            apply_perm_u32(res.read_idx);
            apply_perm_u32(res.ref_idx);
            apply_perm_f32(res.bit_score);
            apply_perm_f32(res.damage_score);
            apply_perm_f32(res.evalue_log10);
            apply_perm_u16(res.dmg_k);
            apply_perm_u16(res.dmg_m);
            apply_perm_f32(res.dmg_ll_a);
            apply_perm_f32(res.dmg_ll_m);
            apply_perm_u16(res.identity_q);
            apply_perm_u16(res.aln_len);
            apply_perm_u16(res.qstart);
            apply_perm_u16(res.qend);
            apply_perm_u16(res.tstart);
            apply_perm_u16(res.tend);
            apply_perm_u16(res.tlen);
            apply_perm_u16(res.qlen);
            apply_perm_u16(res.mismatch);
            apply_perm_u16(res.gapopen);

            if (!config.skip_alignments) {
                std::vector<uint32_t> perm_qaln_offsets(rg_rows + 1, 0);
                std::vector<uint32_t> perm_taln_offsets(rg_rows + 1, 0);
                std::vector<char> perm_qaln_data;
                std::vector<char> perm_taln_data;
                perm_qaln_data.reserve(res.qaln_data.size());
                perm_taln_data.reserve(res.taln_data.size());
                for (uint32_t i = 0; i < rg_rows; ++i) {
                    uint32_t old_idx = perm[i];
                    uint32_t q_start = res.qaln_offsets[old_idx];
                    uint32_t q_end   = res.qaln_offsets[old_idx + 1];
                    uint32_t t_start = res.taln_offsets[old_idx];
                    uint32_t t_end   = res.taln_offsets[old_idx + 1];
                    perm_qaln_offsets[i] = static_cast<uint32_t>(perm_qaln_data.size());
                    perm_qaln_data.insert(perm_qaln_data.end(),
                        res.qaln_data.begin() + q_start, res.qaln_data.begin() + q_end);
                    perm_taln_offsets[i] = static_cast<uint32_t>(perm_taln_data.size());
                    perm_taln_data.insert(perm_taln_data.end(),
                        res.taln_data.begin() + t_start, res.taln_data.begin() + t_end);
                }
                perm_qaln_offsets[rg_rows] = static_cast<uint32_t>(perm_qaln_data.size());
                perm_taln_offsets[rg_rows] = static_cast<uint32_t>(perm_taln_data.size());
                res.qaln_offsets.swap(perm_qaln_offsets);
                res.taln_offsets.swap(perm_taln_offsets);
                res.qaln_data.swap(perm_qaln_data);
                res.taln_data.swap(perm_taln_data);
            }

            if (rg_rows > 0) {
                res.first_read_idx = res.read_idx[0];
                res.last_read_idx  = res.read_idx[rg_rows - 1];
                uint32_t cur = res.read_idx[0];
                uint32_t run = 1;
                for (uint32_t i = 1; i < rg_rows; ++i) {
                    if (res.read_idx[i] == cur) {
                        ++run;
                    } else {
                        res.read_runs.emplace_back(cur, run);
                        cur = res.read_idx[i];
                        run = 1;
                    }
                }
                res.read_runs.emplace_back(cur, run);
            }
            return res;
        };

        // Write a completed RowGroupResult to the output file and update metadata.
        auto flush_result = [&](const RowGroupResult& rg_result) {
            const uint32_t rg_rows = rg_result.rg_rows;
            const size_t rg_i = rg_result.rg_idx;

            for (uint32_t i = 0; i < rg_rows; ++i) {
                if (!has_prev) {
                    has_prev = true;
                } else if (globally_sorted) {
                    if (rg_result.read_idx[i] < prev_read_idx_s ||
                        (rg_result.read_idx[i] == prev_read_idx_s &&
                         rg_result.bit_score[i] > prev_bit_score_s)) {
                        globally_sorted = false;
                    }
                }
                prev_read_idx_s  = rg_result.read_idx[i];
                prev_bit_score_s = rg_result.bit_score[i];
            }

            RowGroupMeta rg_meta{};
            rg_meta.num_rows       = rg_rows;
            rg_meta.first_read_idx = rg_result.first_read_idx;

            RowGroupReadIndex rg_idx_entry{};
            rg_idx_entry.first_read_idx    = rg_result.first_read_idx;
            rg_idx_entry.last_read_idx     = rg_result.last_read_idx;
            rg_idx_entry.start_row         = global_row_cursor;
            rg_idx_entry.end_row_exclusive = global_row_cursor + rg_rows;
            global_row_cursor += rg_rows;

            auto write_chunk = [&](ColumnID col, const void* data, size_t elem_size) {
                auto& chunk = rg_meta.columns[static_cast<size_t>(col)];
                chunk.offset            = out.tellp();
                chunk.uncompressed_size = static_cast<uint32_t>(elem_size * rg_rows);
                chunk.compressed_size   = chunk.uncompressed_size;
                chunk.codec             = Codec::NONE;
                out.write(static_cast<const char*>(data), chunk.uncompressed_size);
            };
            auto write_alignment_chunk = [&](ColumnID col,
                                             const std::vector<uint32_t>& offsets,
                                             const std::vector<char>& data) {
                auto& chunk = rg_meta.columns[static_cast<size_t>(col)];
                chunk.offset = out.tellp();
                const size_t offsets_bytes = static_cast<size_t>(rg_rows + 1) * sizeof(uint32_t);
                const size_t raw_bytes     = offsets_bytes + data.size();
                chunk.uncompressed_size    = static_cast<uint32_t>(raw_bytes);
#ifdef HAVE_ZSTD
                if (config.compress && raw_bytes > 0) {
                    std::vector<char> raw(raw_bytes);
                    std::memcpy(raw.data(), offsets.data(), offsets_bytes);
                    if (!data.empty())
                        std::memcpy(raw.data() + offsets_bytes, data.data(), data.size());
                    auto compressed = compress_zstd(raw.data(), raw.size(), config.compression_level);
                    if (compressed.size() < raw.size()) {
                        chunk.compressed_size = static_cast<uint32_t>(compressed.size());
                        chunk.codec = Codec::ZSTD;
                        out.write(compressed.data(), static_cast<std::streamsize>(compressed.size()));
                        return;
                    }
                }
#endif
                chunk.compressed_size = chunk.uncompressed_size;
                chunk.codec = Codec::NONE;
                out.write(reinterpret_cast<const char*>(offsets.data()),
                          static_cast<std::streamsize>(offsets_bytes));
                if (!data.empty())
                    out.write(data.data(), static_cast<std::streamsize>(data.size()));
            };

            write_chunk(ColumnID::READ_IDX,    rg_result.read_idx.data(),    sizeof(uint32_t));
            write_chunk(ColumnID::REF_IDX,     rg_result.ref_idx.data(),     sizeof(uint32_t));
            write_chunk(ColumnID::BIT_SCORE,   rg_result.bit_score.data(),   sizeof(float));
            write_chunk(ColumnID::DAMAGE_SCORE,rg_result.damage_score.data(),sizeof(float));
            write_chunk(ColumnID::EVALUE_LOG10,rg_result.evalue_log10.data(),sizeof(float));
            write_chunk(ColumnID::DMG_K,       rg_result.dmg_k.data(),       sizeof(uint16_t));
            write_chunk(ColumnID::DMG_M,       rg_result.dmg_m.data(),       sizeof(uint16_t));
            write_chunk(ColumnID::DMG_LL_A,    rg_result.dmg_ll_a.data(),    sizeof(float));
            write_chunk(ColumnID::DMG_LL_M,    rg_result.dmg_ll_m.data(),    sizeof(float));
            write_chunk(ColumnID::IDENTITY_Q,  rg_result.identity_q.data(),  sizeof(uint16_t));
            write_chunk(ColumnID::ALN_LEN,     rg_result.aln_len.data(),     sizeof(uint16_t));
            write_chunk(ColumnID::QSTART,      rg_result.qstart.data(),      sizeof(uint16_t));
            write_chunk(ColumnID::QEND,        rg_result.qend.data(),        sizeof(uint16_t));
            write_chunk(ColumnID::TSTART,      rg_result.tstart.data(),      sizeof(uint16_t));
            write_chunk(ColumnID::TEND,        rg_result.tend.data(),        sizeof(uint16_t));
            write_chunk(ColumnID::TLEN,        rg_result.tlen.data(),        sizeof(uint16_t));
            write_chunk(ColumnID::QLEN,        rg_result.qlen.data(),        sizeof(uint16_t));
            write_chunk(ColumnID::MISMATCH,    rg_result.mismatch.data(),    sizeof(uint16_t));
            write_chunk(ColumnID::GAPOPEN,     rg_result.gapopen.data(),     sizeof(uint16_t));

            if (!config.skip_alignments) {
                write_alignment_chunk(ColumnID::QALN, rg_result.qaln_offsets, rg_result.qaln_data);
                write_alignment_chunk(ColumnID::TALN, rg_result.taln_offsets, rg_result.taln_data);
            }

            row_groups.push_back(rg_meta);
            row_group_read_index_vec.push_back(rg_idx_entry);

            if (config.verbose && (rg_i + 1) % 10 == 0) {
                std::cerr << "\rPass: " << global_row_cursor << " rows, "
                          << (reader.bytes_read() / (1024*1024)) << " MB read...";
            }
        };

        // Main streaming loop: read one line at a time, build dicts inline,
        // accumulate into RowGroupTask, flush synchronously when full.
        size_t cur_rg_idx = 0;
        RowGroupTask current_task = make_row_group_task(0);

        while (true) {
            std::string_view line = reader.next_line();
            if (line.empty()) break;
            ++n_lines;

            const char* p = line.data();
            const char* line_end = p + line.size();
            const char* fields[16];
            size_t lens[16];
            const char* tab1 = simd_find_tab(p, line_end);
            if (tab1 >= line_end) continue;
            fields[0] = p;
            lens[0] = static_cast<size_t>(tab1 - p);

            const char* fp = tab1 + 1;
            for (int f = 1; f < 16; ++f) {
                fields[f] = fp;
                const char* tab = (f == 15) ? line_end : simd_find_tab(fp, line_end);
                lens[f] = static_cast<size_t>(tab - fp);
                fp = (tab < line_end) ? tab + 1 : line_end;
            }
            const std::string_view q(fields[0], lens[0]);
            const std::string_view t(fields[1], lens[1]);

            // Read dict: queries are contiguous in MMseqs2 output; cache miss => new read.
            uint32_t rid;
            if (q == last_read_name) {
                rid = last_read_idx_cached;
            } else {
                rid = static_cast<uint32_t>(read_names_vec.size());
                read_names_vec.emplace_back(q);
                // Lazy damage score lookup — O(1) hash lookup per new read
                uint8_t dmg_q = 0;
                if (damage_reader) {
                    if (const AgdRecord* rec = damage_reader->find(read_names_vec.back()))
                        dmg_q = rec->p_damaged_q;
                }
                read_damage_q.push_back(dmg_q);
                read_counts.push_back(0);
                last_read_name = read_names_vec.back();
                last_read_idx_cached = rid;
            }
            read_counts[rid]++;

            // Ref dict: first-appearance order
            uint32_t tid;
            if (const uint32_t* tid_ptr = ref_lookup.find(t, ref_names_vec)) {
                tid = *tid_ptr;
            } else {
                tid = static_cast<uint32_t>(ref_names_vec.size());
                ref_names_vec.emplace_back(t);
                ref_lookup.emplace(tid, tid, ref_names_vec);
            }

            const float dmg_score = static_cast<float>(read_damage_q[rid]) * (1.0f / 255.0f);
            float fident       = fast_stof({fields[2], lens[2]});
            uint32_t alnlen_v  = fast_stou({fields[3], lens[3]});
            uint32_t mismatch_v = fast_stou({fields[4], lens[4]});
            uint32_t gapopen_v  = fast_stou({fields[5], lens[5]});
            uint32_t qstart_v   = fast_stou({fields[6], lens[6]});
            uint32_t qend_v     = fast_stou({fields[7], lens[7]});
            uint32_t tstart_v   = fast_stou({fields[8], lens[8]});
            uint32_t tend_v     = fast_stou({fields[9], lens[9]});
            float evalue        = fast_stof({fields[10], lens[10]});
            float bits          = fast_stof({fields[11], lens[11]});
            uint32_t qlen_v     = fast_stou({fields[12], lens[12]});
            uint32_t tlen_v     = fast_stou({fields[13], lens[13]});

            const uint16_t qstart_q = static_cast<uint16_t>(std::min(qstart_v, 65535u));
            const uint16_t qlen_q   = static_cast<uint16_t>(std::min(qlen_v, 65535u));

            current_task.read_idx.push_back(rid);
            current_task.ref_idx.push_back(tid);
            current_task.bit_score.push_back(bits);
            current_task.damage_score.push_back(dmg_score);
            current_task.evalue_log10.push_back((evalue > 0) ? std::log10(evalue) : -300.0f);
            current_task.identity_q.push_back(
                static_cast<uint16_t>(std::clamp(fident, 0.0f, 1.0f) * 65535.0f));
            current_task.aln_len.push_back(static_cast<uint16_t>(std::min(alnlen_v, 65535u)));
            current_task.qstart.push_back(qstart_q);
            current_task.qend.push_back(static_cast<uint16_t>(std::min(qend_v, 65535u)));
            current_task.tstart.push_back(static_cast<uint16_t>(std::min(tstart_v, 65535u)));
            current_task.tend.push_back(static_cast<uint16_t>(std::min(tend_v, 65535u)));
            current_task.tlen.push_back(static_cast<uint16_t>(std::min(tlen_v, 65535u)));
            current_task.qlen.push_back(qlen_q);
            current_task.mismatch.push_back(static_cast<uint16_t>(std::min(mismatch_v, 65535u)));
            current_task.gapopen.push_back(static_cast<uint16_t>(std::min(gapopen_v, 65535u)));

            if (!config.skip_alignments) {
                const char* qaln_start = fields[14];
                const char* qaln_end = fields[14] + lens[14];
                const char* taln_start = fields[15];
                const char* taln_end = fields[15] + lens[15];

                const std::string_view qaln_sv(qaln_start, lens[14]);
                const std::string_view taln_sv(taln_start, lens[15]);
                const auto dmg = compute_alignment_damage_evidence(
                    q, qaln_sv, taln_sv, qstart_q, qlen_q, header_d_max, header_lambda);

                current_task.dmg_k.push_back(dmg.k_hits);
                current_task.dmg_m.push_back(dmg.m_sites);
                current_task.dmg_ll_a.push_back(dmg.ll_ancient);
                current_task.dmg_ll_m.push_back(dmg.ll_modern);

                current_task.qaln_offsets.push_back(static_cast<uint32_t>(current_task.qaln_data.size()));
                current_task.qaln_data.insert(current_task.qaln_data.end(), qaln_start, qaln_end);
                current_task.taln_offsets.push_back(static_cast<uint32_t>(current_task.taln_data.size()));
                current_task.taln_data.insert(current_task.taln_data.end(), taln_start, taln_end);
            } else {
                current_task.dmg_k.push_back(0);
                current_task.dmg_m.push_back(0);
                current_task.dmg_ll_a.push_back(0.0f);
                current_task.dmg_ll_m.push_back(0.0f);
            }

            if (current_task.read_idx.size() == rg_size) {
                if (!config.skip_alignments) {
                    current_task.qaln_offsets.push_back(static_cast<uint32_t>(current_task.qaln_data.size()));
                    current_task.taln_offsets.push_back(static_cast<uint32_t>(current_task.taln_data.size()));
                }
                auto result = process_row_group(std::move(current_task));
                flush_result(result);
                current_task = make_row_group_task(++cur_rg_idx);
            }
        }

        // Flush final (possibly partial) row group
        if (!current_task.read_idx.empty()) {
            if (!config.skip_alignments) {
                current_task.qaln_offsets.push_back(static_cast<uint32_t>(current_task.qaln_data.size()));
                current_task.taln_offsets.push_back(static_cast<uint32_t>(current_task.taln_data.size()));
            }
            auto result = process_row_group(std::move(current_task));
            flush_result(result);
        }

        auto pass_end = std::chrono::high_resolution_clock::now();
        auto pass_ms = std::chrono::duration_cast<std::chrono::milliseconds>(
            pass_end - total_start).count();
        if (config.verbose) {
            std::cerr << "\rPass: " << n_lines << " rows, "
                      << read_names_vec.size() << " reads, "
                      << ref_names_vec.size() << " refs"
                      << " (" << dart::log_utils::format_duration_ms(pass_ms) << ")\n";
            if (damage_reader) {
                std::cerr << "Embedded per-read damage scores from "
                          << config.damage_index_path << " into EMI DAMAGE_SCORE\n";
            }
        }

        // =====================================================================
        // Trailer: rg_read_index, csr_offsets, string_dicts (v6), row_group_dir
        // =====================================================================

        if (!config.skip_alignments)    header.flags |= EMI_FLAG_HAS_ALIGNMENT_STRINGS;
        header.flags |= EMI_FLAG_HAS_ROW_GROUP_READ_INDEX;
        header.flags |= EMI_FLAG_64BIT_STRING_DICT;
        if (globally_sorted) {
            header.flags |= EMI_FLAG_SORTED_BY_READ_THEN_SCORE;
        } else if (config.verbose) {
            std::cerr << "Warning: EMI is not globally sorted by (read, bit_score). "
                         "Input is likely not grouped by query.\n";
        }

        // Row group read index
        header.row_group_read_index_offset     = out.tellp();
        header.row_group_read_index_count      = static_cast<uint32_t>(row_group_read_index_vec.size());
        header.row_group_read_index_entry_size = sizeof(RowGroupReadIndex);
        out.write(reinterpret_cast<const char*>(row_group_read_index_vec.data()),
                  row_group_read_index_vec.size() * sizeof(RowGroupReadIndex));

        // CSR offsets (prefix sum of per-read alignment counts)
        header.csr_offsets_offset = out.tellp();
        const uint32_t num_reads  = static_cast<uint32_t>(read_counts.size());
        std::vector<uint64_t> csr_offsets(num_reads + 1);
        csr_offsets[0] = 0;
        for (uint32_t r = 0; r < num_reads; ++r)
            csr_offsets[r + 1] = csr_offsets[r] + read_counts[r];
        out.write(reinterpret_cast<const char*>(csr_offsets.data()),
                  (num_reads + 1) * sizeof(uint64_t));

        // String dictionaries — v6 format: uint64_t offsets + uint64_t size
        header.string_dict_offset = out.tellp();

        // Read names
        {
            std::vector<uint64_t> offsets;
            std::vector<char> blob;
            offsets.reserve(read_names_vec.size());
            for (const auto& s : read_names_vec) {
                offsets.push_back(static_cast<uint64_t>(blob.size()));
                blob.insert(blob.end(), s.begin(), s.end());
                blob.push_back('\0');
            }
            const uint32_t count = static_cast<uint32_t>(offsets.size());
            const uint64_t sz    = blob.size();
            out.write(reinterpret_cast<const char*>(&count), sizeof(count));
            out.write(reinterpret_cast<const char*>(offsets.data()),
                      offsets.size() * sizeof(uint64_t));
            out.write(reinterpret_cast<const char*>(&sz), sizeof(sz));
            out.write(blob.data(), blob.size());
        }

        // Ref names
        {
            std::vector<uint64_t> offsets;
            std::vector<char> blob;
            offsets.reserve(ref_names_vec.size());
            for (const auto& s : ref_names_vec) {
                offsets.push_back(static_cast<uint64_t>(blob.size()));
                blob.insert(blob.end(), s.begin(), s.end());
                blob.push_back('\0');
            }
            const uint32_t count = static_cast<uint32_t>(offsets.size());
            const uint64_t sz    = blob.size();
            out.write(reinterpret_cast<const char*>(&count), sizeof(count));
            out.write(reinterpret_cast<const char*>(offsets.data()),
                      offsets.size() * sizeof(uint64_t));
            out.write(reinterpret_cast<const char*>(&sz), sizeof(sz));
            out.write(blob.data(), blob.size());
        }

        // Row group directory
        header.row_group_dir_offset = out.tellp();
        out.write(reinterpret_cast<const char*>(row_groups.data()),
                  row_groups.size() * sizeof(RowGroupMeta));

        // Finalize header with actual values, seek back and overwrite
        header.num_alignments = n_lines;
        header.num_reads      = num_reads;
        header.num_refs       = static_cast<uint32_t>(ref_names_vec.size());
        header.num_row_groups = static_cast<uint32_t>(row_groups.size());
        header.row_group_size = static_cast<uint32_t>(rg_size);
        out.seekp(0);
        out.write(reinterpret_cast<const char*>(&header), sizeof(header));

        auto total_end = std::chrono::high_resolution_clock::now();
        auto total_ms = std::chrono::duration_cast<std::chrono::milliseconds>(
            total_end - total_start).count();
        if (config.verbose) {
            std::cerr << "Total: "
                      << (n_lines * 1000 / std::max(1LL, static_cast<long long>(total_ms)))
                      << " rows/sec, " << dart::log_utils::format_duration_ms(total_ms) << "\n";
        }

        return n_lines;
    }

};  // class FastColumnarWriter

}  // namespace dart
