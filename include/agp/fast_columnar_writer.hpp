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
 * 6. Two-pass: dictionary build + columnar write
 *
 * Memory bounds:
 * - Input buffer: 128 MB (double-buffered chunks)
 * - Dictionary: configurable cap (default 4 GB, spills to disk)
 * - Row-group buffers: ~50 MB (fixed)
 * - Total: ~4-5 GB regardless of file size
 */

#include "agp/columnar_index.hpp"
#include "agp/fast_tsv.hpp"

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

#ifdef HAVE_ZSTD
#include <zstd.h>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef USE_AVX2
#include <immintrin.h>
#endif

namespace agp {

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

            // Advise sequential access
            posix_fadvise(fd_, 0, 0, POSIX_FADV_SEQUENTIAL);
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
        auto it = in_memory_.find(std::string(sv));
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
    std::pair<std::vector<std::string>, std::unordered_map<std::string, uint32_t>> finalize() {
        std::vector<std::string> names;
        std::unordered_map<std::string, uint32_t> lookup;

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

    std::unordered_map<std::string, uint32_t> in_memory_;
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

}  // namespace agp (compression utilities)

namespace agp {

/**
 * @brief Streaming gzip reader using ring buffer
 *
 * Uses zlib for decompression with a chunked ring buffer,
 * allowing overlap between decompression and parsing.
 */
class GzipStreamReader {
public:
    static constexpr size_t CHUNK_SIZE = 16 * 1024 * 1024;  // 16MB chunks

    explicit GzipStreamReader(const std::string& path) {
        gz_ = gzopen(path.c_str(), "rb");
        if (!gz_) return;

        // Set large buffer for better throughput
        gzbuffer(gz_, 1024 * 1024);  // 1MB internal buffer

        // Pre-allocate ring buffer
        buffer_.resize(CHUNK_SIZE * 2);
    }

    ~GzipStreamReader() {
        if (gz_) gzclose(gz_);
    }

    bool valid() const { return gz_ != nullptr; }

    // Read next chunk, returns (data, size) - data valid until next read
    std::pair<const char*, size_t> read_chunk() {
        if (!gz_ || gzeof(gz_)) return {nullptr, 0};

        int bytes = gzread(gz_, buffer_.data(), static_cast<unsigned>(CHUNK_SIZE));
        if (bytes <= 0) return {nullptr, 0};

        return {buffer_.data(), static_cast<size_t>(bytes)};
    }

    // Read entire file into memory (for smaller files)
    std::vector<char> read_all() {
        std::vector<char> result;
        result.reserve(100 * 1024 * 1024);  // Assume ~100MB uncompressed

        while (true) {
            auto [data, size] = read_chunk();
            if (!data || size == 0) break;
            result.insert(result.end(), data, data + size);
        }

        return result;
    }

private:
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
 * Two-pass streaming algorithm:
 * Pass 1: Stream file to count rows and build dictionaries (with spilling)
 * Pass 2: Stream file again to parse and write row-groups
 *
 * Memory bounded to ~4-5GB regardless of file size.
 */
class FastColumnarWriter {
public:
    struct Config {
        float d_max = 0.0f;
        float lambda = 0.3f;
        bool verbose = false;
        bool compress = true;  // Use ZSTD compression (if available)
        int compression_level = 3;  // ZSTD level (1-19, 3 = good balance)
        bool skip_alignments = false;  // Skip qaln/taln strings (faster, EM doesn't need them)
        size_t memory_limit = 16ULL * 1024 * 1024 * 1024;  // 16GB default memory cap
        bool force_streaming = false;  // Force streaming even for small files
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
     * @param output_path Output .emi2 file
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
        bool is_gz = (tsv_path.size() >= 3 && tsv_path.substr(tsv_path.size() - 3) == ".gz");
        size_t estimated_size = is_gz ? file_size * 5 : file_size;

        // Use mmap if file fits in memory limit, otherwise stream
        if (!config.force_streaming && estimated_size < config.memory_limit) {
            if (config.verbose) {
                std::cerr << "Using mmap mode (file " << (estimated_size / (1024*1024))
                          << " MB < limit " << (config.memory_limit / (1024*1024*1024)) << " GB)\n";
            }
            return convert_mmap(tsv_path, output_path, config);
        } else {
            if (config.verbose) {
                std::cerr << "Using streaming mode (file " << (estimated_size / (1024*1024))
                          << " MB, limit " << (config.memory_limit / (1024*1024*1024)) << " GB)\n";
            }
            return convert_streaming(tsv_path, output_path, config);
        }
    }

    /**
     * @brief Fast mmap-based conversion (for files that fit in memory)
     */
    static uint64_t convert_mmap(
        const std::string& tsv_path,
        const std::string& output_path,
        const Config& config)
    {
        auto total_start = std::chrono::high_resolution_clock::now();

        // Handle gzip vs plain file
        std::vector<char> gz_buffer;
        const char* data = nullptr;
        size_t size = 0;

        if (is_gzip(tsv_path)) {
            if (config.verbose) {
                std::cerr << "Detected gzip compression, decompressing...\n";
            }
            GzipStreamReader gz(tsv_path);
            if (!gz.valid()) {
                throw std::runtime_error("Cannot open gzip TSV: " + tsv_path);
            }
            gz_buffer = gz.read_all();
            data = gz_buffer.data();
            size = gz_buffer.size();
            if (config.verbose) {
                std::cerr << "Decompressed: " << (size / (1024*1024)) << " MB\n";
            }
        }

        // For plain files, use mmap
        std::unique_ptr<MappedTSV> tsv_holder;
        if (!is_gzip(tsv_path)) {
            tsv_holder = std::make_unique<MappedTSV>(tsv_path.c_str());
            if (!tsv_holder->valid()) {
                throw std::runtime_error("Cannot open TSV: " + tsv_path);
            }
            data = tsv_holder->data();
            size = tsv_holder->size();
        }

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
            std::cerr << "Lines: " << n_lines << " (scan: " << scan_ms << " ms), building dictionary...\n";
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

        // Merge into global dictionaries
        constexpr size_t READ_SHARDS = 1024;
        constexpr size_t REF_SHARDS = 256;
        auto read_dict = build_global_dict_from_locals(local_read_uniques, READ_SHARDS);
        auto ref_dict = build_global_dict_from_locals(local_ref_uniques, REF_SHARDS);

        std::vector<std::string_view> read_names = std::move(read_dict.names);
        std::vector<std::string_view> ref_names = std::move(ref_dict.names);

        auto dict_end = std::chrono::high_resolution_clock::now();
        auto dict_ms = std::chrono::duration_cast<std::chrono::milliseconds>(dict_end - pass_start).count();
        if (config.verbose) {
            std::cerr << "  Dictionary: " << dict_ms << " ms (" << read_names.size() << " reads, "
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
        header.d_max = config.d_max;
        header.lambda = config.lambda;

        out.write(reinterpret_cast<const char*>(&header), sizeof(header));

        header.row_group_dir_offset = out.tellp();
        std::vector<RowGroupMeta> row_groups(num_row_groups);
        out.write(reinterpret_cast<const char*>(row_groups.data()), row_groups.size() * sizeof(RowGroupMeta));

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

        // Row-group buffers
        std::vector<uint32_t> rg_read_idx(rg_size);
        std::vector<uint32_t> rg_ref_idx(rg_size);
        std::vector<float> rg_bit_score(rg_size);
        std::vector<float> rg_damage_score(rg_size);
        std::vector<float> rg_evalue_log10(rg_size);
        std::vector<uint16_t> rg_identity_q(rg_size);
        std::vector<uint16_t> rg_aln_len(rg_size);
        std::vector<uint16_t> rg_tstart(rg_size);
        std::vector<uint16_t> rg_tend(rg_size);
        std::vector<uint16_t> rg_tlen(rg_size);
        std::vector<uint16_t> rg_qlen(rg_size);
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
                rg_read_idx[local_idx] = dict_lookup(read_dict, q);
                rg_ref_idx[local_idx] = dict_lookup(ref_dict, t);

                float fident = fast_stof({fields[2], lens[2]});
                uint32_t alnlen = fast_stou({fields[3], lens[3]});
                uint32_t tstart_v = fast_stou({fields[8], lens[8]});
                uint32_t tend_v = fast_stou({fields[9], lens[9]});
                float evalue = fast_stof({fields[10], lens[10]});
                float bits = fast_stof({fields[11], lens[11]});
                uint32_t qlen_v = fast_stou({fields[12], lens[12]});
                uint32_t tlen_v = fast_stou({fields[13], lens[13]});

                rg_bit_score[local_idx] = bits;
                rg_damage_score[local_idx] = 0.0f;
                rg_evalue_log10[local_idx] = (evalue > 0) ? std::log10(evalue) : -300.0f;
                rg_identity_q[local_idx] = static_cast<uint16_t>(std::clamp(fident, 0.0f, 1.0f) * 65535.0f);
                rg_aln_len[local_idx] = static_cast<uint16_t>(std::min(alnlen, 65535u));
                rg_tstart[local_idx] = static_cast<uint16_t>(std::min(tstart_v, 65535u));
                rg_tend[local_idx] = static_cast<uint16_t>(std::min(tend_v, 65535u));
                rg_tlen[local_idx] = static_cast<uint16_t>(std::min(tlen_v, 65535u));
                rg_qlen[local_idx] = static_cast<uint16_t>(std::min(qlen_v, 65535u));
            }

            // Parse string columns (sequential)
            if (!config.skip_alignments) {
                for (uint32_t local_idx = 0; local_idx < rg_rows; ++local_idx) {
                    const char* p = rg_line_starts[local_idx];
                    const char* line_end = rg_line_starts[local_idx + 1];
                    if (line_end > p && *(line_end - 1) == '\n') --line_end;

                    const char* fp = p;
                    for (int f = 0; f < 14 && fp < line_end; ++f) {
                        fp = simd_find_tab(fp, line_end) + 1;
                    }
                    const char* qaln_start = fp;
                    const char* qaln_end = simd_find_tab(fp, line_end);
                    fp = qaln_end + 1;
                    const char* taln_start = fp;
                    const char* taln_end = simd_find_tab(fp, line_end);

                    rg_qaln_offsets[local_idx] = static_cast<uint32_t>(rg_qaln_data.size());
                    rg_qaln_data.insert(rg_qaln_data.end(), qaln_start, qaln_end);
                    rg_taln_offsets[local_idx] = static_cast<uint32_t>(rg_taln_data.size());
                    rg_taln_data.insert(rg_taln_data.end(), taln_start, taln_end);
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

            // Apply permutation to all columns
            auto apply_perm = [&](auto& vec) {
                using T = typename std::remove_reference<decltype(vec)>::type::value_type;
                std::vector<T> temp(rg_rows);
                for (uint32_t i = 0; i < rg_rows; ++i) temp[i] = vec[perm[i]];
                std::copy(temp.begin(), temp.end(), vec.begin());
            };
            apply_perm(rg_read_idx);
            apply_perm(rg_ref_idx);
            apply_perm(rg_bit_score);
            apply_perm(rg_damage_score);
            apply_perm(rg_evalue_log10);
            apply_perm(rg_identity_q);
            apply_perm(rg_aln_len);
            apply_perm(rg_tstart);
            apply_perm(rg_tend);
            apply_perm(rg_tlen);
            apply_perm(rg_qlen);

            // Rebuild qaln/taln with permuted order
            if (!config.skip_alignments) {
                std::vector<std::string> qaln_strs(rg_rows), taln_strs(rg_rows);
                for (uint32_t i = 0; i < rg_rows; ++i) {
                    uint32_t old_i = perm[i];
                    uint32_t q_start = rg_qaln_offsets[old_i];
                    uint32_t q_end = rg_qaln_offsets[old_i + 1];
                    uint32_t t_start = rg_taln_offsets[old_i];
                    uint32_t t_end = rg_taln_offsets[old_i + 1];
                    qaln_strs[i] = std::string(rg_qaln_data.begin() + q_start, rg_qaln_data.begin() + q_end);
                    taln_strs[i] = std::string(rg_taln_data.begin() + t_start, rg_taln_data.begin() + t_end);
                }
                rg_qaln_data.clear();
                rg_taln_data.clear();
                for (uint32_t i = 0; i < rg_rows; ++i) {
                    rg_qaln_offsets[i] = static_cast<uint32_t>(rg_qaln_data.size());
                    rg_qaln_data.insert(rg_qaln_data.end(), qaln_strs[i].begin(), qaln_strs[i].end());
                    rg_taln_offsets[i] = static_cast<uint32_t>(rg_taln_data.size());
                    rg_taln_data.insert(rg_taln_data.end(), taln_strs[i].begin(), taln_strs[i].end());
                }
                rg_qaln_offsets[rg_rows] = static_cast<uint32_t>(rg_qaln_data.size());
                rg_taln_offsets[rg_rows] = static_cast<uint32_t>(rg_taln_data.size());
            }

            // Accumulate CSR counts (now sorted)
            for (uint32_t i = 0; i < rg_rows; ++i) {
                read_counts[rg_read_idx[i]]++;
            }

            // Write row-group
            row_groups[rg].num_rows = rg_rows;
            row_groups[rg].first_read_idx = rg_read_idx[0];

            auto write_chunk = [&](ColumnID col, const void* chunk_data, size_t elem_size) {
                auto& chunk = row_groups[rg].columns[static_cast<size_t>(col)];
                chunk.offset = out.tellp();
                chunk.uncompressed_size = static_cast<uint32_t>(elem_size * rg_rows);
                chunk.compressed_size = chunk.uncompressed_size;
                chunk.codec = Codec::NONE;
                out.write(static_cast<const char*>(chunk_data), chunk.uncompressed_size);
            };

            write_chunk(ColumnID::READ_IDX, rg_read_idx.data(), sizeof(uint32_t));
            write_chunk(ColumnID::REF_IDX, rg_ref_idx.data(), sizeof(uint32_t));
            write_chunk(ColumnID::BIT_SCORE, rg_bit_score.data(), sizeof(float));
            write_chunk(ColumnID::DAMAGE_SCORE, rg_damage_score.data(), sizeof(float));
            write_chunk(ColumnID::EVALUE_LOG10, rg_evalue_log10.data(), sizeof(float));
            write_chunk(ColumnID::IDENTITY_Q, rg_identity_q.data(), sizeof(uint16_t));
            write_chunk(ColumnID::ALN_LEN, rg_aln_len.data(), sizeof(uint16_t));
            write_chunk(ColumnID::TSTART, rg_tstart.data(), sizeof(uint16_t));
            write_chunk(ColumnID::TEND, rg_tend.data(), sizeof(uint16_t));
            write_chunk(ColumnID::TLEN, rg_tlen.data(), sizeof(uint16_t));
            write_chunk(ColumnID::QLEN, rg_qlen.data(), sizeof(uint16_t));

            if (!config.skip_alignments) {
                auto& qaln_chunk = row_groups[rg].columns[static_cast<size_t>(ColumnID::QALN)];
                qaln_chunk.offset = out.tellp();
                out.write(reinterpret_cast<const char*>(rg_qaln_offsets.data()), (rg_rows + 1) * sizeof(uint32_t));
                out.write(rg_qaln_data.data(), rg_qaln_data.size());
                qaln_chunk.uncompressed_size = static_cast<uint32_t>((rg_rows + 1) * sizeof(uint32_t) + rg_qaln_data.size());
                qaln_chunk.compressed_size = qaln_chunk.uncompressed_size;
                qaln_chunk.codec = Codec::NONE;

                auto& taln_chunk = row_groups[rg].columns[static_cast<size_t>(ColumnID::TALN)];
                taln_chunk.offset = out.tellp();
                out.write(reinterpret_cast<const char*>(rg_taln_offsets.data()), (rg_rows + 1) * sizeof(uint32_t));
                out.write(rg_taln_data.data(), rg_taln_data.size());
                taln_chunk.uncompressed_size = static_cast<uint32_t>((rg_rows + 1) * sizeof(uint32_t) + rg_taln_data.size());
                taln_chunk.compressed_size = taln_chunk.uncompressed_size;
                taln_chunk.codec = Codec::NONE;
            }
        }

        // Cleanup dictionaries
        for (auto& shard : read_dict.shard_maps) shard.clear_release();
        for (auto& shard : ref_dict.shard_maps) shard.clear_release();

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
                      << header.num_refs << " refs ("
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

        if (config.verbose) {
            std::cerr << "Streaming mode (memory cap: "
                      << (config.memory_limit / (1024*1024*1024)) << " GB)\n";
        }

        // =====================================================================
        // Pass 1: Stream through file, count lines, build dictionaries
        // =====================================================================
        auto pass1_start = std::chrono::high_resolution_clock::now();

        ChunkedFileReader reader(tsv_path);
        if (!reader.valid()) {
            throw std::runtime_error("Cannot open TSV: " + tsv_path);
        }

        reader.skip_header_if_present();

        // Thread-local dictionary builders
        int n_threads = 1;
#ifdef _OPENMP
        n_threads = omp_get_max_threads();
#endif

        // Use spilling dictionaries for memory-bounded operation
        SpillingStringDict read_dict_builder("reads", config.memory_limit / 2);
        SpillingStringDict ref_dict_builder("refs", config.memory_limit / 4);

        uint64_t n_lines = 0;
        std::mutex dict_mutex;

        // Process file in morsels
        std::vector<Morsel> morsel_pool(n_threads);
        for (auto& m : morsel_pool) m.reserve(MORSEL_SIZE);

        // Read and process morsels
        while (true) {
            // Fill a morsel with lines
            Morsel& morsel = morsel_pool[0];
            morsel.clear();

            for (size_t i = 0; i < MORSEL_SIZE; ++i) {
                std::string_view line = reader.next_line();
                if (line.empty() && reader.bytes_read() > 0) {
                    // Check if truly EOF (empty line at end)
                    break;
                }
                if (line.empty()) break;

                // Parse first two fields for dictionary
                const char* p = line.data();
                const char* end = p + line.size();
                const char* tab1 = simd_find_tab(p, end);
                if (tab1 >= end) continue;
                const char* tab2 = simd_find_tab(tab1 + 1, end);

                std::string_view query(p, static_cast<size_t>(tab1 - p));
                std::string_view target(tab1 + 1, static_cast<size_t>(tab2 - tab1 - 1));

                // Copy to owned strings for dictionary
                morsel.queries.emplace_back(query);
                morsel.targets.emplace_back(target);
            }

            if (morsel.queries.empty()) break;

            n_lines += morsel.queries.size();

            // Insert into dictionaries (thread-safe via mutex for now)
            // Could parallelize with sharded insertion
            {
                std::lock_guard<std::mutex> lock(dict_mutex);
                for (const auto& q : morsel.queries) {
                    read_dict_builder.insert(q);
                }
                for (const auto& t : morsel.targets) {
                    ref_dict_builder.insert(t);
                }
            }

            // Progress
            if (config.verbose && n_lines % 1000000 == 0) {
                std::cerr << "\rPass 1: " << (n_lines / 1000000) << "M lines, "
                          << (reader.bytes_read() / (1024*1024)) << " MB read...";
            }
        }

        // Finalize dictionaries
        auto [read_names_vec_str, read_lookup] = read_dict_builder.finalize();
        auto [ref_names_vec_str, ref_lookup] = ref_dict_builder.finalize();

        auto pass1_end = std::chrono::high_resolution_clock::now();
        auto pass1_ms = std::chrono::duration_cast<std::chrono::milliseconds>(pass1_end - pass1_start).count();

        if (config.verbose) {
            std::cerr << "\rPass 1: " << n_lines << " lines, "
                      << read_names_vec_str.size() << " reads, "
                      << ref_names_vec_str.size() << " refs"
                      << " (" << pass1_ms << " ms)";
            if (read_dict_builder.spill_count() > 0 || ref_dict_builder.spill_count() > 0) {
                std::cerr << " [spilled: " << read_dict_builder.spill_count()
                          << "+" << ref_dict_builder.spill_count() << " files]";
            }
            std::cerr << "\n";
        }

        // =====================================================================
        // Pass 2: Stream through file again, write columnar row-groups
        // =====================================================================
        auto pass2_start = std::chrono::high_resolution_clock::now();

        reader.reset();
        reader.skip_header_if_present();

        // Prepare output file
        std::ofstream out(output_path, std::ios::binary);
        if (!out) {
            throw std::runtime_error("Cannot create file: " + output_path);
        }

        uint64_t total_rows = n_lines;
        const size_t rg_size = ROW_GROUP_SIZE;
        const size_t num_row_groups = (n_lines + rg_size - 1) / rg_size;

        // Build name offset tables
        std::vector<uint32_t> read_name_offsets;
        std::vector<char> read_names_data;
        read_name_offsets.reserve(read_names_vec_str.size());
        for (const auto& s : read_names_vec_str) {
            read_name_offsets.push_back(static_cast<uint32_t>(read_names_data.size()));
            read_names_data.insert(read_names_data.end(), s.begin(), s.end());
            read_names_data.push_back('\0');
        }

        std::vector<uint32_t> ref_name_offsets;
        std::vector<char> ref_names_data;
        ref_name_offsets.reserve(ref_names_vec_str.size());
        for (const auto& s : ref_names_vec_str) {
            ref_name_offsets.push_back(static_cast<uint32_t>(ref_names_data.size()));
            ref_names_data.insert(ref_names_data.end(), s.begin(), s.end());
            ref_names_data.push_back('\0');
        }

        // Write provisional header
        EMIHeader header{};
        header.magic = EMI_MAGIC;
        header.version = EMI_VERSION;
        header.num_alignments = total_rows;
        header.num_reads = static_cast<uint32_t>(read_name_offsets.size());
        header.num_refs = static_cast<uint32_t>(ref_name_offsets.size());
        header.num_row_groups = static_cast<uint32_t>(num_row_groups);
        header.row_group_size = static_cast<uint32_t>(rg_size);
        header.d_max = config.d_max;
        header.lambda = config.lambda;

        out.write(reinterpret_cast<const char*>(&header), sizeof(header));

        // Reserve space for row group directory
        header.row_group_dir_offset = out.tellp();
        std::vector<RowGroupMeta> row_groups(num_row_groups);
        out.write(reinterpret_cast<const char*>(row_groups.data()),
                  row_groups.size() * sizeof(RowGroupMeta));

        // Write string dictionaries
        header.string_dict_offset = out.tellp();
        uint32_t read_names_count = static_cast<uint32_t>(read_name_offsets.size());
        out.write(reinterpret_cast<const char*>(&read_names_count), sizeof(read_names_count));
        out.write(reinterpret_cast<const char*>(read_name_offsets.data()),
                  read_name_offsets.size() * sizeof(uint32_t));
        uint32_t read_names_size = static_cast<uint32_t>(read_names_data.size());
        out.write(reinterpret_cast<const char*>(&read_names_size), sizeof(read_names_size));
        out.write(read_names_data.data(), read_names_data.size());

        uint32_t ref_names_count = static_cast<uint32_t>(ref_name_offsets.size());
        out.write(reinterpret_cast<const char*>(&ref_names_count), sizeof(ref_names_count));
        out.write(reinterpret_cast<const char*>(ref_name_offsets.data()),
                  ref_name_offsets.size() * sizeof(uint32_t));
        uint32_t ref_names_size = static_cast<uint32_t>(ref_names_data.size());
        out.write(reinterpret_cast<const char*>(&ref_names_size), sizeof(ref_names_size));
        out.write(ref_names_data.data(), ref_names_data.size());

        // Release name vectors after writing (keep lookups for row-group processing)
        std::vector<uint32_t>().swap(read_name_offsets);
        std::vector<char>().swap(read_names_data);
        std::vector<uint32_t>().swap(ref_name_offsets);
        std::vector<char>().swap(ref_names_data);
        std::vector<std::string>().swap(read_names_vec_str);
        std::vector<std::string>().swap(ref_names_vec_str);

        // CSR offset building: count alignments per read
        std::vector<uint32_t> read_counts(header.num_reads, 0);

        // Pre-allocate row-group sized buffers (reused for each row-group)
        std::vector<uint32_t> rg_read_idx(rg_size);
        std::vector<uint32_t> rg_ref_idx(rg_size);
        std::vector<float> rg_bit_score(rg_size);
        std::vector<float> rg_damage_score(rg_size);
        std::vector<float> rg_evalue_log10(rg_size);
        std::vector<uint16_t> rg_identity_q(rg_size);
        std::vector<uint16_t> rg_aln_len(rg_size);
        std::vector<uint16_t> rg_tstart(rg_size);
        std::vector<uint16_t> rg_tend(rg_size);
        std::vector<uint16_t> rg_tlen(rg_size);
        std::vector<uint16_t> rg_qlen(rg_size);
        std::vector<uint32_t> rg_qaln_offsets(rg_size + 1);
        std::vector<char> rg_qaln_data;
        rg_qaln_data.reserve(rg_size * 50);
        std::vector<uint32_t> rg_taln_offsets(rg_size + 1);
        std::vector<char> rg_taln_data;
        rg_taln_data.reserve(rg_size * 50);

        // Storage for parsed lines in current row-group (owned copies)
        std::vector<std::string> rg_lines(rg_size);

        // Process each row-group by streaming through file
        uint64_t lines_processed = 0;

        for (size_t rg = 0; rg < num_row_groups; ++rg) {
            size_t rg_start = rg * rg_size;
            size_t rg_end = std::min(rg_start + rg_size, static_cast<size_t>(n_lines));
            uint32_t rg_rows = static_cast<uint32_t>(rg_end - rg_start);

            // Clear string buffers for this row-group
            rg_qaln_data.clear();
            rg_taln_data.clear();

            // Read lines for this row-group (streaming, no mmap)
            for (uint32_t i = 0; i < rg_rows; ++i) {
                std::string_view line = reader.next_line();
                rg_lines[i] = std::string(line);  // Copy to owned storage
            }

            // Fill numeric columns (parallel within row-group)
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
            for (uint32_t local_idx = 0; local_idx < rg_rows; ++local_idx) {
                const std::string& line = rg_lines[local_idx];
                const char* p = line.data();
                const char* line_end = p + line.size();

                // Parse all 16 fields
                const char* fields[16];
                size_t lens[16];
                const char* fp = p;
                for (int f = 0; f < 16; ++f) {
                    fields[f] = fp;
                    const char* tab = (f == 15) ? line_end : simd_find_tab(fp, line_end);
                    lens[f] = static_cast<size_t>(tab - fp);
                    fp = (tab < line_end) ? tab + 1 : line_end;
                }

                std::string_view q(fields[0], lens[0]);
                std::string_view t(fields[1], lens[1]);

                // Lookup in dictionary
                auto it_q = read_lookup.find(std::string(q));
                auto it_t = ref_lookup.find(std::string(t));
                rg_read_idx[local_idx] = (it_q != read_lookup.end()) ? it_q->second : 0;
                rg_ref_idx[local_idx] = (it_t != ref_lookup.end()) ? it_t->second : 0;

                float fident = fast_stof({fields[2], lens[2]});
                uint32_t alnlen = fast_stou({fields[3], lens[3]});
                uint32_t tstart = fast_stou({fields[8], lens[8]});
                uint32_t tend = fast_stou({fields[9], lens[9]});
                float evalue = fast_stof({fields[10], lens[10]});
                float bits = fast_stof({fields[11], lens[11]});
                uint32_t qlen_val = fast_stou({fields[12], lens[12]});
                uint32_t tlen_val = fast_stou({fields[13], lens[13]});

                rg_bit_score[local_idx] = bits;
                rg_damage_score[local_idx] = 0.0f;
                rg_evalue_log10[local_idx] = (evalue > 0) ? std::log10(evalue) : -300.0f;
                rg_identity_q[local_idx] = static_cast<uint16_t>(std::clamp(fident, 0.0f, 1.0f) * 65535.0f);
                rg_aln_len[local_idx] = static_cast<uint16_t>(std::min(alnlen, 65535u));
                rg_tstart[local_idx] = static_cast<uint16_t>(std::min(tstart, 65535u));
                rg_tend[local_idx] = static_cast<uint16_t>(std::min(tend, 65535u));
                rg_tlen[local_idx] = static_cast<uint16_t>(std::min(tlen_val, 65535u));
                rg_qlen[local_idx] = static_cast<uint16_t>(std::min(qlen_val, 65535u));
            }

            // Fill string columns (sequential, uses owned line copies)
            if (!config.skip_alignments) {
                for (uint32_t local_idx = 0; local_idx < rg_rows; ++local_idx) {
                    const std::string& line = rg_lines[local_idx];
                    const char* p = line.data();
                    const char* line_end = p + line.size();

                    const char* fp = p;
                    for (int f = 0; f < 14 && fp < line_end; ++f) {
                        fp = simd_find_tab(fp, line_end) + 1;
                    }
                    const char* qaln_start = fp;
                    const char* qaln_end = simd_find_tab(fp, line_end);
                    fp = (qaln_end < line_end) ? qaln_end + 1 : line_end;
                    const char* taln_start = fp;
                    const char* taln_end = simd_find_tab(fp, line_end);

                    rg_qaln_offsets[local_idx] = static_cast<uint32_t>(rg_qaln_data.size());
                    rg_qaln_data.insert(rg_qaln_data.end(), qaln_start, qaln_end);
                    rg_taln_offsets[local_idx] = static_cast<uint32_t>(rg_taln_data.size());
                    rg_taln_data.insert(rg_taln_data.end(), taln_start, taln_end);
                }
            }
            rg_qaln_offsets[rg_rows] = static_cast<uint32_t>(rg_qaln_data.size());
            rg_taln_offsets[rg_rows] = static_cast<uint32_t>(rg_taln_data.size());

            lines_processed += rg_rows;

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

            // Apply permutation to all numeric columns
            auto apply_perm = [&](auto& vec) {
                using T = typename std::remove_reference<decltype(vec)>::type::value_type;
                std::vector<T> temp(rg_rows);
                for (uint32_t i = 0; i < rg_rows; ++i) temp[i] = vec[perm[i]];
                std::copy(temp.begin(), temp.end(), vec.begin());
            };
            apply_perm(rg_read_idx);
            apply_perm(rg_ref_idx);
            apply_perm(rg_bit_score);
            apply_perm(rg_damage_score);
            apply_perm(rg_evalue_log10);
            apply_perm(rg_identity_q);
            apply_perm(rg_aln_len);
            apply_perm(rg_tstart);
            apply_perm(rg_tend);
            apply_perm(rg_tlen);
            apply_perm(rg_qlen);

            // Rebuild qaln/taln with permuted order
            if (!config.skip_alignments) {
                std::vector<std::string> qaln_strs(rg_rows), taln_strs(rg_rows);
                for (uint32_t i = 0; i < rg_rows; ++i) {
                    uint32_t old_idx = perm[i];
                    uint32_t start = rg_qaln_offsets[old_idx];
                    uint32_t end = rg_qaln_offsets[old_idx + 1];
                    qaln_strs[i] = std::string(rg_qaln_data.begin() + start, rg_qaln_data.begin() + end);

                    start = rg_taln_offsets[old_idx];
                    end = rg_taln_offsets[old_idx + 1];
                    taln_strs[i] = std::string(rg_taln_data.begin() + start, rg_taln_data.begin() + end);
                }
                rg_qaln_data.clear();
                rg_taln_data.clear();
                for (uint32_t i = 0; i < rg_rows; ++i) {
                    rg_qaln_offsets[i] = static_cast<uint32_t>(rg_qaln_data.size());
                    rg_qaln_data.insert(rg_qaln_data.end(), qaln_strs[i].begin(), qaln_strs[i].end());
                    rg_taln_offsets[i] = static_cast<uint32_t>(rg_taln_data.size());
                    rg_taln_data.insert(rg_taln_data.end(), taln_strs[i].begin(), taln_strs[i].end());
                }
                rg_qaln_offsets[rg_rows] = static_cast<uint32_t>(rg_qaln_data.size());
                rg_taln_offsets[rg_rows] = static_cast<uint32_t>(rg_taln_data.size());
            }

            // Accumulate CSR counts
            for (uint32_t i = 0; i < rg_rows; ++i) {
                read_counts[rg_read_idx[i]]++;
            }

            // Write row-group to file
            row_groups[rg].num_rows = rg_rows;
            row_groups[rg].first_read_idx = rg_read_idx[0];

            auto write_chunk = [&](ColumnID col, const void* chunk_data, size_t elem_size) {
                auto& chunk = row_groups[rg].columns[static_cast<size_t>(col)];
                chunk.offset = out.tellp();
                chunk.uncompressed_size = static_cast<uint32_t>(elem_size * rg_rows);
                chunk.compressed_size = chunk.uncompressed_size;
                chunk.codec = Codec::NONE;
                out.write(static_cast<const char*>(chunk_data), chunk.uncompressed_size);
            };

            write_chunk(ColumnID::READ_IDX, rg_read_idx.data(), sizeof(uint32_t));
            write_chunk(ColumnID::REF_IDX, rg_ref_idx.data(), sizeof(uint32_t));
            write_chunk(ColumnID::BIT_SCORE, rg_bit_score.data(), sizeof(float));
            write_chunk(ColumnID::DAMAGE_SCORE, rg_damage_score.data(), sizeof(float));
            write_chunk(ColumnID::EVALUE_LOG10, rg_evalue_log10.data(), sizeof(float));
            write_chunk(ColumnID::IDENTITY_Q, rg_identity_q.data(), sizeof(uint16_t));
            write_chunk(ColumnID::ALN_LEN, rg_aln_len.data(), sizeof(uint16_t));
            write_chunk(ColumnID::TSTART, rg_tstart.data(), sizeof(uint16_t));
            write_chunk(ColumnID::TEND, rg_tend.data(), sizeof(uint16_t));
            write_chunk(ColumnID::TLEN, rg_tlen.data(), sizeof(uint16_t));
            write_chunk(ColumnID::QLEN, rg_qlen.data(), sizeof(uint16_t));

            // Write string columns (combined offsets + data in QALN/TALN columns)
            if (!config.skip_alignments) {
                auto& qaln_chunk = row_groups[rg].columns[static_cast<size_t>(ColumnID::QALN)];
                qaln_chunk.offset = out.tellp();
                out.write(reinterpret_cast<const char*>(rg_qaln_offsets.data()), (rg_rows + 1) * sizeof(uint32_t));
                out.write(rg_qaln_data.data(), rg_qaln_data.size());
                qaln_chunk.uncompressed_size = static_cast<uint32_t>((rg_rows + 1) * sizeof(uint32_t) + rg_qaln_data.size());
                qaln_chunk.compressed_size = qaln_chunk.uncompressed_size;
                qaln_chunk.codec = Codec::NONE;

                auto& taln_chunk = row_groups[rg].columns[static_cast<size_t>(ColumnID::TALN)];
                taln_chunk.offset = out.tellp();
                out.write(reinterpret_cast<const char*>(rg_taln_offsets.data()), (rg_rows + 1) * sizeof(uint32_t));
                out.write(rg_taln_data.data(), rg_taln_data.size());
                taln_chunk.uncompressed_size = static_cast<uint32_t>((rg_rows + 1) * sizeof(uint32_t) + rg_taln_data.size());
                taln_chunk.compressed_size = taln_chunk.uncompressed_size;
                taln_chunk.codec = Codec::NONE;
            }

            // Progress
            if (config.verbose && (rg + 1) % 10 == 0) {
                std::cerr << "\rPass 2: " << lines_processed << "/" << n_lines << " rows...";
            }
        }

        // Build CSR offsets from counts (prefix sum)
        header.csr_offsets_offset = out.tellp();
        uint32_t csr_count = header.num_reads + 1;
        std::vector<uint64_t> csr_offsets(csr_count);
        csr_offsets[0] = 0;
        for (uint32_t r = 0; r < header.num_reads; ++r) {
            csr_offsets[r + 1] = csr_offsets[r] + read_counts[r];
        }
        out.write(reinterpret_cast<const char*>(csr_offsets.data()), csr_count * sizeof(uint64_t));

        // Seek back and update row group directory
        out.seekp(header.row_group_dir_offset);
        out.write(reinterpret_cast<const char*>(row_groups.data()),
                  row_groups.size() * sizeof(RowGroupMeta));

        // Seek back and update header
        out.seekp(0);
        out.write(reinterpret_cast<const char*>(&header), sizeof(header));

        auto pass2_end = std::chrono::high_resolution_clock::now();
        auto pass2_ms = std::chrono::duration_cast<std::chrono::milliseconds>(pass2_end - pass2_start).count();
        auto total_ms = std::chrono::duration_cast<std::chrono::milliseconds>(pass2_end - total_start).count();

        if (config.verbose) {
            std::cerr << "\rPass 2: " << total_rows << " rows written (" << pass2_ms << " ms)\n";
            std::cerr << "Total: " << (total_rows * 1000 / std::max(1LL, static_cast<long long>(total_ms)))
                      << " rows/sec, " << total_ms << " ms\n";
        }

        return total_rows;
    }

};  // class FastColumnarWriter

}  // namespace agp
