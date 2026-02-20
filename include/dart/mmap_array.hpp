#pragma once
/**
 * @file mmap_array.hpp
 * @brief Memory-mapped arrays for handling 10^9+ elements efficiently.
 *
 * Provides MmapArray (vector replacement backed by mmap), CompactAlignment
 * (cache-friendly 24-byte alignment record), and StringTable (deduplicated
 * string storage for read/reference names).
 */

#include <cstdint>
#include <cstring>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>
#include <unordered_map>
#include <vector>

#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

namespace dart {

// ---------------------------------------------------------------------------
// MmapArray: mmap-backed array for large datasets
// ---------------------------------------------------------------------------

template<typename T>
class MmapArray {
    static_assert(std::is_trivially_copyable_v<T>,
                  "MmapArray requires trivially copyable types");

public:
    MmapArray() = default;

    explicit MmapArray(size_t count) {
        alloc_anonymous(count);
    }

    MmapArray(const std::string& path, size_t count) {
        alloc_file_backed(path, count);
    }

    ~MmapArray() { release(); }

    MmapArray(const MmapArray&) = delete;
    MmapArray& operator=(const MmapArray&) = delete;

    MmapArray(MmapArray&& other) noexcept
        : data_(other.data_)
        , size_(other.size_)
        , capacity_(other.capacity_)
        , mapped_bytes_(other.mapped_bytes_)
        , fd_(other.fd_)
        , file_backed_(other.file_backed_) {
        other.data_ = nullptr;
        other.size_ = 0;
        other.capacity_ = 0;
        other.mapped_bytes_ = 0;
        other.fd_ = -1;
        other.file_backed_ = false;
    }

    MmapArray& operator=(MmapArray&& other) noexcept {
        if (this != &other) {
            release();
            data_ = other.data_;
            size_ = other.size_;
            capacity_ = other.capacity_;
            mapped_bytes_ = other.mapped_bytes_;
            fd_ = other.fd_;
            file_backed_ = other.file_backed_;
            other.data_ = nullptr;
            other.size_ = 0;
            other.capacity_ = 0;
            other.mapped_bytes_ = 0;
            other.fd_ = -1;
            other.file_backed_ = false;
        }
        return *this;
    }

    T& operator[](size_t idx) { return data_[idx]; }
    const T& operator[](size_t idx) const { return data_[idx]; }

    T* data() { return data_; }
    const T* data() const { return data_; }
    size_t size() const { return size_; }
    size_t capacity() const { return capacity_; }
    bool empty() const { return size_ == 0; }

    T* begin() { return data_; }
    T* end() { return data_ + size_; }
    const T* begin() const { return data_; }
    const T* end() const { return data_ + size_; }

    void push_back(const T& val) {
        if (size_ >= capacity_) {
            grow(capacity_ ? capacity_ * 2 : 1024);
        }
        data_[size_++] = val;
    }

    void resize(size_t new_size) {
        if (new_size > capacity_) {
            grow(new_size);
        }
        if (new_size > size_) {
            std::memset(data_ + size_, 0, (new_size - size_) * sizeof(T));
        }
        size_ = new_size;
    }

    void clear() { size_ = 0; }

    void advise_sequential() {
        if (data_ && mapped_bytes_) {
            madvise(data_, mapped_bytes_, MADV_SEQUENTIAL);
        }
    }

    void advise_random() {
        if (data_ && mapped_bytes_) {
            madvise(data_, mapped_bytes_, MADV_RANDOM);
        }
    }

    void advise_willneed() {
        if (data_ && mapped_bytes_) {
            madvise(data_, mapped_bytes_, MADV_WILLNEED);
        }
    }

private:
    T* data_ = nullptr;
    size_t size_ = 0;
    size_t capacity_ = 0;
    size_t mapped_bytes_ = 0;
    int fd_ = -1;
    bool file_backed_ = false;

    static size_t page_align(size_t bytes) {
        static const size_t page = static_cast<size_t>(sysconf(_SC_PAGESIZE));
        return (bytes + page - 1) & ~(page - 1);
    }

    void alloc_anonymous(size_t count) {
        if (count == 0) return;
        mapped_bytes_ = page_align(count * sizeof(T));
        void* ptr = mmap(nullptr, mapped_bytes_,
                         PROT_READ | PROT_WRITE,
                         MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
        if (ptr == MAP_FAILED) {
            throw std::runtime_error(
                "MmapArray: anonymous mmap failed for " +
                std::to_string(mapped_bytes_) + " bytes");
        }
        data_ = static_cast<T*>(ptr);
        size_ = count;
        capacity_ = mapped_bytes_ / sizeof(T);
    }

    void alloc_file_backed(const std::string& path, size_t count) {
        if (count == 0) return;
        fd_ = open(path.c_str(), O_RDWR | O_CREAT | O_TRUNC, 0600);
        if (fd_ < 0) {
            throw std::runtime_error(
                "MmapArray: failed to open " + path);
        }
        file_backed_ = true;

        mapped_bytes_ = page_align(count * sizeof(T));
        if (ftruncate(fd_, static_cast<off_t>(mapped_bytes_)) < 0) {
            ::close(fd_);
            fd_ = -1;
            throw std::runtime_error(
                "MmapArray: ftruncate failed for " + path);
        }

        void* ptr = mmap(nullptr, mapped_bytes_,
                         PROT_READ | PROT_WRITE,
                         MAP_SHARED, fd_, 0);
        if (ptr == MAP_FAILED) {
            ::close(fd_);
            fd_ = -1;
            throw std::runtime_error(
                "MmapArray: file-backed mmap failed for " + path);
        }
        data_ = static_cast<T*>(ptr);
        size_ = count;
        capacity_ = mapped_bytes_ / sizeof(T);
    }

    void grow(size_t new_capacity) {
        size_t new_bytes = page_align(new_capacity * sizeof(T));

        if (data_) {
#ifdef __linux__
            void* ptr = mremap(data_, mapped_bytes_, new_bytes, MREMAP_MAYMOVE);
            if (ptr == MAP_FAILED) {
                throw std::runtime_error(
                    "MmapArray: mremap failed for " +
                    std::to_string(new_bytes) + " bytes");
            }
            if (file_backed_ && fd_ >= 0) {
                if (ftruncate(fd_, static_cast<off_t>(new_bytes)) < 0) {
                    throw std::runtime_error("MmapArray: ftruncate failed during grow");
                }
            }
            data_ = static_cast<T*>(ptr);
#else
            // Fallback: allocate new mapping, copy, unmap old
            int flags = file_backed_
                ? MAP_SHARED
                : (MAP_PRIVATE | MAP_ANONYMOUS);
            int map_fd = file_backed_ ? fd_ : -1;

            if (file_backed_ && fd_ >= 0) {
                if (ftruncate(fd_, static_cast<off_t>(new_bytes)) < 0) {
                    throw std::runtime_error("MmapArray: ftruncate failed during grow");
                }
            }

            void* ptr = mmap(nullptr, new_bytes,
                             PROT_READ | PROT_WRITE, flags, map_fd, 0);
            if (ptr == MAP_FAILED) {
                throw std::runtime_error(
                    "MmapArray: mmap failed during grow for " +
                    std::to_string(new_bytes) + " bytes");
            }
            std::memcpy(ptr, data_, size_ * sizeof(T));
            munmap(data_, mapped_bytes_);
            data_ = static_cast<T*>(ptr);
#endif
        } else {
            if (file_backed_) {
                // data_ is null but file_backed_ is true: the backing file path
                // was not retained, so we cannot re-open it. This is an invalid
                // state that should never arise in normal usage.
                throw std::logic_error(
                    "MmapArray::grow: cannot grow a file-backed array that was "
                    "not initially allocated (data_ is null)");
            } else {
                alloc_anonymous(new_capacity);
                size_ = 0;  // alloc_anonymous sets size_ = count
            }
        }
        mapped_bytes_ = new_bytes;
        capacity_ = new_bytes / sizeof(T);
    }

    void release() {
        if (data_) {
            munmap(data_, mapped_bytes_);
            data_ = nullptr;
        }
        if (fd_ >= 0) {
            ::close(fd_);
            fd_ = -1;
        }
        size_ = 0;
        capacity_ = 0;
        mapped_bytes_ = 0;
        file_backed_ = false;
    }
};

// ---------------------------------------------------------------------------
// CompactAlignment: cache-friendly 32-byte alignment record
// ---------------------------------------------------------------------------

struct CompactAlignment {
    uint32_t read_idx;       // Index into read name StringTable
    uint32_t ref_idx;        // Index into reference name StringTable
    float    bit_score;      // Alignment bit score
    float    damage_score;   // Per-read p_damaged prior
    float    damage_ll_a;    // log P(alignment damage evidence | ancient)
    float    damage_ll_m;    // log P(alignment damage evidence | modern)
    uint16_t aln_start;     // Alignment start on reference
    uint16_t aln_end;       // Alignment end on reference
    uint16_t identity_q;    // Identity * 65535 (quantized to [0,1])
    uint16_t flags;          // Reserved for future use
};
static_assert(sizeof(CompactAlignment) == 32,
              "CompactAlignment must be 32 bytes");

// ---------------------------------------------------------------------------
// StringTable: deduplicated string storage
// ---------------------------------------------------------------------------

class StringTable {
public:
    uint32_t insert(std::string_view s) {
        // Use std::string key to avoid dangling pointers when data_ grows
        std::string key(s);
        auto it = index_.find(key);
        if (it != index_.end()) {
            return it->second;
        }
        uint32_t idx = static_cast<uint32_t>(offsets_.size());
        uint32_t offset = static_cast<uint32_t>(data_.size());
        offsets_.push_back(offset);
        data_.insert(data_.end(), s.begin(), s.end());
        data_.push_back('\0');
        index_.emplace(std::move(key), idx);
        return idx;
    }

    std::string_view get(uint32_t idx) const {
        if (idx >= offsets_.size()) {
            throw std::out_of_range(
                "StringTable::get: index " + std::to_string(idx) +
                " >= size " + std::to_string(offsets_.size()));
        }
        const char* start = &data_[offsets_[idx]];
        // Length is distance to next offset (or end of data minus null)
        size_t len;
        if (idx + 1 < offsets_.size()) {
            len = offsets_[idx + 1] - offsets_[idx] - 1;  // -1 for null
        } else {
            len = data_.size() - offsets_[idx] - 1;  // -1 for null
        }
        return {start, len};
    }

    size_t size() const { return offsets_.size(); }

    void reserve(size_t n_strings, size_t avg_len = 32) {
        offsets_.reserve(n_strings);
        data_.reserve(n_strings * (avg_len + 1));
        index_.reserve(n_strings);
    }

private:
    std::vector<char> data_;
    std::vector<uint32_t> offsets_;
    std::unordered_map<std::string, uint32_t> index_;  // Use string, not string_view
};

}  // namespace dart
