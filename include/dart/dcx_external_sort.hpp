#pragma once

#include "dart/dart_columnar.hpp"
#include <algorithm>
#include <atomic>
#include <deque>
#include <cstring>
#include <filesystem>
#include <limits>
#include <mutex>
#include <queue>
#include <stdexcept>
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#ifdef __linux__
#include <sys/vfs.h>
#include <linux/magic.h>
#endif

namespace dcx {

// External sort of .dcx records by (target rank, tstart, source ordinal), the order the
// .dcx is written in. rank is the target's position in byte order of target names, so
// integer keys reproduce the string order. Each thread gets its own Sink on first use
// (local()): it serializes records into its own arena and, when the arena is full, sorts
// and writes it as a run file itself. No lock is taken per record. finish() merges every run through
// mmap. Run files live in a private directory on node-local disk and are removed with
// the ExternalSort.
class ExternalSort {
    struct Header {
        uint64_t src;
        uint32_t rank, start, read, ref, span, query_bytes;
        float gamma, ancient, margin, neff, identity, null_prob, p_read, posterior;
        uint32_t lengths[10];
        uint8_t strand, keep, crank, pad;
    };
    struct Key { uint32_t rank, start; uint64_t src, off; };
    static bool key_less(const Key& a, const Key& b) {
        if (a.rank != b.rank) return a.rank < b.rank;
        if (a.start != b.start) return a.start < b.start;
        return a.src < b.src;
    }

public:
    class Sink {
        friend class ExternalSort;
        ExternalSort* owner_ = nullptr;
        std::vector<char> arena_;
        std::vector<Key> keys_;
    public:
        void append(const Record& r, uint32_t rank) {
            Header h{};
            h.src = r.src_idx; h.rank = rank; h.start = r.tstart;
            h.read = r.read_idx; h.ref = r.ref_idx; h.span = r.span;
            h.query_bytes = static_cast<uint32_t>(r.query_id.size());
            h.gamma = r.gamma; h.ancient = r.gamma_ancient;
            h.margin = r.em_margin; h.neff = r.em_ref_neff;
            h.identity = r.fident; h.null_prob = r.gamma_null;
            h.p_read = r.p_read; h.posterior = r.posterior;
            h.strand = r.strand; h.keep = r.em_keep; h.crank = r.candidate_rank;
            const size_t n[10] = {r.presence.size(), r.ref_aa.size(), r.obs_aa.size(), r.dcls.size(),
                r.obs_codon.size(), r.p_damage.size(), r.dist5.size(), r.dist3.size(),
                r.best_codon.size(), r.best_post.size()};
            for (int i = 0; i < 10; ++i) h.lengths[i] = static_cast<uint32_t>(n[i]);
            keys_.push_back({rank, r.tstart, r.src_idx, arena_.size()});
            put(&h, sizeof(h));
            put(r.query_id.data(), r.query_id.size());
            put(r.presence.data(), r.presence.size()); put(r.ref_aa.data(), r.ref_aa.size());
            put(r.obs_aa.data(), r.obs_aa.size()); put(r.dcls.data(), r.dcls.size());
            put(r.obs_codon.data(), r.obs_codon.size());
            put(r.p_damage.data(), r.p_damage.size() * sizeof(float));
            put(r.dist5.data(), r.dist5.size() * sizeof(uint16_t));
            put(r.dist3.data(), r.dist3.size() * sizeof(uint16_t));
            put(r.best_codon.data(), r.best_codon.size());
            put(r.best_post.data(), r.best_post.size() * sizeof(float));
            if (arena_.size() >= owner_->sink_cap_) spill();
        }
        void spill() {
            if (keys_.empty()) return;
            std::sort(keys_.begin(), keys_.end(), key_less);
            size_t index = 0;
            const std::string path = owner_->new_run(index);
            const int fd = ::open(path.c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0600);
            if (fd < 0) throw std::runtime_error("dcx sort: cannot create run");
            std::vector<char> out;
            out.reserve(std::min<size_t>(arena_.size(), size_t(64) << 20));
            auto flush = [&] {
                const char* p = out.data(); size_t left = out.size();
                while (left) {
                    const ssize_t w = ::write(fd, p, left);
                    if (w <= 0) { ::close(fd); throw std::runtime_error("dcx sort: run write failed"); }
                    p += w; left -= size_t(w);
                }
                out.clear();
            };
            if (index == owner_->fail_run_) { ::close(fd); throw std::runtime_error("dcx sort: injected run failure"); }
            for (size_t i = 0; i < keys_.size(); ++i) {
                const char* rec = arena_.data() + keys_[i].off;
                const size_t len = record_size(rec);
                if (out.size() + len > out.capacity() && !out.empty()) flush();
                out.insert(out.end(), rec, rec + len);
            }
            flush();
            if (::close(fd) != 0) throw std::runtime_error("dcx sort: run close failed");
            arena_.clear(); keys_.clear();
        }
    private:
        void put(const void* p, size_t n) {
            if (n) arena_.insert(arena_.end(), static_cast<const char*>(p), static_cast<const char*>(p) + n);
        }
    };

    // memory_budget is shared by `threads` sinks; each spills at budget / threads.
    ExternalSort(std::string directory, size_t memory_budget, int threads,
                 size_t fail_run = std::numeric_limits<size_t>::max())
        : base_(std::move(directory)),
          sink_cap_(std::max(size_t(1), memory_budget / size_t(std::max(1, threads)))),
          fail_run_(fail_run), id_(next_id().fetch_add(1)) {}
    ~ExternalSort() {
        if (!directory_.empty()) {
            std::error_code ec;
            std::filesystem::remove_all(directory_, ec);
        }
    }
    ExternalSort(const ExternalSort&) = delete;
    ExternalSort& operator=(const ExternalSort&) = delete;

    // The calling thread's sink. Safe from nested OpenMP regions, where
    // omp_get_thread_num() is 0 in every thread.
    Sink& local() {
        thread_local uint64_t owner = 0;
        thread_local Sink* sink = nullptr;
        if (owner != id_) {
            std::lock_guard<std::mutex> lock(dir_mutex_);
            sinks_.emplace_back();
            sinks_.back().owner_ = this;
            sink = &sinks_.back();
            owner = id_;
        }
        return *sink;
    }
    size_t runs() const { return runs_.size(); }

    // Spill what is left, then k-way merge every run in key order. consume() sees each
    // record once, decoded into a reused Record.
    template<class Consume>
    void finish(Consume consume) {
        std::vector<Sink*> left;
        for (auto& s : sinks_) left.push_back(&s);
        std::exception_ptr error;
        #pragma omp parallel for schedule(dynamic, 1)
        for (size_t k = 0; k < left.size(); ++k) {
            try { left[k]->spill(); }
            catch (...) {
                #pragma omp critical(dcx_sort_finish)
                if (!error) error = std::current_exception();
            }
        }
        if (error) std::rethrow_exception(error);
        struct Cursor {
            const char* base = nullptr; size_t size = 0, pos = 0, released = 0;
            Header h{};
        };
        std::vector<Cursor> cur(runs_.size());
        for (size_t i = 0; i < runs_.size(); ++i) {
            const int fd = ::open(runs_[i].c_str(), O_RDONLY);
            if (fd < 0) throw std::runtime_error("dcx sort: cannot open run");
            struct stat st{};
            ::fstat(fd, &st);
            cur[i].size = size_t(st.st_size);
            if (cur[i].size) {
                void* m = ::mmap(nullptr, cur[i].size, PROT_READ, MAP_PRIVATE, fd, 0);
                if (m == MAP_FAILED) { ::close(fd); throw std::runtime_error("dcx sort: cannot map run"); }
                ::madvise(m, cur[i].size, MADV_SEQUENTIAL);
                cur[i].base = static_cast<const char*>(m);
            }
            ::close(fd);
        }
        auto load = [&](Cursor& c) {
            if (c.pos >= c.size) return false;
            std::memcpy(&c.h, c.base + c.pos, sizeof(Header));
            return true;
        };
        auto greater = [&](size_t a, size_t b) {
            const Header& x = cur[a].h; const Header& y = cur[b].h;
            if (x.rank != y.rank) return x.rank > y.rank;
            if (x.start != y.start) return x.start > y.start;
            return x.src > y.src;
        };
        std::priority_queue<size_t, std::vector<size_t>, decltype(greater)> heap(greater);
        for (size_t i = 0; i < cur.size(); ++i) if (load(cur[i])) heap.push(i);
        Record rec;
        constexpr size_t release_step = size_t(64) << 20;
        while (!heap.empty()) {
            const size_t i = heap.top(); heap.pop();
            Cursor& c = cur[i];
            c.pos = decode(c.base + c.pos, rec) - c.base;
            consume(rec, c.h.rank);
            // Drop pages already merged so resident memory stays near the read-ahead window.
            if (c.pos - c.released >= release_step) {
                const size_t page = size_t(::sysconf(_SC_PAGESIZE));
                const size_t upto = c.pos / page * page;
                ::madvise(const_cast<char*>(c.base) + c.released, upto - c.released, MADV_DONTNEED);
                c.released = upto;
            }
            if (load(c)) heap.push(i);
        }
        for (auto& c : cur) if (c.base) ::munmap(const_cast<char*>(c.base), c.size);
    }

private:
    static size_t record_size(const char* p) {
        Header h; std::memcpy(&h, p, sizeof(h));
        return sizeof(Header) + h.query_bytes + h.lengths[0] + h.lengths[1] + h.lengths[2] +
            h.lengths[3] + h.lengths[4] + h.lengths[5] * sizeof(float) +
            (h.lengths[6] + h.lengths[7]) * sizeof(uint16_t) + h.lengths[8] +
            h.lengths[9] * sizeof(float);
    }
    static const char* decode(const char* p, Record& r) {
        Header h; std::memcpy(&h, p, sizeof(h)); p += sizeof(h);
        r.src_idx = h.src; r.tstart = h.start; r.read_idx = h.read; r.ref_idx = h.ref; r.span = h.span;
        r.gamma = h.gamma; r.gamma_ancient = h.ancient; r.em_margin = h.margin; r.em_ref_neff = h.neff;
        r.fident = h.identity; r.gamma_null = h.null_prob; r.p_read = h.p_read; r.posterior = h.posterior;
        r.strand = h.strand; r.em_keep = h.keep; r.candidate_rank = h.crank;
        auto take = [&](auto& v, size_t n) {
            v.resize(n);
            const size_t bytes = n * sizeof(v[0]);
            if (bytes) std::memcpy(v.data(), p, bytes);
            p += bytes;
        };
        take(r.query_id, h.query_bytes);
        take(r.presence, h.lengths[0]); take(r.ref_aa, h.lengths[1]); take(r.obs_aa, h.lengths[2]);
        take(r.dcls, h.lengths[3]); take(r.obs_codon, h.lengths[4]); take(r.p_damage, h.lengths[5]);
        take(r.dist5, h.lengths[6]); take(r.dist3, h.lengths[7]);
        take(r.best_codon, h.lengths[8]); take(r.best_post, h.lengths[9]);
        return p;
    }

    std::string new_run(size_t& index) {
        std::lock_guard<std::mutex> lock(dir_mutex_);
        if (directory_.empty()) {
            if (base_.empty()) throw std::runtime_error("dcx sort: run files need --em-tape-dir or DART_EM_TAPE_DIR");
            if (base_ == "/tmp" || base_.compare(0, 5, "/tmp/") == 0)
                throw std::runtime_error("dcx sort: use node-local disk, not /tmp");
            std::filesystem::create_directories(base_);
#ifdef __linux__
            struct statfs disk{};
            if (::statfs(base_.c_str(), &disk) != 0 || disk.f_type == NFS_SUPER_MAGIC || disk.f_type == TMPFS_MAGIC)
                throw std::runtime_error("dcx sort: run directory must be on node-local disk");
#endif
            std::string pattern = base_ + "/dart_dcx_XXXXXX";
            std::vector<char> path(pattern.begin(), pattern.end()); path.push_back('\0');
            if (!::mkdtemp(path.data())) throw std::runtime_error("dcx sort: cannot create run directory");
            directory_ = path.data();
        }
        index = runs_.size();
        runs_.push_back(directory_ + "/run_" + std::to_string(index));
        return runs_.back();
    }

    static std::atomic<uint64_t>& next_id() { static std::atomic<uint64_t> n{1}; return n; }

    std::string base_, directory_;
    std::deque<Sink> sinks_;
    size_t sink_cap_, fail_run_;
    uint64_t id_;
    std::vector<std::string> runs_;
    std::mutex dir_mutex_;
};

} // namespace dcx
