#pragma once

#include "dart/dart_columnar.hpp"
#include <array>
#include <deque>
#include <filesystem>
#include <fstream>
#include <future>
#include <limits>
#include <queue>
#include <stdexcept>
#include <unistd.h>
#ifdef __linux__
#include <sys/vfs.h>
#include <linux/magic.h>
#endif

namespace dcx {

inline bool record_less(const Record& a, const Record& b) {
    if (a.target_id != b.target_id) return a.target_id < b.target_id;
    if (a.tstart != b.tstart) return a.tstart < b.tstart;
    return a.src_idx < b.src_idx;
}

inline size_t record_memory(const Record& r) {
    // Include vector growth and the payload capacities, not just serialized bytes.
    return 2 * sizeof(Record) + r.query_id.capacity() + r.target_id.capacity() +
        r.presence.capacity() + r.ref_aa.capacity() + r.obs_aa.capacity() +
        r.dcls.capacity() + r.obs_codon.capacity() + r.best_codon.capacity() +
        (r.p_damage.capacity() + r.best_post.capacity()) * sizeof(float) +
        (r.dist5.capacity() + r.dist3.capacity()) * sizeof(uint16_t);
}

// Temporary runs are private to this instance. Appending is serialized by the
// caller; sorting and run writes overlap annotation in a bounded worker window.
class ExternalSort {
    struct Header {
        uint64_t src;
        uint32_t read, ref, start, span, query_bytes, target_bytes;
        float gamma, ancient, margin, neff, identity, null_prob, p_read, posterior;
        uint32_t lengths[10];
        uint8_t strand, keep, rank;
    };
    struct Input {
        std::array<char, 1 << 16> buffer;
        std::ifstream stream;
        Record record;
        explicit Input(const std::string& path) {
            stream.rdbuf()->pubsetbuf(buffer.data(), buffer.size());
            stream.open(path, std::ios::binary);
            if (!stream) throw std::runtime_error("dcx sort: cannot open run");
        }
        bool next() {
            if (stream.peek() == std::char_traits<char>::eof()) {
                if (stream.bad()) throw std::runtime_error("dcx sort: run read failed");
                return false;
            }
            Header h{};
            stream.read(reinterpret_cast<char*>(&h), sizeof(h));
            record.read_idx = h.read; record.ref_idx = h.ref; record.src_idx = h.src;
            record.tstart = h.start; record.span = h.span;
            record.gamma = h.gamma; record.gamma_ancient = h.ancient;
            record.em_margin = h.margin; record.em_ref_neff = h.neff;
            record.fident = h.identity; record.gamma_null = h.null_prob;
            record.p_read = h.p_read; record.posterior = h.posterior;
            record.strand = h.strand; record.em_keep = h.keep; record.candidate_rank = h.rank;
            auto read = [&](auto& v, size_t n) {
                v.resize(n);
                if (n) stream.read(reinterpret_cast<char*>(v.data()), n * sizeof(v[0]));
            };
            if (!stream) throw std::runtime_error("dcx sort: truncated run header");
            read(record.query_id, h.query_bytes); read(record.target_id, h.target_bytes);
            read(record.presence, h.lengths[0]); read(record.ref_aa, h.lengths[1]);
            read(record.obs_aa, h.lengths[2]); read(record.dcls, h.lengths[3]);
            read(record.obs_codon, h.lengths[4]); read(record.p_damage, h.lengths[5]);
            read(record.dist5, h.lengths[6]); read(record.dist3, h.lengths[7]);
            read(record.best_codon, h.lengths[8]); read(record.best_post, h.lengths[9]);
            if (!stream) throw std::runtime_error("dcx sort: truncated run record");
            return true;
        }
    };

    std::string base_, directory_;
    size_t workers_, chunk_cap_, used_ = 0, fail_run_, next_id_ = 0;
    std::vector<Record> buffer_;
    std::vector<std::string> owned_, runs_;
    std::deque<std::future<void>> pending_;
    bool finished_ = false;

    std::string new_path() {
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
        std::string path = directory_ + "/run_" + std::to_string(next_id_++);
        owned_.push_back(path);
        return path;
    }

    static void write(std::ostream& out, const Record& r) {
        Header h{};
        h.src = r.src_idx; h.read = r.read_idx; h.ref = r.ref_idx;
        h.start = r.tstart; h.span = r.span;
        h.query_bytes = static_cast<uint32_t>(r.query_id.size());
        h.target_bytes = static_cast<uint32_t>(r.target_id.size());
        h.gamma = r.gamma; h.ancient = r.gamma_ancient;
        h.margin = r.em_margin; h.neff = r.em_ref_neff;
        h.identity = r.fident; h.null_prob = r.gamma_null;
        h.p_read = r.p_read; h.posterior = r.posterior;
        h.strand = r.strand; h.keep = r.em_keep; h.rank = r.candidate_rank;
        const size_t lengths[] = {r.presence.size(), r.ref_aa.size(), r.obs_aa.size(), r.dcls.size(),
            r.obs_codon.size(), r.p_damage.size(), r.dist5.size(), r.dist3.size(), r.best_codon.size(), r.best_post.size()};
        for (size_t i = 0; i < 10; ++i) h.lengths[i] = static_cast<uint32_t>(lengths[i]);
        out.write(reinterpret_cast<const char*>(&h), sizeof(h));
        auto bytes = [&](const auto& v) {
            if (!v.empty()) out.write(reinterpret_cast<const char*>(v.data()), v.size() * sizeof(v[0]));
        };
        bytes(r.query_id); bytes(r.target_id); bytes(r.presence); bytes(r.ref_aa);
        bytes(r.obs_aa); bytes(r.dcls); bytes(r.obs_codon); bytes(r.p_damage);
        bytes(r.dist5); bytes(r.dist3); bytes(r.best_codon); bytes(r.best_post);
        if (!out) throw std::runtime_error("dcx sort: run write failed");
    }

    template<class Consume>
    static void merge(const std::vector<std::string>& paths, Consume consume) {
        std::vector<std::unique_ptr<Input>> inputs;
        auto greater = [&](size_t a, size_t b) { return record_less(inputs[b]->record, inputs[a]->record); };
        std::priority_queue<size_t, std::vector<size_t>, decltype(greater)> heap(greater);
        for (const auto& path : paths) {
            inputs.push_back(std::make_unique<Input>(path));
            if (inputs.back()->next()) heap.push(inputs.size() - 1);
        }
        while (!heap.empty()) {
            const size_t i = heap.top(); heap.pop();
            consume(inputs[i]->record);
            if (inputs[i]->next()) heap.push(i);
        }
    }

    void spill() {
        if (buffer_.empty()) return;
        if (pending_.size() >= workers_) { pending_.front().get(); pending_.pop_front(); }
        const std::string path = new_path();
        const bool fail = runs_.size() == fail_run_;
        runs_.push_back(path);
        pending_.emplace_back(std::async(std::launch::async,
            [records = std::move(buffer_), path, fail]() mutable {
                std::sort(records.begin(), records.end(), record_less);
                std::array<char, 1 << 20> buf;
                std::ofstream out;
                out.rdbuf()->pubsetbuf(buf.data(), buf.size());
                out.open(path, std::ios::binary | std::ios::trunc);
                if (!out) throw std::runtime_error("dcx sort: cannot create run");
                if (fail) { out.put('\0'); throw std::runtime_error("dcx sort: injected run failure"); }
                for (const auto& rec : records) write(out, rec);
                out.close();
                if (!out) throw std::runtime_error("dcx sort: run close failed");
            }));
        buffer_ = {};
        used_ = 0;
    }

public:
    ExternalSort(std::string directory, size_t memory_budget, size_t workers = 4,
                 size_t fail_run = std::numeric_limits<size_t>::max())
        : base_(std::move(directory)), workers_(std::max(size_t(1), workers)),
          chunk_cap_(std::max(size_t(1), memory_budget / (workers_ + 1))), fail_run_(fail_run) {}
    ~ExternalSort() {
        for (auto& future : pending_) { try { future.get(); } catch (...) {} }
        if (!directory_.empty()) {
            std::error_code ec;
            std::filesystem::remove_all(directory_, ec);
        }
    }
    ExternalSort(const ExternalSort&) = delete;
    ExternalSort& operator=(const ExternalSort&) = delete;

    void append(Record rec) {
        if (finished_) throw std::runtime_error("dcx sort: append after finish");
        used_ += record_memory(rec);
        buffer_.push_back(std::move(rec));
        if (used_ >= chunk_cap_) spill();
    }

    template<class Consume>
    void finish(Consume consume) {
        if (finished_) throw std::runtime_error("dcx sort: already finished");
        finished_ = true;
        if (runs_.empty()) {
            std::sort(buffer_.begin(), buffer_.end(), record_less);
            for (const auto& rec : buffer_) consume(rec);
            std::vector<Record>().swap(buffer_);
            return;
        }
        spill();
        while (!pending_.empty()) { pending_.front().get(); pending_.pop_front(); }
        constexpr size_t fan_in = 64;
        while (runs_.size() > fan_in) {
            std::vector<std::string> next;
            for (size_t i = 0; i < runs_.size(); i += fan_in) {
                const size_t end = std::min(i + fan_in, runs_.size());
                std::vector<std::string> part(runs_.begin() + i, runs_.begin() + end);
                const std::string path = new_path();
                std::ofstream out(path, std::ios::binary);
                if (!out) throw std::runtime_error("dcx sort: cannot create merge run");
                merge(part, [&](const Record& rec) { write(out, rec); });
                out.close();
                if (!out) throw std::runtime_error("dcx sort: merge run close failed");
                next.push_back(path);
                for (const auto& old : part) std::filesystem::remove(old);
            }
            runs_ = std::move(next);
        }
        merge(runs_, consume);
    }
};

} // namespace dcx
