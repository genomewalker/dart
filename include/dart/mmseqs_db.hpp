#ifndef DART_MMSEQS_DB_HPP
#define DART_MMSEQS_DB_HPP
// Minimal, dependency-free reader for mmseqs2 sequence/result DBs plus an alignment-record
// parser byte-compatible with mmseqs2 Matcher::readAlignmentResults for the 10/11-column
// (no-ORF) result format. Replaces the DBReader<unsigned int> + Matcher usage in the native
// streaming merge so dart carries no mmseqs dependency.
//
// Semantics mirror DBReader opened with NOSORT: DBReader::sortIndex reorders the index to
// ascending id whenever the on-disk index isn't already sorted (only HARDNOSORT/SORT_BY_OFFSET
// skip it, DBReader.cpp:296). We replicate that here — after reading .index in file order we
// sort (id,offset) by ascending id — so getData(i)/getDbKey(i) iterate in the SAME order the
// native tool saw. Result DBs are NOT id-sorted on disk; without this sort the read/ref
// first-encounter order (hence every emitted emi byte) diverges from native. Split data files
// (base.0/base.1/...) are mmapped in numeric order and addressed by a cumulative-size offset
// table, exactly as FileUtil::findDatafiles + DBReader::getDataByOffset.
//
// getData(key)-as-position (target/query fetch by dbKey) stays valid because seq DBs are dense
// (id==position after the id sort), matching the native tool's own T.getData((size_t)dbKey).
#include <cstdint>
#include <cstring>
#include <cstdlib>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <limits>
#include <algorithm>
#include <cmath>
#include <charconv>
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>

namespace dart {
namespace mmseqs {

// ---- tokenizer replicas (byte-identical to mmseqs Util) ----
inline size_t skip_ws(const char* d) { size_t c = 0; while (d[c] == ' ' || d[c] == '\t') c++; return c; }
inline size_t skip_nonws(const char* d) {
    size_t c = 0;
    while ((d[c] == ' ' || d[c] == '\t' || d[c] == '\n' || d[c] == '\0') == false) c++;
    return c;
}
inline size_t words_of_line(const char* d, const char** w, size_t maxE) {
    size_t n = 0;
    while (*d != '\n' && *d != '\0') {
        d += skip_ws(d);
        w[n] = d; n++;
        if (n >= maxE) return n;
        d += skip_nonws(d);
    }
    if (n < maxE) w[n] = d;
    return n;
}
template <class T> inline T fast_atoi(const char* s) {
    T val = 0; int sign = 1;
    if (std::numeric_limits<T>::is_signed && *s == '-') { sign = -1; s++; }
    while (*s >= '0' && *s <= '9') val = val * 10 + (*s++ - '0');
    return sign * val;
}
inline char* skip_line(char* d) { while (*d != '\n') d++; return d + 1; }
// mmseqs Matcher::uncompressAlignment: run-length CIGAR -> expanded M/I/D string.
inline std::string uncompress_aln(const char* b, const char* e) {
    std::string bt; bt.reserve(e - b);
    size_t count = 0;
    for (const char* p = b; p < e; ++p) {
        char c = *p;
        if (c >= '0' && c <= '9') count = count * 10 + c - '0';
        else { bt.append(count == 0 ? 1 : count, c); count = 0; }
    }
    return bt;
}
// mmseqs Matcher::computeAlnLength with the -1 -> 0 start adjustment of parseAlignmentRecord.
inline int compute_aln_length(int qStart, int qEnd, int dbStart, int dbEnd) {
    int aqs = (qStart == -1) ? 0 : qStart;
    int ads = (dbStart == -1) ? 0 : dbStart;
    return std::max(std::abs(qEnd - aqs), std::abs(dbEnd - ads)) + 1;
}

// Field names mirror Matcher::result_t so the merge body reads unchanged.
struct AlnRec {
    uint32_t dbKey = 0;
    int score = 0;
    float seqId = 0.0f;
    double eval = 0.0;
    uint32_t alnLength = 0;
    int qStartPos = 0, qEndPos = 0;
    uint32_t qLen = 0;
    int dbStartPos = 0, dbEndPos = 0;
    uint32_t dbLen = 0;
    std::string backtrace;   // expanded (uncompressed)
};

// Parse one DB entry (NUL-terminated concatenation of result lines) into AlnRec rows.
// Mirrors Matcher::readAlignmentResults(result, data, /*readCompressed=*/false).
inline void parse_results(const char* data, std::vector<AlnRec>& out) {
    if (data == nullptr) return;
    char* d = const_cast<char*>(data);
    const char* e[16];
    while (*d != '\0') {
        size_t cols = words_of_line(d, e, 16);
        if (cols >= 10) {
            AlnRec r;
            r.dbKey = fast_atoi<uint32_t>(e[0]);
            r.score = fast_atoi<int>(e[1]);
            double sid = 0.0, ev = 0.0;
            std::from_chars(e[2], e[3] - 1, sid); r.seqId = static_cast<float>(sid);
            std::from_chars(e[3], e[4] - 1, ev);  r.eval = ev;
            r.qStartPos = fast_atoi<int>(e[4]); r.qEndPos = fast_atoi<int>(e[5]); r.qLen = fast_atoi<uint32_t>(e[6]);
            r.dbStartPos = fast_atoi<int>(e[7]); r.dbEndPos = fast_atoi<int>(e[8]); r.dbLen = fast_atoi<uint32_t>(e[9]);
            r.alnLength = static_cast<uint32_t>(compute_aln_length(r.qStartPos, r.qEndPos, r.dbStartPos, r.dbEndPos));
            if (cols >= 11) r.backtrace = uncompress_aln(e[10], e[11]);
            out.push_back(std::move(r));
        }
        d = skip_line(d);
    }
}

// Minimal DBReader (NOSORT / sorted-by-id) over mmseqs data + .index.
class DB {
public:
    explicit DB(std::string path) : base_(std::move(path)) {}
    ~DB() { for (auto& f : files_) if (f.first) ::munmap(const_cast<char*>(f.first), f.second); }
    DB(const DB&) = delete; DB& operator=(const DB&) = delete;

    void open() {
        // data files: base.0, base.1, ... in numeric order (FileUtil::findDatafiles); else base.
        std::vector<std::string> names;
        for (size_t c = 0; ; ++c) {
            std::string n = base_ + "." + std::to_string(c);
            struct stat st; if (::stat(n.c_str(), &st) != 0) break;
            names.push_back(n);
        }
        if (names.empty()) names.push_back(base_);
        cum_.push_back(0);
        for (const auto& n : names) {
            int fd = ::open(n.c_str(), O_RDONLY);
            if (fd < 0) throw std::runtime_error("Cannot open data file " + n);
            struct stat st; if (::fstat(fd, &st) != 0) { ::close(fd); throw std::runtime_error("fstat " + n); }
            size_t sz = static_cast<size_t>(st.st_size);
            const char* p = nullptr;
            if (sz) { void* m = ::mmap(nullptr, sz, PROT_READ, MAP_PRIVATE, fd, 0); if (m == MAP_FAILED) { ::close(fd); throw std::runtime_error("mmap " + n); } p = static_cast<const char*>(m); }
            ::close(fd);
            files_.push_back({p, sz});
            cum_.push_back(cum_.back() + sz);
        }
        read_index(base_ + ".index");
    }

    size_t getSize() const { return id_.size(); }
    uint32_t getDbKey(size_t i) const { return id_[i]; }
    const char* getData(size_t id, int /*thr*/) const { return data_at(off_[id]); }

private:
    void read_index(const std::string& idx) {
        std::ifstream in(idx);
        if (!in) throw std::runtime_error("Cannot open index file " + idx);
        std::string line;
        const char* cols[3];
        while (std::getline(in, line)) {
            if (line.empty()) continue;
            words_of_line(line.c_str(), cols, 3);
            id_.push_back(fast_atoi<uint32_t>(cols[0]));
            off_.push_back(fast_atoi<uint64_t>(cols[1]));
        }
        // DBReader NOSORT: reorder to ascending id if the on-disk index isn't already sorted.
        bool sorted = true;
        for (size_t i = 1; i < id_.size(); ++i) if (id_[i] < id_[i - 1]) { sorted = false; break; }
        if (sorted) return;
        std::vector<uint32_t> perm(id_.size());
        for (uint32_t i = 0; i < perm.size(); ++i) perm[i] = i;
        std::sort(perm.begin(), perm.end(), [this](uint32_t a, uint32_t b) { return id_[a] < id_[b]; });
        std::vector<uint32_t> nid(id_.size());
        std::vector<uint64_t> noff(off_.size());
        for (size_t i = 0; i < perm.size(); ++i) { nid[i] = id_[perm[i]]; noff[i] = off_[perm[i]]; }
        id_.swap(nid);
        off_.swap(noff);
    }
    const char* data_at(uint64_t off) const {
        size_t cnt = 0;
        while ((off >= cum_[cnt] && off < cum_[cnt + 1]) == false) cnt++;
        return files_[cnt].first + (off - cum_[cnt]);
    }

    std::string base_;
    std::vector<std::pair<const char*, size_t>> files_;
    std::vector<uint64_t> cum_;
    std::vector<uint32_t> id_;
    std::vector<uint64_t> off_;
};

// Load a .lookup as key -> name (key in column 0, name in column 1); vector sized n.
inline std::vector<std::string> load_lookup_vec(const std::string& p, size_t n) {
    std::vector<std::string> v(n);
    std::ifstream in(p); std::string l;
    while (std::getline(in, l)) {
        size_t a = l.find('\t'); if (a == std::string::npos) continue;
        size_t b = l.find('\t', a + 1);
        uint32_t k = (uint32_t)strtoul(l.c_str(), 0, 10);
        if (k < n) v[k] = l.substr(a + 1, (b == std::string::npos ? l.size() : b) - a - 1);
    }
    return v;
}

}  // namespace mmseqs
}  // namespace dart

#endif  // DART_MMSEQS_DB_HPP
