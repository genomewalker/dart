#include "dart/dcx_external_sort.hpp"
#include <cstdlib>
#include <iostream>
#include <map>
#include <random>

static void require(bool ok, const char* message) {
    if (!ok) throw std::runtime_error(message);
}

static std::string bytes(const std::string& path) {
    std::ifstream in(path, std::ios::binary);
    require(bool(in), "missing output");
    return {std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>()};
}

static std::vector<dcx::Record> records() {
    std::vector<dcx::Record> out;
    for (uint32_t i = 0; i < 240; ++i) {
        dcx::Record r;
        r.read_idx = i % 31; r.ref_idx = i % 7; r.src_idx = i;
        r.query_id = "read_" + std::to_string(r.read_idx);
        r.target_id = "target_" + std::to_string(r.ref_idx);
        r.tstart = i % 5; r.span = 9 + i % 13;
        r.gamma = (i % 11) / 10.0f; r.gamma_ancient = r.gamma * 0.6f;
        r.posterior = (i % 17) / 16.0f; r.p_read = 0.6f;
        r.fident = 0.95f; r.em_margin = 0.4f; r.em_ref_neff = 1.4f;
        r.gamma_null = 0.1f; r.em_keep = (i % 3 == 0);
        r.strand = i % 2; r.candidate_rank = i % 5;
        r.presence.assign((r.span + 7) / 8, 0);
        for (uint32_t p = 0; p < r.span; ++p) {
            if (p % 3 == 1) continue;
            r.presence[p >> 3] |= uint8_t(1u << (p & 7));
            r.ref_aa.push_back('R'); r.obs_aa.push_back('W'); r.dcls.push_back(p % 4);
            r.obs_codon.push_back(p % 64); r.best_codon.push_back((p + 3) % 64);
            r.p_damage.push_back((p % 7) / 6.0f); r.best_post.push_back(0.9f);
            r.dist5.push_back(p); r.dist3.push_back(r.span - p - 1);
        }
        out.push_back(std::move(r));
    }
    std::mt19937 rng(1847);
    std::shuffle(out.begin(), out.end(), rng);
    return out;
}

static dcx::Batch batch(const std::vector<dcx::Record>& sorted, uint32_t flags) {
    dcx::Batch b; b.col_flags = flags;
    std::map<std::string, uint64_t> qids;
    for (const auto& r : sorted) {
        if (b.target_dict.empty() || b.target_dict.back() != r.target_id) b.target_dict.push_back(r.target_id);
        auto [q, inserted] = qids.try_emplace(r.query_id, b.query_dict.size());
        if (inserted) b.query_dict.push_back(r.query_id);
        b.query_ord.push_back(q->second); b.target_ord.push_back(b.target_dict.size() - 1);
        b.tstart.push_back(r.tstart); b.span.push_back(r.span);
        b.gamma.push_back(r.gamma); b.posterior.push_back(r.posterior);
        b.em_keep.push_back(r.em_keep); b.strand.push_back(r.strand);
        b.gamma_ancient.push_back(r.gamma_ancient); b.em_margin.push_back(r.em_margin);
        b.em_ref_neff.push_back(r.em_ref_neff); b.gamma_null.push_back(r.gamma_null);
        b.p_read.push_back(r.p_read); b.candidate_rank.push_back(r.candidate_rank); b.fident.push_back(r.fident);
        auto add = [](auto& a, const auto& v) { a.insert(a.end(), v.begin(), v.end()); };
        add(b.presence, r.presence); add(b.ref_aa, r.ref_aa); add(b.obs_aa, r.obs_aa);
        add(b.damage_class, r.dcls); add(b.obs_codon, r.obs_codon); add(b.best_codon, r.best_codon);
        add(b.p_damage, r.p_damage); add(b.best_post, r.best_post);
        add(b.dist5, r.dist5); add(b.dist3, r.dist3);
    }
    return b;
}

static void check(const std::string& dir, const std::vector<dcx::Record>& input, uint32_t flags, size_t budget) {
    auto ordered = input;
    std::sort(ordered.begin(), ordered.end(), dcx::record_less);
    const dcx::Batch expected = batch(ordered, flags);
    dcx::DcxWriteOptions opt;
    opt.col_flags = flags; opt.block_raw_cap = 256;
    opt.calib_contract_id[0] = 1; opt.calib_model = "stream-test"; opt.library_label = "same-library";
    require(dcx::write_dcx((dir + "/legacy.dcx").c_str(), expected, opt), "legacy write failed");
    require(dcx::write_dcx_auth((dir + "/legacy.auth").c_str(), expected, opt), "legacy auth failed");
    {
        dcx::ExternalSort sort(dir + "/runs", budget, 2);
        for (const auto& r : input) sort.append(r);
        dcx::StreamWriter full((dir + "/stream.dcx").c_str(), opt);
        dcx::StreamWriter auth((dir + "/stream.auth").c_str(), opt, true);
        std::map<std::string, uint64_t> qids;
        std::vector<std::string> targets, queries;
        size_t n = 0;
        sort.finish([&](const dcx::Record& r) {
            require(!dcx::record_less(r, ordered.at(n)), "merge order mismatch");
            require(!dcx::record_less(ordered.at(n), r), "merge tie mismatch");
            if (targets.empty() || targets.back() != r.target_id) targets.push_back(r.target_id);
            auto [q, inserted] = qids.try_emplace(r.query_id, queries.size());
            if (inserted) queries.push_back(r.query_id);
            full.append(r, targets.size() - 1, q->second, r.posterior);
            auth.append(r, targets.size() - 1, q->second, r.posterior);
            ++n;
        });
        require(n == input.size(), "lost records");
        full.close(targets, queries); auth.close(targets, queries);
    }
    require(bytes(dir + "/legacy.dcx") == bytes(dir + "/stream.dcx"), "dcx bytes differ");
    require(bytes(dir + "/legacy.auth") == bytes(dir + "/stream.auth"), "auth bytes differ");
    if (std::filesystem::exists(dir + "/runs"))
        require(std::filesystem::is_empty(dir + "/runs"), "normal-exit runs leaked");
}

int main() {
    try {
        const char* base = std::getenv("TMPDIR");
        std::string pattern = std::string(base && *base ? base : ".") + "/dcx_stream_XXXXXX";
        std::vector<char> path(pattern.begin(), pattern.end()); path.push_back('\0');
        require(::mkdtemp(path.data()) != nullptr, "cannot create test directory");
        const std::string dir = path.data();
        auto input = records();
        const uint32_t flags = dcx::COL_QID_DICT | dcx::COL_STRAND | dcx::COL_P_READ |
            dcx::COL_GAMMA_ANCIENT | dcx::COL_EM_MARGIN | dcx::COL_EM_REF_NEFF | dcx::COL_GAMMA_NULL |
            dcx::COL_OBS_CODON | dcx::COL_BEST_CODON | dcx::COL_CANDIDATE_RANK | dcx::COL_SOFT_MULTIMAP | dcx::COL_FIDENT;
        for (uint32_t cf : {uint32_t(0), flags, flags | dcx::COL_DIST3_DERIVED}) {
            check(dir, input, cf, size_t(1) << 28);
            check(dir, input, cf, 2048);
            check(dir, {}, cf, 2048);
        }
        bool failed = false;
        try {
            dcx::ExternalSort sort(dir + "/runs", 2048, 2, 1);
            for (const auto& r : input) sort.append(r);
            sort.finish([](const dcx::Record&) {});
        } catch (const std::runtime_error&) { failed = true; }
        require(failed, "failure injection did not fire");
        require(std::filesystem::is_empty(dir + "/runs"), "error-exit runs leaked");
        std::filesystem::remove_all(dir);
        std::cout << "dcx streaming: byte identity, multi-pass merge and cleanup PASS\n";
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n'; return 1;
    }
}
