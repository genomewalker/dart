// EM on pre-indexed alignment data
//
// Uses mmap for zero-copy access to alignment records.
// All parsing done once in em-prepare, then em-solve is pure compute.

#include "agp/em_index.hpp"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <fcntl.h>
#include <fstream>
#include <stdexcept>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace agp {

// ---------------------------------------------------------------------------
// EMIndexWriter
// ---------------------------------------------------------------------------

EMIndexWriter::EMIndexWriter(const std::string& path) : path_(path) {
    read_offsets_.push_back(0);  // First read starts at offset 0
}

EMIndexWriter::~EMIndexWriter() = default;

uint32_t EMIndexWriter::add_ref(std::string_view name) {
    uint32_t idx = num_refs_++;
    ref_name_offsets_.push_back(static_cast<uint32_t>(ref_names_.size()));
    ref_names_.insert(ref_names_.end(), name.begin(), name.end());
    ref_names_.push_back('\0');
    return idx;
}

void EMIndexWriter::add_read(std::string_view name,
                              const std::vector<EMIAlignment>& alignments) {
    // Add read name
    read_names_.insert(read_names_.end(), name.begin(), name.end());
    read_names_.push_back('\0');

    // Add alignments
    alignments_.insert(alignments_.end(), alignments.begin(), alignments.end());

    // Update CSR offset
    read_offsets_.push_back(static_cast<uint32_t>(alignments_.size()));
    num_reads_++;
}

void EMIndexWriter::finalize(float d_max, float lambda) {
    std::ofstream out(path_, std::ios::binary);
    if (!out) {
        throw std::runtime_error("Cannot create EMI file: " + path_);
    }

    // Compute offsets
    uint64_t header_size = sizeof(EMIHeader);
    uint64_t string_table_size = read_names_.size() + ref_names_.size();
    uint64_t alignments_size = alignments_.size() * sizeof(EMIAlignment);
    uint64_t offsets_size = read_offsets_.size() * sizeof(uint64_t);

    uint64_t string_table_offset = header_size;
    uint64_t alignments_offset = string_table_offset + string_table_size;
    uint64_t offsets_offset = alignments_offset + alignments_size;

    // Write header
    EMIHeader header;
    header.magic = EMI_MAGIC;
    header.version = 1;
    header.flags = 0;
    header.num_reads = num_reads_;
    header.num_refs = num_refs_;
    header.num_alignments = alignments_.size();
    header.string_table_offset = string_table_offset;
    header.string_table_size = static_cast<uint32_t>(string_table_size);
    header.alignments_offset = alignments_offset;
    header.offsets_offset = offsets_offset;
    header.d_max = d_max;
    header.lambda = lambda;

    out.write(reinterpret_cast<const char*>(&header), sizeof(header));

    // Write string table (read names, then ref names)
    out.write(read_names_.data(), read_names_.size());
    out.write(ref_names_.data(), ref_names_.size());

    // Write alignments
    out.write(reinterpret_cast<const char*>(alignments_.data()), alignments_size);

    // Write offsets as uint64_t
    std::vector<uint64_t> offsets64(read_offsets_.begin(), read_offsets_.end());
    out.write(reinterpret_cast<const char*>(offsets64.data()), offsets_size);

    out.close();
}

// ---------------------------------------------------------------------------
// EMIndexReader
// ---------------------------------------------------------------------------

EMIndexReader::EMIndexReader(const std::string& path) {
    fd_ = open(path.c_str(), O_RDONLY);
    if (fd_ < 0) {
        throw std::runtime_error("Cannot open EMI file: " + path);
    }

    struct stat st;
    if (fstat(fd_, &st) < 0) {
        close(fd_);
        throw std::runtime_error("Cannot stat EMI file: " + path);
    }
    file_size_ = static_cast<size_t>(st.st_size);

    data_ = mmap(nullptr, file_size_, PROT_READ, MAP_PRIVATE, fd_, 0);
    if (data_ == MAP_FAILED) {
        close(fd_);
        throw std::runtime_error("Cannot mmap EMI file: " + path);
    }

    // Advise sequential access
    madvise(data_, file_size_, MADV_SEQUENTIAL);

    // Parse header
    header_ = static_cast<const EMIHeader*>(data_);
    if (header_->magic != EMI_MAGIC) {
        munmap(data_, file_size_);
        close(fd_);
        throw std::runtime_error("Invalid EMI magic number");
    }

    // Set up pointers
    const char* base = static_cast<const char*>(data_);
    string_table_ = base + header_->string_table_offset;
    alignments_ = reinterpret_cast<const EMIAlignment*>(base + header_->alignments_offset);
    offsets_ = reinterpret_cast<const uint64_t*>(base + header_->offsets_offset);

    // Pre-compute name offsets for O(1) lookup
    read_name_offsets_.reserve(header_->num_reads);
    ref_name_offsets_.reserve(header_->num_refs);

    const char* p = string_table_;
    for (uint32_t i = 0; i < header_->num_reads; ++i) {
        read_name_offsets_.push_back(static_cast<uint32_t>(p - string_table_));
        p += strlen(p) + 1;
    }
    for (uint32_t i = 0; i < header_->num_refs; ++i) {
        ref_name_offsets_.push_back(static_cast<uint32_t>(p - string_table_));
        p += strlen(p) + 1;
    }
}

EMIndexReader::~EMIndexReader() {
    if (data_ && data_ != MAP_FAILED) {
        munmap(data_, file_size_);
    }
    if (fd_ >= 0) {
        close(fd_);
    }
}

std::string_view EMIndexReader::read_name(uint32_t idx) const {
    if (idx >= read_name_offsets_.size()) return {};
    const char* p = string_table_ + read_name_offsets_[idx];
    return {p, strlen(p)};
}

std::string_view EMIndexReader::ref_name(uint32_t idx) const {
    if (idx >= ref_name_offsets_.size()) return {};
    const char* p = string_table_ + ref_name_offsets_[idx];
    return {p, strlen(p)};
}

// ---------------------------------------------------------------------------
// EM solver on indexed data
// ---------------------------------------------------------------------------

// Log-sum-exp for numerical stability
static inline double log_sum_exp(const double* vals, size_t n) {
    if (n == 0) return -std::numeric_limits<double>::infinity();
    if (n == 1) return vals[0];

    double max_val = vals[0];
    for (size_t i = 1; i < n; ++i) {
        if (vals[i] > max_val) max_val = vals[i];
    }

    double sum = 0.0;
    for (size_t i = 0; i < n; ++i) {
        double diff = vals[i] - max_val;
        sum += (diff > -100.0) ? std::exp(diff) : 0.0;
    }
    return max_val + std::log(sum);
}

EMResult em_solve_indexed(
    const EMIndexReader& index,
    double lambda_b,
    uint32_t max_iters,
    double tol,
    bool use_damage)
{
    const uint32_t R = index.num_reads();
    const uint32_t T = index.num_refs();
    const uint64_t A = index.num_alignments();

    EMResult result;
    result.weights.assign(T, 1.0 / static_cast<double>(T));
    result.gamma.resize(A, 0.0);
    if (use_damage) {
        result.gamma_ancient.resize(A, 0.0);
    }

    const double eps = 1e-15;
    const double min_weight = 1e-10;

    // SQUAREM buffers
    std::vector<double> w0(T), w1(T), w2(T), r_buf(T), v_buf(T);

    double prev_ll = -std::numeric_limits<double>::infinity();

    for (uint32_t iter = 0; iter < max_iters; ++iter) {
        // Save current weights
        std::copy(result.weights.begin(), result.weights.end(), w0.begin());

        // --- E-step 1 ---
        double ll1 = 0.0;
        std::fill(w1.begin(), w1.end(), 0.0);

        #pragma omp parallel reduction(+:ll1)
        {
            std::vector<double> log_scores;
            std::vector<double> w_local(T, 0.0);

            #pragma omp for schedule(dynamic, 1024)
            for (uint32_t r = 0; r < R; ++r) {
                const uint64_t start = index.offset(r);
                const uint64_t end = index.offset(r + 1);
                const uint32_t deg = static_cast<uint32_t>(end - start);
                if (deg == 0) continue;

                // Unique mapper fast path
                if (deg == 1) {
                    const auto& aln = index.alignments()[start];
                    result.gamma[start] = 1.0;
                    w_local[aln.ref_idx] += 1.0;
                    ll1 += std::log(std::max(w0[aln.ref_idx], min_weight))
                         + lambda_b * static_cast<double>(aln.bit_score);
                    continue;
                }

                log_scores.resize(deg);
                for (uint32_t j = 0; j < deg; ++j) {
                    const auto& aln = index.alignments()[start + j];
                    log_scores[j] = std::log(std::max(w0[aln.ref_idx], min_weight))
                                  + lambda_b * static_cast<double>(aln.bit_score);
                }

                double lse = log_sum_exp(log_scores.data(), deg);
                ll1 += lse;

                for (uint32_t j = 0; j < deg; ++j) {
                    double diff = log_scores[j] - lse;
                    double g = (diff > -100.0) ? std::exp(diff) : 0.0;
                    result.gamma[start + j] = g;
                    w_local[index.alignments()[start + j].ref_idx] += g;
                }
            }

            #pragma omp critical
            {
                for (uint32_t t = 0; t < T; ++t) {
                    w1[t] += w_local[t];
                }
            }
        }

        // M-step 1: normalize
        double sum1 = 0.0;
        for (uint32_t t = 0; t < T; ++t) {
            w1[t] = std::max(w1[t] / static_cast<double>(R), min_weight);
            sum1 += w1[t];
        }
        for (uint32_t t = 0; t < T; ++t) w1[t] /= sum1;

        // --- E-step 2 ---
        std::fill(w2.begin(), w2.end(), 0.0);

        #pragma omp parallel
        {
            std::vector<double> log_scores;
            std::vector<double> w_local(T, 0.0);

            #pragma omp for schedule(dynamic, 1024)
            for (uint32_t r = 0; r < R; ++r) {
                const uint64_t start = index.offset(r);
                const uint64_t end = index.offset(r + 1);
                const uint32_t deg = static_cast<uint32_t>(end - start);
                if (deg == 0) continue;

                if (deg == 1) {
                    w_local[index.alignments()[start].ref_idx] += 1.0;
                    continue;
                }

                log_scores.resize(deg);
                for (uint32_t j = 0; j < deg; ++j) {
                    const auto& aln = index.alignments()[start + j];
                    log_scores[j] = std::log(std::max(w1[aln.ref_idx], min_weight))
                                  + lambda_b * static_cast<double>(aln.bit_score);
                }

                double lse = log_sum_exp(log_scores.data(), deg);
                for (uint32_t j = 0; j < deg; ++j) {
                    double diff = log_scores[j] - lse;
                    double g = (diff > -100.0) ? std::exp(diff) : 0.0;
                    w_local[index.alignments()[start + j].ref_idx] += g;
                }
            }

            #pragma omp critical
            {
                for (uint32_t t = 0; t < T; ++t) {
                    w2[t] += w_local[t];
                }
            }
        }

        // M-step 2: normalize
        double sum2 = 0.0;
        for (uint32_t t = 0; t < T; ++t) {
            w2[t] = std::max(w2[t] / static_cast<double>(R), min_weight);
            sum2 += w2[t];
        }
        for (uint32_t t = 0; t < T; ++t) w2[t] /= sum2;

        // --- SQUAREM extrapolation ---
        double r_norm_sq = 0.0, v_norm_sq = 0.0;
        #pragma omp simd reduction(+:r_norm_sq, v_norm_sq)
        for (uint32_t t = 0; t < T; ++t) {
            r_buf[t] = w1[t] - w0[t];
            v_buf[t] = (w2[t] - w1[t]) - r_buf[t];
            r_norm_sq += r_buf[t] * r_buf[t];
            v_norm_sq += v_buf[t] * v_buf[t];
        }

        double alpha = -std::sqrt(r_norm_sq / std::max(v_norm_sq, eps));
        alpha = std::clamp(alpha, -1.0, -0.01);

        #pragma omp simd
        for (uint32_t t = 0; t < T; ++t) {
            w0[t] = w0[t] - 2.0 * alpha * r_buf[t] + alpha * alpha * v_buf[t];
            w0[t] = std::max(w0[t], eps);
        }

        // Re-normalize
        double sum0 = 0.0;
        for (uint32_t t = 0; t < T; ++t) sum0 += w0[t];
        for (uint32_t t = 0; t < T; ++t) w0[t] /= sum0;

        // --- Final E-step at accelerated point ---
        double ll = 0.0;

        #pragma omp parallel reduction(+:ll)
        {
            std::vector<double> log_scores;

            #pragma omp for schedule(dynamic, 1024)
            for (uint32_t r = 0; r < R; ++r) {
                const uint64_t start = index.offset(r);
                const uint64_t end = index.offset(r + 1);
                const uint32_t deg = static_cast<uint32_t>(end - start);
                if (deg == 0) continue;

                if (deg == 1) {
                    const auto& aln = index.alignments()[start];
                    result.gamma[start] = 1.0;
                    ll += std::log(std::max(w0[aln.ref_idx], min_weight))
                        + lambda_b * static_cast<double>(aln.bit_score);
                    continue;
                }

                log_scores.resize(deg);
                for (uint32_t j = 0; j < deg; ++j) {
                    const auto& aln = index.alignments()[start + j];
                    log_scores[j] = std::log(std::max(w0[aln.ref_idx], min_weight))
                                  + lambda_b * static_cast<double>(aln.bit_score);
                }

                double lse = log_sum_exp(log_scores.data(), deg);
                ll += lse;

                for (uint32_t j = 0; j < deg; ++j) {
                    double diff = log_scores[j] - lse;
                    result.gamma[start + j] = (diff > -100.0) ? std::exp(diff) : 0.0;
                }
            }
        }

        // Check for worsening (fall back to w2)
        if (ll < prev_ll && iter > 0) {
            std::copy(w2.begin(), w2.end(), result.weights.begin());
        } else {
            std::copy(w0.begin(), w0.end(), result.weights.begin());
        }

        result.log_likelihood = std::max(ll, prev_ll);
        result.iterations = iter + 1;

        // Convergence check
        if (iter > 0) {
            double rel_change = std::abs(ll - prev_ll) / (std::abs(prev_ll) + eps);
            if (rel_change < tol) break;
        }
        prev_ll = ll;
    }

    return result;
}

} // namespace agp
