#include "agp/damage_index_reader.hpp"

#include <cmath>
#include <cstring>
#include <stdexcept>
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>

namespace agp {

DamageIndexReader::DamageIndexReader(const std::string& path) {
    // Open file
    fd_ = open(path.c_str(), O_RDONLY);
    if (fd_ < 0) {
        throw std::runtime_error("Failed to open damage index: " + path);
    }

    // Get file size
    struct stat st;
    if (fstat(fd_, &st) < 0) {
        ::close(fd_);
        fd_ = -1;
        throw std::runtime_error("Failed to stat damage index: " + path);
    }
    file_size_ = st.st_size;

    // Validate minimum size
    if (file_size_ < sizeof(AgdHeader)) {
        ::close(fd_);
        fd_ = -1;
        throw std::runtime_error("Damage index too small: " + path);
    }

    // Memory-map the file
    data_ = mmap(nullptr, file_size_, PROT_READ, MAP_PRIVATE, fd_, 0);
    if (data_ == MAP_FAILED) {
        data_ = nullptr;
        ::close(fd_);
        fd_ = -1;
        throw std::runtime_error("Failed to mmap damage index: " + path);
    }

    // Prefetch entire file into page cache before random access
    // This is critical for NFS performance - avoids per-lookup network round-trips
    madvise(data_, file_size_, MADV_WILLNEED);

    // Then tell kernel we'll do random access (no sequential readahead needed)
    madvise(data_, file_size_, MADV_RANDOM);

    // Parse header
    header_ = reinterpret_cast<const AgdHeader*>(data_);

    // Validate magic and version
    if (header_->magic != AGD_MAGIC) {
        close();
        throw std::runtime_error("Invalid damage index magic: " + path);
    }
    if (header_->version < 2 || header_->version > AGD_VERSION) {
        close();
        throw std::runtime_error("Unsupported damage index version " +
                                 std::to_string(header_->version) + " (expected 2-" +
                                 std::to_string(AGD_VERSION) + "): " + path);
    }

    // Calculate section offsets
    const char* base = reinterpret_cast<const char*>(data_);
    size_t offset = sizeof(AgdHeader);

    // Hash table buckets
    buckets_ = reinterpret_cast<const AgdBucket*>(base + offset);
    offset += header_->num_buckets * sizeof(AgdBucket);

    // Records
    records_ = reinterpret_cast<const AgdRecord*>(base + offset);
    offset += header_->num_records * sizeof(AgdRecord);

    // Chain pointers (for hash collision resolution)
    chain_ = reinterpret_cast<const uint32_t*>(base + offset);

    // Validate file size
    size_t expected_size = sizeof(AgdHeader) +
                           header_->num_buckets * sizeof(AgdBucket) +
                           header_->num_records * sizeof(AgdRecord) +
                           header_->num_records * sizeof(uint32_t);
    if (file_size_ < expected_size) {
        close();
        throw std::runtime_error("Damage index file truncated: " + path);
    }
}

DamageIndexReader::~DamageIndexReader() {
    close();
}

DamageIndexReader::DamageIndexReader(DamageIndexReader&& other) noexcept
    : data_(other.data_)
    , file_size_(other.file_size_)
    , fd_(other.fd_)
    , header_(other.header_)
    , buckets_(other.buckets_)
    , records_(other.records_)
    , chain_(other.chain_) {
    other.data_ = nullptr;
    other.fd_ = -1;
}

DamageIndexReader& DamageIndexReader::operator=(DamageIndexReader&& other) noexcept {
    if (this != &other) {
        close();
        data_ = other.data_;
        file_size_ = other.file_size_;
        fd_ = other.fd_;
        header_ = other.header_;
        buckets_ = other.buckets_;
        records_ = other.records_;
        chain_ = other.chain_;
        other.data_ = nullptr;
        other.fd_ = -1;
    }
    return *this;
}

void DamageIndexReader::close() {
    if (data_) {
        munmap(data_, file_size_);
        data_ = nullptr;
    }
    if (fd_ >= 0) {
        ::close(fd_);
        fd_ = -1;
    }
    header_ = nullptr;
    buckets_ = nullptr;
    records_ = nullptr;
    chain_ = nullptr;
}

const AgdRecord* DamageIndexReader::find(std::string_view read_id) const {
    if (!is_valid() || header_->num_buckets == 0) {
        return nullptr;
    }

    // Strip AGP suffix and compute hash
    std::string_view base_id = strip_agp_suffix(read_id);
    uint64_t hash = fnv1a_hash(base_id);
    uint64_t bucket_idx = hash % header_->num_buckets;

    // Look up in hash table
    uint32_t rec_idx = buckets_[bucket_idx].record_offset;

    while (rec_idx != 0xFFFFFFFF && rec_idx < header_->num_records) {
        const AgdRecord& rec = records_[rec_idx];
        if (rec.id_hash == hash) {
            // Hash match - this is our record
            // (We don't store the full ID, so hash collision would be a false positive,
            //  but with 64-bit hash this is astronomically unlikely for typical datasets)
            return &rec;
        }
        // Follow chain
        rec_idx = chain_[rec_idx];
    }

    return nullptr;
}

const AgdRecord* DamageIndexReader::get_record(size_t idx) const {
    if (!is_valid() || idx >= header_->num_records) {
        return nullptr;
    }
    return &records_[idx];
}

size_t DamageIndexReader::warmup_cache() const {
    if (!is_valid() || file_size_ == 0) return 0;

    // Touch every page to force it into page cache
    // This is critical for NFS where MADV_WILLNEED may be ignored
    const char* ptr = reinterpret_cast<const char*>(data_);
    const long page_size = sysconf(_SC_PAGESIZE);
    if (page_size <= 0) return 0;

    volatile size_t sum = 0;  // volatile prevents optimizer from removing the loop
    for (size_t offset = 0; offset < file_size_; offset += static_cast<size_t>(page_size)) {
        sum += static_cast<unsigned char>(ptr[offset]);
    }
    return sum;
}

// -----------------------------------------------------------------------------
// Synonymous Damage Detection
// -----------------------------------------------------------------------------

namespace {

// Standard genetic code - codon index (0-63) to amino acid
constexpr char CODON_TO_AA[64] = {
    'F', 'F', 'L', 'L',  // TTT, TTC, TTA, TTG
    'S', 'S', 'S', 'S',  // TCT, TCC, TCA, TCG
    'Y', 'Y', '*', '*',  // TAT, TAC, TAA, TAG
    'C', 'C', '*', 'W',  // TGT, TGC, TGA, TGG
    'L', 'L', 'L', 'L',  // CTT, CTC, CTA, CTG
    'P', 'P', 'P', 'P',  // CCT, CCC, CCA, CCG
    'H', 'H', 'Q', 'Q',  // CAT, CAC, CAA, CAG
    'R', 'R', 'R', 'R',  // CGT, CGC, CGA, CGG
    'I', 'I', 'I', 'M',  // ATT, ATC, ATA, ATG
    'T', 'T', 'T', 'T',  // ACT, ACC, ACA, ACG
    'N', 'N', 'K', 'K',  // AAT, AAC, AAA, AAG
    'S', 'S', 'R', 'R',  // AGT, AGC, AGA, AGG
    'V', 'V', 'V', 'V',  // GTT, GTC, GTA, GTG
    'A', 'A', 'A', 'A',  // GCT, GCC, GCA, GCG
    'D', 'D', 'E', 'E',  // GAT, GAC, GAA, GAG
    'G', 'G', 'G', 'G',  // GGT, GGC, GGA, GGG
};

// Check if C→T at position in codon is synonymous
// codon_idx: 0-63 observed codon
// nt_pos: 0, 1, or 2 (position within codon)
// Returns true if changing T back to C would give same amino acid
bool is_ct_synonymous(uint8_t codon_idx, int nt_pos) {
    if (codon_idx > 63) return false;

    // Get observed nucleotide at position
    int shift = 4 - 2 * nt_pos;  // Position 0 is bits 4-5, pos 1 is bits 2-3, pos 2 is bits 0-1
    int observed_nt = (codon_idx >> shift) & 3;  // 0=T, 1=C, 2=A, 3=G

    // If not T, no C→T damage possible
    if (observed_nt != 0) return false;

    // Compute the alternative codon with C instead of T
    int alt_codon_idx = codon_idx ^ (1 << shift);  // T(0) XOR 1 = C(1)

    // Compare amino acids
    return CODON_TO_AA[codon_idx] == CODON_TO_AA[alt_codon_idx];
}

// Check if G→A at position in codon is synonymous
bool is_ga_synonymous(uint8_t codon_idx, int nt_pos) {
    if (codon_idx > 63) return false;

    int shift = 4 - 2 * nt_pos;
    int observed_nt = (codon_idx >> shift) & 3;

    // If not A, no G→A damage possible
    if (observed_nt != 2) return false;

    // Compute alternative with G instead of A
    int alt_codon_idx = codon_idx ^ (1 << shift);  // A(2) XOR 1 = G(3)

    return CODON_TO_AA[codon_idx] == CODON_TO_AA[alt_codon_idx];
}

}  // anonymous namespace

SynonymousDamageResult detect_synonymous_damage(
    const AgdRecord& rec,
    float d_max,
    float lambda) {

    SynonymousDamageResult result;

    // Damage probability threshold (only consider sites with p > 0.05)
    constexpr float P_THRESHOLD = 0.05f;

    // Check 5' terminal codons for C→T damage
    for (int i = 0; i < rec.n_5prime; ++i) {
        uint8_t codon_idx = rec.codons_5prime[i];
        if (codon_idx > 63) continue;

        // Position in read (codon i covers nucleotides 3i to 3i+2)
        // Damage probability decreases exponentially from terminus
        for (int nt_pos = 0; nt_pos < 3; ++nt_pos) {
            int read_pos = i * 3 + nt_pos;
            float p_damage = d_max * std::exp(-lambda * read_pos);

            if (p_damage < P_THRESHOLD) continue;

            // Check for C→T (T observed, could have been C)
            int shift = 4 - 2 * nt_pos;
            int observed_nt = (codon_idx >> shift) & 3;

            if (observed_nt == 0) {  // T observed
                bool synonymous = is_ct_synonymous(codon_idx, nt_pos);

                SynonymousDamageResult::DamageSite site;
                site.codon_idx = i;
                site.nt_position = nt_pos;
                site.observed_nt = 'T';
                site.expected_nt = 'C';
                site.is_synonymous = synonymous;
                result.sites.push_back(site);

                if (synonymous) {
                    result.synonymous_5prime++;
                    result.has_synonymous_damage = true;
                } else {
                    result.nonsynonymous_5prime++;
                }
            }
        }
    }

    // Check 3' terminal codons for G→A damage
    // Note: 3' codons are stored from the end, so codon 0 is the last codon
    for (int i = 0; i < rec.n_3prime; ++i) {
        uint8_t codon_idx = rec.codons_3prime[i];
        if (codon_idx > 63) continue;

        // Position from 3' end
        for (int nt_pos = 0; nt_pos < 3; ++nt_pos) {
            // Distance from 3' terminus
            int dist_from_3prime = i * 3 + (2 - nt_pos);
            float p_damage = d_max * std::exp(-lambda * dist_from_3prime);

            if (p_damage < P_THRESHOLD) continue;

            int shift = 4 - 2 * nt_pos;
            int observed_nt = (codon_idx >> shift) & 3;

            if (observed_nt == 2) {  // A observed
                bool synonymous = is_ga_synonymous(codon_idx, nt_pos);

                SynonymousDamageResult::DamageSite site;
                site.codon_idx = i;
                site.nt_position = nt_pos;
                site.observed_nt = 'A';
                site.expected_nt = 'G';
                site.is_synonymous = synonymous;
                result.sites.push_back(site);

                if (synonymous) {
                    result.synonymous_3prime++;
                    result.has_synonymous_damage = true;
                } else {
                    result.nonsynonymous_3prime++;
                }
            }
        }
    }

    return result;
}

}  // namespace agp
