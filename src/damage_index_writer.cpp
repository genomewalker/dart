#include "dart/damage_index_writer.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace dart {

DamageIndexWriter::DamageIndexWriter(const std::string& path,
                                     const SampleDamageProfile& profile)
    : path_(path) {
    // Initialize header with sample-level parameters
    std::memset(&header_, 0, sizeof(header_));
    header_.magic = AGD_MAGIC;
    header_.version = AGD_VERSION;
    header_.d_max = profile.d_max_combined;
    header_.lambda = (profile.lambda_5prime + profile.lambda_3prime) / 2.0f;

    // Map library type enum to byte
    switch (profile.library_type) {
        case SampleDamageProfile::LibraryType::SINGLE_STRANDED:
            header_.library_type = 1;
            break;
        case SampleDamageProfile::LibraryType::DOUBLE_STRANDED:
            header_.library_type = 2;
            break;
        default:
            header_.library_type = 0;  // unknown
            break;
    }

    // v3: damage informativeness fields
    header_.damage_validated = profile.damage_validated ? 1 : 0;
    header_.damage_artifact = profile.damage_artifact ? 1 : 0;
    header_.channel_b_valid = profile.channel_b_valid ? 1 : 0;
    header_.stop_decay_llr = profile.stop_decay_llr_5prime;
    header_.terminal_shift = profile.terminal_shift_5prime;

    // Reserve space for typical dataset
    records_.reserve(100000);
}

void DamageIndexWriter::add_record(std::string_view read_id, const Gene& gene,
                                    std::string_view dna_sequence) {
    if (finalized_) {
        throw std::runtime_error("Cannot add records after finalize()");
    }

    AgdRecord rec;
    std::memset(&rec, 0, sizeof(rec));

    // Hash the read ID (strip DART suffix first)
    std::string_view base_id = strip_dart_suffix(read_id);
    rec.id_hash = fnv1a_hash(base_id);

    // Basic metadata
    rec.seq_len = static_cast<uint16_t>(std::min<size_t>(dna_sequence.size(), 65535));
    rec.frame_strand = encode_frame_strand(gene.frame, !gene.is_forward);
    rec.damage_pct_q = quantize_damage_pct(gene.damage_score);
    rec.p_damaged_q = quantize_probability(gene.p_read_damaged);

    // Extract terminal codons based on frame and strand
    extract_terminal_codons(dna_sequence, gene.frame, !gene.is_forward, rec);

    // Pack raw nucleotides for synonymous damage detection
    pack_terminal_nucleotides(dna_sequence, !gene.is_forward, rec);

    records_.push_back(rec);
}

void DamageIndexWriter::extract_terminal_codons(std::string_view dna, int frame,
                                                 bool is_reverse, AgdRecord& rec) {
    // Initialize as invalid
    std::memset(rec.codons_5prime, 255, 5);
    std::memset(rec.codons_3prime, 255, 5);
    rec.n_5prime = 0;
    rec.n_3prime = 0;

    const size_t len = dna.size();
    if (len < 3) return;

    // For reverse strand, we work with the reverse complement conceptually,
    // but we extract codons from the original sequence in reverse order.

    if (!is_reverse) {
        // Forward strand: 5' end is at start of sequence
        // First codon starts at position 'frame'
        size_t pos = frame;
        for (int i = 0; i < 5 && pos + 3 <= len; ++i, pos += 3) {
            rec.codons_5prime[i] = encode_codon(dna.data() + pos);
            rec.n_5prime = i + 1;
        }

        // 3' end codons: work backwards from the end
        // Find the last complete codon position
        size_t coding_len = ((len - frame) / 3) * 3;
        if (coding_len >= 3) {
            pos = frame + coding_len - 3;
            for (int i = 0; i < 5 && pos >= frame; ++i, pos -= 3) {
                rec.codons_3prime[i] = encode_codon(dna.data() + pos);
                rec.n_3prime = i + 1;
                if (pos < 3) break;  // Prevent underflow
            }
        }
    } else {
        // Reverse strand: 5' end of protein is at 3' end of DNA
        // Need to consider reverse complement
        // For efficiency, we extract codon indices that would result from
        // reverse-complementing the DNA first

        // Calculate the start position for reverse-strand reading
        // The protein's 5' corresponds to the DNA's 3' end
        size_t coding_len = ((len - frame) / 3) * 3;
        if (coding_len < 3) return;

        // 5' codons (from DNA 3' end, reading backwards)
        size_t pos = len - frame - 3;  // Last codon position
        for (int i = 0; i < 5 && pos + 3 <= len; ++i) {
            // Extract and encode the reverse complement codon
            char codon[3];
            for (int j = 0; j < 3; ++j) {
                char nt = dna[pos + 2 - j];  // Read backwards
                // Complement
                switch (nt) {
                    case 'A': case 'a': codon[j] = 'T'; break;
                    case 'T': case 't': codon[j] = 'A'; break;
                    case 'G': case 'g': codon[j] = 'C'; break;
                    case 'C': case 'c': codon[j] = 'G'; break;
                    default: codon[j] = 'N'; break;
                }
            }
            rec.codons_5prime[i] = encode_codon(codon);
            rec.n_5prime = i + 1;

            if (pos < 3) break;
            pos -= 3;
        }

        // 3' codons (from DNA 5' end)
        pos = frame;
        for (int i = 0; i < 5 && pos + 3 <= len; ++i, pos += 3) {
            char codon[3];
            for (int j = 0; j < 3; ++j) {
                char nt = dna[pos + 2 - j];
                switch (nt) {
                    case 'A': case 'a': codon[j] = 'T'; break;
                    case 'T': case 't': codon[j] = 'A'; break;
                    case 'G': case 'g': codon[j] = 'C'; break;
                    case 'C': case 'c': codon[j] = 'G'; break;
                    default: codon[j] = 'N'; break;
                }
            }
            rec.codons_3prime[i] = encode_codon(codon);
            rec.n_3prime = i + 1;
        }
    }
}

void DamageIndexWriter::pack_terminal_nucleotides(std::string_view dna,
                                                   bool is_reverse,
                                                   AgdRecord& rec) {
    // Pack first 12 nucleotides from each terminus (3 bytes each)
    // This allows detecting exact nucleotide changes for synonymous damage
    std::memset(rec.nt_5prime, 0, 3);
    std::memset(rec.nt_3prime, 0, 3);

    const size_t len = dna.size();

    if (!is_reverse) {
        // 5' nucleotides: positions 0-11
        for (size_t i = 0; i < 12 && i < len; ++i) {
            int val = encode_nucleotide(dna[i]);
            if (val < 0) val = 0;
            size_t byte_idx = i / 4;
            size_t bit_pos = 6 - 2 * (i % 4);
            rec.nt_5prime[byte_idx] |= (val << bit_pos);
        }

        // 3' nucleotides: last 12 positions
        for (size_t i = 0; i < 12 && i < len; ++i) {
            size_t pos = len - 12 + i;
            if (pos >= len) continue;
            int val = encode_nucleotide(dna[pos]);
            if (val < 0) val = 0;
            size_t byte_idx = i / 4;
            size_t bit_pos = 6 - 2 * (i % 4);
            rec.nt_3prime[byte_idx] |= (val << bit_pos);
        }
    } else {
        // For reverse strand, 5' of protein = 3' of DNA (complement)
        // 5' nucleotides: complement of last 12 DNA positions, reversed
        for (size_t i = 0; i < 12 && i < len; ++i) {
            size_t pos = len - 1 - i;
            char nt = dna[pos];
            int val;
            switch (nt) {
                case 'A': case 'a': val = 0; break;  // A->T encoded as T
                case 'T': case 't': val = 2; break;  // T->A encoded as A
                case 'G': case 'g': val = 1; break;  // G->C encoded as C
                case 'C': case 'c': val = 3; break;  // C->G encoded as G
                default: val = 0; break;
            }
            size_t byte_idx = i / 4;
            size_t bit_pos = 6 - 2 * (i % 4);
            rec.nt_5prime[byte_idx] |= (val << bit_pos);
        }

        // 3' nucleotides: complement of first 12 DNA positions, reversed
        for (size_t i = 0; i < 12 && i < len; ++i) {
            size_t pos = 11 - i;
            if (pos >= len) continue;
            char nt = dna[pos];
            int val;
            switch (nt) {
                case 'A': case 'a': val = 0; break;
                case 'T': case 't': val = 2; break;
                case 'G': case 'g': val = 1; break;
                case 'C': case 'c': val = 3; break;
                default: val = 0; break;
            }
            size_t byte_idx = i / 4;
            size_t bit_pos = 6 - 2 * (i % 4);
            rec.nt_3prime[byte_idx] |= (val << bit_pos);
        }
    }
}

void DamageIndexWriter::finalize() {
    if (finalized_) return;
    finalized_ = true;

    if (records_.empty()) {
        // Write empty file with just header
        header_.num_records = 0;
        header_.num_buckets = 0;

        std::ofstream out(path_, std::ios::binary);
        if (!out) {
            throw std::runtime_error("Failed to open output file: " + path_);
        }
        out.write(reinterpret_cast<const char*>(&header_), sizeof(header_));
        return;
    }

    // Calculate number of buckets (load factor ~0.75)
    header_.num_records = records_.size();
    header_.num_buckets = static_cast<uint64_t>(records_.size() * 1.33 + 1);

    // Build hash table with chaining
    std::vector<AgdBucket> buckets(header_.num_buckets);
    for (auto& b : buckets) {
        b.record_offset = 0xFFFFFFFF;
    }

    // We need to handle collisions by chaining
    // Store next pointers in a separate vector (will be merged into records)
    std::vector<uint32_t> next_chain(records_.size(), 0xFFFFFFFF);

    for (size_t i = 0; i < records_.size(); ++i) {
        uint64_t bucket_idx = records_[i].id_hash % header_.num_buckets;
        if (buckets[bucket_idx].record_offset == 0xFFFFFFFF) {
            // Empty bucket, insert directly
            buckets[bucket_idx].record_offset = static_cast<uint32_t>(i);
        } else {
            // Collision: chain at head
            next_chain[i] = buckets[bucket_idx].record_offset;
            buckets[bucket_idx].record_offset = static_cast<uint32_t>(i);
        }
    }

    // Write file
    std::ofstream out(path_, std::ios::binary);
    if (!out) {
        throw std::runtime_error("Failed to open output file: " + path_);
    }

    // Write header
    out.write(reinterpret_cast<const char*>(&header_), sizeof(header_));

    // Write hash table buckets
    out.write(reinterpret_cast<const char*>(buckets.data()),
              buckets.size() * sizeof(AgdBucket));

    // Write records
    out.write(reinterpret_cast<const char*>(records_.data()),
              records_.size() * sizeof(AgdRecord));

    // Write chain pointers (as a separate section after records)
    // This allows the record struct to stay at 32 bytes
    out.write(reinterpret_cast<const char*>(next_chain.data()),
              next_chain.size() * sizeof(uint32_t));

    if (!out) {
        throw std::runtime_error("Failed to write output file: " + path_);
    }
}

}  // namespace dart
