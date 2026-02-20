#pragma once
/**
 * @file damage_index_writer.hpp
 * @brief Writer for .agd binary damage index files.
 *
 * Usage:
 *   DamageIndexWriter writer("output.agd", profile);
 *   for (const auto& gene : genes) {
 *       writer.add_record(read_id, gene, dna_sequence);
 *   }
 *   writer.finalize();  // Writes hash table and closes file
 */

#include "damage_index.hpp"
#include "types.hpp"
#include "frame_selector.hpp"

#include <fstream>
#include <string>
#include <vector>

namespace dart {

/**
 * @brief Writes .agd binary damage index files.
 *
 * Records are buffered in memory, then written with hash table on finalize().
 * Memory usage: ~32 bytes per gene + hash table overhead.
 */
class DamageIndexWriter {
public:
    /**
     * @brief Construct writer for given output path.
     * @param path Output .agd file path
     * @param profile Sample damage profile (for d_max, lambda, library_type)
     */
    DamageIndexWriter(const std::string& path, const SampleDamageProfile& profile);

    /**
     * @brief Add a gene record to the index.
     * @param read_id Original read ID (AGP suffix will be stripped)
     * @param gene Gene with damage metadata
     * @param dna_sequence Original DNA sequence (for extracting terminal codons)
     */
    void add_record(std::string_view read_id, const Gene& gene, std::string_view dna_sequence);

    /**
     * @brief Finalize and write the index file.
     *
     * Builds hash table, writes header + buckets + records.
     * Must be called before destruction.
     */
    void finalize();

    /**
     * @brief Get number of records added.
     */
    size_t record_count() const { return records_.size(); }

private:
    std::string path_;
    AgdHeader header_;
    std::vector<AgdRecord> records_;
    bool finalized_ = false;

    // Extract terminal codons from DNA sequence
    void extract_terminal_codons(std::string_view dna, int frame, bool is_reverse,
                                  AgdRecord& rec);

    // Pack raw nucleotides from terminal region
    void pack_terminal_nucleotides(std::string_view dna, bool is_reverse,
                                    AgdRecord& rec);
};

}  // namespace dart
