#pragma once

#include "agp/frame_selector.hpp"
#include <string>
#include <vector>
#include <array>
#include <cstdint>
#include <cstdio>
#include <cstring>

namespace agp {

/**
 * Parser for metaDMG bdamage binary files
 *
 * The bdamage format stores position-specific nucleotide substitution counts:
 * - Magic header: 'bdam1' (5 bytes, optional)
 * - MAXLENGTH (int32): maximum positions tracked
 * - Per-reference records:
 *   - ID (int32): reference ID
 *   - NREADS (int32): number of reads
 *   - Forward cycles: MAXL × 16 floats (substitution matrix for 5' end)
 *   - Reverse cycles: MAXL × 16 floats (substitution matrix for 3' end)
 *
 * Substitution matrix ordering: AA,AC,AG,AT,CA,CC,CG,CT,GA,GC,GG,GT,TA,TC,TG,TT
 * Index 7 = C→T (for 5' damage)
 * Index 8 = G→A (for 3' damage)
 *
 * Usage:
 *   auto result = load_bdamage("/path/to/sample.bdamage.gz");
 *   if (result.success) {
 *       sample_profile = result.profile;
 *   }
 */

struct BdamageResult {
    SampleDamageProfile profile;
    bool success = false;
    std::string error_message;
    size_t total_reads = 0;
    size_t num_references = 0;
};

/**
 * Load a bdamage file and populate a SampleDamageProfile
 *
 * @param filename Path to bdamage file (.bdamage or .bdamage.gz)
 * @return BdamageResult with profile and status
 */
BdamageResult load_bdamage(const std::string& filename);

/**
 * Load a bdamage file with filtering to specific reference IDs
 *
 * @param filename Path to bdamage file
 * @param ref_ids Set of reference IDs to include (empty = all)
 * @return BdamageResult with aggregated profile
 */
BdamageResult load_bdamage(const std::string& filename,
                           const std::vector<int32_t>& ref_ids);

} // namespace agp
