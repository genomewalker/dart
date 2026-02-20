#pragma once

#include <string>
#include <vector>

namespace dart {

// Optimized functions (no allocations)
std::string translate_sequence(const std::string& seq, int frame);
float calculate_codon_usage_score_direct(const std::string& seq, int frame);
float calculate_dicodon_score(const std::string& seq, int frame);
float calculate_dipeptide_score(const std::string& protein);

// Legacy function (allocates vector)
std::vector<std::string> extract_codons(const std::string& seq, int frame);

} // namespace dart
