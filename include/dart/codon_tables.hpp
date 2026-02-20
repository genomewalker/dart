#pragma once

#include <unordered_map>
#include <string>
#include <vector>

namespace dart {

/**
 * Fast inline character operations
 * Using lookup tables instead of branches/function calls
 */

// Fast uppercase - branchless lookup table
inline char fast_upper(char c) {
    // Uppercase offset: 'a'-'A' = 32
    return (c >= 'a' && c <= 'z') ? (c - 32) : c;
}

// Fast complement - direct lookup
inline char fast_complement(char c) {
    switch (c) {
        case 'A': return 'T';
        case 'T': return 'A';
        case 'C': return 'G';
        case 'G': return 'C';
        case 'a': return 't';
        case 't': return 'a';
        case 'c': return 'g';
        case 'g': return 'c';
        default: return 'N';
    }
}

// Fast base to index (for array lookups)
inline int fast_base_idx(char c) {
    switch (c) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
        default: return -1;
    }
}

/**
 * Static codon and amino acid frequency tables
 * E. coli RSCU values and UniProt SwissProt frequencies
 */

// Fast array-based lookups (preferred for hot paths)
float get_codon_freq_fast(char c1, char c2, char c3);
char translate_codon_fast(char c1, char c2, char c3);
float get_aa_freq_fast(char aa);

/**
 * Canonical reverse complement functions
 * Use these instead of reimplementing in each file.
 */

// Thread-local buffer for cached reverse complement
inline thread_local std::string g_rc_buffer;

// Allocating version - for one-off uses
inline std::string reverse_complement(const std::string& seq) {
    std::string rc;
    rc.reserve(seq.length());
    for (auto it = seq.rbegin(); it != seq.rend(); ++it) {
        rc += fast_complement(*it);
    }
    return rc;
}

// Cached version - for hot paths (returns reference to thread-local buffer)
// WARNING: Only one cached result is valid at a time per thread
inline const std::string& reverse_complement_cached(const std::string& seq) {
    g_rc_buffer.resize(seq.length());
    for (size_t i = 0; i < seq.length(); ++i) {
        g_rc_buffer[i] = fast_complement(seq[seq.length() - 1 - i]);
    }
    return g_rc_buffer;
}

} // namespace dart
