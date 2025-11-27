#!/usr/bin/env python3
"""
Create damage likelihood lookup tables from GTDB hexamer frequencies.

For each hexamer, compute:
- P(hexamer | damaged) / P(hexamer | undamaged)
- This becomes a feature for damage detection

Output: C++ header with lookup tables
"""

import json
import sys

def main():
    with open('/scratch/tmp/gtdb_bacteria_hexamers.json') as f:
        data = json.load(f)

    freq = data['hexamer_frequencies']

    # Build damage likelihood lookup for 5' C→T damage
    # A T-starting hexamer at 5' end is damage evidence if its C-precursor is common
    damage_5prime = {}

    for hex_seq, hex_freq in freq.items():
        if hex_seq.startswith('T'):
            c_precursor = 'C' + hex_seq[1:]
            c_freq = freq.get(c_precursor, 0)

            # Log-likelihood ratio: log(P_damaged / P_undamaged)
            # P(see T-hex | damaged) ≈ P(C-hex) * damage_rate
            # P(see T-hex | undamaged) = P(T-hex)
            # Ratio ≈ P(C-hex) / P(T-hex)

            if hex_freq > 0:
                ratio = c_freq / hex_freq
                if ratio > 1.5:  # Only include if informative
                    damage_5prime[hex_seq] = (hex_freq, c_freq, ratio)

    # Build damage likelihood lookup for 3' G→A damage
    # An A-ending hexamer at 3' end is damage evidence if its G-precursor is common
    damage_3prime = {}

    for hex_seq, hex_freq in freq.items():
        if hex_seq.endswith('A'):
            g_precursor = hex_seq[:-1] + 'G'
            g_freq = freq.get(g_precursor, 0)

            if hex_freq > 0:
                ratio = g_freq / hex_freq
                if ratio > 1.5:
                    damage_3prime[hex_seq] = (hex_freq, g_freq, ratio)

    # Sort by ratio (most informative first)
    sorted_5prime = sorted(damage_5prime.items(), key=lambda x: -x[1][2])
    sorted_3prime = sorted(damage_3prime.items(), key=lambda x: -x[1][2])

    print("=== 5' C→T DAMAGE INDICATORS (top 50) ===")
    print("Hexamer | Undamaged freq | Precursor freq | Ratio")
    for hex_seq, (t_freq, c_freq, ratio) in sorted_5prime[:50]:
        print(f"{hex_seq} | {t_freq:.8f} | {c_freq:.8f} | {ratio:.1f}x")

    print("\n=== 3' G→A DAMAGE INDICATORS (top 50) ===")
    print("Hexamer | Undamaged freq | Precursor freq | Ratio")
    for hex_seq, (a_freq, g_freq, ratio) in sorted_3prime[:50]:
        print(f"{hex_seq} | {a_freq:.8f} | {g_freq:.8f} | {ratio:.1f}x")

    # Generate C++ header
    print("\n\n// Generating C++ header...")

    with open('include/agp/gtdb_damage_likelihood.hpp', 'w') as f:
        f.write("// GTDB-derived damage likelihood ratios\n")
        f.write("// Auto-generated from gtdb_bacteria_hexamers.json\n")
        f.write("// DO NOT EDIT MANUALLY\n\n")
        f.write("#pragma once\n\n")
        f.write("#include <cstdint>\n")
        f.write("#include <cstring>\n\n")
        f.write("namespace agp {\n\n")

        # Hexamer encoding function
        f.write("// Encode hexamer as 12-bit integer\n")
        f.write("inline uint32_t encode_hexamer_dmg(const char* seq) {\n")
        f.write("    static const int8_t base_map[256] = {\n")
        f.write("        -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,\n")
        f.write("        -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,\n")
        f.write("        -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,\n")
        f.write("        -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,\n")
        f.write("        -1, 0,-1, 1,-1,-1,-1, 2,-1,-1,-1,-1,-1,-1,-1,-1,\n")
        f.write("        -1,-1,-1,-1, 3,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,\n")
        f.write("        -1, 0,-1, 1,-1,-1,-1, 2,-1,-1,-1,-1,-1,-1,-1,-1,\n")
        f.write("        -1,-1,-1,-1, 3,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1\n")
        f.write("    };\n")
        f.write("    uint32_t code = 0;\n")
        f.write("    for (int i = 0; i < 6; i++) {\n")
        f.write("        int8_t b = base_map[(uint8_t)seq[i]];\n")
        f.write("        if (b < 0) return UINT32_MAX;\n")
        f.write("        code = (code << 2) | b;\n")
        f.write("    }\n")
        f.write("    return code;\n")
        f.write("}\n\n")

        # 5' damage likelihood table (log-ratio, capped at 10)
        f.write("// 5' C→T damage log-likelihood ratios\n")
        f.write("// Index = encode_hexamer_dmg(seq)\n")
        f.write("// Value = log2(P_C_precursor / P_T_form), capped at 12\n")
        f.write("static constexpr float DAMAGE_LLR_5PRIME[4096] = {\n")

        import math
        llr_5prime = [0.0] * 4096
        for hex_seq, (t_freq, c_freq, ratio) in damage_5prime.items():
            code = encode_hexamer(hex_seq)
            if code < 4096:
                llr_5prime[code] = min(12.0, math.log2(ratio)) if ratio > 1 else 0.0

        for i in range(0, 4096, 16):
            values = [f"{llr_5prime[i+j]:.4f}f" for j in range(16)]
            f.write("    " + ", ".join(values) + ",\n")

        f.write("};\n\n")

        # 3' damage likelihood table
        f.write("// 3' G→A damage log-likelihood ratios\n")
        f.write("static constexpr float DAMAGE_LLR_3PRIME[4096] = {\n")

        llr_3prime = [0.0] * 4096
        for hex_seq, (a_freq, g_freq, ratio) in damage_3prime.items():
            code = encode_hexamer(hex_seq)
            if code < 4096:
                llr_3prime[code] = min(12.0, math.log2(ratio)) if ratio > 1 else 0.0

        for i in range(0, 4096, 16):
            values = [f"{llr_3prime[i+j]:.4f}f" for j in range(16)]
            f.write("    " + ", ".join(values) + ",\n")

        f.write("};\n\n")

        # Helper function to compute damage score
        f.write("// Compute damage score for a sequence using GTDB likelihood ratios\n")
        f.write("inline float compute_gtdb_damage_score(const char* seq, size_t len) {\n")
        f.write("    if (len < 6) return 0.0f;\n")
        f.write("    \n")
        f.write("    float score_5prime = 0.0f;\n")
        f.write("    float score_3prime = 0.0f;\n")
        f.write("    \n")
        f.write("    // Check 5' end (first 3 hexamers, position-weighted)\n")
        f.write("    for (int i = 0; i < 3 && i + 6 <= (int)len; i++) {\n")
        f.write("        uint32_t code = encode_hexamer_dmg(seq + i);\n")
        f.write("        if (code < 4096) {\n")
        f.write("            float weight = 1.0f / (1 << i);  // 1.0, 0.5, 0.25\n")
        f.write("            score_5prime += DAMAGE_LLR_5PRIME[code] * weight;\n")
        f.write("        }\n")
        f.write("    }\n")
        f.write("    \n")
        f.write("    // Check 3' end (last 3 hexamers, position-weighted)\n")
        f.write("    for (int i = 0; i < 3 && i + 6 <= (int)len; i++) {\n")
        f.write("        size_t pos = len - 6 - i;\n")
        f.write("        uint32_t code = encode_hexamer_dmg(seq + pos);\n")
        f.write("        if (code < 4096) {\n")
        f.write("            float weight = 1.0f / (1 << i);\n")
        f.write("            score_3prime += DAMAGE_LLR_3PRIME[code] * weight;\n")
        f.write("        }\n")
        f.write("    }\n")
        f.write("    \n")
        f.write("    // Combine scores (max ~24 if all positions show damage)\n")
        f.write("    return score_5prime + score_3prime;\n")
        f.write("}\n\n")

        f.write("} // namespace agp\n")

    print("Generated include/agp/gtdb_damage_likelihood.hpp")


def encode_hexamer(seq):
    """Encode hexamer to integer code."""
    base_map = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    code = 0
    for c in seq.upper():
        if c not in base_map:
            return 4096  # Invalid
        code = (code << 2) | base_map[c]
    return code


if __name__ == "__main__":
    main()
