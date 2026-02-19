#pragma once

#include <string>
#include <cstdint>
#include <cstring>

#ifdef USE_AVX512
#include <immintrin.h>
#endif

#ifdef USE_AVX2
#include <immintrin.h>
#endif

namespace agp {
namespace simd {

#ifdef USE_AVX512

// AVX-512 reverse complement - processes 64 bytes at a time
inline void reverse_complement_avx512(const char* src, char* dst, size_t len) {
    // Complement lookup: A<->T (0x41<->0x54), C<->G (0x43<->0x47)
    // For uppercase: A=0x41, C=0x43, G=0x47, T=0x54

    const __m512i mask_A = _mm512_set1_epi8('A');
    const __m512i mask_T = _mm512_set1_epi8('T');
    const __m512i mask_C = _mm512_set1_epi8('C');
    const __m512i mask_G = _mm512_set1_epi8('G');
    const __m512i val_T = _mm512_set1_epi8('T');
    const __m512i val_A = _mm512_set1_epi8('A');
    const __m512i val_G = _mm512_set1_epi8('G');
    const __m512i val_C = _mm512_set1_epi8('C');
    const __m512i val_N = _mm512_set1_epi8('N');

    // Reverse index for shuffling
    const __m512i reverse_idx = _mm512_set_epi8(
        0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15,
        16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31,
        32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47,
        48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63
    );

    size_t i = 0;
    size_t dst_pos = len;

    // Process 64-byte chunks from end to beginning
    while (i + 64 <= len) {
        __m512i chunk = _mm512_loadu_si512(reinterpret_cast<const __m512i*>(src + len - i - 64));

        // Reverse the bytes
        chunk = _mm512_permutexvar_epi8(reverse_idx, chunk);

        // Complement: A->T, T->A, C->G, G->C, else->N
        __mmask64 is_A = _mm512_cmpeq_epi8_mask(chunk, mask_A);
        __mmask64 is_T = _mm512_cmpeq_epi8_mask(chunk, mask_T);
        __mmask64 is_C = _mm512_cmpeq_epi8_mask(chunk, mask_C);
        __mmask64 is_G = _mm512_cmpeq_epi8_mask(chunk, mask_G);

        __m512i result = val_N;  // Default to N
        result = _mm512_mask_mov_epi8(result, is_A, val_T);
        result = _mm512_mask_mov_epi8(result, is_T, val_A);
        result = _mm512_mask_mov_epi8(result, is_C, val_G);
        result = _mm512_mask_mov_epi8(result, is_G, val_C);

        dst_pos -= 64;
        _mm512_storeu_si512(reinterpret_cast<__m512i*>(dst + i), result);
        i += 64;
    }

    // Handle remaining bytes with scalar code
    while (i < len) {
        char c = src[len - 1 - i];
        char comp;
        switch (c) {
            case 'A': case 'a': comp = 'T'; break;
            case 'T': case 't': comp = 'A'; break;
            case 'C': case 'c': comp = 'G'; break;
            case 'G': case 'g': comp = 'C'; break;
            default: comp = 'N'; break;
        }
        dst[i] = comp;
        i++;
    }
}

// AVX-512 uppercase conversion - processes 64 bytes at a time
inline void to_uppercase_avx512(char* data, size_t len) {
    const __m512i lower_a = _mm512_set1_epi8('a');
    const __m512i lower_z = _mm512_set1_epi8('z');
    const __m512i diff = _mm512_set1_epi8('a' - 'A');

    size_t i = 0;
    for (; i + 64 <= len; i += 64) {
        __m512i chunk = _mm512_loadu_si512(reinterpret_cast<const __m512i*>(data + i));

        // Find lowercase letters (a-z)
        __mmask64 is_lower = _mm512_cmpge_epi8_mask(chunk, lower_a) &
                            _mm512_cmple_epi8_mask(chunk, lower_z);

        // Subtract 32 to convert to uppercase
        __m512i upper = _mm512_sub_epi8(chunk, diff);
        chunk = _mm512_mask_mov_epi8(chunk, is_lower, upper);

        _mm512_storeu_si512(reinterpret_cast<__m512i*>(data + i), chunk);
    }

    // Scalar remainder
    for (; i < len; i++) {
        if (data[i] >= 'a' && data[i] <= 'z') {
            data[i] -= 32;
        }
    }
}

// AVX-512 base counting
inline void count_bases_avx512(const char* seq, size_t len,
                               size_t& count_A, size_t& count_C,
                               size_t& count_G, size_t& count_T) {
    const __m512i mask_A = _mm512_set1_epi8('A');
    const __m512i mask_C = _mm512_set1_epi8('C');
    const __m512i mask_G = _mm512_set1_epi8('G');
    const __m512i mask_T = _mm512_set1_epi8('T');

    size_t a = 0, c = 0, g = 0, t = 0;
    size_t i = 0;

    for (; i + 64 <= len; i += 64) {
        __m512i chunk = _mm512_loadu_si512(reinterpret_cast<const __m512i*>(seq + i));

        a += _mm_popcnt_u64(_mm512_cmpeq_epi8_mask(chunk, mask_A));
        c += _mm_popcnt_u64(_mm512_cmpeq_epi8_mask(chunk, mask_C));
        g += _mm_popcnt_u64(_mm512_cmpeq_epi8_mask(chunk, mask_G));
        t += _mm_popcnt_u64(_mm512_cmpeq_epi8_mask(chunk, mask_T));
    }

    // Scalar remainder
    for (; i < len; i++) {
        switch (seq[i]) {
            case 'A': case 'a': a++; break;
            case 'C': case 'c': c++; break;
            case 'G': case 'g': g++; break;
            case 'T': case 't': t++; break;
        }
    }

    count_A = a;
    count_C = c;
    count_G = g;
    count_T = t;
}

#elif defined(USE_AVX2)

// AVX2 fallback implementations (32 bytes at a time)
inline void reverse_complement_avx2(const char* src, char* dst, size_t len) {
    // Scalar implementation
    for (size_t i = 0; i < len; i++) {
        char c = src[len - 1 - i];
        char comp;
        switch (c) {
            case 'A': case 'a': comp = 'T'; break;
            case 'T': case 't': comp = 'A'; break;
            case 'C': case 'c': comp = 'G'; break;
            case 'G': case 'g': comp = 'C'; break;
            default: comp = 'N'; break;
        }
        dst[i] = comp;
    }
}

inline void to_uppercase_avx2(char* data, size_t len) {
    for (size_t i = 0; i < len; i++) {
        if (data[i] >= 'a' && data[i] <= 'z') {
            data[i] -= 32;
        }
    }
}

inline void count_bases_avx2(const char* seq, size_t len,
                             size_t& count_A, size_t& count_C,
                             size_t& count_G, size_t& count_T) {
    size_t a = 0, c = 0, g = 0, t = 0;
    for (size_t i = 0; i < len; i++) {
        switch (seq[i]) {
            case 'A': case 'a': a++; break;
            case 'C': case 'c': c++; break;
            case 'G': case 'g': g++; break;
            case 'T': case 't': t++; break;
        }
    }
    count_A = a; count_C = c; count_G = g; count_T = t;
}

#endif

// Dispatch functions that select best implementation
inline std::string reverse_complement(const std::string& seq) {
    std::string result(seq.length(), 'N');

#ifdef USE_AVX512
    reverse_complement_avx512(seq.data(), result.data(), seq.length());
#elif defined(USE_AVX2)
    reverse_complement_avx2(seq.data(), result.data(), seq.length());
#else
    // Scalar fallback
    for (size_t i = 0; i < seq.length(); i++) {
        char c = seq[seq.length() - 1 - i];
        switch (c) {
            case 'A': case 'a': result[i] = 'T'; break;
            case 'T': case 't': result[i] = 'A'; break;
            case 'C': case 'c': result[i] = 'G'; break;
            case 'G': case 'g': result[i] = 'C'; break;
            default: result[i] = 'N'; break;
        }
    }
#endif
    return result;
}

inline void to_uppercase_inplace(std::string& seq) {
#ifdef USE_AVX512
    to_uppercase_avx512(seq.data(), seq.length());
#elif defined(USE_AVX2)
    to_uppercase_avx2(seq.data(), seq.length());
#else
    for (char& c : seq) {
        if (c >= 'a' && c <= 'z') c -= 32;
    }
#endif
}

} // namespace simd
} // namespace agp
