/**
 *******************************************************************************
 * @version 0.1
 * @date 2021-12-16
 * @copyright Copyright © 2021 by University of Luxembourg
 * @author Hao Cheng
 *******************************************************************************
 */

#ifndef _INTRIN_H
#define _INTRIN_H

// AVX-512 header file
#include <immintrin.h>

// packed 64-bit integer arithmetics (AVX-512IFMA)
#define VADD(X, Y)            _mm512_add_epi64(X, Y)
#define VMADD(W, X, Y, Z)     _mm512_mask_add_epi64(W, X, Y, Z)
#define VSUB(X, Y)            _mm512_sub_epi64(X, Y)
#define VMSUB(W, X, Y, Z)     _mm512_mask_sub_epi64(W, X, Y, Z)
#define VMUL(X, Y)            _mm512_mul_epu32(X, Y)
#define VMACLO(X, Y, Z)       _mm512_madd52lo_epu64(X, Y, Z)
#define VMACHI(X, Y, Z)       _mm512_madd52hi_epu64(X, Y, Z)

// bitwise logical operations
#define VAND(X, Y)            _mm512_and_si512(X, Y)
#define VMAND(W, X, Y, Z)     _mm512_mask_and_epi64(W, X, Y, Z)
#define VZAND(X, Y, Z)        _mm512_maskz_and_epi64(X, Y, Z)
#define VOR(X, Y)             _mm512_or_si512(X, Y)
#define VMOR(W, X, Y, Z)      _mm512_mask_or_epi64(W, X, Y, Z)
#define VRDCOR(X)             _mm512_reduce_or_epi64(X)
#define VXOR(X, Y)            _mm512_xor_si512(X, Y)
#define VMXOR(X, Y, Z, W)     _mm512_mask_xor_epi64(X, Y, Z, W)
#define VSHR(X, Y)            _mm512_srli_epi64(X, Y)
#define VSHL(X, Y)            _mm512_slli_epi64(X, Y)
#define VSRA(X, Y)            _mm512_srai_epi64(X, Y)
#define VZSRA(X, Y, Z)        _mm512_maskz_srai_epi64(X, Y, Z)

// broadcast operations
#define VZERO                 _mm512_setzero_si512()
#define VSET1(X)              _mm512_set1_epi64(X)
#define VSET(X7, X6, X5, X4, X3, X2, X1, X0) \
                              _mm512_set_epi64(X7, X6, X5, X4, X3, X2, X1, X0)
#define VEXTR256(X, Y)        _mm512_extracti64x4_epi64(X, Y)
#define VEXTR64(X, Y)         _mm256_extract_epi64(X, Y) 

// permutation & shuffle operations
#define VPERM(X, Y)           _mm512_permutex_epi64(X, Y)     // 256-bit lane
#define VZPERM(X, Y, Z)       _mm512_maskz_permutex_epi64(X, Y, Z) 
#define VPERMV(X, Y)          _mm512_permutexvar_epi64(X, Y)  // 512-bit lane
#define VSHUF(X, Y)           _mm512_shuffle_epi32(X, Y)      // 128-bit lane
#define VZSHUF(X, Y, Z)       _mm512_maskz_shuffle_epi32(X, Y, Z)
#define VMBLEND(X, Y, Z)      _mm512_mask_blend_epi64(X, Y, Z)
#define VALIGNR(X, Y, Z)      _mm512_alignr_epi64(X, Y, Z)

#endif
