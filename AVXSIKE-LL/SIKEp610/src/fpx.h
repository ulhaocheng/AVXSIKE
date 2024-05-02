/**
 *******************************************************************************
 * @version 0.1
 * @date 2021-12-16
 * @copyright Copyright © 2021 by University of Luxembourg
 * @author Hao Cheng
 *******************************************************************************
 */

#ifndef _FPX_H
#define _FPX_H

#include "fp.h"
#include "params.h"

int8_t ct_compare(const uint8_t *a, const uint8_t *b, unsigned int len); 
void ct_cmov(uint8_t *r, const uint8_t *a, unsigned int len, int8_t selector); 

// -----------------------------------------------------------------------------
// (8x1x1)-way Fp2 arithmetic and some (8x1)-way Fp arithmetic.

typedef __m512i  vfelm_t[VNWORDS];
typedef __m512i  vdfelm_t[2*VNWORDS];
typedef vfelm_t  vf2elm_t[2]; 

void mp_add_8x1w(__m512i *r, const __m512i *a, const __m512i *b);
void fpmul_mont_8x1w(vfelm_t r, const vfelm_t a, const vfelm_t b);

void mp2_add_8x1x1w(vf2elm_t r, const vf2elm_t a, const vf2elm_t b);
void mp2_sub_p2_8x1x1w(vf2elm_t r, const vf2elm_t a, const vf2elm_t b);
void fp2add_8x1x1w(vf2elm_t r, const vf2elm_t a, const vf2elm_t b);
void fp2sub_8x1x1w(vf2elm_t r, const vf2elm_t a, const vf2elm_t b);
void fp2sqr_mont_8x1x1w(vf2elm_t r, const vf2elm_t a);
void fp2mul_mont_8x1x1w(vf2elm_t r, const vf2elm_t a, const vf2elm_t b);
void fp2div2_8x1x1w(vf2elm_t r, const vf2elm_t a);

// -----------------------------------------------------------------------------
// (4x2x1)-way Fp2 arithmetic.

// a = a1 | a0 -> r = a0 | a1
static void vec_shuflh_8x1w(vfelm_t r, const vfelm_t a)
{
  __m512i a0 = a[0], a1 = a[1], a2  = a[2],  a3  = a[3];
  __m512i a4 = a[4], a5 = a[5], a6  = a[6],  a7  = a[7];
  __m512i a8 = a[8], a9 = a[9], a10 = a[10], a11 = a[11];
  __m512i r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11;

  r0  = VSHUF(a0,  0x4E); r1  = VSHUF(a1,  0x4E); r2  = VSHUF(a2,  0x4E); 
  r3  = VSHUF(a3,  0x4E); r4  = VSHUF(a4,  0x4E); r5  = VSHUF(a5,  0x4E);
  r6  = VSHUF(a6,  0x4E); r7  = VSHUF(a7,  0x4E); r8  = VSHUF(a8,  0x4E); 
  r9  = VSHUF(a9,  0x4E); r10 = VSHUF(a10, 0x4E); r11 = VSHUF(a11, 0x4E);

  r[0] = r0; r[1] = r1; r[2]  = r2;  r[3]  = r3; 
  r[4] = r4; r[5] = r5; r[6]  = r6;  r[7]  = r7; 
  r[8] = r8; r[9] = r9; r[10] = r10; r[11] = r11;
} 

// a = a1 | a0 -> r = a0 | a0
static void vec_shufll_8x1w(vfelm_t r, const vfelm_t a)
{
  __m512i a0 = a[0], a1 = a[1], a2  = a[2],  a3  = a[3];
  __m512i a4 = a[4], a5 = a[5], a6  = a[6],  a7  = a[7];
  __m512i a8 = a[8], a9 = a[9], a10 = a[10], a11 = a[11];
  __m512i r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11;

  r0  = VSHUF(a0,  0x44); r1  = VSHUF(a1,  0x44); r2  = VSHUF(a2,  0x44); 
  r3  = VSHUF(a3,  0x44); r4  = VSHUF(a4,  0x44); r5  = VSHUF(a5,  0x44);
  r6  = VSHUF(a6,  0x44); r7  = VSHUF(a7,  0x44); r8  = VSHUF(a8,  0x44); 
  r9  = VSHUF(a9,  0x44); r10 = VSHUF(a10, 0x44); r11 = VSHUF(a11, 0x44);

  r[0] = r0; r[1] = r1; r[2]  = r2;  r[3]  = r3; 
  r[4] = r4; r[5] = r5; r[6]  = r6;  r[7]  = r7; 
  r[8] = r8; r[9] = r9; r[10] = r10; r[11] = r11;
}

// a = a1 | a0 -> r = a1 | a1
static void vec_shufhh_8x1w(vfelm_t r, const vfelm_t a)
{
  __m512i a0 = a[0], a1 = a[1], a2  = a[2],  a3  = a[3];
  __m512i a4 = a[4], a5 = a[5], a6  = a[6],  a7  = a[7];
  __m512i a8 = a[8], a9 = a[9], a10 = a[10], a11 = a[11];
  __m512i r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11;

  r0  = VSHUF(a0,  0xEE); r1  = VSHUF(a1,  0xEE); r2  = VSHUF(a2,  0xEE); 
  r3  = VSHUF(a3,  0xEE); r4  = VSHUF(a4,  0xEE); r5  = VSHUF(a5,  0xEE);
  r6  = VSHUF(a6,  0xEE); r7  = VSHUF(a7,  0xEE); r8  = VSHUF(a8,  0xEE); 
  r9  = VSHUF(a9,  0xEE); r10 = VSHUF(a10, 0xEE); r11 = VSHUF(a11, 0xEE);

  r[0] = r0; r[1] = r1; r[2]  = r2;  r[3]  = r3; 
  r[4] = r4; r[5] = r5; r[6]  = r6;  r[7]  = r7; 
  r[8] = r8; r[9] = r9; r[10] = r10; r[11] = r11;
}

// a = a1 | a0 -> r = 0 | a1
static void vec_shufzh_8x1w(vfelm_t r, const vfelm_t a)
{
  __m512i a0 = a[0], a1 = a[1], a2  = a[2],  a3  = a[3];
  __m512i a4 = a[4], a5 = a[5], a6  = a[6],  a7  = a[7];
  __m512i a8 = a[8], a9 = a[9], a10 = a[10], a11 = a[11];
  __m512i r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11;

  r0  = VZSHUF(0x3333, a0,  0xEE); r1  = VZSHUF(0x3333, a1,  0xEE);
  r2  = VZSHUF(0x3333, a2,  0xEE); r3  = VZSHUF(0x3333, a3,  0xEE);
  r4  = VZSHUF(0x3333, a4,  0xEE); r5  = VZSHUF(0x3333, a5,  0xEE);
  r6  = VZSHUF(0x3333, a6,  0xEE); r7  = VZSHUF(0x3333, a7,  0xEE);
  r8  = VZSHUF(0x3333, a8,  0xEE); r9  = VZSHUF(0x3333, a9,  0xEE);
  r10 = VZSHUF(0x3333, a10, 0xEE); r11 = VZSHUF(0x3333, a11, 0xEE);

  r[0] = r0; r[1] = r1; r[2]  = r2;  r[3]  = r3; 
  r[4] = r4; r[5] = r5; r[6]  = r6;  r[7]  = r7; 
  r[8] = r8; r[9] = r9; r[10] = r10; r[11] = r11;
}

void fp2mul_mont_4x2x1w(vfelm_t r, const vfelm_t a, const vfelm_t b);
void fp2sqr_mont_4x2x1w(vfelm_t r, const vfelm_t a);

// -----------------------------------------------------------------------------
// (2x2x2)-way Fp2 arithmetic.

typedef __m512i  vgelm_t[VGWORDS];
typedef __m512i  vdgelm_t[VNWORDS+VGWORDS];
typedef vgelm_t  vg2elm_t[2]; 

void mp_add_4x2w(vgelm_t r, const vgelm_t a, const vgelm_t b);
void fp2mul_mont_2x2x2w(vgelm_t r, const vgelm_t a, const vgelm_t b);
void fp2sqr_mont_2x2x2w(vgelm_t r, const vgelm_t a);

// -----------------------------------------------------------------------------
// 1-way x64 Fp2 arithmetic (radix-2^51).

typedef digit_t felm_r51_t[VNWORDS];
typedef digit_t dfelm_r51_t[2*VNWORDS];                            
typedef felm_r51_t  f2elm_r51_t[2]; 

void fpcopy_1w(felm_r51_t r, const felm_r51_t a);
void fp2copy_1w(f2elm_r51_t r, const f2elm_r51_t a);
void fpzero_1w(felm_r51_t r);
void mp_add_1w(felm_r51_t r, const felm_r51_t a, const felm_r51_t b);
void mp2_add_1w(f2elm_r51_t r, const f2elm_r51_t a, const f2elm_r51_t b);

// -----------------------------------------------------------------------------
// 1-way x64 Fp2 arithmetic (radix-2^64).

typedef digit_t felm_t[NWORDS_FIELD];
typedef digit_t dfelm_t[2*NWORDS_FIELD];                            
typedef felm_t  f2elm_t[2]; 

extern void mp_add610_asm(const digit_t* a, const digit_t* b, digit_t* c);
void decode_to_digits(const unsigned char* x, digit_t* dec, int nbytes, int ndigits);
void mp2_add(const f2elm_t a, const f2elm_t b, f2elm_t c);
void fpmul_mont(const felm_t ma, const felm_t mb, felm_t mc);
void fp2mul_mont(const f2elm_t a, const f2elm_t b, f2elm_t c);
void fpcopy(const felm_t a, felm_t c);
void fp2copy(const f2elm_t a, f2elm_t c);
void fp2inv_mont(f2elm_t a);
void fp2_encode(const f2elm_t x, unsigned char *enc);
void fp2_decode(const unsigned char *x, f2elm_t dec);
void fp2sub(const f2elm_t a, const f2elm_t b, f2elm_t c);
void fp2add(const f2elm_t a, const f2elm_t b, f2elm_t c);
void fp2sqr_mont(const f2elm_t a, f2elm_t c);
void mp2_add(const f2elm_t a, const f2elm_t b, f2elm_t c);

#endif
