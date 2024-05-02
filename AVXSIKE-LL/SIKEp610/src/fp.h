/**
 *******************************************************************************
 * @version 0.1
 * @date 2021-12-16
 * @copyright Copyright © 2021 by University of Luxembourg
 * @author Hao Cheng
 *******************************************************************************
 */

#ifndef _GFPARITH_H
#define _GFPARITH_H

#include <stdint.h>
#include "intrin.h"

// radix-2^51 for the field elements: 12 * 51-bit = 612-bit
#define VNWORDS 12
#define VBRADIX 51
#define VBMASK  0x7FFFFFFFFFFFFULL
#define VGWORDS 6

static const uint64_t vp610[VNWORDS] = {
  0x7FFFFFFFFFFFF, 0x7FFFFFFFFFFFF, 0x7FFFFFFFFFFFF, 0x7FFFFFFFFFFFF,
  0x7FFFFFFFFFFFF, 0x3FFFFFFFFFFFF, 0x22A96AC0B9B80, 0x47FCD5D8BC26F,
  0x52A9AE7BF4504, 0x64AB65F421884, 0x04309479F6232, 0x13DFB53B440C8, };

static const uint64_t vp610p1[VNWORDS] = {
  0x0000000000000, 0x0000000000000, 0x0000000000000, 0x0000000000000, 
  0x0000000000000, 0x4000000000000, 0x22A96AC0B9B80, 0x47FCD5D8BC26F,
  0x52A9AE7BF4504, 0x64AB65F421884, 0x04309479F6232, 0x13DFB53B440C8, };

static const uint64_t vp610x2[VNWORDS] = {
  0x7FFFFFFFFFFFE, 0x7FFFFFFFFFFFF, 0x7FFFFFFFFFFFF, 0x7FFFFFFFFFFFF,
  0x7FFFFFFFFFFFF, 0x7FFFFFFFFFFFF, 0x4552D58173700, 0x0FF9ABB1784DE,
  0x25535CF7E8A09, 0x4956CBE843109, 0x086128F3EC465, 0x27BF6A7688190, };

static const uint64_t vp610x4[VNWORDS] = {
  0x7FFFFFFFFFFFC, 0x7FFFFFFFFFFFF, 0x7FFFFFFFFFFFF, 0x7FFFFFFFFFFFF,
  0x7FFFFFFFFFFFF, 0x7FFFFFFFFFFFF, 0x0AA5AB02E6E01, 0x1FF35762F09BD,
  0x4AA6B9EFD1412, 0x12AD97D086212, 0x10C251E7D88CB, 0x4F7ED4ED10320, };

static const uint64_t vmont_R[VNWORDS] = {
  0x0000000000006, 0x0000000000000, 0x0000000000000, 0x0000000000000, 
  0x0000000000000, 0x0000000000000, 0x30077F7BA5AFD, 0x5012FCEB97164,
  0x1005E918461E4, 0x23FB9C4736CE4, 0x66DC85243B2CF, 0x08C1C09C67B4F, };

static const uint64_t vmont_R2[VNWORDS] = {
  0x163B627392EE7, 0x31BC927BC170B, 0x68315AF05C1E0, 0x56E3FA0CCA068,
  0x26030979EDE54, 0x744E105718847, 0x42E56CDD585D4, 0x277C0450B087A,
  0x6CBEFFAD0601E, 0x7CBD17B195293, 0x2BF3C820BCD5B, 0x07673F03D2224, }; 

// -----------------------------------------------------------------------------
// (8x1)-way fp arithmetic 

void mp_sub_p2_8x1w(__m512i *r, const __m512i *a, const __m512i *b);
void mp_sub_p4_8x1w(__m512i *r, const __m512i *a, const __m512i *b);
void fpadd_8x1w(__m512i *r, const __m512i *a, const __m512i* b);
void fpsub_8x1w(__m512i *r, const __m512i *a, const __m512i* b);
void rdc_mont_8x1w(__m512i *r, const __m512i *a);

void mp_mul_8x1w_v1(__m512i *r, const __m512i *a, const __m512i *b);
void mp_mul_8x1w_v2(__m512i *r, const __m512i *a, const __m512i *b);

#define mp_mul_8x1w mp_mul_8x1w_v2

// -----------------------------------------------------------------------------
// (4x2)-way fp arithmetic 

void mp_sub_p2_4x2w(__m512i *r, const __m512i *a, const __m512i *b);
void mp_sub_p4_4x2w(__m512i *r, const __m512i *a, const __m512i *b);
void fpadd_4x2w(__m512i *r, const __m512i *a, const __m512i* b);
void fpsub_4x2w(__m512i *r, const __m512i *a, const __m512i* b);
void mp_mul_4x2w(__m512i *r, const __m512i *a, const __m512i *b);
void rdc_mont_4x2w(__m512i *r, const __m512i *a);

// -----------------------------------------------------------------------------
// 1-way x64 fp arithmetic 

// radix-2^64 for the field elements: 10 * 64-bit = 640-bit (from PQCrypto-SIDH-3.4)
#define NWORDS_FIELD    10
#define RADIX           64
#define LOG2RADIX       6  
typedef uint64_t        digit_t;        // unsigned 64-bit digit
typedef uint32_t        hdigit_t;       // unsigned 32-bit digit

typedef unsigned uint128_t __attribute__((mode(TI)));

extern const uint64_t p610[NWORDS_FIELD];
extern const uint64_t p610x2[NWORDS_FIELD];
extern const uint64_t p610x4[NWORDS_FIELD];
extern const uint64_t p610p1[NWORDS_FIELD];
extern const uint64_t Montgomery_one[NWORDS_FIELD];
extern const uint64_t Montgomery_R2[NWORDS_FIELD];
                                              
// x64 prototypes

extern void mp_sub610_p4_asm(const digit_t* a, const digit_t* b, digit_t* c); 
extern void fpadd610_asm(const digit_t* a, const digit_t* b, digit_t* c);
extern void fpsub610_asm(const digit_t* a, const digit_t* b, digit_t* c);
extern void rdc610_asm(digit_t* ma, digit_t* mc);
extern void mul610_asm(const digit_t* a, const digit_t* b, digit_t* c);
extern void mp_subadd610x2_asm(const digit_t* a, const digit_t* b, digit_t* c);
extern void mp_dblsub610x2_asm(const digit_t* a, const digit_t* b, digit_t* c);

void mp_mul(const digit_t* a, const digit_t* b, digit_t* c, const unsigned int nwords);
void rdc_mont(digit_t* ma, digit_t* mc);
void fpcorrection(digit_t* a);
void fpadd(const digit_t* a, const digit_t* b, digit_t* c);
void fpneg(digit_t* a);
void fpsub(const digit_t* a, const digit_t* b, digit_t* c);
void mp_sub_p4(const digit_t* a, const digit_t* b, digit_t* c);

// Digit addition with carry
#define ADDC(carryIn, addend1, addend2, carryOut, sumOut)                                         \
    { uint128_t tempReg = (uint128_t)(addend1) + (uint128_t)(addend2) + (uint128_t)(carryIn);     \
    (carryOut) = (digit_t)(tempReg >> RADIX);                                                     \
    (sumOut) = (digit_t)tempReg; }  
    
// Digit subtraction with borrow
#define SUBC(borrowIn, minuend, subtrahend, borrowOut, differenceOut)                             \
    { uint128_t tempReg = (uint128_t)(minuend) - (uint128_t)(subtrahend) - (uint128_t)(borrowIn); \
    (borrowOut) = (digit_t)(tempReg >> (sizeof(uint128_t)*8 - 1));                                \
    (differenceOut) = (digit_t)tempReg; }

#endif
