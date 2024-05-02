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

// radix-2^51 for the field elements: 9 * 51-bit = 459-bit
#define VNWORDS 9
#define VBRADIX 51
#define VBMASK  0x7FFFFFFFFFFFFULL
#define VGWORDS 5

static const uint64_t vp434[VNWORDS] = {
  0x7FFFFFFFFFFFF, 0x7FFFFFFFFFFFF, 0x7FFFFFFFFFFFF, 
  0x7FFFFFFFFFFFF, 0x7DC1767AE2FFF, 0x4B8F062B15D47, 
  0x5A07148159EF1, 0x0BB9A2367E2FE, 0x0000002341F27, };

static const uint64_t vp434p1[VNWORDS] = {
  0x0000000000000, 0x0000000000000, 0x0000000000000, 
  0x0000000000000, 0x7DC1767AE3000, 0x4B8F062B15D47, 
  0x5A07148159EF1, 0x0BB9A2367E2FE, 0x0000002341F27, };

static const uint64_t vp434x2[VNWORDS] = {
  0x7FFFFFFFFFFFE, 0x7FFFFFFFFFFFF, 0x7FFFFFFFFFFFF, 
  0x7FFFFFFFFFFFF, 0x7B82ECF5C5FFF, 0x171E0C562BA8F, 
  0x340E2902B3DE3, 0x1773446CFC5FD, 0x0000004683E4E, };

static const uint64_t vp434x4[VNWORDS] = {
  0x7FFFFFFFFFFFC, 0x7FFFFFFFFFFFF, 0x7FFFFFFFFFFFF, 
  0x7FFFFFFFFFFFF, 0x7705D9EB8BFFF, 0x2E3C18AC5751F, 
  0x681C520567BC6, 0x2EE688D9F8BFA, 0x0000008D07C9C, };

static const uint64_t vmont_R[VNWORDS] = {
  0x0000003A1635C, 0x0000000000000, 0x0000000000000, 
  0x0000000000000, 0x09B6230D6C000, 0x26DCE6EF7EDC0, 
  0x7F1F3F10E0FB2, 0x57CEB5B968572, 0x00000005C9696, };

static const uint64_t vmont_R2[VNWORDS] = {
  0x735A6CC0445BB, 0x4C611472ADB2E, 0x5ACEC73677687, 
  0x4F2E7F0622D11, 0x50B7242C31EAC, 0x18222E43835F8, 
  0x6FE26C5709746, 0x2B1006BB74802, 0x0000001DECB9D, }; 

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

// radix-2^64 for the field elements: 7 * 64-bit = 448-bit (from PQCrypto-SIDH-3.4)
#define NWORDS_FIELD    7
#define RADIX           64
#define LOG2RADIX       6  
typedef uint64_t        digit_t;        // unsigned 64-bit digit
typedef uint32_t        hdigit_t;       // unsigned 32-bit digit

typedef unsigned uint128_t __attribute__((mode(TI)));

extern const uint64_t p434[NWORDS_FIELD];
extern const uint64_t p434x2[NWORDS_FIELD];
extern const uint64_t p434x4[NWORDS_FIELD];
extern const uint64_t p434p1[NWORDS_FIELD];
extern const uint64_t Montgomery_one[NWORDS_FIELD];

// Montgomery constant Montgomery_R2 = (2^448)^2 mod p434
static const uint64_t Montgomery_R2[NWORDS_FIELD] = { 
  0x28E55B65DCD69B30, 0xACEC7367768798C2, 0xAB27973F8311688D, 0x175CC6AF8D6C7C0B,
  0xABCD92BF2DDE347E, 0x69E16A61C7686D9A, 0x000025A89BCDD12A };  

// x64 prototypes

extern void mp_sub434_p4_asm(const digit_t* a, const digit_t* b, digit_t* c); 
extern void fpadd434_asm(const digit_t* a, const digit_t* b, digit_t* c);
extern void fpsub434_asm(const digit_t* a, const digit_t* b, digit_t* c);
extern void rdc434_asm(digit_t* ma, digit_t* mc);
extern void mul434_asm(const digit_t* a, const digit_t* b, digit_t* c);
extern void mp_subadd434x2_asm(const digit_t* a, const digit_t* b, digit_t* c);
extern void mp_dblsub434x2_asm(const digit_t* a, const digit_t* b, digit_t* c);

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
