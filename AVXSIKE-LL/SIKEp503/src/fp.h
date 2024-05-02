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

// radix-2^51 for the field elements: 10 * 51-bit = 510-bit
#define VNWORDS 10
#define VBRADIX 51
#define VBMASK  0x7FFFFFFFFFFFFULL
#define VGWORDS 5

static const uint64_t vp503[VNWORDS] = {
  0x7FFFFFFFFFFFF, 0x7FFFFFFFFFFFF, 0x7FFFFFFFFFFFF, 0x7FFFFFFFFFFFF,
  0x2BFFFFFFFFFFF, 0x0B7B44423CF41, 0x21EDF9F6BC4C2, 0x3BD2680DCDFB6,
  0x61E6045C6BDDA, 0x0080CDEA83023, };

static const uint64_t vp503p1[VNWORDS] = {
  0x0000000000000, 0x0000000000000, 0x0000000000000, 0x0000000000000, 
  0x2C00000000000, 0x0B7B44423CF41, 0x21EDF9F6BC4C2, 0x3BD2680DCDFB6, 
  0x61E6045C6BDDA, 0x0080CDEA83023, };

static const uint64_t vp503x2[VNWORDS] = {
  0x7FFFFFFFFFFFE, 0x7FFFFFFFFFFFF, 0x7FFFFFFFFFFFF, 0x7FFFFFFFFFFFF,
  0x57FFFFFFFFFFF, 0x16F6888479E82, 0x43DBF3ED78984, 0x77A4D01B9BF6C,
  0x43CC08B8D7BB4, 0x01019BD506047, };

static const uint64_t vp503x4[VNWORDS] = {
  0x7FFFFFFFFFFFC, 0x7FFFFFFFFFFFF, 0x7FFFFFFFFFFFF, 0x7FFFFFFFFFFFF,
  0x2FFFFFFFFFFFF, 0x2DED1108F3D05, 0x07B7E7DAF1308, 0x6F49A03737ED9,
  0x07981171AF769, 0x020337AA0C08F, };

static const uint64_t vmont_R[VNWORDS] = {
  0x00000000000FE, 0x0000000000000, 0x0000000000000, 0x0000000000000,  
  0x5800000000000, 0x1BB2464785D2A, 0x55E1FD312C76D, 0x253CC24DA0928, 
  0x5DC7AC4CFA13D, 0x0033B15203C83, };

static const uint64_t vmont_R2[VNWORDS] = {
  0x09A0CF641D011, 0x2E313FDA572A5, 0x3723C5EA6E209, 0x347651D9B2EAC,
  0x56FC57AB6EFF1, 0x1DF6A7AEB3AD1, 0x62BA292CF5EEF, 0x58BBAF42ED687,
  0x5B38EAFFB79B9, 0x0080AEAB77E91, }; 

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

// radix-2^64 for the field elements: 8 * 64-bit = 512-bit (from PQCrypto-SIDH-3.4)
#define NWORDS_FIELD    8
#define RADIX           64
#define LOG2RADIX       6  
typedef uint64_t        digit_t;        // unsigned 64-bit digit
typedef uint32_t        hdigit_t;       // unsigned 32-bit digit

typedef unsigned uint128_t __attribute__((mode(TI)));

extern const uint64_t p503[NWORDS_FIELD];
extern const uint64_t p503x2[NWORDS_FIELD];
extern const uint64_t p503x4[NWORDS_FIELD];
extern const uint64_t p503p1[NWORDS_FIELD];
extern const uint64_t Montgomery_one[NWORDS_FIELD];
extern const uint64_t p503p1x64[NWORDS_FIELD/2];

// Montgomery constant Montgomery_R2 = (2^512)^2 mod p503
static const uint64_t Montgomery_R2[NWORDS_FIELD] = { 
  0x5289A0CF641D011F, 0x9B88257189FED2B9, 0xA3B365D58DC8F17A, 0x5BC57AB6EFF168EC,
  0x9E51998BD84D4423, 0xBF8999CBAC3B5695, 0x46E9127BCE14CDB6, 0x003F6CFCE8B81771, };                                                   

// x64 prototypes

extern void mp_sub503_p4_asm(const digit_t* a, const digit_t* b, digit_t* c); 
extern void fpadd503_asm(const digit_t* a, const digit_t* b, digit_t* c);
extern void fpsub503_asm(const digit_t* a, const digit_t* b, digit_t* c);
extern void rdc503_asm(digit_t* ma, digit_t* mc);
extern void mul503_asm(const digit_t* a, const digit_t* b, digit_t* c);
extern void mp_subadd503x2_asm(const digit_t* a, const digit_t* b, digit_t* c);
extern void mp_dblsub503x2_asm(const digit_t* a, const digit_t* b, digit_t* c);

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
