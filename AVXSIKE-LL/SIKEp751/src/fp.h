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

// radix-2^51 for the field elements: 15 * 51-bit = 765-bit
#define VNWORDS 15
#define VBRADIX 51
#define VBMASK  0x7FFFFFFFFFFFFULL
#define VGWORDS 8

static const uint64_t vp751[VNWORDS] = {
  0x7FFFFFFFFFFFF, 0x7FFFFFFFFFFFF, 0x7FFFFFFFFFFFF, 0x7FFFFFFFFFFFF,
  0x7FFFFFFFFFFFF, 0x7FFFFFFFFFFFF, 0x7FFFFFFFFFFFF, 0x7C3C547757FFF,
  0x476E3EC968549, 0x352B363427EF9, 0x2619F5BAFA1DB, 0x22E592BA40427,
  0x3ADC668562B50, 0x6381C25213F2F, 0x0001BF975507D, };

static const uint64_t vp751p1[VNWORDS] = {
  0x0000000000000, 0x0000000000000, 0x0000000000000, 0x0000000000000, 
  0x0000000000000, 0x0000000000000, 0x0000000000000, 0x7C3C547758000, 
  0x476E3EC968549, 0x352B363427EF9, 0x2619F5BAFA1DB, 0x22E592BA40427, 
  0x3ADC668562B50, 0x6381C25213F2F, 0x0001BF975507D, };

static const uint64_t vp751x2[VNWORDS] = {
  0x7FFFFFFFFFFFE, 0x7FFFFFFFFFFFF, 0x7FFFFFFFFFFFF, 0x7FFFFFFFFFFFF,
  0x7FFFFFFFFFFFF, 0x7FFFFFFFFFFFF, 0x7FFFFFFFFFFFF, 0x7878A8EEAFFFF,
  0x0EDC7D92D0A93, 0x6A566C684FDF3, 0x4C33EB75F43B6, 0x45CB25748084E,
  0x75B8CD0AC56A0, 0x470384A427E5E, 0x00037F2EAA0FB, };

static const uint64_t vp751x4[VNWORDS] = {
  0x7FFFFFFFFFFFC, 0x7FFFFFFFFFFFF, 0x7FFFFFFFFFFFF, 0x7FFFFFFFFFFFF,
  0x7FFFFFFFFFFFF, 0x7FFFFFFFFFFFF, 0x7FFFFFFFFFFFF, 0x70F151DD5FFFF,
  0x1DB8FB25A1527, 0x54ACD8D09FBE6, 0x1867D6EBE876D, 0x0B964AE90109D,
  0x6B719A158AD41, 0x0E0709484FCBD, 0x0006FE5D541F7, };

static const uint64_t vmont_R[VNWORDS] = {
  0x0000000004935, 0x0000000000000, 0x0000000000000, 0x0000000000000,  
  0x0000000000000, 0x0000000000000, 0x0000000000000, 0x136C7B32C8000, 
  0x4645918D44FD5, 0x2B98E7D068C98, 0x358DCEF7BEC40, 0x4FA1831DBEF22, 
  0x7722BD3CF247A, 0x65B95D519629A, 0x00012E6C27835, };

static const uint64_t vmont_R2[VNWORDS] = {
  0x4C1191276B501, 0x0B0D34B229511, 0x72E3FD8EDB010, 0x5CE0CBC6D2828,
  0x47D02FF8820A8, 0x1966544827C3A, 0x71F1EE7FC8149, 0x609FBBDCFFE6B,
  0x4CBAC01F5A7FE, 0x299B961575678, 0x7C1E37997D67B, 0x32A3645EDD616,
  0x2ADB4B9F5ED39, 0x24EB2C434A0E2, 0x0000ABF39813E, }; 

// -----------------------------------------------------------------------------
// (8x1)-way fp arithmetic 

void mp_sub_p2_8x1w(__m512i *r, const __m512i *a, const __m512i *b);
void mp_sub_p4_8x1w(__m512i *r, const __m512i *a, const __m512i *b);
void fpadd_8x1w(__m512i *r, const __m512i *a, const __m512i* b);
void fpsub_8x1w(__m512i *r, const __m512i *a, const __m512i* b);
void mp_mul_8x1w(__m512i *r, const __m512i *a, const __m512i *b);
void rdc_mont_8x1w(__m512i *r, const __m512i *a);

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

// radix-2^64 for the field elements: 12 * 64-bit = 768-bit (from PQCrypto-SIDH-3.4)
#define NWORDS_FIELD    12
#define RADIX           64
#define LOG2RADIX       6  
typedef uint64_t        digit_t;        // unsigned 64-bit digit
typedef uint32_t        hdigit_t;       // unsigned 32-bit digit

typedef unsigned uint128_t __attribute__((mode(TI)));

extern const uint64_t p751[NWORDS_FIELD];
extern const uint64_t p751x2[NWORDS_FIELD];
extern const uint64_t p751x4[NWORDS_FIELD];
extern const uint64_t p751p1[NWORDS_FIELD];
extern const uint64_t Montgomery_one[NWORDS_FIELD];

// Montgomery constant Montgomery_R2 = (2^768)^2 mod p751
static const uint64_t Montgomery_R2[NWORDS_FIELD] = { 
  0x233046449DAD4058, 0xDB010161A696452A, 0x5E36941472E3FD8E, 0xF40BFE2082A2E706, 
  0x4932CCA8904F8751, 0x1F735F1F1EE7FC81, 0xA24F4D80C1048E18, 0xB56C383CCDB607C5, 
  0x441DD47B735F9C90, 0x5673ED2C6A6AC82A, 0x06C905261132294B, 0x000041AD830F1F35, };                                                   

// x64 prototypes

extern void mp_sub751_p4_asm(const digit_t* a, const digit_t* b, digit_t* c); 
extern void fpadd751_asm(const digit_t* a, const digit_t* b, digit_t* c);
extern void fpsub751_asm(const digit_t* a, const digit_t* b, digit_t* c);
extern void rdc751_asm(digit_t* ma, digit_t* mc);
extern void mul751_asm(const digit_t* a, const digit_t* b, digit_t* c);
extern void mp_subadd751x2_asm(const digit_t* a, const digit_t* b, digit_t* c);
extern void mp_dblsub751x2_asm(const digit_t* a, const digit_t* b, digit_t* c);

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
