/**
 *******************************************************************************
 * @version 0.1
 * @date 2021-12-16
 * @copyright Copyright © 2021 by University of Luxembourg
 * @author Hao Cheng
 *******************************************************************************
 */

#include "fpx.h"

int8_t ct_compare(const uint8_t *a, const uint8_t *b, unsigned int len) 
{ // Compare two byte arrays in constant time.
  // Returns 0 if the byte arrays are equal, -1 otherwise.
    uint8_t r = 0;

    for (unsigned int i = 0; i < len; i++)
        r |= a[i] ^ b[i];

    return (-(int8_t)r) >> (8*sizeof(uint8_t)-1);
}

void ct_cmov(uint8_t *r, const uint8_t *a, unsigned int len, int8_t selector) 
{ // Conditional move in constant time.
  // If selector = -1 then load r with a, else if selector = 0 then keep r.

    for (unsigned int i = 0; i < len; i++)
        r[i] ^= selector & (a[i] ^ r[i]);
}

// r = a
void fpcopy(felm_t r, const felm_t a)
{
  int i;

  for (i = 0; i < NWORDS; i++) r[i] = a[i];
}

// r = 0
void fpzero(felm_t r)
{
  int i;

  for (i = 0; i < NWORDS; i++) r[i] = VZERO;
}

// convert from number domain to Montgomery domain
void to_mont(felm_t r, const felm_t a)
{
  felm_t R2;
  int i;

  for (i = 0; i < NWORDS; i++) R2[i] = VSET1(mont_R2[i]);
  fpmul_mont(r, a, R2);
}

// convert from Montgomery domain to number domain 
void from_mont(felm_t r, const felm_t a)
{
  felm_t one;
  int i;

  one[0] = VSET1(1);
  for (i = 1; i < NWORDS; i++) one[i] = VZERO;
  fpmul_mont(r, a, one);
  fpcorrection(r);
}

// Montgomery multplication (to be optimized)
void fpmul_mont_v1(felm_t r, const felm_t a, const felm_t b)
{
  dfelm_t temp;

  mp_mul(temp, a, b);
  rdc_mont(r, temp);
}

// optimized version

// Montgomery squaring (to be optimized)
void fpsqr_mont_v1(felm_t r, const felm_t a)
{
  dfelm_t temp;

  mp_mul(temp, a, a);
  rdc_mont(r, temp);
}

// optimized version
void fpsqr_mont_v2(felm_t r, const felm_t a)
{
  dfelm_t temp;

  mp_sqr(temp, a);
  rdc_mont(r, temp);
}

void fpinv_chain_mont(felm_t r)
{
  felm_t t[27], tt;
  int i, j;
    
  // Precomputed table
  fpsqr_mont(tt, r);
  fpmul_mont(t[0], r, tt);
  fpmul_mont(t[1], t[0], tt);
  fpmul_mont(t[2], t[1], tt);
  fpmul_mont(t[3], t[2], tt); 
  fpmul_mont(t[3], t[3], tt);
  for (i = 3; i <= 8; i++) fpmul_mont(t[i+1], t[i], tt);
  fpmul_mont(t[9], t[9], tt);
  for (i = 9; i <= 20; i++) fpmul_mont(t[i+1], t[i], tt);    
  fpmul_mont(t[21], t[21], tt); 
  for (i = 21; i <= 24; i++) fpmul_mont(t[i+1], t[i], tt); 
  fpmul_mont(t[25], t[25], tt);
  fpmul_mont(t[26], t[25], tt);

  fpcopy(tt, r);
  for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
  fpmul_mont(tt, t[20], tt);
  for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
  fpmul_mont(tt, t[24], tt);
  for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
  fpmul_mont(tt, t[11], tt);
  for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
  fpmul_mont(tt, t[8], tt);
  for (i = 0; i < 8; i++) fpsqr_mont(tt, tt);
  fpmul_mont(tt, t[2], tt);
  for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
  fpmul_mont(tt, t[23], tt);
  for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
  fpmul_mont(tt, t[2], tt);
  for (i = 0; i < 9; i++) fpsqr_mont(tt, tt);
  fpmul_mont(tt, t[2], tt);
  for (i = 0; i < 10; i++) fpsqr_mont(tt, tt);
  fpmul_mont(tt, t[15], tt);
  for (i = 0; i < 8; i++) fpsqr_mont(tt, tt);
  fpmul_mont(tt, t[13], tt);
  for (i = 0; i < 8; i++) fpsqr_mont(tt, tt);
  fpmul_mont(tt, t[26], tt);
  for (i = 0; i < 8; i++) fpsqr_mont(tt, tt);
  fpmul_mont(tt, t[20], tt);
  for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
  fpmul_mont(tt, t[11], tt);
  for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
  fpmul_mont(tt, t[10], tt);
  for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
  fpmul_mont(tt, t[14], tt);
  for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
  fpmul_mont(tt, t[4], tt);
  for (i = 0; i < 10; i++) fpsqr_mont(tt, tt);
  fpmul_mont(tt, t[18], tt);
  for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
  fpmul_mont(tt, t[1], tt);
  for (i = 0; i < 7; i++) fpsqr_mont(tt, tt);
  fpmul_mont(tt, t[22], tt);
  for (i = 0; i < 10; i++) fpsqr_mont(tt, tt);
  fpmul_mont(tt, t[6], tt);
  for (i = 0; i < 7; i++) fpsqr_mont(tt, tt);
  fpmul_mont(tt, t[24], tt);
  for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
  fpmul_mont(tt, t[9], tt);
  for (i = 0; i < 8; i++) fpsqr_mont(tt, tt);
  fpmul_mont(tt, t[18], tt);
  for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
  fpmul_mont(tt, t[17], tt);
  for (i = 0; i < 8; i++) fpsqr_mont(tt, tt);
  fpmul_mont(tt, r, tt);
  for (i = 0; i < 10; i++) fpsqr_mont(tt, tt);
  fpmul_mont(tt, t[16], tt);
  for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
  fpmul_mont(tt, t[7], tt);
  for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
  fpmul_mont(tt, t[0], tt);
  for (i = 0; i < 7; i++) fpsqr_mont(tt, tt);
  fpmul_mont(tt, t[12], tt);
  for (i = 0; i < 7; i++) fpsqr_mont(tt, tt);
  fpmul_mont(tt, t[19], tt);
  for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
  fpmul_mont(tt, t[22], tt);
  for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
  fpmul_mont(tt, t[25], tt);
  for (i = 0; i < 7; i++) fpsqr_mont(tt, tt);
  fpmul_mont(tt, t[2], tt);
  for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
  fpmul_mont(tt, t[10], tt);
  for (i = 0; i < 7; i++) fpsqr_mont(tt, tt);
  fpmul_mont(tt, t[22], tt);
  for (i = 0; i < 8; i++) fpsqr_mont(tt, tt);
  fpmul_mont(tt, t[18], tt);
  for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
  fpmul_mont(tt, t[4], tt);
  for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
  fpmul_mont(tt, t[14], tt);
  for (i = 0; i < 7; i++) fpsqr_mont(tt, tt);
  fpmul_mont(tt, t[13], tt);
  for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
  fpmul_mont(tt, t[5], tt);
  for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
  fpmul_mont(tt, t[23], tt);
  for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
  fpmul_mont(tt, t[21], tt);
  for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
  fpmul_mont(tt, t[2], tt);
  for (i = 0; i < 7; i++) fpsqr_mont(tt, tt);
  fpmul_mont(tt, t[23], tt);
  for (i = 0; i < 8; i++) fpsqr_mont(tt, tt);
  fpmul_mont(tt, t[12], tt);
  for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
  fpmul_mont(tt, t[9], tt);
  for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
  fpmul_mont(tt, t[3], tt);
  for (i = 0; i < 7; i++) fpsqr_mont(tt, tt);
  fpmul_mont(tt, t[13], tt);
  for (i = 0; i < 7; i++) fpsqr_mont(tt, tt);
  fpmul_mont(tt, t[17], tt);
  for (i = 0; i < 8; i++) fpsqr_mont(tt, tt);
  fpmul_mont(tt, t[26], tt);
  for (i = 0; i < 8; i++) fpsqr_mont(tt, tt);
  fpmul_mont(tt, t[5], tt);
  for (i = 0; i < 8; i++) fpsqr_mont(tt, tt);
  fpmul_mont(tt, t[8], tt);
  for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
  fpmul_mont(tt, t[2], tt);
  for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
  fpmul_mont(tt, t[11], tt);
  for (i = 0; i < 7; i++) fpsqr_mont(tt, tt);
  fpmul_mont(tt, t[20], tt);
    for (j = 0; j < 61; j++) {
        for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
        fpmul_mont(tt, t[26], tt);
    }
    fpcopy(r, tt); 
}

void fpinv_mont(felm_t r)
{
  felm_t tt;

  fpcopy(tt, r);
  fpinv_chain_mont(tt);
  fpsqr_mont(tt, tt);
  fpsqr_mont(tt, tt);
  fpmul_mont(r, r, tt);
}

// copy a GF(p^2) element r = a
void fp2copy(f2elm_t r, const f2elm_t a)
{
  fpcopy(r[0], a[0]);
  fpcopy(r[1], a[1]);
}  

// zero a GF(p^2) element r = 0
void fp2zero(f2elm_t r)
{
  fpzero(r[0]);
  fpzero(r[1]);
}

// GF(p^2) negation r = -r in GF(p^2)
void fp2neg(f2elm_t r)
{
  fpneg(r[0]);
  fpneg(r[1]);
}

// GF(p^2) addition r = a + b in GF(p^2)
void fp2add(f2elm_t r, const f2elm_t a, const f2elm_t b)
{
  fpadd(r[0], a[0], b[0]);
  fpadd(r[1], a[1], b[1]);
}

// GF(p^2) subtraction r = a - b in GF(p^2)
void fp2sub(f2elm_t r, const f2elm_t a, const f2elm_t b)
{
  fpsub(r[0], a[0], b[0]);
  fpsub(r[1], a[1], b[1]);
}

// GF(p^2) division by two r = a/2  in GF(p^2)
void fp2div2(f2elm_t r, const f2elm_t a)
{
  fpdiv2(r[0], a[0]);
  fpdiv2(r[1], a[1]);
}

void fp2correction(f2elm_t r)
{
  fpcorrection(r[0]);
  fpcorrection(r[1]);
}

// integer addition r = a + b
void mp_add(__m512i *r, const __m512i *a, const __m512i *b)
{
  __m512i a0  = a[0],  a1  = a[1],  a2  = a[2],  a3  = a[3],  a4 = a[4];
  __m512i a5  = a[5],  a6  = a[6],  a7  = a[7],  a8  = a[8],  a9 = a[9];
  __m512i a10 = a[10], a11 = a[11], a12 = a[12], a13 = a[13], a14 = a[14];
  __m512i b0  = b[0],  b1  = b[1],  b2  = b[2],  b3  = b[3],  b4 = b[4];
  __m512i b5  = b[5],  b6  = b[6],  b7  = b[7],  b8  = b[8],  b9 = b[9];
  __m512i b10 = b[10], b11 = b[11], b12 = b[12], b13 = b[13], b14 = b[14];
  __m512i r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, r12, r13, r14;
  const __m512i vbmask = VSET1(BMASK); 

  // r = a + b
  r0  = VADD(a0,  b0);  r1  = VADD(a1,  b1);  r2  = VADD(a2,  b2); 
  r3  = VADD(a3,  b3);  r4  = VADD(a4,  b4);  r5  = VADD(a5,  b5);
  r6  = VADD(a6,  b6);  r7  = VADD(a7,  b7);  r8  = VADD(a8,  b8);
  r9  = VADD(a9,  b9);  r10 = VADD(a10, b10); r11 = VADD(a11, b11);
  r12 = VADD(a12, b12); r13 = VADD(a13, b13); r14 = VADD(a14, b14);

  // carry propagation 
  r1  = VADD(r1,  VSRA(r0,  BRADIX)); r0  = VAND(r0,  vbmask);
  r2  = VADD(r2,  VSRA(r1,  BRADIX)); r1  = VAND(r1,  vbmask);
  r3  = VADD(r3,  VSRA(r2,  BRADIX)); r2  = VAND(r2,  vbmask);
  r4  = VADD(r4,  VSRA(r3,  BRADIX)); r3  = VAND(r3,  vbmask);
  r5  = VADD(r5,  VSRA(r4,  BRADIX)); r4  = VAND(r4,  vbmask);
  r6  = VADD(r6,  VSRA(r5,  BRADIX)); r5  = VAND(r5,  vbmask);
  r7  = VADD(r7,  VSRA(r6,  BRADIX)); r6  = VAND(r6,  vbmask);
  r8  = VADD(r8,  VSRA(r7,  BRADIX)); r7  = VAND(r7,  vbmask);
  r9  = VADD(r9,  VSRA(r8,  BRADIX)); r8  = VAND(r8,  vbmask);
  r10 = VADD(r10, VSRA(r9,  BRADIX)); r9  = VAND(r9,  vbmask);
  r11 = VADD(r11, VSRA(r10, BRADIX)); r10 = VAND(r10, vbmask);
  r12 = VADD(r12, VSRA(r11, BRADIX)); r11 = VAND(r11, vbmask);
  r13 = VADD(r13, VSRA(r12, BRADIX)); r12 = VAND(r12, vbmask);
  r14 = VADD(r14, VSRA(r13, BRADIX)); r13 = VAND(r13, vbmask);

  r[0]  = r0;  r[1]  = r1;  r[2]  = r2;  r[3]  = r3;  r[4]  = r4; 
  r[5]  = r5;  r[6]  = r6;  r[7]  = r7;  r[8]  = r8;  r[9]  = r9;
  r[10] = r10; r[11] = r11; r[12] = r12; r[13] = r13; r[14] = r14;
}

static void mp_addfast(__m512i *r, const __m512i *a, const __m512i *b)
{
  mp_add(r, a, b);
}

// GF(p^2) addition without correction r = a + b in GF(p^2) 
void mp2_add(f2elm_t r, const f2elm_t a, const f2elm_t b)
{
  mp_addfast(r[0], a[0], b[0]);
  mp_addfast(r[1], a[1], b[1]);
}

// GF(p^2) subtraction with correction with 2p, r = a-b+2p in GF(p^2).
void mp2_sub_p2(f2elm_t r, const f2elm_t a, const f2elm_t b)
{
  mp_sub_p2(r[0], a[0], b[0]);
  mp_sub_p2(r[1], a[1], b[1]);
} 

// GF(p^2) subtraction with correction with 4p, r = a-b+4p in GF(p^2).
static void mp2_sub_p4(f2elm_t r, const f2elm_t a, const f2elm_t b)
{
  mp_sub_p4(r[0], a[0], b[0]);
  mp_sub_p4(r[1], a[1], b[1]);
} 

// GF(p^2) squaring using Montgomery arithmetic r = a^2 in GF(p^2)
void fp2sqr_mont(f2elm_t r, const f2elm_t a)
{
  felm_t t1, t2, t3;

  mp_addfast(t1, a[0], a[1]);           // t1 = a0+a1           in [0, 4p]
  mp_sub_p4(t2, a[0], a[1]);            // t2 = a0-a1           in [0, 4p]
  mp_addfast(t3, a[0], a[0]);           // t3 = 2*a0            in [0, 4p]
  fpmul_mont(r[0], t1, t2);             // r0 = (a0+a1)(a0-a1)  in [0, 2p]
  fpmul_mont(r[1], t3, a[1]);           // r1 = 2a0a1           in [0, 2p]
}

// integer subtraction r = r - a - b where len(r) = len(a) = len(b) = 2*NWORDS
static void mp_dblsubfast(__m512i *r, const __m512i *a, const __m512i *b)
{
  __m512i a0  = a[0],  a1  = a[1],  a2  = a[2],  a3  = a[3],  a4  = a[4];
  __m512i a5  = a[5],  a6  = a[6],  a7  = a[7],  a8  = a[8],  a9  = a[9];
  __m512i a10 = a[10], a11 = a[11], a12 = a[12], a13 = a[13], a14 = a[14];
  __m512i a15 = a[15], a16 = a[16], a17 = a[17], a18 = a[18], a19 = a[19];
  __m512i a20 = a[20], a21 = a[21], a22 = a[22], a23 = a[23], a24 = a[24];
  __m512i a25 = a[25], a26 = a[26], a27 = a[27], a28 = a[28], a29 = a[29];
  __m512i b0  = b[0],  b1  = b[1],  b2  = b[2],  b3  = b[3],  b4  = b[4];
  __m512i b5  = b[5],  b6  = b[6],  b7  = b[7],  b8  = b[8],  b9  = b[9];
  __m512i b10 = b[10], b11 = b[11], b12 = b[12], b13 = b[13], b14 = b[14];
  __m512i b15 = b[15], b16 = b[16], b17 = b[17], b18 = b[18], b19 = b[19];
  __m512i b20 = b[20], b21 = b[21], b22 = b[22], b23 = b[23], b24 = b[24];
  __m512i b25 = b[25], b26 = b[26], b27 = b[27], b28 = b[28], b29 = b[29];
  __m512i r0  = r[0],  r1  = r[1],  r2  = r[2],  r3  = r[3],  r4  = r[4];
  __m512i r5  = r[5],  r6  = r[6],  r7  = r[7],  r8  = r[8],  r9  = r[9];
  __m512i r10 = r[10], r11 = r[11], r12 = r[12], r13 = r[13], r14 = r[14];
  __m512i r15 = r[15], r16 = r[16], r17 = r[17], r18 = r[18], r19 = r[19];
  __m512i r20 = r[20], r21 = r[21], r22 = r[22], r23 = r[23], r24 = r[24];
  __m512i r25 = r[25], r26 = r[26], r27 = r[27], r28 = r[28], r29 = r[29];
  const __m512i vbmask = VSET1(BMASK);

  // r = r - a - b
  r0  = VSUB(r0,  VADD(a0,  b0));  r1  = VSUB(r1,  VADD(a1,  b1));
  r2  = VSUB(r2,  VADD(a2,  b2));  r3  = VSUB(r3,  VADD(a3,  b3));
  r4  = VSUB(r4,  VADD(a4,  b4));  r5  = VSUB(r5,  VADD(a5,  b5));
  r6  = VSUB(r6,  VADD(a6,  b6));  r7  = VSUB(r7,  VADD(a7,  b7));
  r8  = VSUB(r8,  VADD(a8,  b8));  r9  = VSUB(r9,  VADD(a9,  b9));
  r10 = VSUB(r10, VADD(a10, b10)); r11 = VSUB(r11, VADD(a11, b11));
  r12 = VSUB(r12, VADD(a12, b12)); r13 = VSUB(r13, VADD(a13, b13));
  r14 = VSUB(r14, VADD(a14, b14)); r15 = VSUB(r15, VADD(a15, b15));
  r16 = VSUB(r16, VADD(a16, b16)); r17 = VSUB(r17, VADD(a17, b17));
  r18 = VSUB(r18, VADD(a18, b18)); r19 = VSUB(r19, VADD(a19, b19));
  r20 = VSUB(r20, VADD(a20, b20)); r21 = VSUB(r21, VADD(a21, b21));
  r22 = VSUB(r22, VADD(a22, b22)); r23 = VSUB(r23, VADD(a23, b23));
  r24 = VSUB(r24, VADD(a24, b24)); r25 = VSUB(r25, VADD(a25, b25));
  r26 = VSUB(r26, VADD(a26, b26)); r27 = VSUB(r27, VADD(a27, b27));
  r28 = VSUB(r28, VADD(a28, b28)); r29 = VSUB(r29, VADD(a29, b29));

  // carry propagation
  // will be done in reduction
  

  r[0]  = r0;  r[1]  = r1;  r[2]  = r2;  r[3]  = r3;  r[4]  = r4;  
  r[5]  = r5;  r[6]  = r6;  r[7]  = r7;  r[8]  = r8;  r[9]  = r9;  
  r[10] = r10; r[11] = r11; r[12] = r12; r[13] = r13; r[14] = r14; 
  r[15] = r15; r[16] = r16; r[17] = r17; r[18] = r18; r[19] = r19;
  r[20] = r20; r[21] = r21; r[22] = r22; r[23] = r23; r[24] = r24; 
  r[25] = r25; r[26] = r26; r[27] = r27; r[28] = r28; r[29] = r29;
}

// integer subtraction followed by addition with p*765
//  r = a - b + (p*2^765) if a-b<0; otherwise r = a - b
static void mp_subaddfast(__m512i *r, const __m512i *a, const __m512i *b)
{
  __m512i a0  = a[0],  a1  = a[1],  a2  = a[2],  a3  = a[3],  a4  = a[4];
  __m512i a5  = a[5],  a6  = a[6],  a7  = a[7],  a8  = a[8],  a9  = a[9];
  __m512i a10 = a[10], a11 = a[11], a12 = a[12], a13 = a[13], a14 = a[14];
  __m512i a15 = a[15], a16 = a[16], a17 = a[17], a18 = a[18], a19 = a[19];
  __m512i a20 = a[20], a21 = a[21], a22 = a[22], a23 = a[23], a24 = a[24];
  __m512i a25 = a[25], a26 = a[26], a27 = a[27], a28 = a[28], a29 = a[29];
  __m512i b0  = b[0],  b1  = b[1],  b2  = b[2],  b3  = b[3],  b4  = b[4];
  __m512i b5  = b[5],  b6  = b[6],  b7  = b[7],  b8  = b[8],  b9  = b[9];
  __m512i b10 = b[10], b11 = b[11], b12 = b[12], b13 = b[13], b14 = b[14];
  __m512i b15 = b[15], b16 = b[16], b17 = b[17], b18 = b[18], b19 = b[19];
  __m512i b20 = b[20], b21 = b[21], b22 = b[22], b23 = b[23], b24 = b[24];
  __m512i b25 = b[25], b26 = b[26], b27 = b[27], b28 = b[28], b29 = b[29];
  __m512i r0,  r1,  r2,  r3,  r4,  r5,  r6,  r7,  r8,  r9,  r10, r11, r12, r13, r14;
  __m512i r15, r16, r17, r18, r19, r20, r21, r22, r23, r24, r25, r26, r27, r28, r29, smask;
  const __m512i vp0  = VSET1(p751[0]),  vp1  = VSET1(p751[1]),  vp2  = VSET1(p751[2]);
  const __m512i vp3  = VSET1(p751[3]),  vp4  = VSET1(p751[4]),  vp5  = VSET1(p751[5]);
  const __m512i vp6  = VSET1(p751[6]),  vp7  = VSET1(p751[7]),  vp8  = VSET1(p751[8]);
  const __m512i vp9  = VSET1(p751[9]),  vp10 = VSET1(p751[10]), vp11 = VSET1(p751[11]);
  const __m512i vp12 = VSET1(p751[12]), vp13 = VSET1(p751[13]), vp14 = VSET1(p751[14]);
  const __m512i vbmask = VSET1(BMASK);

  // r = a - b
  r0  = VSUB(a0,  b0);  r1  = VSUB(a1,  b1);  r2  = VSUB(a2,  b2);
  r3  = VSUB(a3,  b3);  r4  = VSUB(a4,  b4);  r5  = VSUB(a5,  b5);
  r6  = VSUB(a6,  b6);  r7  = VSUB(a7,  b7);  r8  = VSUB(a8,  b8);
  r9  = VSUB(a9,  b9);  r10 = VSUB(a10, b10); r11 = VSUB(a11, b11);
  r12 = VSUB(a12, b12); r13 = VSUB(a13, b13); r14 = VSUB(a14, b14);
  r15 = VSUB(a15, b15); r16 = VSUB(a16, b16); r17 = VSUB(a17, b17);
  r18 = VSUB(a18, b18); r19 = VSUB(a19, b19); r20 = VSUB(a20, b20);
  r21 = VSUB(a21, b21); r22 = VSUB(a22, b22); r23 = VSUB(a23, b23);
  r24 = VSUB(a24, b24); r25 = VSUB(a25, b25); r26 = VSUB(a26, b26);
  r27 = VSUB(a27, b27); r28 = VSUB(a28, b28); r29 = VSUB(a29, b29);
  
  // carry propagation
  r1  = VADD(r1,  VSRA(r0,  BRADIX)); r0  = VAND(r0,  vbmask);
  r2  = VADD(r2,  VSRA(r1,  BRADIX)); r1  = VAND(r1,  vbmask);
  r3  = VADD(r3,  VSRA(r2,  BRADIX)); r2  = VAND(r2,  vbmask);
  r4  = VADD(r4,  VSRA(r3,  BRADIX)); r3  = VAND(r3,  vbmask);
  r5  = VADD(r5,  VSRA(r4,  BRADIX)); r4  = VAND(r4,  vbmask);
  r6  = VADD(r6,  VSRA(r5,  BRADIX)); r5  = VAND(r5,  vbmask);
  r7  = VADD(r7,  VSRA(r6,  BRADIX)); r6  = VAND(r6,  vbmask);
  r8  = VADD(r8,  VSRA(r7,  BRADIX)); r7  = VAND(r7,  vbmask);
  r9  = VADD(r9,  VSRA(r8,  BRADIX)); r8  = VAND(r8,  vbmask);
  r10 = VADD(r10, VSRA(r9,  BRADIX)); r9  = VAND(r9,  vbmask);
  r11 = VADD(r11, VSRA(r10, BRADIX)); r10 = VAND(r10, vbmask);
  r12 = VADD(r12, VSRA(r11, BRADIX)); r11 = VAND(r11, vbmask);
  r13 = VADD(r13, VSRA(r12, BRADIX)); r12 = VAND(r12, vbmask);
  r14 = VADD(r14, VSRA(r13, BRADIX)); r13 = VAND(r13, vbmask);
  r15 = VADD(r15, VSRA(r14, BRADIX)); r14 = VAND(r14, vbmask);
  r16 = VADD(r16, VSRA(r15, BRADIX)); r15 = VAND(r15, vbmask);
  r17 = VADD(r17, VSRA(r16, BRADIX)); r16 = VAND(r16, vbmask);
  r18 = VADD(r18, VSRA(r17, BRADIX)); r17 = VAND(r17, vbmask);
  r19 = VADD(r19, VSRA(r18, BRADIX)); r18 = VAND(r18, vbmask);
  r20 = VADD(r20, VSRA(r19, BRADIX)); r19 = VAND(r19, vbmask);
  r21 = VADD(r21, VSRA(r20, BRADIX)); r20 = VAND(r20, vbmask);
  r22 = VADD(r22, VSRA(r21, BRADIX)); r21 = VAND(r21, vbmask);
  r23 = VADD(r23, VSRA(r22, BRADIX)); r22 = VAND(r22, vbmask);
  r24 = VADD(r24, VSRA(r23, BRADIX)); r23 = VAND(r23, vbmask);
  r25 = VADD(r25, VSRA(r24, BRADIX)); r24 = VAND(r24, vbmask);
  r26 = VADD(r26, VSRA(r25, BRADIX)); r25 = VAND(r25, vbmask);
  r27 = VADD(r27, VSRA(r26, BRADIX)); r26 = VAND(r26, vbmask);
  r28 = VADD(r28, VSRA(r27, BRADIX)); r27 = VAND(r27, vbmask);
  r29 = VADD(r29, VSRA(r28, BRADIX)); r28 = VAND(r28, vbmask);

  // if r is non-negative, smask = all-0 
  // if r is     negative, smask = all-1
  smask = VSRA(r29, 63);
  // r = r + ((p*2^765) & smask)
  r15 = VADD(r15, VAND(vp0,  smask)); r16 = VADD(r16, VAND(vp1,  smask)); 
  r17 = VADD(r17, VAND(vp2,  smask)); r18 = VADD(r18, VAND(vp3,  smask)); 
  r19 = VADD(r19, VAND(vp4,  smask)); r20 = VADD(r20, VAND(vp5,  smask)); 
  r21 = VADD(r21, VAND(vp6,  smask)); r22 = VADD(r22, VAND(vp7,  smask)); 
  r23 = VADD(r23, VAND(vp8,  smask)); r24 = VADD(r24, VAND(vp9,  smask)); 
  r25 = VADD(r25, VAND(vp10, smask)); r26 = VADD(r26, VAND(vp11, smask)); 
  r27 = VADD(r27, VAND(vp12, smask)); r28 = VADD(r28, VAND(vp13, smask)); 
  r29 = VADD(r29, VAND(vp14, smask)); 

  // carry propagation
  // will be done in reduction 

  r[0]  = r0;  r[1]  = r1;  r[2]  = r2;  r[3]  = r3;  r[4]  = r4;  
  r[5]  = r5;  r[6]  = r6;  r[7]  = r7;  r[8]  = r8;  r[9]  = r9;  
  r[10] = r10; r[11] = r11; r[12] = r12; r[13] = r13; r[14] = r14; 
  r[15] = r15; r[16] = r16; r[17] = r17; r[18] = r18; r[19] = r19;
  r[20] = r20; r[21] = r21; r[22] = r22; r[23] = r23; r[24] = r24; 
  r[25] = r25; r[26] = r26; r[27] = r27; r[28] = r28; r[29] = r29;
}

// GF(p^2) multiplication using Montgomery arithmetic, r = a*b in GF(p^2).
void fp2mul_mont(f2elm_t r, const f2elm_t a, const f2elm_t b)
{
  felm_t t1, t2;
  dfelm_t tt1, tt2, tt3;

  mp_addfast(t1, a[0], a[1]);           // t1 = a0+a1           in [0, 4p]
  mp_addfast(t2, b[0], b[1]);           // t2 = b0+b1           in [0, 4p]
  mp_mul(tt1, a[0], b[0]);              // tt1 = a0*b0          in [0, 4p^2]
  mp_mul(tt2, a[1], b[1]);              // tt2 = a1*b1          in [0, 4p^2]
  mp_mul(tt3, t1, t2);                  // tt3 = (a0+a1)(b0+b1) in [0, 16p^2]
  mp_dblsubfast(tt3, tt1, tt2);         // tt3 = a0b1+a1b0      in [0, 8p^2]    (only be used here)
  mp_subaddfast(tt1, tt1, tt2);         // tt1 = a0b0-a1b1      in [0, p*2^765] (only be used here)
  rdc_mont(r[1], tt3);                  // r1 = a0b1+a1b0 mod 2p                (only be used here)
  rdc_mont(r[0], tt1);                  // r0 = a0b0-a1b1 mod 2p                (only be used here)
 }

// GF(p^2) inversion using Montgomery arithmetic
void fp2inv_mont(f2elm_t r)
{
  f2elm_t t1;

  fpsqr_mont(t1[0], r[0]);  
  fpsqr_mont(t1[1], r[1]);
  fpadd(t1[0], t1[0], t1[1]);
  fpinv_mont(t1[0]);
  fpneg(r[1]);
  fpmul_mont(r[0], r[0], t1[0]);
  fpmul_mont(r[1], r[1], t1[0]);        
}

void to_fp2mont(f2elm_t r, const f2elm_t a)
{
  to_mont(r[0], a[0]);
  to_mont(r[1], a[1]);
}

void from_fp2mont(f2elm_t r, const f2elm_t a)
{
  from_mont(r[0], a[0]);
  from_mont(r[1], a[1]);
}
