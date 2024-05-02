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
void fpmul_mont_v2(felm_t r, const felm_t a, const felm_t b)
{
  __m512i a0 = a[0], a1 = a[1], a2 = a[2];
  __m512i a3 = a[3], a4 = a[4], a5 = a[5];
  __m512i a6 = a[6], a7 = a[7], a8 = a[8];
  __m512i b0 = b[0], b1 = b[1], b2 = b[2];
  __m512i b3 = b[3], b4 = b[4], b5 = b[5];
  __m512i b6 = b[6], b7 = b[7], b8 = b[8];
  __m512i z0  = VZERO, z1  = VZERO, z2  = VZERO;
  __m512i z3  = VZERO, z4  = VZERO, z5  = VZERO;
  __m512i z6  = VZERO, z7  = VZERO, z8  = VZERO;
  __m512i z9  = VZERO, z10 = VZERO, z11 = VZERO;
  __m512i z12 = VZERO, z13 = VZERO, z14 = VZERO;
  __m512i z15 = VZERO, z16 = VZERO, z17 = VZERO;
  __m512i h0  = VZERO, h1  = VZERO, h2  = VZERO;
  __m512i h3  = VZERO, h4  = VZERO, h5  = VZERO;
  __m512i h6  = VZERO, h7  = VZERO, h8  = VZERO;
  __m512i h9  = VZERO, h10 = VZERO, h11 = VZERO;
  __m512i h12 = VZERO, h13 = VZERO, h14 = VZERO;
  __m512i h15 = VZERO, h16 = VZERO, h17 = VZERO;
  __m512i r0, r1, r2, r3, r4, r5, r6, r7, r8, u;
  const __m512i vp4 = VSET1(p434p1[4]), vp5 = VSET1(p434p1[5]);
  const __m512i vp6 = VSET1(p434p1[6]), vp7 = VSET1(p434p1[7]);
  const __m512i vp8 = VSET1(p434p1[8]);
  const __m512i vbmask = VSET1(BMASK);

  // ---------------------------------------------------------------------------
  // integer multiplication 

  z0 = VMACLO(z0, a0, b0);
  h1 = VMACHI(h1, a0, b0);

  z1 = VMACLO(z1, a0, b1); z1 = VMACLO(z1, a1, b0);
  h2 = VMACHI(h2, a0, b1); h2 = VMACHI(h2, a1, b0);

  z2 = VMACLO(z2, a0, b2); z2 = VMACLO(z2, a1, b1); z2 = VMACLO(z2, a2, b0);
  h3 = VMACHI(h3, a0, b2); h3 = VMACHI(h3, a1, b1); h3 = VMACHI(h3, a2, b0);


  z3 = VMACLO(z3, a0, b3); z3 = VMACLO(z3, a1, b2); z3 = VMACLO(z3, a2, b1); 
  z3 = VMACLO(z3, a3, b0);
  h4 = VMACHI(h4, a0, b3); h4 = VMACHI(h4, a1, b2); h4 = VMACHI(h4, a2, b1); 
  h4 = VMACHI(h4, a3, b0);

  z4 = VMACLO(z4, a0, b4); z4 = VMACLO(z4, a1, b3); z4 = VMACLO(z4, a2, b2); 
  z4 = VMACLO(z4, a3, b1); z4 = VMACLO(z4, a4, b0);
  h5 = VMACHI(h5, a0, b4); h5 = VMACHI(h5, a1, b3); h5 = VMACHI(h5, a2, b2); 
  h5 = VMACHI(h5, a3, b1); h5 = VMACHI(h5, a4, b0);

  z5 = VMACLO(z5, a0, b5); z5 = VMACLO(z5, a1, b4); z5 = VMACLO(z5, a2, b3); 
  z5 = VMACLO(z5, a3, b2); z5 = VMACLO(z5, a4, b1); z5 = VMACLO(z5, a5, b0);
  h6 = VMACHI(h6, a0, b5); h6 = VMACHI(h6, a1, b4); h6 = VMACHI(h6, a2, b3); 
  h6 = VMACHI(h6, a3, b2); h6 = VMACHI(h6, a4, b1); h6 = VMACHI(h6, a5, b0);

  z6 = VMACLO(z6, a0, b6); z6 = VMACLO(z6, a1, b5); z6 = VMACLO(z6, a2, b4); 
  z6 = VMACLO(z6, a3, b3); z6 = VMACLO(z6, a4, b2); z6 = VMACLO(z6, a5, b1); 
  z6 = VMACLO(z6, a6, b0);
  h7 = VMACHI(h7, a0, b6); h7 = VMACHI(h7, a1, b5); h7 = VMACHI(h7, a2, b4); 
  h7 = VMACHI(h7, a3, b3); h7 = VMACHI(h7, a4, b2); h7 = VMACHI(h7, a5, b1); 
  h7 = VMACHI(h7, a6, b0);

  z7 = VMACLO(z7, a0, b7); z7 = VMACLO(z7, a1, b6); z7 = VMACLO(z7, a2, b5); 
  z7 = VMACLO(z7, a3, b4); z7 = VMACLO(z7, a4, b3); z7 = VMACLO(z7, a5, b2); 
  z7 = VMACLO(z7, a6, b1); z7 = VMACLO(z7, a7, b0);
  h8 = VMACHI(h8, a0, b7); h8 = VMACHI(h8, a1, b6); h8 = VMACHI(h8, a2, b5); 
  h8 = VMACHI(h8, a3, b4); h8 = VMACHI(h8, a4, b3); h8 = VMACHI(h8, a5, b2); 
  h8 = VMACHI(h8, a6, b1); h8 = VMACHI(h8, a7, b0);

  z8 = VMACLO(z8, a0, b8); z8 = VMACLO(z8, a1, b7); z8 = VMACLO(z8, a2, b6); 
  z8 = VMACLO(z8, a3, b5); z8 = VMACLO(z8, a4, b4); z8 = VMACLO(z8, a5, b3); 
  z8 = VMACLO(z8, a6, b2); z8 = VMACLO(z8, a7, b1); z8 = VMACLO(z8, a8, b0);
  h9 = VMACHI(h9, a0, b8); h9 = VMACHI(h9, a1, b7); h9 = VMACHI(h9, a2, b6); 
  h9 = VMACHI(h9, a3, b5); h9 = VMACHI(h9, a4, b4); h9 = VMACHI(h9, a5, b3); 
  h9 = VMACHI(h9, a6, b2); h9 = VMACHI(h9, a7, b1); h9 = VMACHI(h9, a8, b0);

  // some carry propagations

  z1 = VADD(z1, VADD(h1, h1)); 
  z2 = VADD(z2, VADD(h2, h2));
  z3 = VADD(z3, VADD(h3, h3)); 
  z4 = VADD(z4, VADD(h4, h4));

  z1 = VADD(z1, VSRA(z0, BRADIX)); z0 = VAND(z0, vbmask);
  z2 = VADD(z2, VSRA(z1, BRADIX)); z1 = VAND(z1, vbmask);
  z3 = VADD(z3, VSRA(z2, BRADIX)); z2 = VAND(z2, vbmask);
  z4 = VADD(z4, VSRA(z3, BRADIX)); z3 = VAND(z3, vbmask);

  // 2nd loop of integer multiplication

  u = z0;
  z4 = VMACLO(z4, u, vp4); z5 = VMACLO(z5, u, vp5); z6 = VMACLO(z6, u, vp6); 
  z7 = VMACLO(z7, u, vp7); z8 = VMACLO(z8, u, vp8);
  h5 = VMACHI(h5, u, vp4); h6 = VMACHI(h6, u, vp5); h7 = VMACHI(h7, u, vp6);
  h8 = VMACHI(h8, u, vp7); h9 = VMACHI(h9, u, vp8);
  z5 = VADD(z5, VSRA(z4, BRADIX)); z4 = VAND(z4, vbmask); 
  z5 = VADD(z5, VADD(h5, h5));

  z9 = VMACLO(z9, a1, b8); z9 = VMACLO(z9, a2, b7); z9 = VMACLO(z9, a3, b6); 
  z9 = VMACLO(z9, a4, b5); z9 = VMACLO(z9, a5, b4); z9 = VMACLO(z9, a6, b3); 
  z9 = VMACLO(z9, a7, b2); z9 = VMACLO(z9, a8, b1); 
  h10 = VMACHI(h10, a1, b8); h10 = VMACHI(h10, a2, b7); h10 = VMACHI(h10, a3, b6); 
  h10 = VMACHI(h10, a4, b5); h10 = VMACHI(h10, a5, b4); h10 = VMACHI(h10, a6, b3); 
  h10 = VMACHI(h10, a7, b2); h10 = VMACHI(h10, a8, b1); 

  u = z1; 
  z5 = VMACLO(z5, u, vp4); z6  = VMACLO(z6,  u, vp5); z7 = VMACLO(z7, u, vp6); 
  z8 = VMACLO(z8, u, vp7); z9  = VMACLO(z9,  u, vp8);
  h6 = VMACHI(h6, u, vp4); h7  = VMACHI(h7,  u, vp5); h8 = VMACHI(h8, u, vp6);
  h9 = VMACHI(h9, u, vp7); h10 = VMACHI(h10, u, vp8);
  z6 = VADD(z6, VSRA(z5, BRADIX)); z5 = VAND(z5, vbmask); 
  z6 = VADD(z6, VADD(h6, h6));

  z10 = VMACLO(z10, a2, b8); z10 = VMACLO(z10, a3, b7); z10 = VMACLO(z10, a4, b6);
  z10 = VMACLO(z10, a5, b5); z10 = VMACLO(z10, a6, b4); z10 = VMACLO(z10, a7, b3); 
  z10 = VMACLO(z10, a8, b2); 
  h11 = VMACHI(h11, a2, b8); h11 = VMACHI(h11, a3, b7); h11 = VMACHI(h11, a4, b6); 
  h11 = VMACHI(h11, a5, b5); h11 = VMACHI(h11, a6, b4); h11 = VMACHI(h11, a7, b3); 
  h11 = VMACHI(h11, a8, b2); 

  u = z2; 
  z6  = VMACLO(z6,  u, vp4); z7  = VMACLO(z7,  u, vp5); z8 = VMACLO(z8, u, vp6); 
  z9  = VMACLO(z9,  u, vp7); z10 = VMACLO(z10, u, vp8);
  h7  = VMACHI(h7,  u, vp4); h8  = VMACHI(h8,  u, vp5); h9 = VMACHI(h9, u, vp6);
  h10 = VMACHI(h10, u, vp7); h11 = VMACHI(h11, u, vp8);
  z7 = VADD(z7, VSRA(z6, BRADIX)); z6 = VAND(z6, vbmask); 
  z7 = VADD(z7, VADD(h7, h7));

  z11 = VMACLO(z11, a3, b8); z11 = VMACLO(z11, a4, b7); z11 = VMACLO(z11, a5, b6); 
  z11 = VMACLO(z11, a6, b5); z11 = VMACLO(z11, a7, b4); z11 = VMACLO(z11, a8, b3); 
  h12 = VMACHI(h12, a3, b8); h12 = VMACHI(h12, a4, b7); h12 = VMACHI(h12, a5, b6); 
  h12 = VMACHI(h12, a6, b5); h12 = VMACHI(h12, a7, b4); h12 = VMACHI(h12, a8, b3); 

  u = z3; 
  z7  = VMACLO(z7,  u, vp4); z8  = VMACLO(z8,  u, vp5); z9  = VMACLO(z9,  u, vp6); 
  z10 = VMACLO(z10, u, vp7); z11 = VMACLO(z11, u, vp8);
  h8  = VMACHI(h8,  u, vp4); h9  = VMACHI(h9,  u, vp5); h10 = VMACHI(h10, u, vp6);
  h11 = VMACHI(h11, u, vp7); h12 = VMACHI(h12, u, vp8);
  z8 = VADD(z8, VSRA(z7, BRADIX)); z7 = VAND(z7, vbmask); 
  z8 = VADD(z8, VADD(h8, h8));

  z12 = VMACLO(z12, a4, b8); z12 = VMACLO(z12, a5, b7); z12 = VMACLO(z12, a6, b6); 
  z12 = VMACLO(z12, a7, b5); z12 = VMACLO(z12, a8, b4); 
  h13 = VMACHI(h13, a4, b8); h13 = VMACHI(h13, a5, b7); h13 = VMACHI(h13, a6, b6); 
  h13 = VMACHI(h13, a7, b5); h13 = VMACHI(h13, a8, b4); 

  u = z4; 
  z8  = VMACLO(z8,  u, vp4); z9  = VMACLO(z9,  u, vp5); z10 = VMACLO(z10, u, vp6); 
  z11 = VMACLO(z11, u, vp7); z12 = VMACLO(z12, u, vp8);
  h9  = VMACHI(h9,  u, vp4); h10 = VMACHI(h10, u, vp5); h11 = VMACHI(h11, u, vp6);
  h12 = VMACHI(h12, u, vp7); h13 = VMACHI(h13, u, vp8);
  z9 = VADD(z9, VSRA(z8, BRADIX)); z8 = VAND(z8, vbmask); 
  z9 = VADD(z9, VADD(h9, h9));  

  z13 = VMACLO(z13, a5, b8); z13 = VMACLO(z13, a6, b7); z13 = VMACLO(z13, a7, b6); 
  z13 = VMACLO(z13, a8, b5); 
  h14 = VMACHI(h14, a5, b8); h14 = VMACHI(h14, a6, b7); h14 = VMACHI(h14, a7, b6); 
  h14 = VMACHI(h14, a8, b5); 

  u = z5; 
  z9  = VMACLO(z9,  u, vp4); z10 = VMACLO(z10, u, vp5); z11 = VMACLO(z11, u, vp6); 
  z12 = VMACLO(z12, u, vp7); z13 = VMACLO(z13, u, vp8);
  h10 = VMACHI(h10, u, vp4); h11 = VMACHI(h11, u, vp5); h12 = VMACHI(h12, u, vp6);
  h13 = VMACHI(h13, u, vp7); h14 = VMACHI(h14, u, vp8);
  z10 = VADD(z10, VSRA(z9, BRADIX)); z9 = VAND(z9, vbmask); 
  z10 = VADD(z10, VADD(h10, h10)); 

  z14 = VMACLO(z14, a6, b8); z14 = VMACLO(z14, a7, b7); z14 = VMACLO(z14, a8, b6); 
  h15 = VMACHI(h15, a6, b8); h15 = VMACHI(h15, a7, b7); h15 = VMACHI(h15, a8, b6); 

  u = z6; 
  z10 = VMACLO(z10, u, vp4); z11 = VMACLO(z11, u, vp5); z12 = VMACLO(z12, u, vp6); 
  z13 = VMACLO(z13, u, vp7); z14 = VMACLO(z14, u, vp8);
  h11 = VMACHI(h11, u, vp4); h12 = VMACHI(h12, u, vp5); h13 = VMACHI(h13, u, vp6);
  h14 = VMACHI(h14, u, vp7); h15 = VMACHI(h15, u, vp8);
  z11 = VADD(z11, VSRA(z10, BRADIX)); z10 = VAND(z10, vbmask); 
  z11 = VADD(z11, VADD(h11, h11));

  z15 = VMACLO(z15, a7, b8); z15 = VMACLO(z15, a8, b7); 
  h16 = VMACHI(h16, a7, b8); h16 = VMACHI(h16, a8, b7); 

  u = z7; 
  z11 = VMACLO(z11, u, vp4); z12 = VMACLO(z12, u, vp5); z13 = VMACLO(z13, u, vp6); 
  z14 = VMACLO(z14, u, vp7); z15 = VMACLO(z15, u, vp8);
  h12 = VMACHI(h12, u, vp4); h13 = VMACHI(h13, u, vp5); h14 = VMACHI(h14, u, vp6);
  h15 = VMACHI(h15, u, vp7); h16 = VMACHI(h16, u, vp8);
  z12 = VADD(z12, VSRA(z11, BRADIX)); z11 = VAND(z11, vbmask); 
  z12 = VADD(z12, VADD(h12, h12));

  z16 = VMACLO(z16, a8, b8);
  h17 = VMACHI(h17, a8, b8);

  u = z8; 
  z12 = VMACLO(z12, u, vp4); z13 = VMACLO(z13, u, vp5); z14 = VMACLO(z14, u, vp6); 
  z15 = VMACLO(z15, u, vp7); z16 = VMACLO(z16, u, vp8);
  h13 = VMACHI(h13, u, vp4); h14 = VMACHI(h14, u, vp5); h15 = VMACHI(h15, u, vp6);
  h16 = VMACHI(h16, u, vp7); h17 = VMACHI(h17, u, vp8);
  z13 = VADD(z13, VSRA(z12, BRADIX)); z12 = VAND(z12, vbmask);
  z13 = VADD(z13, VADD(h13, h13)); 

  z14 = VADD(z14, VADD(h14, h14));
  z15 = VADD(z15, VADD(h15, h15));
  z16 = VADD(z16, VADD(h16, h16));
  z17 = VADD(z17, VADD(h17, h17));

  z14 = VADD(z14, VSRA(z13, BRADIX)); z13 = VAND(z13, vbmask);
  z15 = VADD(z15, VSRA(z14, BRADIX)); z14 = VAND(z14, vbmask);
  z16 = VADD(z16, VSRA(z15, BRADIX)); z15 = VAND(z15, vbmask);
  z17 = VADD(z17, VSRA(z16, BRADIX)); z16 = VAND(z16, vbmask);

  r[0] =  z9; r[1] = z10; r[2] = z11; 
  r[3] = z12; r[4] = z13; r[5] = z14; 
  r[6] = z15; r[7] = z16; r[8] = z17; 
}

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
  __m512i a0 = a[0], a1 = a[1], a2 = a[2];
  __m512i a3 = a[3], a4 = a[4], a5 = a[5];
  __m512i a6 = a[6], a7 = a[7], a8 = a[8];
  __m512i d0, d1, d2, d3, d4, d5, d6, d7, d8;
  __m512i z0  = VZERO, z1  = VZERO, z2  = VZERO;
  __m512i z3  = VZERO, z4  = VZERO, z5  = VZERO;
  __m512i z6  = VZERO, z7  = VZERO, z8  = VZERO;
  __m512i z9  = VZERO, z10 = VZERO, z11 = VZERO;
  __m512i z12 = VZERO, z13 = VZERO, z14 = VZERO;
  __m512i z15 = VZERO, z16 = VZERO, z17 = VZERO;
  __m512i h0  = VZERO, h1  = VZERO, h2  = VZERO;
  __m512i h3  = VZERO, h4  = VZERO, h5  = VZERO;
  __m512i h6  = VZERO, h7  = VZERO, h8  = VZERO;
  __m512i h9  = VZERO, h10 = VZERO, h11 = VZERO;
  __m512i h12 = VZERO, h13 = VZERO, h14 = VZERO;
  __m512i h15 = VZERO, h16 = VZERO, h17 = VZERO;
  __m512i r0, r1, r2, r3, r4, r5, r6, r7, r8, u;
  const __m512i vp4 = VSET1(p434p1[4]), vp5 = VSET1(p434p1[5]);
  const __m512i vp6 = VSET1(p434p1[6]), vp7 = VSET1(p434p1[7]);
  const __m512i vp8 = VSET1(p434p1[8]);
  const __m512i vbmask = VSET1(BMASK);

  d0 = VADD(a0, a0); d1 = VADD(a1, a1); d2 = VADD(a2, a2);
  d3 = VADD(a3, a3); d4 = VADD(a4, a4); d5 = VADD(a5, a5);
  d6 = VADD(a6, a6); d7 = VADD(a7, a7); 

  // ---------------------------------------------------------------------------
  // integer multiplication 

  z0 = VMACLO(z0, a0, a0);
  h1 = VMACHI(h1, a0, a0);

  z1 = VMACLO(z1, d0, a1);
  h2 = VMACHI(h2, d0, a1);
   
  z2 = VMACLO(z2, d0, a2); z2 = VMACLO(z2, a1, a1); 
  h3 = VMACHI(h3, d0, a2); h3 = VMACHI(h3, a1, a1);

  z3 = VMACLO(z3, d0, a3); z3 = VMACLO(z3, d1, a2); 
  h4 = VMACHI(h4, d0, a3); h4 = VMACHI(h4, d1, a2); 
 
  z4 = VMACLO(z4, d0, a4); z4 = VMACLO(z4, d1, a3); z4 = VMACLO(z4, a2, a2);
  h5 = VMACHI(h5, d0, a4); h5 = VMACHI(h5, d1, a3); h5 = VMACHI(h5, a2, a2); 

  z5 = VMACLO(z5, d0, a5); z5 = VMACLO(z5, d1, a4); z5 = VMACLO(z5, d2, a3); 
  h6 = VMACHI(h6, d0, a5); h6 = VMACHI(h6, d1, a4); h6 = VMACHI(h6, d2, a3);   

  z6 = VMACLO(z6, d0, a6); z6 = VMACLO(z6, d1, a5); z6 = VMACLO(z6, d2, a4); 
  z6 = VMACLO(z6, a3, a3);
  h7 = VMACHI(h7, d0, a6); h7 = VMACHI(h7, d1, a5); h7 = VMACHI(h7, d2, a4); 
  h7 = VMACHI(h7, a3, a3);

  z7 = VMACLO(z7, d0, a7); z7 = VMACLO(z7, d1, a6); z7 = VMACLO(z7, d2, a5); 
  z7 = VMACLO(z7, d3, a4); 
  h8 = VMACHI(h8, d0, a7); h8 = VMACHI(h8, d1, a6); h8 = VMACHI(h8, d2, a5); 
  h8 = VMACHI(h8, d3, a4); 

  z8 = VMACLO(z8, d0, a8); z8 = VMACLO(z8, d1, a7); z8 = VMACLO(z8, d2, a6); 
  z8 = VMACLO(z8, d3, a5); z8 = VMACLO(z8, a4, a4);
  h9 = VMACHI(h9, d0, a8); h9 = VMACHI(h9, d1, a7); h9 = VMACHI(h9, d2, a6); 
  h9 = VMACHI(h9, d3, a5); h9 = VMACHI(h9, a4, a4);

  // some carry propagations

  z1 = VADD(z1, VADD(h1, h1)); 
  z2 = VADD(z2, VADD(h2, h2));
  z3 = VADD(z3, VADD(h3, h3)); 
  z4 = VADD(z4, VADD(h4, h4));

  z1 = VADD(z1, VSRA(z0, BRADIX)); z0 = VAND(z0, vbmask);
  z2 = VADD(z2, VSRA(z1, BRADIX)); z1 = VAND(z1, vbmask);
  z3 = VADD(z3, VSRA(z2, BRADIX)); z2 = VAND(z2, vbmask);
  z4 = VADD(z4, VSRA(z3, BRADIX)); z3 = VAND(z3, vbmask);

  // 2nd loop of integer squaring

  u = z0;
  z4 = VMACLO(z4, u, vp4); z5 = VMACLO(z5, u, vp5); z6 = VMACLO(z6, u, vp6); 
  z7 = VMACLO(z7, u, vp7); z8 = VMACLO(z8, u, vp8);
  h5 = VMACHI(h5, u, vp4); h6 = VMACHI(h6, u, vp5); h7 = VMACHI(h7, u, vp6);
  h8 = VMACHI(h8, u, vp7); h9 = VMACHI(h9, u, vp8);
  z5 = VADD(z5, VSRA(z4, BRADIX)); z4 = VAND(z4, vbmask); 
  z5 = VADD(z5, VADD(h5, h5));

  z9 = VMACLO(z9, d1, a8); z9 = VMACLO(z9, d2, a7); z9 = VMACLO(z9, d3, a6); 
  z9 = VMACLO(z9, d4, a5); 
  h10 = VMACHI(h10, d1, a8); h10 = VMACHI(h10, d2, a7); h10 = VMACHI(h10, d3, a6); 
  h10 = VMACHI(h10, d4, a5);

  u = z1; 
  z5 = VMACLO(z5, u, vp4); z6  = VMACLO(z6,  u, vp5); z7 = VMACLO(z7, u, vp6); 
  z8 = VMACLO(z8, u, vp7); z9  = VMACLO(z9,  u, vp8);
  h6 = VMACHI(h6, u, vp4); h7  = VMACHI(h7,  u, vp5); h8 = VMACHI(h8, u, vp6);
  h9 = VMACHI(h9, u, vp7); h10 = VMACHI(h10, u, vp8);
  z6 = VADD(z6, VSRA(z5, BRADIX)); z5 = VAND(z5, vbmask); 
  z6 = VADD(z6, VADD(h6, h6));

  z10 = VMACLO(z10, d2, a8); z10 = VMACLO(z10, d3, a7); z10 = VMACLO(z10, d4, a6); 
  z10 = VMACLO(z10, a5, a5); 
  h11 = VMACHI(h11, d2, a8); h11 = VMACHI(h11, d3, a7); h11 = VMACHI(h11, d4, a6);
  h11 = VMACHI(h11, a5, a5); 
 
  u = z2; 
  z6  = VMACLO(z6,  u, vp4); z7  = VMACLO(z7,  u, vp5); z8 = VMACLO(z8, u, vp6); 
  z9  = VMACLO(z9,  u, vp7); z10 = VMACLO(z10, u, vp8);
  h7  = VMACHI(h7,  u, vp4); h8  = VMACHI(h8,  u, vp5); h9 = VMACHI(h9, u, vp6);
  h10 = VMACHI(h10, u, vp7); h11 = VMACHI(h11, u, vp8);
  z7 = VADD(z7, VSRA(z6, BRADIX)); z6 = VAND(z6, vbmask); 
  z7 = VADD(z7, VADD(h7, h7));

  z11 = VMACLO(z11, d3, a8); z11 = VMACLO(z11, d4, a7); z11 = VMACLO(z11, d5, a6); 
  h12 = VMACHI(h12, d3, a8); h12 = VMACHI(h12, d4, a7); h12 = VMACHI(h12, d5, a6);

  u = z3; 
  z7  = VMACLO(z7,  u, vp4); z8  = VMACLO(z8,  u, vp5); z9  = VMACLO(z9,  u, vp6); 
  z10 = VMACLO(z10, u, vp7); z11 = VMACLO(z11, u, vp8);
  h8  = VMACHI(h8,  u, vp4); h9  = VMACHI(h9,  u, vp5); h10 = VMACHI(h10, u, vp6);
  h11 = VMACHI(h11, u, vp7); h12 = VMACHI(h12, u, vp8);
  z8 = VADD(z8, VSRA(z7, BRADIX)); z7 = VAND(z7, vbmask); 
  z8 = VADD(z8, VADD(h8, h8));

  z12 = VMACLO(z12, d4, a8); z12 = VMACLO(z12, d5, a7); z12 = VMACLO(z12, a6, a6); 
  h13 = VMACHI(h13, d4, a8); h13 = VMACHI(h13, d5, a7); h13 = VMACHI(h13, a6, a6); 

  u = z4; 
  z8  = VMACLO(z8,  u, vp4); z9  = VMACLO(z9,  u, vp5); z10 = VMACLO(z10, u, vp6); 
  z11 = VMACLO(z11, u, vp7); z12 = VMACLO(z12, u, vp8);
  h9  = VMACHI(h9,  u, vp4); h10 = VMACHI(h10, u, vp5); h11 = VMACHI(h11, u, vp6);
  h12 = VMACHI(h12, u, vp7); h13 = VMACHI(h13, u, vp8);
  z9 = VADD(z9, VSRA(z8, BRADIX)); z8 = VAND(z8, vbmask); 
  z9 = VADD(z9, VADD(h9, h9));  

  z13 = VMACLO(z13, d5, a8); z13 = VMACLO(z13, d6, a7); 
  h14 = VMACHI(h14, d5, a8); h14 = VMACHI(h14, d6, a7); 

  u = z5; 
  z9  = VMACLO(z9,  u, vp4); z10 = VMACLO(z10, u, vp5); z11 = VMACLO(z11, u, vp6); 
  z12 = VMACLO(z12, u, vp7); z13 = VMACLO(z13, u, vp8);
  h10 = VMACHI(h10, u, vp4); h11 = VMACHI(h11, u, vp5); h12 = VMACHI(h12, u, vp6);
  h13 = VMACHI(h13, u, vp7); h14 = VMACHI(h14, u, vp8);
  z10 = VADD(z10, VSRA(z9, BRADIX)); z9 = VAND(z9, vbmask); 
  z10 = VADD(z10, VADD(h10, h10)); 

  z14 = VMACLO(z14, d6, a8); z14 = VMACLO(z14, a7, a7);
  h15 = VMACHI(h15, d6, a8); h15 = VMACHI(h15, a7, a7); 

  u = z6; 
  z10 = VMACLO(z10, u, vp4); z11 = VMACLO(z11, u, vp5); z12 = VMACLO(z12, u, vp6); 
  z13 = VMACLO(z13, u, vp7); z14 = VMACLO(z14, u, vp8);
  h11 = VMACHI(h11, u, vp4); h12 = VMACHI(h12, u, vp5); h13 = VMACHI(h13, u, vp6);
  h14 = VMACHI(h14, u, vp7); h15 = VMACHI(h15, u, vp8);
  z11 = VADD(z11, VSRA(z10, BRADIX)); z10 = VAND(z10, vbmask); 
  z11 = VADD(z11, VADD(h11, h11));

  z15 = VMACLO(z15, d7, a8); 
  h16 = VMACHI(h16, d7, a8);

  u = z7; 
  z11 = VMACLO(z11, u, vp4); z12 = VMACLO(z12, u, vp5); z13 = VMACLO(z13, u, vp6); 
  z14 = VMACLO(z14, u, vp7); z15 = VMACLO(z15, u, vp8);
  h12 = VMACHI(h12, u, vp4); h13 = VMACHI(h13, u, vp5); h14 = VMACHI(h14, u, vp6);
  h15 = VMACHI(h15, u, vp7); h16 = VMACHI(h16, u, vp8);
  z12 = VADD(z12, VSRA(z11, BRADIX)); z11 = VAND(z11, vbmask); 
  z12 = VADD(z12, VADD(h12, h12));

  z16 = VMACLO(z16, a8, a8);
  h17 = VMACHI(h17, a8, a8); 

  u = z8; 
  z12 = VMACLO(z12, u, vp4); z13 = VMACLO(z13, u, vp5); z14 = VMACLO(z14, u, vp6); 
  z15 = VMACLO(z15, u, vp7); z16 = VMACLO(z16, u, vp8);
  h13 = VMACHI(h13, u, vp4); h14 = VMACHI(h14, u, vp5); h15 = VMACHI(h15, u, vp6);
  h16 = VMACHI(h16, u, vp7); h17 = VMACHI(h17, u, vp8);
  z13 = VADD(z13, VSRA(z12, BRADIX)); z12 = VAND(z12, vbmask);
  z13 = VADD(z13, VADD(h13, h13)); 

  z14 = VADD(z14, VADD(h14, h14));
  z15 = VADD(z15, VADD(h15, h15));
  z16 = VADD(z16, VADD(h16, h16));
  z17 = VADD(z17, VADD(h17, h17));

  z14 = VADD(z14, VSRA(z13, BRADIX)); z13 = VAND(z13, vbmask);
  z15 = VADD(z15, VSRA(z14, BRADIX)); z14 = VAND(z14, vbmask);
  z16 = VADD(z16, VSRA(z15, BRADIX)); z15 = VAND(z15, vbmask);
  z17 = VADD(z17, VSRA(z16, BRADIX)); z16 = VAND(z16, vbmask);

  r[0] =  z9; r[1] = z10; r[2] = z11; 
  r[3] = z12; r[4] = z13; r[5] = z14; 
  r[6] = z15; r[7] = z16; r[8] = z17;  
}

void fpinv_chain_mont(felm_t r)
{
  felm_t t[31], tt;
  int i, j;

  // Precomputed table
  fpsqr_mont(tt, r);
  fpmul_mont(t[0], r, tt);    
  for (i = 0; i <= 29; i++) fpmul_mont(t[i+1], t[i], tt);
  fpcopy(tt, r);
  for (i = 0; i < 7; i++) fpsqr_mont(tt, tt);
  fpmul_mont(tt, t[5], tt);
  for (i = 0; i < 10; i++) fpsqr_mont(tt, tt);
  fpmul_mont(tt, t[14], tt);
  for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
  fpmul_mont(tt, t[3], tt);
  for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
  fpmul_mont(tt, t[23], tt);
  for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
  fpmul_mont(tt, t[13], tt);
  for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
  fpmul_mont(tt, t[24], tt);
  for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
  fpmul_mont(tt, t[7], tt);  
  for (i = 0; i < 8; i++) fpsqr_mont(tt, tt); 
  fpmul_mont(tt, t[12], tt);
  for (i = 0; i < 8; i++) fpsqr_mont(tt, tt);
  fpmul_mont(tt, t[30], tt);
  for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
  fpmul_mont(tt, t[1], tt);
  for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
  fpmul_mont(tt, t[30], tt);    
  for (i = 0; i < 7; i++) fpsqr_mont(tt, tt);
  fpmul_mont(tt, t[21], tt);  
  for (i = 0; i < 9; i++) fpsqr_mont(tt, tt);
  fpmul_mont(tt, t[2], tt);
  for (i = 0; i < 9; i++) fpsqr_mont(tt, tt);
  fpmul_mont(tt, t[19], tt);
  for (i = 0; i < 9; i++) fpsqr_mont(tt, tt);
  fpmul_mont(tt, t[1], tt);
  for (i = 0; i < 7; i++) fpsqr_mont(tt, tt);
  fpmul_mont(tt, t[24], tt);
  for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
  fpmul_mont(tt, t[26], tt);
  for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
  fpmul_mont(tt, t[16], tt);    
  for (i = 0; i < 7; i++) fpsqr_mont(tt, tt);
  fpmul_mont(tt, t[10], tt);
  for (i = 0; i < 7; i++) fpsqr_mont(tt, tt);
  fpmul_mont(tt, t[6], tt);
  for (i = 0; i < 7; i++) fpsqr_mont(tt, tt);
  fpmul_mont(tt, t[0], tt);
  for (i = 0; i < 9; i++) fpsqr_mont(tt, tt);
  fpmul_mont(tt, t[20], tt);
  for (i = 0; i < 8; i++) fpsqr_mont(tt, tt);
  fpmul_mont(tt, t[9], tt);
  for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
  fpmul_mont(tt, t[25], tt);
  for (i = 0; i < 9; i++) fpsqr_mont(tt, tt);
  fpmul_mont(tt, t[30], tt);
  for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
  fpmul_mont(tt, t[26], tt);
  for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
  fpmul_mont(tt, r, tt);
  for (i = 0; i < 7; i++) fpsqr_mont(tt, tt);
  fpmul_mont(tt, t[28], tt);
  for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
  fpmul_mont(tt, t[6], tt);
  for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
  fpmul_mont(tt, t[10], tt);
  for (i = 0; i < 9; i++) fpsqr_mont(tt, tt);
  fpmul_mont(tt, t[22], tt);
  for (j = 0; j < 35; j++) {
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(tt, t[30], tt);
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
  __m512i a0 = a[0], a1 = a[1], a2 = a[2];
  __m512i a3 = a[3], a4 = a[4], a5 = a[5];
  __m512i a6 = a[6], a7 = a[7], a8 = a[8];
  __m512i b0 = b[0], b1 = b[1], b2 = b[2];
  __m512i b3 = b[3], b4 = b[4], b5 = b[5];
  __m512i b6 = b[6], b7 = b[7], b8 = b[8];
  __m512i r0, r1, r2, r3, r4, r5, r6, r7, r8;
  const __m512i vbmask = VSET1(BMASK); 

  // r = a + b
  r0 = VADD(a0, b0); r1 = VADD(a1, b1); r2 = VADD(a2, b2); 
  r3 = VADD(a3, b3); r4 = VADD(a4, b4); r5 = VADD(a5, b5);
  r6 = VADD(a6, b6); r7 = VADD(a7, b7); r8 = VADD(a8, b8);

  // carry propagation 
  r1 = VADD(r1, VSRA(r0, BRADIX)); r0 = VAND(r0, vbmask);
  r2 = VADD(r2, VSRA(r1, BRADIX)); r1 = VAND(r1, vbmask);
  r3 = VADD(r3, VSRA(r2, BRADIX)); r2 = VAND(r2, vbmask);
  r4 = VADD(r4, VSRA(r3, BRADIX)); r3 = VAND(r3, vbmask);
  r5 = VADD(r5, VSRA(r4, BRADIX)); r4 = VAND(r4, vbmask);
  r6 = VADD(r6, VSRA(r5, BRADIX)); r5 = VAND(r5, vbmask);
  r7 = VADD(r7, VSRA(r6, BRADIX)); r6 = VAND(r6, vbmask);
  r8 = VADD(r8, VSRA(r7, BRADIX)); r7 = VAND(r7, vbmask);

  r[0] = r0; r[1] = r1; r[2] = r2;
  r[3] = r3; r[4] = r4; r[5] = r5;
  r[6] = r6; r[7] = r7; r[8] = r8;
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
  __m512i a0  = a[0],  a1  = a[1],  a2  = a[2];
  __m512i a3  = a[3],  a4  = a[4],  a5  = a[5];
  __m512i a6  = a[6],  a7  = a[7],  a8  = a[8];
  __m512i a9  = a[9],  a10 = a[10], a11 = a[11];
  __m512i a12 = a[12], a13 = a[13], a14 = a[14];
  __m512i a15 = a[15], a16 = a[16], a17 = a[17];
  __m512i b0  = b[0],  b1  = b[1],  b2  = b[2];
  __m512i b3  = b[3],  b4  = b[4],  b5  = b[5];
  __m512i b6  = b[6],  b7  = b[7],  b8  = b[8];
  __m512i b9  = b[9],  b10 = b[10], b11 = b[11];
  __m512i b12 = b[12], b13 = b[13], b14 = b[14];
  __m512i b15 = b[15], b16 = b[16], b17 = b[17];
  __m512i r0  = r[0],  r1  = r[1],  r2  = r[2];
  __m512i r3  = r[3],  r4  = r[4],  r5  = r[5];
  __m512i r6  = r[6],  r7  = r[7],  r8  = r[8];
  __m512i r9  = r[9],  r10 = r[10], r11 = r[11];
  __m512i r12 = r[12], r13 = r[13], r14 = r[14];
  __m512i r15 = r[15], r16 = r[16], r17 = r[17];
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

  // carry propagation
  // will be done in reduction

  r[0]  = r0;  r[1]  = r1;  r[2]  = r2;
  r[3]  = r3;  r[4]  = r4;  r[5]  = r5;
  r[6]  = r6;  r[7]  = r7;  r[8]  = r8;
  r[9]  = r9;  r[10] = r10; r[11] = r11;
  r[12] = r12; r[13] = r13; r[14] = r14;
  r[15] = r15; r[16] = r16; r[17] = r17;
}

// integer subtraction followed by addition with p*468
//  r = a - b + (p*2^468) if a-b<0; otherwise r = a - b
static void mp_subaddfast(__m512i *r, const __m512i *a, const __m512i *b)
{
  __m512i a0  = a[0],  a1  = a[1],  a2  = a[2];
  __m512i a3  = a[3],  a4  = a[4],  a5  = a[5];
  __m512i a6  = a[6],  a7  = a[7],  a8  = a[8];
  __m512i a9  = a[9],  a10 = a[10], a11 = a[11];
  __m512i a12 = a[12], a13 = a[13], a14 = a[14];
  __m512i a15 = a[15], a16 = a[16], a17 = a[17];
  __m512i b0  = b[0],  b1  = b[1],  b2  = b[2];
  __m512i b3  = b[3],  b4  = b[4],  b5  = b[5];
  __m512i b6  = b[6],  b7  = b[7],  b8  = b[8];
  __m512i b9  = b[9],  b10 = b[10], b11 = b[11];
  __m512i b12 = b[12], b13 = b[13], b14 = b[14];
  __m512i b15 = b[15], b16 = b[16], b17 = b[17];
  __m512i r0, r1,  r2,  r3,  r4,  r5,  r6,  r7,  r8;
  __m512i r9, r10, r11, r12, r13, r14, r15, r16, r17, smask;
  const __m512i vp0 = VSET1(p434[0]), vp1 = VSET1(p434[1]), vp2 = VSET1(p434[2]);
  const __m512i vp3 = VSET1(p434[3]), vp4 = VSET1(p434[4]), vp5 = VSET1(p434[5]);
  const __m512i vp6 = VSET1(p434[6]), vp7 = VSET1(p434[7]), vp8 = VSET1(p434[8]);
  const __m512i vbmask = VSET1(BMASK);

  // r = a - b
  r0  = VSUB(a0,  b0);  r1  = VSUB(a1,  b1);  r2  = VSUB(a2,  b2);
  r3  = VSUB(a3,  b3);  r4  = VSUB(a4,  b4);  r5  = VSUB(a5,  b5);
  r6  = VSUB(a6,  b6);  r7  = VSUB(a7,  b7);  r8  = VSUB(a8,  b8);
  r9  = VSUB(a9,  b9);  r10 = VSUB(a10, b10); r11 = VSUB(a11, b11);
  r12 = VSUB(a12, b12); r13 = VSUB(a13, b13); r14 = VSUB(a14, b14);
  r15 = VSUB(a15, b15); r16 = VSUB(a16, b16); r17 = VSUB(a17, b17);

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

  // if r is non-negative, smask = all-0 
  // if r is     negative, smask = all-1
  smask = VSRA(r17, 63);
  // r = r + ((p*2^468) & smask)
  r9  = VADD(r9,  VAND(vp0, smask)); r10 = VADD(r10, VAND(vp1, smask));
  r11 = VADD(r11, VAND(vp2, smask)); r12 = VADD(r12, VAND(vp3, smask)); 
  r13 = VADD(r13, VAND(vp4, smask)); r14 = VADD(r14, VAND(vp5, smask)); 
  r15 = VADD(r15, VAND(vp6, smask)); r16 = VADD(r16, VAND(vp7, smask)); 
  r17 = VADD(r17, VAND(vp8, smask)); 

  // carry propagation
  // will be done in reduction 

  r[0]  = r0;  r[1]  = r1;  r[2]  = r2;
  r[3]  = r3;  r[4]  = r4;  r[5]  = r5;
  r[6]  = r6;  r[7]  = r7;  r[8]  = r8;
  r[9]  = r9;  r[10] = r10; r[11] = r11;
  r[12] = r12; r[13] = r13; r[14] = r14;
  r[15] = r15; r[16] = r16; r[17] = r17;
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
  mp_subaddfast(tt1, tt1, tt2);         // tt1 = a0b0-a1b1      in [0, p*2^468] (only be used here)
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
