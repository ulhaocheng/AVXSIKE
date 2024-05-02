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
  __m512i a0 = a[0], a1 = a[1], a2  = a[2],  a3  = a[3];
  __m512i a4 = a[4], a5 = a[5], a6  = a[6],  a7  = a[7];
  __m512i a8 = a[8], a9 = a[9], a10 = a[10], a11 = a[11];
  __m512i d0, d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, d11;
  __m512i r0,  r1,  r2,  r3,  r4,  r5,  r6,  r7,  r8,  r9,  r10, r11, u;
  __m512i z0  = VZERO, z1  = VZERO, z2  = VZERO, z3  = VZERO; 
  __m512i z4  = VZERO, z5  = VZERO, z6  = VZERO, z7  = VZERO; 
  __m512i z8  = VZERO, z9  = VZERO, z10 = VZERO, z11 = VZERO; 
  __m512i z12 = VZERO, z13 = VZERO, z14 = VZERO, z15 = VZERO;
  __m512i z16 = VZERO, z17 = VZERO, z18 = VZERO, z19 = VZERO;
  __m512i z20 = VZERO, z21 = VZERO, z22 = VZERO, z23 = VZERO;
  __m512i h0  = VZERO, h1  = VZERO, h2  = VZERO, h3  = VZERO; 
  __m512i h4  = VZERO, h5  = VZERO, h6  = VZERO, h7  = VZERO; 
  __m512i h8  = VZERO, h9  = VZERO, h10 = VZERO, h11 = VZERO; 
  __m512i h12 = VZERO, h13 = VZERO, h14 = VZERO, h15 = VZERO;
  __m512i h16 = VZERO, h17 = VZERO, h18 = VZERO, h19 = VZERO;
  __m512i h20 = VZERO, h21 = VZERO, h22 = VZERO, h23 = VZERO;
  const __m512i vp5  = VSET1(p610p1[5]),  vp6 = VSET1(p610p1[6]), vp7  = VSET1(p610p1[7]);
  const __m512i vp8  = VSET1(p610p1[8]),  vp9 = VSET1(p610p1[9]), vp10 = VSET1(p610p1[10]);
  const __m512i vp11 = VSET1(p610p1[11]), vbmask = VSET1(BMASK);

  d0 = VADD(a0, a0); d1 = VADD(a1, a1); d2 = VADD(a2, a2);
  d3 = VADD(a3, a3); d4 = VADD(a4, a4); d5 = VADD(a5, a5);
  d6 = VADD(a6, a6); d7 = VADD(a7, a7); d8 = VADD(a8, a8);
  d9 = VADD(a9, a9); d10 = VADD(a10, a10);

  // ---------------------------------------------------------------------------
  // 1st loop of integer squaring

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

  z9 = VMACLO(z9, d0, a9); z9 = VMACLO(z9, d1, a8); z9 = VMACLO(z9, d2, a7); 
  z9 = VMACLO(z9, d3, a6); z9 = VMACLO(z9, d4, a5); 
  h10 = VMACHI(h10, d0, a9); h10 = VMACHI(h10, d1, a8); h10 = VMACHI(h10, d2, a7); 
  h10 = VMACHI(h10, d3, a6); h10 = VMACHI(h10, d4, a5);

  z10 = VMACLO(z10, d0, a10); z10 = VMACLO(z10, d1, a9); z10 = VMACLO(z10, d2, a8); 
  z10 = VMACLO(z10, d3, a7); z10 = VMACLO(z10, d4, a6); z10 = VMACLO(z10, a5, a5); 
  h11 = VMACHI(h11, d0, a10); h11 = VMACHI(h11, d1, a9); h11 = VMACHI(h11, d2, a8); 
  h11 = VMACHI(h11, d3, a7); h11 = VMACHI(h11, d4, a6); h11 = VMACHI(h11, a5, a5); 

  
  z11 = VMACLO(z11, d0, a11); z11 = VMACLO(z11, d1, a10); z11 = VMACLO(z11, d2, a9); 
  z11 = VMACLO(z11, d3, a8); z11 = VMACLO(z11, d4, a7); z11 = VMACLO(z11, d5, a6); 
  h12 = VMACHI(h12, d0, a11); h12 = VMACHI(h12, d1, a10); h12 = VMACHI(h12, d2, a9); 
  h12 = VMACHI(h12, d3, a8); h12 = VMACHI(h12, d4, a7); h12 = VMACHI(h12, d5, a6); 

  // some carry propagations
  z1 = VADD(z1, VADD(h1, h1)); 
  z2 = VADD(z2, VADD(h2, h2));
  z3 = VADD(z3, VADD(h3, h3)); 
  z4 = VADD(z4, VADD(h4, h4));
  z5 = VADD(z5, VADD(h5, h5));

  z1 = VADD(z1, VSRA(z0, BRADIX)); z0 = VAND(z0, vbmask);
  z2 = VADD(z2, VSRA(z1, BRADIX)); z1 = VAND(z1, vbmask);
  z3 = VADD(z3, VSRA(z2, BRADIX)); z2 = VAND(z2, vbmask);
  z4 = VADD(z4, VSRA(z3, BRADIX)); z3 = VAND(z3, vbmask);
  z5 = VADD(z5, VSRA(z4, BRADIX)); z4 = VAND(z4, vbmask);

  // 2nd loop of integer squaring 

  u = z0;
  z5  = VMACLO(z5,  u, vp5); z6  = VMACLO(z6,  u, vp6); z7  = VMACLO(z7,  u, vp7); 
  z8  = VMACLO(z8,  u, vp8); z9  = VMACLO(z9,  u, vp9); z10 = VMACLO(z10, u, vp10);
  z11 = VMACLO(z11, u, vp11);
  h6  = VMACHI(h6,  u, vp5); h7  = VMACHI(h7,  u, vp6); h8  = VMACHI(h8,  u, vp7); 
  h9  = VMACHI(h9,  u, vp8); h10 = VMACHI(h10, u, vp9); h11 = VMACHI(h11, u, vp10);
  h12 = VMACHI(h12, u, vp11);
  z6 = VADD(z6, VSRA(z5, BRADIX)); z5 = VAND(z5, vbmask); 
  z6 = VADD(z6, VADD(h6, h6));
  
  z12 = VMACLO(z12, d1, a11); z12 = VMACLO(z12, d2, a10); z12 = VMACLO(z12, d3, a9); 
  z12 = VMACLO(z12, d4, a8); z12 = VMACLO(z12, d5, a7); z12 = VMACLO(z12, a6, a6); 
  h13 = VMACHI(h13, d1, a11); h13 = VMACHI(h13, d2, a10); h13 = VMACHI(h13, d3, a9); 
  h13 = VMACHI(h13, d4, a8); h13 = VMACHI(h13, d5, a7); h13 = VMACHI(h13, a6, a6); 

  u = z1; 
  z6  = VMACLO(z6,  u, vp5); z7  = VMACLO(z7,  u, vp6); z8  = VMACLO(z8,  u, vp7); 
  z9  = VMACLO(z9,  u, vp8); z10 = VMACLO(z10, u, vp9); z11 = VMACLO(z11, u, vp10);
  z12 = VMACLO(z12, u, vp11);
  h7  = VMACHI(h7,  u, vp5); h8  = VMACHI(h8,  u, vp6); h9  = VMACHI(h9,  u, vp7); 
  h10 = VMACHI(h10, u, vp8); h11 = VMACHI(h11, u, vp9); h12 = VMACHI(h12, u, vp10);
  h13 = VMACHI(h13, u, vp11);
  z7 = VADD(z7, VSRA(z6, BRADIX)); z6 = VAND(z6, vbmask); 
  z7 = VADD(z7, VADD(h7, h7));

  z13 = VMACLO(z13, d2, a11); z13 = VMACLO(z13, d3, a10); z13 = VMACLO(z13, d4, a9); 
  z13 = VMACLO(z13, d5, a8); z13 = VMACLO(z13, d6, a7); 
  h14 = VMACHI(h14, d2, a11); h14 = VMACHI(h14, d3, a10); h14 = VMACHI(h14, d4, a9); 
  h14 = VMACHI(h14, d5, a8); h14 = VMACHI(h14, d6, a7); 

  u = z2; 
  z7  = VMACLO(z7,  u, vp5); z8  = VMACLO(z8,  u, vp6); z9  = VMACLO(z9,  u, vp7); 
  z10 = VMACLO(z10, u, vp8); z11 = VMACLO(z11, u, vp9); z12 = VMACLO(z12, u, vp10);
  z13 = VMACLO(z13, u, vp11);
  h8  = VMACHI(h8,  u, vp5); h9  = VMACHI(h9,  u, vp6); h10 = VMACHI(h10, u, vp7); 
  h11 = VMACHI(h11, u, vp8); h12 = VMACHI(h12, u, vp9); h13 = VMACHI(h13, u, vp10);
  h14 = VMACHI(h14, u, vp11);
  z8 = VADD(z8, VSRA(z7, BRADIX)); z7 = VAND(z7, vbmask); 
  z8 = VADD(z8, VADD(h8, h8));

  z14 = VMACLO(z14, d3, a11); z14 = VMACLO(z14, d4, a10); z14 = VMACLO(z14, d5, a9); 
  z14 = VMACLO(z14, d6, a8); z14 = VMACLO(z14, a7, a7); 
  h15 = VMACHI(h15, d3, a11); h15 = VMACHI(h15, d4, a10); h15 = VMACHI(h15, d5, a9); 
  h15 = VMACHI(h15, d6, a8); h15 = VMACHI(h15, a7, a7); 

  u = z3; 
  z8  = VMACLO(z8,  u, vp5); z9  = VMACLO(z9,  u, vp6); z10 = VMACLO(z10, u, vp7); 
  z11 = VMACLO(z11, u, vp8); z12 = VMACLO(z12, u, vp9); z13 = VMACLO(z13, u, vp10);
  z14 = VMACLO(z14, u, vp11);
  h9  = VMACHI(h9,  u, vp5); h10 = VMACHI(h10, u, vp6); h11 = VMACHI(h11, u, vp7); 
  h12 = VMACHI(h12, u, vp8); h13 = VMACHI(h13, u, vp9); h14 = VMACHI(h14, u, vp10);
  h15 = VMACHI(h15, u, vp11);
  z9 = VADD(z9, VSRA(z8, BRADIX)); z8 = VAND(z8, vbmask); 
  z9 = VADD(z9, VADD(h9, h9));

  z15 = VMACLO(z15, d4, a11); z15 = VMACLO(z15, d5, a10); z15 = VMACLO(z15, d6, a9); 
  z15 = VMACLO(z15, d7, a8); 
  h16 = VMACHI(h16, d4, a11); h16 = VMACHI(h16, d5, a10); h16 = VMACHI(h16, d6, a9); 
  h16 = VMACHI(h16, d7, a8); 

  u = z4; 
  z9  = VMACLO(z9,  u, vp5); z10 = VMACLO(z10, u, vp6); z11 = VMACLO(z11, u, vp7); 
  z12 = VMACLO(z12, u, vp8); z13 = VMACLO(z13, u, vp9); z14 = VMACLO(z14, u, vp10);
  z15 = VMACLO(z15, u, vp11);
  h10 = VMACHI(h10, u, vp5); h11 = VMACHI(h11, u, vp6); h12 = VMACHI(h12, u, vp7); 
  h13 = VMACHI(h13, u, vp8); h14 = VMACHI(h14, u, vp9); h15 = VMACHI(h15, u, vp10);
  h16 = VMACHI(h16, u, vp11);
  z10 = VADD(z10, VSRA(z9, BRADIX)); z9 = VAND(z9, vbmask); 
  z10 = VADD(z10, VADD(h10, h10));

  z16 = VMACLO(z16, d5, a11); z16 = VMACLO(z16, d6, a10); z16 = VMACLO(z16, d7, a9); 
  z16 = VMACLO(z16, a8, a8); 
  h17 = VMACHI(h17, d5, a11); h17 = VMACHI(h17, d6, a10); h17 = VMACHI(h17, d7, a9); 
  h17 = VMACHI(h17, a8, a8); 

  u = z5; 
  z10 = VMACLO(z10, u, vp5); z11 = VMACLO(z11, u, vp6); z12 = VMACLO(z12, u, vp7); 
  z13 = VMACLO(z13, u, vp8); z14 = VMACLO(z14, u, vp9); z15 = VMACLO(z15, u, vp10);
  z16 = VMACLO(z16, u, vp11);
  h11 = VMACHI(h11, u, vp5); h12 = VMACHI(h12, u, vp6); h13 = VMACHI(h13, u, vp7); 
  h14 = VMACHI(h14, u, vp8); h15 = VMACHI(h15, u, vp9); h16 = VMACHI(h16, u, vp10);
  h17 = VMACHI(h17, u, vp11);
  z11 = VADD(z11, VSRA(z10, BRADIX)); z10 = VAND(z10, vbmask); 
  z11 = VADD(z11, VADD(h11, h11));

  z17 = VMACLO(z17, d6, a11); z17 = VMACLO(z17, d7, a10); z17 = VMACLO(z17, d8, a9); 
  h18 = VMACHI(h18, d6, a11); h18 = VMACHI(h18, d7, a10); h18 = VMACHI(h18, d8, a9); 

  u = z6; 
  z11 = VMACLO(z11, u, vp5); z12 = VMACLO(z12, u, vp6); z13 = VMACLO(z13, u, vp7); 
  z14 = VMACLO(z14, u, vp8); z15 = VMACLO(z15, u, vp9); z16 = VMACLO(z16, u, vp10);
  z17 = VMACLO(z17, u, vp11);
  h12 = VMACHI(h12, u, vp5); h13 = VMACHI(h13, u, vp6); h14 = VMACHI(h14, u, vp7); 
  h15 = VMACHI(h15, u, vp8); h16 = VMACHI(h16, u, vp9); h17 = VMACHI(h17, u, vp10);
  h18 = VMACHI(h18, u, vp11);
  z12 = VADD(z12, VSRA(z11, BRADIX)); z11 = VAND(z11, vbmask); 
  z12 = VADD(z12, VADD(h12, h12));

  z18 = VMACLO(z18, d7, a11); z18 = VMACLO(z18, d8, a10); z18 = VMACLO(z18, a9, a9);
  h19 = VMACHI(h19, d7, a11); h19 = VMACHI(h19, d8, a10); h19 = VMACHI(h19, a9, a9);

  u = z7; 
  z12 = VMACLO(z12, u, vp5); z13 = VMACLO(z13, u, vp6); z14 = VMACLO(z14, u, vp7); 
  z15 = VMACLO(z15, u, vp8); z16 = VMACLO(z16, u, vp9); z17 = VMACLO(z17, u, vp10);
  z18 = VMACLO(z18, u, vp11);
  h13 = VMACHI(h13, u, vp5); h14 = VMACHI(h14, u, vp6); h15 = VMACHI(h15, u, vp7); 
  h16 = VMACHI(h16, u, vp8); h17 = VMACHI(h17, u, vp9); h18 = VMACHI(h18, u, vp10);
  h19 = VMACHI(h19, u, vp11);
  z13 = VADD(z13, VSRA(z12, BRADIX)); z12 = VAND(z12, vbmask); 
  z13 = VADD(z13, VADD(h13, h13));

  z19 = VMACLO(z19, d8, a11); z19 = VMACLO(z19, d9, a10);
  h20 = VMACHI(h20, d8, a11); h20 = VMACHI(h20, d9, a10);

  u = z8; 
  z13 = VMACLO(z13, u, vp5); z14 = VMACLO(z14, u, vp6); z15 = VMACLO(z15, u, vp7); 
  z16 = VMACLO(z16, u, vp8); z17 = VMACLO(z17, u, vp9); z18 = VMACLO(z18, u, vp10);
  z19 = VMACLO(z19, u, vp11);
  h14 = VMACHI(h14, u, vp5); h15 = VMACHI(h15, u, vp6); h16 = VMACHI(h16, u, vp7); 
  h17 = VMACHI(h17, u, vp8); h18 = VMACHI(h18, u, vp9); h19 = VMACHI(h19, u, vp10);
  h20 = VMACHI(h20, u, vp11);
  z14 = VADD(z14, VSRA(z13, BRADIX)); z13 = VAND(z13, vbmask);
  z14 = VADD(z14, VADD(h14, h14));

  z20 = VMACLO(z20, d9, a11); z20 = VMACLO(z20, a10, a10);
  h21 = VMACHI(h21, d9, a11); h21 = VMACHI(h21, a10, a10);

  u = z9;
  z14 = VMACLO(z14, u, vp5); z15 = VMACLO(z15, u, vp6); z16 = VMACLO(z16, u, vp7); 
  z17 = VMACLO(z17, u, vp8); z18 = VMACLO(z18, u, vp9); z19 = VMACLO(z19, u, vp10);
  z20 = VMACLO(z20, u, vp11);
  h15 = VMACHI(h15, u, vp5); h16 = VMACHI(h16, u, vp6); h17 = VMACHI(h17, u, vp7); 
  h18 = VMACHI(h18, u, vp8); h19 = VMACHI(h19, u, vp9); h20 = VMACHI(h20, u, vp10);
  h21 = VMACHI(h21, u, vp11);
  z15 = VADD(z15, VSRA(z14, BRADIX)); z14 = VAND(z14, vbmask);
  z15 = VADD(z15, VADD(h15, h15));

  z21 = VMACLO(z21, d10, a11); 
  h22 = VMACHI(h22, d10, a11); 

  u = z10;
  z15 = VMACLO(z15, u, vp5); z16 = VMACLO(z16, u, vp6); z17 = VMACLO(z17, u, vp7); 
  z18 = VMACLO(z18, u, vp8); z19 = VMACLO(z19, u, vp9); z20 = VMACLO(z20, u, vp10);
  z21 = VMACLO(z21, u, vp11);
  h16 = VMACHI(h16, u, vp5); h17 = VMACHI(h17, u, vp6); h18 = VMACHI(h18, u, vp7); 
  h19 = VMACHI(h19, u, vp8); h20 = VMACHI(h20, u, vp9); h21 = VMACHI(h21, u, vp10);
  h22 = VMACHI(h22, u, vp11);
  z16 = VADD(z16, VSRA(z15, BRADIX)); z15 = VAND(z15, vbmask);
  z16 = VADD(z16, VADD(h16, h16));

  z22 = VMACLO(z22, a11, a11); 
  h23 = VMACHI(h23, a11, a11); 

  u = z11;
  z16 = VMACLO(z16, u, vp5); z17 = VMACLO(z17, u, vp6); z18 = VMACLO(z18, u, vp7); 
  z19 = VMACLO(z19, u, vp8); z20 = VMACLO(z20, u, vp9); z21 = VMACLO(z21, u, vp10);
  z22 = VMACLO(z22, u, vp11);
  h17 = VMACHI(h17, u, vp5); h18 = VMACHI(h18, u, vp6); h19 = VMACHI(h19, u, vp7); 
  h20 = VMACHI(h20, u, vp8); h21 = VMACHI(h21, u, vp9); h22 = VMACHI(h22, u, vp10);
  h23 = VMACHI(h23, u, vp11);
  z17 = VADD(z17, VSRA(z16, BRADIX)); z16 = VAND(z16, vbmask);
  z17 = VADD(z17, VADD(h17, h17));

  z18 = VADD(z18, VADD(h18, h18));
  z19 = VADD(z19, VADD(h19, h19));
  z20 = VADD(z20, VADD(h20, h20));
  z21 = VADD(z21, VADD(h21, h21));
  z22 = VADD(z22, VADD(h22, h22));
  z23 = VADD(z23, VADD(h23, h23));

  z18 = VADD(z18, VSRA(z17, BRADIX)); z17 = VAND(z17, vbmask);
  z19 = VADD(z19, VSRA(z18, BRADIX)); z18 = VAND(z18, vbmask);
  z20 = VADD(z20, VSRA(z19, BRADIX)); z19 = VAND(z19, vbmask);
  z21 = VADD(z21, VSRA(z20, BRADIX)); z20 = VAND(z20, vbmask);
  z22 = VADD(z22, VSRA(z21, BRADIX)); z21 = VAND(z21, vbmask);
  z23 = VADD(z23, VSRA(z22, BRADIX)); z22 = VAND(z22, vbmask);

  r[0] = z12; r[1] = z13; r[2]  = z14;  r[3] = z15; 
  r[4] = z16; r[5] = z17; r[6]  = z18;  r[7] = z19; 
  r[8] = z20; r[9] = z21; r[10] = z22; r[11] = z23;
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
  for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
  fpmul_mont(tt, t[6], tt);
  for (i = 0; i < 7; i++) fpsqr_mont(tt, tt);
  fpmul_mont(tt, t[30], tt);
  for (i = 0; i < 7; i++) fpsqr_mont(tt, tt);
  fpmul_mont(tt, t[25], tt);
  for (i = 0; i < 8; i++) fpsqr_mont(tt, tt);
  fpmul_mont(tt, t[28], tt);
  for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
  fpmul_mont(tt, t[7], tt);
  for (i = 0; i < 11; i++) fpsqr_mont(tt, tt);
  fpmul_mont(tt, t[11], tt);
  for (i = 0; i < 8; i++) fpsqr_mont(tt, tt);
  fpmul_mont(tt, r, tt);
  for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
  fpmul_mont(tt, t[0], tt);
  for (i = 0; i < 8; i++) fpsqr_mont(tt, tt);
  fpmul_mont(tt, t[3], tt);
  for (i = 0; i < 7; i++) fpsqr_mont(tt, tt);
  fpmul_mont(tt, t[16], tt);
  for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
  fpmul_mont(tt, t[24], tt);
  for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
  fpmul_mont(tt, t[28], tt);
  for (i = 0; i < 9; i++) fpsqr_mont(tt, tt);
  fpmul_mont(tt, t[16], tt);
  for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
  fpmul_mont(tt, t[4], tt);
  for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
  fpmul_mont(tt, t[3], tt);
  for (i = 0; i < 7; i++) fpsqr_mont(tt, tt);
  fpmul_mont(tt, t[20], tt);
  for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
  fpmul_mont(tt, t[11], tt);
  for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
  fpmul_mont(tt, t[14], tt);
  for (i = 0; i < 7; i++) fpsqr_mont(tt, tt);
  fpmul_mont(tt, t[15], tt);
  for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
  fpmul_mont(tt, t[0], tt);
  for (i = 0; i < 9; i++) fpsqr_mont(tt, tt);
  fpmul_mont(tt, t[15], tt);
  for (i = 0; i < 8; i++) fpsqr_mont(tt, tt);
  fpmul_mont(tt, t[19], tt);
  for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
  fpmul_mont(tt, t[9], tt);
  for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
  fpmul_mont(tt, t[5], tt);
  for (i = 0; i < 7; i++) fpsqr_mont(tt, tt);
  fpmul_mont(tt, t[27], tt);
  for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
  fpmul_mont(tt, t[28], tt);
  for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
  fpmul_mont(tt, t[29], tt);
  for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
  fpmul_mont(tt, t[1], tt);
  for (i = 0; i < 9; i++) fpsqr_mont(tt, tt);
  fpmul_mont(tt, t[3], tt);
  for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
  fpmul_mont(tt, t[2], tt);
  for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
  fpmul_mont(tt, t[30], tt);
  for (i = 0; i < 8; i++) fpsqr_mont(tt, tt);
  fpmul_mont(tt, t[25], tt);
  for (i = 0; i < 7; i++) fpsqr_mont(tt, tt);
  fpmul_mont(tt, t[28], tt);
  for (i = 0; i < 9; i++) fpsqr_mont(tt, tt);
  fpmul_mont(tt, t[22], tt);
  for (i = 0; i < 8; i++) fpsqr_mont(tt, tt);
  fpmul_mont(tt, t[3], tt);
  for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
  fpmul_mont(tt, t[22], tt);
  for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
  fpmul_mont(tt, t[7], tt);
  for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
  fpmul_mont(tt, t[9], tt);
  for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
  fpmul_mont(tt, t[4], tt);
  for (i = 0; i < 7; i++) fpsqr_mont(tt, tt);
  fpmul_mont(tt, t[20], tt);
  for (i = 0; i < 11; i++) fpsqr_mont(tt, tt);
  fpmul_mont(tt, t[10], tt);
  for (i = 0; i < 8; i++) fpsqr_mont(tt, tt);
  fpmul_mont(tt, t[26], tt);
  for (i = 0; i < 11; i++) fpsqr_mont(tt, tt);
  fpmul_mont(tt, t[2], tt);
  for (j = 0; j < 50; j++) {
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
  __m512i a0 = a[0], a1 = a[1], a2  = a[2],  a3  = a[3];
  __m512i a4 = a[4], a5 = a[5], a6  = a[6],  a7  = a[7];
  __m512i a8 = a[8], a9 = a[9], a10 = a[10], a11 = a[11];
  __m512i b0 = b[0], b1 = b[1], b2  = b[2],  b3  = b[3]; 
  __m512i b4 = b[4], b5 = b[5], b6  = b[6],  b7  = b[7];
  __m512i b8 = b[8], b9 = b[9], b10 = b[10], b11 = b[11];
  __m512i r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11;
  const __m512i vbmask = VSET1(BMASK); 

  // r = a + b
  r0 = VADD(a0, b0); r1  = VADD(a1,  b1);  r2  = VADD(a2, b2); 
  r3 = VADD(a3, b3); r4  = VADD(a4,  b4);  r5  = VADD(a5, b5);
  r6 = VADD(a6, b6); r7  = VADD(a7,  b7);  r8  = VADD(a8, b8);
  r9 = VADD(a9, b9); r10 = VADD(a10, b10); r11 = VADD(a11, b11);

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

  r[0] = r0; r[1] = r1; r[2]  = r2;  r[3]  = r3; 
  r[4] = r4; r[5] = r5; r[6]  = r6;  r[7]  = r7; 
  r[8] = r8; r[9] = r9; r[10] = r10; r[11] = r11;
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

// [0, 4p] -> [0, p]
static void fpcorrection_p4(__m512i *r, const __m512i *a)
{
  __m512i a0 = a[0], a1 = a[1], a2  = a[2],  a3  = a[3];
  __m512i a4 = a[4], a5 = a[5], a6  = a[6],  a7  = a[7];
  __m512i a8 = a[8], a9 = a[9], a10 = a[10], a11 = a[11];
  __m512i r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, smask;
  const __m512i vp0  = VSET1(p610x2[0]),  vp1  = VSET1(p610x2[1]);
  const __m512i vp2  = VSET1(p610x2[2]),  vp3  = VSET1(p610x2[3]);
  const __m512i vp4  = VSET1(p610x2[4]),  vp5  = VSET1(p610x2[5]);
  const __m512i vp6  = VSET1(p610x2[6]),  vp7  = VSET1(p610x2[7]);
  const __m512i vp8  = VSET1(p610x2[8]),  vp9  = VSET1(p610x2[9]);
  const __m512i vp10 = VSET1(p610x2[10]), vp11 = VSET1(p610x2[11]);
  const __m512i tp0  = VSET1(p610[0]),  tp1  = VSET1(p610[1]);
  const __m512i tp2  = VSET1(p610[2]),  tp3  = VSET1(p610[3]);
  const __m512i tp4  = VSET1(p610[4]),  tp5  = VSET1(p610[5]);
  const __m512i tp6  = VSET1(p610[6]),  tp7  = VSET1(p610[7]);
  const __m512i tp8  = VSET1(p610[8]),  tp9  = VSET1(p610[9]);
  const __m512i tp10 = VSET1(p610[10]), tp11 = VSET1(p610[11]);
  const __m512i vbmask = VSET1(BMASK); 

  // a - 2p
  r0 = VSUB(a0, vp0); r1  = VSUB(a1,  vp1);  r2  = VSUB(a2, vp2);
  r3 = VSUB(a3, vp3); r4  = VSUB(a4,  vp4);  r5  = VSUB(a5, vp5);
  r6 = VSUB(a6, vp6); r7  = VSUB(a7,  vp7);  r8  = VSUB(a8, vp8);
  r9 = VSUB(a9, vp9); r10 = VSUB(a10, vp10); r11 = VSUB(a11, vp11);

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

  // if r is non-negative, smask = all-0 
  // if r is     negative, smask = all-1
  smask = VSRA(r11, 63);
  // r = r + (2p & smask)
  r0  = VADD(r0,  VAND(vp0,  smask)); r1  = VADD(r1,  VAND(vp1,  smask)); 
  r2  = VADD(r2,  VAND(vp2,  smask)); r3  = VADD(r3,  VAND(vp3,  smask)); 
  r4  = VADD(r4,  VAND(vp4,  smask)); r5  = VADD(r5,  VAND(vp5,  smask)); 
  r6  = VADD(r6,  VAND(vp6,  smask)); r7  = VADD(r7,  VAND(vp7,  smask)); 
  r8  = VADD(r8,  VAND(vp8,  smask)); r9  = VADD(r9,  VAND(vp9,  smask)); 
  r10 = VADD(r10, VAND(vp10, smask)); r11 = VADD(r11, VAND(vp11, smask)); 

  // r - p
  r0 = VSUB(r0, tp0); r1  = VSUB(r1,  tp1);  r2  = VSUB(r2, tp2);
  r3 = VSUB(r3, tp3); r4  = VSUB(r4,  tp4);  r5  = VSUB(r5, tp5);
  r6 = VSUB(r6, tp6); r7  = VSUB(r7,  tp7);  r8  = VSUB(r8, tp8);
  r9 = VSUB(r9, tp9); r10 = VSUB(r10, tp10); r11 = VSUB(r11, tp11);

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

  smask = VSRA(r11, 63);
  // r = r + (p & smask)
  r0  = VADD(r0,  VAND(tp0,  smask)); r1  = VADD(r1,  VAND(tp1,  smask)); 
  r2  = VADD(r2,  VAND(tp2,  smask)); r3  = VADD(r3,  VAND(tp3,  smask)); 
  r4  = VADD(r4,  VAND(tp4,  smask)); r5  = VADD(r5,  VAND(tp5,  smask)); 
  r6  = VADD(r6,  VAND(tp6,  smask)); r7  = VADD(r7,  VAND(tp7,  smask)); 
  r8  = VADD(r8,  VAND(tp8,  smask)); r9  = VADD(r9,  VAND(tp9,  smask)); 
  r10 = VADD(r10, VAND(tp10, smask)); r11 = VADD(r11, VAND(tp11, smask)); 

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

  r[0] = r0; r[1] = r1; r[2]  = r2;  r[3]  = r3; 
  r[4] = r4; r[5] = r5; r[6]  = r6;  r[7]  = r7; 
  r[8] = r8; r[9] = r9; r[10] = r10; r[11] = r11;
}

// GF(p^2) squaring using Montgomery arithmetic r = a^2 in GF(p^2)
void fp2sqr_mont(f2elm_t r, const f2elm_t a)
{
  felm_t t1, t2, t3;
  f2elm_t ta;

  fpcorrection_p4(ta[0], a[0]);
  fpcorrection_p4(ta[1], a[1]);
  
  mp_addfast(t1, ta[0], ta[1]);           // t1 = a0+a1           in [0, 4p]
  // mp_sub_p4(t2, ta[0], ta[1]);            // t2 = a0-a1           in [0, 4p]
  mp_sub_p2(t2, ta[0], ta[1]);
  mp_addfast(t3, ta[0], ta[0]);           // t3 = 2*a0            in [0, 4p]
  fpmul_mont(r[0], t1, t2);             // r0 = (a0+a1)(a0-a1)  in [0, 2p]
  fpmul_mont(r[1], t3, ta[1]);           // r1 = 2a0a1           in [0, 2p]
}

// integer subtraction r = r - a - b where len(r) = len(a) = len(b) = 2*NWORDS
static void mp_dblsubfast(__m512i *r, const __m512i *a, const __m512i *b)
{
  __m512i a0  = a[0],  a1  = a[1],  a2  = a[2],  a3  = a[3];  
  __m512i a4  = a[4],  a5  = a[5],  a6  = a[6],  a7  = a[7];
  __m512i a8  = a[8],  a9  = a[9],  a10 = a[10], a11 = a[11];
  __m512i a12 = a[12], a13 = a[13], a14 = a[14], a15 = a[15];
  __m512i a16 = a[16], a17 = a[17], a18 = a[18], a19 = a[19];
  __m512i a20 = a[20], a21 = a[21], a22 = a[22], a23 = a[23];
  __m512i b0  = b[0],  b1  = b[1],  b2  = b[2],  b3  = b[3];  
  __m512i b4  = b[4],  b5  = b[5],  b6  = b[6],  b7  = b[7];
  __m512i b8  = b[8],  b9  = b[9],  b10 = b[10], b11 = b[11];
  __m512i b12 = b[12], b13 = b[13], b14 = b[14], b15 = b[15];
  __m512i b16 = b[16], b17 = b[17], b18 = b[18], b19 = b[19];
  __m512i b20 = b[20], b21 = b[21], b22 = b[22], b23 = b[23];
  __m512i r0  = r[0],  r1  = r[1],  r2  = r[2],  r3  = r[3];  
  __m512i r4  = r[4],  r5  = r[5],  r6  = r[6],  r7  = r[7];
  __m512i r8  = r[8],  r9  = r[9],  r10 = r[10], r11 = r[11];
  __m512i r12 = r[12], r13 = r[13], r14 = r[14], r15 = r[15];
  __m512i r16 = r[16], r17 = r[17], r18 = r[18], r19 = r[19];
  __m512i r20 = r[20], r21 = r[21], r22 = r[22], r23 = r[23];
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

  // carry propagation
  // will be done in reduction
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


  r[0]  = r0;  r[1]  = r1;  r[2]  = r2;  r[3]  = r3;  
  r[4]  = r4;  r[5]  = r5;  r[6]  = r6;  r[7]  = r7;  
  r[8]  = r8;  r[9]  = r9;  r[10] = r10; r[11] = r11; 
  r[12] = r12; r[13] = r13; r[14] = r14; r[15] = r15; 
  r[16] = r16; r[17] = r17; r[18] = r18; r[19] = r19;
  r[20] = r20; r[21] = r21; r[22] = r22; r[23] = r23;
}

// integer subtraction followed by addition with p*612
//  r = a - b + (p*2^612) if a-b<0; otherwise r = a - b
static void mp_subaddfast(__m512i *r, const __m512i *a, const __m512i *b)
{
  __m512i a0  = a[0],  a1  = a[1],  a2  = a[2],  a3  = a[3];  
  __m512i a4  = a[4],  a5  = a[5],  a6  = a[6],  a7  = a[7];
  __m512i a8  = a[8],  a9  = a[9],  a10 = a[10], a11 = a[11];
  __m512i a12 = a[12], a13 = a[13], a14 = a[14], a15 = a[15];
  __m512i a16 = a[16], a17 = a[17], a18 = a[18], a19 = a[19];
  __m512i a20 = a[20], a21 = a[21], a22 = a[22], a23 = a[23];
  __m512i b0  = b[0],  b1  = b[1],  b2  = b[2],  b3  = b[3];  
  __m512i b4  = b[4],  b5  = b[5],  b6  = b[6],  b7  = b[7];
  __m512i b8  = b[8],  b9  = b[9],  b10 = b[10], b11 = b[11];
  __m512i b12 = b[12], b13 = b[13], b14 = b[14], b15 = b[15];
  __m512i b16 = b[16], b17 = b[17], b18 = b[18], b19 = b[19];
  __m512i b20 = b[20], b21 = b[21], b22 = b[22], b23 = b[23];
  __m512i r0,  r1,  r2,  r3,  r4,  r5,  r6,  r7,  r8,  r9,  r10, r11;
  __m512i r12, r13, r14, r15, r16, r17, r18, r19, r20, r21, r22, r23, smask;
  const __m512i vp0  = VSET1(p610[0]),  vp1  = VSET1(p610[1]);
  const __m512i vp2  = VSET1(p610[2]),  vp3  = VSET1(p610[3]);
  const __m512i vp4  = VSET1(p610[4]),  vp5  = VSET1(p610[5]);
  const __m512i vp6  = VSET1(p610[6]),  vp7  = VSET1(p610[7]);
  const __m512i vp8  = VSET1(p610[8]),  vp9  = VSET1(p610[9]);
  const __m512i vp10 = VSET1(p610[10]), vp11 = VSET1(p610[11]);
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

  // if r is non-negative, smask = all-0 
  // if r is     negative, smask = all-1
  smask = VSRA(r23, 63);
  // r = r + ((p*2^612) & smask)
  r12 = VADD(r12, VAND(vp0,  smask)); r13 = VADD(r13, VAND(vp1,  smask)); 
  r14 = VADD(r14, VAND(vp2,  smask)); r15 = VADD(r15, VAND(vp3,  smask)); 
  r16 = VADD(r16, VAND(vp4,  smask)); r17 = VADD(r17, VAND(vp5,  smask)); 
  r18 = VADD(r18, VAND(vp6,  smask)); r19 = VADD(r19, VAND(vp7,  smask)); 
  r20 = VADD(r20, VAND(vp8,  smask)); r21 = VADD(r21, VAND(vp9,  smask)); 
  r22 = VADD(r22, VAND(vp10, smask)); r23 = VADD(r23, VAND(vp11, smask)); 

  // carry propagation
  // will be done in reduction 
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


  r[0]  = r0;  r[1]  = r1;  r[2]  = r2;  r[3]  = r3;  
  r[4]  = r4;  r[5]  = r5;  r[6]  = r6;  r[7]  = r7;  
  r[8]  = r8;  r[9]  = r9;  r[10] = r10; r[11] = r11; 
  r[12] = r12; r[13] = r13; r[14] = r14; r[15] = r15; 
  r[16] = r16; r[17] = r17; r[18] = r18; r[19] = r19;
  r[20] = r20; r[21] = r21; r[22] = r22; r[23] = r23;
}

// GF(p^2) multiplication using Montgomery arithmetic, r = a*b in GF(p^2).
void fp2mul_mont(f2elm_t r, const f2elm_t a, const f2elm_t b)
{
  felm_t t1, t2;
  dfelm_t tt1, tt2, tt3;
  f2elm_t ta, tb;

  fpcorrection_p4(ta[0], a[0]);
  fpcorrection_p4(ta[1], a[1]);
  fpcorrection_p4(tb[0], b[0]);
  fpcorrection_p4(tb[1], b[1]);

  mp_addfast(t1, ta[0], ta[1]);           // t1 = a0+a1           in [0, 4p]
  mp_addfast(t2, tb[0], tb[1]);           // t2 = b0+b1           in [0, 4p]
  mp_mul(tt1, ta[0], tb[0]);              // tt1 = a0*b0          in [0, 4p^2]
  mp_mul(tt2, ta[1], tb[1]);              // tt2 = a1*b1          in [0, 4p^2]
  mp_mul(tt3, t1, t2);                  // tt3 = (a0+a1)(b0+b1) in [0, 16p^2]
  mp_dblsubfast(tt3, tt1, tt2);         // tt3 = a0b1+a1b0      in [0, 8p^2]    (only be used here)
  mp_subaddfast(tt1, tt1, tt2);         // tt1 = a0b0-a1b1      in [0, p*2^612] (only be used here)
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
