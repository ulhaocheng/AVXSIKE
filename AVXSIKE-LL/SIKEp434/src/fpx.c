/**
 *******************************************************************************
 * @version 0.1
 * @date 2021-12-16
 * @copyright Copyright © 2021 by University of Luxembourg
 * @author Hao Cheng
 *******************************************************************************
 */

#include "fpx.h"
#include <string.h>

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

// -----------------------------------------------------------------------------
// (8x1x1)-way Fp2 arithmetic and some (8x1)-way Fp arithmetic.
// NOTE: (8x1x1)-way means each function performs 8 Fp2 operations, where each 
// Fp2 operation performs 1 Fp operation and each Fp operation uses 1 lane.  

// Montgomery multplication
// optimized version
void fpmul_mont_8x1w(vfelm_t r, const vfelm_t a, const vfelm_t b)
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
  const __m512i vp4 = VSET1(vp434p1[4]), vp5 = VSET1(vp434p1[5]);
  const __m512i vp6 = VSET1(vp434p1[6]), vp7 = VSET1(vp434p1[7]);
  const __m512i vp8 = VSET1(vp434p1[8]);
  const __m512i vbmask = VSET1(VBMASK);

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

  z1 = VADD(z1, VSRA(z0, VBRADIX)); z0 = VAND(z0, vbmask);
  z2 = VADD(z2, VSRA(z1, VBRADIX)); z1 = VAND(z1, vbmask);
  z3 = VADD(z3, VSRA(z2, VBRADIX)); z2 = VAND(z2, vbmask);
  z4 = VADD(z4, VSRA(z3, VBRADIX)); z3 = VAND(z3, vbmask);

  // 2nd loop of integer multiplication

  u = z0;
  z4 = VMACLO(z4, u, vp4); z5 = VMACLO(z5, u, vp5); z6 = VMACLO(z6, u, vp6); 
  z7 = VMACLO(z7, u, vp7); z8 = VMACLO(z8, u, vp8);
  h5 = VMACHI(h5, u, vp4); h6 = VMACHI(h6, u, vp5); h7 = VMACHI(h7, u, vp6);
  h8 = VMACHI(h8, u, vp7); h9 = VMACHI(h9, u, vp8);
  z5 = VADD(z5, VSRA(z4, VBRADIX)); z4 = VAND(z4, vbmask); 
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
  z6 = VADD(z6, VSRA(z5, VBRADIX)); z5 = VAND(z5, vbmask); 
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
  z7 = VADD(z7, VSRA(z6, VBRADIX)); z6 = VAND(z6, vbmask); 
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
  z8 = VADD(z8, VSRA(z7, VBRADIX)); z7 = VAND(z7, vbmask); 
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
  z9 = VADD(z9, VSRA(z8, VBRADIX)); z8 = VAND(z8, vbmask); 
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
  z10 = VADD(z10, VSRA(z9, VBRADIX)); z9 = VAND(z9, vbmask); 
  z10 = VADD(z10, VADD(h10, h10)); 

  z14 = VMACLO(z14, a6, b8); z14 = VMACLO(z14, a7, b7); z14 = VMACLO(z14, a8, b6); 
  h15 = VMACHI(h15, a6, b8); h15 = VMACHI(h15, a7, b7); h15 = VMACHI(h15, a8, b6); 

  u = z6; 
  z10 = VMACLO(z10, u, vp4); z11 = VMACLO(z11, u, vp5); z12 = VMACLO(z12, u, vp6); 
  z13 = VMACLO(z13, u, vp7); z14 = VMACLO(z14, u, vp8);
  h11 = VMACHI(h11, u, vp4); h12 = VMACHI(h12, u, vp5); h13 = VMACHI(h13, u, vp6);
  h14 = VMACHI(h14, u, vp7); h15 = VMACHI(h15, u, vp8);
  z11 = VADD(z11, VSRA(z10, VBRADIX)); z10 = VAND(z10, vbmask); 
  z11 = VADD(z11, VADD(h11, h11));

  z15 = VMACLO(z15, a7, b8); z15 = VMACLO(z15, a8, b7); 
  h16 = VMACHI(h16, a7, b8); h16 = VMACHI(h16, a8, b7); 

  u = z7; 
  z11 = VMACLO(z11, u, vp4); z12 = VMACLO(z12, u, vp5); z13 = VMACLO(z13, u, vp6); 
  z14 = VMACLO(z14, u, vp7); z15 = VMACLO(z15, u, vp8);
  h12 = VMACHI(h12, u, vp4); h13 = VMACHI(h13, u, vp5); h14 = VMACHI(h14, u, vp6);
  h15 = VMACHI(h15, u, vp7); h16 = VMACHI(h16, u, vp8);
  z12 = VADD(z12, VSRA(z11, VBRADIX)); z11 = VAND(z11, vbmask); 
  z12 = VADD(z12, VADD(h12, h12));

  z16 = VMACLO(z16, a8, b8);
  h17 = VMACHI(h17, a8, b8);

  u = z8; 
  z12 = VMACLO(z12, u, vp4); z13 = VMACLO(z13, u, vp5); z14 = VMACLO(z14, u, vp6); 
  z15 = VMACLO(z15, u, vp7); z16 = VMACLO(z16, u, vp8);
  h13 = VMACHI(h13, u, vp4); h14 = VMACHI(h14, u, vp5); h15 = VMACHI(h15, u, vp6);
  h16 = VMACHI(h16, u, vp7); h17 = VMACHI(h17, u, vp8);
  z13 = VADD(z13, VSRA(z12, VBRADIX)); z12 = VAND(z12, vbmask);
  z13 = VADD(z13, VADD(h13, h13)); 

  z14 = VADD(z14, VADD(h14, h14));
  z15 = VADD(z15, VADD(h15, h15));
  z16 = VADD(z16, VADD(h16, h16));
  z17 = VADD(z17, VADD(h17, h17));

  z14 = VADD(z14, VSRA(z13, VBRADIX)); z13 = VAND(z13, vbmask);
  z15 = VADD(z15, VSRA(z14, VBRADIX)); z14 = VAND(z14, vbmask);
  z16 = VADD(z16, VSRA(z15, VBRADIX)); z15 = VAND(z15, vbmask);
  z17 = VADD(z17, VSRA(z16, VBRADIX)); z16 = VAND(z16, vbmask);

  r[0] =  z9; r[1] = z10; r[2] = z11; 
  r[3] = z12; r[4] = z13; r[5] = z14; 
  r[6] = z15; r[7] = z16; r[8] = z17; 
}

// GF(p^2) addition r = a + b in GF(p^2)
void fp2add_8x1x1w(vf2elm_t r, const vf2elm_t a, const vf2elm_t b)
{
  fpadd_8x1w(r[0], a[0], b[0]);
  fpadd_8x1w(r[1], a[1], b[1]);
}

// GF(p^2) subtraction r = a - b in GF(p^2)
void fp2sub_8x1x1w(vf2elm_t r, const vf2elm_t a, const vf2elm_t b)
{
  fpsub_8x1w(r[0], a[0], b[0]);
  fpsub_8x1w(r[1], a[1], b[1]);
}

// integer addition r = a + b
void mp_add_8x1w(__m512i *r, const __m512i *a, const __m512i *b)
{
  __m512i a0 = a[0], a1 = a[1], a2 = a[2];
  __m512i a3 = a[3], a4 = a[4], a5 = a[5];
  __m512i a6 = a[6], a7 = a[7], a8 = a[8];
  __m512i b0 = b[0], b1 = b[1], b2 = b[2];
  __m512i b3 = b[3], b4 = b[4], b5 = b[5];
  __m512i b6 = b[6], b7 = b[7], b8 = b[8];
  __m512i r0, r1, r2, r3, r4, r5, r6, r7, r8;
  const __m512i vbmask = VSET1(VBMASK); 

  // r = a + b
  r0 = VADD(a0, b0); r1 = VADD(a1, b1); r2 = VADD(a2, b2); 
  r3 = VADD(a3, b3); r4 = VADD(a4, b4); r5 = VADD(a5, b5);
  r6 = VADD(a6, b6); r7 = VADD(a7, b7); r8 = VADD(a8, b8);

  // carry propagation 
  r1 = VADD(r1, VSRA(r0, VBRADIX)); r0 = VAND(r0, vbmask);
  r2 = VADD(r2, VSRA(r1, VBRADIX)); r1 = VAND(r1, vbmask);
  r3 = VADD(r3, VSRA(r2, VBRADIX)); r2 = VAND(r2, vbmask);
  r4 = VADD(r4, VSRA(r3, VBRADIX)); r3 = VAND(r3, vbmask);
  r5 = VADD(r5, VSRA(r4, VBRADIX)); r4 = VAND(r4, vbmask);
  r6 = VADD(r6, VSRA(r5, VBRADIX)); r5 = VAND(r5, vbmask);
  r7 = VADD(r7, VSRA(r6, VBRADIX)); r6 = VAND(r6, vbmask);
  r8 = VADD(r8, VSRA(r7, VBRADIX)); r7 = VAND(r7, vbmask);

  r[0] = r0; r[1] = r1; r[2] = r2;
  r[3] = r3; r[4] = r4; r[5] = r5;
  r[6] = r6; r[7] = r7; r[8] = r8;
}

// GF(p^2) addition without correction r = a + b in GF(p^2) 
void mp2_add_8x1x1w(vf2elm_t r, const vf2elm_t a, const vf2elm_t b)
{
  mp_add_8x1w(r[0], a[0], b[0]);
  mp_add_8x1w(r[1], a[1], b[1]);
}

// GF(p^2) subtraction with correction with 2p, r = a-b+2p in GF(p^2).
void mp2_sub_p2_8x1x1w(vf2elm_t r, const vf2elm_t a, const vf2elm_t b)
{
  mp_sub_p2_8x1w(r[0], a[0], b[0]);
  mp_sub_p2_8x1w(r[1], a[1], b[1]);
} 

// GF(p^2) subtraction with correction with 4p, r = a-b+4p in GF(p^2).
static void mp2_sub_p4_8x1x1w(vf2elm_t r, const vf2elm_t a, const vf2elm_t b)
{
  mp_sub_p4_8x1w(r[0], a[0], b[0]);
  mp_sub_p4_8x1w(r[1], a[1], b[1]);
} 

// GF(p^2) squaring using Montgomery arithmetic r = a^2 in GF(p^2)
void fp2sqr_mont_8x1x1w(vf2elm_t r, const vf2elm_t a)
{
  vfelm_t t1, t2, t3;

  mp_add_8x1w(t1, a[0], a[1]);               // t1 = a0+a1           in [0, 8p]
  mp_sub_p4_8x1w(t2, a[0], a[1]);            // t2 = a0-a1           in [0, 8p]
  mp_add_8x1w(t3, a[0], a[0]);               // t3 = 2*a0            in [0, 8p]
  fpmul_mont_8x1w(r[0], t1, t2);             // r0 = (a0+a1)(a0-a1)  in [0, 2p]
  fpmul_mont_8x1w(r[1], t3, a[1]);           // r1 = 2a0a1           in [0, 2p]
}

// integer subtraction r = r - a - b where len(r) = len(a) = len(b) = 2*NWORDS
static void mp_dblsubfast_8x1w(__m512i *r, const __m512i *a, const __m512i *b)
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
  const __m512i vbmask = VSET1(VBMASK);

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

// integer subtraction followed by addition with p*459
//  r = a - b + (p*2^459) if a-b<0; otherwise r = a - b
static void mp_subaddfast_8x1w(__m512i *r, const __m512i *a, const __m512i *b)
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
  const __m512i vp0 = VSET1(vp434[0]), vp1 = VSET1(vp434[1]), vp2 = VSET1(vp434[2]);
  const __m512i vp3 = VSET1(vp434[3]), vp4 = VSET1(vp434[4]), vp5 = VSET1(vp434[5]);
  const __m512i vp6 = VSET1(vp434[6]), vp7 = VSET1(vp434[7]), vp8 = VSET1(vp434[8]);
  const __m512i vbmask = VSET1(VBMASK);

  // r = a - b
  r0  = VSUB(a0,  b0);  r1  = VSUB(a1,  b1);  r2  = VSUB(a2,  b2);
  r3  = VSUB(a3,  b3);  r4  = VSUB(a4,  b4);  r5  = VSUB(a5,  b5);
  r6  = VSUB(a6,  b6);  r7  = VSUB(a7,  b7);  r8  = VSUB(a8,  b8);
  r9  = VSUB(a9,  b9);  r10 = VSUB(a10, b10); r11 = VSUB(a11, b11);
  r12 = VSUB(a12, b12); r13 = VSUB(a13, b13); r14 = VSUB(a14, b14);
  r15 = VSUB(a15, b15); r16 = VSUB(a16, b16); r17 = VSUB(a17, b17);

  // carry propagation
  r1  = VADD(r1,  VSRA(r0,  VBRADIX)); r0  = VAND(r0,  vbmask);
  r2  = VADD(r2,  VSRA(r1,  VBRADIX)); r1  = VAND(r1,  vbmask);
  r3  = VADD(r3,  VSRA(r2,  VBRADIX)); r2  = VAND(r2,  vbmask);
  r4  = VADD(r4,  VSRA(r3,  VBRADIX)); r3  = VAND(r3,  vbmask);
  r5  = VADD(r5,  VSRA(r4,  VBRADIX)); r4  = VAND(r4,  vbmask);
  r6  = VADD(r6,  VSRA(r5,  VBRADIX)); r5  = VAND(r5,  vbmask);
  r7  = VADD(r7,  VSRA(r6,  VBRADIX)); r6  = VAND(r6,  vbmask);
  r8  = VADD(r8,  VSRA(r7,  VBRADIX)); r7  = VAND(r7,  vbmask);
  r9  = VADD(r9,  VSRA(r8,  VBRADIX)); r8  = VAND(r8,  vbmask);
  r10 = VADD(r10, VSRA(r9,  VBRADIX)); r9  = VAND(r9,  vbmask);
  r11 = VADD(r11, VSRA(r10, VBRADIX)); r10 = VAND(r10, vbmask);
  r12 = VADD(r12, VSRA(r11, VBRADIX)); r11 = VAND(r11, vbmask);
  r13 = VADD(r13, VSRA(r12, VBRADIX)); r12 = VAND(r12, vbmask);
  r14 = VADD(r14, VSRA(r13, VBRADIX)); r13 = VAND(r13, vbmask);
  r15 = VADD(r15, VSRA(r14, VBRADIX)); r14 = VAND(r14, vbmask);
  r16 = VADD(r16, VSRA(r15, VBRADIX)); r15 = VAND(r15, vbmask);
  r17 = VADD(r17, VSRA(r16, VBRADIX)); r16 = VAND(r16, vbmask);

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
void fp2mul_mont_8x1x1w(vf2elm_t r, const vf2elm_t a, const vf2elm_t b)
{
  vfelm_t t1, t2;
  vdfelm_t tt1, tt2, tt3;

  mp_add_8x1w(t1, a[0], a[1]);          // t1 = a0+a1           in [0, 4p]
  mp_add_8x1w(t2, b[0], b[1]);          // t2 = b0+b1           in [0, 4p]
  mp_mul_8x1w(tt1, a[0], b[0]);         // tt1 = a0*b0          in [0, 4p^2]
  mp_mul_8x1w(tt2, a[1], b[1]);         // tt2 = a1*b1          in [0, 4p^2]
  mp_mul_8x1w(tt3, t1, t2);             // tt3 = (a0+a1)(b0+b1) in [0, 16p^2]
  mp_dblsubfast_8x1w(tt3, tt1, tt2);    // tt3 = a0b1+a1b0      in [0, 8p^2]    (only be used here)
  mp_subaddfast_8x1w(tt1, tt1, tt2);    // tt1 = a0b0-a1b1      in [0, p*2^510] (only be used here)
  rdc_mont_8x1w(r[1], tt3);             // r1 = a0b1+a1b0 mod 2p                (only be used here)
  rdc_mont_8x1w(r[0], tt1);             // r0 = a0b0-a1b1 mod 2p                (only be used here)
 }

// modular division by two r = a/2 mod p
void fpdiv2_8x1w(__m512i *r, const __m512i *a)
{
  __m512i a0 = a[0], a1 = a[1], a2 = a[2];
  __m512i a3 = a[3], a4 = a[4], a5 = a[5];
  __m512i a6 = a[6], a7 = a[7], a8 = a[8];
  __m512i z0, z1, z2, z3, z4, z5, z6, z7, z8;
  __m512i r0, r1, r2, r3, r4, r5, r6, r7, r8, smask;
  const __m512i vp0 = VSET1(vp434[0]), vp1 = VSET1(vp434[1]), vp2 = VSET1(vp434[2]);
  const __m512i vp3 = VSET1(vp434[3]), vp4 = VSET1(vp434[4]), vp5 = VSET1(vp434[5]);
  const __m512i vp6 = VSET1(vp434[6]), vp7 = VSET1(vp434[7]), vp8 = VSET1(vp434[8]);
  const __m512i vbmask = VSET1(VBMASK), vone = VSET1(1); 

  // if a is even, smask is all-0
  // if a is  odd, smask is all-1
  smask = VSUB(VZERO, VAND(a0, vone));

  // z = a + (p & smask)
  z0 = VADD(a0, VAND(vp0, smask)); z1 = VADD(a1, VAND(vp1, smask)); 
  z2 = VADD(a2, VAND(vp2, smask)); z3 = VADD(a3, VAND(vp3, smask)); 
  z4 = VADD(a4, VAND(vp4, smask)); z5 = VADD(a5, VAND(vp5, smask)); 
  z6 = VADD(a6, VAND(vp6, smask)); z7 = VADD(a7, VAND(vp7, smask)); 
  z8 = VADD(a8, VAND(vp8, smask)); 

  // carry pzopagation 
  z1 = VADD(z1, VSRA(z0, VBRADIX)); z0 = VAND(z0, vbmask);
  z2 = VADD(z2, VSRA(z1, VBRADIX)); z1 = VAND(z1, vbmask);
  z3 = VADD(z3, VSRA(z2, VBRADIX)); z2 = VAND(z2, vbmask);
  z4 = VADD(z4, VSRA(z3, VBRADIX)); z3 = VAND(z3, vbmask);
  z5 = VADD(z5, VSRA(z4, VBRADIX)); z4 = VAND(z4, vbmask);
  z6 = VADD(z6, VSRA(z5, VBRADIX)); z5 = VAND(z5, vbmask);
  z7 = VADD(z7, VSRA(z6, VBRADIX)); z6 = VAND(z6, vbmask);
  z8 = VADD(z8, VSRA(z7, VBRADIX)); z7 = VAND(z7, vbmask);

  // shift 1 bit to the right (to be optimized)
  r0 = VOR(VSHL(VAND(z1, vone), VBRADIX-1), VSRA(z0, 1));
  r1 = VOR(VSHL(VAND(z2, vone), VBRADIX-1), VSRA(z1, 1));
  r2 = VOR(VSHL(VAND(z3, vone), VBRADIX-1), VSRA(z2, 1));
  r3 = VOR(VSHL(VAND(z4, vone), VBRADIX-1), VSRA(z3, 1));
  r4 = VOR(VSHL(VAND(z5, vone), VBRADIX-1), VSRA(z4, 1));
  r5 = VOR(VSHL(VAND(z6, vone), VBRADIX-1), VSRA(z5, 1));
  r6 = VOR(VSHL(VAND(z7, vone), VBRADIX-1), VSRA(z6, 1));
  r7 = VOR(VSHL(VAND(z8, vone), VBRADIX-1), VSRA(z7, 1));
  r8 = VSRA(z8, 1);

  r[0] = r0; r[1] = r1; r[2] = r2;
  r[3] = r3; r[4] = r4; r[5] = r5;
  r[6] = r6; r[7] = r7; r[8] = r8;
}

// GF(p^2) division by two r = a/2  in GF(p^2)
void fp2div2_8x1x1w(vf2elm_t r, const vf2elm_t a)
{
  fpdiv2_8x1w(r[0], a[0]);
  fpdiv2_8x1w(r[1], a[1]);
}

// -----------------------------------------------------------------------------
// (4x2x1)-way Fp2 arithmetic (based on (8x1)-way Fp arithmetic).
// NOTE: (4x2x1)-way means each function performs 4 Fp2 operations, where each 
// Fp2 operation performs 2 Fp operations and each Fp operation uses 1 lane.  

// mixed integer addition and subtration r = <a1+b1 | a0-b0>
void static mp_mixaddsub_8x1w(vdfelm_t r, const vdfelm_t a, const vdfelm_t b)
{
  __m512i a0  = a[0],  a1  = a[1],  a2  = a[2],  a3  = a[3],  a4  = a[4];
  __m512i a5  = a[5],  a6  = a[6],  a7  = a[7],  a8  = a[8],  a9  = a[9];
  __m512i a10 = a[10], a11 = a[11], a12 = a[12], a13 = a[13], a14 = a[14];
  __m512i a15 = a[15], a16 = a[16], a17 = a[17];
  __m512i b0  = b[0],  b1  = b[1],  b2  = b[2],  b3  = b[3],  b4  = b[4];
  __m512i b5  = b[5],  b6  = b[6],  b7  = b[7],  b8  = b[8],  b9  = b[9];
  __m512i b10 = b[10], b11 = b[11], b12 = b[12], b13 = b[13], b14 = b[14];
  __m512i b15 = b[15], b16 = b[16], b17 = b[17];
  __m512i r0,  r1,  r2,  r3,  r4,  r5,  r6,  r7,  r8,  r9;
  __m512i r10, r11, r12, r13, r14, r15, r16, r17, smask;
  __m512i t0,  t1,  t2,  t3,  t4,  t5,  t6,  t7,  t8,  t9;
  __m512i t10, t11, t12, t13, t14, t15, t16, t17;
  const __m512i vp0 = VSET1(vp434[0]), vp1 = VSET1(vp434[1]), vp2 = VSET1(vp434[2]);
  const __m512i vp3 = VSET1(vp434[3]), vp4 = VSET1(vp434[4]), vp5 = VSET1(vp434[5]);
  const __m512i vp6 = VSET1(vp434[6]), vp7 = VSET1(vp434[7]), vp8 = VSET1(vp434[8]);
  const __m512i vbmask = VSET1(VBMASK), zero = VZERO; 

  // r = a1 - 0 | a0 - b0
  r0 = VMSUB(a0, 0x55, a0, b0); r1 = VMSUB(a1, 0x55, a1, b1);
  r2 = VMSUB(a2, 0x55, a2, b2); r3 = VMSUB(a3, 0x55, a3, b3);
  r4 = VMSUB(a4, 0x55, a4, b4); r5 = VMSUB(a5, 0x55, a5, b5);
  r6 = VMSUB(a6, 0x55, a6, b6); r7 = VMSUB(a7, 0x55, a7, b7);
  r8 = VMSUB(a8, 0x55, a8, b8); r9 = VMSUB(a9, 0x55, a9, b9);
  r10 = VMSUB(a10, 0x55, a10, b10); r11 = VMSUB(a11, 0x55, a11, b11);
  r12 = VMSUB(a12, 0x55, a12, b12); r13 = VMSUB(a13, 0x55, a13, b13);
  r14 = VMSUB(a14, 0x55, a14, b14); r15 = VMSUB(a15, 0x55, a15, b15);
  r16 = VMSUB(a16, 0x55, a16, b16); r17 = VMSUB(a17, 0x55, a17, b17);

  // carry propagation
  r1  = VADD(r1,  VSRA(r0,  VBRADIX)); r0  = VAND(r0,  vbmask);
  r2  = VADD(r2,  VSRA(r1,  VBRADIX)); r1  = VAND(r1,  vbmask);
  r3  = VADD(r3,  VSRA(r2,  VBRADIX)); r2  = VAND(r2,  vbmask);
  r4  = VADD(r4,  VSRA(r3,  VBRADIX)); r3  = VAND(r3,  vbmask);
  r5  = VADD(r5,  VSRA(r4,  VBRADIX)); r4  = VAND(r4,  vbmask);
  r6  = VADD(r6,  VSRA(r5,  VBRADIX)); r5  = VAND(r5,  vbmask);
  r7  = VADD(r7,  VSRA(r6,  VBRADIX)); r6  = VAND(r6,  vbmask);
  r8  = VADD(r8,  VSRA(r7,  VBRADIX)); r7  = VAND(r7,  vbmask);
  r9  = VADD(r9,  VSRA(r8,  VBRADIX)); r8  = VAND(r8,  vbmask);
  r10 = VADD(r10, VSRA(r9,  VBRADIX)); r9  = VAND(r9,  vbmask);
  r11 = VADD(r11, VSRA(r10, VBRADIX)); r10 = VAND(r10, vbmask);
  r12 = VADD(r12, VSRA(r11, VBRADIX)); r11 = VAND(r11, vbmask);
  r13 = VADD(r13, VSRA(r12, VBRADIX)); r12 = VAND(r12, vbmask);
  r14 = VADD(r14, VSRA(r13, VBRADIX)); r13 = VAND(r13, vbmask);
  r15 = VADD(r15, VSRA(r14, VBRADIX)); r14 = VAND(r14, vbmask);
  r16 = VADD(r16, VSRA(r15, VBRADIX)); r15 = VAND(r15, vbmask);
  r17 = VADD(r17, VSRA(r16, VBRADIX)); r16 = VAND(r16, vbmask);

  // if r is non-negative, smask = all-0 
  // if r is     negative, smask = all-1
  smask = VSRA(r17, 63);
  // t = b1 | (p*2^448)&smask
  t0 = VMBLEND(0x55, b0, zero); t1 = VMBLEND(0x55, b1, zero); 
  t2 = VMBLEND(0x55, b2, zero); t3 = VMBLEND(0x55, b3, zero); 
  t4 = VMBLEND(0x55, b4, zero); t5 = VMBLEND(0x55, b5, zero); 
  t6 = VMBLEND(0x55, b6, zero); t7 = VMBLEND(0x55, b7, zero); 
  t8 = VMBLEND(0x55, b8, zero); 
  t9  = VMBLEND(0x55, b9,  VAND(vp0, smask)); 
  t10 = VMBLEND(0x55, b10, VAND(vp1, smask)); 
  t11 = VMBLEND(0x55, b11, VAND(vp2, smask));
  t12 = VMBLEND(0x55, b12, VAND(vp3, smask));
  t13 = VMBLEND(0x55, b13, VAND(vp4, smask));
  t14 = VMBLEND(0x55, b14, VAND(vp5, smask));
  t15 = VMBLEND(0x55, b15, VAND(vp6, smask));
  t16 = VMBLEND(0x55, b16, VAND(vp7, smask));
  t17 = VMBLEND(0x55, b17, VAND(vp8, smask));

  // r = r + t = a1 + b1 | a0-b0 + ((p*2^448)&smask)
  r0  = VADD(r0,  t0);  r1  = VADD(r1,  t1);  r2  = VADD(r2,  t2);
  r3  = VADD(r3,  t3);  r4  = VADD(r4,  t4);  r5  = VADD(r5,  t5);
  r6  = VADD(r6,  t6);  r7  = VADD(r7,  t7);  r8  = VADD(r8,  t8);
  r9  = VADD(r9,  t9);  r10 = VADD(r10, t10); r11 = VADD(r11, t11);
  r12 = VADD(r12, t12); r13 = VADD(r13, t13); r14 = VADD(r14, t14);
  r15 = VADD(r15, t15); r16 = VADD(r16, t16); r17 = VADD(r17, t17); 

  // carry propagation
  // will be done in reduction 

  r[0]  = r0;  r[1]  = r1;  r[2]  = r2;  r[3]  = r3;  r[4]  = r4;  
  r[5]  = r5;  r[6]  = r6;  r[7]  = r7;  r[8]  = r8;  r[9]  = r9;  
  r[10] = r10; r[11] = r11; r[12] = r12; r[13] = r13; r[14] = r14;
  r[15] = r15; r[16] = r16; r[17] = r17; 
}

// GF(p^2) multiplication using Montgomery arithmetic, r = a*b in GF(p^2).
// Here uses schoolbook instead of Karatsuba.
// Assume inputs a and b are in the format of <a1 | a0> and <b1 | b0>. 
void fp2mul_mont_4x2x1w(vfelm_t r, const vfelm_t a, const vfelm_t b)
{
  vfelm_t t1, t2, t3;
  vdfelm_t tt1, tt2;

  vec_shufll_8x1w(t1, a);               // t1 = a0 | a0
  vec_shufhh_8x1w(t2, a);               // t2 = a1 | a1 
  vec_shuflh_8x1w(t3, b);               // t3 = b0 | b1
  mp_mul_8x1w(tt1, t1, b);              // tt1 = a0*b1 | a0*b0
  mp_mul_8x1w(tt2, t2, t3);             // tt2 = a1*b0 | a1*b1
  mp_mixaddsub_8x1w(tt1, tt1, tt2);     // tt3 = a0b1 + a1b0 | a0b0 - a1b1
  rdc_mont_8x1w(r, tt1);                // r = r1 | r0 = (a0b1+a1b0) mod 2p | (a0b0-a1b1) mod 2p
}

// GF(p^2) squaring using Montgomery arithmetic r = a^2 in GF(p^2)
void fp2sqr_mont_4x2x1w(vfelm_t r, const vfelm_t a)
{
  vfelm_t t1, t2;

  vec_shufll_8x1w(t1, a);               // t1 = a0 | a0
  vec_shuflh_8x1w(t2, a);               // t2 = a0 | a1
  mp_add_8x1w(t1, t1, t2);              // t1 = 2*a0 | a0+a1
  vec_shufzh_8x1w(t2, a);               // t2 = 0 | a1
  mp_sub_p4_8x1w(t2, a, t2);            // t2 = a1 | a0-a1 
  fpmul_mont_8x1w(r, t1, t2);           // r = r1 | r0 = 2a0*a1 | (a0+a1)*(a0-a1) 
}

// -----------------------------------------------------------------------------
// (2x4x1)-way Fp2 arithmetic (based on (8x1)-way Fp arithmetic).
// NOTE: (2x4x1)-way means each function performs 2 Fp2 operations, where each 
// Fp2 operation performs 4 Fp operations and each Fp operation uses 1 lane.  

// a = a3 | a2 | a1 | a0 -> a0 | a2 | a2 | a0
static void vec_perm0220_8x1w(vfelm_t r, const vfelm_t a)
{
  __m512i a0 = a[0], a1 = a[1], a2 = a[2], a3 = a[3], a4 = a[4];
  __m512i a5 = a[5], a6 = a[6], a7 = a[7], a8 = a[8];
  __m512i r0, r1, r2, r3, r4, r5, r6, r7, r8;

  r0 = VPERM(a0, 0x28); r1 = VPERM(a1, 0x28);
  r2 = VPERM(a2, 0x28); r3 = VPERM(a3, 0x28);
  r4 = VPERM(a4, 0x28); r5 = VPERM(a5, 0x28);
  r6 = VPERM(a6, 0x28); r7 = VPERM(a7, 0x28);
  r8 = VPERM(a8, 0x28);

  r[0] = r0; r[1] = r1; r[2] = r2; r[3] = r3; r[4] = r4; 
  r[5] = r5; r[6] = r6; r[7] = r7; r[8] = r8; 
}

// a = a3 | a2 | a1 | a0 -> 0 | a3 | 0 | a2
static void vec_permz3z2_8x1w(vfelm_t r, const vfelm_t a)
{
  __m512i a0 = a[0], a1 = a[1], a2 = a[2], a3 = a[3], a4 = a[4];
  __m512i a5 = a[5], a6 = a[6], a7 = a[7], a8 = a[8];
  __m512i r0, r1, r2, r3, r4, r5, r6, r7, r8;

  r0 = VZPERM(0x55, a0, 0xFA); r1 = VZPERM(0x55, a1, 0xFA);
  r2 = VZPERM(0x55, a2, 0xFA); r3 = VZPERM(0x55, a3, 0xFA);
  r4 = VZPERM(0x55, a4, 0xFA); r5 = VZPERM(0x55, a5, 0xFA);
  r6 = VZPERM(0x55, a6, 0xFA); r7 = VZPERM(0x55, a7, 0xFA);
  r8 = VZPERM(0x55, a8, 0xFA); 

  r[0] = r0; r[1] = r1; r[2] = r2; r[3] = r3; r[4] = r4; 
  r[5] = r5; r[6] = r6; r[7] = r7; r[8] = r8; 
}

// a = a3 | a2 | a1 | a0 -> 0 | a1 | 0 | a0
static void vec_permz1z0_8x1w(vfelm_t r, const vfelm_t a)
{
  __m512i a0 = a[0], a1 = a[1], a2 = a[2], a3 = a[3], a4 = a[4];
  __m512i a5 = a[5], a6 = a[6], a7 = a[7], a8 = a[8];
  __m512i r0, r1, r2, r3, r4, r5, r6, r7, r8;

  r0 = VZPERM(0x55, a0, 0x50); r1 = VZPERM(0x55, a1, 0x50);
  r2 = VZPERM(0x55, a2, 0x50); r3 = VZPERM(0x55, a3, 0x50);
  r4 = VZPERM(0x55, a4, 0x50); r5 = VZPERM(0x55, a5, 0x50);
  r6 = VZPERM(0x55, a6, 0x50); r7 = VZPERM(0x55, a7, 0x50);
  r8 = VZPERM(0x55, a8, 0x50); 

  r[0] = r0; r[1] = r1; r[2] = r2; r[3] = r3; r[4] = r4; 
  r[5] = r5; r[6] = r6; r[7] = r7; r[8] = r8; 
}

// mixed field addition and subtration r = < 0 | a2+b2 | 0 | a0-b0>
// i.e. r = <0+0-0 | a2+b2-2p | 0+0-0 | a0+0-b0>
static void fp_mixaddsub_8x1w(vfelm_t r, const vfelm_t a, const vfelm_t b)
{
  __m512i a0 = a[0], a1 = a[1], a2 = a[2], a3 = a[3], a4 = a[4];
  __m512i a5 = a[5], a6 = a[6], a7 = a[7], a8 = a[8];
  __m512i b0 = b[0], b1 = b[1], b2 = b[2], b3 = b[3], b4 = b[4];
  __m512i b5 = b[5], b6 = b[6], b7 = b[7], b8 = b[8];
  __m512i r0, r1, r2, r3, r4, r5, r6, r7, r8, smask;
  const __m512i vp0 = VSET(0, vp434x2[0], 0, vp434x2[0], 0, vp434x2[0], 0, vp434x2[0]); 
  const __m512i vp1 = VSET(0, vp434x2[1], 0, vp434x2[1], 0, vp434x2[1], 0, vp434x2[1]); 
  const __m512i vp2 = VSET(0, vp434x2[2], 0, vp434x2[2], 0, vp434x2[2], 0, vp434x2[2]); 
  const __m512i vp3 = VSET(0, vp434x2[3], 0, vp434x2[3], 0, vp434x2[3], 0, vp434x2[3]);
  const __m512i vp4 = VSET(0, vp434x2[4], 0, vp434x2[4], 0, vp434x2[4], 0, vp434x2[4]); 
  const __m512i vp5 = VSET(0, vp434x2[5], 0, vp434x2[5], 0, vp434x2[5], 0, vp434x2[5]);
  const __m512i vp6 = VSET(0, vp434x2[6], 0, vp434x2[6], 0, vp434x2[6], 0, vp434x2[6]); 
  const __m512i vp7 = VSET(0, vp434x2[7], 0, vp434x2[7], 0, vp434x2[7], 0, vp434x2[7]);
  const __m512i vp8 = VSET(0, vp434x2[8], 0, vp434x2[8], 0, vp434x2[8], 0, vp434x2[8]); 
  const __m512i vbmask = VSET1(VBMASK); 
  __m512i t0, t1, t2, t3, t4, t5, t6, t7, t8;

  // r = 0+0 | a2+b2 | 0+0 | a0+0
  r0 = VMADD(a0, 0x44, a0, b0); r1 = VMADD(a1, 0x44, a1, b1);
  r2 = VMADD(a2, 0x44, a2, b2); r3 = VMADD(a3, 0x44, a3, b3);
  r4 = VMADD(a4, 0x44, a4, b4); r5 = VMADD(a5, 0x44, a5, b5);
  r6 = VMADD(a6, 0x44, a6, b6); r7 = VMADD(a7, 0x44, a7, b7);
  r8 = VMADD(a8, 0x44, a8, b8); 

  // t = 0 | 2p | 0 | b0
  t0 = VMBLEND(0x33, vp0, b0); t1 = VMBLEND(0x33, vp1, b1); 
  t2 = VMBLEND(0x33, vp2, b2); t3 = VMBLEND(0x33, vp3, b3); 
  t4 = VMBLEND(0x33, vp4, b4); t5 = VMBLEND(0x33, vp5, b5); 
  t6 = VMBLEND(0x33, vp6, b6); t7 = VMBLEND(0x33, vp7, b7); 
  t8 = VMBLEND(0x33, vp8, b8);  

  // r = r - t
  r0 = VSUB(r0, t0); r1 = VSUB(r1, t1);
  r2 = VSUB(r2, t2); r3 = VSUB(r3, t3);
  r4 = VSUB(r4, t4); r5 = VSUB(r5, t5);
  r6 = VSUB(r6, t6); r7 = VSUB(r7, t7);
  r8 = VSUB(r8, t8); 

  // if r is non-negative, smask = all-0 
  // if r is     negative, smask = all-1
  smask = VSRA(r8, 63);
  // r = r + (2p & smask)
  r0 = VADD(r0, VAND(vp0, smask)); r1 = VADD(r1, VAND(vp1, smask)); 
  r2 = VADD(r2, VAND(vp2, smask)); r3 = VADD(r3, VAND(vp3, smask)); 
  r4 = VADD(r4, VAND(vp4, smask)); r5 = VADD(r5, VAND(vp5, smask)); 
  r6 = VADD(r6, VAND(vp6, smask)); r7 = VADD(r7, VAND(vp7, smask)); 
  r8 = VADD(r8, VAND(vp8, smask)); 

  // carry propagation 
  r1 = VADD(r1, VSRA(r0, VBRADIX)); r0 = VAND(r0, vbmask);
  r2 = VADD(r2, VSRA(r1, VBRADIX)); r1 = VAND(r1, vbmask);
  r3 = VADD(r3, VSRA(r2, VBRADIX)); r2 = VAND(r2, vbmask);
  r4 = VADD(r4, VSRA(r3, VBRADIX)); r3 = VAND(r3, vbmask);
  r5 = VADD(r5, VSRA(r4, VBRADIX)); r4 = VAND(r4, vbmask);
  r6 = VADD(r6, VSRA(r5, VBRADIX)); r5 = VAND(r5, vbmask);
  r7 = VADD(r7, VSRA(r6, VBRADIX)); r6 = VAND(r6, vbmask);
  r8 = VADD(r8, VSRA(r7, VBRADIX)); r7 = VAND(r7, vbmask);

  r[0] = r0; r[1] = r1; r[2] = r2; r[3] = r3; r[4] = r4; 
  r[5] = r5; r[6] = r6; r[7] = r7; r[8] = r8; 

}

// GF(p^2) multiplication using Montgomery arithmetic, r = a*b in GF(p^2).
// Here uses schoolbook instead of Karatsuba.
// Assume inputs a and b are in the format of <0 | a1 | 0 | a0> and <0 | b1 | 0 | b0>. 
void fp2mul_mont_2x4x1w(vfelm_t r, const vfelm_t a, const vfelm_t b)
{
  vfelm_t t1, t2;

  vec_shufll_8x1w(t1, a);               // t1 = a1 | a1 | a0 | a0
  vec_perm0220_8x1w(t2, b);             // t2 = b0 | b1 | b1 | b0
  fpmul_mont_8x1w(t1, t1, t2);          // t1 = a1b0 | a1b1 | a0b1 | a0b0
  vec_permz3z2_8x1w(t2, t1);            // t2 = 0 | a1b0 | 0 | a1b1
  vec_permz1z0_8x1w(t1, t1);            // t1 = 0 | a0b1 | 0 | a0b0
  fp_mixaddsub_8x1w(r, t1, t2);         // r  = 0 | a0b1+a1b0 | 0 | a0b0-a1b1
}

// -----------------------------------------------------------------------------
// (4x1x2)-way Fp2 arithmetic (based on (4x2)-way Fp arithmetic).
// NOTE: (4x2x1)-way means each function performs 4 Fp2 operations, where each 
// Fp2 operation performs 1 Fp operation and each Fp operation uses 2 lanes.  

// a = a3 | a2 | a1 | a0 -> a1 | a0 | a1 | a0
static void vec_perm1010_4x2w(vgelm_t r, const vgelm_t a)
{
  __m512i a0 = a[0], a1 = a[1], a2 = a[2], a3 = a[3], a4 = a[4];
  __m512i r0, r1, r2, r3, r4;

  r0 = VPERM(a0, 0x44); r1 = VPERM(a1, 0x44);
  r2 = VPERM(a2, 0x44); r3 = VPERM(a3, 0x44);
  r4 = VPERM(a4, 0x44); 

  r[0] = r0; r[1] = r1; r[2] = r2; r[3] = r3; r[4] = r4; 
}

// a = a3 | a2 | a1 | a0 -> a3 | a2 | a3 | a2
static void vec_perm3232_4x2w(vgelm_t r, const vgelm_t a)
{
  __m512i a0 = a[0], a1 = a[1], a2 = a[2], a3 = a[3], a4 = a[4];
  __m512i r0, r1, r2, r3, r4;

  r0 = VPERM(a0, 0xEE); r1 = VPERM(a1, 0xEE);
  r2 = VPERM(a2, 0xEE); r3 = VPERM(a3, 0xEE);
  r4 = VPERM(a4, 0xEE); 

  r[0] = r0; r[1] = r1; r[2] = r2; r[3] = r3; r[4] = r4; 
}

// a = a3 | a2 | a1 | a0 -> a1 | a0 | a3 | a2
static void vec_perm1032_4x2w(vgelm_t r, const vgelm_t a)
{
  __m512i a0 = a[0], a1 = a[1], a2 = a[2], a3 = a[3], a4 = a[4];
  __m512i r0, r1, r2, r3, r4;

  r0 = VPERM(a0, 0x4E); r1 = VPERM(a1, 0x4E);
  r2 = VPERM(a2, 0x4E); r3 = VPERM(a3, 0x4E);
  r4 = VPERM(a4, 0x4E); 

  r[0] = r0; r[1] = r1; r[2] = r2; r[3] = r3; r[4] = r4; 
}

// a = a3 | a2 | a1 | a0 -> 0 | 0 | a3 | a2
static void vec_permzz32_4x2w(vgelm_t r, const vgelm_t a)
{
  __m512i a0 = a[0], a1 = a[1], a2 = a[2], a3 = a[3], a4 = a[4];
  __m512i r0, r1, r2, r3, r4;

  r0 = VZPERM(0x33, a0, 0xEE); r1 = VZPERM(0x33, a1, 0xEE);
  r2 = VZPERM(0x33, a2, 0xEE); r3 = VZPERM(0x33, a3, 0xEE);
  r4 = VZPERM(0x33, a4, 0xEE); 

  r[0] = r0; r[1] = r1; r[2] = r2; r[3] = r3; r[4] = r4; 
}

// mixed integer addition and subtration r = <a1+b1 | a0-b0>
static void mp_mixaddsub_4x2w(vdgelm_t r, const vdgelm_t a, const vdgelm_t b)
{
  __m512i a0  = a[0],  a1  = a[1],  a2  = a[2],  a3  = a[3],  a4  = a[4];
  __m512i a5  = a[5],  a6  = a[6],  a7  = a[7],  a8  = a[8],  a9  = a[9];
  __m512i a10 = a[10], a11 = a[11], a12 = a[12], a13 = a[13];
  __m512i b0  = b[0],  b1  = b[1],  b2  = b[2],  b3  = b[3],  b4  = b[4];
  __m512i b5  = b[5],  b6  = b[6],  b7  = b[7],  b8  = b[8],  b9  = b[9];
  __m512i b10 = b[10], b11 = b[11], b12 = b[12], b13 = b[13];
  __m512i r0,  r1,  r2,  r3,  r4,  r5,  r6,  r7,  r8,  r9;
  __m512i r10, r11, r12, r13, smask, c;
  __m512i t0,  t1,  t2,  t3,  t4,  t5,  t6,  t7,  t8,  t9;
  __m512i t10, t11, t12, t13;
  const __m512i vp0 = VSET(vp434[5], vp434[0], vp434[5], vp434[0], vp434[5], vp434[0], vp434[5], vp434[0]);
  const __m512i vp1 = VSET(vp434[6], vp434[1], vp434[6], vp434[1], vp434[6], vp434[1], vp434[6], vp434[1]);
  const __m512i vp2 = VSET(vp434[7], vp434[2], vp434[7], vp434[2], vp434[7], vp434[2], vp434[7], vp434[2]);
  const __m512i vp3 = VSET(vp434[8], vp434[3], vp434[8], vp434[3], vp434[8], vp434[3], vp434[8], vp434[3]);
  const __m512i vp4 = VSET(         0, vp434[4],          0, vp434[4],          0, vp434[4],          0, vp434[4]);
  const __m512i vbmask = VSET1(VBMASK), lbmask = VSET(0, VBMASK, 0, VBMASK, 0, VBMASK, 0, VBMASK), zero = VZERO; 

  // r = a1 - 0 | a0 - b0
  r0  = VMSUB(a0,  0x33, a0,  b0);  r1  = VMSUB(a1,  0x33, a1,  b1);
  r2  = VMSUB(a2,  0x33, a2,  b2);  r3  = VMSUB(a3,  0x33, a3,  b3);
  r4  = VMSUB(a4,  0x33, a4,  b4);  r5  = VMSUB(a5,  0x33, a5,  b5);
  r6  = VMSUB(a6,  0x33, a6,  b6);  r7  = VMSUB(a7,  0x33, a7,  b7);
  r8  = VMSUB(a8,  0x33, a8,  b8);  r9  = VMSUB(a9,  0x33, a9,  b9);
  r10 = VMSUB(a10, 0x33, a10, b10); r11 = VMSUB(a11, 0x33, a11, b11);
  r12 = VMSUB(a12, 0x33, a12, b12); r13 = VMSUB(a13, 0x33, a13, b13);

  // *special* carry propagation
  c = VSRA(r0, VBRADIX);                      // c = r5'>>VBRADIX | r0>>VBRADIX 
  r1 = VMADD(r1, 0x55, r1, c);                // r1 = r6' | r1 + r0>>VBRADIX
  r5 = VMADD(r5, 0x55, r5, VSHUF(r0, 0x4E));  // r5 = r10'| r5 + r5' 
  r0 = VAND(r0, lbmask);                      // r0 =  0  | r0 & VBMASK
  c = VSRA(r1, VBRADIX);                       
  r2 = VMADD(r2, 0x55, r2, c);                
  r6 = VMADD(r6, 0x55, r6, VSHUF(r1, 0x4E));  
  r1 = VAND(r1, lbmask);   
  c = VSRA(r2, VBRADIX);                       
  r3 = VMADD(r3, 0x55, r3, c);                
  r7 = VMADD(r7, 0x55, r7, VSHUF(r2, 0x4E));  
  r2 = VAND(r2, lbmask);                      
  c = VSRA(r3, VBRADIX);                       
  r4 = VMADD(r4, 0x55, r4, c);                
  r8 = VMADD(r8, 0x55, r8, VSHUF(r3, 0x4E));  
  r3 = VAND(r3, lbmask); 
  c = VSRA(r4, VBRADIX);                       
  r5 = VMADD(r5, 0x55, r5, c);                
  r9 = VMADD(r9, 0x55, r9, VSHUF(r4, 0x4E));  
  r4 = VAND(r4, lbmask); 
  c = VSRA(r5, VBRADIX);                       
  r6 = VMADD(r6, 0x55, r6, c);                
  r10 = VMADD(r10, 0x55, r10, VSHUF(r5, 0x4E));  
  r5 = VAND(r5, lbmask); 
  c = VSRA(r6, VBRADIX);                       
  r7 = VMADD(r7, 0x55, r7, c);                
  r11 = VMADD(r11, 0x55, r11, VSHUF(r6, 0x4E));  
  r6 = VAND(r6, lbmask); 
  c = VSRA(r7, VBRADIX);                       
  r8 = VMADD(r8, 0x55, r8, c);                
  r12 = VMADD(r12, 0x55, r12, VSHUF(r7, 0x4E));  
  r7 = VAND(r7, lbmask); 
  c = VSRA(r8, VBRADIX);                       
  r9 = VMADD(r9, 0x55, r9, c);                
  r13 = VMADD(r13, 0x55, r13, VSHUF(r8, 0x4E));  
  r8 = VAND(r8, lbmask); 

  c = VSRA(r9,  VBRADIX);                     // c = r14>>VBRADIX | r9>>VBRADIX
  r10 = VMADD(r10, 0x55, r10, c);             // r10 = r15 | r10 + r9>>BRAIDX
  r9  = VMAND(r9,  0x55, r9,  vbmask);        // r9  = r14 | r9 & VBMASK  
  c = VSRA(r10, VBRADIX);                     // c = r15>>VBRADIX | r10>>VBRADIX
  r11 = VMADD(r11, 0x55, r11, c);             // r11 = r16 | r11 + r10>>BRAIDX
  r10 = VMAND(r10, 0x55, r10, vbmask);        // r10 = r15 | r10 & VBMASK   
  c = VSRA(r11, VBRADIX);                      
  r12 = VMADD(r12, 0x55, r12, c);             
  r11 = VMAND(r11, 0x55, r11, vbmask); 
  c = VSRA(r12, VBRADIX);                      
  r13 = VMADD(r13, 0x55, r13, c);             
  r12 = VMAND(r12, 0x55, r12, vbmask);
  c = VSRA(r13, VBRADIX);   
  r9 = VMADD(r9, 0xAA, r9, VSHUF(c, 0x4E));           
  r13 = VMAND(r13, 0x55, r13, vbmask);        
  c = VSRA(r10, VBRADIX);
  r11 = VMADD(r11, 0xAA, r11, c);
  r10 = VMAND(r10, 0xAA, r10, vbmask);
  c = VSRA(r11, VBRADIX);
  r12 = VMADD(r12, 0xAA, r12, c);
  r11 = VMAND(r11, 0xAA, r11, vbmask);
  // c = VSRA(r12, VBRADIX);
  // r13 = VMADD(r13, 0xAA, r13, c);
  // r12 = VMAND(r12, 0xAA, r12, vbmask);

  smask = VSRA(r12, 63);                      // smask = all-0/-1 X 
  smask = VSHUF(smask, 0xEE);                 // smask = all-0/-1
  // t = b1 | (p*2^459)&smask
  t0  = VMBLEND(0x33, b0,  zero); t1  = VMBLEND(0x33, b1,  zero);
  t2  = VMBLEND(0x33, b2,  zero); t3  = VMBLEND(0x33, b3,  zero);
  t4  = VMBLEND(0x33, b4,  zero); t5  = VMBLEND(0x33, b5,  zero);
  t6  = VMBLEND(0x33, b6,  zero); t7  = VMBLEND(0x33, b7,  zero);
  t8  = VMBLEND(0x33, b8,  zero);
  t9  = VMBLEND(0x33, b9,  VAND(vp0, smask)); 
  t10 = VMBLEND(0x33, b10, VAND(vp1, smask));  
  t11 = VMBLEND(0x33, b11, VAND(vp2, smask)); 
  t12 = VMBLEND(0x33, b12, VAND(vp3, smask));  
  t13 = VMBLEND(0x33, b13, VAND(vp4, smask));

  // t8  = VMBLEND(0x33, b8,  VAND(vp1, smask));
  // t9  = VMBLEND(0x33, b9,  VAND(vp1, smask)); 
  // t10 = VMBLEND(0x33, b10, VAND(vp2, smask));  
  // t11 = VMBLEND(0x33, b11, VAND(vp3, smask)); 
  // t12 = VMBLEND(0x33, b12, VAND(vp4, smask));  
  // t13 = VMBLEND(0x33, b13, zero);


  // r = r + t = a1 + b1 | a0-b0 + ((p*2^459)&smask)
  r0  = VADD(r0,  t0);  r1  = VADD(r1,  t1);  r2  = VADD(r2,  t2);
  r3  = VADD(r3,  t3);  r4  = VADD(r4,  t4);  r5  = VADD(r5,  t5);
  r6  = VADD(r6,  t6);  r7  = VADD(r7,  t7);  r8  = VADD(r8,  t8);
  r9  = VADD(r9,  t9);  r10 = VADD(r10, t10); r11 = VADD(r11, t11);
  r12 = VADD(r12, t12); r13 = VADD(r13, t13);

  // carry propagation
  // will be done in reduction 

  r[0]  = r0;  r[1]  = r1;  r[2]  = r2;  r[3]  = r3;  r[4]  = r4;  
  r[5]  = r5;  r[6]  = r6;  r[7]  = r7;  r[8]  = r8;  r[9]  = r9;  
  r[10] = r10; r[11] = r11; r[12] = r12; r[13] = r13;
}

// Montgomery multplication 
void fpmul_mont_4x2w(__m512i *r, const __m512i *a, const __m512i *b)
{
  vdgelm_t tt;

  mp_mul_4x2w(tt, a, b);
  rdc_mont_4x2w(r, tt);
}

// integer addition r = a + b
void mp_add_4x2w(__m512i *r, const __m512i *a, const __m512i *b)
{
  __m512i a0 = a[0], a1 = a[1], a2 = a[2], a3 = a[3], a4 = a[4];
  __m512i b0 = b[0], b1 = b[1], b2 = b[2], b3 = b[3], b4 = b[4];
  __m512i r0, r1, r2, r3, r4, c;
  const __m512i vbmask = VSET1(VBMASK);

  // r = a +b
  r0 = VADD(a0, b0); r1 = VADD(a1, b1); r2 = VADD(a2, b2); 
  r3 = VADD(a3, b3); r4 = VADD(a4, b4);

  // *simple* carry propagation
  c = VSRA(r0, VBRADIX); r0 = VAND(r0, vbmask); r1 = VADD(r1, c);
  c = VSRA(r1, VBRADIX); r1 = VAND(r1, vbmask); r2 = VADD(r2, c);
  c = VSRA(r2, VBRADIX); r2 = VAND(r2, vbmask); r3 = VADD(r3, c);
  c = VSRA(r3, VBRADIX); r3 = VAND(r3, vbmask); r4 = VADD(r4, c);
  c = VSRA(r4, VBRADIX); r4 = VAND(r4, vbmask);
  r0 = VMADD(r0, 0xAA, r0, VSHUF(c, 0x4E));

  r[0] = r0; r[1] = r1; r[2] = r2; r[3] = r3; r[4] = r4;
}

// GF(p^2) multiplication using Montgomery arithmetic, r = a*b in GF(p^2).
// Here uses schoolbook instead of Karatsuba.
// Assume inputs a and b are in the format of <a1 | a0> and <b1 | b0>. 
void fp2mul_mont_2x2x2w(vgelm_t r, const vgelm_t a, const vgelm_t b)
{
  vgelm_t t1, t2, t3;
  vdgelm_t tt1, tt2;

  vec_perm1010_4x2w(t1, a);             // t1 = a0 | a0
  vec_perm3232_4x2w(t2, a);             // t2 = a1 | a1
  vec_perm1032_4x2w(t3, b);             // t3 = b0 | b1 
  mp_mul_4x2w(tt1, t1, b);              // tt1 = a0*b1 | a0*b0
  mp_mul_4x2w(tt2, t2, t3);             // tt2 = a1*b0 | a1*b1
  mp_mixaddsub_4x2w(tt1, tt1, tt2);     // tt3 = a0b1 + a1b0 | a0b0 - a1b1
  rdc_mont_4x2w(r, tt1);                // r = r1 | r0 = (a0b1+a1b0) mod 2p | (a0b0-a1b1) mod 2p
}

// GF(p^2) squaring using Montgomery arithmetic r = a^2 in GF(p^2).
void fp2sqr_mont_2x2x2w(vgelm_t r, const vgelm_t a)
{
  vgelm_t t1, t2;

  vec_perm1010_4x2w(t1, a);             // t1 = a0 | a0
  vec_perm1032_4x2w(t2, a);             // t2 = a0 | a1
  mp_add_4x2w(t1, t1, t2);              // t1 = 2*a0 | a0+a1
  vec_permzz32_4x2w(t2, a);             // t2 = 0 | a1  
  mp_sub_p4_4x2w(t2, a, t2);            // t2 = a1 | a0-a1
  fpmul_mont_4x2w(r, t1, t2);           // r = r1 | r0 = 2a0*a1 | (a0+a1)*(a0-a1) 
}

// -----------------------------------------------------------------------------
// 1-way x64 implementation for radix-2^51

void fpcopy_1w(felm_r51_t r, const felm_r51_t a)
{
  int i;

  for (i = 0; i < VNWORDS; i++) r[i] = a[i];
}

void fp2copy_1w(f2elm_r51_t r, const f2elm_r51_t a)
{
  fpcopy_1w(r[0], a[0]);
  fpcopy_1w(r[1], a[1]);
}

void fpzero_1w(felm_r51_t r)
{
  int i;

  for (i = 0; i < VNWORDS; i++) r[i] = 0;
}

void mp_add_1w(felm_r51_t r, const felm_r51_t a, const felm_r51_t b)
{
  uint64_t a0 = a[0], a1 = a[1], a2 = a[2], a3 = a[3], a4 = a[4];
  uint64_t a5 = a[5], a6 = a[6], a7 = a[7], a8 = a[8];
  uint64_t b0 = b[0], b1 = b[1], b2 = b[2], b3 = b[3], b4 = b[4];
  uint64_t b5 = b[5], b6 = b[6], b7 = b[7], b8 = b[8];
  uint64_t r0, r1, r2, r3, r4, r5, r6, r7, r8;

  // r = a + b
  r0 = a0 + b0; r1 = a1 + b1; r2 = a2 + b2; 
  r3 = a3 + b3; r4 = a4 + b4; r5 = a5 + b5; 
  r6 = a6 + b6; r7 = a7 + b7; r8 = a8 + b8; 

  // carry propagation
  r1 += r0>>VBRADIX; r0 &= VBMASK;
  r2 += r1>>VBRADIX; r1 &= VBMASK;
  r3 += r2>>VBRADIX; r2 &= VBMASK;
  r4 += r3>>VBRADIX; r3 &= VBMASK;
  r5 += r4>>VBRADIX; r4 &= VBMASK;
  r6 += r5>>VBRADIX; r5 &= VBMASK;
  r7 += r6>>VBRADIX; r6 &= VBMASK;
  r8 += r7>>VBRADIX; r7 &= VBMASK;

  r[0] = r0; r[1] = r1; r[2] = r2; r[3] = r3; r[4] = r4; 
  r[5] = r5; r[6] = r6; r[7] = r7; r[8] = r8; 
}

void mp2_add_1w(f2elm_r51_t r, const f2elm_r51_t a, const f2elm_r51_t b)
{
  mp_add_1w(r[0], a[0], b[0]);
  mp_add_1w(r[1], a[1], b[1]);
}

// -----------------------------------------------------------------------------
// 1-way x64 implementation from PQCrypto-SIDH-3.4

void decode_to_digits(const unsigned char* x, digit_t* dec, int nbytes, int ndigits)
{ // Decoding bytes to digits according to endianness

    dec[ndigits - 1] = 0;
    memcpy((unsigned char*)dec, x, nbytes);
}

void fpmul_mont(const felm_t ma, const felm_t mb, felm_t mc)
{ // Multiprecision multiplication, c = a*b mod p.
    dfelm_t temp = {0};

    mp_mul(ma, mb, temp, NWORDS_FIELD);
    rdc_mont(temp, mc);
}

void from_mont(const felm_t ma, felm_t c)
{ // Conversion from Montgomery representation to standard representation,
  // c = ma*R^(-1) mod p = a mod p, where ma in [0, p-1].
    digit_t one[NWORDS_FIELD] = {0};
    
    one[0] = 1;
    fpmul_mont(ma, one, c);
    fpcorrection(c);
}

static void mp_addfast(const digit_t* a, const digit_t* b, digit_t* c)
{ // Multiprecision addition, c = a+b.        
    mp_add434_asm(a, b, c);    
}

static void mp_subaddfast(const digit_t* a, const digit_t* b, digit_t* c)
{ // Multiprecision subtraction followed by addition with p*2^MAXBITS_FIELD, c = a-b+(p*2^MAXBITS_FIELD) if a-b < 0, otherwise c=a-b. 
    mp_subadd434x2_asm(a, b, c);     
}

static void mp_dblsubfast(const digit_t* a, const digit_t* b, digit_t* c)
{ // Multiprecision subtraction, c = c-a-b, where lng(a) = lng(b) = 2*NWORDS_FIELD.
    mp_dblsub434x2_asm(a, b, c);
}

void fp2mul_mont(const f2elm_t a, const f2elm_t b, f2elm_t c)
{ // GF(p^2) multiplication using Montgomery arithmetic, c = a*b in GF(p^2).
  // Inputs: a = a0+a1*i and b = b0+b1*i, where a0, a1, b0, b1 are in [0, 2*p-1] 
  // Output: c = c0+c1*i, where c0, c1 are in [0, 2*p-1] 
    felm_t t1, t2;
    dfelm_t tt1, tt2, tt3; 
    
    mp_addfast(a[0], a[1], t1);                      // t1 = a0+a1
    mp_addfast(b[0], b[1], t2);                      // t2 = b0+b1
    mp_mul(a[0], b[0], tt1, NWORDS_FIELD);           // tt1 = a0*b0
    mp_mul(a[1], b[1], tt2, NWORDS_FIELD);           // tt2 = a1*b1
    mp_mul(t1, t2, tt3, NWORDS_FIELD);               // tt3 = (a0+a1)*(b0+b1)
    mp_dblsubfast(tt1, tt2, tt3);                    // tt3 = (a0+a1)*(b0+b1) - a0*b0 - a1*b1
    mp_subaddfast(tt1, tt2, tt1);                    // tt1 = a0*b0 - a1*b1 + p*2^MAXBITS_FIELD if a0*b0 - a1*b1 < 0, else tt1 = a0*b0 - a1*b1
    rdc_mont(tt3, c[1]);                             // c[1] = (a0+a1)*(b0+b1) - a0*b0 - a1*b1 
    rdc_mont(tt1, c[0]);                             // c[0] = a0*b0 - a1*b1
}

void fpcopy(const felm_t a, felm_t c)
{ // Copy a field element, c = a.
    unsigned int i;

    for (i = 0; i < NWORDS_FIELD; i++)
        c[i] = a[i];
}

void fp2copy(const f2elm_t a, f2elm_t c)
{ // Copy a GF(p^2) element, c = a.
    fpcopy(a[0], c[0]);
    fpcopy(a[1], c[1]);
}

void fpsqr_mont(const felm_t ma, felm_t mc)
{ // Multiprecision squaring, c = a^2 mod p.
    dfelm_t temp = {0};

    mp_mul(ma, ma, temp, NWORDS_FIELD);
    rdc_mont(temp, mc);
}

void fpinv_chain_mont(felm_t a)
{ // Chain to compute a^(p-3)/4 using Montgomery arithmetic.
    unsigned int i, j;

    felm_t t[31], tt;

    // Precomputed table
    fpsqr_mont(a, tt);
    fpmul_mont(a, tt, t[0]);
    for (i = 0; i <= 29; i++) fpmul_mont(t[i], tt, t[i+1]);

    fpcopy(a, tt);
    for (i = 0; i < 7; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[5], tt, tt);
    for (i = 0; i < 10; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[14], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[3], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[23], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[13], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[24], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[7], tt, tt);
    for (i = 0; i < 8; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[12], tt, tt);
    for (i = 0; i < 8; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[30], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[1], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[30], tt, tt);
    for (i = 0; i < 7; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[21], tt, tt);
    for (i = 0; i < 9; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[2], tt, tt);
    for (i = 0; i < 9; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[19], tt, tt);
    for (i = 0; i < 9; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[1], tt, tt);
    for (i = 0; i < 7; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[24], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[26], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[16], tt, tt);
    for (i = 0; i < 7; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[10], tt, tt);
    for (i = 0; i < 7; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[6], tt, tt);
    for (i = 0; i < 7; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[0], tt, tt);
    for (i = 0; i < 9; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[20], tt, tt);
    for (i = 0; i < 8; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[9], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[25], tt, tt);
    for (i = 0; i < 9; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[30], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[26], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(a, tt, tt);
    for (i = 0; i < 7; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[28], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[6], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[10], tt, tt);
    for (i = 0; i < 9; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[22], tt, tt);
    for (j = 0; j < 35; j++) {
        for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
        fpmul_mont(t[30], tt, tt);
    }
    fpcopy(tt, a);   
}

void fpinv_mont(felm_t a)
{ // Field inversion using Montgomery arithmetic, a = a^(-1)*R mod p.
    felm_t tt;

    fpcopy(a, tt);
    fpinv_chain_mont(tt);
    fpsqr_mont(tt, tt);
    fpsqr_mont(tt, tt);
    fpmul_mont(a, tt, a);
}

void fp2inv_mont(f2elm_t a)
{// GF(p^2) inversion using Montgomery arithmetic, a = (a0-i*a1)/(a0^2+a1^2).
    f2elm_t t1;

    fpsqr_mont(a[0], t1[0]);                         // t10 = a0^2
    fpsqr_mont(a[1], t1[1]);                         // t11 = a1^2
    fpadd(t1[0], t1[1], t1[0]);                      // t10 = a0^2+a1^2
    fpinv_mont(t1[0]);                               // t10 = (a0^2+a1^2)^-1
    fpneg(a[1]);                                     // a = a0-i*a1
    fpmul_mont(a[0], t1[0], a[0]);
    fpmul_mont(a[1], t1[0], a[1]);                   // a = (a0-i*a1)*(a0^2+a1^2)^-1
}

static void encode_to_bytes(const digit_t* x, unsigned char* enc, int nbytes)
{ // Encoding digits to bytes according to endianness  
    memcpy(enc, (const unsigned char*)x, nbytes);
}

void from_fp2mont(const f2elm_t ma, f2elm_t c)
{ // Conversion of a GF(p^2) element from Montgomery representation to standard representation,
  // c_i = ma_i*R^(-1) = a_i in GF(p^2).

    from_mont(ma[0], c[0]);
    from_mont(ma[1], c[1]);
}

void fp2_encode(const f2elm_t x, unsigned char *enc)
{ // Conversion of GF(p^2) element from Montgomery to standard representation, and encoding by removing leading 0 bytes
    f2elm_t t;

    from_fp2mont(x, t);
    encode_to_bytes(t[0], enc, FP2_ENCODED_BYTES / 2);
    encode_to_bytes(t[1], enc + FP2_ENCODED_BYTES / 2, FP2_ENCODED_BYTES / 2);
}

void to_mont(const felm_t a, felm_t mc)
{ // Conversion to Montgomery representation,
  // mc = a*R^2*R^(-1) mod p = a*R mod p, where a in [0, p-1].
  // The Montgomery constant R^2 mod p is the global value "Montgomery_R2". 

    fpmul_mont(a, (digit_t*)&Montgomery_R2, mc);
}

void to_fp2mont(const f2elm_t a, f2elm_t mc)
{ // Conversion of a GF(p^2) element to Montgomery representation,
  // mc_i = a_i*R^2*R^(-1) = a_i*R in GF(p^2). 

    to_mont(a[0], mc[0]);
    to_mont(a[1], mc[1]);
}

void fp2_decode(const unsigned char *x, f2elm_t dec)
{ // Parse byte sequence back into GF(p^2) element, and conversion to Montgomery representation

    decode_to_digits(x, dec[0], FP2_ENCODED_BYTES / 2, NWORDS_FIELD);
    decode_to_digits(x + FP2_ENCODED_BYTES / 2, dec[1], FP2_ENCODED_BYTES / 2, NWORDS_FIELD);
    to_fp2mont(dec, dec);
}

void fp2sub(const f2elm_t a, const f2elm_t b, f2elm_t c)          
{ // GF(p^2) subtraction, c = a-b in GF(p^2).
    fpsub(a[0], b[0], c[0]);
    fpsub(a[1], b[1], c[1]);
}

void fp2add(const f2elm_t a, const f2elm_t b, f2elm_t c)           
{ // GF(p^2) addition, c = a+b in GF(p^2).
    fpadd(a[0], b[0], c[0]);
    fpadd(a[1], b[1], c[1]);
}

void fp2sqr_mont(const f2elm_t a, f2elm_t c)
{ // GF(p^2) squaring using Montgomery arithmetic, c = a^2 in GF(p^2).
  // Inputs: a = a0+a1*i, where a0, a1 are in [0, 2*p-1] 
  // Output: c = c0+c1*i, where c0, c1 are in [0, 2*p-1] 
    felm_t t1, t2, t3;
    
    mp_addfast(a[0], a[1], t1);                      // t1 = a0+a1 
    mp_sub_p4(a[0], a[1], t2);                       // t2 = a0-a1
    mp_addfast(a[0], a[0], t3);                      // t3 = 2a0
    fpmul_mont(t1, t2, c[0]);                        // c0 = (a0+a1)(a0-a1)
    fpmul_mont(t3, a[1], c[1]);                      // c1 = 2a0*a1
}

void mp2_add(const f2elm_t a, const f2elm_t b, f2elm_t c)       
{ // GF(p^2) addition without correction, c = a+b in GF(p^2). 
    mp_addfast(a[0], b[0], c[0]);
    mp_addfast(a[1], b[1], c[1]);
}
