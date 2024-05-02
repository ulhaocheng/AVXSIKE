/**
 *******************************************************************************
 * @version 0.1
 * @date 2021-12-16
 * @copyright Copyright © 2021 by University of Luxembourg
 * @author Hao Cheng
 *******************************************************************************
 */

#include "fp.h"

// -----------------------------------------------------------------------------
// (8x1)-way fp arithmetic 

// integer subtraction r = a - b + 2p
void mp_sub_p2_8x1w(__m512i *r, const __m512i *a, const __m512i *b)
{
  __m512i a0 = a[0], a1 = a[1], a2  = a[2],  a3  = a[3];
  __m512i a4 = a[4], a5 = a[5], a6  = a[6],  a7  = a[7];
  __m512i a8 = a[8], a9 = a[9], a10 = a[10], a11 = a[11];
  __m512i b0 = b[0], b1 = b[1], b2  = b[2],  b3  = b[3]; 
  __m512i b4 = b[4], b5 = b[5], b6  = b[6],  b7  = b[7];
  __m512i b8 = b[8], b9 = b[9], b10 = b[10], b11 = b[11];
  __m512i r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11;
  const __m512i vp0  = VSET1(vp610x2[0]),  vp1  = VSET1(vp610x2[1]);
  const __m512i vp2  = VSET1(vp610x2[2]),  vp3  = VSET1(vp610x2[3]);
  const __m512i vp4  = VSET1(vp610x2[4]),  vp5  = VSET1(vp610x2[5]);
  const __m512i vp6  = VSET1(vp610x2[6]),  vp7  = VSET1(vp610x2[7]);
  const __m512i vp8  = VSET1(vp610x2[8]),  vp9  = VSET1(vp610x2[9]);
  const __m512i vp10 = VSET1(vp610x2[10]), vp11 = VSET1(vp610x2[11]);
  const __m512i vbmask = VSET1(VBMASK); 

  // r = a + 2p 
  r0 = VADD(a0, vp0); r1  = VADD(a1,  vp1);  r2  = VADD(a2,  vp2); 
  r3 = VADD(a3, vp3); r4  = VADD(a4,  vp4);  r5  = VADD(a5,  vp5);
  r6 = VADD(a6, vp6); r7  = VADD(a7,  vp7);  r8  = VADD(a8,  vp8);
  r9 = VADD(a9, vp9); r10 = VADD(a10, vp10); r11 = VADD(a11, vp11);

  // r = r - b
  r0 = VSUB(r0, b0); r1  = VSUB(r1,  b1);  r2  = VSUB(r2, b2);
  r3 = VSUB(r3, b3); r4  = VSUB(r4,  b4);  r5  = VSUB(r5, b5);
  r6 = VSUB(r6, b6); r7  = VSUB(r7,  b7);  r8  = VSUB(r8, b8);
  r9 = VSUB(r9, b9); r10 = VSUB(r10, b10); r11 = VSUB(r11, b11);

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

  r[0] = r0; r[1] = r1; r[2]  = r2;  r[3]  = r3; 
  r[4] = r4; r[5] = r5; r[6]  = r6;  r[7]  = r7; 
  r[8] = r8; r[9] = r9; r[10] = r10; r[11] = r11;
}

// integer subtraction r = a - b + 4p
void mp_sub_p4_8x1w(__m512i *r, const __m512i *a, const __m512i *b)
{
  __m512i a0 = a[0], a1 = a[1], a2  = a[2],  a3  = a[3];
  __m512i a4 = a[4], a5 = a[5], a6  = a[6],  a7  = a[7];
  __m512i a8 = a[8], a9 = a[9], a10 = a[10], a11 = a[11];
  __m512i b0 = b[0], b1 = b[1], b2  = b[2],  b3  = b[3]; 
  __m512i b4 = b[4], b5 = b[5], b6  = b[6],  b7  = b[7];
  __m512i b8 = b[8], b9 = b[9], b10 = b[10], b11 = b[11];
  __m512i r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11;
  const __m512i vp0  = VSET1(vp610x4[0]),  vp1  = VSET1(vp610x4[1]);
  const __m512i vp2  = VSET1(vp610x4[2]),  vp3  = VSET1(vp610x4[3]);
  const __m512i vp4  = VSET1(vp610x4[4]),  vp5  = VSET1(vp610x4[5]);
  const __m512i vp6  = VSET1(vp610x4[6]),  vp7  = VSET1(vp610x4[7]);
  const __m512i vp8  = VSET1(vp610x4[8]),  vp9  = VSET1(vp610x4[9]);
  const __m512i vp10 = VSET1(vp610x4[10]), vp11 = VSET1(vp610x4[11]);
  const __m512i vbmask = VSET1(VBMASK); 

  // r = a + 4p 
  r0 = VADD(a0, vp0); r1  = VADD(a1,  vp1);  r2  = VADD(a2,  vp2); 
  r3 = VADD(a3, vp3); r4  = VADD(a4,  vp4);  r5  = VADD(a5,  vp5);
  r6 = VADD(a6, vp6); r7  = VADD(a7,  vp7);  r8  = VADD(a8,  vp8);
  r9 = VADD(a9, vp9); r10 = VADD(a10, vp10); r11 = VADD(a11, vp11);

  // r = r - b
  r0 = VSUB(r0, b0); r1  = VSUB(r1,  b1);  r2  = VSUB(r2, b2);
  r3 = VSUB(r3, b3); r4  = VSUB(r4,  b4);  r5  = VSUB(r5, b5);
  r6 = VSUB(r6, b6); r7  = VSUB(r7,  b7);  r8  = VSUB(r8, b8);
  r9 = VSUB(r9, b9); r10 = VSUB(r10, b10); r11 = VSUB(r11, b11);

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

  r[0] = r0; r[1] = r1; r[2]  = r2;  r[3]  = r3; 
  r[4] = r4; r[5] = r5; r[6]  = r6;  r[7]  = r7; 
  r[8] = r8; r[9] = r9; r[10] = r10; r[11] = r11;
}

// modular addition r = a + b mod 2p
void fpadd_8x1w(__m512i *r, const __m512i *a, const __m512i* b)
{
  __m512i a0 = a[0], a1 = a[1], a2  = a[2],  a3  = a[3];
  __m512i a4 = a[4], a5 = a[5], a6  = a[6],  a7  = a[7];
  __m512i a8 = a[8], a9 = a[9], a10 = a[10], a11 = a[11];
  __m512i b0 = b[0], b1 = b[1], b2  = b[2],  b3  = b[3]; 
  __m512i b4 = b[4], b5 = b[5], b6  = b[6],  b7  = b[7];
  __m512i b8 = b[8], b9 = b[9], b10 = b[10], b11 = b[11];
  __m512i r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, smask;
  const __m512i vp0  = VSET1(vp610x2[0]),  vp1  = VSET1(vp610x2[1]);
  const __m512i vp2  = VSET1(vp610x2[2]),  vp3  = VSET1(vp610x2[3]);
  const __m512i vp4  = VSET1(vp610x2[4]),  vp5  = VSET1(vp610x2[5]);
  const __m512i vp6  = VSET1(vp610x2[6]),  vp7  = VSET1(vp610x2[7]);
  const __m512i vp8  = VSET1(vp610x2[8]),  vp9  = VSET1(vp610x2[9]);
  const __m512i vp10 = VSET1(vp610x2[10]), vp11 = VSET1(vp610x2[11]);
  const __m512i vbmask = VSET1(VBMASK); 

  // r = a + b
  r0 = VADD(a0, b0); r1  = VADD(a1,  b1);  r2  = VADD(a2, b2); 
  r3 = VADD(a3, b3); r4  = VADD(a4,  b4);  r5  = VADD(a5, b5);
  r6 = VADD(a6, b6); r7  = VADD(a7,  b7);  r8  = VADD(a8, b8);
  r9 = VADD(a9, b9); r10 = VADD(a10, b10); r11 = VADD(a11, b11);

  // r = r - 2p
  r0 = VSUB(r0, vp0); r1  = VSUB(r1,  vp1);  r2  = VSUB(r2, vp2);
  r3 = VSUB(r3, vp3); r4  = VSUB(r4,  vp4);  r5  = VSUB(r5, vp5);
  r6 = VSUB(r6, vp6); r7  = VSUB(r7,  vp7);  r8  = VSUB(r8, vp8);
  r9 = VSUB(r9, vp9); r10 = VSUB(r10, vp10); r11 = VSUB(r11, vp11);

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

  r[0] = r0; r[1] = r1; r[2]  = r2;  r[3]  = r3; 
  r[4] = r4; r[5] = r5; r[6]  = r6;  r[7]  = r7; 
  r[8] = r8; r[9] = r9; r[10] = r10; r[11] = r11;
}

// modular subtraction r = a - b mod 2p
void fpsub_8x1w(__m512i *r, const __m512i *a, const __m512i* b)
{
  __m512i a0 = a[0], a1 = a[1], a2  = a[2],  a3  = a[3];
  __m512i a4 = a[4], a5 = a[5], a6  = a[6],  a7  = a[7];
  __m512i a8 = a[8], a9 = a[9], a10 = a[10], a11 = a[11];
  __m512i b0 = b[0], b1 = b[1], b2  = b[2],  b3  = b[3]; 
  __m512i b4 = b[4], b5 = b[5], b6  = b[6],  b7  = b[7];
  __m512i b8 = b[8], b9 = b[9], b10 = b[10], b11 = b[11];
  __m512i r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, smask;
  const __m512i vp0  = VSET1(vp610x2[0]),  vp1  = VSET1(vp610x2[1]);
  const __m512i vp2  = VSET1(vp610x2[2]),  vp3  = VSET1(vp610x2[3]);
  const __m512i vp4  = VSET1(vp610x2[4]),  vp5  = VSET1(vp610x2[5]);
  const __m512i vp6  = VSET1(vp610x2[6]),  vp7  = VSET1(vp610x2[7]);
  const __m512i vp8  = VSET1(vp610x2[8]),  vp9  = VSET1(vp610x2[9]);
  const __m512i vp10 = VSET1(vp610x2[10]), vp11 = VSET1(vp610x2[11]);
  const __m512i vbmask = VSET1(VBMASK); 

  // r = a - b
  r0 = VSUB(a0, b0); r1  = VSUB(a1,  b1);  r2  = VSUB(a2, b2); 
  r3 = VSUB(a3, b3); r4  = VSUB(a4,  b4);  r5  = VSUB(a5, b5);
  r6 = VSUB(a6, b6); r7  = VSUB(a7,  b7);  r8  = VSUB(a8, b8);
  r9 = VSUB(a9, b9); r10 = VSUB(a10, b10); r11 = VSUB(a11, b11);

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

  r[0] = r0; r[1] = r1; r[2]  = r2;  r[3]  = r3; 
  r[4] = r4; r[5] = r5; r[6]  = r6;  r[7]  = r7; 
  r[8] = r8; r[9] = r9; r[10] = r10; r[11] = r11; 
}

// integer multiplication r = a * b (product-scanning)
// len(r) = 2*NWORDS; len(a) = len(b) = NOWRDS
// the carry propagation is "interleaved with" the multiplication 
void mp_mul_8x1w_v1(__m512i *r, const __m512i *a, const __m512i *b)
{
  __m512i a0 = a[0], a1 = a[1], a2  = a[2],  a3  = a[3];
  __m512i a4 = a[4], a5 = a[5], a6  = a[6],  a7  = a[7];
  __m512i a8 = a[8], a9 = a[9], a10 = a[10], a11 = a[11];
  __m512i b0 = b[0], b1 = b[1], b2  = b[2],  b3  = b[3]; 
  __m512i b4 = b[4], b5 = b[5], b6  = b[6],  b7  = b[7];
  __m512i b8 = b[8], b9 = b[9], b10 = b[10], b11 = b[11];
  __m512i r0,  r1,  r2,  r3,  r4,  r5,  r6,  r7,  r8,  r9,  r10, r11;
  __m512i r12, r13, r14, r15, r16, r17, r18, r19, r20, r21, r22, r23;
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
  const __m512i vbmask = VSET1(VBMASK);

  // ---------------------------------------------------------------------------
  // integer multiplication 

  z0 = VMACLO(z0, a0, b0);
  h1 = VMACHI(h1, a0, b0); h1 = VADD(h1, h1);
  r0 = z0; z1 = VADD(z1, h1);

  z1 = VMACLO(z1, a0, b1); z1 = VMACLO(z1, a1, b0); 
  h2 = VMACHI(h2, a0, b1); h2 = VMACHI(h2, a1, b0); h2 = VADD(h2, h2);
  r1 = VAND(z1, vbmask); z2 = VADD(h2, VSRA(z1, VBRADIX));
  
  z2 = VMACLO(z2, a0, b2); z2 = VMACLO(z2, a1, b1); z2 = VMACLO(z2, a2, b0);
  h3 = VMACHI(h3, a0, b2); h3 = VMACHI(h3, a1, b1); h3 = VMACHI(h3, a2, b0);
  h3 = VADD(h3, h3);
  r2 = VAND(z2, vbmask); z3 = VADD(h3, VSRA(z2, VBRADIX));

  z3 = VMACLO(z3, a0, b3); z3 = VMACLO(z3, a1, b2); z3 = VMACLO(z3, a2, b1); 
  z3 = VMACLO(z3, a3, b0);
  h4 = VMACHI(h4, a0, b3); h4 = VMACHI(h4, a1, b2); h4 = VMACHI(h4, a2, b1); 
  h4 = VMACHI(h4, a3, b0); h4 = VADD(h4, h4);
  r3 = VAND(z3, vbmask); z4 = VADD(h4, VSRA(z3, VBRADIX));

  z4 = VMACLO(z4, a0, b4); z4 = VMACLO(z4, a1, b3); z4 = VMACLO(z4, a2, b2); 
  z4 = VMACLO(z4, a3, b1); z4 = VMACLO(z4, a4, b0);
  h5 = VMACHI(h5, a0, b4); h5 = VMACHI(h5, a1, b3); h5 = VMACHI(h5, a2, b2); 
  h5 = VMACHI(h5, a3, b1); h5 = VMACHI(h5, a4, b0); h5 = VADD(h5, h5);
  r4 = VAND(z4, vbmask); z5 = VADD(h5, VSRA(z4, VBRADIX));

  z5 = VMACLO(z5, a0, b5); z5 = VMACLO(z5, a1, b4); z5 = VMACLO(z5, a2, b3); 
  z5 = VMACLO(z5, a3, b2); z5 = VMACLO(z5, a4, b1); z5 = VMACLO(z5, a5, b0);
  h6 = VMACHI(h6, a0, b5); h6 = VMACHI(h6, a1, b4); h6 = VMACHI(h6, a2, b3); 
  h6 = VMACHI(h6, a3, b2); h6 = VMACHI(h6, a4, b1); h6 = VMACHI(h6, a5, b0);
  h6 = VADD(h6, h6);
  r5 = VAND(z5, vbmask); z6 = VADD(h6, VSRA(z5, VBRADIX));

  z6 = VMACLO(z6, a0, b6); z6 = VMACLO(z6, a1, b5); z6 = VMACLO(z6, a2, b4); 
  z6 = VMACLO(z6, a3, b3); z6 = VMACLO(z6, a4, b2); z6 = VMACLO(z6, a5, b1); 
  z6 = VMACLO(z6, a6, b0);
  h7 = VMACHI(h7, a0, b6); h7 = VMACHI(h7, a1, b5); h7 = VMACHI(h7, a2, b4); 
  h7 = VMACHI(h7, a3, b3); h7 = VMACHI(h7, a4, b2); h7 = VMACHI(h7, a5, b1); 
  h7 = VMACHI(h7, a6, b0); h7 = VADD(h7, h7);
  r6 = VAND(z6, vbmask); z7 = VADD(h7, VSRA(z6, VBRADIX));

  z7 = VMACLO(z7, a0, b7); z7 = VMACLO(z7, a1, b6); z7 = VMACLO(z7, a2, b5); 
  z7 = VMACLO(z7, a3, b4); z7 = VMACLO(z7, a4, b3); z7 = VMACLO(z7, a5, b2); 
  z7 = VMACLO(z7, a6, b1); z7 = VMACLO(z7, a7, b0);
  h8 = VMACHI(h8, a0, b7); h8 = VMACHI(h8, a1, b6); h8 = VMACHI(h8, a2, b5); 
  h8 = VMACHI(h8, a3, b4); h8 = VMACHI(h8, a4, b3); h8 = VMACHI(h8, a5, b2); 
  h8 = VMACHI(h8, a6, b1); h8 = VMACHI(h8, a7, b0); h8 = VADD(h8, h8);
  r7 = VAND(z7, vbmask); z8 = VADD(h8, VSRA(z7, VBRADIX));

  z8 = VMACLO(z8, a0, b8); z8 = VMACLO(z8, a1, b7); z8 = VMACLO(z8, a2, b6); 
  z8 = VMACLO(z8, a3, b5); z8 = VMACLO(z8, a4, b4); z8 = VMACLO(z8, a5, b3); 
  z8 = VMACLO(z8, a6, b2); z8 = VMACLO(z8, a7, b1); z8 = VMACLO(z8, a8, b0);
  h9 = VMACHI(h9, a0, b8); h9 = VMACHI(h9, a1, b7); h9 = VMACHI(h9, a2, b6); 
  h9 = VMACHI(h9, a3, b5); h9 = VMACHI(h9, a4, b4); h9 = VMACHI(h9, a5, b3); 
  h9 = VMACHI(h9, a6, b2); h9 = VMACHI(h9, a7, b1); h9 = VMACHI(h9, a8, b0);
  h9 = VADD(h9, h9);
  r8 = VAND(z8, vbmask); z9 = VADD(h9, VSRA(z8, VBRADIX));

  z9 = VMACLO(z9, a0, b9); z9 = VMACLO(z9, a1, b8); z9 = VMACLO(z9, a2, b7); 
  z9 = VMACLO(z9, a3, b6); z9 = VMACLO(z9, a4, b5); z9 = VMACLO(z9, a5, b4); 
  z9 = VMACLO(z9, a6, b3); z9 = VMACLO(z9, a7, b2); z9 = VMACLO(z9, a8, b1); 
  z9 = VMACLO(z9, a9, b0);
  h10 = VMACHI(h10, a0, b9); h10 = VMACHI(h10, a1, b8); h10 = VMACHI(h10, a2, b7); 
  h10 = VMACHI(h10, a3, b6); h10 = VMACHI(h10, a4, b5); h10 = VMACHI(h10, a5, b4); 
  h10 = VMACHI(h10, a6, b3); h10 = VMACHI(h10, a7, b2); h10 = VMACHI(h10, a8, b1); 
  h10 = VMACHI(h10, a9, b0); h10 = VADD(h10, h10);
  r9 = VAND(z9, vbmask); z10 = VADD(h10, VSRA(z9, VBRADIX));

  z10 = VMACLO(z10, a1, b9); z10 = VMACLO(z10, a2, b8); z10 = VMACLO(z10, a3, b7); 
  z10 = VMACLO(z10, a4, b6); z10 = VMACLO(z10, a5, b5); z10 = VMACLO(z10, a6, b4); 
  z10 = VMACLO(z10, a7, b3); z10 = VMACLO(z10, a8, b2); z10 = VMACLO(z10, a9, b1);
  z10 = VMACLO(z10, a0, b10); z10 = VMACLO(z10, a10, b0);
  h11 = VMACHI(h11, a1, b9); h11 = VMACHI(h11, a2, b8); h11 = VMACHI(h11, a3, b7); 
  h11 = VMACHI(h11, a4, b6); h11 = VMACHI(h11, a5, b5); h11 = VMACHI(h11, a6, b4); 
  h11 = VMACHI(h11, a7, b3); h11 = VMACHI(h11, a8, b2); h11 = VMACHI(h11, a9, b1); 
  h11 = VMACHI(h11, a0, b10); h11 = VMACHI(h11, a10, b0);
  h11 = VADD(h11, h11);
  r10 = VAND(z10, vbmask); z11 = VADD(h11, VSRA(z10, VBRADIX));

  z11 = VMACLO(z11, a2, b9); z11 = VMACLO(z11, a3, b8); z11 = VMACLO(z11, a4, b7); 
  z11 = VMACLO(z11, a5, b6); z11 = VMACLO(z11, a6, b5); z11 = VMACLO(z11, a7, b4); 
  z11 = VMACLO(z11, a8, b3); z11 = VMACLO(z11, a9, b2);
  z11 = VMACLO(z11, a0, b11); z11 = VMACLO(z11, a11, b0);
  z11 = VMACLO(z11, a1, b10); z11 = VMACLO(z11, a10, b1);
  h12 = VMACHI(h12, a2, b9); h12 = VMACHI(h12, a3, b8); h12 = VMACHI(h12, a4, b7); 
  h12 = VMACHI(h12, a5, b6); h12 = VMACHI(h12, a6, b5); h12 = VMACHI(h12, a7, b4); 
  h12 = VMACHI(h12, a8, b3); h12 = VMACHI(h12, a9, b2); 
  h12 = VMACHI(h12, a0, b11); h12 = VMACHI(h12, a11, b0);
  h12 = VMACHI(h12, a1, b10); h12 = VMACHI(h12, a10, b1);
  h12 = VADD(h12, h12);
  r11 = VAND(z11, vbmask); z12 = VADD(h12, VSRA(z11, VBRADIX));

  z12 = VMACLO(z12, a3, b9); z12 = VMACLO(z12, a4, b8); z12 = VMACLO(z12, a5, b7); 
  z12 = VMACLO(z12, a6, b6); z12 = VMACLO(z12, a7, b5); z12 = VMACLO(z12, a8, b4); 
  z12 = VMACLO(z12, a9, b3);
  z12 = VMACLO(z12, a1, b11); z12 = VMACLO(z12, a2, b10);
  z12 = VMACLO(z12, a11, b1); z12 = VMACLO(z12, a10, b2);
  h13 = VMACHI(h13, a3, b9); h13 = VMACHI(h13, a4, b8); h13 = VMACHI(h13, a5, b7); 
  h13 = VMACHI(h13, a6, b6); h13 = VMACHI(h13, a7, b5); h13 = VMACHI(h13, a8, b4); 
  h13 = VMACHI(h13, a9, b3); 
  h13 = VMACHI(h13, a1, b11); h13 = VMACHI(h13, a2, b10);
  h13 = VMACHI(h13, a11, b1); h13 = VMACHI(h13, a10, b2);
  h13 = VADD(h13, h13);
  r12 = VAND(z12, vbmask); z13 = VADD(h13, VSRA(z12, VBRADIX)); 

  z13 = VMACLO(z13, a4, b9); z13 = VMACLO(z13, a5, b8); z13 = VMACLO(z13, a6, b7); 
  z13 = VMACLO(z13, a7, b6); z13 = VMACLO(z13, a8, b5); z13 = VMACLO(z13, a9, b4);
  z13 = VMACLO(z13, a2, b11); z13 = VMACLO(z13, a3, b10);
  z13 = VMACLO(z13, a11, b2); z13 = VMACLO(z13, a10, b3);
  h14 = VMACHI(h14, a4, b9); h14 = VMACHI(h14, a5, b8); h14 = VMACHI(h14, a6, b7); 
  h14 = VMACHI(h14, a7, b6); h14 = VMACHI(h14, a8, b5); h14 = VMACHI(h14, a9, b4);
  h14 = VMACHI(h14, a2, b11); h14 = VMACHI(h14, a3, b10);
  h14 = VMACHI(h14, a11, b2); h14 = VMACHI(h14, a10, b3);
  h14 = VADD(h14, h14);
  r13 = VAND(z13, vbmask); z14 = VADD(h14, VSRA(z13, VBRADIX));

  z14 = VMACLO(z14, a5, b9); z14 = VMACLO(z14, a6, b8); z14 = VMACLO(z14, a7, b7); 
  z14 = VMACLO(z14, a8, b6); z14 = VMACLO(z14, a9, b5);
  z14 = VMACLO(z14, a3, b11); z14 = VMACLO(z14, a4, b10);
  z14 = VMACLO(z14, a11, b3); z14 = VMACLO(z14, a10, b4);
  h15 = VMACHI(h15, a5, b9); h15 = VMACHI(h15, a6, b8); h15 = VMACHI(h15, a7, b7); 
  h15 = VMACHI(h15, a8, b6); h15 = VMACHI(h15, a9, b5); 
  h15 = VMACHI(h15, a3, b11); h15 = VMACHI(h15, a4, b10);
  h15 = VMACHI(h15, a11, b3); h15 = VMACHI(h15, a10, b4);
  h15 = VADD(h15, h15);
  r14 = VAND(z14, vbmask); z15 = VADD(h15, VSRA(z14, VBRADIX));


  z15 = VMACLO(z15, a6, b9); z15 = VMACLO(z15, a7, b8); z15 = VMACLO(z15, a8, b7); 
  z15 = VMACLO(z15, a9, b6);
  z15 = VMACLO(z15, a4, b11); z15 = VMACLO(z15, a5, b10);
  z15 = VMACLO(z15, a11, b4); z15 = VMACLO(z15, a10, b5);
  h16 = VMACHI(h16, a6, b9); h16 = VMACHI(h16, a7, b8); h16 = VMACHI(h16, a8, b7); 
  h16 = VMACHI(h16, a9, b6); 
  h16 = VMACHI(h16, a4, b11); h16 = VMACHI(h16, a5, b10);
  h16 = VMACHI(h16, a11, b4); h16 = VMACHI(h16, a10, b5);
  h16 = VADD(h16, h16);
  r15 = VAND(z15, vbmask); z16 = VADD(h16, VSRA(z15, VBRADIX));

  z16 = VMACLO(z16, a7, b9); z16 = VMACLO(z16, a8, b8); z16 = VMACLO(z16, a9, b7);
  z16 = VMACLO(z16, a5, b11); z16 = VMACLO(z16, a6, b10);
  z16 = VMACLO(z16, a11, b5); z16 = VMACLO(z16, a10, b6);
  h17 = VMACHI(h17, a7, b9); h17 = VMACHI(h17, a8, b8); h17 = VMACHI(h17, a9, b7);
  h17 = VMACHI(h17, a5, b11); h17 = VMACHI(h17, a6, b10);
  h17 = VMACHI(h17, a11, b5); h17 = VMACHI(h17, a10, b6);
  h17 = VADD(h17, h17);
  r16 = VAND(z16, vbmask); z17 = VADD(h17, VSRA(z16, VBRADIX));
 
  z17 = VMACLO(z17, a8, b9); z17 = VMACLO(z17, a9, b8);
  z17 = VMACLO(z17, a6, b11); z17 = VMACLO(z17, a7, b10);
  z17 = VMACLO(z17, a11, b6); z17 = VMACLO(z17, a10, b7);
  h18 = VMACHI(h18, a8, b9); h18 = VMACHI(h18, a9, b8); 
  h18 = VMACHI(h18, a6, b11); h18 = VMACHI(h18, a7, b10);
  h18 = VMACHI(h18, a11, b6); h18 = VMACHI(h18, a10, b7);
  h18 = VADD(h18, h18);
  r17 = VAND(z17, vbmask); z18 = VADD(h18, VSRA(z17, VBRADIX));

  z18 = VMACLO(z18, a9, b9);
  z18 = VMACLO(z18, a7, b11); z18 = VMACLO(z18, a8, b10);
  z18 = VMACLO(z18, a11, b7); z18 = VMACLO(z18, a10, b8);
  h19 = VMACHI(h19, a9, b9); 
  h19 = VMACHI(h19, a7, b11); h19 = VMACHI(h19, a8, b10);
  h19 = VMACHI(h19, a11, b7); h19 = VMACHI(h19, a10, b8);
  h19 = VADD(h19, h19);
  r18 = VAND(z18, vbmask); z19 = VADD(h19, VSRA(z18, VBRADIX));

  z19 = VMACLO(z19, a8, b11); z19 = VMACLO(z19, a9, b10);
  z19 = VMACLO(z19, a11, b8); z19 = VMACLO(z19, a10, b9);
  h20 = VMACHI(h20, a8, b11); h20 = VMACHI(h20, a9, b10);
  h20 = VMACHI(h20, a11, b8); h20 = VMACHI(h20, a10, b9);
  h20 = VADD(h20, h20);
  r19 = VAND(z19, vbmask); z20 = VADD(h20, VSRA(z19, VBRADIX));

  z20 = VMACLO(z20, a9, b11); z20 = VMACLO(z20, a10, b10);
  z20 = VMACLO(z20, a11, b9); 
  h21 = VMACHI(h21, a9, b11); h21 = VMACHI(h21, a10, b10);
  h21 = VMACHI(h21, a11, b9); 
  h21 = VADD(h21, h21);
  r20 = VAND(z20, vbmask); z21 = VADD(h21, VSRA(z20, VBRADIX));

  z21 = VMACLO(z21, a10, b11); z21 = VMACLO(z21, a11, b10);
  h22 = VMACHI(h22, a10, b11); h22 = VMACHI(h22, a11, b10);
  h22 = VADD(h22, h22);
  r21 = VAND(z21, vbmask); z22 = VADD(h22, VSRA(z21, VBRADIX));

  z22 = VMACLO(z22, a11, b11); 
  h23 = VMACHI(h23, a11, b11); 
  h23 = VADD(h23, h23);
  r22 = VAND(z22, vbmask); z23 = VADD(h23, VSRA(z22, VBRADIX));

  r23 = z23;

  // ---------------------------------------------------------------------------
  r[0]  = r0;  r[1]  = r1;  r[2]  = r2;  r[3]  = r3;  
  r[4]  = r4;  r[5]  = r5;  r[6]  = r6;  r[7]  = r7;  
  r[8]  = r8;  r[9]  = r9;  r[10] = r10; r[11] = r11; 
  r[12] = r12; r[13] = r13; r[14] = r14; r[15] = r15; 
  r[16] = r16; r[17] = r17; r[18] = r18; r[19] = r19;
  r[20] = r20; r[21] = r21; r[22] = r22; r[23] = r23;
}

// Karatsuba
void mp_mul_8x1w_v2(__m512i *r, const __m512i *a, const __m512i *b)
{
  __m512i a0 = a[0], a1 = a[1], a2  = a[2],  a3  = a[3];
  __m512i a4 = a[4], a5 = a[5], a6  = a[6],  a7  = a[7];
  __m512i a8 = a[8], a9 = a[9], a10 = a[10], a11 = a[11];
  __m512i b0 = b[0], b1 = b[1], b2  = b[2],  b3  = b[3]; 
  __m512i b4 = b[4], b5 = b[5], b6  = b[6],  b7  = b[7];
  __m512i b8 = b[8], b9 = b[9], b10 = b[10], b11 = b[11];
  __m512i ta0, ta1, ta2, ta3, ta4, ta5;
  __m512i tb0, tb1, tb2, tb3, tb4, tb5;
  __m512i r0,  r1,  r2,  r3,  r4,  r5,  r6,  r7,  r8,  r9,  r10, r11;
  __m512i r12, r13, r14, r15, r16, r17, r18, r19, r20, r21, r22, r23;
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
  __m512i m0  = VZERO, m1  = VZERO, m2  = VZERO, m3  = VZERO;
  __m512i m4  = VZERO, m5  = VZERO, m6  = VZERO, m7  = VZERO;
  __m512i m8  = VZERO, m9  = VZERO, m10 = VZERO, m11 = VZERO;
  const __m512i vbmask = VSET1(VBMASK);

  // ---------------------------------------------------------------------------
  // compute zL(z0~z11) by aL(a0~a5) * bL(b0~b5)
  z0 = VMACLO(z0, a0, b0);
  h1 = VMACHI(h1, a0, b0); h1 = VADD(h1, h1);

  z1 = VMACLO(h1, a0, b1); z1 = VMACLO(z1, a1, b0); 
  h2 = VMACHI(h2, a0, b1); h2 = VMACHI(h2, a1, b0); h2 = VADD(h2, h2);

  z2 = VMACLO(h2, a0, b2); z2 = VMACLO(z2, a1, b1); z2 = VMACLO(z2, a2, b0);
  h3 = VMACHI(h3, a0, b2); h3 = VMACHI(h3, a1, b1); h3 = VMACHI(h3, a2, b0);
  h3 = VADD(h3, h3);

  z3 = VMACLO(h3, a0, b3); z3 = VMACLO(z3, a1, b2); z3 = VMACLO(z3, a2, b1); 
  z3 = VMACLO(z3, a3, b0);
  h4 = VMACHI(h4, a0, b3); h4 = VMACHI(h4, a1, b2); h4 = VMACHI(h4, a2, b1); 
  h4 = VMACHI(h4, a3, b0); h4 = VADD(h4, h4);

  z4 = VMACLO(h4, a0, b4); z4 = VMACLO(z4, a1, b3); z4 = VMACLO(z4, a2, b2); 
  z4 = VMACLO(z4, a3, b1); z4 = VMACLO(z4, a4, b0);
  h5 = VMACHI(h5, a0, b4); h5 = VMACHI(h5, a1, b3); h5 = VMACHI(h5, a2, b2); 
  h5 = VMACHI(h5, a3, b1); h5 = VMACHI(h5, a4, b0); h5 = VADD(h5, h5);

  z5 = VMACLO(h5, a0, b5); z5 = VMACLO(z5, a1, b4); z5 = VMACLO(z5, a2, b3); 
  z5 = VMACLO(z5, a3, b2); z5 = VMACLO(z5, a4, b1); z5 = VMACLO(z5, a5, b0);
  h6 = VMACHI(h6, a0, b5); h6 = VMACHI(h6, a1, b4); h6 = VMACHI(h6, a2, b3); 
  h6 = VMACHI(h6, a3, b2); h6 = VMACHI(h6, a4, b1); h6 = VMACHI(h6, a5, b0);
  h6 = VADD(h6, h6);

  z6 = VMACLO(h6, a1, b5); z6 = VMACLO(z6, a2, b4); z6 = VMACLO(z6, a3, b3); 
  z6 = VMACLO(z6, a4, b2); z6 = VMACLO(z6, a5, b1);
  h7 = VMACHI(h7, a1, b5); h7 = VMACHI(h7, a2, b4); h7 = VMACHI(h7, a3, b3); 
  h7 = VMACHI(h7, a4, b2); h7 = VMACHI(h7, a5, b1); h7 = VADD(h7, h7);

  z7 = VMACLO(h7, a2, b5); z7 = VMACLO(z7, a3, b4); z7 = VMACLO(z7, a4, b3); 
  z7 = VMACLO(z7, a5, b2);
  h8 = VMACHI(h8, a2, b5); h8 = VMACHI(h8, a3, b4); h8 = VMACHI(h8, a4, b3); 
  h8 = VMACHI(h8, a5, b2); h8 = VADD(h8, h8);

  z8 = VMACLO(h8, a3, b5); z8 = VMACLO(z8, a4, b4); z8 = VMACLO(z8, a5, b3);
  h9 = VMACHI(h9, a3, b5); h9 = VMACHI(h9, a4, b4); h9 = VMACHI(h9, a5, b3);
  h9 = VADD(h9, h9); 

  z9 = VMACLO(h9, a4, b5); z9 = VMACLO(z9, a5, b4);
  h10 = VMACHI(h10, a4, b5); h10 = VMACHI(h10, a5, b4); h10 = VADD(h10, h10); 

  z10 = VMACLO(h10, a5, b5);
  h11 = VMACHI(h11, a5, b5); h11 = VADD(h11, h11);

  z11 = h11;

  // ---------------------------------------------------------------------------
  // compute zH(z12~z23) by aH(a6~a11) * bH(b6~b11)
  
  z12 = VMACLO(h12, a6, b6); 
  h13 = VMACHI(h13, a6, b6); h13 = VADD(h13, h13);

  z13 = VMACLO(h13, a6, b7); z13 = VMACLO(z13, a7, b6); 
  h14 = VMACHI(h14, a6, b7); h14 = VMACHI(h14, a7, b6); h14 = VADD(h14, h14);

  z14 = VMACLO(h14, a6, b8); z14 = VMACLO(z14, a7, b7); z14 = VMACLO(z14, a8, b6); 
  h15 = VMACHI(h15, a6, b8); h15 = VMACHI(h15, a7, b7); h15 = VMACHI(h15, a8, b6); 
  h15 = VADD(h15, h15);

  z15 = VMACLO(h15, a6, b9); z15 = VMACLO(z15, a7, b8); z15 = VMACLO(z15, a8, b7); 
  z15 = VMACLO(z15, a9, b6);
  h16 = VMACHI(h16, a6, b9); h16 = VMACHI(h16, a7, b8); h16 = VMACHI(h16, a8, b7); 
  h16 = VMACHI(h16, a9, b6); h16 = VADD(h16, h16);

  z16 = VMACLO(h16, a6, b10); z16 = VMACLO(z16, a7,  b9); z16 = VMACLO(z16, a8, b8); 
  z16 = VMACLO(z16, a9, b7);  z16 = VMACLO(z16, a10, b6);
  h17 = VMACHI(h17, a6, b10); h17 = VMACHI(h17, a7,  b9); h17 = VMACHI(h17, a8, b8); 
  h17 = VMACHI(h17, a9, b7);  h17 = VMACHI(h17, a10, b6); h17 = VADD(h17, h17);

  z17 = VMACLO(h17, a6, b11); z17 = VMACLO(z17, a7,  b10); z17 = VMACLO(z17, a8,  b9); 
  z17 = VMACLO(z17, a9, b8);  z17 = VMACLO(z17, a10, b7);  z17 = VMACLO(z17, a11, b6);
  h18 = VMACHI(h18, a6, b11); h18 = VMACHI(h18, a7,  b10); h18 = VMACHI(h18, a8,  b9); 
  h18 = VMACHI(h18, a9, b8);  h18 = VMACHI(h18, a10, b7);  h18 = VMACHI(h18, a11, b6);
  h18 = VADD(h18, h18);

  z18 = VMACLO(h18, a7, b11); z18 = VMACLO(z18, a8, b10); z18 = VMACLO(z18, a9, b9);
  z18 = VMACLO(z18, a10, b8); z18 = VMACLO(z18, a11, b7);
  h19 = VMACHI(h19, a7, b11); h19 = VMACHI(h19, a8, b10); h19 = VMACHI(h19, a9, b9);
  h19 = VMACHI(h19, a10, b8); h19 = VMACHI(h19, a11, b7); h19 = VADD(h19, h19);

  z19 = VMACLO(h19, a8, b11); z19 = VMACLO(z19, a9, b10); z19 = VMACLO(z19, a10, b9);
  z19 = VMACLO(z19, a11, b8);
  h20 = VMACHI(h20, a8, b11); h20 = VMACHI(h20, a9, b10); h20 = VMACHI(h20, a10, b9);
  h20 = VMACHI(h20, a11, b8); h20 = VADD(h20, h20);

  z20 = VMACLO(h20, a9, b11); z20 = VMACLO(z20, a10, b10); z20 = VMACLO(z20, a11, b9);
  h21 = VMACHI(h21, a9, b11); h21 = VMACHI(h21, a10, b10); h21 = VMACHI(h21, a11, b9);
  h21 = VADD(h21, h21);

  z21 = VMACLO(h21, a10, b11); z21 = VMACLO(z21, a11, b10);
  h22 = VMACHI(h22, a10, b11); h22 = VMACHI(h22, a11, b10); h22 = VADD(h22, h22);

  z22 = VMACLO(h22, a11, b11); 
  h23 = VMACHI(h23, a11, b11); h23 = VADD(h23, h23);

  z23 = h23;

  // ---------------------------------------------------------------------------
  // ta(ta0~ta5) = aL(a0~a5) + aH(a6~a11) 
  ta0 = VADD(a0, a6); ta1 = VADD(a1, a7);  ta2 = VADD(a2, a8);
  ta3 = VADD(a3, a9); ta4 = VADD(a4, a10); ta5 = VADD(a5, a11);

  // tb(tb0~tb5) = bL(b0~b5) + bH(b6~b11) 
  tb0 = VADD(b0, b6); tb1 = VADD(b1, b7);  tb2 = VADD(b2, b8);
  tb3 = VADD(b3, b9); tb4 = VADD(b4, b10); tb5 = VADD(b5, b11);

  // ---------------------------------------------------------------------------
  // zM = ta * tb - zL - zH

  h0  = h1  = h2  = h3  = h4  = h5  = VZERO;
  h6  = h7  = h8  = h9  = h10 = h11 = VZERO;

  m0 = VMACLO(m0, ta0, tb0);
  h1 = VMACHI(h1, ta0, tb0); h1 = VADD(h1, h1);

  m1 = VMACLO(h1, ta0, tb1); m1 = VMACLO(m1, ta1, tb0); 
  h2 = VMACHI(h2, ta0, tb1); h2 = VMACHI(h2, ta1, tb0); h2 = VADD(h2, h2);

  m2 = VMACLO(h2, ta0, tb2); m2 = VMACLO(m2, ta1, tb1); m2 = VMACLO(m2, ta2, tb0);
  h3 = VMACHI(h3, ta0, tb2); h3 = VMACHI(h3, ta1, tb1); h3 = VMACHI(h3, ta2, tb0);
  h3 = VADD(h3, h3);

  m3 = VMACLO(h3, ta0, tb3); m3 = VMACLO(m3, ta1, tb2); m3 = VMACLO(m3, ta2, tb1); 
  m3 = VMACLO(m3, ta3, tb0);
  h4 = VMACHI(h4, ta0, tb3); h4 = VMACHI(h4, ta1, tb2); h4 = VMACHI(h4, ta2, tb1); 
  h4 = VMACHI(h4, ta3, tb0); h4 = VADD(h4, h4);

  m4 = VMACLO(h4, ta0, tb4); m4 = VMACLO(m4, ta1, tb3); m4 = VMACLO(m4, ta2, tb2); 
  m4 = VMACLO(m4, ta3, tb1); m4 = VMACLO(m4, ta4, tb0);
  h5 = VMACHI(h5, ta0, tb4); h5 = VMACHI(h5, ta1, tb3); h5 = VMACHI(h5, ta2, tb2); 
  h5 = VMACHI(h5, ta3, tb1); h5 = VMACHI(h5, ta4, tb0); h5 = VADD(h5, h5);

  m5 = VMACLO(h5, ta0, tb5); m5 = VMACLO(m5, ta1, tb4); m5 = VMACLO(m5, ta2, tb3); 
  m5 = VMACLO(m5, ta3, tb2); m5 = VMACLO(m5, ta4, tb1); m5 = VMACLO(m5, ta5, tb0);
  h6 = VMACHI(h6, ta0, tb5); h6 = VMACHI(h6, ta1, tb4); h6 = VMACHI(h6, ta2, tb3); 
  h6 = VMACHI(h6, ta3, tb2); h6 = VMACHI(h6, ta4, tb1); h6 = VMACHI(h6, ta5, tb0);
  h6 = VADD(h6, h6);

  m6 = VMACLO(h6, ta1, tb5); m6 = VMACLO(m6, ta2, tb4); m6 = VMACLO(m6, ta3, tb3); 
  m6 = VMACLO(m6, ta4, tb2); m6 = VMACLO(m6, ta5, tb1);
  h7 = VMACHI(h7, ta1, tb5); h7 = VMACHI(h7, ta2, tb4); h7 = VMACHI(h7, ta3, tb3); 
  h7 = VMACHI(h7, ta4, tb2); h7 = VMACHI(h7, ta5, tb1); h7 = VADD(h7, h7);

  m7 = VMACLO(h7, ta2, tb5); m7 = VMACLO(m7, ta3, tb4); m7 = VMACLO(m7, ta4, tb3); 
  m7 = VMACLO(m7, ta5, tb2);
  h8 = VMACHI(h8, ta2, tb5); h8 = VMACHI(h8, ta3, tb4); h8 = VMACHI(h8, ta4, tb3); 
  h8 = VMACHI(h8, ta5, tb2); h8 = VADD(h8, h8);

  m8 = VMACLO(h8, ta3, tb5); m8 = VMACLO(m8, ta4, tb4); m8 = VMACLO(m8, ta5, tb3);
  h9 = VMACHI(h9, ta3, tb5); h9 = VMACHI(h9, ta4, tb4); h9 = VMACHI(h9, ta5, tb3);
  h9 = VADD(h9, h9); 

  m9 = VMACLO(h9, ta4, tb5); m9 = VMACLO(m9, ta5, tb4);
  h10 = VMACHI(h10, ta4, tb5); h10 = VMACHI(h10, ta5, tb4); h10 = VADD(h10, h10); 

  m10 = VMACLO(h10, ta5, tb5);
  h11 = VMACHI(h11, ta5, tb5); h11 = VADD(h11, h11);

  m11 = h11;

  m0  = VSUB(m0,  VADD(z0,  z12)); m1  = VSUB(m1,  VADD(z1, z13));
  m2  = VSUB(m2,  VADD(z2,  z14)); m3  = VSUB(m3,  VADD(z3, z15));
  m4  = VSUB(m4,  VADD(z4,  z16)); m5  = VSUB(m5,  VADD(z5, z17));
  m6  = VSUB(m6,  VADD(z6,  z18)); m7  = VSUB(m7,  VADD(z7, z19));
  m8  = VSUB(m8,  VADD(z8,  z20)); m9  = VSUB(m9,  VADD(z9, z21));
  m10 = VSUB(m10, VADD(z10, z22)); m11 = VSUB(m11, VADD(z11, z23));

  // z = z + zM
  z6  = VADD(z6,  m0);  z7  = VADD(z7,  m1); 
  z8  = VADD(z8,  m2);  z9  = VADD(z9,  m3); 
  z10 = VADD(z10, m4);  z11 = VADD(z11, m5); 
  z12 = VADD(z12, m6);  z13 = VADD(z13, m7); 
  z14 = VADD(z14, m8);  z15 = VADD(z15, m9); 
  z16 = VADD(z16, m10); z17 = VADD(z17, m11); 

  // ---------------------------------------------------------------------------
  r[0]  = z0;  r[1]  = z1;  r[2]  = z2;  r[3]  = z3;  r[4]  = z4;  
  r[5]  = z5;  r[6]  = z6;  r[7]  = z7;  r[8]  = z8;  r[9]  = z9;  
  r[10] = z10; r[11] = z11; r[12] = z12; r[13] = z13; r[14] = z14; 
  r[15] = z15; r[16] = z16; r[17] = z17; r[18] = z18; r[19] = z19;
  r[20] = z20; r[21] = z21; r[22] = z22; r[23] = z23;
}

// Montgomery reduction r = a * R^-1 mod 2p, where R = 2^468 (operand-scanning)
void rdc_mont_8x1w(__m512i *r, const __m512i *a)
{
  __m512i a0  = a[0],  a1  = a[1],  a2  = a[2],  a3  = a[3];  
  __m512i a4  = a[4],  a5  = a[5],  a6  = a[6],  a7  = a[7];
  __m512i a8  = a[8],  a9  = a[9],  a10 = a[10], a11 = a[11];
  __m512i a12 = a[12], a13 = a[13], a14 = a[14], a15 = a[15];
  __m512i a16 = a[16], a17 = a[17], a18 = a[18], a19 = a[19];
  __m512i a20 = a[20], a21 = a[21], a22 = a[22], a23 = a[23];
  __m512i h6  = VZERO, h7  = VZERO, h8  = VZERO, h9  = VZERO;
  __m512i h10 = VZERO, h11 = VZERO, h12 = VZERO, h13 = VZERO; 
  __m512i h14 = VZERO, h15 = VZERO, h16 = VZERO, h17 = VZERO;
  __m512i h18 = VZERO, h19 = VZERO, h20 = VZERO, h21 = VZERO;
  __m512i h22 = VZERO, h23 = VZERO, u; 
  const __m512i vp5  = VSET1(vp610p1[5]),  vp6 = VSET1(vp610p1[6]), vp7  = VSET1(vp610p1[7]);
  const __m512i vp8  = VSET1(vp610p1[8]),  vp9 = VSET1(vp610p1[9]), vp10 = VSET1(vp610p1[10]);
  const __m512i vp11 = VSET1(vp610p1[11]), vbmask = VSET1(VBMASK);

  // some pre carry propagations
  a1 = VADD(a1, VSRA(a0, VBRADIX)); a0 = VAND(a0, vbmask);
  a2 = VADD(a2, VSRA(a1, VBRADIX)); a1 = VAND(a1, vbmask);
  a3 = VADD(a3, VSRA(a2, VBRADIX)); a2 = VAND(a2, vbmask);
  a4 = VADD(a4, VSRA(a3, VBRADIX)); a3 = VAND(a3, vbmask);
  a5 = VADD(a5, VSRA(a4, VBRADIX)); a4 = VAND(a4, vbmask);

  u = a0;
  a5  = VMACLO(a5,  u, vp5); a6  = VMACLO(a6,  u, vp6); a7  = VMACLO(a7,  u, vp7); 
  a8  = VMACLO(a8,  u, vp8); a9  = VMACLO(a9,  u, vp9); a10 = VMACLO(a10, u, vp10);
  a11 = VMACLO(a11, u, vp11);
  h6  = VMACHI(h6,  u, vp5); h7  = VMACHI(h7,  u, vp6); h8  = VMACHI(h8,  u, vp7); 
  h9  = VMACHI(h9,  u, vp8); h10 = VMACHI(h10, u, vp9); h11 = VMACHI(h11, u, vp10);
  h12 = VMACHI(h12, u, vp11);
  a6 = VADD(a6, VSRA(a5, VBRADIX)); a5 = VAND(a5, vbmask); 
  a6 = VADD(a6, VADD(h6, h6));

  u = a1; 
  a6  = VMACLO(a6,  u, vp5); a7  = VMACLO(a7,  u, vp6); a8  = VMACLO(a8,  u, vp7); 
  a9  = VMACLO(a9,  u, vp8); a10 = VMACLO(a10, u, vp9); a11 = VMACLO(a11, u, vp10);
  a12 = VMACLO(a12, u, vp11);
  h7  = VMACHI(h7,  u, vp5); h8  = VMACHI(h8,  u, vp6); h9  = VMACHI(h9,  u, vp7); 
  h10 = VMACHI(h10, u, vp8); h11 = VMACHI(h11, u, vp9); h12 = VMACHI(h12, u, vp10);
  h13 = VMACHI(h13, u, vp11);
  a7 = VADD(a7, VSRA(a6, VBRADIX)); a6 = VAND(a6, vbmask); 
  a7 = VADD(a7, VADD(h7, h7));

  u = a2; 
  a7  = VMACLO(a7,  u, vp5); a8  = VMACLO(a8,  u, vp6); a9  = VMACLO(a9,  u, vp7); 
  a10 = VMACLO(a10, u, vp8); a11 = VMACLO(a11, u, vp9); a12 = VMACLO(a12, u, vp10);
  a13 = VMACLO(a13, u, vp11);
  h8  = VMACHI(h8,  u, vp5); h9  = VMACHI(h9,  u, vp6); h10 = VMACHI(h10, u, vp7); 
  h11 = VMACHI(h11, u, vp8); h12 = VMACHI(h12, u, vp9); h13 = VMACHI(h13, u, vp10);
  h14 = VMACHI(h14, u, vp11);
  a8 = VADD(a8, VSRA(a7, VBRADIX)); a7 = VAND(a7, vbmask); 
  a8 = VADD(a8, VADD(h8, h8));

  u = a3; 
  a8  = VMACLO(a8,  u, vp5); a9  = VMACLO(a9,  u, vp6); a10 = VMACLO(a10, u, vp7); 
  a11 = VMACLO(a11, u, vp8); a12 = VMACLO(a12, u, vp9); a13 = VMACLO(a13, u, vp10);
  a14 = VMACLO(a14, u, vp11);
  h9  = VMACHI(h9,  u, vp5); h10 = VMACHI(h10, u, vp6); h11 = VMACHI(h11, u, vp7); 
  h12 = VMACHI(h12, u, vp8); h13 = VMACHI(h13, u, vp9); h14 = VMACHI(h14, u, vp10);
  h15 = VMACHI(h15, u, vp11);
  a9 = VADD(a9, VSRA(a8, VBRADIX)); a8 = VAND(a8, vbmask); 
  a9 = VADD(a9, VADD(h9, h9));

  u = a4; 
  a9  = VMACLO(a9,  u, vp5); a10 = VMACLO(a10, u, vp6); a11 = VMACLO(a11, u, vp7); 
  a12 = VMACLO(a12, u, vp8); a13 = VMACLO(a13, u, vp9); a14 = VMACLO(a14, u, vp10);
  a15 = VMACLO(a15, u, vp11);
  h10 = VMACHI(h10, u, vp5); h11 = VMACHI(h11, u, vp6); h12 = VMACHI(h12, u, vp7); 
  h13 = VMACHI(h13, u, vp8); h14 = VMACHI(h14, u, vp9); h15 = VMACHI(h15, u, vp10);
  h16 = VMACHI(h16, u, vp11);
  a10 = VADD(a10, VSRA(a9, VBRADIX)); a9 = VAND(a9, vbmask); 
  a10 = VADD(a10, VADD(h10, h10));

  u = a5; 
  a10 = VMACLO(a10, u, vp5); a11 = VMACLO(a11, u, vp6); a12 = VMACLO(a12, u, vp7); 
  a13 = VMACLO(a13, u, vp8); a14 = VMACLO(a14, u, vp9); a15 = VMACLO(a15, u, vp10);
  a16 = VMACLO(a16, u, vp11);
  h11 = VMACHI(h11, u, vp5); h12 = VMACHI(h12, u, vp6); h13 = VMACHI(h13, u, vp7); 
  h14 = VMACHI(h14, u, vp8); h15 = VMACHI(h15, u, vp9); h16 = VMACHI(h16, u, vp10);
  h17 = VMACHI(h17, u, vp11);
  a11 = VADD(a11, VSRA(a10, VBRADIX)); a10 = VAND(a10, vbmask); 
  a11 = VADD(a11, VADD(h11, h11));

  u = a6; 
  a11 = VMACLO(a11, u, vp5); a12 = VMACLO(a12, u, vp6); a13 = VMACLO(a13, u, vp7); 
  a14 = VMACLO(a14, u, vp8); a15 = VMACLO(a15, u, vp9); a16 = VMACLO(a16, u, vp10);
  a17 = VMACLO(a17, u, vp11);
  h12 = VMACHI(h12, u, vp5); h13 = VMACHI(h13, u, vp6); h14 = VMACHI(h14, u, vp7); 
  h15 = VMACHI(h15, u, vp8); h16 = VMACHI(h16, u, vp9); h17 = VMACHI(h17, u, vp10);
  h18 = VMACHI(h18, u, vp11);
  a12 = VADD(a12, VSRA(a11, VBRADIX)); a11 = VAND(a11, vbmask); 
  a12 = VADD(a12, VADD(h12, h12));

  u = a7; 
  a12 = VMACLO(a12, u, vp5); a13 = VMACLO(a13, u, vp6); a14 = VMACLO(a14, u, vp7); 
  a15 = VMACLO(a15, u, vp8); a16 = VMACLO(a16, u, vp9); a17 = VMACLO(a17, u, vp10);
  a18 = VMACLO(a18, u, vp11);
  h13 = VMACHI(h13, u, vp5); h14 = VMACHI(h14, u, vp6); h15 = VMACHI(h15, u, vp7); 
  h16 = VMACHI(h16, u, vp8); h17 = VMACHI(h17, u, vp9); h18 = VMACHI(h18, u, vp10);
  h19 = VMACHI(h19, u, vp11);
  a13 = VADD(a13, VSRA(a12, VBRADIX)); a12 = VAND(a12, vbmask); 
  a13 = VADD(a13, VADD(h13, h13));

  u = a8; 
  a13 = VMACLO(a13, u, vp5); a14 = VMACLO(a14, u, vp6); a15 = VMACLO(a15, u, vp7); 
  a16 = VMACLO(a16, u, vp8); a17 = VMACLO(a17, u, vp9); a18 = VMACLO(a18, u, vp10);
  a19 = VMACLO(a19, u, vp11);
  h14 = VMACHI(h14, u, vp5); h15 = VMACHI(h15, u, vp6); h16 = VMACHI(h16, u, vp7); 
  h17 = VMACHI(h17, u, vp8); h18 = VMACHI(h18, u, vp9); h19 = VMACHI(h19, u, vp10);
  h20 = VMACHI(h20, u, vp11);
  a14 = VADD(a14, VSRA(a13, VBRADIX)); a13 = VAND(a13, vbmask);
  a14 = VADD(a14, VADD(h14, h14));

  u = a9;
  a14 = VMACLO(a14, u, vp5); a15 = VMACLO(a15, u, vp6); a16 = VMACLO(a16, u, vp7); 
  a17 = VMACLO(a17, u, vp8); a18 = VMACLO(a18, u, vp9); a19 = VMACLO(a19, u, vp10);
  a20 = VMACLO(a20, u, vp11);
  h15 = VMACHI(h15, u, vp5); h16 = VMACHI(h16, u, vp6); h17 = VMACHI(h17, u, vp7); 
  h18 = VMACHI(h18, u, vp8); h19 = VMACHI(h19, u, vp9); h20 = VMACHI(h20, u, vp10);
  h21 = VMACHI(h21, u, vp11);
  a15 = VADD(a15, VSRA(a14, VBRADIX)); a14 = VAND(a14, vbmask);
  a15 = VADD(a15, VADD(h15, h15));

  u = a10;
  a15 = VMACLO(a15, u, vp5); a16 = VMACLO(a16, u, vp6); a17 = VMACLO(a17, u, vp7); 
  a18 = VMACLO(a18, u, vp8); a19 = VMACLO(a19, u, vp9); a20 = VMACLO(a20, u, vp10);
  a21 = VMACLO(a21, u, vp11);
  h16 = VMACHI(h16, u, vp5); h17 = VMACHI(h17, u, vp6); h18 = VMACHI(h18, u, vp7); 
  h19 = VMACHI(h19, u, vp8); h20 = VMACHI(h20, u, vp9); h21 = VMACHI(h21, u, vp10);
  h22 = VMACHI(h22, u, vp11);
  a16 = VADD(a16, VSRA(a15, VBRADIX)); a15 = VAND(a15, vbmask);
  a16 = VADD(a16, VADD(h16, h16));

  u = a11;
  a16 = VMACLO(a16, u, vp5); a17 = VMACLO(a17, u, vp6); a18 = VMACLO(a18, u, vp7); 
  a19 = VMACLO(a19, u, vp8); a20 = VMACLO(a20, u, vp9); a21 = VMACLO(a21, u, vp10);
  a22 = VMACLO(a22, u, vp11);
  h17 = VMACHI(h17, u, vp5); h18 = VMACHI(h18, u, vp6); h19 = VMACHI(h19, u, vp7); 
  h20 = VMACHI(h20, u, vp8); h21 = VMACHI(h21, u, vp9); h22 = VMACHI(h22, u, vp10);
  h23 = VMACHI(h23, u, vp11);
  a17 = VADD(a17, VSRA(a16, VBRADIX)); a16 = VAND(a16, vbmask);
  a17 = VADD(a17, VADD(h17, h17));

  a18 = VADD(a18, VADD(h18, h18));
  a19 = VADD(a19, VADD(h19, h19));
  a20 = VADD(a20, VADD(h20, h20));
  a21 = VADD(a21, VADD(h21, h21));
  a22 = VADD(a22, VADD(h22, h22));
  a23 = VADD(a23, VADD(h23, h23));

  a18 = VADD(a18, VSRA(a17, VBRADIX)); a17 = VAND(a17, vbmask);
  a19 = VADD(a19, VSRA(a18, VBRADIX)); a18 = VAND(a18, vbmask);
  a20 = VADD(a20, VSRA(a19, VBRADIX)); a19 = VAND(a19, vbmask);
  a21 = VADD(a21, VSRA(a20, VBRADIX)); a20 = VAND(a20, vbmask);
  a22 = VADD(a22, VSRA(a21, VBRADIX)); a21 = VAND(a21, vbmask);
  a23 = VADD(a23, VSRA(a22, VBRADIX)); a22 = VAND(a22, vbmask);
  
  r[0] = a12; r[1] = a13; r[2]  = a14;  r[3] = a15; 
  r[4] = a16; r[5] = a17; r[6]  = a18;  r[7] = a19; 
  r[8] = a20; r[9] = a21; r[10] = a22; r[11] = a23;
}

// -----------------------------------------------------------------------------
// (4x2)-way Fp arithmetic 

// integer subtraction r = a - b + 2p
void mp_sub_p2_4x2w(__m512i *r, const __m512i *a, const __m512i *b)
{
  __m512i a0 = a[0], a1 = a[1], a2 = a[2], a3 = a[3], a4 = a[4], a5 = a[5];
  __m512i b0 = b[0], b1 = b[1], b2 = b[2], b3 = b[3], b4 = b[4], b5 = b[5];
  __m512i r0, r1, r2, r3, r4, r5, c;
  const __m512i vp0 = VSET(vp610x2[6],  vp610x2[0], vp610x2[6],  vp610x2[0], vp610x2[6],  vp610x2[0], vp610x2[6],  vp610x2[0]);
  const __m512i vp1 = VSET(vp610x2[7],  vp610x2[1], vp610x2[7],  vp610x2[1], vp610x2[7],  vp610x2[1], vp610x2[7],  vp610x2[1]);
  const __m512i vp2 = VSET(vp610x2[8],  vp610x2[2], vp610x2[8],  vp610x2[2], vp610x2[8],  vp610x2[2], vp610x2[8],  vp610x2[2]);
  const __m512i vp3 = VSET(vp610x2[9],  vp610x2[3], vp610x2[9],  vp610x2[3], vp610x2[9],  vp610x2[3], vp610x2[9],  vp610x2[3]);
  const __m512i vp4 = VSET(vp610x2[10], vp610x2[4], vp610x2[10], vp610x2[4], vp610x2[10], vp610x2[4], vp610x2[10], vp610x2[4]);
  const __m512i vp5 = VSET(vp610x2[11], vp610x2[5], vp610x2[11], vp610x2[5], vp610x2[11], vp610x2[5], vp610x2[11], vp610x2[5]);
  const __m512i vbmask = VSET1(VBMASK); 

  // r = a + 2p
  r0 = VADD(a0, vp0); r1 = VADD(a1, vp1); r2 = VADD(a2, vp2); 
  r3 = VADD(a3, vp3); r4 = VADD(a4, vp4); r5 = VADD(a5, vp5); 

  // r = r - b
  r0 = VSUB(r0, b0); r1 = VSUB(r1, b1); r2 = VSUB(r2, b2); 
  r3 = VSUB(r3, b3); r4 = VSUB(r4, b4); r5 = VSUB(r5, b5); 

  // *simple* carry propagation
  // some limbs are finally 52-bit not 51-bit 

  c = VSRA(r0, VBRADIX); r0 = VAND(r0, vbmask); r1 = VADD(r1, c);
  c = VSRA(r1, VBRADIX); r1 = VAND(r1, vbmask); r2 = VADD(r2, c);
  c = VSRA(r2, VBRADIX); r2 = VAND(r2, vbmask); r3 = VADD(r3, c);
  c = VSRA(r3, VBRADIX); r3 = VAND(r3, vbmask); r4 = VADD(r4, c);
  c = VSRA(r4, VBRADIX); r4 = VAND(r4, vbmask); r5 = VADD(r5, c);
  c = VSRA(r5, VBRADIX); r5 = VAND(r5, vbmask); 
  r0 = VMADD(r0, 0xAA, r0, VSHUF(c, 0x4E));

  r[0] = r0; r[1] = r1; r[2] = r2; r[3] = r3; r[4] = r4; r[5] = r5; 
}

// integer subtraction r = a - b + 4p
void mp_sub_p4_4x2w(__m512i *r, const __m512i *a, const __m512i *b)
{
  __m512i a0 = a[0], a1 = a[1], a2 = a[2], a3 = a[3], a4 = a[4], a5 = a[5];
  __m512i b0 = b[0], b1 = b[1], b2 = b[2], b3 = b[3], b4 = b[4], b5 = b[5];
  __m512i r0, r1, r2, r3, r4, r5, c;
  const __m512i vp0 = VSET(vp610x4[6],  vp610x4[0], vp610x4[6],  vp610x4[0], vp610x4[6],  vp610x4[0], vp610x4[6],  vp610x4[0]);
  const __m512i vp1 = VSET(vp610x4[7],  vp610x4[1], vp610x4[7],  vp610x4[1], vp610x4[7],  vp610x4[1], vp610x4[7],  vp610x4[1]);
  const __m512i vp2 = VSET(vp610x4[8],  vp610x4[2], vp610x4[8],  vp610x4[2], vp610x4[8],  vp610x4[2], vp610x4[8],  vp610x4[2]);
  const __m512i vp3 = VSET(vp610x4[9],  vp610x4[3], vp610x4[9],  vp610x4[3], vp610x4[9],  vp610x4[3], vp610x4[9],  vp610x4[3]);
  const __m512i vp4 = VSET(vp610x4[10], vp610x4[4], vp610x4[10], vp610x4[4], vp610x4[10], vp610x4[4], vp610x4[10], vp610x4[4]);
  const __m512i vp5 = VSET(vp610x4[11], vp610x4[5], vp610x4[11], vp610x4[5], vp610x4[11], vp610x4[5], vp610x4[11], vp610x4[5]);

  const __m512i vbmask = VSET1(VBMASK); 

  // r = a + 4p
  r0 = VADD(a0, vp0); r1 = VADD(a1, vp1); r2 = VADD(a2, vp2); 
  r3 = VADD(a3, vp3); r4 = VADD(a4, vp4); r5 = VADD(a5, vp5); 

  // r = r - b
  r0 = VSUB(r0, b0); r1 = VSUB(r1, b1); r2 = VSUB(r2, b2); 
  r3 = VSUB(r3, b3); r4 = VSUB(r4, b4); r5 = VSUB(r5, b5); 

  // *simple* carry propagation
  // some limbs are finally 52-bit not 51-bit 

  c = VSRA(r0, VBRADIX); r0 = VAND(r0, vbmask); r1 = VADD(r1, c);
  c = VSRA(r1, VBRADIX); r1 = VAND(r1, vbmask); r2 = VADD(r2, c);
  c = VSRA(r2, VBRADIX); r2 = VAND(r2, vbmask); r3 = VADD(r3, c);
  c = VSRA(r3, VBRADIX); r3 = VAND(r3, vbmask); r4 = VADD(r4, c);
  c = VSRA(r4, VBRADIX); r4 = VAND(r4, vbmask); r5 = VADD(r5, c);
  c = VSRA(r5, VBRADIX); r5 = VAND(r5, vbmask); 
  r0 = VMADD(r0, 0xAA, r0, VSHUF(c, 0x4E));

  r[0] = r0; r[1] = r1; r[2] = r2; r[3] = r3; r[4] = r4; r[5] = r5; 
}

// modular addition r = a + b mod 2p
void fpadd_4x2w(__m512i *r, const __m512i *a, const __m512i* b)
{
  __m512i a0 = a[0], a1 = a[1], a2 = a[2], a3 = a[3], a4 = a[4], a5 = a[5];
  __m512i b0 = b[0], b1 = b[1], b2 = b[2], b3 = b[3], b4 = b[4], b5 = b[5];
  __m512i r0, r1, r2, r3, r4, r5, c, carry, smask;
  const __m512i vp0 = VSET(vp610x2[6],  vp610x2[0], vp610x2[6],  vp610x2[0], vp610x2[6],  vp610x2[0], vp610x2[6],  vp610x2[0]);
  const __m512i vp1 = VSET(vp610x2[7],  vp610x2[1], vp610x2[7],  vp610x2[1], vp610x2[7],  vp610x2[1], vp610x2[7],  vp610x2[1]);
  const __m512i vp2 = VSET(vp610x2[8],  vp610x2[2], vp610x2[8],  vp610x2[2], vp610x2[8],  vp610x2[2], vp610x2[8],  vp610x2[2]);
  const __m512i vp3 = VSET(vp610x2[9],  vp610x2[3], vp610x2[9],  vp610x2[3], vp610x2[9],  vp610x2[3], vp610x2[9],  vp610x2[3]);
  const __m512i vp4 = VSET(vp610x2[10], vp610x2[4], vp610x2[10], vp610x2[4], vp610x2[10], vp610x2[4], vp610x2[10], vp610x2[4]);
  const __m512i vp5 = VSET(vp610x2[11], vp610x2[5], vp610x2[11], vp610x2[5], vp610x2[11], vp610x2[5], vp610x2[11], vp610x2[5]);
  const __m512i vbmask = VSET1(VBMASK); 

  // r = a + b
  r0 = VADD(a0, b0); r1 = VADD(a1, b1); r2 = VADD(a2, b2); 
  r3 = VADD(a3, b3); r4 = VADD(a4, b4); r5 = VADD(a5, b5); 

  // r = r - 2p
  r0 = VSUB(r0, vp0); r1 = VSUB(r1, vp1); r2 = VSUB(r2, vp2); 
  r3 = VSUB(r3, vp3); r4 = VSUB(r4, vp4); r5 = VSUB(r5, vp5); 

  // carry propagation 
  carry = VSRA(r0, VBRADIX); r0 = VMAND(r0, 0x55, r0, vbmask); r1 = VMADD(r1, 0x55, r1, carry);
  carry = VSRA(r1, VBRADIX); r1 = VMAND(r1, 0x55, r1, vbmask); r2 = VMADD(r2, 0x55, r2, carry);
  carry = VSRA(r2, VBRADIX); r2 = VMAND(r2, 0x55, r2, vbmask); r3 = VMADD(r3, 0x55, r3, carry);
  carry = VSRA(r3, VBRADIX); r3 = VMAND(r3, 0x55, r3, vbmask); r4 = VMADD(r4, 0x55, r4, carry);
  carry = VSRA(r4, VBRADIX); r4 = VMAND(r4, 0x55, r4, vbmask); r5 = VMADD(r5, 0x55, r5, carry);
  carry = VSRA(r5, VBRADIX); r5 = VMAND(r5, 0x55, r5, vbmask); 
  carry = VZSHUF(0xCCCC, carry, 0x4E); r0 = VMADD(r0, 0xAA, r0, carry); 
  carry = VSRA(r0, VBRADIX); r0 = VMAND(r0, 0xAA, r0, vbmask); r1 = VMADD(r1, 0xAA, r1, carry);
  carry = VSRA(r1, VBRADIX); r1 = VMAND(r1, 0xAA, r1, vbmask); r2 = VMADD(r2, 0xAA, r2, carry);
  carry = VSRA(r2, VBRADIX); r2 = VMAND(r2, 0xAA, r2, vbmask); r3 = VMADD(r3, 0xAA, r3, carry);
  carry = VSRA(r3, VBRADIX); r3 = VMAND(r3, 0xAA, r3, vbmask); r4 = VMADD(r4, 0xAA, r4, carry);
  carry = VSRA(r4, VBRADIX); r4 = VMAND(r4, 0xAA, r4, vbmask); r5 = VMADD(r5, 0xAA, r5, carry);

  // if r is non-negative, smask = all-0 
  // if r is     negative, smask = all-1
  smask = VSRA(r5, 63);                 // smask = all-0/-1 X | all-0/-1 X | all-0/-1 X | all-0/-1 X | 
  smask = VSHUF(smask, 0xEE);           // smask = all-0/-1 | all-0/-1 | all-0/-1 | all-0/-1
  // r = r + (2p & smask)
  r0 = VADD(r0, VAND(vp0, smask)); r1 = VADD(r1, VAND(vp1, smask));
  r2 = VADD(r2, VAND(vp2, smask)); r3 = VADD(r3, VAND(vp3, smask));
  r4 = VADD(r4, VAND(vp4, smask)); r5 = VADD(r5, VAND(vp5, smask)); 

  // *simple* carry propagation
  // some limbs are finally 52-bit not 51-bit 

  c = VSRA(r0, VBRADIX); r0 = VAND(r0, vbmask); r1 = VADD(r1, c);
  c = VSRA(r1, VBRADIX); r1 = VAND(r1, vbmask); r2 = VADD(r2, c);
  c = VSRA(r2, VBRADIX); r2 = VAND(r2, vbmask); r3 = VADD(r3, c);
  c = VSRA(r3, VBRADIX); r3 = VAND(r3, vbmask); r4 = VADD(r4, c);
  c = VSRA(r4, VBRADIX); r4 = VAND(r4, vbmask); r5 = VADD(r5, c);
  c = VSRA(r5, VBRADIX); r5 = VAND(r5, vbmask); 
  r0 = VMADD(r0, 0xAA, r0, VSHUF(c, 0x4E));

  r[0] = r0; r[1] = r1; r[2] = r2; r[3] = r3; r[4] = r4; r[5] = r5; 
}

// modular subtraction r = a + b mod 2p
void fpsub_4x2w(__m512i *r, const __m512i *a, const __m512i* b)
{
  __m512i a0 = a[0], a1 = a[1], a2 = a[2], a3 = a[3], a4 = a[4], a5 = a[5];
  __m512i b0 = b[0], b1 = b[1], b2 = b[2], b3 = b[3], b4 = b[4], b5 = b[5];
  __m512i r0, r1, r2, r3, r4, r5, c, carry, smask;
  const __m512i vp0 = VSET(vp610x2[6],  vp610x2[0], vp610x2[6],  vp610x2[0], vp610x2[6],  vp610x2[0], vp610x2[6],  vp610x2[0]);
  const __m512i vp1 = VSET(vp610x2[7],  vp610x2[1], vp610x2[7],  vp610x2[1], vp610x2[7],  vp610x2[1], vp610x2[7],  vp610x2[1]);
  const __m512i vp2 = VSET(vp610x2[8],  vp610x2[2], vp610x2[8],  vp610x2[2], vp610x2[8],  vp610x2[2], vp610x2[8],  vp610x2[2]);
  const __m512i vp3 = VSET(vp610x2[9],  vp610x2[3], vp610x2[9],  vp610x2[3], vp610x2[9],  vp610x2[3], vp610x2[9],  vp610x2[3]);
  const __m512i vp4 = VSET(vp610x2[10], vp610x2[4], vp610x2[10], vp610x2[4], vp610x2[10], vp610x2[4], vp610x2[10], vp610x2[4]);
  const __m512i vp5 = VSET(vp610x2[11], vp610x2[5], vp610x2[11], vp610x2[5], vp610x2[11], vp610x2[5], vp610x2[11], vp610x2[5]);
  const __m512i vbmask = VSET1(VBMASK); 

  // r = a - b
  r0 = VSUB(a0, b0); r1 = VSUB(a1, b1); r2 = VSUB(a2, b2); 
  r3 = VSUB(a3, b3); r4 = VSUB(a4, b4); r5 = VSUB(a5, b5); 

  // carry propagation 
  carry = VSRA(r0, VBRADIX); r0 = VMAND(r0, 0x55, r0, vbmask); r1 = VMADD(r1, 0x55, r1, carry);
  carry = VSRA(r1, VBRADIX); r1 = VMAND(r1, 0x55, r1, vbmask); r2 = VMADD(r2, 0x55, r2, carry);
  carry = VSRA(r2, VBRADIX); r2 = VMAND(r2, 0x55, r2, vbmask); r3 = VMADD(r3, 0x55, r3, carry);
  carry = VSRA(r3, VBRADIX); r3 = VMAND(r3, 0x55, r3, vbmask); r4 = VMADD(r4, 0x55, r4, carry);
  carry = VSRA(r4, VBRADIX); r4 = VMAND(r4, 0x55, r4, vbmask); r5 = VMADD(r5, 0x55, r5, carry);
  carry = VSRA(r5, VBRADIX); r5 = VMAND(r5, 0x55, r5, vbmask); 
  carry = VZSHUF(0xCCCC, carry, 0x4E); r0 = VMADD(r0, 0xAA, r0, carry); 
  carry = VSRA(r0, VBRADIX); r0 = VMAND(r0, 0xAA, r0, vbmask); r1 = VMADD(r1, 0xAA, r1, carry);
  carry = VSRA(r1, VBRADIX); r1 = VMAND(r1, 0xAA, r1, vbmask); r2 = VMADD(r2, 0xAA, r2, carry);
  carry = VSRA(r2, VBRADIX); r2 = VMAND(r2, 0xAA, r2, vbmask); r3 = VMADD(r3, 0xAA, r3, carry);
  carry = VSRA(r3, VBRADIX); r3 = VMAND(r3, 0xAA, r3, vbmask); r4 = VMADD(r4, 0xAA, r4, carry);
  carry = VSRA(r4, VBRADIX); r4 = VMAND(r4, 0xAA, r4, vbmask); r5 = VMADD(r5, 0xAA, r5, carry);

  // if r is non-negative, smask = all-0 
  // if r is     negative, smask = all-1
  smask = VSRA(r5, 63);                 // smask = all-0/-1 X | all-0/-1 X | all-0/-1 X | all-0/-1 X | 
  smask = VSHUF(smask, 0xEE);           // smask = all-0/-1 | all-0/-1 | all-0/-1 | all-0/-1
  // r = r + (2p & smask)
  r0 = VADD(r0, VAND(vp0, smask)); r1 = VADD(r1, VAND(vp1, smask));
  r2 = VADD(r2, VAND(vp2, smask)); r3 = VADD(r3, VAND(vp3, smask));
  r4 = VADD(r4, VAND(vp4, smask)); r5 = VADD(r5, VAND(vp5, smask)); 

  // *simple* carry propagation
  // some limbs are finally 52-bit not 51-bit 

  c = VSRA(r0, VBRADIX); r0 = VAND(r0, vbmask); r1 = VADD(r1, c);
  c = VSRA(r1, VBRADIX); r1 = VAND(r1, vbmask); r2 = VADD(r2, c);
  c = VSRA(r2, VBRADIX); r2 = VAND(r2, vbmask); r3 = VADD(r3, c);
  c = VSRA(r3, VBRADIX); r3 = VAND(r3, vbmask); r4 = VADD(r4, c);
  c = VSRA(r4, VBRADIX); r4 = VAND(r4, vbmask); r5 = VADD(r5, c);
  c = VSRA(r5, VBRADIX); r5 = VAND(r5, vbmask); 
  r0 = VMADD(r0, 0xAA, r0, VSHUF(c, 0x4E));

  r[0] = r0; r[1] = r1; r[2] = r2; r[3] = r3; r[4] = r4; r[5] = r5; 
}

// integer multiplication r = a * b (operand-scanning)
// len(r) = 2*VNWORDS; len(a) = len(b) = NOWRDS
// no carry propagation 
void mp_mul_4x2w(__m512i *r, const __m512i *a, const __m512i *b)
{
  __m512i a0 = a[0], a1 = a[1], a2 = a[2], a3 = a[3], a4 = a[4], a5 = a[5];
  __m512i b0 = b[0], b1 = b[1], b2 = b[2], b3 = b[3], b4 = b[4], b5 = b[5];
  __m512i r0,  r1,  r2,  r3,  r4,  r5,  r6,  r7,  r8,  r9;
  __m512i r10, r11, r12, r13, r14, r15, r16, r17, tb;
  __m512i z0  = VZERO, z1  = VZERO, z2  = VZERO, z3  = VZERO, z4  = VZERO;
  __m512i z5  = VZERO, z6  = VZERO, z7  = VZERO, z8  = VZERO, z9  = VZERO;
  __m512i z10 = VZERO, z11 = VZERO, z12 = VZERO, z13 = VZERO, z14 = VZERO;
  __m512i z15 = VZERO, z16 = VZERO, z17 = VZERO;
  __m512i h0  = VZERO, h1  = VZERO, h2  = VZERO, h3  = VZERO, h4  = VZERO;
  __m512i h5  = VZERO, h6  = VZERO, h7  = VZERO, h8  = VZERO, h9  = VZERO;
  __m512i h10 = VZERO, h11 = VZERO, h12 = VZERO, h13 = VZERO, h14 = VZERO;
  __m512i h15 = VZERO, h16 = VZERO, h17 = VZERO;
  const __m512i vbmask = VSET1(VBMASK);

  // integer multiplication 
  tb = VSHUF(b0, 0x44); 
  z0 = VMACLO(z0, tb, a0); z1 = VMACLO(z1, tb, a1); z2 = VMACLO(z2, tb, a2); 
  z3 = VMACLO(z3, tb, a3); z4 = VMACLO(z4, tb, a4); z5 = VMACLO(z5, tb, a5);
  h0 = VMACHI(h0, tb, a0); h1 = VMACHI(h1, tb, a1); h2 = VMACHI(h2, tb, a2);
  h3 = VMACHI(h3, tb, a3); h4 = VMACHI(h4, tb, a4); h5 = VMACHI(h5, tb, a5);

  tb = VSHUF(b1, 0x44);
  z1 = VMACLO(z1, tb, a0); z2 = VMACLO(z2, tb, a1); z3 = VMACLO(z3, tb, a2);
  z4 = VMACLO(z4, tb, a3); z5 = VMACLO(z5, tb, a4); z6 = VMACLO(z6, tb, a5);
  h1 = VMACHI(h1, tb, a0); h2 = VMACHI(h2, tb, a1); h3 = VMACHI(h3, tb, a2);
  h4 = VMACHI(h4, tb, a3); h5 = VMACHI(h5, tb, a4); h6 = VMACHI(h6, tb, a5);

  tb = VSHUF(b2, 0x44);
  z2 = VMACLO(z2, tb, a0); z3 = VMACLO(z3, tb, a1); z4 = VMACLO(z4, tb, a2);
  z5 = VMACLO(z5, tb, a3); z6 = VMACLO(z6, tb, a4); z7 = VMACLO(z7, tb, a5);
  h2 = VMACHI(h2, tb, a0); h3 = VMACHI(h3, tb, a1); h4 = VMACHI(h4, tb, a2);
  h5 = VMACHI(h5, tb, a3); h6 = VMACHI(h6, tb, a4); h7 = VMACHI(h7, tb, a5);

  tb = VSHUF(b3, 0x44);
  z3 = VMACLO(z3, tb, a0); z4  = VMACLO(z4,  tb, a1); z5 = VMACLO(z5, tb, a2);
  z6 = VMACLO(z6, tb, a3); z7  = VMACLO(z7,  tb, a4); z8 = VMACLO(z8, tb, a5);
  h3 = VMACHI(h3, tb, a0); h4  = VMACHI(h4,  tb, a1); h5 = VMACHI(h5, tb, a2);
  h6 = VMACHI(h6, tb, a3); h7  = VMACHI(h7,  tb, a4); h8 = VMACHI(h8, tb, a5);

  tb = VSHUF(b4, 0x44);
  z4  = VMACLO(z4,  tb, a0); z5  = VMACLO(z5,  tb, a1); z6 = VMACLO(z6, tb, a2);
  z7  = VMACLO(z7,  tb, a3); z8  = VMACLO(z8,  tb, a4); z9 = VMACLO(z9, tb, a5);
  h4  = VMACHI(h4,  tb, a0); h5  = VMACHI(h5,  tb, a1); h6 = VMACHI(h6, tb, a2);
  h7  = VMACHI(h7,  tb, a3); h8  = VMACHI(h8,  tb, a4); h9 = VMACHI(h9, tb, a5);

  tb = VSHUF(b5, 0x44);
  z5  = VMACLO(z5,  tb, a0); z6  = VMACLO(z6,  tb, a1); z7  = VMACLO(z7,  tb, a2);
  z8  = VMACLO(z8,  tb, a3); z9  = VMACLO(z9,  tb, a4); z10 = VMACLO(z10, tb, a5);
  h5  = VMACHI(h5,  tb, a0); h6  = VMACHI(h6,  tb, a1); h7  = VMACHI(h7,  tb, a2);
  h8  = VMACHI(h8,  tb, a3); h9  = VMACHI(h9,  tb, a4); h10 = VMACHI(h10, tb, a5);

  tb = VSHUF(b0, 0xEE);
  z6  = VMACLO(z6,  tb, a0); z7  = VMACLO(z7,  tb, a1); z8  = VMACLO(z8,  tb, a2);
  z9  = VMACLO(z9,  tb, a3); z10 = VMACLO(z10, tb, a4); z11 = VMACLO(z11, tb, a5);
  h6  = VMACHI(h6,  tb, a0); h7  = VMACHI(h7,  tb, a1); h8  = VMACHI(h8,  tb, a2);
  h9  = VMACHI(h9,  tb, a3); h10 = VMACHI(h10, tb, a4); h11 = VMACHI(h11, tb, a5);

  tb = VSHUF(b1, 0xEE);
  z7  = VMACLO(z7,  tb, a0); z8  = VMACLO(z8,  tb, a1); z9  = VMACLO(z9,  tb, a2);
  z10 = VMACLO(z10, tb, a3); z11 = VMACLO(z11, tb, a4); z12 = VMACLO(z12, tb, a5);
  h7  = VMACHI(h7,  tb, a0); h8  = VMACHI(h8,  tb, a1); h9  = VMACHI(h9,  tb, a2);
  h10 = VMACHI(h10, tb, a3); h11 = VMACHI(h11, tb, a4); h12 = VMACHI(h12, tb, a5);

  tb = VSHUF(b2, 0xEE);
  z8  = VMACLO(z8,  tb, a0); z9  = VMACLO(z9,  tb, a1); z10 = VMACLO(z10, tb, a2);
  z11 = VMACLO(z11, tb, a3); z12 = VMACLO(z12, tb, a4); z13 = VMACLO(z13, tb, a5);
  h8  = VMACHI(h8,  tb, a0); h9  = VMACHI(h9,  tb, a1); h10 = VMACHI(h10, tb, a2);
  h11 = VMACHI(h11, tb, a3); h12 = VMACHI(h12, tb, a4); h13 = VMACHI(h13, tb, a5);

  tb = VSHUF(b3, 0xEE);
  z9  = VMACLO(z9,  tb, a0); z10 = VMACLO(z10, tb, a1); z11 = VMACLO(z11, tb, a2);
  z12 = VMACLO(z12, tb, a3); z13 = VMACLO(z13, tb, a4); z14 = VMACLO(z14, tb, a5);
  h9  = VMACHI(h9,  tb, a0); h10 = VMACHI(h10, tb, a1); h11 = VMACHI(h11, tb, a2);
  h12 = VMACHI(h12, tb, a3); h13 = VMACHI(h13, tb, a4); h14 = VMACHI(h14, tb, a5);

  tb = VSHUF(b4, 0xEE);
  z10 = VMACLO(z10, tb, a0); z11 = VMACLO(z11, tb, a1); z12 = VMACLO(z12, tb, a2);
  z13 = VMACLO(z13, tb, a3); z14 = VMACLO(z14, tb, a4); z15 = VMACLO(z15, tb, a5);
  h10 = VMACHI(h10, tb, a0); h11 = VMACHI(h11, tb, a1); h12 = VMACHI(h12, tb, a2);
  h13 = VMACHI(h13, tb, a3); h14 = VMACHI(h14, tb, a4); h15 = VMACHI(h15, tb, a5);

  tb = VSHUF(b5, 0xEE);
  z11 = VMACLO(z11, tb, a0); z12 = VMACLO(z12, tb, a1); z13 = VMACLO(z13, tb, a2);
  z14 = VMACLO(z14, tb, a3); z15 = VMACLO(z15, tb, a4); z16 = VMACLO(z16, tb, a5);
  h11 = VMACHI(h11, tb, a0); h12 = VMACHI(h12, tb, a1); h13 = VMACHI(h13, tb, a2);
  h14 = VMACHI(h14, tb, a3); h15 = VMACHI(h15, tb, a4); h16 = VMACHI(h16, tb, a5);
   
  // merge Z and H 
  z0  = z0;
  z1  = VADD(z1,  VADD(h0,  h0));  z2  = VADD(z2,  VADD(h1,  h1));
  z3  = VADD(z3,  VADD(h2,  h2));  z4  = VADD(z4,  VADD(h3,  h3));
  z5  = VADD(z5,  VADD(h4,  h4));  z6  = VADD(z6,  VADD(h5,  h5));
  z7  = VADD(z7,  VADD(h6,  h6));  z8  = VADD(z8,  VADD(h7,  h7));
  z9  = VADD(z9,  VADD(h8,  h8));  z10 = VADD(z10, VADD(h9,  h9));
  z11 = VADD(z11, VADD(h10, h10)); z12 = VADD(z12, VADD(h11, h11));
  z13 = VADD(z13, VADD(h12, h12)); z14 = VADD(z14, VADD(h13, h13));
  z15 = VADD(z15, VADD(h14, h14)); z16 = VADD(z16, VADD(h15, h15));
  z17 = VADD(z17, VADD(h16, h16)); 

  r[0]  = z0;  r[1]  = z1;  r[2]  = z2;  r[3]  = z3;  r[4]  = z4;
  r[5]  = z5;  r[6]  = z6;  r[7]  = z7;  r[8]  = z8;  r[9]  = z9;
  r[10] = z10; r[11] = z11; r[12] = z12; r[13] = z13; r[14] = z14;
  r[15] = z15; r[16] = z16; r[17] = z17; 
}

// Montgomery reduction r = a * R^-1 mod 2p, where R = 2^612 (operand-scanning)
void rdc_mont_4x2w(__m512i *r, const __m512i *a)
{
  __m512i z0  = a[0],  z1  = a[1],  z2  = a[2],  z3  = a[3],  z4  = a[4];
  __m512i z5  = a[5],  z6  = a[6],  z7  = a[7],  z8  = a[8],  z9  = a[9];
  __m512i z10 = a[10], z11 = a[11], z12 = a[12], z13 = a[13], z14 = a[14];
  __m512i z15 = a[15], z16 = a[16], z17 = a[17];
  __m512i h0  = VZERO, h1  = VZERO, h2  = VZERO, h3  = VZERO, h4  = VZERO;
  __m512i h5  = VZERO, h6  = VZERO, h7  = VZERO, h8  = VZERO, h9  = VZERO;
  __m512i h10 = VZERO, h11 = VZERO, h12 = VZERO, h13 = VZERO, h14 = VZERO; 
  __m512i h15 = VZERO, h16 = VZERO, h17 = VZERO;
  const __m512i vp5  = VSET(vp610p1[11], vp610p1[5],  vp610p1[11], vp610p1[5],  vp610p1[11], vp610p1[5],  vp610p1[11], vp610p1[5]);
  const __m512i vp6  = VSET(vp610p1[7],  vp610p1[6],  vp610p1[7],  vp610p1[6],  vp610p1[7],  vp610p1[6],  vp610p1[7],  vp610p1[6]);
  const __m512i vp8  = VSET(vp610p1[9],  vp610p1[8],  vp610p1[9],  vp610p1[8],  vp610p1[9],  vp610p1[8],  vp610p1[9],  vp610p1[8]);
  const __m512i vp10 = VSET(         0,  vp610p1[10],          0,  vp610p1[10],          0,  vp610p1[10],          0,  vp610p1[10]);
  const __m512i _vp10 = VSET(vp610p1[10],          0,  vp610p1[10],          0,  vp610p1[10],          0,  vp610p1[10],          0);
  const __m512i vbmask = VSET1(VBMASK), vzero = VZERO;
  __m512i r0, r1, r2, r3, r4, r5;
  __m512i c, u, t0, t1, t2, t3;

  c = VSRA(z0, VBRADIX); z0 = VAND(z0, vbmask); z1 = VADD(z1, c);
  u = VSHUF(z0, 0x44); 
  z5 = VMACLO(z5, u, vp5); 
  h5 = VMACHI(h5, u, vp5);
  t0 = VMACLO(vzero, u, vp6);  // t0 = lo(u*p7 | u*p6)
  t1 = VMACHI(vzero, u, vp6);  // t1 = hi(u*p7 | u*p6)
  t2 = VMACLO(vzero, u, vp8);  // t2 = lo(u*p9 | u*p8)
  t3 = VMACHI(vzero, u, vp8);  // t3 = hi(u*p9 | u*p8)
  z10 = VMACLO(z10, u, vp10);
  h10 = VMACHI(h10, u, vp10);
  z6 = VMADD(z6, 0x55, z6, t0);               // z6 = z12 | z6+u*p6
  z7 = VMADD(z7, 0x55, z7, VSHUF(t0, 0x4E));  // z7 = z13 | z7+u*p7
  h6 = VMADD(h6, 0x55, h6, t1);               // h6 = h12 | h6+u*p6
  h7 = VMADD(h7, 0x55, h7, VSHUF(t1, 0x4E));  // h7 = h13 | h7+u*p7
  z2 = VMADD(z2, 0xAA, z2, VSHUF(t2, 0x4E));  // z2 = z8+u*p8 | z2
  z3 = VMADD(z3, 0xAA, z3, t2);               // z3 = z9+u*p9 | z3
  h2 = VMADD(h2, 0xAA, h2, VSHUF(t3, 0x4E));  // h2 = h8+u*p8 | h2
  h3 = VMADD(h3, 0xAA, h3, t3);               // h3 = h9+u*p9 | h3
  z6 = VMADD(z6, 0x55, z6, VSHUF(z0, 0x4E));  

  c = VSRA(z1, VBRADIX); z1 = VAND(z1, vbmask); z2 = VADD(z2, c);
  u = VSHUF(z1, 0x44);
  z6 = VMACLO(z6, u, vp5); h6 = VMACHI(h6, u, vp5);
  t0 = VMACLO(vzero, u, vp6);  t1 = VMACHI(vzero, u, vp6);
  t2 = VMACLO(vzero, u, vp8);  t3 = VMACHI(vzero, u, vp8);
  z11 = VMACLO(z11, u, vp10); h11 = VMACHI(h11, u, vp10);
  z7 = VMADD(z7, 0x55, z7, t0); z8 = VMADD(z8, 0x55, z8, VSHUF(t0, 0x4E));
  h7 = VMADD(h7, 0x55, h7, t1); h8 = VMADD(h8, 0x55, h8, VSHUF(t1, 0x4E));
  z3 = VMADD(z3, 0xAA, z3, VSHUF(t2, 0x4E)); z4 = VMADD(z4, 0xAA, z4, t2);
  h3 = VMADD(h3, 0xAA, h3, VSHUF(t3, 0x4E)); h4 = VMADD(h4, 0xAA, h4, t3); z7 = VMADD(z7, 0x55, z7, VSHUF(z1, 0x4E));

  c = VSRA(z2, VBRADIX); z2 = VAND(z2, vbmask); z3 = VADD(z3, c);
  u = VSHUF(z2, 0x44);
  z7 = VMACLO(z7, u, vp5); h7 = VMACHI(h7, u, vp5);
  t0 = VMACLO(vzero, u, vp6);  t1 = VMACHI(vzero, u, vp6);
  t2 = VMACLO(vzero, u, vp8);  t3 = VMACHI(vzero, u, vp8);
  z12 = VMACLO(z12, u, vp10); h12 = VMACHI(h12, u, vp10);
  z8 = VMADD(z8, 0x55, z8, t0); z9 = VMADD(z9, 0x55, z9, VSHUF(t0, 0x4E));
  h8 = VMADD(h8, 0x55, h8, t1); h9 = VMADD(h9, 0x55, h9, VSHUF(t1, 0x4E));
  z4 = VMADD(z4, 0xAA, z4, VSHUF(t2, 0x4E)); z5 = VMADD(z5, 0xAA, z5, t2);
  h4 = VMADD(h4, 0xAA, h4, VSHUF(t3, 0x4E)); h5 = VMADD(h5, 0xAA, h5, t3); z8 = VMADD(z8, 0x55, z8, VSHUF(z2, 0x4E));
  z3 = VADD(z3, VADD(h2, h2));

  c = VSRA(z3, VBRADIX); z3 = VAND(z3, vbmask); z4 = VADD(z4, c);
  u = VSHUF(z3, 0x44);
  z8 = VMACLO(z8, u, vp5); h8 = VMACHI(h8, u, vp5);
  t0 = VMACLO(vzero, u, vp6);  t1 = VMACHI(vzero, u, vp6);
  t2 = VMACLO(vzero, u, vp8);  t3 = VMACHI(vzero, u, vp8);
  z13 = VMACLO(z13, u, vp10); h13 = VMACHI(h13, u, vp10);
  z9 = VMADD(z9, 0x55, z9, t0); z10 = VMADD(z10, 0x55, z10, VSHUF(t0, 0x4E));
  h9 = VMADD(h9, 0x55, h9, t1); h10 = VMADD(h10, 0x55, h10, VSHUF(t1, 0x4E));
  z5 = VMADD(z5, 0xAA, z5, VSHUF(t2, 0x4E)); z6 = VMADD(z6, 0xAA, z6, t2);
  h5 = VMADD(h5, 0xAA, h5, VSHUF(t3, 0x4E)); h6 = VMADD(h6, 0xAA, h6, t3); z9 = VMADD(z9, 0x55, z9, VSHUF(z3, 0x4E));
  z4 = VADD(z4, VADD(h3, h3));

  c = VSRA(z4, VBRADIX); z4 = VAND(z4, vbmask); z5 = VADD(z5, c);
  u = VSHUF(z4, 0x44);
  z9 = VMACLO(z9, u, vp5); h9 = VMACHI(h9, u, vp5);
  t0 = VMACLO(vzero, u, vp6);  t1 = VMACHI(vzero, u, vp6);
  t2 = VMACLO(vzero, u, vp8);  t3 = VMACHI(vzero, u, vp8);
  z14 = VMACLO(z14, u, vp10); h14 = VMACHI(h14, u, vp10);
  z10 = VMADD(z10, 0x55, z10, t0); z11 = VMADD(z11, 0x55, z11, VSHUF(t0, 0x4E));
  h10 = VMADD(h10, 0x55, h10, t1); h11 = VMADD(h11, 0x55, h11, VSHUF(t1, 0x4E));
  z6 = VMADD(z6, 0xAA, z6, VSHUF(t2, 0x4E)); z7 = VMADD(z7, 0xAA, z7, t2);
  h6 = VMADD(h6, 0xAA, h6, VSHUF(t3, 0x4E)); h7 = VMADD(h7, 0xAA, h7, t3); z10 = VMADD(z10, 0x55, z10, VSHUF(z4, 0x4E));
  z5 = VADD(z5, VADD(h4, h4));

  c = VSRA(z5, VBRADIX); z5 = VAND(z5, vbmask); z6 = VADD(z6, c);
  u = VSHUF(z5, 0x44);
  z10 = VMACLO(z10, u, vp5); h10 = VMACHI(h10, u, vp5);
  t0 = VMACLO(vzero, u, vp6);  t1 = VMACHI(vzero, u, vp6);
  t2 = VMACLO(vzero, u, vp8);  t3 = VMACHI(vzero, u, vp8);
  z15 = VMACLO(z15, u, vp10); h15 = VMACHI(h15, u, vp10);
  z11 = VMADD(z11, 0x55, z11, t0); z12 = VMADD(z12, 0x55, z12, VSHUF(t0, 0x4E));
  h11 = VMADD(h11, 0x55, h11, t1); h12 = VMADD(h12, 0x55, h12, VSHUF(t1, 0x4E));
  z7 = VMADD(z7, 0xAA, z7, VSHUF(t2, 0x4E)); z8 = VMADD(z8, 0xAA, z8, t2);
  h7 = VMADD(h7, 0xAA, h7, VSHUF(t3, 0x4E)); h8 = VMADD(h8, 0xAA, h8, t3); z11 = VMADD(z11, 0x55, z11, VSHUF(z5, 0x4E));
  z6 = VADD(z6, VADD(h5, h5));

  c = VSRA(z6, VBRADIX); z6 = VAND(z6, vbmask); z7 = VADD(z7, c);
  u = VSHUF(z6, 0x44);
  z11 = VMACLO(z11, u, vp5); h11 = VMACHI(h11, u, vp5);
  t0 = VMACLO(vzero, u, vp6);  t1 = VMACHI(vzero, u, vp6);
  t2 = VMACLO(vzero, u, vp8);  t3 = VMACHI(vzero, u, vp8);
  z16 = VMACLO(z16, u, vp10); h16 = VMACHI(h16, u, vp10);
  z12 = VMADD(z12, 0x55, z12, t0); z13 = VMADD(z13, 0x55, z13, VSHUF(t0, 0x4E));
  h12 = VMADD(h12, 0x55, h12, t1); h13 = VMADD(h13, 0x55, h13, VSHUF(t1, 0x4E));
  z8 = VMADD(z8, 0xAA, z8, VSHUF(t2, 0x4E)); z9 = VMADD(z9, 0xAA, z9, t2);
  h8 = VMADD(h8, 0xAA, h8, VSHUF(t3, 0x4E)); h9 = VMADD(h9, 0xAA, h9, t3); z12 = VMADD(z12, 0x55, z12, VSHUF(z6, 0x4E));
  z7 = VADD(z7, VADD(h6, h6));

  c = VSRA(z7, VBRADIX); z7 = VAND(z7, vbmask); z8 = VADD(z8, c);
  u = VSHUF(z7, 0x44);
  z12 = VMACLO(z12, u, vp5); h12 = VMACHI(h12, u, vp5);
  t0 = VMACLO(vzero, u, vp6);  t1 = VMACHI(vzero, u, vp6);
  t2 = VMACLO(vzero, u, vp8);  t3 = VMACHI(vzero, u, vp8);
  z17 = VMACLO(z17, u, vp10); h11 = VMACHI(h11, u, _vp10);
  z13 = VMADD(z13, 0x55, z13, t0); z14 = VMADD(z14, 0x55, z14, VSHUF(t0, 0x4E));
  h13 = VMADD(h13, 0x55, h13, t1); h14 = VMADD(h14, 0x55, h14, VSHUF(t1, 0x4E));
  z9 = VMADD(z9, 0xAA, z9, VSHUF(t2, 0x4E)); z10 = VMADD(z10, 0xAA, z10, t2);
  h9 = VMADD(h9, 0xAA, h9, VSHUF(t3, 0x4E)); h10 = VMADD(h10, 0xAA, h10, t3); z13 = VMADD(z13, 0x55, z13, VSHUF(z7, 0x4E));
  z8 = VADD(z8, VADD(h7, h7));

  c = VSRA(z8, VBRADIX); z8 = VAND(z8, vbmask); z9 = VADD(z9, c);
  u = VSHUF(z8, 0x44);
  z13 = VMACLO(z13, u, vp5); h13 = VMACHI(h13, u, vp5);
  t0 = VMACLO(vzero, u, vp6);  t1 = VMACHI(vzero, u, vp6);
  t2 = VMACLO(vzero, u, vp8);  t3 = VMACHI(vzero, u, vp8);
  z12 = VMACLO(z12, u, _vp10); h12 = VMACHI(h12, u, _vp10);
  z14 = VMADD(z14, 0x55, z14, t0); z15 = VMADD(z15, 0x55, z15, VSHUF(t0, 0x4E));
  h14 = VMADD(h14, 0x55, h14, t1); h15 = VMADD(h15, 0x55, h15, VSHUF(t1, 0x4E));
  z10 = VMADD(z10, 0xAA, z10, VSHUF(t2, 0x4E)); z11 = VMADD(z11, 0xAA, z11, t2);
  h10 = VMADD(h10, 0xAA, h10, VSHUF(t3, 0x4E)); h11 = VMADD(h11, 0xAA, h11, t3); z14 = VMADD(z14, 0x55, z14, VSHUF(z8, 0x4E));
  z9 = VADD(z9, VADD(h8, h8));

  c = VSRA(z9, VBRADIX); z9 = VAND(z9, vbmask); z10 = VADD(z10, c);
  u = VSHUF(z9, 0x44);
  z14 = VMACLO(z14, u, vp5); h14 = VMACHI(h14, u, vp5);
  t0 = VMACLO(vzero, u, vp6);  t1 = VMACHI(vzero, u, vp6);
  t2 = VMACLO(vzero, u, vp8);  t3 = VMACHI(vzero, u, vp8);
  z13 = VMACLO(z13, u, _vp10); h13 = VMACHI(h13, u, _vp10);
  z15 = VMADD(z15, 0x55, z15, t0); z16 = VMADD(z16, 0x55, z16, VSHUF(t0, 0x4E));
  h15 = VMADD(h15, 0x55, h15, t1); h16 = VMADD(h16, 0x55, h16, VSHUF(t1, 0x4E));
  z11 = VMADD(z11, 0xAA, z11, VSHUF(t2, 0x4E)); z12 = VMADD(z12, 0xAA, z12, t2);
  h11 = VMADD(h11, 0xAA, h11, VSHUF(t3, 0x4E)); h12 = VMADD(h12, 0xAA, h12, t3); z15 = VMADD(z15, 0x55, z15, VSHUF(z9, 0x4E));
  z10 = VADD(z10, VADD(h9, h9));

  c = VSRA(z10, VBRADIX); z10 = VAND(z10, vbmask); z11 = VADD(z11, c);
  u = VSHUF(z10, 0x44);
  z15 = VMACLO(z15, u, vp5); h15 = VMACHI(h15, u, vp5);
  t0 = VMACLO(vzero, u, vp6);  t1 = VMACHI(vzero, u, vp6);
  t2 = VMACLO(vzero, u, vp8);  t3 = VMACHI(vzero, u, vp8);
  z14 = VMACLO(z14, u, _vp10); h14 = VMACHI(h14, u, _vp10);
  z16 = VMADD(z16, 0x55, z16, t0); z17 = VMADD(z17, 0x55, z17, VSHUF(t0, 0x4E));
  h16 = VMADD(h16, 0x55, h16, t1); h11 = VMADD(h11, 0xAA, h11, t1);
  z12 = VMADD(z12, 0xAA, z12, VSHUF(t2, 0x4E)); z13 = VMADD(z13, 0xAA, z13, t2);
  h12 = VMADD(h12, 0xAA, h12, VSHUF(t3, 0x4E)); h13 = VMADD(h13, 0xAA, h13, t3); z16 = VMADD(z16, 0x55, z16, VSHUF(z10, 0x4E));
  z11 = VADD(z11, VADD(h10, h10));

  c = VSRA(z11, VBRADIX); z11 = VAND(z11, vbmask); z12 = VADD(z12, c);
  u = VSHUF(z11, 0x44);
  z16 = VMACLO(z16, u, vp5); h16 = VMACHI(h16, u, vp5);
  t0 = VMACLO(vzero, u, vp6);  t1 = VMACHI(vzero, u, vp6);
  t2 = VMACLO(vzero, u, vp8);  t3 = VMACHI(vzero, u, vp8);
  z15 = VMACLO(z15, u, _vp10); h15 = VMACHI(h15, u, _vp10);
  z17 = VMADD(z17, 0x55, z17, t0); z12 = VMADD(z12, 0xAA, z12, t0);
  h11 = VMADD(h11, 0xAA, h11, VSHUF(t1, 0x4E)); h12 = VMADD(h12, 0xAA, h12, t1);
  z13 = VMADD(z13, 0xAA, z13, VSHUF(t2, 0x4E)); z14 = VMADD(z14, 0xAA, z14, t2);
  h13 = VMADD(h13, 0xAA, h13, VSHUF(t3, 0x4E)); h14 = VMADD(h14, 0xAA, h14, t3); z17 = VMADD(z17, 0x55, z17, VSHUF(z11, 0x4E));
  z12 = VADD(z12, VADD(h11, h11));

  z13 = VADD(z13, VADD(h12, h12));
  z14 = VADD(z14, VADD(h13, h13));
  z15 = VADD(z15, VADD(h14, h14));
  z16 = VADD(z16, VADD(h15, h15));
  z17 = VADD(z17, VADD(h16, h16));

  // ---------------------------------------------------------------------------
  // *simple* carry propagation
  // some limbs are finally 52-bit not 51-bit 

  r0 = z12; r1 = z13; r2 = z14; r3 = z15; r4 = z16; r5 = z17;

  c = VSRA(r0, VBRADIX); r0 = VAND(r0, vbmask); r1 = VADD(r1, c);
  c = VSRA(r1, VBRADIX); r1 = VAND(r1, vbmask); r2 = VADD(r2, c);
  c = VSRA(r2, VBRADIX); r2 = VAND(r2, vbmask); r3 = VADD(r3, c);
  c = VSRA(r3, VBRADIX); r3 = VAND(r3, vbmask); r4 = VADD(r4, c);
  c = VSRA(r4, VBRADIX); r4 = VAND(r4, vbmask); r5 = VADD(r5, c);
  c = VSRA(r5, VBRADIX); r5 = VAND(r5, vbmask); 
  r0 = VMADD(r0, 0xAA, r0, VSHUF(c, 0x4E));

  r[0] = r0; r[1] = r1; r[2] = r2; r[3] = r3; r[4] = r4; r[5] = r5; 
}

// -----------------------------------------------------------------------------
// 1-way x64 Fp arithmetic (from PQCrypto-SIDH-3.4) 

void mp_mul(const digit_t* a, const digit_t* b, digit_t* c, const unsigned int nwords)
{ // Multiprecision multiply, c = a*b, where lng(a) = lng(b) = nwords.
    
    mul610_asm(a, b, c);
}

void rdc_mont(digit_t* ma, digit_t* mc)
{ // Montgomery reduction exploiting special form of the prime.
  // mc = ma*R^-1 mod p503x2, where R = 2^512.
  // If ma < 2^512*p503, the output mc is in the range [0, 2*p503-1].
  // ma is assumed to be in Montgomery representation.
   
    rdc610_asm(ma, mc);    
}

void fpcorrection(digit_t* a)
{ // Modular correction to reduce field element a in [0, 2*p503-1] to [0, p503-1].
    unsigned int i, borrow = 0;
    digit_t mask;

    for (i = 0; i < NWORDS_FIELD; i++) {
        SUBC(borrow, a[i], ((digit_t*)p610)[i], borrow, a[i]); 
    }
    mask = 0 - (digit_t)borrow;

    borrow = 0;
    for (i = 0; i < NWORDS_FIELD; i++) {
        ADDC(borrow, a[i], ((digit_t*)p610)[i] & mask, borrow, a[i]); 
    }
}

void fpadd(const digit_t* a, const digit_t* b, digit_t* c)
{ // Modular addition, c = a+b mod p503.
  // Inputs: a, b in [0, 2*p503-1] 
  // Output: c in [0, 2*p503-1]               
    
    fpadd610_asm(a, b, c);    
} 

void fpneg(digit_t* a)
{ // Modular negation, a = -a mod p503.
  // Input/output: a in [0, 2*p503-1] 
    unsigned int i, borrow = 0;
    
    for (i = 0; i < NWORDS_FIELD; i++) {
        SUBC(borrow, ((digit_t*)p610x2)[i], a[i], borrow, a[i]); 
    }
}

void fpsub(const digit_t* a, const digit_t* b, digit_t* c)
{ // Modular subtraction, c = a-b mod p503.
  // Inputs: a, b in [0, 2*p503-1] 
  // Output: c in [0, 2*p503-1] 
    
    fpsub610_asm(a, b, c);    
}

void mp_sub_p4(const digit_t* a, const digit_t* b, digit_t* c)
{ // Multiprecision subtraction with correction with 4*p, c = a-b+4p.                   
    
    mp_sub610_p4_asm(a, b, c);    
} 
