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
  __m512i a0  = a[0],  a1  = a[1],  a2  = a[2],  a3  = a[3],  a4  = a[4];
  __m512i a5  = a[5],  a6  = a[6],  a7  = a[7],  a8  = a[8],  a9  = a[9];
  __m512i a10 = a[10], a11 = a[11], a12 = a[12], a13 = a[13], a14 = a[14];
  __m512i b0  = b[0],  b1  = b[1],  b2  = b[2],  b3  = b[3],  b4  = b[4];
  __m512i b5  = b[5],  b6  = b[6],  b7  = b[7],  b8  = b[8],  b9  = b[9];
  __m512i b10 = b[10], b11 = b[11], b12 = b[12], b13 = b[13], b14 = b[14];
  __m512i r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, r12, r13, r14;
  const __m512i vp0  = VSET1(vp751x2[0]),  vp1  = VSET1(vp751x2[1]),  vp2  = VSET1(vp751x2[2]);
  const __m512i vp3  = VSET1(vp751x2[3]),  vp4  = VSET1(vp751x2[4]),  vp5  = VSET1(vp751x2[5]);
  const __m512i vp6  = VSET1(vp751x2[6]),  vp7  = VSET1(vp751x2[7]),  vp8  = VSET1(vp751x2[8]);
  const __m512i vp9  = VSET1(vp751x2[9]),  vp10 = VSET1(vp751x2[10]), vp11 = VSET1(vp751x2[11]);
  const __m512i vp12 = VSET1(vp751x2[12]), vp13 = VSET1(vp751x2[13]), vp14 = VSET1(vp751x2[14]);
  const __m512i vbmask = VSET1(VBMASK); 

  // r = a + 2p 
  r0  = VADD(a0,  vp0);  r1  = VADD(a1,  vp1);  r2  = VADD(a2,  vp2); 
  r3  = VADD(a3,  vp3);  r4  = VADD(a4,  vp4);  r5  = VADD(a5,  vp5);
  r6  = VADD(a6,  vp6);  r7  = VADD(a7,  vp7);  r8  = VADD(a8,  vp8);
  r9  = VADD(a9,  vp9);  r10 = VADD(a10, vp10); r11 = VADD(a11, vp11);
  r12 = VADD(a12, vp12); r13 = VADD(a13, vp13); r14 = VADD(a14, vp14);

  // r = r - b
  r0  = VSUB(r0,  b0);  r1  = VSUB(r1,  b1);  r2  = VSUB(r2,  b2);
  r3  = VSUB(r3,  b3);  r4  = VSUB(r4,  b4);  r5  = VSUB(r5,  b5);
  r6  = VSUB(r6,  b6);  r7  = VSUB(r7,  b7);  r8  = VSUB(r8,  b8);
  r9  = VSUB(r9,  b9);  r10 = VSUB(r10, b10); r11 = VSUB(r11, b11);
  r12 = VSUB(r12, b12); r13 = VSUB(r13, b13); r14 = VSUB(r14, b14);

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

  r[0]  = r0;  r[1]  = r1;  r[2]  = r2;  r[3]  = r3;  r[4]  = r4; 
  r[5]  = r5;  r[6]  = r6;  r[7]  = r7;  r[8]  = r8;  r[9]  = r9;
  r[10] = r10; r[11] = r11; r[12] = r12; r[13] = r13; r[14] = r14;
}

// integer subtraction r = a - b + 4p
void mp_sub_p4_8x1w(__m512i *r, const __m512i *a, const __m512i *b)
{
  __m512i a0  = a[0],  a1  = a[1],  a2  = a[2],  a3  = a[3],  a4  = a[4];
  __m512i a5  = a[5],  a6  = a[6],  a7  = a[7],  a8  = a[8],  a9  = a[9];
  __m512i a10 = a[10], a11 = a[11], a12 = a[12], a13 = a[13], a14 = a[14];
  __m512i b0  = b[0],  b1  = b[1],  b2  = b[2],  b3  = b[3],  b4  = b[4];
  __m512i b5  = b[5],  b6  = b[6],  b7  = b[7],  b8  = b[8],  b9  = b[9];
  __m512i b10 = b[10], b11 = b[11], b12 = b[12], b13 = b[13], b14 = b[14];
  __m512i r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, r12, r13, r14;
  const __m512i vp0  = VSET1(vp751x4[0]),  vp1  = VSET1(vp751x4[1]),  vp2  = VSET1(vp751x4[2]);
  const __m512i vp3  = VSET1(vp751x4[3]),  vp4  = VSET1(vp751x4[4]),  vp5  = VSET1(vp751x4[5]);
  const __m512i vp6  = VSET1(vp751x4[6]),  vp7  = VSET1(vp751x4[7]),  vp8  = VSET1(vp751x4[8]);
  const __m512i vp9  = VSET1(vp751x4[9]),  vp10 = VSET1(vp751x4[10]), vp11 = VSET1(vp751x4[11]);
  const __m512i vp12 = VSET1(vp751x4[12]), vp13 = VSET1(vp751x4[13]), vp14 = VSET1(vp751x4[14]);
  const __m512i vbmask = VSET1(VBMASK); 

  // r = a + 4p 
  r0  = VADD(a0,  vp0);  r1  = VADD(a1,  vp1);  r2  = VADD(a2,  vp2); 
  r3  = VADD(a3,  vp3);  r4  = VADD(a4,  vp4);  r5  = VADD(a5,  vp5);
  r6  = VADD(a6,  vp6);  r7  = VADD(a7,  vp7);  r8  = VADD(a8,  vp8);
  r9  = VADD(a9,  vp9);  r10 = VADD(a10, vp10); r11 = VADD(a11, vp11);
  r12 = VADD(a12, vp12); r13 = VADD(a13, vp13); r14 = VADD(a14, vp14);

  // r = r - b
  r0  = VSUB(r0,  b0);  r1  = VSUB(r1,  b1);  r2  = VSUB(r2,  b2);
  r3  = VSUB(r3,  b3);  r4  = VSUB(r4,  b4);  r5  = VSUB(r5,  b5);
  r6  = VSUB(r6,  b6);  r7  = VSUB(r7,  b7);  r8  = VSUB(r8,  b8);
  r9  = VSUB(r9,  b9);  r10 = VSUB(r10, b10); r11 = VSUB(r11, b11);
  r12 = VSUB(r12, b12); r13 = VSUB(r13, b13); r14 = VSUB(r14, b14);

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

  r[0]  = r0;  r[1]  = r1;  r[2]  = r2;  r[3]  = r3;  r[4]  = r4; 
  r[5]  = r5;  r[6]  = r6;  r[7]  = r7;  r[8]  = r8;  r[9]  = r9;
  r[10] = r10; r[11] = r11; r[12] = r12; r[13] = r13; r[14] = r14;
}

// modular addition r = a + b mod 2p
void fpadd_8x1w(__m512i *r, const __m512i *a, const __m512i* b)
{
  __m512i a0  = a[0],  a1  = a[1],  a2  = a[2],  a3  = a[3],  a4 = a[4];
  __m512i a5  = a[5],  a6  = a[6],  a7  = a[7],  a8  = a[8],  a9 = a[9];
  __m512i a10 = a[10], a11 = a[11], a12 = a[12], a13 = a[13], a14 = a[14];
  __m512i b0  = b[0],  b1  = b[1],  b2  = b[2],  b3  = b[3],  b4 = b[4];
  __m512i b5  = b[5],  b6  = b[6],  b7  = b[7],  b8  = b[8],  b9 = b[9];
  __m512i b10 = b[10], b11 = b[11], b12 = b[12], b13 = b[13], b14 = b[14];
  __m512i r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, r12, r13, r14, smask;
  const __m512i vp0  = VSET1(vp751x2[0]),  vp1  = VSET1(vp751x2[1]),  vp2  = VSET1(vp751x2[2]);
  const __m512i vp3  = VSET1(vp751x2[3]),  vp4  = VSET1(vp751x2[4]),  vp5  = VSET1(vp751x2[5]);
  const __m512i vp6  = VSET1(vp751x2[6]),  vp7  = VSET1(vp751x2[7]),  vp8  = VSET1(vp751x2[8]);
  const __m512i vp9  = VSET1(vp751x2[9]),  vp10 = VSET1(vp751x2[10]), vp11 = VSET1(vp751x2[11]);
  const __m512i vp12 = VSET1(vp751x2[12]), vp13 = VSET1(vp751x2[13]), vp14 = VSET1(vp751x2[14]);
  const __m512i vbmask = VSET1(VBMASK); 

  // r = a + b
  r0  = VADD(a0,  b0);  r1  = VADD(a1,  b1);  r2  = VADD(a2,  b2); 
  r3  = VADD(a3,  b3);  r4  = VADD(a4,  b4);  r5  = VADD(a5,  b5);
  r6  = VADD(a6,  b6);  r7  = VADD(a7,  b7);  r8  = VADD(a8,  b8);
  r9  = VADD(a9,  b9);  r10 = VADD(a10, b10); r11 = VADD(a11, b11);
  r12 = VADD(a12, b12); r13 = VADD(a13, b13); r14 = VADD(a14, b14);

  // r = r - 2p
  r0  = VSUB(r0,  vp0);  r1  = VSUB(r1,  vp1);  r2  = VSUB(r2,  vp2);
  r3  = VSUB(r3,  vp3);  r4  = VSUB(r4,  vp4);  r5  = VSUB(r5,  vp5);
  r6  = VSUB(r6,  vp6);  r7  = VSUB(r7,  vp7);  r8  = VSUB(r8,  vp8);
  r9  = VSUB(r9,  vp9);  r10 = VSUB(r10, vp10); r11 = VSUB(r11, vp11);
  r12 = VSUB(r12, vp12); r13 = VSUB(r13, vp13); r14 = VSUB(r14, vp14);

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

  // if r is non-negative, smask = all-0 
  // if r is     negative, smask = all-1
  smask = VSRA(r14, 63);
  // r = r + (2p & smask)
  r0  = VADD(r0,  VAND(vp0,  smask)); r1  = VADD(r1,  VAND(vp1,  smask)); 
  r2  = VADD(r2,  VAND(vp2,  smask)); r3  = VADD(r3,  VAND(vp3,  smask)); 
  r4  = VADD(r4,  VAND(vp4,  smask)); r5  = VADD(r5,  VAND(vp5,  smask)); 
  r6  = VADD(r6,  VAND(vp6,  smask)); r7  = VADD(r7,  VAND(vp7,  smask)); 
  r8  = VADD(r8,  VAND(vp8,  smask)); r9  = VADD(r9,  VAND(vp9,  smask)); 
  r10 = VADD(r10, VAND(vp10, smask)); r11 = VADD(r11, VAND(vp11, smask)); 
  r12 = VADD(r12, VAND(vp12, smask)); r13 = VADD(r13, VAND(vp13, smask)); 
  r14 = VADD(r14, VAND(vp14, smask)); 

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

  r[0]  = r0;  r[1]  = r1;  r[2]  = r2;  r[3]  = r3;  r[4] = r4; 
  r[5]  = r5;  r[6]  = r6;  r[7]  = r7;  r[8]  = r8;  r[9] = r9;
  r[10] = r10; r[11] = r11; r[12] = r12; r[13] = r13; r[14] = r14;
}

// modular subtraction r = a - b mod 2p
void fpsub_8x1w(__m512i *r, const __m512i *a, const __m512i* b)
{
  __m512i a0  = a[0],  a1  = a[1],  a2  = a[2],  a3  = a[3],  a4 = a[4];
  __m512i a5  = a[5],  a6  = a[6],  a7  = a[7],  a8  = a[8],  a9 = a[9];
  __m512i a10 = a[10], a11 = a[11], a12 = a[12], a13 = a[13], a14 = a[14];
  __m512i b0  = b[0],  b1  = b[1],  b2  = b[2],  b3  = b[3],  b4 = b[4];
  __m512i b5  = b[5],  b6  = b[6],  b7  = b[7],  b8  = b[8],  b9 = b[9];
  __m512i b10 = b[10], b11 = b[11], b12 = b[12], b13 = b[13], b14 = b[14];
  __m512i r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, r12, r13, r14, smask;
  const __m512i vp0  = VSET1(vp751x2[0]),  vp1  = VSET1(vp751x2[1]),  vp2  = VSET1(vp751x2[2]);
  const __m512i vp3  = VSET1(vp751x2[3]),  vp4  = VSET1(vp751x2[4]),  vp5  = VSET1(vp751x2[5]);
  const __m512i vp6  = VSET1(vp751x2[6]),  vp7  = VSET1(vp751x2[7]),  vp8  = VSET1(vp751x2[8]);
  const __m512i vp9  = VSET1(vp751x2[9]),  vp10 = VSET1(vp751x2[10]), vp11 = VSET1(vp751x2[11]);
  const __m512i vp12 = VSET1(vp751x2[12]), vp13 = VSET1(vp751x2[13]), vp14 = VSET1(vp751x2[14]);
  const __m512i vbmask = VSET1(VBMASK); 

  // r = a - b
  r0  = VSUB(a0,  b0);  r1  = VSUB(a1,  b1);  r2  = VSUB(a2,  b2); 
  r3  = VSUB(a3,  b3);  r4  = VSUB(a4,  b4);  r5  = VSUB(a5,  b5);
  r6  = VSUB(a6,  b6);  r7  = VSUB(a7,  b7);  r8  = VSUB(a8,  b8);
  r9  = VSUB(a9,  b9);  r10 = VSUB(a10, b10); r11 = VSUB(a11, b11);
  r12 = VSUB(a12, b12); r13 = VSUB(a13, b13); r14 = VSUB(a14, b14);

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

  // if r is non-negative, smask = all-0 
  // if r is     negative, smask = all-1
  smask = VSRA(r14, 63);
  // r = r + (2p & smask)
  r0  = VADD(r0,  VAND(vp0,  smask)); r1  = VADD(r1,  VAND(vp1,  smask)); 
  r2  = VADD(r2,  VAND(vp2,  smask)); r3  = VADD(r3,  VAND(vp3,  smask)); 
  r4  = VADD(r4,  VAND(vp4,  smask)); r5  = VADD(r5,  VAND(vp5,  smask)); 
  r6  = VADD(r6,  VAND(vp6,  smask)); r7  = VADD(r7,  VAND(vp7,  smask)); 
  r8  = VADD(r8,  VAND(vp8,  smask)); r9  = VADD(r9,  VAND(vp9,  smask)); 
  r10 = VADD(r10, VAND(vp10, smask)); r11 = VADD(r11, VAND(vp11, smask)); 
  r12 = VADD(r12, VAND(vp12, smask)); r13 = VADD(r13, VAND(vp13, smask)); 
  r14 = VADD(r14, VAND(vp14, smask)); 

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

  r[0]  = r0;  r[1]  = r1;  r[2]  = r2;  r[3]  = r3;  r[4] = r4; 
  r[5]  = r5;  r[6]  = r6;  r[7]  = r7;  r[8]  = r8;  r[9] = r9;
  r[10] = r10; r[11] = r11; r[12] = r12; r[13] = r13; r[14] = r14;
}

// Karatsuba 
void mp_mul_8x1w(__m512i *r, const __m512i *a, const __m512i *b)
{
  __m512i a0  = a[0],  a1  = a[1],  a2  = a[2],  a3  = a[3],  a4 = a[4];
  __m512i a5  = a[5],  a6  = a[6],  a7  = a[7],  a8  = a[8],  a9 = a[9];
  __m512i a10 = a[10], a11 = a[11], a12 = a[12], a13 = a[13], a14 = a[14];
  __m512i b0  = b[0],  b1  = b[1],  b2  = b[2],  b3  = b[3],  b4 = b[4];
  __m512i b5  = b[5],  b6  = b[6],  b7  = b[7],  b8  = b[8],  b9 = b[9];
  __m512i b10 = b[10], b11 = b[11], b12 = b[12], b13 = b[13], b14 = b[14];
  __m512i ta0, ta1, ta2, ta3, ta4, ta5, ta6, ta7;
  __m512i tb0, tb1, tb2, tb3, tb4, tb5, tb6, tb7;
  __m512i r0,  r1,  r2,  r3,  r4,  r5,  r6,  r7,  r8,  r9,  r10, r11, r12, r13, r14;
  __m512i r15, r16, r17, r18, r19, r20, r21, r22, r23, r24, r25, r26, r27, r28, r29;
  __m512i z0  = VZERO, z1  = VZERO, z2  = VZERO, z3  = VZERO, z4  = VZERO;
  __m512i z5  = VZERO, z6  = VZERO, z7  = VZERO, z8  = VZERO, z9  = VZERO;
  __m512i z10 = VZERO, z11 = VZERO, z12 = VZERO, z13 = VZERO, z14 = VZERO;
  __m512i z15 = VZERO, z16 = VZERO, z17 = VZERO, z18 = VZERO, z19 = VZERO;
  __m512i z20 = VZERO, z21 = VZERO, z22 = VZERO, z23 = VZERO, z24 = VZERO;
  __m512i z25 = VZERO, z26 = VZERO, z27 = VZERO, z28 = VZERO, z29 = VZERO;
  __m512i h0  = VZERO, h1  = VZERO, h2  = VZERO, h3  = VZERO, h4  = VZERO;
  __m512i h5  = VZERO, h6  = VZERO, h7  = VZERO, h8  = VZERO, h9  = VZERO;
  __m512i h10 = VZERO, h11 = VZERO, h12 = VZERO, h13 = VZERO, h14 = VZERO;
  __m512i h15 = VZERO, h16 = VZERO, h17 = VZERO, h18 = VZERO, h19 = VZERO;
  __m512i h20 = VZERO, h21 = VZERO, h22 = VZERO, h23 = VZERO, h24 = VZERO;
  __m512i h25 = VZERO, h26 = VZERO, h27 = VZERO, h28 = VZERO, h29 = VZERO;
  __m512i m0  = VZERO, m1  = VZERO, m2  = VZERO, m3  = VZERO, m4  = VZERO;
  __m512i m5  = VZERO, m6  = VZERO, m7  = VZERO, m8  = VZERO, m9  = VZERO;
  __m512i m10 = VZERO, m11 = VZERO, m12 = VZERO, m13 = VZERO, m14 = VZERO;
  __m512i m15 = VZERO;
  const __m512i vbmask = VSET1(VBMASK);

  // ---------------------------------------------------------------------------
  // compute zL(z0~z14) by aL(a0~a7) * bL(b0~b7)
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

  z6 = VMACLO(h6, a0, b6); z6 = VMACLO(z6, a1, b5); z6 = VMACLO(z6, a2, b4); 
  z6 = VMACLO(z6, a3, b3); z6 = VMACLO(z6, a4, b2); z6 = VMACLO(z6, a5, b1);
  z6 = VMACLO(z6, a6, b0);
  h7 = VMACHI(h7, a0, b6); h7 = VMACHI(h7, a1, b5); h7 = VMACHI(h7, a2, b4); 
  h7 = VMACHI(h7, a3, b3); h7 = VMACHI(h7, a4, b2); h7 = VMACHI(h7, a5, b1); 
  h7 = VMACHI(h7, a6, b0); h7 = VADD(h7, h7);

  z7 = VMACLO(h7, a0, b7); z7 = VMACLO(z7, a1, b6); z7 = VMACLO(z7, a2, b5); 
  z7 = VMACLO(z7, a3, b4); z7 = VMACLO(z7, a4, b3); z7 = VMACLO(z7, a5, b2);
  z7 = VMACLO(z7, a6, b1); z7 = VMACLO(z7, a7, b0);
  h8 = VMACHI(h8, a0, b7); h8 = VMACHI(h8, a1, b6); h8 = VMACHI(h8, a2, b5); 
  h8 = VMACHI(h8, a3, b4); h8 = VMACHI(h8, a4, b3); h8 = VMACHI(h8, a5, b2); 
  h8 = VMACHI(h8, a6, b1); h8 = VMACHI(h8, a7, b0); h8 = VADD(h8, h8);

  z8 = VMACLO(h8, a1, b7); z8 = VMACLO(z8, a2, b6); z8 = VMACLO(z8, a3, b5); 
  z8 = VMACLO(z8, a4, b4); z8 = VMACLO(z8, a5, b3); z8 = VMACLO(z8, a6, b2); 
  z8 = VMACLO(z8, a7, b1); 
  h9 = VMACHI(h9, a1, b7); h9 = VMACHI(h9, a2, b6); h9 = VMACHI(h9, a3, b5); 
  h9 = VMACHI(h9, a4, b4); h9 = VMACHI(h9, a5, b3); h9 = VMACHI(h9, a6, b2);
  h9 = VMACHI(h9, a7, b1); h9 = VADD(h9, h9); 

  z9 = VMACLO(h9, a2, b7); z9 = VMACLO(z9, a3, b6); z9 = VMACLO(z9, a4, b5); 
  z9 = VMACLO(z9, a5, b4); z9 = VMACLO(z9, a6, b3); z9 = VMACLO(z9, a7, b2);
  h10 = VMACHI(h10, a2, b7); h10 = VMACHI(h10, a3, b6); h10 = VMACHI(h10, a4, b5); 
  h10 = VMACHI(h10, a5, b4); h10 = VMACHI(h10, a6, b3); h10 = VMACHI(h10, a7, b2);
  h10 = VADD(h10, h10); 

  z10 = VMACLO(h10, a3, b7); z10 = VMACLO(z10, a4, b6); z10 = VMACLO(z10, a5, b5);
  z10 = VMACLO(z10, a6, b4); z10 = VMACLO(z10, a7, b3);
  h11 = VMACHI(h11, a3, b7); h11 = VMACHI(h11, a4, b6); h11 = VMACHI(h11, a5, b5);
  h11 = VMACHI(h11, a6, b4); h11 = VMACHI(h11, a7, b3); h11 = VADD(h11, h11);

  z11 = VMACLO(h11, a4, b7); z11 = VMACLO(z11, a5, b6); z11 = VMACLO(z11, a6, b5);
  z11 = VMACLO(z11, a7, b4);
  h12 = VMACHI(h12, a4, b7); h12 = VMACHI(h12, a5, b6); h12 = VMACHI(h12, a6, b5);
  h12 = VMACHI(h12, a7, b4); h12 = VADD(h12, h12);

  z12 = VMACLO(h12, a5, b7); z12 = VMACLO(z12, a6, b6); z12 = VMACLO(z12, a7, b5);
  h13 = VMACHI(h13, a5, b7); h13 = VMACHI(h13, a6, b6); h13 = VMACHI(h13, a7, b5);
  h13 = VADD(h13, h13);

  z13 = VMACLO(h13, a6, b7); z13 = VMACLO(z13, a7, b6);
  h14 = VMACHI(h14, a6, b7); h14 = VMACHI(h14, a7, b6); h14 = VADD(h14, h14);

  z14 = VMACLO(h14, a7, b7); 
  h15 = VMACHI(h15, a7, b7); h15 = VADD(h15, h15);

  z15 = h15;

  // ---------------------------------------------------------------------------
  // compute zH(z16~z29) by aH(a8~a15) * bH(b8~b15)
  z16 = VMACLO(h16, a8, b8);
  h17 = VMACHI(h17, a8, b8); h17 = VADD(h17, h17);

  z17 = VMACLO(h17, a8, b9); z17 = VMACLO(z17, a9, b8); 
  h18 = VMACHI(h18, a8, b9); h18 = VMACHI(h18, a9, b8); h18 = VADD(h18, h18);

  z18 = VMACLO(h18, a8, b10); z18 = VMACLO(z18, a9, b9); z18 = VMACLO(z18, a10, b8);
  h19 = VMACHI(h19, a8, b10); h19 = VMACHI(h19, a9, b9); h19 = VMACHI(h19, a10, b8);
  h19 = VADD(h19, h19);

  z19 = VMACLO(h19, a8, b11); z19 = VMACLO(z19, a9, b10); z19 = VMACLO(z19, a10, b9); 
  z19 = VMACLO(z19, a11, b8);
  h20 = VMACHI(h20, a8, b11); h20 = VMACHI(h20, a9, b10); h20 = VMACHI(h20, a10, b9); 
  h20 = VMACHI(h20, a11, b8); h20 = VADD(h20, h20);

  z20 = VMACLO(h20, a8, b12); z20 = VMACLO(z20, a9, b11); z20 = VMACLO(z20, a10, b10); 
  z20 = VMACLO(z20, a11, b9); z20 = VMACLO(z20, a12, b8);
  h21 = VMACHI(h21, a8, b12); h21 = VMACHI(h21, a9, b11); h21 = VMACHI(h21, a10, b10); 
  h21 = VMACHI(h21, a11, b9); h21 = VMACHI(h21, a12, b8); h21 = VADD(h21, h21);

  z21 = VMACLO(h21, a8,  b13); z21 = VMACLO(z21, a9,  b12); z21 = VMACLO(z21, a10, b11); 
  z21 = VMACLO(z21, a11, b10); z21 = VMACLO(z21, a12, b9);  z21 = VMACLO(z21, a13, b8);
  h22 = VMACHI(h22, a8,  b13); h22 = VMACHI(h22, a9,  b12); h22 = VMACHI(h22, a10, b11); 
  h22 = VMACHI(h22, a11, b10); h22 = VMACHI(h22, a12, b9);  h22 = VMACHI(h22, a13, b8);
  h22 = VADD(h22, h22);

  z22 = VMACLO(h22, a8, b14); z22 = VMACLO(z22, a9, b13); z22 = VMACLO(z22, a10, b12); 
  z22 = VMACLO(z22, a11, b11); z22 = VMACLO(z22, a12, b10); z22 = VMACLO(z22, a13, b9);
  z22 = VMACLO(z22, a14, b8);
  h23 = VMACHI(h23, a8, b14); h23 = VMACHI(h23, a9, b13); h23 = VMACHI(h23, a10, b12); 
  h23 = VMACHI(h23, a11, b11); h23 = VMACHI(h23, a12, b10); h23 = VMACHI(h23, a13, b9);
  h23 = VMACHI(h23, a14, b8); h23 = VADD(h23, h23);

  z23 = VMACLO(h23, a9, b14); z23 = VMACLO(z23, a10, b13); z23 = VMACLO(z23, a11, b12); 
  z23 = VMACLO(z23, a12, b11); z23 = VMACLO(z23, a13, b10); z23 = VMACLO(z23, a14, b9);
  h24 = VMACHI(h24, a9, b14); h24 = VMACHI(h24, a10, b13); h24 = VMACHI(h24, a11, b12); 
  h24 = VMACHI(h24, a12, b11); h24 = VMACHI(h24, a13, b10); h24 = VMACHI(h24, a14, b9);
  h24 = VADD(h24, h24);

  z24 = VMACLO(h24, a10, b14); z24 = VMACLO(z24, a11, b13); z24 = VMACLO(z24, a12, b12); 
  z24 = VMACLO(z24, a13, b11); z24 = VMACLO(z24, a14, b10); 
  h25 = VMACHI(h25, a10, b14); h25 = VMACHI(h25, a11, b13); h25 = VMACHI(h25, a12, b12); 
  h25 = VMACHI(h25, a13, b11); h25 = VMACHI(h25, a14, b10); h25 = VADD(h25, h25); 

  z25 = VMACLO(h25, a11, b14); z25 = VMACLO(z25, a12, b13); z25 = VMACLO(z25, a13, b12); 
  z25 = VMACLO(z25, a14, b11); 
  h26 = VMACHI(h26, a11, b14); h26 = VMACHI(h26, a12, b13); h26 = VMACHI(h26, a13, b12); 
  h26 = VMACHI(h26, a14, b11); h26 = VADD(h26, h26); 

  z26 = VMACLO(h26, a12, b14); z26 = VMACLO(z26, a13, b13); z26 = VMACLO(z26, a14, b12); 
  h27 = VMACHI(h27, a12, b14); h27 = VMACHI(h27, a13, b13); h27 = VMACHI(h27, a14, b12); 
  h27 = VADD(h27, h27);

  z27 = VMACLO(h27, a13, b14); z27 = VMACLO(z27, a14, b13);
  h28 = VMACHI(h28, a13, b14); h28 = VMACHI(h28, a14, b13); h28 = VADD(h28, h28);

  z28 = VMACLO(h28, a14, b14); 
  h29 = VMACHI(h29, a14, b14); h29 = VADD(h29, h29);

  z29 = h29;

  // ---------------------------------------------------------------------------
  // ta(ta0~ta7) = aL(a0~a7) + aH(a8~a14) 
  ta0 = VADD(a0, a8); ta1 = VADD(a1, a9); ta2 = VADD(a2, a10);
  ta3 = VADD(a3, a11); ta4 = VADD(a4, a12); ta5 = VADD(a5, a13);
  ta6 = VADD(a6, a14); ta7 = a7;

  // tb(tb0~tb7) = bL(b0~b7) + bH(b8~b14)
  tb0 = VADD(b0, b8); tb1 = VADD(b1, b9); tb2 = VADD(b2, b10);
  tb3 = VADD(b3, b11); tb4 = VADD(b4, b12); tb5 = VADD(b5, b13);
  tb6 = VADD(b6, b14); tb7 = b7;

  // ---------------------------------------------------------------------------
  // zM = ta * tb - zL - zH

  h0  = h1  = h2  = h3  = h4  = h5  = h6  = h7  = VZERO;
  h8  = h9  = h10 = h11 = h12 = h13 = h14 = h15 = VZERO;

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

  m6 = VMACLO(h6, ta0, tb6); m6 = VMACLO(m6, ta1, tb5); m6 = VMACLO(m6, ta2, tb4); 
  m6 = VMACLO(m6, ta3, tb3); m6 = VMACLO(m6, ta4, tb2); m6 = VMACLO(m6, ta5, tb1);
  m6 = VMACLO(m6, ta6, tb0);
  h7 = VMACHI(h7, ta0, tb6); h7 = VMACHI(h7, ta1, tb5); h7 = VMACHI(h7, ta2, tb4); 
  h7 = VMACHI(h7, ta3, tb3); h7 = VMACHI(h7, ta4, tb2); h7 = VMACHI(h7, ta5, tb1); 
  h7 = VMACHI(h7, ta6, tb0); h7 = VADD(h7, h7);

  m7 = VMACLO(h7, ta0, tb7); m7 = VMACLO(m7, ta1, tb6); m7 = VMACLO(m7, ta2, tb5); 
  m7 = VMACLO(m7, ta3, tb4); m7 = VMACLO(m7, ta4, tb3); m7 = VMACLO(m7, ta5, tb2);
  m7 = VMACLO(m7, ta6, tb1); m7 = VMACLO(m7, ta7, tb0);
  h8 = VMACHI(h8, ta0, tb7); h8 = VMACHI(h8, ta1, tb6); h8 = VMACHI(h8, ta2, tb5); 
  h8 = VMACHI(h8, ta3, tb4); h8 = VMACHI(h8, ta4, tb3); h8 = VMACHI(h8, ta5, tb2); 
  h8 = VMACHI(h8, ta6, tb1); h8 = VMACHI(h8, ta7, tb0); h8 = VADD(h8, h8);

  m8 = VMACLO(h8, ta1, tb7); m8 = VMACLO(m8, ta2, tb6); m8 = VMACLO(m8, ta3, tb5); 
  m8 = VMACLO(m8, ta4, tb4); m8 = VMACLO(m8, ta5, tb3); m8 = VMACLO(m8, ta6, tb2); 
  m8 = VMACLO(m8, ta7, tb1); 
  h9 = VMACHI(h9, ta1, tb7); h9 = VMACHI(h9, ta2, tb6); h9 = VMACHI(h9, ta3, tb5); 
  h9 = VMACHI(h9, ta4, tb4); h9 = VMACHI(h9, ta5, tb3); h9 = VMACHI(h9, ta6, tb2);
  h9 = VMACHI(h9, ta7, tb1); h9 = VADD(h9, h9); 

  m9 = VMACLO(h9, ta2, tb7); m9 = VMACLO(m9, ta3, tb6); m9 = VMACLO(m9, ta4, tb5); 
  m9 = VMACLO(m9, ta5, tb4); m9 = VMACLO(m9, ta6, tb3); m9 = VMACLO(m9, ta7, tb2);
  h10 = VMACHI(h10, ta2, tb7); h10 = VMACHI(h10, ta3, tb6); h10 = VMACHI(h10, ta4, tb5); 
  h10 = VMACHI(h10, ta5, tb4); h10 = VMACHI(h10, ta6, tb3); h10 = VMACHI(h10, ta7, tb2);
  h10 = VADD(h10, h10); 

  m10 = VMACLO(h10, ta3, tb7); m10 = VMACLO(m10, ta4, tb6); m10 = VMACLO(m10, ta5, tb5);
  m10 = VMACLO(m10, ta6, tb4); m10 = VMACLO(m10, ta7, tb3);
  h11 = VMACHI(h11, ta3, tb7); h11 = VMACHI(h11, ta4, tb6); h11 = VMACHI(h11, ta5, tb5);
  h11 = VMACHI(h11, ta6, tb4); h11 = VMACHI(h11, ta7, tb3); h11 = VADD(h11, h11);

  m11 = VMACLO(h11, ta4, tb7); m11 = VMACLO(m11, ta5, tb6); m11 = VMACLO(m11, ta6, tb5);
  m11 = VMACLO(m11, ta7, tb4);
  h12 = VMACHI(h12, ta4, tb7); h12 = VMACHI(h12, ta5, tb6); h12 = VMACHI(h12, ta6, tb5);
  h12 = VMACHI(h12, ta7, tb4); h12 = VADD(h12, h12);

  m12 = VMACLO(h12, ta5, tb7); m12 = VMACLO(m12, ta6, tb6); m12 = VMACLO(m12, ta7, tb5);
  h13 = VMACHI(h13, ta5, tb7); h13 = VMACHI(h13, ta6, tb6); h13 = VMACHI(h13, ta7, tb5);
  h13 = VADD(h13, h13);

  m13 = VMACLO(h13, ta6, tb7); m13 = VMACLO(m13, ta7, tb6);
  h14 = VMACHI(h14, ta6, tb7); h14 = VMACHI(h14, ta7, tb6); h14 = VADD(h14, h14);

  m14 = VMACLO(h14, ta7, tb7); 
  h15 = VMACHI(h15, ta7, tb7); h15 = VADD(h15, h15);

  m15 = h15;

  m0  = VSUB(m0,  VADD(z0,  z16)); m1  = VSUB(m1,  VADD(z1, z17));
  m2  = VSUB(m2,  VADD(z2,  z18)); m3  = VSUB(m3,  VADD(z3, z19));
  m4  = VSUB(m4,  VADD(z4,  z20)); m5  = VSUB(m5,  VADD(z5, z21));
  m6  = VSUB(m6,  VADD(z6,  z22)); m7  = VSUB(m7,  VADD(z7, z23));
  m8  = VSUB(m8,  VADD(z8,  z24)); m9  = VSUB(m9,  VADD(z9, z25));
  m10 = VSUB(m10, VADD(z10, z26)); m11 = VSUB(m11, VADD(z11, z27));
  m12 = VSUB(m12, VADD(z12, z28)); m13 = VSUB(m13, VADD(z13, z29));
  m14 = VSUB(m14, z14);            m15 = VSUB(m15, z15);


  // z = z + zM
  z8  = VADD(z8,  m0);  z9  = VADD(z9,  m1); 
  z10 = VADD(z10, m2);  z11 = VADD(z11, m3); 
  z12 = VADD(z12, m4);  z13 = VADD(z13, m5); 
  z14 = VADD(z14, m6);  z15 = VADD(z15, m7); 
  z16 = VADD(z16, m8);  z17 = VADD(z17, m9);
  z18 = VADD(z18, m10); z19 = VADD(z19, m11);
  z20 = VADD(z20, m12); z21 = VADD(z21, m13);
  z22 = VADD(z22, m14); z23 = VADD(z23, m15);

  // ---------------------------------------------------------------------------
  r[0]  = z0;  r[1]  = z1;  r[2]  = z2;  r[3]  = z3;  r[4]  = z4;  
  r[5]  = z5;  r[6]  = z6;  r[7]  = z7;  r[8]  = z8;  r[9]  = z9;  
  r[10] = z10; r[11] = z11; r[12] = z12; r[13] = z13; r[14] = z14; 
  r[15] = z15; r[16] = z16; r[17] = z17; r[18] = z18; r[19] = z19;
  r[20] = z20; r[21] = z21; r[22] = z22; r[23] = z23; r[24] = z24; 
  r[25] = z25; r[26] = z26; r[27] = z27; r[28] = z28; r[29] = z29;
}

// Montgomery reduction r = a * R^-1 mod 2p, where R = 2^765 (operand-scanning)
void rdc_mont_8x1w(__m512i *r, const __m512i *a)
{
  __m512i a0  = a[0],  a1  = a[1],  a2  = a[2],  a3  = a[3],  a4  = a[4];
  __m512i a5  = a[5],  a6  = a[6],  a7  = a[7],  a8  = a[8],  a9  = a[9];
  __m512i a10 = a[10], a11 = a[11], a12 = a[12], a13 = a[13], a14 = a[14];
  __m512i a15 = a[15], a16 = a[16], a17 = a[17], a18 = a[18], a19 = a[19];
  __m512i a20 = a[20], a21 = a[21], a22 = a[22], a23 = a[23], a24 = a[24];
  __m512i a25 = a[25], a26 = a[26], a27 = a[27], a28 = a[28], a29 = a[29];
  __m512i h8  = VZERO, h9  = VZERO, h10 = VZERO, h11 = VZERO, h12 = VZERO;
  __m512i h13 = VZERO, h14 = VZERO, h15 = VZERO, h16 = VZERO, h17 = VZERO;
  __m512i h18 = VZERO, h19 = VZERO, h20 = VZERO, h21 = VZERO, h22 = VZERO;
  __m512i h23 = VZERO, h24 = VZERO, h25 = VZERO, h26 = VZERO, h27 = VZERO;
  __m512i h28 = VZERO, h29 = VZERO, u;
  const __m512i vp7  = VSET1(vp751p1[7]),  vp8  = VSET1(vp751p1[8]),  vp9  = VSET1(vp751p1[9]);
  const __m512i vp10 = VSET1(vp751p1[10]), vp11 = VSET1(vp751p1[11]), vp12 = VSET1(vp751p1[12]);
  const __m512i vp13 = VSET1(vp751p1[13]), vp14 = VSET1(vp751p1[14]);
  const __m512i vbmask = VSET1(VBMASK); 

  // some pre carry propagations

  a1 = VADD(a1, VSRA(a0, VBRADIX)); a0 = VAND(a0, vbmask);
  a2 = VADD(a2, VSRA(a1, VBRADIX)); a1 = VAND(a1, vbmask);
  a3 = VADD(a3, VSRA(a2, VBRADIX)); a2 = VAND(a2, vbmask);
  a4 = VADD(a4, VSRA(a3, VBRADIX)); a3 = VAND(a3, vbmask);
  a5 = VADD(a5, VSRA(a4, VBRADIX)); a4 = VAND(a4, vbmask);
  a6 = VADD(a6, VSRA(a5, VBRADIX)); a5 = VAND(a5, vbmask);
  a7 = VADD(a7, VSRA(a6, VBRADIX)); a6 = VAND(a6, vbmask);

  u = a0;
  a7  = VMACLO(a7,  u, vp7);  a8  = VMACLO(a8,  u, vp8);  a9  = VMACLO(a9,  u, vp9); 
  a10 = VMACLO(a10, u, vp10); a11 = VMACLO(a11, u, vp11); a12 = VMACLO(a12, u, vp12);
  a13 = VMACLO(a13, u, vp13); a14 = VMACLO(a14, u, vp14);
  h8  = VMACHI(h8,  u, vp7);  h9  = VMACHI(h9,  u, vp8);  h10 = VMACHI(h10, u, vp9);
  h11 = VMACHI(h11, u, vp10); h12 = VMACHI(h12, u, vp11); h13 = VMACHI(h13, u, vp12);
  h14 = VMACHI(h14, u, vp13); h15 = VMACHI(h15, u, vp14); 
  a8 = VADD(a8, VSRA(a7, VBRADIX)); a7 = VAND(a7, vbmask); 
  a8 = VADD(a8, VADD(h8, h8));

  u = a1;
  a8  = VMACLO(a8,  u, vp7);  a9  = VMACLO(a9,  u, vp8);  a10 = VMACLO(a10, u, vp9);
  a11 = VMACLO(a11, u, vp10); a12 = VMACLO(a12, u, vp11); a13 = VMACLO(a13, u, vp12);
  a14 = VMACLO(a14, u, vp13); a15 = VMACLO(a15, u, vp14);
  h9  = VMACHI(h9,  u, vp7);  h10 = VMACHI(h10, u, vp8);  h11 = VMACHI(h11, u, vp9);
  h12 = VMACHI(h12, u, vp10); h13 = VMACHI(h13, u, vp11); h14 = VMACHI(h14, u, vp12);
  h15 = VMACHI(h15, u, vp13); h16 = VMACHI(h16, u, vp14);
  a9 = VADD(a9, VSRA(a8, VBRADIX)); a8 = VAND(a8, vbmask);
  a9 = VADD(a9, VADD(h9, h9));

  u = a2;
  a9  = VMACLO(a9,  u, vp7);  a10 = VMACLO(a10, u, vp8);  a11 = VMACLO(a11, u, vp9);
  a12 = VMACLO(a12, u, vp10); a13 = VMACLO(a13, u, vp11); a14 = VMACLO(a14, u, vp12);
  a15 = VMACLO(a15, u, vp13); a16 = VMACLO(a16, u, vp14);
  h10 = VMACHI(h10, u, vp7);  h11 = VMACHI(h11, u, vp8);  h12 = VMACHI(h12, u, vp9);
  h13 = VMACHI(h13, u, vp10); h14 = VMACHI(h14, u, vp11); h15 = VMACHI(h15, u, vp12);
  h16 = VMACHI(h16, u, vp13); h17 = VMACHI(h17, u, vp14);
  a10 = VADD(a10, VSRA(a9, VBRADIX)); a9 = VAND(a9, vbmask);
  a10 = VADD(a10, VADD(h10, h10));

  u = a3;
  a10 = VMACLO(a10, u, vp7);  a11 = VMACLO(a11, u, vp8);  a12 = VMACLO(a12, u, vp9);
  a13 = VMACLO(a13, u, vp10); a14 = VMACLO(a14, u, vp11); a15 = VMACLO(a15, u, vp12);
  a16 = VMACLO(a16, u, vp13); a17 = VMACLO(a17, u, vp14);
  h11 = VMACHI(h11, u, vp7);  h12 = VMACHI(h12, u, vp8);  h13 = VMACHI(h13, u, vp9);
  h14 = VMACHI(h14, u, vp10); h15 = VMACHI(h15, u, vp11); h16 = VMACHI(h16, u, vp12);
  h17 = VMACHI(h17, u, vp13); h18 = VMACHI(h18, u, vp14);
  a11 = VADD(a11, VSRA(a10, VBRADIX)); a10 = VAND(a10, vbmask);
  a11 = VADD(a11, VADD(h11, h11));

  u = a4;
  a11 = VMACLO(a11, u, vp7);  a12 = VMACLO(a12, u, vp8);  a13 = VMACLO(a13, u, vp9);
  a14 = VMACLO(a14, u, vp10); a15 = VMACLO(a15, u, vp11); a16 = VMACLO(a16, u, vp12);
  a17 = VMACLO(a17, u, vp13); a18 = VMACLO(a18, u, vp14);
  h12 = VMACHI(h12, u, vp7);  h13 = VMACHI(h13, u, vp8);  h14 = VMACHI(h14, u, vp9);
  h15 = VMACHI(h15, u, vp10); h16 = VMACHI(h16, u, vp11); h17 = VMACHI(h17, u, vp12);
  h18 = VMACHI(h18, u, vp13); h19 = VMACHI(h19, u, vp14);
  a12 = VADD(a12, VSRA(a11, VBRADIX)); a11 = VAND(a11, vbmask);
  a12 = VADD(a12, VADD(h12, h12));

  u = a5;
  a12 = VMACLO(a12, u, vp7);  a13 = VMACLO(a13, u, vp8);  a14 = VMACLO(a14, u, vp9);
  a15 = VMACLO(a15, u, vp10); a16 = VMACLO(a16, u, vp11); a17 = VMACLO(a17, u, vp12);
  a18 = VMACLO(a18, u, vp13); a19 = VMACLO(a19, u, vp14);
  h13 = VMACHI(h13, u, vp7);  h14 = VMACHI(h14, u, vp8);  h15 = VMACHI(h15, u, vp9);
  h16 = VMACHI(h16, u, vp10); h17 = VMACHI(h17, u, vp11); h18 = VMACHI(h18, u, vp12);
  h19 = VMACHI(h19, u, vp13); h20 = VMACHI(h20, u, vp14);
  a13 = VADD(a13, VSRA(a12, VBRADIX)); a12 = VAND(a12, vbmask);
  a13 = VADD(a13, VADD(h13, h13));

  u = a6;
  a13 = VMACLO(a13, u, vp7);  a14 = VMACLO(a14, u, vp8);  a15 = VMACLO(a15, u, vp9);
  a16 = VMACLO(a16, u, vp10); a17 = VMACLO(a17, u, vp11); a18 = VMACLO(a18, u, vp12);
  a19 = VMACLO(a19, u, vp13); a20 = VMACLO(a20, u, vp14);
  h14 = VMACHI(h14, u, vp7);  h15 = VMACHI(h15, u, vp8);  h16 = VMACHI(h16, u, vp9);
  h17 = VMACHI(h17, u, vp10); h18 = VMACHI(h18, u, vp11); h19 = VMACHI(h19, u, vp12);
  h20 = VMACHI(h20, u, vp13); h21 = VMACHI(h21, u, vp14);
  a14 = VADD(a14, VSRA(a13, VBRADIX)); a13 = VAND(a13, vbmask);
  a14 = VADD(a14, VADD(h14, h14));

  u = a7;
  a14 = VMACLO(a14, u, vp7);  a15 = VMACLO(a15, u, vp8);  a16 = VMACLO(a16, u, vp9);
  a17 = VMACLO(a17, u, vp10); a18 = VMACLO(a18, u, vp11); a19 = VMACLO(a19, u, vp12);
  a20 = VMACLO(a20, u, vp13); a21 = VMACLO(a21, u, vp14);
  h15 = VMACHI(h15, u, vp7);  h16 = VMACHI(h16, u, vp8);  h17 = VMACHI(h17, u, vp9);
  h18 = VMACHI(h18, u, vp10); h19 = VMACHI(h19, u, vp11); h20 = VMACHI(h20, u, vp12);
  h21 = VMACHI(h21, u, vp13); h22 = VMACHI(h22, u, vp14);
  a15 = VADD(a15, VSRA(a14, VBRADIX)); a14 = VAND(a14, vbmask);
  a15 = VADD(a15, VADD(h15, h15));

  u = a8;
  a15 = VMACLO(a15, u, vp7);  a16 = VMACLO(a16, u, vp8);  a17 = VMACLO(a17, u, vp9);
  a18 = VMACLO(a18, u, vp10); a19 = VMACLO(a19, u, vp11); a20 = VMACLO(a20, u, vp12);
  a21 = VMACLO(a21, u, vp13); a22 = VMACLO(a22, u, vp14);
  h16 = VMACHI(h16, u, vp7);  h17 = VMACHI(h17, u, vp8);  h18 = VMACHI(h18, u, vp9);
  h19 = VMACHI(h19, u, vp10); h20 = VMACHI(h20, u, vp11); h21 = VMACHI(h21, u, vp12);
  h22 = VMACHI(h22, u, vp13); h23 = VMACHI(h23, u, vp14);
  a16 = VADD(a16, VSRA(a15, VBRADIX)); a15 = VAND(a15, vbmask);
  a16 = VADD(a16, VADD(h16, h16));

  u = a9;
  a16 = VMACLO(a16, u, vp7);  a17 = VMACLO(a17, u, vp8);  a18 = VMACLO(a18, u, vp9);
  a19 = VMACLO(a19, u, vp10); a20 = VMACLO(a20, u, vp11); a21 = VMACLO(a21, u, vp12);
  a22 = VMACLO(a22, u, vp13); a23 = VMACLO(a23, u, vp14);
  h17 = VMACHI(h17, u, vp7);  h18 = VMACHI(h18, u, vp8);  h19 = VMACHI(h19, u, vp9);
  h20 = VMACHI(h20, u, vp10); h21 = VMACHI(h21, u, vp11); h22 = VMACHI(h22, u, vp12);
  h23 = VMACHI(h23, u, vp13); h24 = VMACHI(h24, u, vp14);
  a17 = VADD(a17, VSRA(a16, VBRADIX)); a16 = VAND(a16, vbmask);
  a17 = VADD(a17, VADD(h17, h17));

  u = a10;
  a17 = VMACLO(a17, u, vp7);  a18 = VMACLO(a18, u, vp8);  a19 = VMACLO(a19, u, vp9);
  a20 = VMACLO(a20, u, vp10); a21 = VMACLO(a21, u, vp11); a22 = VMACLO(a22, u, vp12);
  a23 = VMACLO(a23, u, vp13); a24 = VMACLO(a24, u, vp14);
  h18 = VMACHI(h18, u, vp7);  h19 = VMACHI(h19, u, vp8);  h20 = VMACHI(h20, u, vp9);
  h21 = VMACHI(h21, u, vp10); h22 = VMACHI(h22, u, vp11); h23 = VMACHI(h23, u, vp12);
  h24 = VMACHI(h24, u, vp13); h25 = VMACHI(h25, u, vp14);
  a18 = VADD(a18, VSRA(a17, VBRADIX)); a17 = VAND(a17, vbmask);
  a18 = VADD(a18, VADD(h18, h18));

  u = a11;
  a18 = VMACLO(a18, u, vp7);  a19 = VMACLO(a19, u, vp8);  a20 = VMACLO(a20, u, vp9);
  a21 = VMACLO(a21, u, vp10); a22 = VMACLO(a22, u, vp11); a23 = VMACLO(a23, u, vp12);
  a24 = VMACLO(a24, u, vp13); a25 = VMACLO(a25, u, vp14);
  h19 = VMACHI(h19, u, vp7);  h20 = VMACHI(h20, u, vp8);  h21 = VMACHI(h21, u, vp9);
  h22 = VMACHI(h22, u, vp10); h23 = VMACHI(h23, u, vp11); h24 = VMACHI(h24, u, vp12);
  h25 = VMACHI(h25, u, vp13); h26 = VMACHI(h26, u, vp14);
  a19 = VADD(a19, VSRA(a18, VBRADIX)); a18 = VAND(a18, vbmask);
  a19 = VADD(a19, VADD(h19, h19));

  u = a12;
  a19 = VMACLO(a19, u, vp7);  a20 = VMACLO(a20, u, vp8);  a21 = VMACLO(a21, u, vp9);
  a22 = VMACLO(a22, u, vp10); a23 = VMACLO(a23, u, vp11); a24 = VMACLO(a24, u, vp12);
  a25 = VMACLO(a25, u, vp13); a26 = VMACLO(a26, u, vp14);
  h20 = VMACHI(h20, u, vp7);  h21 = VMACHI(h21, u, vp8);  h22 = VMACHI(h22, u, vp9);
  h23 = VMACHI(h23, u, vp10); h24 = VMACHI(h24, u, vp11); h25 = VMACHI(h25, u, vp12);
  h26 = VMACHI(h26, u, vp13); h27 = VMACHI(h27, u, vp14);
  a20 = VADD(a20, VSRA(a19, VBRADIX)); a19 = VAND(a19, vbmask);
  a20 = VADD(a20, VADD(h20, h20));

  u = a13;
  a20 = VMACLO(a20, u, vp7);  a21 = VMACLO(a21, u, vp8);  a22 = VMACLO(a22, u, vp9);
  a23 = VMACLO(a23, u, vp10); a24 = VMACLO(a24, u, vp11); a25 = VMACLO(a25, u, vp12);
  a26 = VMACLO(a26, u, vp13); a27 = VMACLO(a27, u, vp14);
  h21 = VMACHI(h21, u, vp7);  h22 = VMACHI(h22, u, vp8);  h23 = VMACHI(h23, u, vp9);
  h24 = VMACHI(h24, u, vp10); h25 = VMACHI(h25, u, vp11); h26 = VMACHI(h26, u, vp12);
  h27 = VMACHI(h27, u, vp13); h28 = VMACHI(h28, u, vp14);
  a21 = VADD(a21, VSRA(a20, VBRADIX)); a20 = VAND(a20, vbmask);
  a21 = VADD(a21, VADD(h21, h21));

  u = a14;
  a21 = VMACLO(a21, u, vp7);  a22 = VMACLO(a22, u, vp8);  a23 = VMACLO(a23, u, vp9);
  a24 = VMACLO(a24, u, vp10); a25 = VMACLO(a25, u, vp11); a26 = VMACLO(a26, u, vp12);
  a27 = VMACLO(a27, u, vp13); a28 = VMACLO(a28, u, vp14);
  h22 = VMACHI(h22, u, vp7);  h23 = VMACHI(h23, u, vp8);  h24 = VMACHI(h24, u, vp9);
  h25 = VMACHI(h25, u, vp10); h26 = VMACHI(h26, u, vp11); h27 = VMACHI(h27, u, vp12);
  h28 = VMACHI(h28, u, vp13); h29 = VMACHI(h29, u, vp14);
  a22 = VADD(a22, VSRA(a21, VBRADIX)); a21 = VAND(a21, vbmask);
  a22 = VADD(a22, VADD(h22, h22));

  a23 = VADD(a23, VADD(h23, h23));
  a24 = VADD(a24, VADD(h24, h24));
  a25 = VADD(a25, VADD(h25, h25));
  a26 = VADD(a26, VADD(h26, h26));
  a27 = VADD(a27, VADD(h27, h27));
  a28 = VADD(a28, VADD(h28, h28));
  a29 = VADD(a29, VADD(h29, h29));

  a23 = VADD(a23, VSRA(a22, VBRADIX)); a22 = VAND(a22, vbmask);
  a24 = VADD(a24, VSRA(a23, VBRADIX)); a23 = VAND(a23, vbmask);
  a25 = VADD(a25, VSRA(a24, VBRADIX)); a24 = VAND(a24, vbmask);
  a26 = VADD(a26, VSRA(a25, VBRADIX)); a25 = VAND(a25, vbmask);
  a27 = VADD(a27, VSRA(a26, VBRADIX)); a26 = VAND(a26, vbmask);
  a28 = VADD(a28, VSRA(a27, VBRADIX)); a27 = VAND(a27, vbmask);
  a29 = VADD(a29, VSRA(a28, VBRADIX)); a28 = VAND(a28, vbmask);

  r[0]  = a15; r[1]  = a16; r[2]  = a17; r[3]  = a18; r[4]  = a19;
  r[5]  = a20; r[6]  = a21; r[7]  = a22; r[8]  = a23; r[9]  = a24;
  r[10] = a25; r[11] = a26; r[12] = a27; r[13] = a28; r[14] = a29;
}

// -----------------------------------------------------------------------------
// (4x2)-way Fp arithmetic 

// integer subtraction r = a - b + 2p
void mp_sub_p2_4x2w(__m512i *r, const __m512i *a, const __m512i *b)
{
  __m512i a0 = a[0], a1 = a[1], a2 = a[2], a3 = a[3];
  __m512i a4 = a[4], a5 = a[5], a6 = a[6], a7 = a[7];
  __m512i b0 = b[0], b1 = b[1], b2 = b[2], b3 = b[3];
  __m512i b4 = b[4], b5 = b[5], b6 = b[6], b7 = b[7];
  __m512i r0, r1, r2, r3, r4, r5, r6, r7, c;
  const __m512i vp0 = VSET(vp751x2[8],  vp751x2[0], vp751x2[8],  vp751x2[0], vp751x2[8],  vp751x2[0], vp751x2[8],  vp751x2[0]);
  const __m512i vp1 = VSET(vp751x2[9],  vp751x2[1], vp751x2[9],  vp751x2[1], vp751x2[9],  vp751x2[1], vp751x2[9],  vp751x2[1]);
  const __m512i vp2 = VSET(vp751x2[10], vp751x2[2], vp751x2[10], vp751x2[2], vp751x2[10], vp751x2[2], vp751x2[10], vp751x2[2]);
  const __m512i vp3 = VSET(vp751x2[11], vp751x2[3], vp751x2[11], vp751x2[3], vp751x2[11], vp751x2[3], vp751x2[11], vp751x2[3]);
  const __m512i vp4 = VSET(vp751x2[12], vp751x2[4], vp751x2[12], vp751x2[4], vp751x2[12], vp751x2[4], vp751x2[12], vp751x2[4]);
  const __m512i vp5 = VSET(vp751x2[13], vp751x2[5], vp751x2[13], vp751x2[5], vp751x2[13], vp751x2[5], vp751x2[13], vp751x2[5]);
  const __m512i vp6 = VSET(vp751x2[14], vp751x2[6], vp751x2[14], vp751x2[6], vp751x2[14], vp751x2[6], vp751x2[14], vp751x2[6]);
  const __m512i vp7 = VSET(          0, vp751x2[7],           0, vp751x2[7],           0, vp751x2[7],           0, vp751x2[7]);   
  const __m512i vbmask = VSET1(VBMASK); 

  // r = a + 2p
  r0 = VADD(a0, vp0); r1 = VADD(a1, vp1); r2 = VADD(a2, vp2); r3 = VADD(a3, vp3); 
  r4 = VADD(a4, vp4); r5 = VADD(a5, vp5); r6 = VADD(a6, vp6); r7 = VADD(a7, vp7);

  // r = r - b
  r0 = VSUB(r0, b0); r1 = VSUB(r1, b1); r2 = VSUB(r2, b2); r3 = VSUB(r3, b3); 
  r4 = VSUB(r4, b4); r5 = VSUB(r5, b5); r6 = VSUB(r6, b6); r7 = VSUB(r7, b7);

  // *simple* carry propagation
  // some limbs are finally 52-bit not 51-bit 

  c = VSRA(r0, VBRADIX); r0 = VAND(r0, vbmask); r1 = VADD(r1, c);
  c = VSRA(r1, VBRADIX); r1 = VAND(r1, vbmask); r2 = VADD(r2, c);
  c = VSRA(r2, VBRADIX); r2 = VAND(r2, vbmask); r3 = VADD(r3, c);
  c = VSRA(r3, VBRADIX); r3 = VAND(r3, vbmask); r4 = VADD(r4, c);
  c = VSRA(r4, VBRADIX); r4 = VAND(r4, vbmask); r5 = VADD(r5, c);
  c = VSRA(r5, VBRADIX); r5 = VAND(r5, vbmask); r6 = VADD(r6, c);  
  c = VSRA(r6, VBRADIX); r6 = VAND(r6, vbmask); r7 = VADD(r7, c);  
  c = VSRA(r7, VBRADIX); r7 = VAND(r7, vbmask);
  r0 = VMADD(r0, 0xAA, r0, VSHUF(c, 0x4E));

  r[0] = r0; r[1] = r1; r[2] = r2; r[3] = r3; 
  r[4] = r4; r[5] = r5; r[6] = r6; r[7] = r7;
}

// integer subtraction r = a - b + 4p
void mp_sub_p4_4x2w(__m512i *r, const __m512i *a, const __m512i *b)
{
  __m512i a0 = a[0], a1 = a[1], a2 = a[2], a3 = a[3];
  __m512i a4 = a[4], a5 = a[5], a6 = a[6], a7 = a[7];
  __m512i b0 = b[0], b1 = b[1], b2 = b[2], b3 = b[3];
  __m512i b4 = b[4], b5 = b[5], b6 = b[6], b7 = b[7];
  __m512i r0, r1, r2, r3, r4, r5, r6, r7, c;
  const __m512i vp0 = VSET(vp751x4[8],  vp751x4[0], vp751x4[8],  vp751x4[0], vp751x4[8],  vp751x4[0], vp751x4[8],  vp751x4[0]);
  const __m512i vp1 = VSET(vp751x4[9],  vp751x4[1], vp751x4[9],  vp751x4[1], vp751x4[9],  vp751x4[1], vp751x4[9],  vp751x4[1]);
  const __m512i vp2 = VSET(vp751x4[10], vp751x4[2], vp751x4[10], vp751x4[2], vp751x4[10], vp751x4[2], vp751x4[10], vp751x4[2]);
  const __m512i vp3 = VSET(vp751x4[11], vp751x4[3], vp751x4[11], vp751x4[3], vp751x4[11], vp751x4[3], vp751x4[11], vp751x4[3]);
  const __m512i vp4 = VSET(vp751x4[12], vp751x4[4], vp751x4[12], vp751x4[4], vp751x4[12], vp751x4[4], vp751x4[12], vp751x4[4]);
  const __m512i vp5 = VSET(vp751x4[13], vp751x4[5], vp751x4[13], vp751x4[5], vp751x4[13], vp751x4[5], vp751x4[13], vp751x4[5]);
  const __m512i vp6 = VSET(vp751x4[14], vp751x4[6], vp751x4[14], vp751x4[6], vp751x4[14], vp751x4[6], vp751x4[14], vp751x4[6]);
  const __m512i vp7 = VSET(          0, vp751x4[7],           0, vp751x4[7],           0, vp751x4[7],           0, vp751x4[7]);   
  const __m512i vbmask = VSET1(VBMASK); 

  // r = a + 4p
  r0 = VADD(a0, vp0); r1 = VADD(a1, vp1); r2 = VADD(a2, vp2); r3 = VADD(a3, vp3); 
  r4 = VADD(a4, vp4); r5 = VADD(a5, vp5); r6 = VADD(a6, vp6); r7 = VADD(a7, vp7);

  // r = r - b
  r0 = VSUB(r0, b0); r1 = VSUB(r1, b1); r2 = VSUB(r2, b2); r3 = VSUB(r3, b3); 
  r4 = VSUB(r4, b4); r5 = VSUB(r5, b5); r6 = VSUB(r6, b6); r7 = VSUB(r7, b7);

  // *simple* carry propagation
  // some limbs are finally 52-bit not 51-bit 

  c = VSRA(r0, VBRADIX); r0 = VAND(r0, vbmask); r1 = VADD(r1, c);
  c = VSRA(r1, VBRADIX); r1 = VAND(r1, vbmask); r2 = VADD(r2, c);
  c = VSRA(r2, VBRADIX); r2 = VAND(r2, vbmask); r3 = VADD(r3, c);
  c = VSRA(r3, VBRADIX); r3 = VAND(r3, vbmask); r4 = VADD(r4, c);
  c = VSRA(r4, VBRADIX); r4 = VAND(r4, vbmask); r5 = VADD(r5, c);
  c = VSRA(r5, VBRADIX); r5 = VAND(r5, vbmask); r6 = VADD(r6, c);  
  c = VSRA(r6, VBRADIX); r6 = VAND(r6, vbmask); r7 = VADD(r7, c);  
  c = VSRA(r7, VBRADIX); r7 = VAND(r7, vbmask);
  r0 = VMADD(r0, 0xAA, r0, VSHUF(c, 0x4E));

  r[0] = r0; r[1] = r1; r[2] = r2; r[3] = r3; 
  r[4] = r4; r[5] = r5; r[6] = r6; r[7] = r7;
}

// modular addition r = a + b mod 2p
void fpadd_4x2w(__m512i *r, const __m512i *a, const __m512i* b)
{
  __m512i a0 = a[0], a1 = a[1], a2 = a[2], a3 = a[3];
  __m512i a4 = a[4], a5 = a[5], a6 = a[6], a7 = a[7];
  __m512i b0 = b[0], b1 = b[1], b2 = b[2], b3 = b[3];
  __m512i b4 = b[4], b5 = b[5], b6 = b[6], b7 = b[7];
  __m512i r0, r1, r2, r3, r4, r5, r6, r7, carry, smask, c;
  const __m512i vp0 = VSET(vp751x2[8],  vp751x2[0], vp751x2[8],  vp751x2[0], vp751x2[8],  vp751x2[0], vp751x2[8],  vp751x2[0]);
  const __m512i vp1 = VSET(vp751x2[9],  vp751x2[1], vp751x2[9],  vp751x2[1], vp751x2[9],  vp751x2[1], vp751x2[9],  vp751x2[1]);
  const __m512i vp2 = VSET(vp751x2[10], vp751x2[2], vp751x2[10], vp751x2[2], vp751x2[10], vp751x2[2], vp751x2[10], vp751x2[2]);
  const __m512i vp3 = VSET(vp751x2[11], vp751x2[3], vp751x2[11], vp751x2[3], vp751x2[11], vp751x2[3], vp751x2[11], vp751x2[3]);
  const __m512i vp4 = VSET(vp751x2[12], vp751x2[4], vp751x2[12], vp751x2[4], vp751x2[12], vp751x2[4], vp751x2[12], vp751x2[4]);
  const __m512i vp5 = VSET(vp751x2[13], vp751x2[5], vp751x2[13], vp751x2[5], vp751x2[13], vp751x2[5], vp751x2[13], vp751x2[5]);
  const __m512i vp6 = VSET(vp751x2[14], vp751x2[6], vp751x2[14], vp751x2[6], vp751x2[14], vp751x2[6], vp751x2[14], vp751x2[6]);
  const __m512i vp7 = VSET(          0, vp751x2[7],           0, vp751x2[7],           0, vp751x2[7],           0, vp751x2[7]);   
  const __m512i vbmask = VSET1(VBMASK); 

  // r = a + b
  r0 = VADD(a0, b0); r1 = VADD(a1, b1); r2 = VADD(a2, b2); r3 = VADD(a3, b3); 
  r4 = VADD(a4, b4); r5 = VADD(a5, b5); r6 = VADD(a6, b6); r7 = VADD(a7, b7);

  // r = r - 2p
  r0 = VSUB(r0, vp0); r1 = VSUB(r1, vp1); r2 = VSUB(r2, vp2); r3 = VSUB(r3, vp3); 
  r4 = VSUB(r4, vp4); r5 = VSUB(r5, vp5); r6 = VSUB(r6, vp6); r7 = VSUB(r7, vp7);

  // carry propagation 
  carry = VSRA(r0, VBRADIX); r0 = VMAND(r0, 0x55, r0, vbmask); r1 = VMADD(r1, 0x55, r1, carry);
  carry = VSRA(r1, VBRADIX); r1 = VMAND(r1, 0x55, r1, vbmask); r2 = VMADD(r2, 0x55, r2, carry);
  carry = VSRA(r2, VBRADIX); r2 = VMAND(r2, 0x55, r2, vbmask); r3 = VMADD(r3, 0x55, r3, carry);
  carry = VSRA(r3, VBRADIX); r3 = VMAND(r3, 0x55, r3, vbmask); r4 = VMADD(r4, 0x55, r4, carry);
  carry = VSRA(r4, VBRADIX); r4 = VMAND(r4, 0x55, r4, vbmask); r5 = VMADD(r5, 0x55, r5, carry);
  carry = VSRA(r5, VBRADIX); r5 = VMAND(r5, 0x55, r5, vbmask); r6 = VMADD(r6, 0x55, r6, carry);
  carry = VSRA(r6, VBRADIX); r6 = VMAND(r6, 0x55, r6, vbmask); r7 = VMADD(r7, 0x55, r7, carry);
  carry = VSRA(r7, VBRADIX); r7 = VMAND(r7, 0x55, r7, vbmask); 
  carry = VZSHUF(0xCCCC, carry, 0x4E); r0 = VMADD(r0, 0xAA, r0, carry); 
  carry = VSRA(r0, VBRADIX); r0 = VMAND(r0, 0xAA, r0, vbmask); r1 = VMADD(r1, 0xAA, r1, carry);
  carry = VSRA(r1, VBRADIX); r1 = VMAND(r1, 0xAA, r1, vbmask); r2 = VMADD(r2, 0xAA, r2, carry);
  carry = VSRA(r2, VBRADIX); r2 = VMAND(r2, 0xAA, r2, vbmask); r3 = VMADD(r3, 0xAA, r3, carry);
  carry = VSRA(r3, VBRADIX); r3 = VMAND(r3, 0xAA, r3, vbmask); r4 = VMADD(r4, 0xAA, r4, carry);
  carry = VSRA(r4, VBRADIX); r4 = VMAND(r4, 0xAA, r4, vbmask); r5 = VMADD(r5, 0xAA, r5, carry);
  carry = VSRA(r5, VBRADIX); r5 = VMAND(r5, 0xAA, r5, vbmask); r6 = VMADD(r6, 0xAA, r6, carry);

  // if r is non-negative, smask = all-0 
  // if r is     negative, smask = all-1
  smask = VSRA(r6, 63);                 // smask = all-0/-1 X | all-0/-1 X | all-0/-1 X | all-0/-1 X | 
  smask = VSHUF(smask, 0xEE);           // smask = all-0/-1 | all-0/-1 | all-0/-1 | all-0/-1
  // r = r + (2p & smask)
  r0 = VADD(r0, VAND(vp0, smask)); r1 = VADD(r1, VAND(vp1, smask));
  r2 = VADD(r2, VAND(vp2, smask)); r3 = VADD(r3, VAND(vp3, smask));
  r4 = VADD(r4, VAND(vp4, smask)); r5 = VADD(r5, VAND(vp5, smask)); 
  r6 = VADD(r6, VAND(vp6, smask)); r7 = VADD(r7, VAND(vp7, smask)); 

  // *simple* carry propagation
  // some limbs are finally 52-bit not 51-bit 

  c = VSRA(r0, VBRADIX); r0 = VAND(r0, vbmask); r1 = VADD(r1, c);
  c = VSRA(r1, VBRADIX); r1 = VAND(r1, vbmask); r2 = VADD(r2, c);
  c = VSRA(r2, VBRADIX); r2 = VAND(r2, vbmask); r3 = VADD(r3, c);
  c = VSRA(r3, VBRADIX); r3 = VAND(r3, vbmask); r4 = VADD(r4, c);
  c = VSRA(r4, VBRADIX); r4 = VAND(r4, vbmask); r5 = VADD(r5, c);
  c = VSRA(r5, VBRADIX); r5 = VAND(r5, vbmask); r6 = VADD(r6, c);  
  c = VSRA(r6, VBRADIX); r6 = VAND(r6, vbmask); r7 = VADD(r7, c);  
  c = VSRA(r7, VBRADIX); r7 = VAND(r7, vbmask);
  r0 = VMADD(r0, 0xAA, r0, VSHUF(c, 0x4E));

  r[0] = r0; r[1] = r1; r[2] = r2; r[3] = r3; 
  r[4] = r4; r[5] = r5; r[6] = r6; r[7] = r7;
}

// modular subtraction r = a + b mod 2p
void fpsub_4x2w(__m512i *r, const __m512i *a, const __m512i* b)
{
  __m512i a0 = a[0], a1 = a[1], a2 = a[2], a3 = a[3];
  __m512i a4 = a[4], a5 = a[5], a6 = a[6], a7 = a[7];
  __m512i b0 = b[0], b1 = b[1], b2 = b[2], b3 = b[3];
  __m512i b4 = b[4], b5 = b[5], b6 = b[6], b7 = b[7];
  __m512i r0, r1, r2, r3, r4, r5, r6, r7, carry, smask, c;
  const __m512i vp0 = VSET(vp751x2[8],  vp751x2[0], vp751x2[8],  vp751x2[0], vp751x2[8],  vp751x2[0], vp751x2[8],  vp751x2[0]);
  const __m512i vp1 = VSET(vp751x2[9],  vp751x2[1], vp751x2[9],  vp751x2[1], vp751x2[9],  vp751x2[1], vp751x2[9],  vp751x2[1]);
  const __m512i vp2 = VSET(vp751x2[10], vp751x2[2], vp751x2[10], vp751x2[2], vp751x2[10], vp751x2[2], vp751x2[10], vp751x2[2]);
  const __m512i vp3 = VSET(vp751x2[11], vp751x2[3], vp751x2[11], vp751x2[3], vp751x2[11], vp751x2[3], vp751x2[11], vp751x2[3]);
  const __m512i vp4 = VSET(vp751x2[12], vp751x2[4], vp751x2[12], vp751x2[4], vp751x2[12], vp751x2[4], vp751x2[12], vp751x2[4]);
  const __m512i vp5 = VSET(vp751x2[13], vp751x2[5], vp751x2[13], vp751x2[5], vp751x2[13], vp751x2[5], vp751x2[13], vp751x2[5]);
  const __m512i vp6 = VSET(vp751x2[14], vp751x2[6], vp751x2[14], vp751x2[6], vp751x2[14], vp751x2[6], vp751x2[14], vp751x2[6]);
  const __m512i vp7 = VSET(          0, vp751x2[7],           0, vp751x2[7],           0, vp751x2[7],           0, vp751x2[7]);   
  const __m512i vbmask = VSET1(VBMASK); 

  // r = a - b
  r0 = VSUB(a0, b0); r1 = VSUB(a1, b1); r2 = VSUB(a2, b2); r3 = VSUB(a3, b3); 
  r4 = VSUB(a4, b4); r5 = VSUB(a5, b5); r6 = VSUB(a6, b6); r7 = VSUB(a7, b7);

  // carry propagation 
  carry = VSRA(r0, VBRADIX); r0 = VMAND(r0, 0x55, r0, vbmask); r1 = VMADD(r1, 0x55, r1, carry);
  carry = VSRA(r1, VBRADIX); r1 = VMAND(r1, 0x55, r1, vbmask); r2 = VMADD(r2, 0x55, r2, carry);
  carry = VSRA(r2, VBRADIX); r2 = VMAND(r2, 0x55, r2, vbmask); r3 = VMADD(r3, 0x55, r3, carry);
  carry = VSRA(r3, VBRADIX); r3 = VMAND(r3, 0x55, r3, vbmask); r4 = VMADD(r4, 0x55, r4, carry);
  carry = VSRA(r4, VBRADIX); r4 = VMAND(r4, 0x55, r4, vbmask); r5 = VMADD(r5, 0x55, r5, carry);
  carry = VSRA(r5, VBRADIX); r5 = VMAND(r5, 0x55, r5, vbmask); r6 = VMADD(r6, 0x55, r6, carry);
  carry = VSRA(r6, VBRADIX); r6 = VMAND(r6, 0x55, r6, vbmask); r7 = VMADD(r7, 0x55, r7, carry);
  carry = VSRA(r7, VBRADIX); r7 = VMAND(r7, 0x55, r7, vbmask); 
  carry = VZSHUF(0xCCCC, carry, 0x4E); r0 = VMADD(r0, 0xAA, r0, carry); 
  carry = VSRA(r0, VBRADIX); r0 = VMAND(r0, 0xAA, r0, vbmask); r1 = VMADD(r1, 0xAA, r1, carry);
  carry = VSRA(r1, VBRADIX); r1 = VMAND(r1, 0xAA, r1, vbmask); r2 = VMADD(r2, 0xAA, r2, carry);
  carry = VSRA(r2, VBRADIX); r2 = VMAND(r2, 0xAA, r2, vbmask); r3 = VMADD(r3, 0xAA, r3, carry);
  carry = VSRA(r3, VBRADIX); r3 = VMAND(r3, 0xAA, r3, vbmask); r4 = VMADD(r4, 0xAA, r4, carry);
  carry = VSRA(r4, VBRADIX); r4 = VMAND(r4, 0xAA, r4, vbmask); r5 = VMADD(r5, 0xAA, r5, carry);
  carry = VSRA(r5, VBRADIX); r5 = VMAND(r5, 0xAA, r5, vbmask); r6 = VMADD(r6, 0xAA, r6, carry);

  // if r is non-negative, smask = all-0 
  // if r is     negative, smask = all-1
  smask = VSRA(r6, 63);                 // smask = all-0/-1 X | all-0/-1 X | all-0/-1 X | all-0/-1 X | 
  smask = VSHUF(smask, 0xEE);           // smask = all-0/-1 | all-0/-1 | all-0/-1 | all-0/-1
  // r = r + (2p & smask)
  r0 = VADD(r0, VAND(vp0, smask)); r1 = VADD(r1, VAND(vp1, smask));
  r2 = VADD(r2, VAND(vp2, smask)); r3 = VADD(r3, VAND(vp3, smask));
  r4 = VADD(r4, VAND(vp4, smask)); r5 = VADD(r5, VAND(vp5, smask)); 
  r6 = VADD(r6, VAND(vp6, smask)); r7 = VADD(r7, VAND(vp7, smask)); 

  // *simple* carry propagation
  // some limbs are finally 52-bit not 51-bit 

  c = VSRA(r0, VBRADIX); r0 = VAND(r0, vbmask); r1 = VADD(r1, c);
  c = VSRA(r1, VBRADIX); r1 = VAND(r1, vbmask); r2 = VADD(r2, c);
  c = VSRA(r2, VBRADIX); r2 = VAND(r2, vbmask); r3 = VADD(r3, c);
  c = VSRA(r3, VBRADIX); r3 = VAND(r3, vbmask); r4 = VADD(r4, c);
  c = VSRA(r4, VBRADIX); r4 = VAND(r4, vbmask); r5 = VADD(r5, c);
  c = VSRA(r5, VBRADIX); r5 = VAND(r5, vbmask); r6 = VADD(r6, c);  
  c = VSRA(r6, VBRADIX); r6 = VAND(r6, vbmask); r7 = VADD(r7, c);  
  c = VSRA(r7, VBRADIX); r7 = VAND(r7, vbmask);
  r0 = VMADD(r0, 0xAA, r0, VSHUF(c, 0x4E));

  r[0] = r0; r[1] = r1; r[2] = r2; r[3] = r3; 
  r[4] = r4; r[5] = r5; r[6] = r6; r[7] = r7;
}

// integer multiplication r = a * b (operand-scanning)
// len(r) = 2*VNWORDS; len(a) = len(b) = NOWRDS
// no carry propagation 
void mp_mul_4x2w(__m512i *r, const __m512i *a, const __m512i *b)
{
  __m512i a0 = a[0], a1 = a[1], a2 = a[2], a3 = a[3];
  __m512i a4 = a[4], a5 = a[5], a6 = a[6], a7 = a[7];
  __m512i b0 = b[0], b1 = b[1], b2 = b[2], b3 = b[3];
  __m512i b4 = b[4], b5 = b[5], b6 = b[6], b7 = b[7];
  __m512i r0,  r1,  r2,  r3,  r4,  r5,  r6,  r7,  r8,  r9;
  __m512i r10, r11, r12, r13, r14, r15, r16, r17, r18, r19;
  __m512i r20, r21, r22, tb;
  __m512i z0  = VZERO, z1  = VZERO, z2  = VZERO, z3  = VZERO, z4  = VZERO;
  __m512i z5  = VZERO, z6  = VZERO, z7  = VZERO, z8  = VZERO, z9  = VZERO;
  __m512i z10 = VZERO, z11 = VZERO, z12 = VZERO, z13 = VZERO, z14 = VZERO;
  __m512i z15 = VZERO, z16 = VZERO, z17 = VZERO, z18 = VZERO, z19 = VZERO;
  __m512i z20 = VZERO, z21 = VZERO, z22 = VZERO;
  __m512i h0  = VZERO, h1  = VZERO, h2  = VZERO, h3  = VZERO, h4  = VZERO;
  __m512i h5  = VZERO, h6  = VZERO, h7  = VZERO, h8  = VZERO, h9  = VZERO;
  __m512i h10 = VZERO, h11 = VZERO, h12 = VZERO, h13 = VZERO, h14 = VZERO;
  __m512i h15 = VZERO, h16 = VZERO, h17 = VZERO, h18 = VZERO, h19 = VZERO;
  __m512i h20 = VZERO, h21 = VZERO, h22 = VZERO;
  const __m512i vbmask = VSET1(VBMASK);

  // integer multiplication 
  tb = VSHUF(b0, 0x44); 
  z0 = VMACLO(z0, tb, a0); z1 = VMACLO(z1, tb, a1); z2 = VMACLO(z2, tb, a2); 
  z3 = VMACLO(z3, tb, a3); z4 = VMACLO(z4, tb, a4); z5 = VMACLO(z5, tb, a5);
  z6 = VMACLO(z6, tb, a6); z7 = VMACLO(z7, tb, a7);
  h0 = VMACHI(h0, tb, a0); h1 = VMACHI(h1, tb, a1); h2 = VMACHI(h2, tb, a2);
  h3 = VMACHI(h3, tb, a3); h4 = VMACHI(h4, tb, a4); h5 = VMACHI(h5, tb, a5);
  h6 = VMACHI(h6, tb, a6); h7 = VMACHI(h7, tb, a7);

  tb = VSHUF(b1, 0x44);
  z1 = VMACLO(z1, tb, a0); z2 = VMACLO(z2, tb, a1); z3 = VMACLO(z3, tb, a2);
  z4 = VMACLO(z4, tb, a3); z5 = VMACLO(z5, tb, a4); z6 = VMACLO(z6, tb, a5);
  z7 = VMACLO(z7, tb, a6); z8 = VMACLO(z8, tb, a7); 
  h1 = VMACHI(h1, tb, a0); h2 = VMACHI(h2, tb, a1); h3 = VMACHI(h3, tb, a2);
  h4 = VMACHI(h4, tb, a3); h5 = VMACHI(h5, tb, a4); h6 = VMACHI(h6, tb, a5);
  h7 = VMACHI(h7, tb, a6); h8 = VMACHI(h8, tb, a7);

  tb = VSHUF(b2, 0x44);
  z2 = VMACLO(z2, tb, a0); z3 = VMACLO(z3, tb, a1); z4 = VMACLO(z4, tb, a2);
  z5 = VMACLO(z5, tb, a3); z6 = VMACLO(z6, tb, a4); z7 = VMACLO(z7, tb, a5);
  z8 = VMACLO(z8, tb, a6); z9 = VMACLO(z9, tb, a7);
  h2 = VMACHI(h2, tb, a0); h3 = VMACHI(h3, tb, a1); h4 = VMACHI(h4, tb, a2);
  h5 = VMACHI(h5, tb, a3); h6 = VMACHI(h6, tb, a4); h7 = VMACHI(h7, tb, a5);
  h8 = VMACHI(h8, tb, a6); h9 = VMACHI(h9, tb, a7);

  tb = VSHUF(b3, 0x44);
  z3 = VMACLO(z3, tb, a0); z4  = VMACLO(z4,  tb, a1); z5 = VMACLO(z5, tb, a2);
  z6 = VMACLO(z6, tb, a3); z7  = VMACLO(z7,  tb, a4); z8 = VMACLO(z8, tb, a5);
  z9 = VMACLO(z9, tb, a6); z10 = VMACLO(z10, tb, a7);
  h3 = VMACHI(h3, tb, a0); h4  = VMACHI(h4,  tb, a1); h5 = VMACHI(h5, tb, a2);
  h6 = VMACHI(h6, tb, a3); h7  = VMACHI(h7,  tb, a4); h8 = VMACHI(h8, tb, a5);
  h9 = VMACHI(h9, tb, a6); h10 = VMACHI(h10, tb, a7);

  tb = VSHUF(b4, 0x44);
  z4  = VMACLO(z4,  tb, a0); z5  = VMACLO(z5,  tb, a1); z6 = VMACLO(z6, tb, a2);
  z7  = VMACLO(z7,  tb, a3); z8  = VMACLO(z8,  tb, a4); z9 = VMACLO(z9, tb, a5);
  z10 = VMACLO(z10, tb, a6); z11 = VMACLO(z11, tb, a7);
  h4  = VMACHI(h4,  tb, a0); h5  = VMACHI(h5,  tb, a1); h6 = VMACHI(h6, tb, a2);
  h7  = VMACHI(h7,  tb, a3); h8  = VMACHI(h8,  tb, a4); h9 = VMACHI(h9, tb, a5);
  h10 = VMACHI(h10, tb, a6); h11 = VMACHI(h11, tb, a7);

  tb = VSHUF(b5, 0x44);
  z5  = VMACLO(z5,  tb, a0); z6  = VMACLO(z6,  tb, a1); z7  = VMACLO(z7,  tb, a2);
  z8  = VMACLO(z8,  tb, a3); z9  = VMACLO(z9,  tb, a4); z10 = VMACLO(z10, tb, a5);
  z11 = VMACLO(z11, tb, a6); z12 = VMACLO(z12, tb, a7);
  h5  = VMACHI(h5,  tb, a0); h6  = VMACHI(h6,  tb, a1); h7  = VMACHI(h7,  tb, a2);
  h8  = VMACHI(h8,  tb, a3); h9  = VMACHI(h9,  tb, a4); h10 = VMACHI(h10, tb, a5);
  h11 = VMACHI(h11, tb, a6); h12 = VMACHI(h12, tb, a7);

  tb = VSHUF(b6, 0x44);
  z6  = VMACLO(z6,  tb, a0); z7  = VMACLO(z7,  tb, a1); z8  = VMACLO(z8,  tb, a2);
  z9  = VMACLO(z9,  tb, a3); z10 = VMACLO(z10, tb, a4); z11 = VMACLO(z11, tb, a5);
  z12 = VMACLO(z12, tb, a6); z13 = VMACLO(z13, tb, a7);
  h6  = VMACHI(h6,  tb, a0); h7  = VMACHI(h7,  tb, a1); h8  = VMACHI(h8,  tb, a2);
  h9  = VMACHI(h9,  tb, a3); h10 = VMACHI(h10, tb, a4); h11 = VMACHI(h11, tb, a5);
  h12 = VMACHI(h12, tb, a6); h13 = VMACHI(h13, tb, a7);

  tb = VSHUF(b7, 0x44);
  z7  = VMACLO(z7,  tb, a0); z8  = VMACLO(z8,  tb, a1); z9  = VMACLO(z9,  tb, a2);
  z10 = VMACLO(z10, tb, a3); z11 = VMACLO(z11, tb, a4); z12 = VMACLO(z12, tb, a5);
  z13 = VMACLO(z13, tb, a6); z14 = VMACLO(z14, tb, a7);
  h7  = VMACHI(h7,  tb, a0); h8  = VMACHI(h8,  tb, a1); h9  = VMACHI(h9,  tb, a2);
  h10 = VMACHI(h10, tb, a3); h11 = VMACHI(h11, tb, a4); h12 = VMACHI(h12, tb, a5);
  h13 = VMACHI(h13, tb, a6); h14 = VMACHI(h14, tb, a7); 

  tb = VSHUF(b0, 0xEE);
  z8  = VMACLO(z8,  tb, a0); z9  = VMACLO(z9,  tb, a1); z10 = VMACLO(z10, tb, a2);
  z11 = VMACLO(z11, tb, a3); z12 = VMACLO(z12, tb, a4); z13 = VMACLO(z13, tb, a5);
  z14 = VMACLO(z14, tb, a6); z15 = VMACLO(z15, tb, a7);
  h8  = VMACHI(h8,  tb, a0); h9  = VMACHI(h9,  tb, a1); h10 = VMACHI(h10, tb, a2);
  h11 = VMACHI(h11, tb, a3); h12 = VMACHI(h12, tb, a4); h13 = VMACHI(h13, tb, a5);
  h14 = VMACHI(h14, tb, a6); h15 = VMACHI(h15, tb, a7); 

  tb = VSHUF(b1, 0xEE);
  z9  = VMACLO(z9,  tb, a0); z10 = VMACLO(z10, tb, a1); z11 = VMACLO(z11, tb, a2);
  z12 = VMACLO(z12, tb, a3); z13 = VMACLO(z13, tb, a4); z14 = VMACLO(z14, tb, a5);
  z15 = VMACLO(z15, tb, a6); z16 = VMACLO(z16, tb, a7);
  h9  = VMACHI(h9,  tb, a0); h10 = VMACHI(h10, tb, a1); h11 = VMACHI(h11, tb, a2);
  h12 = VMACHI(h12, tb, a3); h13 = VMACHI(h13, tb, a4); h14 = VMACHI(h14, tb, a5);
  h15 = VMACHI(h15, tb, a6); h16 = VMACHI(h16, tb, a7); 

  tb = VSHUF(b2, 0xEE);
  z10 = VMACLO(z10, tb, a0); z11 = VMACLO(z11, tb, a1); z12 = VMACLO(z12, tb, a2);
  z13 = VMACLO(z13, tb, a3); z14 = VMACLO(z14, tb, a4); z15 = VMACLO(z15, tb, a5);
  z16 = VMACLO(z16, tb, a6); z17 = VMACLO(z17, tb, a7);
  h10 = VMACHI(h10, tb, a0); h11 = VMACHI(h11, tb, a1); h12 = VMACHI(h12, tb, a2);
  h13 = VMACHI(h13, tb, a3); h14 = VMACHI(h14, tb, a4); h15 = VMACHI(h15, tb, a5);
  h16 = VMACHI(h16, tb, a6); h17 = VMACHI(h17, tb, a7); 

  tb = VSHUF(b3, 0xEE);
  z11 = VMACLO(z11, tb, a0); z12 = VMACLO(z12, tb, a1); z13 = VMACLO(z13, tb, a2);
  z14 = VMACLO(z14, tb, a3); z15 = VMACLO(z15, tb, a4); z16 = VMACLO(z16, tb, a5);
  z17 = VMACLO(z17, tb, a6); z18 = VMACLO(z18, tb, a7);
  h11 = VMACHI(h11, tb, a0); h12 = VMACHI(h12, tb, a1); h13 = VMACHI(h13, tb, a2);
  h14 = VMACHI(h14, tb, a3); h15 = VMACHI(h15, tb, a4); h16 = VMACHI(h16, tb, a5);
  h17 = VMACHI(h17, tb, a6); h18 = VMACHI(h18, tb, a7); 

  tb = VSHUF(b4, 0xEE);
  z12 = VMACLO(z12, tb, a0); z13 = VMACLO(z13, tb, a1); z14 = VMACLO(z14, tb, a2);
  z15 = VMACLO(z15, tb, a3); z16 = VMACLO(z16, tb, a4); z17 = VMACLO(z17, tb, a5);
  z18 = VMACLO(z18, tb, a6); z19 = VMACLO(z19, tb, a7);
  h12 = VMACHI(h12, tb, a0); h13 = VMACHI(h13, tb, a1); h14 = VMACHI(h14, tb, a2);
  h15 = VMACHI(h15, tb, a3); h16 = VMACHI(h16, tb, a4); h17 = VMACHI(h17, tb, a5);
  h18 = VMACHI(h18, tb, a6); h19 = VMACHI(h19, tb, a7);

  tb = VSHUF(b5, 0xEE);
  z13 = VMACLO(z13, tb, a0); z14 = VMACLO(z14, tb, a1); z15 = VMACLO(z15, tb, a2);
  z16 = VMACLO(z16, tb, a3); z17 = VMACLO(z17, tb, a4); z18 = VMACLO(z18, tb, a5);
  z19 = VMACLO(z19, tb, a6); z20 = VMACLO(z20, tb, a7);
  h13 = VMACHI(h13, tb, a0); h14 = VMACHI(h14, tb, a1); h15 = VMACHI(h15, tb, a2);
  h16 = VMACHI(h16, tb, a3); h17 = VMACHI(h17, tb, a4); h18 = VMACHI(h18, tb, a5);
  h19 = VMACHI(h19, tb, a6); h20 = VMACHI(h20, tb, a7);

  tb = VSHUF(b6, 0xEE);
  z14 = VMACLO(z14, tb, a0); z15 = VMACLO(z15, tb, a1); z16 = VMACLO(z16, tb, a2);
  z17 = VMACLO(z17, tb, a3); z18 = VMACLO(z18, tb, a4); z19 = VMACLO(z19, tb, a5);
  z20 = VMACLO(z20, tb, a6); z21 = VMACLO(z21, tb, a7);
  h14 = VMACHI(h14, tb, a0); h15 = VMACHI(h15, tb, a1); h16 = VMACHI(h16, tb, a2);
  h17 = VMACHI(h17, tb, a3); h18 = VMACHI(h18, tb, a4); h19 = VMACHI(h19, tb, a5);
  h20 = VMACHI(h20, tb, a6); h21 = VMACHI(h21, tb, a7);    
  
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
  z17 = VADD(z17, VADD(h16, h16)); z18 = VADD(z18, VADD(h17, h17));
  z19 = VADD(z19, VADD(h18, h18)); z20 = VADD(z20, VADD(h19, h19));
  z21 = VADD(z21, VADD(h20, h20)); z22 = VADD(z22, VADD(h21, h21));

  r[0]  = z0;  r[1]  = z1;  r[2]  = z2;  r[3]  = z3;  r[4]  = z4;
  r[5]  = z5;  r[6]  = z6;  r[7]  = z7;  r[8]  = z8;  r[9]  = z9;
  r[10] = z10; r[11] = z11; r[12] = z12; r[13] = z13; r[14] = z14;
  r[15] = z15; r[16] = z16; r[17] = z17; r[18] = z18; r[19] = z19;
  r[20] = z20; r[21] = z21; r[22] = z22;
}

// Montgomery reduction r = a * R^-1 mod 2p, where R = 2^765 (operand-scanning)
void rdc_mont_4x2w(__m512i *r, const __m512i *a)
{
  __m512i z0  = a[0],  z1  = a[1],  z2  = a[2],  z3  = a[3],  z4  = a[4];
  __m512i z5  = a[5],  z6  = a[6],  z7  = a[7],  z8  = a[8],  z9  = a[9];
  __m512i z10 = a[10], z11 = a[11], z12 = a[12], z13 = a[13], z14 = a[14];
  __m512i z15 = a[15], z16 = a[16], z17 = a[17], z18 = a[18], z19 = a[19];
  __m512i z20 = a[20], z21 = a[21], z22 = a[22];
  __m512i h0  = VZERO, h1  = VZERO, h2  = VZERO, h3  = VZERO, h4  = VZERO;
  __m512i h5  = VZERO, h6  = VZERO, h7  = VZERO, h8  = VZERO, h9  = VZERO;
  __m512i h10 = VZERO, h11 = VZERO, h12 = VZERO, h13 = VZERO, h14 = VZERO; 
  __m512i h15 = VZERO, h16 = VZERO, h17 = VZERO, h18 = VZERO, h19 = VZERO;
  __m512i h20 = VZERO, h21 = VZERO, h22 = VZERO;
  const __m512i vp7  = VSET(vp751p1[8],  vp751p1[7],  vp751p1[8],  vp751p1[7],  vp751p1[8],  vp751p1[7],  vp751p1[8],  vp751p1[7]);
  const __m512i vp9  = VSET(vp751p1[10], vp751p1[9],  vp751p1[10], vp751p1[9],  vp751p1[10], vp751p1[9],  vp751p1[10], vp751p1[9]);
  const __m512i vp11 = VSET(vp751p1[12], vp751p1[11], vp751p1[12], vp751p1[11], vp751p1[12], vp751p1[11], vp751p1[12], vp751p1[11]);
  const __m512i vp13 = VSET(vp751p1[14], vp751p1[13], vp751p1[14], vp751p1[13], vp751p1[14], vp751p1[13], vp751p1[14], vp751p1[13]);
  const __m512i vbmask = VSET1(VBMASK), vzero = VZERO;
  __m512i r0, r1, r2, r3, r4, r5, r6, r7;
  __m512i c, u, t0, t1, t2, t3, t4, t5, t6, t7;

  c = VSRA(z0, VBRADIX); z0 = VAND(z0, vbmask); z1 = VADD(z1, c);
  u = VSHUF(z0, 0x44);
  t0 = VMACLO(vzero, u, vp7);  // t0 = lo(u*p8  | u*p7)
  t1 = VMACHI(vzero, u, vp7);  // t1 = hi(u*p8  | u*p7)
  t2 = VMACLO(vzero, u, vp9);  // t2 = lo(u*p10 | u*p9)
  t3 = VMACHI(vzero, u, vp9);  // t3 = hi(u*p10 | u*p9)
  t4 = VMACLO(vzero, u, vp11); // t4 = lo(u*p12 | u*p11)
  t5 = VMACHI(vzero, u, vp11); // t5 = hi(u*p12 | u*p11)
  t6 = VMACLO(vzero, u, vp13); // t6 = lo(u*p14 | u*p13)
  t7 = VMACHI(vzero, u, vp13); // t7 = hi(u*p14 | u*p13)
  z7  = VMADD(z7,  0x55, z7, t0);               // z7  = z15 | z7+u*p7
  z8  = VMADD(z8,  0x55, z8, VSHUF(t0, 0x4E));  // z8  = z16 | z8+u*p8
  h7  = VMADD(h7,  0x55, h7, t1);               // h7  = h15 | h7+u*p7
  h8  = VMADD(h8,  0x55, h8, VSHUF(t1, 0x4E));  // h8  = h16 | h8+u*p8
  z9  = VMADD(z9,  0x55, z9, t2);               // z9  = z17 | z9+u*p9
  z10 = VMADD(z10, 0x55, z10, VSHUF(t2, 0x4E)); // z10 = z18 | z10+u*p10
  h9  = VMADD(h9,  0x55, h9, t3);               // h9  = h17 | h9+u*p9
  h10 = VMADD(h10, 0x55, h10, VSHUF(t3, 0x4E)); // h10 = h18 | h10+u*p10
  z3  = VMADD(z3,  0xAA, z3, VSHUF(t4, 0x4E));  // z3  = z11+u*p11 | z3
  z4  = VMADD(z4,  0xAA, z4, t4);               // z4  = z12+u*p12 | z4
  h3  = VMADD(h3,  0xAA, h3, VSHUF(t5, 0x4E));  // h3  = h11+u*p11 | h3 
  h4  = VMADD(h4,  0xAA, h4, t5);               // h4  = h12+u*p12 | h4
  z5  = VMADD(z5,  0xAA, z5, VSHUF(t6, 0x4E));  // z5  = z13+u*p13 | z5
  z6  = VMADD(z6,  0xAA, z6, t6);               // z6  = z14+u*p14 | z6 
  h5  = VMADD(h5,  0xAA, h5, VSHUF(t7, 0x4E));  // h5  = h13+u*p11 | h5 
  h6  = VMADD(h6,  0xAA, h6, t7);               // h6  = h14+u*p12 | h6
  z8  = VMADD(z8,  0x55, z8, VSHUF(z0, 0x4E));

  c = VSRA(z1, VBRADIX); z1 = VAND(z1, vbmask); z2 = VADD(z2, c);
  u = VSHUF(z1, 0x44);
  t0 = VMACLO(vzero, u, vp7);  t1 = VMACHI(vzero, u, vp7);
  t2 = VMACLO(vzero, u, vp9);  t3 = VMACHI(vzero, u, vp9);
  t4 = VMACLO(vzero, u, vp11); t5 = VMACHI(vzero, u, vp11);
  t6 = VMACLO(vzero, u, vp13); t7 = VMACHI(vzero, u, vp13);
  z8 = VMADD(z8, 0x55, z8, t0); z9 = VMADD(z9, 0x55, z9, VSHUF(t0, 0x4E));
  h8 = VMADD(h8, 0x55, h8, t1); h9 = VMADD(h9, 0x55, h9, VSHUF(t1, 0x4E));
  z10 = VMADD(z10, 0x55, z10, t2); z11 = VMADD(z11, 0x55, z11, VSHUF(t2, 0x4E));
  h10 = VMADD(h10, 0x55, h10, t3); h11 = VMADD(h11, 0x55, h11, VSHUF(t3, 0x4E));
  z4 = VMADD(z4, 0xAA, z4, VSHUF(t4, 0x4E)); z5 = VMADD(z5, 0xAA, z5, t4);
  h4 = VMADD(h4, 0xAA, h4, VSHUF(t5, 0x4E)); h5 = VMADD(h5, 0xAA, h5, t5);
  z6 = VMADD(z6, 0xAA, z6, VSHUF(t6, 0x4E)); z7 = VMADD(z7, 0xAA, z7, t6);
  h6 = VMADD(h6, 0xAA, h6, VSHUF(t7, 0x4E)); h7 = VMADD(h7, 0xAA, h7, t7);
  z9 = VMADD(z9, 0x55, z9, VSHUF(z1, 0x4E));

  c = VSRA(z2, VBRADIX); z2 = VAND(z2, vbmask); z3 = VADD(z3, c);
  u = VSHUF(z2, 0x44);
  t0 = VMACLO(vzero, u, vp7);  t1 = VMACHI(vzero, u, vp7);
  t2 = VMACLO(vzero, u, vp9);  t3 = VMACHI(vzero, u, vp9);
  t4 = VMACLO(vzero, u, vp11); t5 = VMACHI(vzero, u, vp11);
  t6 = VMACLO(vzero, u, vp13); t7 = VMACHI(vzero, u, vp13);
  z9 = VMADD(z9, 0x55, z9, t0); z10 = VMADD(z10, 0x55, z10, VSHUF(t0, 0x4E));
  h9 = VMADD(h9, 0x55, h9, t1); h10 = VMADD(h10, 0x55, h10, VSHUF(t1, 0x4E));
  z11 = VMADD(z11, 0x55, z11, t2); z12 = VMADD(z12, 0x55, z12, VSHUF(t2, 0x4E));
  h11 = VMADD(h11, 0x55, h11, t3); h12 = VMADD(h12, 0x55, h12, VSHUF(t3, 0x4E));
  z5 = VMADD(z5, 0xAA, z5, VSHUF(t4, 0x4E)); z6 = VMADD(z6, 0xAA, z6, t4);
  h5 = VMADD(h5, 0xAA, h5, VSHUF(t5, 0x4E)); h6 = VMADD(h6, 0xAA, h6, t5);
  z7 = VMADD(z7, 0xAA, z7, VSHUF(t6, 0x4E)); z8 = VMADD(z8, 0xAA, z8, t6);
  h7 = VMADD(h7, 0xAA, h7, VSHUF(t7, 0x4E)); h8 = VMADD(h8, 0xAA, h8, t7);
  z10 = VMADD(z10, 0x55, z10, VSHUF(z2, 0x4E));

  c = VSRA(z3, VBRADIX); z3 = VAND(z3, vbmask); z4 = VADD(z4, c);
  u = VSHUF(z3, 0x44);
  t0 = VMACLO(vzero, u, vp7);  t1 = VMACHI(vzero, u, vp7);
  t2 = VMACLO(vzero, u, vp9);  t3 = VMACHI(vzero, u, vp9);
  t4 = VMACLO(vzero, u, vp11); t5 = VMACHI(vzero, u, vp11);
  t6 = VMACLO(vzero, u, vp13); t7 = VMACHI(vzero, u, vp13);
  z10 = VMADD(z10, 0x55, z10, t0); z11 = VMADD(z11, 0x55, z11, VSHUF(t0, 0x4E));
  h10 = VMADD(h10, 0x55, h10, t1); h11 = VMADD(h11, 0x55, h11, VSHUF(t1, 0x4E));
  z12 = VMADD(z12, 0x55, z12, t2); z13 = VMADD(z13, 0x55, z13, VSHUF(t2, 0x4E));
  h12 = VMADD(h12, 0x55, h12, t3); h13 = VMADD(h13, 0x55, h13, VSHUF(t3, 0x4E));
  z6 = VMADD(z6, 0xAA, z6, VSHUF(t4, 0x4E)); z7 = VMADD(z7, 0xAA, z7, t4);
  h6 = VMADD(h6, 0xAA, h6, VSHUF(t5, 0x4E)); h7 = VMADD(h7, 0xAA, h7, t5);
  z8 = VMADD(z8, 0xAA, z8, VSHUF(t6, 0x4E)); z9 = VMADD(z9, 0xAA, z9, t6);
  h8 = VMADD(h8, 0xAA, h8, VSHUF(t7, 0x4E)); h9 = VMADD(h9, 0xAA, h9, t7);
  z11 = VMADD(z11, 0x55, z11, VSHUF(z3, 0x4E));
  z4 = VADD(z4, VADD(h3, h3));

  c = VSRA(z4, VBRADIX); z4 = VAND(z4, vbmask); z5 = VADD(z5, c);
  u = VSHUF(z4, 0x44);
  t0 = VMACLO(vzero, u, vp7);  t1 = VMACHI(vzero, u, vp7);
  t2 = VMACLO(vzero, u, vp9);  t3 = VMACHI(vzero, u, vp9);
  t4 = VMACLO(vzero, u, vp11); t5 = VMACHI(vzero, u, vp11);
  t6 = VMACLO(vzero, u, vp13); t7 = VMACHI(vzero, u, vp13);
  z11 = VMADD(z11, 0x55, z11, t0); z12 = VMADD(z12, 0x55, z12, VSHUF(t0, 0x4E));
  h11 = VMADD(h11, 0x55, h11, t1); h12 = VMADD(h12, 0x55, h12, VSHUF(t1, 0x4E));
  z13 = VMADD(z13, 0x55, z13, t2); z14 = VMADD(z14, 0x55, z14, VSHUF(t2, 0x4E));
  h13 = VMADD(h13, 0x55, h13, t3); h14 = VMADD(h14, 0x55, h14, VSHUF(t3, 0x4E));
  z7 = VMADD(z7, 0xAA, z7, VSHUF(t4, 0x4E)); z8 = VMADD(z8, 0xAA, z8, t4);
  h7 = VMADD(h7, 0xAA, h7, VSHUF(t5, 0x4E)); h8 = VMADD(h8, 0xAA, h8, t5);
  z9 = VMADD(z9, 0xAA, z9, VSHUF(t6, 0x4E)); z10 = VMADD(z10, 0xAA, z10, t6);
  h9 = VMADD(h9, 0xAA, h9, VSHUF(t7, 0x4E)); h10 = VMADD(h10, 0xAA, h10, t7);
  z12 = VMADD(z12, 0x55, z12, VSHUF(z4, 0x4E));
  z5 = VADD(z5, VADD(h4, h4));

  c = VSRA(z5, VBRADIX); z5 = VAND(z5, vbmask); z6 = VADD(z6, c);
  u = VSHUF(z5, 0x44);
  t0 = VMACLO(vzero, u, vp7);  t1 = VMACHI(vzero, u, vp7);
  t2 = VMACLO(vzero, u, vp9);  t3 = VMACHI(vzero, u, vp9);
  t4 = VMACLO(vzero, u, vp11); t5 = VMACHI(vzero, u, vp11);
  t6 = VMACLO(vzero, u, vp13); t7 = VMACHI(vzero, u, vp13);
  z12 = VMADD(z12, 0x55, z12, t0); z13 = VMADD(z13, 0x55, z13, VSHUF(t0, 0x4E));
  h12 = VMADD(h12, 0x55, h12, t1); h13 = VMADD(h13, 0x55, h13, VSHUF(t1, 0x4E));
  z14 = VMADD(z14, 0x55, z14, t2); z15 = VMADD(z15, 0x55, z15, VSHUF(t2, 0x4E));
  h14 = VMADD(h14, 0x55, h14, t3); h15 = VMADD(h15, 0x55, h15, VSHUF(t3, 0x4E));
  z8 = VMADD(z8, 0xAA, z8, VSHUF(t4, 0x4E)); z9 = VMADD(z9, 0xAA, z9, t4);
  h8 = VMADD(h8, 0xAA, h8, VSHUF(t5, 0x4E)); h9 = VMADD(h9, 0xAA, h9, t5);
  z10 = VMADD(z10, 0xAA, z10, VSHUF(t6, 0x4E)); z11 = VMADD(z11, 0xAA, z11, t6);
  h10 = VMADD(h10, 0xAA, h10, VSHUF(t7, 0x4E)); h11 = VMADD(h11, 0xAA, h11, t7);
  z13 = VMADD(z13, 0x55, z13, VSHUF(z5, 0x4E));
  z6 = VADD(z6, VADD(h5, h5));

  c = VSRA(z6, VBRADIX); z6 = VAND(z6, vbmask); z7 = VADD(z7, c);
  u = VSHUF(z6, 0x44);
  t0 = VMACLO(vzero, u, vp7);  t1 = VMACHI(vzero, u, vp7);
  t2 = VMACLO(vzero, u, vp9);  t3 = VMACHI(vzero, u, vp9);
  t4 = VMACLO(vzero, u, vp11); t5 = VMACHI(vzero, u, vp11);
  t6 = VMACLO(vzero, u, vp13); t7 = VMACHI(vzero, u, vp13);
  z13 = VMADD(z13, 0x55, z13, t0); z14 = VMADD(z14, 0x55, z14, VSHUF(t0, 0x4E));
  h13 = VMADD(h13, 0x55, h13, t1); h14 = VMADD(h14, 0x55, h14, VSHUF(t1, 0x4E));
  z15 = VMADD(z15, 0x55, z15, t2); z16 = VMADD(z16, 0x55, z16, VSHUF(t2, 0x4E));
  h15 = VMADD(h15, 0x55, h15, t3); h16 = VMADD(h16, 0x55, h16, VSHUF(t3, 0x4E));
  z9 = VMADD(z9, 0xAA, z9, VSHUF(t4, 0x4E)); z10 = VMADD(z10, 0xAA, z10, t4);
  h9 = VMADD(h9, 0xAA, h9, VSHUF(t5, 0x4E)); h10 = VMADD(h10, 0xAA, h10, t5);
  z11 = VMADD(z11, 0xAA, z11, VSHUF(t6, 0x4E)); z12 = VMADD(z12, 0xAA, z12, t6);
  h11 = VMADD(h11, 0xAA, h11, VSHUF(t7, 0x4E)); h12 = VMADD(h12, 0xAA, h12, t7);
  z14 = VMADD(z14, 0x55, z14, VSHUF(z6, 0x4E));
  z7 = VADD(z7, VADD(h6, h6));

  c = VSRA(z7, VBRADIX); z7 = VAND(z7, vbmask); z8 = VADD(z8, c);
  u = VSHUF(z7, 0x44);
  t0 = VMACLO(vzero, u, vp7);  t1 = VMACHI(vzero, u, vp7);
  t2 = VMACLO(vzero, u, vp9);  t3 = VMACHI(vzero, u, vp9);
  t4 = VMACLO(vzero, u, vp11); t5 = VMACHI(vzero, u, vp11);
  t6 = VMACLO(vzero, u, vp13); t7 = VMACHI(vzero, u, vp13);
  z14 = VMADD(z14, 0x55, z14, t0); z15 = VMADD(z15, 0x55, z15, VSHUF(t0, 0x4E));
  h14 = VMADD(h14, 0x55, h14, t1); h15 = VMADD(h15, 0x55, h15, VSHUF(t1, 0x4E));
  z16 = VMADD(z16, 0x55, z16, t2); z17 = VMADD(z17, 0x55, z17, VSHUF(t2, 0x4E));
  h16 = VMADD(h16, 0x55, h16, t3); h17 = VMADD(h17, 0x55, h17, VSHUF(t3, 0x4E));
  z10 = VMADD(z10, 0xAA, z10, VSHUF(t4, 0x4E)); z11 = VMADD(z11, 0xAA, z11, t4);
  h10 = VMADD(h10, 0xAA, h10, VSHUF(t5, 0x4E)); h11 = VMADD(h11, 0xAA, h11, t5);
  z12 = VMADD(z12, 0xAA, z12, VSHUF(t6, 0x4E)); z13 = VMADD(z13, 0xAA, z13, t6);
  h12 = VMADD(h12, 0xAA, h12, VSHUF(t7, 0x4E)); h13 = VMADD(h13, 0xAA, h13, t7);
  z15 = VMADD(z15, 0x55, z15, VSHUF(z7, 0x4E));
  z8 = VADD(z8, VADD(h7, h7));

  c = VSRA(z8, VBRADIX); z8 = VAND(z8, vbmask); z9 = VADD(z9, c);
  u = VSHUF(z8, 0x44);
  t0 = VMACLO(vzero, u, vp7);  t1 = VMACHI(vzero, u, vp7);
  t2 = VMACLO(vzero, u, vp9);  t3 = VMACHI(vzero, u, vp9);
  t4 = VMACLO(vzero, u, vp11); t5 = VMACHI(vzero, u, vp11);
  t6 = VMACLO(vzero, u, vp13); t7 = VMACHI(vzero, u, vp13);
  z15 = VMADD(z15, 0x55, z15, t0); z16 = VMADD(z16, 0x55, z16, VSHUF(t0, 0x4E));
  h15 = VMADD(h15, 0x55, h15, t1); h16 = VMADD(h16, 0x55, h16, VSHUF(t1, 0x4E));
  z17 = VMADD(z17, 0x55, z17, t2); z18 = VMADD(z18, 0x55, z18, VSHUF(t2, 0x4E));
  h17 = VMADD(h17, 0x55, h17, t3); h18 = VMADD(h18, 0x55, h18, VSHUF(t3, 0x4E));
  z11 = VMADD(z11, 0xAA, z11, VSHUF(t4, 0x4E)); z12 = VMADD(z12, 0xAA, z12, t4);
  h11 = VMADD(h11, 0xAA, h11, VSHUF(t5, 0x4E)); h12 = VMADD(h12, 0xAA, h12, t5);
  z13 = VMADD(z13, 0xAA, z13, VSHUF(t6, 0x4E)); z14 = VMADD(z14, 0xAA, z14, t6);
  h13 = VMADD(h13, 0xAA, h13, VSHUF(t7, 0x4E)); h14 = VMADD(h14, 0xAA, h14, t7);
  z16 = VMADD(z16, 0x55, z16, VSHUF(z8, 0x4E));
  z9 = VADD(z9, VADD(h8, h8));

  c = VSRA(z9, VBRADIX); z9 = VAND(z9, vbmask); z10 = VADD(z10, c);
  u = VSHUF(z9, 0x44);
  t0 = VMACLO(vzero, u, vp7);  t1 = VMACHI(vzero, u, vp7);
  t2 = VMACLO(vzero, u, vp9);  t3 = VMACHI(vzero, u, vp9);
  t4 = VMACLO(vzero, u, vp11); t5 = VMACHI(vzero, u, vp11);
  t6 = VMACLO(vzero, u, vp13); t7 = VMACHI(vzero, u, vp13);
  z16 = VMADD(z16, 0x55, z16, t0); z17 = VMADD(z17, 0x55, z17, VSHUF(t0, 0x4E));
  h16 = VMADD(h16, 0x55, h16, t1); h17 = VMADD(h17, 0x55, h17, VSHUF(t1, 0x4E));
  z18 = VMADD(z18, 0x55, z18, t2); z19 = VMADD(z19, 0x55, z19, VSHUF(t2, 0x4E));
  h18 = VMADD(h18, 0x55, h18, t3); h19 = VMADD(h19, 0x55, h19, VSHUF(t3, 0x4E));
  z12 = VMADD(z12, 0xAA, z12, VSHUF(t4, 0x4E)); z13 = VMADD(z13, 0xAA, z13, t4);
  h12 = VMADD(h12, 0xAA, h12, VSHUF(t5, 0x4E)); h13 = VMADD(h13, 0xAA, h13, t5);
  z14 = VMADD(z14, 0xAA, z14, VSHUF(t6, 0x4E)); z15 = VMADD(z15, 0xAA, z15, t6);
  h14 = VMADD(h14, 0xAA, h14, VSHUF(t7, 0x4E)); h15 = VMADD(h15, 0xAA, h15, t7);
  z17 = VMADD(z17, 0x55, z17, VSHUF(z9, 0x4E));
  z10 = VADD(z10, VADD(h9, h9));

  c = VSRA(z10, VBRADIX); z10 = VAND(z10, vbmask); z11 = VADD(z11, c);
  u = VSHUF(z10, 0x44);
  t0 = VMACLO(vzero, u, vp7);  t1 = VMACHI(vzero, u, vp7);
  t2 = VMACLO(vzero, u, vp9);  t3 = VMACHI(vzero, u, vp9);
  t4 = VMACLO(vzero, u, vp11); t5 = VMACHI(vzero, u, vp11);
  t6 = VMACLO(vzero, u, vp13); t7 = VMACHI(vzero, u, vp13);
  z17 = VMADD(z17, 0x55, z17, t0); z18 = VMADD(z18, 0x55, z18, VSHUF(t0, 0x4E));
  h17 = VMADD(h17, 0x55, h17, t1); h18 = VMADD(h18, 0x55, h18, VSHUF(t1, 0x4E));
  z19 = VMADD(z19, 0x55, z19, t2); z20 = VMADD(z20, 0x55, z20, VSHUF(t2, 0x4E));
  h19 = VMADD(h19, 0x55, h19, t3); h20 = VMADD(h20, 0x55, h20, VSHUF(t3, 0x4E));
  z13 = VMADD(z13, 0xAA, z13, VSHUF(t4, 0x4E)); z14 = VMADD(z14, 0xAA, z14, t4);
  h13 = VMADD(h13, 0xAA, h13, VSHUF(t5, 0x4E)); h14 = VMADD(h14, 0xAA, h14, t5);
  z15 = VMADD(z15, 0xAA, z15, VSHUF(t6, 0x4E)); z16 = VMADD(z16, 0xAA, z16, t6);
  h15 = VMADD(h15, 0xAA, h15, VSHUF(t7, 0x4E)); h16 = VMADD(h16, 0xAA, h16, t7);
  z18 = VMADD(z18, 0x55, z18, VSHUF(z10, 0x4E));
  z11 = VADD(z11, VADD(h10, h10));

  c = VSRA(z11, VBRADIX); z11 = VAND(z11, vbmask); z12 = VADD(z12, c);
  u = VSHUF(z11, 0x44);
  t0 = VMACLO(vzero, u, vp7);  t1 = VMACHI(vzero, u, vp7);
  t2 = VMACLO(vzero, u, vp9);  t3 = VMACHI(vzero, u, vp9);
  t4 = VMACLO(vzero, u, vp11); t5 = VMACHI(vzero, u, vp11);
  t6 = VMACLO(vzero, u, vp13); t7 = VMACHI(vzero, u, vp13);
  z18 = VMADD(z18, 0x55, z18, t0); z19 = VMADD(z19, 0x55, z19, VSHUF(t0, 0x4E));
  h18 = VMADD(h18, 0x55, h18, t1); h19 = VMADD(h19, 0x55, h19, VSHUF(t1, 0x4E));
  z20 = VMADD(z20, 0x55, z20, t2); z21 = VMADD(z21, 0x55, z21, VSHUF(t2, 0x4E));
  h20 = VMADD(h20, 0x55, h20, t3); h21 = VMADD(h21, 0x55, h21, VSHUF(t3, 0x4E));
  z14 = VMADD(z14, 0xAA, z14, VSHUF(t4, 0x4E)); z15 = VMADD(z15, 0xAA, z15, t4);
  h14 = VMADD(h14, 0xAA, h14, VSHUF(t5, 0x4E)); h15 = VMADD(h15, 0xAA, h15, t5);
  z16 = VMADD(z16, 0xAA, z16, VSHUF(t6, 0x4E)); z17 = VMADD(z17, 0xAA, z17, t6);
  h16 = VMADD(h16, 0xAA, h16, VSHUF(t7, 0x4E)); h17 = VMADD(h17, 0xAA, h17, t7);
  z19 = VMADD(z19, 0x55, z19, VSHUF(z11, 0x4E));
  z12 = VADD(z12, VADD(h11, h11));

  c = VSRA(z12, VBRADIX); z12 = VAND(z12, vbmask); z13 = VADD(z13, c);
  u = VSHUF(z12, 0x44);
  t0 = VMACLO(vzero, u, vp7);  t1 = VMACHI(vzero, u, vp7);
  t2 = VMACLO(vzero, u, vp9);  t3 = VMACHI(vzero, u, vp9);
  t4 = VMACLO(vzero, u, vp11); t5 = VMACHI(vzero, u, vp11);
  t6 = VMACLO(vzero, u, vp13); t7 = VMACHI(vzero, u, vp13);
  z19 = VMADD(z19, 0x55, z19, t0); z20 = VMADD(z20, 0x55, z20, VSHUF(t0, 0x4E));
  h19 = VMADD(h19, 0x55, h19, t1); h20 = VMADD(h20, 0x55, h20, VSHUF(t1, 0x4E));
  z21 = VMADD(z21, 0x55, z21, t2); z22 = VMADD(z22, 0x55, z22, VSHUF(t2, 0x4E));
  h21 = VMADD(h21, 0x55, h21, t3); h14 = VMADD(h14, 0xAA, h14, t3);
  z15 = VMADD(z15, 0xAA, z15, VSHUF(t4, 0x4E)); z16 = VMADD(z16, 0xAA, z16, t4);
  h15 = VMADD(h15, 0xAA, h15, VSHUF(t5, 0x4E)); h16 = VMADD(h16, 0xAA, h16, t5);
  z17 = VMADD(z17, 0xAA, z17, VSHUF(t6, 0x4E)); z18 = VMADD(z18, 0xAA, z18, t6);
  h17 = VMADD(h17, 0xAA, h17, VSHUF(t7, 0x4E)); h18 = VMADD(h18, 0xAA, h18, t7);
  z20 = VMADD(z20, 0x55, z20, VSHUF(z12, 0x4E));
  z13 = VADD(z13, VADD(h12, h12));

  c = VSRA(z13, VBRADIX); z13 = VAND(z13, vbmask); z14 = VADD(z14, c);
  u = VSHUF(z13, 0x44);
  t0 = VMACLO(vzero, u, vp7);  t1 = VMACHI(vzero, u, vp7);
  t2 = VMACLO(vzero, u, vp9);  t3 = VMACHI(vzero, u, vp9);
  t4 = VMACLO(vzero, u, vp11); t5 = VMACHI(vzero, u, vp11);
  t6 = VMACLO(vzero, u, vp13); t7 = VMACHI(vzero, u, vp13);
  z20 = VMADD(z20, 0x55, z20, t0); z21 = VMADD(z21, 0x55, z21, VSHUF(t0, 0x4E));
  h20 = VMADD(h20, 0x55, h20, t1); h21 = VMADD(h21, 0x55, h21, VSHUF(t1, 0x4E));
  z22 = VMADD(z22, 0x55, z22, t2); z15 = VMADD(z15, 0xAA, z15, t2);
  h14 = VMADD(h14, 0xAA, h14, VSHUF(t3, 0x4E)); h15 = VMADD(h15, 0xAA, h15, t3);
  z16 = VMADD(z16, 0xAA, z16, VSHUF(t4, 0x4E)); z17 = VMADD(z17, 0xAA, z17, t4);
  h16 = VMADD(h16, 0xAA, h16, VSHUF(t5, 0x4E)); h17 = VMADD(h17, 0xAA, h17, t5);
  z18 = VMADD(z18, 0xAA, z18, VSHUF(t6, 0x4E)); z19 = VMADD(z19, 0xAA, z19, t6);
  h18 = VMADD(h18, 0xAA, h18, VSHUF(t7, 0x4E)); h19 = VMADD(h19, 0xAA, h19, t7);
  z21 = VMADD(z21, 0x55, z21, VSHUF(z13, 0x4E));
  z14 = VADD(z14, VADD(h13, h13));

  c = VSRA(z14, VBRADIX); z14 = VAND(z14, vbmask); z15 = VADD(z15, c);
  u = VSHUF(z14, 0x44);
  t0 = VMACLO(vzero, u, vp7);  t1 = VMACHI(vzero, u, vp7);
  t2 = VMACLO(vzero, u, vp9);  t3 = VMACHI(vzero, u, vp9);
  t4 = VMACLO(vzero, u, vp11); t5 = VMACHI(vzero, u, vp11);
  t6 = VMACLO(vzero, u, vp13); t7 = VMACHI(vzero, u, vp13);
  z21 = VMADD(z21, 0x55, z21, t0); z22 = VMADD(z22, 0x55, z22, VSHUF(t0, 0x4E));
  h21 = VMADD(h21, 0x55, h21, t1); h14 = VMADD(h14, 0xAA, h14, t1);
  z15 = VMADD(z15, 0xAA, z15, VSHUF(t2, 0x4E)); z16 = VMADD(z16, 0xAA, z16, t2);
  h15 = VMADD(h15, 0xAA, h15, VSHUF(t3, 0x4E)); h16 = VMADD(h16, 0xAA, h16, t3);
  z17 = VMADD(z17, 0xAA, z17, VSHUF(t4, 0x4E)); z18 = VMADD(z18, 0xAA, z18, t4);
  h17 = VMADD(h17, 0xAA, h17, VSHUF(t5, 0x4E)); h18 = VMADD(h18, 0xAA, h18, t5);
  z19 = VMADD(z19, 0xAA, z19, VSHUF(t6, 0x4E)); z20 = VMADD(z20, 0xAA, z20, t6);
  h19 = VMADD(h19, 0xAA, h19, VSHUF(t7, 0x4E)); h20 = VMADD(h20, 0xAA, h20, t7);
  z22 = VMADD(z22, 0x55, z22, VSHUF(z14, 0x4E));
  z15 = VADD(z15, VADD(h14, h14));

  z16 = VADD(z16, VADD(h15, h15));
  z17 = VADD(z17, VADD(h16, h16));
  z18 = VADD(z18, VADD(h17, h17));
  z19 = VADD(z19, VADD(h18, h18));
  z20 = VADD(z20, VADD(h19, h19));
  z21 = VADD(z21, VADD(h20, h20));
  z22 = VADD(z22, VADD(h21, h21));  

  // ---------------------------------------------------------------------------
  // *simple* carry propagation
  // some limbs are finally 52-bit not 51-bit 

  r0 = z15; r1 = z16; r2 = z17; r3 = z18; 
  r4 = z19; r5 = z20; r6 = z21; r7 = z22;

  c = VSRA(r0, VBRADIX); r0 = VAND(r0, vbmask); r1 = VADD(r1, c);
  c = VSRA(r1, VBRADIX); r1 = VAND(r1, vbmask); r2 = VADD(r2, c);
  c = VSRA(r2, VBRADIX); r2 = VAND(r2, vbmask); r3 = VADD(r3, c);
  c = VSRA(r3, VBRADIX); r3 = VAND(r3, vbmask); r4 = VADD(r4, c);
  c = VSRA(r4, VBRADIX); r4 = VAND(r4, vbmask); r5 = VADD(r5, c);
  c = VSRA(r5, VBRADIX); r5 = VAND(r5, vbmask); r6 = VADD(r6, c);  
  c = VSRA(r6, VBRADIX); r6 = VAND(r6, vbmask); r7 = VADD(r7, c);  
  c = VSRA(r7, VBRADIX); r7 = VAND(r7, vbmask);
  r0 = VMADD(r0, 0xAA, r0, VSHUF(c, 0x4E));

  r[0] = r0; r[1] = r1; r[2] = r2; r[3] = r3; 
  r[4] = r4; r[5] = r5; r[6] = r6; r[7] = r7;
}

// -----------------------------------------------------------------------------
// 1-way x64 Fp arithmetic (from PQCrypto-SIDH-3.4) 

void mp_mul(const digit_t* a, const digit_t* b, digit_t* c, const unsigned int nwords)
{ // Multiprecision multiply, c = a*b, where lng(a) = lng(b) = nwords.
    
    mul751_asm(a, b, c);
}

void rdc_mont(digit_t* ma, digit_t* mc)
{ // Montgomery reduction exploiting special form of the prime.
  // mc = ma*R^-1 mod p503x2, where R = 2^512.
  // If ma < 2^512*p503, the output mc is in the range [0, 2*p503-1].
  // ma is assumed to be in Montgomery representation.
   
    rdc751_asm(ma, mc);    
}

void fpcorrection(digit_t* a)
{ // Modular correction to reduce field element a in [0, 2*p503-1] to [0, p503-1].
    unsigned int i, borrow = 0;
    digit_t mask;

    for (i = 0; i < NWORDS_FIELD; i++) {
        SUBC(borrow, a[i], ((digit_t*)p751)[i], borrow, a[i]); 
    }
    mask = 0 - (digit_t)borrow;

    borrow = 0;
    for (i = 0; i < NWORDS_FIELD; i++) {
        ADDC(borrow, a[i], ((digit_t*)p751)[i] & mask, borrow, a[i]); 
    }
}

void fpadd(const digit_t* a, const digit_t* b, digit_t* c)
{ // Modular addition, c = a+b mod p503.
  // Inputs: a, b in [0, 2*p503-1] 
  // Output: c in [0, 2*p503-1]               
    
    fpadd751_asm(a, b, c);    
} 

void fpneg(digit_t* a)
{ // Modular negation, a = -a mod p503.
  // Input/output: a in [0, 2*p503-1] 
    unsigned int i, borrow = 0;
    
    for (i = 0; i < NWORDS_FIELD; i++) {
        SUBC(borrow, ((digit_t*)p751x2)[i], a[i], borrow, a[i]); 
    }
}

void fpsub(const digit_t* a, const digit_t* b, digit_t* c)
{ // Modular subtraction, c = a-b mod p503.
  // Inputs: a, b in [0, 2*p503-1] 
  // Output: c in [0, 2*p503-1] 
    
    fpsub751_asm(a, b, c);    
}

void mp_sub_p4(const digit_t* a, const digit_t* b, digit_t* c)
{ // Multiprecision subtraction with correction with 4*p, c = a-b+4p.                   
    
    mp_sub751_p4_asm(a, b, c);    
} 
