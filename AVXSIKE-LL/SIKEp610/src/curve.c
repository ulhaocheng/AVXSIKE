/**
 *******************************************************************************
 * @version 0.1
 * @date 2021-12-16
 * @copyright Copyright © 2021 by University of Luxembourg
 * @author Hao Cheng
 *******************************************************************************
 */

#include "curve.h"
#include <stdio.h>
#include "utils.h"

// -----------------------------------------------------------------------------
// (8x1x1x1)-way curve/isogeny arithmetic (based on (8x1x1)-way Fp2 arithmetic).
// NOTE: (8x1x1x1)-way means each function performs 8 curve/isogeny operations, 
// where each curve/isogeny operation performs 1 Fp2 operation, and each Fp2 
// operation performs 1 Fp operation and each Fp operation uses 1 lane. 

void xDBL_8x1x1x1w(const vpoint_proj_t P, vpoint_proj_t Q, const vf2elm_t A24plus, const vf2elm_t C24)
{
  vf2elm_t t0, t1;

  mp2_sub_p2_8x1x1w(t0, P->X, P->Z);           // t0 = X1-Z1          
  mp2_add_8x1x1w(t1, P->X, P->Z);              // t1 = X1+Z1          
  fp2sqr_mont_8x1x1w(t0, t0);                  // t0 = (X1-Z1)^2      
  fp2sqr_mont_8x1x1w(t1, t1);                  // t1 = (X1+Z1)^2     
  fp2mul_mont_8x1x1w(Q->Z, C24, t0);           // Z2 = C24*(X1-Z1)^2  
  fp2mul_mont_8x1x1w(Q->X, t1, Q->Z);          // X2 = C24*(X1-Z1)^2*(X1+Z1)^2 
  mp2_sub_p2_8x1x1w(t1, t1, t0);               // t1 = (X1+Z1)^2-(X1-Z1)^2  
  fp2mul_mont_8x1x1w(t0, A24plus, t1);         // t0 = A24plus*[(X1+Z1)^2-(X1-Z1)^2]   
  mp2_add_8x1x1w(Q->Z, Q->Z, t0);              // Z2 = A24plus*[(X1+Z1)^2-(X1-Z1)^2] + C24*(X1-Z1)^2
  fp2mul_mont_8x1x1w(Q->Z, Q->Z, t1);          // Z2 = [A24plus*[(X1+Z1)^2-(X1-Z1)^2] + C24*(X1-Z1)^2]*[(X1+Z1)^2-(X1-Z1)^2]
}

void get_4_isog_8x1x1x1w(const vpoint_proj_t P, vf2elm_t A24plus, vf2elm_t C24, vf2elm_t *coeff)
{
  mp2_sub_p2_8x1x1w(coeff[1], P->X, P->Z);               // coeff[1] = X4-Z4
  mp2_add_8x1x1w(coeff[2], P->X, P->Z);                  // coeff[2] = X4+Z4
  fp2sqr_mont_8x1x1w(coeff[0], P->Z);                    // coeff[0] = Z4^2
  mp2_add_8x1x1w(coeff[0], coeff[0], coeff[0]);          // coeff[0] = 2*Z4^2
  fp2sqr_mont_8x1x1w(C24, coeff[0]);                     // C24 = 4*Z4^4
  mp2_add_8x1x1w(coeff[0], coeff[0], coeff[0]);          // coeff[0] = 4*Z4^2
  fp2sqr_mont_8x1x1w(A24plus, P->X);                     // A24plus = X4^2
  mp2_add_8x1x1w(A24plus, A24plus, A24plus);             // A24plus = 2*X4^2
  fp2sqr_mont_8x1x1w(A24plus, A24plus);                  // A24plus = 4*X4^4
}

void eval_4_isog_8x1x1x1w(vpoint_proj_t P, vf2elm_t *coeff)
{
  vf2elm_t t0, t1;

  mp2_add_8x1x1w(t0, P->X, P->Z);                        // t0 = X+Z
  mp2_sub_p2_8x1x1w(t1, P->X, P->Z);                     // t1 = X-Z
  fp2mul_mont_8x1x1w(P->X, t0, coeff[1]);                // X = (X+Z)*coeff[1]
  fp2mul_mont_8x1x1w(P->Z, t1, coeff[2]);                // Z = (X-Z)*coeff[2]
  fp2mul_mont_8x1x1w(t0, t0, t1);                        // t0 = (X+Z)*(X-Z)
  fp2mul_mont_8x1x1w(t0, coeff[0], t0);                  // t0 = coeff[0]*(X+Z)*(X-Z)
  mp2_add_8x1x1w(t1, P->X, P->Z);                        // t1 = (X-Z)*coeff[2] + (X+Z)*coeff[1]
  mp2_sub_p2_8x1x1w(P->Z, P->X, P->Z);                   // Z = (X-Z)*coeff[2] - (X+Z)*coeff[1]
  fp2sqr_mont_8x1x1w(t1, t1);                            // t1 = [(X-Z)*coeff[2] + (X+Z)*coeff[1]]^2
  fp2sqr_mont_8x1x1w(P->Z, P->Z);                        // Z = [(X-Z)*coeff[2] - (X+Z)*coeff[1]]^2
  mp2_add_8x1x1w(P->X, t1, t0);                          // X = coeff[0]*(X+Z)*(X-Z) + [(X-Z)*coeff[2] + (X+Z)*coeff[1]]^2
  mp2_sub_p2_8x1x1w(t0, P->Z, t0);                       // t0 = [(X-Z)*coeff[2] - (X+Z)*coeff[1]]^2 - coeff[0]*(X+Z)*(X-Z) 
  fp2mul_mont_8x1x1w(P->X, P->X, t1);                    // Xfinal
  fp2mul_mont_8x1x1w(P->Z, P->Z, t0);                    // Zfinal
}

void xTPL_8x1x1x1w(const vpoint_proj_t P, vpoint_proj_t Q, const vf2elm_t A24minus, const vf2elm_t A24plus)  
{
  vf2elm_t t0, t1, t2, t3, t4, t5, t6;

  mp2_sub_p2_8x1x1w(t0, P->X, P->Z);                     // t0 = X-Z 
  fp2sqr_mont_8x1x1w(t2, t0);                            // t2 = (X-Z)^2           
  mp2_add_8x1x1w(t1, P->X, P->Z);                        // t1 = X+Z 
  fp2sqr_mont_8x1x1w(t3, t1);                            // t3 = (X+Z)^2
  mp2_add_8x1x1w(t4, P->X, P->X);                        // t4 = 2*X
  mp2_add_8x1x1w(t0, P->Z, P->Z);                        // t0 = 2*Z 
  fp2sqr_mont_8x1x1w(t1, t4);                            // t1 = 4*X^2
  mp2_sub_p2_8x1x1w(t1, t1, t3);                         // t1 = 4*X^2 - (X+Z)^2 
  mp2_sub_p2_8x1x1w(t1, t1, t2);                         // t1 = 4*X^2 - (X+Z)^2 - (X-Z)^2
  fp2mul_mont_8x1x1w(t5, A24plus, t3);                   // t5 = A24plus*(X+Z)^2 
  fp2mul_mont_8x1x1w(t3, t3, t5);                        // t3 = A24plus*(X+Z)^4
  fp2mul_mont_8x1x1w(t6, A24minus, t2);                  // t6 = A24minus*(X-Z)^2
  fp2mul_mont_8x1x1w(t2, t2, t6);                        // t2 = A24minus*(X-Z)^4
  mp2_sub_p2_8x1x1w(t3, t2, t3);                         // t3 = A24minus*(X-Z)^4 - A24plus*(X+Z)^4
  mp2_sub_p2_8x1x1w(t2, t5, t6);                         // t2 = A24plus*(X+Z)^2 - A24minus*(X-Z)^2
  fp2mul_mont_8x1x1w(t1, t1, t2);                        // t1 = [4*X^2 - (X+Z)^2 - (X-Z)^2]*[A24plus*(X+Z)^2 - A24minus*(X-Z)^2]
  fp2add_8x1x1w(t2, t3, t1);                             // t2 = [4*X^2 - (X+Z)^2 - (X-Z)^2]*[A24plus*(X+Z)^2 - A24minus*(X-Z)^2] + A24minus*(X-Z)^4 - A24plus*(X+Z)^4
  fp2sqr_mont_8x1x1w(t2, t2);                            // t2 = t2^2
  fp2mul_mont_8x1x1w(Q->X, t4, t2);                      // X3 = 2*X*t2
  fp2sub_8x1x1w(t1, t3, t1);                             // t1 = A24minus*(X-Z)^4 - A24plus*(X+Z)^4 - [4*X^2 - (X+Z)^2 - (X-Z)^2]*[A24plus*(X+Z)^2 - A24minus*(X-Z)^2]
  fp2sqr_mont_8x1x1w(t1, t1);                            // t1 = t1^2
  fp2mul_mont_8x1x1w(Q->Z, t0, t1);                      // Z3 = 2*Z*t1
}

void get_3_isog_8x1x1x1w(const vpoint_proj_t P, vf2elm_t A24minus, vf2elm_t A24plus, vf2elm_t* coeff)
{
  vf2elm_t t0, t1, t2, t3, t4;

  mp2_sub_p2_8x1x1w(coeff[0], P->X, P->Z);               // coeff0 = X-Z
  fp2sqr_mont_8x1x1w(t0, coeff[0]);                      // t0 = (X-Z)^2
  mp2_add_8x1x1w(coeff[1], P->X, P->Z);                  // coeff1 = X+Z
  fp2sqr_mont_8x1x1w(t1, coeff[1]);                      // t1 = (X+Z)^2
  mp2_add_8x1x1w(t3, P->X, P->X);                        // t3 = 2*X
  fp2sqr_mont_8x1x1w(t3, t3);                            // t3 = 4*X^2 
  fp2sub_8x1x1w(t2, t3, t0);                             // t2 = 4*X^2 - (X-Z)^2 
  fp2sub_8x1x1w(t3, t3, t1);                             // t3 = 4*X^2 - (X+Z)^2
  mp2_add_8x1x1w(t4, t0, t3);                            // t4 = 4*X^2 - (X+Z)^2 + (X-Z)^2 
  mp2_add_8x1x1w(t4, t4, t4);                            // t4 = 2(4*X^2 - (X+Z)^2 + (X-Z)^2) 
  mp2_add_8x1x1w(t4, t1, t4);                            // t4 = 8*X^2 - (X+Z)^2 + 2*(X-Z)^2
  fp2mul_mont_8x1x1w(A24minus, t2, t4);                  // A24minus = [4*X^2 - (X-Z)^2]*[8*X^2 - (X+Z)^2 + 2*(X-Z)^2]
  mp2_add_8x1x1w(t4, t1, t2);                            // t4 = 4*X^2 + (X+Z)^2 - (X-Z)^2
  mp2_add_8x1x1w(t4, t4, t4);                            // t4 = 2(4*X^2 + (X+Z)^2 - (X-Z)^2) 
  mp2_add_8x1x1w(t4, t0, t4);                            // t4 = 8*X^2 + 2*(X+Z)^2 - (X-Z)^2
  fp2mul_mont_8x1x1w(A24plus, t3, t4);                   // A24plus = [4*X^2 - (X+Z)^2]*[8*X^2 + 2*(X+Z)^2 - (X-Z)^2]
}

void eval_3_isog_8x1x1x1w(vpoint_proj_t Q, const vf2elm_t *coeff)
{
  vf2elm_t t0, t1, t2;

  mp2_add_8x1x1w(t0, Q->X, Q->Z);                      // t0 = X+Z
  mp2_sub_p2_8x1x1w(t1, Q->X, Q->Z);                   // t1 = X-Z
  fp2mul_mont_8x1x1w(t0, coeff[0], t0);                // t0 = coeff0*(X+Z)
  fp2mul_mont_8x1x1w(t1, coeff[1], t1);                // t1 = coeff1*(X-Z)
  mp2_add_8x1x1w(t2, t0, t1);                          // t2 = coeff0*(X+Z) + coeff1*(X-Z)
  mp2_sub_p2_8x1x1w(t0, t1, t0);                       // t0 = coeff1*(X-Z) - coeff0*(X+Z)
  fp2sqr_mont_8x1x1w(t2, t2);                          // t2 = [coeff0*(X+Z) + coeff1*(X-Z)]^2
  fp2sqr_mont_8x1x1w(t0, t0);                          // t0 = [coeff1*(X-Z) - coeff0*(X+Z)]^2
  fp2mul_mont_8x1x1w(Q->X, Q->X, t2);                  // X3final = X*[coeff0*(X+Z) + coeff1*(X-Z)]^2        
  fp2mul_mont_8x1x1w(Q->Z, Q->Z, t0);                  // Z3final = Z*[coeff1*(X-Z) - coeff0*(X+Z)]^2
}

void xDBLADD_8x1x1x1w(vpoint_proj_t P, vpoint_proj_t Q, const vf2elm_t XPQ, const vf2elm_t ZPQ, const vf2elm_t A24)
{
  vf2elm_t t0, t1, t2;

  mp2_add_8x1x1w(t0, P->X, P->Z);                        // t0 = XP+ZP
  mp2_sub_p2_8x1x1w(t1, P->X, P->Z);                     // t1 = XP-ZP
  fp2sqr_mont_8x1x1w(P->X, t0);                          // XP = (XP+ZP)^2
  mp2_sub_p2_8x1x1w(t2, Q->X, Q->Z);                     // t2 = XQ-ZQ
  mp2_add_8x1x1w(Q->X, Q->X, Q->Z);                      // XQ = XQ+ZQ
  fp2mul_mont_8x1x1w(t0, t0, t2);                        // t0 = (XP+ZP)*(XQ-ZQ)
  fp2sqr_mont_8x1x1w(P->Z, t1);                          // ZP = (XP-ZP)^2
  fp2mul_mont_8x1x1w(t1, t1, Q->X);                      // t1 = (XP-ZP)*(XQ+ZQ)
  mp2_sub_p2_8x1x1w(t2, P->X, P->Z);                     // t2 = (XP+ZP)^2-(XP-ZP)^2
  fp2mul_mont_8x1x1w(P->X, P->X, P->Z);                  // XP = (XP+ZP)^2*(XP-ZP)^2
  fp2mul_mont_8x1x1w(Q->X, A24, t2);                     // XQ = A24*[(XP+ZP)^2-(XP-ZP)^2]
  mp2_sub_p2_8x1x1w(Q->Z, t0, t1);                       // ZQ = (XP+ZP)*(XQ-ZQ)-(XP-ZP)*(XQ+ZQ)
  mp2_add_8x1x1w(P->Z, Q->X, P->Z);                      // ZP = A24*[(XP+ZP)^2-(XP-ZP)^2]+(XP-ZP)^2
  mp2_add_8x1x1w(Q->X, t0, t1);                          // XQ = (XP+ZP)*(XQ-ZQ)+(XP-ZP)*(XQ+ZQ)
  fp2mul_mont_8x1x1w(P->Z, P->Z, t2);                    // ZP = [A24*[(XP+ZP)^2-(XP-ZP)^2]+(XP-ZP)^2]*[(XP+ZP)^2-(XP-ZP)^2]
  fp2sqr_mont_8x1x1w(Q->Z, Q->Z);                        // ZQ = [(XP+ZP)*(XQ-ZQ)-(XP-ZP)*(XQ+ZQ)]^2
  fp2sqr_mont_8x1x1w(Q->X, Q->X);                        // XQ = [(XP+ZP)*(XQ-ZQ)+(XP-ZP)*(XQ+ZQ)]^2
  fp2mul_mont_8x1x1w(Q->Z, Q->Z, XPQ);                   // ZQ = xPQ*[(XP+ZP)*(XQ-ZQ)-(XP-ZP)*(XQ+ZQ)]^2
  fp2mul_mont_8x1x1w(Q->X, Q->X, ZPQ);                   // XQ = ZPQ*[(XP+ZP)*(XQ-ZQ)+(XP-ZP)*(XQ+ZQ)]^2 
}

// -----------------------------------------------------------------------------
// (1x4x2x1)-way curve/isogeny arithmetic (based on (4x2x1)-way Fp2 arithmetic).
// NOTE: (1x4x2x1)-way means each function performs 1 curve/isogeny operation, 
// where each curve/isogeny operation performs 4 Fp2 operations, and each Fp2 
// operation performs 2 Fp operations and each Fp operation uses 1 lane. 

// INPUT  a = < D1 | D0 | C1 | C0 | B1 | B0 | A1 | A0 >                         all in [0, 2p)
// OUTPUT r = < D1+C1 | D0+C0 | D1-C1 | D0-C0 | B1+A1 | B0+A0 | B1-A1 | B0-A0>  all in [0, 4p)
static void mp2_hadamard_4x2x1w(vfelm_t r, const vfelm_t a)
{
  __m512i a0 = a[0], a1 = a[1], a2  = a[2],  a3  = a[3];
  __m512i a4 = a[4], a5 = a[5], a6  = a[6],  a7  = a[7];
  __m512i a8 = a[8], a9 = a[9], a10 = a[10], a11 = a[11];
  __m512i r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11;
  __m512i t0, t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11;
  const __m512i vp0  = VSET1(vp610x2[0]),  vp1  = VSET1(vp610x2[1]);
  const __m512i vp2  = VSET1(vp610x2[2]),  vp3  = VSET1(vp610x2[3]);
  const __m512i vp4  = VSET1(vp610x2[4]),  vp5  = VSET1(vp610x2[5]);
  const __m512i vp6  = VSET1(vp610x2[6]),  vp7  = VSET1(vp610x2[7]);
  const __m512i vp8  = VSET1(vp610x2[8]),  vp9  = VSET1(vp610x2[9]);
  const __m512i vp10 = VSET1(vp610x2[10]), vp11 = VSET1(vp610x2[11]);
  const __m512i vbmask = VSET1(VBMASK); 

  // t = C1 | C0 | D1 | D0 | A1 | A0 | B1 | B0
  t0 = VPERM(a0, 0x4E); t1 = VPERM(a1, 0x4E);
  t2 = VPERM(a2, 0x4E); t3 = VPERM(a3, 0x4E);
  t4 = VPERM(a4, 0x4E); t5 = VPERM(a5, 0x4E);
  t6 = VPERM(a6, 0x4E); t7 = VPERM(a7, 0x4E);
  t8 = VPERM(a8, 0x4E); t9 = VPERM(a9, 0x4E);
  t10 = VPERM(a10, 0x4E); t11 = VPERM(a11, 0x4E);
           
  // r = D1 | D0 | 2p | 2p | B1 | B0 | 2p | 2p
  r0 = VMBLEND(0x33, a0, vp0); r1 = VMBLEND(0x33, a1, vp1);
  r2 = VMBLEND(0x33, a2, vp2); r3 = VMBLEND(0x33, a3, vp3);
  r4 = VMBLEND(0x33, a4, vp4); r5 = VMBLEND(0x33, a5, vp5);
  r6 = VMBLEND(0x33, a6, vp6); r7 = VMBLEND(0x33, a7, vp7);
  r8 = VMBLEND(0x33, a8, vp8); r9 = VMBLEND(0x33, a9, vp9);
  r10 = VMBLEND(0x33, a10, vp10); r11 = VMBLEND(0x33, a11, vp11);

  // r =  D1+C1 | D0+C0 | 2p+D1 | 2p+D0 | B1+A1 | B0+A0 | 2p+B1 | 2p+B0
  r0 = VADD(r0, t0); r1 = VADD(r1, t1);
  r2 = VADD(r2, t2); r3 = VADD(r3, t3);
  r4 = VADD(r4, t4); r5 = VADD(r5, t5);
  r6 = VADD(r6, t6); r7 = VADD(r7, t7);
  r8 = VADD(r8, t8); r9 = VADD(r9, t9);
  r10 = VADD(r10, t10); r11 = VADD(r11, t11);

  // r = D1+C1 | D0+C0 | 2p+D1-C1 | 2p+D0-C0 | B1+A1 | B0+A0 | 2p+B1-A1 | 2p+B0-A0
  r0 = VMSUB(r0, 0x33, r0, a0); r1 = VMSUB(r1, 0x33, r1, a1);
  r2 = VMSUB(r2, 0x33, r2, a2); r3 = VMSUB(r3, 0x33, r3, a3);
  r4 = VMSUB(r4, 0x33, r4, a4); r5 = VMSUB(r5, 0x33, r5, a5);
  r6 = VMSUB(r6, 0x33, r6, a6); r7 = VMSUB(r7, 0x33, r7, a7);
  r8 = VMSUB(r8, 0x33, r8, a8); r9 = VMSUB(r9, 0x33, r9, a9);
  r10 = VMSUB(r10, 0x33, r10, a10); r11 = VMSUB(r11, 0x33, r11, a11);

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

// INPUT  a = < D1 | D0 | C1 | C0 | B1 | B0 | A1 | A0 >
// OUTPUT r = < A1 | A0 | B1 | B0 | B1 | B0 | A1 | A0 >
static void fp2align_4x2x1w(vfelm_t r, vfelm_t a)
{
  __m512i a0 = a[0], a1 = a[1], a2  = a[2],  a3  = a[3];
  __m512i a4 = a[4], a5 = a[5], a6  = a[6],  a7  = a[7];
  __m512i a8 = a[8], a9 = a[9], a10 = a[10], a11 = a[11];
  __m512i r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11;
  __m512i mask = VSET(1, 0, 3, 2, 3, 2, 1, 0); 

  r0 = VPERMV(mask, a0); r1 = VPERMV(mask, a1);
  r2 = VPERMV(mask, a2); r3 = VPERMV(mask, a3);
  r4 = VPERMV(mask, a4); r5 = VPERMV(mask, a5);
  r6 = VPERMV(mask, a6); r7 = VPERMV(mask, a7);
  r8 = VPERMV(mask, a8); r9 = VPERMV(mask, a9);
  r10 = VPERMV(mask, a10); r11 = VPERMV(mask, a11);

  r[0] = r0; r[1] = r1; r[2]  = r2;  r[3]  = r3; 
  r[4] = r4; r[5] = r5; r[6]  = r6;  r[7]  = r7; 
  r[8] = r8; r[9] = r9; r[10] = r10; r[11] = r11;
}

// INPUT  a = < D1 | D0 | C1 | C0 | B1 | B0 | A1 | A0 >
//        b = < H1 | H0 | G1 | G0 | F1 | F0 | E1 | E0 >
//        c = < L1 | L0 | K1 | K0 | J1 | J0 | I1 | I0 >
// OUTPUT r = < D1 | D0 | C1 | C0 | F1 | F0 | A1 | A0 >
//        s = < L1 | L0 | K1 | K0 | E1 | E0 | I1 | I0 >
static void fp2premix_4x2x1w(vfelm_t r, vfelm_t s, const vfelm_t a, const vfelm_t b, const vfelm_t c)
{
  __m512i a0 = a[0], a1 = a[1], a2  = a[2],  a3  = a[3];
  __m512i a4 = a[4], a5 = a[5], a6  = a[6],  a7  = a[7];
  __m512i a8 = a[8], a9 = a[9], a10 = a[10], a11 = a[11];
  __m512i b0 = b[0], b1 = b[1], b2  = b[2],  b3  = b[3]; 
  __m512i b4 = b[4], b5 = b[5], b6  = b[6],  b7  = b[7];
  __m512i b8 = b[8], b9 = b[9], b10 = b[10], b11 = b[11];
  __m512i c0 = c[0], c1 = c[1], c2  = c[2],  c3  = c[3]; 
  __m512i c4 = c[4], c5 = c[5], c6  = c[6],  c7  = c[7];
  __m512i c8 = c[8], c9 = c[9], c10 = c[10], c11 = c[11];
  __m512i r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11;
  __m512i s0, s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, s11;

  // r = < D1 | D0 | C1 | C0 | F1 | F0 | A1 | A0 >
  r0 = VMBLEND(0x0C, a0, b0); r1 = VMBLEND(0x0C, a1, b1);
  r2 = VMBLEND(0x0C, a2, b2); r3 = VMBLEND(0x0C, a3, b3);
  r4 = VMBLEND(0x0C, a4, b4); r5 = VMBLEND(0x0C, a5, b5);
  r6 = VMBLEND(0x0C, a6, b6); r7 = VMBLEND(0x0C, a7, b7);
  r8 = VMBLEND(0x0C, a8, b8); r9 = VMBLEND(0x0C, a9, b9);
  r10 = VMBLEND(0x0C, a10, b10); r11 = VMBLEND(0x0C, a11, b11);

  // b = < G1 | G0 | H1 | H0 | E1 | E0 | F1 | F0 >
  b0 = VPERM(b0, 0x4E); b1 = VPERM(b1, 0x4E);
  b2 = VPERM(b2, 0x4E); b3 = VPERM(b3, 0x4E);
  b4 = VPERM(b4, 0x4E); b5 = VPERM(b5, 0x4E);
  b6 = VPERM(b6, 0x4E); b7 = VPERM(b7, 0x4E);
  b8 = VPERM(b8, 0x4E); b9 = VPERM(b9, 0x4E);
  b10 = VPERM(b10, 0x4E); b11 = VPERM(b11, 0x4E);

  // s = < L1 | L0 | K1 | K0 | E1 | E0 | I1 | I0 >
  s0 = VMBLEND(0x0C, c0, b0); s1 = VMBLEND(0x0C, c1, b1);
  s2 = VMBLEND(0x0C, c2, b2); s3 = VMBLEND(0x0C, c3, b3);
  s4 = VMBLEND(0x0C, c4, b4); s5 = VMBLEND(0x0C, c5, b5);
  s6 = VMBLEND(0x0C, c6, b6); s7 = VMBLEND(0x0C, c7, b7);
  s8 = VMBLEND(0x0C, c8, b8); s9 = VMBLEND(0x0C, c9, b9);
  s10 = VMBLEND(0x0C, c10, b10); s11 = VMBLEND(0x0C, c11, b11);

  r[0] = r0; r[1] = r1; r[2]  = r2;  r[3]  = r3; 
  r[4] = r4; r[5] = r5; r[6]  = r6;  r[7]  = r7; 
  r[8] = r8; r[9] = r9; r[10] = r10; r[11] = r11;

  s[0] = s0; s[1] = s1; s[2]  = s2;  s[3]  = s3; 
  s[4] = s4; s[5] = s5; s[6]  = s6;  s[7]  = s7; 
  s[8] = s8; s[9] = s9; s[10] = s10; s[11] = s11;
}

// INPUT  a = < D1 | D0 | C1 | C0 | B1 | B0 | A1 | A0 >               all in [0, 2p)
//        b = < H1 | H0 | G1 | G0 | F1 | F0 | E1 | E0 >               all in [0, 2p)
// OUTPUT r = < D1 | D0 | C1 | C0 | F1-E1 | F0-E0 | 2B1+A1 | 2B0+A0 > all in [0, 2p)
static void fp2postmix_4x2x1w(vfelm_t r, const vfelm_t a, const vfelm_t b)
{
  __m512i a0 = a[0], a1 = a[1], a2  = a[2],  a3  = a[3];
  __m512i a4 = a[4], a5 = a[5], a6  = a[6],  a7  = a[7];
  __m512i a8 = a[8], a9 = a[9], a10 = a[10], a11 = a[11];
  __m512i b0 = b[0], b1 = b[1], b2  = b[2],  b3  = b[3]; 
  __m512i b4 = b[4], b5 = b[5], b6  = b[6],  b7  = b[7];
  __m512i b8 = b[8], b9 = b[9], b10 = b[10], b11 = b[11];
  __m512i r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, smask;
  __m512i t0, t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11;
  __m512i s0, s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, s11;
  const __m512i vp0  = VSET1(vp610x2[0]),  vp1  = VSET1(vp610x2[1]);
  const __m512i vp2  = VSET1(vp610x2[2]),  vp3  = VSET1(vp610x2[3]);
  const __m512i vp4  = VSET1(vp610x2[4]),  vp5  = VSET1(vp610x2[5]);
  const __m512i vp6  = VSET1(vp610x2[6]),  vp7  = VSET1(vp610x2[7]);
  const __m512i vp8  = VSET1(vp610x2[8]),  vp9  = VSET1(vp610x2[9]);
  const __m512i vp10 = VSET1(vp610x2[10]), vp11 = VSET1(vp610x2[11]);
  const __m512i vbmask = VSET1(VBMASK); 

  // t = < 2D1 | 2D0 | 2C1 | 2C0 | 2B1 | 2B0 | 2A1 | 2A0 >
  t0 = VADD(a0, a0); t1 = VADD(a1, a1); 
  t2 = VADD(a2, a2); t3 = VADD(a3, a3); 
  t4 = VADD(a4, a4); t5 = VADD(a5, a5); 
  t6 = VADD(a6, a6); t7 = VADD(a7, a7); 
  t8 = VADD(a8, a8); t9 = VADD(a9, a9); 
  t10 = VADD(a10, a10); t11 = VADD(a11, a11); 

  // t = < 2C1 | 2C0 | 2D1 | 2D0 | 2A1 | 2A0 | 2B1 | 2B0 >
  t0 = VPERM(t0, 0x4E); t1 = VPERM(t1, 0x4E);
  t2 = VPERM(t2, 0x4E); t3 = VPERM(t3, 0x4E);
  t4 = VPERM(t4, 0x4E); t5 = VPERM(t5, 0x4E);
  t6 = VPERM(t6, 0x4E); t7 = VPERM(t7, 0x4E);
  t8 = VPERM(t8, 0x4E); t9 = VPERM(t9, 0x4E);
  t10 = VPERM(t10, 0x4E); t11 = VPERM(t11, 0x4E);

  // s = < G1 | G0 | H1 | H0 | E1 | E0 | F1 | F0 >
  s0 = VPERM(b0, 0x4E); s1 = VPERM(b1, 0x4E);
  s2 = VPERM(b2, 0x4E); s3 = VPERM(b3, 0x4E);
  s4 = VPERM(b4, 0x4E); s5 = VPERM(b5, 0x4E);
  s6 = VPERM(b6, 0x4E); s7 = VPERM(b7, 0x4E);
  s8 = VPERM(b8, 0x4E); s9 = VPERM(b9, 0x4E);
  s10 = VPERM(b10, 0x4E); s11 = VPERM(b11, 0x4E);

  // t = t + a = < 2C1+D1 | 2C0+D0 | 2D1+C1 | 2D0+C0 | 2A1+B1 | 2A0+B0 | 2B1+A1 | 2B0+A0 >
  t0 = VADD(t0, a0); t1 = VADD(t1, a1);
  t2 = VADD(t2, a2); t3 = VADD(t3, a3);
  t4 = VADD(t4, a4); t5 = VADD(t5, a5);
  t6 = VADD(t6, a6); t7 = VADD(t7, a7);
  t8 = VADD(t8, a8); t9 = VADD(t9, a9);
  t10 = VADD(t10, a10); t11 = VADD(t11, a11);

  // s = b + 2p - s = < H1+2p-G1 | H0+2p-G0 | G1+2p-H1 | G0+2p-H0 | F1+2p-E1 | F0+2p-E0 | E1+2p-F1 | E0+2p-F0 >
  s0 = VSUB(VADD(b0, vp0), s0); s1 = VSUB(VADD(b1, vp1), s1); 
  s2 = VSUB(VADD(b2, vp2), s2); s3 = VSUB(VADD(b3, vp3), s3); 
  s4 = VSUB(VADD(b4, vp4), s4); s5 = VSUB(VADD(b5, vp5), s5); 
  s6 = VSUB(VADD(b6, vp6), s6); s7 = VSUB(VADD(b7, vp7), s7); 
  s8 = VSUB(VADD(b8, vp8), s8); s9 = VSUB(VADD(b9, vp9), s9); 
  s10 = VSUB(VADD(b10, vp10), s10); s11 = VSUB(VADD(b11, vp11), s11); 

  // r = < D1 | D0 | C1 | C0 | F1-E1+2p | F0-E0+2p | A1 | A0 > 
  r0 = VMBLEND(0x0C, a0, s0); r1 = VMBLEND(0x0C, a1, s1); 
  r2 = VMBLEND(0x0C, a2, s2); r3 = VMBLEND(0x0C, a3, s3); 
  r4 = VMBLEND(0x0C, a4, s4); r5 = VMBLEND(0x0C, a5, s5); 
  r6 = VMBLEND(0x0C, a6, s6); r7 = VMBLEND(0x0C, a7, s7); 
  r8 = VMBLEND(0x0C, a8, s8); r9 = VMBLEND(0x0C, a9, s9);
  r10 = VMBLEND(0x0C, a10, s10); r11 = VMBLEND(0x0C, a11, s11); 

  // r =  < D1 | D0 | C1 | C0 | F1-E1+2p | F0-E0+2p | 2B1+A1 | 2B0+A0 >
  //        2p | 2p | 2p | 2p | 4p       | 4p       | 6p     | 6p 
  r0 = VMBLEND(0x03, r0, t0); r1 = VMBLEND(0x03, r1, t1); 
  r2 = VMBLEND(0x03, r2, t2); r3 = VMBLEND(0x03, r3, t3); 
  r4 = VMBLEND(0x03, r4, t4); r5 = VMBLEND(0x03, r5, t5); 
  r6 = VMBLEND(0x03, r6, t6); r7 = VMBLEND(0x03, r7, t7); 
  r8 = VMBLEND(0x03, r8, t8); r9 = VMBLEND(0x03, r9, t9);
  r10 = VMBLEND(0x03, r10, t10); r11 = VMBLEND(0x03, r11, t11); 

  // TODO: try to optimze the following part (3 carry propagations)

  // r = r - 2p
  r0  = VSUB(r0,  vp0);  r1  = VSUB(r1,  vp1);  r2  = VSUB(r2,  vp2);
  r3  = VSUB(r3,  vp3);  r4  = VSUB(r4,  vp4);  r5  = VSUB(r5,  vp5);
  r6  = VSUB(r6,  vp6);  r7  = VSUB(r7,  vp7);  r8  = VSUB(r8,  vp8);
  r9  = VSUB(r9,  vp9);  r10 = VSUB(r10, vp10); r11 = VSUB(r11, vp11);

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
  // r =  < D1 | D0 | C1 | C0 | F1-E1 | F0-E0 | 2B1+A1 | 2B0+A0 >
  //        2p | 2p | 2p | 2p | 2p    | 2p    | 4p     | 4p 

  // r = r - 2p
  r0  = VSUB(r0,  vp0);  r1  = VSUB(r1,  vp1);  r2  = VSUB(r2,  vp2);
  r3  = VSUB(r3,  vp3);  r4  = VSUB(r4,  vp4);  r5  = VSUB(r5,  vp5);
  r6  = VSUB(r6,  vp6);  r7  = VSUB(r7,  vp7);  r8  = VSUB(r8,  vp8);
  r9  = VSUB(r9,  vp9);  r10 = VSUB(r10, vp10); r11 = VSUB(r11, vp11);

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
  // r =  < D1 | D0 | C1 | C0 | F1-E1 | F0-E0 | 2B1+A1 | 2B0+A0 >
  //        2p | 2p | 2p | 2p | 2p    | 2p    | 2p     | 2p

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

// Here we use the algorithm presented in [HEY20] and try to keep our implementation 
// as same as [HEY20] implementation. 
// INPUT 
// x3z3x2z2: <  XQ1 |  XQ0 |  ZQ1 |  ZQ0 | XP1 | XP0 | ZP1 | ZP0 >
// z1x1_A  : < ZPQ1 | ZPQ0 | XPQ1 | XPQ0 |  0  |  0  |  A1 | A0  >
// OUTPUT   P = [2]P, Q = P+Q
// x3z3x2z2: <  XQ1 |  XQ0 |  ZQ1 |  ZQ0 | XP1 | XP0 | ZP1 | ZP0 >
void xDBLADD_1x4x2x1w(vfelm_t x3z3x2z2, vfelm_t z1x1_A)
{
  vfelm_t z2x2x2z2, x5z5x4z4, z1x1x2A, x5z5z2z4, x7z7x6z6;

  mp2_hadamard_4x2x1w(x3z3x2z2, x3z3x2z2);            // x3z3x2z2 = x3+z3 | x3-z3 | x2+z2 | x2-z2
  fp2align_4x2x1w(z2x2x2z2, x3z3x2z2);                // z2x2x2z2 = x2-z2 | x2+z2 | x2+z2 | x2-z2
  fp2mul_mont_4x2x1w(x3z3x2z2, x3z3x2z2, z2x2x2z2);   // x3z3x2z2 = (x3+z3)(x2-z2) | (x3-z3)(x2+z2) | (x2+z2)^2 | (x2-z2)^2
  mp2_hadamard_4x2x1w(x3z3x2z2, x3z3x2z2);            // x3z3x2z2 = 2(x2x3-z2z3) | 2(x2z3-x3z2) | 2(x2^2+z2^2) | 4x2z2
  fp2sqr_mont_4x2x1w(x5z5x4z4, x3z3x2z2);             // x5z5x4z4 = 4(x2x3-z2z3)^2 | 4(x2z3-x3z2)^2 | 4(x2^2+z2^2)^2 | 16 x2^2 z2^2
  fp2premix_4x2x1w(z1x1x2A, x5z5z2z4, z1x1_A, x3z3x2z2, x5z5x4z4); // z1x1x2A = z1 | x1 | 2(x2^2+z2^2) | A, x5z5z2z4 = 4(x2x3-z2z3)^2 | 4(x2z3-x3z2)^2 | 4x2z2 | 16 x2^2 z2^2
  fp2mul_mont_4x2x1w(x7z7x6z6, z1x1x2A, x5z5z2z4);    // x7z7x6z6 = z1(x2x3-z2z3)^2 | x1(x2z3-x3z2)^2 | 2x2z2(x2^2+z2^2) | 4A x2^2 z2^2
  fp2postmix_4x2x1w(x3z3x2z2, x7z7x6z6, x5z5x4z4);    // x3z3x2z2 = z1(x2x3-z2z3)^2 | x1(x2z3-x3z2)^2 | 4(x2^2-z2^2)^2 | 4x2z2(x2^2+z2^2+A x2^2 z2^2)
}

// -----------------------------------------------------------------------------
// (2x4x1x1)-way curve/isogeny arithmetic (based on (8x1x1)-way Fp2 arithmetic).
// NOTE: (2x4x1x1)-way means each function performs 2 curve/isogeny operation, 
// where each curve/isogeny operation performs 4 Fp2 operations, and each Fp2 
// operation performs 1 Fp operations and each Fp operation uses 1 lane. 

// INPUT  a = < H | G | F | E | D | C | B | A >                              all in [0, 2p)
// OUTPUT r = < H+G | H-G+2p | F+E | F-E+2p | D+C | D-C+2p | B+A | B-A+2p >  all in [0, 4p)
static void mp_hadamard_8x1w(vfelm_t r, const vfelm_t a)
{
  __m512i a0 = a[0], a1 = a[1], a2  = a[2],  a3  = a[3];
  __m512i a4 = a[4], a5 = a[5], a6  = a[6],  a7  = a[7];
  __m512i a8 = a[8], a9 = a[9], a10 = a[10], a11 = a[11];
  __m512i r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11;
  __m512i t0, t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11;
  const __m512i vp0  = VSET1(vp610x2[0]),  vp1  = VSET1(vp610x2[1]);
  const __m512i vp2  = VSET1(vp610x2[2]),  vp3  = VSET1(vp610x2[3]);
  const __m512i vp4  = VSET1(vp610x2[4]),  vp5  = VSET1(vp610x2[5]);
  const __m512i vp6  = VSET1(vp610x2[6]),  vp7  = VSET1(vp610x2[7]);
  const __m512i vp8  = VSET1(vp610x2[8]),  vp9  = VSET1(vp610x2[9]);
  const __m512i vp10 = VSET1(vp610x2[10]), vp11 = VSET1(vp610x2[11]);
  const __m512i vbmask = VSET1(VBMASK); 

  // t = < Z | X | Z | X | Z | X | Z | X >
  t0 = VSHUF(a0, 0x4E); t1 = VSHUF(a1, 0x4E);
  t2 = VSHUF(a2, 0x4E); t3 = VSHUF(a3, 0x4E);
  t4 = VSHUF(a4, 0x4E); t5 = VSHUF(a5, 0x4E);
  t6 = VSHUF(a6, 0x4E); t7 = VSHUF(a7, 0x4E);
  t8 = VSHUF(a8, 0x4E); t9 = VSHUF(a9, 0x4E);
  t10 = VSHUF(a10, 0x4E); t11 = VSHUF(a11, 0x4E);

  // r = < X | 2p | X | 2p | X | 2p | X | 2p >
  r0 = VMBLEND(0x55, a0, vp0); r1 = VMBLEND(0x55, a1, vp1);
  r2 = VMBLEND(0x55, a2, vp2); r3 = VMBLEND(0x55, a3, vp3);
  r4 = VMBLEND(0x55, a4, vp4); r5 = VMBLEND(0x55, a5, vp5);
  r6 = VMBLEND(0x55, a6, vp6); r7 = VMBLEND(0x55, a7, vp7);
  r8 = VMBLEND(0x55, a8, vp8); r9 = VMBLEND(0x55, a9, vp9);
  r10 = VMBLEND(0x55, a10, vp10); r11 = VMBLEND(0x55, a11, vp11);

  // r = < X+Z | X+2p | X+Z | X+2p | X+Z | X+2p | X+Z | X+2p >
  r0 = VADD(r0, t0); r1 = VADD(r1, t1);
  r2 = VADD(r2, t2); r3 = VADD(r3, t3);
  r4 = VADD(r4, t4); r5 = VADD(r5, t5);
  r6 = VADD(r6, t6); r7 = VADD(r7, t7);
  r8 = VADD(r8, t8); r9 = VADD(r9, t9);
  r10 = VADD(r10, t10); r11 = VADD(r11, t11);

  // r = < X+Z | X+2p-Z | X+Z | X+2p-Z | X+Z | X+2p-Z | X+Z | X+2p-Z >
  r0 = VMSUB(r0, 0x55, r0, a0); r1 = VMSUB(r1, 0x55, r1, a1);
  r2 = VMSUB(r2, 0x55, r2, a2); r3 = VMSUB(r3, 0x55, r3, a3);
  r4 = VMSUB(r4, 0x55, r4, a4); r5 = VMSUB(r5, 0x55, r5, a5);
  r6 = VMSUB(r6, 0x55, r6, a6); r7 = VMSUB(r7, 0x55, r7, a7);
  r8 = VMSUB(r8, 0x55, r8, a8); r9 = VMSUB(r9, 0x55, r9, a9);
  r10 = VMSUB(r10, 0x55, r10, a10); r11 = VMSUB(r11, 0x55, r11, a11);

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

static void mp2_hadamard_8x1x1w(vf2elm_t r, const vf2elm_t a)
{
  mp_hadamard_8x1w(r[0], a[0]);
  mp_hadamard_8x1w(r[1], a[1]);
}

// INPUT  a = < H | G | F | E | D | C | B | A > 
// OUTPUT r = < D | E | E | D | A | B | B | A >
static void fpalign2_8x1w(vfelm_t r, const vfelm_t a)
{
  __m512i a0 = a[0], a1 = a[1], a2  = a[2],  a3  = a[3];
  __m512i a4 = a[4], a5 = a[5], a6  = a[6],  a7  = a[7];
  __m512i a8 = a[8], a9 = a[9], a10 = a[10], a11 = a[11];
  __m512i r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11;

  r0 = VPERM(a0, 0x14); r1 = VPERM(a1, 0x14);
  r2 = VPERM(a2, 0x14); r3 = VPERM(a3, 0x14);
  r4 = VPERM(a4, 0x14); r5 = VPERM(a5, 0x14);
  r6 = VPERM(a6, 0x14); r7 = VPERM(a7, 0x14);
  r8 = VPERM(a8, 0x14); r9 = VPERM(a9, 0x14);
  r10 = VPERM(a10, 0x14); r11 = VPERM(a11, 0x14);

  r[0] = r0; r[1] = r1; r[2]  = r2;  r[3]  = r3; 
  r[4] = r4; r[5] = r5; r[6]  = r6;  r[7]  = r7; 
  r[8] = r8; r[9] = r9; r[10] = r10; r[11] = r11;
}

static void fp2align2_8x1x1w(vf2elm_t r, const vf2elm_t a)
{
  fpalign2_8x1w(r[0], a[0]);
  fpalign2_8x1w(r[1], a[1]);
}

// INPUT  a = < D | C | B | A >
//        b = < H | G | F | E >
//        c = < L | K | J | I >
// OUTPUT r = < D | C | F | A >
//        s = < L | K | E | I >
static void fppremix_8x1w(vfelm_t r, vfelm_t s, const vfelm_t a, const vfelm_t b, const vfelm_t c)
{
  __m512i a0 = a[0], a1 = a[1], a2  = a[2],  a3  = a[3];
  __m512i a4 = a[4], a5 = a[5], a6  = a[6],  a7  = a[7];
  __m512i a8 = a[8], a9 = a[9], a10 = a[10], a11 = a[11];
  __m512i b0 = b[0], b1 = b[1], b2  = b[2],  b3  = b[3]; 
  __m512i b4 = b[4], b5 = b[5], b6  = b[6],  b7  = b[7];
  __m512i b8 = b[8], b9 = b[9], b10 = b[10], b11 = b[11];
  __m512i c0 = c[0], c1 = c[1], c2  = c[2],  c3  = c[3];
  __m512i c4 = c[4], c5 = c[5], c6  = c[6],  c7  = c[7];
  __m512i c8 = c[8], c9 = c[9], c10 = c[10], c11 = c[11];
  __m512i r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11;
  __m512i s0, s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, s11;

  // r = < D | C | F | A >
  r0 = VMBLEND(0x22, a0, b0); r1 = VMBLEND(0x22, a1, b1);
  r2 = VMBLEND(0x22, a2, b2); r3 = VMBLEND(0x22, a3, b3);
  r4 = VMBLEND(0x22, a4, b4); r5 = VMBLEND(0x22, a5, b5);
  r6 = VMBLEND(0x22, a6, b6); r7 = VMBLEND(0x22, a7, b7);
  r8 = VMBLEND(0x22, a8, b8); r9 = VMBLEND(0x22, a9, b9);
  r10 = VMBLEND(0x22, a10, b10); r11 = VMBLEND(0x22, a11, b11);

  // b = < G | H | E | F >
  b0 = VSHUF(b0, 0x4E); b1 = VSHUF(b1, 0x4E);
  b2 = VSHUF(b2, 0x4E); b3 = VSHUF(b3, 0x4E);
  b4 = VSHUF(b4, 0x4E); b5 = VSHUF(b5, 0x4E);
  b6 = VSHUF(b6, 0x4E); b7 = VSHUF(b7, 0x4E);
  b8 = VSHUF(b8, 0x4E); b9 = VSHUF(b9, 0x4E);
  b10 = VSHUF(b10, 0x4E); b11 = VSHUF(b11, 0x4E);

  // s = < L | K | E | I >
  s0 = VMBLEND(0x22, c0, b0); s1 = VMBLEND(0x22, c1, b1);
  s2 = VMBLEND(0x22, c2, b2); s3 = VMBLEND(0x22, c3, b3);
  s4 = VMBLEND(0x22, c4, b4); s5 = VMBLEND(0x22, c5, b5);
  s6 = VMBLEND(0x22, c6, b6); s7 = VMBLEND(0x22, c7, b7);
  s8 = VMBLEND(0x22, c8, b8); s9 = VMBLEND(0x22, c9, b9);
  s10 = VMBLEND(0x22, c10, b10); s11 = VMBLEND(0x22, c11, b11);

  r[0] = r0; r[1] = r1; r[2]  = r2;  r[3]  = r3; 
  r[4] = r4; r[5] = r5; r[6]  = r6;  r[7]  = r7; 
  r[8] = r8; r[9] = r9; r[10] = r10; r[11] = r11;

  s[0] = s0; s[1] = s1; s[2]  = s2;  s[3]  = s3; 
  s[4] = s4; s[5] = s5; s[6]  = s6;  s[7]  = s7; 
  s[8] = s8; s[9] = s9; s[10] = s10; s[11] = s11;
}

static void fp2premix_8x1x1w(vf2elm_t r, vf2elm_t s, const vf2elm_t a, const vf2elm_t b, const vf2elm_t c)
{
  fppremix_8x1w(r[0], s[0], a[0], b[0], c[0]);
  fppremix_8x1w(r[1], s[1], a[1], b[1], c[1]);  
}

// INPUT  a = < D | C | B | A >               all in [0, 2p)
//        b = < H | G | F | E >               all in [0, 2p)
// OUTPUT r = < D | C | F-E | 2B+A >          all in [0, 2p)
static void fppostmix_8x1w(vfelm_t r, const vfelm_t a, const vfelm_t b)
{
  __m512i a0 = a[0], a1 = a[1], a2  = a[2],  a3  = a[3];
  __m512i a4 = a[4], a5 = a[5], a6  = a[6],  a7  = a[7];
  __m512i a8 = a[8], a9 = a[9], a10 = a[10], a11 = a[11];
  __m512i b0 = b[0], b1 = b[1], b2  = b[2],  b3  = b[3]; 
  __m512i b4 = b[4], b5 = b[5], b6  = b[6],  b7  = b[7];
  __m512i b8 = b[8], b9 = b[9], b10 = b[10], b11 = b[11];
  __m512i r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, smask;
  __m512i t0, t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11;
  __m512i s0, s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, s11;
  const __m512i vp0  = VSET1(vp610x2[0]),  vp1  = VSET1(vp610x2[1]);
  const __m512i vp2  = VSET1(vp610x2[2]),  vp3  = VSET1(vp610x2[3]);
  const __m512i vp4  = VSET1(vp610x2[4]),  vp5  = VSET1(vp610x2[5]);
  const __m512i vp6  = VSET1(vp610x2[6]),  vp7  = VSET1(vp610x2[7]);
  const __m512i vp8  = VSET1(vp610x2[8]),  vp9  = VSET1(vp610x2[9]);
  const __m512i vp10 = VSET1(vp610x2[10]), vp11 = VSET1(vp610x2[11]);
  const __m512i vbmask = VSET1(VBMASK); 

  // t = < 2D | 2C | 2B | 2A >
  t0 = VADD(a0, a0); t1 = VADD(a1, a1); 
  t2 = VADD(a2, a2); t3 = VADD(a3, a3); 
  t4 = VADD(a4, a4); t5 = VADD(a5, a5); 
  t6 = VADD(a6, a6); t7 = VADD(a7, a7); 
  t8 = VADD(a8, a8); t9 = VADD(a9, a9); 
  t10 = VADD(a10, a10); t11 = VADD(a11, a11); 

  // t = < 2C | 2D | 2A | 2B >
  t0 = VSHUF(t0, 0x4E); t1 = VSHUF(t1, 0x4E);
  t2 = VSHUF(t2, 0x4E); t3 = VSHUF(t3, 0x4E);
  t4 = VSHUF(t4, 0x4E); t5 = VSHUF(t5, 0x4E);
  t6 = VSHUF(t6, 0x4E); t7 = VSHUF(t7, 0x4E);
  t8 = VSHUF(t8, 0x4E); t9 = VSHUF(t9, 0x4E);
  t10 = VSHUF(t10, 0x4E); t11 = VSHUF(t11, 0x4E);

  // s = < G | H | E | F >
  s0 = VSHUF(b0, 0x4E); s1 = VSHUF(b1, 0x4E);
  s2 = VSHUF(b2, 0x4E); s3 = VSHUF(b3, 0x4E);
  s4 = VSHUF(b4, 0x4E); s5 = VSHUF(b5, 0x4E);
  s6 = VSHUF(b6, 0x4E); s7 = VSHUF(b7, 0x4E);
  s8 = VSHUF(b8, 0x4E); s9 = VSHUF(b9, 0x4E);
  s10 = VSHUF(b10, 0x4E); s11 = VSHUF(b11, 0x4E);

  // t = t + a = < 2C+D | 2D+C | 2A+B | 2B+A >
  t0 = VADD(t0, a0); t1 = VADD(t1, a1);
  t2 = VADD(t2, a2); t3 = VADD(t3, a3);
  t4 = VADD(t4, a4); t5 = VADD(t5, a5);
  t6 = VADD(t6, a6); t7 = VADD(t7, a7);
  t8 = VADD(t8, a8); t9 = VADD(t9, a9);
  t10 = VADD(t10, a10); t11 = VADD(t11, a11);

  // s = b + 2p - s = < H+2p-G | G+2p-H | F+2p-E | E+2p-F >
  s0 = VSUB(VADD(b0, vp0), s0); s1 = VSUB(VADD(b1, vp1), s1); 
  s2 = VSUB(VADD(b2, vp2), s2); s3 = VSUB(VADD(b3, vp3), s3); 
  s4 = VSUB(VADD(b4, vp4), s4); s5 = VSUB(VADD(b5, vp5), s5); 
  s6 = VSUB(VADD(b6, vp6), s6); s7 = VSUB(VADD(b7, vp7), s7); 
  s8 = VSUB(VADD(b8, vp8), s8); s9 = VSUB(VADD(b9, vp9), s9); 
  s10 = VSUB(VADD(b10, vp10), s10); s11 = VSUB(VADD(b11, vp11), s11); 

  // r = < D | C | F-E+2p | A > 
  r0 = VMBLEND(0x22, a0, s0); r1 = VMBLEND(0x22, a1, s1); 
  r2 = VMBLEND(0x22, a2, s2); r3 = VMBLEND(0x22, a3, s3); 
  r4 = VMBLEND(0x22, a4, s4); r5 = VMBLEND(0x22, a5, s5); 
  r6 = VMBLEND(0x22, a6, s6); r7 = VMBLEND(0x22, a7, s7); 
  r8 = VMBLEND(0x22, a8, s8); r9 = VMBLEND(0x22, a9, s9);
  r10 = VMBLEND(0x22, a10, s10); r11 = VMBLEND(0x22, a11, s11); 

  // r =  < D  | C  | F-E+2p | 2B+A >
  //        2p | 2p | 4p     | 6p 
  r0 = VMBLEND(0x11, r0, t0); r1 = VMBLEND(0x11, r1, t1); 
  r2 = VMBLEND(0x11, r2, t2); r3 = VMBLEND(0x11, r3, t3); 
  r4 = VMBLEND(0x11, r4, t4); r5 = VMBLEND(0x11, r5, t5); 
  r6 = VMBLEND(0x11, r6, t6); r7 = VMBLEND(0x11, r7, t7); 
  r8 = VMBLEND(0x11, r8, t8); r9 = VMBLEND(0x11, r9, t9);
  r10 = VMBLEND(0x11, r10, t10); r11 = VMBLEND(0x11, r11, t11); 

  // TODO: try to optimze the following part (3 carry propagations)

  // r = r - 2p
  r0 = VSUB(r0, vp0); r1 = VSUB(r1, vp1);
  r2 = VSUB(r2, vp2); r3 = VSUB(r3, vp3);
  r4 = VSUB(r4, vp4); r5 = VSUB(r5, vp5);
  r6 = VSUB(r6, vp6); r7 = VSUB(r7, vp7);
  r8 = VSUB(r8, vp8); r9 = VSUB(r9, vp9);
  r10 = VSUB(r10, vp10); r11 = VSUB(r11, vp11);

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
  r0 = VADD(r0, VAND(vp0, smask)); r1 = VADD(r1, VAND(vp1, smask)); 
  r2 = VADD(r2, VAND(vp2, smask)); r3 = VADD(r3, VAND(vp3, smask)); 
  r4 = VADD(r4, VAND(vp4, smask)); r5 = VADD(r5, VAND(vp5, smask)); 
  r6 = VADD(r6, VAND(vp6, smask)); r7 = VADD(r7, VAND(vp7, smask)); 
  r8 = VADD(r8, VAND(vp8, smask)); r9 = VADD(r9, VAND(vp9, smask)); 
  r10 = VADD(r10, VAND(vp10, smask)); r11 = VADD(r11, VAND(vp11, smask)); 
  // r =  < D  | C  | F-E | 2B+A >
  //        2p | 2p | 2p  | 4p    

  // r = r - 2p
  r0 = VSUB(r0, vp0); r1 = VSUB(r1, vp1);
  r2 = VSUB(r2, vp2); r3 = VSUB(r3, vp3);
  r4 = VSUB(r4, vp4); r5 = VSUB(r5, vp5);
  r6 = VSUB(r6, vp6); r7 = VSUB(r7, vp7);
  r8 = VSUB(r8, vp8); r9 = VSUB(r9, vp9);
  r10 = VSUB(r10, vp10); r11 = VSUB(r11, vp11);

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
  r0 = VADD(r0, VAND(vp0, smask)); r1 = VADD(r1, VAND(vp1, smask)); 
  r2 = VADD(r2, VAND(vp2, smask)); r3 = VADD(r3, VAND(vp3, smask)); 
  r4 = VADD(r4, VAND(vp4, smask)); r5 = VADD(r5, VAND(vp5, smask)); 
  r6 = VADD(r6, VAND(vp6, smask)); r7 = VADD(r7, VAND(vp7, smask)); 
  r8 = VADD(r8, VAND(vp8, smask)); r9 = VADD(r9, VAND(vp9, smask)); 
  r10 = VADD(r10, VAND(vp10, smask)); r11 = VADD(r11, VAND(vp11, smask)); 
  // r =  < D  | C  | F-E | 2B+A >
  //        2p | 2p | 2p  | 2p 

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

static void fp2postmix_8x1x1w(vf2elm_t r, const vf2elm_t a, const vf2elm_t b)
{
  fppostmix_8x1w(r[0], a[0], b[0]);
  fppostmix_8x1w(r[1], a[1], b[1]);
}

// Here we use the algorithm presented in [HEY20] and try to keep our implementation 
// as same as [HEY20] implementation. 
void xDBLADD_2x4x1x1w(vf2elm_t x3z3x2z2, vf2elm_t z1x1_A)
{
  vf2elm_t z2x2x2z2, x5z5x4z4, z1x1x2A, x5z5z2z4, x7z7x6z6;

  mp2_hadamard_8x1x1w(x3z3x2z2, x3z3x2z2);            // x3z3x2z2 = x3+z3 | x3-z3 | x2+z2 | x2-z2
  fp2align2_8x1x1w(z2x2x2z2, x3z3x2z2);               // z2x2x2z2 = x2-z2 | x2+z2 | x2+z2 | x2-z2
  fp2mul_mont_8x1x1w(x3z3x2z2, x3z3x2z2, z2x2x2z2);   // x3z3x2z2 = (x3+z3)(x2-z2) | (x3-z3)(x2+z2) | (x2+z2)^2 | (x2-z2)^2
  mp2_hadamard_8x1x1w(x3z3x2z2, x3z3x2z2);            // x3z3x2z2 = 2(x2x3-z2z3) | 2(x2z3-x3z2) | 2(x2^2+z2^2) | 4x2z2
  fp2sqr_mont_8x1x1w(x5z5x4z4, x3z3x2z2);             // x5z5x4z4 = 4(x2x3-z2z3)^2 | 4(x2z3-x3z2)^2 | 4(x2^2+z2^2)^2 | 16 x2^2 z2^2
  fp2premix_8x1x1w(z1x1x2A, x5z5z2z4, z1x1_A, x3z3x2z2, x5z5x4z4); // z1x1x2A = z1 | x1 | 2(x2^2+z2^2) | A, x5z5z2z4 = 4(x2x3-z2z3)^2 | 4(x2z3-x3z2)^2 | 4x2z2 | 16 x2^2 z2^2
  fp2mul_mont_8x1x1w(x7z7x6z6, z1x1x2A, x5z5z2z4);    // x7z7x6z6 = z1(x2x3-z2z3)^2 | x1(x2z3-x3z2)^2 | 2x2z2(x2^2+z2^2) | 4A x2^2 z2^2
  fp2postmix_8x1x1w(x3z3x2z2, x7z7x6z6, x5z5x4z4);    // x3z3x2z2 = z1(x2x3-z2z3)^2 | x1(x2z3-x3z2)^2 | 4(x2^2-z2^2)^2 | 4x2z2(x2^2+z2^2+A x2^2 z2^2)
}

// -----------------------------------------------------------------------------
// (1x2x2x2)-way curve/isogeny arithmetic (based on (2x2x2)-way Fp2 arithmetic).
// NOTE: (1x2x2x2)-way means each function performs 1 curve/isogeny operation, 
// where each curve/isogeny operation performs 2 Fp2 operations, and each Fp2 
// operation performs 2 Fp operations and each Fp operation uses 2 lanes. 

// INPUT  a = < X1'     | X1    | X0'     | X0    | Z1'     | Z1    | Z0'     | Z0    >  all in [0, 2p)
// OUTPUT r = < X1'+Z1' | X1+Z1 | X0'+Z0' | X0+Z0 | X1'-Z1' | X1-Z1 | X0'-Z0' | X0-Z0 >  all in [0, 4p)
static void mp2_hadamard_2x2x2w(vgelm_t r, const vgelm_t a)
{
  __m512i a0 = a[0], a1 = a[1], a2 = a[2], a3 = a[3], a4 = a[4], a5 = a[5];
  __m512i r0, r1, r2, r3, r4, r5, c, carry;
  __m512i t0, t1, t2, t3, t4, t5;
  const __m512i vp0 = VSET(vp610x2[6],  vp610x2[0], vp610x2[6],  vp610x2[0], vp610x2[6],  vp610x2[0], vp610x2[6],  vp610x2[0]);
  const __m512i vp1 = VSET(vp610x2[7],  vp610x2[1], vp610x2[7],  vp610x2[1], vp610x2[7],  vp610x2[1], vp610x2[7],  vp610x2[1]);
  const __m512i vp2 = VSET(vp610x2[8],  vp610x2[2], vp610x2[8],  vp610x2[2], vp610x2[8],  vp610x2[2], vp610x2[8],  vp610x2[2]);
  const __m512i vp3 = VSET(vp610x2[9],  vp610x2[3], vp610x2[9],  vp610x2[3], vp610x2[9],  vp610x2[3], vp610x2[9],  vp610x2[3]);
  const __m512i vp4 = VSET(vp610x2[10], vp610x2[4], vp610x2[10], vp610x2[4], vp610x2[10], vp610x2[4], vp610x2[10], vp610x2[4]);
  const __m512i vp5 = VSET(vp610x2[11], vp610x2[5], vp610x2[11], vp610x2[5], vp610x2[11], vp610x2[5], vp610x2[11], vp610x2[5]);
  const __m512i vbmask = VSET1(VBMASK); 

  // t = <Z1' | Z1 | Z0' | Z0 | X1' | X1 | X0' | X0>
  t0 = VALIGNR(a0, a0, 4); t1 = VALIGNR(a1, a1, 4);
  t2 = VALIGNR(a2, a2, 4); t3 = VALIGNR(a3, a3, 4);
  t4 = VALIGNR(a4, a4, 4); t5 = VALIGNR(a5, a5, 4);

  // r = <X1' | X1 | X0' | X0 | 2p' | 2p | 2p' | 2p>
  r0 = VMBLEND(0x0F, a0, vp0); r1 = VMBLEND(0x0F, a1, vp1);
  r2 = VMBLEND(0x0F, a2, vp2); r3 = VMBLEND(0x0F, a3, vp3);
  r4 = VMBLEND(0x0F, a4, vp4); r5 = VMBLEND(0x0F, a5, vp5);

  // r = <X1'+Z1' | X1+Z1 | X0'+Z0' | X0+Z0 | X1'+2p' | X1+2p | X0'+2p' | X0+2p>
  r0 = VADD(r0, t0); r1 = VADD(r1, t1);
  r2 = VADD(r2, t2); r3 = VADD(r3, t3);
  r4 = VADD(r4, t4); r5 = VADD(r5, t5);

  // r = <X1'+Z1' | X1+Z1 | X0'+Z0' | X0+Z0 | X1'-Z1' | X1-Z1 | X0'-Z0' | X0-Z0>
  r0 = VMSUB(r0, 0x0F, r0, a0); r1 = VMSUB(r1, 0x0F, r1, a1); 
  r2 = VMSUB(r2, 0x0F, r2, a2); r3 = VMSUB(r3, 0x0F, r3, a3); 
  r4 = VMSUB(r4, 0x0F, r4, a4); r5 = VMSUB(r5, 0x0F, r5, a5); 

  // *simple* carry propagation
  // some limbs are finally 52-bit not 51-bit 

  // c = VSRA(r0, VBRADIX); r0 = VAND(r0, vbmask); r1 = VADD(r1, c);
  // c = VSRA(r1, VBRADIX); r1 = VAND(r1, vbmask); r2 = VADD(r2, c);
  // c = VSRA(r2, VBRADIX); r2 = VAND(r2, vbmask); r3 = VADD(r3, c);
  // c = VSRA(r3, VBRADIX); r3 = VAND(r3, vbmask); r4 = VADD(r4, c);
  // c = VSRA(r4, VBRADIX); r4 = VAND(r4, vbmask); r5 = VADD(r5, c);
  // c = VSRA(r5, VBRADIX); r5 = VAND(r5, vbmask); 
  // r0 = VMADD(r0, 0xAA, r0, VSHUF(c, 0x4E));

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

  r[0] = r0; r[1] = r1; r[2] = r2; r[3] = r3; r[4] = r4; r[5] = r5; 
}

// INPUT  a = < C1' | C1 | C0' | C0 | A1' | A1 | A0' | A0 >
//        b = < D1' | D1 | D0' | D0 | B1' | B1 | B0' | B0 >
// OUTPUT r = < B1' | B1 | B0' | B0 | A1' | A1 | A0' | A0 >
static void fp2mix1_2x2x2w(vgelm_t r, const vgelm_t a, const vgelm_t b)
{
  __m512i a0 = a[0], a1 = a[1], a2 = a[2], a3 = a[3], a4 = a[4], a5 = a[5];
  __m512i b0 = b[0], b1 = b[1], b2 = b[2], b3 = b[3], b4 = b[4], b5 = b[5];
  __m512i r0, r1, r2, r3, r4, r5;
  __m512i t0, t1, t2, t3, t4, t5;
  
  // t = <B1' | B1 | B0' | B0 | D1' | D1 | D0' | D0>
  t0 = VALIGNR(b0, b0, 4); t1 = VALIGNR(b1, b1, 4);
  t2 = VALIGNR(b2, b2, 4); t3 = VALIGNR(b3, b3, 4);
  t4 = VALIGNR(b4, b4, 4); t5 = VALIGNR(b5, b5, 4);

  // r = < B1' | B1 | B0' | B0 | A1' | A1 | A0' | A0 >
  r0 = VMBLEND(0x0F, t0, a0); r1 = VMBLEND(0x0F, t1, a1);
  r2 = VMBLEND(0x0F, t2, a2); r3 = VMBLEND(0x0F, t3, a3);
  r4 = VMBLEND(0x0F, t4, a4); r5 = VMBLEND(0x0F, t5, a5);

  r[0] = r0; r[1] = r1; r[2] = r2; r[3] = r3; r[4] = r4; r[5] = r5; 
}

// INPUT  a = < E1' | E1 | E0' | E0 | A1' | A1 | A0' | A0 >
//        b = < F1' | F1 | F0' | F0 | B1' | B1 | B0' | B0 >
//        c = < G1' | G1 | G0' | G0 | C1' | C1 | C0' | C0 >
//        d = < H1' | H1 | H0' | H0 | D1' | D1 | D0' | D0 >
// OUTPUT r = < E1' | E1 | E0' | E0 | B1' | B1 | B0' | B0 >
//        s = < G1' | G1 | G0' | G0 | H1' | H1 | H0' | H0 >
static void fp2mix2_2x2x2w(vgelm_t r, vgelm_t s, const vgelm_t a, const vgelm_t b, const vgelm_t c, const vgelm_t d)
{
  __m512i a0 = a[0], a1 = a[1], a2 = a[2], a3 = a[3], a4 = a[4], a5 = a[5];
  __m512i b0 = b[0], b1 = b[1], b2 = b[2], b3 = b[3], b4 = b[4], b5 = b[5];
  __m512i c0 = c[0], c1 = c[1], c2 = c[2], c3 = c[3], c4 = c[4], c5 = c[5];
  __m512i d0 = d[0], d1 = d[1], d2 = d[2], d3 = d[3], d4 = d[4], d5 = d[5];
  __m512i r0, r1, r2, r3, r4, r5;
  __m512i s0, s1, s2, s3, s4, s5;

  // r = < E1' | E1 | E0' | E0 | B1' | B1 | B0' | B0 >
  r0 = VMBLEND(0x0F, a0, b0); r1 = VMBLEND(0x0F, a1, b1);
  r2 = VMBLEND(0x0F, a2, b2); r3 = VMBLEND(0x0F, a3, b3);
  r4 = VMBLEND(0x0F, a4, b4); r5 = VMBLEND(0x0F, a5, b5);

  // s = < D1' | D1 | D0' | D0 | H1' | H1 | H0' | H0 >
  s0 = VALIGNR(d0, d0, 4); s1 = VALIGNR(d1, d1, 4);
  s2 = VALIGNR(d2, d2, 4); s3 = VALIGNR(d3, d3, 4);
  s4 = VALIGNR(d4, d4, 4); s5 = VALIGNR(d5, d5, 4);

  // s = < G1' | G1 | G0' | G0 | H1' | H1 | H0' | H0 >
  s0 = VMBLEND(0x0F, c0, s0); s1 = VMBLEND(0x0F, c1, s1); 
  s2 = VMBLEND(0x0F, c2, s2); s3 = VMBLEND(0x0F, c3, s3); 
  s4 = VMBLEND(0x0F, c4, s4); s5 = VMBLEND(0x0F, c5, s5); 

  r[0] = r0; r[1] = r1; r[2] = r2; r[3] = r3; r[4] = r4; r[5] = r5; 
  s[0] = s0; s[1] = s1; s[2] = s2; s[3] = s3; s[4] = s4; s[5] = s5;
}

// INPUT
// P           : < XP1' | XP1 | XP0' | XP0 | ZP1' | ZP1 | ZP0' | ZP0 >
// C24_A24plus : < C1'  | C1  | C0'  | C0  | A1'  | A1  | A0'  | A0  >
// OUTPUT 
// Q           : < XQ1' | XQ1 | XQ0' | XQ0 | ZQ1' | ZQ1 | ZQ0' | ZQ0 >
void xDBL_1x2x2x2w(const vgelm_t P, vgelm_t Q, const vgelm_t C24_A24plus)
{
  vgelm_t t1, t2, t3, t4;

  mp2_hadamard_2x2x2w(t1, P);             // t1 = X+Z               | X-Z                   [4p]
  fp2sqr_mont_2x2x2w(t2, t1);             // t2 = (X+Z)^2           | (X-Z)^2               [2p]
  mp2_hadamard_2x2x2w(t1, t2);            // t1 = ((X+Z)^2+(X-Z)^2) | 4*X*Z                 [4p]
  fp2mix1_2x2x2w(t3, t1, t2);             // t3 = (X-Z)^2           | 4XZ                   [2p|4p]        
  fp2mul_mont_2x2x2w(t3, C24_A24plus, t3);// t3 = C(X-Z)^2          | 4AXZ                  [2p]
  mp2_hadamard_2x2x2w(t4, t3);            // t4 = C(X-Z)^2+4AXZ     | C(X-Z)^2-4AXZ         [4p]
  fp2mix2_2x2x2w(t1, t2, t2, t1, t3, t4); // t1 = (X+Z)^2           | 4XZ                   [2p|4p]
                                          // t2 = C(X-Z)^2          | C(X-Z)^2+4AXZ         [2p|4p]
  fp2mul_mont_2x2x2w(Q, t1, t2);          // Q  = C(X-Z)^2(X+Z)^2   | (C(X-Z)^2+4AXZ)*4XZ   [2p]
}

// INPUT  a = < C1' | C1 | C0' | C0 | A1' | A1 | A0' | A0 >               all in [0, 2p)
//        b = < D1' | D1 | D0' | D0 | B1' | B1 | B0' | B0 >               all in [0, 2p)
// OUTPUT r = < C1+D1' | C1+D1 | C0+D0' | C0+D0 | A1-B1' | A1-B1 | A0-B0' | A0-B0 >   all in [0, 4p)
static void mp2_addsub_2x2x2w(vgelm_t r, const vgelm_t a, const vgelm_t b)
{
  __m512i a0 = a[0], a1 = a[1], a2 = a[2], a3 = a[3], a4 = a[4], a5 = a[5];
  __m512i b0 = b[0], b1 = b[1], b2 = b[2], b3 = b[3], b4 = b[4], b5 = b[5];
  __m512i r0, r1, r2, r3, r4, r5, c, carry;
  __m512i t0, t1, t2, t3, t4, t5;
  const __m512i vp0 = VSET(vp610x2[6],  vp610x2[0], vp610x2[6],  vp610x2[0], vp610x2[6],  vp610x2[0], vp610x2[6],  vp610x2[0]);
  const __m512i vp1 = VSET(vp610x2[7],  vp610x2[1], vp610x2[7],  vp610x2[1], vp610x2[7],  vp610x2[1], vp610x2[7],  vp610x2[1]);
  const __m512i vp2 = VSET(vp610x2[8],  vp610x2[2], vp610x2[8],  vp610x2[2], vp610x2[8],  vp610x2[2], vp610x2[8],  vp610x2[2]);
  const __m512i vp3 = VSET(vp610x2[9],  vp610x2[3], vp610x2[9],  vp610x2[3], vp610x2[9],  vp610x2[3], vp610x2[9],  vp610x2[3]);
  const __m512i vp4 = VSET(vp610x2[10], vp610x2[4], vp610x2[10], vp610x2[4], vp610x2[10], vp610x2[4], vp610x2[10], vp610x2[4]);
  const __m512i vp5 = VSET(vp610x2[11], vp610x2[5], vp610x2[11], vp610x2[5], vp610x2[11], vp610x2[5], vp610x2[11], vp610x2[5]);
  const __m512i vbmask = VSET1(VBMASK); 

  // t = < D1' | D1 | D0' | D0 | 2p' | 2p | 2p' | 2p >
  t0 = VMBLEND(0x0F, b0, vp0); t1 = VMBLEND(0x0F, b1, vp1);
  t2 = VMBLEND(0x0F, b2, vp2); t3 = VMBLEND(0x0F, b3, vp3);
  t4 = VMBLEND(0x0F, b4, vp4); t5 = VMBLEND(0x0F, b5, vp5);

  // r = < C1+D1' | C1+D1 | C0+D0' | C0+D0 | A1+2p' | A1+2p | A0+2p' | A0+2p >
  r0 = VADD(a0, t0); r1 = VADD(a1, t1);
  r2 = VADD(a2, t2); r3 = VADD(a3, t3);
  r4 = VADD(a4, t4); r5 = VADD(a5, t5);

  // r = < C1+D1' | C1+D1 | C0+D0' | C0+D0 | A1+2p-B1' | A1+2p-B1 | A0+2p-B0' | A0+2p-B0 >
  r0 = VMSUB(r0, 0x0F, r0, b0); r1 = VMSUB(r1, 0x0F, r1, b1);
  r2 = VMSUB(r2, 0x0F, r2, b2); r3 = VMSUB(r3, 0x0F, r3, b3);
  r4 = VMSUB(r4, 0x0F, r4, b4); r5 = VMSUB(r5, 0x0F, r5, b5);

  // *simple* carry propagation
  // some limbs are finally 52-bit not 51-bit 

  // c = VSRA(r0, VBRADIX); r0 = VAND(r0, vbmask); r1 = VADD(r1, c);
  // c = VSRA(r1, VBRADIX); r1 = VAND(r1, vbmask); r2 = VADD(r2, c);
  // c = VSRA(r2, VBRADIX); r2 = VAND(r2, vbmask); r3 = VADD(r3, c);
  // c = VSRA(r3, VBRADIX); r3 = VAND(r3, vbmask); r4 = VADD(r4, c);
  // c = VSRA(r4, VBRADIX); r4 = VAND(r4, vbmask); r5 = VADD(r5, c);
  // c = VSRA(r5, VBRADIX); r5 = VAND(r5, vbmask); 
  // r0 = VMADD(r0, 0xAA, r0, VSHUF(c, 0x4E));

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

  r[0] = r0; r[1] = r1; r[2] = r2; r[3] = r3; r[4] = r4; r[5] = r5; 
}

// INPUT  a = < C1' | C1 | C0' | C0 | A1' | A1 | A0' | A0 >               all in [0, 2p)
//        b = < D1' | D1 | D0' | D0 | B1' | B1 | B0' | B0 >               all in [0, 2p)
// OUTPUT r = < C1+D1' | C1+D1 | C0+D0' | C0+D0 | A1-B1' | A1-B1 | A0-B0' | A0-B0 >   all in [0, 2p)
static void fp2addsub_2x2x2w(vgelm_t r, const vgelm_t a, const vgelm_t b)
{
  __m512i a0 = a[0], a1 = a[1], a2 = a[2], a3 = a[3], a4 = a[4], a5 = a[5];
  __m512i b0 = b[0], b1 = b[1], b2 = b[2], b3 = b[3], b4 = b[4], b5 = b[5];
  __m512i r0, r1, r2, r3, r4, r5, c, smask, carry;
  __m512i t0, t1, t2, t3, t4, t5;
  const __m512i vp0 = VSET(vp610x2[6],  vp610x2[0], vp610x2[6],  vp610x2[0], vp610x2[6],  vp610x2[0], vp610x2[6],  vp610x2[0]);
  const __m512i vp1 = VSET(vp610x2[7],  vp610x2[1], vp610x2[7],  vp610x2[1], vp610x2[7],  vp610x2[1], vp610x2[7],  vp610x2[1]);
  const __m512i vp2 = VSET(vp610x2[8],  vp610x2[2], vp610x2[8],  vp610x2[2], vp610x2[8],  vp610x2[2], vp610x2[8],  vp610x2[2]);
  const __m512i vp3 = VSET(vp610x2[9],  vp610x2[3], vp610x2[9],  vp610x2[3], vp610x2[9],  vp610x2[3], vp610x2[9],  vp610x2[3]);
  const __m512i vp4 = VSET(vp610x2[10], vp610x2[4], vp610x2[10], vp610x2[4], vp610x2[10], vp610x2[4], vp610x2[10], vp610x2[4]);
  const __m512i vp5 = VSET(vp610x2[11], vp610x2[5], vp610x2[11], vp610x2[5], vp610x2[11], vp610x2[5], vp610x2[11], vp610x2[5]);
  const __m512i vbmask = VSET1(VBMASK); 

  // r = < C1+D1' | C1+D1 | C0+D0' | C0+D0 | A1+0' | A1+0 | A0+0' | A0+0 >
  r0 = VMADD(a0, 0xF0, a0, b0); r1 = VMADD(a1, 0xF0, a1, b1);
  r2 = VMADD(a2, 0xF0, a2, b2); r3 = VMADD(a3, 0xF0, a3, b3);
  r4 = VMADD(a4, 0xF0, a4, b4); r5 = VMADD(a5, 0xF0, a5, b5);

  // t = < 2p' | 2p | 2p' | 2p | B1' | B1 | B0' | B0 >
  t0 = VMBLEND(0xF0, b0, vp0); t1 = VMBLEND(0xF0, b1, vp1);
  t2 = VMBLEND(0xF0, b2, vp2); t3 = VMBLEND(0xF0, b3, vp3);
  t4 = VMBLEND(0xF0, b4, vp4); t5 = VMBLEND(0xF0, b5, vp5);

  // r = < C1+D1-2p' | C1+D1-2p | C0+D0-2p' | C0+D0-2p | A1-B1' | A1-B1 | A0-B0' | A0-B0 >
  r0 = VSUB(r0, t0); r1 = VSUB(r1, t1);
  r2 = VSUB(r2, t2); r3 = VSUB(r3, t3);
  r4 = VSUB(r4, t4); r5 = VSUB(r5, t5);

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

  // c = VSRA(r0, VBRADIX); r0 = VAND(r0, vbmask); r1 = VADD(r1, c);
  // c = VSRA(r1, VBRADIX); r1 = VAND(r1, vbmask); r2 = VADD(r2, c);
  // c = VSRA(r2, VBRADIX); r2 = VAND(r2, vbmask); r3 = VADD(r3, c);
  // c = VSRA(r3, VBRADIX); r3 = VAND(r3, vbmask); r4 = VADD(r4, c);
  // c = VSRA(r4, VBRADIX); r4 = VAND(r4, vbmask); r5 = VADD(r5, c);
  // c = VSRA(r5, VBRADIX); r5 = VAND(r5, vbmask); 
  // r0 = VMADD(r0, 0xAA, r0, VSHUF(c, 0x4E));

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

  r[0] = r0; r[1] = r1; r[2] = r2; r[3] = r3; r[4] = r4; r[5] = r5; 
}

// INPUT  a = < C1' | C1 | C0' | C0 | A1' | A1 | A0' | A0 >
//        b = < D1' | D1 | D0' | D0 | B1' | B1 | B0' | B0 >
// OUTPUT r = < C1' | C1 | C0' | C0 | B1' | B1 | B0' | B0 >
//        s = < C1' | C1 | C0' | C0 | D1' | D1 | D0' | D0 >
static void fp2mix3_2x2x2w(vgelm_t r, vgelm_t s, const vgelm_t a, const vgelm_t b)
{
  __m512i a0 = a[0], a1 = a[1], a2 = a[2], a3 = a[3], a4 = a[4], a5 = a[5];
  __m512i b0 = b[0], b1 = b[1], b2 = b[2], b3 = b[3], b4 = b[4], b5 = b[5];
  __m512i r0, r1, r2, r3, r4, r5;
  __m512i s0, s1, s2, s3, s4, s5;

  // r = < C1' | C1 | C0' | C0 | B1' | B1 | B0' | B0 >
  r0 = VMBLEND(0x0F, a0, b0); r1 = VMBLEND(0x0F, a1, b1); 
  r2 = VMBLEND(0x0F, a2, b2); r3 = VMBLEND(0x0F, a3, b3); 
  r4 = VMBLEND(0x0F, a4, b4); r5 = VMBLEND(0x0F, a5, b5); 

  // s = < B1' | B1 | B0' | B0 | D1' | D1 | D0' | D0 >
  s0 = VALIGNR(b0, b0, 4); s1 = VALIGNR(b1, b1, 4);
  s2 = VALIGNR(b2, b2, 4); s3 = VALIGNR(b3, b3, 4);
  s4 = VALIGNR(b4, b4, 4); s5 = VALIGNR(b5, b5, 4);

  // s = < C1' | C1 | C0' | C0 | D1' | D1 | D0' | D0 >
  s0 = VMBLEND(0x0F, a0, s0); s1 = VMBLEND(0x0F, a1, s1); 
  s2 = VMBLEND(0x0F, a2, s2); s3 = VMBLEND(0x0F, a3, s3); 
  s4 = VMBLEND(0x0F, a4, s4); s5 = VMBLEND(0x0F, a5, s5); 

  r[0] = r0; r[1] = r1; r[2] = r2; r[3] = r3; r[4] = r4; r[5] = r5; 
  s[0] = s0; s[1] = s1; s[2] = s2; s[3] = s3; s[4] = s4; s[5] = s5;
}

// INPUT  a = < C1' | C1 | C0' | C0 | A1' | A1 | A0' | A0 >
//        b = < D1' | D1 | D0' | D0 | B1' | B1 | B0' | B0 >
// OUTPUT r = < C1' | C1 | C0' | C0 | B1' | B1 | B0' | B0 >
static void vec_blend_4x2w(vfelm_t r, const vfelm_t a, const vfelm_t b, const __mmask8 mask)
{
  __m512i a0 = a[0], a1 = a[1], a2 = a[2], a3 = a[3], a4 = a[4], a5 = a[5];
  __m512i b0 = b[0], b1 = b[1], b2 = b[2], b3 = b[3], b4 = b[4], b5 = b[5];
  __m512i r0, r1, r2, r3, r4, r5;

  r0 = VMBLEND(mask, a0, b0); r1 = VMBLEND(mask, a1, b1); 
  r2 = VMBLEND(mask, a2, b2); r3 = VMBLEND(mask, a3, b3); 
  r4 = VMBLEND(mask, a4, b4); r5 = VMBLEND(mask, a5, b5); 

  r[0] = r0; r[1] = r1; r[2] = r2; r[3] = r3; r[4] = r4; r[5] = r5; 
}

// INPUT  a = < D1' | D1 | D0' | D0 | A1' | A1 | A0' | A0 >
//        b = < E1' | E1 | E0' | E0 | B1' | B1 | B0' | B0 >
//        c = < F1' | F1 | F0' | F0 | C1' | C1 | C0' | C0 >
// OUTPUT r = < D1' | D1 | D0' | D0 | F1' | F1 | F0' | F0 >
//        s = < E1' | E1 | E0' | E0 | F1' | F1 | F0' | F0 >
static void fp2mix4_2x2x2w(vgelm_t r, vgelm_t s, const vgelm_t a, const vgelm_t b, const vgelm_t c)
{
  __m512i a0 = a[0], a1 = a[1], a2 = a[2], a3 = a[3], a4 = a[4], a5 = a[5];
  __m512i b0 = b[0], b1 = b[1], b2 = b[2], b3 = b[3], b4 = b[4], b5 = b[5];
  __m512i c0 = c[0], c1 = c[1], c2 = c[2], c3 = c[3], c4 = c[4], c5 = c[5];
  __m512i r0, r1, r2, r3, r4, r5;
  __m512i s0, s1, s2, s3, s4, s5;

  // c = < C1' | C1 | C0' | C0 | F1' | F1 | F0' | F0 >
  c0 = VALIGNR(c0, c0, 4); c1 = VALIGNR(c1, c1, 4);
  c2 = VALIGNR(c2, c2, 4); c3 = VALIGNR(c3, c3, 4);
  c4 = VALIGNR(c4, c4, 4); c5 = VALIGNR(c5, c5, 4);

  // r = < D1' | D1 | D0' | D0 | F1' | F1 | F0' | F0 >
  r0 = VMBLEND(0x0F, a0, c0); r1 = VMBLEND(0x0F, a1, c1); 
  r2 = VMBLEND(0x0F, a2, c2); r3 = VMBLEND(0x0F, a3, c3); 
  r4 = VMBLEND(0x0F, a4, c4); r5 = VMBLEND(0x0F, a5, c5); 

  // s = < E1' | E1 | E0' | E0 | F1' | F1 | F0' | F0 >
  s0 = VMBLEND(0x0F, b0, c0); s1 = VMBLEND(0x0F, b1, c1); 
  s2 = VMBLEND(0x0F, b2, c2); s3 = VMBLEND(0x0F, b3, c3); 
  s4 = VMBLEND(0x0F, b4, c4); s5 = VMBLEND(0x0F, b5, c5); 

  r[0] = r0; r[1] = r1; r[2] = r2; r[3] = r3; r[4] = r4; r[5] = r5; 
  s[0] = s0; s[1] = s1; s[2] = s2; s[3] = s3; s[4] = s4; s[5] = s5;
}

// INPUT  a = < C1' | C1 | C0' | C0 | A1' | A1 | A0' | A0 >
//        b = < D1' | D1 | D0' | D0 | B1' | B1 | B0' | B0 >
// OUTPUT r = < C1' | C1 | C0' | C0 | D1' | D1 | D0' | D0 >
//        s = < A1' | A1 | A0' | A0 | B1' | B1 | B0' | B0 >
static void fp2mix5_2x2x2w(vgelm_t r, vgelm_t s, const vgelm_t a, const vgelm_t b)
{
  __m512i a0 = a[0], a1 = a[1], a2 = a[2], a3 = a[3], a4 = a[4], a5 = a[5];
  __m512i b0 = b[0], b1 = b[1], b2 = b[2], b3 = b[3], b4 = b[4], b5 = b[5];
  __m512i r0, r1, r2, r3, r4, r5;
  __m512i s0, s1, s2, s3, s4, s5;

  // r = < B1' | B1 | B0' | B0 | D1' | D1 | D0' | D0 >
  r0 = VALIGNR(b0, b0, 4); r1 = VALIGNR(b1, b1, 4);
  r2 = VALIGNR(b2, b2, 4); r3 = VALIGNR(b3, b3, 4);
  r4 = VALIGNR(b4, b4, 4); r5 = VALIGNR(b5, b5, 4);

  // r = < C1' | C1 | C0' | C0 | D1' | D1 | D0' | D0 >
  r0 = VMBLEND(0x0F, a0, r0); r1 = VMBLEND(0x0F, a1, r1); 
  r2 = VMBLEND(0x0F, a2, r2); r3 = VMBLEND(0x0F, a3, r3); 
  r4 = VMBLEND(0x0F, a4, r4); r5 = VMBLEND(0x0F, a5, r5); 

  // s = < A1' | A1 | A0' | A0 | C1' | C1 | C0' | C0 >
  s0 = VALIGNR(a0, a0, 4); s1 = VALIGNR(a1, a1, 4);
  s2 = VALIGNR(a2, a2, 4); s3 = VALIGNR(a3, a3, 4);
  s4 = VALIGNR(a4, a4, 4); s5 = VALIGNR(a5, a5, 4);

  // r = < A1' | A1 | A0' | A0 | B1' | B1 | B0' | B0 >
  s0 = VMBLEND(0x0F, s0, b0); s1 = VMBLEND(0x0F, s1, b1); 
  s2 = VMBLEND(0x0F, s2, b2); s3 = VMBLEND(0x0F, s3, b3); 
  s4 = VMBLEND(0x0F, s4, b4); s5 = VMBLEND(0x0F, s5, b5); 

  r[0] = r0; r[1] = r1; r[2] = r2; r[3] = r3; r[4] = r4; r[5] = r5; 
  s[0] = s0; s[1] = s1; s[2] = s2; s[3] = s3; s[4] = s4; s[5] = s5;
}

// INPUT  a = < C1' | C1 | C0' | C0 | A1' | A1 | A0' | A0 >
//        b = < D1' | D1 | D0' | D0 | B1' | B1 | B0' | B0 >
// OUTPUT r = < C1' | C1 | C0' | C0 | C1' | C1 | C0' | C0 >
//        s = < A1' | A1 | A0' | A0 | B1' | B1 | B0' | B0 >
static void fp2mix6_2x2x2w(vgelm_t r, vgelm_t s, const vgelm_t a, const vgelm_t b)
{
  __m512i a0 = a[0], a1 = a[1], a2 = a[2], a3 = a[3], a4 = a[4], a5 = a[5];
  __m512i b0 = b[0], b1 = b[1], b2 = b[2], b3 = b[3], b4 = b[4], b5 = b[5];
  __m512i r0, r1, r2, r3, r4, r5;
  __m512i s0, s1, s2, s3, s4, s5;
  __m512i mask = VSET(7, 6, 5, 4, 7, 6, 5, 4); 

  // s = < A1' | A1 | A0' | A0 | C1' | C1 | C0' | C0 >
  s0 = VALIGNR(a0, a0, 4); s1 = VALIGNR(a1, a1, 4);
  s2 = VALIGNR(a2, a2, 4); s3 = VALIGNR(a3, a3, 4);
  s4 = VALIGNR(a4, a4, 4); s5 = VALIGNR(a5, a5, 4);

  // s = < A1' | A1 | A0' | A0 | B1' | B1 | B0' | B0 >
  s0 = VMBLEND(0x0F, s0, b0); s1 = VMBLEND(0x0F, s1, b1); 
  s2 = VMBLEND(0x0F, s2, b2); s3 = VMBLEND(0x0F, s3, b3); 
  s4 = VMBLEND(0x0F, s4, b4); s5 = VMBLEND(0x0F, s5, b5); 

  // r = < C1' | C1 | C0' | C0 | C1' | C1 | C0' | C0 >
  r0 = VPERMV(mask, a0); r1 = VPERMV(mask, a1);
  r2 = VPERMV(mask, a2); r3 = VPERMV(mask, a3);
  r4 = VPERMV(mask, a4); r5 = VPERMV(mask, a5);

  r[0] = r0; r[1] = r1; r[2] = r2; r[3] = r3; r[4] = r4; r[5] = r5; 
  s[0] = s0; s[1] = s1; s[2] = s2; s[3] = s3; s[4] = s4; s[5] = s5;
}

// INPUT  a = < B1' | B1 | B0' | B0 | A1' | A1 | A0' | A0 >
// OUTPUT r = < B1' | B1 | B0' | B0 | B1' | B1 | B0' | B0 >
static void fp2align_2x2x2w(vgelm_t r, const vgelm_t a)
{
  __m512i a0 = a[0], a1 = a[1], a2 = a[2], a3 = a[3], a4 = a[4], a5 = a[5];
  __m512i r0, r1, r2, r3, r4, r5;
  __m512i mask = VSET(7, 6, 5, 4, 7, 6, 5, 4); 

  // r = < B1' | B1 | B0' | B0 | B1' | B1 | B0' | B0 >
  r0 = VPERMV(mask, a0); r1 = VPERMV(mask, a1);
  r2 = VPERMV(mask, a2); r3 = VPERMV(mask, a3);
  r4 = VPERMV(mask, a4); r5 = VPERMV(mask, a5);

  r[0] = r0; r[1] = r1; r[2] = r2; r[3] = r3; r[4] = r4; r[5] = r5; 
}

// INPUT
// P           : < XP1' | XP1 | XP0' | XP0 | ZP1' | ZP1 | ZP0' | ZP0 >
// C24_A24plus : <  C1' |  C1 |  C0' |  C0 |  A1' |  A1 |  A0' |  A0 >
// OUTPUT 
// Q           : < XQ1' | XQ1 | XQ0' | XQ0 | ZQ1' | ZQ1 | ZQ0' | ZQ0 >
void xTPL_1x2x2x2w(const vgelm_t P, vgelm_t Q, const vgelm_t C24_A24plus)
{
  vgelm_t t1, t2, t3, t4, t5;

  mp2_hadamard_2x2x2w(t1, P);             //  t1 = X+Z                | X-Z                 [4p]
  fp2sqr_mont_2x2x2w(t2, t1);             //  t2 = (X+Z)^2            | (X-Z)^2             [2p]
  fp2mix3_2x2x2w(t1, t3, P, t2);          //  t1 = X                  | (X-Z)^2             [2p]
                                          //  t3 = X                  | (X+Z)^2             [2p]
  fp2addsub_2x2x2w(t1, t3, t1);           //  t1 = 2X                 | 4XZ                 [2p]
  vec_blend_4x2w(t3, t2, t1, 0x0F);       //  t3 = (X+Z)^2            | 4XZ                 [2p]
  fpadd_4x2w(t3, t3, t3);                 //  t3 = 2(X+Z)^2           | 8XZ                 [2p]
  fp2mix1_2x2x2w(t4, t1, t2);             //  t4 = (X-Z)^2            | 4XZ                 [2p]
  fp2mul_mont_2x2x2w(t4, C24_A24plus, t4);//  t4 = C(X-Z)^2           | 4AXZ                [2p]
  fp2mix4_2x2x2w(t1, t2, t4, t2, t1);     //  t1 = C(X-Z)^2           | 2X                  [2p]
                                          //  t2 = (X+Z)^2            | 2X                  [2p]
  fp2mul_mont_2x2x2w(t1, t1, t2);         // !t1 = C(X^2-Z^2)^2       | 4X^2                [2p]
  vec_blend_4x2w(t2, t3, t1, 0x0F);       //  t2 = 2(X+Z)^2           | 4X^2                [2p]
  fpadd_4x2w(t2, t2, t3);                 //  t2 = 4(X+Z)^2           | 4X^2+8XZ            [2p]
  fp2mix5_2x2x2w(t3, t2, t4, t2);         //  t3 = C(X-Z)^2           | 4(X+Z)^2            [2p]
                                          //  t2 = 4AXZ               | 4X^2+8XZ            [2p]
  mp2_addsub_2x2x2w(t4, t3, t2);          // !t4 = C(X-Z)^2+4AXZ      | 4Z^2                [4p]
  fp2mix6_2x2x2w(t5, t2, t4, t1);         //  t5 = C(X-Z)^2+4AXZ      | C(X-Z)^2+4AXZ       [4p]
                                          // !t2 = 4Z^2               | 4X^2                [4p|2p]
  fp2mul_mont_2x2x2w(t5, t5, t2);         //  t5 = 4Z^2(C(X-Z)^2+4AXZ)| 4X^2(C(X-Z)^2+4AXZ) [2p]
  fp2align_2x2x2w(t3, t1);                // !t3 = C(X^2-Z^2)^2       | C(X^2-Z^2)^2        [2p]
  mp_sub_p2_4x2w(t5, t3, t5);             //  t5 = C(X^2-Z^2)^2-4Z^2(C(X-Z)^2+4AXZ) | C(X^2-Z^2)^2-4X^2(C(X-Z)^2+4AXZ)  [4p]
  fp2sqr_mont_2x2x2w(t5, t5);             // !t5 = (C(X^2-Z^2)^2-4Z^2(C(X-Z)^2+4AXZ))^2 | (C(X^2-Z^2)^2-4X^2(C(X-Z)^2+4AXZ))^2   [2p]
  fp2mul_mont_2x2x2w(Q, t5, P);           //   Q = X(C(X^2-Z^2)^2-4Z^2(C(X-Z)^2+4AXZ))^2| Z(C(X^2-Z^2)^2-4X^2(C(X-Z)^2+4AXZ))^2  [2p]
}

// INPUT  a = < B1' | B1 | B0' | B0 | A1' | A1 | A0' | A0 >
// OUTPUT r = < A1' | A1 | A0' | A0 | B1' | B1 | B0' | B0 >
static void vec_alignr4_4x2w(vgelm_t r, const vgelm_t a)
{  
  __m512i a0 = a[0], a1 = a[1], a2 = a[2], a3 = a[3], a4 = a[4], a5 = a[5];
  __m512i r0, r1, r2, r3, r4, r5;

  // r = < A1' | A1 | A0' | A0 | B1' | B1 | B0' | B0 >
  r0 = VALIGNR(a0, a0, 4); r1 = VALIGNR(a1, a1, 4);
  r2 = VALIGNR(a2, a2, 4); r3 = VALIGNR(a3, a3, 4);
  r4 = VALIGNR(a4, a4, 4); r5 = VALIGNR(a5, a5, 4);

  r[0] = r0; r[1] = r1; r[2] = r2; r[3] = r3; r[4] = r4; r[5] = r5; 
}

// INPUT  a = < B1' | B1 | B0' | B0 | A1' | A1 | A0' | A0 >
// OUTPUT r = < A1' | A1 | A0' | A0 | A1' | A1 | A0' | A0 >
static void fp2align2_2x2x2w(vgelm_t r, const vgelm_t a)
{
  __m512i a0 = a[0], a1 = a[1], a2 = a[2], a3 = a[3], a4 = a[4], a5 = a[5];
  __m512i r0, r1, r2, r3, r4, r5;
  __m512i mask = VSET(3, 2, 1, 0, 3, 2, 1, 0); 

  // r = < A1' | A1 | A0' | A0 | A1' | A1 | A0' | A0 >
  r0 = VPERMV(mask, a0); r1 = VPERMV(mask, a1);
  r2 = VPERMV(mask, a2); r3 = VPERMV(mask, a3);
  r4 = VPERMV(mask, a4); r5 = VPERMV(mask, a5);

  r[0] = r0; r[1] = r1; r[2] = r2; r[3] = r3; r[4] = r4; r[5] = r5; 
}

// INPUT
// P           : < XP1' | XP1 | XP0' | XP0 |   ZP1' |  ZP1 |  ZP0' |  ZP0 >
// OUTPUT 
// C24_A24plus : <   C1' |   C1 |   C0' |   C0 |    A1' |   A1 |   A0' |   A0 >
// coeff0      : <    %  |    % |    %  |    % |  cf01' | cf01 | cf00' | cf00 >
// coeff2_1    : < cf21' | cf21 | cf20' | cf20 |  cf11' | cf11 | cf10' | cf10 >
void get_4_isog_1x2x2x2w(const vgelm_t P, vgelm_t C24_A24plus, vgelm_t coeff__0, vgelm_t coeff2_1 )
{
  vgelm_t t1, t2;

  mp2_hadamard_2x2x2w(coeff2_1, P);     //    coeff2_1 = X+Z  | X-Z
  fp2sqr_mont_2x2x2w(t1, P);            //          t1 = X^2  | Z^2
  mp_add_4x2w(t1, t1, t1);              //          t1 = 2X^2 | 2Z^2
  fp2sqr_mont_2x2x2w(t2, t1);           //          t2 = 4X^4 | 4Z^4
  vec_alignr4_4x2w(C24_A24plus, t2);    // C24_A24plus = 4Z^4 | 4X^4
  mp_add_4x2w(coeff__0, t1, t1);        //      coeff0 =    % | 4Z^2 
}

// INPUT  a = < D1' | D1 | D0' | D0 | A1' | A1 | A0' | A0 >
//        b = < E1' | E1 | E0' | E0 | B1' | B1 | B0' | B0 >
//        c = < F1' | F1 | F0' | F0 | C1' | C1 | C0' | C0 >
// OUTPUT r = < E1' | E1 | E0' | E0 | A1' | A1 | A0' | A0 >
//        s = < C1' | C1 | C0' | C0 | D1' | D1 | D0' | D0 >
static void fp2mix8_2x2x2w(vgelm_t r, vgelm_t s, const vgelm_t a, const vgelm_t b, const vgelm_t c)
{
  __m512i a0 = a[0], a1 = a[1], a2 = a[2], a3 = a[3], a4 = a[4], a5 = a[5];
  __m512i b0 = b[0], b1 = b[1], b2 = b[2], b3 = b[3], b4 = b[4], b5 = b[5];
  __m512i c0 = c[0], c1 = c[1], c2 = c[2], c3 = c[3], c4 = c[4], c5 = c[5];
  __m512i r0, r1, r2, r3, r4, r5;
  __m512i s0, s1, s2, s3, s4, s5;

  // c = < C1' | C1 | C0' | C0 | F1' | F1 | F0' | F0 >
  c0 = VALIGNR(c0, c0, 4); c1 = VALIGNR(c1, c1, 4);
  c2 = VALIGNR(c2, c2, 4); c3 = VALIGNR(c3, c3, 4);
  c4 = VALIGNR(c4, c4, 4); c5 = VALIGNR(c5, c5, 4);

  // r = < E1' | E1 | E0' | E0 | A1' | A1 | A0' | A0 >
  r0 = VMBLEND(0x0F, b0, a0); r1 = VMBLEND(0x0F, b1, a1); 
  r2 = VMBLEND(0x0F, b2, a2); r3 = VMBLEND(0x0F, b3, a3); 
  r4 = VMBLEND(0x0F, b4, a4); r5 = VMBLEND(0x0F, b5, a5); 

  // a = < A1' | A1 | A0' | A0 | D1' | D1 | D0' | D0 >
  a0 = VALIGNR(a0, a0, 4); a1 = VALIGNR(a1, a1, 4);
  a2 = VALIGNR(a2, a2, 4); a3 = VALIGNR(a3, a3, 4);
  a4 = VALIGNR(a4, a4, 4); a5 = VALIGNR(a5, a5, 4);

  // s = < C1' | C1 | C0' | C0 | D1' | D1 | D0' | D0 >
  s0 = VMBLEND(0x0F, c0, a0); s1 = VMBLEND(0x0F, c1, a1); 
  s2 = VMBLEND(0x0F, c2, a2); s3 = VMBLEND(0x0F, c3, a3); 
  s4 = VMBLEND(0x0F, c4, a4); s5 = VMBLEND(0x0F, c5, a5); 

  r[0] = r0; r[1] = r1; r[2] = r2; r[3] = r3; r[4] = r4; r[5] = r5; 
  s[0] = s0; s[1] = s1; s[2] = s2; s[3] = s3; s[4] = s4; s[5] = s5;
}

// INPUT
// P           : < XP1' | XP1 | XP0' | XP0 |   ZP1' |  ZP1 |  ZP0' |  ZP0 >
// OUTPUT 
// C24_A24plus : <   C1' |   C1 |   C0' |   C0 |    A1' |   A1 |   A0' |   A0 >
// coeff1_0    : < cf11' | cf11 | cf10' | cf10 |  cf01' | cf01 | cf00' | cf00 >
void get_3_isog_1x2x2x2w(const vgelm_t P, vgelm_t C24_A24plus, vgelm_t coeff1_0)
{
  vgelm_t t1, t2, t3, t4, t5;
  vgelm_t zero = {0};

  mp2_hadamard_2x2x2w(coeff1_0, P);     //    coeff1_0 = X+Z                | X-Z           [4p]
  fpadd_4x2w(coeff1_0, coeff1_0, zero);
  fp2sqr_mont_2x2x2w(t2, coeff1_0);     //          t2 = (X+Z)^2            | (X-Z)^2       [2p]
  fp2mix3_2x2x2w(t1, t3, P, t2);        //          t1 = X                  | (X-Z)^2       [2p]
                                        //          t3 = X                  | (X+Z)^2       [2p]
  fp2addsub_2x2x2w(t1, t1, t3);         //          t1 = 2X                 | -4XZ          [2p] 
  vec_blend_4x2w(t3, t1, P, 0x0F);      //          t3 = 2X                 | Z             [2p]
  fp2sqr_mont_2x2x2w(t3, t3);           //          t3 = 4X^2               | Z^2           [2p]
  vec_blend_4x2w(t1, t3, t1, 0x0F);     //          t1 = 4X^2               | -4XZ          [2p]
  fpadd_4x2w(t1, t1, t1);               //          t1 = 8X^2               | -8XZ          [2p]
  vec_alignr4_4x2w(t3, t3);             //         !t3 = Z^2                | 4X^2          [2p]
  fp2mix3_2x2x2w(t1, t4, t3, t1);       //         !t1 = Z^2                | -8XZ          [2p]
                                        //          t4 = Z^2                | 8X^2          [2p]
  fp2addsub_2x2x2w(t4, t4, t1);         //         !t4 = 2Z^2               | 8X^2+8XZ      [2p]
  vec_alignr4_4x2w(t2, t2);             //         !t2 = (X-Z)^2            | (X+Z)^2       [2p]
  vec_alignr4_4x2w(t5, t4);             //         !t5 = 8X^2+8XZ           | 2Z^2          [2p]
  vec_blend_4x2w(t5, t5, t3, 0x0F);     //         !t5 = 8X^2+8XZ           | 4X^2          [2p]
  mp2_addsub_2x2x2w(t5, t5, t2);        //          t5 = (3X+Z)^2           | 3X^2-2XZ-Z^2  [4p]
  fp2mix8_2x2x2w(t1, t2, t5, t4, t1);   //          t1 = 2Z^2               | 3X^2-2XZ-Z^2  [2p|4p] 
                                        //          t2 = -8XZ               | (3X+Z)^2      [2p|4p]
  fp2mul_mont_2x2x2w(C24_A24plus, t1, t2); // C24_A24plus = -16XZ^3         | (3X^2-2XZ-Z^2)*(3X+Z)^2 [2p]
}

// INPUT  a = < C1' | C1 | C0' | C0 | A1' | A1 | A0' | A0 >
//        b = < D1' | D1 | D0' | D0 | B1' | B1 | B0' | B0 >
// OUTPUT r = < C1' | C1 | C0' | C0 | D1' | D1 | D0' | D0 >
static void fp2mix7_2x2x2w(vgelm_t r, const vgelm_t a, const vgelm_t b)
{
  __m512i a0 = a[0], a1 = a[1], a2 = a[2], a3 = a[3], a4 = a[4], a5 = a[5];
  __m512i b0 = b[0], b1 = b[1], b2 = b[2], b3 = b[3], b4 = b[4], b5 = b[5];
  __m512i r0, r1, r2, r3, r4, r5;

  // r = < B1' | B1 | B0' | B0 | D1' | D1 | D0' | D0 >
  r0 = VALIGNR(b0, b0, 4); r1 = VALIGNR(b1, b1, 4);
  r2 = VALIGNR(b2, b2, 4); r3 = VALIGNR(b3, b3, 4);
  r4 = VALIGNR(b4, b4, 4); r5 = VALIGNR(b5, b5, 4);

  // r = < C1' | C1 | C0' | C0 | D1' | D1 | D0' | D0 >
  r0 = VMBLEND(0x0F, a0, r0); r1 = VMBLEND(0x0F, a1, r1); 
  r2 = VMBLEND(0x0F, a2, r2); r3 = VMBLEND(0x0F, a3, r3); 
  r4 = VMBLEND(0x0F, a4, r4); r5 = VMBLEND(0x0F, a5, r5); 

  r[0] = r0; r[1] = r1; r[2] = r2; r[3] = r3; r[4] = r4; r[5] = r5; 
}

// INPUT        
// P      : <  XP1' |  XP1 |  XP0' |  XP0 |   ZP1' |  ZP1 |  ZP0' |   ZP0 >
// coeff0 : <    0  |    0 |    0  |    0 |  cf01' | cf01 | cf00' |  cf00 >
// coeff1 : < cf11' | cf11 | cf10' | cf10 |     0  |    0 |    0  |     0 >
// coeff2 : < cf21' | cf21 | cf20' | cf20 |     0  |    0 |    0  |     0 >
// OUTPUT 
// P      : <  XP1' |  XP1 |  XP0' |  XP0 |   ZP1' |  ZP1 |  ZP0' |   ZP0 >
void eval_4_isog_1x2x2x2w(vgelm_t P, const vgelm_t *coeff)
{
  vgelm_t t1, t2, t3;

  mp2_hadamard_2x2x2w(t1, P);           // t1 = X+Z                   | X-Z                   [4p]
  fp2mix7_2x2x2w(t2, coeff[1], t1);     // t2 = cf1                   | X+Z                   [2p|4p]
  fp2mul_mont_2x2x2w(t3, t2, t1);       // t3 = cf1(X+Z)              | (X+Z)(X-Z)            [2p]
  fp2mix1_2x2x2w(t1, t3, t1);           // t1 = X-Z                   | (X+Z)(X-Z)            [4p|2p]
  vec_blend_4x2w(t2, coeff[2], coeff[0], 0x0F);   // t2 = cf2         | cf0                   [2p]
  fp2mul_mont_2x2x2w(t1, t2, t1);       // t1 = cf2(X-Z)              | cf0(X+Z)(X-Z)         [2p]
  fp2mix7_2x2x2w(t2, t1, t3);           // t2 = cf2(X-Z)              | cf1(X+Z)              [2p]
  mp2_hadamard_2x2x2w(t2, t2);          // t2 = cf2(X-Z)+cf1(X+Z)     | cf2(X-Z)-cf1(X+Z)     [4p]
  fp2sqr_mont_2x2x2w(t2, t2);           // t2 = [cf2(X-Z)+cf1(X+Z)]^2 | [cf2(X-Z)-cf1(X+Z)]^2 [2p]
  fp2align2_2x2x2w(t1, t1);             // t1 = cf0(X+Z)(X-Z)         | cf0(X+Z)(X-Z)         [2p]
  mp2_addsub_2x2x2w(t1, t2, t1);        // t1 = [cf2(X-Z)+cf1(X+Z)]^2+cf0(X+Z)(X-Z) |         [4p]
                                        //      [cf2(X-Z)-cf1(X+Z)]^2-cf0(X+Z)(X-Z)           [4p]
  fp2mul_mont_2x2x2w(P, t1, t2);        //  P = Xfinal                | Zfinal                [2p]
}

// INPUT        
// Q       : <  XQ1' |  XQ1 |  XQ0' |  XQ0  |  ZQ1' |  ZQ1 |  ZQ0' |  ZQ0 >
// coeff10 : <  cf11' | cf11 | cf10' | cf10 | cf01' | cf01 | cf00' |  cf00>
// OUTPUT 
// Q       : <  XQ1' |  XQ1 |  XQ0' |  XQ0  |  ZQ1' |  ZQ1 |  ZQ0' |  ZQ0 >
void eval_3_isog_1x2x2x2w(vgelm_t Q, vgelm_t coeff0_1)
{
  vgelm_t t1;

  mp2_hadamard_2x2x2w(t1, Q);           // t1 = X+Z                   | X-Z                   [4p]
  fp2mul_mont_2x2x2w(t1, coeff0_1, t1); // t1 = cf0(X+Z)              | cf1(X-Z)              [2p]
  mp2_hadamard_2x2x2w(t1, t1);          // t1 = cf0(X+Z)+cf1(X-Z)     | cf0(X+Z)-cf1(X-Z)     [4p]
  fp2sqr_mont_2x2x2w(t1, t1);           // t1 = [cf0(X+Z)+cf1(X-Z)]^2 | [cf0(X+Z)-cf1(X-Z)]^2 [2p]
  fp2mul_mont_2x2x2w(Q, t1, Q);         //  Q = X[cf0(X+Z)+cf1(X-Z)]^2| Z[cf0(X+Z)-cf1(X-Z)]^2[2p]     
}

// INPUT
// P           : < XP1' | XP1 | XP0' | XP0 |   ZP1' |  ZP1 |  ZP0' |  ZP0 >
// OUTPUT 
// C24_A24plus : <   C1' |   C1 |   C0' |   C0 |    A1' |   A1 |   A0' |   A0 >
void get_2_isog_1x2x2x2w(const vgelm_t P, vgelm_t C24_A24plus)
{
  vgelm_t t1, t2;

  fp2sqr_mont_2x2x2w(t1, P);            //          t1 = X^2      | Z^2
  vec_alignr4_4x2w(t2, t1);             //          t2 = Z^2      | X^2
  fpsub_4x2w(t1, t1, t2);               //          t1 = X^2-Z^2  | Z^2-X^2
  vec_blend_4x2w(C24_A24plus, t2, t1, 0x0F); // C24_A24plus = Z^2 | Z^2-X^2 
}

void eval_2_isog_1x2x2x2w(vgelm_t P, vgelm_t Q)
{
  vgelm_t t1, t2;

  mp2_hadamard_2x2x2w(t1, Q);           // t1 = X2+Z2         | X2-Z2         [4p]
  mp2_hadamard_2x2x2w(t2, P);           // t2 = X+Z           | X-Z           [4p]
  vec_alignr4_4x2w(t2, t2);             // t2 = X-Z           | X+Z           [4p]
  fp2mul_mont_2x2x2w(t1, t1, t2);       // t1 = (X2+Z2)*(X-Z) | (X2-Z2)*(X+Z) [2p]
  mp2_hadamard_2x2x2w(t1, t1);          // t1 = (X2+Z2)(X-Z)+(X2-Z2)(X+Z) | (X2+Z2)(X-Z)-(X2-Z2)(X+Z)        [4p]
  fp2mul_mont_2x2x2w(P, t1, P);         //  P =  X[(X2+Z2)(X-Z)+(X2-Z2)(X+Z)] | Z[(X2+Z2)(X-Z)-(X2-Z2)(X+Z)] [2p]
}

// -----------------------------------------------------------------------------
// (4x2x1x1)-way curve/isogeny arithmetic (based on (8x1x1)-way Fp2 arithmetic).
// NOTE: (4x2x1x1)-way means each function performs 4 curve/isogeny operations, 
// where each curve/isogeny operation performs 2 Fp2 operations, and each Fp2 
// operation performs 1 Fp operation and each Fp operation uses 1 lane. 

// INPUT  a = < C | A | C | A | C | A | C | A >
//        b = < D | B | D | B | D | B | D | B >
// OUTPUT r = < C | D | C | D | C | D | C | D >
static void fpmix1_8x1w(vfelm_t r, const vfelm_t a, const vfelm_t b)
{
  __m512i a0 = a[0], a1 = a[1], a2  = a[2],  a3  = a[3];
  __m512i a4 = a[4], a5 = a[5], a6  = a[6],  a7  = a[7];
  __m512i a8 = a[8], a9 = a[9], a10 = a[10], a11 = a[11];
  __m512i b0 = b[0], b1 = b[1], b2  = b[2],  b3  = b[3]; 
  __m512i b4 = b[4], b5 = b[5], b6  = b[6],  b7  = b[7];
  __m512i b8 = b[8], b9 = b[9], b10 = b[10], b11 = b[11];
  __m512i r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11;

  // b = < B | D | B | D | B | D | B | D > 
  b0 = VSHUF(b0, 0x4E); b1 = VSHUF(b1, 0x4E); 
  b2 = VSHUF(b2, 0x4E); b3 = VSHUF(b3, 0x4E); 
  b4 = VSHUF(b4, 0x4E); b5 = VSHUF(b5, 0x4E); 
  b6 = VSHUF(b6, 0x4E); b7 = VSHUF(b7, 0x4E); 
  b8 = VSHUF(b8, 0x4E); b9 = VSHUF(b9, 0x4E);
  b10 = VSHUF(b10, 0x4E); b11 = VSHUF(b11, 0x4E); 

  // r = < C | D | C | D | C | D | C | D >
  r0 = VMBLEND(0x55, a0, b0); r1 = VMBLEND(0x55, a1, b1);
  r2 = VMBLEND(0x55, a2, b2); r3 = VMBLEND(0x55, a3, b3);
  r4 = VMBLEND(0x55, a4, b4); r5 = VMBLEND(0x55, a5, b5);
  r6 = VMBLEND(0x55, a6, b6); r7 = VMBLEND(0x55, a7, b7);
  r8 = VMBLEND(0x55, a8, b8); r9 = VMBLEND(0x55, a9, b9);
  r10 = VMBLEND(0x55, a10, b10); r11 = VMBLEND(0x55, a11, b11);

  r[0] = r0; r[1] = r1; r[2]  = r2;  r[3]  = r3; 
  r[4] = r4; r[5] = r5; r[6]  = r6;  r[7]  = r7; 
  r[8] = r8; r[9] = r9; r[10] = r10; r[11] = r11;
}

static void fp2mix1_8x1x1w(vf2elm_t r, const vf2elm_t a, const vf2elm_t b)
{
  fpmix1_8x1w(r[0], a[0], b[0]);
  fpmix1_8x1w(r[1], a[1], b[1]);
}

// INPUT  a = < C | A | C | A | C | A | C | A >
//        b = < D | B | D | B | D | B | D | B >
// OUTPUT r = < B | A | B | A | B | A | B | A >
static void fpmix2_8x1w(vfelm_t r, const vfelm_t a, const vfelm_t b)
{
  __m512i a0 = a[0], a1 = a[1], a2  = a[2],  a3  = a[3];
  __m512i a4 = a[4], a5 = a[5], a6  = a[6],  a7  = a[7];
  __m512i a8 = a[8], a9 = a[9], a10 = a[10], a11 = a[11];
  __m512i b0 = b[0], b1 = b[1], b2  = b[2],  b3  = b[3]; 
  __m512i b4 = b[4], b5 = b[5], b6  = b[6],  b7  = b[7];
  __m512i b8 = b[8], b9 = b[9], b10 = b[10], b11 = b[11];
  __m512i r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11;

  // b = < B | D | B | D | B | D | B | D > 
  b0 = VSHUF(b0, 0x4E); b1 = VSHUF(b1, 0x4E); 
  b2 = VSHUF(b2, 0x4E); b3 = VSHUF(b3, 0x4E); 
  b4 = VSHUF(b4, 0x4E); b5 = VSHUF(b5, 0x4E); 
  b6 = VSHUF(b6, 0x4E); b7 = VSHUF(b7, 0x4E); 
  b8 = VSHUF(b8, 0x4E); b9 = VSHUF(b9, 0x4E);
  b10 = VSHUF(b10, 0x4E); b11 = VSHUF(b11, 0x4E); 

  // r = < B | A | B | A | B | A | B | A >
  r0 = VMBLEND(0xAA, a0, b0); r1 = VMBLEND(0xAA, a1, b1);
  r2 = VMBLEND(0xAA, a2, b2); r3 = VMBLEND(0xAA, a3, b3);
  r4 = VMBLEND(0xAA, a4, b4); r5 = VMBLEND(0xAA, a5, b5);
  r6 = VMBLEND(0xAA, a6, b6); r7 = VMBLEND(0xAA, a7, b7);
  r8 = VMBLEND(0xAA, a8, b8); r9 = VMBLEND(0xAA, a9, b9);
  r10 = VMBLEND(0xAA, a10, b10); r11 = VMBLEND(0xAA, a11, b11);

  r[0] = r0; r[1] = r1; r[2]  = r2;  r[3]  = r3; 
  r[4] = r4; r[5] = r5; r[6]  = r6;  r[7]  = r7; 
  r[8] = r8; r[9] = r9; r[10] = r10; r[11] = r11;
}

static void fp2mix2_8x1x1w(vf2elm_t r, const vf2elm_t a, const vf2elm_t b)
{
  fpmix2_8x1w(r[0], a[0], b[0]);
  fpmix2_8x1w(r[1], a[1], b[1]);
}

// INPUT  a = < 0 | C1 | 0 | C0 | 0 | A1 | 0 | A0 >
//        b = < 0 | D1 | 0 | D0 | 0 | B1 | 0 | B0 >
// OUTPUT r = < 0 | C1 | 0 | C0 | 0 | B1 | 0 | B0 >
static void vec_blend_8x1w(vfelm_t r, const vfelm_t a, const vfelm_t b, const __mmask8 mask)
{
  __m512i a0 = a[0], a1 = a[1], a2  = a[2],  a3  = a[3];
  __m512i a4 = a[4], a5 = a[5], a6  = a[6],  a7  = a[7];
  __m512i a8 = a[8], a9 = a[9], a10 = a[10], a11 = a[11];
  __m512i b0 = b[0], b1 = b[1], b2  = b[2],  b3  = b[3]; 
  __m512i b4 = b[4], b5 = b[5], b6  = b[6],  b7  = b[7];
  __m512i b8 = b[8], b9 = b[9], b10 = b[10], b11 = b[11];
  __m512i r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11;

  r0 = VMBLEND(mask, a0, b0); r1 = VMBLEND(mask, a1, b1); 
  r2 = VMBLEND(mask, a2, b2); r3 = VMBLEND(mask, a3, b3); 
  r4 = VMBLEND(mask, a4, b4); r5 = VMBLEND(mask, a5, b5); 
  r6 = VMBLEND(mask, a6, b6); r7 = VMBLEND(mask, a7, b7); 
  r8 = VMBLEND(mask, a8, b8); r9 = VMBLEND(mask, a9, b9); 
  r10 = VMBLEND(mask, a10, b10); r11 = VMBLEND(mask, a11, b11); 

  r[0] = r0; r[1] = r1; r[2]  = r2;  r[3]  = r3; 
  r[4] = r4; r[5] = r5; r[6]  = r6;  r[7]  = r7; 
  r[8] = r8; r[9] = r9; r[10] = r10; r[11] = r11;
}

static void fp2blend_8x1x1w(vf2elm_t r, const vf2elm_t a, const vf2elm_t b, const __mmask8 mask)
{
  vec_blend_8x1w(r[0], a[0], b[0], mask);
  vec_blend_8x1w(r[1], a[1], b[1], mask);
}

// INPUT  a = < B | A | B | A | B | A | B | A >
// OUTPUT r = < A | A | A | A | A | A | A | A >
static void fpalign_8x1w(vfelm_t r, const vfelm_t a)
{
  __m512i a0 = a[0], a1 = a[1], a2  = a[2],  a3  = a[3];
  __m512i a4 = a[4], a5 = a[5], a6  = a[6],  a7  = a[7];
  __m512i a8 = a[8], a9 = a[9], a10 = a[10], a11 = a[11];
  __m512i r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11;

  // r = < A | A | A | A | A | A | A | A >
  r0 = VSHUF(a0, 0x44); r1 = VSHUF(a1, 0x44);
  r2 = VSHUF(a2, 0x44); r3 = VSHUF(a3, 0x44);
  r4 = VSHUF(a4, 0x44); r5 = VSHUF(a5, 0x44);
  r6 = VSHUF(a6, 0x44); r7 = VSHUF(a7, 0x44);
  r8 = VSHUF(a8, 0x44); r9 = VSHUF(a9, 0x44);
  r10 = VSHUF(a10, 0x44); r11 = VSHUF(a11, 0x44);

  r[0] = r0; r[1] = r1; r[2]  = r2;  r[3]  = r3; 
  r[4] = r4; r[5] = r5; r[6]  = r6;  r[7]  = r7; 
  r[8] = r8; r[9] = r9; r[10] = r10; r[11] = r11;
}

static void fp2align_8x1x1w(vf2elm_t r, const vf2elm_t a)
{
  fpalign_8x1w(r[0], a[0]);
  fpalign_8x1w(r[1], a[1]);
}

// INPUT  a = < C | A | C | A | C | A | C | A >                   all in [0, 2p)
//        b = < D | B | D | B | D | B | D | B >                   all in [0, 2p)
// OUTPUT r = < C+D | A-B | C+D | A-B | C+D | A-B | C+D | A-B >   all in [0, 4p)
static void mp_addsub_8x1w(vfelm_t r, const vfelm_t a, const vfelm_t b)
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

  // r = < D | 2p | D | 2p | D | 2p | D | 2p >
  r0 = VMBLEND(0x55, b0, vp0); r1 = VMBLEND(0x55, b1, vp1);
  r2 = VMBLEND(0x55, b2, vp2); r3 = VMBLEND(0x55, b3, vp3);
  r4 = VMBLEND(0x55, b4, vp4); r5 = VMBLEND(0x55, b5, vp5);
  r6 = VMBLEND(0x55, b6, vp6); r7 = VMBLEND(0x55, b7, vp7);
  r8 = VMBLEND(0x55, b8, vp8); r9 = VMBLEND(0x55, b9, vp9);
  r10 = VMBLEND(0x55, b10, vp10); r11 = VMBLEND(0x55, b11, vp11);

  // r = < C+D | A+2p | C+D | A+2p | C+D | A+2p | C+D | A+2p >
  r0 = VADD(r0, a0); r1 = VADD(r1, a1);
  r2 = VADD(r2, a2); r3 = VADD(r3, a3);
  r4 = VADD(r4, a4); r5 = VADD(r5, a5);
  r6 = VADD(r6, a6); r7 = VADD(r7, a7);
  r8 = VADD(r8, a8); r9 = VADD(r9, a9);
  r10 = VADD(r10, a10); r11 = VADD(r11, a11);

  // r = < C+D | A+2p-B | C+D | A+2p-B | C+D | A+2p-B | C+D | A+2p-B > 
  r0 = VMSUB(r0, 0x55, r0, b0); r1 = VMSUB(r1, 0x55, r1, b1);
  r2 = VMSUB(r2, 0x55, r2, b2); r3 = VMSUB(r3, 0x55, r3, b3);
  r4 = VMSUB(r4, 0x55, r4, b4); r5 = VMSUB(r5, 0x55, r5, b5);
  r6 = VMSUB(r6, 0x55, r6, b6); r7 = VMSUB(r7, 0x55, r7, b7);
  r8 = VMSUB(r8, 0x55, r8, b8); r9 = VMSUB(r9, 0x55, r9, b9);
  r10 = VMSUB(r10, 0x55, r10, b10); r11 = VMSUB(r11, 0x55, r11, b11);

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

static void mp2_addsub_8x1x1w(vf2elm_t r, const vf2elm_t a, const vf2elm_t b)
{
  mp_addsub_8x1w(r[0], a[0], b[0]);
  mp_addsub_8x1w(r[1], a[1], b[1]);
}

// INPUT       
// P0      : < XP0 | ZP0 | XP0 | ZP0 | XP0 | ZP0 | XP0 | ZP0>
// P1      : < XP1 | ZP1 | XP1 | ZP1 | XP1 | ZP1 | XP1 | ZP1>
// coeff00 : < 0 | cf00 | 0 | cf00 | 0 | cf00 | 0 | cf00 >
// coeff01 : < 0 | cf01 | 0 | cf01 | 0 | cf01 | 0 | cf01 >
// coeff10 : < cf10 | 0 | cf10 | 0 | cf10 | 0 | cf10 | 0 >
// coeff11 : < cf11 | 0 | cf11 | 0 | cf11 | 0 | cf11 | 0 >
// coeff20 : < cf20 | 0 | cf20 | 0 | cf20 | 0 | cf20 | 0 >
// coeff21 : < cf21 | 0 | cf21 | 0 | cf21 | 0 | cf21 | 0 >
// OUTPUT  
// P0      : < XP0 | ZP0 | XP0 | ZP0 | XP0 | ZP0 | XP0 | ZP0>
// P1      : < XP1 | ZP1 | XP1 | ZP1 | XP1 | ZP1 | XP1 | ZP1>
void eval_4_isog_4x2x1x1w(vf2elm_t P, vf2elm_t *coeff)
{
  vf2elm_t t1, t2, t3;

  mp2_hadamard_8x1x1w(t1, P);           // t1 = X+Z                   | X-Z                   [4p]
  fp2mix1_8x1x1w(t2, coeff[1], t1);     // t2 = cf1                   | X+Z                   [2p|4p]
  fp2mul_mont_8x1x1w(t3, t2, t1);       // t3 = cf1(X+Z)              | (X+Z)(X-Z)            [2p]
  fp2mix2_8x1x1w(t1, t3, t1);           // t1 = X-Z                   | (X+Z)(X-Z)            [4p|2p]
  fp2blend_8x1x1w(t2, coeff[2], coeff[0], 0x55);  // t2 = cf2         | cf0                   [2p]
  fp2mul_mont_8x1x1w(t1, t2, t1);       // t1 = cf2(X-Z)              | cf0(X+Z)(X-Z)         [2p]
  fp2mix1_8x1x1w(t2, t1, t3);           // t2 = cf2(X-Z)              | cf1(X+Z)              [2p]
  mp2_hadamard_8x1x1w(t2, t2);          // t2 = cf2(X-Z)+cf1(X+Z)     | cf2(X-Z)-cf1(X+Z)     [4p]
  fp2sqr_mont_8x1x1w(t2, t2);           // t2 = [cf2(X-Z)+cf1(X+Z)]^2 | [cf2(X-Z)-cf1(X+Z)]^2 [2p]
  fp2align_8x1x1w(t1, t1);              // t1 = cf0(X+Z)(X-Z)         | cf0(X+Z)(X-Z)         [2p]
  mp2_addsub_8x1x1w(t1, t2, t1);        // t1 = [cf2(X-Z)+cf1(X+Z)]^2+cf0(X+Z)(X-Z) |         [4p]
                                        //      [cf2(X-Z)-cf1(X+Z)]^2-cf0(X+Z)(X-Z)           [4p]
  fp2mul_mont_8x1x1w(P, t1, t2);        //  P = Xfinal                | Zfinal                [2p]
}

// INPUT       
// Q0       : < XQ0 | ZQ0 | XQ0 | ZQ0 | XQ0 | ZQ0 | XQ0 | ZP0>
// Q1       : < XQ1 | ZQ1 | XQ1 | ZQ1 | XQ1 | ZQ1 | XQ1 | ZP1>
// coeff010 : < cf00 | cf10 | cf00 | cf10 | cf00 | cf10 | cf00 | cf10 >
// coeff011 : < cf01 | cf11 | cf01 | cf11 | cf01 | cf11 | cf01 | cf11 >
// OUTPUT  
// Q0       : < XQ0 | ZQ0 | XQ0 | ZQ0 | XP0 | ZP0 | XP0 | ZP0>
// Q1       : < XQ1 | ZQ1 | XQ1 | ZQ1 | XP1 | ZP1 | XP1 | ZP1>
void eval_3_isog_4x2x1x1w(vf2elm_t Q, const vf2elm_t coeff0_1)
{
  vf2elm_t t1;

  mp2_hadamard_8x1x1w(t1, Q);           // t1 = X+Z                   | X-Z                   [4p]
  fp2mul_mont_8x1x1w(t1, coeff0_1, t1); // t1 = cf0(X+Z)              | cf1(X-Z)              [2p]
  mp2_hadamard_8x1x1w(t1, t1);          // t1 = cf0(X+Z)+cf1(X-Z)     | cf0(X+Z)-cf1(X-Z)     [4p]
  fp2sqr_mont_8x1x1w(t1, t1);           // t1 = [cf0(X+Z)+cf1(X-Z)]^2 | [cf0(X+Z)-cf1(X-Z)]   [2p]
  fp2mul_mont_8x1x1w(Q, t1, Q);         //  Q = X[cf0(X+Z)+cf1(X-Z)]^2| Z[cf0(X+Z)-cf1(X-Z)]^2[2p]            
}

// -----------------------------------------------------------------------------
// (2x2x2x1)-way curve/isogeny arithmetic (based on (4x2x1)-way Fp2 arithmetic).
// NOTE: (2x2x2x1)-way means each function performs 2 curve/isogeny operations, 
// where each curve/isogeny operation performs 2 Fp2 operations, and each Fp2 
// operation performs 2 Fp operations and each Fp operation uses 1 lane. 

// INPUT  a = < C1 | C0 | A1 | A0 | C1 | C0 | A1 | A0 >  
//        b = < D1 | D0 | B1 | B0 | D1 | D0 | B1 | B0 >   
// OUTPUT r = < C1 | C0 | D1 | D0 | C1 | C0 | D1 | D0 >  
static void fp2mix1_4x2x1w(vfelm_t r, const vfelm_t a, const vfelm_t b)
{
  __m512i a0 = a[0], a1 = a[1], a2  = a[2],  a3  = a[3];
  __m512i a4 = a[4], a5 = a[5], a6  = a[6],  a7  = a[7];
  __m512i a8 = a[8], a9 = a[9], a10 = a[10], a11 = a[11];
  __m512i b0 = b[0], b1 = b[1], b2  = b[2],  b3  = b[3]; 
  __m512i b4 = b[4], b5 = b[5], b6  = b[6],  b7  = b[7];
  __m512i b8 = b[8], b9 = b[9], b10 = b[10], b11 = b[11];
  __m512i r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11;

  // b = < B1 | B0 | D1 | D0 | B1 | B0 | D1 | D0 > 
  b0 = VPERM(b0, 0x4E); b1 = VPERM(b1, 0x4E); 
  b2 = VPERM(b2, 0x4E); b3 = VPERM(b3, 0x4E); 
  b4 = VPERM(b4, 0x4E); b5 = VPERM(b5, 0x4E); 
  b6 = VPERM(b6, 0x4E); b7 = VPERM(b7, 0x4E); 
  b8 = VPERM(b8, 0x4E); b9 = VPERM(b9, 0x4E);
  b10 = VPERM(b10, 0x4E); b11 = VPERM(b11, 0x4E);

  // r = < C1 | C0 | D1 | D0 | C1 | C0 | D1 | D0 >  
  r0 = VMBLEND(0x33, a0, b0); r1 = VMBLEND(0x33, a1, b1);
  r2 = VMBLEND(0x33, a2, b2); r3 = VMBLEND(0x33, a3, b3);
  r4 = VMBLEND(0x33, a4, b4); r5 = VMBLEND(0x33, a5, b5);
  r6 = VMBLEND(0x33, a6, b6); r7 = VMBLEND(0x33, a7, b7);
  r8 = VMBLEND(0x33, a8, b8); r9 = VMBLEND(0x33, a9, b9);
  r10 = VMBLEND(0x33, a10, b10); r11 = VMBLEND(0x33, a11, b11);

  r[0] = r0; r[1] = r1; r[2]  = r2;  r[3]  = r3; 
  r[4] = r4; r[5] = r5; r[6]  = r6;  r[7]  = r7; 
  r[8] = r8; r[9] = r9; r[10] = r10; r[11] = r11;
}

// INPUT  a = < C1 | C0 | A1 | A0 | C1 | C0 | A1 | A0 >  
//        b = < D1 | D0 | B1 | B0 | D1 | D0 | B1 | B0 > 
// OUTPUT r = < B1 | B0 | A1 | A0 | B1 | B0 | A1 | A0 > 
static void fp2mix2_4x2x1w(vfelm_t r, const vfelm_t a, const vfelm_t b)
{
  __m512i a0 = a[0], a1 = a[1], a2  = a[2],  a3  = a[3];
  __m512i a4 = a[4], a5 = a[5], a6  = a[6],  a7  = a[7];
  __m512i a8 = a[8], a9 = a[9], a10 = a[10], a11 = a[11];
  __m512i b0 = b[0], b1 = b[1], b2  = b[2],  b3  = b[3]; 
  __m512i b4 = b[4], b5 = b[5], b6  = b[6],  b7  = b[7];
  __m512i b8 = b[8], b9 = b[9], b10 = b[10], b11 = b[11];
  __m512i r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11;

  // b = < B1 | B0 | D1 | D0 | B1 | B0 | D1 | D0 > 
  b0 = VPERM(b0, 0x4E); b1 = VPERM(b1, 0x4E); 
  b2 = VPERM(b2, 0x4E); b3 = VPERM(b3, 0x4E); 
  b4 = VPERM(b4, 0x4E); b5 = VPERM(b5, 0x4E); 
  b6 = VPERM(b6, 0x4E); b7 = VPERM(b7, 0x4E); 
  b8 = VPERM(b8, 0x4E); b9 = VPERM(b9, 0x4E);
  b10 = VPERM(b10, 0x4E); b11 = VPERM(b11, 0x4E);

  // r = < B1 | B0 | A1 | A0 | B1 | B0 | A1 | A0 > 
  r0 = VMBLEND(0xCC, a0, b0); r1 = VMBLEND(0xCC, a1, b1);
  r2 = VMBLEND(0xCC, a2, b2); r3 = VMBLEND(0xCC, a3, b3);
  r4 = VMBLEND(0xCC, a4, b4); r5 = VMBLEND(0xCC, a5, b5);
  r6 = VMBLEND(0xCC, a6, b6); r7 = VMBLEND(0xCC, a7, b7);
  r8 = VMBLEND(0xCC, a8, b8); r9 = VMBLEND(0xCC, a9, b9);
  r10 = VMBLEND(0xCC, a10, b10); r11 = VMBLEND(0xCC, a11, b11);

  r[0] = r0; r[1] = r1; r[2]  = r2;  r[3]  = r3; 
  r[4] = r4; r[5] = r5; r[6]  = r6;  r[7]  = r7; 
  r[8] = r8; r[9] = r9; r[10] = r10; r[11] = r11;
}

// INPUT  a = < B1 | B0 | A1 | A0 | B1 | B0 | A1 | A0 >
// OUTPUT r = < A1 | A0 | A1 | A0 | A1 | A0 | A1 | A0 >
static void fp2align2_4x2x1w(vfelm_t r, const vfelm_t a)
{
  __m512i a0 = a[0], a1 = a[1], a2  = a[2],  a3  = a[3];
  __m512i a4 = a[4], a5 = a[5], a6  = a[6],  a7  = a[7];
  __m512i a8 = a[8], a9 = a[9], a10 = a[10], a11 = a[11];
  __m512i r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11;

  // r = < A1 | A0 | A1 | A0 | A1 | A0 | A1 | A0 >
  r0 = VPERM(a0, 0x44); r1 = VPERM(a1, 0x44);
  r2 = VPERM(a2, 0x44); r3 = VPERM(a3, 0x44);
  r4 = VPERM(a4, 0x44); r5 = VPERM(a5, 0x44);
  r6 = VPERM(a6, 0x44); r7 = VPERM(a7, 0x44);
  r8 = VPERM(a8, 0x44); r9 = VPERM(a9, 0x44);
  r10 = VPERM(a10, 0x44); r11 = VPERM(a11, 0x44);

  r[0] = r0; r[1] = r1; r[2]  = r2;  r[3]  = r3; 
  r[4] = r4; r[5] = r5; r[6]  = r6;  r[7]  = r7; 
  r[8] = r8; r[9] = r9; r[10] = r10; r[11] = r11;
}

// INPUT  a = < C1 | C0 | A1 | A0 | C1 | C0 | A1 | A0 >                         all in [0, 2p)
//        b = < D1 | D0 | B1 | B0 | D1 | D0 | B1 | B0 >                         all in [0, 2p)
// OUTPUT r = < C1+D1 | C0+D0 | A1-B1 | A0-B0 | C1+D1 | C0+D0 | A1-B1 | A0-B0 > all in [0, 4p)
static void mp2_addsub_4x2x1w(vfelm_t r, const vfelm_t a, const vfelm_t b)
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

  // r = < D1 | D0 | 2p | 2p | D1 | D0 | 2p | 2p >
  r0 = VMBLEND(0x33, b0, vp0); r1 = VMBLEND(0x33, b1, vp1);
  r2 = VMBLEND(0x33, b2, vp2); r3 = VMBLEND(0x33, b3, vp3);
  r4 = VMBLEND(0x33, b4, vp4); r5 = VMBLEND(0x33, b5, vp5);
  r6 = VMBLEND(0x33, b6, vp6); r7 = VMBLEND(0x33, b7, vp7);
  r8 = VMBLEND(0x33, b8, vp8); r9 = VMBLEND(0x33, b9, vp9);
  r10 = VMBLEND(0x33, b10, vp10); r11 = VMBLEND(0x33, b11, vp11);

  // r = < C1+D1 | C0+D0 | A1+2p | A0+2p | C1+D1 | C0+D0 | A1+2p | A0+2p >
  r0 = VADD(r0, a0); r1 = VADD(r1, a1);
  r2 = VADD(r2, a2); r3 = VADD(r3, a3);
  r4 = VADD(r4, a4); r5 = VADD(r5, a5);
  r6 = VADD(r6, a6); r7 = VADD(r7, a7);
  r8 = VADD(r8, a8); r9 = VADD(r9, a9);
  r10 = VADD(r10, a10); r11 = VADD(r11, a11);

  // r = < C1+D1 | C0+D0 | A1+2p-B1 | A0+2p-B0 | C1+D1 | C0+D0 | A1+2p-B1 | A0+2p-B0 > 
  r0 = VMSUB(r0, 0x33, r0, b0); r1 = VMSUB(r1, 0x33, r1, b1);
  r2 = VMSUB(r2, 0x33, r2, b2); r3 = VMSUB(r3, 0x33, r3, b3);
  r4 = VMSUB(r4, 0x33, r4, b4); r5 = VMSUB(r5, 0x33, r5, b5);
  r6 = VMSUB(r6, 0x33, r6, b6); r7 = VMSUB(r7, 0x33, r7, b7);
  r8 = VMSUB(r8, 0x33, r8, b8); r9 = VMSUB(r9, 0x33, r9, b9);
  r10 = VMSUB(r10, 0x33, r10, b10); r11 = VMSUB(r11, 0x33, r11, b11);

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

// INPUT        
// P      : <  XP1 |  XP0 |  ZP1 |  ZP0 |  XP1 |  XP0 |  ZP1 |  ZP0 >
// coeff0 : <    0 |    0 | cf01 | cf00 |    0 |    0 | cf01 | cf00 >
// coeff1 : < cf11 | cf10 |    0 |    0 | cf11 | cf10 |    0 |    0 >
// coeff2 : < cf21 | cf20 |    0 |    0 | cf21 | cf20 |    0 |    0 >
// OUTPUT 
// P      : <  XP1 |  XP0 |  ZP1 |  ZP0 |  XP1 |  XP0 |  ZP1 |  ZP0 >
void eval_4_isog_2x2x2x1w(vfelm_t P, vfelm_t *coeff)
{
  vfelm_t t1, t2, t3;

  mp2_hadamard_4x2x1w(t1, P);           // t1 = X+Z                   | X-Z                   [4p]
  fp2mix1_4x2x1w(t2, coeff[1], t1);     // t2 = cf1                   | X+Z                   [2p|4p]
  fp2mul_mont_4x2x1w(t3, t2, t1);       // t3 = cf1(X+Z)              | (X+Z)(X-Z)            [2p] 
  fp2mix2_4x2x1w(t1, t3, t1);           // t1 = X-Z                   | (X+Z)(X-Z)            [4p|2p]
  vec_blend_8x1w(t2, coeff[2], coeff[0], 0x33);  // t2 = cf2          | cf0                   [2p]
  fp2mul_mont_4x2x1w(t1, t2, t1);       // t1 = cf2(X-Z)              | cf0(X+Z)(X-Z)         [2p]  
  fp2mix1_4x2x1w(t2, t1, t3);           // t2 = cf2(X-Z)              | cf1(X+Z)              [2p]
  mp2_hadamard_4x2x1w(t2, t2);          // t2 = cf2(X-Z)+cf1(X+Z)     | cf2(X-Z)-cf1(X+Z)     [4p]
  fp2sqr_mont_4x2x1w(t2, t2);           // t2 = [cf2(X-Z)+cf1(X+Z)]^2 | [cf2(X-Z)-cf1(X+Z)]^2 [2p]
  fp2align2_4x2x1w(t1, t1);             // t1 = cf0(X+Z)(X-Z)         | cf0(X+Z)(X-Z)         [2p]
  mp2_addsub_4x2x1w(t1, t2, t1);        // t1 = [cf2(X-Z)+cf1(X+Z)]^2+cf0(X+Z)(X-Z) |         [4p]
                                        //      [cf2(X-Z)-cf1(X+Z)]^2-cf0(X+Z)(X-Z)           [4p]
  fp2mul_mont_4x2x1w(P, t1, t2);        //  P = Xfinal                | Zfinal                [2p]
}

// INPUT        
// Q       : <  XQ1 |  XQ0 |  ZQ1 |  ZQ0 |  XQ1 |  XQ0 |  ZQ1 |  ZQ0 >
// coeff01 : < cf01 | cf00 | cf11 | cf10 | cf01 | cf00 | cf11 | cf10 >
// OUTPUT 
// Q       : <  XQ1 |  XQ0 |  ZQ1 |  ZQ0 |  XQ1 |  XQ0 |  ZQ1 |  ZQ0 >
void eval_3_isog_2x2x2x1w(vfelm_t Q, vfelm_t coeff0_1)
{
  vfelm_t t1;

  mp2_hadamard_4x2x1w(t1, Q);           // t1 = X+Z                   | X-Z                   [4p]
  fp2mul_mont_4x2x1w(t1, coeff0_1, t1); // t1 = cf0(X+Z)              | cf1(X-Z)              [2p]
  mp2_hadamard_4x2x1w(t1, t1);          // t1 = cf0(X+Z)+cf1(X-Z)     | cf0(X+Z)-cf1(X-Z)     [4p]
  fp2sqr_mont_4x2x1w(t1, t1);           // t1 = [cf0(X+Z)+cf1(X-Z)]^2 | [cf0(X+Z)-cf1(X-Z)]   [2p]
  fp2mul_mont_4x2x1w(Q, t1, Q);         //  Q = X[cf0(X+Z)+cf1(X-Z)]^2| Z[cf0(X+Z)-cf1(X-Z)]^2[2p]     
}

// INPUT  a = < E1 | E0 | A1 | A0 | E1 | E0 | A1 | A0 >
//        b = < F1 | F0 | B1 | B0 | F1 | F0 | B1 | B0 >
//        c = < G1 | G0 | C1 | C0 | G1 | G0 | C1 | C0 >
//        d = < H1 | H0 | D1 | D0 | H1 | H0 | D1 | D0 >
// OUTPUT r = < E1 | E0 | B1 | B0 | E1 | E0 | B1 | B0 >
//        s = < G1 | G0 | H1 | H0 | G1 | G0 | H1 | H0 >
static void fp2mix3_4x2x1w(vfelm_t r, vfelm_t s, const vfelm_t a, const vfelm_t b, const vfelm_t c, const vfelm_t d)
{
  __m512i a0 = a[0], a1 = a[1], a2  = a[2],  a3  = a[3];
  __m512i a4 = a[4], a5 = a[5], a6  = a[6],  a7  = a[7];
  __m512i a8 = a[8], a9 = a[9], a10 = a[10], a11 = a[11];
  __m512i b0 = b[0], b1 = b[1], b2  = b[2],  b3  = b[3]; 
  __m512i b4 = b[4], b5 = b[5], b6  = b[6],  b7  = b[7];
  __m512i b8 = b[8], b9 = b[9], b10 = b[10], b11 = b[11];
  __m512i c0 = c[0], c1 = c[1], c2  = c[2],  c3  = c[3];
  __m512i c4 = c[4], c5 = c[5], c6  = c[6],  c7  = c[7];
  __m512i c8 = c[8], c9 = c[9], c10 = c[10], c11 = c[11];
  __m512i d0 = d[0], d1 = d[1], d2  = d[2],  d3  = d[3];
  __m512i d4 = d[4], d5 = d[5], d6  = d[6],  d7  = d[7];
  __m512i d8 = d[8], d9 = d[9], d10 = d[10], d11 = d[11];
  __m512i r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11;
  __m512i s0, s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, s11;

  // r = < E1 | E0 | B1 | B0 | E1 |  E0 | B1 | B0 >
  r0 = VMBLEND(0x33, a0, b0); r1 = VMBLEND(0x33, a1, b1);
  r2 = VMBLEND(0x33, a2, b2); r3 = VMBLEND(0x33, a3, b3);
  r4 = VMBLEND(0x33, a4, b4); r5 = VMBLEND(0x33, a5, b5);
  r6 = VMBLEND(0x33, a6, b6); r7 = VMBLEND(0x33, a7, b7);
  r8 = VMBLEND(0x33, a8, b8); r9 = VMBLEND(0x33, a9, b9);
  r10 = VMBLEND(0x33, a10, b10); r11 = VMBLEND(0x33, a11, b11);

  // s = < D1 | D0 | H1 | H0 | D1 | D0 | H1 | H0 >
  s0 = VPERM(d0, 0x4E); s1 = VPERM(d1, 0x4E);
  s2 = VPERM(d2, 0x4E); s3 = VPERM(d3, 0x4E);
  s4 = VPERM(d4, 0x4E); s5 = VPERM(d5, 0x4E);
  s6 = VPERM(d6, 0x4E); s7 = VPERM(d7, 0x4E);
  s8 = VPERM(d8, 0x4E); s9 = VPERM(d9, 0x4E);
  s10 = VPERM(d10, 0x4E); s11 = VPERM(d11, 0x4E);

  // s = < G1 | G0 | H1 | H0 | G1 | G0 | H1 | H0 >
  s0 = VMBLEND(0x33, c0, s0); s1 = VMBLEND(0x33, c1, s1); 
  s2 = VMBLEND(0x33, c2, s2); s3 = VMBLEND(0x33, c3, s3); 
  s4 = VMBLEND(0x33, c4, s4); s5 = VMBLEND(0x33, c5, s5); 
  s6 = VMBLEND(0x33, c6, s6); s7 = VMBLEND(0x33, c7, s7); 
  s8 = VMBLEND(0x33, c8, s8); s9 = VMBLEND(0x33, c9, s9);
  s10 = VMBLEND(0x33, c10, s10); s11 = VMBLEND(0x33, c11, s11); 

  r[0] = r0; r[1] = r1; r[2]  = r2;  r[3]  = r3; 
  r[4] = r4; r[5] = r5; r[6]  = r6;  r[7]  = r7; 
  r[8] = r8; r[9] = r9; r[10] = r10; r[11] = r11;

  s[0] = s0; s[1] = s1; s[2]  = s2;  s[3]  = s3; 
  s[4] = s4; s[5] = s5; s[6]  = s6;  s[7]  = s7; 
  s[8] = s8; s[9] = s9; s[10] = s10; s[11] = s11;
}

// INPUT        
// P          : < XP1 |  XP0 |   ZP1 |  ZP0 |  XP1 |  XP0 |  ZP1 |  ZP0 >
// C24_A24plus: <  C1 |   C0 |    A1 |   A0 |   C1 |   C0 |   A1 |   A0 >
// OUTPUT 
// Q          : <  XQ1 |  XQ0 |  ZQ1 |  ZQ0 |  XQ1 |  XQ0 |  ZQ1 |  ZQ0 >
void xDBL_2x2x2x1w(const vfelm_t P, vfelm_t Q, const vfelm_t C24_A24plus)
{
  vfelm_t t1, t2, t3, t4;

  mp2_hadamard_4x2x1w(t1, P);             // t1 = X+Z               | X-Z                   [4p]
  fp2sqr_mont_4x2x1w(t2, t1);             // t2 = (X+Z)^2           | (X-Z)^2               [2p]
  // fp2sqr_mont_4x2x1w(Q, t1);
  mp2_hadamard_4x2x1w(t1, t2);            // t1 = ((X+Z)^2+(X-Z)^2) | 4*X*Z                 [4p]
  // mp2_hadamard_4x2x1w(Q, t2);
  fp2mix2_4x2x1w(t3, t1, t2);             // t3 = (X-Z)^2           | 4XZ                   [2p|4p]
  fp2mul_mont_4x2x1w(t3, C24_A24plus, t3);// t3 = C(X-Z)^2          | 4AXZ                  [2p]
  mp2_hadamard_4x2x1w(t4, t3);            // t4 = C(X-Z)^2+4AXZ     | C(X-Z)^2-4AXZ         [4p]
  fp2mix3_4x2x1w(t1, t2, t2, t1, t3, t4); // t1 = (X+Z)^2           | 4XZ                   [2p|4p]
                                          // t2 = C(X-Z)^2          | C(X-Z)^2+4AXZ         [2p|4p]
  fp2mul_mont_4x2x1w(Q, t1, t2);          // Q  = C(X-Z)^2(X+Z)^2   | (C(X-Z)^2+4AXZ)*4XZ   [2p]
}

// INPUT  a = < D1 | D0 | C1 | C0 | B1 | B0 | A1 | A0 >
// OUTPUT r = < C1 | C0 | D1 | D0 | A1 | A0 | B1 | B0 >
static void fp2align3_4x2x1w(vfelm_t r, const vfelm_t a)
{
  __m512i a0 = a[0], a1 = a[1], a2  = a[2],  a3  = a[3];
  __m512i a4 = a[4], a5 = a[5], a6  = a[6],  a7  = a[7];
  __m512i a8 = a[8], a9 = a[9], a10 = a[10], a11 = a[11];
  __m512i r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11;

  // r = C1 | C0 | D1 | D0 | A1 | A0 | B1 | B0
  r0 = VPERM(a0, 0x4E); r1 = VPERM(a1, 0x4E);
  r2 = VPERM(a2, 0x4E); r3 = VPERM(a3, 0x4E);
  r4 = VPERM(a4, 0x4E); r5 = VPERM(a5, 0x4E);
  r6 = VPERM(a6, 0x4E); r7 = VPERM(a7, 0x4E);
  r8 = VPERM(a8, 0x4E); r9 = VPERM(a9, 0x4E);
  r10 = VPERM(a10, 0x4E); r11 = VPERM(a11, 0x4E);

  r[0] = r0; r[1] = r1; r[2]  = r2;  r[3]  = r3; 
  r[4] = r4; r[5] = r5; r[6]  = r6;  r[7]  = r7; 
  r[8] = r8; r[9] = r9; r[10] = r10; r[11] = r11;
}

// INPUT
// P           : < XP1 | XP0 | ZP1 | ZP0 | XP1 | XP0 | ZP1 | ZP0 >
// OUTPUT 
// C24_A24plus : <   C1 |   C0 |   A1 |   A0 |   C1 |   C0 |   A1 |   A0 >
// coeff0      : <    % |    % | cf01 | cf00 |    % |    % | cf01 | cf00 >
// coeff2_1    : < cf21 | cf20 | cf11 | cf10 | cf21 | cf20 | cf11 | cf10 >
void get_4_isog_2x2x2x1w(const vfelm_t P, vfelm_t C24_A24plus, vfelm_t coeff__0, vfelm_t coeff2_1)
{
  vfelm_t t1, t2;

  mp2_hadamard_4x2x1w(coeff2_1, P);     //    coeff2_1 = X+Z  | X-Z
  fp2sqr_mont_4x2x1w(t1, P);            //          t1 = X^2  | Z^2
  mp_add_8x1w(t1, t1, t1);              //          t1 = 2X^2 | 2Z^2
  fp2sqr_mont_4x2x1w(t2, t1);           //          t2 = 4X^4 | 4Z^4
  fp2align3_4x2x1w(C24_A24plus, t2);    // C24_A24plus = 4Z^4 | 4X^4
  mp_add_8x1w(coeff__0, t1, t1);        //      coeff0 =    % | 4Z^2 
}

void get_2_isog_2x2x2x1w(const vfelm_t P, vfelm_t C24_A24plus)
{
  vfelm_t t1, t2;

  fp2sqr_mont_4x2x1w(t1, P);            //          t1 = X^2      | Z^2
  fp2align3_4x2x1w(t2, t1);             //          t2 = Z^2      | X^2
  fpsub_8x1w(t1, t1, t2);               //          t1 = X^2-Z^2  | Z^2-X^2
  vec_blend_8x1w(C24_A24plus, t2, t1, 0x33); // C24_A24plus = Z^2 | Z^2-X^2
}

void eval_2_isog_2x2x2x1w(vfelm_t P, vfelm_t Q)
{
  vfelm_t t1, t2;

  mp2_hadamard_4x2x1w(t1, Q);           // t1 = X2+Z2         | X2-Z2         [4p]
  mp2_hadamard_4x2x1w(t2, P);           // t2 = X+Z           | X-Z           [4p]
  fp2align3_4x2x1w(t2, t2);             // t2 = X-Z           | X+Z           [4p]
  fp2mul_mont_4x2x1w(t1, t1, t2);       // t1 = (X2+Z2)*(X-Z) | (X2-Z2)*(X+Z) [2p]
  mp2_hadamard_4x2x1w(t1, t1);          // t1 = (X2+Z2)(X-Z)+(X2-Z2)(X+Z) | (X2+Z2)(X-Z)-(X2-Z2)(X+Z)        [4p]
  fp2mul_mont_4x2x1w(P, t1, P);         //  P =  X[(X2+Z2)(X-Z)+(X2-Z2)(X+Z)] | Z[(X2+Z2)(X-Z)-(X2-Z2)(X+Z)] [2p]
}


// -----------------------------------------------------------------------------
// 1-way or mixed vectorized curve/isogeny arithmetic. 

void pointcopy_1w(point_proj_r51_t Q, const point_proj_r51_t P)
{
  fp2copy_1w(Q->X, P->X);
  fp2copy_1w(Q->Z, P->Z);
}

// if swap == 0, then x3z3x2z2, z1x1_A -> x3z3x2z2, z1x1_A
// if swap == 1, then x3z3x2z2, z1x1_A -> x1z1x2z2, z3x3_A
void swap_points_1x4x2x1w(vfelm_t x3z3x2z2, vfelm_t z1x1_A ,int swap)
{
  vfelm_t x1z1_A, z3x3x2z2;
  uint8_t mask = ((uint8_t)(-swap)) & 0xF0; // mask = 0x00 or 0xF0
  const __m512i perm = VSET(5, 4, 7, 6, 3, 2, 1, 0);
  int i;

  for (i = 0; i < VNWORDS; i++) {
    z3x3x2z2[i] = VPERMV(perm, x3z3x2z2[i]);
    x1z1_A[i]   = VPERMV(perm, z1x1_A[i]);
  }

  for (i = 0; i < VNWORDS; i++) {
    x3z3x2z2[i] = VMBLEND(mask, x3z3x2z2[i], x1z1_A[i]);
    z1x1_A[i]   = VMBLEND(mask, z1x1_A[i], z3x3x2z2[i]);
  }
}

void LADDER3PT_1x4x2x1w(const f2elm_r51_t xP, const f2elm_r51_t xQ, const f2elm_r51_t xPQ, const digit_t* m, const unsigned int AliceOrBob, point_proj_r51_t R, const f2elm_r51_t A)
{
  vfelm_t x3z3x2z2, z1x1_A;
  int i, nbits, bit, swap, prevbit = 0;

  if (AliceOrBob == ALICE) nbits = OALICE_BITS;
  else                     nbits = OBOB_BITS - 1;

  // Initializing AVX-512 vectors 
  // x3z3x2z2: < x31 | x30 | z31 | z30 | x21 | x20 | z21 | z20 >
  // z1x1_A  : < z11 | z10 | x11 | x10 |  0  |  0  |  A1 | A0  >
  for (i = 0; i < VNWORDS; i++) {
    x3z3x2z2[i] = _SET(xPQ[1][i], xPQ[0][i], 0, vmont_R[i], xQ[1][i], xQ[0][i], 0, vmont_R[i]); 
    z1x1_A[i]   = _SET(0, vmont_R[i], xP[1][i], xP[0][i], 0, 0, A[1][i], A[0][i]);
  }

  // Main loop
  for (i = 0; i < nbits; i++) {
    bit = (m[i >> LOG2RADIX] >> (i & (RADIX-1))) & 1;
    swap = bit ^ prevbit;
    prevbit = bit;

    swap_points_1x4x2x1w(x3z3x2z2, z1x1_A, swap);
    xDBLADD_1x4x2x1w(x3z3x2z2, z1x1_A);
  }
  swap_points_1x4x2x1w(x3z3x2z2, z1x1_A, prevbit);

  // extract the point R
  get_channel_8x1w(R->X[0], z1x1_A, 4);
  get_channel_8x1w(R->X[1], z1x1_A, 5);
  get_channel_8x1w(R->Z[0], z1x1_A, 6);
  get_channel_8x1w(R->Z[1], z1x1_A, 7);
}

// if swap == 0, then x3z3x2z2, z1x1_A -> x3z3x2z2, z1x1_A
// if swap == 1, then x3z3x2z2, z1x1_A -> x1z1x2z2, z3x3_A
void swap_points_2x4x1x1w(vf2elm_t x3z3x2z2, vf2elm_t z1x1_A ,int swap)
{
  vf2elm_t x1z1_A, z3x3x2z2;
  uint8_t mask = ((uint8_t)(-swap)) & 0xCC; // mask = 0x00 or 0xCC
  int i;

  for (i = 0; i < VNWORDS; i++) {
    z3x3x2z2[0][i] = VPERM(x3z3x2z2[0][i], 0xB4);
    z3x3x2z2[1][i] = VPERM(x3z3x2z2[1][i], 0xB4);
    x1z1_A[0][i]   = VPERM(z1x1_A[0][i], 0xB4);
    x1z1_A[1][i]   = VPERM(z1x1_A[1][i], 0xB4);
  }

  for (i = 0; i < VNWORDS; i++) {
    x3z3x2z2[0][i] = VMBLEND(mask, x3z3x2z2[0][i], x1z1_A[0][i]);
    x3z3x2z2[1][i] = VMBLEND(mask, x3z3x2z2[1][i], x1z1_A[1][i]);
    z1x1_A[0][i]   = VMBLEND(mask, z1x1_A[0][i], z3x3x2z2[0][i]);
    z1x1_A[1][i]   = VMBLEND(mask, z1x1_A[1][i], z3x3x2z2[1][i]);
  }
}

void LADDER3PT_2x4x1x1w(const f2elm_r51_t xP, const f2elm_r51_t xQ, const f2elm_r51_t xPQ, \
                        const f2elm_r51_t _xP, const f2elm_r51_t _xQ, const f2elm_r51_t _xPQ, \
                        const digit_t* m, const unsigned int AliceOrBob, \
                        point_proj_r51_t R, const f2elm_r51_t A, \
                        point_proj_r51_t _R, const f2elm_r51_t _A)
{
  vf2elm_t x3z3x2z2, z1x1_A;
  int i, nbits, bit, swap, prevbit = 0;

  if (AliceOrBob == ALICE) nbits = OALICE_BITS;
  else                     nbits = OBOB_BITS - 1;

  // Initializing AVX-512 vectors 
  // x3z3x2z2: < x3 | z3 | x2 | z2 | _x3 | _z3 | _x2 | _z2 >
  // z1x1_A  : < z1 | x1 | 0  | A  | _z1 | _x1 |  0  | _A  >
  for (i = 0; i < VNWORDS; i++) {
    x3z3x2z2[0][i] = _SET( xPQ[0][i], vmont_R[i],  xQ[0][i], vmont_R[i], \
                          _xPQ[0][i], vmont_R[i], _xQ[0][i], vmont_R[i]); 
    x3z3x2z2[1][i] = _SET( xPQ[1][i], 0,  xQ[1][i], 0, \
                          _xPQ[1][i], 0, _xQ[1][i], 0);
    z1x1_A[0][i]   = _SET(vmont_R[i],  xP[0][i], 0,  A[0][i], \
                          vmont_R[i], _xP[0][i], 0, _A[0][i]);
    z1x1_A[1][i]   = _SET(0,  xP[1][i], 0,  A[1][i], \
                          0, _xP[1][i], 0, _A[1][i]);
  }

  // Main loop
  for (i = 0; i < nbits; i++) {
    bit = (m[i >> LOG2RADIX] >> (i & (RADIX-1))) & 1;
    swap = bit ^ prevbit;
    prevbit = bit;

    swap_points_2x4x1x1w(x3z3x2z2, z1x1_A, swap);
    xDBLADD_2x4x1x1w(x3z3x2z2, z1x1_A);
  }
  swap_points_2x4x1x1w(x3z3x2z2, z1x1_A, prevbit);

  // extract the point R
  get_channel_8x1w(R->X[0], z1x1_A[0], 6);
  get_channel_8x1w(R->X[1], z1x1_A[1], 6);
  get_channel_8x1w(R->Z[0], z1x1_A[0], 7);
  get_channel_8x1w(R->Z[1], z1x1_A[1], 7);
  get_channel_8x1w(_R->X[0], z1x1_A[0], 2);
  get_channel_8x1w(_R->X[1], z1x1_A[1], 2);
  get_channel_8x1w(_R->Z[0], z1x1_A[0], 3);
  get_channel_8x1w(_R->Z[1], z1x1_A[1], 3);
}

void carryp_1w(uint64_t *a)
{
  int i;

  for (i = 0; i < VNWORDS-1; i++) {
    a[i+1] += a[i]>>VBRADIX;
    a[i] &= VBMASK;
  }
}

void xDBLe_1x2x2x2w(vgelm_t vP, point_proj_r51_t Q, const vgelm_t C24_A24plus, const int e)
{
  int i;

  for (i = 0; i < e; i++) xDBL_1x2x2x2w(vP, vP, C24_A24plus);

  // extract the point Q from vP
  get_channel_4x2w(Q->X[0], vP, 4); get_channel_4x2w(Q->X[1], vP, 6);
  get_channel_4x2w(Q->Z[0], vP, 0); get_channel_4x2w(Q->Z[1], vP, 2);

  // final carry propagation making Q strictly in radix-2^51
  carryp_1w(Q->X[0]); carryp_1w(Q->X[1]);
  carryp_1w(Q->Z[0]); carryp_1w(Q->Z[1]);
}

void xDBLe_2x2x2x1w(vfelm_t vP, point_proj_r51_t Q, point_proj_r51_t _Q, const vfelm_t C24_A24plus, const int e)
{
  int i;

  for (i = 0; i < e; i++) xDBL_2x2x2x1w(vP, vP, C24_A24plus);

  // extract the points Q and _Q from vP
  get_channel_8x1w(Q->X[0], vP, 6); get_channel_8x1w(Q->X[1], vP, 7);
  get_channel_8x1w(Q->Z[0], vP, 4); get_channel_8x1w(Q->Z[1], vP, 5);
  get_channel_8x1w(_Q->X[0], vP, 2); get_channel_8x1w(_Q->X[1], vP, 3);
  get_channel_8x1w(_Q->Z[0], vP, 0); get_channel_8x1w(_Q->Z[1], vP, 1);
}

void eval_4_isog_parallel_kg(point_proj_r51_t *pts, point_proj_r51_t phiP, point_proj_r51_t phiQ, point_proj_r51_t phiR, const vgelm_t coeff__0, const vgelm_t coeff2_1, const int num)
{
  felm_r51_t cf00, cf01, cf10, cf11, cf20, cf21;
  int i;

  // extract constants
  get_channel_4x2w(cf00, coeff__0, 0); get_channel_4x2w(cf01, coeff__0, 2);
  get_channel_4x2w(cf10, coeff2_1, 0); get_channel_4x2w(cf11, coeff2_1, 2);
  get_channel_4x2w(cf20, coeff2_1, 4); get_channel_4x2w(cf21, coeff2_1, 6);

  // carry propagation
  carryp_1w(cf00); carryp_1w(cf01);
  carryp_1w(cf10); carryp_1w(cf11);
  carryp_1w(cf20); carryp_1w(cf21);

  // depends on public info
  if (num == 3 | num == 4 | num == 5 | num == 6) {
    // 4x2x1x1w < pts[0] | phiP | phiQ | phiR >
    vf2elm_t P, coeff[3];

    // initialize vectors 
    for (i = 0; i < VNWORDS; i++) {
      P[0][i] = _SET(pts[0]->X[0][i], pts[0]->Z[0][i], phiP->X[0][i], phiP->Z[0][i], \
                       phiQ->X[0][i],   phiQ->Z[0][i], phiR->X[0][i], phiR->Z[0][i]);
      P[1][i] = _SET(pts[0]->X[1][i], pts[0]->Z[1][i], phiP->X[1][i], phiP->Z[1][i], \
                       phiQ->X[1][i],   phiQ->Z[1][i], phiR->X[1][i], phiR->Z[1][i]);
      coeff[0][0][i] = _SET(0, cf00[i], 0, cf00[i], 0, cf00[i], 0, cf00[i]);
      coeff[0][1][i] = _SET(0, cf01[i], 0, cf01[i], 0, cf01[i], 0, cf01[i]);
      coeff[1][0][i] = _SET(cf10[i], 0, cf10[i], 0, cf10[i], 0, cf10[i], 0);
      coeff[1][1][i] = _SET(cf11[i], 0, cf11[i], 0, cf11[i], 0, cf11[i], 0);
      coeff[2][0][i] = _SET(cf20[i], 0, cf20[i], 0, cf20[i], 0, cf20[i], 0);
      coeff[2][1][i] = _SET(cf21[i], 0, cf21[i], 0, cf21[i], 0, cf21[i], 0);
    }

    eval_4_isog_4x2x1x1w(P, coeff);

    // extract results
    get_channel_8x1w(phiP->X[0], P[0], 5);   get_channel_8x1w(phiP->X[1], P[1], 5);
    get_channel_8x1w(phiP->Z[0], P[0], 4);   get_channel_8x1w(phiP->Z[1], P[1], 4);
    get_channel_8x1w(phiQ->X[0], P[0], 3);   get_channel_8x1w(phiQ->X[1], P[1], 3);
    get_channel_8x1w(phiQ->Z[0], P[0], 2);   get_channel_8x1w(phiQ->Z[1], P[1], 2);
    get_channel_8x1w(phiR->X[0], P[0], 1);   get_channel_8x1w(phiR->X[1], P[1], 1);
    get_channel_8x1w(phiR->Z[0], P[0], 0);   get_channel_8x1w(phiR->Z[1], P[1], 0);
  
    // depends on public info
    if (num >= 4) {
      get_channel_8x1w(pts[0]->X[0], P[0], 7); get_channel_8x1w(pts[0]->X[1], P[1], 7);
      get_channel_8x1w(pts[0]->Z[0], P[0], 6); get_channel_8x1w(pts[0]->Z[1], P[1], 6);
    }

    // depends on public info
    if (num == 5) {
      // 1x2x2x2w < pts[1] >
      vgelm_t _P, _coeff[3];

      // initialize vectors 
      for (i = 0; i < VGWORDS; i++) {
        _P[i] = _SET(pts[1]->X[1][i+VGWORDS], pts[1]->X[1][i], pts[1]->X[0][i+VGWORDS], pts[1]->X[0][i],\
                     pts[1]->Z[1][i+VGWORDS], pts[1]->Z[1][i], pts[1]->Z[0][i+VGWORDS], pts[1]->Z[0][i]);
        _coeff[0][i] = _SET(0, 0, 0, 0, cf01[i+VGWORDS], cf01[i], cf00[i+VGWORDS], cf00[i]); 
        _coeff[1][i] = _SET(cf11[i+VGWORDS], cf11[i], cf10[i+VGWORDS], cf10[i], 0, 0, 0, 0); 
        _coeff[2][i] = _SET(cf21[i+VGWORDS], cf21[i], cf20[i+VGWORDS], cf20[i], 0, 0, 0, 0); 
      }

      eval_4_isog_1x2x2x2w(_P, _coeff);

      // extract results
      get_channel_4x2w(pts[1]->X[0], _P, 4); get_channel_4x2w(pts[1]->X[1], _P, 6);
      get_channel_4x2w(pts[1]->Z[0], _P, 0); get_channel_4x2w(pts[1]->Z[1], _P, 2);

      // carry propagation 
      carryp_1w(pts[1]->X[0]); carryp_1w(pts[1]->X[1]); 
      carryp_1w(pts[1]->Z[0]); carryp_1w(pts[1]->Z[1]);  
    }

    // depends on public info
    if (num == 6) {
      // 2x2x2x1w < pts[1] | pts[2] >
      vfelm_t _P, _coeff[3];

      for (i = 0; i < VNWORDS; i++) {
        _P[i] = _SET(pts[1]->X[1][i], pts[1]->X[0][i], pts[1]->Z[1][i], pts[1]->Z[0][i], \
                     pts[2]->X[1][i], pts[2]->X[0][i], pts[2]->Z[1][i], pts[2]->Z[0][i]);
        _coeff[0][i] = _SET(0, 0, cf01[i], cf00[i], 0, 0, cf01[i], cf00[i]);
        _coeff[1][i] = _SET(cf11[i], cf10[i], 0, 0, cf11[i], cf10[i], 0, 0);
        _coeff[2][i] = _SET(cf21[i], cf20[i], 0, 0, cf21[i], cf20[i], 0, 0);
      }

      eval_4_isog_2x2x2x1w(_P, _coeff);

      // extract results
      get_channel_8x1w(pts[1]->X[0], _P, 6); get_channel_8x1w(pts[1]->X[1], _P, 7);
      get_channel_8x1w(pts[1]->Z[0], _P, 4); get_channel_8x1w(pts[1]->Z[1], _P, 5);
      get_channel_8x1w(pts[2]->X[0], _P, 2); get_channel_8x1w(pts[2]->X[1], _P, 3);
      get_channel_8x1w(pts[2]->Z[0], _P, 0); get_channel_8x1w(pts[2]->Z[1], _P, 1);
    }
  }
  // depends on public info
  else if (num == 7 | num == 8 | num == 9 | num == 10 | num == 11) {
    // 8x1x1x1w < pts[0] | pts[1] | pts[2] | pts[3] | pts[4] | phiP | phiQ | phiR >
    vpoint_proj_t P;
    vf2elm_t coeff[3];

    // initialize vectors 
    for (i = 0; i < VNWORDS; i++) {
      P->X[0][i] = _SET(pts[0]->X[0][i], pts[1]->X[0][i], pts[2]->X[0][i], pts[3]->X[0][i], \
                        pts[4]->X[0][i],   phiP->X[0][i],   phiQ->X[0][i],   phiR->X[0][i]);
      P->X[1][i] = _SET(pts[0]->X[1][i], pts[1]->X[1][i], pts[2]->X[1][i], pts[3]->X[1][i], \
                        pts[4]->X[1][i],   phiP->X[1][i],   phiQ->X[1][i],   phiR->X[1][i]);
      P->Z[0][i] = _SET(pts[0]->Z[0][i], pts[1]->Z[0][i], pts[2]->Z[0][i], pts[3]->Z[0][i], \
                        pts[4]->Z[0][i],   phiP->Z[0][i],   phiQ->Z[0][i],   phiR->Z[0][i]);
      P->Z[1][i] = _SET(pts[0]->Z[1][i], pts[1]->Z[1][i], pts[2]->Z[1][i], pts[3]->Z[1][i], \
                        pts[4]->Z[1][i],   phiP->Z[1][i],   phiQ->Z[1][i],   phiR->Z[1][i]);
      coeff[0][0][i] = VSET1(cf00[i]); coeff[0][1][i] = VSET1(cf01[i]);
      coeff[1][0][i] = VSET1(cf10[i]); coeff[1][1][i] = VSET1(cf11[i]);
      coeff[2][0][i] = VSET1(cf20[i]); coeff[2][1][i] = VSET1(cf21[i]);
    }

    eval_4_isog_8x1x1x1w(P, coeff);

    // extract results
    get_channel_8x1w(pts[0]->X[0], P->X[0], 7); get_channel_8x1w(pts[0]->X[1], P->X[1], 7);
    get_channel_8x1w(pts[0]->Z[0], P->Z[0], 7); get_channel_8x1w(pts[0]->Z[1], P->Z[1], 7);
    get_channel_8x1w(pts[1]->X[0], P->X[0], 6); get_channel_8x1w(pts[1]->X[1], P->X[1], 6);
    get_channel_8x1w(pts[1]->Z[0], P->Z[0], 6); get_channel_8x1w(pts[1]->Z[1], P->Z[1], 6);
    get_channel_8x1w(pts[2]->X[0], P->X[0], 5); get_channel_8x1w(pts[2]->X[1], P->X[1], 5);
    get_channel_8x1w(pts[2]->Z[0], P->Z[0], 5); get_channel_8x1w(pts[2]->Z[1], P->Z[1], 5);
    get_channel_8x1w(pts[3]->X[0], P->X[0], 4); get_channel_8x1w(pts[3]->X[1], P->X[1], 4);
    get_channel_8x1w(pts[3]->Z[0], P->Z[0], 4); get_channel_8x1w(pts[3]->Z[1], P->Z[1], 4);
    get_channel_8x1w(phiP->X[0], P->X[0], 2);   get_channel_8x1w(phiP->X[1], P->X[1], 2);
    get_channel_8x1w(phiP->Z[0], P->Z[0], 2);   get_channel_8x1w(phiP->Z[1], P->Z[1], 2);
    get_channel_8x1w(phiQ->X[0], P->X[0], 1);   get_channel_8x1w(phiQ->X[1], P->X[1], 1);
    get_channel_8x1w(phiQ->Z[0], P->Z[0], 1);   get_channel_8x1w(phiQ->Z[1], P->Z[1], 1);
    get_channel_8x1w(phiR->X[0], P->X[0], 0);   get_channel_8x1w(phiR->X[1], P->X[1], 0);
    get_channel_8x1w(phiR->Z[0], P->Z[0], 0);   get_channel_8x1w(phiR->Z[1], P->Z[1], 0);

    // depends on public info
    if (num >= 8) {
      get_channel_8x1w(pts[4]->X[0], P->X[0], 3); get_channel_8x1w(pts[4]->X[1], P->X[1], 3);
      get_channel_8x1w(pts[4]->Z[0], P->Z[0], 3); get_channel_8x1w(pts[4]->Z[1], P->Z[1], 3);
    }

    // depends on public info
    if (num == 9) {
      // 1x2x2x2w < pts[5] >
      vgelm_t _P, _coeff[3];

      // initialize vectors 
      for (i = 0; i < VGWORDS; i++) {
        _P[i] = _SET(pts[5]->X[1][i+VGWORDS], pts[5]->X[1][i], pts[5]->X[0][i+VGWORDS], pts[5]->X[0][i],\
                     pts[5]->Z[1][i+VGWORDS], pts[5]->Z[1][i], pts[5]->Z[0][i+VGWORDS], pts[5]->Z[0][i]);
        _coeff[0][i] = _SET(0, 0, 0, 0, cf01[i+VGWORDS], cf01[i], cf00[i+VGWORDS], cf00[i]); 
        _coeff[1][i] = _SET(cf11[i+VGWORDS], cf11[i], cf10[i+VGWORDS], cf10[i], 0, 0, 0, 0); 
        _coeff[2][i] = _SET(cf21[i+VGWORDS], cf21[i], cf20[i+VGWORDS], cf20[i], 0, 0, 0, 0); 
      }

      eval_4_isog_1x2x2x2w(_P, _coeff);

      // extract results
      get_channel_4x2w(pts[5]->X[0], _P, 4); get_channel_4x2w(pts[5]->X[1], _P, 6);
      get_channel_4x2w(pts[5]->Z[0], _P, 0); get_channel_4x2w(pts[5]->Z[1], _P, 2);

      // carry propagation 
      carryp_1w(pts[5]->X[0]); carryp_1w(pts[5]->X[1]); 
      carryp_1w(pts[5]->Z[0]); carryp_1w(pts[5]->Z[1]); 
    }

    // depends on public info
    if (num == 10) {
      // 2x2x2x1w < pts[5] | pts[6] >
      vfelm_t _P, _coeff[3];

      for (i = 0; i < VNWORDS; i++) {
        _P[i] = _SET(pts[5]->X[1][i], pts[5]->X[0][i], pts[5]->Z[1][i], pts[5]->Z[0][i], \
                     pts[6]->X[1][i], pts[6]->X[0][i], pts[6]->Z[1][i], pts[6]->Z[0][i]);
        _coeff[0][i] = _SET(0, 0, cf01[i], cf00[i], 0, 0, cf01[i], cf00[i]);
        _coeff[1][i] = _SET(cf11[i], cf10[i], 0, 0, cf11[i], cf10[i], 0, 0);
        _coeff[2][i] = _SET(cf21[i], cf20[i], 0, 0, cf21[i], cf20[i], 0, 0);
      }

      eval_4_isog_2x2x2x1w(_P, _coeff);

      // extract results
      get_channel_8x1w(pts[5]->X[0], _P, 6); get_channel_8x1w(pts[5]->X[1], _P, 7);
      get_channel_8x1w(pts[5]->Z[0], _P, 4); get_channel_8x1w(pts[5]->Z[1], _P, 5);
      get_channel_8x1w(pts[6]->X[0], _P, 2); get_channel_8x1w(pts[6]->X[1], _P, 3);
      get_channel_8x1w(pts[6]->Z[0], _P, 0); get_channel_8x1w(pts[6]->Z[1], _P, 1);
    }

    // depends on public info
    if (num == 11) {
      // 4x2x1x1w < pts[5] | pts[6] | pts[7] | % >
      vf2elm_t _P, _coeff[3];

      // initialize vectors 
      for (i = 0; i < VNWORDS; i++) {
        _P[0][i] = _SET(pts[5]->X[0][i], pts[5]->Z[0][i], pts[6]->X[0][i], pts[6]->Z[0][i], \
                        pts[7]->X[0][i], pts[7]->Z[0][i],               0,               0);
        _P[1][i] = _SET(pts[5]->X[1][i], pts[5]->Z[1][i], pts[6]->X[1][i], pts[6]->Z[1][i], \
                        pts[7]->X[1][i], pts[7]->Z[1][i],               0,               0);
        _coeff[0][0][i] = _SET(0, cf00[i], 0, cf00[i], 0, cf00[i], 0, cf00[i]);
        _coeff[0][1][i] = _SET(0, cf01[i], 0, cf01[i], 0, cf01[i], 0, cf01[i]);
        _coeff[1][0][i] = _SET(cf10[i], 0, cf10[i], 0, cf10[i], 0, cf10[i], 0);
        _coeff[1][1][i] = _SET(cf11[i], 0, cf11[i], 0, cf11[i], 0, cf11[i], 0);
        _coeff[2][0][i] = _SET(cf20[i], 0, cf20[i], 0, cf20[i], 0, cf20[i], 0);
        _coeff[2][1][i] = _SET(cf21[i], 0, cf21[i], 0, cf21[i], 0, cf21[i], 0);
      }

      eval_4_isog_4x2x1x1w(_P, _coeff);

      // extract results
      get_channel_8x1w(pts[5]->X[0], _P[0], 7); get_channel_8x1w(pts[5]->X[1], _P[1], 7);
      get_channel_8x1w(pts[5]->Z[0], _P[0], 6); get_channel_8x1w(pts[5]->Z[1], _P[1], 6);
      get_channel_8x1w(pts[6]->X[0], _P[0], 5); get_channel_8x1w(pts[6]->X[1], _P[1], 5);
      get_channel_8x1w(pts[6]->Z[0], _P[0], 4); get_channel_8x1w(pts[6]->Z[1], _P[1], 4);
      get_channel_8x1w(pts[7]->X[0], _P[0], 3); get_channel_8x1w(pts[7]->X[1], _P[1], 3);
      get_channel_8x1w(pts[7]->Z[0], _P[0], 2); get_channel_8x1w(pts[7]->Z[1], _P[1], 2);
    }
  }
  // depends on public info
  else {
    puts("Error!");
  }
}

void eval_4_isog_parallel_ss(point_proj_r51_t *pts, const vgelm_t coeff__0, const vgelm_t coeff2_1, const int num)
{
  felm_r51_t cf00, cf01, cf10, cf11, cf20, cf21;
  int i;

  // extract constants
  get_channel_4x2w(cf00, coeff__0, 0); get_channel_4x2w(cf01, coeff__0, 2);
  get_channel_4x2w(cf10, coeff2_1, 0); get_channel_4x2w(cf11, coeff2_1, 2);
  get_channel_4x2w(cf20, coeff2_1, 4); get_channel_4x2w(cf21, coeff2_1, 6);

  // carry propagation
  carryp_1w(cf00); carryp_1w(cf01);
  carryp_1w(cf10); carryp_1w(cf11);
  carryp_1w(cf20); carryp_1w(cf21);

  // depends on public info
  if (num == 1) {
    // 1x2x2x2w < pts[0] >
    vgelm_t _P, _coeff[3];

    // initialize vectors 
    for (i = 0; i < VGWORDS; i++) {
      _P[i] = _SET(pts[0]->X[1][i+VGWORDS], pts[0]->X[1][i], pts[0]->X[0][i+VGWORDS], pts[0]->X[0][i],\
                   pts[0]->Z[1][i+VGWORDS], pts[0]->Z[1][i], pts[0]->Z[0][i+VGWORDS], pts[0]->Z[0][i]);
      _coeff[0][i] = _SET(0, 0, 0, 0, cf01[i+VGWORDS], cf01[i], cf00[i+VGWORDS], cf00[i]); 
      _coeff[1][i] = _SET(cf11[i+VGWORDS], cf11[i], cf10[i+VGWORDS], cf10[i], 0, 0, 0, 0); 
      _coeff[2][i] = _SET(cf21[i+VGWORDS], cf21[i], cf20[i+VGWORDS], cf20[i], 0, 0, 0, 0); 
    }

    eval_4_isog_1x2x2x2w(_P, _coeff);

    // extract results
    get_channel_4x2w(pts[0]->X[0], _P, 4); get_channel_4x2w(pts[0]->X[1], _P, 6);
    get_channel_4x2w(pts[0]->Z[0], _P, 0); get_channel_4x2w(pts[0]->Z[1], _P, 2);

    // carry propagation 
    carryp_1w(pts[0]->X[0]); carryp_1w(pts[0]->X[1]); 
    carryp_1w(pts[0]->Z[0]); carryp_1w(pts[0]->Z[1]);  
  }
  // depends on public info
  else if (num == 2) {
    // 2x2x2x1w < pts[0] | pts[1] >
    vfelm_t _P, _coeff[3];

    for (i = 0; i < VNWORDS; i++) {
      _P[i] = _SET(pts[0]->X[1][i], pts[0]->X[0][i], pts[0]->Z[1][i], pts[0]->Z[0][i], \
                   pts[1]->X[1][i], pts[1]->X[0][i], pts[1]->Z[1][i], pts[1]->Z[0][i]);
      _coeff[0][i] = _SET(0, 0, cf01[i], cf00[i], 0, 0, cf01[i], cf00[i]);
      _coeff[1][i] = _SET(cf11[i], cf10[i], 0, 0, cf11[i], cf10[i], 0, 0);
      _coeff[2][i] = _SET(cf21[i], cf20[i], 0, 0, cf21[i], cf20[i], 0, 0);
    }

    eval_4_isog_2x2x2x1w(_P, _coeff);

    // extract results
    get_channel_8x1w(pts[0]->X[0], _P, 6); get_channel_8x1w(pts[0]->X[1], _P, 7);
    get_channel_8x1w(pts[0]->Z[0], _P, 4); get_channel_8x1w(pts[0]->Z[1], _P, 5);
    get_channel_8x1w(pts[1]->X[0], _P, 2); get_channel_8x1w(pts[1]->X[1], _P, 3);
    get_channel_8x1w(pts[1]->Z[0], _P, 0); get_channel_8x1w(pts[1]->Z[1], _P, 1);
  }
  // depends on public info
  else if (num == 3 | num == 4 | num == 5 | num == 6) {
    // 4x2x1x1w < pts[0] | pts[1] | pts[2] | pts[3] >
    vf2elm_t P, coeff[3];

    // initialize vectors 
    for (i = 0; i < VNWORDS; i++) {
      P[0][i] = _SET(pts[0]->X[0][i], pts[0]->Z[0][i], pts[1]->X[0][i], pts[1]->Z[0][i], \
                     pts[2]->X[0][i], pts[2]->Z[0][i], pts[3]->X[0][i], pts[3]->Z[0][i]);
      P[1][i] = _SET(pts[0]->X[1][i], pts[0]->Z[1][i], pts[1]->X[1][i], pts[1]->Z[1][i], \
                     pts[2]->X[1][i], pts[2]->Z[1][i], pts[3]->X[1][i], pts[3]->Z[1][i]);
      coeff[0][0][i] = _SET(0, cf00[i], 0, cf00[i], 0, cf00[i], 0, cf00[i]);
      coeff[0][1][i] = _SET(0, cf01[i], 0, cf01[i], 0, cf01[i], 0, cf01[i]);
      coeff[1][0][i] = _SET(cf10[i], 0, cf10[i], 0, cf10[i], 0, cf10[i], 0);
      coeff[1][1][i] = _SET(cf11[i], 0, cf11[i], 0, cf11[i], 0, cf11[i], 0);
      coeff[2][0][i] = _SET(cf20[i], 0, cf20[i], 0, cf20[i], 0, cf20[i], 0);
      coeff[2][1][i] = _SET(cf21[i], 0, cf21[i], 0, cf21[i], 0, cf21[i], 0);
    }

    eval_4_isog_4x2x1x1w(P, coeff);

    // extract results
    get_channel_8x1w(pts[0]->X[0], P[0], 7); get_channel_8x1w(pts[0]->X[1], P[1], 7);
    get_channel_8x1w(pts[0]->Z[0], P[0], 6); get_channel_8x1w(pts[0]->Z[1], P[1], 6);
    get_channel_8x1w(pts[1]->X[0], P[0], 5); get_channel_8x1w(pts[1]->X[1], P[1], 5);
    get_channel_8x1w(pts[1]->Z[0], P[0], 4); get_channel_8x1w(pts[1]->Z[1], P[1], 4);
    get_channel_8x1w(pts[2]->X[0], P[0], 3); get_channel_8x1w(pts[2]->X[1], P[1], 3);
    get_channel_8x1w(pts[2]->Z[0], P[0], 2); get_channel_8x1w(pts[2]->Z[1], P[1], 2);

    // depends on public info
    if (num >= 4) {
      get_channel_8x1w(pts[3]->X[0], P[0], 1); get_channel_8x1w(pts[3]->X[1], P[1], 1);
      get_channel_8x1w(pts[3]->Z[0], P[0], 0); get_channel_8x1w(pts[3]->Z[1], P[1], 0);
    }

    // depends on public info
    if (num == 5) {
      // 1x2x2x2w < pts[4] >
      vgelm_t _P, _coeff[3];

      // initialize vectors 
      for (i = 0; i < VGWORDS; i++) {
        _P[i] = _SET(pts[4]->X[1][i+VGWORDS], pts[4]->X[1][i], pts[4]->X[0][i+VGWORDS], pts[4]->X[0][i],\
                     pts[4]->Z[1][i+VGWORDS], pts[4]->Z[1][i], pts[4]->Z[0][i+VGWORDS], pts[4]->Z[0][i]);
        _coeff[0][i] = _SET(0, 0, 0, 0, cf01[i+VGWORDS], cf01[i], cf00[i+VGWORDS], cf00[i]); 
        _coeff[1][i] = _SET(cf11[i+VGWORDS], cf11[i], cf10[i+VGWORDS], cf10[i], 0, 0, 0, 0); 
        _coeff[2][i] = _SET(cf21[i+VGWORDS], cf21[i], cf20[i+VGWORDS], cf20[i], 0, 0, 0, 0); 
      }

      eval_4_isog_1x2x2x2w(_P, _coeff);

      // extract results
      get_channel_4x2w(pts[4]->X[0], _P, 4); get_channel_4x2w(pts[4]->X[1], _P, 6);
      get_channel_4x2w(pts[4]->Z[0], _P, 0); get_channel_4x2w(pts[4]->Z[1], _P, 2);

      // carry propagation 
      carryp_1w(pts[4]->X[0]); carryp_1w(pts[4]->X[1]); 
      carryp_1w(pts[4]->Z[0]); carryp_1w(pts[4]->Z[1]);  
    }

    // depends on public info
    if (num == 6) {
      // 2x2x2x1w < pts[4] | pts[5] >
      vfelm_t _P, _coeff[3];

      for (i = 0; i < VNWORDS; i++) {
        _P[i] = _SET(pts[4]->X[1][i], pts[4]->X[0][i], pts[4]->Z[1][i], pts[4]->Z[0][i], \
                     pts[5]->X[1][i], pts[5]->X[0][i], pts[5]->Z[1][i], pts[5]->Z[0][i]);
        _coeff[0][i] = _SET(0, 0, cf01[i], cf00[i], 0, 0, cf01[i], cf00[i]);
        _coeff[1][i] = _SET(cf11[i], cf10[i], 0, 0, cf11[i], cf10[i], 0, 0);
        _coeff[2][i] = _SET(cf21[i], cf20[i], 0, 0, cf21[i], cf20[i], 0, 0);
      }

      eval_4_isog_2x2x2x1w(_P, _coeff);

      // extract results
      get_channel_8x1w(pts[4]->X[0], _P, 6); get_channel_8x1w(pts[4]->X[1], _P, 7);
      get_channel_8x1w(pts[4]->Z[0], _P, 4); get_channel_8x1w(pts[4]->Z[1], _P, 5);
      get_channel_8x1w(pts[5]->X[0], _P, 2); get_channel_8x1w(pts[5]->X[1], _P, 3);
      get_channel_8x1w(pts[5]->Z[0], _P, 0); get_channel_8x1w(pts[5]->Z[1], _P, 1);
    }
  }
  // depends on public info
  else if (num == 7 | num == 8) {
    // 8x1x1x1w < pts[0] | pts[1] | pts[2] | pts[3] | pts[4] | pts[5] | pts[6] | pts[7] >
    vpoint_proj_t P;
    vf2elm_t coeff[3];

    // initialize vectors 
    for (i = 0; i < VNWORDS; i++) {
      P->X[0][i] = _SET(pts[0]->X[0][i], pts[1]->X[0][i], pts[2]->X[0][i], pts[3]->X[0][i], \
                        pts[4]->X[0][i], pts[5]->X[0][i], pts[6]->X[0][i], pts[7]->X[0][i]);
      P->X[1][i] = _SET(pts[0]->X[1][i], pts[1]->X[1][i], pts[2]->X[1][i], pts[3]->X[1][i], \
                        pts[4]->X[1][i], pts[5]->X[1][i], pts[6]->X[1][i], pts[7]->X[1][i]);
      P->Z[0][i] = _SET(pts[0]->Z[0][i], pts[1]->Z[0][i], pts[2]->Z[0][i], pts[3]->Z[0][i], \
                        pts[4]->Z[0][i], pts[5]->Z[0][i], pts[6]->Z[0][i], pts[7]->Z[0][i]);
      P->Z[1][i] = _SET(pts[0]->Z[1][i], pts[1]->Z[1][i], pts[2]->Z[1][i], pts[3]->Z[1][i], \
                        pts[4]->Z[1][i], pts[5]->Z[1][i], pts[6]->Z[1][i], pts[7]->Z[1][i]);
      coeff[0][0][i] = VSET1(cf00[i]); coeff[0][1][i] = VSET1(cf01[i]);
      coeff[1][0][i] = VSET1(cf10[i]); coeff[1][1][i] = VSET1(cf11[i]);
      coeff[2][0][i] = VSET1(cf20[i]); coeff[2][1][i] = VSET1(cf21[i]);
    }

    eval_4_isog_8x1x1x1w(P, coeff);

    // extract results
    get_channel_8x1w(pts[0]->X[0], P->X[0], 7); get_channel_8x1w(pts[0]->X[1], P->X[1], 7);
    get_channel_8x1w(pts[0]->Z[0], P->Z[0], 7); get_channel_8x1w(pts[0]->Z[1], P->Z[1], 7);
    get_channel_8x1w(pts[1]->X[0], P->X[0], 6); get_channel_8x1w(pts[1]->X[1], P->X[1], 6);
    get_channel_8x1w(pts[1]->Z[0], P->Z[0], 6); get_channel_8x1w(pts[1]->Z[1], P->Z[1], 6);
    get_channel_8x1w(pts[2]->X[0], P->X[0], 5); get_channel_8x1w(pts[2]->X[1], P->X[1], 5);
    get_channel_8x1w(pts[2]->Z[0], P->Z[0], 5); get_channel_8x1w(pts[2]->Z[1], P->Z[1], 5);
    get_channel_8x1w(pts[3]->X[0], P->X[0], 4); get_channel_8x1w(pts[3]->X[1], P->X[1], 4);
    get_channel_8x1w(pts[3]->Z[0], P->Z[0], 4); get_channel_8x1w(pts[3]->Z[1], P->Z[1], 4);
    get_channel_8x1w(pts[4]->X[0], P->X[0], 3); get_channel_8x1w(pts[4]->X[1], P->X[1], 3);
    get_channel_8x1w(pts[4]->Z[0], P->Z[0], 3); get_channel_8x1w(pts[4]->Z[1], P->Z[1], 3);
    get_channel_8x1w(pts[5]->X[0], P->X[0], 2); get_channel_8x1w(pts[5]->X[1], P->X[1], 2);
    get_channel_8x1w(pts[5]->Z[0], P->Z[0], 2); get_channel_8x1w(pts[5]->Z[1], P->Z[1], 2);
    get_channel_8x1w(pts[6]->X[0], P->X[0], 1); get_channel_8x1w(pts[6]->X[1], P->X[1], 1);
    get_channel_8x1w(pts[6]->Z[0], P->Z[0], 1); get_channel_8x1w(pts[6]->Z[1], P->Z[1], 1);

    // depends on public info
    if (num == 8) {
      get_channel_8x1w(pts[7]->X[0], P->X[0], 0); get_channel_8x1w(pts[7]->X[1], P->X[1], 0);
      get_channel_8x1w(pts[7]->Z[0], P->Z[0], 0); get_channel_8x1w(pts[7]->Z[1], P->Z[1], 0);
    }
  }
  // depends on public info
  else {
    puts("Error!");    
  }
}

void eval_4_isog_parallel_kgss(point_proj_r51_t *pts, point_proj_r51_t *_pts, point_proj_r51_t phiP, point_proj_r51_t phiQ, point_proj_r51_t phiR, const vfelm_t coeff__0, const vfelm_t coeff2_1, const int num)
{
  felm_r51_t cf00, cf01, cf10, cf11, cf20, cf21;
  felm_r51_t _cf00, _cf01, _cf10, _cf11, _cf20, _cf21;
  int i;

  // extract constants
  get_channel_8x1w(cf00, coeff__0, 4); get_channel_8x1w(cf01, coeff__0, 5);
  get_channel_8x1w(cf10, coeff2_1, 4); get_channel_8x1w(cf11, coeff2_1, 5);
  get_channel_8x1w(cf20, coeff2_1, 6); get_channel_8x1w(cf21, coeff2_1, 7);
  get_channel_8x1w(_cf00, coeff__0, 0); get_channel_8x1w(_cf01, coeff__0, 1);
  get_channel_8x1w(_cf10, coeff2_1, 0); get_channel_8x1w(_cf11, coeff2_1, 1);
  get_channel_8x1w(_cf20, coeff2_1, 2); get_channel_8x1w(_cf21, coeff2_1, 3);

  // num = 2*npts+3 could be 3, 5, 7, 9, 11, 13, 15, 17, 19
  // depends on public info
  if (num == 3 | num == 5) {
    // 4x2x1x1w < pts[0] | phiP | phiQ | phiR >
    vf2elm_t P, coeff[3];

    // initialize vectors 
    for (i = 0; i < VNWORDS; i++) {
      P[0][i] = _SET(pts[0]->X[0][i], pts[0]->Z[0][i], phiP->X[0][i], phiP->Z[0][i], \
                       phiQ->X[0][i],   phiQ->Z[0][i], phiR->X[0][i], phiR->Z[0][i]);
      P[1][i] = _SET(pts[0]->X[1][i], pts[0]->Z[1][i], phiP->X[1][i], phiP->Z[1][i], \
                       phiQ->X[1][i],   phiQ->Z[1][i], phiR->X[1][i], phiR->Z[1][i]);
      coeff[0][0][i] = _SET(0, cf00[i], 0, cf00[i], 0, cf00[i], 0, cf00[i]);
      coeff[0][1][i] = _SET(0, cf01[i], 0, cf01[i], 0, cf01[i], 0, cf01[i]);
      coeff[1][0][i] = _SET(cf10[i], 0, cf10[i], 0, cf10[i], 0, cf10[i], 0);
      coeff[1][1][i] = _SET(cf11[i], 0, cf11[i], 0, cf11[i], 0, cf11[i], 0);
      coeff[2][0][i] = _SET(cf20[i], 0, cf20[i], 0, cf20[i], 0, cf20[i], 0);
      coeff[2][1][i] = _SET(cf21[i], 0, cf21[i], 0, cf21[i], 0, cf21[i], 0);
    }

    eval_4_isog_4x2x1x1w(P, coeff);

    // extract results
    get_channel_8x1w(phiP->X[0], P[0], 5);   get_channel_8x1w(phiP->X[1], P[1], 5);
    get_channel_8x1w(phiP->Z[0], P[0], 4);   get_channel_8x1w(phiP->Z[1], P[1], 4);
    get_channel_8x1w(phiQ->X[0], P[0], 3);   get_channel_8x1w(phiQ->X[1], P[1], 3);
    get_channel_8x1w(phiQ->Z[0], P[0], 2);   get_channel_8x1w(phiQ->Z[1], P[1], 2);
    get_channel_8x1w(phiR->X[0], P[0], 1);   get_channel_8x1w(phiR->X[1], P[1], 1);
    get_channel_8x1w(phiR->Z[0], P[0], 0);   get_channel_8x1w(phiR->Z[1], P[1], 0);

    // depends on public info
    if (num == 5) {
      get_channel_8x1w(pts[0]->X[0], P[0], 7); get_channel_8x1w(pts[0]->X[1], P[1], 7);
      get_channel_8x1w(pts[0]->Z[0], P[0], 6); get_channel_8x1w(pts[0]->Z[1], P[1], 6);

      // 1x2x2x2w < _pts[0] >
      vgelm_t _P, _coeff[3];

      // initialize vectors 
      for (i = 0; i < VGWORDS; i++) {
        _P[i] = _SET(_pts[0]->X[1][i+VGWORDS], _pts[0]->X[1][i], _pts[0]->X[0][i+VGWORDS], _pts[0]->X[0][i],\
                     _pts[0]->Z[1][i+VGWORDS], _pts[0]->Z[1][i], _pts[0]->Z[0][i+VGWORDS], _pts[0]->Z[0][i]);
        _coeff[0][i] = _SET(0, 0, 0, 0, _cf01[i+VGWORDS], _cf01[i], _cf00[i+VGWORDS], _cf00[i]); 
        _coeff[1][i] = _SET(_cf11[i+VGWORDS], _cf11[i], _cf10[i+VGWORDS], _cf10[i], 0, 0, 0, 0); 
        _coeff[2][i] = _SET(_cf21[i+VGWORDS], _cf21[i], _cf20[i+VGWORDS], _cf20[i], 0, 0, 0, 0); 
      }

      eval_4_isog_1x2x2x2w(_P, _coeff);

      // extract results
      get_channel_4x2w(_pts[0]->X[0], _P, 4); get_channel_4x2w(_pts[0]->X[1], _P, 6);
      get_channel_4x2w(_pts[0]->Z[0], _P, 0); get_channel_4x2w(_pts[0]->Z[1], _P, 2);

      // carry propagation 
      carryp_1w(_pts[0]->X[0]); carryp_1w(_pts[0]->X[1]); 
      carryp_1w(_pts[0]->Z[0]); carryp_1w(_pts[0]->Z[1]);  
    }
  }
  // depends on public info
  else if (num == 7 | num == 9 | num == 11 | num == 13 | num == 15 | num == 17 | num == 19) {
    // 8x1x1x1w < _pts[0] | _pts[1] | pts[0] | pts[1] | pts[2] | phiP | phiQ | phiR >
    vpoint_proj_t P;
    vf2elm_t coeff[3];

    // initialize vectors 
    for (i = 0; i < VNWORDS; i++) {
      P->X[0][i] = _SET(_pts[0]->X[0][i], _pts[1]->X[0][i], pts[0]->X[0][i], pts[1]->X[0][i], \
                         pts[2]->X[0][i],    phiP->X[0][i],   phiQ->X[0][i],   phiR->X[0][i]);
      P->X[1][i] = _SET(_pts[0]->X[1][i], _pts[1]->X[1][i], pts[0]->X[1][i], pts[1]->X[1][i], \
                         pts[2]->X[1][i],    phiP->X[1][i],   phiQ->X[1][i],   phiR->X[1][i]);
      P->Z[0][i] = _SET(_pts[0]->Z[0][i], _pts[1]->Z[0][i], pts[0]->Z[0][i], pts[1]->Z[0][i], \
                         pts[2]->Z[0][i],    phiP->Z[0][i],   phiQ->Z[0][i],   phiR->Z[0][i]);
      P->Z[1][i] = _SET(_pts[0]->Z[1][i], _pts[1]->Z[1][i], pts[0]->Z[1][i], pts[1]->Z[1][i], \
                         pts[2]->Z[1][i],    phiP->Z[1][i],   phiQ->Z[1][i],   phiR->Z[1][i]);
      coeff[0][0][i] = _SET(_cf00[i], _cf00[i], cf00[i], cf00[i], cf00[i], cf00[i], cf00[i], cf00[i]); 
      coeff[0][1][i] = _SET(_cf01[i], _cf01[i], cf01[i], cf01[i], cf01[i], cf01[i], cf01[i], cf01[i]);
      coeff[1][0][i] = _SET(_cf10[i], _cf10[i], cf10[i], cf10[i], cf10[i], cf10[i], cf10[i], cf10[i]); 
      coeff[1][1][i] = _SET(_cf11[i], _cf11[i], cf11[i], cf11[i], cf11[i], cf11[i], cf11[i], cf11[i]);
      coeff[2][0][i] = _SET(_cf20[i], _cf20[i], cf20[i], cf20[i], cf20[i], cf20[i], cf20[i], cf20[i]); 
      coeff[2][1][i] = _SET(_cf21[i], _cf21[i], cf21[i], cf21[i], cf21[i], cf21[i], cf21[i], cf21[i]);
    }

    eval_4_isog_8x1x1x1w(P, coeff);

    // extract results
    get_channel_8x1w(_pts[0]->X[0], P->X[0], 7); get_channel_8x1w(_pts[0]->X[1], P->X[1], 7);
    get_channel_8x1w(_pts[0]->Z[0], P->Z[0], 7); get_channel_8x1w(_pts[0]->Z[1], P->Z[1], 7);
    get_channel_8x1w(_pts[1]->X[0], P->X[0], 6); get_channel_8x1w(_pts[1]->X[1], P->X[1], 6);
    get_channel_8x1w(_pts[1]->Z[0], P->Z[0], 6); get_channel_8x1w(_pts[1]->Z[1], P->Z[1], 6);
    get_channel_8x1w(pts[0]->X[0], P->X[0], 5); get_channel_8x1w(pts[0]->X[1], P->X[1], 5);
    get_channel_8x1w(pts[0]->Z[0], P->Z[0], 5); get_channel_8x1w(pts[0]->Z[1], P->Z[1], 5);
    get_channel_8x1w(pts[1]->X[0], P->X[0], 4); get_channel_8x1w(pts[1]->X[1], P->X[1], 4);
    get_channel_8x1w(pts[1]->Z[0], P->Z[0], 4); get_channel_8x1w(pts[1]->Z[1], P->Z[1], 4);
    get_channel_8x1w(phiP->X[0], P->X[0], 2);   get_channel_8x1w(phiP->X[1], P->X[1], 2);
    get_channel_8x1w(phiP->Z[0], P->Z[0], 2);   get_channel_8x1w(phiP->Z[1], P->Z[1], 2);
    get_channel_8x1w(phiQ->X[0], P->X[0], 1);   get_channel_8x1w(phiQ->X[1], P->X[1], 1);
    get_channel_8x1w(phiQ->Z[0], P->Z[0], 1);   get_channel_8x1w(phiQ->Z[1], P->Z[1], 1);
    get_channel_8x1w(phiR->X[0], P->X[0], 0);   get_channel_8x1w(phiR->X[1], P->X[1], 0);
    get_channel_8x1w(phiR->Z[0], P->Z[0], 0);   get_channel_8x1w(phiR->Z[1], P->Z[1], 0);

    // depends on public info
    if (num >= 9) {
      get_channel_8x1w(pts[2]->X[0], P->X[0], 3); get_channel_8x1w(pts[2]->X[1], P->X[1], 3);
      get_channel_8x1w(pts[2]->Z[0], P->Z[0], 3); get_channel_8x1w(pts[2]->Z[1], P->Z[1], 3);
    }

    // depends on public info
    if (num == 9) {
      // 1x2x2x2w < _pts[2] >       
      vgelm_t _P, _coeff[3];

      // initialize vectors 
      for (i = 0; i < VGWORDS; i++) {
        _P[i] = _SET(_pts[2]->X[1][i+VGWORDS], _pts[2]->X[1][i], _pts[2]->X[0][i+VGWORDS], _pts[2]->X[0][i],\
                     _pts[2]->Z[1][i+VGWORDS], _pts[2]->Z[1][i], _pts[2]->Z[0][i+VGWORDS], _pts[2]->Z[0][i]);
        _coeff[0][i] = _SET(0, 0, 0, 0, _cf01[i+VGWORDS], _cf01[i], _cf00[i+VGWORDS], _cf00[i]); 
        _coeff[1][i] = _SET(_cf11[i+VGWORDS], _cf11[i], _cf10[i+VGWORDS], _cf10[i], 0, 0, 0, 0); 
        _coeff[2][i] = _SET(_cf21[i+VGWORDS], _cf21[i], _cf20[i+VGWORDS], _cf20[i], 0, 0, 0, 0); 
      }

      eval_4_isog_1x2x2x2w(_P, _coeff);    

      // extract results
      get_channel_4x2w(_pts[2]->X[0], _P, 4); get_channel_4x2w(_pts[2]->X[1], _P, 6);
      get_channel_4x2w(_pts[2]->Z[0], _P, 0); get_channel_4x2w(_pts[2]->Z[1], _P, 2);

      // carry propagation 
      carryp_1w(_pts[2]->X[0]); carryp_1w(_pts[2]->X[1]); 
      carryp_1w(_pts[2]->Z[0]); carryp_1w(_pts[2]->Z[1]); 
    }

    // depends on public info
    if (num == 11 | num == 13) {
      // 4x2x1x1w < _pts[2] | _pts[3] | pts[3] | pts[4] >
      vf2elm_t _P, _coeff[3];

      // initialize vectors 
      for (i = 0; i < VNWORDS; i++) {
        _P[0][i] = _SET(_pts[2]->X[0][i], _pts[2]->Z[0][i], _pts[3]->X[0][i], _pts[3]->Z[0][i], \
                         pts[3]->X[0][i],  pts[3]->Z[0][i],  pts[4]->X[0][i],  pts[4]->Z[0][i]);
        _P[1][i] = _SET(_pts[2]->X[1][i], _pts[2]->Z[1][i], _pts[3]->X[1][i], _pts[3]->Z[1][i], \
                         pts[3]->X[1][i],  pts[3]->Z[1][i],  pts[4]->X[1][i],  pts[4]->Z[1][i]);
        _coeff[0][0][i] = _SET(0, _cf00[i], 0, _cf00[i], 0, cf00[i], 0, cf00[i]);
        _coeff[0][1][i] = _SET(0, _cf01[i], 0, _cf01[i], 0, cf01[i], 0, cf01[i]);
        _coeff[1][0][i] = _SET(_cf10[i], 0, _cf10[i], 0, cf10[i], 0, cf10[i], 0);
        _coeff[1][1][i] = _SET(_cf11[i], 0, _cf11[i], 0, cf11[i], 0, cf11[i], 0);
        _coeff[2][0][i] = _SET(_cf20[i], 0, _cf20[i], 0, cf20[i], 0, cf20[i], 0);
        _coeff[2][1][i] = _SET(_cf21[i], 0, _cf21[i], 0, cf21[i], 0, cf21[i], 0);
      }

      eval_4_isog_4x2x1x1w(_P, _coeff);

      // extract results
      get_channel_8x1w(_pts[2]->X[0], _P[0], 7); get_channel_8x1w(_pts[2]->X[1], _P[1], 7);
      get_channel_8x1w(_pts[2]->Z[0], _P[0], 6); get_channel_8x1w(_pts[2]->Z[1], _P[1], 6);
      get_channel_8x1w(_pts[3]->X[0], _P[0], 5); get_channel_8x1w(_pts[3]->X[1], _P[1], 5);
      get_channel_8x1w(_pts[3]->Z[0], _P[0], 4); get_channel_8x1w(_pts[3]->Z[1], _P[1], 4);
      get_channel_8x1w(pts[3]->X[0], _P[0], 3); get_channel_8x1w(pts[3]->X[1], _P[1], 3);
      get_channel_8x1w(pts[3]->Z[0], _P[0], 2); get_channel_8x1w(pts[3]->Z[1], _P[1], 2);

      // depends on public info
      if (num == 13) {
        get_channel_8x1w(pts[4]->X[0], _P[0], 1); get_channel_8x1w(pts[4]->X[1], _P[1], 1);
        get_channel_8x1w(pts[4]->Z[0], _P[0], 0); get_channel_8x1w(pts[4]->Z[1], _P[1], 0);

        // 1x2x2x2w < _pts[4] >
        vgelm_t __P, __coeff[3];

        // initialize vectors 
        for (i = 0; i < VGWORDS; i++) {
          __P[i] = _SET(_pts[4]->X[1][i+VGWORDS], _pts[4]->X[1][i], _pts[4]->X[0][i+VGWORDS], _pts[4]->X[0][i],\
                        _pts[4]->Z[1][i+VGWORDS], _pts[4]->Z[1][i], _pts[4]->Z[0][i+VGWORDS], _pts[4]->Z[0][i]);
          __coeff[0][i] = _SET(0, 0, 0, 0, _cf01[i+VGWORDS], _cf01[i], _cf00[i+VGWORDS], _cf00[i]); 
          __coeff[1][i] = _SET(_cf11[i+VGWORDS], _cf11[i], _cf10[i+VGWORDS], _cf10[i], 0, 0, 0, 0); 
          __coeff[2][i] = _SET(_cf21[i+VGWORDS], _cf21[i], _cf20[i+VGWORDS], _cf20[i], 0, 0, 0, 0); 
        }

        eval_4_isog_1x2x2x2w(__P, __coeff);

        // extract results
        get_channel_4x2w(_pts[4]->X[0], __P, 4); get_channel_4x2w(_pts[4]->X[1], __P, 6);
        get_channel_4x2w(_pts[4]->Z[0], __P, 0); get_channel_4x2w(_pts[4]->Z[1], __P, 2);

        // carry propagation 
        carryp_1w(_pts[4]->X[0]); carryp_1w(_pts[4]->X[1]); 
        carryp_1w(_pts[4]->Z[0]); carryp_1w(_pts[4]->Z[1]); 
      }
    }

    // depends on public info
    if (num == 15 | num == 17 | num == 19) {
      // 8x1x1x1w < _pts[2] | _pts[3] | _pts[4] | _pts[5] | pts[6] | pts[3] | pts[4] | pts[5] >
      vpoint_proj_t _P;
      vf2elm_t _coeff[3];

      // initialize vectors 
      for (i = 0; i < VNWORDS; i++) {
        _P->X[0][i] = _SET(_pts[2]->X[0][i], _pts[3]->X[0][i], _pts[4]->X[0][i], _pts[5]->X[0][i], \
                            pts[6]->X[0][i],  pts[3]->X[0][i],  pts[4]->X[0][i],  pts[5]->X[0][i]);
        _P->X[1][i] = _SET(_pts[2]->X[1][i], _pts[3]->X[1][i], _pts[4]->X[1][i], _pts[5]->X[1][i], \
                            pts[6]->X[1][i],  pts[3]->X[1][i],  pts[4]->X[1][i],  pts[5]->X[1][i]);
        _P->Z[0][i] = _SET(_pts[2]->Z[0][i], _pts[3]->Z[0][i], _pts[4]->Z[0][i], _pts[5]->Z[0][i], \
                            pts[6]->Z[0][i],  pts[3]->Z[0][i],  pts[4]->Z[0][i],  pts[5]->Z[0][i]);
        _P->Z[1][i] = _SET(_pts[2]->Z[1][i], _pts[3]->Z[1][i], _pts[4]->Z[1][i], _pts[5]->Z[1][i], \
                            pts[6]->Z[1][i],  pts[3]->Z[1][i],  pts[4]->Z[1][i],  pts[5]->Z[1][i]);
        _coeff[0][0][i] = _SET(_cf00[i], _cf00[i], _cf00[i], _cf00[i], cf00[i], cf00[i], cf00[i], cf00[i]); 
        _coeff[0][1][i] = _SET(_cf01[i], _cf01[i], _cf01[i], _cf01[i], cf01[i], cf01[i], cf01[i], cf01[i]);
        _coeff[1][0][i] = _SET(_cf10[i], _cf10[i], _cf10[i], _cf10[i], cf10[i], cf10[i], cf10[i], cf10[i]); 
        _coeff[1][1][i] = _SET(_cf11[i], _cf11[i], _cf11[i], _cf11[i], cf11[i], cf11[i], cf11[i], cf11[i]);
        _coeff[2][0][i] = _SET(_cf20[i], _cf20[i], _cf20[i], _cf20[i], cf20[i], cf20[i], cf20[i], cf20[i]); 
        _coeff[2][1][i] = _SET(_cf21[i], _cf21[i], _cf21[i], _cf21[i], cf21[i], cf21[i], cf21[i], cf21[i]);
      }

      eval_4_isog_8x1x1x1w(_P, _coeff);

      // extract results
      get_channel_8x1w(_pts[2]->X[0], _P->X[0], 7); get_channel_8x1w(_pts[2]->X[1], _P->X[1], 7);
      get_channel_8x1w(_pts[2]->Z[0], _P->Z[0], 7); get_channel_8x1w(_pts[2]->Z[1], _P->Z[1], 7);
      get_channel_8x1w(_pts[3]->X[0], _P->X[0], 6); get_channel_8x1w(_pts[3]->X[1], _P->X[1], 6);
      get_channel_8x1w(_pts[3]->Z[0], _P->Z[0], 6); get_channel_8x1w(_pts[3]->Z[1], _P->Z[1], 6);
      get_channel_8x1w(_pts[4]->X[0], _P->X[0], 5); get_channel_8x1w(_pts[4]->X[1], _P->X[1], 5);
      get_channel_8x1w(_pts[4]->Z[0], _P->Z[0], 5); get_channel_8x1w(_pts[4]->Z[1], _P->Z[1], 5);
      get_channel_8x1w(_pts[5]->X[0], _P->X[0], 4); get_channel_8x1w(_pts[5]->X[1], _P->X[1], 4);
      get_channel_8x1w(_pts[5]->Z[0], _P->Z[0], 4); get_channel_8x1w(_pts[5]->Z[1], _P->Z[1], 4);
      get_channel_8x1w(pts[3]->X[0], _P->X[0], 2); get_channel_8x1w(pts[3]->X[1], _P->X[1], 2);
      get_channel_8x1w(pts[3]->Z[0], _P->Z[0], 2); get_channel_8x1w(pts[3]->Z[1], _P->Z[1], 2);
      get_channel_8x1w(pts[4]->X[0], _P->X[0], 1); get_channel_8x1w(pts[4]->X[1], _P->X[1], 1);
      get_channel_8x1w(pts[4]->Z[0], _P->Z[0], 1); get_channel_8x1w(pts[4]->Z[1], _P->Z[1], 1);
      get_channel_8x1w(pts[5]->X[0], _P->X[0], 0); get_channel_8x1w(pts[5]->X[1], _P->X[1], 0);
      get_channel_8x1w(pts[5]->Z[0], _P->Z[0], 0); get_channel_8x1w(pts[5]->Z[1], _P->Z[1], 0);

      // depends on public info
      if (num >= 17) {
        get_channel_8x1w(pts[6]->X[0], _P->X[0], 3); get_channel_8x1w(pts[6]->X[1], _P->X[1], 3);
        get_channel_8x1w(pts[6]->Z[0], _P->Z[0], 3); get_channel_8x1w(pts[6]->Z[1], _P->Z[1], 3);
      }

      // depends on public info
      if (num == 17) {
        // 1x2x2x2w < _pts[6] >
        vgelm_t __P, __coeff[3];

        // initialize vectors 
        for (i = 0; i < VGWORDS; i++) {
          __P[i] = _SET(_pts[6]->X[1][i+VGWORDS], _pts[6]->X[1][i], _pts[6]->X[0][i+VGWORDS], _pts[6]->X[0][i],\
                        _pts[6]->Z[1][i+VGWORDS], _pts[6]->Z[1][i], _pts[6]->Z[0][i+VGWORDS], _pts[6]->Z[0][i]);
          __coeff[0][i] = _SET(0, 0, 0, 0, _cf01[i+VGWORDS], _cf01[i], _cf00[i+VGWORDS], _cf00[i]); 
          __coeff[1][i] = _SET(_cf11[i+VGWORDS], _cf11[i], _cf10[i+VGWORDS], _cf10[i], 0, 0, 0, 0); 
          __coeff[2][i] = _SET(_cf21[i+VGWORDS], _cf21[i], _cf20[i+VGWORDS], _cf20[i], 0, 0, 0, 0); 
        }

        eval_4_isog_1x2x2x2w(__P, __coeff);

        // extract results
        get_channel_4x2w(_pts[6]->X[0], __P, 4); get_channel_4x2w(_pts[6]->X[1], __P, 6);
        get_channel_4x2w(_pts[6]->Z[0], __P, 0); get_channel_4x2w(_pts[6]->Z[1], __P, 2);

        // carry propagation 
        carryp_1w(_pts[6]->X[0]); carryp_1w(_pts[6]->X[1]); 
        carryp_1w(_pts[6]->Z[0]); carryp_1w(_pts[6]->Z[1]); 
      }

      // depends on public info
      if (num == 19) {
        // 4x2x1x1w < _pts[6] | _pts[7] | pts[7] | pts[7] >
        vf2elm_t __P, __coeff[3];

        // initialize vectors 
        for (i = 0; i < VNWORDS; i++) {
          __P[0][i] = _SET(_pts[6]->X[0][i], _pts[6]->Z[0][i], _pts[7]->X[0][i], _pts[7]->Z[0][i], \
                            pts[7]->X[0][i],  pts[7]->Z[0][i],  pts[7]->X[0][i],  pts[7]->Z[0][i]);
          __P[1][i] = _SET(_pts[6]->X[1][i], _pts[6]->Z[1][i], _pts[7]->X[1][i], _pts[7]->Z[1][i], \
                            pts[7]->X[1][i],  pts[7]->Z[1][i],  pts[7]->X[1][i],  pts[7]->Z[1][i]);
          __coeff[0][0][i] = _SET(0, _cf00[i], 0, _cf00[i], 0, cf00[i], 0, cf00[i]);
          __coeff[0][1][i] = _SET(0, _cf01[i], 0, _cf01[i], 0, cf01[i], 0, cf01[i]);
          __coeff[1][0][i] = _SET(_cf10[i], 0, _cf10[i], 0, cf10[i], 0, cf10[i], 0);
          __coeff[1][1][i] = _SET(_cf11[i], 0, _cf11[i], 0, cf11[i], 0, cf11[i], 0);
          __coeff[2][0][i] = _SET(_cf20[i], 0, _cf20[i], 0, cf20[i], 0, cf20[i], 0);
          __coeff[2][1][i] = _SET(_cf21[i], 0, _cf21[i], 0, cf21[i], 0, cf21[i], 0);
        }

        eval_4_isog_4x2x1x1w(__P, __coeff);

        // extract results
        get_channel_8x1w(_pts[6]->X[0], __P[0], 7); get_channel_8x1w(_pts[6]->X[1], __P[1], 7);
        get_channel_8x1w(_pts[6]->Z[0], __P[0], 6); get_channel_8x1w(_pts[6]->Z[1], __P[1], 6);
        get_channel_8x1w(_pts[7]->X[0], __P[0], 5); get_channel_8x1w(_pts[7]->X[1], __P[1], 5);
        get_channel_8x1w(_pts[7]->Z[0], __P[0], 4); get_channel_8x1w(_pts[7]->Z[1], __P[1], 4);
        get_channel_8x1w(pts[7]->X[0], __P[0], 3); get_channel_8x1w(pts[7]->X[1], __P[1], 3);
        get_channel_8x1w(pts[7]->Z[0], __P[0], 2); get_channel_8x1w(pts[7]->Z[1], __P[1], 2);        
      }
    }
  }
  // depends on public info
  else {
    puts("Error!");
  }
}

void xTPLe_1x2x2x2w(vgelm_t vP, point_proj_r51_t Q, const vgelm_t C24_A24plus, const int e)
{
  int i;

  for (i = 0; i < e; i++) xTPL_1x2x2x2w(vP, vP, C24_A24plus);

  // extract the point Q from vP
  get_channel_4x2w(Q->X[0], vP, 4); get_channel_4x2w(Q->X[1], vP, 6);
  get_channel_4x2w(Q->Z[0], vP, 0); get_channel_4x2w(Q->Z[1], vP, 2);

  // final carry propagation making Q strictly in radix-2^51
  carryp_1w(Q->X[0]); carryp_1w(Q->X[1]);
  carryp_1w(Q->Z[0]); carryp_1w(Q->Z[1]);
}

void eval_3_isog_parallel_kg(point_proj_r51_t *pts, point_proj_r51_t phiP, point_proj_r51_t phiQ, point_proj_r51_t phiR, const vgelm_t coeff1_0, const int num)
{
  felm_r51_t cf00, cf01, cf10, cf11;
  int i;

  // extract constants
  get_channel_4x2w(cf00, coeff1_0, 0); get_channel_4x2w(cf01, coeff1_0, 2);
  get_channel_4x2w(cf10, coeff1_0, 4); get_channel_4x2w(cf11, coeff1_0, 6);

  // carry propagation
  carryp_1w(cf00); carryp_1w(cf01);
  carryp_1w(cf10); carryp_1w(cf11);

  // depends on public info
  if (num == 3 | num == 4 | num == 5 | num == 6) {
    // 4x2x1x1w < pts[0] | phiP | phiQ | phiR >
    vf2elm_t P, coeff0_1;

    // initialize vectors 
    for (i = 0; i < VNWORDS; i++) {
      P[0][i] = _SET(pts[0]->X[0][i], pts[0]->Z[0][i], phiP->X[0][i], phiP->Z[0][i], \
                       phiQ->X[0][i],   phiQ->Z[0][i], phiR->X[0][i], phiR->Z[0][i]);
      P[1][i] = _SET(pts[0]->X[1][i], pts[0]->Z[1][i], phiP->X[1][i], phiP->Z[1][i], \
                       phiQ->X[1][i],   phiQ->Z[1][i], phiR->X[1][i], phiR->Z[1][i]);
      coeff0_1[0][i] = _SET(cf00[i], cf10[i], cf00[i], cf10[i], cf00[i], cf10[i], cf00[i], cf10[i]);
      coeff0_1[1][i] = _SET(cf01[i], cf11[i], cf01[i], cf11[i], cf01[i], cf11[i], cf01[i], cf11[i]);
    }

    eval_3_isog_4x2x1x1w(P, coeff0_1);

    // extract results
    get_channel_8x1w(phiP->X[0], P[0], 5);   get_channel_8x1w(phiP->X[1], P[1], 5);
    get_channel_8x1w(phiP->Z[0], P[0], 4);   get_channel_8x1w(phiP->Z[1], P[1], 4);
    get_channel_8x1w(phiQ->X[0], P[0], 3);   get_channel_8x1w(phiQ->X[1], P[1], 3);
    get_channel_8x1w(phiQ->Z[0], P[0], 2);   get_channel_8x1w(phiQ->Z[1], P[1], 2);
    get_channel_8x1w(phiR->X[0], P[0], 1);   get_channel_8x1w(phiR->X[1], P[1], 1);
    get_channel_8x1w(phiR->Z[0], P[0], 0);   get_channel_8x1w(phiR->Z[1], P[1], 0);
  
    // depends on public info
    if (num >= 4) {
      get_channel_8x1w(pts[0]->X[0], P[0], 7); get_channel_8x1w(pts[0]->X[1], P[1], 7);
      get_channel_8x1w(pts[0]->Z[0], P[0], 6); get_channel_8x1w(pts[0]->Z[1], P[1], 6);
    }

    // depends on public info
    if (num == 5) {
      // 1x2x2x2w < pts[1] >
      vgelm_t _P, _coeff0_1;

      // initialize vectors 
      for (i = 0; i < VGWORDS; i++) {
        _P[i] = _SET(pts[1]->X[1][i+VGWORDS], pts[1]->X[1][i], pts[1]->X[0][i+VGWORDS], pts[1]->X[0][i],\
                     pts[1]->Z[1][i+VGWORDS], pts[1]->Z[1][i], pts[1]->Z[0][i+VGWORDS], pts[1]->Z[0][i]);
        _coeff0_1[i] = _SET(cf01[i+VGWORDS], cf01[i], cf00[i+VGWORDS], cf00[i], 
                            cf11[i+VGWORDS], cf11[i], cf10[i+VGWORDS], cf10[i] ); 
      }

      eval_3_isog_1x2x2x2w(_P, _coeff0_1);

      // extract results
      get_channel_4x2w(pts[1]->X[0], _P, 4); get_channel_4x2w(pts[1]->X[1], _P, 6);
      get_channel_4x2w(pts[1]->Z[0], _P, 0); get_channel_4x2w(pts[1]->Z[1], _P, 2);

      // carry propagation 
      carryp_1w(pts[1]->X[0]); carryp_1w(pts[1]->X[1]); 
      carryp_1w(pts[1]->Z[0]); carryp_1w(pts[1]->Z[1]);  
    }

    // depends on public info
    if (num == 6) {
      // 2x2x2x1w < pts[1] | pts[2] >
      vfelm_t _P, _coeff0_1;

      for (i = 0; i < VNWORDS; i++) {
        _P[i] = _SET(pts[1]->X[1][i], pts[1]->X[0][i], pts[1]->Z[1][i], pts[1]->Z[0][i], \
                     pts[2]->X[1][i], pts[2]->X[0][i], pts[2]->Z[1][i], pts[2]->Z[0][i]);
        _coeff0_1[i] = _SET(cf01[i], cf00[i], cf11[i], cf10[i], cf01[i], cf00[i], cf11[i], cf10[i]);
      }

      eval_3_isog_2x2x2x1w(_P, _coeff0_1);

      // extract results
      get_channel_8x1w(pts[1]->X[0], _P, 6); get_channel_8x1w(pts[1]->X[1], _P, 7);
      get_channel_8x1w(pts[1]->Z[0], _P, 4); get_channel_8x1w(pts[1]->Z[1], _P, 5);
      get_channel_8x1w(pts[2]->X[0], _P, 2); get_channel_8x1w(pts[2]->X[1], _P, 3);
      get_channel_8x1w(pts[2]->Z[0], _P, 0); get_channel_8x1w(pts[2]->Z[1], _P, 1);
    }
  }
  // depends on public info
  else if (num == 7 | num == 8 | num == 9 | num == 10 | num == 11 | num == 12 | num == 13) {
    // 8x1x1x1w < pts[0] | pts[1] | pts[2] | pts[3] | pts[4] | phiP | phiQ | phiR >
    vpoint_proj_t P;
    vf2elm_t coeff[2];

    // initialize vectors 
    for (i = 0; i < VNWORDS; i++) {
      P->X[0][i] = _SET(pts[0]->X[0][i], pts[1]->X[0][i], pts[2]->X[0][i], pts[3]->X[0][i], \
                        pts[4]->X[0][i],   phiP->X[0][i],   phiQ->X[0][i],   phiR->X[0][i]);
      P->X[1][i] = _SET(pts[0]->X[1][i], pts[1]->X[1][i], pts[2]->X[1][i], pts[3]->X[1][i], \
                        pts[4]->X[1][i],   phiP->X[1][i],   phiQ->X[1][i],   phiR->X[1][i]);
      P->Z[0][i] = _SET(pts[0]->Z[0][i], pts[1]->Z[0][i], pts[2]->Z[0][i], pts[3]->Z[0][i], \
                        pts[4]->Z[0][i],   phiP->Z[0][i],   phiQ->Z[0][i],   phiR->Z[0][i]);
      P->Z[1][i] = _SET(pts[0]->Z[1][i], pts[1]->Z[1][i], pts[2]->Z[1][i], pts[3]->Z[1][i], \
                        pts[4]->Z[1][i],   phiP->Z[1][i],   phiQ->Z[1][i],   phiR->Z[1][i]);
      coeff[0][0][i] = VSET1(cf00[i]); coeff[0][1][i] = VSET1(cf01[i]);
      coeff[1][0][i] = VSET1(cf10[i]); coeff[1][1][i] = VSET1(cf11[i]);
    }

    eval_3_isog_8x1x1x1w(P, coeff);

    // extract results
    get_channel_8x1w(pts[0]->X[0], P->X[0], 7); get_channel_8x1w(pts[0]->X[1], P->X[1], 7);
    get_channel_8x1w(pts[0]->Z[0], P->Z[0], 7); get_channel_8x1w(pts[0]->Z[1], P->Z[1], 7);
    get_channel_8x1w(pts[1]->X[0], P->X[0], 6); get_channel_8x1w(pts[1]->X[1], P->X[1], 6);
    get_channel_8x1w(pts[1]->Z[0], P->Z[0], 6); get_channel_8x1w(pts[1]->Z[1], P->Z[1], 6);
    get_channel_8x1w(pts[2]->X[0], P->X[0], 5); get_channel_8x1w(pts[2]->X[1], P->X[1], 5);
    get_channel_8x1w(pts[2]->Z[0], P->Z[0], 5); get_channel_8x1w(pts[2]->Z[1], P->Z[1], 5);
    get_channel_8x1w(pts[3]->X[0], P->X[0], 4); get_channel_8x1w(pts[3]->X[1], P->X[1], 4);
    get_channel_8x1w(pts[3]->Z[0], P->Z[0], 4); get_channel_8x1w(pts[3]->Z[1], P->Z[1], 4);
    get_channel_8x1w(phiP->X[0], P->X[0], 2);   get_channel_8x1w(phiP->X[1], P->X[1], 2);
    get_channel_8x1w(phiP->Z[0], P->Z[0], 2);   get_channel_8x1w(phiP->Z[1], P->Z[1], 2);
    get_channel_8x1w(phiQ->X[0], P->X[0], 1);   get_channel_8x1w(phiQ->X[1], P->X[1], 1);
    get_channel_8x1w(phiQ->Z[0], P->Z[0], 1);   get_channel_8x1w(phiQ->Z[1], P->Z[1], 1);
    get_channel_8x1w(phiR->X[0], P->X[0], 0);   get_channel_8x1w(phiR->X[1], P->X[1], 0);
    get_channel_8x1w(phiR->Z[0], P->Z[0], 0);   get_channel_8x1w(phiR->Z[1], P->Z[1], 0);

    // depends on public info
    if (num >= 8) {
      get_channel_8x1w(pts[4]->X[0], P->X[0], 3); get_channel_8x1w(pts[4]->X[1], P->X[1], 3);
      get_channel_8x1w(pts[4]->Z[0], P->Z[0], 3); get_channel_8x1w(pts[4]->Z[1], P->Z[1], 3);
    }

    // depends on public info
    if (num == 9) {
      // 1x2x2x2w < pts[5] >
      vgelm_t _P, _coeff0_1;

      // initialize vectors 
      for (i = 0; i < VGWORDS; i++) {
        _P[i] = _SET(pts[5]->X[1][i+VGWORDS], pts[5]->X[1][i], pts[5]->X[0][i+VGWORDS], pts[5]->X[0][i],\
                     pts[5]->Z[1][i+VGWORDS], pts[5]->Z[1][i], pts[5]->Z[0][i+VGWORDS], pts[5]->Z[0][i]);
        _coeff0_1[i] = _SET(cf01[i+VGWORDS], cf01[i], cf00[i+VGWORDS], cf00[i], 
                            cf11[i+VGWORDS], cf11[i], cf10[i+VGWORDS], cf10[i] ); 
      }

      eval_3_isog_1x2x2x2w(_P, _coeff0_1);

      // extract results
      get_channel_4x2w(pts[5]->X[0], _P, 4); get_channel_4x2w(pts[5]->X[1], _P, 6);
      get_channel_4x2w(pts[5]->Z[0], _P, 0); get_channel_4x2w(pts[5]->Z[1], _P, 2);

      // carry propagation 
      carryp_1w(pts[5]->X[0]); carryp_1w(pts[5]->X[1]); 
      carryp_1w(pts[5]->Z[0]); carryp_1w(pts[5]->Z[1]); 
    }

    // depends on public info
    if (num == 10) {
      // 2x2x2x1w < pts[5] | pts[6] >
      vfelm_t _P, _coeff0_1;

      for (i = 0; i < VNWORDS; i++) {
        _P[i] = _SET(pts[5]->X[1][i], pts[5]->X[0][i], pts[5]->Z[1][i], pts[5]->Z[0][i], \
                     pts[6]->X[1][i], pts[6]->X[0][i], pts[6]->Z[1][i], pts[6]->Z[0][i]);
        _coeff0_1[i] = _SET(cf01[i], cf00[i], cf11[i], cf10[i], cf01[i], cf00[i], cf11[i], cf10[i]);
      }

      eval_3_isog_2x2x2x1w(_P, _coeff0_1);

      // extract results
      get_channel_8x1w(pts[5]->X[0], _P, 6); get_channel_8x1w(pts[5]->X[1], _P, 7);
      get_channel_8x1w(pts[5]->Z[0], _P, 4); get_channel_8x1w(pts[5]->Z[1], _P, 5);
      get_channel_8x1w(pts[6]->X[0], _P, 2); get_channel_8x1w(pts[6]->X[1], _P, 3);
      get_channel_8x1w(pts[6]->Z[0], _P, 0); get_channel_8x1w(pts[6]->Z[1], _P, 1);
    }

    // depends on public info
    if (num == 11 | num == 12 | num == 13) {
      // 4x2x1x1w < pts[5] | pts[6] | pts[7] | pts[8] >
      vf2elm_t _P, _coeff0_1;

      // initialize vectors 
      for (i = 0; i < VNWORDS; i++) {
        _P[0][i] = _SET(pts[5]->X[0][i], pts[5]->Z[0][i], pts[6]->X[0][i], pts[6]->Z[0][i], \
                        pts[7]->X[0][i], pts[7]->Z[0][i], pts[8]->X[0][i], pts[8]->Z[0][i]);
        _P[1][i] = _SET(pts[5]->X[1][i], pts[5]->Z[1][i], pts[6]->X[1][i], pts[6]->Z[1][i], \
                        pts[7]->X[1][i], pts[7]->Z[1][i], pts[8]->X[1][i], pts[8]->Z[1][i]);
        _coeff0_1[0][i] = _SET(cf00[i], cf10[i], cf00[i], cf10[i], cf00[i], cf10[i], cf00[i], cf10[i]);
        _coeff0_1[1][i] = _SET(cf01[i], cf11[i], cf01[i], cf11[i], cf01[i], cf11[i], cf01[i], cf11[i]);
      }

      eval_3_isog_4x2x1x1w(_P, _coeff0_1);

      // extract results
      get_channel_8x1w(pts[5]->X[0], _P[0], 7); get_channel_8x1w(pts[5]->X[1], _P[1], 7);
      get_channel_8x1w(pts[5]->Z[0], _P[0], 6); get_channel_8x1w(pts[5]->Z[1], _P[1], 6);
      get_channel_8x1w(pts[6]->X[0], _P[0], 5); get_channel_8x1w(pts[6]->X[1], _P[1], 5);
      get_channel_8x1w(pts[6]->Z[0], _P[0], 4); get_channel_8x1w(pts[6]->Z[1], _P[1], 4);
      get_channel_8x1w(pts[7]->X[0], _P[0], 3); get_channel_8x1w(pts[7]->X[1], _P[1], 3);
      get_channel_8x1w(pts[7]->Z[0], _P[0], 2); get_channel_8x1w(pts[7]->Z[1], _P[1], 2);

      // depends on public info
      if (num >= 12) {
        get_channel_8x1w(pts[8]->X[0], _P[0], 1); get_channel_8x1w(pts[8]->X[1], _P[1], 1);
        get_channel_8x1w(pts[8]->Z[0], _P[0], 0); get_channel_8x1w(pts[8]->Z[1], _P[1], 0);
      }

      // depends on public info
      if (num == 13) {
        // 1x2x2x2w < pts[9] >
        vgelm_t __P, __coeff0_1;

        // initialize vectors 
        for (i = 0; i < VGWORDS; i++) {
          __P[i] = _SET(pts[9]->X[1][i+VGWORDS], pts[9]->X[1][i], pts[9]->X[0][i+VGWORDS], pts[9]->X[0][i],\
                        pts[9]->Z[1][i+VGWORDS], pts[9]->Z[1][i], pts[9]->Z[0][i+VGWORDS], pts[9]->Z[0][i]);
          __coeff0_1[i] = _SET(cf01[i+VGWORDS], cf01[i], cf00[i+VGWORDS], cf00[i], 
                               cf11[i+VGWORDS], cf11[i], cf10[i+VGWORDS], cf10[i] ); 
        }

        eval_3_isog_1x2x2x2w(__P, __coeff0_1);

        // extract results
        get_channel_4x2w(pts[9]->X[0], __P, 4); get_channel_4x2w(pts[9]->X[1], __P, 6);
        get_channel_4x2w(pts[9]->Z[0], __P, 0); get_channel_4x2w(pts[9]->Z[1], __P, 2);

        // carry propagation 
        carryp_1w(pts[9]->X[0]); carryp_1w(pts[9]->X[1]); 
        carryp_1w(pts[9]->Z[0]); carryp_1w(pts[9]->Z[1]); 
      }
    }
  }
  // depends on public info
  else {
    puts("Error!");
  }
}

void eval_3_isog_parallel_ss(point_proj_r51_t *pts, const vgelm_t coeff1_0, const int num)
{
  felm_r51_t cf00, cf01, cf10, cf11;
  int i;

  // extract constants
  get_channel_4x2w(cf00, coeff1_0, 0); get_channel_4x2w(cf01, coeff1_0, 2);
  get_channel_4x2w(cf10, coeff1_0, 4); get_channel_4x2w(cf11, coeff1_0, 6);

  // carry propagation
  carryp_1w(cf00); carryp_1w(cf01);
  carryp_1w(cf10); carryp_1w(cf11);

  // depends on public info
  if (num == 1) {
    // 1x2x2x2w < pts[0] >
    vgelm_t _P, _coeff0_1;

    // initialize vectors 
    for (i = 0; i < VGWORDS; i++) {
      _P[i] = _SET(pts[0]->X[1][i+VGWORDS], pts[0]->X[1][i], pts[0]->X[0][i+VGWORDS], pts[0]->X[0][i],\
                   pts[0]->Z[1][i+VGWORDS], pts[0]->Z[1][i], pts[0]->Z[0][i+VGWORDS], pts[0]->Z[0][i]);
      _coeff0_1[i] = _SET(cf01[i+VGWORDS], cf01[i], cf00[i+VGWORDS], cf00[i], 
                          cf11[i+VGWORDS], cf11[i], cf10[i+VGWORDS], cf10[i] ); 
    }

    eval_3_isog_1x2x2x2w(_P, _coeff0_1);

    // extract results
    get_channel_4x2w(pts[0]->X[0], _P, 4); get_channel_4x2w(pts[0]->X[1], _P, 6);
    get_channel_4x2w(pts[0]->Z[0], _P, 0); get_channel_4x2w(pts[0]->Z[1], _P, 2);

    // carry propagation 
    carryp_1w(pts[0]->X[0]); carryp_1w(pts[0]->X[1]); 
    carryp_1w(pts[0]->Z[0]); carryp_1w(pts[0]->Z[1]); 
  }
  // depends on public info
  else if (num == 2) {
    // 2x2x2x1w < pts[0] | pts[1] >
    vfelm_t _P, _coeff0_1;

    for (i = 0; i < VNWORDS; i++) {
      _P[i] = _SET(pts[0]->X[1][i], pts[0]->X[0][i], pts[0]->Z[1][i], pts[0]->Z[0][i], \
                   pts[1]->X[1][i], pts[1]->X[0][i], pts[1]->Z[1][i], pts[1]->Z[0][i]);
      _coeff0_1[i] = _SET(cf01[i], cf00[i], cf11[i], cf10[i], cf01[i], cf00[i], cf11[i], cf10[i]);
    }

    eval_3_isog_2x2x2x1w(_P, _coeff0_1);

    // extract results
    get_channel_8x1w(pts[0]->X[0], _P, 6); get_channel_8x1w(pts[0]->X[1], _P, 7);
    get_channel_8x1w(pts[0]->Z[0], _P, 4); get_channel_8x1w(pts[0]->Z[1], _P, 5);
    get_channel_8x1w(pts[1]->X[0], _P, 2); get_channel_8x1w(pts[1]->X[1], _P, 3);
    get_channel_8x1w(pts[1]->Z[0], _P, 0); get_channel_8x1w(pts[1]->Z[1], _P, 1);
  }
  // depends on public info
  else if (num == 3 | num == 4 | num == 5 | num == 6) {
    // 4x2x1x1w < pts[0] | pts[1] | pts[2] | pts[3] >
    vf2elm_t P, coeff0_1;

    // initialize vectors 
    for (i = 0; i < VNWORDS; i++) {
      P[0][i] = _SET(pts[0]->X[0][i], pts[0]->Z[0][i], pts[1]->X[0][i], pts[1]->Z[0][i], \
                     pts[2]->X[0][i], pts[2]->Z[0][i], pts[3]->X[0][i], pts[3]->Z[0][i]);
      P[1][i] = _SET(pts[0]->X[1][i], pts[0]->Z[1][i], pts[1]->X[1][i], pts[1]->Z[1][i], \
                     pts[2]->X[1][i], pts[2]->Z[1][i], pts[3]->X[1][i], pts[3]->Z[1][i]);
      coeff0_1[0][i] = _SET(cf00[i], cf10[i], cf00[i], cf10[i], cf00[i], cf10[i], cf00[i], cf10[i]);
      coeff0_1[1][i] = _SET(cf01[i], cf11[i], cf01[i], cf11[i], cf01[i], cf11[i], cf01[i], cf11[i]);
    }

    eval_3_isog_4x2x1x1w(P, coeff0_1);

    // extract results
    get_channel_8x1w(pts[0]->X[0], P[0], 7); get_channel_8x1w(pts[0]->X[1], P[1], 7);
    get_channel_8x1w(pts[0]->Z[0], P[0], 6); get_channel_8x1w(pts[0]->Z[1], P[1], 6);
    get_channel_8x1w(pts[1]->X[0], P[0], 5); get_channel_8x1w(pts[1]->X[1], P[1], 5);
    get_channel_8x1w(pts[1]->Z[0], P[0], 4); get_channel_8x1w(pts[1]->Z[1], P[1], 4);
    get_channel_8x1w(pts[2]->X[0], P[0], 3); get_channel_8x1w(pts[2]->X[1], P[1], 3);
    get_channel_8x1w(pts[2]->Z[0], P[0], 2); get_channel_8x1w(pts[2]->Z[1], P[1], 2);

    // depends on public info
    if (num >= 4) {
      get_channel_8x1w(pts[3]->X[0], P[0], 1); get_channel_8x1w(pts[3]->X[1], P[1], 1);
      get_channel_8x1w(pts[3]->Z[0], P[0], 0); get_channel_8x1w(pts[3]->Z[1], P[1], 0);
    }

    // depends on public info
    if (num == 5) {
      // 1x2x2x2w < pts[4] >
      vgelm_t _P, _coeff0_1;

      // initialize vectors 
      for (i = 0; i < VGWORDS; i++) {
        _P[i] = _SET(pts[4]->X[1][i+VGWORDS], pts[4]->X[1][i], pts[4]->X[0][i+VGWORDS], pts[4]->X[0][i],\
                     pts[4]->Z[1][i+VGWORDS], pts[4]->Z[1][i], pts[4]->Z[0][i+VGWORDS], pts[4]->Z[0][i]);
        _coeff0_1[i] = _SET(cf01[i+VGWORDS], cf01[i], cf00[i+VGWORDS], cf00[i], 
                            cf11[i+VGWORDS], cf11[i], cf10[i+VGWORDS], cf10[i] ); 
      }

      eval_3_isog_1x2x2x2w(_P, _coeff0_1);

      // extract results
      get_channel_4x2w(pts[4]->X[0], _P, 4); get_channel_4x2w(pts[4]->X[1], _P, 6);
      get_channel_4x2w(pts[4]->Z[0], _P, 0); get_channel_4x2w(pts[4]->Z[1], _P, 2);

      // carry propagation 
      carryp_1w(pts[4]->X[0]); carryp_1w(pts[4]->X[1]); 
      carryp_1w(pts[4]->Z[0]); carryp_1w(pts[4]->Z[1]);  
    }

    // depends on public info
    if (num == 6) {
      // 2x2x2x1w < pts[4] | pts[5] >
      vfelm_t _P, _coeff0_1;

      for (i = 0; i < VNWORDS; i++) {
        _P[i] = _SET(pts[4]->X[1][i], pts[4]->X[0][i], pts[4]->Z[1][i], pts[4]->Z[0][i], \
                     pts[5]->X[1][i], pts[5]->X[0][i], pts[5]->Z[1][i], pts[5]->Z[0][i]);
        _coeff0_1[i] = _SET(cf01[i], cf00[i], cf11[i], cf10[i], cf01[i], cf00[i], cf11[i], cf10[i]);
      }

      eval_3_isog_2x2x2x1w(_P, _coeff0_1);

      // extract results
      get_channel_8x1w(pts[4]->X[0], _P, 6); get_channel_8x1w(pts[4]->X[1], _P, 7);
      get_channel_8x1w(pts[4]->Z[0], _P, 4); get_channel_8x1w(pts[4]->Z[1], _P, 5);
      get_channel_8x1w(pts[5]->X[0], _P, 2); get_channel_8x1w(pts[5]->X[1], _P, 3);
      get_channel_8x1w(pts[5]->Z[0], _P, 0); get_channel_8x1w(pts[5]->Z[1], _P, 1);
    }
  }
  // depends on public info
  else if (num == 7 | num == 8 | num == 9 | num == 10) {
    // 8x1x1x1w < pts[0] | pts[1] | pts[2] | pts[3] | pts[4] | pts[5] | pts[6] | pts[7] >
    vpoint_proj_t P;
    vf2elm_t coeff[2];

    // initialize vectors 
    for (i = 0; i < VNWORDS; i++) {
      P->X[0][i] = _SET(pts[0]->X[0][i], pts[1]->X[0][i], pts[2]->X[0][i], pts[3]->X[0][i], \
                        pts[4]->X[0][i], pts[5]->X[0][i], pts[6]->X[0][i], pts[7]->X[0][i]);
      P->X[1][i] = _SET(pts[0]->X[1][i], pts[1]->X[1][i], pts[2]->X[1][i], pts[3]->X[1][i], \
                        pts[4]->X[1][i], pts[5]->X[1][i], pts[6]->X[1][i], pts[7]->X[1][i]);
      P->Z[0][i] = _SET(pts[0]->Z[0][i], pts[1]->Z[0][i], pts[2]->Z[0][i], pts[3]->Z[0][i], \
                        pts[4]->Z[0][i], pts[5]->Z[0][i], pts[6]->Z[0][i], pts[7]->Z[0][i]);
      P->Z[1][i] = _SET(pts[0]->Z[1][i], pts[1]->Z[1][i], pts[2]->Z[1][i], pts[3]->Z[1][i], \
                        pts[4]->Z[1][i], pts[5]->Z[1][i], pts[6]->Z[1][i], pts[7]->Z[1][i]);
      coeff[0][0][i] = VSET1(cf00[i]); coeff[0][1][i] = VSET1(cf01[i]);
      coeff[1][0][i] = VSET1(cf10[i]); coeff[1][1][i] = VSET1(cf11[i]);
    }

    eval_3_isog_8x1x1x1w(P, coeff);

    // extract results
    get_channel_8x1w(pts[0]->X[0], P->X[0], 7); get_channel_8x1w(pts[0]->X[1], P->X[1], 7);
    get_channel_8x1w(pts[0]->Z[0], P->Z[0], 7); get_channel_8x1w(pts[0]->Z[1], P->Z[1], 7);
    get_channel_8x1w(pts[1]->X[0], P->X[0], 6); get_channel_8x1w(pts[1]->X[1], P->X[1], 6);
    get_channel_8x1w(pts[1]->Z[0], P->Z[0], 6); get_channel_8x1w(pts[1]->Z[1], P->Z[1], 6);
    get_channel_8x1w(pts[2]->X[0], P->X[0], 5); get_channel_8x1w(pts[2]->X[1], P->X[1], 5);
    get_channel_8x1w(pts[2]->Z[0], P->Z[0], 5); get_channel_8x1w(pts[2]->Z[1], P->Z[1], 5);
    get_channel_8x1w(pts[3]->X[0], P->X[0], 4); get_channel_8x1w(pts[3]->X[1], P->X[1], 4);
    get_channel_8x1w(pts[3]->Z[0], P->Z[0], 4); get_channel_8x1w(pts[3]->Z[1], P->Z[1], 4);
    get_channel_8x1w(pts[4]->X[0], P->X[0], 3); get_channel_8x1w(pts[4]->X[1], P->X[1], 3);
    get_channel_8x1w(pts[4]->Z[0], P->Z[0], 3); get_channel_8x1w(pts[4]->Z[1], P->Z[1], 3);
    get_channel_8x1w(pts[5]->X[0], P->X[0], 2); get_channel_8x1w(pts[5]->X[1], P->X[1], 2);
    get_channel_8x1w(pts[5]->Z[0], P->Z[0], 2); get_channel_8x1w(pts[5]->Z[1], P->Z[1], 2);
    get_channel_8x1w(pts[6]->X[0], P->X[0], 1); get_channel_8x1w(pts[6]->X[1], P->X[1], 1);
    get_channel_8x1w(pts[6]->Z[0], P->Z[0], 1); get_channel_8x1w(pts[6]->Z[1], P->Z[1], 1);

    // depends on public info
    if (num >= 8) {
      get_channel_8x1w(pts[7]->X[0], P->X[0], 0); get_channel_8x1w(pts[7]->X[1], P->X[1], 0);
      get_channel_8x1w(pts[7]->Z[0], P->Z[0], 0); get_channel_8x1w(pts[7]->Z[1], P->Z[1], 0);
    }

    // depends on public info
    if (num == 9) {
      // 1x2x2x2w < pts[8] >
      vgelm_t _P, _coeff0_1;

      // initialize vectors 
      for (i = 0; i < VGWORDS; i++) {
        _P[i] = _SET(pts[8]->X[1][i+VGWORDS], pts[8]->X[1][i], pts[8]->X[0][i+VGWORDS], pts[8]->X[0][i],\
                     pts[8]->Z[1][i+VGWORDS], pts[8]->Z[1][i], pts[8]->Z[0][i+VGWORDS], pts[8]->Z[0][i]);
        _coeff0_1[i] = _SET(cf01[i+VGWORDS], cf01[i], cf00[i+VGWORDS], cf00[i], 
                            cf11[i+VGWORDS], cf11[i], cf10[i+VGWORDS], cf10[i] ); 
      }

      eval_3_isog_1x2x2x2w(_P, _coeff0_1);

      // extract results
      get_channel_4x2w(pts[8]->X[0], _P, 4); get_channel_4x2w(pts[8]->X[1], _P, 6);
      get_channel_4x2w(pts[8]->Z[0], _P, 0); get_channel_4x2w(pts[8]->Z[1], _P, 2);

      // carry propagation 
      carryp_1w(pts[8]->X[0]); carryp_1w(pts[8]->X[1]); 
      carryp_1w(pts[8]->Z[0]); carryp_1w(pts[8]->Z[1]);
    }

    // depends on public info
    if (num == 10) {
      // 2x2x2x1w < pts[8] | pts[9] >
      vfelm_t _P, _coeff0_1;

      for (i = 0; i < VNWORDS; i++) {
        _P[i] = _SET(pts[8]->X[1][i], pts[8]->X[0][i], pts[8]->Z[1][i], pts[8]->Z[0][i], \
                     pts[9]->X[1][i], pts[9]->X[0][i], pts[9]->Z[1][i], pts[9]->Z[0][i]);
        _coeff0_1[i] = _SET(cf01[i], cf00[i], cf11[i], cf10[i], cf01[i], cf00[i], cf11[i], cf10[i]);
      }

      eval_3_isog_2x2x2x1w(_P, _coeff0_1);

      // extract results
      get_channel_8x1w(pts[8]->X[0], _P, 6); get_channel_8x1w(pts[8]->X[1], _P, 7);
      get_channel_8x1w(pts[8]->Z[0], _P, 4); get_channel_8x1w(pts[8]->Z[1], _P, 5);
      get_channel_8x1w(pts[9]->X[0], _P, 2); get_channel_8x1w(pts[9]->X[1], _P, 3);
      get_channel_8x1w(pts[9]->Z[0], _P, 0); get_channel_8x1w(pts[9]->Z[1], _P, 1);
    }

  }
  // depends on public info
  else {
    puts("Error!");    
  }
}

void inv_3_way(f2elm_t z1, f2elm_t z2, f2elm_t z3)
{ // 3-way simultaneous inversion
  // Input:  z1,z2,z3
  // Output: 1/z1,1/z2,1/z3 (override inputs).
    f2elm_t t0, t1, t2, t3;

    fp2mul_mont(z1, z2, t0);                      // t0 = z1*z2
    fp2mul_mont(z3, t0, t1);                      // t1 = z1*z2*z3
    fp2inv_mont(t1);                              // t1 = 1/(z1*z2*z3)
    fp2mul_mont(z3, t1, t2);                      // t2 = 1/(z1*z2) 
    fp2mul_mont(t2, z2, t3);                      // t3 = 1/z1
    fp2mul_mont(t2, z1, z2);                      // z2 = 1/z2
    fp2mul_mont(t0, t1, z3);                      // z3 = 1/z3
    fp2copy(t3, z1);                              // z1 = 1/z1
}

void get_A(const f2elm_t xP, const f2elm_t xQ, const f2elm_t xR, f2elm_t A)
{ // Given the x-coordinates of P, Q, and R, returns the value A corresponding to the Montgomery curve E_A: y^2=x^3+A*x^2+x such that R=Q-P on E_A.
  // Input:  the x-coordinates xP, xQ, and xR of the points P, Q and R.
  // Output: the coefficient A corresponding to the curve E_A: y^2=x^3+A*x^2+x.
    f2elm_t t0, t1, one = {0};
    
    fpcopy((digit_t*)&Montgomery_one, one[0]);
    fp2add(xP, xQ, t1);                           // t1 = xP+xQ
    fp2mul_mont(xP, xQ, t0);                      // t0 = xP*xQ
    fp2mul_mont(xR, t1, A);                       // A = xR*t1
    fp2add(t0, A, A);                             // A = A+t0
    fp2mul_mont(t0, xR, t0);                      // t0 = t0*xR
    fp2sub(A, one, A);                            // A = A-1
    fp2add(t0, t0, t0);                           // t0 = t0+t0
    fp2add(t1, xR, t1);                           // t1 = t1+xR
    fp2add(t0, t0, t0);                           // t0 = t0+t0
    fp2sqr_mont(A, A);                            // A = A^2
    fp2inv_mont(t0);                              // t0 = 1/t0
    fp2mul_mont(A, t0, A);                        // A = A*t0
    fp2sub(A, t1, A);                             // Afinal = A-t1
}

void j_inv(const f2elm_t A, const f2elm_t C, f2elm_t jinv)
{ // Computes the j-invariant of a Montgomery curve with projective constant.
  // Input: A,C in GF(p^2).
  // Output: j=256*(A^2-3*C^2)^3/(C^4*(A^2-4*C^2)), which is the j-invariant of the Montgomery curve B*y^2=x^3+(A/C)*x^2+x or (equivalently) j-invariant of B'*y^2=C*x^3+A*x^2+C*x.
    f2elm_t t0, t1;
    
    fp2sqr_mont(A, jinv);                           // jinv = A^2        
    fp2sqr_mont(C, t1);                             // t1 = C^2
    fp2add(t1, t1, t0);                             // t0 = t1+t1
    fp2sub(jinv, t0, t0);                           // t0 = jinv-t0
    fp2sub(t0, t1, t0);                             // t0 = t0-t1
    fp2sub(t0, t1, jinv);                           // jinv = t0-t1
    fp2sqr_mont(t1, t1);                            // t1 = t1^2
    fp2mul_mont(jinv, t1, jinv);                    // jinv = jinv*t1
    fp2add(t0, t0, t0);                             // t0 = t0+t0
    fp2add(t0, t0, t0);                             // t0 = t0+t0
    fp2sqr_mont(t0, t1);                            // t1 = t0^2
    fp2mul_mont(t0, t1, t0);                        // t0 = t0*t1
    fp2add(t0, t0, t0);                             // t0 = t0+t0
    fp2add(t0, t0, t0);                             // t0 = t0+t0
    fp2inv_mont(jinv);                              // jinv = 1/jinv 
    fp2mul_mont(jinv, t0, jinv);                    // jinv = t0*jinv
}
