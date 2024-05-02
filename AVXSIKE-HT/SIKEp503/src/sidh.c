/**
 *******************************************************************************
 * @version 0.1
 * @date 2021-12-16
 * @copyright Copyright © 2021 by University of Luxembourg
 * @author Hao Cheng
 *******************************************************************************
 */

#include "sidh.h"

static void init_basis(const uint64_t *gen, f2elm_t XP, f2elm_t XQ, f2elm_t XR)
{
  int i;

  for (i = 0; i < NWORDS; i++) {
    XP[0][i] = VSET1(gen[i          ]);
    XP[1][i] = VSET1(gen[i + NWORDS]);
    XQ[0][i] = VSET1(gen[i+2*NWORDS]);
    XQ[1][i] = VSET1(gen[i+3*NWORDS]);
    XR[0][i] = VSET1(gen[i+4*NWORDS]);
    XR[1][i] = VSET1(gen[i+5*NWORDS]);
  }
}

// generation of Bob's secret key
void random_mod_order_B(__m512i *random_digits, uint8_t *sk)
{
  uint8_t s0[SK_B_VECTS*8] = { 0 }, s1[SK_B_VECTS*8] = { 0 };
  uint8_t s2[SK_B_VECTS*8] = { 0 }, s3[SK_B_VECTS*8] = { 0 };
  uint8_t s4[SK_B_VECTS*8] = { 0 }, s5[SK_B_VECTS*8] = { 0 };
  uint8_t s6[SK_B_VECTS*8] = { 0 }, s7[SK_B_VECTS*8] = { 0 };
  uint64_t *s64_0 = (uint64_t *)s0, *s64_1 = (uint64_t *)s1; 
  uint64_t *s64_2 = (uint64_t *)s2, *s64_3 = (uint64_t *)s3; 
  uint64_t *s64_4 = (uint64_t *)s4, *s64_5 = (uint64_t *)s5; 
  uint64_t *s64_6 = (uint64_t *)s6, *s64_7 = (uint64_t *)s7; 
  int i;

  randombytes(s0, SECRETKEY_B_BYTES); s0[SECRETKEY_B_BYTES-1] &= MASK_BOB;
  randombytes(s1, SECRETKEY_B_BYTES); s1[SECRETKEY_B_BYTES-1] &= MASK_BOB;
  randombytes(s2, SECRETKEY_B_BYTES); s2[SECRETKEY_B_BYTES-1] &= MASK_BOB;
  randombytes(s3, SECRETKEY_B_BYTES); s3[SECRETKEY_B_BYTES-1] &= MASK_BOB;
  randombytes(s4, SECRETKEY_B_BYTES); s4[SECRETKEY_B_BYTES-1] &= MASK_BOB;
  randombytes(s5, SECRETKEY_B_BYTES); s5[SECRETKEY_B_BYTES-1] &= MASK_BOB;
  randombytes(s6, SECRETKEY_B_BYTES); s6[SECRETKEY_B_BYTES-1] &= MASK_BOB;
  randombytes(s7, SECRETKEY_B_BYTES); s7[SECRETKEY_B_BYTES-1] &= MASK_BOB;

  // form strings in vectors
  for (i = 0; i < SK_B_VECTS; i++) 
    random_digits[i] = set_vector(s64_7[i], s64_6[i], s64_5[i], s64_4[i], s64_3[i], s64_2[i], s64_1[i], s64_0[i]);

  for (i = 0; i < SECRETKEY_B_BYTES; i++) {
    sk[                        MSG_BYTES+i] = s0[i];
    sk[  CRYPTO_SECRETKEYBYTES+MSG_BYTES+i] = s1[i];
    sk[2*CRYPTO_SECRETKEYBYTES+MSG_BYTES+i] = s2[i];
    sk[3*CRYPTO_SECRETKEYBYTES+MSG_BYTES+i] = s3[i];
    sk[4*CRYPTO_SECRETKEYBYTES+MSG_BYTES+i] = s4[i];
    sk[5*CRYPTO_SECRETKEYBYTES+MSG_BYTES+i] = s5[i];
    sk[6*CRYPTO_SECRETKEYBYTES+MSG_BYTES+i] = s6[i];
    sk[7*CRYPTO_SECRETKEYBYTES+MSG_BYTES+i] = s7[i];
  }
}

void EphemeralKeyGeneration_A(const __m512i *PrivateKeyA, __m512i *PublicKeyA)
{
  point_proj_t R, phiP, phiQ, phiR, pts[MAX_INT_POINTS_ALICE];
  f2elm_t XPA, XQA, XRA, coeff[3], A24plus, C24, A;
  __m512i SecretKeyA[SK_A_VECTS];
  int i, row, m, index = 0, pts_index[MAX_INT_POINTS_ALICE], npts = 0, ii = 0;

  // initialize basis points 
  init_basis(A_gen, XPA, XQA, XRA);
  init_basis(B_gen, phiP->X, phiQ->X, phiR->X);
  for (i = 0; i < NWORDS; i++) {
    phiP->Z[0][i] = phiQ->Z[0][i] = phiR->Z[0][i] = VSET1(mont_R[i]);
    phiP->Z[1][i] = phiQ->Z[1][i] = phiR->Z[1][i] = VZERO;
  } 

  // initialize constants: A24plus = A+2C, C24 = 4C, where A = 6, C = 1
  for (i = 0; i < NWORDS; i++) {
    A24plus[0][i] = VSET1(mont_R[i]);
    A24plus[1][i] = VZERO;
  }

  mp2_add(A24plus, A24plus, A24plus);   // A24plus = 2
  mp2_add(C24, A24plus, A24plus);       // C24 = 4
  mp2_add(A, A24plus, C24);             // A = 6
  mp2_add(A24plus, C24, C24);           // A24plus = 8

  // retrieve kernel point
  for (i = 0; i < SK_A_VECTS; i++) SecretKeyA[i] = PrivateKeyA[i];
  LADDER3PT(XPA, XQA, XRA, SecretKeyA, ALICE, R, A);

  // traverse tree
  index = 0;
  // row = 1;
  for (row = 1; row < MAX_Alice; row++) {
    while (index < MAX_Alice-row) {
      pointcopy(pts[npts], R);
      pts_index[npts++] = index;
      m = strat_Alice[ii++];
      xDBLe(R, R, A24plus, C24, 2*m);
      index += m;
    }

    get_4_isog(R, A24plus, C24, coeff);

    for (i = 0; i < npts; i++) eval_4_isog(pts[i], coeff);
    eval_4_isog(phiP, coeff);
    eval_4_isog(phiQ, coeff);
    eval_4_isog(phiR, coeff);   

    pointcopy(R, pts[npts-1]);
    index = pts_index[npts-1];
    npts -= 1;
  }

  get_4_isog(R, A24plus, C24, coeff);
  eval_4_isog(phiP, coeff);
  eval_4_isog(phiQ, coeff);
  eval_4_isog(phiR, coeff);

  inv_3_way(phiP->Z, phiQ->Z, phiR->Z);
  fp2mul_mont(phiP->X, phiP->X, phiP->Z);
  fp2mul_mont(phiQ->X, phiQ->X, phiQ->Z);
  fp2mul_mont(phiR->X, phiR->X, phiR->Z);

  // form public key 
  from_fp2mont(phiP->X, phiP->X);
  from_fp2mont(phiQ->X, phiQ->X);
  from_fp2mont(phiR->X, phiR->X);
  for (i = 0; i < NWORDS; i++) {
    PublicKeyA[i]            = phiP->X[0][i];
    PublicKeyA[i +   NWORDS] = phiP->X[1][i];
    PublicKeyA[i + 2*NWORDS] = phiQ->X[0][i];
    PublicKeyA[i + 3*NWORDS] = phiQ->X[1][i];
    PublicKeyA[i + 4*NWORDS] = phiR->X[0][i];
    PublicKeyA[i + 5*NWORDS] = phiR->X[1][i];
  }
}

void EphemeralKeyGeneration_B(const __m512i *PrivateKeyB, __m512i *PublicKeyB)
{
  point_proj_t R, phiP, phiQ, phiR, pts[MAX_INT_POINTS_BOB];
  f2elm_t XPB, XQB, XRB, coeff[3], A24plus, A24minus, A;
  __m512i SecretKeyB[SK_B_VECTS];
  int i, row, m, index = 0, pts_index[MAX_INT_POINTS_BOB], npts = 0, ii = 0;

  // initialize basis points 
  init_basis(B_gen, XPB, XQB, XRB);
  init_basis(A_gen, phiP->X, phiQ->X, phiR->X);
  for (i = 0; i < NWORDS; i++) {
    phiP->Z[0][i] = phiQ->Z[0][i] = phiR->Z[0][i] = VSET1(mont_R[i]);
    phiP->Z[1][i] = phiQ->Z[1][i] = phiR->Z[1][i] = VZERO;
  } 

  // initialize constants: A24plus = A+2C, A24minus = A-2C, where A = 6, C = 1
  for (i = 0; i < NWORDS; i++) {
    A24plus[0][i] = VSET1(mont_R[i]);
    A24plus[1][i] = VZERO;
  }

  mp2_add(A24plus, A24plus, A24plus);   // A24plus = 2
  mp2_add(A24minus, A24plus, A24plus);  // A24minus = 4
  mp2_add(A, A24plus, A24minus);        // A = 6
  mp2_add(A24plus, A24minus, A24minus); // A24plus = 8

  // retrieve kernel point
  for (i = 0; i < SK_B_VECTS; i++) SecretKeyB[i] = PrivateKeyB[i];
  LADDER3PT(XPB, XQB, XRB, SecretKeyB, BOB, R, A);

  // traverse tree
  index = 0;
  for (row = 1; row < MAX_Bob; row++) {
    while (index < MAX_Bob-row) {
      pointcopy(pts[npts], R);
      pts_index[npts++] = index;
      m = strat_Bob[ii++];
      xTPLe(R, R, A24minus, A24plus, m);
      index += m;
    }
    get_3_isog(R, A24minus, A24plus, coeff);

    for (i = 0; i < npts; i++) eval_3_isog(pts[i], coeff);
    eval_3_isog(phiP, coeff);
    eval_3_isog(phiQ, coeff);
    eval_3_isog(phiR, coeff);

    pointcopy(R, pts[npts-1]);
    index = pts_index[npts-1];
    npts -= 1;
  }

  get_3_isog(R, A24minus, A24plus, coeff);
  eval_3_isog(phiP, coeff);
  eval_3_isog(phiQ, coeff);
  eval_3_isog(phiR, coeff);

  inv_3_way(phiP->Z, phiQ->Z, phiR->Z);
  fp2mul_mont(phiP->X, phiP->X, phiP->Z);
  fp2mul_mont(phiQ->X, phiQ->X, phiQ->Z);
  fp2mul_mont(phiR->X, phiR->X, phiR->Z);

  // form public key 
  from_fp2mont(phiP->X, phiP->X);
  from_fp2mont(phiQ->X, phiQ->X);
  from_fp2mont(phiR->X, phiR->X);
  for (i = 0; i < NWORDS; i++) {
    PublicKeyB[i]            = phiP->X[0][i];
    PublicKeyB[i +   NWORDS] = phiP->X[1][i];
    PublicKeyB[i + 2*NWORDS] = phiQ->X[0][i];
    PublicKeyB[i + 3*NWORDS] = phiQ->X[1][i];
    PublicKeyB[i + 4*NWORDS] = phiR->X[0][i];
    PublicKeyB[i + 5*NWORDS] = phiR->X[1][i];
  }
}

void EphemeralSecretAgreement_A(const __m512i *PrivateKeyA, const __m512i *PublicKeyB, __m512i *SharedSecretA)
{
  point_proj_t R, pts[MAX_INT_POINTS_ALICE];
  f2elm_t coeff[3], PKB[3], jinv;
  f2elm_t A24plus, C24, A;
  __m512i SecretKeyA[SK_A_VECTS] = {0};
  int i, row, m, index = 0, pts_index[MAX_INT_POINTS_ALICE], npts = 0, ii = 0;

  // initialize images of Bob's basis
  for (i = 0; i < NWORDS; i++) {
    PKB[0][0][i] = PublicKeyB[i           ];
    PKB[0][1][i] = PublicKeyB[i +   NWORDS];
    PKB[1][0][i] = PublicKeyB[i + 2*NWORDS];
    PKB[1][1][i] = PublicKeyB[i + 3*NWORDS];
    PKB[2][0][i] = PublicKeyB[i + 4*NWORDS];
    PKB[2][1][i] = PublicKeyB[i + 5*NWORDS];
  }
  to_fp2mont(PKB[0], PKB[0]);
  to_fp2mont(PKB[1], PKB[1]);
  to_fp2mont(PKB[2], PKB[2]);

  // initialize constants: A24plus = A+2C, C24 = 4C, where C=1
  get_A(PKB[0], PKB[1], PKB[2], A);
  for (i = 0; i < NWORDS; i++) {
    A24plus[0][i] = A24plus[1][i] = C24[1][i] = VZERO;
    C24[0][i] = VSET1(mont_R[i]);
  }
  mp2_add(C24, C24, C24);               // C24 = 2 
  mp2_add(A24plus, A, C24);             // A24plus = A + 2
  mp2_add(C24, C24, C24);               // C24 = 4 
  
  // retrieve kernel point
  for (i = 0; i < SK_A_VECTS; i++) SecretKeyA[i] = PrivateKeyA[i];
  LADDER3PT(PKB[0], PKB[1], PKB[2], SecretKeyA, ALICE, R, A); 

  // traverse tree
  index = 0;
  for (row = 1; row < MAX_Alice; row++) {
    while (index < MAX_Alice-row) {
      pointcopy(pts[npts], R);
      pts_index[npts++] = index;
      m = strat_Alice[ii++];
      xDBLe(R, R, A24plus, C24, 2*m);
      index += m;
    }
    get_4_isog(R, A24plus, C24, coeff);

    for (i = 0; i < npts; i++) eval_4_isog(pts[i], coeff);

    pointcopy(R, pts[npts-1]);
    index = pts_index[npts-1];
    npts -= 1;
  }

  get_4_isog(R, A24plus, C24, coeff); 
  mp2_add(A24plus, A24plus, A24plus);
  fp2sub(A24plus, A24plus, C24);
  fp2add(A24plus, A24plus, A24plus);
  j_inv(A24plus, C24, jinv);

  // encode ss
  from_fp2mont(jinv, jinv);
  for (i = 0; i < NWORDS; i++) {
    SharedSecretA[i]        = jinv[0][i];
    SharedSecretA[i+NWORDS] = jinv[1][i];
  }
}

void EphemeralSecretAgreement_B(const __m512i *PrivateKeyB, const __m512i *PublicKeyA, __m512i *SharedSecretB)
{
  point_proj_t R, pts[MAX_INT_POINTS_BOB];
  f2elm_t coeff[3], PKB[3], jinv;
  f2elm_t A24plus, A24minus, A;
  __m512i SecretKeyB[SK_B_VECTS];
  int i, row, m, index = 0, pts_index[MAX_INT_POINTS_ALICE], npts = 0, ii = 0;

  // initialize images of Alice's basis
  for (i = 0; i < NWORDS; i++) {
    PKB[0][0][i] = PublicKeyA[i           ];
    PKB[0][1][i] = PublicKeyA[i +   NWORDS];
    PKB[1][0][i] = PublicKeyA[i + 2*NWORDS];
    PKB[1][1][i] = PublicKeyA[i + 3*NWORDS];
    PKB[2][0][i] = PublicKeyA[i + 4*NWORDS];
    PKB[2][1][i] = PublicKeyA[i + 5*NWORDS];
  }
  to_fp2mont(PKB[0], PKB[0]);
  to_fp2mont(PKB[1], PKB[1]);
  to_fp2mont(PKB[2], PKB[2]);

  // initialize constants: A24plus = A+2C, A24minus = A-2C, where C=1
  get_A(PKB[0], PKB[1], PKB[2], A);
  for (i = 0; i < NWORDS; i++) {
    A24plus[0][i] = A24plus[1][i] = A24minus[1][i] = VZERO;
    A24minus[0][i] = VSET1(mont_R[i]);
  }
  mp2_add(A24minus, A24minus, A24minus);// A24minus = 2 
  mp2_add(A24plus, A, A24minus);        // A24plus = A + 2
  mp2_sub_p2(A24minus, A, A24minus);       // A24minus = A - 2 
  
  // retrieve kernel point
  for (i = 0; i < SK_B_VECTS; i++) SecretKeyB[i] = PrivateKeyB[i];
  LADDER3PT(PKB[0], PKB[1], PKB[2], SecretKeyB, BOB, R, A);  

  // traverse tree
  index = 0;
  for (row = 1; row < MAX_Bob; row++) {
    while (index < MAX_Bob-row) {
      pointcopy(pts[npts], R);
      pts_index[npts++] = index;
      m = strat_Bob[ii++];
      xTPLe(R, R, A24minus, A24plus, m);
      index += m;
    }
    get_3_isog(R, A24minus, A24plus, coeff);

    for (i = 0; i < npts; i++) eval_3_isog(pts[i], coeff);

    pointcopy(R, pts[npts-1]);
    index = pts_index[npts-1];
    npts -= 1;
  }

  get_3_isog(R, A24minus, A24plus, coeff);
  fp2add(A, A24plus, A24minus); 
  fp2add(A, A, A); 
  fp2sub(A24plus, A24plus, A24minus);      
  j_inv(A, A24plus, jinv);

  // encode ss
  from_fp2mont(jinv, jinv);
  for (i = 0; i < NWORDS; i++) {
    SharedSecretB[i]        = jinv[0][i];
    SharedSecretB[i+NWORDS] = jinv[1][i];
  }
}
