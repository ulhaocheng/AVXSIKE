/**
 *******************************************************************************
 * @version 0.1
 * @date 2021-12-16
 * @copyright Copyright © 2021 by University of Luxembourg
 * @author Hao Cheng
 *******************************************************************************
 */

#include "sidh.h"
#include "random.h"

static void init_basis_1w(uint64_t *gen, f2elm_r51_t XP, f2elm_r51_t XQ, f2elm_r51_t XR)
{
  int i;

  for (i = 0; i < VNWORDS; i++) {
    XP[0][i] = gen[i];
    XP[1][i] = gen[i+VNWORDS];
    XQ[0][i] = gen[i+2*VNWORDS];
    XQ[1][i] = gen[i+3*VNWORDS];
    XR[0][i] = gen[i+4*VNWORDS];
    XR[1][i] = gen[i+5*VNWORDS];
  }
}

void random_mod_order_B(unsigned char* random_digits)
{  // Generation of Bob's secret key  
   // Outputs random value in [0, 2^Floor(Log(2, oB)) - 1]

    randombytes(random_digits, SECRETKEY_B_BYTES);
    random_digits[SECRETKEY_B_BYTES-1] &= MASK_BOB;     // Masking last byte 
}

// The vectorized Alice's ephemeral public key generation
// PKA in radix-2^64
void EphemeralKeyGeneration_A(const unsigned char* PrivateKeyA, unsigned char* PublicKeyA)
{
  point_proj_r51_t R, phiP = { 0 }, phiQ = { 0 }, phiR = { 0 }, pts[MAX_INT_POINTS_ALICE];
  f2elm_r51_t XPA, XQA, XRA, A24plus = { 0 }, C24 = { 0 }, A = { 0 };
  digit_t SecretKeyA[NWORDS_ORDER] = { 0 };
  int i, row, m, index = 0, pts_index[MAX_INT_POINTS_ALICE], npts = 0, ii = 0;
  vgelm_t vR, C24_A24plus, coeff__0, coeff2_1;
  point_proj_t phiP64 = { 0 }, phiQ64 = { 0 }, phiR64 = { 0 };

  // Initialize basis points
  init_basis_1w((digit_t*)vA_gen, XPA, XQA, XRA); 
  init_basis_1w((digit_t*)vB_gen, phiP->X, phiQ->X, phiR->X);
  fpcopy_1w(phiP->Z[0], vmont_R);
  fpcopy_1w(phiQ->Z[0], vmont_R);
  fpcopy_1w(phiR->Z[0], vmont_R);

  // Initialize constants: A24plus = A+2C, C24 = 4C, where A=6, C=1
  fpcopy_1w(A24plus[0], vmont_R);       // A24plus = 1
  mp2_add_1w(A24plus, A24plus, A24plus);// A24plus = 2
  mp2_add_1w(C24, A24plus, A24plus);    // C24 = 4
  mp2_add_1w(A, A24plus, C24);          // A = 6 
  mp2_add_1w(A24plus, C24, C24);        // A24plus = 8

  // Retrieve kernel point
  decode_to_digits(PrivateKeyA, SecretKeyA, SECRETKEY_A_BYTES, NWORDS_ORDER);
  // vectorized Montgomery ladder (ladder step is 1x4x2x1w)
  // the ouput R is in radix-2^51
  LADDER3PT_1x4x2x1w(XPA, XQA, XRA, SecretKeyA, ALICE, R, A); 

  // intialize C24_A24plus : < C1' | C1 | C0' | C0 | A1' | A1 | A0' | A0 >
  for (i = 0; i < VGWORDS-1; i++) 
    C24_A24plus[i] = _SET(C24[1][i+VGWORDS], C24[1][i], \
                          C24[0][i+VGWORDS], C24[0][i], \
                          A24plus[1][i+VGWORDS], A24plus[1][i], \
                          A24plus[0][i+VGWORDS], A24plus[0][i]);
  C24_A24plus[i] = _SET(0,     C24[1][i], 0, C24[0][i], \
                        0, A24plus[1][i], 0, A24plus[0][i]);
  
  // Traverse tree
  index = 0;
  for (row = 1; row < MAX_Alice; row++) {
    // update vR : < XR1' | XR1 | XR0' | XR0 | ZR1' | ZR1 | ZR0' | ZR0 >
    for (i = 0; i < VGWORDS-1; i++) 
      vR[i] = _SET(R->X[1][i+VGWORDS], R->X[1][i], R->X[0][i+VGWORDS], R->X[0][i], \
                  R->Z[1][i+VGWORDS], R->Z[1][i], R->Z[0][i+VGWORDS], R->Z[0][i]);
    vR[i] = _SET(0, R->X[1][i], 0, R->X[0][i], \
                 0, R->Z[1][i], 0, R->Z[0][i]);
    
    while (index < MAX_Alice-row) {
      pointcopy_1w(pts[npts], R);
      pts_index[npts++] = index;
      m = strat_Alice[ii++];
      xDBLe_1x2x2x2w(vR, R, C24_A24plus, 2*m);  // both vR and R get updated here
      index += m;
    }
    
    get_4_isog_1x2x2x2w(vR, C24_A24plus, coeff__0, coeff2_1);

    // eval_4_isog could be parallelized in 1/2/4/8-way, depending on npts. 
    eval_4_isog_parallel_kg(pts, phiP, phiQ, phiR, coeff__0, coeff2_1, npts+3);

    pointcopy_1w(R, pts[npts-1]);       // update R 
    index = pts_index[npts-1];
    npts -= 1;
  }

  // update vR : < XR1' | XR1 | XR0' | XR0 | ZR1' | ZR1 | ZR0' | ZR0 >
  for (i = 0; i < VGWORDS-1; i++) 
    vR[i] = _SET(R->X[1][i+VGWORDS], R->X[1][i], R->X[0][i+VGWORDS], R->X[0][i], \
                 R->Z[1][i+VGWORDS], R->Z[1][i], R->Z[0][i+VGWORDS], R->Z[0][i]);
  vR[i] = _SET(0, R->X[1][i], 0, R->X[0][i], \
               0, R->Z[1][i], 0, R->Z[0][i]);

  get_4_isog_1x2x2x2w(vR, C24_A24plus, coeff__0, coeff2_1);
  eval_4_isog_parallel_kg(pts, phiP, phiQ, phiR, coeff__0, coeff2_1, 3);             

  // projective -> affine for phiP, phiQ, phiR
  // note that phiP, phiQ, phiR are now in radix-2^51 Montgomery domain (montR = 2^765)
  
  // 1. convert to radix-2^64
  mpi_conv_51to64(phiP64->X[0], phiP->X[0], NWORDS_FIELD, VNWORDS);
  mpi_conv_51to64(phiP64->X[1], phiP->X[1], NWORDS_FIELD, VNWORDS);
  mpi_conv_51to64(phiP64->Z[0], phiP->Z[0], NWORDS_FIELD, VNWORDS);
  mpi_conv_51to64(phiP64->Z[1], phiP->Z[1], NWORDS_FIELD, VNWORDS);
  mpi_conv_51to64(phiQ64->X[0], phiQ->X[0], NWORDS_FIELD, VNWORDS);
  mpi_conv_51to64(phiQ64->X[1], phiQ->X[1], NWORDS_FIELD, VNWORDS);
  mpi_conv_51to64(phiQ64->Z[0], phiQ->Z[0], NWORDS_FIELD, VNWORDS);
  mpi_conv_51to64(phiQ64->Z[1], phiQ->Z[1], NWORDS_FIELD, VNWORDS);
  mpi_conv_51to64(phiR64->X[0], phiR->X[0], NWORDS_FIELD, VNWORDS);
  mpi_conv_51to64(phiR64->X[1], phiR->X[1], NWORDS_FIELD, VNWORDS);
  mpi_conv_51to64(phiR64->Z[0], phiR->Z[0], NWORDS_FIELD, VNWORDS);
  mpi_conv_51to64(phiR64->Z[1], phiR->Z[1], NWORDS_FIELD, VNWORDS);

  // 2. convert phiP, phiQ, phiR to radix-2^64 Montgomery domain (montR = 2^768)
  // MontMul64(a', 2^771) = [a * 2^765] * 2^771 * 2^(-768) = a * 2^768 = a''
  fpmul_mont(phiP64->X[0], montRx8, phiP64->X[0]);
  fpmul_mont(phiP64->X[1], montRx8, phiP64->X[1]);
  fpmul_mont(phiP64->Z[0], montRx8, phiP64->Z[0]);
  fpmul_mont(phiP64->Z[1], montRx8, phiP64->Z[1]);
  fpmul_mont(phiQ64->X[0], montRx8, phiQ64->X[0]);
  fpmul_mont(phiQ64->X[1], montRx8, phiQ64->X[1]);
  fpmul_mont(phiQ64->Z[0], montRx8, phiQ64->Z[0]);
  fpmul_mont(phiQ64->Z[1], montRx8, phiQ64->Z[1]);
  fpmul_mont(phiR64->X[0], montRx8, phiR64->X[0]);
  fpmul_mont(phiR64->X[1], montRx8, phiR64->X[1]);
  fpmul_mont(phiR64->Z[0], montRx8, phiR64->Z[0]);
  fpmul_mont(phiR64->Z[1], montRx8, phiR64->Z[1]);

  // 3. the rest is the same as PQCrypto-SIDH code
  inv_3_way(phiP64->Z, phiQ64->Z, phiR64->Z);
  fp2mul_mont(phiP64->X, phiP64->Z, phiP64->X);
  fp2mul_mont(phiQ64->X, phiQ64->Z, phiQ64->X);
  fp2mul_mont(phiR64->X, phiR64->Z, phiR64->X);
                
  // Format public key                   
  fp2_encode(phiP64->X, PublicKeyA);
  fp2_encode(phiQ64->X, PublicKeyA + FP2_ENCODED_BYTES);
  fp2_encode(phiR64->X, PublicKeyA + 2*FP2_ENCODED_BYTES);
}

void EphemeralSecretAgreement_A(const unsigned char* PrivateKeyA, const unsigned char* PublicKeyB, unsigned char* SharedSecretA)
{
  point_proj_r51_t R, pts[MAX_INT_POINTS_ALICE];
  f2elm_r51_t PKB[3];
  f2elm_r51_t A24plus = { 0 }, C24 = { 0 }, A = { 0 };
  digit_t SecretKeyA[NWORDS_ORDER] = { 0 };
  int i, row, m, index = 0, pts_index[MAX_INT_POINTS_ALICE], npts = 0, ii = 0;
  vgelm_t vR, C24_A24plus, coeff__0, coeff2_1;
  f2elm_t PKB64[3], jinv64, A64 = { 0 }, C2464 = { 0 }, A24plus64 = { 0 };

  // Initialize images of Bob's basis
  fp2_decode(PublicKeyB, PKB64[0]);
  fp2_decode(PublicKeyB + FP2_ENCODED_BYTES, PKB64[1]);
  fp2_decode(PublicKeyB + 2*FP2_ENCODED_BYTES, PKB64[2]);

  // Initialize constants: A24plus = A+2C, C24 = 4C, where C=1
  get_A(PKB64[0], PKB64[1], PKB64[2], A64);

  // convert PKB[3] to radix-2^51
  // MontMul(PKB, 2^765) = [PKB * 2^768] * 2^765 * 2^(-768) = PKB * 2^765 = PKB''
  fpmul_mont(PKB64[0][0], montRdiv8, PKB64[0][0]);
  fpmul_mont(PKB64[0][1], montRdiv8, PKB64[0][1]);
  fpmul_mont(PKB64[1][0], montRdiv8, PKB64[1][0]);
  fpmul_mont(PKB64[1][1], montRdiv8, PKB64[1][1]);
  fpmul_mont(PKB64[2][0], montRdiv8, PKB64[2][0]);
  fpmul_mont(PKB64[2][1], montRdiv8, PKB64[2][1]);
  mpi_conv_64to51(PKB[0][0], PKB64[0][0], VNWORDS, NWORDS_FIELD);
  mpi_conv_64to51(PKB[0][1], PKB64[0][1], VNWORDS, NWORDS_FIELD);
  mpi_conv_64to51(PKB[1][0], PKB64[1][0], VNWORDS, NWORDS_FIELD);
  mpi_conv_64to51(PKB[1][1], PKB64[1][1], VNWORDS, NWORDS_FIELD);
  mpi_conv_64to51(PKB[2][0], PKB64[2][0], VNWORDS, NWORDS_FIELD);
  mpi_conv_64to51(PKB[2][1], PKB64[2][1], VNWORDS, NWORDS_FIELD);

  // convert A to radix-2^51
  // MontMul(A', 2^765) = [A * 2^768] * 2^765 * 2^(-768) = A * 2^768 = A''
  fpmul_mont(A64[0], montRdiv8, A64[0]);
  fpmul_mont(A64[1], montRdiv8, A64[1]);
  mpi_conv_64to51(A[0], A64[0], VNWORDS, NWORDS_FIELD);
  mpi_conv_64to51(A[1], A64[1], VNWORDS, NWORDS_FIELD);
  mp_add_1w(C24[0], vmont_R, vmont_R);       // C24 = 2
  mp2_add_1w(A24plus, A, C24);               // A24plus = A+2C
  mp_add_1w(C24[0], C24[0], C24[0]);         // C24 = 4

  // Retrieve kernel point
  decode_to_digits(PrivateKeyA, SecretKeyA, SECRETKEY_A_BYTES, NWORDS_ORDER);
  LADDER3PT_1x4x2x1w(PKB[0], PKB[1], PKB[2], SecretKeyA, ALICE, R, A);   

  // intialize C24_A24plus : < C1' | C1 | C0' | C0 | A1' | A1 | A0' | A0 >
  for (i = 0; i < VGWORDS-1; i++) 
    C24_A24plus[i] = _SET(C24[1][i+VGWORDS], C24[1][i], \
                          C24[0][i+VGWORDS], C24[0][i], \
                          A24plus[1][i+VGWORDS], A24plus[1][i], \
                          A24plus[0][i+VGWORDS], A24plus[0][i]);
  C24_A24plus[i] = _SET(0,     C24[1][i], 0, C24[0][i], \
                        0, A24plus[1][i], 0, A24plus[0][i]);

  // Traverse tree
  index = 0;
  for (row = 1; row < MAX_Alice; row++) {
    // update vR : < XR1' | XR1 | XR0' | XR0 | ZR1' | ZR1 | ZR0' | ZR0 >
    for (i = 0; i < VGWORDS-1; i++) 
      vR[i] = _SET(R->X[1][i+VGWORDS], R->X[1][i], R->X[0][i+VGWORDS], R->X[0][i], \
                  R->Z[1][i+VGWORDS], R->Z[1][i], R->Z[0][i+VGWORDS], R->Z[0][i]);
    vR[i] = _SET(0, R->X[1][i], 0, R->X[0][i], \
                0, R->Z[1][i], 0, R->Z[0][i]);
    
    while (index < MAX_Alice-row) {
      pointcopy_1w(pts[npts], R);
      pts_index[npts++] = index;
      m = strat_Alice[ii++];
      xDBLe_1x2x2x2w(vR, R, C24_A24plus, 2*m);  // both vR and R get updated here
      index += m;
    }
    
    get_4_isog_1x2x2x2w(vR, C24_A24plus, coeff__0, coeff2_1);

    // eval_4_isog could be parallelized in 1/2/4/8-way, depending on npts. 
    eval_4_isog_parallel_ss(pts, coeff__0, coeff2_1, npts);

    pointcopy_1w(R, pts[npts-1]);       // update R 
    index = pts_index[npts-1];
    npts -= 1;
  }

  // update vR : < XR1' | XR1 | XR0' | XR0 | ZR1' | ZR1 | ZR0' | ZR0 >
  for (i = 0; i < VGWORDS-1; i++) 
    vR[i] = _SET(R->X[1][i+VGWORDS], R->X[1][i], R->X[0][i+VGWORDS], R->X[0][i], \
                 R->Z[1][i+VGWORDS], R->Z[1][i], R->Z[0][i+VGWORDS], R->Z[0][i]);
  vR[i] = _SET(0, R->X[1][i], 0, R->X[0][i], \
               0, R->Z[1][i], 0, R->Z[0][i]);

  get_4_isog_1x2x2x2w(vR, C24_A24plus, coeff__0, coeff2_1);

  // extract constants
  get_channel_4x2w(A24plus[0], C24_A24plus, 0);
  get_channel_4x2w(A24plus[1], C24_A24plus, 2);
  get_channel_4x2w(C24[0], C24_A24plus, 4);
  get_channel_4x2w(C24[1], C24_A24plus, 6);

  // carry propagation
  carryp_1w(A24plus[0]); 
  carryp_1w(A24plus[1]);
  carryp_1w(C24[0]); 
  carryp_1w(C24[1]);

  // convert A24plus, C24 to radix-2^64
  mpi_conv_51to64(A24plus64[0], A24plus[0], NWORDS_FIELD, VNWORDS);
  mpi_conv_51to64(A24plus64[1], A24plus[1], NWORDS_FIELD, VNWORDS);
  mpi_conv_51to64(C2464[0], C24[0], NWORDS_FIELD, VNWORDS);
  mpi_conv_51to64(C2464[1], C24[1], NWORDS_FIELD, VNWORDS);
  // MontMul64(a', 2^771) = [a * 2^765] * 2^771 * 2^(-768) = a * 2^768 = a''
  fpmul_mont(A24plus64[0], montRx8, A24plus64[0]);
  fpmul_mont(A24plus64[1], montRx8, A24plus64[1]);
  fpmul_mont(C2464[0], montRx8, C2464[0]);
  fpmul_mont(C2464[1], montRx8, C2464[1]);

  fp2add(A24plus64, A24plus64, A24plus64);    // A24plus = 2A+4C
  fp2sub(A24plus64, C2464, A24plus64);        // A24plus = 2A
  fp2add(A24plus64, A24plus64, A24plus64);    // A24plus = 4A
  j_inv(A24plus64, C2464, jinv64);  
  fp2_encode(jinv64, SharedSecretA);          // Format shared secret 
}

void EphemeralKeyGeneration_B(const unsigned char* PrivateKeyB, unsigned char* PublicKeyB)
{
  point_proj_r51_t R, phiP = { 0 }, phiQ = { 0 }, phiR = { 0 }, pts[MAX_INT_POINTS_BOB];
  f2elm_r51_t XPB, XQB, XRB, A24plus = { 0 }, C24 = { 0 }, A = { 0 };
  digit_t SecretKeyB[NWORDS_ORDER] = { 0 };
  int i, row, m, index = 0, pts_index[MAX_INT_POINTS_BOB], npts = 0, ii = 0;
  vgelm_t vR, C24_A24plus, coeff1_0;
  point_proj_t phiP64 = { 0 }, phiQ64 = { 0 }, phiR64 = { 0 };

  // Initialize basis points
  init_basis_1w((digit_t*)vB_gen, XPB, XQB, XRB); 
  init_basis_1w((digit_t*)vA_gen, phiP->X, phiQ->X, phiR->X);
  fpcopy_1w(phiP->Z[0], vmont_R);
  fpcopy_1w(phiQ->Z[0], vmont_R);
  fpcopy_1w(phiR->Z[0], vmont_R);

  // Initialize constants: A24plus = A+2C, C24 = 4C, where A=6, C=1
  fpcopy_1w(A24plus[0], vmont_R);       // A24plus = 1
  mp2_add_1w(A24plus, A24plus, A24plus);// A24plus = 2
  mp2_add_1w(C24, A24plus, A24plus);    // C24 = 4
  mp2_add_1w(A, A24plus, C24);          // A = 6 
  mp2_add_1w(A24plus, C24, C24);        // A24plus = 8

  // Retrieve kernel point
  decode_to_digits(PrivateKeyB, SecretKeyB, SECRETKEY_B_BYTES, NWORDS_ORDER);
  // vectorized Montgomery ladder (ladder step is 1x4x2x1w)
  // the ouput R is in radix-2^51
  LADDER3PT_1x4x2x1w(XPB, XQB, XRB, SecretKeyB, BOB, R, A); 

  // intialize C24_A24plus : < C1' | C1 | C0' | C0 | A1' | A1 | A0' | A0 >
  for (i = 0; i < VGWORDS-1; i++) 
    C24_A24plus[i] = _SET(C24[1][i+VGWORDS], C24[1][i], \
                          C24[0][i+VGWORDS], C24[0][i], \
                          A24plus[1][i+VGWORDS], A24plus[1][i], \
                          A24plus[0][i+VGWORDS], A24plus[0][i]);
  C24_A24plus[i] = _SET(0,     C24[1][i], 0, C24[0][i], \
                        0, A24plus[1][i], 0, A24plus[0][i]);

  // Traverse tree 
  index = 0;
  for (row = 1; row < MAX_Bob; row++) {
    // update vR : < XR1' | XR1 | XR0' | XR0 | ZR1' | ZR1 | ZR0' | ZR0 >
    for (i = 0; i < VGWORDS-1; i++) 
      vR[i] = _SET(R->X[1][i+VGWORDS], R->X[1][i], R->X[0][i+VGWORDS], R->X[0][i], \
                  R->Z[1][i+VGWORDS], R->Z[1][i], R->Z[0][i+VGWORDS], R->Z[0][i]);
    vR[i] = _SET(0, R->X[1][i], 0, R->X[0][i], \
                 0, R->Z[1][i], 0, R->Z[0][i]);
    
    while (index < MAX_Bob-row) {
      pointcopy_1w(pts[npts], R);
      pts_index[npts++] = index;
      m = strat_Bob[ii++];
      xTPLe_1x2x2x2w(vR, R, C24_A24plus, m);    // both vR and R get updated here
      index += m;
    } 

    get_3_isog_1x2x2x2w(vR, C24_A24plus, coeff1_0);

    // eval_3_isog could be parallelized in 1/2/4/8-way, depending on npts. 
    eval_3_isog_parallel_kg(pts, phiP, phiQ, phiR, coeff1_0, npts+3);

    pointcopy_1w(R, pts[npts-1]);       // update R 
    index = pts_index[npts-1];
    npts -= 1;
  }

  // update vR : < XR1' | XR1 | XR0' | XR0 | ZR1' | ZR1 | ZR0' | ZR0 >
  for (i = 0; i < VGWORDS-1; i++) 
    vR[i] = _SET(R->X[1][i+VGWORDS], R->X[1][i], R->X[0][i+VGWORDS], R->X[0][i], \
                 R->Z[1][i+VGWORDS], R->Z[1][i], R->Z[0][i+VGWORDS], R->Z[0][i]);
  vR[i] = _SET(0, R->X[1][i], 0, R->X[0][i], \
               0, R->Z[1][i], 0, R->Z[0][i]);

  get_3_isog_1x2x2x2w(vR, C24_A24plus, coeff1_0);
  eval_3_isog_parallel_kg(pts, phiP, phiQ, phiR, coeff1_0, 3);           

  // projective -> affine for phiP, phiQ, phiR
  // note that phiP, phiQ, phiR are now in radix-2^51 Montgomery domain (montR = 2^765)
  
  // 1. convert to radix-2^64
  mpi_conv_51to64(phiP64->X[0], phiP->X[0], NWORDS_FIELD, VNWORDS);
  mpi_conv_51to64(phiP64->X[1], phiP->X[1], NWORDS_FIELD, VNWORDS);
  mpi_conv_51to64(phiP64->Z[0], phiP->Z[0], NWORDS_FIELD, VNWORDS);
  mpi_conv_51to64(phiP64->Z[1], phiP->Z[1], NWORDS_FIELD, VNWORDS);
  mpi_conv_51to64(phiQ64->X[0], phiQ->X[0], NWORDS_FIELD, VNWORDS);
  mpi_conv_51to64(phiQ64->X[1], phiQ->X[1], NWORDS_FIELD, VNWORDS);
  mpi_conv_51to64(phiQ64->Z[0], phiQ->Z[0], NWORDS_FIELD, VNWORDS);
  mpi_conv_51to64(phiQ64->Z[1], phiQ->Z[1], NWORDS_FIELD, VNWORDS);
  mpi_conv_51to64(phiR64->X[0], phiR->X[0], NWORDS_FIELD, VNWORDS);
  mpi_conv_51to64(phiR64->X[1], phiR->X[1], NWORDS_FIELD, VNWORDS);
  mpi_conv_51to64(phiR64->Z[0], phiR->Z[0], NWORDS_FIELD, VNWORDS);
  mpi_conv_51to64(phiR64->Z[1], phiR->Z[1], NWORDS_FIELD, VNWORDS);

  // 2. convert phiP, phiQ, phiR to radix-2^64 Montgomery domain (montR = 2^768)
  // MontMul64(a', 2^771) = [a * 2^765] * 2^771 * 2^(-768) = a * 2^768 = a''
  fpmul_mont(phiP64->X[0], montRx8, phiP64->X[0]);
  fpmul_mont(phiP64->X[1], montRx8, phiP64->X[1]);
  fpmul_mont(phiP64->Z[0], montRx8, phiP64->Z[0]);
  fpmul_mont(phiP64->Z[1], montRx8, phiP64->Z[1]);
  fpmul_mont(phiQ64->X[0], montRx8, phiQ64->X[0]);
  fpmul_mont(phiQ64->X[1], montRx8, phiQ64->X[1]);
  fpmul_mont(phiQ64->Z[0], montRx8, phiQ64->Z[0]);
  fpmul_mont(phiQ64->Z[1], montRx8, phiQ64->Z[1]);
  fpmul_mont(phiR64->X[0], montRx8, phiR64->X[0]);
  fpmul_mont(phiR64->X[1], montRx8, phiR64->X[1]);
  fpmul_mont(phiR64->Z[0], montRx8, phiR64->Z[0]);
  fpmul_mont(phiR64->Z[1], montRx8, phiR64->Z[1]);

  // 3. the rest is the same as PQCrypto-SIDH code
  inv_3_way(phiP64->Z, phiQ64->Z, phiR64->Z);
  fp2mul_mont(phiP64->X, phiP64->Z, phiP64->X);
  fp2mul_mont(phiQ64->X, phiQ64->Z, phiQ64->X);
  fp2mul_mont(phiR64->X, phiR64->Z, phiR64->X);
                
  // Format public key                   
  fp2_encode(phiP64->X, PublicKeyB);
  fp2_encode(phiQ64->X, PublicKeyB + FP2_ENCODED_BYTES);
  fp2_encode(phiR64->X, PublicKeyB + 2*FP2_ENCODED_BYTES);
}

void EphemeralSecretAgreement_B(const unsigned char* PrivateKeyB, const unsigned char* PublicKeyA, unsigned char* SharedSecretB)
{
  point_proj_r51_t R, pts[MAX_INT_POINTS_BOB];
  f2elm_r51_t PKB[3];
  f2elm_r51_t A24plus = { 0 }, C24 = { 0 }, A = { 0 };
  digit_t SecretKeyB[NWORDS_ORDER] = { 0 };
  int i, row, m, index = 0, pts_index[MAX_INT_POINTS_ALICE], npts = 0, ii = 0;
  vgelm_t vR, C24_A24plus, coeff1_0;
  f2elm_t PKB64[3], jinv64, A64 = { 0 }, C2464 = { 0 }, A24plus64 = { 0 };

  // Initialize images of Alice's basis
  fp2_decode(PublicKeyA, PKB64[0]);
  fp2_decode(PublicKeyA + FP2_ENCODED_BYTES, PKB64[1]);
  fp2_decode(PublicKeyA + 2*FP2_ENCODED_BYTES, PKB64[2]);

  // Initialize constants: A24plus = A+2C, C24 = 4C, where C=1
  get_A(PKB64[0], PKB64[1], PKB64[2], A64);

  // convert PKB[3] to radix-2^51
  // MontMul(PKB, 2^765) = [PKB * 2^768] * 2^765 * 2^(-768) = PKB * 2^765 = PKB''
  fpmul_mont(PKB64[0][0], montRdiv8, PKB64[0][0]);
  fpmul_mont(PKB64[0][1], montRdiv8, PKB64[0][1]);
  fpmul_mont(PKB64[1][0], montRdiv8, PKB64[1][0]);
  fpmul_mont(PKB64[1][1], montRdiv8, PKB64[1][1]);
  fpmul_mont(PKB64[2][0], montRdiv8, PKB64[2][0]);
  fpmul_mont(PKB64[2][1], montRdiv8, PKB64[2][1]);
  mpi_conv_64to51(PKB[0][0], PKB64[0][0], VNWORDS, NWORDS_FIELD);
  mpi_conv_64to51(PKB[0][1], PKB64[0][1], VNWORDS, NWORDS_FIELD);
  mpi_conv_64to51(PKB[1][0], PKB64[1][0], VNWORDS, NWORDS_FIELD);
  mpi_conv_64to51(PKB[1][1], PKB64[1][1], VNWORDS, NWORDS_FIELD);
  mpi_conv_64to51(PKB[2][0], PKB64[2][0], VNWORDS, NWORDS_FIELD);
  mpi_conv_64to51(PKB[2][1], PKB64[2][1], VNWORDS, NWORDS_FIELD);

  // convert A to radix-2^51
  // MontMul(A', 2^765) = [A * 2^768] * 2^765 * 2^(-768) = A * 2^768 = A''
  fpmul_mont(A64[0], montRdiv8, A64[0]);
  fpmul_mont(A64[1], montRdiv8, A64[1]);
  mpi_conv_64to51(A[0], A64[0], VNWORDS, NWORDS_FIELD);
  mpi_conv_64to51(A[1], A64[1], VNWORDS, NWORDS_FIELD);
  mp_add_1w(C24[0], vmont_R, vmont_R);       // C24 = 2
  mp2_add_1w(A24plus, A, C24);               // A24plus = A+2C
  mp_add_1w(C24[0], C24[0], C24[0]);         // C24 = 4

  // Retrieve kernel point
  decode_to_digits(PrivateKeyB, SecretKeyB, SECRETKEY_B_BYTES, NWORDS_ORDER);
  LADDER3PT_1x4x2x1w(PKB[0], PKB[1], PKB[2], SecretKeyB, BOB, R, A);   

  // intialize C24_A24plus : < C1' | C1 | C0' | C0 | A1' | A1 | A0' | A0 >
  for (i = 0; i < VGWORDS-1; i++) 
    C24_A24plus[i] = _SET(C24[1][i+VGWORDS], C24[1][i], \
                          C24[0][i+VGWORDS], C24[0][i], \
                          A24plus[1][i+VGWORDS], A24plus[1][i], \
                          A24plus[0][i+VGWORDS], A24plus[0][i]);
  C24_A24plus[i] = _SET(0,     C24[1][i], 0, C24[0][i], \
                        0, A24plus[1][i], 0, A24plus[0][i]);

  // Traverse tree 
  index = 0;
  for (row = 1; row < MAX_Bob; row++) {
    // update vR : < XR1' | XR1 | XR0' | XR0 | ZR1' | ZR1 | ZR0' | ZR0 >
    for (i = 0; i < VGWORDS-1; i++) 
      vR[i] = _SET(R->X[1][i+VGWORDS], R->X[1][i], R->X[0][i+VGWORDS], R->X[0][i], \
                  R->Z[1][i+VGWORDS], R->Z[1][i], R->Z[0][i+VGWORDS], R->Z[0][i]);
    vR[i] = _SET(0, R->X[1][i], 0, R->X[0][i], \
                 0, R->Z[1][i], 0, R->Z[0][i]);
    
    while (index < MAX_Bob-row) {
      pointcopy_1w(pts[npts], R);
      pts_index[npts++] = index;
      m = strat_Bob[ii++];
      xTPLe_1x2x2x2w(vR, R, C24_A24plus, m);    // both vR and R get updated here
      index += m;
    } 

    get_3_isog_1x2x2x2w(vR, C24_A24plus, coeff1_0);

    // eval_3_isog could be parallelized in 1/2/4/8-way, depending on npts. 
    eval_3_isog_parallel_ss(pts, coeff1_0, npts);

    pointcopy_1w(R, pts[npts-1]);       // update R 
    index = pts_index[npts-1];
    npts -= 1;
  }

  // update vR : < XR1' | XR1 | XR0' | XR0 | ZR1' | ZR1 | ZR0' | ZR0 >
  for (i = 0; i < VGWORDS-1; i++) 
    vR[i] = _SET(R->X[1][i+VGWORDS], R->X[1][i], R->X[0][i+VGWORDS], R->X[0][i], \
                 R->Z[1][i+VGWORDS], R->Z[1][i], R->Z[0][i+VGWORDS], R->Z[0][i]);
  vR[i] = _SET(0, R->X[1][i], 0, R->X[0][i], \
               0, R->Z[1][i], 0, R->Z[0][i]);

  get_3_isog_1x2x2x2w(vR, C24_A24plus, coeff1_0);

  // extract constants
  get_channel_4x2w(A24plus[0], C24_A24plus, 0);
  get_channel_4x2w(A24plus[1], C24_A24plus, 2);
  get_channel_4x2w(C24[0], C24_A24plus, 4);
  get_channel_4x2w(C24[1], C24_A24plus, 6);

  // carry propagation
  carryp_1w(A24plus[0]); 
  carryp_1w(A24plus[1]);
  carryp_1w(C24[0]); 
  carryp_1w(C24[1]);

  // convert A24plus, C24 to radix-2^64
  mpi_conv_51to64(A24plus64[0], A24plus[0], NWORDS_FIELD, VNWORDS);
  mpi_conv_51to64(A24plus64[1], A24plus[1], NWORDS_FIELD, VNWORDS);
  mpi_conv_51to64(C2464[0], C24[0], NWORDS_FIELD, VNWORDS);
  mpi_conv_51to64(C2464[1], C24[1], NWORDS_FIELD, VNWORDS);
  // MontMul64(a', 2^771) = [a * 2^765] * 2^771 * 2^(-768) = a * 2^768 = a''
  fpmul_mont(A24plus64[0], montRx8, A24plus64[0]);
  fpmul_mont(A24plus64[1], montRx8, A24plus64[1]);
  fpmul_mont(C2464[0], montRx8, C2464[0]);
  fpmul_mont(C2464[1], montRx8, C2464[1]);

  fp2add(A24plus64, A24plus64, A24plus64);    // A24plus = 2A+4C
  fp2sub(A24plus64, C2464, A24plus64);        // A24plus = 2A
  fp2add(A24plus64, A24plus64, A24plus64);    // A24plus = 4A
  j_inv(A24plus64, C2464, jinv64);  
  fp2_encode(jinv64, SharedSecretB);          // Format shared secret 
}

void EphemeralKeyGenSecAgr_A_parallel(const unsigned char* PrivateKeyA, unsigned char* PublicKeyA, const unsigned char* PublicKeyB, unsigned char* SharedSecretA)
{
  point_proj_r51_t R, phiP = { 0 }, phiQ = { 0 }, phiR = { 0 }, pts[MAX_INT_POINTS_ALICE];
  f2elm_r51_t XPA, XQA, XRA, A24plus = { 0 }, C24 = { 0 }, A = { 0 };
  digit_t SecretKeyA[NWORDS_ORDER] = { 0 };
  int i, row, m, index = 0, pts_index[MAX_INT_POINTS_ALICE], npts = 0, ii = 0;
  point_proj_t phiP64 = { 0 }, phiQ64 = { 0 }, phiR64 = { 0 };

  point_proj_r51_t _R, _pts[MAX_INT_POINTS_ALICE];
  f2elm_r51_t _PKB[3];
  f2elm_r51_t _A24plus = { 0 }, _C24 = { 0 }, _A = { 0 };
  f2elm_t _PKB64[3], _jinv64, _A64 = { 0 }, _C2464 = { 0 }, _A24plus64 = { 0 };

  vfelm_t vR, C24_A24plus, coeff__0, coeff2_1;

  // ---------------------------------------------------------------------------
  // KeyGen part

  // Initialize basis points
  init_basis_1w((digit_t*)vA_gen, XPA, XQA, XRA); 
  init_basis_1w((digit_t*)vB_gen, phiP->X, phiQ->X, phiR->X);
  fpcopy_1w(phiP->Z[0], vmont_R);
  fpcopy_1w(phiQ->Z[0], vmont_R);
  fpcopy_1w(phiR->Z[0], vmont_R);

  // Initialize constants: A24plus = A+2C, C24 = 4C, where A=6, C=1
  fpcopy_1w(A24plus[0], vmont_R);       // A24plus = 1
  mp2_add_1w(A24plus, A24plus, A24plus);// A24plus = 2
  mp2_add_1w(C24, A24plus, A24plus);    // C24 = 4
  mp2_add_1w(A, A24plus, C24);          // A = 6 
  mp2_add_1w(A24plus, C24, C24);        // A24plus = 8

  // ---------------------------------------------------------------------------
  // SecAgr part

  // Initialize images of Bob's basis
  fp2_decode(PublicKeyB, _PKB64[0]);
  fp2_decode(PublicKeyB + FP2_ENCODED_BYTES, _PKB64[1]);
  fp2_decode(PublicKeyB + 2*FP2_ENCODED_BYTES, _PKB64[2]);

  // Initialize constants: A24plus = A+2C, C24 = 4C, where C=1
  get_A(_PKB64[0], _PKB64[1], _PKB64[2], _A64);

  // convert PKB[3] to radix-2^51
  // MontMul(PKB, 2^765) = [PKB * 2^768] * 2^765 * 2^(-768) = PKB * 2^765 = PKB''
  fpmul_mont(_PKB64[0][0], montRdiv8, _PKB64[0][0]);
  fpmul_mont(_PKB64[0][1], montRdiv8, _PKB64[0][1]);
  fpmul_mont(_PKB64[1][0], montRdiv8, _PKB64[1][0]);
  fpmul_mont(_PKB64[1][1], montRdiv8, _PKB64[1][1]);
  fpmul_mont(_PKB64[2][0], montRdiv8, _PKB64[2][0]);
  fpmul_mont(_PKB64[2][1], montRdiv8, _PKB64[2][1]);
  mpi_conv_64to51(_PKB[0][0], _PKB64[0][0], VNWORDS, NWORDS_FIELD);
  mpi_conv_64to51(_PKB[0][1], _PKB64[0][1], VNWORDS, NWORDS_FIELD);
  mpi_conv_64to51(_PKB[1][0], _PKB64[1][0], VNWORDS, NWORDS_FIELD);
  mpi_conv_64to51(_PKB[1][1], _PKB64[1][1], VNWORDS, NWORDS_FIELD);
  mpi_conv_64to51(_PKB[2][0], _PKB64[2][0], VNWORDS, NWORDS_FIELD);
  mpi_conv_64to51(_PKB[2][1], _PKB64[2][1], VNWORDS, NWORDS_FIELD);

  // convert A to radix-2^51
  // MontMul(A', 2^765) = [A * 2^768] * 2^765 * 2^(-768) = A * 2^768 = A''
  fpmul_mont(_A64[0], montRdiv8, _A64[0]);
  fpmul_mont(_A64[1], montRdiv8, _A64[1]);
  mpi_conv_64to51(_A[0], _A64[0], VNWORDS, NWORDS_FIELD);
  mpi_conv_64to51(_A[1], _A64[1], VNWORDS, NWORDS_FIELD);
  mp_add_1w(_C24[0], vmont_R, vmont_R);         // C24 = 2
  mp2_add_1w(_A24plus, _A, _C24);               // A24plus = A+2C
  mp_add_1w(_C24[0], _C24[0], _C24[0]);         // C24 = 4

  // ---------------------------------------------------------------------------
  // Parallel part

  // Retrieve kernel point
  decode_to_digits(PrivateKeyA, SecretKeyA, SECRETKEY_A_BYTES, NWORDS_ORDER);

  // LADDER3PT_1x4x2x1w(XPA, XQA, XRA, SecretKeyA, ALICE, R, A); 
  // LADDER3PT_1x4x2x1w(_PKB[0], _PKB[1], _PKB[2], SecretKeyA, ALICE, _R, _A); 
  LADDER3PT_2x4x1x1w(XPA, XQA, XRA, _PKB[0], _PKB[1], _PKB[2], SecretKeyA, ALICE, R, A, _R, _A);  

  // intialize C24_A24plus :  <  C1 | C0 | A1 | A0 | _C1 | _C0 | _A1 | _A0 >
  for (i = 0; i < VNWORDS; i++) 
    C24_A24plus[i] = _SET( C24[1][i],  C24[0][i],  A24plus[1][i],  A24plus[0][i], \
                          _C24[1][i], _C24[0][i], _A24plus[1][i], _A24plus[0][i]);

  // Traverse tree 
  // Use 2x2x2x1w isog instead of 1x2x2x2w isog -> 8x1w fp instead of 4x2w fp (big difference)
  index = 0;
  for (row = 1; row < MAX_Alice; row++) {
    // update vR : < XR1 |  XR0 | ZR1 | ZR0 | _XR1 | _XR0 | _ZR1 |  _ZR0 >
    for (i = 0; i < VNWORDS; i++) 
      vR[i] = _SET( R->X[1][i],  R->X[0][i],  R->Z[1][i],  R->Z[0][i], \
                   _R->X[1][i], _R->X[0][i], _R->Z[1][i], _R->Z[0][i]);

    while (index < MAX_Alice-row) {
      pointcopy_1w(pts[npts], R);
      pointcopy_1w(_pts[npts], _R);
      pts_index[npts++] = index;
      m = strat_Alice[ii++];
      xDBLe_2x2x2x1w(vR, R, _R, C24_A24plus, 2*m);  // all vR, R, _R get updated here
      index += m;
    }

    get_4_isog_2x2x2x1w(vR, C24_A24plus, coeff__0, coeff2_1);

    // eval_4_isog could be parallelized in 1/2/4/8-way, depending on npts. 
    eval_4_isog_parallel_kgss(pts, _pts, phiP, phiQ, phiR, coeff__0, coeff2_1, 2*npts+3);

    pointcopy_1w(R, pts[npts-1]);       // update R 
    pointcopy_1w(_R, _pts[npts-1]);     // update _R 
    index = pts_index[npts-1];
    npts -= 1;
  }

  // update vR : < XR1 |  XR0 | ZR1 | ZR0 | _XR1 | _XR0 | _ZR1 |  _ZR0 >
  for (i = 0; i < VNWORDS; i++) 
    vR[i] = _SET( R->X[1][i],  R->X[0][i],  R->Z[1][i],  R->Z[0][i], \
                 _R->X[1][i], _R->X[0][i], _R->Z[1][i], _R->Z[0][i]);

  get_4_isog_2x2x2x1w(vR, C24_A24plus, coeff__0, coeff2_1);

  eval_4_isog_parallel_kgss(pts, _pts, phiP, phiQ, phiR, coeff__0, coeff2_1, 3); 

  // ---------------------------------------------------------------------------
  // KeyGen part
  // projective -> affine for phiP, phiQ, phiR
  // note that phiP, phiQ, phiR are now in radix-2^51 Montgomery domain (montR = 2^765)
  
  // 1. convert to radix-2^64
  mpi_conv_51to64(phiP64->X[0], phiP->X[0], NWORDS_FIELD, VNWORDS);
  mpi_conv_51to64(phiP64->X[1], phiP->X[1], NWORDS_FIELD, VNWORDS);
  mpi_conv_51to64(phiP64->Z[0], phiP->Z[0], NWORDS_FIELD, VNWORDS);
  mpi_conv_51to64(phiP64->Z[1], phiP->Z[1], NWORDS_FIELD, VNWORDS);
  mpi_conv_51to64(phiQ64->X[0], phiQ->X[0], NWORDS_FIELD, VNWORDS);
  mpi_conv_51to64(phiQ64->X[1], phiQ->X[1], NWORDS_FIELD, VNWORDS);
  mpi_conv_51to64(phiQ64->Z[0], phiQ->Z[0], NWORDS_FIELD, VNWORDS);
  mpi_conv_51to64(phiQ64->Z[1], phiQ->Z[1], NWORDS_FIELD, VNWORDS);
  mpi_conv_51to64(phiR64->X[0], phiR->X[0], NWORDS_FIELD, VNWORDS);
  mpi_conv_51to64(phiR64->X[1], phiR->X[1], NWORDS_FIELD, VNWORDS);
  mpi_conv_51to64(phiR64->Z[0], phiR->Z[0], NWORDS_FIELD, VNWORDS);
  mpi_conv_51to64(phiR64->Z[1], phiR->Z[1], NWORDS_FIELD, VNWORDS);

  // 2. convert phiP, phiQ, phiR to radix-2^64 Montgomery domain (montR = 2^768)
  // MontMul64(a', 2^771) = [a * 2^765] * 2^771 * 2^(-768) = a * 2^768 = a''
  fpmul_mont(phiP64->X[0], montRx8, phiP64->X[0]);
  fpmul_mont(phiP64->X[1], montRx8, phiP64->X[1]);
  fpmul_mont(phiP64->Z[0], montRx8, phiP64->Z[0]);
  fpmul_mont(phiP64->Z[1], montRx8, phiP64->Z[1]);
  fpmul_mont(phiQ64->X[0], montRx8, phiQ64->X[0]);
  fpmul_mont(phiQ64->X[1], montRx8, phiQ64->X[1]);
  fpmul_mont(phiQ64->Z[0], montRx8, phiQ64->Z[0]);
  fpmul_mont(phiQ64->Z[1], montRx8, phiQ64->Z[1]);
  fpmul_mont(phiR64->X[0], montRx8, phiR64->X[0]);
  fpmul_mont(phiR64->X[1], montRx8, phiR64->X[1]);
  fpmul_mont(phiR64->Z[0], montRx8, phiR64->Z[0]);
  fpmul_mont(phiR64->Z[1], montRx8, phiR64->Z[1]);

  // 3. the rest is the same as PQCrypto-SIDH code
  inv_3_way(phiP64->Z, phiQ64->Z, phiR64->Z);
  fp2mul_mont(phiP64->X, phiP64->Z, phiP64->X);
  fp2mul_mont(phiQ64->X, phiQ64->Z, phiQ64->X);
  fp2mul_mont(phiR64->X, phiR64->Z, phiR64->X);
                
  // Format public key                   
  fp2_encode(phiP64->X, PublicKeyA);
  fp2_encode(phiQ64->X, PublicKeyA + FP2_ENCODED_BYTES);
  fp2_encode(phiR64->X, PublicKeyA + 2*FP2_ENCODED_BYTES);

  // ---------------------------------------------------------------------------
  // SecAgr part

  // extract constants
  get_channel_8x1w(_A24plus[0], C24_A24plus, 0);
  get_channel_8x1w(_A24plus[1], C24_A24plus, 1);
  get_channel_8x1w(_C24[0], C24_A24plus, 2);
  get_channel_8x1w(_C24[1], C24_A24plus, 3);

  // convert A24plus, C24 to radix-2^64
  mpi_conv_51to64(_A24plus64[0], _A24plus[0], NWORDS_FIELD, VNWORDS);
  mpi_conv_51to64(_A24plus64[1], _A24plus[1], NWORDS_FIELD, VNWORDS);
  mpi_conv_51to64(_C2464[0], _C24[0], NWORDS_FIELD, VNWORDS);
  mpi_conv_51to64(_C2464[1], _C24[1], NWORDS_FIELD, VNWORDS);
  // MontMul64(a', 2^771) = [a * 2^765] * 2^771 * 2^(-768) = a * 2^768 = a''
  fpmul_mont(_A24plus64[0], montRx8, _A24plus64[0]);
  fpmul_mont(_A24plus64[1], montRx8, _A24plus64[1]);
  fpmul_mont(_C2464[0], montRx8, _C2464[0]);
  fpmul_mont(_C2464[1], montRx8, _C2464[1]);

  fp2add(_A24plus64, _A24plus64, _A24plus64);    // A24plus = 2A+4C
  fp2sub(_A24plus64, _C2464, _A24plus64);        // A24plus = 2A
  fp2add(_A24plus64, _A24plus64, _A24plus64);    // A24plus = 4A
  j_inv(_A24plus64, _C2464, _jinv64);  
  fp2_encode(_jinv64, SharedSecretA);          // Format shared secret 
}

