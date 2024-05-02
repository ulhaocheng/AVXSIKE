/**
 *******************************************************************************
 * @version 0.1
 * @date 2021-12-16
 * @copyright Copyright © 2021 by University of Luxembourg
 * @author Hao Cheng
 *******************************************************************************
 */

#include "fp.h"
#include "fpx.h"
#include "curve.h"
#include "sidh.h"
#include "sike.h"
#include "utils.h"
#include "testvec.h"
#include <string.h>

// the function of measuring CPU cycles 
extern uint64_t read_tsc();

// MACROs for benchmarking 

#define ITER_L 10000
#define ITER_M 1000
#define ITER_S 100

#define LOAD_CACHE(X, ITER) for (i = 0; i < (ITER); i++) (X)

#define MEASURE_TIME(X, ITER)                         \
  start_cycles = read_tsc();                          \
  for (i = 0; i < (ITER); i++) (X);                   \
  end_cycles = read_tsc();                            \
  diff_cycles = (end_cycles-start_cycles)/(ITER)

static void carryp(uint64_t *a)
{
  int i;

  for (i = 0; i < VNWORDS-1; i++) {
    a[i+1] += a[i]>>VBRADIX;
    a[i] &= VBMASK;
  }
}

void test_fp()
{
  uint64_t a64[NWORDS_FIELD] = { TV_A1 }, b64[NWORDS_FIELD] = { TV_B1 }, r64[NWORDS_FIELD] = { 0 };
  uint64_t a51[VNWORDS], b51[VNWORDS], r51[VNWORDS];
  __m512i a_8x1w[VNWORDS], b_8x1w[VNWORDS], r_8x1w[VNWORDS], z_8x1w[2*VNWORDS];
  __m512i a_4x2w[VGWORDS], b_4x2w[VGWORDS], r_4x2w[VGWORDS], z_4x2w[3*VGWORDS];
  int i;

  mpi_conv_64to51(a51, a64, VNWORDS, NWORDS_FIELD);
  mpi_conv_64to51(b51, b64, VNWORDS, NWORDS_FIELD);

  for (i = 0; i < VNWORDS; i++) {
    a_8x1w[i] = VSET1(a51[i]);
    b_8x1w[i] = VSET1(b51[i]); 
  }

  for (i = 0; i < VGWORDS; i++) {
    a_4x2w[i] = VSET(0, 0, 0, 0, 0, 0, a51[i+VGWORDS], a51[i]);
    b_4x2w[i] = VSET(0, 0, 0, 0, 0, 0, b51[i+VGWORDS], b51[i]);
  }

  // 8x1w Mont mul 
  mp_mul_8x1w(z_8x1w, a_8x1w, b_8x1w);
  rdc_mont_8x1w(r_8x1w, z_8x1w);

  get_channel_8x1w(r51, r_8x1w, 0);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 8x1w MONT MUL R: ", r64, NWORDS_FIELD);

  // 4x2w Mont mul
  mp_mul_4x2w(z_4x2w, a_4x2w, b_4x2w);
  rdc_mont_4x2w(r_4x2w, z_4x2w);

  get_channel_4x2w(r51, r_4x2w, 0);
  carryp(r51);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 4x2w MONT MUL R: ", r64, NWORDS_FIELD);


  // ---------------------------------------------------------------------------
  // measure timings 
  uint64_t start_cycles, end_cycles, diff_cycles;

  LOAD_CACHE(mp_mul_8x1w(z_8x1w, a_8x1w, b_8x1w), ITER_M);
  MEASURE_TIME(mp_mul_8x1w(z_8x1w, a_8x1w, b_8x1w), ITER_L);
  printf("* 8x1w INT MUL: %ld\n", diff_cycles);

  LOAD_CACHE(mp_mul_4x2w(z_4x2w, a_4x2w, b_4x2w), ITER_M);
  MEASURE_TIME(mp_mul_4x2w(z_4x2w, a_4x2w, b_4x2w), ITER_L);
  printf("* 4x2w INT MUL: %ld\n", diff_cycles);
}

void test_fpx()
{
  vf2elm_t a_8x1x1w, b_8x1x1w, r_8x1x1w;
  vfelm_t  a_4x2x1w, b_4x2x1w, r_4x2x1w;
  vg2elm_t a_4x1x2w, b_4x1x2w, r_4x1x2w;
  vfelm_t  a_2x4x1w, b_2x4x1w, r_2x4x1w; 
  vgelm_t  a_2x2x2w, b_2x2x2w, r_2x2x2w;
  uint64_t a64[NWORDS_FIELD] = { TV_A0 }, b64[NWORDS_FIELD] = { TV_B0 }, r64[NWORDS_FIELD] = { 0 };
  uint64_t a51[VNWORDS] = { 0 }, b51[VNWORDS] = { 0 }, r51[VNWORDS];
  int i;

  mpi_conv_64to51(a51, a64, VNWORDS, NWORDS_FIELD);
  mpi_conv_64to51(b51, b64, VNWORDS, NWORDS_FIELD);

  for (i = 0; i < VNWORDS; i++) {
    a_8x1x1w[0][i] = VSET1(a51[i]);
    a_8x1x1w[1][i] = VSET1(b51[i]);
    a_4x2x1w[i] = VSET(0, 0, 0, 0, 0, 0, b51[i], a51[i]);
    a_2x4x1w[i] = VSET(0, 0, 0, 0, 0, b51[i], 0, a51[i]);
  }
  for (i = 0; i < VGWORDS; i++) {
    a_4x1x2w[0][i] = VSET(0, 0, 0, 0, 0, 0, a51[i+VGWORDS], a51[i]);
    a_4x1x2w[1][i] = VSET(0, 0, 0, 0, 0, 0, b51[i+VGWORDS], b51[i]);
    a_2x2x2w[i] = VSET(0, 0, 0, 0, b51[i+VGWORDS], b51[i], a51[i+VGWORDS], a51[i]);
  }

  puts("\nCORRECTNESS TEST: ");

  // 8x1x1w Fp2 mul 
  fp2mul_mont_8x1x1w(r_8x1x1w, a_8x1x1w, a_8x1x1w);
  
  get_channel_8x1w(r51, r_8x1x1w[0], 0);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 8x1x1w FPX MUL R[0]: ", r64, NWORDS_FIELD);
  get_channel_8x1w(r51, r_8x1x1w[1], 0);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 8x1x1w FPX MUL R[1]: ", r64, NWORDS_FIELD);

  // 8x1x1w Fp2 sqr
  fp2sqr_mont_8x1x1w(r_8x1x1w, a_8x1x1w);
  
  get_channel_8x1w(r51, r_8x1x1w[0], 0);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 8x1x1w FPX SQR R[0]: ", r64, NWORDS_FIELD);
  get_channel_8x1w(r51, r_8x1x1w[1], 0);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 8x1x1w FPX SQR R[1]: ", r64, NWORDS_FIELD);

  // 4x2x1w Fp2 mul
  fp2mul_mont_4x2x1w(r_4x2x1w, a_4x2x1w, a_4x2x1w);

  get_channel_8x1w(r51, r_4x2x1w, 0);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 4x2x1w FPX MUL R[0]: ", r64, NWORDS_FIELD);
  get_channel_8x1w(r51, r_4x2x1w, 1);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 4x2x1w FPX MUL R[1]: ", r64, NWORDS_FIELD);

  // 4x2x1w Fp2 sqr
  fp2sqr_mont_4x2x1w(r_4x2x1w, a_4x2x1w);

  get_channel_8x1w(r51, r_4x2x1w, 0);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 4x2x1w FPX SQR R[0]: ", r64, NWORDS_FIELD);
  get_channel_8x1w(r51, r_4x2x1w, 1);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 4x2x1w FPX SQR R[1]: ", r64, NWORDS_FIELD);

  // 4x1x2w Fp2 mul
  fp2mul_mont_4x1x2w(r_4x1x2w, a_4x1x2w, a_4x1x2w);

  get_channel_4x2w(r51, r_4x1x2w[0], 0);
  carryp(r51);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 4x1x2w FPX MUL R[0]: ", r64, NWORDS_FIELD);
  get_channel_4x2w(r51, r_4x1x2w[1], 0);
  carryp(r51);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 4x1x2w FPX MUL R[1]: ", r64, NWORDS_FIELD);

  // 4x1x2w Fp2 sqr
  fp2sqr_mont_4x1x2w(r_4x1x2w, a_4x1x2w);

  get_channel_4x2w(r51, r_4x1x2w[0], 0);
  carryp(r51);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 4x1x2w FPX SQR R[0]: ", r64, NWORDS_FIELD);
  get_channel_4x2w(r51, r_4x1x2w[1], 0);
  carryp(r51);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 4x1x2w FPX SQR R[1]: ", r64, NWORDS_FIELD);

  // 2x4x1w Fp2 mul 
  fp2mul_mont_2x4x1w(r_2x4x1w, a_2x4x1w, a_2x4x1w);

  get_channel_8x1w(r51, r_2x4x1w, 0);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 2x4x1w FPX MUL R[0]: ", r64, NWORDS_FIELD); 
  get_channel_8x1w(r51, r_2x4x1w, 2);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 2x4x1w FPX MUL R[1]: ", r64, NWORDS_FIELD);   

  // 2x2x2w Fp2 mul
  fp2mul_mont_2x2x2w(r_2x2x2w, a_2x2x2w, a_2x2x2w);

  get_channel_4x2w(r51, r_2x2x2w, 0);
  carryp(r51);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 2x2x2w FPX MUL R[0]: ", r64, NWORDS_FIELD);
  get_channel_4x2w(r51, r_2x2x2w, 2);
  carryp(r51);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 2x2x2w FPX MUL R[1]: ", r64, NWORDS_FIELD);

  // 2x2x2w Fp2 sqr
  fp2sqr_mont_2x2x2w(r_2x2x2w, a_2x2x2w);

  get_channel_4x2w(r51, r_2x2x2w, 0);
  carryp(r51);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 2x2x2w FPX SQR R[0]: ", r64, NWORDS_FIELD);
  get_channel_4x2w(r51, r_2x2x2w, 2);
  carryp(r51);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 2x2x2w FPX SQR R[1]: ", r64, NWORDS_FIELD);

  // ---------------------------------------------------------------------------
  // measure timings 
  uint64_t start_cycles, end_cycles, diff_cycles;

  puts("\nTIMINGS: ");

  LOAD_CACHE(fp2mul_mont_8x1x1w(r_8x1x1w, a_8x1x1w, a_8x1x1w), ITER_M);
  MEASURE_TIME(fp2mul_mont_8x1x1w(r_8x1x1w, a_8x1x1w, a_8x1x1w), ITER_L);
  printf("* 8x1x1w FPX MUL: %ld\n", diff_cycles);

  LOAD_CACHE(fp2mul_mont_4x2x1w(r_4x2x1w, a_4x2x1w, a_4x2x1w), ITER_M);
  MEASURE_TIME(fp2mul_mont_4x2x1w(r_4x2x1w, a_4x2x1w, a_4x2x1w), ITER_L);
  printf("* 4x2x1w FPX MUL: %ld\n", diff_cycles);

  LOAD_CACHE(fp2mul_mont_4x1x2w(r_4x1x2w, a_4x1x2w, a_4x1x2w), ITER_M);
  MEASURE_TIME(fp2mul_mont_4x1x2w(r_4x1x2w, a_4x1x2w, a_4x1x2w), ITER_L);
  printf("* 4x1x2w FPX MUL: %ld\n", diff_cycles);

  LOAD_CACHE(fp2mul_mont_2x4x1w(r_2x4x1w, a_2x4x1w, a_2x4x1w), ITER_M);
  MEASURE_TIME(fp2mul_mont_2x4x1w(r_2x4x1w, a_2x4x1w, a_2x4x1w), ITER_L);
  printf("* 2x4x1w FPX MUL: %ld\n", diff_cycles);

  LOAD_CACHE(fp2mul_mont_2x2x2w(r_2x2x2w, a_2x2x2w, a_2x2x2w), ITER_M);
  MEASURE_TIME(fp2mul_mont_2x2x2w(r_2x2x2w, a_2x2x2w, a_2x2x2w), ITER_L);
  printf("* 2x2x2w FPX MUL: %ld\n", diff_cycles);

  puts("");

  LOAD_CACHE(fp2sqr_mont_8x1x1w(r_8x1x1w, a_8x1x1w), ITER_M);
  MEASURE_TIME(fp2sqr_mont_8x1x1w(r_8x1x1w, a_8x1x1w), ITER_L);
  printf("* 8x1x1w FPX SQR: %ld\n", diff_cycles);

  LOAD_CACHE(fp2sqr_mont_4x2x1w(r_4x2x1w, a_4x2x1w), ITER_M);
  MEASURE_TIME(fp2sqr_mont_4x2x1w(r_4x2x1w, a_4x2x1w), ITER_L);
  printf("* 4x2x1w FPX SQR: %ld\n", diff_cycles);

  LOAD_CACHE(fp2sqr_mont_4x1x2w(r_4x1x2w, a_4x1x2w), ITER_M);
  MEASURE_TIME(fp2sqr_mont_4x1x2w(r_4x1x2w, a_4x1x2w), ITER_L);
  printf("* 4x1x2w FPX SQR: %ld\n", diff_cycles);

  LOAD_CACHE(fp2sqr_mont_2x2x2w(r_2x2x2w, a_2x2x2w), ITER_M);
  MEASURE_TIME(fp2sqr_mont_2x2x2w(r_2x2x2w, a_2x2x2w), ITER_L);
  printf("* 2x2x2w FPX SQR: %ld\n", diff_cycles);
}

void test_curve()
{
  vpoint_proj P_8x1x1x1w, Q_8x1x1x1w, PQ_8x1x1x1w;
  vf2elm_t A_8x1x1x1w, A24_8x1x1x1w, C24_8x1x1x1w, A1_8x1x1x1w, A2_8x1x1x1w;
  vfelm_t x3z3x2z2, z1x1_A;
  vf2elm_t x3z3x2z2_2x4x1x1w, z1x1_A_2x4x1x1w;
  vfelm_t P_1x2x4x1w, Q_1x2x4x1w, C24_A1_1x2x4x1w;
  vfelm_t P_2x2x2x1w, Q_2x2x2x1w, C24_A1_2x2x2x1w;
  vgelm_t P_1x2x2x2w, Q_1x2x2x2w, C24_A1_1x2x2x2w;
  uint64_t xP0[VNWORDS], xP1[VNWORDS], zP0[VNWORDS], zP1[VNWORDS];
  uint64_t xQ0[VNWORDS], xQ1[VNWORDS], zQ0[VNWORDS], zQ1[VNWORDS];
  uint64_t xPQ0[VNWORDS], xPQ1[VNWORDS], zPQ0[VNWORDS], zPQ1[VNWORDS];
  uint64_t r51[VNWORDS], r64[VNWORDS];
  int i;

  uint64_t vmont_Rx4[VNWORDS] = {
    0x00000000003F8, 0x0000000000000, 0x0000000000000, 0x0000000000000, 
    0x6000000000000, 0x6EC9191E174AA, 0x5787F4C4B1DB4, 0x14F30936824A2, 
    0x771EB133E84F5, 0x00CEC5480F20E, };

  uint64_t vmont_Rx6[VNWORDS] = {
    0x00000000005F6, 0x0000000000000, 0x0000000000000, 0x0000000000000, 
    0x3800000000000, 0x0F371D28A907D, 0x3F6FFB399230B, 0x67C7BDB627787, 
    0x6EE2011504BBA, 0x00348C1710ACE, };

  puts("\nCORRECTNESS TEST: ");

  // ---------------------------------------------------------------------------
  // 8x1x1x1w xDBLADD

  for (i = 0; i < VNWORDS; i++) {
    xP0[i]  = vA_gen[i];
    xP1[i]  = vA_gen[i+VNWORDS];
    xQ0[i]  = vA_gen[i+2*VNWORDS];
    xQ1[i]  = vA_gen[i+3*VNWORDS];
    xPQ0[i] = vA_gen[i+4*VNWORDS];
    xPQ1[i] = vA_gen[i+5*VNWORDS];
    zP0[i] = zQ0[i] = zPQ0[i] = vmont_R[i];
    zP1[i] = zQ1[i] = zPQ1[i] = 0;
  }

  // initialize vectors for 8x1x1x1w 
  for (i = 0; i < VNWORDS; i++) {
    P_8x1x1x1w.X[0][i] = VSET1(xP0[i]);
    P_8x1x1x1w.X[1][i] = VSET1(xP1[i]);
    P_8x1x1x1w.Z[0][i] = VSET1(zP0[i]);
    P_8x1x1x1w.Z[1][i] = VSET1(zP1[i]);
    Q_8x1x1x1w.X[0][i] = VSET1(xQ0[i]);
    Q_8x1x1x1w.X[1][i] = VSET1(xQ1[i]);
    Q_8x1x1x1w.Z[0][i] = VSET1(zQ0[i]);
    Q_8x1x1x1w.Z[1][i] = VSET1(zQ1[i]);
    PQ_8x1x1x1w.X[0][i] = VSET1(xPQ0[i]);
    PQ_8x1x1x1w.X[1][i] = VSET1(xPQ1[i]);
    PQ_8x1x1x1w.Z[0][i] = VSET1(zPQ0[i]);
    PQ_8x1x1x1w.Z[1][i] = VSET1(zPQ1[i]);
    A_8x1x1x1w[0][i] = VZERO;
    A_8x1x1x1w[1][i] = VZERO;
  }

  // initializing A24
  for (i = 0; i < VNWORDS; i++) {
    A1_8x1x1x1w[0][i] = A24_8x1x1x1w[0][i] = VSET1(vmont_R[i]);
    A1_8x1x1x1w[1][i] = A24_8x1x1x1w[1][i] = VZERO;
  }
  fp2add_8x1x1w(A24_8x1x1x1w, A24_8x1x1x1w, A24_8x1x1x1w);       // A24 = 2
  fp2mul_mont_8x1x1w(C24_8x1x1x1w, A24_8x1x1x1w, A24_8x1x1x1w);  // C24 = 4;
  fp2add_8x1x1w(A24_8x1x1x1w, A_8x1x1x1w, A24_8x1x1x1w);         // A24 = A+2
  fp2div2_8x1x1w(A24_8x1x1x1w, A24_8x1x1x1w);                    // A24 = (A+2)/2
  fp2div2_8x1x1w(A24_8x1x1x1w, A24_8x1x1x1w);                    // A24 = (A+2)/4

  xDBLADD_8x1x1x1w(&P_8x1x1x1w, &Q_8x1x1x1w, PQ_8x1x1x1w.X, PQ_8x1x1x1w.Z, A24_8x1x1x1w);

  puts("\n8x1x1x1w xDBLADD: ");
  get_channel_8x1w(r51, P_8x1x1x1w.X[0], 0);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 8x1x1x1w XP0: ", r64, NWORDS_FIELD);
  get_channel_8x1w(r51, P_8x1x1x1w.X[1], 0);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 8x1x1x1w XP1: ", r64, NWORDS_FIELD);
  get_channel_8x1w(r51, P_8x1x1x1w.Z[0], 0);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 8x1x1x1w ZP0: ", r64, NWORDS_FIELD);
  get_channel_8x1w(r51, P_8x1x1x1w.Z[1], 0);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 8x1x1x1w ZP1: ", r64, NWORDS_FIELD);
  get_channel_8x1w(r51, Q_8x1x1x1w.X[0], 0);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 8x1x1x1w XQ0: ", r64, NWORDS_FIELD);
  get_channel_8x1w(r51, Q_8x1x1x1w.X[1], 0);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 8x1x1x1w XQ1: ", r64, NWORDS_FIELD);
  get_channel_8x1w(r51, Q_8x1x1x1w.Z[0], 0);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 8x1x1x1w ZQ0: ", r64, NWORDS_FIELD);
  get_channel_8x1w(r51, Q_8x1x1x1w.Z[1], 0);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 8x1x1x1w ZQ1: ", r64, NWORDS_FIELD);

  // 2x4x1x1w xDBLADD
  // initialize vectors for 2x4x1x1w 
  for (i = 0; i < VNWORDS; i++) {
    x3z3x2z2_2x4x1x1w[0][i] = VSET(xQ0[i], zQ0[i], xP0[i], zP0[i], xQ0[i], zQ0[i], xP0[i], zP0[i]);
    x3z3x2z2_2x4x1x1w[1][i] = VSET(xQ1[i], zQ1[i], xP1[i], zP1[i], xQ1[i], zQ1[i], xP1[i], zP1[i]);
    z1x1_A_2x4x1x1w[0][i] = VSET(zPQ0[i], xPQ0[i], 0, 0, zPQ0[i], xPQ0[i], 0, 0);
    z1x1_A_2x4x1x1w[1][i] = VSET(zPQ1[i], xPQ1[i], 0, 0, zPQ1[i], xPQ1[i], 0, 0);
  }

  xDBLADD_2x4x1x1w(x3z3x2z2_2x4x1x1w, z1x1_A_2x4x1x1w);

  puts("\n2x4x1x1w xDBLADD: ");
  get_channel_8x1w(r51, x3z3x2z2_2x4x1x1w[0], 1);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 1x4x2x1w XP0: ", r64, NWORDS_FIELD);
  get_channel_8x1w(r51, x3z3x2z2_2x4x1x1w[1], 1);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 1x4x2x1w XP1: ", r64, NWORDS_FIELD);
  get_channel_8x1w(r51, x3z3x2z2_2x4x1x1w[0], 0);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 1x4x2x1w ZP0: ", r64, NWORDS_FIELD);
  get_channel_8x1w(r51, x3z3x2z2_2x4x1x1w[1], 0);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 1x4x2x1w ZP1: ", r64, NWORDS_FIELD);
  get_channel_8x1w(r51, x3z3x2z2_2x4x1x1w[0], 3);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 1x4x2x1w XQ0: ", r64, NWORDS_FIELD);
  get_channel_8x1w(r51, x3z3x2z2_2x4x1x1w[1], 3);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 1x4x2x1w XQ1: ", r64, NWORDS_FIELD);
  get_channel_8x1w(r51, x3z3x2z2_2x4x1x1w[0], 2);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 1x4x2x1w ZQ0: ", r64, NWORDS_FIELD);
  get_channel_8x1w(r51, x3z3x2z2_2x4x1x1w[1], 2);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 1x4x2x1w ZQ1: ", r64, NWORDS_FIELD);

  // 1x4x2x1w xDBLADD

  // initialize vectors for 1x4x2x1w 
  for (i = 0; i < VNWORDS; i++) {
    x3z3x2z2[i] = VSET(xQ1[i], xQ0[i], zQ1[i], zQ0[i], xP1[i], xP0[i], zP1[i], zP0[i]);
    z1x1_A[i] = VSET(zPQ1[i], zPQ0[i], xPQ1[i], xPQ0[i], 0, 0, 0, 0);
  }

  xDBLADD_1x4x2x1w(x3z3x2z2, z1x1_A);

  puts("\n1x4x2x1w xDBLADD: ");
  get_channel_8x1w(r51, x3z3x2z2, 2);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 1x4x2x1w XP0: ", r64, NWORDS_FIELD);
  get_channel_8x1w(r51, x3z3x2z2, 3);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 1x4x2x1w XP1: ", r64, NWORDS_FIELD);
  get_channel_8x1w(r51, x3z3x2z2, 0);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 1x4x2x1w ZP0: ", r64, NWORDS_FIELD);
  get_channel_8x1w(r51, x3z3x2z2, 1);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 1x4x2x1w ZP1: ", r64, NWORDS_FIELD);
  get_channel_8x1w(r51, x3z3x2z2, 6);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 1x4x2x1w XQ0: ", r64, NWORDS_FIELD);
  get_channel_8x1w(r51, x3z3x2z2, 7);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 1x4x2x1w XQ1: ", r64, NWORDS_FIELD);
  get_channel_8x1w(r51, x3z3x2z2, 4);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 1x4x2x1w ZQ0: ", r64, NWORDS_FIELD);
  get_channel_8x1w(r51, x3z3x2z2, 5);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 1x4x2x1w ZQ1: ", r64, NWORDS_FIELD);

  // ---------------------------------------------------------------------------
  // xDBL 8x1x1x1w 

  // initialize vectors for 8x1x1x1w 
  for (i = 0; i < VNWORDS; i++) {
    P_8x1x1x1w.X[0][i] = VSET1(xP0[i]);
    P_8x1x1x1w.X[1][i] = VSET1(xP1[i]);
    P_8x1x1x1w.Z[0][i] = VSET1(zP0[i]);
    P_8x1x1x1w.Z[1][i] = VSET1(zP1[i]);
    Q_8x1x1x1w.X[0][i] = VZERO;
    Q_8x1x1x1w.X[1][i] = VZERO;
    Q_8x1x1x1w.Z[0][i] = VZERO;
    Q_8x1x1x1w.Z[1][i] = VZERO;
  }

  xDBL_8x1x1x1w(&P_8x1x1x1w, &Q_8x1x1x1w, A1_8x1x1x1w, C24_8x1x1x1w);      // A1 = 1, C24 = 4

  puts("\n8x1x1x1w xDBL: ");
  get_channel_8x1w(r51, Q_8x1x1x1w.X[0], 0);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 8x1x1x1w XQ0: ", r64, NWORDS_FIELD);
  get_channel_8x1w(r51, Q_8x1x1x1w.X[1], 0);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 8x1x1x1w XQ1: ", r64, NWORDS_FIELD);
  get_channel_8x1w(r51, Q_8x1x1x1w.Z[0], 0);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 8x1x1x1w ZQ0: ", r64, NWORDS_FIELD);
  get_channel_8x1w(r51, Q_8x1x1x1w.Z[1], 0);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 8x1x1x1w ZQ1: ", r64, NWORDS_FIELD);

  // xDBL_2x2x2x1w 
  // initialize vectors for 2x2x2x1w
  for (i = 0; i < VNWORDS; i++) {
    P_2x2x2x1w[i] = VSET(xP1[i], xP0[i], zP1[i], zP0[i], xP1[i], xP0[i], zP1[i], zP0[i]);
    C24_A1_2x2x2x1w[i] = VSET(0, vmont_Rx4[i], 0, vmont_R[i], 0, vmont_Rx4[i], 0, vmont_R[i]);
  }

  xDBL_2x2x2x1w(P_2x2x2x1w, Q_2x2x2x1w, C24_A1_2x2x2x1w);

  puts("\n2x2x2x1w xDBL: ");  
  get_channel_8x1w(r51, Q_2x2x2x1w, 2);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 2x2x2x1w XQ0: ", r64, NWORDS_FIELD);
  get_channel_8x1w(r51, Q_2x2x2x1w, 3);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 2x2x2x1w XQ1: ", r64, NWORDS_FIELD);
  get_channel_8x1w(r51, Q_2x2x2x1w, 0);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 2x2x2x1w ZQ0: ", r64, NWORDS_FIELD);
  get_channel_8x1w(r51, Q_2x2x2x1w, 1);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 2x2x2x1w ZQ1: ", r64, NWORDS_FIELD);

  // xDBL 1x2x4x1w

  // initialize vectors for 1x4x2x1w 
  for (i = 0; i < VNWORDS; i++) {
    P_1x2x4x1w[i]      = VSET(0, xP1[i], 0, xP0[i], 0, zP1[i], 0, zP0[i]);
    Q_1x2x4x1w[i]      = VZERO;
    C24_A1_1x2x4x1w[i] = VSET(0, 0, 0, vmont_Rx4[i], 0, 0, 0, vmont_R[i]);
  }

  xDBL_1x2x4x1w(P_1x2x4x1w, Q_1x2x4x1w, C24_A1_1x2x4x1w);

  puts("\n1x2x4x1w xDBL: ");
  get_channel_8x1w(r51, Q_1x2x4x1w, 4);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 1x2x4x1w XQ0: ", r64, NWORDS_FIELD);
  get_channel_8x1w(r51, Q_1x2x4x1w, 6);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 1x2x4x1w XQ1: ", r64, NWORDS_FIELD);
  get_channel_8x1w(r51, Q_1x2x4x1w, 0);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 1x2x4x1w ZQ0: ", r64, NWORDS_FIELD);
  get_channel_8x1w(r51, Q_1x2x4x1w, 2);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 1x2x4x1w ZQ1: ", r64, NWORDS_FIELD);

  // xDBL_1x2x2x2w

  // initialize vectors for 1x2x2x2w
  for (i = 0; i < VGWORDS; i++) {
    P_1x2x2x2w[i]      = VSET(xP1[i+VGWORDS], xP1[i], xP0[i+VGWORDS], xP0[i], zP1[i+VGWORDS], zP1[i], zP0[i+VGWORDS], zP0[i]);
    Q_1x2x2x2w[i]      = VZERO;
    C24_A1_1x2x2x2w[i] = VSET(0, 0, vmont_Rx4[i+VGWORDS], vmont_Rx4[i], 0, 0, vmont_R[i+VGWORDS], vmont_R[i]);
  } 

  xDBL_1x2x2x2w(P_1x2x2x2w, Q_1x2x2x2w, C24_A1_1x2x2x2w);

  puts("\n1x2x2x2w xDBL: ");
  get_channel_4x2w(r51, Q_1x2x2x2w, 4);
  carryp(r51);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 1x2x2x2w XQ0: ", r64, NWORDS_FIELD);
  get_channel_4x2w(r51, Q_1x2x2x2w, 6);
  carryp(r51);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 1x2x2x2w XQ1: ", r64, NWORDS_FIELD);
  get_channel_4x2w(r51, Q_1x2x2x2w, 0);
  carryp(r51);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 1x2x2x2w ZQ0: ", r64, NWORDS_FIELD);
  get_channel_4x2w(r51, Q_1x2x2x2w, 2);
  carryp(r51);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 1x2x2x2w ZQ1: ", r64, NWORDS_FIELD);

  // ---------------------------------------------------------------------------
  // xTPL 8x1x1x1w 

  // initialize vectors for 8x1x1x1w 
  for (i = 0; i < VNWORDS; i++) {
    P_8x1x1x1w.X[0][i] = VSET1(xP0[i]);
    P_8x1x1x1w.X[1][i] = VSET1(xP1[i]);
    P_8x1x1x1w.Z[0][i] = VSET1(zP0[i]);
    P_8x1x1x1w.Z[1][i] = VSET1(zP1[i]);
    Q_8x1x1x1w.X[0][i] = VZERO;
    Q_8x1x1x1w.X[1][i] = VZERO;
    Q_8x1x1x1w.Z[0][i] = VZERO;
    Q_8x1x1x1w.Z[1][i] = VZERO;
  }

  // initializing A24
  for (i = 0; i < VNWORDS; i++) {
    A1_8x1x1x1w[0][i] = A2_8x1x1x1w[0][i] = VSET1(vmont_R[i]);
    A1_8x1x1x1w[1][i] = A2_8x1x1x1w[1][i] = VZERO;
  }
  fp2add_8x1x1w(A2_8x1x1x1w, A2_8x1x1x1w, A2_8x1x1x1w);         // A2 = 2
  fp2add_8x1x1w(A1_8x1x1x1w, A1_8x1x1x1w, A1_8x1x1x1w);         // A1 = 2
  fp2add_8x1x1w(A1_8x1x1x1w, A1_8x1x1x1w, A1_8x1x1x1w);         // A1 = 4
  fp2add_8x1x1w(A1_8x1x1x1w, A1_8x1x1x1w, A2_8x1x1x1w);         // A1 = 6

  xTPL_8x1x1x1w(&P_8x1x1x1w, &Q_8x1x1x1w, A2_8x1x1x1w, A1_8x1x1x1w);// A2 = 2, A1 = 6, C24 = 4

  puts("\n8x1x1x1w xTPL: ");
  get_channel_8x1w(r51, Q_8x1x1x1w.X[0], 0);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 8x1x1x1w XQ0: ", r64, NWORDS_FIELD);
  get_channel_8x1w(r51, Q_8x1x1x1w.X[1], 0);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 8x1x1x1w XQ1: ", r64, NWORDS_FIELD);
  get_channel_8x1w(r51, Q_8x1x1x1w.Z[0], 0);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 8x1x1x1w ZQ0: ", r64, NWORDS_FIELD);
  get_channel_8x1w(r51, Q_8x1x1x1w.Z[1], 0);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 8x1x1x1w ZQ1: ", r64, NWORDS_FIELD);

  // xTPL 1x2x4x1w

  // initialize vectors for 1x2x4x1w 
  for (i = 0; i < VNWORDS; i++) {
    P_1x2x4x1w[i]      = VSET(0, xP1[i], 0, xP0[i], 0, zP1[i], 0, zP0[i]);
    Q_1x2x4x1w[i]      = VZERO;
    C24_A1_1x2x4x1w[i] = VSET(0, 0, 0, vmont_Rx4[i], 0, 0, 0, vmont_Rx6[i]);
  }

  xTPL_1x2x4x1w(P_1x2x4x1w, Q_1x2x4x1w, C24_A1_1x2x4x1w);

  puts("\n1x2x4x1w xTPL: ");
  get_channel_8x1w(r51, Q_1x2x4x1w, 4);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 1x2x4x1w XQ0: ", r64, NWORDS_FIELD);
  get_channel_8x1w(r51, Q_1x2x4x1w, 6);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 1x2x4x1w XQ1: ", r64, NWORDS_FIELD);
  get_channel_8x1w(r51, Q_1x2x4x1w, 0);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 1x2x4x1w ZQ0: ", r64, NWORDS_FIELD);
  get_channel_8x1w(r51, Q_1x2x4x1w, 2);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 1x2x4x1w ZQ1: ", r64, NWORDS_FIELD);

  // xTPL_1x2x2x2w

  // initialize vectors for 1x2x2x2w
  for (i = 0; i < VGWORDS; i++) {
    P_1x2x2x2w[i]      = VSET(xP1[i+VGWORDS], xP1[i], xP0[i+VGWORDS], xP0[i], zP1[i+VGWORDS], zP1[i], zP0[i+VGWORDS], zP0[i]);
    Q_1x2x2x2w[i]      = VZERO;
    C24_A1_1x2x2x2w[i] = VSET(0, 0, vmont_Rx4[i+VGWORDS], vmont_Rx4[i], 0, 0, vmont_Rx6[i+VGWORDS], vmont_Rx6[i]);
  } 

  xTPL_1x2x2x2w(P_1x2x2x2w, Q_1x2x2x2w, C24_A1_1x2x2x2w);

  puts("\n1x2x2x2w xTPL: ");
  get_channel_4x2w(r51, Q_1x2x2x2w, 4);
  carryp(r51);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 1x2x2x2w XQ0: ", r64, NWORDS_FIELD);
  get_channel_4x2w(r51, Q_1x2x2x2w, 6);
  carryp(r51);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 1x2x2x2w XQ1: ", r64, NWORDS_FIELD);
  get_channel_4x2w(r51, Q_1x2x2x2w, 0);
  carryp(r51);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 1x2x2x2w ZQ0: ", r64, NWORDS_FIELD);
  get_channel_4x2w(r51, Q_1x2x2x2w, 2);
  carryp(r51);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 1x2x2x2w ZQ1: ", r64, NWORDS_FIELD);


  // ---------------------------------------------------------------------------
  // measure timings 
  uint64_t start_cycles, end_cycles, diff_cycles;

  puts("\nTIMINGS: ");

  LOAD_CACHE(xDBLADD_8x1x1x1w(&P_8x1x1x1w, &Q_8x1x1x1w, PQ_8x1x1x1w.X, PQ_8x1x1x1w.Z, A24_8x1x1x1w), ITER_S);
  MEASURE_TIME(xDBLADD_8x1x1x1w(&P_8x1x1x1w, &Q_8x1x1x1w, PQ_8x1x1x1w.X, PQ_8x1x1x1w.Z, A24_8x1x1x1w), ITER_L);
  printf("* 8x1x1x1w xDBLADD: %ld\n", diff_cycles);

  LOAD_CACHE(xDBLADD_2x4x1x1w(x3z3x2z2_2x4x1x1w, z1x1_A_2x4x1x1w), ITER_S);
  MEASURE_TIME(xDBLADD_2x4x1x1w(x3z3x2z2_2x4x1x1w, z1x1_A_2x4x1x1w), ITER_L);
  printf("* 2x4x1x1w xDBLADD: %ld\n", diff_cycles);

  LOAD_CACHE(xDBLADD_1x4x2x1w(x3z3x2z2, z1x1_A), ITER_S);
  MEASURE_TIME(xDBLADD_1x4x2x1w(x3z3x2z2, z1x1_A), ITER_L);
  printf("* 1x4x2x1w xDBLADD: %ld\n", diff_cycles);

  puts("");

  LOAD_CACHE(xDBL_8x1x1x1w(&P_8x1x1x1w, &Q_8x1x1x1w, A1_8x1x1x1w, C24_8x1x1x1w), ITER_S);
  MEASURE_TIME(xDBL_8x1x1x1w(&P_8x1x1x1w, &Q_8x1x1x1w, A1_8x1x1x1w, C24_8x1x1x1w), ITER_L);
  printf("* 8x1x1x1w xDBL: %ld\n", diff_cycles);

  LOAD_CACHE(xDBL_2x2x2x1w(P_2x2x2x1w, Q_2x2x2x1w, C24_A1_2x2x2x1w), ITER_S);
  MEASURE_TIME(xDBL_2x2x2x1w(P_2x2x2x1w, Q_2x2x2x1w, C24_A1_2x2x2x1w), ITER_L);
  printf("* 2x2x2x1w xDBL: %ld\n", diff_cycles);

  LOAD_CACHE(xDBL_1x2x4x1w(P_1x2x4x1w, Q_1x2x4x1w, C24_A1_1x2x4x1w), ITER_S);
  MEASURE_TIME(xDBL_1x2x4x1w(P_1x2x4x1w, Q_1x2x4x1w, C24_A1_1x2x4x1w), ITER_L);
  printf("* 1x2x4x1w xDBL: %ld\n", diff_cycles);

  LOAD_CACHE(xDBL_1x2x2x2w(P_1x2x2x2w, Q_1x2x2x2w, C24_A1_1x2x2x2w), ITER_S);
  MEASURE_TIME(xDBL_1x2x2x2w(P_1x2x2x2w, Q_1x2x2x2w, C24_A1_1x2x2x2w), ITER_L);
  printf("* 1x2x2x2w xDBL: %ld\n", diff_cycles);

  puts("");

  LOAD_CACHE(xTPL_8x1x1x1w(&P_8x1x1x1w, &Q_8x1x1x1w, A1_8x1x1x1w, A2_8x1x1x1w), ITER_S);
  MEASURE_TIME(xTPL_8x1x1x1w(&P_8x1x1x1w, &Q_8x1x1x1w, A1_8x1x1x1w, A2_8x1x1x1w), ITER_L);
  printf("* 8x1x1x1w xTPL: %ld\n", diff_cycles);

  LOAD_CACHE(xTPL_1x2x4x1w(P_1x2x4x1w, Q_1x2x4x1w, C24_A1_1x2x4x1w), ITER_S);
  MEASURE_TIME(xTPL_1x2x4x1w(P_1x2x4x1w, Q_1x2x4x1w, C24_A1_1x2x4x1w), ITER_L);
  printf("* 1x2x4x1w xTPL: %ld\n", diff_cycles);

  LOAD_CACHE(xTPL_1x2x2x2w(P_1x2x2x2w, Q_1x2x2x2w, C24_A1_1x2x2x2w), ITER_S);
  MEASURE_TIME(xTPL_1x2x2x2w(P_1x2x2x2w, Q_1x2x2x2w, C24_A1_1x2x2x2w), ITER_L);
  printf("* 1x2x2x2w xTPL: %ld\n", diff_cycles);
}

void test_isog()
{
  vpoint_proj P_8x1x1x1w, Q_8x1x1x1w;
  vf2elm_t A24plus_8x1x1x1w, A24minus_8x1x1x1w, C24_8x1x1x1w, coeff_8x1x1x1w[3];
  vf2elm_t P_4x2x1x1w, coeff_4x2x1x1w[3], Q_4x2x1x1w, coeff0_1_4x2x1x1w;
  vfelm_t P_2x2x2x1w, coeff_2x2x2x1w[3], Q_2x2x2x1w, coeff0_1_2x2x2x1w;
  vfelm_t C24_A24plus_2x2x2x1w, coeff__0_2x2x2x1w, coeff2_1_2x2x2x1w;
  vgelm_t P_1x2x2x2w, C24_A24plus_1x2x2x2w, coeff__0_1x2x2x2w, coeff2_1_1x2x2x2w, coeff_1x2x2x2w[3];
  vgelm_t Q_1x2x2x2w, coeff0_1_1x2x2x2w, coeff1_0_1x2x2x2w;
  uint64_t xP0[VNWORDS], xP1[VNWORDS], zP0[VNWORDS], zP1[VNWORDS], yP0[VNWORDS], yP1[VNWORDS];
  uint64_t r51[VNWORDS], r64[VNWORDS];
  int i;

  puts("\nCORRECTNESS TEST: ");

  for (i = 0; i < VNWORDS; i++) {
    xP0[i] = vA_gen[i];
    xP1[i] = vA_gen[i+VNWORDS];
    zP0[i] = vA_gen[i+2*VNWORDS];
    zP1[i] = vA_gen[i+3*VNWORDS];
    yP0[i] = vA_gen[i+4*VNWORDS];
    yP1[i] = vA_gen[i+5*VNWORDS];
  }

  // ---------------------------------------------------------------------------
  // get_4_isog 8x1x1x1w 

  // initialize vectors for 8x1x1x1w 
  for (i = 0; i < VNWORDS; i++) {
    P_8x1x1x1w.X[0][i] = VSET1(xP0[i]);
    P_8x1x1x1w.X[1][i] = VSET1(xP1[i]);
    P_8x1x1x1w.Z[0][i] = VSET1(zP0[i]);
    P_8x1x1x1w.Z[1][i] = VSET1(zP1[i]);
    A24plus_8x1x1x1w[0][i] = VZERO;
    A24plus_8x1x1x1w[1][i] = VZERO;
    C24_8x1x1x1w[0][i] = VZERO;
    C24_8x1x1x1w[1][i] = VZERO;
  }

  get_4_isog_8x1x1x1w(&P_8x1x1x1w, A24plus_8x1x1x1w, C24_8x1x1x1w, coeff_8x1x1x1w);

  puts("\n8x1x1x1w get_4_isog: ");
  get_channel_8x1w(r51, A24plus_8x1x1x1w[0], 0);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 8x1x1x1w A24+   0: ", r64, NWORDS_FIELD);
  get_channel_8x1w(r51, A24plus_8x1x1x1w[1], 0);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 8x1x1x1w A24+   1: ", r64, NWORDS_FIELD);
  get_channel_8x1w(r51, C24_8x1x1x1w[0], 0);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 8x1x1x1w C24    0: ", r64, NWORDS_FIELD);
  get_channel_8x1w(r51, C24_8x1x1x1w[1], 0);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 8x1x1x1w C24    1: ", r64, NWORDS_FIELD);
  get_channel_8x1w(r51, coeff_8x1x1x1w[0][0], 0);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 8x1x1x1w COEFF 00: ", r64, NWORDS_FIELD);
  get_channel_8x1w(r51, coeff_8x1x1x1w[0][1], 0);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 8x1x1x1w COEFF 01: ", r64, NWORDS_FIELD);
  get_channel_8x1w(r51, coeff_8x1x1x1w[1][0], 0);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 8x1x1x1w COEFF 10: ", r64, NWORDS_FIELD);
  get_channel_8x1w(r51, coeff_8x1x1x1w[1][1], 0);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 8x1x1x1w COEFF 11: ", r64, NWORDS_FIELD);
  get_channel_8x1w(r51, coeff_8x1x1x1w[2][0], 0);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 8x1x1x1w COEFF 20: ", r64, NWORDS_FIELD);
  get_channel_8x1w(r51, coeff_8x1x1x1w[2][1], 0);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 8x1x1x1w COEFF 21: ", r64, NWORDS_FIELD);

  // get_4_isog 2x2x2x1w 

  // initialize vectors for 2x2x2x1w
  for (i = 0; i < VNWORDS; i++) {
    P_2x2x2x1w[i] = VSET(xP1[i], xP0[i], zP1[i], zP0[i], xP1[i], xP0[i], zP1[i], zP0[i]);
    coeff__0_2x2x2x1w[i] = VZERO;
    coeff2_1_2x2x2x1w[i] = VZERO;
  }

  get_4_isog_2x2x2x1w(P_2x2x2x1w, C24_A24plus_2x2x2x1w, coeff__0_2x2x2x1w, coeff2_1_2x2x2x1w);

  puts("\n2x2x2x1w get_4_isog: ");  
  get_channel_8x1w(r51, C24_A24plus_2x2x2x1w, 0);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 2x2x2x1w A24+   0: ", r64, NWORDS_FIELD);
  get_channel_8x1w(r51, C24_A24plus_2x2x2x1w, 1);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 2x2x2x1w A24+   1: ", r64, NWORDS_FIELD);
  get_channel_8x1w(r51, C24_A24plus_2x2x2x1w, 2);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 2x2x2x1w C24    0: ", r64, NWORDS_FIELD);
  get_channel_8x1w(r51, C24_A24plus_2x2x2x1w, 3);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 2x2x2x1w C24    1: ", r64, NWORDS_FIELD);
  get_channel_8x1w(r51, coeff__0_2x2x2x1w, 0);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 2x2x2x1w COEFF 00: ", r64, NWORDS_FIELD);
  get_channel_8x1w(r51, coeff__0_2x2x2x1w, 1);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 2x2x2x1w COEFF 01: ", r64, NWORDS_FIELD);
  get_channel_8x1w(r51, coeff2_1_2x2x2x1w, 0);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 2x2x2x1w COEFF 10: ", r64, NWORDS_FIELD);
  get_channel_8x1w(r51, coeff2_1_2x2x2x1w, 1);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 2x2x2x1w COEFF 11: ", r64, NWORDS_FIELD);
  get_channel_8x1w(r51, coeff2_1_2x2x2x1w, 2);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 2x2x2x1w COEFF 20: ", r64, NWORDS_FIELD);
  get_channel_8x1w(r51, coeff2_1_2x2x2x1w, 3);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 2x2x2x1w COEFF 21: ", r64, NWORDS_FIELD);

  // get_4_isog 1x2x2x2w 

  // initialize vectors for 1x2x2x2w
  for (i = 0; i < VGWORDS; i++) {
    P_1x2x2x2w[i]           = VSET(xP1[i+VGWORDS], xP1[i], xP0[i+VGWORDS], xP0[i], zP1[i+VGWORDS], zP1[i], zP0[i+VGWORDS], zP0[i]);
    C24_A24plus_1x2x2x2w[i] = VZERO;
    coeff__0_1x2x2x2w[i]    = VZERO;
    coeff2_1_1x2x2x2w[i]    = VZERO;
  } 

  get_4_isog_1x2x2x2w(P_1x2x2x2w, C24_A24plus_1x2x2x2w, coeff__0_1x2x2x2w, coeff2_1_1x2x2x2w);

  puts("\n1x2x2x2w get_4_isog: ");
  get_channel_4x2w(r51, C24_A24plus_1x2x2x2w, 0);
  carryp(r51);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 1x2x2x2w A24+   0: ", r64, NWORDS_FIELD);
  get_channel_4x2w(r51, C24_A24plus_1x2x2x2w, 2);
  carryp(r51);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 1x2x2x2w A24+   1: ", r64, NWORDS_FIELD);
  get_channel_4x2w(r51, C24_A24plus_1x2x2x2w, 4);
  carryp(r51);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 1x2x2x2w C24    0: ", r64, NWORDS_FIELD);
  get_channel_4x2w(r51, C24_A24plus_1x2x2x2w, 6);
  carryp(r51);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 1x2x2x2w C24    1: ", r64, NWORDS_FIELD);
  get_channel_4x2w(r51, coeff__0_1x2x2x2w, 0);
  carryp(r51);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 1x2x2x2w COEFF 00: ", r64, NWORDS_FIELD);
  get_channel_4x2w(r51, coeff__0_1x2x2x2w, 2);
  carryp(r51);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 1x2x2x2w COEFF 01: ", r64, NWORDS_FIELD);
  get_channel_4x2w(r51, coeff2_1_1x2x2x2w, 0);
  carryp(r51);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 1x2x2x2w COEFF 10: ", r64, NWORDS_FIELD);
  get_channel_4x2w(r51, coeff2_1_1x2x2x2w, 2);
  carryp(r51);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 1x2x2x2w COEFF 11: ", r64, NWORDS_FIELD);
  get_channel_4x2w(r51, coeff2_1_1x2x2x2w, 4);
  carryp(r51);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 1x2x2x2w COEFF 20: ", r64, NWORDS_FIELD);
  get_channel_4x2w(r51, coeff2_1_1x2x2x2w, 6);
  carryp(r51);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 1x2x2x2w COEFF 21: ", r64, NWORDS_FIELD);

  // ---------------------------------------------------------------------------
  // eval_4_isog 8x1x1x1w 

  // initialize vectors for 8x1x1x1w 
  for (i = 0; i < VNWORDS; i++) {
    P_8x1x1x1w.X[0][i] = VSET1(xP0[i]);
    P_8x1x1x1w.X[1][i] = VSET1(xP1[i]);
    P_8x1x1x1w.Z[0][i] = VSET1(zP0[i]);
    P_8x1x1x1w.Z[1][i] = VSET1(zP1[i]);
    coeff_8x1x1x1w[0][0][i] = VSET1(xP0[i]);
    coeff_8x1x1x1w[0][1][i] = VSET1(xP1[i]);
    coeff_8x1x1x1w[1][0][i] = VSET1(zP0[i]);
    coeff_8x1x1x1w[1][1][i] = VSET1(zP1[i]);
    coeff_8x1x1x1w[2][0][i] = VSET1(yP0[i]);
    coeff_8x1x1x1w[2][1][i] = VSET1(yP1[i]);
  }

  eval_4_isog_8x1x1x1w(&P_8x1x1x1w, coeff_8x1x1x1w);

  puts("\n8x1x1x1w eval_4_isog: ");
  get_channel_8x1w(r51, P_8x1x1x1w.X[0], 0);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 8x1x1x1w XP0: ", r64, NWORDS_FIELD);
  get_channel_8x1w(r51, P_8x1x1x1w.X[1], 0);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 8x1x1x1w XP1: ", r64, NWORDS_FIELD);
  get_channel_8x1w(r51, P_8x1x1x1w.Z[0], 0);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 8x1x1x1w ZP0: ", r64, NWORDS_FIELD);
  get_channel_8x1w(r51, P_8x1x1x1w.Z[1], 0);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 8x1x1x1w ZP1: ", r64, NWORDS_FIELD);

  // eval_4_isog 4x2x1x1w

  // initialize vectors for 4x2x1x1w
  for (i = 0; i < VNWORDS; i++) {
    P_4x2x1x1w[0][i] = VSET(xP0[i], zP0[i], xP0[i], zP0[i], xP0[i], zP0[i], xP0[i], zP0[i]);
    P_4x2x1x1w[1][i] = VSET(xP1[i], zP1[i], xP1[i], zP1[i], xP1[i], zP1[i], xP1[i], zP1[i]);
    coeff_4x2x1x1w[0][0][i] = VSET(0, xP0[i], 0, xP0[i], 0, xP0[i], 0, xP0[i]);
    coeff_4x2x1x1w[0][1][i] = VSET(0, xP1[i], 0, xP1[i], 0, xP1[i], 0, xP1[i]);
    coeff_4x2x1x1w[1][0][i] = VSET(zP0[i], 0, zP0[i], 0, zP0[i], 0, zP0[i], 0);
    coeff_4x2x1x1w[1][1][i] = VSET(zP1[i], 0, zP1[i], 0, zP1[i], 0, zP1[i], 0);
    coeff_4x2x1x1w[2][0][i] = VSET(yP0[i], 0, yP0[i], 0, yP0[i], 0, yP0[i], 0);
    coeff_4x2x1x1w[2][1][i] = VSET(yP1[i], 0, yP1[i], 0, yP1[i], 0, yP1[i], 0);
  }

  eval_4_isog_4x2x1x1w(P_4x2x1x1w, coeff_4x2x1x1w);

  puts("\n4x2x1x1w eval_4_isog: ");  
  get_channel_8x1w(r51, P_4x2x1x1w[0], 1);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 4x2x1x1w XP0: ", r64, NWORDS_FIELD);
  get_channel_8x1w(r51, P_4x2x1x1w[1], 1);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 4x2x1x1w XP1: ", r64, NWORDS_FIELD);
  get_channel_8x1w(r51, P_4x2x1x1w[0], 0);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 4x2x1x1w ZP0: ", r64, NWORDS_FIELD);
  get_channel_8x1w(r51, P_4x2x1x1w[1], 0);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 4x2x1x1w ZP1: ", r64, NWORDS_FIELD);

  // eval_4_isog 2x2x2x1w 

  // initialize vectors for 2x2x2x1w
  for (i = 0; i < VNWORDS; i++) {
    P_2x2x2x1w[i] = VSET(xP1[i], xP0[i], zP1[i], zP0[i], xP1[i], xP0[i], zP1[i], zP0[i]);
    coeff_2x2x2x1w[0][i] = VSET(0, 0, xP1[i], xP0[i], 0, 0, xP1[i], xP0[i]);
    coeff_2x2x2x1w[1][i] = VSET(zP1[i], zP0[i], 0, 0, zP1[i], zP0[i], 0, 0);
    coeff_2x2x2x1w[2][i] = VSET(yP1[i], yP0[i], 0, 0, yP1[i], yP0[i], 0, 0);
  }

  eval_4_isog_2x2x2x1w(P_2x2x2x1w, coeff_2x2x2x1w);

  puts("\n2x2x2x1w eval_4_isog: ");  
  get_channel_8x1w(r51, P_2x2x2x1w, 2);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 2x2x2x1w XP0: ", r64, NWORDS_FIELD);
  get_channel_8x1w(r51, P_2x2x2x1w, 3);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 2x2x2x1w XP1: ", r64, NWORDS_FIELD);
  get_channel_8x1w(r51, P_2x2x2x1w, 0);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 2x2x2x1w ZP0: ", r64, NWORDS_FIELD);
  get_channel_8x1w(r51, P_2x2x2x1w, 1);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 2x2x2x1w ZP1: ", r64, NWORDS_FIELD);

  // eval_4_isog 1x2x2x2w 
  for (i = 0; i < VGWORDS; i++) {
    P_1x2x2x2w[i] = VSET(xP1[i+VGWORDS], xP1[i], xP0[i+VGWORDS], xP0[i], zP1[i+VGWORDS], zP1[i], zP0[i+VGWORDS], zP0[i]);
    coeff_1x2x2x2w[0][i] = VSET(0, 0, 0, 0, xP1[i+VGWORDS], xP1[i], xP0[i+VGWORDS], xP0[i]);
    coeff_1x2x2x2w[1][i] = VSET(zP1[i+VGWORDS], zP1[i], zP0[i+VGWORDS], zP0[i], 0, 0, 0, 0);
    coeff_1x2x2x2w[2][i] = VSET(yP1[i+VGWORDS], yP1[i], yP0[i+VGWORDS], yP0[i], 0, 0, 0, 0);
  }

  eval_4_isog_1x2x2x2w(P_1x2x2x2w, coeff_1x2x2x2w);

  puts("\n1x2x2x2w eval_4_isog: ");  
  get_channel_4x2w(r51, P_1x2x2x2w, 4);
  carryp(r51);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 1x2x2x2w XP0: ", r64, NWORDS_FIELD);
  get_channel_4x2w(r51, P_1x2x2x2w, 6);
  carryp(r51);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 1x2x2x2w XP1: ", r64, NWORDS_FIELD);
  get_channel_4x2w(r51, P_1x2x2x2w, 0);
  carryp(r51);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 1x2x2x2w ZP0: ", r64, NWORDS_FIELD);
  get_channel_4x2w(r51, P_1x2x2x2w, 2);
  carryp(r51);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 1x2x2x2w ZP1: ", r64, NWORDS_FIELD);

  // ---------------------------------------------------------------------------
  // eval_3_isog 8x1x1x1w 

  // initialize vectors for 8x1x1x1w 
  for (i = 0; i < VNWORDS; i++) {
    Q_8x1x1x1w.X[0][i] = VSET1(xP0[i]);
    Q_8x1x1x1w.X[1][i] = VSET1(xP1[i]);
    Q_8x1x1x1w.Z[0][i] = VSET1(zP0[i]);
    Q_8x1x1x1w.Z[1][i] = VSET1(zP1[i]);
    coeff_8x1x1x1w[0][0][i] = VSET1(yP0[i]);
    coeff_8x1x1x1w[0][1][i] = VSET1(yP1[i]);
    coeff_8x1x1x1w[1][0][i] = VSET1(xP0[i]);
    coeff_8x1x1x1w[1][1][i] = VSET1(xP1[i]);
  }

  eval_3_isog_8x1x1x1w(&Q_8x1x1x1w, coeff_8x1x1x1w);

  puts("\n8x1x1x1w eval_3_isog: ");
  get_channel_8x1w(r51, Q_8x1x1x1w.X[0], 0);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 8x1x1x1w XQ0: ", r64, NWORDS_FIELD);
  get_channel_8x1w(r51, Q_8x1x1x1w.X[1], 0);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 8x1x1x1w XQ1: ", r64, NWORDS_FIELD);
  get_channel_8x1w(r51, Q_8x1x1x1w.Z[0], 0);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 8x1x1x1w ZQ0: ", r64, NWORDS_FIELD);
  get_channel_8x1w(r51, Q_8x1x1x1w.Z[1], 0);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 8x1x1x1w ZQ1: ", r64, NWORDS_FIELD);

  // eval_3_isog 4x2x1x1w

  // initialize vectors for 4x2x1x1w
  for (i = 0; i < VNWORDS; i++) {
    Q_4x2x1x1w[0][i] = VSET(xP0[i], zP0[i], xP0[i], zP0[i], xP0[i], zP0[i], xP0[i], zP0[i]);
    Q_4x2x1x1w[1][i] = VSET(xP1[i], zP1[i], xP1[i], zP1[i], xP1[i], zP1[i], xP1[i], zP1[i]);
    coeff0_1_4x2x1x1w[0][i] = VSET(yP0[i], xP0[i], yP0[i], xP0[i], yP0[i], xP0[i], yP0[i], xP0[i]);
    coeff0_1_4x2x1x1w[1][i] = VSET(yP1[i], xP1[i], yP1[i], xP1[i], yP1[i], xP1[i], yP1[i], xP1[i]);
  }

  eval_3_isog_4x2x1x1w(Q_4x2x1x1w, coeff0_1_4x2x1x1w);

  puts("\n4x2x1x1w eval_3_isog: ");  
  get_channel_8x1w(r51, Q_4x2x1x1w[0], 1); 
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 4x2x1x1w XQ0: ", r64, NWORDS_FIELD);
  get_channel_8x1w(r51, Q_4x2x1x1w[1], 1); 
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 4x2x1x1w XQ1: ", r64, NWORDS_FIELD);
  get_channel_8x1w(r51, Q_4x2x1x1w[0], 0); 
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 4x2x1x1w ZQ0: ", r64, NWORDS_FIELD);
  get_channel_8x1w(r51, Q_4x2x1x1w[1], 0); 
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 4x2x1x1w ZQ1: ", r64, NWORDS_FIELD);

  // eval_3_isog 2x2x2x1w 

  // initialize vectors for 2x2x2x1w
  for (i = 0; i < VNWORDS; i++) {
    Q_2x2x2x1w[i] = VSET(xP1[i], xP0[i], zP1[i], zP0[i], xP1[i], xP0[i], zP1[i], zP0[i]);
    coeff0_1_2x2x2x1w[i] = VSET(yP1[i], yP0[i], xP1[i], xP0[i], yP1[i], yP0[i], xP1[i], xP0[i]);
  }

  eval_3_isog_2x2x2x1w(Q_2x2x2x1w, coeff0_1_2x2x2x1w);

  puts("\n2x2x2x1w eval_3_isog: ");  
  get_channel_8x1w(r51, Q_2x2x2x1w, 2);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 2x2x2x1w XQ0: ", r64, NWORDS_FIELD);
  get_channel_8x1w(r51, Q_2x2x2x1w, 3);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 2x2x2x1w XQ1: ", r64, NWORDS_FIELD);
  get_channel_8x1w(r51, Q_2x2x2x1w, 0);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 2x2x2x1w ZQ0: ", r64, NWORDS_FIELD);
  get_channel_8x1w(r51, Q_2x2x2x1w, 1);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 2x2x2x1w ZQ1: ", r64, NWORDS_FIELD);

  // eval_3_isog 1x2x2x2w 
  for (i = 0; i < VGWORDS; i++) {
    Q_1x2x2x2w[i] = VSET(xP1[i+VGWORDS], xP1[i], xP0[i+VGWORDS], xP0[i], zP1[i+VGWORDS], zP1[i], zP0[i+VGWORDS], zP0[i]);
    coeff0_1_1x2x2x2w[i] = VSET(yP1[i+VGWORDS], yP1[i], yP0[i+VGWORDS], yP0[i], xP1[i+VGWORDS], xP1[i], xP0[i+VGWORDS], xP0[i]);
  }

  eval_3_isog_1x2x2x2w(Q_1x2x2x2w, coeff0_1_1x2x2x2w);

  puts("\n1x2x2x2w eval_3_isog: ");  
  get_channel_4x2w(r51, Q_1x2x2x2w, 4);
  carryp(r51);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 1x2x2x2w XQ0: ", r64, NWORDS_FIELD);
  get_channel_4x2w(r51, Q_1x2x2x2w, 6);
  carryp(r51);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 1x2x2x2w XQ1: ", r64, NWORDS_FIELD);
  get_channel_4x2w(r51, Q_1x2x2x2w, 0);
  carryp(r51);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 1x2x2x2w ZQ0: ", r64, NWORDS_FIELD);
  get_channel_4x2w(r51, Q_1x2x2x2w, 2);
  carryp(r51);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 1x2x2x2w ZQ1: ", r64, NWORDS_FIELD);

  // ---------------------------------------------------------------------------
  // get_3_isog 8x1x1x1w 

  // initialize vectors for 8x1x1x1w 
  for (i = 0; i < VNWORDS; i++) {
    P_8x1x1x1w.X[0][i] = VSET1(xP0[i]);
    P_8x1x1x1w.X[1][i] = VSET1(xP1[i]);
    P_8x1x1x1w.Z[0][i] = VSET1(zP0[i]);
    P_8x1x1x1w.Z[1][i] = VSET1(zP1[i]);
    A24plus_8x1x1x1w[0][i] = VZERO;
    A24plus_8x1x1x1w[1][i] = VZERO;
    C24_8x1x1x1w[0][i] = VZERO;
    C24_8x1x1x1w[1][i] = VZERO;
  }

  get_3_isog_8x1x1x1w(&P_8x1x1x1w, C24_8x1x1x1w, A24plus_8x1x1x1w, coeff_8x1x1x1w);

  fp2sub_8x1x1w(C24_8x1x1x1w, A24plus_8x1x1x1w, C24_8x1x1x1w);  // C24 = A24plus - A24minus

  puts("\n8x1x1x1w get_3_isog: ");
  get_channel_8x1w(r51, A24plus_8x1x1x1w[0], 0);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 8x1x1x1w A24+   0: ", r64, NWORDS_FIELD);
  get_channel_8x1w(r51, A24plus_8x1x1x1w[1], 0);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 8x1x1x1w A24+   1: ", r64, NWORDS_FIELD);
  get_channel_8x1w(r51, C24_8x1x1x1w[0], 0);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 8x1x1x1w C24    0: ", r64, NWORDS_FIELD);
  get_channel_8x1w(r51, C24_8x1x1x1w[1], 0);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 8x1x1x1w C24    1: ", r64, NWORDS_FIELD);
  get_channel_8x1w(r51, coeff_8x1x1x1w[0][0], 0);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 8x1x1x1w COEFF 00: ", r64, NWORDS_FIELD);
  get_channel_8x1w(r51, coeff_8x1x1x1w[0][1], 0);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 8x1x1x1w COEFF 01: ", r64, NWORDS_FIELD);
  get_channel_8x1w(r51, coeff_8x1x1x1w[1][0], 0);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 8x1x1x1w COEFF 10: ", r64, NWORDS_FIELD);
  get_channel_8x1w(r51, coeff_8x1x1x1w[1][1], 0);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 8x1x1x1w COEFF 11: ", r64, NWORDS_FIELD);

  // ---------------------------------------------------------------------------
  // get_3_isog 1x2x2x2w 

  // initialize vectors for 1x2x2x2w
  for (i = 0; i < VGWORDS; i++) {
    P_1x2x2x2w[i]           = VSET(xP1[i+VGWORDS], xP1[i], xP0[i+VGWORDS], xP0[i], zP1[i+VGWORDS], zP1[i], zP0[i+VGWORDS], zP0[i]);
    C24_A24plus_1x2x2x2w[i] = VZERO;
    coeff1_0_1x2x2x2w[i]    = VZERO;
  } 

  get_3_isog_1x2x2x2w(P_1x2x2x2w, C24_A24plus_1x2x2x2w, coeff1_0_1x2x2x2w);

  puts("\n1x2x2x2w get_3_isog: ");
  get_channel_4x2w(r51, C24_A24plus_1x2x2x2w, 0);
  carryp(r51);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 1x2x2x2w A24+   0: ", r64, NWORDS_FIELD);
  get_channel_4x2w(r51, C24_A24plus_1x2x2x2w, 2);
  carryp(r51);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 1x2x2x2w A24+   1: ", r64, NWORDS_FIELD);
  get_channel_4x2w(r51, C24_A24plus_1x2x2x2w, 4);
  carryp(r51);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 1x2x2x2w C24    0: ", r64, NWORDS_FIELD);
  get_channel_4x2w(r51, C24_A24plus_1x2x2x2w, 6);
  carryp(r51);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 1x2x2x2w C24    1: ", r64, NWORDS_FIELD);
  get_channel_4x2w(r51, coeff1_0_1x2x2x2w, 0);
  carryp(r51);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 1x2x2x2w COEFF 00: ", r64, NWORDS_FIELD);
  get_channel_4x2w(r51, coeff1_0_1x2x2x2w, 2);
  carryp(r51);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 1x2x2x2w COEFF 01: ", r64, NWORDS_FIELD);
  get_channel_4x2w(r51, coeff1_0_1x2x2x2w, 4);
  carryp(r51);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 1x2x2x2w COEFF 10: ", r64, NWORDS_FIELD);
  get_channel_4x2w(r51, coeff1_0_1x2x2x2w, 6);
  carryp(r51);
  mpi_conv_51to64(r64, r51, NWORDS_FIELD, VNWORDS);
  mpi_print("* 1x2x2x2w COEFF 11: ", r64, NWORDS_FIELD);

  // ---------------------------------------------------------------------------
  // measure timings 
  uint64_t start_cycles, end_cycles, diff_cycles;

  puts("\nTIMINGS: ");

  LOAD_CACHE(get_4_isog_8x1x1x1w(&P_8x1x1x1w, A24plus_8x1x1x1w, C24_8x1x1x1w, coeff_8x1x1x1w), ITER_S);
  MEASURE_TIME(get_4_isog_8x1x1x1w(&P_8x1x1x1w, A24plus_8x1x1x1w, C24_8x1x1x1w, coeff_8x1x1x1w), ITER_L);
  printf("* 8x1x1x1w get_4_isog: %ld\n", diff_cycles);

  LOAD_CACHE(get_4_isog_2x2x2x1w(P_2x2x2x1w, C24_A24plus_2x2x2x1w, coeff__0_2x2x2x1w, coeff2_1_2x2x2x1w), ITER_S);
  MEASURE_TIME(get_4_isog_2x2x2x1w(P_2x2x2x1w, C24_A24plus_2x2x2x1w, coeff__0_2x2x2x1w, coeff2_1_2x2x2x1w), ITER_L);
  printf("* 2x2x2x1w get_4_isog: %ld\n", diff_cycles);

  LOAD_CACHE(get_4_isog_1x2x2x2w(P_1x2x2x2w, C24_A24plus_1x2x2x2w, coeff__0_1x2x2x2w, coeff2_1_1x2x2x2w), ITER_S);
  MEASURE_TIME(get_4_isog_1x2x2x2w(P_1x2x2x2w, C24_A24plus_1x2x2x2w, coeff__0_1x2x2x2w, coeff2_1_1x2x2x2w), ITER_L);
  printf("* 1x2x2x2w get_4_isog: %ld\n", diff_cycles);

  puts("");

  LOAD_CACHE(eval_4_isog_8x1x1x1w(&P_8x1x1x1w, coeff_8x1x1x1w), ITER_S);
  MEASURE_TIME(eval_4_isog_8x1x1x1w(&P_8x1x1x1w, coeff_8x1x1x1w), ITER_L);
  printf("* 8x1x1x1w eval_4_isog: %ld\n", diff_cycles);

  LOAD_CACHE(eval_4_isog_4x2x1x1w(P_4x2x1x1w, coeff_4x2x1x1w), ITER_S);
  MEASURE_TIME(eval_4_isog_4x2x1x1w(P_4x2x1x1w, coeff_4x2x1x1w), ITER_L);
  printf("* 4x2x1x1w eval_4_isog: %ld\n", diff_cycles);

  LOAD_CACHE(eval_4_isog_2x2x2x1w(P_2x2x2x1w, coeff_2x2x2x1w), ITER_S);
  MEASURE_TIME(eval_4_isog_2x2x2x1w(P_2x2x2x1w, coeff_2x2x2x1w), ITER_L);
  printf("* 2x2x2x1w eval_4_isog: %ld\n", diff_cycles);

  LOAD_CACHE(eval_4_isog_1x2x2x2w(P_1x2x2x2w, coeff_1x2x2x2w), ITER_S);
  MEASURE_TIME(eval_4_isog_1x2x2x2w(P_1x2x2x2w, coeff_1x2x2x2w), ITER_L);
  printf("* 1x2x2x2w eval_4_isog: %ld\n", diff_cycles);

  puts("");

  LOAD_CACHE(eval_3_isog_8x1x1x1w(&Q_8x1x1x1w, coeff_8x1x1x1w), ITER_S);
  MEASURE_TIME(eval_3_isog_8x1x1x1w(&Q_8x1x1x1w, coeff_8x1x1x1w), ITER_L);
  printf("* 8x1x1x1w eval_3_isog: %ld\n", diff_cycles);

  LOAD_CACHE(eval_3_isog_4x2x1x1w(Q_4x2x1x1w, coeff0_1_4x2x1x1w), ITER_S);
  MEASURE_TIME(eval_3_isog_4x2x1x1w(Q_4x2x1x1w, coeff0_1_4x2x1x1w), ITER_L);
  printf("* 4x2x1x1w eval_3_isog: %ld\n", diff_cycles);

  LOAD_CACHE(eval_3_isog_2x2x2x1w(Q_2x2x2x1w, coeff0_1_2x2x2x1w), ITER_S);
  MEASURE_TIME(eval_3_isog_2x2x2x1w(Q_2x2x2x1w, coeff0_1_2x2x2x1w), ITER_L);
  printf("* 2x2x2x1w eval_3_isog: %ld\n", diff_cycles);

  LOAD_CACHE(eval_3_isog_1x2x2x2w(Q_1x2x2x2w, coeff0_1_1x2x2x2w), ITER_S);
  MEASURE_TIME(eval_3_isog_1x2x2x2w(Q_1x2x2x2w, coeff0_1_1x2x2x2w), ITER_L);
  printf("* 1x2x2x2w eval_3_isog: %ld\n", diff_cycles);

  puts("");

  LOAD_CACHE(get_3_isog_8x1x1x1w(&P_8x1x1x1w, A24plus_8x1x1x1w, C24_8x1x1x1w, coeff_8x1x1x1w), ITER_S);
  MEASURE_TIME(get_3_isog_8x1x1x1w(&P_8x1x1x1w, A24plus_8x1x1x1w, C24_8x1x1x1w, coeff_8x1x1x1w), ITER_L);
  printf("* 8x1x1x1w get_3_isog: %ld\n", diff_cycles);

  LOAD_CACHE(get_3_isog_1x2x2x2w(P_1x2x2x2w, C24_A24plus_1x2x2x2w, coeff1_0_1x2x2x2w), ITER_S);
  MEASURE_TIME(get_3_isog_1x2x2x2w(P_1x2x2x2w, C24_A24plus_1x2x2x2w, coeff1_0_1x2x2x2w), ITER_L);
  printf("* 1x2x2x2w get_3_isog: %ld\n", diff_cycles);
}

void test_sidh()
{
  uint8_t skA[SECRETKEY_A_BYTES] = { 123 }, pkA[CRYPTO_PUBLICKEYBYTES];
  uint8_t skB[SECRETKEY_B_BYTES] = { 234 }, pkB[CRYPTO_PUBLICKEYBYTES];
  uint8_t ssA[FP2_ENCODED_BYTES], ssB[FP2_ENCODED_BYTES];

  int i;

  puts("\nCORRECTNESS TEST: ");

  randombytes(skA, SECRETKEY_A_BYTES);
  randombytes(skB, SECRETKEY_B_BYTES);

  EphemeralKeyGeneration_A(skA, pkA);
  EphemeralKeyGeneration_B(skB, pkB);

  puts("PKA: ");
  for (i = 0; i < CRYPTO_PUBLICKEYBYTES; i++) printf("%02X", pkA[i]);
  puts("\n");

  puts("PKB: ");
  for (i = 0; i < CRYPTO_PUBLICKEYBYTES; i++) printf("%02X", pkB[i]);
  puts("\n");

  EphemeralSecretAgreement_A(skA, pkB, ssA);
  EphemeralSecretAgreement_B(skB, pkA, ssB);

  puts("SSA: ");
  for (i = 0; i < FP2_ENCODED_BYTES; i++) printf("%02X", ssA[i]);
  puts("\n");

  puts("SSB: ");
  for (i = 0; i < FP2_ENCODED_BYTES; i++) printf("%02X", ssB[i]);
  puts("\n");

  // ---------------------------------------------------------------------------
  // measure timings 
  uint64_t start_cycles, end_cycles, diff_cycles;

  // vf2elm_t xP_8x1x1x1w, xQ_8x1x1x1w, xPQ_8x1x1x1w, A_8x1x1x1w;
  // __m512i m_8x1x1x1w[4];
  // vpoint_proj_t R_8x1x1x1w;
  // f2elm_r51_t xP, xQ, xPQ, A;
  // digit_t m[NWORDS_FIELD];
  // point_proj_r51_t R;

  puts("\nTIMINGS: ");

  // LOAD_CACHE(LADDER3PT_8x1x1x1w(xP_8x1x1x1w, xQ_8x1x1x1w, xPQ_8x1x1x1w, m_8x1x1x1w, ALICE, R_8x1x1x1w, A_8x1x1x1w), ITER_S);
  // MEASURE_TIME(LADDER3PT_8x1x1x1w(xP_8x1x1x1w, xQ_8x1x1x1w, xPQ_8x1x1x1w, m_8x1x1x1w, ALICE, R_8x1x1x1w, A_8x1x1x1w), ITER_M);
  // printf("* 8x1x1x1w LADDER3PT: %ld\n", diff_cycles);

  // LOAD_CACHE(LADDER3PT_1x4x2x1w(xP, xQ, xPQ, m, ALICE, R, A), ITER_S);
  // MEASURE_TIME(LADDER3PT_1x4x2x1w(xP, xQ, xPQ, m, ALICE, R, A), ITER_M);
  // printf("* 1x4x2x1w LADDER3PT: %ld\n", diff_cycles);

  LOAD_CACHE(EphemeralKeyGeneration_A(skA, pkA), ITER_S);
  MEASURE_TIME(EphemeralKeyGeneration_A(skA, pkA), ITER_M);
  printf("* KeyGenA: %ld\n", diff_cycles);

  LOAD_CACHE(EphemeralKeyGeneration_B(skB, pkB), ITER_S);
  MEASURE_TIME(EphemeralKeyGeneration_B(skB, pkB), ITER_M);
  printf("* KeyGenB: %ld\n", diff_cycles);

  LOAD_CACHE(EphemeralSecretAgreement_A(skA, pkB, ssA), ITER_S);
  MEASURE_TIME(EphemeralSecretAgreement_A(skA, pkB, ssA), ITER_M);
  printf("* SecAgrA: %ld\n", diff_cycles);

  LOAD_CACHE(EphemeralSecretAgreement_B(skB, pkA, ssB), ITER_S);
  MEASURE_TIME(EphemeralSecretAgreement_B(skB, pkA, ssB), ITER_M);
  printf("* SecAgrB: %ld\n", diff_cycles);
}

void test_sike()
{
  int i, wrong = 0;
  unsigned char sk[CRYPTO_SECRETKEYBYTES] = { 0 };
  unsigned char pk[CRYPTO_PUBLICKEYBYTES] = { 0 };
  unsigned char ct[CRYPTO_CIPHERTEXTBYTES] = { 0 };
  unsigned char ss[CRYPTO_BYTES] = { 0 };
  unsigned char ss_[CRYPTO_BYTES] = { 0 };

  crypto_kem_keypair(pk, sk);
  // crypto_kem_enc(ct, ss, pk);
  crypto_kem_enc_opt(ct, ss, pk);
  crypto_kem_dec(ss_, ct, sk);

  printf("SSA: ");
  for (i = 0; i < CRYPTO_BYTES; i++) printf("%02X", ss[i]);
  puts("");

  printf("SSB: ");
  for (i = 0; i < CRYPTO_BYTES; i++) printf("%02X", ss_[i]);
  puts("");

  wrong = memcmp(ss, ss_, CRYPTO_BYTES);

  if (wrong == 0) printf("\x1b[32m EQUAL!\x1b[0m\n");
  else            printf("\x1b[31m NOT EQUAL!\x1b[0m\n");

  // ---------------------------------------------------------------------------
  // measure timings 
  uint64_t start_cycles, end_cycles, diff_cycles;

  puts("\nTIMINGS: ");

  LOAD_CACHE(crypto_kem_keypair(pk, sk), ITER_S);
  MEASURE_TIME(crypto_kem_keypair(pk, sk), ITER_M);
  printf("* KeyGen: %ld\n", diff_cycles);

  // LOAD_CACHE(crypto_kem_enc(ct, ss, pk), ITER_S);
  // MEASURE_TIME(crypto_kem_enc(ct, ss, pk), ITER_M);
  // printf("* Encaps (NOT optimized)  : %ld\n", diff_cycles);

  LOAD_CACHE(crypto_kem_enc_opt(ct, ss, pk), ITER_S);
  MEASURE_TIME(crypto_kem_enc_opt(ct, ss, pk), ITER_M);
  printf("* Encaps: %ld\n", diff_cycles);

  LOAD_CACHE(crypto_kem_dec(ss_, ct, sk), ITER_S);
  MEASURE_TIME(crypto_kem_dec(ss_, ct, sk), ITER_M);
  printf("* Decaps: %ld\n", diff_cycles);
}

int main()
{
  // test_fp();
  // test_fpx();
  // test_curve();
  // test_isog();
  // test_sidh();
  test_sike();

  return 0;
}
