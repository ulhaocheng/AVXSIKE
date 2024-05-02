/**
 *******************************************************************************
 * @version 0.1
 * @date 2021-12-16
 * @copyright Copyright © 2021 by University of Luxembourg
 * @author Hao Cheng
 *******************************************************************************
 */

#include "sidh.h"
#include "sike.h"
#include "utils.h"
#include <time.h>
#include <string.h>


// the function to measure CPU cycles 
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


// simple test for sidh (just test one lane)
void test_sidh()
{
  __m512i skA[SK_A_VECTS] = { 0 }, pkA[6*NWORDS]; 
  __m512i skB[SK_B_VECTS] = { 0 }, pkB[6*NWORDS];
  __m512i ssA[2*NWORDS], ssB[2*NWORDS];
  uint64_t pkA51[NWORDS], pkA64[NWORDS] = { 0 };
  uint64_t pkB51[NWORDS], pkB64[NWORDS] = { 0 };
  uint64_t ssA51[2*NWORDS], ssA64[2*NWORDS] = { 0 };
  uint64_t ssB51[2*NWORDS], ssB64[2*NWORDS] = { 0 };
  int i, seed;

  seed = (int)time(NULL);
	srandom(seed);

  // // initialize skA 
  skA[0] = VSET1(random()); skA[1] = VSET1(random());
  skA[2] = VSET1(random()); skA[3] = VSET1(random());
  skA[4] = VSET1(random()); skA[5] = VSET1(random());
  // // initialize skB
  skB[0] = VSET1(random()); skB[1] = VSET1(random());
  skB[2] = VSET1(random()); skB[3] = VSET1(random());
  skB[4] = VSET1(random()); skB[5] = VSET1(random()); 

  EphemeralKeyGeneration_A(skA, pkA);
  EphemeralKeyGeneration_B(skB, pkB);
  
  EphemeralSecretAgreement_A(skA, pkB, ssA);
  EphemeralSecretAgreement_B(skB, pkA, ssB);

  get_channel(ssA51, ssA, 0);
  get_channel(&ssA51[NWORDS], &ssA[NWORDS], 0);
  get_channel(ssB51, ssB, 0);
  get_channel(&ssB51[NWORDS], &ssB[NWORDS], 0);

  mpi_conv_51to64(ssA64, ssA51, 2*NWORDS, 2*NWORDS);
  mpi_print(" ssA: ", ssA64, 2*NWORDS);
  mpi_conv_51to64(ssB64, ssB51, 2*NWORDS, 2*NWORDS);
  mpi_print(" ssB: ", ssB64, 2*NWORDS);
}

void timing_sidh()
{
  __m512i skA[SK_A_VECTS], pkA[6*NWORDS]; 
  __m512i skB[SK_B_VECTS], pkB[6*NWORDS];
  __m512i ssA[2*NWORDS], ssB[2*NWORDS];

  uint64_t start_cycles, end_cycles, diff_cycles;
  int i;

  LOAD_CACHE(EphemeralKeyGeneration_A(skA, pkA), ITER_S);
  MEASURE_TIME(EphemeralKeyGeneration_A(skA, pkA), ITER_S);
  printf("* KEYGEN A: %ld\n", diff_cycles);

  LOAD_CACHE(EphemeralKeyGeneration_B(skB, pkB), ITER_S);
  MEASURE_TIME(EphemeralKeyGeneration_B(skB, pkB), ITER_S);
  printf("* KEYGEN B: %ld\n", diff_cycles);

  LOAD_CACHE(EphemeralSecretAgreement_A(skA, pkB, ssA), ITER_S);
  MEASURE_TIME(EphemeralSecretAgreement_A(skA, pkB, ssA), ITER_S);
  printf("* SHRSEC A: %ld\n", diff_cycles);

  LOAD_CACHE(EphemeralSecretAgreement_B(skB, pkA, ssB), ITER_S);
  MEASURE_TIME(EphemeralSecretAgreement_B(skB, pkA, ssB), ITER_S);
  printf("* SHRSEC B: %ld\n", diff_cycles);
}

// 8-lane test for sike
void test_sike()
{
  // 8 instances
  uint8_t sk[INSTANCES][CRYPTO_SECRETKEYBYTES] = { 0 };
  uint8_t pk[INSTANCES][CRYPTO_PUBLICKEYBYTES] = { 0 };
  uint8_t ct[INSTANCES][CRYPTO_CIPHERTEXTBYTES] = { 0 };
  uint8_t ssa[INSTANCES][CRYPTO_BYTES] = { 0 };
  uint8_t ssb[INSTANCES][CRYPTO_BYTES] = { 0 };
  int i, k, wrong = 0;

  crypto_kem_keypair((uint8_t *)pk, (uint8_t *)sk);
  crypto_kem_enc((uint8_t *)ct, (uint8_t *)ssa, (uint8_t *)pk);
  crypto_kem_dec((uint8_t *)ssb, (uint8_t *)ct, (uint8_t *)sk);

  puts("\n*******************************************************************");
  puts("CORRECTNESS TEST:");
  puts("-------------------------------------------------------------------");
  puts("EIGHT instances in each execution of the program.\n");

  for (k = 0; k < INSTANCES; k++) {
    printf("* %d-instance\n", k);
    printf("- SSA: ");
    for (i = 0; i < CRYPTO_BYTES; i++) printf("%02X", ssa[k][i]);
    puts("");
    printf("- SSB: ");
    for (i = 0; i < CRYPTO_BYTES; i++) printf("%02X", ssb[k][i]);
    puts("\n");
    wrong |= memcmp(ssa[k], ssb[k], CRYPTO_BYTES);
  }

  if (wrong == 0) printf("\x1b[32mALL EQUAL!\x1b[0m\n");
  else            printf("\x1b[31mNOT EQUAL!\x1b[0m\n");
}

void timing_sike()
{
  // 8 instances
  uint8_t sk[INSTANCES][CRYPTO_SECRETKEYBYTES] = { 0 };
  uint8_t pk[INSTANCES][CRYPTO_PUBLICKEYBYTES] = { 0 };
  uint8_t ct[INSTANCES][CRYPTO_CIPHERTEXTBYTES] = { 0 };
  uint8_t ssa[INSTANCES][CRYPTO_BYTES] = { 0 };
  uint8_t ssb[INSTANCES][CRYPTO_BYTES] = { 0 };

  uint64_t start_cycles, end_cycles, diff_cycles;
  int i;

  puts("\n*******************************************************************");
  puts("Execution Time:\n");

  LOAD_CACHE(crypto_kem_keypair((uint8_t *)pk, (uint8_t *)sk), ITER_S);
  MEASURE_TIME(crypto_kem_keypair((uint8_t *)pk, (uint8_t *)sk), ITER_M);
  printf("* KEYGEN : %ld\n", diff_cycles);

  LOAD_CACHE(crypto_kem_enc((uint8_t *)ct, (uint8_t *)ssa, (uint8_t *)pk), ITER_S);
  MEASURE_TIME(crypto_kem_enc((uint8_t *)ct, (uint8_t *)ssa, (uint8_t *)pk), ITER_M);
  printf("* ENCAPS : %ld\n", diff_cycles);

  LOAD_CACHE(crypto_kem_dec((uint8_t *)ssb, (uint8_t *)ct, (uint8_t *)sk), ITER_S);
  MEASURE_TIME(crypto_kem_dec((uint8_t *)ssb, (uint8_t *)ct, (uint8_t *)sk), ITER_M);
  printf("* DECAPS : %ld\n", diff_cycles);
}

int main()
{
  // test_sidh();
  // timing_sidh();
  test_sike();
  timing_sike();
  return 0;
}
