/**
 *******************************************************************************
 * @version 0.1
 * @date 2021-12-16
 * @copyright Copyright © 2021 by University of Luxembourg
 * @author Hao Cheng
 *******************************************************************************
 */

#include "sike.h"
#include <string.h>

// we let SIDH primitives input/output radix-51 and vectors
// but SIKE primitives input/output radix-64 and strings

static void append_pk_to_sk(uint8_t *sk, const uint8_t *pk)
{
  int k;

  for (k = 0; k < INSTANCES; k++) 
    memcpy(&sk[k*CRYPTO_SECRETKEYBYTES+MSG_BYTES+SECRETKEY_B_BYTES], &pk[k*CRYPTO_PUBLICKEYBYTES], CRYPTO_PUBLICKEYBYTES);
}

// convert radix-51 vector to eight radix-64 strings 
static void vec_to_str(uint8_t *str, const int slen, const __m512i *vec, const int vnum)
{
  uint64_t t51[INSTANCES][vnum][NWORDS];
  uint64_t t64[vnum][NWORDS];
  int i, j, k;

  for (j = 0; j < vnum; j++) {
    get_channel(t51[0][j], &vec[j*NWORDS], 0);
    get_channel(t51[1][j], &vec[j*NWORDS], 1);
    get_channel(t51[2][j], &vec[j*NWORDS], 2);
    get_channel(t51[3][j], &vec[j*NWORDS], 3);
    get_channel(t51[4][j], &vec[j*NWORDS], 4);
    get_channel(t51[5][j], &vec[j*NWORDS], 5);
    get_channel(t51[6][j], &vec[j*NWORDS], 6);
    get_channel(t51[7][j], &vec[j*NWORDS], 7);
  }

  for (k = 0; k < INSTANCES; k++) {
    for (j = 0; j < vnum; j++) {
      mpi_conv_51to64(t64[j], t51[k][j], NWORDS, NWORDS); 
      memcpy(&str[k*slen+j*GFP_BYTES], t64[j], GFP_BYTES);
    }
  }
}

static void str_to_vec(__m512i *vec, const int vnum, const uint8_t *str, const int slen)
{
  uint64_t t51[INSTANCES][vnum][NWORDS];
  uint64_t t64[vnum][NWORDS];
  int i, j, k;

  memset(t64, 0, vnum*NWORDS*8);

  for (k = 0; k < INSTANCES; k++) {
    for (j = 0; j < vnum; j++)  {
      memcpy(t64[j], &str[k*slen+j*GFP_BYTES], GFP_BYTES);
      mpi_conv_64to51(t51[k][j], t64[j], NWORDS, NWORDS);
    }
  }

  for (j = 0; j < vnum; j++) {
    for (i = 0; i < NWORDS; i++) 
      vec[j*NWORDS+i] = set_vector(t51[7][j][i], t51[6][j][i], t51[5][j][i], t51[4][j][i], \
                                   t51[3][j][i], t51[2][j][i], t51[1][j][i], t51[0][j][i]);
  }
}

// SIKE key generation 
// single sk = s || SK || PK (16+28+55*6)-bytes <-> (2+4+6*9)-vectors
// single pk = PK            (55*6)-bytes       <-> (6*9)-vectors
void crypto_kem_keypair(uint8_t *pk, uint8_t *sk)
{
  __m512i vsk[SK_B_VECTS], vpk[6*NWORDS];

  // generate random strings for 8 instances  
  randombytes(sk                        , MSG_BYTES);
  randombytes(sk+  CRYPTO_SECRETKEYBYTES, MSG_BYTES);
  randombytes(sk+2*CRYPTO_SECRETKEYBYTES, MSG_BYTES);
  randombytes(sk+3*CRYPTO_SECRETKEYBYTES, MSG_BYTES);
  randombytes(sk+4*CRYPTO_SECRETKEYBYTES, MSG_BYTES);
  randombytes(sk+5*CRYPTO_SECRETKEYBYTES, MSG_BYTES);
  randombytes(sk+6*CRYPTO_SECRETKEYBYTES, MSG_BYTES);
  randombytes(sk+7*CRYPTO_SECRETKEYBYTES, MSG_BYTES);

  // generate random private keys for 8 instances 
  random_mod_order_B(vsk, sk);

  // generate vectorized public key (radix-51)
  EphemeralKeyGeneration_B(vsk, vpk);

  // vectorized pk -> pk strings (radix-51 -> radix-64)
  vec_to_str(pk, CRYPTO_PUBLICKEYBYTES, vpk, 6);

  // append pk strings to sk strings (radix-64)
  append_pk_to_sk(sk, pk);
}

// SIKE encapsulation
void crypto_kem_enc(uint8_t *ct, uint8_t *ss, const uint8_t *pk)
{
  __m512i vsk[SK_A_VECTS], vpk[6*NWORDS], vct[6*NWORDS], vjinv[2*NWORDS];
  uint8_t sk[INSTANCES][8*SK_A_VECTS] = { 0 };
  uint64_t *sk64 = (uint64_t *)sk;
  uint8_t temp[INSTANCES][CRYPTO_CIPHERTEXTBYTES+MSG_BYTES] = { 0 };
  uint8_t jinv[INSTANCES][2*GFP_BYTES];
  uint8_t h[INSTANCES][MSG_BYTES];
  int i, k;

  for (k = 0; k < INSTANCES; k++) {
    randombytes(temp[k], MSG_BYTES);
    memcpy(&temp[k][MSG_BYTES], &pk[k*CRYPTO_PUBLICKEYBYTES], CRYPTO_PUBLICKEYBYTES);
  }  
  shake256_8x1w(sk[0], sk[1], sk[2], sk[3], sk[4], sk[5], sk[6], sk[7], SECRETKEY_A_BYTES,
                temp[0], temp[1], temp[2], temp[3], temp[4], temp[5], temp[6], temp[7], CRYPTO_PUBLICKEYBYTES+MSG_BYTES);
  for (k = 0; k < INSTANCES; k++) sk[k][SECRETKEY_A_BYTES-1] &= MASK_ALICE;

  // form vectorized private key
  for (i = 0; i < SK_A_VECTS; i++) 
    vsk[i] = set_vector(sk64[7*SK_A_VECTS+i], sk64[6*SK_A_VECTS+i], \
                        sk64[5*SK_A_VECTS+i], sk64[4*SK_A_VECTS+i], \
                        sk64[3*SK_A_VECTS+i], sk64[2*SK_A_VECTS+i], \
                        sk64[  SK_A_VECTS+i], sk64[i]);

  // pk strings -> vectorized pk (radix-64 -> radix-51)
  str_to_vec(vpk, 6, pk, CRYPTO_PUBLICKEYBYTES);

  // Encrypt
  EphemeralKeyGeneration_A(vsk, vct);             
  EphemeralSecretAgreement_A(vsk, vpk, vjinv);    

  // vectorized ct and jinv -> ct and jinv strings (radix-51 -> radix-64)
  vec_to_str(ct, CRYPTO_CIPHERTEXTBYTES, vct, 6);                          
  vec_to_str((uint8_t *)jinv, 2*GFP_BYTES, vjinv, 2);

  shake256_8x1w(h[0], h[1], h[2], h[3], h[4], h[5], h[6], h[7], MSG_BYTES, \
                jinv[0], jinv[1], jinv[2], jinv[3], jinv[4], jinv[5], jinv[6], jinv[7], 2*GFP_BYTES);
  for (k = 0; k < INSTANCES; k++) {
    for (i = 0; i < MSG_BYTES; i++) 
      ct[i+CRYPTO_PUBLICKEYBYTES+k*CRYPTO_CIPHERTEXTBYTES] = temp[k][i] ^ h[k][i];
    memcpy(&temp[k][MSG_BYTES], &ct[k*CRYPTO_CIPHERTEXTBYTES], CRYPTO_CIPHERTEXTBYTES);
  }
  shake256_8x1w(ss, &ss[CRYPTO_BYTES], &ss[2*CRYPTO_BYTES], &ss[3*CRYPTO_BYTES], \
                &ss[4*CRYPTO_BYTES], &ss[5*CRYPTO_BYTES], &ss[6*CRYPTO_BYTES], &ss[7*CRYPTO_BYTES], CRYPTO_BYTES,\
                temp[0], temp[1], temp[2], temp[3], temp[4], temp[5], temp[6], temp[7], CRYPTO_CIPHERTEXTBYTES+MSG_BYTES);
}

// SIKE decapsulation
void crypto_kem_dec(uint8_t *ss, const uint8_t *ct, const uint8_t *sk)
{
  __m512i vct[6*NWORDS], vjinv[2*NWORDS], vsk[SK_A_VECTS], vc0[6*NWORDS];
  uint8_t tsk[INSTANCES][8*SK_A_VECTS] = { 0 };
  uint64_t *sk64 = (uint64_t *)tsk;
  uint8_t jinv[INSTANCES][2*GFP_BYTES];
  uint8_t h[INSTANCES][MSG_BYTES];
  uint8_t temp[INSTANCES][CRYPTO_CIPHERTEXTBYTES+MSG_BYTES];
  uint8_t c0[INSTANCES][CRYPTO_PUBLICKEYBYTES];
  int i, k;
  int8_t selector;

  // ct strings -> vectorized ct (radix-64 -> radix-51)
  str_to_vec(vct, 6, ct, CRYPTO_CIPHERTEXTBYTES);

  // form vectorized tsk
  for (k = 0; k < INSTANCES; k++) memcpy(tsk[k], sk+k*CRYPTO_SECRETKEYBYTES+MSG_BYTES, SECRETKEY_B_BYTES);
  for (i = 0; i < SK_B_VECTS; i++)
    vsk[i] = set_vector(sk64[7*SK_B_VECTS+i], sk64[6*SK_B_VECTS+i], \
                        sk64[5*SK_B_VECTS+i], sk64[4*SK_B_VECTS+i], \
                        sk64[3*SK_B_VECTS+i], sk64[2*SK_B_VECTS+i], \
                        sk64[  SK_B_VECTS+i], sk64[i]);

  // decrypt 
  EphemeralSecretAgreement_B(vsk, vct, vjinv);

  // vectorized jinv -> jinv strings (radix-51 -> radix-64)
  vec_to_str((uint8_t *)jinv, 2*GFP_BYTES, vjinv, 2);

  // clear tsk array and use it as esk
  memset(tsk, 0, INSTANCES*8*SK_A_VECTS);

  shake256_8x1w(h[0], h[1], h[2], h[3], h[4], h[5], h[6], h[7], MSG_BYTES, \
                jinv[0], jinv[1], jinv[2], jinv[3], jinv[4], jinv[5], jinv[6], jinv[7], 2*GFP_BYTES);
  for (k = 0; k < INSTANCES; k++) {
    for (i = 0; i < MSG_BYTES; i++) 
      temp[k][i] = ct[i+CRYPTO_PUBLICKEYBYTES+k*CRYPTO_CIPHERTEXTBYTES] ^ h[k][i];
    memcpy(&temp[k][MSG_BYTES], &sk[MSG_BYTES+SECRETKEY_B_BYTES+k*CRYPTO_SECRETKEYBYTES], CRYPTO_PUBLICKEYBYTES);
  }
  shake256_8x1w(tsk[0], tsk[1], tsk[2], tsk[3], tsk[4], tsk[5], tsk[6], tsk[7], SECRETKEY_A_BYTES, \
                temp[0], temp[1], temp[2], temp[3], temp[4], temp[5], temp[6], temp[7], CRYPTO_PUBLICKEYBYTES+MSG_BYTES);
  for (k = 0; k < INSTANCES; k++) tsk[k][SECRETKEY_A_BYTES-1] &= MASK_ALICE;

  // form vectorized esk
  for (i = 0; i < SK_A_VECTS; i++)
    vsk[i] = set_vector(sk64[7*SK_A_VECTS+i], sk64[6*SK_A_VECTS+i], \
                        sk64[5*SK_A_VECTS+i], sk64[4*SK_A_VECTS+i], \
                        sk64[3*SK_A_VECTS+i], sk64[2*SK_A_VECTS+i], \
                        sk64[  SK_A_VECTS+i], sk64[i]);

  // generate shared secret ss <- H(m||ct), or output ss <- H(s||ct) in case of ct verification failure
  EphemeralKeyGeneration_A(vsk, vc0);

  // vectorized c0 -> c0 strings (radix-51 -> radix-64)
  vec_to_str((uint8_t *)c0, CRYPTO_PUBLICKEYBYTES, vc0, 6);

  for (k = 0; k < INSTANCES; k++) {
    selector = ct_compare(c0[k], &ct[k*CRYPTO_CIPHERTEXTBYTES], CRYPTO_PUBLICKEYBYTES);
    ct_cmov(&temp[k][0], &sk[k*CRYPTO_SECRETKEYBYTES], MSG_BYTES, selector);
    memcpy(&temp[k][MSG_BYTES], &ct[k*CRYPTO_CIPHERTEXTBYTES], CRYPTO_CIPHERTEXTBYTES);
  }
  shake256_8x1w(ss, &ss[CRYPTO_BYTES], &ss[2*CRYPTO_BYTES], &ss[3*CRYPTO_BYTES], \
                &ss[4*CRYPTO_BYTES], &ss[5*CRYPTO_BYTES], &ss[6*CRYPTO_BYTES], &ss[7*CRYPTO_BYTES], CRYPTO_BYTES,\
                temp[0], temp[1], temp[2], temp[3], temp[4], temp[5], temp[6], temp[7], CRYPTO_CIPHERTEXTBYTES+MSG_BYTES);
}
