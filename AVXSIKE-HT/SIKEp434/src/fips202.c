/********************************************************************************************
* SHA3-derived function SHAKE
*
* Based on the public domain implementation in crypto_hash/keccakc512/simple/ 
* from http://bench.cr.yp.to/supercop.html by Ronny Van Keer 
* and the public domain "TweetFips202" implementation from https://twitter.com/tweetfips202 
* by Gilles Van Assche, Daniel J. Bernstein, and Peter Schwabe
*
* See NIST Special Publication 800-185 for more information:
* http://nvlpubs.nist.gov/nistpubs/SpecialPublications/NIST.SP.800-185.pdf
*
*********************************************************************************************/  

#include <stdint.h>
#include <assert.h>
#include "fips202.h"
#include "utils.h"

// Host to little endian, little endian to host
#define HTOLE_64(i) (i)
#define LETOH_64(i) (i)

static uint64_t load64(const uint8_t *x)
{
  return LETOH_64(*((uint64_t*)x));
}

static void store64(uint8_t *x, uint64_t u)
{
  *(uint64_t*)x = HTOLE_64(u);
}

// AVX-512 implementation 
#include "intrin.h"

extern void  KeccakP1600times8_PermuteAll_24rounds(void *states);

static void keccak_absorb_8x1w(__m512i *s, unsigned int r, const uint8_t *m0, const uint8_t *m1, 
                               const uint8_t *m2, const uint8_t *m3, const uint8_t *m4, const uint8_t *m5,
                               const uint8_t *m6, const uint8_t *m7, unsigned long long int mlen, unsigned char p) 
{
  unsigned long long i;
  unsigned char t0[200], t1[200], t2[200], t3[200], t4[200], t5[200], t6[200], t7[200];
  __m512i a;

  while (mlen >= r)
  {
    for (i = 0; i < r/8; i++) {
      a = set_vector(load64(m7+8*i), load64(m6+8*i), load64(m5+8*i), load64(m4+8*i), \
                     load64(m3+8*i), load64(m2+8*i), load64(m1+8*i), load64(m0+8*i));
      s[i] = VXOR(s[i], a);
    }

    KeccakP1600times8_PermuteAll_24rounds(s);
    mlen -= r;
    m0 += r; m1 += r; m2 += r; m3 += r;
    m4 += r; m5 += r; m6 += r; m7 += r;
  }

  for (i = 0; i < r; i++) {
    t0[i] = t1[i] = t2[i] = t3[i] = 0;
    t4[i] = t5[i] = t6[i] = t7[i] = 0;
  }

  for (i = 0; i < mlen; i++) {
    t0[i] = m0[i]; t1[i] = m1[i]; t2[i] = m2[i]; t3[i] = m3[i];
    t4[i] = m4[i]; t5[i] = m5[i]; t6[i] = m6[i]; t7[i] = m7[i];
  }

  t0[i] = t1[i] = t2[i] = t3[i] = p;
  t4[i] = t5[i] = t6[i] = t7[i] = p;

  t0[r-1] |= 128; t1[r-1] |= 128; t2[r-1] |= 128; t3[r-1] |= 128; 
  t4[r-1] |= 128; t5[r-1] |= 128; t6[r-1] |= 128; t7[r-1] |= 128;  

  for (i = 0; i < r/8; i++) {
    a = set_vector(load64(t7+8*i), load64(t6+8*i), load64(t5+8*i), load64(t4+8*i), \
                   load64(t3+8*i), load64(t2+8*i), load64(t1+8*i), load64(t0+8*i));
    s[i] = VXOR(s[i], a);
  }
}

static void keccak_squeezeblocks_8x1w(uint8_t *h0, uint8_t *h1, uint8_t *h2, uint8_t *h3, 
                                      uint8_t *h4, uint8_t *h5, uint8_t *h6, uint8_t *h7, 
                                      unsigned long long int nblocks, __m512i *s, unsigned int r)
{
  unsigned int i;

  while (nblocks > 0) {
    KeccakP1600times8_PermuteAll_24rounds(s);
    for (i = 0; i < (r>>3); i++) {
      store64(h0+8*i, ((uint64_t *)&s[i])[0]);
      store64(h1+8*i, ((uint64_t *)&s[i])[1]);
      store64(h2+8*i, ((uint64_t *)&s[i])[2]);
      store64(h3+8*i, ((uint64_t *)&s[i])[3]);
      store64(h4+8*i, ((uint64_t *)&s[i])[4]);
      store64(h5+8*i, ((uint64_t *)&s[i])[5]);
      store64(h6+8*i, ((uint64_t *)&s[i])[6]);
      store64(h7+8*i, ((uint64_t *)&s[i])[7]);
    }
    h0 += r; h1 += r; h2 += r; h3 += r;
    h4 += r; h5 += r; h6 += r; h7 += r;
    nblocks--;
  }
}

void shake256_8x1w(uint8_t *out0, uint8_t *out1, uint8_t *out2, uint8_t *out3, 
                   uint8_t *out4, uint8_t *out5, uint8_t *out6, uint8_t *out7,
                   unsigned long long outlen,
                   const uint8_t *in0, const uint8_t *in1, const uint8_t *in2, const uint8_t *in3,
                   const uint8_t *in4, const uint8_t *in5, const uint8_t *in6, const uint8_t *in7,
                   unsigned long long inlen)
{
  __m512i s[25];
  uint8_t t0[SHAKE256_RATE], t1[SHAKE256_RATE], t2[SHAKE256_RATE], t3[SHAKE256_RATE];
  uint8_t t4[SHAKE256_RATE], t5[SHAKE256_RATE], t6[SHAKE256_RATE], t7[SHAKE256_RATE];
  unsigned long long nblocks = outlen/SHAKE256_RATE;
  unsigned int i;

  for (i = 0; i < 25; i++) s[i] = VZERO;

  /* Absorb input */
  keccak_absorb_8x1w(s, SHAKE256_RATE, in0, in1, in2, in3, in4, in5, in6, in7, inlen, 0x1F);

  /* Squeeze output */
  keccak_squeezeblocks_8x1w(out0, out1, out2, out3, out4, out5, out6, out7, nblocks, s, SHAKE256_RATE);  

  out0 += nblocks*SHAKE256_RATE;  out1 += nblocks*SHAKE256_RATE;
  out2 += nblocks*SHAKE256_RATE;  out3 += nblocks*SHAKE256_RATE;
  out4 += nblocks*SHAKE256_RATE;  out5 += nblocks*SHAKE256_RATE;
  out6 += nblocks*SHAKE256_RATE;  out7 += nblocks*SHAKE256_RATE;
  outlen -= nblocks*SHAKE256_RATE;

  if (outlen) {
    keccak_squeezeblocks_8x1w(t0, t1, t2, t3, t4, t5, t6, t7, 1, s, SHAKE256_RATE);

    for (i = 0; i < outlen; i++) {
      out0[i] = t0[i]; out1[i] = t1[i];
      out2[i] = t2[i]; out3[i] = t3[i];
      out4[i] = t4[i]; out5[i] = t5[i];
      out6[i] = t6[i]; out7[i] = t7[i];
    }
  }
}
