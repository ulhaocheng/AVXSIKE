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

#define NROUNDS 24
#define ROL(a, offset) ((a << offset) ^ (a >> (64-offset)))

// Host to little endian, little endian to host
#define HTOLE_64(i) (i)
#define LETOH_64(i) (i)

static uint64_t load64(const unsigned char *x)
{
  return LETOH_64(*((uint64_t*)x));
}


static void store64(uint8_t *x, uint64_t u)
{
  *(uint64_t*)x = HTOLE_64(u);
}


static const uint64_t KeccakF_RoundConstants[NROUNDS] = 
{
    (uint64_t)0x0000000000000001ULL,
    (uint64_t)0x0000000000008082ULL,
    (uint64_t)0x800000000000808aULL,
    (uint64_t)0x8000000080008000ULL,
    (uint64_t)0x000000000000808bULL,
    (uint64_t)0x0000000080000001ULL,
    (uint64_t)0x8000000080008081ULL,
    (uint64_t)0x8000000000008009ULL,
    (uint64_t)0x000000000000008aULL,
    (uint64_t)0x0000000000000088ULL,
    (uint64_t)0x0000000080008009ULL,
    (uint64_t)0x000000008000000aULL,
    (uint64_t)0x000000008000808bULL,
    (uint64_t)0x800000000000008bULL,
    (uint64_t)0x8000000000008089ULL,
    (uint64_t)0x8000000000008003ULL,
    (uint64_t)0x8000000000008002ULL,
    (uint64_t)0x8000000000000080ULL,
    (uint64_t)0x000000000000800aULL,
    (uint64_t)0x800000008000000aULL,
    (uint64_t)0x8000000080008081ULL,
    (uint64_t)0x8000000000008080ULL,
    (uint64_t)0x0000000080000001ULL,
    (uint64_t)0x8000000080008008ULL
};


void KeccakF1600_StatePermute(uint64_t * state)
{
  int round;

        uint64_t Aba, Abe, Abi, Abo, Abu;
        uint64_t Aga, Age, Agi, Ago, Agu;
        uint64_t Aka, Ake, Aki, Ako, Aku;
        uint64_t Ama, Ame, Ami, Amo, Amu;
        uint64_t Asa, Ase, Asi, Aso, Asu;
        uint64_t BCa, BCe, BCi, BCo, BCu;
        uint64_t Da, De, Di, Do, Du;
        uint64_t Eba, Ebe, Ebi, Ebo, Ebu;
        uint64_t Ega, Ege, Egi, Ego, Egu;
        uint64_t Eka, Eke, Eki, Eko, Eku;
        uint64_t Ema, Eme, Emi, Emo, Emu;
        uint64_t Esa, Ese, Esi, Eso, Esu;

        //copyFromState(A, state)
        Aba = state[ 0];
        Abe = state[ 1];
        Abi = state[ 2];
        Abo = state[ 3];
        Abu = state[ 4];
        Aga = state[ 5];
        Age = state[ 6];
        Agi = state[ 7];
        Ago = state[ 8];
        Agu = state[ 9];
        Aka = state[10];
        Ake = state[11];
        Aki = state[12];
        Ako = state[13];
        Aku = state[14];
        Ama = state[15];
        Ame = state[16];
        Ami = state[17];
        Amo = state[18];
        Amu = state[19];
        Asa = state[20];
        Ase = state[21];
        Asi = state[22];
        Aso = state[23];
        Asu = state[24];

        for( round = 0; round < NROUNDS; round += 2 )
        {
            //    prepareTheta
            BCa = Aba^Aga^Aka^Ama^Asa;
            BCe = Abe^Age^Ake^Ame^Ase;
            BCi = Abi^Agi^Aki^Ami^Asi;
            BCo = Abo^Ago^Ako^Amo^Aso;
            BCu = Abu^Agu^Aku^Amu^Asu;

            //thetaRhoPiChiIotaPrepareTheta(round  , A, E)
            Da = BCu^ROL(BCe, 1);
            De = BCa^ROL(BCi, 1);
            Di = BCe^ROL(BCo, 1);
            Do = BCi^ROL(BCu, 1);
            Du = BCo^ROL(BCa, 1);

            Aba ^= Da;
            BCa = Aba;
            Age ^= De;
            BCe = ROL(Age, 44);
            Aki ^= Di;
            BCi = ROL(Aki, 43);
            Amo ^= Do;
            BCo = ROL(Amo, 21);
            Asu ^= Du;
            BCu = ROL(Asu, 14);
            Eba =   BCa ^((~BCe)&  BCi );
            Eba ^= (uint64_t)KeccakF_RoundConstants[round];
            Ebe =   BCe ^((~BCi)&  BCo );
            Ebi =   BCi ^((~BCo)&  BCu );
            Ebo =   BCo ^((~BCu)&  BCa );
            Ebu =   BCu ^((~BCa)&  BCe );

            Abo ^= Do;
            BCa = ROL(Abo, 28);
            Agu ^= Du;
            BCe = ROL(Agu, 20);
            Aka ^= Da;
            BCi = ROL(Aka,  3);
            Ame ^= De;
            BCo = ROL(Ame, 45);
            Asi ^= Di;
            BCu = ROL(Asi, 61);
            Ega =   BCa ^((~BCe)&  BCi );
            Ege =   BCe ^((~BCi)&  BCo );
            Egi =   BCi ^((~BCo)&  BCu );
            Ego =   BCo ^((~BCu)&  BCa );
            Egu =   BCu ^((~BCa)&  BCe );

            Abe ^= De;
            BCa = ROL(Abe,  1);
            Agi ^= Di;
            BCe = ROL(Agi,  6);
            Ako ^= Do;
            BCi = ROL(Ako, 25);
            Amu ^= Du;
            BCo = ROL(Amu,  8);
            Asa ^= Da;
            BCu = ROL(Asa, 18);
            Eka =   BCa ^((~BCe)&  BCi );
            Eke =   BCe ^((~BCi)&  BCo );
            Eki =   BCi ^((~BCo)&  BCu );
            Eko =   BCo ^((~BCu)&  BCa );
            Eku =   BCu ^((~BCa)&  BCe );

            Abu ^= Du;
            BCa = ROL(Abu, 27);
            Aga ^= Da;
            BCe = ROL(Aga, 36);
            Ake ^= De;
            BCi = ROL(Ake, 10);
            Ami ^= Di;
            BCo = ROL(Ami, 15);
            Aso ^= Do;
            BCu = ROL(Aso, 56);
            Ema =   BCa ^((~BCe)&  BCi );
            Eme =   BCe ^((~BCi)&  BCo );
            Emi =   BCi ^((~BCo)&  BCu );
            Emo =   BCo ^((~BCu)&  BCa );
            Emu =   BCu ^((~BCa)&  BCe );

            Abi ^= Di;
            BCa = ROL(Abi, 62);
            Ago ^= Do;
            BCe = ROL(Ago, 55);
            Aku ^= Du;
            BCi = ROL(Aku, 39);
            Ama ^= Da;
            BCo = ROL(Ama, 41);
            Ase ^= De;
            BCu = ROL(Ase,  2);
            Esa =   BCa ^((~BCe)&  BCi );
            Ese =   BCe ^((~BCi)&  BCo );
            Esi =   BCi ^((~BCo)&  BCu );
            Eso =   BCo ^((~BCu)&  BCa );
            Esu =   BCu ^((~BCa)&  BCe );

            //    prepareTheta
            BCa = Eba^Ega^Eka^Ema^Esa;
            BCe = Ebe^Ege^Eke^Eme^Ese;
            BCi = Ebi^Egi^Eki^Emi^Esi;
            BCo = Ebo^Ego^Eko^Emo^Eso;
            BCu = Ebu^Egu^Eku^Emu^Esu;

            //thetaRhoPiChiIotaPrepareTheta(round+1, E, A)
            Da = BCu^ROL(BCe, 1);
            De = BCa^ROL(BCi, 1);
            Di = BCe^ROL(BCo, 1);
            Do = BCi^ROL(BCu, 1);
            Du = BCo^ROL(BCa, 1);

            Eba ^= Da;
            BCa = Eba;
            Ege ^= De;
            BCe = ROL(Ege, 44);
            Eki ^= Di;
            BCi = ROL(Eki, 43);
            Emo ^= Do;
            BCo = ROL(Emo, 21);
            Esu ^= Du;
            BCu = ROL(Esu, 14);
            Aba =   BCa ^((~BCe)&  BCi );
            Aba ^= (uint64_t)KeccakF_RoundConstants[round+1];
            Abe =   BCe ^((~BCi)&  BCo );
            Abi =   BCi ^((~BCo)&  BCu );
            Abo =   BCo ^((~BCu)&  BCa );
            Abu =   BCu ^((~BCa)&  BCe );

            Ebo ^= Do;
            BCa = ROL(Ebo, 28);
            Egu ^= Du;
            BCe = ROL(Egu, 20);
            Eka ^= Da;
            BCi = ROL(Eka, 3);
            Eme ^= De;
            BCo = ROL(Eme, 45);
            Esi ^= Di;
            BCu = ROL(Esi, 61);
            Aga =   BCa ^((~BCe)&  BCi );
            Age =   BCe ^((~BCi)&  BCo );
            Agi =   BCi ^((~BCo)&  BCu );
            Ago =   BCo ^((~BCu)&  BCa );
            Agu =   BCu ^((~BCa)&  BCe );

            Ebe ^= De;
            BCa = ROL(Ebe, 1);
            Egi ^= Di;
            BCe = ROL(Egi, 6);
            Eko ^= Do;
            BCi = ROL(Eko, 25);
            Emu ^= Du;
            BCo = ROL(Emu, 8);
            Esa ^= Da;
            BCu = ROL(Esa, 18);
            Aka =   BCa ^((~BCe)&  BCi );
            Ake =   BCe ^((~BCi)&  BCo );
            Aki =   BCi ^((~BCo)&  BCu );
            Ako =   BCo ^((~BCu)&  BCa );
            Aku =   BCu ^((~BCa)&  BCe );

            Ebu ^= Du;
            BCa = ROL(Ebu, 27);
            Ega ^= Da;
            BCe = ROL(Ega, 36);
            Eke ^= De;
            BCi = ROL(Eke, 10);
            Emi ^= Di;
            BCo = ROL(Emi, 15);
            Eso ^= Do;
            BCu = ROL(Eso, 56);
            Ama =   BCa ^((~BCe)&  BCi );
            Ame =   BCe ^((~BCi)&  BCo );
            Ami =   BCi ^((~BCo)&  BCu );
            Amo =   BCo ^((~BCu)&  BCa );
            Amu =   BCu ^((~BCa)&  BCe );

            Ebi ^= Di;
            BCa = ROL(Ebi, 62);
            Ego ^= Do;
            BCe = ROL(Ego, 55);
            Eku ^= Du;
            BCi = ROL(Eku, 39);
            Ema ^= Da;
            BCo = ROL(Ema, 41);
            Ese ^= De;
            BCu = ROL(Ese, 2);
            Asa =   BCa ^((~BCe)&  BCi );
            Ase =   BCe ^((~BCi)&  BCo );
            Asi =   BCi ^((~BCo)&  BCu );
            Aso =   BCo ^((~BCu)&  BCa );
            Asu =   BCu ^((~BCa)&  BCe );
        }

        //copyToState(state, A)
        state[ 0] = Aba;
        state[ 1] = Abe;
        state[ 2] = Abi;
        state[ 3] = Abo;
        state[ 4] = Abu;
        state[ 5] = Aga;
        state[ 6] = Age;
        state[ 7] = Agi;
        state[ 8] = Ago;
        state[ 9] = Agu;
        state[10] = Aka;
        state[11] = Ake;
        state[12] = Aki;
        state[13] = Ako;
        state[14] = Aku;
        state[15] = Ama;
        state[16] = Ame;
        state[17] = Ami;
        state[18] = Amo;
        state[19] = Amu;
        state[20] = Asa;
        state[21] = Ase;
        state[22] = Asi;
        state[23] = Aso;
        state[24] = Asu;

        #undef    round
}



#include <string.h>
#define MIN(a, b) ((a) < (b) ? (a) : (b))


static void keccak_absorb(uint64_t *s, unsigned int r, const unsigned char *m, unsigned long long int mlen, unsigned char p)
{
  unsigned long long i;
  unsigned char t[200];
 
  while (mlen >= r) 
  {
    for (i = 0; i < r / 8; ++i)
      s[i] ^= load64(m + 8 * i);
    
    KeccakF1600_StatePermute(s);
    mlen -= r;
    m += r;
  }

  for (i = 0; i < r; ++i)
    t[i] = 0;
  for (i = 0; i < mlen; ++i)
    t[i] = m[i];
  t[i] = p;
  t[r - 1] |= 128;
  for (i = 0; i < r / 8; ++i)
    s[i] ^= load64(t + 8 * i);
}


static void keccak_squeezeblocks(unsigned char *h, unsigned long long int nblocks, uint64_t *s, unsigned int r)
{
  unsigned int i;

  while(nblocks > 0) 
  {
    KeccakF1600_StatePermute(s);
    for (i = 0; i < (r>>3); i++)
    {
      store64(h+8*i, s[i]);
    }
    h += r;
    nblocks--;
  }
}

/********** SHAKE256 ***********/

void shake256_absorb(uint64_t *s, const unsigned char *input, unsigned int inputByteLen)
{
	keccak_absorb(s, SHAKE256_RATE, input, inputByteLen, 0x1F);
}


void shake256_squeezeblocks(unsigned char *output, unsigned long long nblocks, uint64_t *s)
{
	keccak_squeezeblocks(output, nblocks, s, SHAKE256_RATE);
}


void shake256(unsigned char *output, unsigned long long outlen, const unsigned char *input,  unsigned long long inlen)
{
  uint64_t s[25];
  unsigned char t[SHAKE256_RATE];
  unsigned long long nblocks = outlen/SHAKE256_RATE;
  size_t i;

  for (i = 0; i < 25; ++i)
    s[i] = 0;
  
  /* Absorb input */
  keccak_absorb(s, SHAKE256_RATE, input, inlen, 0x1F);

  /* Squeeze output */
  keccak_squeezeblocks(output, nblocks, s, SHAKE256_RATE);

  output += nblocks*SHAKE256_RATE;
  outlen -= nblocks*SHAKE256_RATE;

  if (outlen) 
  {
    keccak_squeezeblocks(t, 1, s, SHAKE256_RATE);
    for (i = 0; i < outlen; i++)
      output[i] = t[i];
  }
}


// AVX-512 implementation 
#include "intrin.h"

extern void  KeccakP1600times8_PermuteAll_24rounds(void *states);
#define KeccakF1600_StatePermutex8 KeccakP1600times8_PermuteAll_24rounds

static void keccak_absorbx8(__m512i *s, unsigned int r, const uint8_t *m0, const uint8_t *m1, 
                            const uint8_t *m2, const uint8_t *m3, const uint8_t *m4, const uint8_t *m5,
                            const uint8_t *m6, const uint8_t *m7, unsigned long long int mlen, unsigned char p) 
{
  unsigned long long i;
  unsigned char t0[200], t1[200], t2[200], t3[200], t4[200], t5[200], t6[200], t7[200];
  unsigned long long *ss = (unsigned long long *)s;

  while (mlen >= r)
  {
    for (i = 0; i < r/8; i++) {
      ss[8*i  ] ^= load64(m0+8*i);  ss[8*i+1] ^= load64(m1+8*i);
      ss[8*i+2] ^= load64(m2+8*i);  ss[8*i+3] ^= load64(m3+8*i);
      ss[8*i+4] ^= load64(m4+8*i);  ss[8*i+5] ^= load64(m5+8*i);
      ss[8*i+6] ^= load64(m6+8*i);  ss[8*i+7] ^= load64(m7+8*i);
    }

    KeccakF1600_StatePermutex8(s);
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
    ss[8*i  ] ^= load64(t0+8*i);  ss[8*i+1] ^= load64(t1+8*i);
    ss[8*i+2] ^= load64(t2+8*i);  ss[8*i+3] ^= load64(t3+8*i);
    ss[8*i+4] ^= load64(t4+8*i);  ss[8*i+5] ^= load64(t5+8*i);
    ss[8*i+6] ^= load64(t6+8*i);  ss[8*i+7] ^= load64(t7+8*i);
  }
}

static void keccak_squeezeblocksx8(uint8_t *h0, uint8_t *h1, uint8_t *h2, uint8_t *h3, 
                                   uint8_t *h4, uint8_t *h5, uint8_t *h6, uint8_t *h7, 
                                   unsigned long long int nblocks, __m512i *s, unsigned int r)
{
  unsigned int i;
  unsigned long long *ss = (unsigned long long *)s;

  while (nblocks > 0) {
    KeccakF1600_StatePermutex8(s);
    for (i = 0; i < (r>>3); i++) {
      store64(h0+8*i, ss[8*i  ]);  store64(h1+8*i, ss[8*i+1]);
      store64(h2+8*i, ss[8*i+2]);  store64(h3+8*i, ss[8*i+3]);
      store64(h4+8*i, ss[8*i+4]);  store64(h5+8*i, ss[8*i+5]);
      store64(h6+8*i, ss[8*i+6]);  store64(h7+8*i, ss[8*i+7]);
    }
    h0 += r; h1 += r; h2 += r; h3 += r;
    h4 += r; h5 += r; h6 += r; h7 += r;
    nblocks--;
  }
}
