#ifndef FIPS202_H
#define FIPS202_H

#include <stdint.h>

#define SHAKE256_RATE 136

void shake256_absorb(uint64_t *s, const unsigned char *input, unsigned int inputByteLen);
void shake256_squeezeblocks(unsigned char *output, unsigned long long nblocks, uint64_t *s);
void shake256(unsigned char *output, unsigned long long outlen, const unsigned char *input,  unsigned long long inlen);

#endif
