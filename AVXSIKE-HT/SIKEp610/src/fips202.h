#ifndef FIPS202_H
#define FIPS202_H

#include <stdint.h>

#define SHAKE256_RATE 136

void shake256_8x1w(uint8_t *out0, uint8_t *out1, uint8_t *out2, uint8_t *out3, 
                   uint8_t *out4, uint8_t *out5, uint8_t *out6, uint8_t *out7,
                   unsigned long long outlen,
                   const uint8_t *in0, const uint8_t *in1, const uint8_t *in2, const uint8_t *in3,
                   const uint8_t *in4, const uint8_t *in5, const uint8_t *in6, const uint8_t *in7,
                   unsigned long long inlen);

#endif
