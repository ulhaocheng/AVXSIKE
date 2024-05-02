/**
 *******************************************************************************
 * @version 0.1
 * @date 2021-12-16
 * @copyright Copyright © 2021 by University of Luxembourg
 * @author Hao Cheng
 *******************************************************************************
 */

#ifndef _UTILS_H
#define _UTILS_H

#include <stdio.h>
#include <stdint.h>
#include "fp.h"

void mpi_print(const char *c, const uint64_t *a, int len);
void mpi_conv_64to51(uint64_t *r, const uint64_t *a, int rlen, int alen);
void mpi_conv_51to64(uint64_t *r, const uint64_t *a, int rlen, int alen);

static __m512i set_vector(const uint64_t a7, const uint64_t a6, const uint64_t a5, const uint64_t a4, 
                          const uint64_t a3, const uint64_t a2, const uint64_t a1, const uint64_t a0)
{
  __m512i r;

  ((uint64_t *)&r)[0] = a0; ((uint64_t *)&r)[1] = a1;
  ((uint64_t *)&r)[2] = a2; ((uint64_t *)&r)[3] = a3;
  ((uint64_t *)&r)[4] = a4; ((uint64_t *)&r)[5] = a5;
  ((uint64_t *)&r)[6] = a6; ((uint64_t *)&r)[7] = a7;

  return r;
}

#define _SET set_vector

#endif
