/**
 *******************************************************************************
 * @version 0.1
 * @date 2021-12-16
 * @copyright Copyright © 2021 by University of Luxembourg
 * @author Hao Cheng
 *******************************************************************************
 */

#ifndef _GFPARIT_H
#define _GFPARIT_H

#include <stdint.h>
#include "intrin.h"

// radix-2^51 for the field elements: 10 * 51-bit = 510-bit
#define NWORDS 10
#define BRADIX 51
#define BMASK 0x7FFFFFFFFFFFFULL

static const uint64_t p503[NWORDS] = {
  0x7FFFFFFFFFFFF, 0x7FFFFFFFFFFFF, 0x7FFFFFFFFFFFF, 0x7FFFFFFFFFFFF,
  0x2BFFFFFFFFFFF, 0x0B7B44423CF41, 0x21EDF9F6BC4C2, 0x3BD2680DCDFB6,
  0x61E6045C6BDDA, 0x0080CDEA83023, };

static const uint64_t p503p1[NWORDS] = {
  0x0000000000000, 0x0000000000000, 0x0000000000000, 0x0000000000000, 
  0x2C00000000000, 0x0B7B44423CF41, 0x21EDF9F6BC4C2, 0x3BD2680DCDFB6, 
  0x61E6045C6BDDA, 0x0080CDEA83023, };

static const uint64_t p503x2[NWORDS] = {
  0x7FFFFFFFFFFFE, 0x7FFFFFFFFFFFF, 0x7FFFFFFFFFFFF, 0x7FFFFFFFFFFFF,
  0x57FFFFFFFFFFF, 0x16F6888479E82, 0x43DBF3ED78984, 0x77A4D01B9BF6C,
  0x43CC08B8D7BB4, 0x01019BD506047, };

static const uint64_t p503x4[NWORDS] = {
  0x7FFFFFFFFFFFC, 0x7FFFFFFFFFFFF, 0x7FFFFFFFFFFFF, 0x7FFFFFFFFFFFF,
  0x2FFFFFFFFFFFF, 0x2DED1108F3D05, 0x07B7E7DAF1308, 0x6F49A03737ED9,
  0x07981171AF769, 0x020337AA0C08F, };

static const uint64_t mont_R[NWORDS] = {
  0x00000000000FE, 0x0000000000000, 0x0000000000000, 0x0000000000000,  
  0x5800000000000, 0x1BB2464785D2A, 0x55E1FD312C76D, 0x253CC24DA0928, 
  0x5DC7AC4CFA13D, 0x0033B15203C83, };

static const uint64_t mont_R2[NWORDS] = {
  0x09A0CF641D011, 0x2E313FDA572A5, 0x3723C5EA6E209, 0x347651D9B2EAC,
  0x56FC57AB6EFF1, 0x1DF6A7AEB3AD1, 0x62BA292CF5EEF, 0x58BBAF42ED687,
  0x5B38EAFFB79B9, 0x0080AEAB77E91, }; 

void mp_sub_p2(__m512i *r, const __m512i *a, const __m512i *b);
void mp_sub_p4(__m512i *r, const __m512i *a, const __m512i *b);
void fpadd(__m512i *r, const __m512i *a, const __m512i* b);
void fpsub(__m512i *r, const __m512i *a, const __m512i* b);
void fpneg(__m512i *r);
void fpdiv2(__m512i *r, const __m512i *a);
void fpcorrection(__m512i *r);
void rdc_mont(__m512i *r, const __m512i *a);

void mp_mul_v1(__m512i *r, const __m512i *a, const __m512i *b);
void mp_mul_v2(__m512i *r, const __m512i *a, const __m512i *b);

#define mp_mul mp_mul_v2

#endif