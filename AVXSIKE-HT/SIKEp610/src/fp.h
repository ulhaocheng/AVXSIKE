/**
 *******************************************************************************
 * @version 0.1
 * @date 2021-12-16
 * @copyright Copyright © 2021 by University of Luxembourg
 * @author Hao Cheng
 *******************************************************************************
 */

#ifndef _GFPARITH_H
#define _GFPARITH_H

#include <stdint.h>
#include "intrin.h"

// radix-2^51 for the field elements: 12 * 51-bit = 612-bit
#define NWORDS 12
#define BRADIX 51
#define BMASK 0x7FFFFFFFFFFFFULL

static const uint64_t p610[NWORDS] = {
  0x7FFFFFFFFFFFF, 0x7FFFFFFFFFFFF, 0x7FFFFFFFFFFFF, 0x7FFFFFFFFFFFF,
  0x7FFFFFFFFFFFF, 0x3FFFFFFFFFFFF, 0x22A96AC0B9B80, 0x47FCD5D8BC26F,
  0x52A9AE7BF4504, 0x64AB65F421884, 0x04309479F6232, 0x13DFB53B440C8, };

static const uint64_t p610p1[NWORDS] = {
  0x0000000000000, 0x0000000000000, 0x0000000000000, 0x0000000000000, 
  0x0000000000000, 0x4000000000000, 0x22A96AC0B9B80, 0x47FCD5D8BC26F,
  0x52A9AE7BF4504, 0x64AB65F421884, 0x04309479F6232, 0x13DFB53B440C8, };

static const uint64_t p610x2[NWORDS] = {
  0x7FFFFFFFFFFFE, 0x7FFFFFFFFFFFF, 0x7FFFFFFFFFFFF, 0x7FFFFFFFFFFFF,
  0x7FFFFFFFFFFFF, 0x7FFFFFFFFFFFF, 0x4552D58173700, 0x0FF9ABB1784DE,
  0x25535CF7E8A09, 0x4956CBE843109, 0x086128F3EC465, 0x27BF6A7688190, };

static const uint64_t p610x4[NWORDS] = {
  0x7FFFFFFFFFFFC, 0x7FFFFFFFFFFFF, 0x7FFFFFFFFFFFF, 0x7FFFFFFFFFFFF,
  0x7FFFFFFFFFFFF, 0x7FFFFFFFFFFFF, 0x0AA5AB02E6E01, 0x1FF35762F09BD,
  0x4AA6B9EFD1412, 0x12AD97D086212, 0x10C251E7D88CB, 0x4F7ED4ED10320, };

static const uint64_t mont_R[NWORDS] = {
  0x0000000000006, 0x0000000000000, 0x0000000000000, 0x0000000000000, 
  0x0000000000000, 0x0000000000000, 0x30077F7BA5AFD, 0x5012FCEB97164,
  0x1005E918461E4, 0x23FB9C4736CE4, 0x66DC85243B2CF, 0x08C1C09C67B4F, };

static const uint64_t mont_R2[NWORDS] = {
  0x163B627392EE7, 0x31BC927BC170B, 0x68315AF05C1E0, 0x56E3FA0CCA068,
  0x26030979EDE54, 0x744E105718847, 0x42E56CDD585D4, 0x277C0450B087A,
  0x6CBEFFAD0601E, 0x7CBD17B195293, 0x2BF3C820BCD5B, 0x07673F03D2224, }; 

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
