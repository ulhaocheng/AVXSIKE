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

// radix-2^51 for the field elements: 9 * 51-bit = 459-bit
#define NWORDS 9
#define BRADIX 51
#define BMASK  0x7FFFFFFFFFFFFULL

static const uint64_t p434[NWORDS] = {
  0x7FFFFFFFFFFFF, 0x7FFFFFFFFFFFF, 0x7FFFFFFFFFFFF, 
  0x7FFFFFFFFFFFF, 0x7DC1767AE2FFF, 0x4B8F062B15D47, 
  0x5A07148159EF1, 0x0BB9A2367E2FE, 0x0000002341F27, };

static const uint64_t p434p1[NWORDS] = {
  0x0000000000000, 0x0000000000000, 0x0000000000000, 
  0x0000000000000, 0x7DC1767AE3000, 0x4B8F062B15D47, 
  0x5A07148159EF1, 0x0BB9A2367E2FE, 0x0000002341F27, };

static const uint64_t p434x2[NWORDS] = {
  0x7FFFFFFFFFFFE, 0x7FFFFFFFFFFFF, 0x7FFFFFFFFFFFF, 
  0x7FFFFFFFFFFFF, 0x7B82ECF5C5FFF, 0x171E0C562BA8F, 
  0x340E2902B3DE3, 0x1773446CFC5FD, 0x0000004683E4E, };

static const uint64_t p434x4[NWORDS] = {
  0x7FFFFFFFFFFFC, 0x7FFFFFFFFFFFF, 0x7FFFFFFFFFFFF, 
  0x7FFFFFFFFFFFF, 0x7705D9EB8BFFF, 0x2E3C18AC5751F, 
  0x681C520567BC6, 0x2EE688D9F8BFA, 0x0000008D07C9C, };

static const uint64_t mont_R[NWORDS] = {
  0x0000003A1635C, 0x0000000000000, 0x0000000000000, 
  0x0000000000000, 0x09B6230D6C000, 0x26DCE6EF7EDC0, 
  0x7F1F3F10E0FB2, 0x57CEB5B968572, 0x00000005C9696, };

static const uint64_t mont_R2[NWORDS] = {
  0x735A6CC0445BB, 0x4C611472ADB2E, 0x5ACEC73677687, 
  0x4F2E7F0622D11, 0x50B7242C31EAC, 0x18222E43835F8, 
  0x6FE26C5709746, 0x2B1006BB74802, 0x0000001DECB9D, }; 

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
