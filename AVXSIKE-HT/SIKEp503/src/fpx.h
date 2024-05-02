/**
 *******************************************************************************
 * @version 0.1
 * @date 2021-12-16
 * @copyright Copyright © 2021 by University of Luxembourg
 * @author Hao Cheng
 *******************************************************************************
 */

#ifndef _FPX_H
#define _FPX_H

#include "fp.h"

typedef __m512i  felm_t[NWORDS];
typedef __m512i  dfelm_t[2*NWORDS];
typedef felm_t   f2elm_t[2]; 

int8_t ct_compare(const uint8_t *a, const uint8_t *b, unsigned int len); 
void ct_cmov(uint8_t *r, const uint8_t *a, unsigned int len, int8_t selector); 

void fpcopy(felm_t r, const felm_t a);
void fpzero(felm_t r);

void to_mont(felm_t r, const felm_t a);
void from_mont(felm_t r, const felm_t a);

#define fpmul_mont fpmul_mont_v1
#define fpsqr_mont fpsqr_mont_v2

void fpmul_mont_v1(felm_t r, const felm_t a, const felm_t b);
void fpmul_mont_v2(felm_t r, const felm_t a, const felm_t b);
void fpsqr_mont_v1(felm_t r, const felm_t a);
void fpsqr_mont_v2(felm_t r, const felm_t a);
void fpinv_chain_mont(felm_t r);
void fpinv_mont(felm_t r);

void mp2_add(f2elm_t r, const felm_t *a, const felm_t *b);
void mp2_sub_p2(f2elm_t r, const felm_t *a, const felm_t *b);

void fp2copy(f2elm_t r, const f2elm_t a);
void fp2zero(f2elm_t r);
void fp2neg(f2elm_t r);
void fp2add(f2elm_t r, const f2elm_t a, const f2elm_t b);
void fp2sub(f2elm_t r, const f2elm_t a, const f2elm_t b);
void fp2div2(f2elm_t r, const f2elm_t a);
void fp2correction(f2elm_t r);
void fp2sqr_mont(f2elm_t r, const f2elm_t a);
void fp2mul_mont(f2elm_t r, const f2elm_t a, const f2elm_t b);
void fp2inv_mont(f2elm_t r);
void to_fp2mont(f2elm_t r, const f2elm_t a);
void from_fp2mont(f2elm_t r, const f2elm_t a);

#endif