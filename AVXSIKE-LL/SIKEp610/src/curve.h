/**
 *******************************************************************************
 * @version 0.1
 * @date 2021-12-16
 * @copyright Copyright © 2021 by University of Luxembourg
 * @author Hao Cheng
 *******************************************************************************
 */

#ifndef _CURVE_H
#define _CURVE_H

#include "fpx.h"

// -----------------------------------------------------------------------------
// (8x1x1x1)-way curve/isogeny arithmetic.

typedef struct { vf2elm_t X; vf2elm_t Z; } vpoint_proj;
typedef vpoint_proj vpoint_proj_t[1]; 

void xDBL_8x1x1x1w(const vpoint_proj_t P, vpoint_proj_t Q, const vf2elm_t A24plus, const vf2elm_t C24);
void get_4_isog_8x1x1x1w(const vpoint_proj_t P, vf2elm_t A24plus, vf2elm_t C24, vf2elm_t *coeff);
void eval_4_isog_8x1x1x1w(vpoint_proj_t P, vf2elm_t *coeff);

void xTPL_8x1x1x1w(const vpoint_proj_t P, vpoint_proj_t Q, const vf2elm_t A24minus, const vf2elm_t A24plus);
void get_3_isog_8x1x1x1w(const vpoint_proj_t P, vf2elm_t A24minus, vf2elm_t A24plus, vf2elm_t* coeff);
void eval_3_isog_8x1x1x1w(vpoint_proj_t Q, const vf2elm_t *coeff);

void xDBLADD_8x1x1x1w(vpoint_proj_t P, vpoint_proj_t Q, const vf2elm_t XPQ, const vf2elm_t ZPQ, const vf2elm_t A24);

// -----------------------------------------------------------------------------
// (1x4x2x1)-way curve/isogeny arithmetic.

void xDBLADD_1x4x2x1w(vfelm_t x3z3x2z2, vfelm_t z1x1_A);

// -----------------------------------------------------------------------------
// (2x4x1x1)-way curve/isogeny arithmetic.

void xDBLADD_2x4x1x1w(vf2elm_t x3z3x2z2, vf2elm_t z1x1_A);

// -----------------------------------------------------------------------------
// (1x2x2x2)-way curve/isogeny arithmetic.

void xDBL_1x2x2x2w(const vgelm_t P, vgelm_t Q, const vgelm_t C24_A24plus);
void xTPL_1x2x2x2w(const vgelm_t P, vgelm_t Q, const vgelm_t C24_A24plus);
void get_4_isog_1x2x2x2w(const vgelm_t P, vgelm_t C24_A24plus, vgelm_t coeff__0, vgelm_t coeff2_1);
void get_3_isog_1x2x2x2w(const vgelm_t P, vgelm_t C24_A24plus, vgelm_t coeff1_0);
void eval_4_isog_1x2x2x2w(vgelm_t P, const vgelm_t *coeff);
void eval_3_isog_1x2x2x2w(vgelm_t Q, vgelm_t coeff0_1);

void get_2_isog_1x2x2x2w(const vgelm_t P, vgelm_t C24_A24plus);
void eval_2_isog_1x2x2x2w(vgelm_t P, vgelm_t Q);

// -----------------------------------------------------------------------------
// (4x2x1x1)-way curve/isogeny arithmetic.

void eval_4_isog_4x2x1x1w(vf2elm_t P, vf2elm_t *coeff);
void eval_3_isog_4x2x1x1w(vf2elm_t Q, const vf2elm_t coeff0_1);

// -----------------------------------------------------------------------------
// (2x2x2x1)-way curve/isogeny arithmetic.

void eval_4_isog_2x2x2x1w(vfelm_t P, vfelm_t *coeff);
void eval_3_isog_2x2x2x1w(vfelm_t Q, vfelm_t coeff0_1);
void xDBL_2x2x2x1w(const vfelm_t P, vfelm_t Q, const vfelm_t C24_A24plus);
void get_4_isog_2x2x2x1w(const vfelm_t P, vfelm_t C24_A24plus, vfelm_t coeff__0, vfelm_t coeff2_1);

void get_2_isog_2x2x2x1w(const vfelm_t P, vfelm_t C24_A24plus);
void eval_2_isog_2x2x2x1w(vfelm_t P, vfelm_t Q);

// -----------------------------------------------------------------------------
// 1-way x64 curve/isogeny arithmetic.

typedef struct { f2elm_t X; f2elm_t Z; } point_proj;               
typedef point_proj point_proj_t[1]; 

typedef struct { f2elm_r51_t X; f2elm_r51_t Z; } point_proj_r51;               
typedef point_proj_r51 point_proj_r51_t[1]; 

void pointcopy_1w(point_proj_r51_t Q, const point_proj_r51_t P);
void LADDER3PT_1x4x2x1w(const f2elm_r51_t xP, const f2elm_r51_t xQ, const f2elm_r51_t xPQ, const digit_t* m, const unsigned int AliceOrBob, point_proj_r51_t R, const f2elm_r51_t A);
void carryp_1w(uint64_t *a);

void LADDER3PT_2x4x1x1w(const f2elm_r51_t xP, const f2elm_r51_t xQ, const f2elm_r51_t xPQ, \
                        const f2elm_r51_t _xP, const f2elm_r51_t _xQ, const f2elm_r51_t _xPQ, \
                        const digit_t* m, const unsigned int AliceOrBob, \
                        point_proj_r51_t R, const f2elm_r51_t A, \
                        point_proj_r51_t _R, const f2elm_r51_t _A);

void xDBLe_1x2x2x2w(vgelm_t vP, point_proj_r51_t Q, const vgelm_t C24_A24plus, const int e);
void xDBLe_2x2x2x1w(vfelm_t vP, point_proj_r51_t Q, point_proj_r51_t _Q, const vfelm_t C24_A24plus, const int e);
void eval_4_isog_parallel_kg(point_proj_r51_t *pts, point_proj_r51_t phiP, point_proj_r51_t phiQ, point_proj_r51_t phiR, const vgelm_t coeff__0, const vgelm_t coeff2_1, const int num);
void eval_4_isog_parallel_ss(point_proj_r51_t *pts, const vgelm_t coeff__0, const vgelm_t coeff2_1, const int num);
void eval_4_isog_parallel_kgss(point_proj_r51_t *pts, point_proj_r51_t *_pts, point_proj_r51_t phiP, point_proj_r51_t phiQ, point_proj_r51_t phiR, const vfelm_t coeff__0, const vfelm_t coeff2_1, const int num);

void xTPLe_1x2x2x2w(vgelm_t vP, point_proj_r51_t Q, const vgelm_t C24_A24plus, const int e);
void eval_3_isog_parallel_kg(point_proj_r51_t *pts, point_proj_r51_t phiP, point_proj_r51_t phiQ, point_proj_r51_t phiR, const vgelm_t coeff1_0, const int num);
void eval_3_isog_parallel_ss(point_proj_r51_t *pts, const vgelm_t coeff1_0, const int num);

void inv_3_way(f2elm_t z1, f2elm_t z2, f2elm_t z3);
void get_A(const f2elm_t xP, const f2elm_t xQ, const f2elm_t xR, f2elm_t A);
void j_inv(const f2elm_t A, const f2elm_t C, f2elm_t jinv);

static void get_channel_8x1w(felm_r51_t r, const vfelm_t a, const int ch) 
{
  int i;

  for(i = 0; i < VNWORDS; i++) {
    r[i] = ((uint64_t *)&a[i])[ch];
  }
}

static void get_channel_4x2w(felm_r51_t r, const vgelm_t a, const int ch) 
{
  int i;

  for(i = 0; i < VGWORDS; i++) {
    r[i] = ((uint64_t *)&a[i])[ch];
    r[i+VGWORDS] = ((uint64_t *)&a[i])[ch+1];
  }
}

#endif
