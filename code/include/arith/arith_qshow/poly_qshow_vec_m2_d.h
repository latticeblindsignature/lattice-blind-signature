#ifndef POLY_QSHOW_VEC_M2_D_H
#define POLY_QSHOW_VEC_M2_D_H

#include "poly_qshow.h"

typedef struct {
	poly_qshow entries[PARAM_M2_D_SHOW];
} __poly_qshow_vec_m2_d;

typedef __poly_qshow_vec_m2_d poly_qshow_vec_m2_d[1];

void poly_qshow_vec_m2_d_setup(void);
void poly_qshow_vec_m2_d_teardown(void);

void poly_qshow_vec_m2_d_init(poly_qshow_vec_m2_d arg);
void poly_qshow_vec_m2_d_clear(poly_qshow_vec_m2_d arg);
void poly_qshow_vec_m2_d_zero(poly_qshow_vec_m2_d arg);
void poly_qshow_vec_m2_d_set(poly_qshow_vec_m2_d res, const poly_qshow_vec_m2_d arg);
void poly_qshow_vec_m2_d_get_poly(poly_qshow res, const poly_qshow_vec_m2_d arg, size_t pos);
void poly_qshow_vec_m2_d_set_poly(poly_qshow_vec_m2_d res, const poly_qshow arg, size_t pos);
void poly_qshow_vec_m2_d_neg(poly_qshow_vec_m2_d res, const poly_qshow_vec_m2_d arg);
void poly_qshow_vec_m2_d_add(poly_qshow_vec_m2_d res, const poly_qshow_vec_m2_d lhs, const poly_qshow_vec_m2_d rhs);
void poly_qshow_vec_m2_d_sub(poly_qshow_vec_m2_d res, const poly_qshow_vec_m2_d lhs, const poly_qshow_vec_m2_d rhs);
void poly_qshow_vec_m2_d_mul_scalar(poly_qshow_vec_m2_d res, const poly_qshow_vec_m2_d arg, const coeff_qshow fac);
void poly_qshow_vec_m2_d_mul_poly_qshow(poly_qshow_vec_m2_d res, const poly_qshow_vec_m2_d lhs, const poly_qshow rhs);
void poly_qshow_vec_m2_d_mul_inner(poly_qshow res, const poly_qshow_vec_m2_d lhs, const poly_qshow_vec_m2_d rhs);
uint128 poly_qshow_vec_m2_d_norm2(const poly_qshow_vec_m2_d arg);
int poly_qshow_vec_m2_d_equal(const poly_qshow_vec_m2_d lhs, const poly_qshow_vec_m2_d rhs);
void poly_qshow_vec_m2_d_dump(const poly_qshow_vec_m2_d arg);

#endif /* POLY_QSHOW_VEC_M2_D_H */
