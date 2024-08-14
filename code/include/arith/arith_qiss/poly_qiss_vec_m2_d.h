#ifndef POLY_QISS_VEC_M2_D_H
#define POLY_QISS_VEC_M2_D_H

#include "poly_qiss.h"

typedef struct {
	poly_qiss entries[PARAM_M2_D_ISS];
} __poly_qiss_vec_m2_d;

typedef __poly_qiss_vec_m2_d poly_qiss_vec_m2_d[1];

void poly_qiss_vec_m2_d_setup(void);
void poly_qiss_vec_m2_d_teardown(void);

void poly_qiss_vec_m2_d_init(poly_qiss_vec_m2_d arg);
void poly_qiss_vec_m2_d_clear(poly_qiss_vec_m2_d arg);
void poly_qiss_vec_m2_d_zero(poly_qiss_vec_m2_d arg);
void poly_qiss_vec_m2_d_set(poly_qiss_vec_m2_d res, const poly_qiss_vec_m2_d arg);
void poly_qiss_vec_m2_d_get_poly(poly_qiss res, const poly_qiss_vec_m2_d arg, size_t pos);
void poly_qiss_vec_m2_d_set_poly(poly_qiss_vec_m2_d res, const poly_qiss arg, size_t pos);
void poly_qiss_vec_m2_d_neg(poly_qiss_vec_m2_d res, const poly_qiss_vec_m2_d arg);
void poly_qiss_vec_m2_d_add(poly_qiss_vec_m2_d res, const poly_qiss_vec_m2_d lhs, const poly_qiss_vec_m2_d rhs);
void poly_qiss_vec_m2_d_sub(poly_qiss_vec_m2_d res, const poly_qiss_vec_m2_d lhs, const poly_qiss_vec_m2_d rhs);
void poly_qiss_vec_m2_d_mul_scalar(poly_qiss_vec_m2_d res, const poly_qiss_vec_m2_d arg, const coeff_qiss fac);
void poly_qiss_vec_m2_d_mul_poly_qiss(poly_qiss_vec_m2_d res, const poly_qiss_vec_m2_d lhs, const poly_qiss rhs);
void poly_qiss_vec_m2_d_mul_inner(poly_qiss res, const poly_qiss_vec_m2_d lhs, const poly_qiss_vec_m2_d rhs);
uint128 poly_qiss_vec_m2_d_norm2(const poly_qiss_vec_m2_d arg);
int poly_qiss_vec_m2_d_equal(const poly_qiss_vec_m2_d lhs, const poly_qiss_vec_m2_d rhs);
void poly_qiss_vec_m2_d_dump(const poly_qiss_vec_m2_d arg);

#endif /* POLY_QISS_VEC_M2_D_H */
