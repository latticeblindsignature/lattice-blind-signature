#ifndef POLY_QISS_VEC_M1_H
#define POLY_QISS_VEC_M1_H

#include "poly_qiss.h"

typedef struct {
	poly_qiss entries[PARAM_M1_ISS];
} __poly_qiss_vec_m1;

typedef __poly_qiss_vec_m1 poly_qiss_vec_m1[1];

void poly_qiss_vec_m1_setup(void);
void poly_qiss_vec_m1_teardown(void);

void poly_qiss_vec_m1_init(poly_qiss_vec_m1 arg);
void poly_qiss_vec_m1_clear(poly_qiss_vec_m1 arg);
void poly_qiss_vec_m1_zero(poly_qiss_vec_m1 arg);
void poly_qiss_vec_m1_set(poly_qiss_vec_m1 res, const poly_qiss_vec_m1 arg);
void poly_qiss_vec_m1_get_poly(poly_qiss res, const poly_qiss_vec_m1 arg, size_t pos);
void poly_qiss_vec_m1_set_poly(poly_qiss_vec_m1 res, const poly_qiss arg, size_t pos);
void poly_qiss_vec_m1_neg(poly_qiss_vec_m1 res, const poly_qiss_vec_m1 arg);
void poly_qiss_vec_m1_add(poly_qiss_vec_m1 res, const poly_qiss_vec_m1 lhs, const poly_qiss_vec_m1 rhs);
void poly_qiss_vec_m1_sub(poly_qiss_vec_m1 res, const poly_qiss_vec_m1 lhs, const poly_qiss_vec_m1 rhs);
void poly_qiss_vec_m1_mul_scalar(poly_qiss_vec_m1 res, const poly_qiss_vec_m1 arg, const coeff_qiss fac);
void poly_qiss_vec_m1_mul_poly_qiss(poly_qiss_vec_m1 res, const poly_qiss_vec_m1 lhs, const poly_qiss rhs);
void poly_qiss_vec_m1_mul_inner(poly_qiss res, const poly_qiss_vec_m1 lhs, const poly_qiss_vec_m1 rhs);
void poly_qiss_vec_m1_conjugate(poly_qiss_vec_m1 res, const poly_qiss_vec_m1 arg);
uint128 poly_qiss_vec_m1_norm2(const poly_qiss_vec_m1 arg);
int poly_qiss_vec_m1_equal(const poly_qiss_vec_m1 lhs, const poly_qiss_vec_m1 rhs);
void poly_qiss_vec_m1_dump(const poly_qiss_vec_m1 arg);

#endif /* POLY_QISS_VEC_M1_H */
