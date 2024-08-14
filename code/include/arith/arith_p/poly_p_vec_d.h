#ifndef POLY_P_VEC_D_H
#define POLY_P_VEC_D_H

#include "poly_p.h"

typedef struct {
	poly_p entries[PARAM_DE];
} __poly_p_vec_d;

typedef __poly_p_vec_d poly_p_vec_d[1];

void poly_p_vec_d_setup(void);
void poly_p_vec_d_teardown(void);

void poly_p_vec_d_init(poly_p_vec_d arg);
void poly_p_vec_d_clear(poly_p_vec_d arg);
void poly_p_vec_d_zero(poly_p_vec_d arg);
void poly_p_vec_d_set(poly_p_vec_d res, const poly_p_vec_d arg);
void poly_p_vec_d_get_poly(poly_p res, const poly_p_vec_d arg, size_t pos);
void poly_p_vec_d_set_poly(poly_p_vec_d res, const poly_p arg, size_t pos);
void poly_p_vec_d_neg(poly_p_vec_d res, const poly_p_vec_d arg);
void poly_p_vec_d_add(poly_p_vec_d res, const poly_p_vec_d lhs, const poly_p_vec_d rhs);
void poly_p_vec_d_sub(poly_p_vec_d res, const poly_p_vec_d lhs, const poly_p_vec_d rhs);
void poly_p_vec_d_mul_scalar(poly_p_vec_d res, const poly_p_vec_d arg, coeff_p fac);
void poly_p_vec_d_mul_poly(poly_p_vec_d res, const poly_p_vec_d arg, const poly_p arg2);
void poly_p_vec_d_mul_inner(poly_p res, const poly_p_vec_d lhs, const poly_p_vec_d rhs);
int poly_p_vec_d_equal(const poly_p_vec_d lhs, const poly_p_vec_d rhs);
void poly_p_vec_d_dump(const poly_p_vec_d arg);

#endif /* POLY_P_VEC_D_H */
