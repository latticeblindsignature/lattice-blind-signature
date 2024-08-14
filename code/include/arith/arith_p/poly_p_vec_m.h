#ifndef POLY_P_VEC_M_H
#define POLY_P_VEC_M_H

#include "poly_p.h"

typedef struct {
	poly_p entries[PARAM_ME];
} __poly_p_vec_m;

typedef __poly_p_vec_m poly_p_vec_m[1];

void poly_p_vec_m_setup(void);
void poly_p_vec_m_teardown(void);

void poly_p_vec_m_init(poly_p_vec_m arg);
void poly_p_vec_m_clear(poly_p_vec_m arg);
void poly_p_vec_m_zero(poly_p_vec_m arg);
void poly_p_vec_m_set(poly_p_vec_m res, const poly_p_vec_m arg);
void poly_p_vec_m_get_poly(poly_p res, const poly_p_vec_m arg, size_t pos);
void poly_p_vec_m_set_poly(poly_p_vec_m res, const poly_p arg, size_t pos);
void poly_p_vec_m_neg(poly_p_vec_m res, const poly_p_vec_m arg);
void poly_p_vec_m_add(poly_p_vec_m res, const poly_p_vec_m lhs, const poly_p_vec_m rhs);
void poly_p_vec_m_sub(poly_p_vec_m res, const poly_p_vec_m lhs, const poly_p_vec_m rhs);
void poly_p_vec_m_mul_scalar(poly_p_vec_m res, const poly_p_vec_m arg, const poly_p fac);
void poly_p_vec_m_mul_inner(poly_p res, const poly_p_vec_m lhs, const poly_p_vec_m rhs);
int poly_p_vec_m_equal(const poly_p_vec_m lhs, const poly_p_vec_m rhs);
uint64_t poly_p_vec_m_norm2(const poly_p_vec_m arg);
void poly_p_vec_m_dump(const poly_p_vec_m arg);

#endif /* POLY_P_VEC_D_H */
