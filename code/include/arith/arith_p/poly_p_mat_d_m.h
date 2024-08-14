#ifndef POLY_P_MAT_D_M_H
#define POLY_P_MAT_D_M_H

#include "arith_p.h"

typedef struct {
    poly_p_vec_m rows[PARAM_DE];
} __poly_p_mat_d_m;

typedef __poly_p_mat_d_m poly_p_mat_d_m[1];

void poly_p_mat_d_m_setup(void);
void poly_p_mat_d_m_teardown(void);

void poly_p_mat_d_m_init(poly_p_mat_d_m res);
void poly_p_mat_d_m_clear(poly_p_mat_d_m res);
void poly_p_mat_d_m_zero(poly_p_mat_d_m res);
void poly_p_mat_d_m_set(poly_p_mat_d_m res, const poly_p_mat_d_m arg);
void poly_p_mat_d_m_neg(poly_p_mat_d_m res, const poly_p_mat_d_m arg);
void poly_p_mat_d_m_add(poly_p_mat_d_m res, const poly_p_mat_d_m lhs, const poly_p_mat_d_m rhs);
void poly_p_mat_d_m_sub(poly_p_mat_d_m res, const poly_p_mat_d_m lhs, const poly_p_mat_d_m rhs);
void poly_p_mat_d_m_mul_vec_m(poly_p_vec_d res, const poly_p_mat_d_m lhs, const poly_p_vec_m rhs);
int poly_p_mat_d_m_equal(const poly_p_mat_d_m lhs, const poly_p_mat_d_m rhs);
void poly_p_mat_d_m_dump(const poly_p_mat_d_m arg);

#endif /* POLY_P_MAT_D_M_H */
