#ifndef POLY_QISS_MAT_D_M1_H
#define POLY_QISS_MAT_D_M1_H

#include "arith_qiss.h"

typedef struct {
    poly_qiss_vec_m1 rows[PARAM_D_ISS];
} __poly_qiss_mat_d_m1;

typedef __poly_qiss_mat_d_m1 poly_qiss_mat_d_m1[1];

void poly_qiss_mat_d_m1_setup(void);
void poly_qiss_mat_d_m1_teardown(void);

void poly_qiss_mat_d_m1_init(poly_qiss_mat_d_m1 res);
void poly_qiss_mat_d_m1_clear(poly_qiss_mat_d_m1 res);
void poly_qiss_mat_d_m1_zero(poly_qiss_mat_d_m1 res);
void poly_qiss_mat_d_m1_set(poly_qiss_mat_d_m1 res, const poly_qiss_mat_d_m1 arg);
void poly_qiss_mat_d_m1_neg(poly_qiss_mat_d_m1 res, const poly_qiss_mat_d_m1 arg);
void poly_qiss_mat_d_m1_add(poly_qiss_mat_d_m1 res, const poly_qiss_mat_d_m1 lhs, const poly_qiss_mat_d_m1 rhs);
void poly_qiss_mat_d_m1_sub(poly_qiss_mat_d_m1 res, const poly_qiss_mat_d_m1 lhs, const poly_qiss_mat_d_m1 rhs);
void poly_qiss_mat_d_m1_mul_vec_m1(poly_qiss_vec_d res, const poly_qiss_mat_d_m1 lhs, const poly_qiss_vec_m1 rhs);
int poly_qiss_mat_d_m1_equal(const poly_qiss_mat_d_m1 lhs, const poly_qiss_mat_d_m1 rhs);
void poly_qiss_mat_d_m1_dump(const poly_qiss_mat_d_m1 arg);

#endif /* POLY_QISS_MAT_D_M1_H */
