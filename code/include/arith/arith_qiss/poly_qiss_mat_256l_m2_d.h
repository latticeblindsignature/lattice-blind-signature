#ifndef POLY_QISS_MAT_256L_M2_D_H
#define POLY_QISS_MAT_256L_M2_D_H

#include "arith_qiss.h"

typedef struct {
    poly_qiss_vec_m2_d rows[PARAM_ARP_DIV_N_L_ISS];
} __poly_qiss_mat_256l_m2_d;

typedef __poly_qiss_mat_256l_m2_d poly_qiss_mat_256l_m2_d[1];

void poly_qiss_mat_256l_m2_d_setup(void);
void poly_qiss_mat_256l_m2_d_teardown(void);

void poly_qiss_mat_256l_m2_d_init(poly_qiss_mat_256l_m2_d res);
void poly_qiss_mat_256l_m2_d_clear(poly_qiss_mat_256l_m2_d res);
void poly_qiss_mat_256l_m2_d_zero(poly_qiss_mat_256l_m2_d res);
void poly_qiss_mat_256l_m2_d_set(poly_qiss_mat_256l_m2_d res, const poly_qiss_mat_256l_m2_d arg);
void poly_qiss_mat_256l_m2_d_neg(poly_qiss_mat_256l_m2_d res, const poly_qiss_mat_256l_m2_d arg);
void poly_qiss_mat_256l_m2_d_add(poly_qiss_mat_256l_m2_d res, const poly_qiss_mat_256l_m2_d lhs, const poly_qiss_mat_256l_m2_d rhs);
void poly_qiss_mat_256l_m2_d_sub(poly_qiss_mat_256l_m2_d res, const poly_qiss_mat_256l_m2_d lhs, const poly_qiss_mat_256l_m2_d rhs);
void poly_qiss_mat_256l_m2_d_mul_vec_m2_d(poly_qiss_vec_256_l res, const poly_qiss_mat_256l_m2_d lhs, const poly_qiss_vec_m2_d rhs);
int poly_qiss_mat_256l_m2_d_equal(const poly_qiss_mat_256l_m2_d lhs, const poly_qiss_mat_256l_m2_d rhs);
void poly_qiss_mat_256l_m2_d_dump(const poly_qiss_mat_256l_m2_d arg);

#endif /* POLY_QISS_MAT_256L_M2_D_H */
