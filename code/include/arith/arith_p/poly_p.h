#ifndef POLY_P_H
#define POLY_P_H

#include <assert.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include <flint/nmod_poly.h>

#include "params.h"

typedef nmod_poly_t poly_p;

typedef int64_t coeff_p;

void arith_p_setup(void);
void arith_p_teardown(void);
void poly_p_init(poly_p res);
void poly_p_clear(poly_p arg);
void poly_p_zero(poly_p res);
void poly_p_set(poly_p res, const poly_p arg);
coeff_p poly_p_get_coeff(const poly_p arg, size_t n);
coeff_p poly_p_get_coeff_centered(const poly_p arg, size_t n);
void poly_p_set_coeff(poly_p arg, size_t n, coeff_p c);
void poly_p_from_bits(poly_p arg, const uint8_t coeffs[PARAM_N / 8]);
void poly_p_neg(poly_p res, const poly_p arg);
void poly_p_add(poly_p res, const poly_p lhs, const poly_p rhs);
void poly_p_sub(poly_p res, const poly_p lhs, const poly_p rhs);
void poly_p_mul(poly_p res, const poly_p lhs, const poly_p rhs);
void poly_p_mul_scalar(poly_p res, const poly_p arg, const coeff_p fac);
int poly_p_equal(const poly_p lhs, const poly_p rhs);
void poly_p_dump(const poly_p arg);
uint64_t poly_p_sq_norm2(const poly_p arg);

#endif /* POLY_P_H */
