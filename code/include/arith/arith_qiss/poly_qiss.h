#ifndef POLY_QISS_H
#define POLY_QISS_H

#include <assert.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include <flint/nmod_poly.h>

#include "params.h"

#ifdef __SIZEOF_INT128__
__extension__ typedef __int128 int128;
__extension__ typedef unsigned __int128 uint128;
#endif

#define COEFFQISS_PACKEDBYTES 8
#define POLYQISS_PACKEDBYTES (PARAM_N_ISS * COEFFQISS_PACKEDBYTES)

typedef nmod_poly_t poly_qiss;

typedef int64_t coeff_qiss;

void arith_qiss_setup(void);
void arith_qiss_teardown(void);

void poly_qiss_init(poly_qiss res);
void poly_qiss_clear(poly_qiss arg);
void poly_qiss_zero(poly_qiss res);
void poly_qiss_set(poly_qiss res, const poly_qiss arg);
coeff_qiss poly_qiss_get_coeff(const poly_qiss arg, size_t n);
coeff_qiss poly_qiss_get_coeff_centered(const poly_qiss arg, size_t n);
void poly_qiss_set_coeff(poly_qiss arg, size_t n, coeff_qiss c);
void poly_qiss_neg(poly_qiss res, const poly_qiss arg);
void poly_qiss_add(poly_qiss res, const poly_qiss lhs, const poly_qiss rhs);
void poly_qiss_neg(poly_qiss res, const poly_qiss arg);
void poly_qiss_sub(poly_qiss res, const poly_qiss lhs, const poly_qiss rhs);
void poly_qiss_mul(poly_qiss res, const poly_qiss lhs, const poly_qiss rhs);
void poly_qiss_mul_x(poly_qiss res, const poly_qiss arg);
void poly_qiss_mul_xj(poly_qiss res, const poly_qiss arg, const size_t j);
void poly_qiss_mul_scalar(poly_qiss out, const poly_qiss lhs, const coeff_qiss rhs);
void poly_qiss_muladd_constant(poly_qiss arg, const coeff_qiss c0_lhs, const coeff_qiss c0_rhs);
void poly_qiss_shift_left(poly_qiss res, const poly_qiss arg, size_t n);
void poly_qiss_conjugate(poly_qiss out, const poly_qiss arg);
void poly_qiss_decompose(poly_qiss high, poly_qiss low, const poly_qiss arg);
void poly_qiss_power2round(poly_qiss high, poly_qiss low, const poly_qiss arg);
void poly_qiss_makeGhint(poly_qiss hint, const poly_qiss z, const poly_qiss r);
void poly_qiss_useGhint(poly_qiss high_zr, const poly_qiss hint, const poly_qiss r);
int poly_qiss_equal(const poly_qiss lhs, const poly_qiss rhs);
void poly_qiss_dump(const poly_qiss arg);
uint128 poly_qiss_sq_norm2(const poly_qiss arg);
int64_t poly_qiss_norm_inf(const poly_qiss arg);
void poly_qiss_pack(uint8_t buf[POLYQISS_PACKEDBYTES], const poly_qiss arg);
void coeff_qiss_pack(uint8_t buf[COEFFQISS_PACKEDBYTES], const coeff_qiss arg);

#endif /* POLY_QISS_H */
