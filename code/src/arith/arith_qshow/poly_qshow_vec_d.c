#include "arith_qshow.h"
#include "macros.h"

static poly_qshow TMP;

/*************************************************
* Name:        poly_qshow_vec_d_setup
*
* Description: Initialize and setup the backend for arithmetic 
*              mod PARAM_Q_SHOW integer vectors with PARAM_D_SHOW entries. 
* 			   		 This is strictly required and must be called once 
* 			   		 before any other function from here is used.
**************************************************/
void poly_qshow_vec_d_setup(void) {
  poly_qshow_init(TMP);
}

/*************************************************
* Name:        poly_qshow_vec_d_teardown
*
* Description: Clean up and teardown the backend for arithmetic 
*              mod PARAM_Q_SHOW integer vectors with PARAM_D_SHOW entries. 
* 			   		 This is strictly required and must be called once 
* 			   		 at the very end to release any resources.
**************************************************/
void poly_qshow_vec_d_teardown(void) {
  poly_qshow_clear(TMP);
}

/*************************************************
* Name:        poly_qshow_vec_d_init
*
* Description: Initialize polynomial vector with PARAM_D_SHOW entries.
*              This is strictly required before any operations 
*              are done with/on the vector.
* 
* Arguments:   - poly_qshow_vec_d arg: polynomial vector to be initialized
**************************************************/
void poly_qshow_vec_d_init(poly_qshow_vec_d arg) {
	for (size_t i = 0; i < PARAM_D_SHOW; ++i) {
		poly_qshow_init(arg->entries[i]);
	}
}

/*************************************************
* Name:        poly_qshow_vec_d_clear
*
* Description: Clear polynomial vector with PARAM_D_SHOW entries.
*              This is strictly required to avoid memory leaks and the 
*              polynomial vector must not be used again (unless reinitialized).
* 
* Arguments:   - poly_qshow_vec_d arg: polynomial vector to be cleared
**************************************************/
void poly_qshow_vec_d_clear(poly_qshow_vec_d arg) {
	for (size_t i = 0; i < PARAM_D_SHOW; ++i) {
		poly_qshow_clear(arg->entries[i]);
	}
}

/*************************************************
* Name:        poly_qshow_vec_d_zero
*
* Description: Set an initialized polynomial vector with PARAM_D_SHOW entries to zero
* 
* Arguments:   - poly_qshow_vec_d arg: polynomial vector to be zeroized (initialized)
**************************************************/
void poly_qshow_vec_d_zero(poly_qshow_vec_d arg) {
	for (size_t i = 0; i < PARAM_D_SHOW; ++i) {
		poly_qshow_zero(arg->entries[i]);
	}
}

/*************************************************
* Name:        poly_qshow_vec_d_set
*
* Description: Set a polynomial vector with PARAM_D_SHOW entries equal to another polynomial vector
* 
* Arguments:   - poly_qshow_vec_d res: polynomial vector to be set (initialized)
* 			   		 - const poly_qshow_vec_d arg: polynomial vector to be read
**************************************************/
void poly_qshow_vec_d_set(poly_qshow_vec_d res, const poly_qshow_vec_d arg) {
	for (size_t i = 0; i < PARAM_D_SHOW; ++i) {
		poly_qshow_set(res->entries[i], arg->entries[i]);
	}
}

/*************************************************
* Name:        poly_qshow_vec_d_get_poly
*
* Description: Get pos-th entry polynomial of the vector
*              condition: [0 <= pos < PARAM_D_SHOW]
* 
* Arguments:   - poly_qshow res: polynomial to host the pos-th entry (initialized)
* 			   		 - const poly_qshow_vec_d arg: polynomial vector to be read
* 			   		 - size_t pos: position to get in the vector
**************************************************/
void poly_qshow_vec_d_get_poly(poly_qshow res, const poly_qshow_vec_d arg, size_t pos) {
	ASSERT_DEBUG(pos < PARAM_D_SHOW, "Illegal argument: cannot get entry of vector at given position.");
	poly_qshow_set(res, arg->entries[pos]);
}

/*************************************************
* Name:        poly_qshow_vec_d_set_poly
*
* Description: Get pos-th entry polynomial of the vector
*              condition: [0 <= pos < PARAM_D_SHOW]
* 
* Arguments:   - poly_qshow_vec_d res: polynomial vector to be set (initialized)
* 			   		 - const poly_qshow arg: polynomial to set pos-th entry
* 			   		 - size_t pos: position to get in the vector
**************************************************/
void poly_qshow_vec_d_set_poly(poly_qshow_vec_d res, const poly_qshow arg, size_t pos) {
	ASSERT_DEBUG(pos < PARAM_D_SHOW, "Illegal argument: cannot set entry of vector at given position.");
	poly_qshow_set(res->entries[pos], arg);
}

/*************************************************
* Name:        poly_qshow_vec_d_neg
*
* Description: Negate a polynomial vector with PARAM_D_SHOW entries
* 
* Arguments:   - poly_qshow_vec_d res: polynomial vector to host the negation (initialized)
* 			   		 - const poly_qshow_vec_d arg: polynomial vector to be negated
**************************************************/
void poly_qshow_vec_d_neg(poly_qshow_vec_d res, const poly_qshow_vec_d arg) {
  for (size_t i = 0; i < PARAM_D_SHOW; ++i) {
    poly_qshow_neg(res->entries[i], arg->entries[i]);
  }
}

/*************************************************
* Name:        poly_qshow_vec_d_add
*
* Description: Add two polynomial vectors with PARAM_D_SHOW entries
* 
* Arguments:   - poly_qshow_vec_d res: polynomial vector to host the sum (initialized)
* 			   		 - const poly_qshow_vec_d lhs: first polynomial vector summand
* 			   		 - const poly_qshow_vec_d rhs: second polynomial vector summand
**************************************************/
void poly_qshow_vec_d_add(poly_qshow_vec_d res, const poly_qshow_vec_d lhs, const poly_qshow_vec_d rhs) {
	for (size_t i = 0; i < PARAM_D_SHOW; ++i) {
		poly_qshow_add(res->entries[i], lhs->entries[i], rhs->entries[i]);
	}
}

/*************************************************
* Name:        poly_qshow_vec_d_sub
*
* Description: Substract two polynomial vectors with PARAM_D_SHOW entries
* 
* Arguments:   - poly_qshow_vec_d res: polynomial vector to host the difference (initialized)
* 			   		 - const poly_qshow_vec_d lhs: first polynomial vector term
* 			   		 - const poly_qshow_vec_d rhs: second polynomial vector term
**************************************************/
void poly_qshow_vec_d_sub(poly_qshow_vec_d res, const poly_qshow_vec_d lhs, const poly_qshow_vec_d rhs) {
	for (size_t i = 0; i < PARAM_D_SHOW; ++i) {
		poly_qshow_sub(res->entries[i], lhs->entries[i], rhs->entries[i]);
	}
}

/*************************************************
* Name:        poly_qshow_vec_d_mul_scalar
*
* Description: Multiplication of a polynomial vector with PARAM_D_SHOW entries by a integer scalar
* 
* Arguments:   - poly_qshow_vec_d res: polynomial vector to host the multiplication (initialized)
* 			   		 - const poly_qshow_vec_d arg: polynomial vector factor
* 			   		 - coeff_qshow fac: integer factor
**************************************************/
void poly_qshow_vec_d_mul_scalar(poly_qshow_vec_d res, const poly_qshow_vec_d arg, const coeff_qshow fac) {
	for (size_t i = 0; i < PARAM_D_SHOW; ++i) {
		poly_qshow_mul_scalar(res->entries[i], arg->entries[i], fac);
	}
}

/*************************************************
* Name:        poly_qshow_vec_d_mul_poly_qshow
*
* Description: Multiplication of a polynomial vector with PARAM_D_SHOW entries by a polynomial
* 
* Arguments:   - poly_qshow_vec_d res: polynomial vector to host the multiplication (initialized)
* 			   		 - const poly_qshow_vec_d lhs: first polynomial vector factor
* 			   		 - const poly_qshow rhs: second polynomial factor
**************************************************/
void poly_qshow_vec_d_mul_poly_qshow(poly_qshow_vec_d res, const poly_qshow_vec_d lhs, const poly_qshow rhs) {
	for (size_t i = 0; i < PARAM_D_SHOW; i++) {
		poly_qshow_mul(res->entries[i], lhs->entries[i], rhs);
	}
}

/*************************************************
* Name:        poly_qshow_vec_d_mul_inner
*
* Description: Inner product of two polynomial vectors with PARAM_D_SHOW entries
* 
* Arguments:   - poly_qshow res: polynomial to host the inner product (initialized)
* 			   		 - const poly_qshow_vec_d lhs: first polynomial vector
* 			   		 - const poly_qshow_vec_d rhs: second polynomial vector
**************************************************/
void poly_qshow_vec_d_mul_inner(poly_qshow res, const poly_qshow_vec_d lhs, const poly_qshow_vec_d rhs) {
	poly_qshow_zero(res);
	for (size_t i = 0; i < PARAM_D_SHOW; ++i) {
		poly_qshow_mul(TMP, lhs->entries[i], rhs->entries[i]);
		poly_qshow_add(res, res, TMP);
	}
}

/*************************************************
* Name:        poly_qshow_vec_d_decompose
*
* Description: Compute the decomposition x_L + PARAM_GAMMA_SHOW.x_H of a polynomial vector x
*              with x_L with coefficients in (-PARAM_GAMMA_SHOW/2, PARAM_GAMMA_SHOW/2]
*              except if coefficients of x_H = (PARAM_Q_SHOW-1)/PARAM_GAMMA_SHOW where we set 
*              the coefficient of x_H to 0 and that of x_L to x mod+ PARAM_Q_SHOW
* 
* Arguments:   - poly_qshow_vec_d high: polynomial vector to host the high order decomposition (initialized)
*              - poly_qshow_vec_d low: polynomial vector to host the low order decomposition (initialized)
*              - const poly_qshow_vec_d arg: polynomial vector to be decomposed
**************************************************/
void poly_qshow_vec_d_decompose(poly_qshow_vec_d high, poly_qshow_vec_d low, const poly_qshow_vec_d arg) {
  for (size_t i = 0; i < PARAM_D_SHOW; i++) {
  	poly_qshow_decompose(high->entries[i], low->entries[i], arg->entries[i]);
  }
}

/*************************************************
* Name:        poly_qshow_vec_d_power2round
*
* Description: Compute the decomposition x_L + 2^PARAM_D_ROUND_SHOW.x_H of a polynomial vector x
*              with x_L with coefficients in (-2^(PARAM_D_ROUND_SHOW-1), 2^(PARAM_D_ROUND_SHOW-1)]
* 
* Arguments:   - poly_qshow_vec_d high: polynomial vector to host the high order decomposition (initialized)
*              - poly_qshow_vec_d low: polynomial vector to host the low order decomposition (initialized)
*              - const poly_qshow_vec_d arg: polynomial vector to be decomposed
**************************************************/
void poly_qshow_vec_d_power2round(poly_qshow_vec_d high, poly_qshow_vec_d low, const poly_qshow_vec_d arg) {
  for (size_t i = 0; i < PARAM_D_SHOW; i++) {
  	poly_qshow_power2round(high->entries[i], low->entries[i], arg->entries[i]);
  }
}

/*************************************************
* Name:        poly_qshow_vec_d_makeGhint
*
* Description: 
* 
* Arguments:   - poly_qshow_vec_d hint: polynomial vector to host the hint polynomial (initialized)
*              - const poly_qshow_vec_d z: polynomial vector 
*              - const poly_qshow_vec_d r: polynomial vector 
**************************************************/
void poly_qshow_vec_d_makeGhint(poly_qshow_vec_d hint, const poly_qshow_vec_d z, const poly_qshow_vec_d r) {
  for (size_t i = 0; i < PARAM_D_SHOW; i++) {
  	poly_qshow_makeGhint(hint->entries[i], z->entries[i], r->entries[i]);
  }
}

/*************************************************
* Name:        poly_qshow_vec_d_useGhint
*
* Description: 
* 
* Arguments:   - poly_qshow_vec_d high_zr: polynomial vector to host the hint polynomial (initialized)
*              - const poly_qshow_vec_d hint: polynomial vector 
*              - const poly_qshow_vec_d r: polynomial vector 
**************************************************/
void poly_qshow_vec_d_useGhint(poly_qshow_vec_d high_zr, const poly_qshow_vec_d hint, const poly_qshow_vec_d r) {
  for (size_t i = 0; i < PARAM_D_SHOW; i++) {
  	poly_qshow_useGhint(high_zr->entries[i], hint->entries[i], r->entries[i]);
  }
}

/*************************************************
* Name:        poly_qshow_vec_d_norm2
*
* Description: Compute the square l2 norm of a polynomial vector
* 
* Arguments:   - const poly_qshow_vec_d arg: the polynomial vector
* 
* Returns an unsigned 64-bit integer with the square l2 norm
**************************************************/
uint128 poly_qshow_vec_d_norm2(const poly_qshow_vec_d arg) {
  uint128 sq_norm2 = 0;
  for (size_t i = 0; i < PARAM_D_SHOW; i++) {
  	sq_norm2 += poly_qshow_sq_norm2(arg->entries[i]);
  }
  return sq_norm2;
}

/*************************************************
* Name:        poly_qshow_vec_d_infnorm
*
* Description: Compute the infinity norm of a polynomial vector
*
* Arguments:   - const poly_qshow_vec_d arg: the polynomial vector
* 
* Returns an 64-bit integer with the linf norm
**************************************************/
int64_t poly_qshow_vec_d_norm_inf(const poly_qshow_vec_d arg) {
 	int64_t max, tmpmax;
  max = poly_qshow_norm_inf(arg->entries[0]);
  for (size_t i = 1; i < PARAM_D_SHOW; i++) {
    tmpmax = poly_qshow_norm_inf(arg->entries[i]);
    max = (tmpmax > max) ? tmpmax : max;
  }
  return max;
}

/*************************************************
* Name:        poly_qshow_vec_d_equal
*
* Description: Equality test between two polynomial vectors with PARAM_D_SHOW entries
* 
* Arguments:   - const poly_qshow_vec_d lhs: first polynomial vector
* 			   		 - const poly_qshow_vec_d rhs: second polynomial vector
* 
* Returns 1 if the polynomial vectors are equal, 0 otherwise
**************************************************/
int poly_qshow_vec_d_equal(const poly_qshow_vec_d lhs, const poly_qshow_vec_d rhs) {
  for (size_t i = 0; i < PARAM_D_SHOW; ++i) {
    if (!poly_qshow_equal(lhs->entries[i], rhs->entries[i])) {
      return 0;
    }
  }
  return 1;
}

/*************************************************
* Name:        poly_qshow_vec_d_dump
*
* Description: Print a polynomial vector with PARAM_D_SHOW entries
* 
* Arguments:   - const poly_qshow_vec_d arg: polynomial vector to be printed
**************************************************/
void poly_qshow_vec_d_dump(const poly_qshow_vec_d arg) {
	printf("[");
	for (size_t i = 0; i < PARAM_D_SHOW - 1; ++i) {
		poly_qshow_dump(arg->entries[i]);
		printf(", ");
	}
	poly_qshow_dump(arg->entries[PARAM_D_SHOW - 1]);
	printf("]");
}

/*************************************************
* Name:        poly_qshow_vec_d_pack
*
* Description: Pack a polynomial vector mod PARAM_Q_SHOW into a byte array
* 
* Arguments:   - uint8_t buf: output byte array (allocated POLYQSHOW_VECD_PACKEDBYTES bytes)
*              - const poly_qshow arg: the polynomial vector to be packed
**************************************************/
void poly_qshow_vec_d_pack(uint8_t buf[POLYQSHOW_VECD_PACKEDBYTES], const poly_qshow_vec_d arg) {
	for (size_t i = 0; i < PARAM_D_SHOW; i++) {
		poly_qshow_pack(&buf[i * POLYQSHOW_PACKEDBYTES], arg->entries[i]);
	}
}
