#include "arith_p.h"
#include "random.h"
#include "macros.h"

static poly_p TMP;

/*************************************************
* Name:        poly_p_vec_d_setup
*
* Description: Initialize and setup the backend for arithmetic 
*              mod PARAM_Q integer vectors with PARAM_DE entries. 
* 			   		 This is strictly required and must be called once 
* 			   		 before any other function from here is used.
**************************************************/
void poly_p_vec_d_setup(void) {
  poly_p_init(TMP);
}

/*************************************************
* Name:        poly_p_vec_d_teardown
*
* Description: Clean up and teardown the backend for arithmetic 
*              mod PARAM_Q integer vectors with PARAM_DE entries. 
* 			   		 This is strictly required and must be called once 
* 			   		 at the very end to release any resources.
**************************************************/
void poly_p_vec_d_teardown(void) {
  poly_p_clear(TMP);
}

/*************************************************
* Name:        poly_p_vec_d_init
*
* Description: Initialize polynomial vector with PARAM_DE entries.
*              This is strictly required before any operations 
*              are done with/on the vector.
* 
* Arguments:   - poly_p_vec_d arg: polynomial vector to be initialized
**************************************************/
void poly_p_vec_d_init(poly_p_vec_d arg) {
	for (size_t i = 0; i < PARAM_DE; ++i) {
		poly_p_init(arg->entries[i]);
	}
}

/*************************************************
* Name:        poly_p_vec_d_clear
*
* Description: Clear polynomial vector with PARAM_DE entries.
*              This is strictly required to avoid memory leaks and the 
*              polynomial vector must not be used again (unless reinitialized).
* 
* Arguments:   - poly_p_vec_d arg: polynomial vector to be cleared
**************************************************/
void poly_p_vec_d_clear(poly_p_vec_d arg) {
	for (size_t i = 0; i < PARAM_DE; ++i) {
		poly_p_clear(arg->entries[i]);
	}
}

/*************************************************
* Name:        poly_p_vec_d_zero
*
* Description: Set an initialized polynomial vector with PARAM_DE entries to zero
* 
* Arguments:   - poly_p_vec_d arg: polynomial vector to be zeroized (initialized)
**************************************************/
void poly_p_vec_d_zero(poly_p_vec_d arg) {
	for (size_t i = 0; i < PARAM_DE; ++i) {
		poly_p_zero(arg->entries[i]);
	}
}

/*************************************************
* Name:        poly_p_vec_d_set
*
* Description: Set a polynomial vector with PARAM_DE entries equal to another polynomial vector
* 
* Arguments:   - poly_p_vec_d res: polynomial vector to be set (initialized)
* 			   		 - const poly_p_vec_d arg: polynomial vector to be read
**************************************************/
void poly_p_vec_d_set(poly_p_vec_d res, const poly_p_vec_d arg) {
	for (size_t i = 0; i < PARAM_DE; ++i) {
		poly_p_set(res->entries[i], arg->entries[i]);
	}
}

/*************************************************
* Name:        poly_p_vec_d_get_poly
*
* Description: Get pos-th entry polynomial of the vector
*              condition: [0 <= pos < PARAM_DE]
* 
* Arguments:   - poly_p res: polynomial to host the pos-th entry (initialized)
* 			   		 - const poly_p_vec_d arg: polynomial vector to be read
* 			   		 - size_t pos: position to get in the vector
**************************************************/
void poly_p_vec_d_get_poly(poly_p res, const poly_p_vec_d arg, size_t pos) {
  ASSERT_DEBUG(pos < PARAM_DE, "Illegal argument: cannot get entry of vector at given position.");
	poly_p_set(res, arg->entries[pos]);
}

/*************************************************
* Name:        poly_p_vec_d_set_poly
*
* Description: Get pos-th entry polynomial of the vector
*              condition: [0 <= pos < PARAM_DE]
* 
* Arguments:   - poly_p_vec_d res: polynomial vector to be set (initialized)
* 			   		 - const poly_p arg: polynomial to set pos-th entry
* 			   		 - size_t pos: position to get in the vector
**************************************************/
void poly_p_vec_d_set_poly(poly_p_vec_d res, const poly_p arg, size_t pos) {
  ASSERT_DEBUG(pos < PARAM_DE, "Illegal argument: cannot set entry of vector at given position.");
	poly_p_set(res->entries[pos], arg);
}

/*************************************************
* Name:        poly_p_vec_d_neg
*
* Description: Negate a polynomial vector with PARAM_DE entries
* 
* Arguments:   - poly_p_vec_d res: polynomial vector to host the negation (initialized)
* 			   		 - const poly_p_vec_d arg: polynomial vector to be negated
**************************************************/
void poly_p_vec_d_neg(poly_p_vec_d res, const poly_p_vec_d arg) {
  for (size_t i = 0; i < PARAM_DE; ++i) {
  	poly_p_neg(res->entries[i], arg->entries[i]);
  }
}

/*************************************************
* Name:        poly_p_vec_d_add
*
* Description: Add two polynomial vectors with PARAM_DE entries
* 
* Arguments:   - poly_p_vec_d res: polynomial vector to host the sum (initialized)
* 			   		 - const poly_p_vec_d lhs: first polynomial vector summand
* 			   		 - const poly_p_vec_d rhs: second polynomial vector summand
**************************************************/
void poly_p_vec_d_add(poly_p_vec_d res, const poly_p_vec_d lhs, const poly_p_vec_d rhs) {
	for (size_t i = 0; i < PARAM_DE; ++i) {
		poly_p_add(res->entries[i], lhs->entries[i], rhs->entries[i]);
	}
}

/*************************************************
* Name:        poly_p_vec_d_sub
*
* Description: Substract two polynomial vectors with PARAM_DE entries
* 
* Arguments:   - poly_p_vec_d res: polynomial vector to host the difference (initialized)
* 			   		 - const poly_p_vec_d lhs: first polynomial vector term
* 			   		 - const poly_p_vec_d rhs: second polynomial vector term
**************************************************/
void poly_p_vec_d_sub(poly_p_vec_d res, const poly_p_vec_d lhs, const poly_p_vec_d rhs) {
	for (size_t i = 0; i < PARAM_DE; ++i) {
		poly_p_sub(res->entries[i], lhs->entries[i], rhs->entries[i]);
	}
}

/*************************************************
* Name:        poly_p_vec_d_mul_scalar
*
* Description: Multiplication of a polynomial vector with PARAM_DE entries by a integer scalar
* 
* Arguments:   - poly_p_vec_d res: polynomial vector to host the multiplication (initialized)
* 			   		 - const poly_p_vec_d arg: polynomial vector factor
* 			   		 - coeff_q fac: integer factor
**************************************************/
void poly_p_vec_d_mul_scalar(poly_p_vec_d res, const poly_p_vec_d arg, coeff_p fac) {
	for (size_t i = 0; i < PARAM_DE; ++i) {
		poly_p_mul_scalar(res->entries[i], arg->entries[i], fac);
	}
}

/*************************************************
* Name:        poly_p_vec_d_mul_poly
*
* Description: Multiplication of a polynomial vector with PARAM_DE entries by a polynomial
* 
* Arguments:   - poly_p_vec_d res: polynomial vector to host the multiplication (initialized)
* 			   		 - const poly_p_vec_d arg: first polynomial vector factor
* 			   		 - const poly_p arg2: second polynomial factor
**************************************************/
void poly_p_vec_d_mul_poly(poly_p_vec_d res, const poly_p_vec_d arg, const poly_p arg2) {
  for (size_t i = 0; i < PARAM_DE; i++) {
    poly_p_mul(res->entries[i], arg->entries[i], arg2);
  }
}

/*************************************************
* Name:        poly_p_vec_d_mul_inner
*
* Description: Inner product of two polynomial vectors with PARAM_DE entries
* 
* Arguments:   - poly_p res: polynomial to host the inner product (initialized)
* 			   		 - const poly_p_vec_d lhs: first polynomial vector
* 			   		 - const poly_p_vec_d rhs: second polynomial vector
**************************************************/
void poly_p_vec_d_mul_inner(poly_p res, const poly_p_vec_d lhs, const poly_p_vec_d rhs) {
	poly_p_zero(res);
	for (size_t i = 0; i < PARAM_DE; ++i) {
		poly_p_mul(TMP, lhs->entries[i], rhs->entries[i]);
		poly_p_add(res, res, TMP);
	}
}

/*************************************************
* Name:        poly_p_vec_d_equal
*
* Description: Equality test between two polynomial vectors with PARAM_DE entries
* 
* Arguments:   - const poly_p_vec_d lhs: first polynomial vector
* 			   		 - const poly_p_vec_d rhs: second polynomial vector
* 
* Returns 1 if the polynomial vectors are equal, 0 otherwise
**************************************************/
int poly_p_vec_d_equal(const poly_p_vec_d lhs, const poly_p_vec_d rhs) {
  for (size_t i = 0; i < PARAM_DE; ++i) {
    if (!poly_p_equal(lhs->entries[i], rhs->entries[i])) {
        return 0;
    }
  }
  return 1;
}

/*************************************************
* Name:        poly_p_vec_d_dump
*
* Description: Print a polynomial vector with PARAM_DE entries
* 
* Arguments:   - const poly_p_vec_d arg: polynomial vector to be printed
**************************************************/
void poly_p_vec_d_dump(const poly_p_vec_d arg) {
	printf("[");
	for (size_t i = 0; i < PARAM_DE - 1; ++i) {
		poly_p_dump(arg->entries[i]);
		printf(", ");
	}
	poly_p_dump(arg->entries[PARAM_DE - 1]);
	printf("]");
}
