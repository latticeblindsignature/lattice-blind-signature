#include "arith_qshow.h"
#include "macros.h"

static poly_qshow TMP;

/*************************************************
* Name:        poly_qshow_vec_256_setup
*
* Description: Initialize and setup the backend for arithmetic 
*              mod PARAM_Q_SHOW integer vectors with PARAM_ARP_DIV_N_SHOW entries. 
* 			   		 This is strictly required and must be called once 
* 			   		 before any other function from here is used.
**************************************************/
void poly_qshow_vec_256_setup(void) {
  poly_qshow_init(TMP);
}

/*************************************************
* Name:        poly_qshow_vec_256_teardown
*
* Description: Clean up and teardown the backend for arithmetic 
*              mod PARAM_Q_SHOW integer vectors with PARAM_ARP_DIV_N_SHOW entries. 
* 			   		 This is strictly required and must be called once 
* 			   		 at the very end to release any resources.
**************************************************/
void poly_qshow_vec_256_teardown(void) {
  poly_qshow_clear(TMP);
}

/*************************************************
* Name:        poly_qshow_vec_256_init
*
* Description: Initialize polynomial vector with PARAM_ARP_DIV_N_SHOW entries.
*              This is strictly required before any operations 
*              are done with/on the vector.
* 
* Arguments:   - poly_qshow_vec_256 arg: polynomial vector to be initialized
**************************************************/
void poly_qshow_vec_256_init(poly_qshow_vec_256 arg) {
	for (size_t i = 0; i < PARAM_ARP_DIV_N_SHOW; ++i) {
		poly_qshow_init(arg->entries[i]);
	}
}

/*************************************************
* Name:        poly_qshow_vec_256_clear
*
* Description: Clear polynomial vector with PARAM_ARP_DIV_N_SHOW entries.
*              This is strictly required to avoid memory leaks and the 
*              polynomial vector must not be used again (unless reinitialized).
* 
* Arguments:   - poly_qshow_vec_256 arg: polynomial vector to be cleared
**************************************************/
void poly_qshow_vec_256_clear(poly_qshow_vec_256 arg) {
	for (size_t i = 0; i < PARAM_ARP_DIV_N_SHOW; ++i) {
		poly_qshow_clear(arg->entries[i]);
	}
}

/*************************************************
* Name:        poly_qshow_vec_256_zero
*
* Description: Set an initialized polynomial vector with PARAM_ARP_DIV_N_SHOW entries to zero
* 
* Arguments:   - poly_qshow_vec_256 arg: polynomial vector to be zeroized (initialized)
**************************************************/
void poly_qshow_vec_256_zero(poly_qshow_vec_256 arg) {
	for (size_t i = 0; i < PARAM_ARP_DIV_N_SHOW; ++i) {
		poly_qshow_zero(arg->entries[i]);
	}
}

/*************************************************
* Name:        poly_qshow_vec_256_set
*
* Description: Set a polynomial vector with PARAM_ARP_DIV_N_SHOW entries equal to another polynomial vector
* 
* Arguments:   - poly_qshow_vec_256 res: polynomial vector to be set (initialized)
* 			   		 - const poly_qshow_vec_256 arg: polynomial vector to be read
**************************************************/
void poly_qshow_vec_256_set(poly_qshow_vec_256 res, const poly_qshow_vec_256 arg) {
	for (size_t i = 0; i < PARAM_ARP_DIV_N_SHOW; ++i) {
		poly_qshow_set(res->entries[i], arg->entries[i]);
	}
}

/*************************************************
* Name:        poly_qshow_vec_256_get_poly
*
* Description: Get pos-th entry polynomial of the vector
*              condition: [0 <= pos < PARAM_ARP_DIV_N_SHOW]
* 
* Arguments:   - poly_qshow res: polynomial to host the pos-th entry (initialized)
* 			   		 - const poly_qshow_vec_256 arg: polynomial vector to be read
* 			   		 - size_t pos: position to get in the vector
**************************************************/
void poly_qshow_vec_256_get_poly(poly_qshow res, const poly_qshow_vec_256 arg, size_t pos) {
	ASSERT_DEBUG(pos < PARAM_ARP_DIV_N_SHOW, "Illegal argument: cannot get entry of vector at given position.");
	poly_qshow_set(res, arg->entries[pos]);
}

/*************************************************
* Name:        poly_qshow_vec_256_set_poly
*
* Description: Get pos-th entry polynomial of the vector
*              condition: [0 <= pos < PARAM_ARP_DIV_N_SHOW]
* 
* Arguments:   - poly_qshow_vec_256 res: polynomial vector to be set (initialized)
* 			   		 - const poly_qshow arg: polynomial to set pos-th entry
* 			   		 - size_t pos: position to get in the vector
**************************************************/
void poly_qshow_vec_256_set_poly(poly_qshow_vec_256 res, const poly_qshow arg, size_t pos) {
	ASSERT_DEBUG(pos < PARAM_ARP_DIV_N_SHOW, "Illegal argument: cannot set entry of vector at given position.");
	poly_qshow_set(res->entries[pos], arg);
}

/*************************************************
* Name:        poly_qshow_vec_256_neg
*
* Description: Negate a polynomial vector with PARAM_ARP_DIV_N_SHOW entries
* 
* Arguments:   - poly_qshow_vec_256 res: polynomial vector to host the negation (initialized)
* 			   		 - const poly_qshow_vec_256 arg: polynomial vector to be negated
**************************************************/
void poly_qshow_vec_256_neg(poly_qshow_vec_256 res, const poly_qshow_vec_256 arg) {
  for (size_t i = 0; i < PARAM_ARP_DIV_N_SHOW; ++i) {
    poly_qshow_neg(res->entries[i], arg->entries[i]);
  }
}

/*************************************************
* Name:        poly_qshow_vec_256_add
*
* Description: Add two polynomial vectors with PARAM_ARP_DIV_N_SHOW entries
* 
* Arguments:   - poly_qshow_vec_256 res: polynomial vector to host the sum (initialized)
* 			   		 - const poly_qshow_vec_256 lhs: first polynomial vector summand
* 			   		 - const poly_qshow_vec_256 rhs: second polynomial vector summand
**************************************************/
void poly_qshow_vec_256_add(poly_qshow_vec_256 res, const poly_qshow_vec_256 lhs, const poly_qshow_vec_256 rhs) {
	for (size_t i = 0; i < PARAM_ARP_DIV_N_SHOW; ++i) {
		poly_qshow_add(res->entries[i], lhs->entries[i], rhs->entries[i]);
	}
}

/*************************************************
* Name:        poly_qshow_vec_256_sub
*
* Description: Substract two polynomial vectors with PARAM_ARP_DIV_N_SHOW entries
* 
* Arguments:   - poly_qshow_vec_256 res: polynomial vector to host the difference (initialized)
* 			   		 - const poly_qshow_vec_256 lhs: first polynomial vector term
* 			   		 - const poly_qshow_vec_256 rhs: second polynomial vector term
**************************************************/
void poly_qshow_vec_256_sub(poly_qshow_vec_256 res, const poly_qshow_vec_256 lhs, const poly_qshow_vec_256 rhs) {
	for (size_t i = 0; i < PARAM_ARP_DIV_N_SHOW; ++i) {
		poly_qshow_sub(res->entries[i], lhs->entries[i], rhs->entries[i]);
	}
}

/*************************************************
* Name:        poly_qshow_vec_256_mul_scalar
*
* Description: Multiplication of a polynomial vector with PARAM_ARP_DIV_N_SHOW entries by a integer scalar
* 
* Arguments:   - poly_qshow_vec_256 res: polynomial vector to host the multiplication (initialized)
* 			   		 - const poly_qshow_vec_256 arg: polynomial vector factor
* 			   		 - coeff_qshow fac: integer factor
**************************************************/
void poly_qshow_vec_256_mul_scalar(poly_qshow_vec_256 res, const poly_qshow_vec_256 arg, const coeff_qshow fac) {
	for (size_t i = 0; i < PARAM_ARP_DIV_N_SHOW; ++i) {
		poly_qshow_mul_scalar(res->entries[i], arg->entries[i], fac);
	}
}

/*************************************************
* Name:        poly_qshow_vec_256_mul_inner
*
* Description: Inner product of two polynomial vectors with PARAM_ARP_DIV_N_SHOW entries
* 
* Arguments:   - poly_qshow res: polynomial to host the inner product (initialized)
* 			   		 - const poly_qshow_vec_256 lhs: first polynomial vector
* 			   		 - const poly_qshow_vec_256 rhs: second polynomial vector
**************************************************/
void poly_qshow_vec_256_mul_inner(poly_qshow res, const poly_qshow_vec_256 lhs, const poly_qshow_vec_256 rhs) {
	poly_qshow_zero(res);
	for (size_t i = 0; i < PARAM_ARP_DIV_N_SHOW; ++i) {
		poly_qshow_mul(TMP, lhs->entries[i], rhs->entries[i]);
		poly_qshow_add(res, res, TMP);
	}
}

/*************************************************
* Name:        poly_qshow_vec_256_equal
*
* Description: Equality test between two polynomial vectors with PARAM_ARP_DIV_N_SHOW entries
* 
* Arguments:   - const poly_qshow_vec_256 lhs: first polynomial vector
* 			   		 - const poly_qshow_vec_256 rhs: second polynomial vector
* 
* Returns 1 if the polynomial vectors are equal, 0 otherwise
**************************************************/
int poly_qshow_vec_256_equal(const poly_qshow_vec_256 lhs, const poly_qshow_vec_256 rhs) {
  for (size_t i = 0; i < PARAM_ARP_DIV_N_SHOW; ++i) {
    if (!poly_qshow_equal(lhs->entries[i], rhs->entries[i])) {
      return 0;
    }
  }
  return 1;
}

/*************************************************
* Name:        poly_qshow_vec_256_dump
*
* Description: Print a polynomial vector with PARAM_ARP_DIV_N_SHOW entries
* 
* Arguments:   - const poly_qshow_vec_256 arg: polynomial vector to be printed
**************************************************/
void poly_qshow_vec_256_dump(const poly_qshow_vec_256 arg) {
	printf("[");
	for (size_t i = 0; i < PARAM_ARP_DIV_N_SHOW - 1; ++i) {
		poly_qshow_dump(arg->entries[i]);
		printf(", ");
	}
	poly_qshow_dump(arg->entries[PARAM_ARP_DIV_N_SHOW - 1]);
	printf("]");
}
