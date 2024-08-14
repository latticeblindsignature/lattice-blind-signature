#include "arith_p.h"
#include "macros.h"

static poly_p TMP;

/*************************************************
* Name:        poly_p_vec_m_setup
*
* Description: Initialize and setup the backend for arithmetic 
*              mod PARAM_Q integer vectors with PARAM_ME entries. 
* 			   		 This is strictly required and must be called once 
* 			   		 before any other function from here is used.
**************************************************/
void poly_p_vec_m_setup(void) {
  poly_p_init(TMP);
}

/*************************************************
* Name:        poly_p_vec_m_teardown
*
* Description: Clean up and teardown the backend for arithmetic 
*              mod PARAM_Q integer vectors with PARAM_ME entries. 
* 			   		 This is strictly required and must be called once 
* 			   		 at the very end to release any resources.
**************************************************/
void poly_p_vec_m_teardown(void) {
  poly_p_clear(TMP);
}

/*************************************************
* Name:        poly_p_vec_m_init
*
* Description: Initialize polynomial vector with PARAM_ME entries.
*              This is strictly required before any operations 
*              are done with/on the vector.
* 
* Arguments:   - poly_p_vec_m arg: polynomial vector to be initialized
**************************************************/
void poly_p_vec_m_init(poly_p_vec_m arg) {
	for (size_t i = 0; i < PARAM_ME; ++i) {
		poly_p_init(arg->entries[i]);
	}
}

/*************************************************
* Name:        poly_p_vec_m_clear
*
* Description: Clear polynomial vector with PARAM_ME entries.
*              This is strictly required to avoid memory leaks and the 
*              polynomial vector must not be used again (unless reinitialized).
* 
* Arguments:   - poly_p_vec_m arg: polynomial vector to be cleared
**************************************************/
void poly_p_vec_m_clear(poly_p_vec_m arg) {
	for (size_t i = 0; i < PARAM_ME; ++i) {
		poly_p_clear(arg->entries[i]);
	}
}

/*************************************************
* Name:        poly_p_vec_m_zero
*
* Description: Set an initialized polynomial vector with PARAM_ME entries to zero
* 
* Arguments:   - poly_p_vec_m arg: polynomial vector to be zeroized (initialized)
**************************************************/
void poly_p_vec_m_zero(poly_p_vec_m arg) {
	for (size_t i = 0; i < PARAM_ME; ++i) {
		poly_p_zero(arg->entries[i]);
	}
}

/*************************************************
* Name:        poly_p_vec_m_set
*
* Description: Set a polynomial vector with PARAM_ME entries equal to another polynomial vector
* 
* Arguments:   - poly_p_vec_m res: polynomial vector to be set (initialized)
* 			   		 - const poly_p_vec_m arg: polynomial vector to be read
**************************************************/
void poly_p_vec_m_set(poly_p_vec_m res, const poly_p_vec_m arg) {
	for (size_t i = 0; i < PARAM_ME; ++i) {
		poly_p_set(res->entries[i], arg->entries[i]);
	}
}

/*************************************************
* Name:        poly_p_vec_m_get_poly
*
* Description: Get pos-th entry polynomial of the vector
*              condition: [0 <= pos < PARAM_ME]
* 
* Arguments:   - poly_p res: polynomial to host the pos-th entry (initialized)
* 			   		 - const poly_p_vec_m arg: polynomial vector to be read
* 			   		 - size_t pos: position to get in the vector
**************************************************/
void poly_p_vec_m_get_poly(poly_p res, const poly_p_vec_m arg, size_t pos) {
  ASSERT_DEBUG(pos < PARAM_ME, "Illegal argument: cannot get entry of vector at given position.");
	poly_p_set(res, arg->entries[pos]);
}

/*************************************************
* Name:        poly_p_vec_m_set_poly
*
* Description: Get pos-th entry polynomial of the vector
*              condition: [0 <= pos < PARAM_ME]
* 
* Arguments:   - poly_p_vec_m res: polynomial vector to be set (initialized)
* 			   		 - const poly_p arg: polynomial to set pos-th entry
* 			   		 - size_t pos: position to get in the vector
**************************************************/
void poly_p_vec_m_set_poly(poly_p_vec_m res, const poly_p arg, size_t pos) {
  ASSERT_DEBUG(pos < PARAM_ME, "Illegal argument: cannot set entry of vector at given position.");
	poly_p_set(res->entries[pos], arg);
}

/*************************************************
* Name:        poly_p_vec_m_neg
*
* Description: Negate a polynomial vector with PARAM_ME entries
* 
* Arguments:   - poly_p_vec_m res: polynomial vector to host the negation (initialized)
* 			   		 - const poly_p_vec_m arg: polynomial vector to be negated
**************************************************/
void poly_p_vec_m_neg(poly_p_vec_m res, const poly_p_vec_m arg) {
  for (size_t i = 0; i < PARAM_ME; ++i) {
    poly_p_neg(res->entries[i], arg->entries[i]);
  }
}

/*************************************************
* Name:        poly_p_vec_m_add
*
* Description: Add two polynomial vectors with PARAM_ME entries
* 
* Arguments:   - poly_p_vec_m res: polynomial vector to host the sum (initialized)
* 			   		 - const poly_p_vec_m lhs: first polynomial vector summand
* 			   		 - const poly_p_vec_m rhs: second polynomial vector summand
**************************************************/
void poly_p_vec_m_add(poly_p_vec_m res, const poly_p_vec_m lhs, const poly_p_vec_m rhs) {
	for (size_t i = 0; i < PARAM_ME; ++i) {
		poly_p_add(res->entries[i], lhs->entries[i], rhs->entries[i]);
	}
}

/*************************************************
* Name:        poly_p_vec_m_sub
*
* Description: Substract two polynomial vectors with PARAM_ME entries
* 
* Arguments:   - poly_p_vec_m res: polynomial vector to host the difference (initialized)
* 			   		 - const poly_p_vec_m lhs: first polynomial vector term
* 			   		 - const poly_p_vec_m rhs: second polynomial vector term
**************************************************/
void poly_p_vec_m_sub(poly_p_vec_m res, const poly_p_vec_m lhs, const poly_p_vec_m rhs) {
	for (size_t i = 0; i < PARAM_ME; ++i) {
		poly_p_sub(res->entries[i], lhs->entries[i], rhs->entries[i]);
	}
}

/*************************************************
* Name:        poly_p_vec_m_mul_scalar
*
* Description: Multiplication of a polynomial vector with PARAM_ME entries by a integer scalar
* 
* Arguments:   - poly_p_vec_m res: polynomial vector to host the multiplication (initialized)
* 			   		 - const poly_p_vec_m arg: polynomial vector factor
* 			   		 - coeff_q fac: integer factor
**************************************************/
void poly_p_vec_m_mul_scalar(poly_p_vec_m res, const poly_p_vec_m arg, const poly_p fac) {
	for (size_t i = 0; i < PARAM_ME; ++i) {
		poly_p_mul(res->entries[i], arg->entries[i], fac);
	}
}

/*************************************************
* Name:        poly_p_vec_m_mul_inner
*
* Description: Inner product of two polynomial vectors with PARAM_ME entries
* 
* Arguments:   - poly_p res: polynomial to host the inner product (initialized)
* 			   		 - const poly_p_vec_m lhs: first polynomial vector
* 			   		 - const poly_p_vec_m rhs: second polynomial vector
**************************************************/
void poly_p_vec_m_mul_inner(poly_p res, const poly_p_vec_m lhs, const poly_p_vec_m rhs) {
	poly_p_zero(res);
	for (size_t i = 0; i < PARAM_ME; ++i) {
		poly_p_mul(TMP, lhs->entries[i], rhs->entries[i]);
		poly_p_add(res, res, TMP);
	}
}

/*************************************************
* Name:        poly_p_vec_m_equal
*
* Description: Equality test between two polynomial vectors with PARAM_ME entries
* 
* Arguments:   - const poly_p_vec_m lhs: first polynomial vector
* 			   		 - const poly_p_vec_m rhs: second polynomial vector
* 
* Returns 1 if the polynomial vectors are equal, 0 otherwise
**************************************************/
int poly_p_vec_m_equal(const poly_p_vec_m lhs, const poly_p_vec_m rhs) {
  for (size_t i = 0; i < PARAM_ME; ++i) {
    if (!poly_p_equal(lhs->entries[i], rhs->entries[i])) {
      return 0;
    }
  }
  return 1;
}

/*************************************************
* Name:        poly_p_vec_m_norm2
*
* Description: Compute the square l2 norm of a polynomial vector
* 
* Arguments:   - const poly_p_vec_m arg: the polynomial vector
* 
* Returns an unsigned 64-bit integer with the square l2 norm
**************************************************/
uint64_t poly_p_vec_m_norm2(const poly_p_vec_m arg) {
  uint64_t sq_norm2, tmp;
  size_t i;
  sq_norm2 = 0;
  for (i = 0; i < PARAM_ME; i++) {
  	tmp = poly_p_sq_norm2(arg->entries[i]);
		CHK_UI_OVF_ADDITION(sq_norm2, tmp);
  }
  return sq_norm2;
}

/*************************************************
* Name:        poly_p_vec_m_dump
*
* Description: Print a polynomial vector with PARAM_ME entries
* 
* Arguments:   - const poly_p_vec_m arg: polynomial vector to be printed
**************************************************/
void poly_p_vec_m_dump(const poly_p_vec_m arg) {
	printf("[");
	for (size_t i = 0; i < PARAM_ME - 1; ++i) {
		poly_p_dump(arg->entries[i]);
		printf(", ");
	}
	poly_p_dump(arg->entries[PARAM_ME - 1]);
	printf("]");
}
