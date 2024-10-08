#include "arith_p.h"
#include "macros.h"

static poly_p TMP;

/*************************************************
* Name:        poly_p_mat_d_m_setup
*
* Description: Initialize and setup the backend for arithmetic 
*              mod PARAM_Q integer matrices with PARAM_DE x PARAM_ME entries. 
*              This is strictly required and must be called once 
*              before any other function from here is used.
**************************************************/
void poly_p_mat_d_m_setup(void) {
  poly_p_init(TMP);
}

/*************************************************
* Name:        poly_p_mat_d_m_teardown
*
* Description: Clean up and teardown the backend for arithmetic 
*              mod PARAM_Q integer matrices with PARAM_DE x PARAM_ME entries. 
*              This is strictly required and must be called once 
*              at the very end to release any resources.
**************************************************/
void poly_p_mat_d_m_teardown(void) {
  poly_p_clear(TMP);
}

/*************************************************
* Name:        poly_p_mat_d_m_init
*
* Description: Initialize polynomial matrix with PARAM_DE x PARAM_ME entries.
*              This is strictly required before any operations 
*              are done with/on the matrix.
* 
* Arguments:   - poly_p_mat_d_m res: polynomial matrix to be initialized
**************************************************/
void poly_p_mat_d_m_init(poly_p_mat_d_m res) {
  for (size_t i = 0; i < PARAM_DE; ++i) {
    poly_p_vec_m_init(res->rows[i]);
  }
}

/*************************************************
* Name:        poly_p_mat_d_m_clear
*
* Description: Clear polynomial matrix with PARAM_DE x PARAM_ME entries.
*              This is strictly required to avoid memory leaks and the 
*              polynomial matrix must not be used again (unless reinitialized).
* 
* Arguments:   - poly_p_mat_d_m res: polynomial matrix to be cleared
**************************************************/
void poly_p_mat_d_m_clear(poly_p_mat_d_m res) {
  for (size_t i = 0; i < PARAM_DE; ++i) {
    poly_p_vec_m_clear(res->rows[i]);
  }
}

/*************************************************
* Name:        poly_p_mat_d_m_zero
*
* Description: Set an initialized polynomial matrix with PARAM_DE x PARAM_ME entries to zero
* 
* Arguments:   - poly_p_mat_d_m res: polynomial matrix to be zeroized (initialized)
**************************************************/
void poly_p_mat_d_m_zero(poly_p_mat_d_m res) {
  for (size_t i = 0; i < PARAM_DE; ++i) {
    poly_p_vec_m_zero(res->rows[i]);
  }
}

/*************************************************
* Name:        poly_p_mat_d_m_set
*
* Description: Set a polynomial matrix with PARAM_DE x PARAM_ME entries equal to another polynomial matrix
* 
* Arguments:   - poly_p_mat_d_m res: polynomial matrix to be set (initialized)
*              - const poly_p_mat_d_m arg: polynomial matrix to be read
**************************************************/
void poly_p_mat_d_m_set(poly_p_mat_d_m res, const poly_p_mat_d_m arg) {
  for (size_t i = 0; i < PARAM_DE; ++i) {
    poly_p_vec_m_set(res->rows[i], arg->rows[i]);
  }
}

/*************************************************
* Name:        poly_p_mat_d_m_neg
*
* Description: Negate a polynomial matrix with PARAM_DE x PARAM_ME entries
* 
* Arguments:   - poly_p_mat_d_m res: polynomial matrix to host the negation (initialized)
*              - const poly_p_mat_d_m arg: polynomial matrix to be negated
**************************************************/
void poly_p_mat_d_m_neg(poly_p_mat_d_m res, const poly_p_mat_d_m arg) {
  for (size_t i = 0; i < PARAM_DE; ++i) {
    poly_p_vec_m_neg(res->rows[i], arg->rows[i]);
  }
}

/*************************************************
* Name:        poly_p_mat_d_m_add
*
* Description: Add two polynomial matrices with PARAM_DE x PARAM_ME entries
* 
* Arguments:   - poly_p_mat_d_m res: polynomial matrix to host the sum (initialized)
*              - const poly_p_mat_d_m lhs: first polynomial matrix summand
*              - const poly_p_mat_d_m rhs: second polynomial matrix summand
**************************************************/
void poly_p_mat_d_m_add(poly_p_mat_d_m res, const poly_p_mat_d_m lhs, const poly_p_mat_d_m rhs) {
  for (size_t i = 0; i < PARAM_DE; ++i) {
    poly_p_vec_m_add(res->rows[i], lhs->rows[i], rhs->rows[i]);
  }
}

/*************************************************
* Name:        poly_p_mat_d_m_sub
*
* Description: Substract two polynomial matrices with PARAM_DE x PARAM_ME entries
* 
* Arguments:   - poly_p_mat_d_m res: polynomial matrix to host the difference (initialized)
*              - const poly_p_mat_d_m lhs: first polynomial matrix term
*              - const poly_p_mat_d_m rhs: second polynomial matrix term
**************************************************/
void poly_p_mat_d_m_sub(poly_p_mat_d_m res, const poly_p_mat_d_m lhs, const poly_p_mat_d_m rhs) {
  for (size_t i = 0; i < PARAM_DE; ++i) {
    poly_p_vec_m_sub(res->rows[i], lhs->rows[i], rhs->rows[i]);
  }
}

/*************************************************
* Name:        poly_p_mat_d_m_mul_vec_m
*
* Description: Product of a polynomial matrix with PARAM_DE x PARAM_ME entries
*              with a polynomial vector with PARAM_ME entries
* 
* Arguments:   - poly_p_vec_d res: polynomial vector to host the multiplication (initialized)
*              - const poly_p_mat_d_m lhs: polynomial matrix to multiply
*              - const poly_p_vec_k rhs: polynomial vector to multiply
**************************************************/
void poly_p_mat_d_m_mul_vec_m(poly_p_vec_d res, const poly_p_mat_d_m lhs, const poly_p_vec_m rhs) {
  for (size_t i = 0; i < PARAM_DE; ++i) {
    poly_p_vec_m_mul_inner(TMP, lhs->rows[i], rhs);
    poly_p_vec_d_set_poly(res, TMP, i);
  }
}

/*************************************************
* Name:        poly_p_mat_d_m_equal
*
* Description: Equality test between two polynomial matrices with PARAM_DE x PARAM_ME entries
* 
* Arguments:   - const poly_p_mat_d_m lhs: first polynomial matrix
*              - const poly_p_mat_d_m rhs: second polynomial matrix
* 
* Returns 1 if the polynomial matrices are equal, 0 otherwise
**************************************************/
int poly_p_mat_d_m_equal(const poly_p_mat_d_m lhs, const poly_p_mat_d_m rhs) {
  for (size_t i = 0; i < PARAM_DE; ++i) {
    if (!poly_p_vec_m_equal(lhs->rows[i], rhs->rows[i])) {
      return 0;
    }
  }
  return 1;
}

/*************************************************
* Name:        poly_p_mat_d_m_dump
*
* Description: Print a polynomial matrix with PARAM_DE x PARAM_ME entries
* 
* Arguments:   - const poly_p_mat_d_m arg: polynomial matrix to be printed
**************************************************/
void poly_p_mat_d_m_dump(const poly_p_mat_d_m arg) {
	printf("[");
	for (size_t i = 0; i < PARAM_DE - 1; ++i) {
		poly_p_vec_m_dump(arg->rows[i]);
		printf(", ");
	}
	poly_p_vec_m_dump(arg->rows[PARAM_DE - 1]);
	printf("]");
}
