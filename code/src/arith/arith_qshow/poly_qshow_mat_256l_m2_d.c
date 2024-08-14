#include "arith_qshow.h"

static poly_qshow TMP;

/*************************************************
* Name:        poly_qshow_mat_256l_m2_d_setup
*
* Description: Initialize and setup the backend for arithmetic 
*              mod PARAM_Q_SHOW integer matrices with PARAM_ARP_DIV_N_L_SHOW x PARAM_M2_D_SHOW entries. 
*              This is strictly required and must be called once 
*              before any other function from here is used.
**************************************************/
void poly_qshow_mat_256l_m2_d_setup(void) {
  poly_qshow_init(TMP);
}

/*************************************************
* Name:        poly_qshow_mat_256l_m2_d_teardown
*
* Description: Clean up and teardown the backend for arithmetic 
*              mod PARAM_Q_SHOW integer matrices with PARAM_ARP_DIV_N_L_SHOW x PARAM_M2_D_SHOW entries. 
*              This is strictly required and must be called once 
*              at the very end to release any resources.
**************************************************/
void poly_qshow_mat_256l_m2_d_teardown(void) {
  poly_qshow_clear(TMP);
}

/*************************************************
* Name:        poly_qshow_mat_256l_m2_d_init
*
* Description: Initialize polynomial matrix with PARAM_ARP_DIV_N_L_SHOW x PARAM_M2_D_SHOW entries.
*              This is strictly required before any operations 
*              are done with/on the matrix.
* 
* Arguments:   - poly_qshow_mat_256l_m2_d res: polynomial matrix to be initialized
**************************************************/
void poly_qshow_mat_256l_m2_d_init(poly_qshow_mat_256l_m2_d res) {
  for (size_t i = 0; i < PARAM_ARP_DIV_N_L_SHOW; ++i) {
    poly_qshow_vec_m2_d_init(res->rows[i]);
  }
}

/*************************************************
* Name:        poly_qshow_mat_256l_m2_d_clear
*
* Description: Clear polynomial matrix with PARAM_ARP_DIV_N_L_SHOW x PARAM_M2_D_SHOW entries.
*              This is strictly required to avoid memory leaks and the 
*              polynomial matrix must not be used again (unless reinitialized).
* 
* Arguments:   - poly_qshow_mat_256l_m2_d res: polynomial matrix to be cleared
**************************************************/
void poly_qshow_mat_256l_m2_d_clear(poly_qshow_mat_256l_m2_d res) {
  for (size_t i = 0; i < PARAM_ARP_DIV_N_L_SHOW; ++i) {
    poly_qshow_vec_m2_d_clear(res->rows[i]);
  }
}

/*************************************************
* Name:        poly_qshow_mat_256l_m2_d_zero
*
* Description: Set an initialized polynomial matrix with PARAM_ARP_DIV_N_L_SHOW x PARAM_M2_D_SHOW entries to zero
* 
* Arguments:   - poly_qshow_mat_256l_m2_d res: polynomial matrix to be zeroized (initialized)
**************************************************/
void poly_qshow_mat_256l_m2_d_zero(poly_qshow_mat_256l_m2_d res) {
  for (size_t i = 0; i < PARAM_ARP_DIV_N_L_SHOW; ++i) {
    poly_qshow_vec_m2_d_zero(res->rows[i]);
  }
}

/*************************************************
* Name:        poly_qshow_mat_256l_m2_d_set
*
* Description: Set a polynomial matrix with PARAM_ARP_DIV_N_L_SHOW x PARAM_M2_D_SHOW entries equal to another polynomial matrix
* 
* Arguments:   - poly_qshow_mat_256l_m2_d res: polynomial matrix to be set (initialized)
*              - const poly_qshow_mat_256l_m2_d arg: polynomial matrix to be read
**************************************************/
void poly_qshow_mat_256l_m2_d_set(poly_qshow_mat_256l_m2_d res, const poly_qshow_mat_256l_m2_d arg) {
  for (size_t i = 0; i < PARAM_ARP_DIV_N_L_SHOW; ++i) {
    poly_qshow_vec_m2_d_set(res->rows[i], arg->rows[i]);
  }
}

/*************************************************
* Name:        poly_qshow_mat_256l_m2_d_neg
*
* Description: Negate a polynomial matrix with PARAM_ARP_DIV_N_L_SHOW x PARAM_M2_D_SHOW entries
* 
* Arguments:   - poly_qshow_mat_256l_m2_d res: polynomial matrix to host the negation (initialized)
*              - const poly_qshow_mat_256l_m2_d arg: polynomial matrix to be negated
**************************************************/
void poly_qshow_mat_256l_m2_d_neg(poly_qshow_mat_256l_m2_d res, const poly_qshow_mat_256l_m2_d arg) {
  for (size_t i = 0; i < PARAM_ARP_DIV_N_L_SHOW; ++i) {
    poly_qshow_vec_m2_d_neg(res->rows[i], arg->rows[i]);
  }
}

/*************************************************
* Name:        poly_qshow_mat_256l_m2_d_add
*
* Description: Add two polynomial matrices with PARAM_ARP_DIV_N_L_SHOW x PARAM_M2_D_SHOW entries
* 
* Arguments:   - poly_qshow_mat_256l_m2_d res: polynomial matrix to host the sum (initialized)
*              - const poly_qshow_mat_256l_m2_d lhs: first polynomial matrix summand
*              - const poly_qshow_mat_256l_m2_d rhs: second polynomial matrix summand
**************************************************/
void poly_qshow_mat_256l_m2_d_add(poly_qshow_mat_256l_m2_d res, const poly_qshow_mat_256l_m2_d lhs, const poly_qshow_mat_256l_m2_d rhs) {
  for (size_t i = 0; i < PARAM_ARP_DIV_N_L_SHOW; ++i) {
    poly_qshow_vec_m2_d_add(res->rows[i], lhs->rows[i], rhs->rows[i]);
  }
}

/*************************************************
* Name:        poly_qshow_mat_256l_m2_d_sub
*
* Description: Substract two polynomial matrices with PARAM_ARP_DIV_N_L_SHOW x PARAM_M2_D_SHOW entries
* 
* Arguments:   - poly_qshow_mat_256l_m2_d res: polynomial matrix to host the difference (initialized)
*              - const poly_qshow_mat_256l_m2_d lhs: first polynomial matrix term
*              - const poly_qshow_mat_256l_m2_d rhs: second polynomial matrix term
**************************************************/
void poly_qshow_mat_256l_m2_d_sub(poly_qshow_mat_256l_m2_d res, const poly_qshow_mat_256l_m2_d lhs, const poly_qshow_mat_256l_m2_d rhs) {
  for (size_t i = 0; i < PARAM_ARP_DIV_N_L_SHOW; ++i) {
    poly_qshow_vec_m2_d_sub(res->rows[i], lhs->rows[i], rhs->rows[i]);
  }
}

/*************************************************
* Name:        poly_qshow_mat_256l_m2_d_mul_vec_m2_d
*
* Description: Product of a polynomial matrix with PARAM_ARP_DIV_N_L_SHOW x PARAM_M2_D_SHOW entries
*              with a polynomial vector with PARAM_M2_D_SHOW entries
* 
* Arguments:   - poly_qshow_vec_256_l res: polynomial vector to host the multiplication (initialized)
*              - const poly_qshow_mat_256l_m2_d lhs: polynomial matrix to multiply
*              - const poly_qshow_vec_m2_d rhs: polynomial vector to multiply
**************************************************/
void poly_qshow_mat_256l_m2_d_mul_vec_m2_d(poly_qshow_vec_256_l res, const poly_qshow_mat_256l_m2_d lhs, const poly_qshow_vec_m2_d rhs) {
  for (size_t i = 0; i < PARAM_ARP_DIV_N_L_SHOW; ++i) {
    poly_qshow_vec_m2_d_mul_inner(TMP, lhs->rows[i], rhs);
    poly_qshow_vec_256_l_set_poly(res, TMP, i);
  }
}

/*************************************************
* Name:        poly_qshow_mat_256l_m2_d_equal
*
* Description: Equality test between two polynomial matrices with PARAM_ARP_DIV_N_L_SHOW x PARAM_M2_D_SHOW entries
* 
* Arguments:   - const poly_qshow_mat_256l_m2_d lhs: first polynomial matrix
*              - const poly_qshow_mat_256l_m2_d rhs: second polynomial matrix
* 
* Returns 1 if the polynomial matrices are equal, 0 otherwise
**************************************************/
int poly_qshow_mat_256l_m2_d_equal(const poly_qshow_mat_256l_m2_d lhs, const poly_qshow_mat_256l_m2_d rhs) {
  for (size_t i = 0; i < PARAM_ARP_DIV_N_L_SHOW; ++i) {
    if (!poly_qshow_vec_m2_d_equal(lhs->rows[i], rhs->rows[i])) {
      return 0;
    }
  }
  return 1;
}

/*************************************************
* Name:        poly_qshow_mat_256l_m2_d_dump
*
* Description: Print a polynomial matrix with PARAM_ARP_DIV_N_L_SHOW x PARAM_M2_D_SHOW entries
* 
* Arguments:   - const poly_qshow_mat_256l_m2_d arg: polynomial matrix to be printed
**************************************************/
void poly_qshow_mat_256l_m2_d_dump(const poly_qshow_mat_256l_m2_d arg) {
	printf("[");
	for (size_t i = 0; i < PARAM_ARP_DIV_N_L_SHOW - 1; ++i) {
		poly_qshow_vec_m2_d_dump(arg->rows[i]);
		printf(", ");
	}
	poly_qshow_vec_m2_d_dump(arg->rows[PARAM_ARP_DIV_N_L_SHOW - 1]);
	printf("]");
}
