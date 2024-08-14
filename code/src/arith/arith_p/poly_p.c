#include "arith_p.h"
#include "macros.h"

static nmod_poly_t POLY_F;
static nmod_poly_t POLY_F_REV_INV;

static poly_p TMP;
static poly_p AUX;

/*************************************************
* Name:        nmod_poly_invert
*
* Description: Invert an nmod_poly to later compute faster multiplication
*              using the precomputed inverse of reverse(x^PARAM_N + 1)
*
* Arguments:   - nmod_poly_t res: polynomial to host the inverse
*              - const nmod_poly_t arg: polynomial to invert
*              - const nmod_poly_t mod: polynomial to reduce with
*              - int64_t p: modulus
**************************************************/
static void nmod_poly_invert(nmod_poly_t res, const nmod_poly_t arg, const nmod_poly_t mod, int64_t p) {
  nmod_poly_t G, S;
  nmod_poly_init(G, p);
  nmod_poly_init(S, p);
  nmod_poly_xgcd(G, S, res, mod, arg);
  ASSERT_DEBUG(nmod_poly_is_one(G), "GDC != 1");
  nmod_poly_mulmod(G, res, arg, mod);
  ASSERT_DEBUG(nmod_poly_is_one(G), "arg * res % mod != 1");
  nmod_poly_clear(G);
  nmod_poly_clear(S);
}

/*************************************************
* Name:        arith_q_setup
*
* Description: Initialize and setup the backend for arithmetic
*              modulo PARAM_P. This is strictly required
*              and must be called once before any other function
*              from here is used.
**************************************************/
void arith_p_setup(void) {
  ASSERT_DEBUG(PARAM_N % 8 == 0, "Illegal parameter: ring degree must be divisible by 8 to allow `poly_from_bits` function.");
  nmod_poly_init(POLY_F, PARAM_P);
  nmod_poly_set_coeff_ui(POLY_F, 0, 1);
  nmod_poly_set_coeff_ui(POLY_F, PARAM_N, 1);

  nmod_poly_t f_len;
  nmod_poly_init(f_len, PARAM_P);
  nmod_poly_set_coeff_ui(f_len, PARAM_N + 1, 1);

  nmod_poly_t f_rev;
  nmod_poly_init(f_rev, PARAM_P);

  nmod_poly_init(POLY_F_REV_INV, PARAM_P);
  nmod_poly_reverse(f_rev, POLY_F, PARAM_N + 1);
  nmod_poly_invert(POLY_F_REV_INV, f_rev, f_len, PARAM_P);

  nmod_poly_clear(f_len);
  nmod_poly_clear(f_rev);

  poly_p_init(TMP);
  poly_p_init(AUX);
}

/*************************************************
* Name:        arith_q_teardown
*
* Description: Clean up and teardown the backend for arithmetic
*              modulo PARAM_P. This is strictly required
*              and must be called once at the very end to release
*              any resources.
**************************************************/
void arith_p_teardown(void) {
  poly_p_clear(TMP);
  poly_p_clear(AUX);

  nmod_poly_clear(POLY_F);
  nmod_poly_clear(POLY_F_REV_INV);
}

/*************************************************
* Name:        poly_p_init
*
* Description: Initialize polynomial and set it to zero
*              This is strictly required before any operations
*              are done with/on the polynomial.
*
* Arguments:   - poly_p res: polynomial to be initialized
**************************************************/
void poly_p_init(poly_p res) {
	nmod_poly_init(res, PARAM_P);
}

/*************************************************
* Name:        poly_p_clear
*
* Description: Clears a polynomial and releases all associated memory.
*              This is strictly required to avoid memory leaks and the
*              polynomial must not be used again (unless reinitialized).
*
* Arguments:   - poly_p arg: polynomial to be cleared
**************************************************/
void poly_p_clear(poly_p arg) {
	nmod_poly_clear(arg);
}

/*************************************************
* Name:        poly_p_zero
*
* Description: Set an initialized polynomial to zero
*
* Arguments:   - poly_p res: polynomial to be zeroized (initialized)
**************************************************/
void poly_p_zero(poly_p res) {
	nmod_poly_zero(res);
}

/*************************************************
* Name:        poly_p_set
*
* Description: Set a polynomial equal to another polynomial.
*              Coefficients are reduced mod PARAM_P
*
* Arguments:   - poly_p res: polynomial to be set (initialized)
*              - const poly_p arg: polynomial to be read
**************************************************/
void poly_p_set(poly_p res, const poly_p arg) {
	nmod_poly_set(res, arg);
}

/*************************************************
* Name:        poly_p_get_coeff
*
* Description: Get coefficient of x^n of a polynomial
*              condition: [0 <= n < PARAM_N]
*
* Arguments:   - const poly_p arg: polynomial to be read
*              - size_t n: degree of the coefficient to be read
*
* Returns the coefficients of x^n of arg
**************************************************/
coeff_p poly_p_get_coeff(const poly_p arg, size_t n) {
  ASSERT_DEBUG(n < PARAM_N, "Illegal argument: cannot get coefficient of poly at given position.");
  return nmod_poly_get_coeff_ui(arg, n);
}

/*************************************************
* Name:        poly_p_get_coeff_centered
*
* Description: Get coefficient of x^n of a polynomial in centered representation
*              condition: [0 <= n < PARAM_N]
*
* Arguments:   - const poly_p arg: polynomial to be read
*              - size_t n: degree of the coefficient to be read
*
* Returns the coefficients of x^n of arg in [-PARAM_P/2, PARAM_P/2]
**************************************************/
coeff_p poly_p_get_coeff_centered(const poly_p arg, size_t n) {
  coeff_p tmp = poly_p_get_coeff(arg, n);
  return tmp - ((~((tmp - PARAM_P/2) >> (sizeof(coeff_p)*8-1))) & PARAM_P);
}

/*************************************************
* Name:        poly_p_set_coeff
*
* Description: Set coefficient of x^n of a polynomial
*              condition: [0 <= n < PARAM_N]
*              Coefficient is reduced mod PARAM_P
*
* Arguments:   - poly_p arg: polynomial whose n-th coefficient is set (initialized)
*              - size_t n: degree of the coefficient to be set
*              - coeff_p c: the new coefficient
**************************************************/
void poly_p_set_coeff(poly_p arg, size_t n, coeff_p c) {
  ASSERT_DEBUG(n < PARAM_N, "Illegal argument: cannot set coefficient of poly at given position.");
  c %= PARAM_P;
  c += (c >> (sizeof(coeff_p)*8-1)) & PARAM_P;
  ASSERT_DEBUG(c >= 0, "Don't give me too much negativity!");
  nmod_poly_set_coeff_ui(arg, n, c);
}

/*************************************************
* Name:        poly_p_from_bits
*
* Description: Set polynomial with {0,1} coefficients from byte-array
*              Example: Given an array { 0xB2, 0x08 } where 0xB2 = 0b1011 0010
*                   and 0x08 = 0b0000 1000 and a polynomial p with precisely 16
*                   coefficients, this functions sets p = x^15 + x^13 + x^12 + x^9 + x^3.
*
* Arguments:   - poly_p arg: polynomial to be set from byte-array (initialized)
*              - const uint8_t *coeffs: byte array for coefficients (allocated PARAM_N/8 bytes)
**************************************************/
void poly_p_from_bits(poly_p arg, const uint8_t coeffs[PARAM_N / 8]) {
  poly_p_zero(arg);
  for (size_t i = 0; i < PARAM_N; ++i) {
    coeff_p c = coeffs[i/8];
    c = (c >> (i % 8)) & 1;
    poly_p_set_coeff(arg, i, c);
  }
}

/*************************************************
* Name:        poly_p_neg
*
* Description: Negate a polynomial coefficient-wise
*              Coefficients are reduced mod PARAM_P
*
* Arguments:   - poly_p res: polynomial to host the negation (initialized)
*              - const poly_p arg: polynomial to be negated
**************************************************/
void poly_p_neg(poly_p res, const poly_p arg) {
  nmod_poly_neg(res, arg);
}

/*************************************************
* Name:        poly_p_add
*
* Description: Add two polynomials. Coefficients are reduced mod PARAM_P
*
* Arguments:   - poly_p res: polynomial to host the sum (initialized)
*              - const poly_p lhs: first polynomial summand
*              - const poly_p rhs: second polynomial summand
**************************************************/
void poly_p_add(poly_p res, const poly_p lhs, const poly_p rhs) {
	nmod_poly_add(res, lhs, rhs);
}

/*************************************************
* Name:        poly_p_sub
*
* Description: Substract two polynomials. Coefficients are reduced mod PARAM_P
*
* Arguments:   - poly_p res: polynomial to host the difference (initialized)
*              - const poly_p lhs: first polynomial term
*              - const poly_p rhs: second polynomial term
**************************************************/
void poly_p_sub(poly_p res, const poly_p lhs, const poly_p rhs) {
  nmod_poly_sub(res, lhs, rhs);
}

/*************************************************
* Name:        poly_p_mul
*
* Description: Multiplication of two polynomials reduced
*              modulo x^PARAM_N + 1.
*              Coefficients are reduced mod PARAM_P
*
* Arguments:   - poly_p res: polynomial to host the multiplication (initialized)
*              - const poly_p lhs: first polynomial factor
*              - const poly_p rhs: second polynomial factor
**************************************************/
void poly_p_mul(poly_p res, const poly_p lhs, const poly_p rhs) {
  ASSERT_DEBUG(nmod_poly_degree(lhs) < PARAM_N, "Argument to `poly_p_mul` must already be reduced for FLINT.");
  ASSERT_DEBUG(nmod_poly_degree(rhs) < PARAM_N, "Argument to `poly_p_mul` must already be reduced for FLINT.");
	nmod_poly_mulmod_preinv(res, lhs, rhs, POLY_F, POLY_F_REV_INV);
}

/*************************************************
* Name:        poly_p_mul_scalar
*
* Description: Multiplication of a polynomials by a scalar
*
* Arguments:   - poly_p res: polynomial to host the multiplication (initialized)
*              - const poly_p arg: first polynomial factor
*              - coeff_p fac: second scalar factor
**************************************************/
void poly_p_mul_scalar(poly_p res, const poly_p arg, coeff_p fac) {
  ASSERT_DEBUG(fac < PARAM_P, "Scalar must not exceed q for scalar polynomial multiplication.");
  nmod_poly_scalar_mul_nmod(res, arg, fac);
}

/*************************************************
* Name:        poly_p_equal
*
* Description: Equality test between two polynomials
*
* Arguments:   - const poly_p lhs: first polynomial
*              - const poly_p rhs: second polynomial
*
* Returns 1 if the polynomials are equal, 0 otherwise
**************************************************/
int poly_p_equal(const poly_p lhs, const poly_p rhs) {
  return nmod_poly_equal(lhs, rhs);
}

/*************************************************
* Name:        poly_p_dump
*
* Description: Print a polynomial
*
* Arguments:   - const poly_p arg: the polynomial to be printed
**************************************************/
void poly_p_dump(const poly_p arg) {
  printf("[");
  for (size_t i = 0; i < PARAM_N; i++) {
    printf("%ld, ", poly_p_get_coeff_centered(arg, i));
  }
  printf("]");
}

/*************************************************
* Name:        poly_p_sq_norm2
*
* Description: Compute the square l2 norm of a polynomial
*
* Arguments:   - const poly_p arg: the polynomial
*
* Returns an unsigned 64-bit integer with the square l2 norm
**************************************************/
uint64_t poly_p_sq_norm2(const poly_p arg) {
	uint64_t sq_norm2 = 0;
	coeff_p c;
	for (size_t i = 0; i < PARAM_N; i++) {
		c = poly_p_get_coeff_centered(arg, i);
		CHK_UI_OVF_ADDITION(sq_norm2, (uint64_t) (c * c));
	}
	return sq_norm2;
}