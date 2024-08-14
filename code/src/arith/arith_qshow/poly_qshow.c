#include "arith_qshow.h"
#include "macros.h"

static nmod_poly_t POLY_F;
static nmod_poly_t POLY_F_REV_INV;

/*************************************************
* Name:        ct_abs
*
* Description: Absolute value
*
* Arguments:   - coeff_qshow c: coefficient 
* 
* Returns an 64-bit integer with the absolute value of c
**************************************************/
static int64_t ct_abs(coeff_qshow c) {
  int64_t const mask = c >> (sizeof(coeff_qshow)*8 - 1);
  return (c + mask) ^ mask;
}

/*************************************************
* Name:        nmod_poly_invert
*
* Description: Invert an nmod_poly to later compute faster multiplication
*              using the precomputed inverse of reverse(x^PARAM_N + 1)
* 
* Arguments:   - nmod_poly_t res: polynomial to host the inverse
*              - const nmod_poly_t arg: polynomial to invert
*              - const nmod_poly_t mod: polynomial to reduce with
*              - int64_t q: modulus
**************************************************/
static void nmod_poly_invert(nmod_poly_t res, const nmod_poly_t arg, const nmod_poly_t mod, int64_t q) {
    nmod_poly_t G, S;
    nmod_poly_init(G, q);
    nmod_poly_init(S, q);
    nmod_poly_xgcd(G, S, res, mod, arg);
    ASSERT_DEBUG(nmod_poly_is_one(G), "GDC != 1");
    nmod_poly_mulmod(G, res, arg, mod);
    ASSERT_DEBUG(nmod_poly_is_one(G), "arg * res % mod != 1");
    nmod_poly_clear(G);
    nmod_poly_clear(S);
}

/*************************************************
* Name:        decompose_show
*
* Description: Decompose x into x_L + PARAM_GAMMA_SHOW.x_H with x_L 
*              in (-PARAM_GAMMA_SHOW/2, PARAM_GAMMA_SHOW/2] except if 
*              x_H = (PARAM_Q_SHOW-1)/PARAM_GAMMA_SHOW where we set 
*              x_H = 0 and x_L = x mod+ PARAM_Q_SHOW
*
* Arguments:   - coeff_qshow *low_c: coefficient to host low part of c
*              - const coeff_qshow c: coefficient to be decomposed
* 
* Returns high_c the high part of c
**************************************************/
static coeff_qshow decompose_show(coeff_qshow *low_c, const coeff_qshow c) {
  coeff_qshow high_c;
  high_c = (c + (PARAM_GAMMA_SHOW>>1) - 1) / PARAM_GAMMA_SHOW; // floor division here
  high_c -= (((PARAM_Q_GAMMA_SHOW - high_c - 1) >> (sizeof(coeff_qshow)*8-1)) & high_c); // putting high_c to 0 if high_c = (PARAM_Q_SHOW-1)/PARAM_GAMMA_SHOW
  *low_c = c - PARAM_GAMMA_SHOW * high_c;
  return high_c;
}

/*************************************************
* Name:        highbits_show
*
* Description: Decompose x into x_L + PARAM_GAMMA_SHOW.x_H with x_L 
*              in (-PARAM_GAMMA_SHOW/2, PARAM_GAMMA_SHOW/2] except if 
*              x_H = (PARAM_Q_SHOW-1)/PARAM_GAMMA_SHOW where we set 
*              x_H = 0 and x_L = x mod+ PARAM_Q_SHOW. Return only x_H
*
* Arguments:   - const coeff_qshow c: coefficient to be decomposed
* 
* Returns high_c the high part of c
**************************************************/
static coeff_qshow highbits_show(const coeff_qshow c) {
  coeff_qshow high_c = (c + PARAM_Q_SHOW) % PARAM_Q_SHOW;
  high_c = (high_c + (PARAM_GAMMA_SHOW>>1) - 1) / PARAM_GAMMA_SHOW; // floor division here
  return high_c - (((PARAM_Q_GAMMA_SHOW - high_c - 1) >> (sizeof(coeff_qshow)*8-1)) & high_c); // putting high_c to 0 if high_c = (PARAM_Q_SHOW-1)/PARAM_GAMMA_SHOW;
}

/*************************************************
* Name:        power2round_show
*
* Description: Decompose x into x_L + 2^PARAM_D_ROUND_SHOW.x_H with x_L 
*              in (-2^(PARAM_D_ROUND_SHOW-1), 2^(PARAM_D_ROUND_SHOW-1)]
*
* Arguments:   - coeff_qshow *low_c: coefficient to host low part of c
*              - const coeff_qshow c: coefficient to be decomposed
* 
* Returns high_c the high part of c
**************************************************/
static coeff_qshow power2round_show(coeff_qshow *low_c, const coeff_qshow c) {
  coeff_qshow high_c;
  high_c = (c + (1 << (PARAM_D_ROUND_SHOW-1)) - 1) >> PARAM_D_ROUND_SHOW;
  *low_c = c - (high_c << PARAM_D_ROUND_SHOW); // high_c positive so no runtime error when left shift
  return high_c;
}

/*************************************************
* Name:        makeGhint_show
*
* Description: 
*
* Arguments:   - const coeff_qshow z: 
*              - const coeff_qshow r: 
* 
* Returns h_c the hint coeff on z
**************************************************/
static coeff_qshow makeGhint_show(const coeff_qshow z, const coeff_qshow r) {
  coeff_qshow h_c;
  h_c = (highbits_show(z+r) - highbits_show(r)) % PARAM_Q_GAMMA_SHOW;
  return h_c - ((~((h_c - (PARAM_Q_GAMMA_SHOW>>1) - 1) >> (sizeof(coeff_qshow)*8-1))) & PARAM_Q_GAMMA_SHOW);
}

/*************************************************
* Name:        useGhint_show
*
* Description: 
*
* Arguments:   - const coeff_qshow h_c: 
*              - const coeff_qshow r: 
* 
* Returns high_zr high part of z+r
**************************************************/
static coeff_qshow useGhint_show(const coeff_qshow hint_c, const coeff_qshow r) {
  return (hint_c + highbits_show(r)) % PARAM_Q_GAMMA_SHOW;
}

/*************************************************
* Name:        arith_qshow_setup
*
* Description: Initialize and setup the backend for arithmetic 
*              modulo PARAM_Q_SHOW. This is strictly required 
*              and must be called once before any other function 
*              from here is used.
**************************************************/
void arith_qshow_setup(void) {
  nmod_poly_init(POLY_F, PARAM_Q_SHOW);

	nmod_poly_set_coeff_ui(POLY_F, 0, 1);
  nmod_poly_set_coeff_ui(POLY_F, PARAM_N_SHOW, 1);

	nmod_poly_t f_len;
	nmod_poly_init(f_len, PARAM_Q_SHOW);
	nmod_poly_set_coeff_ui(f_len, PARAM_N_SHOW + 1, 1);

	nmod_poly_t f_rev;
	nmod_poly_init(f_rev, PARAM_Q_SHOW);

	nmod_poly_init(POLY_F_REV_INV, PARAM_Q_SHOW);
	nmod_poly_reverse(f_rev, POLY_F, PARAM_N_SHOW + 1);
	nmod_poly_invert(POLY_F_REV_INV, f_rev, f_len, PARAM_Q_SHOW);

	nmod_poly_clear(f_len);
	nmod_poly_clear(f_rev);
}

/*************************************************
* Name:        arith_qshow_teardown
*
* Description: Clean up and teardown the backend for arithmetic 
*              modulo PARAM_Q_SHOW. This is strictly required 
*              and must be called once at the very end to release 
*              any resources.
**************************************************/
void arith_qshow_teardown(void) {
  nmod_poly_clear(POLY_F);
	nmod_poly_clear(POLY_F_REV_INV);
}

/*************************************************
* Name:        poly_qshow_init
*
* Description: Initialize polynomial and set it to zero
*              This is strictly required before any operations 
*              are done with/on the polynomial.
* 
* Arguments:   - poly_qshow res: polynomial to be initialized
**************************************************/
void poly_qshow_init(poly_qshow res) {
  nmod_poly_init(res, PARAM_Q_SHOW);
}

/*************************************************
* Name:        poly_qshow_clear
*
* Description: Clears a polynomial and releases all associated memory. 
*              This is strictly required to avoid memory leaks and the 
*              polynomial must not be used again (unless reinitialized).
* 
* Arguments:   - poly_qshow arg: polynomial to be cleared
**************************************************/
void poly_qshow_clear(poly_qshow arg) {
  nmod_poly_clear(arg);
}

/*************************************************
* Name:        poly_qshow_zero
*
* Description: Set an initialized polynomial to zero
* 
* Arguments:   - poly_qshow res: polynomial to be zeroized (initialized)
**************************************************/
void poly_qshow_zero(poly_qshow res) {
  nmod_poly_zero(res);
}

/*************************************************
* Name:        poly_qshow_set
*
* Description: Set a polynomial equal to another polynomial.
*              Coefficients are reduced mod PARAM_Q_SHOW
* 
* Arguments:   - poly_qshow res: polynomial to be set (initialized)
*              - const poly_qshow arg: polynomial to be read
**************************************************/
void poly_qshow_set(poly_qshow res, const poly_qshow arg) {
	nmod_poly_set(res, arg);
}

/*************************************************
* Name:        poly_qshow_get_coeff
*
* Description: Get coefficient of x^n of a polynomial
*              condition: [0 <= n < PARAM_N_SHOW]
* 
* Arguments:   - const poly_qshow arg: polynomial to be read
*              - size_t n: degree of the coefficient to be read
* 
* Returns the coefficients of x^n of arg
**************************************************/
coeff_qshow poly_qshow_get_coeff(const poly_qshow arg, size_t n) {
	ASSERT_DEBUG(n < PARAM_N_SHOW, "Illegal argument: cannot get coefficient of poly at given position.");
	return nmod_poly_get_coeff_ui(arg, n);
}

/*************************************************
* Name:        poly_qshow_get_coeff_centered
*
* Description: Get coefficient of x^n of a polynomial in centered representation
*              condition: [0 <= n < PARAM_N]
* 
* Arguments:   - const poly_qshow arg: polynomial to be read
*              - size_t n: degree of the coefficient to be read
* 
* Returns the coefficients of x^n of arg in [-PARAM_Q_SHOW/2, PARAM_Q_SHOW/2]
**************************************************/
coeff_qshow poly_qshow_get_coeff_centered(const poly_qshow arg, size_t n) {
	coeff_qshow tmp = poly_qshow_get_coeff(arg, n);
	return tmp - ((~((tmp - PARAM_Q_SHOW/2) >> (sizeof(coeff_qshow)*8-1))) & PARAM_Q_SHOW);
}

/*************************************************
* Name:        poly_qshow_set_coeff
*
* Description: Set coefficient of x^n of a polynomial
*              condition: [0 <= n < PARAM_N_SHOW]
*              Coefficient is reduced mod PARAM_Q_SHOW
* 
* Arguments:   - poly_qshow arg: polynomial whose n-th coefficient is set (initialized)
*              - size_t n: degree of the coefficient to be set
*              - coeff_qshow c: the new coefficient
**************************************************/
void poly_qshow_set_coeff(poly_qshow arg, size_t n, coeff_qshow c) {
	ASSERT_DEBUG(n < PARAM_N_SHOW, "Illegal argument: cannot set coefficient of poly at given position.");
  c += PARAM_Q_SHOW;
  c %= PARAM_Q_SHOW;
  c += (c >> (sizeof(coeff_qshow)*8 - 1)) & PARAM_Q_SHOW;
	ASSERT_DEBUG(c >= 0, "Coefficients of polys must not be negative.");
	nmod_poly_set_coeff_ui(arg, n, c);
}

/*************************************************
* Name:        poly_qshow_neg
*
* Description: Negate a polynomial coefficient-wise
*              Coefficients are reduced mod PARAM_Q_SHOW
* 
* Arguments:   - poly_qshow res: polynomial to host the negation (initialized)
*              - const poly_qshow arg: polynomial to be negated
**************************************************/
void poly_qshow_neg(poly_qshow res, const poly_qshow arg) {
  nmod_poly_neg(res, arg);
}

/*************************************************
* Name:        poly_qshow_add
*
* Description: Add two polynomials. Coefficients are reduced mod PARAM_Q_SHOW
* 
* Arguments:   - poly_qshow res: polynomial to host the sum (initialized)
*              - const poly_qshow lhs: first polynomial summand
*              - const poly_qshow rhs: second polynomial summand
**************************************************/
void poly_qshow_add(poly_qshow res, const poly_qshow lhs, const poly_qshow rhs) {
	nmod_poly_add(res, lhs, rhs);
}

/*************************************************
* Name:        poly_qshow_sub
*
* Description: Substract two polynomials. Coefficients are reduced mod PARAM_Q_SHOW
* 
* Arguments:   - poly_qshow res: polynomial to host the difference (initialized)
*              - const poly_qshow lhs: first polynomial term
*              - const poly_qshow rhs: second polynomial term
**************************************************/
void poly_qshow_sub(poly_qshow res, const poly_qshow lhs, const poly_qshow rhs) {
	nmod_poly_sub(res, lhs, rhs);
}

/*************************************************
* Name:        poly_qshow_mul
*
* Description: Multiplication of two polynomials reduced
*              modulo x^PARAM_N_SHOW + 1.
*              Coefficients are reduced mod PARAM_Q_SHOW
* 
* Arguments:   - poly_qshow res: polynomial to host the multiplication (initialized)
*              - const poly_qshow lhs: first polynomial factor
*              - const poly_qshow rhs: second polynomial factor
**************************************************/
void poly_qshow_mul(poly_qshow res, const poly_qshow lhs, const poly_qshow rhs) {
	ASSERT_DEBUG(nmod_poly_degree(lhs) < PARAM_N_SHOW, "Argument to `poly_qshow_mul` must already be reduced for FLINT.");
	ASSERT_DEBUG(nmod_poly_degree(rhs) < PARAM_N_SHOW, "Argument to `poly_qshow_mul` must already be reduced for FLINT.");
	nmod_poly_mulmod_preinv(res, lhs, rhs, POLY_F, POLY_F_REV_INV);
}

/*************************************************
* Name:        poly_qshow_mul_x
*
* Description: Multiplication of a polynomial by x and reduced
*              modulo x^PARAM_N_SHOW + 1.
*              Coefficients are reduced mod PARAM_Q_SHOW
* 
* Arguments:   - poly_qshow res: polynomial to host the multiplication (initialized)
*              - const poly_qshow arg: polynomial to be multiplied by x
**************************************************/
void poly_qshow_mul_x(poly_qshow res, const poly_qshow arg) {
  nmod_poly_shift_left(res, arg, 1);
  nmod_poly_set_coeff_ui(res, 0, PARAM_Q_SHOW-nmod_poly_get_coeff_ui(res, PARAM_N_SHOW));
  nmod_poly_set_coeff_ui(res, PARAM_N_SHOW, 0);
}

/*************************************************
* Name:        poly_qshow_mul_xj
*
* Description: Multiplication of a polynomial by x^j and reduced
*              modulo x^PARAM_N_SHOW + 1.
*              Coefficients are reduced mod PARAM_Q_SHOW
* 
* Arguments:   - poly_qiss res: polynomial to host the multiplication (initialized)
*              - const poly_qiss arg: polynomial to be multiplied by x
*              - const size_t j: exponent to multiply by x^j
**************************************************/
void poly_qshow_mul_xj(poly_qshow res, const poly_qshow arg, const size_t j) {
  for (size_t i = 0; i < j; i++) {
    poly_qshow_set_coeff(res, i, -poly_qshow_get_coeff(arg, PARAM_N_SHOW + i - j));
  }
  for (size_t i = j; i < PARAM_N_SHOW; i++) {
    poly_qshow_set_coeff(res, i, poly_qshow_get_coeff(arg, i - j));
  }
}

/*************************************************
* Name:        poly_qshow_mul_scalar
*
* Description: Multiplication of a polynomials by a scalar
* 
* Arguments:   - poly_qshow res: polynomial to host the multiplication (initialized)
*              - const poly_qshow lhs: first polynomial factor
*              - coeff_qshow rhs: second scalar factor
**************************************************/
void poly_qshow_mul_scalar(poly_qshow out, const poly_qshow lhs, const coeff_qshow rhs) {
  coeff_qshow f = rhs % PARAM_Q_SHOW;
  f += (f >> (sizeof(coeff_qshow)*8 - 1)) & PARAM_Q_SHOW;
  ASSERT_DEBUG(f >= 0, "Factors for scalar multiplication must not be negative.");
  nmod_poly_scalar_mul_nmod(out, lhs, f);
}

/*************************************************
* Name:        poly_qshow_muladd_constant
*
* Description: Increment constant coefficient of a polynomial
*              by the product of two scalar integers
* 
* Arguments:   - poly_qshow arg: polynomial to be incremented (initialized)
*              - const coeff_qshow c0_lhs: first scalar factor
*              - const coeff_qshow c0_rhs: second scalar factor
**************************************************/
void poly_qshow_muladd_constant(poly_qshow arg, const coeff_qshow c0_lhs, const coeff_qshow c0_rhs) {
  int128 tmp = ((int128)c0_lhs * (int128)c0_rhs) % PARAM_Q_SHOW;
  tmp += PARAM_Q_SHOW;
  tmp %= PARAM_Q_SHOW;
  assert(tmp >= 0);
  tmp += nmod_poly_get_coeff_ui(arg, 0);
  tmp %= PARAM_Q_SHOW;
  nmod_poly_set_coeff_ui(arg, 0, tmp);
}

/*************************************************
* Name:        poly_qshow_shift_left
*
* Description: Shift the coefficients of a polynomial to
*              the left by n places. Corresponds to a 
*              multiplication by x^n (NOT REDUCED MOD X^PARAM_N_SHOW + 1)
*              Trailing zeros are inserted
* 
* Arguments:   - poly_qshow res: polynomial to host the shift (initialized)
*              - const poly_qshow arg: polynomial to be shifted
*              - size_t n: amount of shift
**************************************************/
void poly_qshow_shift_left(poly_qshow res, const poly_qshow arg, size_t n) {
  nmod_poly_shift_left(res, arg, n);
}

/*************************************************
* Name:        poly_qshow_conjugate
*
* Description: Compute the conjugate of a polynomial (evaluation at x^-1)
* 
* Arguments:   - poly_qshow res: polynomial to host the conjugate (initialized)
*              - const poly_qshow arg: polynomial to be conjugated
**************************************************/
void poly_qshow_conjugate(poly_qshow res, const poly_qshow arg) {
  coeff_qshow c = poly_qshow_get_coeff(arg, 0), d = poly_qshow_get_coeff(arg, PARAM_N_SHOW/2);
  poly_qshow_set_coeff(res, 0, c);
  poly_qshow_set_coeff(res, PARAM_N_SHOW/2, PARAM_Q_SHOW-d);
  for (size_t i = 1; i < PARAM_N_SHOW/2; ++i) {
    c = poly_qshow_get_coeff(arg, PARAM_N_SHOW - i);
    d = poly_qshow_get_coeff(arg, i);
    poly_qshow_set_coeff(res,                i, PARAM_Q_SHOW-c);
    poly_qshow_set_coeff(res, PARAM_N_SHOW - i, PARAM_Q_SHOW-d);
  }
}

/*************************************************
* Name:        poly_qshow_decompose
*
* Description: Compute the decomposition x_L + PARAM_GAMMA_SHOW.x_H of a polynomial x
*              with x_L with coefficients in (-PARAM_GAMMA_SHOW/2, PARAM_GAMMA_SHOW/2]
*              except if coefficients of x_H = (PARAM_Q_SHOW-1)/PARAM_GAMMA_SHOW where we set 
*              the coefficient of x_H to 0 and that of x_L to x mod+ PARAM_Q_SHOW
* 
* Arguments:   - poly_qshow high: polynomial to host the high order decomposition (initialized)
*              - poly_qshow low: polynomial to host the low order decomposition (initialized)
*              - const poly_qshow arg: polynomial to be decomposed
**************************************************/
void poly_qshow_decompose(poly_qshow high, poly_qshow low, const poly_qshow arg) {
  coeff_qshow c, low_c;
  for (size_t i = 0; i < PARAM_N_SHOW; i++) {
    c = poly_qshow_get_coeff(arg, i); // must be on uncentered representation
    poly_qshow_set_coeff(high, i, decompose_show(&low_c, c));
    poly_qshow_set_coeff(low, i, low_c);
  }
}

/*************************************************
* Name:        poly_qshow_power2round
*
* Description: Compute the decomposition x_L + 2^PARAM_D_ROUND_SHOW.x_H of a polynomial x
*              with x_L with coefficients in (-2^(PARAM_D_ROUND_SHOW-1), 2^(PARAM_D_ROUND_SHOW-1)]
* 
* Arguments:   - poly_qshow high: polynomial to host the high order decomposition (initialized)
*              - poly_qshow low: polynomial to host the low order decomposition (initialized)
*              - const poly_qshow arg: polynomial to be decomposed
**************************************************/
void poly_qshow_power2round(poly_qshow high, poly_qshow low, const poly_qshow arg) {
  coeff_qshow c, low_c;
  for (size_t i = 0; i < PARAM_N_SHOW; i++) {
    c = poly_qshow_get_coeff(arg, i); // must be on uncentered representation
    poly_qshow_set_coeff(high, i, power2round_show(&low_c, c));
    poly_qshow_set_coeff(low, i, low_c);
  }
}

/*************************************************
* Name:        poly_qshow_makeGhint
*
* Description: 
* 
* Arguments:   - poly_qshow hint: polynomial to host the hint polynomial (initialized)
*              - const poly_qshow z: polynomial 
*              - const poly_qshow r: polynomial 
**************************************************/
void poly_qshow_makeGhint(poly_qshow hint, const poly_qshow z, const poly_qshow r) {
  coeff_qshow cz, cr;
  for (size_t i = 0; i < PARAM_N_SHOW; i++) {
    cz = poly_qshow_get_coeff(z, i); // must be on uncentered representation for Highbits
    cr = poly_qshow_get_coeff(r, i); // must be on uncentered representation for Hightbits
    poly_qshow_set_coeff(hint, i, makeGhint_show(cz,cr));
  }
}

/*************************************************
* Name:        poly_qshow_useGhint
*
* Description: 
* 
* Arguments:   - poly_qshow high_zr: polynomial to host the hint polynomial (initialized)
*              - const poly_qshow hint: polynomial 
*              - const poly_qshow r: polynomial 
**************************************************/
void poly_qshow_useGhint(poly_qshow high_zr, const poly_qshow hint, const poly_qshow r) {
  coeff_qshow chint, cr;
  for (size_t i = 0; i < PARAM_N_SHOW; i++) {
    chint = poly_qshow_get_coeff_centered(hint, i); // must be on centered representation
    cr = poly_qshow_get_coeff(r, i); // must be on uncentered representation for Highbits
    poly_qshow_set_coeff(high_zr, i, useGhint_show(chint,cr));
  }
}

/*************************************************
* Name:        poly_qshow_equal
*
* Description: Equality test between two polynomials
* 
* Arguments:   - const poly_qshow lhs: first polynomial
*              - const poly_qshow rhs: second polynomial
* 
* Returns 1 if the polynomials are equal, 0 otherwise
**************************************************/
int poly_qshow_equal(const poly_qshow lhs, const poly_qshow rhs) {
  return nmod_poly_equal(lhs, rhs);
}

/*************************************************
* Name:        poly_qshow_dump
*
* Description: Print a polynomial
* 
* Arguments:   - const poly_qshow arg: the polynomial to be printed
**************************************************/
void poly_qshow_dump(const poly_qshow arg) {
  nmod_poly_print(arg);
}

/*************************************************
* Name:        poly_qshow_sq_norm2
*
* Description: Compute the square l2 norm of a polynomial
* 
* Arguments:   - const poly_qshow arg: the polynomial
* 
* Returns an unsigned 128-bit integer with the square l2 norm
**************************************************/
uint128 poly_qshow_sq_norm2(const poly_qshow arg) {
	uint128 sq_norm2 = 0;
	coeff_qshow c;
	for (size_t i = 0; i < PARAM_N_SHOW; i++) {
		c = poly_qshow_get_coeff_centered(arg, i);
    sq_norm2 += (uint128)c * (uint128)c;
	}
	return sq_norm2;
}

/*************************************************
* Name:        poly_qshow_infnorm
*
* Description: Compute the infinity norm of a polynomial
*
* Arguments:   - const poly_qshow arg: the polynomial
* 
* Returns an 64-bit integer with the linf norm
**************************************************/
int64_t poly_qshow_norm_inf(const poly_qshow arg) {
  int64_t max, tmpmax;
  max = ct_abs(poly_qshow_get_coeff_centered(arg, 0));
  for (size_t i = 1; i < PARAM_N_ISS; i++) {
    tmpmax = ct_abs(poly_qshow_get_coeff_centered(arg, i));
    max = (tmpmax > max) ? tmpmax : max;
  }
  return max;
}


/*************************************************
* Name:        poly_qshow_pack
*
* Description: Pack a polynomial mod PARAM_Q_SHOW into a byte array
* 
* Arguments:   - uint8_t buf: output byte array (allocated POLYQSHOW_PACKEDBYTES bytes)
*              - const poly_qshow arg: the polynomial to be packed
**************************************************/
void poly_qshow_pack(uint8_t buf[POLYQSHOW_PACKEDBYTES], const poly_qshow arg) {
  uint64_t x;
  for (size_t i = 0; i < PARAM_N_SHOW; i++) {
    x = poly_qshow_get_coeff(arg, i);
#if PARAM_Q_SHOW_BITLEN > 64
#error "poly_qshow_pack assumes that PARAM_Q_SHOW_BITLEN <= 64"
#endif
    buf[8*i + 0] = x >>  0;
    buf[8*i + 1] = x >>  8;
    buf[8*i + 2] = x >> 16;
    buf[8*i + 3] = x >> 24;
    buf[8*i + 4] = x >> 32;
    buf[8*i + 5] = x >> 40;
    buf[8*i + 6] = x >> 48;
    buf[8*i + 7] = x >> 56;
  }
}

/*************************************************
* Name:        coeff_qshow_pack
*
* Description: Pack a polynomial coefficient mod PARAM_Q_SHOW into a byte array
* 
* Arguments:   - uint8_t buf: output byte array (allocated COEFFQSHOW_PACKEDBYTES bytes)
*              - const coeff_qshow arg: the coefficient to be packed
**************************************************/
void coeff_qshow_pack(uint8_t buf[COEFFQSHOW_PACKEDBYTES], const coeff_qshow arg) {
  uint64_t x;
  coeff_qshow tmp = arg + ((arg >> (sizeof(coeff_qshow)*8-1)) & PARAM_Q_SHOW);
  x = tmp % PARAM_Q_SHOW;
#if PARAM_Q_SHOW_BITLEN > 64
#error "coeff_qshow_pack assumes that PARAM_Q_SHOW_BITLEN <= 64"
#endif
  buf[0] = x >>  0;
  buf[1] = x >>  8;
  buf[2] = x >> 16;
  buf[3] = x >> 24;
  buf[4] = x >> 32;
  buf[5] = x >> 40;
  buf[6] = x >> 48;
  buf[7] = x >> 56;
}
