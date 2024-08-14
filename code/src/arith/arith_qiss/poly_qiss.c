#include "arith_q.h"
#include "arith_qiss.h"
#include "macros.h"

static nmod_poly_t POLY_F;
static nmod_poly_t POLY_F_REV_INV;

/*************************************************
* Name:        ct_abs
*
* Description: Absolute value
*
* Arguments:   - coeff_qiss c: coefficient 
* 
* Returns an 64-bit integer with the absolute value of c
**************************************************/
static int64_t ct_abs(coeff_qiss c) {
  int64_t const mask = c >> (sizeof(coeff_qiss)*8 - 1);
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
* Name:        decompose_iss
*
* Description: Decompose x into x_L + PARAM_GAMMA_ISS.x_H with x_L 
*              in (-PARAM_GAMMA_ISS/2, PARAM_GAMMA_ISS/2] except if 
*              x_H = (PARAM_Q_ISS-1)/PARAM_GAMMA_ISS where we set 
*              x_H = 0 and x_L = x mod+ PARAM_Q_ISS
*
* Arguments:   - coeff_qiss *low_c: coefficient to host low part of c
*              - const coeff_qiss c: coefficient to be decomposed
* 
* Returns high_c the high part of c
**************************************************/
static coeff_qiss decompose_iss(coeff_qiss *low_c, const coeff_qiss c) {
  coeff_qiss high_c;
  high_c = (c + (PARAM_GAMMA_ISS>>1) - 1) / PARAM_GAMMA_ISS; // floor division here
  high_c -= (((PARAM_Q_GAMMA_ISS - high_c - 1) >> (sizeof(coeff_qiss)*8-1)) & high_c); // putting high_c to 0 if high_c = (PARAM_Q_ISS-1)/PARAM_GAMMA_ISS
  *low_c = c - PARAM_GAMMA_ISS * high_c;
  return high_c;
}

/*************************************************
* Name:        highbits_iss
*
* Description: Decompose x into x_L + PARAM_GAMMA_ISS.x_H with x_L 
*              in (-PARAM_GAMMA_ISS/2, PARAM_GAMMA_ISS/2] except if 
*              x_H = (PARAM_Q_ISS-1)/PARAM_GAMMA_ISS where we set 
*              x_H = 0 and x_L = x mod+ PARAM_Q_ISS. Return only x_H
*
* Arguments:   - const coeff_qiss c: coefficient to be decomposed
* 
* Returns high_c the high part of c
**************************************************/
static coeff_qiss highbits_iss(const coeff_qiss c) {
  coeff_qiss high_c = (c + PARAM_Q_ISS) % PARAM_Q_ISS;
  high_c = (high_c + (PARAM_GAMMA_ISS>>1) - 1) / PARAM_GAMMA_ISS; // floor division here
  return high_c - ((((coeff_qiss)PARAM_Q_GAMMA_ISS - high_c - 1) >> (sizeof(coeff_qiss)*8-1)) & high_c); // putting high_c to 0 if high_c = (PARAM_Q_ISS-1)/PARAM_GAMMA_ISS;
}

/*************************************************
* Name:        power2round_iss
*
* Description: Decompose x into x_L + 2^PARAM_D_ROUND_ISS.x_H with x_L 
*              in (-2^(PARAM_D_ROUND_ISS-1), 2^(PARAM_D_ROUND_ISS-1)]
*
* Arguments:   - coeff_qiss *low_c: coefficient to host low part of c
*              - const coeff_qiss c: coefficient to be decomposed
* 
* Returns high_c the high part of c
**************************************************/
static coeff_qiss power2round_iss(coeff_qiss *low_c, const coeff_qiss c) {
  coeff_qiss high_c;
  high_c = (c + (1 << (PARAM_D_ROUND_ISS-1)) - 1) >> PARAM_D_ROUND_ISS;
  *low_c = c - (high_c << PARAM_D_ROUND_ISS); // high_c positive so no runtime error when left shift
  return high_c;
}

/*************************************************
* Name:        makeGhint_iss
*
* Description: 
*
* Arguments:   - const coeff_qiss z: 
*              - const coeff_qiss r: 
* 
* Returns h_c the hint coeff 
**************************************************/
static coeff_qiss makeGhint_iss(const coeff_qiss z, const coeff_qiss r) {
  coeff_qiss h_c;
  h_c = (highbits_iss(z+r) - highbits_iss(r)) % PARAM_Q_GAMMA_ISS;
  return h_c - ((~((h_c - (PARAM_Q_GAMMA_ISS>>1) - 1) >> (sizeof(coeff_qiss)*8-1))) & PARAM_Q_GAMMA_ISS);
}

/*************************************************
* Name:        useGhint_iss
*
* Description: 
*
* Arguments:   - const coeff_qiss h_c: 
*              - const coeff_qiss r: 
* 
* Returns high_zr high part of z+r
**************************************************/
static coeff_qiss useGhint_iss(const coeff_qiss hint_c, const coeff_qiss r) {
  return (hint_c + highbits_iss(r)) % PARAM_Q_GAMMA_ISS;
}

/*************************************************
* Name:        arith_qiss_setup
*
* Description: Initialize and setup the backend for arithmetic 
*              modulo PARAM_Q_ISS. This is strictly required 
*              and must be called once before any other function 
*              from here is used.
**************************************************/
void arith_qiss_setup(void) {
	nmod_poly_init(POLY_F, PARAM_Q_ISS);
	nmod_poly_set_coeff_ui(POLY_F, 0, 1);
	nmod_poly_set_coeff_ui(POLY_F, PARAM_N_ISS, 1);

  nmod_poly_t f_len;
  nmod_poly_init(f_len, PARAM_Q_ISS);
  nmod_poly_set_coeff_ui(f_len, PARAM_N_ISS + 1, 1);

  nmod_poly_t f_rev;
  nmod_poly_init(f_rev, PARAM_Q_ISS);

  nmod_poly_init(POLY_F_REV_INV, PARAM_Q_ISS);
  nmod_poly_reverse(f_rev, POLY_F, PARAM_N_ISS + 1);
  nmod_poly_invert(POLY_F_REV_INV, f_rev, f_len, PARAM_Q_ISS);

  nmod_poly_clear(f_len);
  nmod_poly_clear(f_rev);
}

/*************************************************
* Name:        arith_qiss_teardown
*
* Description: Clean up and teardown the backend for arithmetic 
*              modulo PARAM_Q_ISS. This is strictly required 
*              and must be called once at the very end to release 
*              any resources.
**************************************************/
void arith_qiss_teardown(void) {
  nmod_poly_clear(POLY_F);
  nmod_poly_clear(POLY_F_REV_INV);
}

/*************************************************
* Name:        poly_qiss_init
*
* Description: Initialize polynomial and set it to zero
*              This is strictly required before any operations 
*              are done with/on the polynomial.
* 
* Arguments:   - poly_qiss res: polynomial to be initialized
**************************************************/
void poly_qiss_init(poly_qiss res) {
  nmod_poly_init(res, PARAM_Q_ISS);
}

/*************************************************
* Name:        poly_qiss_clear
*
* Description: Clears a polynomial and releases all associated memory. 
*              This is strictly required to avoid memory leaks and the 
*              polynomial must not be used again (unless reinitialized).
* 
* Arguments:   - poly_qiss arg: polynomial to be cleared
**************************************************/
void poly_qiss_clear(poly_qiss arg) {
  nmod_poly_clear(arg);
}

/*************************************************
* Name:        poly_qiss_zero
*
* Description: Set an initialized polynomial to zero
* 
* Arguments:   - poly_qiss res: polynomial to be zeroized (initialized)
**************************************************/
void poly_qiss_zero(poly_qiss res) {
  nmod_poly_zero(res);
}

/*************************************************
* Name:        poly_qiss_set
*
* Description: Set a polynomial equal to another polynomial.
*              Coefficients are reduced mod PARAM_Q_ISS
* 
* Arguments:   - poly_qiss res: polynomial to be set (initialized)
*              - const poly_qiss arg: polynomial to be read
**************************************************/
void poly_qiss_set(poly_qiss res, const poly_qiss arg) {
	nmod_poly_set(res, arg);
}

/*************************************************
* Name:        poly_qiss_get_coeff
*
* Description: Get coefficient of x^n of a polynomial
*              condition: [0 <= n < PARAM_N_ISS]
* 
* Arguments:   - const poly_qiss arg: polynomial to be read
*              - size_t n: degree of the coefficient to be read
* 
* Returns the coefficients of x^n of arg
**************************************************/
coeff_qiss poly_qiss_get_coeff(const poly_qiss arg, size_t n) {
  ASSERT_DEBUG(n < PARAM_N_ISS, "Illegal argument: cannot get coefficient of poly at given position.");
	return nmod_poly_get_coeff_ui(arg, n);
}

/*************************************************
* Name:        poly_qiss_get_coeff_centered
*
* Description: Get coefficient of x^n of a polynomial in centered representation
*              condition: [0 <= n < PARAM_N]
* 
* Arguments:   - const poly_qiss arg: polynomial to be read
*              - size_t n: degree of the coefficient to be read
* 
* Returns the coefficients of x^n of arg in [-PARAM_Q_ISS/2, PARAM_Q_ISS/2]
**************************************************/
coeff_qiss poly_qiss_get_coeff_centered(const poly_qiss arg, size_t n) {
	coeff_qiss tmp = poly_qiss_get_coeff(arg, n);
	return tmp - ((~((tmp - PARAM_Q_ISS/2) >> (sizeof(coeff_qiss)*8-1))) & PARAM_Q_ISS);
}

/*************************************************
* Name:        poly_qiss_set_coeff
*
* Description: Set coefficient of x^n of a polynomial
*              condition: [0 <= n < PARAM_N_ISS]
*              Coefficient is reduced mod PARAM_Q_ISS
* 
* Arguments:   - poly_qiss arg: polynomial whose n-th coefficient is set (initialized)
*              - size_t n: degree of the coefficient to be set
*              - coeff_qiss c: the new coefficient
**************************************************/
void poly_qiss_set_coeff(poly_qiss arg, size_t n, coeff_qiss c) {
  ASSERT_DEBUG(n < PARAM_N_ISS, "Illegal argument: cannot set coefficient of poly at given position.");
	c += PARAM_Q_ISS;
	c %= PARAM_Q_ISS;
  c += (c >> (sizeof(coeff_qiss)*8 - 1)) & PARAM_Q_ISS;
  ASSERT_DEBUG(c >= 0, "Coefficients of polys must not be negative.");
  nmod_poly_set_coeff_ui(arg, n, c);
}

/*************************************************
* Name:        poly_qiss_neg
*
* Description: Negate a polynomial coefficient-wise
*              Coefficients are reduced mod PARAM_Q_ISS
* 
* Arguments:   - poly_qiss res: polynomial to host the negation (initialized)
*              - const poly_qiss arg: polynomial to be negated
**************************************************/
void poly_qiss_neg(poly_qiss res, const poly_qiss arg) {
  nmod_poly_neg(res, arg);
}

/*************************************************
* Name:        poly_qiss_add
*
* Description: Add two polynomials. Coefficients are reduced mod PARAM_Q_ISS
* 
* Arguments:   - poly_qiss res: polynomial to host the sum (initialized)
*              - const poly_qiss lhs: first polynomial summand
*              - const poly_qiss rhs: second polynomial summand
**************************************************/
void poly_qiss_add(poly_qiss res, const poly_qiss lhs, const poly_qiss rhs) {
	nmod_poly_add(res, lhs, rhs);
}

/*************************************************
* Name:        poly_qiss_sub
*
* Description: Substract two polynomials. Coefficients are reduced mod PARAM_Q_ISS
* 
* Arguments:   - poly_qiss res: polynomial to host the difference (initialized)
*              - const poly_qiss lhs: first polynomial term
*              - const poly_qiss rhs: second polynomial term
**************************************************/
void poly_qiss_sub(poly_qiss res, const poly_qiss lhs, const poly_qiss rhs) {
	nmod_poly_sub(res, lhs, rhs);
}

/*************************************************
* Name:        poly_qiss_mul
*
* Description: Multiplication of two polynomials reduced
*              modulo x^PARAM_N_ISS + 1.
*              Coefficients are reduced mod PARAM_Q_ISS
* 
* Arguments:   - poly_qiss res: polynomial to host the multiplication (initialized)
*              - const poly_qiss lhs: first polynomial factor
*              - const poly_qiss rhs: second polynomial factor
**************************************************/
void poly_qiss_mul(poly_qiss res, const poly_qiss lhs, const poly_qiss rhs) {
  ASSERT_DEBUG(nmod_poly_degree(lhs) < PARAM_N_ISS, "Argument to `poly_qiss_mul` must already be reduced for FLINT.");
  ASSERT_DEBUG(nmod_poly_degree(rhs) < PARAM_N_ISS, "Argument to `poly_qiss_mul` must already be reduced for FLINT.");
	nmod_poly_mulmod_preinv(res, lhs, rhs, POLY_F, POLY_F_REV_INV);
}

/*************************************************
* Name:        poly_qiss_mul_x
*
* Description: Multiplication of a polynomial by x and reduced
*              modulo x^PARAM_N_ISS + 1.
*              Coefficients are reduced mod PARAM_Q_ISS
* 
* Arguments:   - poly_qiss res: polynomial to host the multiplication (initialized)
*              - const poly_qiss arg: polynomial to be multiplied by x
**************************************************/
void poly_qiss_mul_x(poly_qiss res, const poly_qiss arg) {
  coeff_qiss c = poly_qiss_get_coeff(arg, PARAM_N_ISS - 1);
  poly_qiss_shift_left(res, arg, 1);
  poly_qiss_set_coeff(res, 0, -c);
  // removing leading coefficient manually as poly_qiss_set_coeff only authorizes exponent < PARAM_N_ISS
  nmod_poly_set_coeff_ui(res, PARAM_N_ISS, 0);
}

/*************************************************
* Name:        poly_qiss_mul_xj
*
* Description: Multiplication of a polynomial by x^j and reduced
*              modulo x^PARAM_N_ISS + 1.
*              Coefficients are reduced mod PARAM_Q_ISS
* 
* Arguments:   - poly_qiss res: polynomial to host the multiplication (initialized)
*              - const poly_qiss arg: polynomial to be multiplied by x
*              - const size_t j: exponent to multiply by x^j
**************************************************/
void poly_qiss_mul_xj(poly_qiss res, const poly_qiss arg, const size_t j) {
  for (size_t i = 0; i < j; i++) {
    poly_qiss_set_coeff(res, i, -poly_qiss_get_coeff(arg, PARAM_N_ISS + i - j));
  }
  for (size_t i = j; i < PARAM_N_ISS; i++) {
    poly_qiss_set_coeff(res, i, poly_qiss_get_coeff(arg, i - j));
  }
}

/*************************************************
* Name:        poly_qiss_mul_scalar
*
* Description: Multiplication of a polynomials by a scalar
* 
* Arguments:   - poly_qiss res: polynomial to host the multiplication (initialized)
*              - const poly_qiss arg: first polynomial factor
*              - coeff_qiss fac: second scalar factor
**************************************************/
void poly_qiss_mul_scalar(poly_qiss out, const poly_qiss lhs, const coeff_qiss rhs) {
  coeff_qiss f = rhs % PARAM_Q_ISS;
  f += (f >> (sizeof(coeff_qiss)*8 - 1)) & PARAM_Q_ISS;
  ASSERT_DEBUG(f >= 0, "Factors for scalar multiplication must not be negative.");
  nmod_poly_scalar_mul_nmod(out, lhs, f);
}

/*************************************************
* Name:        poly_qiss_muladd_constant
*
* Description: Increment constant coefficient of a polynomial
*              by the product of two scalar integers
* 
* Arguments:   - poly_qiss arg: polynomial to be incremented (initialized)
*              - const coeff_qiss c0_lhs: first scalar factor
*              - const coeff_qiss c0_rhs: second scalar factor
**************************************************/
void poly_qiss_muladd_constant(poly_qiss arg, const coeff_qiss c0_lhs, const coeff_qiss c0_rhs) {
  int128 tmp = ((int128)c0_lhs * (int128)c0_rhs) % PARAM_Q_ISS;
  tmp += PARAM_Q_ISS;
  tmp %= PARAM_Q_ISS;
  assert(tmp >= 0);
  tmp += nmod_poly_get_coeff_ui(arg, 0);
  tmp %= PARAM_Q_ISS;
  nmod_poly_set_coeff_ui(arg, 0, tmp);
}

/*************************************************
* Name:        poly_qiss_shift_left
*
* Description: Shift the coefficients of a polynomial to
*              the left by n places. Corresponds to a 
*              multiplication by x^n (NOT REDUCED MOD X^PARAM_N_ISS + 1)
*              Trailing zeros are inserted
* 
* Arguments:   - poly_qiss res: polynomial to host the shift (initialized)
*              - const poly_qiss arg: polynomial to be shifted
*              - size_t n: amount of shift
**************************************************/
void poly_qiss_shift_left(poly_qiss res, const poly_qiss arg, size_t n) {
  nmod_poly_shift_left(res, arg, (signed long) n);
}

/*************************************************
* Name:        poly_qiss_conjugate
*
* Description: Compute the conjugate of a polynomial (evaluation at x^-1)
* 
* Arguments:   - poly_qiss res: polynomial to host the conjugate (initialized)
*              - const poly_qiss arg: polynomial to be conjugated
**************************************************/
void poly_qiss_conjugate(poly_qiss res, const poly_qiss arg) {
  ASSERT_DEBUG(res != arg, "Input and output of conjugation must not be idential.");
  coeff_qiss c = poly_qiss_get_coeff(arg, 0);
  poly_qiss_set_coeff(res, 0, c);
  for (size_t i = 1; i < PARAM_N_ISS; ++i) {
    c = poly_qiss_get_coeff_centered(arg, PARAM_N_ISS - i);
    poly_qiss_set_coeff(res, i, -c);
  }
}

/*************************************************
* Name:        poly_qiss_decompose
*
* Description: Compute the decomposition x_L + PARAM_GAMMA_ISS.x_H of a polynomial x
*              with x_L with coefficients in (-PARAM_GAMMA_ISS/2, PARAM_GAMMA_ISS/2]
*              except if coefficients of x_H = (PARAM_Q_ISS-1)/PARAM_GAMMA_ISS where we set 
*              the coefficient of x_H to 0 and that of x_L to x mod+ PARAM_Q_ISS
* 
* Arguments:   - poly_qiss high: polynomial to host the high order decomposition (initialized)
*              - poly_qiss low: polynomial to host the low order decomposition (initialized)
*              - const poly_qiss arg: polynomial to be decomposed
**************************************************/
void poly_qiss_decompose(poly_qiss high, poly_qiss low, const poly_qiss arg) {
  coeff_qiss c, low_c;
  for (size_t i = 0; i < PARAM_N_ISS; i++) {
    c = poly_qiss_get_coeff(arg, i); // must be on uncentered representation
    poly_qiss_set_coeff(high, i, decompose_iss(&low_c, c));
    poly_qiss_set_coeff(low, i, low_c);
  }
}

/*************************************************
* Name:        poly_qiss_power2round
*
* Description: Compute the decomposition x_L + 2^PARAM_D_ROUND_ISS.x_H of a polynomial x
*              with x_L with coefficients in (-2^(PARAM_D_ROUND_ISS-1), 2^(PARAM_D_ROUND_ISS-1)]
* 
* Arguments:   - poly_qiss high: polynomial to host the high order decomposition (initialized)
*              - poly_qiss low: polynomial to host the low order decomposition (initialized)
*              - const poly_qiss arg: polynomial to be decomposed
**************************************************/
void poly_qiss_power2round(poly_qiss high, poly_qiss low, const poly_qiss arg) {
  coeff_qiss c, low_c;
  for (size_t i = 0; i < PARAM_N_ISS; i++) {
    c = poly_qiss_get_coeff(arg, i); // must be on uncentered representation
    poly_qiss_set_coeff(high, i, power2round_iss(&low_c, c));
    poly_qiss_set_coeff(low, i, low_c);
  }
}

/*************************************************
* Name:        poly_qiss_makeGhint
*
* Description: 
* 
* Arguments:   - poly_qiss hint: polynomial to host the hint polynomial (initialized)
*              - const poly_qiss z: polynomial 
*              - const poly_qiss r: polynomial 
**************************************************/
void poly_qiss_makeGhint(poly_qiss hint, const poly_qiss z, const poly_qiss r) {
  coeff_qiss cz, cr;
  for (size_t i = 0; i < PARAM_N_ISS; i++) {
    cz = poly_qiss_get_coeff(z, i); // must be on uncentered representation for Highbits
    cr = poly_qiss_get_coeff(r, i); // must be on uncentered representation for Hightbits
    poly_qiss_set_coeff(hint, i, makeGhint_iss(cz,cr));
  }
}

/*************************************************
* Name:        poly_qiss_useGhint
*
* Description: 
* 
* Arguments:   - poly_qiss high_zr: polynomial to host the hint polynomial (initialized)
*              - const poly_qiss hint: polynomial 
*              - const poly_qiss r: polynomial 
**************************************************/
void poly_qiss_useGhint(poly_qiss high_zr, const poly_qiss hint, const poly_qiss r) {
  coeff_qiss chint, cr;
  for (size_t i = 0; i < PARAM_N_ISS; i++) {
    chint = poly_qiss_get_coeff_centered(hint, i); // must be on centered representation
    cr = poly_qiss_get_coeff(r, i); // must be on uncentered representation for Highbits
    poly_qiss_set_coeff(high_zr, i, useGhint_iss(chint,cr));
  }
}

/*************************************************
* Name:        poly_qiss_equal
*
* Description: Equality test between two polynomials
* 
* Arguments:   - const poly_qiss lhs: first polynomial
*              - const poly_qiss rhs: second polynomial
* 
* Returns 1 if the polynomials are equal, 0 otherwise
**************************************************/
int poly_qiss_equal(const poly_qiss lhs, const poly_qiss rhs) {
  return nmod_poly_equal(lhs, rhs);
}

/*************************************************
* Name:        poly_qiss_dump
*
* Description: Print a polynomial
* 
* Arguments:   - const poly_qiss arg: the polynomial to be printed
**************************************************/
void poly_qiss_dump(const poly_qiss arg) {
	nmod_poly_print(arg);
}

/*************************************************
* Name:        poly_qiss_sq_norm2
*
* Description: Compute the square l2 norm of a polynomial
* 
* Arguments:   - const poly_qiss arg: the polynomial
* 
* Returns an unsigned 128-bit integer with the square l2 norm
**************************************************/
uint128 poly_qiss_sq_norm2(const poly_qiss arg) {
	uint128 sq_norm2 = 0;
	coeff_qiss c;
	for (size_t i = 0; i < PARAM_N_ISS; i++) {
		c = poly_qiss_get_coeff_centered(arg, i);
    sq_norm2 += (uint128)c * (uint128)c;
	}
	return sq_norm2;
}

/*************************************************
* Name:        poly_qiss_infnorm
*
* Description: Compute the infinity norm of a polynomial
*
* Arguments:   - const poly_qiss arg: the polynomial
* 
* Returns an 64-bit integer with the linf norm
**************************************************/
int64_t poly_qiss_norm_inf(const poly_qiss arg) {
  int64_t max, tmpmax;
  max = ct_abs(poly_qiss_get_coeff_centered(arg, 0));
  for (size_t i = 1; i < PARAM_N_ISS; i++) {
    tmpmax = ct_abs(poly_qiss_get_coeff_centered(arg, i));
    max = (tmpmax > max) ? tmpmax : max;
  }
  return max;
}

/*************************************************
* Name:        poly_qiss_pack
*
* Description: Pack a polynomial mod PARAM_Q_ISS into a byte array
* 
* Arguments:   - uint8_t buf: output byte array (allocated POLYQISS_PACKEDBYTES bytes)
*              - const poly_qiss arg: the polynomial to be packed
**************************************************/
void poly_qiss_pack(uint8_t buf[POLYQISS_PACKEDBYTES], const poly_qiss arg) {
  uint64_t x;
  for (size_t i = 0; i < PARAM_N_ISS; i++) {
	  x = poly_qiss_get_coeff(arg, i);
#if PARAM_Q_ISS_BITLEN > 64
#error "poly_qiss_pack assumes that PARAM_Q_ISS_BITLEN <= 64"
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
* Name:        coeff_qiss_pack
*
* Description: Pack a polynomial coefficient mod PARAM_Q_ISS into a byte array
* 
* Arguments:   - uint8_t buf: output byte array (allocated COEFFQISS_PACKEDBYTES bytes)
*              - const coeff_qiss arg: the coefficient to be packed
**************************************************/
void coeff_qiss_pack(uint8_t buf[COEFFQISS_PACKEDBYTES], const coeff_qiss arg) {
  uint64_t x;
  coeff_qiss tmp = arg + ((arg >> (sizeof(coeff_qiss)*8-1)) & PARAM_Q_ISS);
  x = tmp % PARAM_Q_ISS;
#if PARAM_Q_ISS_BITLEN > 64
#error "poly_qiss_pack assumes that PARAM_Q_ISS_BITLEN <= 64"
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
