#include "poly_p_sampling.h"
#include "fips202.h"
#include "arith.h"
#include "random.h"

/*************************************************
* Name:        poly_p_uniform
*
* Description: Sample a uniformly random polynomial modulo
*              PARAM_P deterministically from a seed.
* 
* Arguments:   - poly_p pout: output uniform polynomial (initialized)
*              - const uint8_t *seed: pointer to byte array containing seed (allocated SEED_BYTES bytes)
*              - uint32_t domain_separator: domain separator for XOF
*              - size_t i: index domain separator for XOF
*              - size_t j: index domain separator for XOF
**************************************************/
static void poly_p_uniform(poly_p pout, const uint8_t seed[SEED_BYTES], uint32_t domain_separator, size_t i, size_t j) {
  uint8_t output[SHAKE128_RATE * 2];
  keccak_state state;
  size_t k,cnt,off,bytecnt;
  shake128_init(&state);
  shake128_absorb(&state, seed, SEED_BYTES);
  shake128_absorb(&state, (const uint8_t*)&domain_separator, sizeof(uint32_t));
  shake128_absorb(&state, (const uint8_t*) &i, sizeof(i));
  shake128_absorb(&state, (const uint8_t*) &j, sizeof(j));
  shake128_finalize(&state);
  shake128_squeezeblocks(output, 2, &state);
  bytecnt = 2*SHAKE128_RATE;

  cnt = 0;
  off = 0;
  while (cnt < PARAM_N) {
#if PARAM_P_BITLEN > 13
#error "PARAM_P_BITLEN too big for uniform sampling."
#else
    // idea: take 5 byte, divide them into three partitions each of 13 bits, potentially ignore the MSBs, perform rejection sampling
    if (bytecnt < 5) {
      for (k = 0; k < bytecnt; k++) {
        output[k] = output[off++];
      }
      shake128_squeezeblocks(&output[bytecnt], 1, &state);
      off = 0;
      bytecnt += SHAKE128_RATE;
    }
    int64_t tmp5byte = (int64_t) (output[off] | ((uint64_t)output[off+1] << 8) | ((uint64_t)output[off+2] << 16) | ((uint64_t)output[off+3] << 24) | ((uint64_t)output[off+4] << 32));
    tmp5byte = tmp5byte & ((1UL<<40)-1);
    int64_t tmp = tmp5byte & ((1<<PARAM_P_BITLEN)-1);
    if (tmp < (coeff_p)PARAM_P) {
      poly_p_set_coeff(pout, cnt++, tmp);
      if (cnt == PARAM_N) {
        break;
      }
    }
    tmp5byte >>= PARAM_P_BITLEN;
    tmp = tmp5byte & ((1<<PARAM_P_BITLEN)-1);
    if (tmp < (coeff_p)PARAM_P) {
      poly_p_set_coeff(pout, cnt++, tmp);
      if (cnt == PARAM_N) {
        break;
      }
    }
    tmp5byte >>= PARAM_P_BITLEN;
    tmp = tmp5byte & ((1<<PARAM_P_BITLEN)-1);
    if (tmp < (coeff_p)PARAM_P) {
      poly_p_set_coeff(pout, cnt++, tmp);
    }

    off += 5;
    bytecnt -= 5;
#if PARAM_P_BITLEN < 12
#warning "PARAM_P_BITLEN maybe unsuitable for efficient uniform sampling."
#endif
#endif
  }
}

/*************************************************
* Name:        poly_p_mat_d_m_uniform
*
* Description: Sample a uniformly random polynomial matrix of
*              size PARAM_DE x PARAM_ME modulo PARAM_P deterministically 
*              from a seed.
* 
* Arguments:   - poly_p_mat_d_m mat: output uniform polynomial matrix (initialized)
*              - const uint8_t *seed: pointer to byte array containing seed (allocated SEED_BYTES bytes)
*              - uint32_t domain_separator: domain separator for XOF
**************************************************/
void poly_p_mat_d_m_uniform(poly_p_mat_d_m mat, const uint8_t seed[SEED_BYTES], uint32_t domain_separator) {
  size_t i,j;
  for (i = 0; i < PARAM_DE; i++) {
    for (j = 0; j < PARAM_ME; j++) {
      poly_p_uniform(mat->rows[i]->entries[j], seed, domain_separator, i, j);
    }
  }
}

/*************************************************
* Name:        poly_p_vec_m_uniform
*
* Description: Sample a uniformly random polynomial vector of
*              size PARAM_ME modulo PARAM_P deterministically 
*              from a seed.
* 
* Arguments:   - poly_p_vec_m vec: output uniform polynomial matrix (initialized)
*              - const uint8_t *seed: pointer to byte array containing seed (allocated SEED_BYTES bytes)
*              - uint32_t domain_separator: domain separator for XOF
**************************************************/
void poly_p_vec_m_uniform(poly_p_vec_m vec, const uint8_t seed[SEED_BYTES], uint32_t domain_separator) {
  size_t i;
  for (i = 0; i < PARAM_ME; i++) {
    poly_p_uniform(vec->entries[i], seed, domain_separator, i, 0);
  }
}

/*************************************************
* Name:        poly_p_vec_m_binomial
*
* Description: Sample a centered binomial polynomial vector of
*              size PARAM_ME with binomial parameter 1 
*              deterministically from a seed.
* 
* Arguments:   - poly_p_vec_m res: output binomial polynomial vector (initialized)
*              - const uint8_t *seed: pointer to byte array containing seed (allocated SEED_BYTES bytes)
*              - uint32_t cnt: repetition domain separator for XOF
*              - uint32_t domain_separator: domain separator for XOF
**************************************************/
void poly_p_vec_m_binomial(poly_p_vec_m res, const uint8_t buf[SEED_BYTES], const uint32_t cnt, const uint32_t domain_separator) {
#if (PARAM_N%64) != 0
#error "PARAM_N must be divisible by 64"
#endif
  uint64_t output[PARAM_ME*PARAM_N*2/64]; // 2 bits per coefficient
  uint64_t coef_lsb[PARAM_N/64];
  uint64_t coef_sign[PARAM_N/64];
  keccak_state state;
  size_t k,l;
  shake256_init(&state);
  shake256_absorb(&state, buf, SEED_BYTES);
  shake256_absorb(&state, (const uint8_t*)&domain_separator, sizeof(uint32_t));
  shake256_absorb(&state, (const uint8_t*)&cnt, sizeof(uint32_t));
  shake256_finalize(&state);
  shake256_squeeze((uint8_t*)output, PARAM_ME*PARAM_N*2/8, &state);
  for (k = 0; k < PARAM_ME; k++) {
    for (l = 0; l < PARAM_N/64; l++) {
      coef_lsb[l]  = output[2*l+2*PARAM_N/64*k] ^ output[2*l+1+2*PARAM_N/64*k];
      coef_sign[l] = output[2*l+2*PARAM_N/64*k] & output[2*l+1+2*PARAM_N/64*k];
    }
    for (l = 0; l < PARAM_N; l++) {
      poly_p_set_coeff(res->entries[k], l, (int32_t)((coef_lsb[l/64] >> (l%64))&1) + (int32_t)(((coef_sign[l/64] >> ((l%64))) << 1)&2) - 1);
      // we have for sign||lsb either 00 (->-1) or 01 (->0) or 10 (->1), so we reconstruct this and subtract one
    }
  }
}