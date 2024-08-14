#include "poly_qiss_sampling.h"
#include "poly_q_sampling.h"
#include "fips202.h"
#include "random.h"

/*************************************************
* Name:        poly_qiss_uniform
*
* Description: Sample a uniformly random polynomial modulo
*              PARAM_Q_ISS deterministically from a seed.
* 
* Arguments:   - poly_qiss pout: output uniform polynomial (initialized)
*              - const uint8_t *seed: pointer to byte array containing seed (allocated SEED_BYTES bytes)
*              - uint32_t domain_separator: domain separator for XOF
*              - size_t i: index domain separator for XOF
*              - size_t j: index domain separator for XOF
**************************************************/
static void poly_qiss_uniform(poly_qiss pout, const uint8_t seed[SEED_BYTES], uint32_t domain_separator, size_t i, size_t j) {
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
  while (cnt < PARAM_N_ISS) {
#if PARAM_Q_ISS_BITLEN > 60
#error "PARAM_Q_ISS_BITLEN too big for uniform sampling."
#else
    // idea: take 15 bytes, ignore MSBs
    if (bytecnt < 15) {
      for (k = 0; k < bytecnt; k++) {
        output[k] = output[off++];
      }
      shake128_squeezeblocks(&output[bytecnt], 1, &state);
      off = 0;
      bytecnt += SHAKE128_RATE;
    }
    uint64_t tmp8byte = output[off] | ((uint64_t)output[off+1] << 8) | ((uint64_t)output[off+2] << 16) | ((uint64_t)output[off+3] << 24) | ((uint64_t)output[off+4] << 32) | ((uint64_t)output[off+5] << 40) | ((uint64_t)output[off+6] << 48) | ((uint64_t)output[off+7] << 56);
    uint64_t tmp = tmp8byte & ((1UL<<PARAM_Q_ISS_BITLEN)-1);
    if (tmp < PARAM_Q_ISS) {
      poly_qshow_set_coeff(pout, cnt++, tmp);
    }
    if (cnt >= PARAM_N_ISS) {
      break;
    }
    tmp8byte >>= PARAM_Q_ISS_BITLEN;
    tmp8byte |= (output[off+8] | ((uint64_t)output[off+9] << 8) | ((uint64_t)output[off+10] << 16) | ((uint64_t)output[off+11] << 24) | ((uint64_t)output[off+12] << 32) | ((uint64_t)output[off+13] << 40) | ((uint64_t)output[off+14] << 48)) << 7;
    tmp = tmp8byte & ((1UL<<PARAM_Q_ISS_BITLEN)-1);
    if (tmp < PARAM_Q_ISS) {
      poly_qshow_set_coeff(pout, cnt++, tmp);
    }
    off += 15;
    bytecnt -= 15;
#if PARAM_Q_ISS_BITLEN < 57
#warning "PARAM_Q_ISS_BITLEN maybe unsuitable for efficient uniform sampling."
#endif
#endif
  }
}

/*************************************************
* Name:        vec_qiss_uniform
*
* Description: Sample a uniformly random vector of integer modulo
*              PARAM_Q_ISS of size PARAM_ARP_ISS + 1 deterministically 
*              from an input buffer.
* 
* Arguments:   - coeff_qiss *out: output uniform vector (allocated PARAM_ARP_ISS+1 coeff_qiss)
*              - const uint8_t *buf: pointer to byte buffer containing the XOF input
*              - const uint32_t domain_separator: domain separator for XOF
*              - const uint32_t counter: rejection domain separator for XOF
*              - size_t buflen: length of the input buffer
**************************************************/
void vec_qiss_uniform(coeff_qiss out[PARAM_ARP_ISS + 4 + PARAM_N_ISS - 1], const uint8_t *buf, const uint32_t domain_separator, const uint32_t counter, size_t buflen) {
  uint8_t output[SHAKE128_RATE * 2];
  keccak_state state;
  size_t k,cnt,off,bytecnt;
  shake128_init(&state);
  shake128_absorb(&state, buf, buflen);
  shake128_absorb(&state, (const uint8_t*)&domain_separator, sizeof(uint32_t));
  shake128_absorb(&state, (const uint8_t*)&counter, sizeof(uint32_t));
  shake128_finalize(&state);
  shake128_squeezeblocks(output, 2, &state);
  bytecnt = 2*SHAKE128_RATE;

  cnt = 0;
  off = 0;
  while (cnt < PARAM_ARP_ISS + 4 + PARAM_N_ISS - 1) {
#if PARAM_Q_ISS_BITLEN > 60
#error "PARAM_Q_ISS_BITLEN too big for uniform sampling."
#else
    // idea: take 15 bytes, ignore MSBs
    if (bytecnt < 15) {
      for (k = 0; k < bytecnt; k++) {
        output[k] = output[off++];
      }
      shake128_squeezeblocks(&output[bytecnt], 1, &state);
      off = 0;
      bytecnt += SHAKE128_RATE;
    }
    uint64_t tmp8byte = output[off] | ((uint64_t)output[off+1] << 8) | ((uint64_t)output[off+2] << 16) | ((uint64_t)output[off+3] << 24) | ((uint64_t)output[off+4] << 32) | ((uint64_t)output[off+5] << 40) | ((uint64_t)output[off+6] << 48) | ((uint64_t)output[off+7] << 56);
    uint64_t tmp = tmp8byte & ((1UL<<PARAM_Q_ISS_BITLEN)-1);
    if (tmp < PARAM_Q_ISS) {
      out[cnt++] = tmp;
    }
    if (cnt >= PARAM_ARP_ISS + 4 + PARAM_N_ISS - 1) {
      break;
    }
    tmp8byte >>= PARAM_Q_ISS_BITLEN;
    tmp8byte |= (output[off+8] | ((uint64_t)output[off+9] << 8) | ((uint64_t)output[off+10] << 16) | ((uint64_t)output[off+11] << 24) | ((uint64_t)output[off+12] << 32) | ((uint64_t)output[off+13] << 40) | ((uint64_t)output[off+14] << 48)) << 7;
    tmp = tmp8byte & ((1UL<<PARAM_Q_ISS_BITLEN)-1);
    if (tmp < PARAM_Q_ISS) {
      out[cnt++] = tmp;
    }
    off += 15;
    bytecnt -= 15;
#if PARAM_Q_ISS_BITLEN < 57
#warning "PARAM_Q_ISS_BITLEN maybe unsuitable for efficient uniform sampling."
#endif
#endif
  }
}

/*************************************************
* Name:        poly_qiss_uniform_but_zero_half
*
* Description: Sample a uniformly random polynomial modulo
*              PARAM_Q_ISS with zero constant coefficient,
*              and zero n/2-th coefficient deterministically 
*              from a seed.
* 
* Arguments:   - poly_qiss out: output uniform polynomial (initialized)
*              - const uint8_t *seed: pointer to byte array containing seed (allocated SEED_BYTES bytes)
*              - uint32_t kappa: rejection domain separator for XOF
*              - uint32_t domain_separator: domain separator for XOF
**************************************************/
void poly_qiss_uniform_but_zero_half(poly_qiss out, const uint8_t seed[SEED_BYTES], uint32_t kappa, uint32_t domain_separator) {
  // TODO implement dedicated function?
  poly_qiss_uniform(out, seed, domain_separator, kappa, 0);
  poly_qiss_set_coeff(out, 0, 0);
  poly_qiss_set_coeff(out, PARAM_N_ISS/2, 0);
}

/*************************************************
* Name:        poly_qiss_mat_d_m1_uniform
*
* Description: Sample a uniformly random polynomial matrix of
*              size PARAM_D_ISS x PARAM_M1_ISS modulo PARAM_Q_ISS 
*              deterministically from a seed.
* 
* Arguments:   - poly_qiss_mat_d_m1 mat: output uniform polynomial matrix (initialized)
*              - const uint8_t *seed: pointer to byte array containing seed (allocated SEED_BYTES bytes)
*              - uint32_t domain_separator: domain separator for XOF
**************************************************/
void poly_qiss_mat_d_m1_uniform(poly_qiss_mat_d_m1 mat, const uint8_t seed[SEED_BYTES], uint32_t domain_separator) {
  size_t i,j;
  for (i = 0; i < PARAM_D_ISS; i++) {
    for (j = 0; j < PARAM_M1_ISS; j++) {
      poly_qiss_uniform(mat->rows[i]->entries[j], seed, domain_separator, i, j);
    }
  }
}

/*************************************************
* Name:        poly_qiss_mat_d_m2_d_uniform
*
* Description: Sample a uniformly random polynomial matrix of
*              size PARAM_D_ISS x PARAM_M2_D_ISS modulo PARAM_Q_ISS 
*              deterministically from a seed.
* 
* Arguments:   - poly_qiss_mat_d_m2_d mat: output uniform polynomial matrix (initialized)
*              - const uint8_t *seed: pointer to byte array containing seed (allocated SEED_BYTES bytes)
*              - uint32_t domain_separator: domain separator for XOF
**************************************************/
void poly_qiss_mat_d_m2_d_uniform(poly_qiss_mat_d_m2_d mat, const uint8_t seed[SEED_BYTES], uint32_t domain_separator) {
  size_t i,j;
  for (i = 0; i < PARAM_D_ISS; i++) {
    for (j = 0; j < PARAM_M2_D_ISS; j++) {
      poly_qiss_uniform(mat->rows[i]->entries[j], seed, domain_separator, i, j);
    }
  }
}

/*************************************************
* Name:        poly_qiss_mat_256l_m2_d_uniform
*
* Description: Sample a uniformly random polynomial matrix of
*              size PARAM_ARP_DIV_N_L_ISS x PARAM_M2_D_ISS modulo PARAM_Q_ISS 
*              deterministically from a seed.
* 
* Arguments:   - poly_qiss_mat_256l_m2_d mat: output uniform polynomial matrix (initialized)
*              - const uint8_t *seed: pointer to byte array containing seed (allocated SEED_BYTES bytes)
*              - uint32_t domain_separator: domain separator for XOF
**************************************************/
void poly_qiss_mat_256l_m2_d_uniform(poly_qiss_mat_256l_m2_d mat, const uint8_t seed[SEED_BYTES], uint32_t domain_separator) {
  size_t i,j;
  for (i = 0; i < PARAM_ARP_DIV_N_L_ISS; i++) {
    for (j = 0; j < PARAM_M2_D_ISS; j++) {
      poly_qiss_uniform(mat->rows[i]->entries[j], seed, domain_separator, i, j);
    }
  }
}

/*************************************************
* Name:        poly_qiss_vec_m2_d_uniform
*
* Description: Sample a uniformly random polynomial vector of
*              size PARAM_M2_D_ISS modulo PARAM_Q_ISS 
*              deterministically from a seed.
* 
* Arguments:   - poly_qiss_vec_m2 vec: output uniform polynomial matrix (initialized)
*              - const uint8_t *seed: pointer to byte array containing seed (allocated SEED_BYTES bytes)
*              - uint32_t domain_separator: domain separator for XOF
**************************************************/
void poly_qiss_vec_m2_d_uniform(poly_qiss_vec_m2_d vec, const uint8_t seed[SEED_BYTES], uint32_t domain_separator) {
  size_t i;
  for (i = 0; i < PARAM_M2_D_ISS; i++) {
    poly_qiss_uniform(vec->entries[i], seed, domain_separator, i, 0);
  }
}

/*************************************************
* Name:        poly_qiss_vec_l_1_uniform
*
* Description: Sample a uniformly random polynomial vector of
*              size PARAM_L_ISS modulo PARAM_Q_ISS and extra polynomial
*              deterministically from an input buffer.
* 
* Arguments:   - poly_qiss_vec_l vec: output uniform polynomial vector (initialized)
*              - poly_qiss pout: output uniform polynomial (initialized)
*              - const uint8_t *buf: pointer to byte buffer containing the XOF input
*              - uint32_t domain_separator: domain separator for XOF
*              - size_t buflen: length of the input buffer
**************************************************/
void poly_qiss_vec_l_1_uniform(poly_qiss_vec_l vec, poly_qiss pout, const uint8_t *buf, uint32_t domain_separator, size_t buflen) {
  size_t i;
  uint8_t lazybuf[SEED_BYTES];
  sha3_256(lazybuf, buf, buflen); // TODO remove this and use buf directly for sampling input
  for (i = 0; i < PARAM_L_ISS; i++) {
    poly_qiss_uniform(vec->entries[i], lazybuf, domain_separator, i, 0);
  }
  poly_qiss_uniform(pout, lazybuf, domain_separator, PARAM_L_ISS, 0);
}

/*************************************************
* Name:        poly_qiss_vec_k_uniform
*
* Description: Sample a uniformly random polynomial vector of
*              size PARAM_K_ISS modulo PARAM_Q_ISS 
*              deterministically from an input buffer.
* 
* Arguments:   - poly_qiss_vec_k vec: output uniform polynomial matrix (initialized)
*              - const uint8_t *buf: pointer to byte buffer containing the XOF input
*              - uint32_t domain_separator: domain separator for XOF
*              - uint32_t cnt: rejection domain separator for XOF
*              - size_t buflen: length of the input buffer
**************************************************/
void poly_qiss_vec_k_uniform(poly_qiss_vec_k vec, const uint8_t *buf, uint32_t domain_separator, uint32_t cnt, size_t buflen) {
  size_t i;
  uint8_t lazybuf[SEED_BYTES];
  sha3_256(lazybuf, buf, buflen); // TODO remove this and use buf directly for sampling input
  for (i = 0; i < PARAM_K_ISS; i++) {
    poly_qiss_uniform(vec->entries[i], lazybuf, domain_separator, i, cnt);
  }
}

// TODO streamline the binomial sampling function signatures
/*************************************************
* Name:        poly_qiss_vec_m1_vec_k_binomial
*
* Description: Sample a centered binomial polynomial vector of
*              size PARAM_M1_ISS + PARAM_K_ISS*(PARAM_DE+1) with binomial parameter 1 
*              deterministically from a buffer.
* 
* Arguments:   - poly_qiss_vec_m1 res: output binomial polynomial vector (initialized)
*              - const uint8_t *buf: pointer to byte buffer containing the XOF input
*              - uint32_t domain_separator: domain separator for XOF
*              - size_t i: index domain separator for XOF
*              - size_t j: index domain separator for XOF
*              - size_t inlen: length of the input buffer
**************************************************/
void poly_qiss_vec_m1_vec_k_binomial(poly_qiss_vec_m1 res_1, poly_qiss_vec_k res_2[PARAM_DE+1], const uint8_t *buf, uint32_t domain_separator, uint32_t i, size_t inlen) {
#if (PARAM_N_ISS%64) != 0
#error "PARAM_N_ISS must be divisible by 64"
#endif
  uint64_t output[(PARAM_M1_ISS+PARAM_K_ISS*(PARAM_DE+1))*2*PARAM_N_ISS/64]; // 2 bits per coefficient
  uint64_t coef_lsb[PARAM_N_ISS/64];
  uint64_t coef_sign[PARAM_N_ISS/64];
  keccak_state state;
  size_t j,k,l;
  shake256_init(&state);
  shake256_absorb(&state, buf, inlen);
  shake256_absorb(&state, (const uint8_t*)&domain_separator, sizeof(uint32_t));
  shake256_absorb(&state, (const uint8_t*) &i, sizeof(uint32_t));
  shake256_finalize(&state);
  shake256_squeeze((uint8_t*)output, (PARAM_M1_ISS+PARAM_K_ISS*(PARAM_DE+1))*2*PARAM_N_ISS/8, &state);
  for (k = 0; k < PARAM_M1_ISS; k++) {
    for (l = 0; l < PARAM_N_ISS/64; l++) {
      coef_lsb[l]  = output[2*l+2*PARAM_N_ISS/64*k] ^ output[2*l+1+2*PARAM_N_ISS/64*k];
      coef_sign[l] = output[2*l+2*PARAM_N_ISS/64*k] & output[2*l+1+2*PARAM_N_ISS/64*k];
    }
    for (l = 0; l < PARAM_N_ISS; l++) {
      poly_qiss_set_coeff(res_1->entries[k], l, (int32_t)((coef_lsb[l/64] >> (l%64))&1) + (int32_t)(((coef_sign[l/64] >> ((l%64))) << 1)&2) - 1);
      // we have for sign||lsb either 00 (->-1) or 01 (->0) or 10 (->1), so we reconstruct this and subtract one
    }
  }

  for (j = 0; j < PARAM_DE+1; j++) {
    for (k = 0; k < PARAM_K_ISS; k++) {
      for (l = 0; l < PARAM_N_ISS/64; l++) {
        coef_lsb[l]  = output[2*PARAM_N_ISS/64*(PARAM_M1_ISS + PARAM_K_ISS*j) + 2*l+2*PARAM_N_ISS/64*k] ^ output[2*PARAM_N_ISS/64*(PARAM_M1_ISS + PARAM_K_ISS*j) + 2*l+1+2*PARAM_N_ISS/64*k];
        coef_sign[l] = output[2*PARAM_N_ISS/64*(PARAM_M1_ISS + PARAM_K_ISS*j) + 2*l+2*PARAM_N_ISS/64*k] & output[2*PARAM_N_ISS/64*(PARAM_M1_ISS + PARAM_K_ISS*j) + 2*l+1+2*PARAM_N_ISS/64*k];
      }
      for (l = 0; l < PARAM_N_ISS; l++) {
        poly_qiss_set_coeff(res_2[j]->entries[k], l, (int32_t)((coef_lsb[l/64] >> (l%64))&1) + (int32_t)(((coef_sign[l/64] >> ((l%64))) << 1)&2) - 1);
        // we have for sign||lsb either 00 (->-1) or 01 (->0) or 10 (->1), so we reconstruct this and subtract one
      }
    }
  }
}

/*************************************************
* Name:        poly_qiss_vec_m2_d_binomial
*
* Description: Sample a centered binomial polynomial vector of
*              size PARAM_M2_D_ISS with binomial parameter 1 
*              deterministically from a buffer.
* 
* Arguments:   - poly_qiss_vec_m2_d res_1: first output binomial polynomial vector (initialized)
*              - poly_qiss_vec_d res_2: second output binomial polynomial vector (initialized)
*              - const uint8_t *buf: pointer to byte buffer containing the XOF input (allocated SEED_BYTES bytes)
*              - const uint32_t cnt: rejection domain separator for XOF
*              - const uint32_t domain_separator: domain separator for XOF
**************************************************/
void poly_qiss_vec_m2_d_binomial(poly_qiss_vec_m2_d res_1, poly_qiss_vec_d res_2, const uint8_t buf[SEED_BYTES], const uint32_t cnt, const uint32_t domain_separator) {
#if (PARAM_N_ISS%64) != 0
#error "PARAM_N_ISS must be divisible by 64"
#endif
  uint64_t output[2*(PARAM_M2_D_ISS + PARAM_D_ISS)*PARAM_N_ISS/64]; // 2 bits per coefficient
  uint64_t coef_lsb[PARAM_N_ISS/64];
  uint64_t coef_sign[PARAM_N_ISS/64];
  keccak_state state;
  size_t k,l;
  shake256_init(&state);
  shake256_absorb(&state, buf, SEED_BYTES);
  shake256_absorb(&state, (const uint8_t*)&domain_separator, sizeof(uint32_t));
  shake256_absorb(&state, (const uint8_t*)&cnt, sizeof(uint32_t));
  shake256_finalize(&state);
  shake256_squeeze((uint8_t*)output, 2*(PARAM_M2_D_ISS + PARAM_D_ISS)*PARAM_N_ISS/8, &state);
  for (k = 0; k < PARAM_M2_D_ISS; k++) {
    for (l = 0; l < PARAM_N_ISS/64; l++) {
      coef_lsb[l]  = output[2*l+2*PARAM_N_ISS/64*k] ^ output[2*l+1+2*PARAM_N_ISS/64*k];
      coef_sign[l] = output[2*l+2*PARAM_N_ISS/64*k] & output[2*l+1+2*PARAM_N_ISS/64*k];
    }
    for (l = 0; l < PARAM_N_ISS; l++) {
      poly_qiss_set_coeff(res_1->entries[k], l, (int32_t)((coef_lsb[l/64] >> (l%64))&1) + (int32_t)(((coef_sign[l/64] >> ((l%64))) << 1)&2) - 1);
      // we have for sign||lsb either 00 (->-1) or 01 (->0) or 10 (->1), so we reconstruct this and subtract one
    }
  }
  for (k = 0; k < PARAM_D_ISS; k++) {
    for (l = 0; l < PARAM_N_ISS/64; l++) {
      coef_lsb[l]  = output[2*PARAM_N_ISS/64*PARAM_M2_D_ISS  +  2*l+2*PARAM_N_ISS/64*k] ^ output[2*PARAM_N_ISS/64*PARAM_M2_D_ISS  +  2*l+1+2*PARAM_N_ISS/64*k];
      coef_sign[l] = output[2*PARAM_N_ISS/64*PARAM_M2_D_ISS  +  2*l+2*PARAM_N_ISS/64*k] & output[2*PARAM_N_ISS/64*PARAM_M2_D_ISS  +  2*l+1+2*PARAM_N_ISS/64*k];
    }
    for (l = 0; l < PARAM_N_ISS; l++) {
      poly_qiss_set_coeff(res_2->entries[k], l, (int32_t)((coef_lsb[l/64] >> (l%64))&1) + (int32_t)(((coef_sign[l/64] >> ((l%64))) << 1)&2) - 1);
      // we have for sign||lsb either 00 (->-1) or 01 (->0) or 10 (->1), so we reconstruct this and subtract one
    }
  }
}

/*************************************************
* Name:        poly_qiss_vec_m1_sample_gaussian_s1
*
* Description: Sample a polynomial vector with PARAM_M1_ISS entries
*              from the centered spherical Gaussian with parameter
*              PARAM_S1_ISS
* 
* Arguments:   - poly_qiss_vec_m1 res: the polynomial to host the Gaussian sample
**************************************************/
void poly_qiss_vec_m1_sample_gaussian_s1(poly_qiss_vec_m1 res) {
  for (size_t i = 0; i < PARAM_M1_ISS; i++) {
    for (size_t j = 0; j < PARAM_N_ISS; j++) {
      poly_qiss_set_coeff(res->entries[i], j, SampleZ(0, PARAM_S1_ISS));
    }
  }
}

/*************************************************
* Name:        poly_qiss_vec_m2_d_sample_gaussian_s2
*
* Description: Sample a polynomial vector with PARAM_M2_ISS entries
*              from the centered spherical Gaussian with parameter
*              PARAM_S2_ISS
* 
* Arguments:   - poly_qiss_vec_m2_d res_1: first polynomial to host the Gaussian sample
*              - poly_qiss_vec_d res_2: second polynomial to host the Gaussian sample
**************************************************/
void poly_qiss_vec_m2_d_sample_gaussian_s2(poly_qiss_vec_m2_d res_1, poly_qiss_vec_d res_2) {
  for (size_t i = 0; i < PARAM_M2_D_ISS; i++) {
    for (size_t j = 0; j < PARAM_N_ISS; j++) {
      poly_qiss_set_coeff(res_1->entries[i], j, SampleZ(0, PARAM_S2_ISS));
    }
  }
  for (size_t i = 0; i < PARAM_D_ISS; i++) {
    for (size_t j = 0; j < PARAM_N_ISS; j++) {
      poly_qiss_set_coeff(res_2->entries[i], j, SampleZ(0, PARAM_S2_ISS));
    }
  }
}

/*************************************************
* Name:        poly_qiss_sample_challenge
*
* Description: Sample a uniformly random polynomial with 
*              coefficients between [-PARAM_RHO_ISS, PARAM_RHO_ISS]
*              and self-adjoint deterministically from an input buffer.
* 
* Arguments:   - poly_qiss out: output uniform bounded self-adjoint polynomial (initialized)
*              - const uint8_t *buf: pointer to byte buffer containing the XOF input
*              - const uint32_t domain_separator: domain separator for XOF
*              - const uint32_t counter: rejection domain separator for XOF
*              - size_t buflen: length of the input buffer
**************************************************/
void poly_qiss_sample_challenge(poly_qiss out, const uint8_t *buf, const uint32_t domain_separator, const uint32_t counter, size_t buflen) {
  uint8_t outbuf[SHAKE256_RATE];
  size_t outcnt = 0, bytecnt, pos = 0;
  keccak_state state;
  shake256_init(&state);
  shake256_absorb(&state, buf, buflen);
  shake256_absorb(&state, (const uint8_t*)&domain_separator, sizeof(uint32_t));
  shake256_absorb(&state, (const uint8_t*)&counter, sizeof(uint32_t));
  shake256_finalize(&state);
  shake256_squeezeblocks(outbuf, 1, &state);
  bytecnt = SHAKE256_RATE;

#if PARAM_RHO_ISS != 8
#error "poly_qiss_sample_challenge is implemented specifically for PARAM_RHO_ISS = 8"
#endif
  // the idea: for each coefficient sample a byte and reject the byte iff it is 255
  // if accepted, the coefficient is (byte%17) - 8
  while (outcnt < PARAM_N_ISS/2) {
    if (bytecnt == 0) {
      shake256_squeezeblocks(outbuf, 1, &state);
      bytecnt = SHAKE256_RATE;
      pos = 0;
    }
    if (outbuf[pos] < 255) {
      coeff_qiss tmp = ((coeff_qiss)(outbuf[pos]%17)) - 8;
      poly_qiss_set_coeff(out, outcnt, tmp);
      if (outcnt > 0) {
        poly_qiss_set_coeff(out, PARAM_N_ISS-outcnt, -tmp);
      }
      outcnt += 1;
    }
    bytecnt -= 1;
    pos += 1;
  }
}
