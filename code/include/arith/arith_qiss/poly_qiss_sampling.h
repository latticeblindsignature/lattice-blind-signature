#ifndef POLY_QISS_SAMPLING_H
#define POLY_QISS_SAMPLING_H

#include <stdint.h>
#include "arith.h"

void vec_qiss_uniform(coeff_qiss out[PARAM_ARP_ISS + 4 + PARAM_N_ISS - 1], const uint8_t *buf, const uint32_t domain_separator, const uint32_t counter, size_t buflen);
void poly_qiss_uniform_but_zero_half(poly_qiss out, const uint8_t seed[SEED_BYTES], uint32_t kappa, uint32_t domain_separator);
void poly_qiss_mat_d_m1_uniform(poly_qiss_mat_d_m1 mat, const uint8_t seed[SEED_BYTES], uint32_t domain_separator);
void poly_qiss_mat_d_m2_d_uniform(poly_qiss_mat_d_m2_d mat, const uint8_t seed[SEED_BYTES], uint32_t domain_separator);
void poly_qiss_mat_256l_m2_d_uniform(poly_qiss_mat_256l_m2_d mat, const uint8_t seed[SEED_BYTES], uint32_t domain_separator);
void poly_qiss_vec_m2_d_uniform(poly_qiss_vec_m2_d vec, const uint8_t seed[SEED_BYTES], uint32_t domain_separator);
void poly_qiss_vec_l_1_uniform(poly_qiss_vec_l vec, poly_qiss pout, const uint8_t *buf, uint32_t domain_separator, size_t buflen);
void poly_qiss_vec_k_uniform(poly_qiss_vec_k vec, const uint8_t *buf, uint32_t domain_separator, uint32_t cnt, size_t buflen);
// TODO streamline the binomial sampling function signatures
void poly_qiss_vec_m1_vec_k_binomial(poly_qiss_vec_m1 res_1, poly_qiss_vec_k res_2[PARAM_DE+1], const uint8_t *buf, uint32_t domain_separator, uint32_t i, size_t inlen);
void poly_qiss_vec_m2_d_binomial(poly_qiss_vec_m2_d res_1, poly_qiss_vec_d res_2, const uint8_t buf[SEED_BYTES], const uint32_t cnt, const uint32_t domain_separator);

void poly_qiss_vec_m1_sample_gaussian_s1(poly_qiss_vec_m1 res);
void poly_qiss_vec_m2_d_sample_gaussian_s2(poly_qiss_vec_m2_d res_1, poly_qiss_vec_d res_2);
void poly_qiss_sample_challenge(poly_qiss out, const uint8_t *buf, const uint32_t domain_separator, const uint32_t counter, size_t buflen);

#endif /* POLY_QISS_SAMPLING_H */
