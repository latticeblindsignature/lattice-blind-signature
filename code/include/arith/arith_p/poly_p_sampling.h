#ifndef SAMPLING_P_H
#define SAMPLING_P_H

#include <stdint.h>
#include "params.h"
#include "arith.h"

void poly_p_mat_d_m_uniform(poly_p_mat_d_m mat, const uint8_t seed[SEED_BYTES], uint32_t domain_separator);
void poly_p_vec_m_uniform(poly_p_vec_m vec, const uint8_t seed[SEED_BYTES], uint32_t domain_separator);
void poly_p_vec_m_binomial(poly_p_vec_m mat, const uint8_t seed[SEED_BYTES], uint32_t cnt, uint32_t domain_separator);

#endif /* SAMPLING_P_H */
