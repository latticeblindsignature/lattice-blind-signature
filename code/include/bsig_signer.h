#ifndef BSIG_SIGNER_H
#define BSIG_SIGNER_H

#include <stdint.h>
#include "params.h"
#include "arith.h"

typedef struct {
  poly_q_mat_d_d R[2][PARAM_K];
  poly_real_mat_2d_2d S;
} sk_t;

typedef struct {
  poly_q_mat_d_d B[PARAM_K];
  uint8_t seed[SEED_BYTES];
} pk_t;

typedef struct {
  poly_q_vec_d v12;
  poly_q_vec_d v2[PARAM_K];
  poly_q_vec_k v3;
} pre_sig_t;

void keys_init(pk_t *pk, sk_t *sk);
void keys_clear(pk_t *pk, sk_t *sk);
void pre_sig_init(pre_sig_t *pre_sig);
void pre_sig_clear(pre_sig_t *pre_sig);

void keygen(pk_t *pk, sk_t *sk);
void tag_gen(poly_q tag, uint8_t state[STATE_BYTES]);
void pre_sign_commitment(pre_sig_t *pre_sig, const sk_t *sk, const pk_t *pk, const poly_q_vec_d cmt, const poly_q tag);
int pre_sig_verify_from_commitment(poly_q_vec_d v11, const poly_q tag, const pre_sig_t *pre_sig, const poly_q_vec_d cmt, const pk_t *pk);

#endif /* BSIG_SIGNER_H */
