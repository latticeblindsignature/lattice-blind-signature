#include "bench_bsig_signer.h"

#include "bsig_signer.h"
#include "bsig_user.h"
#include "randombytes.h"
#include "random.h"

double tag_gen_bench(timer* t) {
  double time;
  uint8_t state[STATE_BYTES];
  poly_q tag;
  poly_q_init(tag);
  randombytes(state, STATE_BYTES);
  /* ----------- BEGIN: Code under measurement --------- */
  start_timer(t);
  tag_gen(tag, state);
  time = stop_timer(t);
  /* ----------- END: Code under measurement ----------- */
  poly_q_clear(tag);
  return time;
}

double keygen_bench(timer* t) {
  double time;
  sk_t sk;
  pk_t pk;
  /* ----------- BEGIN: Code under measurement --------- */
  start_timer(t);
  keys_init(&pk, &sk);
  keygen(&pk, &sk);
  time = stop_timer(t);
  /* ----------- END: Code under measurement ----------- */
  keys_clear(&pk, &sk);
  return time;
}

double pre_sign_commitment_bench(timer* t) {
  double time;
  sk_t sk;
  pk_t pk;
  pre_sig_t pre_sig;
  rand_t rand;
  poly_q tag;
  poly_q_vec_d cmt, v11;
  uint8_t state[STATE_BYTES], msg[PARAM_N/8];

  keys_init(&pk, &sk);
  rand_init(&rand);
  poly_q_vec_d_init(cmt);
  poly_q_vec_d_init(v11);
  poly_q_init(tag);

  keygen(&pk, &sk);
  randombytes(msg, PARAM_N/8);
  randombytes(state, STATE_BYTES);
  tag_gen(tag, state);
  commit(&rand, cmt, tag, msg, &pk);
  /* ----------- BEGIN: Code under measurement --------- */
  start_timer(t);
  pre_sig_init(&pre_sig);
  pre_sign_commitment(&pre_sig, &sk, &pk, cmt, tag);
  time = stop_timer(t);
  /* ----------- END: Code under measurement ----------- */
  if (!pre_sig_verify_from_commitment(v11, tag, &pre_sig, cmt, &pk)) {
    printf("FATAL ERROR: benchmarked pre signature is not valid\n");
  }
  pre_sig_clear(&pre_sig);
  keys_clear(&pk, &sk);
  rand_clear(&rand);
  poly_q_vec_d_clear(cmt);
  poly_q_vec_d_clear(v11);
  poly_q_clear(tag);
  return time;
}

double pre_sig_verify_from_commitment_valid_bench(timer* t) {
  double time;
  sk_t sk;
  pk_t pk;
  pre_sig_t pre_sig;
  rand_t rand;
  poly_q tag;
  poly_q_vec_d cmt, v11;
  uint8_t state[STATE_BYTES], msg[PARAM_N/8];

  pre_sig_init(&pre_sig);
  keys_init(&pk, &sk);
  rand_init(&rand);
  poly_q_vec_d_init(cmt);
  poly_q_vec_d_init(v11);
  poly_q_init(tag);

  keygen(&pk, &sk);
  randombytes(msg, PARAM_N/8);
  randombytes(state, STATE_BYTES);
  tag_gen(tag, state);
  commit(&rand, cmt, tag, msg, &pk);
  pre_sign_commitment(&pre_sig, &sk, &pk, cmt, tag);
  /* ----------- BEGIN: Code under measurement --------- */
  start_timer(t);
  int is_valid = pre_sig_verify_from_commitment(v11, tag, &pre_sig, cmt, &pk);
  time = stop_timer(t);
  /* ----------- END: Code under measurement ----------- */
  if (!is_valid) {
    printf("FATAL ERROR: benchmarked pre signature is not valid\n");
  }
  pre_sig_clear(&pre_sig);
  keys_clear(&pk, &sk);
  rand_clear(&rand);
  poly_q_vec_d_clear(cmt);
  poly_q_vec_d_clear(v11);
  poly_q_clear(tag);
  return time;
}

double pre_sig_verify_from_commitment_invalid_bench(timer* t) {
  double time;
  sk_t sk;
  pk_t pk;
  pre_sig_t pre_sig;
  rand_t rand;
  poly_q tag;
  poly_q_vec_d cmt, v11;
  uint8_t state[STATE_BYTES], msg[PARAM_N/8];

  pre_sig_init(&pre_sig);
  keys_init(&pk, &sk);
  rand_init(&rand);
  poly_q_vec_d_init(cmt);
  poly_q_vec_d_init(v11);
  poly_q_init(tag);

  keygen(&pk, &sk);
  randombytes(msg, PARAM_N/8);
  randombytes(state, STATE_BYTES);
  tag_gen(tag, state);
  commit(&rand, cmt, tag, msg, &pk);
  pre_sign_commitment(&pre_sig, &sk, &pk, cmt, tag);
  msg[0] ^= 1;
  commit(&rand, cmt, tag, msg, &pk);
  /* ----------- BEGIN: Code under measurement --------- */
  start_timer(t);
  int is_valid = pre_sig_verify_from_commitment(v11, tag, &pre_sig, cmt, &pk);
  time = stop_timer(t);
  /* ----------- END: Code under measurement ----------- */
  if (is_valid) {
    printf("FATAL ERROR: benchmarked pre signature is valid\n");
  }
  pre_sig_clear(&pre_sig);
  keys_clear(&pk, &sk);
  rand_clear(&rand);
  poly_q_vec_d_clear(cmt);
  poly_q_vec_d_clear(v11);
  poly_q_clear(tag);
  return time;
}

