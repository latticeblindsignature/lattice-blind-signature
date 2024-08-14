#include "bench_bsig_user.h"

#include "bsig_user.h"
#include "bsig_signer.h"
#include "randombytes.h"
#include "random.h"

double tag_verify_bench(timer* t) {
  double time;
  uint8_t state[STATE_BYTES];
  poly_q tag;
  poly_q_init(tag);
  randombytes(state, STATE_BYTES);
  tag_gen(tag, state);
  /* ----------- BEGIN: Code under measurement --------- */
  start_timer(t);
  int64_t weight = poly_q_weight(tag);
  if (weight != PARAM_W) {
    printf("FATAL ERROR: benchmarked tag is not valid\n");
  }
  time = stop_timer(t);
  /* ----------- END: Code under measurement ----------- */
  poly_q_clear(tag);
  return time;
}

double commit_bench(timer* t) {
  double time;
  uint8_t state[STATE_BYTES], msg[PARAM_N/8];
  sk_t sk;
  pk_t pk;
  rand_t rand;
  poly_q tag;
  poly_q_vec_d cmt;

  keys_init(&pk, &sk);
  poly_q_init(tag);
  poly_q_vec_d_init(cmt);
  randombytes(state, STATE_BYTES);
  randombytes(msg, PARAM_N/8);
  tag_gen(tag, state);
  keygen(&pk, &sk);
  /* ----------- BEGIN: Code under measurement --------- */
  start_timer(t);
  rand_init(&rand);
  commit(&rand, cmt, tag, msg, &pk);
  time = stop_timer(t);
  /* ----------- END: Code under measurement ----------- */
  keys_clear(&pk, &sk);
  rand_clear(&rand);
  poly_q_clear(tag);
  poly_q_vec_d_clear(cmt);
  return time;
}

double encrypt_bench(timer* t) {
  double time;
  ct_t ct;
  sk_t sk;
  pk_t pk;
  poly_p_vec_m r_e;
  uint8_t msg[PARAM_N/8];

  keys_init(&pk, &sk);
  keygen(&pk, &sk);
  randombytes(msg, PARAM_N/8);
  /* ----------- BEGIN: Code under measurement --------- */
  start_timer(t);
  ct_init(&ct);
  poly_p_vec_m_init(r_e);
  encrypt(&ct, r_e, msg, &pk);
  time = stop_timer(t);
  /* ----------- END: Code under measurement ----------- */
  keys_clear(&pk, &sk);
  ct_clear(&ct);
  poly_p_vec_m_clear(r_e);
  return time;
}

double embed_1_bench(timer *t) {
  double time;
  uint8_t state[STATE_BYTES], msg[PARAM_N/8];
  size_t i,j;
  sk_t sk;
  pk_t pk;
  rand_t rand;
  ct_t ct;
  poly_q tag;
  poly_q_vec_d cmt;
  poly_p_vec_m r_e;
  poly_qiss_mat_k_k A_embed[PARAM_D][PARAM_D], B_embed[PARAM_D][PARAM_K*(PARAM_D+1)], D_embed[PARAM_D], Ae_be_embed[PARAM_DE+1][PARAM_ME];
  poly_qiss_vec_k cmt_embed[PARAM_D], ct_embed[PARAM_DE+1];
  poly_qiss_vec_m1 s1;

  keys_init(&pk, &sk);
  rand_init(&rand);
  ct_init(&ct);
  poly_q_init(tag);
  poly_q_vec_d_init(cmt);
  poly_p_vec_m_init(r_e);

  randombytes(state, STATE_BYTES);
  keygen(&pk, &sk);
  tag_gen(tag, state);
  randombytes(msg, PARAM_N/8);
  commit(&rand, cmt, tag, msg, &pk);
  encrypt(&ct, r_e, msg, &pk);

  /* ----------- BEGIN: Code under measurement --------- */
  start_timer(t);
  // init
  for (i = 0; i < PARAM_D; i++) {
    for (j = 0; j < PARAM_D; j++) {
      poly_qiss_mat_k_k_init(A_embed[i][j]);
    }
    for (j = 0; j < PARAM_K*(PARAM_D+1); j++) {
      poly_qiss_mat_k_k_init(B_embed[i][j]);
    }
    poly_qiss_mat_k_k_init(D_embed[i]);
    poly_qiss_vec_k_init(cmt_embed[i]);
  }
  for (i = 0; i < PARAM_DE + 1; i++) {
    for (j = 0; j < PARAM_ME; j++) {
      poly_qiss_mat_k_k_init(Ae_be_embed[i][j]);
    }
    poly_qiss_vec_k_init(ct_embed[i]);
  }
  poly_qiss_vec_m1_init(s1);
  // embed
  embed_1(A_embed, B_embed, D_embed, Ae_be_embed, cmt_embed, ct_embed, s1, &pk, cmt, &ct, tag, &rand, r_e, msg);
  time = stop_timer(t);
  /* ----------- END: Code under measurement ----------- */

  keys_clear(&pk, &sk);
  rand_clear(&rand);
  ct_clear(&ct);
  poly_q_clear(tag);
  poly_q_vec_d_clear(cmt);
  poly_p_vec_m_clear(r_e);
  for (i = 0; i < PARAM_D; i++) {
    for (j = 0; j < PARAM_D; j++) {
      poly_qiss_mat_k_k_clear(A_embed[i][j]);
    }
    for (j = 0; j < PARAM_K*(PARAM_D+1); j++) {
      poly_qiss_mat_k_k_clear(B_embed[i][j]);
    }
    poly_qiss_mat_k_k_clear(D_embed[i]);
    poly_qiss_vec_k_clear(cmt_embed[i]);
  }
  for (i = 0; i < PARAM_DE + 1; i++) {
    for (j = 0; j < PARAM_ME; j++) {
      poly_qiss_mat_k_k_clear(Ae_be_embed[i][j]);
    }
    poly_qiss_vec_k_clear(ct_embed[i]);
  }
  poly_qiss_vec_m1_clear(s1);
  return time;
}

double embed_1_verifier_bench(timer *t) {
  double time;
  uint8_t state[STATE_BYTES], msg[PARAM_N/8];
  size_t i,j;
  sk_t sk;
  pk_t pk;
  rand_t rand;
  ct_t ct;
  poly_q tag;
  poly_q_vec_d cmt;
  poly_p_vec_m r_e;
  poly_qiss_mat_k_k A_embed[PARAM_D][PARAM_D], B_embed[PARAM_D][PARAM_K*(PARAM_D+1)], D_embed[PARAM_D], Ae_be_embed[PARAM_DE+1][PARAM_ME];
  poly_qiss_vec_k cmt_embed[PARAM_D], ct_embed[PARAM_DE+1];

  keys_init(&pk, &sk);
  rand_init(&rand);
  ct_init(&ct);
  poly_q_init(tag);
  poly_q_vec_d_init(cmt);
  poly_p_vec_m_init(r_e);

  randombytes(state, STATE_BYTES);
  keygen(&pk, &sk);
  tag_gen(tag, state);
  randombytes(msg, PARAM_N/8);
  commit(&rand, cmt, tag, msg, &pk);
  encrypt(&ct, r_e, msg, &pk);

  /* ----------- BEGIN: Code under measurement --------- */
  start_timer(t);
  // init
  for (i = 0; i < PARAM_D; i++) {
    for (j = 0; j < PARAM_D; j++) {
      poly_qiss_mat_k_k_init(A_embed[i][j]);
    }
    for (j = 0; j < PARAM_K*(PARAM_D+1); j++) {
      poly_qiss_mat_k_k_init(B_embed[i][j]);
    }
    poly_qiss_mat_k_k_init(D_embed[i]);
    poly_qiss_vec_k_init(cmt_embed[i]);
  }
  for (i = 0; i < PARAM_DE + 1; i++) {
    for (j = 0; j < PARAM_ME; j++) {
      poly_qiss_mat_k_k_init(Ae_be_embed[i][j]);
    }
    poly_qiss_vec_k_init(ct_embed[i]);
  }
  // embed
  embed_1_verifier(A_embed, B_embed, D_embed, Ae_be_embed, cmt_embed, ct_embed, &pk, cmt, &ct, tag);
  time = stop_timer(t);
  /* ----------- END: Code under measurement ----------- */

  keys_clear(&pk, &sk);
  rand_clear(&rand);
  ct_clear(&ct);
  poly_q_clear(tag);
  poly_q_vec_d_clear(cmt);
  poly_p_vec_m_clear(r_e);
  for (i = 0; i < PARAM_D; i++) {
    for (j = 0; j < PARAM_D; j++) {
      poly_qiss_mat_k_k_clear(A_embed[i][j]);
    }
    for (j = 0; j < PARAM_K*(PARAM_D+1); j++) {
      poly_qiss_mat_k_k_clear(B_embed[i][j]);
    }
    poly_qiss_mat_k_k_clear(D_embed[i]);
    poly_qiss_vec_k_clear(cmt_embed[i]);
  }
  for (i = 0; i < PARAM_DE + 1; i++) {
    for (j = 0; j < PARAM_ME; j++) {
      poly_qiss_mat_k_k_clear(Ae_be_embed[i][j]);
    }
    poly_qiss_vec_k_clear(ct_embed[i]);
  }
  return time;
}

double complete_decompose_bench(timer* t) {
  double time;
  sk_t sk;
  pk_t pk;
  pre_sig_t pre_sig;
  rand_t rand;
  bsig_t bsig;
  poly_q tag;
  poly_q_vec_d cmt, v11, w1H[2], w2H[PARAM_K];
  poly_q_vec_k w3H;
  uint8_t state[STATE_BYTES], msg[PARAM_N/8];

  keys_init(&pk, &sk);
  pre_sig_init(&pre_sig);
  rand_init(&rand);
  poly_q_vec_d_init(cmt);
  poly_q_vec_d_init(v11);
  poly_q_init(tag);

  randombytes(state, STATE_BYTES);
  keygen(&pk, &sk);
  tag_gen(tag, state);
  randombytes(msg, PARAM_N/8);
  commit(&rand, cmt, tag, msg, &pk);
  pre_sign_commitment(&pre_sig, &sk, &pk, cmt, tag);
  int is_valid = pre_sig_verify_from_commitment(v11, tag, &pre_sig, cmt, &pk);
  if (!is_valid) {
    printf("FATAL ERROR: benchmarked pre signature is not valid\n");
  }
  /* ----------- BEGIN: Code under measurement --------- */
  start_timer(t);
  bsig_init(&bsig);
  poly_q_vec_d_init(w1H[0]);
  poly_q_vec_d_init(w1H[1]);
  for (int k = 0; k < PARAM_K; k++) {
    poly_q_vec_d_init(w2H[k]);
  }
  poly_q_vec_k_init(w3H);
  complete_decompose(&bsig, w1H, w2H, w3H, v11, &pre_sig, &rand);
  time = stop_timer(t);
  /* ----------- END: Code under measurement ----------- */
  pre_sig_clear(&pre_sig);
  keys_clear(&pk, &sk);
  rand_clear(&rand);
  bsig_clear(&bsig);
  poly_q_vec_d_clear(cmt);
  poly_q_vec_d_clear(v11);
  poly_q_vec_d_clear(w1H[0]);
  poly_q_vec_d_clear(w1H[1]);
  for (int k = 0; k < PARAM_K; k++) {
    poly_q_vec_d_clear(w2H[k]);
  }
  poly_q_vec_k_clear(w3H);
  poly_q_clear(tag);
  return time;
}

double embed_2_bench(timer *t) {
  double time;
  uint8_t state[STATE_BYTES], msg[PARAM_N/8];
  size_t i,j;
  sk_t sk;
  pk_t pk;
  pre_sig_t pre_sig;
  rand_t rand;
  bsig_t bsig;
  poly_q tag;
  poly_q_vec_d cmt, v11, w1H[2], w2H[PARAM_K];
  poly_q_vec_k w3H;
  poly_qshow_mat_k_k A_embed[PARAM_D][PARAM_D], B_embed[PARAM_D][PARAM_K*PARAM_D], A3_embed[PARAM_D][PARAM_K], G_embed[PARAM_D];
  poly_qshow_vec_k u_embed[PARAM_D];
  poly_qshow_vec_m1 s1;

  keys_init(&pk, &sk);
  pre_sig_init(&pre_sig);
  rand_init(&rand);
  bsig_init(&bsig);
  poly_q_vec_d_init(cmt);
  poly_q_vec_d_init(v11);
  poly_q_vec_d_init(w1H[0]);
  poly_q_vec_d_init(w1H[1]);
  for (int k = 0; k < PARAM_K; k++) {
    poly_q_vec_d_init(w2H[k]);
  }
  poly_q_vec_k_init(w3H);
  poly_q_init(tag);
  
  randombytes(state, STATE_BYTES);
  keygen(&pk, &sk);
  tag_gen(tag, state);
  randombytes(msg, PARAM_N/8);
  commit(&rand, cmt, tag, msg, &pk);
  pre_sign_commitment(&pre_sig, &sk, &pk, cmt, tag);
  int is_valid = pre_sig_verify_from_commitment(v11, tag, &pre_sig, cmt, &pk);
  if (!is_valid) {
    printf("FATAL ERROR: benchmarked pre signature is not valid\n");
  }
  complete_decompose(&bsig, w1H, w2H, w3H, v11, &pre_sig, &rand);

  /* ----------- BEGIN: Code under measurement --------- */
  start_timer(t);
  // init
  for (i = 0; i < PARAM_D; i++) {
    for (j = 0; j < PARAM_D; j++) {
      poly_qshow_mat_k_k_init(A_embed[i][j]);
    }
    for (j = 0; j < PARAM_K*PARAM_D; j++) {
      poly_qshow_mat_k_k_init(B_embed[i][j]);
    }
    for (j = 0; j < PARAM_K; j++) {
      poly_qshow_mat_k_k_init(A3_embed[i][j]);
    }
    poly_qshow_mat_k_k_init(G_embed[i]);
    poly_qshow_vec_k_init(u_embed[i]);
  }
  poly_qshow_vec_m1_init(s1);
  // embed
  embed_2(A_embed, B_embed, A3_embed, G_embed, u_embed, s1, &pk, tag, &bsig, w1H, w2H, w3H, msg);
  time = stop_timer(t);
  /* ----------- END: Code under measurement ----------- */

  keys_clear(&pk, &sk);
  pre_sig_clear(&pre_sig);
  rand_clear(&rand);
  bsig_clear(&bsig);
  poly_q_vec_d_clear(cmt);
  poly_q_vec_d_clear(v11);
  poly_q_vec_d_clear(w1H[0]);
  poly_q_vec_d_clear(w1H[1]);
  for (int k = 0; k < PARAM_K; k++) {
    poly_q_vec_d_clear(w2H[k]);
  }
  poly_q_vec_k_clear(w3H);
  poly_q_clear(tag);
  for (i = 0; i < PARAM_D; i++) {
    for (j = 0; j < PARAM_D; j++) {
      poly_qshow_mat_k_k_clear(A_embed[i][j]);
    }
    for (j = 0; j < PARAM_K*PARAM_D; j++) {
      poly_qshow_mat_k_k_clear(B_embed[i][j]);
    }
    for (j = 0; j < PARAM_K; j++) {
      poly_qshow_mat_k_k_clear(A3_embed[i][j]);
    }
    poly_qshow_mat_k_k_clear(G_embed[i]);
    poly_qshow_vec_k_clear(u_embed[i]);
  }
  poly_qshow_vec_m1_clear(s1);
  return time;
}

double embed_2_verifier_bench(timer *t) {
  double time;
  uint8_t state[STATE_BYTES], msg[PARAM_N/8];
  size_t i,j;
  sk_t sk;
  pk_t pk;
  pre_sig_t pre_sig;
  rand_t rand;
  bsig_t bsig;
  poly_q tag;
  poly_q_vec_d cmt, v11, w1H[2], w2H[PARAM_K];
  poly_q_vec_k w3H;
  poly_qshow_mat_k_k A_embed[PARAM_D][PARAM_D], B_embed[PARAM_D][PARAM_K*PARAM_D], A3_embed[PARAM_D][PARAM_K], G_embed[PARAM_D];
  poly_qshow_vec_k u_embed[PARAM_D];

  keys_init(&pk, &sk);
  pre_sig_init(&pre_sig);
  rand_init(&rand);
  bsig_init(&bsig);
  poly_q_vec_d_init(cmt);
  poly_q_vec_d_init(v11);
  poly_q_vec_d_init(w1H[0]);
  poly_q_vec_d_init(w1H[1]);
  for (int k = 0; k < PARAM_K; k++) {
    poly_q_vec_d_init(w2H[k]);
  }
  poly_q_vec_k_init(w3H);
  poly_q_init(tag);
  
  randombytes(state, STATE_BYTES);
  keygen(&pk, &sk);
  tag_gen(tag, state);
  randombytes(msg, PARAM_N/8);
  commit(&rand, cmt, tag, msg, &pk);
  pre_sign_commitment(&pre_sig, &sk, &pk, cmt, tag);
  int is_valid = pre_sig_verify_from_commitment(v11, tag, &pre_sig, cmt, &pk);
  if (!is_valid) {
    printf("FATAL ERROR: benchmarked pre signature is not valid\n");
  }
  complete_decompose(&bsig, w1H, w2H, w3H, v11, &pre_sig, &rand);

  /* ----------- BEGIN: Code under measurement --------- */
  start_timer(t);
  // init
  for (i = 0; i < PARAM_D; i++) {
    for (j = 0; j < PARAM_D; j++) {
      poly_qshow_mat_k_k_init(A_embed[i][j]);
    }
    for (j = 0; j < PARAM_K*PARAM_D; j++) {
      poly_qshow_mat_k_k_init(B_embed[i][j]);
    }
    for (j = 0; j < PARAM_K; j++) {
      poly_qshow_mat_k_k_init(A3_embed[i][j]);
    }
    poly_qshow_mat_k_k_init(G_embed[i]);
    poly_qshow_vec_k_init(u_embed[i]);
  }
  // embed
  embed_2_verifier(A_embed, B_embed, A3_embed, G_embed, u_embed, &pk, &bsig, msg);
  time = stop_timer(t);
  /* ----------- END: Code under measurement ----------- */

  keys_clear(&pk, &sk);
  pre_sig_clear(&pre_sig);
  rand_clear(&rand);
  bsig_clear(&bsig);
  poly_q_vec_d_clear(cmt);
  poly_q_vec_d_clear(v11);
  poly_q_vec_d_clear(w1H[0]);
  poly_q_vec_d_clear(w1H[1]);
  for (int k = 0; k < PARAM_K; k++) {
    poly_q_vec_d_clear(w2H[k]);
  }
  poly_q_vec_k_clear(w3H);
  poly_q_clear(tag);
  for (i = 0; i < PARAM_D; i++) {
    for (j = 0; j < PARAM_D; j++) {
      poly_qshow_mat_k_k_clear(A_embed[i][j]);
    }
    for (j = 0; j < PARAM_K*PARAM_D; j++) {
      poly_qshow_mat_k_k_clear(B_embed[i][j]);
    }
    for (j = 0; j < PARAM_K; j++) {
      poly_qshow_mat_k_k_clear(A3_embed[i][j]);
    }
    poly_qshow_mat_k_k_clear(G_embed[i]);
    poly_qshow_vec_k_clear(u_embed[i]);
  }
  return time;
}