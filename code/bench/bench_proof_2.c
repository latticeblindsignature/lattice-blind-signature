#include "bench_proof_2.h"

#include "bsig_user.h"
#include "bsig_signer.h"
#include "randombytes.h"
#include "random.h"

double prove_2_bench(timer* t) {
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

  randombytes(state, STATE_BYTES);
  keygen(&pk, &sk);
  randombytes(bsig.crs_seed, CRS_SEED_BYTES);
  tag_gen(tag, state);
  randombytes(msg, PARAM_N/8);
  commit(&rand, cmt, tag, msg, &pk);
  pre_sign_commitment(&pre_sig, &sk, &pk, cmt, tag);
  if (!pre_sig_verify_from_commitment(v11, tag, &pre_sig, cmt, &pk)) {
    printf("FATAL ERROR: pre signature is not valid\n");
  }
  complete_decompose(&bsig, w1H, w2H, w3H, v11, &pre_sig, &rand);
  embed_2(A_embed, B_embed, A3_embed, G_embed, u_embed, s1, &pk, tag, &bsig, w1H, w2H, w3H, msg);
 
  /* ----------- BEGIN: Code under measurement --------- */
  start_timer(t);
  prove_2(&bsig.proof_2, A_embed, B_embed, A3_embed, G_embed, u_embed, s1, bsig.crs_seed, pk.seed);
  time = stop_timer(t);
  /* ----------- END: Code under measurement ----------- */

  if (!verify_2(&bsig.proof_2, A_embed, B_embed, A3_embed, G_embed, u_embed, bsig.crs_seed, pk.seed)) {
    printf("verify_2 returned zero for a valid proof.\n");
  }

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

double verify_2_valid_bench(timer* t) {
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

  randombytes(state, STATE_BYTES);
  keygen(&pk, &sk);
  randombytes(bsig.crs_seed, CRS_SEED_BYTES);
  tag_gen(tag, state);
  randombytes(msg, PARAM_N/8);
  commit(&rand, cmt, tag, msg, &pk);
  pre_sign_commitment(&pre_sig, &sk, &pk, cmt, tag);
  if (!pre_sig_verify_from_commitment(v11, tag, &pre_sig, cmt, &pk)) {
    printf("FATAL ERROR: pre signature is not valid\n");
  }
  complete_decompose(&bsig, w1H, w2H, w3H, v11, &pre_sig, &rand);
  embed_2(A_embed, B_embed, A3_embed, G_embed, u_embed, s1, &pk, tag, &bsig, w1H, w2H, w3H, msg);
  prove_2(&bsig.proof_2, A_embed, B_embed, A3_embed, G_embed, u_embed, s1, bsig.crs_seed, pk.seed);
  /* ----------- BEGIN: Code under measurement --------- */
  start_timer(t);
  int is_valid = verify_2(&bsig.proof_2, A_embed, B_embed, A3_embed, G_embed, u_embed, bsig.crs_seed, pk.seed);
  time = stop_timer(t);
  /* ----------- END: Code under measurement ----------- */

  if (!is_valid) {
    printf("FATAL ERROR: benchmarked proof_2 is not valid.\n");
  }

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

double verify_2_invalid_c_bench(timer* t) {
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

  randombytes(state, STATE_BYTES);
  keygen(&pk, &sk);
  randombytes(bsig.crs_seed, CRS_SEED_BYTES);
  tag_gen(tag, state);
  randombytes(msg, PARAM_N/8);
  commit(&rand, cmt, tag, msg, &pk);
  pre_sign_commitment(&pre_sig, &sk, &pk, cmt, tag);
  if (!pre_sig_verify_from_commitment(v11, tag, &pre_sig, cmt, &pk)) {
    printf("FATAL ERROR: pre signature is not valid\n");
  }
  complete_decompose(&bsig, w1H, w2H, w3H, v11, &pre_sig, &rand);
  embed_2(A_embed, B_embed, A3_embed, G_embed, u_embed, s1, &pk, tag, &bsig, w1H, w2H, w3H, msg);
  prove_2(&bsig.proof_2, A_embed, B_embed, A3_embed, G_embed, u_embed, s1, bsig.crs_seed, pk.seed);
  poly_qshow_muladd_constant(bsig.proof_2.c, 1, 1);
  /* ----------- BEGIN: Code under measurement --------- */
  start_timer(t);
  int is_valid = verify_2(&bsig.proof_2, A_embed, B_embed, A3_embed, G_embed, u_embed, bsig.crs_seed, pk.seed);
  time = stop_timer(t);
  /* ----------- END: Code under measurement ----------- */

  if (is_valid) {
    printf("FATAL ERROR: benchmarked proof_2 is valid for wrong c.\n");
  }

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

double verify_2_invalid_f_bench(timer* t) {
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

  randombytes(state, STATE_BYTES);
  keygen(&pk, &sk);
  randombytes(bsig.crs_seed, CRS_SEED_BYTES);
  tag_gen(tag, state);
  randombytes(msg, PARAM_N/8);
  commit(&rand, cmt, tag, msg, &pk);
  pre_sign_commitment(&pre_sig, &sk, &pk, cmt, tag);
  if (!pre_sig_verify_from_commitment(v11, tag, &pre_sig, cmt, &pk)) {
    printf("FATAL ERROR: pre signature is not valid\n");
  }
  complete_decompose(&bsig, w1H, w2H, w3H, v11, &pre_sig, &rand);
  embed_2(A_embed, B_embed, A3_embed, G_embed, u_embed, s1, &pk, tag, &bsig, w1H, w2H, w3H, msg);
  prove_2(&bsig.proof_2, A_embed, B_embed, A3_embed, G_embed, u_embed, s1, bsig.crs_seed, pk.seed);
  poly_qshow_muladd_constant(bsig.proof_2.f->entries[0], 1, 1);
  /* ----------- BEGIN: Code under measurement --------- */
  start_timer(t);
  int is_valid = verify_2(&bsig.proof_2, A_embed, B_embed, A3_embed, G_embed, u_embed, bsig.crs_seed, pk.seed);
  time = stop_timer(t);
  /* ----------- END: Code under measurement ----------- */

  if (is_valid) {
    printf("FATAL ERROR: benchmarked proof_2 is valid for wrong f.\n");
  }

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