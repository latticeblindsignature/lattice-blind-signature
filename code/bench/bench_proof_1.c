#include "bench_proof_1.h"

#include "bsig_user.h"
#include "bsig_signer.h"
#include "randombytes.h"
#include "random.h"

double prove_1_bench(timer* t) {
  double time;
  uint8_t state[STATE_BYTES], msg[PARAM_N/8], crs_seed[CRS_SEED_BYTES];
  size_t i,j;
  sk_t sk;
  pk_t pk;
  rand_t rand;
  ct_t ct;
  proof_1_t proof_1;
  poly_q tag;
  poly_q_vec_d cmt;
  poly_p_vec_m r_e;
  poly_qiss_mat_k_k A_embed[PARAM_D][PARAM_D], B_embed[PARAM_D][PARAM_K*(PARAM_D+1)], D_embed[PARAM_D], Ae_be_embed[PARAM_DE+1][PARAM_ME];
  poly_qiss_vec_k cmt_embed[PARAM_D], ct_embed[PARAM_DE+1];
  poly_qiss_vec_m1 s1;

  keys_init(&pk, &sk);
  rand_init(&rand);
  ct_init(&ct);
  proof_1_init(&proof_1);
  poly_q_init(tag);
  poly_q_vec_d_init(cmt);
  poly_p_vec_m_init(r_e);
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

  randombytes(state, STATE_BYTES);
  keygen(&pk, &sk);
  randombytes(crs_seed, CRS_SEED_BYTES);
  tag_gen(tag, state);
  randombytes(msg, PARAM_N/8);
  commit(&rand, cmt, tag, msg, &pk);
  encrypt(&ct, r_e, msg, &pk);
  embed_1(A_embed, B_embed, D_embed, Ae_be_embed, cmt_embed, ct_embed, s1, &pk, cmt, &ct, tag, &rand, r_e, msg);
  /* ----------- BEGIN: Code under measurement --------- */
  start_timer(t);
  prove_1(&proof_1, A_embed, B_embed, D_embed, Ae_be_embed, cmt_embed, ct_embed, s1, crs_seed, pk.seed);
  time = stop_timer(t);
  /* ----------- END: Code under measurement ----------- */
  keys_clear(&pk, &sk);
  rand_clear(&rand);
  ct_clear(&ct);
  proof_1_clear(&proof_1);
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

double verify_1_valid_bench(timer* t) {
  double time;
  uint8_t state[STATE_BYTES], msg[PARAM_N/8], crs_seed[CRS_SEED_BYTES];
  size_t i,j;
  sk_t sk;
  pk_t pk;
  rand_t rand;
  ct_t ct;
  proof_1_t proof_1;
  poly_q tag;
  poly_q_vec_d cmt;
  poly_p_vec_m r_e;
  poly_qiss_mat_k_k A_embed[PARAM_D][PARAM_D], B_embed[PARAM_D][PARAM_K*(PARAM_D+1)], D_embed[PARAM_D], Ae_be_embed[PARAM_DE+1][PARAM_ME];
  poly_qiss_vec_k cmt_embed[PARAM_D], ct_embed[PARAM_DE+1];
  poly_qiss_vec_m1 s1;

  keys_init(&pk, &sk);
  rand_init(&rand);
  ct_init(&ct);
  proof_1_init(&proof_1);
  poly_q_init(tag);
  poly_q_vec_d_init(cmt);
  poly_p_vec_m_init(r_e);
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

  randombytes(state, STATE_BYTES);
  keygen(&pk, &sk);
  randombytes(crs_seed, CRS_SEED_BYTES);
  tag_gen(tag, state);
  randombytes(msg, PARAM_N/8);
  commit(&rand, cmt, tag, msg, &pk);
  encrypt(&ct, r_e, msg, &pk);
  embed_1(A_embed, B_embed, D_embed, Ae_be_embed, cmt_embed, ct_embed, s1, &pk, cmt, &ct, tag, &rand, r_e, msg);
  prove_1(&proof_1, A_embed, B_embed, D_embed, Ae_be_embed, cmt_embed, ct_embed, s1, crs_seed, pk.seed);
  /* ----------- BEGIN: Code under measurement --------- */
  start_timer(t);
  int is_valid = verify_1(&proof_1, A_embed, B_embed, D_embed, Ae_be_embed, cmt_embed, ct_embed, crs_seed, pk.seed);
  time = stop_timer(t);
  /* ----------- END: Code under measurement ----------- */
  if (!is_valid) {
    printf("FATAL ERROR: benchmarked proof_1 is not valid\n");
  }
  keys_clear(&pk, &sk);
  rand_clear(&rand);
  ct_clear(&ct);
  proof_1_clear(&proof_1);
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

double verify_1_invalid_c_bench(timer* t) {
  double time;
  uint8_t state[STATE_BYTES], msg[PARAM_N/8], crs_seed[CRS_SEED_BYTES];
  size_t i,j;
  sk_t sk;
  pk_t pk;
  rand_t rand;
  ct_t ct;
  proof_1_t proof_1;
  poly_q tag;
  poly_q_vec_d cmt;
  poly_p_vec_m r_e;
  poly_qiss_mat_k_k A_embed[PARAM_D][PARAM_D], B_embed[PARAM_D][PARAM_K*(PARAM_D+1)], D_embed[PARAM_D], Ae_be_embed[PARAM_DE+1][PARAM_ME];
  poly_qiss_vec_k cmt_embed[PARAM_D], ct_embed[PARAM_DE+1];
  poly_qiss_vec_m1 s1;

  keys_init(&pk, &sk);
  rand_init(&rand);
  ct_init(&ct);
  proof_1_init(&proof_1);
  poly_q_init(tag);
  poly_q_vec_d_init(cmt);
  poly_p_vec_m_init(r_e);
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

  randombytes(state, STATE_BYTES);
  keygen(&pk, &sk);
  randombytes(crs_seed, CRS_SEED_BYTES);
  tag_gen(tag, state);
  randombytes(msg, PARAM_N/8);
  commit(&rand, cmt, tag, msg, &pk);
  encrypt(&ct, r_e, msg, &pk);
  embed_1(A_embed, B_embed, D_embed, Ae_be_embed, cmt_embed, ct_embed, s1, &pk, cmt, &ct, tag, &rand, r_e, msg);
  prove_1(&proof_1, A_embed, B_embed, D_embed, Ae_be_embed, cmt_embed, ct_embed, s1, crs_seed, pk.seed);
  poly_qiss_muladd_constant(proof_1.c, 1, 1);
  /* ----------- BEGIN: Code under measurement --------- */
  start_timer(t);
  int is_valid = verify_1(&proof_1, A_embed, B_embed, D_embed, Ae_be_embed, cmt_embed, ct_embed, crs_seed, pk.seed);
  time = stop_timer(t);
  /* ----------- END: Code under measurement ----------- */
  if (is_valid) {
    printf("FATAL ERROR: benchmarked proof_1 is valid\n");
  }
  keys_clear(&pk, &sk);
  rand_clear(&rand);
  ct_clear(&ct);
  proof_1_clear(&proof_1);
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

double verify_1_invalid_f_bench(timer* t) {
  double time;
  uint8_t state[STATE_BYTES], msg[PARAM_N/8], crs_seed[CRS_SEED_BYTES];
  size_t i,j;
  sk_t sk;
  pk_t pk;
  rand_t rand;
  ct_t ct;
  proof_1_t proof_1;
  poly_q tag;
  poly_q_vec_d cmt;
  poly_p_vec_m r_e;
  poly_qiss_mat_k_k A_embed[PARAM_D][PARAM_D], B_embed[PARAM_D][PARAM_K*(PARAM_D+1)], D_embed[PARAM_D], Ae_be_embed[PARAM_DE+1][PARAM_ME];
  poly_qiss_vec_k cmt_embed[PARAM_D], ct_embed[PARAM_DE+1];
  poly_qiss_vec_m1 s1;

  keys_init(&pk, &sk);
  rand_init(&rand);
  ct_init(&ct);
  proof_1_init(&proof_1);
  poly_q_init(tag);
  poly_q_vec_d_init(cmt);
  poly_p_vec_m_init(r_e);
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

  randombytes(state, STATE_BYTES);
  keygen(&pk, &sk);
  randombytes(crs_seed, CRS_SEED_BYTES);
  tag_gen(tag, state);
  randombytes(msg, PARAM_N/8);
  commit(&rand, cmt, tag, msg, &pk);
  encrypt(&ct, r_e, msg, &pk);
  embed_1(A_embed, B_embed, D_embed, Ae_be_embed, cmt_embed, ct_embed, s1, &pk, cmt, &ct, tag, &rand, r_e, msg);
  prove_1(&proof_1, A_embed, B_embed, D_embed, Ae_be_embed, cmt_embed, ct_embed, s1, crs_seed, pk.seed);
  poly_qiss_muladd_constant(proof_1.f->entries[0], 1, 1);
  /* ----------- BEGIN: Code under measurement --------- */
  start_timer(t);
  int is_valid = verify_1(&proof_1, A_embed, B_embed, D_embed, Ae_be_embed, cmt_embed, ct_embed, crs_seed, pk.seed);
  time = stop_timer(t);
  /* ----------- END: Code under measurement ----------- */
  if (is_valid) {
    printf("FATAL ERROR: benchmarked proof_1 is valid\n");
  }
  keys_clear(&pk, &sk);
  rand_clear(&rand);
  ct_clear(&ct);
  proof_1_clear(&proof_1);
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