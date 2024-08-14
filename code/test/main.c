#include <stdio.h>
#include "poly_p_sampling.h"
#include "poly_q_sampling.h"
#include "poly_qiss_sampling.h"
#include "arith.h"
#include "bsig_user.h"
#include "bsig_signer.h"
#include "bsig_verify.h"
#include "covariance.h"
#include "randombytes.h"
#include "random.h"

#define NTESTS 1
#define NSUBTESTS 5

static int tag_gen_test(void)
{
  int rval = 1;
  uint8_t state[STATE_BYTES];
  poly_q prev_tag, tag;
  int64_t tag_weight;

  printf("\ntag_gen_test\n");
  poly_q_init(prev_tag);
  poly_q_init(tag);

  for (int i = 0; i < NSUBTESTS; i++)
  {
    randombytes(state, STATE_BYTES);
    poly_q_zero(prev_tag);
    for (int j = 0; j < NSUBTESTS; j++)
    {
      tag_gen(tag, state);
      tag_weight = poly_q_weight(tag);

      if (tag_weight == -1) {
        printf("tag_gen generated non binary tag.\n");
        rval = 0;
        goto tag_gen_test_cleanup;
      }

      if (tag_weight != PARAM_W) {
        printf("tag_gen generated binary tag with wrong hamming weight.\n");
        rval = 0;
        goto tag_gen_test_cleanup;
      }

      if (poly_q_equal(tag, prev_tag)) {
        printf("tag_gen generated same tag twice.\n");
        rval = 0;
        goto tag_gen_test_cleanup;
      }
      printf(":");
      fflush(stdout);
    }
  }

tag_gen_test_cleanup:
  poly_q_clear(prev_tag);
  poly_q_clear(tag);

  return rval;
}

static int keygen_test(void)
{
  int rval = 1;
  sk_t sk;
  pk_t pk;
  poly_q_mat_d_d A, RRstar[2][2];

  printf("\nkeygen_test\n");

  keys_init(&pk, &sk);
  poly_q_mat_d_d_init(A);
  poly_q_mat_d_d_init(RRstar[0][0]);
  poly_q_mat_d_d_init(RRstar[0][1]);
  poly_q_mat_d_d_init(RRstar[1][0]);
  poly_q_mat_d_d_init(RRstar[1][1]);

  for (int i = 0; i < NSUBTESTS; i++)
  {
    keygen(&pk, &sk);
    if (sk_sq_spectral_norm(RRstar, sk.R) > PARAM_R_MAX_SQ_SPECTRAL_NORM) {
      printf("keygen generated keys with large spectral norm.\n");
      rval = 0;
      goto keygen_test_cleanup;
    }

    poly_q_mat_d_d_uniform(A, pk.seed, DOMAIN_SEPARATOR_A, 0);
    for (size_t j = 0; j < PARAM_K; j++) {
      poly_q_mat_d_d_mul_mat_d_d(RRstar[0][0], A, sk.R[1][j]);
      poly_q_mat_d_d_add(RRstar[0][0], RRstar[0][0], sk.R[0][j]);
      if (!poly_q_mat_d_d_equal(RRstar[0][0], pk.B[j])) {
        printf("keygen generated keys with B != AR.\n");
        rval = 0;
        goto keygen_test_cleanup;
      }
    }

    poly_q_set_coeff(sk.R[1][0]->rows[0]->entries[0], PARAM_N/2, 3); // changing one coefficient of R
    poly_q_mat_d_d_mul_mat_d_d(RRstar[0][0], A, sk.R[1][0]);
    poly_q_mat_d_d_add(RRstar[0][0], RRstar[0][0], sk.R[0][0]);
    if (poly_q_mat_d_d_equal(RRstar[0][0], pk.B[0])) {
      printf("found R' such that AR' = B.\n");
      rval = 0;
      goto keygen_test_cleanup;
    }

    printf(":");
    fflush(stdout);
  }

keygen_test_cleanup:
  keys_clear(&pk, &sk);
  poly_q_mat_d_d_clear(A);
  poly_q_mat_d_d_clear(RRstar[0][0]);
  poly_q_mat_d_d_clear(RRstar[0][1]);
  poly_q_mat_d_d_clear(RRstar[1][0]);
  poly_q_mat_d_d_clear(RRstar[1][1]);

  return rval;
}

static int commit_test(void)
{
  int rval = 1;
  uint8_t state[STATE_BYTES], msg[PARAM_N/8];
  sk_t sk;
  pk_t pk;
  rand_t rand;
  poly_q tag;
  poly_q_vec_d cmt, tmp;
  int64_t min, max, tmpmin, tmpmax;

  printf("\ncommit_test\n");
  keys_init(&pk, &sk);
  rand_init(&rand);
  poly_q_init(tag);
  poly_q_vec_d_init(cmt);
  poly_q_vec_d_init(tmp);

  for (int i = 0; i < NSUBTESTS; i++)
  {
    randombytes(state, STATE_BYTES);
    keygen(&pk, &sk);
    for (int j = 0; j < NSUBTESTS; j++)
    {
      tag_gen(tag, state);
      randombytes(msg, PARAM_N/8);

      commit(&rand, cmt, tag, msg, &pk);

      msg[0] ^= 1;

      commit(&rand, tmp, tag, msg, &pk);

      if (poly_q_vec_d_equal(cmt, tmp)) {
        printf("commitment not binding.\n");
        rval = 0;
        goto commit_test_cleanup;
      }

      max = poly_q_vec_d_minmax(&min, rand.r1[0]);
      tmpmax = poly_q_vec_d_minmax(&tmpmin, rand.r1[1]);
      max = (tmpmax > max) ? tmpmax : max;
      min = (tmpmin < min) ? tmpmin : min;
      if ((min < -2*PARAM_B1) | (max > 2*PARAM_B1-1)) {
        printf("commitment randomness r1 not in correct interval.\n");
        rval = 0;
        goto commit_test_cleanup;
      }

      max = poly_q_vec_k_minmax(&min, rand.r3);
      for (int k = 0; k < PARAM_K; k++) {
        tmpmax = poly_q_vec_d_minmax(&tmpmin, rand.r2[k]);
        max = (tmpmax > max) ? tmpmax : max;
        min = (tmpmin < min) ? tmpmin : min;
      }
      if ((min < -PARAM_B2) || (max >= PARAM_B2)) {
        printf("commitment randomness [r2|r3] not in correct interval.\n");
        rval = 0;
        goto commit_test_cleanup;
      }
      printf(":");
      fflush(stdout);
    }
  }

commit_test_cleanup:
  keys_clear(&pk, &sk);
  rand_clear(&rand);
  poly_q_clear(tag);
  poly_q_vec_d_clear(cmt);
  poly_q_vec_d_clear(tmp);

  return rval;
}

static int encrypt_test(void)
{
  int rval = 1;
  ct_t ct;
  sk_t sk;
  pk_t pk;
  poly_p_vec_m r_e;
  uint8_t msg[PARAM_N/8];

  printf("\nencrypt_test\n");

  keys_init(&pk, &sk);
  ct_init(&ct);
  poly_p_vec_m_init(r_e);

  for (int i = 0; i < NSUBTESTS; i++)
  {
    keygen(&pk, &sk);
    for (int j = 0; j < NSUBTESTS; j++)
    {
      randombytes(msg, PARAM_N/8);
      // no real test because encryption to the sky. Just to check for runtime errors
      encrypt(&ct, r_e, msg, &pk);
      printf(":");
      fflush(stdout);
    }
  }

  keys_clear(&pk, &sk);
  ct_clear(&ct);
  poly_p_vec_m_clear(r_e);
  return rval;
}

static int pre_sig_test(void)
{
  int rval = 1;
  sk_t sk;
  pk_t pk;
  pre_sig_t pre_sig;
  rand_t rand;
  poly_q tag;
  poly_q_vec_d cmt, v11;
  uint8_t state[STATE_BYTES], msg[PARAM_N/8];

  printf("\npre_sig_test\n");

  keys_init(&pk, &sk);
  pre_sig_init(&pre_sig);
  rand_init(&rand);
  poly_q_vec_d_init(cmt);
  poly_q_vec_d_init(v11);
  poly_q_init(tag);

  for (int i = 0; i < NSUBTESTS; i++)
  {
    randombytes(state, STATE_BYTES);
    keygen(&pk, &sk);
    for (int j = 0; j < NSUBTESTS; j++)
    {
      tag_gen(tag, state);
      randombytes(msg, PARAM_N/8);
      commit(&rand, cmt, tag, msg, &pk);

      pre_sign_commitment(&pre_sig, &sk, &pk, cmt, tag);
      if (!pre_sig_verify_from_commitment(v11, tag, &pre_sig, cmt, &pk))
      {
        printf("pre_sig_verify_from_commitment returned zero for a valid pre signature.\n");
        rval = 0;
        goto pre_sig_test_cleanup;
      }
      
      msg[0] ^= 1;
      commit(&rand, cmt, tag, msg, &pk);

      if (pre_sig_verify_from_commitment(v11, tag, &pre_sig, cmt, &pk))
      {
        printf("pre_sig_verify_from_commitment returned non-zero for a valid pre signature on a commitment to a wrong message.\n");
        rval = 0;
        goto pre_sig_test_cleanup;
      }
      printf(":");
      fflush(stdout);
    }
  }

pre_sig_test_cleanup:
  pre_sig_clear(&pre_sig);
  keys_clear(&pk, &sk);
  rand_clear(&rand);
  poly_q_vec_d_clear(cmt);
  poly_q_vec_d_clear(v11);
  poly_q_clear(tag);
  return rval;
}

static int complete_decompose_test(void)
{
  int rval = 1;
  sk_t sk;
  pk_t pk;
  pre_sig_t pre_sig;
  rand_t rand;
  bsig_t bsig;
  poly_q tag;
  poly_q_vec_d cmt, v11, w1H[2], w2H[PARAM_K];
  poly_q_vec_k w3H;
  uint8_t state[STATE_BYTES], msg[PARAM_N/8];
  int64_t min, max, tmpmin, tmpmax;
  uint64_t sq_norm;

  printf("\ncomplete_decompose_test\n");

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

  for (int i = 0; i < NSUBTESTS; i++)
  {
    randombytes(state, STATE_BYTES);
    keygen(&pk, &sk);
    for (int j = 0; j < NSUBTESTS; j++)
    {
      tag_gen(tag, state);
      randombytes(msg, PARAM_N/8);
      commit(&rand, cmt, tag, msg, &pk);

      pre_sign_commitment(&pre_sig, &sk, &pk, cmt, tag);
      if (!pre_sig_verify_from_commitment(v11, tag, &pre_sig, cmt, &pk))
      {
        printf("pre_sig_verify_from_commitment returned zero for a valid pre signature.\n");
        rval = 0;
        goto complete_decompose_test_cleanup;
      }
      
      complete_decompose(&bsig, w1H, w2H, w3H, v11, &pre_sig, &rand);

      max = poly_q_vec_d_minmax(&min, bsig.w1L[0]);
      tmpmax = poly_q_vec_d_minmax(&tmpmin, bsig.w1L[1]);
      max = (tmpmax > max) ? tmpmax : max;
      min = (tmpmin < min) ? tmpmin : min;
      if ((min < -PARAM_B1) || (max > PARAM_B1-1)) {
        printf("decompose term w1L not in correct interval.\n");
        rval = 0;
        goto complete_decompose_test_cleanup;
      }

      max = poly_q_vec_k_minmax(&min, bsig.w3L);
      for (int k = 0; k < PARAM_K; k++) {
        tmpmax = poly_q_vec_d_minmax(&tmpmin, bsig.w2L[k]);
        max = (tmpmax > max) ? tmpmax : max;
        min = (tmpmin < min) ? tmpmin : min;
      }
      if ((min < -PARAM_B2) || (max >= PARAM_B2)) {
        printf("decompose term [w2L|w3L] not in correct interval.\n");
        rval = 0;
        goto complete_decompose_test_cleanup;
      }

      sq_norm = poly_q_vec_d_norm2(w1H[0]);
      sq_norm += poly_q_vec_d_norm2(w1H[1]);
      if (sq_norm > PARAM_B1_PRIME_SQ) {
        printf("decompose term w1H exceeds bound.\n");
        rval = 0;
        goto complete_decompose_test_cleanup;
      }

      sq_norm = poly_q_vec_k_norm2(w3H);
      for (int k = 0; k < PARAM_K; k++) {
        sq_norm += poly_q_vec_d_norm2(w2H[k]);
      }
      if (sq_norm > PARAM_B2_PRIME_SQ) {
        printf("decompose term [w2H|w3H] exceeds bound.\n");
        rval = 0;
        goto complete_decompose_test_cleanup;
      }
      
      printf(":");
      fflush(stdout);
    }
  }

complete_decompose_test_cleanup:
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
  return rval;
}

static int embed_1_test(void)
{
  int rval = 1;
  uint8_t state[STATE_BYTES], msg[PARAM_N/8];
  size_t i,j,k,l;
  sk_t sk;
  pk_t pk;
  rand_t rand;
  ct_t ct;
  poly_q tag;
  poly_q_vec_d cmt;
  poly_p_vec_m r_e;
  poly_qiss_mat_k_k A_embed[PARAM_D][PARAM_D], B_embed[PARAM_D][PARAM_K*(PARAM_D+1)], D_embed[PARAM_D], Ae_be_embed[PARAM_DE+1][PARAM_ME];
  poly_qiss_vec_k cmt_embed[PARAM_D], ct_embed[PARAM_DE+1];
  poly_qiss tmp, tmp2;
  poly_qiss_vec_m1 s1;
  coeff_qiss tmp_coeff;

  printf("\nembed_1_test\n");
  keys_init(&pk, &sk);
  rand_init(&rand);
  ct_init(&ct);
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
  poly_qiss_init(tmp);
  poly_qiss_init(tmp2);

  for (i = 0; i < NSUBTESTS; i++)
  {
    randombytes(state, STATE_BYTES);
    keygen(&pk, &sk);
    for (j = 0; j < NSUBTESTS; j++)
    {
      tag_gen(tag, state);
      randombytes(msg, PARAM_N/8);
      commit(&rand, cmt, tag, msg, &pk);
      encrypt(&ct, r_e, msg, &pk);

      embed_1(A_embed, B_embed, D_embed, Ae_be_embed, cmt_embed, ct_embed, s1, &pk, cmt, &ct, tag, &rand, r_e, msg);

      for (k = 0; k < PARAM_D*PARAM_K_ISS; k++) {
        // A_theta.s_{1,1}
        poly_qiss_mul_scalar(tmp, s1->entries[k], PARAM_Q1_ISS);
        for (l = 0; l < PARAM_D*PARAM_K_ISS; l++) {
          poly_qiss_mul(tmp2, A_embed[k / PARAM_K_ISS][l / PARAM_K_ISS]->rows[k % PARAM_K_ISS]->entries[l % PARAM_K_ISS], s1->entries[IDX_R12_ISS + l]);
          poly_qiss_add(tmp, tmp, tmp2);
        }
        // + B_theta.s_{1,23}
        for (l = 0; l < PARAM_K*(PARAM_D+1)*PARAM_K_ISS; l++) {
          poly_qiss_mul(tmp2, B_embed[k / PARAM_K_ISS][l / PARAM_K_ISS]->rows[k % PARAM_K_ISS]->entries[l % PARAM_K_ISS], s1->entries[IDX_R23_ISS + l]);
          poly_qiss_add(tmp, tmp, tmp2);
        }
        // + D_theta.s_{1,h}
        for (l = 0; l < PARAM_K_ISS; l++) {
          poly_qiss_mul(tmp2, D_embed[k / PARAM_K_ISS]->rows[k % PARAM_K_ISS]->entries[l], s1->entries[IDX_MSG_ISS + l]);
          poly_qiss_add(tmp, tmp, tmp2);
        }
        // = s_{1,b}.u
        poly_qiss_mul(tmp2, cmt_embed[k / PARAM_K_ISS]->entries[k % PARAM_K_ISS], s1->entries[IDX_B_ISS]);

        if (!poly_qiss_equal(tmp, tmp2))
        {
          printf("embed_1 fails to embed the correct mod q equation.\n");
          rval = 0;
          goto embed_1_test_cleanup;
        }
      }
      
      for (k = 0; k < (PARAM_DE+1)*PARAM_K_ISS; k++) {
        // A_e'.s_{1,e}
        if (k < PARAM_DE*PARAM_K_ISS) {
          poly_qiss_zero(tmp);
        } else {
          poly_qiss_mul_scalar(tmp, s1->entries[IDX_MSG_ISS + (k % PARAM_K_ISS)], PARAM_HALF_P);
        }
        for (l = 0; l < PARAM_ME*PARAM_K_ISS; l++) {
          poly_qiss_mul(tmp2, Ae_be_embed[k / PARAM_K_ISS][l / PARAM_K_ISS]->rows[k % PARAM_K_ISS]->entries[l % PARAM_K_ISS], s1->entries[IDX_RE_ISS + l]);
          poly_qiss_add(tmp, tmp, tmp2);
        }
        
        // - s_{1,b}.ct'
        poly_qiss_mul(tmp2, ct_embed[k / PARAM_K_ISS]->entries[k % PARAM_K_ISS], s1->entries[IDX_B_ISS]);
        poly_qiss_sub(tmp, tmp, tmp2);

        // manual mod p
        for (l = 0; l < PARAM_N_ISS; l++) {
          tmp_coeff = poly_qiss_get_coeff_centered(tmp, l);
          if ((tmp_coeff % PARAM_P) != 0) {
            printf("embed_1 fails to embed the correct mod p equation.\n");
            rval = 0;
            goto embed_1_test_cleanup;
          }
        }
      }

      printf(":");
      fflush(stdout);
    }
  }

embed_1_test_cleanup:
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
  poly_qiss_clear(tmp);
  poly_qiss_clear(tmp2);

  return rval;
}

static int prove_1_test(void)
{
  int rval = 1;
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

  printf("\nprove_1_test\n");
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

  for (i = 0; i < NSUBTESTS; i++)
  {
    randombytes(state, STATE_BYTES);
    keygen(&pk, &sk);
    randombytes(crs_seed, CRS_SEED_BYTES);
    for (j = 0; j < NSUBTESTS; j++)
    {
      tag_gen(tag, state);
      randombytes(msg, PARAM_N/8);
      commit(&rand, cmt, tag, msg, &pk);
      encrypt(&ct, r_e, msg, &pk);

      embed_1(A_embed, B_embed, D_embed, Ae_be_embed, cmt_embed, ct_embed, s1, &pk, cmt, &ct, tag, &rand, r_e, msg);
      prove_1(&proof_1, A_embed, B_embed, D_embed, Ae_be_embed, cmt_embed, ct_embed, s1, crs_seed, pk.seed);

      if (!verify_1(&proof_1, A_embed, B_embed, D_embed, Ae_be_embed, cmt_embed, ct_embed, crs_seed, pk.seed)) {
        printf("verify_1 returned zero for a valid proof.\n");
        rval = 0;
        goto prove_1_test_cleanup;
      }

      poly_qiss_muladd_constant(proof_1.c, 1, 1);

      if (verify_1(&proof_1, A_embed, B_embed, D_embed, Ae_be_embed, cmt_embed, ct_embed, crs_seed, pk.seed)) {
        printf("verify_1 returned one for an invalid proof.\n");
        rval = 0;
        goto prove_1_test_cleanup;
      }
      
      printf(":");
      fflush(stdout);
    }
  }

prove_1_test_cleanup:
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

  return rval;
}

static int prove_2_test(void)
{
  int rval = 1;
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

  printf("\nprove_2_test\n");
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

  for (i = 0; i < NSUBTESTS; i++)
  {
    randombytes(state, STATE_BYTES);
    keygen(&pk, &sk);
    randombytes(bsig.crs_seed, CRS_SEED_BYTES);
    for (j = 0; j < NSUBTESTS; j++)
    {
      tag_gen(tag, state);
      randombytes(msg, PARAM_N/8);
      commit(&rand, cmt, tag, msg, &pk);
      pre_sign_commitment(&pre_sig, &sk, &pk, cmt, tag);
      if (!pre_sig_verify_from_commitment(v11, tag, &pre_sig, cmt, &pk)) {
        printf("FATAL ERROR: pre signature is not valid\n");
        rval = 0;
        goto prove_2_test_cleanup;
      }
      complete_decompose(&bsig, w1H, w2H, w3H, v11, &pre_sig, &rand);
      embed_2(A_embed, B_embed, A3_embed, G_embed, u_embed, s1, &pk, tag, &bsig, w1H, w2H, w3H, msg);
      
      prove_2(&bsig.proof_2, A_embed, B_embed, A3_embed, G_embed, u_embed, s1, bsig.crs_seed, pk.seed);

      if (!verify_2(&bsig.proof_2, A_embed, B_embed, A3_embed, G_embed, u_embed, bsig.crs_seed, pk.seed)) {
        printf("verify_2 returned zero for a valid proof.\n");
        rval = 0;
        goto prove_2_test_cleanup;
      }

      poly_qshow_muladd_constant(bsig.proof_2.c, 1, 1);

      if (verify_2(&bsig.proof_2, A_embed, B_embed, A3_embed, G_embed, u_embed, bsig.crs_seed, pk.seed)) {
        printf("verify_2 returned one for an invalid proof.\n");
        rval = 0;
        goto prove_2_test_cleanup;
      }
      
      printf(":");
      fflush(stdout);
    }
  }

prove_2_test_cleanup:
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

  return rval;
}

static int full_bsig_test(void)
{
  int rval = 1;
  uint8_t state[STATE_BYTES], msg[PARAM_N/8];
  size_t i,j;
  sk_t sk;
  pk_t pk;
  pre_sig_t pre_sig;
  rand_t rand;
  ct_t ct;
  proof_1_t proof_1;
  bsig_t bsig;
  poly_q tag;
  poly_p_vec_m r_e;
  poly_q_vec_d cmt, v11, w1H[2], w2H[PARAM_K];
  poly_q_vec_k w3H;
  poly_qiss_mat_k_k A_embed_1[PARAM_D][PARAM_D], B_embed_1[PARAM_D][PARAM_K*(PARAM_D+1)], D_embed_1[PARAM_D], Ae_be_embed_1[PARAM_DE+1][PARAM_ME];
  poly_qiss_vec_k cmt_embed_1[PARAM_D], ct_embed_1[PARAM_DE+1];
  poly_qiss_vec_m1 s1_1;
  poly_qshow_mat_k_k A_embed_2[PARAM_D][PARAM_D], B_embed_2[PARAM_D][PARAM_K*PARAM_D], A3_embed_2[PARAM_D][PARAM_K], G_embed_2[PARAM_D];
  poly_qshow_vec_k u_embed_2[PARAM_D];
  poly_qshow_vec_m1 s1_2;

  printf("\nbsig_test\n");
  keys_init(&pk, &sk);
  pre_sig_init(&pre_sig);
  rand_init(&rand);
  ct_init(&ct);
  proof_1_init(&proof_1);
  bsig_init(&bsig);
  poly_p_vec_m_init(r_e);
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
      poly_qiss_mat_k_k_init(A_embed_1[i][j]);
      poly_qshow_mat_k_k_init(A_embed_2[i][j]);
    }
    for (j = 0; j < PARAM_K*(PARAM_D+1); j++) {
      poly_qiss_mat_k_k_init(B_embed_1[i][j]);
    }
    for (j = 0; j < PARAM_K*PARAM_D; j++) {
      poly_qshow_mat_k_k_init(B_embed_2[i][j]);
    }
    for (j = 0; j < PARAM_K; j++) {
      poly_qshow_mat_k_k_init(A3_embed_2[i][j]);
    }
    poly_qiss_mat_k_k_init(D_embed_1[i]);
    poly_qiss_vec_k_init(cmt_embed_1[i]);
    poly_qshow_mat_k_k_init(G_embed_2[i]);
    poly_qshow_vec_k_init(u_embed_2[i]);
  }
  for (i = 0; i < PARAM_DE + 1; i++) {
    for (j = 0; j < PARAM_ME; j++) {
      poly_qiss_mat_k_k_init(Ae_be_embed_1[i][j]);
    }
    poly_qiss_vec_k_init(ct_embed_1[i]);
  }
  poly_qiss_vec_m1_init(s1_1);
  poly_qshow_vec_m1_init(s1_2);

  for (i = 0; i < NSUBTESTS; i++)
  {
    randombytes(state, STATE_BYTES);
    keygen(&pk, &sk);
    randombytes(bsig.crs_seed, CRS_SEED_BYTES);
    tag_gen(tag, state);
    randombytes(msg, PARAM_N/8);
    commit(&rand, cmt, tag, msg, &pk);
    encrypt(&ct, r_e, msg, &pk);
    embed_1(A_embed_1, B_embed_1, D_embed_1, Ae_be_embed_1, cmt_embed_1, ct_embed_1, s1_1, &pk, cmt, &ct, tag, &rand, r_e, msg);
    prove_1(&proof_1, A_embed_1, B_embed_1, D_embed_1, Ae_be_embed_1, cmt_embed_1, ct_embed_1, s1_1, bsig.crs_seed, pk.seed);
    if (!verify_1(&proof_1, A_embed_1, B_embed_1, D_embed_1, Ae_be_embed_1, cmt_embed_1, ct_embed_1, bsig.crs_seed, pk.seed)) {
      printf("verify_1 returned zero for a valid proof.\n");
      rval = 0;
      goto bsig_test_cleanup;
    }
    pre_sign_commitment(&pre_sig, &sk, &pk, cmt, tag);
    if (!pre_sig_verify_from_commitment(v11, tag, &pre_sig, cmt, &pk)) {
      printf("FATAL ERROR: pre signature is not valid\n");
      rval = 0;
      goto bsig_test_cleanup;
    }
    complete_decompose(&bsig, w1H, w2H, w3H, v11, &pre_sig, &rand);
    embed_2(A_embed_2, B_embed_2, A3_embed_2, G_embed_2, u_embed_2, s1_2, &pk, tag, &bsig, w1H, w2H, w3H, msg);
    prove_2(&bsig.proof_2, A_embed_2, B_embed_2, A3_embed_2, G_embed_2, u_embed_2, s1_2, bsig.crs_seed, pk.seed);

    if (!bsig_verify(&bsig, &pk, msg)) {
      printf("bsig_verify returned zero for a valid blind signature.\n");
      rval = 0;
      goto bsig_test_cleanup;
    }

    poly_qshow_muladd_constant(bsig.proof_2.c, 1, 1);

    if (bsig_verify(&bsig, &pk, msg)) {
      printf("bsig_verify returned one for a invalid blind signature.\n");
      rval = 0;
      goto bsig_test_cleanup;
    }

    printf(":");
    fflush(stdout);
  }

bsig_test_cleanup:
  keys_clear(&pk, &sk);
  pre_sig_clear(&pre_sig);
  rand_clear(&rand);
  ct_clear(&ct);
  proof_1_clear(&proof_1);
  bsig_clear(&bsig);
  poly_p_vec_m_clear(r_e);
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
      poly_qiss_mat_k_k_clear(A_embed_1[i][j]);
      poly_qshow_mat_k_k_clear(A_embed_2[i][j]);
    }
    for (j = 0; j < PARAM_K*(PARAM_D+1); j++) {
      poly_qiss_mat_k_k_clear(B_embed_1[i][j]);
    }
    for (j = 0; j < PARAM_K*PARAM_D; j++) {
      poly_qshow_mat_k_k_clear(B_embed_2[i][j]);
    }
    for (j = 0; j < PARAM_K; j++) {
      poly_qshow_mat_k_k_clear(A3_embed_2[i][j]);
    }
    poly_qiss_mat_k_k_clear(D_embed_1[i]);
    poly_qiss_vec_k_clear(cmt_embed_1[i]);
    poly_qshow_mat_k_k_clear(G_embed_2[i]);
    poly_qshow_vec_k_clear(u_embed_2[i]);
  }
  for (i = 0; i < PARAM_DE + 1; i++) {
    for (j = 0; j < PARAM_ME; j++) {
      poly_qiss_mat_k_k_clear(Ae_be_embed_1[i][j]);
    }
    poly_qiss_vec_k_clear(ct_embed_1[i]);
  }
  poly_qiss_vec_m1_clear(s1_1);
  poly_qshow_vec_m1_clear(s1_2);

  return rval;
}

int main(void) {
  int pass = 1;
  arith_setup();
  random_init();
  printf("Hello from the unit tests.\n");
  for (int i = 0; i < NTESTS; i++)
  {
    pass &= tag_gen_test();
    pass &= keygen_test();
    pass &= commit_test();
    pass &= encrypt_test(); // just check for runtime errors
    pass &= embed_1_test();
    pass &= prove_1_test();
    pass &= pre_sig_test();
    pass &= complete_decompose_test();
    pass &= prove_2_test();
    pass &= full_bsig_test();

    if (!pass)
    {
      printf("FAILED!\n");
      break;
    } else {
      printf(".");
    }
  }
  if (pass)
  {
    printf("\npassed.\n");
  }
  arith_teardown();
  return 0;
}
