#include "arith.h"
#include "randombytes.h"
#include "poly_q_sampling.h"
#include "poly_p_sampling.h"
#include "bsig_signer.h"
#include "bsig_user.h"
#include "four_squares.h"

/*************************************************
* Name:        rand_init
*
* Description: Initialize structure to host the commitment randomness
*              by calling Flint initialization
*
* Arguments:   - rand_t *rand: pointer to randomness structure
**************************************************/
void rand_init(rand_t *rand) {
  size_t i;
  for (i = 0; i < PARAM_K; i++) {
    poly_q_vec_d_init(rand->r2[i]);
  }
  poly_q_vec_d_init(rand->r1[0]);
  poly_q_vec_d_init(rand->r1[1]);
  poly_q_vec_k_init(rand->r3);
}

/*************************************************
* Name:        rand_clear
*
* Description: Clear structure that hosts the commitment randomness
*              by calling Flint clean up
*
* Arguments:   - rand_t *sig: pointer to randomness structure
**************************************************/
void rand_clear(rand_t *rand) {
  size_t i;
  for (i = 0; i < PARAM_K; i++) {
    poly_q_vec_d_clear(rand->r2[i]);
  }
  poly_q_vec_d_clear(rand->r1[0]);
  poly_q_vec_d_clear(rand->r1[1]);
  poly_q_vec_k_clear(rand->r3);
}

/*************************************************
* Name:        ct_init
*
* Description: Initialize structure to host the ciphertext (encryption to the sky)
*              by calling Flint initialization
*
* Arguments:   - ct_t *ct: pointer to ciphertext structure
**************************************************/
void ct_init(ct_t *ct) {
  poly_p_vec_d_init(ct->ct0);
  poly_p_init(ct->ct1);
}

/*************************************************
* Name:        ct_clear
*
* Description: Clear structure that hosts the ciphertext (encryption to the sky)
*              by calling Flint clean up
*
* Arguments:   - ct_t *ct: pointer to ciphertext structure
**************************************************/
void ct_clear(ct_t *ct) {
  poly_p_vec_d_clear(ct->ct0);
  poly_p_clear(ct->ct1);
}

/*************************************************
* Name:        proof_1_init
*
* Description: Initialize structure to host the issuance proof
*              by calling Flint initialization
*
* Arguments:   - proof_1_t *proof_1: pointer to issuance proof structure
**************************************************/
void proof_1_init(proof_1_t *proof_1) {
  poly_qiss_vec_m1_init(proof_1->z1);
  poly_qiss_vec_d_init(proof_1->tA_H);
  poly_qiss_vec_d_init(proof_1->hint);
  poly_qiss_vec_256_l_init(proof_1->tB);
  poly_qiss_vec_l_init(proof_1->f);
  poly_qiss_init(proof_1->t1);
  poly_qiss_init(proof_1->c);
  poly_qiss_vec_m2_d_init(proof_1->z2_1);
}

/*************************************************
* Name:        proof_1_clear
*
* Description: Clear structure that hosts the issuance proof
*              by calling Flint clean up
*
* Arguments:   - proof_1_t *proof_1: pointer to issuance proof structure
**************************************************/
void proof_1_clear(proof_1_t *proof_1) {
  poly_qiss_vec_m1_clear(proof_1->z1);
  poly_qiss_vec_d_clear(proof_1->tA_H);
  poly_qiss_vec_d_clear(proof_1->hint);
  poly_qiss_vec_256_l_clear(proof_1->tB);
  poly_qiss_vec_l_clear(proof_1->f);
  poly_qiss_clear(proof_1->t1);
  poly_qiss_clear(proof_1->c);
  poly_qiss_vec_m2_d_clear(proof_1->z2_1);
}

/*************************************************
* Name:        proof_2_init
*
* Description: Initialize structure to host the show proof
*              by calling Flint initialization
*
* Arguments:   - proof_2_t *proof: pointer to show proof structure
**************************************************/
void proof_2_init(proof_2_t *proof_2) {
  poly_qshow_vec_m1_init(proof_2->z1);
  poly_qshow_vec_d_init(proof_2->tA_H);
  poly_qshow_vec_d_init(proof_2->hint);
  poly_qshow_vec_256_l_init(proof_2->tB);
  poly_qshow_vec_l_init(proof_2->f);
  poly_qshow_init(proof_2->t1);
  poly_qshow_init(proof_2->c);
  poly_qshow_vec_m2_d_init(proof_2->z2_1);
}

/*************************************************
* Name:        proof_2_init
*
* Description: Clear structure that hosts the show proof
*              by calling Flint clean up
*
* Arguments:   - proof_2_t *proof_2: pointer to show proof structure
**************************************************/
void proof_2_clear(proof_2_t *proof_2) {
  poly_qshow_vec_m1_clear(proof_2->z1);
  poly_qshow_vec_d_clear(proof_2->tA_H);
  poly_qshow_vec_d_clear(proof_2->hint);
  poly_qshow_vec_256_l_clear(proof_2->tB);
  poly_qshow_vec_l_clear(proof_2->f);
  poly_qshow_clear(proof_2->t1);
  poly_qshow_clear(proof_2->c);
  poly_qshow_vec_m2_d_clear(proof_2->z2_1);
}

/*************************************************
* Name:        bsig_init
*
* Description: Initialize structure to host the blind signature
*              by calling Flint initialization
*
* Arguments:   - bsig_t *bsig: pointer to blind signature structure
**************************************************/
void bsig_init(bsig_t *bsig) {
  poly_q_vec_d_init(bsig->w1L[0]);
  poly_q_vec_d_init(bsig->w1L[1]);
  for (size_t i = 0; i < PARAM_K; i++) {
    poly_q_vec_d_init(bsig->w2L[i]);
  }
  poly_q_vec_k_init(bsig->w3L);
  proof_2_init(&bsig->proof_2);
}

/*************************************************
* Name:        bsig_clear
*
* Description: Clear structure that hosts the blind signature
*              by calling Flint clean up
*
* Arguments:   - bsig_t *bsig: pointer to blind signature structure
**************************************************/
void bsig_clear(bsig_t *bsig) {
  poly_q_vec_d_clear(bsig->w1L[0]);
  poly_q_vec_d_clear(bsig->w1L[1]);
  for (size_t i = 0; i < PARAM_K; i++) {
    poly_q_vec_d_clear(bsig->w2L[i]);
  }
  poly_q_vec_k_clear(bsig->w3L);
  proof_2_clear(&bsig->proof_2);
}

/*************************************************
* Name:        commit
*
* Description: Compute commitment to m by cmt = (I|A')r1 + (tG - B)r2 + A3.r3 + d.msg
*
* Arguments:   - rand_t *rand: pointer to user commitment randomness structure (initialized)
*              - poly_q_vec_d cmt: polynomial vector for user commitment (initialized)
*              - const poly_q tag: tag polynomial
*              - const uint8_t *msg: pointer to input message byte array (allocated PARAM_N/8 bytes)
*              - pk_t *pk: pointer to signer public key structure
**************************************************/
void commit(rand_t *rand, poly_q_vec_d cmt, const poly_q tag, const uint8_t msg[PARAM_N/8], const pk_t *pk) {
  size_t i,j;
  uint8_t randomness_seed[SEED_BYTES];
  int64_t bexpi;
  poly_q_mat_d_d A;
  poly_q_mat_d_k A3;
  poly_q_vec_d d;
  poly_q m;

  // init matrices and vectors
  poly_q_mat_d_d_init(A);
  poly_q_mat_d_k_init(A3);
  poly_q_vec_d_init(d);
  poly_q_init(m);

  // generate random secret seed for commitment randomness
  randombytes(randomness_seed, SEED_BYTES);

  // expand uniform A', D from seed
  poly_q_mat_d_d_uniform(A, pk->seed, DOMAIN_SEPARATOR_A, 0);
  poly_q_mat_d_k_uniform(A3, pk->seed, DOMAIN_SEPARATOR_A3);
  poly_q_vec_d_uniform(d, pk->seed, DOMAIN_SEPARATOR_D);

  // storing message as polynomial vector
  poly_q_from_bits(m, msg);

  // sample r1, r2, r3
  poly_q_vec_d_sample_r1(rand->r1[0], randomness_seed, DOMAIN_SEPARATOR_R1, 0);
  poly_q_vec_d_sample_r1(rand->r1[1], randomness_seed, DOMAIN_SEPARATOR_R1, 1);
  for (i = 0; i < PARAM_K; i++) {
    poly_q_vec_d_uniform_b2(rand->r2[i], randomness_seed, DOMAIN_SEPARATOR_R2, i);
  }
  poly_q_vec_k_uniform_b2(rand->r3, randomness_seed, DOMAIN_SEPARATOR_R3);

  // cmt = r1[0] + A'.r1[1] + (tag.G - B).r2 + A3.r3 + d.m
  poly_q_mat_d_d_mul_vec_d(cmt, A, rand->r1[1]);
  poly_q_vec_d_add(cmt, cmt, rand->r1[0]);
  poly_q_vec_d_mul_poly(d, d, m);
  poly_q_vec_d_add(cmt, cmt, d); // A, d and m can now be used as temporary variable
  poly_q_mat_d_k_mul_vec_k(d, A3, rand->r3);
  poly_q_vec_d_add(cmt, cmt, d);
  bexpi = 1;
  for (i = 0; i < PARAM_K; i++) {
    // copy -B[i] to A
    poly_q_mat_d_d_neg(A, pk->B[i]);
    for (j = 0; j < PARAM_D; j++) {
      if (i == 0) {
        poly_q_add(A->rows[j]->entries[j], A->rows[j]->entries[j], tag);
      } else {
        poly_q_mul_scalar(m, tag, bexpi);
        poly_q_add(A->rows[j]->entries[j], A->rows[j]->entries[j], m);
      }
    }
    poly_q_mat_d_d_muladd_vec_d(cmt, A, rand->r2[i]);
    bexpi *= PARAM_B;
  }

  // clean up matrices and vectors
  poly_q_mat_d_d_clear(A);
  poly_q_mat_d_k_clear(A3);
  poly_q_vec_d_clear(d);
  poly_q_clear(m);
}

/*************************************************
* Name:        encrypt
*
* Description: Compute encryption to the sky of msg by 
*              (ct0, ct1) = (A_e^T.r_e, b_e^T.r_e + round(p/2).msg)
*
* Arguments:   - ct_t *ct: pointer to ciphertext structure (initialized)
*              - poly_p_vec_m r_e: polynomial vector for encryption randomness (initialized)
*              - const uint8_t *msg: pointer to input message byte array (allocated PARAM_N/8 bytes)
*              - const sep_pk_t *pk: pointer to signer public key structure
**************************************************/
void encrypt(ct_t *ct, poly_p_vec_m r_e, const uint8_t msg[PARAM_N/8], const pk_t *pk) {
  uint8_t randomness_seed[SEED_BYTES];
  poly_p_mat_d_m A_e_T;
  poly_p_vec_m b_e;
  poly_p m;
  uint32_t kappa;
  
  // init matrices and vectors
  poly_p_mat_d_m_init(A_e_T);
  poly_p_vec_m_init(b_e);
  poly_p_init(m);

  // generate random secret seed for commitment randomness
  randombytes(randomness_seed, SEED_BYTES);
  
  // expand uniform A_e^T, b_e from seed
  poly_p_mat_d_m_uniform(A_e_T, pk->seed, DOMAIN_SEPARATOR_AE);
  poly_p_vec_m_uniform(b_e, pk->seed, DOMAIN_SEPARATOR_BE);

  // storing message as polynomial vector and multiplying by round(p/2)
  poly_p_from_bits(m, msg);
  poly_p_mul_scalar(m, m, (coeff_p)PARAM_HALF_P);

  // sample r_e from B_1^{m_e} with bound B_re
  kappa = 0;
  do {
    poly_p_vec_m_binomial(r_e, randomness_seed, kappa++, DOMAIN_SEPARATOR_RE);
  } while (poly_p_vec_m_norm2(r_e) > PARAM_B_RE_SQ);

  // ct0 = A_e^T.r_e    and    ct1 = b_e^T.r_e + round(p/2).m
  poly_p_mat_d_m_mul_vec_m(ct->ct0, A_e_T, r_e);
  poly_p_vec_m_mul_inner(ct->ct1, b_e, r_e);
  poly_p_add(ct->ct1, ct->ct1, m);

  // clean up matrices and vectors
  poly_p_mat_d_m_clear(A_e_T);
  poly_p_vec_m_clear(b_e);
  poly_p_clear(m);
}

/*************************************************
* Name:        embed_1
*
* Description: Embedding the user relation for the issuance proof
*              of commitment opening and verifiable encryption
*
* Arguments:   - poly_qiss_mat_k_k *A_embed: array of polynomial matrices to host subring embedding of q1.A'
*              - poly_qiss_mat_k_k *B_embed: array of polynomial matrices to host subring embedding of q1.[tG - B|A3]
*              - poly_qiss_mat_k_k *D_embed: array of polynomial matrices to host subring embedding of q1.d
*              - poly_qiss_mat_k_k *Ae_be_embed: array of polynomial matrices to host subring embedding of [A_e | b_e]^T
*              - poly_qiss_vec_k *cmt_embed: array of polynomial vectors to host subring embedding of q1.cmt
*              - poly_qiss_vec_k *ct_embed: array of polynomial vectors to host subring embedding of ct
*              - poly_qiss_vec_k *s1: array of polynomial vectors to host subring embedding of witness
*              - const pk_t *pk: pointer to signer public key structure
*              - const poly_q_vec_d cmt: polynomial vector for commitment
*              - const ct_t *ct: pointer to ciphertext structure
*              - const poly_q tag: polynomial for tag
*              - const rand_t *rand: pointer to commitment randomness structure
*              - const poly_p_vec_m r_e: polynomial vector for encryption randomness
*              - const uint8_t *msg: pointer to input message byte array (allocated PARAM_N/8 bytes)
**************************************************/
void embed_1(
    poly_qiss_mat_k_k A_embed[PARAM_D][PARAM_D], 
    poly_qiss_mat_k_k B_embed[PARAM_D][PARAM_K*(PARAM_D+1)], 
    poly_qiss_mat_k_k D_embed[PARAM_D], 
    poly_qiss_mat_k_k Ae_be_embed[PARAM_DE+1][PARAM_ME],
    poly_qiss_vec_k cmt_embed[PARAM_D], 
    poly_qiss_vec_k ct_embed[PARAM_DE+1], 
    poly_qiss_vec_m1 s1, 
    const pk_t *pk,
    const poly_q_vec_d cmt, 
    const ct_t *ct,
    const poly_q tag,
    const rand_t *rand, 
    const poly_p_vec_m r_e, 
    const uint8_t msg[PARAM_N/8]) {
  uint8_t random_sign[1];
  size_t i,j,k;
  uint64_t norm2sq, four_squares_res[4];
  coeff_qiss tmp_coeff;
  poly_q tmp;
  poly_q_mat_d_d A;
  poly_q_mat_d_k A3;
  poly_q_vec_d d;
  poly_p_mat_d_m A_e_T;
  poly_p_vec_m b_e;
  poly_qiss_vec_k tmp_vec_k;

  // init matrices and polynomials
  poly_q_init(tmp);
  poly_q_mat_d_d_init(A);
  poly_q_mat_d_k_init(A3);
  poly_q_vec_d_init(d);
  poly_p_mat_d_m_init(A_e_T);
  poly_p_vec_m_init(b_e);
  poly_qiss_vec_k_init(tmp_vec_k);

  // computing square norms and four-square decompositions
  // B_{r,1}^2 - |r1|^2
  norm2sq = poly_q_vec_d_norm2(rand->r1[0]);
  norm2sq += poly_q_vec_d_norm2(rand->r1[1]);
  assert(norm2sq <= PARAM_B_R1_SQ);
  four_squares(four_squares_res, PARAM_B_R1_SQ - norm2sq);
  for (i = 0; i < 4; i++) {
    poly_qiss_set_coeff(s1->entries[IDX_R23_ISS - 1], i, (coeff_qiss)(four_squares_res[i]));
  }
  // B_{r,2}^2 - |[r2,r3]|^2
  norm2sq = poly_q_vec_k_norm2(rand->r3);
  for (i = 0; i < PARAM_K; i++) {
    norm2sq += poly_q_vec_d_norm2(rand->r2[i]);
  }
  assert(norm2sq <= PARAM_B_R2_SQ);
  four_squares(four_squares_res, PARAM_B_R2_SQ - norm2sq);
  for (i = 0; i < 4; i++) {
    poly_qiss_set_coeff(s1->entries[IDX_RE_ISS - 1], i, (coeff_qiss)(four_squares_res[i]));
  }
  // B_{r,e}^2 - |re|^2
  norm2sq = poly_p_vec_m_norm2(r_e);
  assert(norm2sq <= PARAM_B_RE_SQ);
  four_squares(four_squares_res, PARAM_B_RE_SQ - norm2sq);
  for (i = 0; i < 4; i++) {
    poly_qiss_set_coeff(s1->entries[IDX_MSG_ISS - 1], i, (coeff_qiss)(four_squares_res[i]));
  }

  // embedding witness vector s1 = [theta(r1)|a1 | theta(r23)|a23 | theta(re)|ae | theta(h) | 1] and multiplying by random sign b
  // r1
  for (i = 0; i < 2*PARAM_D; i++) {
    poly_q_qiss_subring_embed_vec_k(tmp_vec_k, rand->r1[i / PARAM_D]->entries[i % PARAM_D], 1);
    for (j = 0; j < PARAM_K_ISS; j++) {
      poly_qiss_set(s1->entries[IDX_R1_ISS + i*PARAM_K_ISS + j], tmp_vec_k->entries[j]);
    }
  }
  // r2
  for (i = 0; i < PARAM_D*PARAM_K; i++) {
    poly_q_qiss_subring_embed_vec_k(tmp_vec_k, rand->r2[i / PARAM_D]->entries[i % PARAM_D], 1);
    for (j = 0; j < PARAM_K_ISS; j++) {
      poly_qiss_set(s1->entries[IDX_R23_ISS + i*PARAM_K_ISS + j], tmp_vec_k->entries[j]);
    }
  } 
  // r3
  for (i = 0; i < PARAM_K; i++) {
    poly_q_qiss_subring_embed_vec_k(tmp_vec_k, rand->r3->entries[i], 1);
    for (j = 0; j < PARAM_K_ISS; j++) {
      poly_qiss_set(s1->entries[IDX_R3_ISS + i*PARAM_K_ISS + j], tmp_vec_k->entries[j]);
    }
  }
  // re
  for (i = 0; i < PARAM_ME; i++) {
    poly_p_qiss_subring_embed_vec_k(tmp_vec_k, r_e->entries[i]);
    for (j = 0; j < PARAM_K_ISS; j++) {
      poly_qiss_set(s1->entries[IDX_RE_ISS + i*PARAM_K_ISS + j], tmp_vec_k->entries[j]);
    }
  }
  // msg
  poly_q_from_bits(tmp, msg);
  poly_q_qiss_subring_embed_vec_k(tmp_vec_k, tmp, 1);
  for (j = 0; j < PARAM_K_ISS; j++) {
    poly_qiss_set(s1->entries[IDX_MSG_ISS + j], tmp_vec_k->entries[j]);
  }
  // b
  poly_qiss_zero(s1->entries[IDX_B_ISS]);
  poly_qiss_set_coeff(s1->entries[IDX_B_ISS], 0, 1);

  // sampling random sign b
  randombytes(random_sign, 1);
  tmp_coeff = (coeff_qiss)(((uint64_t)random_sign[0] >> 7) << 1) - 1;
  poly_qiss_vec_m1_mul_scalar(s1, s1, tmp_coeff); // multiply by b


  // embedding commitment u = q1 * theta(cmt)
  for (i = 0; i < PARAM_D; i++) {
    poly_q_qiss_subring_embed_vec_k(cmt_embed[i], cmt->entries[i], PARAM_Q1_ISS);
  }

  // embedding ciphertext ct' = theta([ct0 | ct1])
  for (i = 0; i < PARAM_DE; i++) {
    poly_p_qiss_subring_embed_vec_k(ct_embed[i], ct->ct0->entries[i]);
  }
  poly_p_qiss_subring_embed_vec_k(ct_embed[PARAM_DE], ct->ct1);

  // expanding uniform public parameters
  poly_q_mat_d_d_uniform(A, pk->seed, DOMAIN_SEPARATOR_A, 0);
  poly_q_mat_d_k_uniform(A3, pk->seed, DOMAIN_SEPARATOR_A3);
  poly_q_vec_d_uniform(d, pk->seed, DOMAIN_SEPARATOR_D);
  poly_p_mat_d_m_uniform(A_e_T, pk->seed, DOMAIN_SEPARATOR_AE);
  poly_p_vec_m_uniform(b_e, pk->seed, DOMAIN_SEPARATOR_BE);

  // embedding matrices
  // A, d
  for (i = 0; i < PARAM_D; i++) {
    for (j = 0; j < PARAM_D; j++) {
      poly_q_qiss_subring_embed_mat_k_k(A_embed[i][j], A->rows[i]->entries[j], PARAM_Q1_ISS); // A'
    }
    poly_q_qiss_subring_embed_mat_k_k(D_embed[i], d->entries[i], PARAM_Q1_ISS); // d
  } // A and m can now be used as temporary variable
  // tG - B
  tmp_coeff = 1;
  for (i = 0; i < PARAM_K; i++) {
    // copy -B[i] to A
    poly_q_mat_d_d_neg(A, pk->B[i]);
    for (j = 0; j < PARAM_D; j++) {
      if (i == 0) {
        poly_q_add(A->rows[j]->entries[j], A->rows[j]->entries[j], tag);
      } else {
        poly_q_mul_scalar(tmp, tag, tmp_coeff);
        poly_q_add(A->rows[j]->entries[j], A->rows[j]->entries[j], tmp);
      }
    }
    // A contains the i-th dxd block of tG-B
    for (j = 0; j < PARAM_D; j++) {
      for (k = 0; k < PARAM_D; k++) {
        poly_q_qiss_subring_embed_mat_k_k(B_embed[j][k + i*PARAM_D], A->rows[j]->entries[k], PARAM_Q1_ISS); 
      }
    }
    tmp_coeff *= PARAM_B;
  }
  // A3
  for (i = 0; i < PARAM_D; i++) {
    for (j = 0; j < PARAM_K; j++) {
      poly_q_qiss_subring_embed_mat_k_k(B_embed[i][PARAM_D*PARAM_K + j], A3->rows[i]->entries[j], PARAM_Q1_ISS); 
    }
  }
  // A_e^T
  for (j = 0; j < PARAM_ME; j++) {
    for (i = 0; i < PARAM_DE; i++) {
      poly_p_qiss_subring_embed_mat_k_k(Ae_be_embed[i][j], A_e_T->rows[i]->entries[j]); 
    }
    poly_p_qiss_subring_embed_mat_k_k(Ae_be_embed[PARAM_DE][j], b_e->entries[j]); 
  }

  // clean up matrices and polynomials
  poly_q_clear(tmp);
  poly_q_mat_d_d_clear(A);
  poly_q_mat_d_k_clear(A3);
  poly_q_vec_d_clear(d);
  poly_p_mat_d_m_clear(A_e_T);
  poly_p_vec_m_clear(b_e);
  poly_qiss_vec_k_clear(tmp_vec_k);
}

/*************************************************
* Name:        embed_1_verifier
*
* Description: Embedding the user relation for the issuance proof
*              of commitment opening and verifiable encryption on the
*              verifier side
*
* Arguments:   - poly_qiss_mat_k_k *A_embed: array of polynomial matrices to host subring embedding of q1.A'
*              - poly_qiss_mat_k_k *B_embed: array of polynomial matrices to host subring embedding of q1.[tG - B|A3]
*              - poly_qiss_mat_k_k *D_embed: array of polynomial matrices to host subring embedding of q1.d
*              - poly_qiss_mat_k_k *Ae_be_embed: array of polynomial matrices to host subring embedding of [A_e | b_e]^T
*              - poly_qiss_vec_k *cmt_embed: array of polynomial vectors to host subring embedding of q1.cmt
*              - poly_qiss_vec_k *ct_embed: array of polynomial vectors to host subring embedding of ct
*              - const pk_t *pk: pointer to signer public key structure
*              - const poly_q_vec_d cmt: polynomial vector for commitment
*              - const ct_t *ct: pointer to ciphertext structure
*              - const poly_q tag: polynomial for tag
**************************************************/
void embed_1_verifier(
    poly_qiss_mat_k_k A_embed[PARAM_D][PARAM_D], 
    poly_qiss_mat_k_k B_embed[PARAM_D][PARAM_K*(PARAM_D+1)], 
    poly_qiss_mat_k_k D_embed[PARAM_D], 
    poly_qiss_mat_k_k Ae_be_embed[PARAM_DE+1][PARAM_ME],
    poly_qiss_vec_k cmt_embed[PARAM_D], 
    poly_qiss_vec_k ct_embed[PARAM_DE+1], 
    const pk_t *pk,
    const poly_q_vec_d cmt, 
    const ct_t *ct,
    const poly_q tag) {
  size_t i,j,k;
  coeff_qiss tmp_coeff;
  poly_q tmp;
  poly_q_mat_d_d A;
  poly_q_mat_d_k A3;
  poly_q_vec_d d;
  poly_p_mat_d_m A_e_T;
  poly_p_vec_m b_e;

  // init matrices and polynomials
  poly_q_init(tmp);
  poly_q_mat_d_d_init(A);
  poly_q_mat_d_k_init(A3);
  poly_q_vec_d_init(d);
  poly_p_mat_d_m_init(A_e_T);
  poly_p_vec_m_init(b_e);

  // embedding commitment u = q1 * theta(cmt)
  for (i = 0; i < PARAM_D; i++) {
    poly_q_qiss_subring_embed_vec_k(cmt_embed[i], cmt->entries[i], PARAM_Q1_ISS);
  }

  // embedding ciphertext ct' = theta([ct0 | ct1])
  for (i = 0; i < PARAM_DE; i++) {
    poly_p_qiss_subring_embed_vec_k(ct_embed[i], ct->ct0->entries[i]);
  }
  poly_p_qiss_subring_embed_vec_k(ct_embed[PARAM_DE], ct->ct1);

  // expanding uniform public parameters
  poly_q_mat_d_d_uniform(A, pk->seed, DOMAIN_SEPARATOR_A, 0);
  poly_q_mat_d_k_uniform(A3, pk->seed, DOMAIN_SEPARATOR_A3);
  poly_q_vec_d_uniform(d, pk->seed, DOMAIN_SEPARATOR_D);
  poly_p_mat_d_m_uniform(A_e_T, pk->seed, DOMAIN_SEPARATOR_AE);
  poly_p_vec_m_uniform(b_e, pk->seed, DOMAIN_SEPARATOR_BE);

  // embedding matrices
  // A, d
  for (i = 0; i < PARAM_D; i++) {
    for (j = 0; j < PARAM_D; j++) {
      poly_q_qiss_subring_embed_mat_k_k(A_embed[i][j], A->rows[i]->entries[j], PARAM_Q1_ISS); // A'
    }
    poly_q_qiss_subring_embed_mat_k_k(D_embed[i], d->entries[i], PARAM_Q1_ISS); // d
  } // A and m can now be used as temporary variable
  // tG - B
  tmp_coeff = 1;
  for (i = 0; i < PARAM_K; i++) {
    // copy -B[i] to A
    poly_q_mat_d_d_neg(A, pk->B[i]);
    for (j = 0; j < PARAM_D; j++) {
      if (i == 0) {
        poly_q_add(A->rows[j]->entries[j], A->rows[j]->entries[j], tag);
      } else {
        poly_q_mul_scalar(tmp, tag, tmp_coeff);
        poly_q_add(A->rows[j]->entries[j], A->rows[j]->entries[j], tmp);
      }
    }
    // A contains the i-th dxd block of tG-B
    for (j = 0; j < PARAM_D; j++) {
      for (k = 0; k < PARAM_D; k++) {
        poly_q_qiss_subring_embed_mat_k_k(B_embed[j][k + i*PARAM_D], A->rows[j]->entries[k], PARAM_Q1_ISS); 
      }
    }
    tmp_coeff *= PARAM_B;
  }
  // A3
  for (i = 0; i < PARAM_D; i++) {
    for (j = 0; j < PARAM_K; j++) {
      poly_q_qiss_subring_embed_mat_k_k(B_embed[i][PARAM_D*PARAM_K + j], A3->rows[i]->entries[j], PARAM_Q1_ISS); 
    }
  }
  // A_e^T
  for (j = 0; j < PARAM_ME; j++) {
    for (i = 0; i < PARAM_DE; i++) {
      poly_p_qiss_subring_embed_mat_k_k(Ae_be_embed[i][j], A_e_T->rows[i]->entries[j]); 
    }
    poly_p_qiss_subring_embed_mat_k_k(Ae_be_embed[PARAM_DE][j], b_e->entries[j]); 
  }

  // clean up matrices and polynomials
  poly_q_clear(tmp);
  poly_q_mat_d_d_clear(A);
  poly_q_mat_d_k_clear(A3);
  poly_q_vec_d_clear(d);
  poly_p_mat_d_m_clear(A_e_T);
  poly_p_vec_m_clear(b_e);
}

/*************************************************
* Name:        pre_sig_verify
*
* Description: User verification of the obtained pre signature
*
* Arguments:   - poly_q_vec_d v11: polynomial vector to host the reconstructed v11 (initialized)
*              - const poly_q tag: polynomial with the tag used in this issuance
*              - const pre_sig_t *pre_sig: pointer to the pre signature structure 
*              - const poly_q_vec_d cmt: polynomial vector with the signed commitment
*              - const sep_pk_t *pk: pointer to signer public key structure
* 
* Returns 1 if signature could be verified correctly and 0 otherwise
**************************************************/
int pre_sig_verify(poly_q_vec_d v11, const poly_q tag, const pre_sig_t *pre_sig, const poly_q_vec_d cmt, const pk_t *pk) {
  return pre_sig_verify_from_commitment(v11, tag, pre_sig, cmt, pk);
}

/*************************************************
* Name:        complete_decompose
*
* Description: Merge commitment randomness with pre signature and decompose
*
* Arguments:   - bsig_t *bsig: pointer to blind signature structure (to host w_{i,L}) (initialized)
*              - poly_q_vec_d *w1H: array of polynomial vectors to host w_{1,H} (initialized)
*              - poly_q_vec_d *w2H: array of polynomial vectors to host w_{2,H} (initialized)
*              - poly_q_vec_k w3H: polynomial vector to host w_{3,H} (initialized)
*              - const poly_q_vec_d v11: polynomial vector with the reconstructed v11
*              - const pre_sig_t *pre_sig: pointer to pre signature structure 
*              - const rand_t *rand: pointer to commitment randomness structure
**************************************************/
void complete_decompose(
    bsig_t *bsig,
    poly_q_vec_d w1H[2],
    poly_q_vec_d w2H[PARAM_K],
    poly_q_vec_k w3H,
    const poly_q_vec_d v11,
    const pre_sig_t *pre_sig, 
    const rand_t *rand) {
  size_t i;
  poly_q_vec_d tmp_L, tmp_H;
  // init vectors
  poly_q_vec_d_init(tmp_L);
  poly_q_vec_d_init(tmp_H);
  
  // computing w1L[0] and w1H[0]
  poly_q_vec_d_decompose_b1(tmp_H, tmp_L, rand->r1[0]);
  poly_q_vec_d_sub(tmp_L, v11, tmp_L);
  poly_q_vec_d_decompose_b1(w1H[0], bsig->w1L[0], tmp_L);
  poly_q_vec_d_sub(w1H[0], w1H[0], tmp_H);
  // computing w1L[1] and w1H[1]
  poly_q_vec_d_decompose_b1(tmp_H, tmp_L, rand->r1[1]);
  poly_q_vec_d_sub(tmp_L, pre_sig->v12, tmp_L);
  poly_q_vec_d_decompose_b1(w1H[1], bsig->w1L[1], tmp_L);
  poly_q_vec_d_sub(w1H[1], w1H[1], tmp_H);

  // computing w2L and w2H
  for (i = 0; i < PARAM_K; i++) {
    poly_q_vec_d_sub(tmp_L, pre_sig->v2[i], rand->r2[i]);
    poly_q_vec_d_decompose_b2(w2H[i], bsig->w2L[i], tmp_L);
  }

  // computing w3L and w3H
  poly_q_vec_k_sub(w3H, pre_sig->v3, rand->r3);
  poly_q_vec_k_decompose_b2(w3H, bsig->w3L, w3H);
  
  // clean up vectors
  poly_q_vec_d_clear(tmp_L);
  poly_q_vec_d_clear(tmp_H);
}

/*************************************************
* Name:        embed_2
*
* Description: Embedding the user relation for the issuance proof
*              of commitment opening and verifiable encryption
*
* Arguments:   - poly_qshow_mat_k_k *A_embed: array of polynomial matrices to host subring embedding of q1.b1.A'
*              - poly_qshow_mat_k_k *B_embed: array of polynomial matrices to host subring embedding of q1.b2.B
*              - poly_qshow_mat_k_k *A3_embed: array of polynomial matrices to host subring embedding of q1.b2.A3
*              - poly_qshow_mat_k_k *G_embed: array of polynomial matrices to host subring embedding of q1.G.w_{2,L}
*              - poly_qshow_vec_k *u_embed: array of polynomial vectors to host subring embedding of q1.(u + d.H(m) - [I|A'].w_{1,L} + B.w_{2,L} - A3.w_{3,L})
*              - poly_qshow_vec_m1 s1: polynomial vector to host subring embedding of witness
*              - const pk_t *pk: pointer to signer public key structure
*              - const poly_q tag: polynomial for tag
*              - const bsig_t *bsig: pointer to blind signature structure
*              - const poly_q_vec_d *w1H: array of polynomial vectors for w_{1,H}
*              - const poly_q_vec_d *w2H: array of polynomial vectors for w_{2,H}
*              - const poly_q_vec_k w3H: polynomial vector for w_{3,H}
*              - const uint8_t *msg: pointer to input message byte array (allocated PARAM_N/8 bytes)
**************************************************/
void embed_2(
    poly_qshow_mat_k_k A_embed[PARAM_D][PARAM_D], 
    poly_qshow_mat_k_k B_embed[PARAM_D][PARAM_K*PARAM_D], 
    poly_qshow_mat_k_k A3_embed[PARAM_D][PARAM_K], 
    poly_qshow_mat_k_k G_embed[PARAM_D], 
    poly_qshow_vec_k u_embed[PARAM_D], 
    poly_qshow_vec_m1 s1, 
    const pk_t *pk,
    const poly_q tag,
    const bsig_t *bsig, 
    const poly_q_vec_d w1H[2],
    const poly_q_vec_d w2H[PARAM_K],
    const poly_q_vec_k w3H, 
    const uint8_t msg[PARAM_N/8]) {
  uint8_t random_sign[1];
  size_t i,j;
  uint64_t norm2sq, four_squares_res[4];
  coeff_q tmp_coeff;
  poly_q tmp;
  poly_q_mat_d_d A;
  poly_q_mat_d_k A3;
  poly_q_vec_d d, u;
  poly_qshow_vec_k tmp_vec_k;

  // init matrices and polynomials
  poly_q_init(tmp);
  poly_q_mat_d_d_init(A);
  poly_q_mat_d_k_init(A3);
  poly_q_vec_d_init(d);
  poly_q_vec_d_init(u);
  poly_qshow_vec_k_init(tmp_vec_k);

  // computing square norms and four-square decompositions
  // B_1'^2 - |w1H|^2
  norm2sq = poly_q_vec_d_norm2(w1H[0]);
  norm2sq += poly_q_vec_d_norm2(w1H[1]);
  assert(norm2sq <= PARAM_B1_PRIME_SQ);
  four_squares(four_squares_res, PARAM_B1_PRIME_SQ - norm2sq);
  for (i = 0; i < 4; i++) {
    poly_qshow_set_coeff(s1->entries[IDX_W23_SHOW - 1], i, (coeff_qshow)(four_squares_res[i]));
  }
  // B_2'^2 - |[w2H,w3H]|^2
  norm2sq = poly_q_vec_k_norm2(w3H);
  for (i = 0; i < PARAM_K; i++) {
    norm2sq += poly_q_vec_d_norm2(w2H[i]);
  }
  assert(norm2sq <= PARAM_B2_PRIME_SQ);
  four_squares(four_squares_res, PARAM_B2_PRIME_SQ - norm2sq);
  for (i = 0; i < 4; i++) {
    poly_qshow_set_coeff(s1->entries[IDX_TAG_SHOW - 1], i, (coeff_qshow)(four_squares_res[i]));
  }

  // embedding witness vector s1 = [theta(w1H)|a1 | theta(wH23)|a23 | theta(tag) | 1] and multiplying by random sign b
  // w1H
  for (i = 0; i < 2*PARAM_D; i++) {
    poly_qshow_subring_embed_vec_k(tmp_vec_k, w1H[i / PARAM_D]->entries[i % PARAM_D], 1);
    for (j = 0; j < PARAM_K_SHOW; j++) {
      poly_qshow_set(s1->entries[IDX_W1_SHOW + i*PARAM_K_SHOW + j], tmp_vec_k->entries[j]);
    }
  }
  // w2H
  for (i = 0; i < PARAM_D*PARAM_K; i++) {
    poly_qshow_subring_embed_vec_k(tmp_vec_k, w2H[i / PARAM_D]->entries[i % PARAM_D], 1);
    for (j = 0; j < PARAM_K_SHOW; j++) {
      poly_qshow_set(s1->entries[IDX_W23_SHOW + i*PARAM_K_SHOW + j], tmp_vec_k->entries[j]);
    }
  } 
  // w3H
  for (i = 0; i < PARAM_K; i++) {
    poly_qshow_subring_embed_vec_k(tmp_vec_k, w3H->entries[i], 1);
    for (j = 0; j < PARAM_K_SHOW; j++) {
      poly_qshow_set(s1->entries[IDX_W3_SHOW + i*PARAM_K_SHOW + j], tmp_vec_k->entries[j]);
    }
  }
  // tag
  poly_qshow_subring_embed_vec_k(tmp_vec_k, tag, 1);
  for (j = 0; j < PARAM_K_SHOW; j++) {
    poly_qshow_set(s1->entries[IDX_TAG_SHOW + j], tmp_vec_k->entries[j]);
  }
  // b
  poly_qshow_zero(s1->entries[IDX_B_SHOW]);
  poly_qshow_set_coeff(s1->entries[IDX_B_SHOW], 0, 1);

  // sampling random sign b
  randombytes(random_sign, 1);
  tmp_coeff = (coeff_qshow)(((uint64_t)random_sign[0] >> 7) << 1) - 1;
  poly_qshow_vec_m1_mul_scalar(s1, s1, tmp_coeff); // multiply by b

  // expanding uniform public parameters
  poly_q_mat_d_d_uniform(A, pk->seed, DOMAIN_SEPARATOR_A, 0);
  poly_q_mat_d_k_uniform(A3, pk->seed, DOMAIN_SEPARATOR_A3);
  poly_q_vec_d_uniform(d, pk->seed, DOMAIN_SEPARATOR_D);
  poly_q_vec_d_uniform(u, pk->seed, DOMAIN_SEPARATOR_U);

  // embedding commitment u = q1 * theta(u + d.H(m) - [I|A'].w_{1,L} + B.w_{2,L} - A3.w_{3,L})
  poly_q_from_bits(tmp, msg);
  poly_q_vec_d_mul_poly(d, d, tmp);
  poly_q_vec_d_add(u, u, d); // d can be used as temp variable
  poly_q_vec_d_sub(u, u, bsig->w1L[0]);
  poly_q_mat_d_d_mul_vec_d(d, A, bsig->w1L[1]);
  poly_q_vec_d_sub(u, u, d); // u contains u + d.H(m) - [I | A'].w_{1,L}
  for (i = 0; i < PARAM_K; i++) {
    poly_q_mat_d_d_mul_vec_d(d, pk->B[i], bsig->w2L[i]);
    poly_q_vec_d_add(u, u, d);
  }
  poly_q_mat_d_k_mul_vec_k(d, A3, bsig->w3L);
  poly_q_vec_d_sub(u, u, d);
  for (i = 0; i < PARAM_D; i++) {
    poly_qshow_subring_embed_vec_k(u_embed[i], u->entries[i], PARAM_Q1_SHOW);
  }

  // embedding matrices
  // computing G.w_{2,L}
  poly_q_vec_d_set(d, bsig->w2L[0]);
  tmp_coeff = PARAM_B;
  for (i = 1; i < PARAM_K; i++) {
    poly_q_vec_d_mul_scalar(u, bsig->w2L[i], tmp_coeff);
    poly_q_vec_d_add(d, d, u);
    tmp_coeff *= PARAM_B;
  }
  // A, G.w_{2,L}
  for (i = 0; i < PARAM_D; i++) {
    for (j = 0; j < PARAM_D; j++) {
      poly_qshow_subring_embed_mat_k_k(A_embed[i][j], A->rows[i]->entries[j], PARAM_B1*PARAM_Q1_SHOW); // A'
    }
    poly_qshow_subring_embed_mat_k_k(G_embed[i], d->entries[i], PARAM_Q1_SHOW); // G.w_{2,L}
  }
  // B
  for (i = 0; i < PARAM_D; i++) {
    for (j = 0; j < PARAM_D*PARAM_K; j++) {
      poly_qshow_subring_embed_mat_k_k(B_embed[i][j], pk->B[j / PARAM_D]->rows[i]->entries[j % PARAM_D], PARAM_B2*PARAM_Q1_SHOW); // B
    }
  }
  // A3
  for (i = 0; i < PARAM_D; i++) {
    for (j = 0; j < PARAM_K; j++) {
      poly_qshow_subring_embed_mat_k_k(A3_embed[i][j], A3->rows[i]->entries[j], PARAM_B2*PARAM_Q1_SHOW); 
    }
  }

  // clean up matrices and polynomials
  poly_q_clear(tmp);
  poly_q_mat_d_d_clear(A);
  poly_q_mat_d_k_clear(A3);
  poly_q_vec_d_clear(d);
  poly_q_vec_d_clear(u);
  poly_qshow_vec_k_clear(tmp_vec_k);
}

/*************************************************
* Name:        embed_2_verifier
*
* Description: Embedding the user relation for the issuance proof
*              of commitment opening and verifiable encryption
*
* Arguments:   - poly_qshow_mat_k_k *A_embed: array of polynomial matrices to host subring embedding of q1.b1.A'
*              - poly_qshow_mat_k_k *B_embed: array of polynomial matrices to host subring embedding of q1.b2.B
*              - poly_qshow_mat_k_k *A3_embed: array of polynomial matrices to host subring embedding of q1.b2.A3
*              - poly_qshow_mat_k_k *G_embed: array of polynomial matrices to host subring embedding of q1.G.w_{2,L}
*              - poly_qshow_vec_k *u_embed: array of polynomial vectors to host subring embedding of q1.(u + d.H(m) - [I|A'].w_{1,L} + B.w_{2,L} - A3.w_{3,L})
*              - const pk_t *pk: pointer to signer public key structure
*              - const bsig_t *bsig: pointer to blind signature structure
*              - const uint8_t *msg: pointer to input message byte array (allocated PARAM_N/8 bytes)
**************************************************/
void embed_2_verifier(
    poly_qshow_mat_k_k A_embed[PARAM_D][PARAM_D], 
    poly_qshow_mat_k_k B_embed[PARAM_D][PARAM_K*PARAM_D], 
    poly_qshow_mat_k_k A3_embed[PARAM_D][PARAM_K], 
    poly_qshow_mat_k_k G_embed[PARAM_D], 
    poly_qshow_vec_k u_embed[PARAM_D], 
    const pk_t *pk,
    const bsig_t *bsig, 
    const uint8_t msg[PARAM_N/8]) {
  size_t i,j;
  coeff_q tmp_coeff;
  poly_q tmp;
  poly_q_mat_d_d A;
  poly_q_mat_d_k A3;
  poly_q_vec_d d, u;
  poly_qshow_vec_k tmp_vec_k;

  // init matrices and polynomials
  poly_q_init(tmp);
  poly_q_mat_d_d_init(A);
  poly_q_mat_d_k_init(A3);
  poly_q_vec_d_init(d);
  poly_q_vec_d_init(u);
  poly_qshow_vec_k_init(tmp_vec_k);

  // expanding uniform public parameters
  poly_q_mat_d_d_uniform(A, pk->seed, DOMAIN_SEPARATOR_A, 0);
  poly_q_mat_d_k_uniform(A3, pk->seed, DOMAIN_SEPARATOR_A3);
  poly_q_vec_d_uniform(d, pk->seed, DOMAIN_SEPARATOR_D);
  poly_q_vec_d_uniform(u, pk->seed, DOMAIN_SEPARATOR_U);

  // embedding commitment u = q1 * theta(u + d.H(m) - [I|A'].w_{1,L} + B.w_{2,L} - A3.w_{3,L})
  poly_q_from_bits(tmp, msg);
  poly_q_vec_d_mul_poly(d, d, tmp);
  poly_q_vec_d_add(u, u, d); // d can be used as temp variable
  poly_q_vec_d_sub(u, u, bsig->w1L[0]);
  poly_q_mat_d_d_mul_vec_d(d, A, bsig->w1L[1]);
  poly_q_vec_d_sub(u, u, d); // u contains u + d.H(m) - [I | A'].w_{1,L}
  for (i = 0; i < PARAM_K; i++) {
    poly_q_mat_d_d_mul_vec_d(d, pk->B[i], bsig->w2L[i]);
    poly_q_vec_d_add(u, u, d);
  }
  poly_q_mat_d_k_mul_vec_k(d, A3, bsig->w3L);
  poly_q_vec_d_sub(u, u, d);
  for (i = 0; i < PARAM_D; i++) {
    poly_qshow_subring_embed_vec_k(u_embed[i], u->entries[i], PARAM_Q1_SHOW);
  }

  // embedding matrices
  // computing G.w_{2,L}
  poly_q_vec_d_set(d, bsig->w2L[0]);
  tmp_coeff = PARAM_B;
  for (i = 1; i < PARAM_K; i++) {
    poly_q_vec_d_mul_scalar(u, bsig->w2L[i], tmp_coeff);
    poly_q_vec_d_add(d, d, u);
    tmp_coeff *= PARAM_B;
  }
  // A, G.w_{2,L}
  for (i = 0; i < PARAM_D; i++) {
    for (j = 0; j < PARAM_D; j++) {
      poly_qshow_subring_embed_mat_k_k(A_embed[i][j], A->rows[i]->entries[j], PARAM_B1*PARAM_Q1_SHOW); // A'
    }
    poly_qshow_subring_embed_mat_k_k(G_embed[i], d->entries[i], PARAM_Q1_SHOW); // G.w_{2,L}
  }
  // B
  for (i = 0; i < PARAM_D; i++) {
    for (j = 0; j < PARAM_D*PARAM_K; j++) {
      poly_qshow_subring_embed_mat_k_k(B_embed[i][j], pk->B[j / PARAM_D]->rows[i]->entries[j % PARAM_D], PARAM_B2*PARAM_Q1_SHOW); // B
    }
  }
  // A3
  for (i = 0; i < PARAM_D; i++) {
    for (j = 0; j < PARAM_K; j++) {
      poly_qshow_subring_embed_mat_k_k(A3_embed[i][j], A3->rows[i]->entries[j], PARAM_B2*PARAM_Q1_SHOW); 
    }
  }

  // clean up matrices and polynomials
  poly_q_clear(tmp);
  poly_q_mat_d_d_clear(A);
  poly_q_mat_d_k_clear(A3);
  poly_q_vec_d_clear(d);
  poly_q_vec_d_clear(u);
  poly_qshow_vec_k_clear(tmp_vec_k);
}