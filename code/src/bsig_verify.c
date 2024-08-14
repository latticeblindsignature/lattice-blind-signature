#include "arith.h"
#include "randombytes.h"
#include "poly_qshow_sampling.h"
#include "bsig_signer.h"
#include "bsig_user.h"
#include "bsig_verify.h"
#include "macros.h"
#include "fips202.h"

/*************************************************
* Name:        bsig_verify
*
* Description: Blind signature verification
*
* Arguments:   - const bsig_t *bsig: pointer to blind signature structure
*              - const pk_t *pk: pointer to signature public key structure
*              - const uint8_t *msg: pointer to input message byte array (allocated PARAM_N/8 bytes)
**************************************************/
int bsig_verify(const bsig_t *bsig, const pk_t *pk, const uint8_t msg[PARAM_N/8]) {
  int is_valid;
  size_t i,j;
  int64_t min, max, tmpmin, tmpmax;
  poly_qshow_mat_k_k A_embed[PARAM_D][PARAM_D], B_embed[PARAM_D][PARAM_K*PARAM_D], A3_embed[PARAM_D][PARAM_K], G_embed[PARAM_D];
  poly_qshow_vec_k u_embed[PARAM_D];

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

  // checking w1L
  max = poly_q_vec_d_minmax(&min, bsig->w1L[0]);
  tmpmax = poly_q_vec_d_minmax(&tmpmin, bsig->w1L[1]);
  max = (tmpmax > max) ? tmpmax : max;
  min = (tmpmin < min) ? tmpmin : min;
  is_valid = (min >= -PARAM_B1) && (max < PARAM_B1);
  
  // checking [w2L | w3L]
  max = poly_q_vec_k_minmax(&min, bsig->w3L);
  for (int k = 0; k < PARAM_K; k++) {
    tmpmax = poly_q_vec_d_minmax(&tmpmin, bsig->w2L[k]);
    max = (tmpmax > max) ? tmpmax : max;
    min = (tmpmin < min) ? tmpmin : min;
  }
  is_valid = is_valid && (min >= -PARAM_B2) && (max < PARAM_B2);

  embed_2_verifier(A_embed, B_embed, A3_embed, G_embed, u_embed, pk, bsig, msg);
  is_valid = is_valid && verify_2(&bsig->proof_2, A_embed, B_embed, A3_embed, G_embed, u_embed, bsig->crs_seed, pk->seed);

  // clean up
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

  return is_valid;
}