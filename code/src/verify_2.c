#include "arith.h"
#include "randombytes.h"
#include "poly_qshow_sampling.h"
#include "bsig_signer.h"
#include "bsig_user.h"
#include "macros.h"
#include "fips202.h"

/*************************************************
* Name:        verify_2
*
* Description: Verification of the zero-knowledge proof included in the blind signature
*
* Arguments:   - const proof_2_t *proof_2: pointer to show proof structure
*              - const poly_qshow_mat_k_k *A_embed: array of polynomial matrices to host subring embedding of q1.b1.A'
*              - const poly_qshow_mat_k_k *B_embed: array of polynomial matrices to host subring embedding of q1.b2.B
*              - const poly_qshow_mat_k_k *A3_embed: array of polynomial matrices to host subring embedding of q1.b2.A3
*              - const poly_qshow_mat_k_k *G_embed: array of polynomial matrices to host subring embedding of q1.G.w_{2,L}
*              - const poly_qshow_vec_k *u_embed: array of polynomial vectors to host subring embedding of q1.(u + d.H(m) - [I|A'].w_{1,L} + B.w_{2,L} - A3.w_{3,L})
*              - const uint8_t *crs_seed: pointer to byte array containing the CRS seed (allocated SEED_BYTES bytes)
*              - const uint8_t *seed: pointer to byte array containing the seed for public parameters (allocated SEED_BYTES bytes)
* 
* Returns 1 if proof_2 could be verified correctly and 0 otherwise
**************************************************/
int verify_2(
    const proof_2_t          *proof_2, 
    const poly_qshow_mat_k_k A_embed[PARAM_D][PARAM_D], 
    const poly_qshow_mat_k_k B_embed[PARAM_D][PARAM_K*PARAM_D], 
    const poly_qshow_mat_k_k A3_embed[PARAM_D][PARAM_K], 
    const poly_qshow_mat_k_k G_embed[PARAM_D], 
    const poly_qshow_vec_k   u_embed[PARAM_D], 
    const uint8_t            crs_seed[CRS_SEED_BYTES], 
    const uint8_t            seed[SEED_BYTES]) {
  size_t i,j,k,i_k_quot,i_k_rem;
  int is_valid = 1;
  uint8_t buf[CHAL4_SHOW_INPUT_BYTES] = {0}, challenge_seed[SEED_BYTES];
  uint128 sq_norm_z2;
  uint64_t sq_norm_z1, sq_norm_z3;
  coeff_qshow tmp_coeff;
  poly_qshow tmp_poly, t0, sum_mu_gamma[4], c, z1s_z1, z1i_star, sum_gamma_X_star_ij_1, sum_gamma_X_star_ij_2;
  poly_qshow_vec_m2_d b;
  poly_qshow_vec_d tmp_vec_d, w_H;
  poly_qshow_vec_l c_i, c_i_prime, d_i;
  poly_qshow_vec_256_l tmp_vec_256_l, Z;
  poly_qshow_vec_k Gz1_w2, tmp_vec_k, sum_vec_k;
  poly_qshow_mat_k_k chal_3_quad_matrix;
  poly_qshow_mat_d_m1 A1;
  poly_qshow_mat_d_m2_d A2;
  poly_qshow_mat_256l_m2_d Byg;
  // challenges
  poly_qshow_vec_m1 chal_1[PARAM_ARP_SHOW];
  coeff_qshow chal_2[2*PARAM_L_SHOW][PARAM_ARP_SHOW + 4 + PARAM_N_SHOW - 1];
  poly_qshow_vec_l chal_3_l;
  poly_qshow_vec_k chal_3_dk[PARAM_D];
  poly_qshow chal_3_1;

  // init
  // init polynomials
  poly_qshow_init(tmp_poly);
  poly_qshow_init(t0);
  poly_qshow_init(c);
  poly_qshow_init(z1s_z1);
  poly_qshow_init(z1i_star);
  poly_qshow_init(sum_gamma_X_star_ij_1);
  poly_qshow_init(sum_gamma_X_star_ij_2);
  for (i = 0; i < 4; i++) {
    poly_qshow_init(sum_mu_gamma[i]);
    poly_qshow_zero(sum_mu_gamma[i]);
  }
  // init vectors and matrices
  poly_qshow_vec_m2_d_init(b);
  poly_qshow_vec_d_init(tmp_vec_d);
  poly_qshow_vec_d_init(w_H);
  poly_qshow_vec_l_init(c_i);
  poly_qshow_vec_l_init(c_i_prime);
  poly_qshow_vec_l_init(d_i);
  poly_qshow_vec_256_l_init(tmp_vec_256_l);
  poly_qshow_vec_256_l_init(Z);
  poly_qshow_vec_k_init(Gz1_w2);
  poly_qshow_vec_k_init(tmp_vec_k);
  poly_qshow_vec_k_init(sum_vec_k);
  poly_qshow_mat_k_k_init(chal_3_quad_matrix);
  poly_qshow_mat_d_m1_init(A1);
  poly_qshow_mat_d_m2_d_init(A2);
  poly_qshow_mat_256l_m2_d_init(Byg);
   // init challenges
  for (i = 0; i < PARAM_ARP_SHOW; i++) {
    poly_qshow_vec_m1_init(chal_1[i]);
  }
  poly_qshow_vec_l_init(chal_3_l);
  for (i = 0; i < PARAM_D; i++) {
    poly_qshow_vec_k_init(chal_3_dk[i]);
  }
  poly_qshow_init(chal_3_1);

  // expanding CRS
  poly_qshow_mat_d_m1_uniform(A1, crs_seed, DOMAIN_SEPARATOR_A1_SHOW);
  poly_qshow_mat_d_m2_d_uniform(A2, crs_seed, DOMAIN_SEPARATOR_A2_SHOW);
  poly_qshow_mat_256l_m2_d_uniform(Byg, crs_seed, DOMAIN_SEPARATOR_BYG_SHOW);
  poly_qshow_vec_m2_d_uniform(b, crs_seed, DOMAIN_SEPARATOR_B_SHOW);

  // byte-packing the statement
  // public parameters are derived from seed
  for (i = 0; i < CRS_SEED_BYTES; i++) {
    buf[1 + i] = crs_seed[i];
  }
  for (i = 0; i < SEED_BYTES; i++) {
    buf[1 + CRS_SEED_BYTES + i] = seed[i];
  }
  for (i = 0; i < PARAM_D; i++) {
    poly_qshow_vec_k_pack(buf + 1 + CRS_SEED_BYTES + SEED_BYTES + i*POLYQSHOW_VECK_PACKEDBYTES, u_embed[i]);
  }
  ////////////////////////////////////// NO INITIALIZATIONS BELOW THIS LINE

  /****** reconstructing first message ******/
  // recomputing tmp = A1.z1 + A2.z2_1 - c.2^D.tA_H 
  poly_qshow_mat_d_m1_mul_vec_m1(tmp_vec_d, A1, proof_2->z1);
  poly_qshow_mat_d_m2_d_mul_vec_m2_d(w_H, A2, proof_2->z2_1);  // using w_H as temp variable
  poly_qshow_vec_d_add(tmp_vec_d, tmp_vec_d, w_H);
  poly_qshow_vec_d_mul_poly_qshow(w_H, proof_2->tA_H, proof_2->c); // using w_H as temp variable
  poly_qshow_vec_d_mul_scalar(w_H, w_H, (coeff_qshow)(1 << PARAM_D_ROUND_SHOW));
  poly_qshow_vec_d_sub(tmp_vec_d, tmp_vec_d, w_H);
  // using hint
  poly_qshow_vec_d_useGhint(w_H, proof_2->hint, tmp_vec_d); 

  // computing first challenge with w_H so that we can use w_H as temp variable afterwards
  buf[0] = 1;
  poly_qshow_vec_d_pack(buf + SHOW_CHALLENGE_BASE_BYTES, proof_2->tA_H);
  poly_qshow_vec_d_pack(buf + SHOW_CHALLENGE_BASE_BYTES + POLYQSHOW_VECD_PACKEDBYTES, w_H);
  poly_qshow_vec_256_l_pack(buf + SHOW_CHALLENGE_BASE_BYTES + 2*POLYQSHOW_VECD_PACKEDBYTES, proof_2->tB);
  // recomputing first challenge
  shake256(challenge_seed, SEED_BYTES, buf, CHAL1_SHOW_INPUT_BYTES);
  for (i = 0; i < PARAM_ARP_SHOW; i++) {
    poly_qshow_vec_m1_binomial(chal_1[i], challenge_seed, DOMAIN_SEPARATOR_CHAL1_SHOW, i, SEED_BYTES);
  }

  /****** performing all checks except c=c' for early reject ******/
  poly_qshow_vec_d_mul_scalar(w_H, w_H, PARAM_GAMMA_SHOW);
  poly_qshow_vec_d_sub(w_H, w_H, tmp_vec_d); // w_H contains z2_2
  sq_norm_z1 = poly_qshow_vec_m1_norm2(proof_2->z1);
  sq_norm_z2 = poly_qshow_vec_m2_d_norm2(proof_2->z2_1) + poly_qshow_vec_d_norm2(w_H);
  sq_norm_z3 = 0;
  for (i = 0; i < PARAM_ARP_SHOW; i++) {
    tmp_coeff = proof_2->z3[i];
    CHK_UI_OVF_ADDITION(sq_norm_z3, tmp_coeff * tmp_coeff);
  }
  tmp_coeff = poly_qshow_vec_d_norm_inf(proof_2->hint);
  is_valid = is_valid && (sq_norm_z1 <= PARAM_B1SQ_SHOW) && (sq_norm_z2 <= ((uint128)PARAM_B2SQ_SHOW_LOW64 + (((uint128)PARAM_B2SQ_SHOW_HIGH64) << 64))) && (sq_norm_z3 <= PARAM_B3SQ_SHOW);
  is_valid = is_valid && (tmp_coeff <= PARAM_Q_GAMMA_SHOW/2);
  for (i = 0; i < PARAM_L_SHOW; i++) {
    is_valid = is_valid && (poly_qshow_get_coeff(proof_2->f->entries[i], 0) == 0) && (poly_qshow_get_coeff(proof_2->f->entries[i], PARAM_N_SHOW/2) == 0);
  }
  if (!is_valid) {
    goto verify_2_cleanup; 
  }

  /****** reconstructing second message ******/
  buf[0] = 2;
  for (i = 0; i < PARAM_ARP_SHOW; i++) {
    coeff_qshow_pack(buf + CHAL1_SHOW_INPUT_BYTES + i*COEFFQSHOW_PACKEDBYTES, proof_2->z3[i]);
  }
  // recomputing second challenge
  shake256(challenge_seed, SEED_BYTES, buf, CHAL2_SHOW_INPUT_BYTES);
  for (i = 0; i < 2*PARAM_L_SHOW; i++) {
    vec_qshow_uniform(chal_2[i], challenge_seed, DOMAIN_SEPARATOR_CHAL2_SHOW, i, SEED_BYTES); // writes PARAM_ARP_SHOW + 4 + PARAM_N_SHOW - 1 uniform numbers to the first argument
  }

  /****** reconstructing third message ******/
  // computing third challenge
  buf[0] = 3;
  poly_qshow_vec_l_pack(buf + CHAL2_SHOW_INPUT_BYTES, proof_2->f);
  // recomputing third challenge
  shake256(challenge_seed, SEED_BYTES, buf, CHAL3_SHOW_INPUT_BYTES);
  poly_qshow_vec_l_1_uniform(chal_3_l, chal_3_1, buf, DOMAIN_SEPARATOR_CHAL3_SHOW, SEED_BYTES);
  for (i = 0; i < PARAM_D; i++) {
    poly_qshow_vec_k_uniform(chal_3_dk[i], challenge_seed, DOMAIN_SEPARATOR_CHAL3_SHOW, i+PARAM_L_SHOW+1, SEED_BYTES);
  }

  /****** reconstructing fourth message ******/
  /**********************************************
  * Recomputing commitment t0
  **********************************************/
  // computing <b, z2> - c.t1
  poly_qshow_vec_m2_d_mul_inner(t0, b, proof_2->z2_1);
  poly_qshow_mul(tmp_poly, proof_2->c, proof_2->t1);
  poly_qshow_sub(t0, t0, tmp_poly);

  // Computing sum_i mu_i (gamma_{2i,256+j} + x^(n_iss/2).gamma_{2i+1,256+j}), c_i, c_i' and d_i
  poly_qshow_vec_l_zero(d_i);
  for (i = 0; i < PARAM_L_SHOW; i++) {
    poly_qshow_mul_xj(c, chal_3_l->entries[i], PARAM_N_SHOW/2); // using c as temp variable
    for (j = 0; j < 4; j++) {
      poly_qshow_mul_scalar(tmp_poly, chal_3_l->entries[i], chal_2[2*i][PARAM_ARP_SHOW + j]); // mu_i gamma_{2i,256+j}
      poly_qshow_add(sum_mu_gamma[j], sum_mu_gamma[j], tmp_poly);
      poly_qshow_mul_scalar(tmp_poly, c, chal_2[2*i+1][PARAM_ARP_SHOW + j]); // mu_i.x^(n_iss/2).gamma_{2i+1,256+j}
      poly_qshow_add(sum_mu_gamma[j], sum_mu_gamma[j], tmp_poly);
    }
    // c_i and c_i'
    poly_qshow_set_coeff(c_i->entries[i], 0, -chal_2[2*i+1][PARAM_ARP_SHOW+3+PARAM_N_SHOW/2]);
    poly_qshow_set_coeff(c_i_prime->entries[i], 0, chal_2[2*i+1][PARAM_ARP_SHOW+3+PARAM_N_SHOW/2]);
    for (j = 1; j < PARAM_N_SHOW/2; j++) {
      poly_qshow_set_coeff(c_i->entries[i], j, chal_2[2*i][PARAM_ARP_SHOW+3+j] - chal_2[2*i+1][PARAM_ARP_SHOW+3+j+PARAM_N_SHOW/2]);
      poly_qshow_set_coeff(c_i_prime->entries[i], j, chal_2[2*i+1][PARAM_ARP_SHOW+3+PARAM_N_SHOW/2-j] - chal_2[2*i][PARAM_ARP_SHOW+3+PARAM_N_SHOW-j]);
    }
    poly_qshow_set_coeff(c_i->entries[i], PARAM_N_SHOW/2, chal_2[2*i][PARAM_ARP_SHOW+3+PARAM_N_SHOW/2]);
    poly_qshow_set_coeff(c_i_prime->entries[i], PARAM_N_SHOW/2, -chal_2[2*i][PARAM_ARP_SHOW+3+PARAM_N_SHOW/2]);
    for (j = PARAM_N_SHOW/2 + 1; j < PARAM_N_SHOW; j++) {
      poly_qshow_set_coeff(c_i->entries[i], j, chal_2[2*i][PARAM_ARP_SHOW+3+j] + chal_2[2*i+1][PARAM_ARP_SHOW+3+j-PARAM_N_SHOW/2]);
      poly_qshow_set_coeff(c_i_prime->entries[i], j, -(chal_2[2*i][PARAM_ARP_SHOW+3+PARAM_N_SHOW-j] + chal_2[2*i+1][PARAM_ARP_SHOW+3+3*PARAM_N_SHOW/2-j]));
    }
    poly_qshow_mul_scalar(c_i->entries[i], c_i->entries[i], PARAM_TWO_INVMOD_Q_SHOW);
    poly_qshow_mul_scalar(c_i_prime->entries[i], c_i_prime->entries[i], PARAM_TWO_INVMOD_Q_SHOW);
    // d_i
    poly_qshow_zero(tmp_poly);
    poly_qshow_muladd_constant(tmp_poly, chal_2[2*i+1][PARAM_ARP_SHOW    ], PARAM_B1_PRIME_SQ);
    poly_qshow_muladd_constant(tmp_poly, chal_2[2*i+1][PARAM_ARP_SHOW + 1], PARAM_B2_PRIME_SQ);
    poly_qshow_muladd_constant(tmp_poly, chal_2[2*i+1][PARAM_ARP_SHOW + 2], PARAM_W);
    for (j = 0; j < PARAM_ARP_SHOW; j++) {
      poly_qshow_muladd_constant(tmp_poly, chal_2[2*i+1][j], proof_2->z3[j]);
    }
    poly_qshow_mul_xj(d_i->entries[i], tmp_poly, PARAM_N_SHOW/2);
    poly_qshow_muladd_constant(d_i->entries[i], chal_2[2*i][PARAM_ARP_SHOW    ], PARAM_B1_PRIME_SQ);
    poly_qshow_muladd_constant(d_i->entries[i], chal_2[2*i][PARAM_ARP_SHOW + 1], PARAM_B2_PRIME_SQ);
    poly_qshow_muladd_constant(d_i->entries[i], chal_2[2*i][PARAM_ARP_SHOW + 2], PARAM_W);
    for (j = 0; j < PARAM_ARP_SHOW; j++) {
      poly_qshow_muladd_constant(d_i->entries[i], chal_2[2*i][j], proof_2->z3[j]);
    }
  }
  poly_qshow_add(sum_mu_gamma[2], sum_mu_gamma[2], sum_mu_gamma[3]); // sum_mu_gamma[2] = sum_i mu_i(gamma_{2i,258}+gamma_{2i,259} + x^{n_show/2}.(gamma_{2i+1,258}+gamma_{2i+1,259}))
  poly_qshow_mul_scalar(sum_mu_gamma[3], sum_mu_gamma[3], PARAM_TWO_INVMOD_Q_SHOW);

  /**********************************************
  * Quadratic terms in t0:
  *  sum_mu_gamma[0].<z_{1,1}'*, z_{1,1}'> + sum_mu_gamma[1].<z_{1,23}'*, z_{1,23}'> 
  *    + sum_mu_gamma[2].<z_{1,t}*, z_{1,t}> - sum_mu_gamma[3].(<z_{1,t}*, z_{1,b}.one> + <z_{1,t}*, z_{1,b}.one>*)
  *    + z_{1,t}.chal_3_quad_matrix.z_{1,2} + sum_{i < dk} mu_{l+i}.z_{1,b}.([A_theta.z_{1,1} - B_theta.z_{1,2} + A3_theta.z_{1,3} + G_theta.z_{1,t}]_i)
  *    + mu_{l+dk+1}.z_{1,b}^2
  **********************************************/
  // {1,1}
  poly_qshow_zero(z1s_z1);
  for (i = IDX_W1_SHOW; i < IDX_W23_SHOW; i++) {
    poly_qshow_conjugate(z1i_star, proof_2->z1->entries[i]);
    poly_qshow_mul(tmp_poly, z1i_star, proof_2->z1->entries[i]);
    poly_qshow_add(z1s_z1, z1s_z1, tmp_poly);
  }
  poly_qshow_mul(tmp_poly, sum_mu_gamma[0], z1s_z1);
  poly_qshow_add(t0, t0, tmp_poly);
  // <z_{1,23}'*,z_{1,23}'>
  poly_qshow_zero(z1s_z1);
  for (i = IDX_W23_SHOW; i < IDX_TAG_SHOW; i++) {
    poly_qshow_conjugate(z1i_star, proof_2->z1->entries[i]);
    poly_qshow_mul(tmp_poly, z1i_star, proof_2->z1->entries[i]);
    poly_qshow_add(z1s_z1, z1s_z1, tmp_poly);
  }
  poly_qshow_mul(tmp_poly, z1s_z1, sum_mu_gamma[1]);
  poly_qshow_add(t0, t0, tmp_poly);
  // <z_{1,t}*,z_{1,t}>
  poly_qshow_zero(z1s_z1);
  for (i = IDX_TAG_SHOW; i < IDX_B_SHOW; i++) {
    poly_qshow_conjugate(z1i_star, proof_2->z1->entries[i]);
    poly_qshow_mul(tmp_poly, z1i_star, proof_2->z1->entries[i]);
    poly_qshow_add(z1s_z1, z1s_z1, tmp_poly);
  }
  poly_qshow_mul(tmp_poly, z1s_z1, sum_mu_gamma[2]);
  poly_qshow_add(t0, t0, tmp_poly);
  // one term
  for (i = 0; i < PARAM_N_SHOW; i++) {
    poly_qshow_set_coeff(c, i, 1); // poly with all one // using c as temp variable
  }
  poly_qshow_zero(z1s_z1);
  poly_qshow_mul(c, proof_2->z1->entries[IDX_B_SHOW], c); // using c as temp variable
  for (i = IDX_TAG_SHOW; i < IDX_B_SHOW; i++) {
    poly_qshow_conjugate(z1i_star, proof_2->z1->entries[i]);
    poly_qshow_mul(tmp_poly, z1i_star, c);
    poly_qshow_add(z1s_z1, z1s_z1, tmp_poly);
  }
  poly_qshow_conjugate(tmp_poly, z1s_z1);
  poly_qshow_add(z1s_z1, z1s_z1, tmp_poly);
  poly_qshow_mul(tmp_poly, z1s_z1, sum_mu_gamma[3]);
  poly_qshow_sub(t0, t0, tmp_poly);
  // {1,b}
  poly_qshow_mul(tmp_poly, chal_3_1, proof_2->z1->entries[IDX_B_SHOW]);
  poly_qshow_mul(tmp_poly, tmp_poly, proof_2->z1->entries[IDX_B_SHOW]); 
  poly_qshow_add(t0, t0, tmp_poly);

  // quadratic part depending on chal_3_quad_matrix
  // compute gadget quadratic matrix necessary to compute sum_i mu_{l + i} G_i"
  poly_qshow_vec_k_zero(sum_vec_k);
  for (i = 0; i < PARAM_D; i++) {
    poly_qshow_mat_k_k_chal_3_embed(chal_3_quad_matrix, chal_3_dk[i]);
    for (j = 0; j < PARAM_K_SHOW; j++) {
      poly_qshow_zero(Gz1_w2->entries[j]);
      tmp_coeff = PARAM_B2*PARAM_Q1_SHOW;
      for (k = 0; k < PARAM_K; k++) {
        poly_qshow_mul_scalar(tmp_poly, proof_2->z1->entries[IDX_W23_SHOW + k*PARAM_D*PARAM_K_SHOW + i*PARAM_K_SHOW + j], tmp_coeff);
        poly_qshow_add(Gz1_w2->entries[j], Gz1_w2->entries[j], tmp_poly);
        tmp_coeff *= PARAM_B;
      }
    }
    poly_qshow_mat_k_k_mul_vec_k(tmp_vec_k, chal_3_quad_matrix, Gz1_w2);
    poly_qshow_vec_k_add(sum_vec_k, sum_vec_k, tmp_vec_k); // sum_vec_k = (sum_i mu_{l+i}G_i").z1_{v2}
  }
  // <z1_t, (sum_i mu_{l+i}G_i").z1_{v2}>
  for (i = 0; i < PARAM_K_SHOW; i++) {
    poly_qshow_mul(tmp_poly, proof_2->z1->entries[IDX_TAG_SHOW + i], sum_vec_k->entries[i]);
    poly_qshow_add(t0, t0, tmp_poly);
  }

  // part depending on embedded matrices A_embed, B_embed, A3_embed, G_embed
  // storing sum_{i < dk} mu_{l+i}.[A_theta.y_{1,1} - B_theta.y_{1,2} + A3_theta.y_{1,3} + G_theta.y_{1,t}]_i in z1s_z1 before multiplying by z_{1,b}
  poly_qshow_zero(z1s_z1);
  for (i = 0; i < PARAM_D*PARAM_K_SHOW; i++) {
    i_k_quot = i / PARAM_K_SHOW;
    i_k_rem = i % PARAM_K_SHOW;
    // z1i_star used as temp variable to host [A_theta.z_{1,1} - B_theta.z_{1,2} + A3_theta.z_{1,3} + G_theta.z_{1,t}]_i
    poly_qshow_mul_scalar(z1i_star, proof_2->z1->entries[i], PARAM_B1*PARAM_Q1_SHOW);
    for (j = 0; j < PARAM_D*PARAM_K_SHOW; j++) {
      poly_qshow_mul(tmp_poly, A_embed[i_k_quot][j / PARAM_K_SHOW]->rows[i_k_rem]->entries[j % PARAM_K_SHOW], proof_2->z1->entries[IDX_W12_SHOW + j]);
      poly_qshow_add(z1i_star, z1i_star, tmp_poly);
    }
    for (j = 0; j < PARAM_D*PARAM_K*PARAM_K_SHOW; j++) {
      poly_qshow_mul(tmp_poly, B_embed[i_k_quot][j / PARAM_K_SHOW]->rows[i_k_rem]->entries[j % PARAM_K_SHOW], proof_2->z1->entries[IDX_W23_SHOW + j]);
      poly_qshow_sub(z1i_star, z1i_star, tmp_poly); // substraction
    }
    for (j = 0; j < PARAM_K*PARAM_K_SHOW; j++) {
      poly_qshow_mul(tmp_poly, A3_embed[i_k_quot][j / PARAM_K_SHOW]->rows[i_k_rem]->entries[j % PARAM_K_SHOW], proof_2->z1->entries[IDX_W3_SHOW + j]);
      poly_qshow_add(z1i_star, z1i_star, tmp_poly);
    }
    for (j = 0; j < PARAM_K_SHOW; j++) {
      poly_qshow_mul(tmp_poly, G_embed[i_k_quot]->rows[i_k_rem]->entries[j], proof_2->z1->entries[IDX_TAG_SHOW + j]);
      poly_qshow_add(z1i_star, z1i_star, tmp_poly); 
    }
    poly_qshow_mul(z1i_star, z1i_star, chal_3_dk[i_k_quot]->entries[i_k_rem]);
    poly_qshow_add(z1s_z1, z1s_z1, z1i_star);
  }
  poly_qshow_mul(tmp_poly, z1s_z1, proof_2->z1->entries[IDX_B_SHOW]);
  poly_qshow_add(t0, t0, tmp_poly);

  /**********************************************
  * Linear terms in t0:
  *   Z = c.tB - Byg.z2
  * 
  *   c.sum_{i < l} mu_i.(Z_i + 2^{-1}.(SE_{2i} + x^{n_iss/2}.SE_{2i+1}).Z*_{:256/n_iss}
  *                    + 2^{-1}.(SE_{2i}* + x^{n_iss/2}.SE_{2i+1}*).Z_{:256/n_iss}
  *                    + 2^{-1}.(SR_{2i} + x^{n_iss/2}.SR_{2i+1}).z1* + 2^{-1}.(SR_{2i}* + x^{n_iss/2}.SR_{2i+1}*).z1
  *                    + c_i[i].z_{1,b}* + c_i_prime[i].z_{1,b})
  **********************************************/
  // computing Z = c.tB - Byg.z2
  poly_qshow_mat_256l_m2_d_mul_vec_m2_d(tmp_vec_256_l, Byg, proof_2->z2_1);
  poly_qshow_vec_256_l_mul_poly_qshow(Z, proof_2->tB, proof_2->c);
  poly_qshow_vec_256_l_sub(Z, Z, tmp_vec_256_l);

  poly_qshow_zero(sum_mu_gamma[0]); // using sum_mu_gamma[0] as temp variable to store sum before multiplying by c
  for (i = 0; i < PARAM_L_SHOW; i++) {
    // using z1s_z1 to store sum before multiplying by mu_i
    poly_qshow_mul(z1s_z1, c_i_prime->entries[i], proof_2->z1->entries[IDX_B_SHOW]);
    poly_qshow_conjugate(tmp_poly, proof_2->z1->entries[IDX_B_SHOW]);
    poly_qshow_mul(tmp_poly, tmp_poly, c_i->entries[i]);
    poly_qshow_add(z1s_z1, z1s_z1, tmp_poly);

    poly_qshow_add(z1s_z1, z1s_z1, Z->entries[PARAM_ARP_DIV_N_SHOW + i]);

    poly_qshow_zero(c); // c used as temp variable
    // term in SE_i
    for (j = 0; j < PARAM_ARP_DIV_N_SHOW; j++) {
      // [SE_{2i}*]_j
      poly_qshow_set_coeff(sum_gamma_X_star_ij_1, 0, chal_2[2*i][j * PARAM_N_SHOW]);
      for (k = 1; k < PARAM_N_SHOW; k++) {
        poly_qshow_set_coeff(sum_gamma_X_star_ij_1, k, - chal_2[2*i][(j + 1) * PARAM_N_SHOW - k]); // set conjugate directly
      }
      // [SE_{2i+1}*]_j
      poly_qshow_set_coeff(sum_gamma_X_star_ij_2, 0, chal_2[2*i+1][j * PARAM_N_SHOW]);
      for (k = 1; k < PARAM_N_SHOW; k++) {
        poly_qshow_set_coeff(sum_gamma_X_star_ij_2, k, - chal_2[2*i+1][(j + 1) * PARAM_N_SHOW - k]); // set conjugate directly
      }
      poly_qshow_mul_xj(tmp_poly, sum_gamma_X_star_ij_2, PARAM_N_SHOW/2); 
      poly_qshow_add(z1i_star, sum_gamma_X_star_ij_1, tmp_poly); // SE_{2i}* + x^{n_iss/2}.SE_{2i+1}*
      poly_qshow_sub(tmp_poly, sum_gamma_X_star_ij_1, tmp_poly); // SE_{2i}* - x^{n_iss/2}.SE_{2i+1}*
      poly_qshow_mul(z1i_star, z1i_star, Z->entries[j]); // (SE_{2i}* + x^{n_iss/2}.SE_{2i+1}*)_j . Z_j
      poly_qshow_add(c, c, z1i_star);
      poly_qshow_mul(tmp_poly, tmp_poly, Z->entries[j]); // (SE_{2i}* - x^{n_iss/2}.SE_{2i+1}*)_j . Z_j
      poly_qshow_conjugate(z1i_star, tmp_poly); // (SE_{2i} + x^{n_iss/2}.SE_{2i+1})_j . Z_j*
      poly_qshow_add(c, c, z1i_star);
    }
    // term in SR_i
    for (j = 0; j < PARAM_M1_SHOW; j++) {
      poly_qshow_zero(sum_gamma_X_star_ij_1);
      poly_qshow_zero(sum_gamma_X_star_ij_2);
      for (k = 0; k < PARAM_ARP_SHOW; k++) {
        poly_qshow_conjugate(z1i_star, chal_1[k]->entries[j]); // z1i_star used as temp variable
        // [SR_{2i,1}*]_j
        poly_qshow_mul_scalar(tmp_poly, z1i_star, chal_2[2*i][k]);
        poly_qshow_add(sum_gamma_X_star_ij_1, sum_gamma_X_star_ij_1, tmp_poly);
        // [SR_{2i+1,1}*]_j
        poly_qshow_mul_scalar(tmp_poly, z1i_star, chal_2[2*i+1][k]);
        poly_qshow_add(sum_gamma_X_star_ij_2, sum_gamma_X_star_ij_2, tmp_poly);
      }
      poly_qshow_mul_xj(tmp_poly, sum_gamma_X_star_ij_2, PARAM_N_SHOW/2); 
      poly_qshow_add(z1i_star, sum_gamma_X_star_ij_1, tmp_poly); // SR_{2i,1}* + x^{n_iss/2}.SR_{2i+1,1}*
      poly_qshow_sub(tmp_poly, sum_gamma_X_star_ij_1, tmp_poly); // SR_{2i,1}* - x^{n_iss/2}.SR_{2i+1,1}*
      poly_qshow_mul(z1i_star, z1i_star, proof_2->z1->entries[j]); // (SR_{2i,1}* + x^{n_iss/2}.SR_{2i+1,1}*)_j . [z_1]_j
      poly_qshow_add(c, c, z1i_star);
      poly_qshow_mul(tmp_poly, tmp_poly, proof_2->z1->entries[j]); // (SR_{2i,1}* - x^{n_iss/2}.SR_{2i+1,1}*)_j . [z_1]_j
      poly_qshow_conjugate(z1i_star, tmp_poly); // (SR_{2i,1} + x^{n_iss/2}.SR_{2i+1,1})_j . [z_1*]_j
      poly_qshow_add(c, c, z1i_star);
    }
    poly_qshow_mul_scalar(c, c, PARAM_TWO_INVMOD_Q_SHOW);
    poly_qshow_add(z1s_z1, z1s_z1, c);

    poly_qshow_mul(tmp_poly, z1s_z1, chal_3_l->entries[i]); // multiply by mu_i
    poly_qshow_add(sum_mu_gamma[0], sum_mu_gamma[0], tmp_poly); // add to sum_mu_gamma[0]
  }
  poly_qshow_mul(sum_mu_gamma[0], sum_mu_gamma[0], proof_2->c);
  poly_qshow_add(t0, t0, sum_mu_gamma[0]);
  
  /**********************************************
  * Constant terms in t0:
  *   -c^2.(mu_{l+dk+1} + sum_{i < l} mu_i.(d_i + f_i) + sum_{i < dk} mu_{l+i}.u_embed[i])
  **********************************************/
  poly_qshow_set(tmp_poly, chal_3_1);
  for (i = 0; i < PARAM_L_SHOW; i++) {
    poly_qshow_add(c, d_i->entries[i], proof_2->f->entries[i]); // c used as temp variable
    poly_qshow_mul(c, c, chal_3_l->entries[i]);
    poly_qshow_add(tmp_poly, tmp_poly, c);
  }
  for (i = 0; i < PARAM_D; i++) {
    poly_qshow_vec_k_mul_inner(c, chal_3_dk[i], u_embed[i]);
    poly_qshow_add(tmp_poly, tmp_poly, c);
  }
  poly_qshow_mul(c, proof_2->c, proof_2->c);
  poly_qshow_mul(tmp_poly, tmp_poly, c);
  poly_qshow_sub(t0, t0, tmp_poly);

  // computing fourth challenge
  poly_qshow_zero(c);
  buf[0] = 4;
  poly_qshow_pack(buf + CHAL3_SHOW_INPUT_BYTES, t0);
  poly_qshow_pack(buf + CHAL3_SHOW_INPUT_BYTES + POLYQSHOW_PACKEDBYTES, proof_2->t1);
  shake256(challenge_seed, SEED_BYTES, buf, CHAL4_SHOW_INPUT_BYTES);
  poly_qshow_sample_challenge(c, challenge_seed, DOMAIN_SEPARATOR_CHAL4_SHOW, proof_2->ctr_c, SEED_BYTES);

  /****** last check c=c' ******/
  is_valid = is_valid && poly_qshow_equal(proof_2->c, c);

verify_2_cleanup:
  // clean up
  // clean up polynomials
  poly_qshow_clear(tmp_poly);
  poly_qshow_clear(t0);
  poly_qshow_clear(c);
  poly_qshow_clear(z1s_z1);
  poly_qshow_clear(z1i_star);
  poly_qshow_clear(sum_gamma_X_star_ij_1);
  poly_qshow_clear(sum_gamma_X_star_ij_2);
  for (i = 0; i < 4; i++) {
    poly_qshow_clear(sum_mu_gamma[i]);
    poly_qshow_zero(sum_mu_gamma[i]);
  }
  // clean up vectors and matrices
  poly_qshow_vec_m2_d_clear(b);
  poly_qshow_vec_d_clear(tmp_vec_d);
  poly_qshow_vec_d_clear(w_H);
  poly_qshow_vec_l_clear(c_i);
  poly_qshow_vec_l_clear(c_i_prime);
  poly_qshow_vec_l_clear(d_i);
  poly_qshow_vec_256_l_clear(tmp_vec_256_l);
  poly_qshow_vec_256_l_clear(Z);
  poly_qshow_vec_k_clear(Gz1_w2);
  poly_qshow_vec_k_clear(tmp_vec_k);
  poly_qshow_vec_k_clear(sum_vec_k);
  poly_qshow_mat_k_k_clear(chal_3_quad_matrix);
  poly_qshow_mat_d_m1_clear(A1);
  poly_qshow_mat_d_m2_d_clear(A2);
  poly_qshow_mat_256l_m2_d_clear(Byg);
   // clean up challenges
  for (i = 0; i < PARAM_ARP_SHOW; i++) {
    poly_qshow_vec_m1_clear(chal_1[i]);
  }
  poly_qshow_vec_l_clear(chal_3_l);
  for (i = 0; i < PARAM_D; i++) {
    poly_qshow_vec_k_clear(chal_3_dk[i]);
  }
  poly_qshow_clear(chal_3_1);

  return is_valid;
}