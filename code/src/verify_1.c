#include "arith.h"
#include "randombytes.h"
#include "poly_qiss_sampling.h"
#include "bsig_signer.h"
#include "bsig_user.h"
#include "macros.h"
#include "fips202.h"

/*************************************************
* Name:        verify_1
*
* Description: Signer verification of the zero-knowledge proof of commitment opening
*              and verifiable encryption
*
* Arguments:   - const proof_1_t *proof_1: pointer to issuance proof structure
*              - const poly_qiss_mat_k_k *A_embed: array of polynomial matrices, subring embedding of q1.A'
*              - const poly_qiss_mat_k_k *B_embed: array of polynomial matrices to host subring embedding of q1.[tG - B|A3]
*              - const poly_qiss_mat_k_k *D_embed: array of polynomial matrices to host subring embedding of q1.d
*              - const poly_qiss_mat_k_k *Ae_be_embed: array of polynomial matrices to host subring embedding of [A_e | b_e]^T
*              - const poly_qiss_vec_k *cmt_embed: array of polynomial vectors to host subring embedding of q1.cmt
*              - const poly_qiss_vec_k *ct_embed: array of polynomial vectors to host subring embedding of ct
*              - const uint8_t *crs_seed: pointer to byte array containing the CRS seed (allocated SEED_BYTES bytes)
*              - const uint8_t *seed: pointer to byte array containing the seed for public parameters (allocated SEED_BYTES bytes)
* 
* Returns 1 if proof_1 could be verified correctly and 0 otherwise
**************************************************/
int verify_1(
    const proof_1_t         *proof_1, 
    const poly_qiss_mat_k_k A_embed[PARAM_D][PARAM_D], 
    const poly_qiss_mat_k_k B_embed[PARAM_D][PARAM_K*(PARAM_D+1)], 
    const poly_qiss_mat_k_k D_embed[PARAM_D], 
    const poly_qiss_mat_k_k Ae_be_embed[PARAM_DE+1][PARAM_ME],
    const poly_qiss_vec_k   cmt_embed[PARAM_D], 
    const poly_qiss_vec_k   ct_embed[PARAM_DE+1],  
    const uint8_t           crs_seed[CRS_SEED_BYTES], 
    const uint8_t           seed[SEED_BYTES]) {
  size_t i,j,k;
  int is_valid = 1;
  uint8_t buf[CHAL4_ISS_INPUT_BYTES] = {0}, challenge_seed[SEED_BYTES];
  uint128 sq_norm_z1, sq_norm_z2;
  uint64_t sq_norm_z3;
  coeff_qiss tmp_coeff;
  poly_qiss tmp_poly, t0, sum_mu_gamma[4], c, z1s_z1, z1i_star, sum_gamma_X_star_ij_1, sum_gamma_X_star_ij_2;
  poly_qiss_vec_m2_d b;
  poly_qiss_vec_d tmp_vec_d, w_H;
  poly_qiss_vec_l c_i, c_i_prime, d_i;
  poly_qiss_vec_256_l tmp_vec_256_l, Z;
  poly_qiss_mat_d_m1 A1;
  poly_qiss_mat_d_m2_d A2;
  poly_qiss_mat_256l_m2_d Byg;
  // challenges
  poly_qiss_vec_m1 chal_11[PARAM_ARP_ISS];
  poly_qiss_vec_k chal_12[PARAM_ARP_ISS][PARAM_DE+1], enc_arp_z1[PARAM_DE+1];
  coeff_qiss chal_2[2*PARAM_L_ISS][PARAM_ARP_ISS + 4 + PARAM_N_ISS - 1];
  poly_qiss_vec_l chal_3_l;
  poly_qiss_vec_k chal_3_dk[PARAM_D];
  poly_qiss chal_3_1;

  // init
  // init polynomials
  poly_qiss_init(tmp_poly);
  poly_qiss_init(t0);
  poly_qiss_init(c);
  poly_qiss_init(z1s_z1);
  poly_qiss_init(z1i_star);
  poly_qiss_init(sum_gamma_X_star_ij_1);
  poly_qiss_init(sum_gamma_X_star_ij_2);
  for (i = 0; i < 4; i++) {
    poly_qiss_init(sum_mu_gamma[i]);
    poly_qiss_zero(sum_mu_gamma[i]);
  }
  // init vectors and matrices
  poly_qiss_vec_m2_d_init(b);
  poly_qiss_vec_d_init(tmp_vec_d);
  poly_qiss_vec_d_init(w_H);
  poly_qiss_vec_l_init(c_i);
  poly_qiss_vec_l_init(c_i_prime);
  poly_qiss_vec_l_init(d_i);
  poly_qiss_vec_256_l_init(tmp_vec_256_l);
  poly_qiss_vec_256_l_init(Z);
  for (i = 0; i < PARAM_DE+1; i++) {
    poly_qiss_vec_k_init(enc_arp_z1[i]);
  }
  poly_qiss_mat_d_m1_init(A1);
  poly_qiss_mat_d_m2_d_init(A2);
  poly_qiss_mat_256l_m2_d_init(Byg);
   // init challenges
  for (i = 0; i < PARAM_ARP_ISS; i++) {
    poly_qiss_vec_m1_init(chal_11[i]);
    for (j = 0; j < PARAM_DE+1; j++) {
      poly_qiss_vec_k_init(chal_12[i][j]);
    }
  }
  poly_qiss_vec_l_init(chal_3_l);
  for (i = 0; i < PARAM_D; i++) {
    poly_qiss_vec_k_init(chal_3_dk[i]);
  }
  poly_qiss_init(chal_3_1);

  // expanding CRS
  poly_qiss_mat_d_m1_uniform(A1, crs_seed, DOMAIN_SEPARATOR_A1_ISS);
  poly_qiss_mat_d_m2_d_uniform(A2, crs_seed, DOMAIN_SEPARATOR_A2_ISS);
  poly_qiss_mat_256l_m2_d_uniform(Byg, crs_seed, DOMAIN_SEPARATOR_BYG_ISS);
  poly_qiss_vec_m2_d_uniform(b, crs_seed, DOMAIN_SEPARATOR_B_ISS);

  // byte-packing the statement
  // the matrices A_embed, B_embed, D_embed, Ae_be_embed, are derived from seed
  // the statement is thus (crs_seed || seed || byte-packing(cmt_embed) || byte-packing(ct_embed)).
  for (i = 0; i < CRS_SEED_BYTES; i++) {
    buf[1 + i] = crs_seed[i];
  }
  for (i = 0; i < SEED_BYTES; i++) {
    buf[1 + CRS_SEED_BYTES + i] = seed[i];
  }
  for (i = 0; i < PARAM_D; i++) {
    poly_qiss_vec_k_pack(buf + 1 + CRS_SEED_BYTES + SEED_BYTES + i*POLYQISS_VECK_PACKEDBYTES, cmt_embed[i]);
  }
  for (i = 0; i < PARAM_DE+1; i++) {
    poly_qiss_vec_k_pack(buf + 1 + CRS_SEED_BYTES + SEED_BYTES + PARAM_D*POLYQISS_VECK_PACKEDBYTES + i*POLYQISS_VECK_PACKEDBYTES, ct_embed[i]);
  }

  ////////////////////////////////////// NO INITIALIZATIONS BELOW THIS LINE

  /****** reconstructing first message ******/
  // recomputing tmp = A1.z1 + A2.z2_1 - c.2^D.tA_H 
  poly_qiss_mat_d_m1_mul_vec_m1(tmp_vec_d, A1, proof_1->z1);
  poly_qiss_mat_d_m2_d_mul_vec_m2_d(w_H, A2, proof_1->z2_1);  // using w_H as temp variable
  poly_qiss_vec_d_add(tmp_vec_d, tmp_vec_d, w_H);
  poly_qiss_vec_d_mul_poly_qiss(w_H, proof_1->tA_H, proof_1->c); // using w_H as temp variable
  poly_qiss_vec_d_mul_scalar(w_H, w_H, (coeff_qiss)(1 << PARAM_D_ROUND_ISS));
  poly_qiss_vec_d_sub(tmp_vec_d, tmp_vec_d, w_H);
  // using hint
  poly_qiss_vec_d_useGhint(w_H, proof_1->hint, tmp_vec_d); 

  // computing first challenge with w_H so that we can use w_H as temp variable afterwards
  // recomputing first challenge
  buf[0] = 1;
  poly_qiss_vec_d_pack(buf + ISS_CHALLENGE_BASE_BYTES, proof_1->tA_H);
  poly_qiss_vec_d_pack(buf + ISS_CHALLENGE_BASE_BYTES + POLYQISS_VECD_PACKEDBYTES, w_H);
  poly_qiss_vec_256_l_pack(buf + ISS_CHALLENGE_BASE_BYTES + 2*POLYQISS_VECD_PACKEDBYTES, proof_1->tB);
  shake256(challenge_seed, SEED_BYTES, buf, CHAL1_ISS_INPUT_BYTES);
  for (i = 0; i < PARAM_ARP_ISS; i++) {
    poly_qiss_vec_m1_vec_k_binomial(chal_11[i], chal_12[i], challenge_seed, DOMAIN_SEPARATOR_CHAL1_ISS, i, SEED_BYTES);
  }

  /****** performing all checks except c=c' for early reject ******/
  poly_qiss_vec_d_mul_scalar(w_H, w_H, PARAM_GAMMA_ISS);
  poly_qiss_vec_d_sub(w_H, w_H, tmp_vec_d); // w_H contains z2_2
  sq_norm_z1 = poly_qiss_vec_m1_norm2(proof_1->z1);
  sq_norm_z2 = poly_qiss_vec_m2_d_norm2(proof_1->z2_1) + poly_qiss_vec_d_norm2(w_H);
  sq_norm_z3 = 0;
  for (i = 0; i < PARAM_ARP_ISS; i++) {
    tmp_coeff = proof_1->z3[i];
    CHK_UI_OVF_ADDITION(sq_norm_z3, tmp_coeff * tmp_coeff);
  }
  tmp_coeff = poly_qiss_vec_d_norm_inf(proof_1->hint);
  is_valid = is_valid && (sq_norm_z1 <= PARAM_B1SQ_ISS) && (sq_norm_z2 <= ((uint128)PARAM_B2SQ_ISS_LOW64 + (((uint128)PARAM_B2SQ_ISS_HIGH64) << 64))) && (sq_norm_z3 <= PARAM_B3SQ_ISS);
  is_valid = is_valid && (tmp_coeff <= PARAM_Q_GAMMA_ISS/2);
  for (i = 0; i < PARAM_L_ISS; i++) {
    is_valid = is_valid && (poly_qiss_get_coeff(proof_1->f->entries[i], 0) == 0) && (poly_qiss_get_coeff(proof_1->f->entries[i], PARAM_N_ISS/2) == 0);
  }
  if (!is_valid) {
    goto verify_1_cleanup; 
  }

  /****** reconstructing second message ******/
  buf[0] = 2;
  for (i = 0; i < PARAM_ARP_ISS; i++) {
    coeff_qiss_pack(buf + CHAL1_ISS_INPUT_BYTES + i*COEFFQISS_PACKEDBYTES, proof_1->z3[i]);
  }
  // recomputing second challenge
  shake256(challenge_seed, SEED_BYTES, buf, CHAL2_ISS_INPUT_BYTES);
  for (i = 0; i < 2*PARAM_L_ISS; i++) {
    vec_qiss_uniform(chal_2[i], challenge_seed, DOMAIN_SEPARATOR_CHAL2_ISS, i, SEED_BYTES); // writes PARAM_ARP_ISS + 4 + PARAM_N_ISS - 1 uniform numbers to the first argument
  }

  /****** reconstructing third message ******/
  buf[0] = 3;
  poly_qiss_vec_l_pack(buf + CHAL2_ISS_INPUT_BYTES, proof_1->f);
  // recomputing third challenge
  shake256(challenge_seed, SEED_BYTES, buf, CHAL3_ISS_INPUT_BYTES);
  poly_qiss_vec_l_1_uniform(chal_3_l, chal_3_1, challenge_seed, DOMAIN_SEPARATOR_CHAL3_ISS, SEED_BYTES);
  for (i = 0; i < PARAM_D; i++) {
    poly_qiss_vec_k_uniform(chal_3_dk[i], challenge_seed, DOMAIN_SEPARATOR_CHAL3_ISS, i+PARAM_L_ISS+1, SEED_BYTES);
  }
  
  /****** reconstructing fourth message ******/
  /**********************************************
  * Recomputing commitment t0
  **********************************************/
  // computing <b, z2> - c.t1
  poly_qiss_vec_m2_d_mul_inner(t0, b, proof_1->z2_1);
  poly_qiss_mul(tmp_poly, proof_1->c, proof_1->t1);
  poly_qiss_sub(t0, t0, tmp_poly);

  // Computing sum_i mu_i (gamma_{2i,256+j} + x^(n_iss/2).gamma_{2i+1,256+j}), c_i, c_i' and d_i
  poly_qiss_vec_l_zero(d_i);
  for (i = 0; i < PARAM_L_ISS; i++) {
    poly_qiss_mul_xj(c, chal_3_l->entries[i], PARAM_N_ISS/2); // using c as temp variable
    for (j = 0; j < 4; j++) {
      poly_qiss_mul_scalar(tmp_poly, chal_3_l->entries[i], chal_2[2*i][PARAM_ARP_ISS + j]); // mu_i gamma_{2i,256+j}
      poly_qiss_add(sum_mu_gamma[j], sum_mu_gamma[j], tmp_poly);
      poly_qiss_mul_scalar(tmp_poly, c, chal_2[2*i+1][PARAM_ARP_ISS + j]); // mu_i.x^(n_iss/2).gamma_{2i+1,256+j}
      poly_qiss_add(sum_mu_gamma[j], sum_mu_gamma[j], tmp_poly);
    }
    // c_i and c_i'
    poly_qiss_set_coeff(c_i->entries[i], 0, -chal_2[2*i+1][PARAM_ARP_ISS+3+PARAM_N_ISS/2]);
    poly_qiss_set_coeff(c_i_prime->entries[i], 0, chal_2[2*i+1][PARAM_ARP_ISS+3+PARAM_N_ISS/2]);
    for (j = 1; j < PARAM_N_ISS/2; j++) {
      poly_qiss_set_coeff(c_i->entries[i], j, chal_2[2*i][PARAM_ARP_ISS+3+j] - chal_2[2*i+1][PARAM_ARP_ISS+3+j+PARAM_N_ISS/2]);
      poly_qiss_set_coeff(c_i_prime->entries[i], j, chal_2[2*i+1][PARAM_ARP_ISS+3+PARAM_N_ISS/2-j] - chal_2[2*i][PARAM_ARP_ISS+3+PARAM_N_ISS-j]);
    }
    poly_qiss_set_coeff(c_i->entries[i], PARAM_N_ISS/2, chal_2[2*i][PARAM_ARP_ISS+3+PARAM_N_ISS/2]);
    poly_qiss_set_coeff(c_i_prime->entries[i], PARAM_N_ISS/2, -chal_2[2*i][PARAM_ARP_ISS+3+PARAM_N_ISS/2]);
    for (j = PARAM_N_ISS/2 + 1; j < PARAM_N_ISS; j++) {
      poly_qiss_set_coeff(c_i->entries[i], j, chal_2[2*i][PARAM_ARP_ISS+3+j] + chal_2[2*i+1][PARAM_ARP_ISS+3+j-PARAM_N_ISS/2]);
      poly_qiss_set_coeff(c_i_prime->entries[i], j, -(chal_2[2*i][PARAM_ARP_ISS+3+PARAM_N_ISS-j] + chal_2[2*i+1][PARAM_ARP_ISS+3+3*PARAM_N_ISS/2-j]));
    }
    poly_qiss_mul_scalar(c_i->entries[i], c_i->entries[i], PARAM_TWO_INVMOD_Q_ISS);
    poly_qiss_mul_scalar(c_i_prime->entries[i], c_i_prime->entries[i], PARAM_TWO_INVMOD_Q_ISS);
    // d_i
    poly_qiss_zero(tmp_poly);
    poly_qiss_muladd_constant(tmp_poly, chal_2[2*i+1][PARAM_ARP_ISS    ], PARAM_B_R1_SQ);
    poly_qiss_muladd_constant(tmp_poly, chal_2[2*i+1][PARAM_ARP_ISS + 1], PARAM_B_R2_SQ);
    poly_qiss_muladd_constant(tmp_poly, chal_2[2*i+1][PARAM_ARP_ISS + 2], PARAM_B_RE_SQ);
    for (j = 0; j < PARAM_ARP_ISS; j++) {
      poly_qiss_muladd_constant(tmp_poly, chal_2[2*i+1][j], proof_1->z3[j]);
    }
    poly_qiss_mul_xj(d_i->entries[i], tmp_poly, PARAM_N_ISS/2);
    poly_qiss_muladd_constant(d_i->entries[i], chal_2[2*i][PARAM_ARP_ISS    ], PARAM_B_R1_SQ);
    poly_qiss_muladd_constant(d_i->entries[i], chal_2[2*i][PARAM_ARP_ISS + 1], PARAM_B_R2_SQ);
    poly_qiss_muladd_constant(d_i->entries[i], chal_2[2*i][PARAM_ARP_ISS + 2], PARAM_B_RE_SQ);
    for (j = 0; j < PARAM_ARP_ISS; j++) {
      poly_qiss_muladd_constant(d_i->entries[i], chal_2[2*i][j], proof_1->z3[j]);
    }
  }
  poly_qiss_mul_scalar(sum_mu_gamma[3], sum_mu_gamma[3], PARAM_TWO_INVMOD_Q_ISS);

  // p^{-1}.(z_{1,b}.ct_embed - Ae_be_embed.z_{1,e} - [0 | round(p/2).z_{1,h}]^T) // z1i_star used as temp variable
  for (i = 0; i < (PARAM_DE+1)*PARAM_K_ISS; i++) {
    if (i < PARAM_DE*PARAM_K_ISS) {
      poly_qiss_zero(tmp_poly);
    } else {
      poly_qiss_mul_scalar(tmp_poly, proof_1->z1->entries[IDX_MSG_ISS + (i % PARAM_K_ISS)], PARAM_HALF_P); // round(p/2).z_{1,h}
    }
    for (j = 0; j < PARAM_ME*PARAM_K_ISS; j++) {
      poly_qiss_mul(z1i_star, Ae_be_embed[i / PARAM_K_ISS][j / PARAM_K_ISS]->rows[i % PARAM_K_ISS]->entries[j % PARAM_K_ISS], proof_1->z1->entries[IDX_RE_ISS + j]); 
      poly_qiss_add(tmp_poly, tmp_poly, z1i_star);
    }
    // - z_{1,b}.ct'
    poly_qiss_mul(z1i_star, ct_embed[i / PARAM_K_ISS]->entries[i % PARAM_K_ISS], proof_1->z1->entries[IDX_B_ISS]);
    poly_qiss_sub(tmp_poly, tmp_poly, z1i_star);
    // .-p^{-1}
    poly_qiss_mul_scalar(tmp_poly, tmp_poly, PARAM_P_INVMOD_Q_ISS_NEG);

    poly_qiss_set(enc_arp_z1[i / PARAM_K_ISS]->entries[i % PARAM_K_ISS], tmp_poly);
  }

  /**********************************************
  * Quadratic terms in t0:
  *  sum_mu_gamma[0].<z_{1,1}'*, z_{1,1}'> + sum_mu_gamma[1].<z_{1,23}'*, z_{1,23}'> 
  *    + sum_mu_gamma[2].<z_{1,e}'*, z_{1,e}'> + sum_mu_gamma[3].(<z_{1,h}*, z_{1,h} - z_{1,b}.one> + <z_{1,h}*, z_{1,h} - z_{1,b}.one>*)
  *    + mu_{l+dk+1}.z_{1,b}^2
  **********************************************/
  // <z_{1,1}'*,z_{1,1}'>
  poly_qiss_zero(z1s_z1);
  for (i = IDX_R1_ISS; i < IDX_R23_ISS; i++) {
    poly_qiss_conjugate(z1i_star, proof_1->z1->entries[i]);
    poly_qiss_mul(tmp_poly, z1i_star, proof_1->z1->entries[i]);
    poly_qiss_add(z1s_z1, z1s_z1, tmp_poly);
  }
  poly_qiss_mul(tmp_poly, z1s_z1, sum_mu_gamma[0]);
  poly_qiss_add(t0, t0, tmp_poly);
  // <z_{1,23}'*,z_{1,23}'>
  poly_qiss_zero(z1s_z1);
  for (i = IDX_R23_ISS; i < IDX_RE_ISS; i++) {
    poly_qiss_conjugate(z1i_star, proof_1->z1->entries[i]);
    poly_qiss_mul(tmp_poly, z1i_star, proof_1->z1->entries[i]);
    poly_qiss_add(z1s_z1, z1s_z1, tmp_poly);
  }
  poly_qiss_mul(tmp_poly, z1s_z1, sum_mu_gamma[1]);
  poly_qiss_add(t0, t0, tmp_poly);
  // <z_{1,e}'*,z_{1,e}'>
  poly_qiss_zero(z1s_z1);
  for (i = IDX_RE_ISS; i < IDX_MSG_ISS; i++) {
    poly_qiss_conjugate(z1i_star, proof_1->z1->entries[i]);
    poly_qiss_mul(tmp_poly, z1i_star, proof_1->z1->entries[i]);
    poly_qiss_add(z1s_z1, z1s_z1, tmp_poly);
  }
  poly_qiss_mul(tmp_poly, z1s_z1, sum_mu_gamma[2]);
  poly_qiss_add(t0, t0, tmp_poly);
  // one term
  for (i = 0; i < PARAM_N_ISS; i++) {
    poly_qiss_set_coeff(c, i, 1); // poly with all one // using c as temp variable
  }
  poly_qiss_zero(z1s_z1);
  poly_qiss_mul(c, proof_1->z1->entries[IDX_B_ISS], c); // using c as temp variable
  for (i = IDX_MSG_ISS; i < IDX_B_ISS; i++) {
    poly_qiss_conjugate(z1i_star, proof_1->z1->entries[i]);
    poly_qiss_sub(tmp_poly, proof_1->z1->entries[i], c);
    poly_qiss_mul(tmp_poly, z1i_star, tmp_poly);
    poly_qiss_add(z1s_z1, z1s_z1, tmp_poly);
  }
  poly_qiss_conjugate(tmp_poly, z1s_z1);
  poly_qiss_add(z1s_z1, z1s_z1, tmp_poly);
  poly_qiss_mul(tmp_poly, z1s_z1, sum_mu_gamma[3]);
  poly_qiss_add(t0, t0, tmp_poly);
  // mu_{l+dk+1}.z_{1,b}^2
  poly_qiss_mul(tmp_poly, chal_3_1, proof_1->z1->entries[IDX_B_ISS]);
  poly_qiss_mul(tmp_poly, tmp_poly, proof_1->z1->entries[IDX_B_ISS]);
  poly_qiss_add(t0, t0, tmp_poly);

  /**********************************************
  * Linear terms in t0:
  *   Z = c.tB - Byg.z2
  * 
  *   c.sum_{i < l} mu_i.(Z_i + 2^{-1}.(SE_{2i} + x^{n_iss/2}.SE_{2i+1}).Z*_{:256/n_iss}
  *                    + 2^{-1}.(SE_{2i}* + x^{n_iss/2}.SE_{2i+1}*).Z_{:256/n_iss}
  *                    + 2^{-1}.(SR_{2i,1} + x^{n_iss/2}.SR_{2i+1,1}).z1* + 2^{-1}.(SR_{2i,1}* + x^{n_iss/2}.SR_{2i+1,1}*).z1
  *                    + 2^{-1}.(SR_{2i,2} + x^{n_iss/2}.SR_{2i+1,2}).enc_arp_z1* + 2^{-1}.(SR_{2i,2}* + x^{n_iss/2}.SR_{2i+1,2}*).enc_arp_z1
  *                    + c_i[i].z_{1,b}* + c_i_prime[i].z_{1,b})
  *   + c.sum_{i < dk} mu_{l+i}.[A_theta.z_{1,1} + B_theta.z_{1,23} + D.z_{1,h} - z_{1,b}.cmt_embed]_i
  **********************************************/
  // computing Z = c.tB - Byg.z2
  poly_qiss_mat_256l_m2_d_mul_vec_m2_d(tmp_vec_256_l, Byg, proof_1->z2_1);
  poly_qiss_vec_256_l_mul_poly_qiss(Z, proof_1->tB, proof_1->c);
  poly_qiss_vec_256_l_sub(Z, Z, tmp_vec_256_l);

  poly_qiss_zero(sum_mu_gamma[0]); // using sum_mu_gamma[0] as temp variable to store sum before multiplying by c
  for (i = 0; i < PARAM_L_ISS; i++) {
    // using z1s_z1 to store sum before multiplying by mu_i
    poly_qiss_mul(z1s_z1, c_i_prime->entries[i], proof_1->z1->entries[IDX_B_ISS]);
    poly_qiss_conjugate(tmp_poly, proof_1->z1->entries[IDX_B_ISS]);
    poly_qiss_mul(tmp_poly, tmp_poly, c_i->entries[i]);
    poly_qiss_add(z1s_z1, z1s_z1, tmp_poly);

    poly_qiss_add(z1s_z1, z1s_z1, Z->entries[PARAM_ARP_DIV_N_ISS + i]);

    poly_qiss_zero(c); // c used as temp variable
    // term in SE_i
    for (j = 0; j < PARAM_ARP_DIV_N_ISS; j++) {
      // [SE_{2i}*]_j
      poly_qiss_set_coeff(sum_gamma_X_star_ij_1, 0, chal_2[2*i][j * PARAM_N_ISS]);
      for (k = 1; k < PARAM_N_ISS; k++) {
        poly_qiss_set_coeff(sum_gamma_X_star_ij_1, k, - chal_2[2*i][(j + 1) * PARAM_N_ISS - k]); // set conjugate directly
      }
      // [SE_{2i+1}*]_j
      poly_qiss_set_coeff(sum_gamma_X_star_ij_2, 0, chal_2[2*i+1][j * PARAM_N_ISS]);
      for (k = 1; k < PARAM_N_ISS; k++) {
        poly_qiss_set_coeff(sum_gamma_X_star_ij_2, k, - chal_2[2*i+1][(j + 1) * PARAM_N_ISS - k]); // set conjugate directly
      }
      poly_qiss_mul_xj(tmp_poly, sum_gamma_X_star_ij_2, PARAM_N_ISS/2); 
      poly_qiss_add(z1i_star, sum_gamma_X_star_ij_1, tmp_poly); // SE_{2i}* + x^{n_iss/2}.SE_{2i+1}*
      poly_qiss_sub(tmp_poly, sum_gamma_X_star_ij_1, tmp_poly); // SE_{2i}* - x^{n_iss/2}.SE_{2i+1}*
      poly_qiss_mul(z1i_star, z1i_star, Z->entries[j]); // (SE_{2i}* + x^{n_iss/2}.SE_{2i+1}*)_j . Z_j
      poly_qiss_add(c, c, z1i_star);
      poly_qiss_mul(tmp_poly, tmp_poly, Z->entries[j]); // (SE_{2i}* - x^{n_iss/2}.SE_{2i+1}*)_j . Z_j
      poly_qiss_conjugate(z1i_star, tmp_poly); // (SE_{2i} + x^{n_iss/2}.SE_{2i+1})_j . Z_j*
      poly_qiss_add(c, c, z1i_star);
    }
    // term in SR_{i,1}
    for (j = 0; j < PARAM_M1_ISS; j++) {
      poly_qiss_zero(sum_gamma_X_star_ij_1);
      poly_qiss_zero(sum_gamma_X_star_ij_2);
      for (k = 0; k < PARAM_ARP_ISS; k++) {
        poly_qiss_conjugate(z1i_star, chal_11[k]->entries[j]); // z1i_star used as temp variable
        // [SR_{2i,1}*]_j
        poly_qiss_mul_scalar(tmp_poly, z1i_star, chal_2[2*i][k]);
        poly_qiss_add(sum_gamma_X_star_ij_1, sum_gamma_X_star_ij_1, tmp_poly);
        // [SR_{2i+1,1}*]_j
        poly_qiss_mul_scalar(tmp_poly, z1i_star, chal_2[2*i+1][k]);
        poly_qiss_add(sum_gamma_X_star_ij_2, sum_gamma_X_star_ij_2, tmp_poly);
      }
      poly_qiss_mul_xj(tmp_poly, sum_gamma_X_star_ij_2, PARAM_N_ISS/2); 
      poly_qiss_add(z1i_star, sum_gamma_X_star_ij_1, tmp_poly); // SR_{2i,1}* + x^{n_iss/2}.SR_{2i+1,1}*
      poly_qiss_sub(tmp_poly, sum_gamma_X_star_ij_1, tmp_poly); // SR_{2i,1}* - x^{n_iss/2}.SR_{2i+1,1}*
      poly_qiss_mul(z1i_star, z1i_star, proof_1->z1->entries[j]); // (SR_{2i,1}* + x^{n_iss/2}.SR_{2i+1,1}*)_j . [z_1]_j
      poly_qiss_add(c, c, z1i_star);
      poly_qiss_mul(tmp_poly, tmp_poly, proof_1->z1->entries[j]); // (SR_{2i,1}* - x^{n_iss/2}.SR_{2i+1,1}*)_j . [z_1]_j
      poly_qiss_conjugate(z1i_star, tmp_poly); // (SR_{2i,1} + x^{n_iss/2}.SR_{2i+1,1})_j . [z_1*]_j
      poly_qiss_add(c, c, z1i_star);
    }
    // term in SR_{i,2}
    for (j = 0; j < (PARAM_DE+1)*PARAM_K_ISS; j++) {
      poly_qiss_zero(sum_gamma_X_star_ij_1);
      poly_qiss_zero(sum_gamma_X_star_ij_2);
      for (k = 0; k < PARAM_ARP_ISS; k++) {
        poly_qiss_conjugate(z1i_star, chal_12[k][j / PARAM_K_ISS]->entries[j % PARAM_K_ISS]); // z1i_star used as temp variable
        // [SR_{2i,1}*]_j
        poly_qiss_mul_scalar(tmp_poly, z1i_star, chal_2[2*i][k]);
        poly_qiss_add(sum_gamma_X_star_ij_1, sum_gamma_X_star_ij_1, tmp_poly);
        // [SR_{2i+1,1}*]_j
        poly_qiss_mul_scalar(tmp_poly, z1i_star, chal_2[2*i+1][k]);
        poly_qiss_add(sum_gamma_X_star_ij_2, sum_gamma_X_star_ij_2, tmp_poly);
      }
      poly_qiss_mul_xj(tmp_poly, sum_gamma_X_star_ij_2, PARAM_N_ISS/2); 
      poly_qiss_add(z1i_star, sum_gamma_X_star_ij_1, tmp_poly); // SR_{2i,2}* + x^{n_iss/2}.SR_{2i+1,2}*
      poly_qiss_sub(tmp_poly, sum_gamma_X_star_ij_1, tmp_poly); // SR_{2i,2}* - x^{n_iss/2}.SR_{2i+1,2}*
      poly_qiss_mul(z1i_star, z1i_star, enc_arp_z1[j / PARAM_K_ISS]->entries[j % PARAM_K_ISS]); // (SR_{2i,2}* + x^{n_iss/2}.SR_{2i+1,2}*)_j . [enc_arp_z1]_j
      poly_qiss_add(c, c, z1i_star);
      poly_qiss_mul(tmp_poly, tmp_poly, enc_arp_z1[j / PARAM_K_ISS]->entries[j % PARAM_K_ISS]); // (SR_{2i,2}* - x^{n_iss/2}.SR_{2i+1,2}*)_j . [enc_arp_z1]_j
      poly_qiss_conjugate(z1i_star, tmp_poly); // (SR_{2i,2} + x^{n_iss/2}.SR_{2i+1,2})_j . [enc_arp_z1*]_j
      poly_qiss_add(c, c, z1i_star);
    }
    poly_qiss_mul_scalar(c, c, PARAM_TWO_INVMOD_Q_ISS);
    poly_qiss_add(z1s_z1, z1s_z1, c);

    poly_qiss_mul(tmp_poly, z1s_z1, chal_3_l->entries[i]); // multiply by mu_i
    poly_qiss_add(sum_mu_gamma[0], sum_mu_gamma[0], tmp_poly); // add to sum_mu_gamma[0]
  }
  for (i = 0; i < PARAM_D*PARAM_K_ISS; i++) {
    // c temp variable
    // A_theta.z_{1,1}
    poly_qiss_mul_scalar(tmp_poly, proof_1->z1->entries[i], PARAM_Q1_ISS);
    for (j = 0; j < PARAM_D*PARAM_K_ISS; j++) {
      poly_qiss_mul(c, A_embed[i / PARAM_K_ISS][j / PARAM_K_ISS]->rows[i % PARAM_K_ISS]->entries[j % PARAM_K_ISS], proof_1->z1->entries[IDX_R12_ISS + j]); 
      poly_qiss_add(tmp_poly, tmp_poly, c);
    }
    // + B_theta.z_{1,23}
    for (j = 0; j < PARAM_K*(PARAM_D+1)*PARAM_K_ISS; j++) {
      poly_qiss_mul(c, B_embed[i / PARAM_K_ISS][j / PARAM_K_ISS]->rows[i % PARAM_K_ISS]->entries[j % PARAM_K_ISS], proof_1->z1->entries[IDX_R23_ISS + j]);
      poly_qiss_add(tmp_poly, tmp_poly, c);
    }
    // + D_theta.z_{1,h}
    for (j = 0; j < PARAM_K_ISS; j++) {
      poly_qiss_mul(c, D_embed[i / PARAM_K_ISS]->rows[i % PARAM_K_ISS]->entries[j], proof_1->z1->entries[IDX_MSG_ISS + j]);
      poly_qiss_add(tmp_poly, tmp_poly, c);
    }
    // - z_{1,b}.u
    poly_qiss_mul(c, cmt_embed[i / PARAM_K_ISS]->entries[i % PARAM_K_ISS], proof_1->z1->entries[IDX_B_ISS]);
    poly_qiss_sub(tmp_poly, tmp_poly, c);
    // multiply by mu_{l+i}
    poly_qiss_mul(tmp_poly, tmp_poly, chal_3_dk[i / PARAM_K_ISS]->entries[i % PARAM_K_ISS]);
    poly_qiss_add(sum_mu_gamma[0], sum_mu_gamma[0], tmp_poly);
  }
  poly_qiss_mul(sum_mu_gamma[0], sum_mu_gamma[0], proof_1->c);
  poly_qiss_add(t0, t0, sum_mu_gamma[0]);

  /**********************************************
  * Constant terms in t0:
  *   -c^2.(mu_{l+dk+1} + sum_{i < l} mu_i.(d_i + f_i))
  **********************************************/
  poly_qiss_set(tmp_poly, chal_3_1);
  for (i = 0; i < PARAM_L_ISS; i++) {
    poly_qiss_add(c, d_i->entries[i], proof_1->f->entries[i]); // c used as temp variable
    poly_qiss_mul(c, c, chal_3_l->entries[i]);
    poly_qiss_add(tmp_poly, tmp_poly, c);
  }
  poly_qiss_mul(c, proof_1->c, proof_1->c);
  poly_qiss_mul(tmp_poly, tmp_poly, c);
  poly_qiss_sub(t0, t0, tmp_poly);

  // computing fourth challenge
  poly_qiss_zero(c); // resetting c
  buf[0] = 4;
  poly_qiss_pack(buf + CHAL3_ISS_INPUT_BYTES, t0);
  poly_qiss_pack(buf + CHAL3_ISS_INPUT_BYTES + POLYQISS_PACKEDBYTES, proof_1->t1);
  shake256(challenge_seed, SEED_BYTES, buf, CHAL4_ISS_INPUT_BYTES);
  poly_qiss_sample_challenge(c, challenge_seed, DOMAIN_SEPARATOR_CHAL4_ISS, proof_1->ctr_c, SEED_BYTES);

  /****** last check c=c' ******/
  is_valid = is_valid && poly_qiss_equal(proof_1->c, c);

verify_1_cleanup:
  // clean up
  // clean up polynomials
  poly_qiss_clear(tmp_poly);
  poly_qiss_clear(t0);
  poly_qiss_clear(c);
  poly_qiss_clear(z1s_z1);
  poly_qiss_clear(z1i_star);
  poly_qiss_clear(sum_gamma_X_star_ij_1);
  poly_qiss_clear(sum_gamma_X_star_ij_2);
  for (i = 0; i < 4; i++) {
    poly_qiss_clear(sum_mu_gamma[i]);
    poly_qiss_zero(sum_mu_gamma[i]);
  }
  // clean up vectors and matrices
  poly_qiss_vec_m2_d_clear(b);
  poly_qiss_vec_d_clear(tmp_vec_d);
  poly_qiss_vec_d_clear(w_H);
  poly_qiss_vec_l_clear(c_i);
  poly_qiss_vec_l_clear(c_i_prime);
  poly_qiss_vec_l_clear(d_i);
  poly_qiss_vec_256_l_clear(tmp_vec_256_l);
  poly_qiss_vec_256_l_clear(Z);
  for (i = 0; i < PARAM_DE+1; i++) {
    poly_qiss_vec_k_clear(enc_arp_z1[i]);
  }
  poly_qiss_mat_d_m1_clear(A1);
  poly_qiss_mat_d_m2_d_clear(A2);
  poly_qiss_mat_256l_m2_d_clear(Byg);
   // clean up challenges
  for (i = 0; i < PARAM_ARP_ISS; i++) {
    poly_qiss_vec_m1_clear(chal_11[i]);
    for (j = 0; j < PARAM_DE+1; j++) {
      poly_qiss_vec_k_clear(chal_12[i][j]);
    }
  }
  poly_qiss_vec_l_clear(chal_3_l);
  for (i = 0; i < PARAM_D; i++) {
    poly_qiss_vec_k_clear(chal_3_dk[i]);
  }
  poly_qiss_clear(chal_3_1);

  return is_valid;
}