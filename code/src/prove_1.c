#include "arith.h"
#include "randombytes.h"
#include "poly_qiss_sampling.h"
#include "bsig_signer.h"
#include "bsig_user.h"
#include "macros.h"
#include "fips202.h"

#include <math.h>

/*************************************************
* Name:        _reject_exp [static]
*
* Description: Rejection sampling given a probability threshold. Accepts
*              with probability threshold
*
* Arguments:   - long double threshold: acceptance probability threshold
* 
* Returns 1 if the sample should be rejected (ie if udbl > threshold), 0 otherwise
**************************************************/
static int _reject_exp(long double threshold) {
  uint64_t u = 0;
  long double udbl = 0;
  randombytes((uint8_t*)&u, 6); // 48 random bits
  udbl = (long double)u / (long double)(1ul << 48);
  if (udbl > threshold){
    return 1;
  }
  return 0;
}

/*************************************************
* Name:        prove_1_round1 [static]
*
* Description: Compute round 1 of zero-knowledge proof of commitment opening
*              and verifiable encryption
*
* Arguments:   - proof_1_t *proof_1: pointer to issuance proof structure
*              - poly_qiss_vec_m1 *chal_11: array of polynomial vectors to host first challenge (R = R0 - R1) (initialized)
*              - poly_qiss_vec_k *chal_12: array of polynomial vectors to host first challenge (R = R0 - R1) (initialized)
*              - uint8_t *buf: array for XOF input (allocated CHAL1_ISS_INPUT_BYTES bytes)
*              - poly_qiss_vec_m2_d s2_1: polynomial vectors to host ABDLOP commitment randomness (initialized)
*              - poly_qiss_vec_d s2_2: polynomial vectors to host ABDLOP commitment randomness (initialized)
*              - poly_qiss_vec_m1 y1: polynomial vectors to host mask for c.s1 (initialized)
*              - poly_qiss_vec_m2_d y2_1: polynomial vectors to host mask for c.s2 (initialized)
*              - poly_qiss_vec_d y2_2: polynomial vectors to host mask for c.s2 (initialized)
*              - poly_qiss_vec_256_l y3_g: polynomial vectors to host mask for R.s1 and automorphism eq. (initialized)
*              - poly_qiss_vec_d tA_L: polynomial vector to host the lower-order bits of witness commitment tA (initialized)
*              - poly_qiss_vec_d w_L: polynomial vector to host the lower-order bits of mask commitment w (initialized)
*              - poly_qiss_vec_d w_H: polynomial vector to host the high-order bits of mask commitment w (initialized)
*              - const poly_qiss_vec_m1 s1: polynomial vectors, witness
*              - const poly_qiss_mat_d_m1 A1: polynomial matrix for A1 (CRS)
*              - const poly_qiss_mat_d_m2_d A2: polynomial matrix for A2 (CRS)
*              - const poly_qiss_mat_256l_m2_d Byg: polynomial matrix for B_{y,g} (CRS)
*              - const uint8_t *randomness_seed: seed to extract randomness from (allocated SEED_BYTES bytes)
*              - const uint32_t kappa: XOF domain separator due to rejections 
**************************************************/
static void prove_1_round1(
    proof_1_t                   *proof_1,
    poly_qiss_vec_m1            chal_11[PARAM_ARP_ISS],
    poly_qiss_vec_k             chal_12[PARAM_ARP_ISS][PARAM_DE+1],
    uint8_t                     buf[CHAL1_ISS_INPUT_BYTES], 
    poly_qiss_vec_m2_d          s2_1, 
    poly_qiss_vec_d             s2_2, 
    poly_qiss_vec_m1            y1, 
    poly_qiss_vec_m2_d          y2_1, 
    poly_qiss_vec_d             y2_2, 
    poly_qiss_vec_256_l         y3_g,
    poly_qiss_vec_d             tA_L,
    poly_qiss_vec_d             w_H,
    poly_qiss_vec_d             w_L,
    const poly_qiss_vec_m1      s1,
    const poly_qiss_mat_d_m1    A1,
    const poly_qiss_mat_d_m2_d  A2,
    const poly_qiss_mat_256l_m2_d Byg,
    const uint8_t               randomness_seed[SEED_BYTES],
    const uint32_t              kappa) {
  size_t i,j;
  uint32_t kpp = kappa;
  poly_qiss_vec_d tmp_vec_d;
  uint8_t challenge_seed[SEED_BYTES];

  // init vectors
  poly_qiss_vec_d_init(tmp_vec_d);

  // sampling ABDLOP commitment randomness
  poly_qiss_vec_m2_d_binomial(s2_1, s2_2, randomness_seed, kpp++, DOMAIN_SEPARATOR_RAND_S2_ISS);

  // sampling Gaussian masks for c.s1, c.s2, R.s1
  poly_qiss_vec_m1_sample_gaussian_s1(y1);
  poly_qiss_vec_m2_d_sample_gaussian_s2(y2_1, y2_2);
  for (i = 0; i < PARAM_ARP_DIV_N_ISS; i++) {
    for (j = 0; j < PARAM_N_ISS; j++) {
      poly_qiss_set_coeff(y3_g->entries[i], j, SampleZ(0, PARAM_S3_ISS));
    }
  }

  // sampling uniform mask for equations with automorphisms
  for (i = PARAM_ARP_DIV_N_ISS; i < PARAM_ARP_DIV_N_L_ISS; i++) {
    poly_qiss_uniform_but_zero_half(y3_g->entries[i], randomness_seed, kpp++, DOMAIN_SEPARATOR_RAND_G_ISS);
  }

  // computing commitments tA = A1.s1 + A2.s2_1 + s2_2, w = A1.y1 + A2.y2_1 + y2_2, tB = Byg.s2_2 + y3_g
  // tA
  poly_qiss_mat_d_m1_mul_vec_m1(proof_1->tA_H, A1, s1);
  poly_qiss_mat_d_m2_d_mul_vec_m2_d(tmp_vec_d, A2, s2_1);
  poly_qiss_vec_d_add(proof_1->tA_H, proof_1->tA_H, tmp_vec_d);
  poly_qiss_vec_d_add(proof_1->tA_H, proof_1->tA_H, s2_2);
  // w
  poly_qiss_mat_d_m1_mul_vec_m1(w_H, A1, y1);
  poly_qiss_mat_d_m2_d_mul_vec_m2_d(tmp_vec_d, A2, y2_1);
  poly_qiss_vec_d_add(w_H, w_H, tmp_vec_d);
  poly_qiss_vec_d_add(w_H, w_H, y2_2);
  // tB
  poly_qiss_mat_256l_m2_d_mul_vec_m2_d(proof_1->tB, Byg, s2_1); 
  poly_qiss_vec_256_l_add(proof_1->tB, proof_1->tB, y3_g);

  // rounding
  poly_qiss_vec_d_power2round(proof_1->tA_H, tA_L, proof_1->tA_H);
  poly_qiss_vec_d_decompose(w_H, w_L, w_H);

  // computing first challenge
  buf[0] = 1;
  poly_qiss_vec_d_pack(buf + ISS_CHALLENGE_BASE_BYTES, proof_1->tA_H);
  poly_qiss_vec_d_pack(buf + ISS_CHALLENGE_BASE_BYTES + POLYQISS_VECD_PACKEDBYTES, w_H);
  poly_qiss_vec_256_l_pack(buf + ISS_CHALLENGE_BASE_BYTES + 2*POLYQISS_VECD_PACKEDBYTES, proof_1->tB);
  shake256(challenge_seed, SEED_BYTES, buf, CHAL1_ISS_INPUT_BYTES);
  for (i = 0; i < PARAM_ARP_ISS; i++) {
    poly_qiss_vec_m1_vec_k_binomial(chal_11[i], chal_12[i], challenge_seed, DOMAIN_SEPARATOR_CHAL1_ISS, i, SEED_BYTES);
  }

  // clean up vectors
  poly_qiss_vec_d_clear(tmp_vec_d);
}

/*************************************************
* Name:        prove_1_round2 [static]
*
* Description: Compute round 2 of zero-knowledge proof of commitment opening
*              and verifiable encryption
*
* Arguments:   - proof_1_t *proof_1: pointer to issuance proof structure
*              - coeff_qiss *chal_2: array of coeff_qiss to host second challenge (gamma_{i,j}) (allocated)
*              - int64_t *z3_Rx_inner: pointer to coeff_qiss to host <z_3, Rx> (allocated)
*              - uint8_t *buf: array for XOF input (allocated CHAL2_ISS_INPUT_BYTES bytes)
*              - const poly_qiss_vec_k *s1: array of polynomial vectors, witness
*              - const poly_qiss_vec_k *chal_1: array of polynomial vectors, first challenge R = R0 - R1
*              - const poly_qiss_vec_256_l y3_g: polynomial vector, mask for R.s1 (ARP)
* 
* Returns the 64-bit unsigned integer corresponding to ||z_3 - y_3||^2
**************************************************/
static uint64_t prove_1_round2(
    proof_1_t                 *proof_1,
    coeff_qiss                chal_2[2*PARAM_L_ISS][PARAM_ARP_ISS + 4 + PARAM_N_ISS - 1],
    int64_t                   *z3_Rx_inner,
    uint8_t                   buf[CHAL2_ISS_INPUT_BYTES], 
    const poly_qiss_vec_m1    s1,
    const poly_qiss_vec_k     enc_arp_embed[PARAM_DE+1],
    const poly_qiss_vec_m1    chal_11[PARAM_ARP_ISS],
    const poly_qiss_vec_k     chal_12[PARAM_ARP_ISS][PARAM_DE+1],
    const poly_qiss_vec_256_l y3_g) {
  size_t i,j,k,l;
  uint64_t sq_norm_arp = 0;
  coeff_qiss tmp_coeff;
  uint8_t challenge_seed[SEED_BYTES];

  // computing z3 (in Z not in ring) and square norms
  *z3_Rx_inner = 0;
  for (i = 0; i < PARAM_ARP_ISS; i++) {
    tmp_coeff = 0;
    for (j = 0; j < PARAM_M1_ISS; j++) {
      for (k = 0; k < PARAM_N_ISS; k++) {
        tmp_coeff += poly_qiss_get_coeff_centered(chal_11[i]->entries[j], k) * poly_qiss_get_coeff_centered(s1->entries[j], k);
      }
    }
    for (j = 0; j < PARAM_DE+1; j++) {
      for (k = 0; k < PARAM_K_ISS; k++) {
        for (l = 0; l < PARAM_N_ISS; l++) {
          tmp_coeff += poly_qiss_get_coeff_centered(chal_12[i][j]->entries[k], l) * poly_qiss_get_coeff_centered(enc_arp_embed[j]->entries[k], l);
        }
      }
    }
    CHK_UI_OVF_ADDITION(sq_norm_arp, (uint64_t)(tmp_coeff * tmp_coeff));
    proof_1->z3[i] = poly_qiss_get_coeff_centered(y3_g->entries[i / PARAM_N_ISS], i % PARAM_N_ISS) + tmp_coeff;
    *z3_Rx_inner += (int64_t)(proof_1->z3[i] * tmp_coeff); // <z_3, Rx> should not overflow 63 bits
  }

  // computing second challenge
  buf[0] = 2;
  for (i = 0; i < PARAM_ARP_ISS; i++) {
    coeff_qiss_pack(buf + CHAL1_ISS_INPUT_BYTES + i*COEFFQISS_PACKEDBYTES, proof_1->z3[i]);
  }
  shake256(challenge_seed, SEED_BYTES, buf, CHAL2_ISS_INPUT_BYTES);
  for (i = 0; i < 2*PARAM_L_ISS; i++) {
    vec_qiss_uniform(chal_2[i], challenge_seed, DOMAIN_SEPARATOR_CHAL2_ISS, i, SEED_BYTES); // writes PARAM_ARP_ISS + 4 + PARAM_N_ISS - 1 uniform numbers to the first argument
  }
  return sq_norm_arp; // holds ||z_3 - y_3||^2
}

/*************************************************
* Name:        prove_1_round3 [static]
*
* Description: Compute round 3 of zero-knowledge proof of commitment opening
*              and verifiable encryption
*
* Arguments:   - proof_1_t *proof_1: pointer to issuance proof structure
*              - poly_qiss_vec_l chal_3_l: polynomial vector to host third challenge (mu_i)_{i < l} (initialized)
*              - poly_qiss_vec_k *chal_3_dk: array of polynomial vectors to host third challenge (mu_i)_{l <= i < dk+l} (initialized)
*              - poly_qiss chal_3_1: polynomial to host third challenge mu_{dk+l+1} (initialized)
*              - uint8_t *buf: array for XOF input (allocated CHAL3_ISS_INPUT_BYTES bytes)
*              - poly_qiss_vec_256 *sum_gamma_e_star: array polynomial vectors to host Sum_j gamma_{ij}e_j* (initialized)
*              - poly_qiss_vec_m1 *sum_gamma_r1_star: array polynomial vectors to host Sum_j gamma_{ij}r_{j,1}* (initialized)
*              - poly_qiss_vec_k *sum_gamma_r2_star: array polynomial vectors to host Sum_j gamma_{ij}r_{j,2}* (initialized)
*              - const poly_qiss_vec_m1 s1: polynomial vectors, witness
*              - const poly_qiss_vec_m1 *chal_11: array of polynomial vectors, first challenge (R = R0 - R1) 
*              - const poly_qiss_vec_k *chal_12: array of polynomial vectors, first challenge (R = R0 - R1) 
*              - const coeff_qiss *chal_2: array of coeff_qiss, second challenge (gamma_{i,j})
*              - const poly_qiss_vec_256_l y3_g: polynomial vectors, mask for arp and automorphism eq.
*              - const poly_qiss_vec_k *enc_arp_embed: array of polynomial vectors, precomputation of p^{-1}.(s_{1,b}.ct - [A_e | b_e]^T.s_{1,e} - [0 | round(p/2).s_{1,h}]^T)
*              - const poly_qiss *quadratic_precomp: array of polynomials, precomputations quadratic norm terms
**************************************************/
static void prove_1_round3(
    proof_1_t                 *proof_1,
    poly_qiss_vec_l           chal_3_l,
    poly_qiss_vec_k           chal_3_dk[PARAM_D],
    poly_qiss                 chal_3_1,
    uint8_t                   buf[CHAL3_ISS_INPUT_BYTES], 
    poly_qiss_vec_256         sum_gamma_e_star[2*PARAM_L_ISS],
    poly_qiss_vec_m1          sum_gamma_r1_star[2*PARAM_L_ISS],
    poly_qiss_vec_k           sum_gamma_r2_star[2*PARAM_L_ISS][PARAM_DE+1],
    const poly_qiss_vec_m1    s1,
    const poly_qiss_vec_m1    chal_11[PARAM_ARP_ISS],
    const poly_qiss_vec_k     chal_12[PARAM_ARP_ISS][PARAM_DE+1],
    const coeff_qiss          chal_2[2*PARAM_L_ISS][PARAM_ARP_ISS + 4 + PARAM_N_ISS - 1],
    const poly_qiss_vec_256_l y3_g,
    const poly_qiss_vec_k     enc_arp_embed[PARAM_DE+1],
    const poly_qiss           quadratic_precomp[4]) {
  size_t i,j,k;
  poly_qiss tmp_poly, tmp_i[2*PARAM_L_ISS];
  poly_qiss_vec_m1 tmp_vec_m1;
  poly_qiss_vec_k tmp_vec_k;
  uint8_t challenge_seed[SEED_BYTES];

  // init vectors and polynomials
  poly_qiss_init(tmp_poly);
  poly_qiss_vec_m1_init(tmp_vec_m1);
  poly_qiss_vec_k_init(tmp_vec_k);

  for (i = 0; i < 2*PARAM_L_ISS; i++) {
    poly_qiss_init(tmp_i[i]);
    poly_qiss_zero(tmp_i[i]);

    // sum of -gamma_{ij}z3_j
    for (j = 0; j < PARAM_ARP_ISS; j++) {
      poly_qiss_muladd_constant(tmp_i[i], chal_2[i][j], -proof_1->z3[j]);
    }

    // sum of gamma_{ij}r_{j,1}*
    poly_qiss_vec_m1_conjugate(sum_gamma_r1_star[i], chal_11[0]);
    poly_qiss_vec_m1_mul_scalar(sum_gamma_r1_star[i], sum_gamma_r1_star[i], chal_2[i][0]);
    for (j = 1; j < PARAM_ARP_ISS; j++) {
      poly_qiss_vec_m1_conjugate(tmp_vec_m1, chal_11[j]);
      poly_qiss_vec_m1_mul_scalar(tmp_vec_m1, tmp_vec_m1, chal_2[i][j]);
      poly_qiss_vec_m1_add(sum_gamma_r1_star[i], sum_gamma_r1_star[i], tmp_vec_m1);
    }

    // sum of gamma_{ij}r_{j,2}*
    for (k = 0; k < PARAM_DE+1; k++) {
      poly_qiss_vec_k_conjugate(sum_gamma_r2_star[i][k], chal_12[0][k]);
      poly_qiss_vec_k_mul_scalar(sum_gamma_r2_star[i][k], sum_gamma_r2_star[i][k], chal_2[i][0]);
      for (j = 1; j < PARAM_ARP_ISS; j++) {
        poly_qiss_vec_k_conjugate(tmp_vec_k, chal_12[j][k]);
        poly_qiss_vec_k_mul_scalar(tmp_vec_k, tmp_vec_k, chal_2[i][j]);
        poly_qiss_vec_k_add(sum_gamma_r2_star[i][k], sum_gamma_r2_star[i][k], tmp_vec_k);
      }
    }

    // sum of gamma_{ij}e_j* = conjugate(tau^-1([gamma_{i,0} | ... | gamma_{i,256}]))
    for (j = 0; j < PARAM_ARP_DIV_N_ISS; j++) {
      poly_qiss_set_coeff(sum_gamma_e_star[i]->entries[j], 0, chal_2[i][j * PARAM_N_ISS]);
      for (k = 1; k < PARAM_N_ISS; k++) {
        poly_qiss_set_coeff(sum_gamma_e_star[i]->entries[j], k, - chal_2[i][(j + 1) * PARAM_N_ISS - k]); // set conjugate directly
      }
    }

    /**********************************************
    * Computing
    *   tmp_i = - sum_j gamma_{ij}.z3_j
    *         + sum_j gamma_{ij}.e_j*.y3 
    *         + sum_j gamma_{ij}.r_{j,1}*.s1 
    *         + sum_j gamma_{ij}.r_{j,2}*.enc_arp_embed 
    *         + gamma_{i,256}.(<s_{1,1}'*, s_{1,1}'> - B_{r,1}^2)
    *         + gamma_{i,257}.(<s_{1,23}'*, s_{1,23}'> - B_{r,2}^2)
    *         + gamma_{i,258}.(<s_{1,e}'*, s_{1,e}'> - B_{r,e}^2)
    *         + gamma_{i,259}.<s_{1,h}*, s_{1,h} - s_{1,b}.one>
    *         - sum_{0 < j < n_iss} gamma_{i,259+j}.s_{1,b}.x^{n_iss - j}
    **********************************************/
    for (j = 0; j < PARAM_ARP_DIV_N_ISS; j++) { // + <sum_e_gamma_star,y3>
      poly_qiss_mul(tmp_poly, sum_gamma_e_star[i]->entries[j], y3_g->entries[j]);
      poly_qiss_add(tmp_i[i], tmp_i[i], tmp_poly);
    }
    poly_qiss_vec_m1_mul_inner(tmp_poly, sum_gamma_r1_star[i], s1); // + <sum_r1_gamma_star,s1>
    poly_qiss_add(tmp_i[i], tmp_i[i], tmp_poly);
    for (j = 0; j < PARAM_DE+1; j++) { // + <sum_r2_gamma_star,enc_arp_embed>
      poly_qiss_vec_k_mul_inner(tmp_poly, sum_gamma_r2_star[i][j], enc_arp_embed[j]);
      poly_qiss_add(tmp_i[i], tmp_i[i], tmp_poly);
    }
    for (j = 0; j < 4; j++){
      poly_qiss_mul_scalar(tmp_poly, quadratic_precomp[j], chal_2[i][PARAM_ARP_ISS + j]); // + gamma_{i,256+j}.quadratic_precomp[j]
      poly_qiss_add(tmp_i[i], tmp_i[i], tmp_poly);
    }
    for (j = 1; j < PARAM_N_ISS; j++) { // - sum_{0 < j < n_iss} gamma_{i,259+j}.s_{1,b}.x^{n_iss - j}
      poly_qiss_mul_xj(tmp_poly, s1->entries[IDX_B_ISS], PARAM_N_ISS - j);
      poly_qiss_mul_scalar(tmp_poly, tmp_poly, chal_2[i][PARAM_ARP_ISS+3+j]);
      poly_qiss_sub(tmp_i[i], tmp_i[i], tmp_poly);
    }
  }

  // computing f_i = g_i + 2^{-1}(tmp_{2i-1} + tmp_{2i-1}*) + 2^{-1}.x^{n_iss/2}.(tmp_{2i} + tmp_{2i}*) 
  for (i = 0; i < PARAM_L_ISS; i++) {
    poly_qiss_conjugate(tmp_poly, tmp_i[2*i+1]);
    poly_qiss_add(tmp_poly, tmp_poly, tmp_i[2*i+1]);
    poly_qiss_mul_xj(proof_1->f->entries[i], tmp_poly, PARAM_N_ISS/2);
    poly_qiss_conjugate(tmp_poly, tmp_i[2*i]);
    poly_qiss_add(tmp_poly, tmp_poly, tmp_i[2*i]);
    poly_qiss_add(proof_1->f->entries[i], proof_1->f->entries[i], tmp_poly);
    poly_qiss_mul_scalar(proof_1->f->entries[i], proof_1->f->entries[i], PARAM_TWO_INVMOD_Q_ISS);
    poly_qiss_add(proof_1->f->entries[i], proof_1->f->entries[i], y3_g->entries[PARAM_ARP_DIV_N_ISS + i]);
  }

  // computing third challenge
  buf[0] = 3;
  poly_qiss_vec_l_pack(buf + CHAL2_ISS_INPUT_BYTES, proof_1->f);
  shake256(challenge_seed, SEED_BYTES, buf, CHAL3_ISS_INPUT_BYTES);
  poly_qiss_vec_l_1_uniform(chal_3_l, chal_3_1, challenge_seed, DOMAIN_SEPARATOR_CHAL3_ISS, SEED_BYTES);
  for (i = 0; i < PARAM_D; i++) {
    poly_qiss_vec_k_uniform(chal_3_dk[i], challenge_seed, DOMAIN_SEPARATOR_CHAL3_ISS, i+PARAM_L_ISS+1, SEED_BYTES);
  }

  // clean up vectors and polynomials
  poly_qiss_clear(tmp_poly);
  poly_qiss_vec_k_clear(tmp_vec_k);
  poly_qiss_vec_m1_clear(tmp_vec_m1);
  for (i = 0; i < 2*PARAM_L_ISS; i++) {
    poly_qiss_clear(tmp_i[i]);
  }
}

/*************************************************
* Name:        prove_1_round4 [static]
*
* Description: Compute round 4 of zero-knowledge proof of commitment opening
*              and verifiable encryption
*
* Arguments:   - proof_1_t *proof_1: pointer to issuance proof structure
*              - uint8_t *buf: array for XOF input (allocated CHAL4_ISS_INPUT_BYTES bytes)
*              - const poly_qiss_vec_m1 s1: polynomial vectors, witness
*              - const poly_qiss_vec_m1 s1_star: polynomial vectors, witness conjugate
*              - const poly_qiss_vec_m2_d s2_1: polynomial vector, ABDLOP commitment randomness
*              - const poly_qiss_mat_256l_m2_d Byg: polynomial matrix for B_{y,g} (CRS)
*              - const poly_qiss_vec_m2_d b: polynomial vector for b (CRS)
*              - const poly_qiss_mat_k_k *A_embed: array of polynomial matrices, subring embedding of q1.A'
*              - const poly_qiss_mat_k_k *B_embed: array of polynomial matrices to host subring embedding of q1.[tG - B|A3]
*              - const poly_qiss_mat_k_k *D_embed: array of polynomial matrices to host subring embedding of q1.d
*              - const poly_qiss_mat_k_k *Ae_be_embed: array of polynomial matrices to host subring embedding of [A_e | b_e]^T
*              - const poly_qiss_vec_k *cmt_embed: array of polynomial vectors to host subring embedding of q1.cmt
*              - const poly_qiss_vec_k *ct_embed: array of polynomial vectors to host subring embedding of ct
*              - const poly_qiss_vec_256 *sum_gamma_e_star: array polynomial vectors, Sum_j gamma_{ij}e_j*
*              - const poly_qiss_vec_m1 *sum_gamma_r1_star: array polynomial vectors, Sum_j gamma_{ij}r_{j,1}*
*              - const poly_qiss_vec_m1 *sum_gamma_r2_star: array polynomial vectors, Sum_j gamma_{ij}r_{j,2}*
*              - const poly_qiss_vec_m1 y1: polynomial vectors, mask for c.s1
*              - const poly_qiss_vec_m2_d y2_1: polynomial vector, mask for c.s2
*              - const coeff_qiss *chal_2: array of coeff_qiss, second challenge (gamma_{i,j})
*              - const poly_qiss_vec_l chal_3_l: polynomial vector, third challenge (mu_i)_{i < l} 
*              - const poly_qiss_vec_k *chal_3_dk: array of polynomial vectors, third challenge (mu_i)_{l <= i < dk+l}
*              - const poly_qiss chal_3_1: array of polynomial vectors, third challenge mu_{dk+l+1}
*              - const poly_qiss s1b_one: polynomial s_{1,b}.one
*              - const poly_qiss one: polynomial with all ones
**************************************************/
static void prove_1_round4(
    proof_1_t                     *proof_1,
    uint8_t                       buf[CHAL4_ISS_INPUT_BYTES],
    const poly_qiss_vec_m1        s1,
    const poly_qiss_vec_m1        s1_star,
    const poly_qiss_vec_m2_d      s2_1,
    const poly_qiss_mat_256l_m2_d Byg, 
    const poly_qiss_vec_m2_d      b,
    const poly_qiss_mat_k_k       A_embed[PARAM_D][PARAM_D], 
    const poly_qiss_mat_k_k       B_embed[PARAM_D][PARAM_K*(PARAM_D+1)], 
    const poly_qiss_mat_k_k       D_embed[PARAM_D], 
    const poly_qiss_mat_k_k       Ae_be_embed[PARAM_DE+1][PARAM_ME],
    const poly_qiss_vec_k         cmt_embed[PARAM_D], 
    const poly_qiss_vec_k         ct_embed[PARAM_DE+1], 
    const poly_qiss_vec_256       sum_gamma_e_star[2*PARAM_L_ISS],
    const poly_qiss_vec_m1        sum_gamma_r1_star[2*PARAM_L_ISS],
    const poly_qiss_vec_k         sum_gamma_r2_star[2*PARAM_L_ISS][PARAM_DE+1],
    const poly_qiss_vec_m1        y1,
    const poly_qiss_vec_m2_d      y2_1,
    const coeff_qiss              chal_2[2*PARAM_L_ISS][PARAM_ARP_ISS + 4 + PARAM_N_ISS - 1],
    const poly_qiss_vec_l         chal_3_l,
    const poly_qiss_vec_k         chal_3_dk[PARAM_D],
    const poly_qiss               chal_3_1,
    const poly_qiss               s1b_one,
    const poly_qiss               one) {
  size_t i,j;
  uint32_t kappa_c;
  poly_qiss sum_mu_gamma[4], tmp_poly, e0, e1, y1s_y1, y1s_s1, y1i_star, y1b_one, t0;
  poly_qiss_vec_k enc_arp_y1[PARAM_DE+1];
  poly_qiss_vec_256_l tmp_vec_256_l;
  poly_qiss_vec_l c_i, c_i_prime;
  uint8_t challenge_seed[SEED_BYTES];
  
  // init vectors and polynomials
  for (i = 0; i < PARAM_DE+1; i++) {
    poly_qiss_vec_k_init(enc_arp_y1[i]);
  }
  poly_qiss_vec_256_l_init(tmp_vec_256_l);
  poly_qiss_vec_l_init(c_i);
  poly_qiss_vec_l_init(c_i_prime);
  for (i = 0; i < 4; i++) {
    poly_qiss_init(sum_mu_gamma[i]);
    poly_qiss_zero(sum_mu_gamma[i]);
  }
  poly_qiss_init(tmp_poly);
  poly_qiss_init(e0);
  poly_qiss_init(e1);
  poly_qiss_init(y1s_y1);
  poly_qiss_init(y1s_s1);
  poly_qiss_init(y1i_star);
  poly_qiss_init(y1b_one);
  poly_qiss_init(t0);

  // Computing sum_i mu_i (gamma_{2i,256+j} + x^(n_iss/2).gamma_{2i+1,256+j}) and c_i, c_i'
  for (i = 0; i < PARAM_L_ISS; i++) {
    poly_qiss_mul_xj(e0, chal_3_l->entries[i], PARAM_N_ISS/2); // using e0 as tmp variable
    for (j = 0; j < 4; j++) {
      poly_qiss_mul_scalar(tmp_poly, chal_3_l->entries[i], chal_2[2*i][PARAM_ARP_ISS + j]); // mu_i gamma_{2i,256+j}
      poly_qiss_add(sum_mu_gamma[j], sum_mu_gamma[j], tmp_poly);
      poly_qiss_mul_scalar(tmp_poly, e0, chal_2[2*i+1][PARAM_ARP_ISS + j]); // mu_i.x^(n_iss/2).gamma_{2i+1,256+j}
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
  }
  poly_qiss_mul_scalar(sum_mu_gamma[3], sum_mu_gamma[3], PARAM_TWO_INVMOD_Q_ISS);

  // p^{-1}.(y_{1,b}.ct_embed - Ae_be_embed.y_{1,e} - [0 | round(p/2).y_{1,h}]^T) // y1i_star used as temp variable
  for (i = 0; i < (PARAM_DE+1)*PARAM_K_ISS; i++) {
    if (i < PARAM_DE*PARAM_K_ISS) {
      poly_qiss_zero(tmp_poly);
    } else {
      poly_qiss_mul_scalar(tmp_poly, y1->entries[IDX_MSG_ISS + (i % PARAM_K_ISS)], PARAM_HALF_P); // round(p/2)*y_{1,h}
    }
    for (j = 0; j < PARAM_ME*PARAM_K_ISS; j++) {
      poly_qiss_mul(y1i_star, Ae_be_embed[i / PARAM_K_ISS][j / PARAM_K_ISS]->rows[i % PARAM_K_ISS]->entries[j % PARAM_K_ISS], y1->entries[IDX_RE_ISS + j]); 
      poly_qiss_add(tmp_poly, tmp_poly, y1i_star);
    }
    // - y_{1,b}.ct'
    poly_qiss_mul(y1i_star, ct_embed[i / PARAM_K_ISS]->entries[i % PARAM_K_ISS], y1->entries[IDX_B_ISS]);
    poly_qiss_sub(tmp_poly, tmp_poly, y1i_star);
    // .-p^{-1}
    poly_qiss_mul_scalar(tmp_poly, tmp_poly, PARAM_P_INVMOD_Q_ISS_NEG);

    poly_qiss_set(enc_arp_y1[i / PARAM_K_ISS]->entries[i % PARAM_K_ISS], tmp_poly);
  }

  /**********************************************
  * Computing quadratic terms in e0 and e1
  *  For e0: sum_mu_gamma[0].<y_{1,1}'*, y_{1,1}'> + sum_mu_gamma[1].<y_{1,23}'*, y_{1,23}'> 
  *          + sum_mu_gamma[2].<y_{1,e}'*, y_{1,e}'> + sum_mu_gamma[3].(<y_{1,h}*, y_{1,h} - y_{1,b}.one> + <y_{1,h}*, y_{1,h} - y_{1,b}.one>*)
  *          + mu_{l+dk+1}.y_{1,b}^2
  *   
  *  For e1: sum_mu_gamma[0].(<y_{1,1}'*, s_{1,1}'> + <y_{1,1}'*, s_{1,1}'>*)
  *          + sum_mu_gamma[1].(<y_{1,23}'*, s_{1,23}'> + <y_{1,23}'*, s_{1,23}'>*)
  *          + sum_mu_gamma[2].(<y_{1,e}'*, s_{1,e}'> + <y_{1,e}'*, s_{1,e}'>*)
  *          + sum_mu_gamma[3].(<y_{1,h}*, s_{1,h} - s_{1,b}.one> + <y_{1,h}*, s_{1,h} - s_{1,b}.one>* + <s_{1,h}*, y_{1,h} - y_{1,b}.one> + <s_{1,h}*, y_{1,h} - y_{1,b}.one>*)
  *          + 2.mu_{l+dk+1}.s_{1,b}.y_{1.b}
  **********************************************/
  // {1,1}
  poly_qiss_zero(y1s_y1);
  poly_qiss_zero(y1s_s1);
  for (i = IDX_R1_ISS; i < IDX_R23_ISS; i++) {
    poly_qiss_conjugate(y1i_star, y1->entries[i]);
    poly_qiss_mul(tmp_poly, y1i_star, y1->entries[i]);
    poly_qiss_add(y1s_y1, y1s_y1, tmp_poly);
    poly_qiss_mul(tmp_poly, y1i_star, s1->entries[i]);
    poly_qiss_add(y1s_s1, y1s_s1, tmp_poly);
  }
  poly_qiss_mul(e0, sum_mu_gamma[0], y1s_y1);
  poly_qiss_conjugate(tmp_poly, y1s_s1);
  poly_qiss_add(y1s_s1, y1s_s1, tmp_poly); // <y_{1,1}'*,s_{1,1}'> + <y_{1,1}'*,s_{1,1}'>*
  poly_qiss_mul(e1, sum_mu_gamma[0], y1s_s1);
  // {1,23}
  poly_qiss_zero(y1s_y1);
  poly_qiss_zero(y1s_s1);
  for (i = IDX_R23_ISS; i < IDX_RE_ISS; i++) {
    poly_qiss_conjugate(y1i_star, y1->entries[i]);
    poly_qiss_mul(tmp_poly, y1i_star, y1->entries[i]);
    poly_qiss_add(y1s_y1, y1s_y1, tmp_poly);
    poly_qiss_mul(tmp_poly, y1i_star, s1->entries[i]);
    poly_qiss_add(y1s_s1, y1s_s1, tmp_poly);
  }
  poly_qiss_mul(tmp_poly, sum_mu_gamma[1], y1s_y1);
  poly_qiss_add(e0, e0, tmp_poly); 
  poly_qiss_conjugate(tmp_poly, y1s_s1);
  poly_qiss_add(y1s_s1, y1s_s1, tmp_poly); // <y_{1,23}'*,s_{1,23}'> + <y_{1,23}'*,s_{1,23}'>*
  poly_qiss_mul(tmp_poly, sum_mu_gamma[1], y1s_s1);
  poly_qiss_add(e1, e1, tmp_poly);
  // {1,e}
  poly_qiss_zero(y1s_y1);
  poly_qiss_zero(y1s_s1);
  for (i = IDX_RE_ISS; i < IDX_MSG_ISS; i++) {
    poly_qiss_conjugate(y1i_star, y1->entries[i]);
    poly_qiss_mul(tmp_poly, y1i_star, y1->entries[i]);
    poly_qiss_add(y1s_y1, y1s_y1, tmp_poly);
    poly_qiss_mul(tmp_poly, y1i_star, s1->entries[i]);
    poly_qiss_add(y1s_s1, y1s_s1, tmp_poly);
  }
  poly_qiss_mul(tmp_poly, sum_mu_gamma[2], y1s_y1);
  poly_qiss_add(e0, e0, tmp_poly); 
  poly_qiss_conjugate(tmp_poly, y1s_s1);
  poly_qiss_add(y1s_s1, y1s_s1, tmp_poly); // <y_{1,e}'*,s_{1,e}'> + <y_{1,e}'*,s_{1,e}'>*
  poly_qiss_mul(tmp_poly, sum_mu_gamma[2], y1s_s1);
  poly_qiss_add(e1, e1, tmp_poly);
  // one term
  poly_qiss_zero(y1s_y1);
  poly_qiss_zero(y1s_s1);
  poly_qiss_mul(y1b_one, y1->entries[IDX_B_ISS], one);
  for (i = IDX_MSG_ISS; i < IDX_B_ISS; i++) {
    poly_qiss_conjugate(y1i_star, y1->entries[i]);
    poly_qiss_sub(tmp_poly, s1->entries[i], s1b_one);
    poly_qiss_mul(tmp_poly, tmp_poly, y1i_star);
    poly_qiss_add(y1s_s1, y1s_s1, tmp_poly); // <y_{1,h}*, s_{1,h} - s_{1,b}.one>

    poly_qiss_sub(tmp_poly, y1->entries[i], y1b_one);
    poly_qiss_mul(y1i_star, y1i_star, tmp_poly); // y1i_star used as tmp variable
    poly_qiss_add(y1s_y1, y1s_y1, y1i_star); // <y_{1,h}*, y_{1,h} - y_{1,b}.one>

    poly_qiss_mul(tmp_poly, s1_star->entries[i], tmp_poly);
    poly_qiss_add(y1s_s1, y1s_s1, tmp_poly); // + <s_{1,h}*, y_{1,h} - y_{1,b}.one>
  }
  poly_qiss_conjugate(tmp_poly, y1s_y1);
  poly_qiss_add(y1s_y1, y1s_y1, tmp_poly); // (<y_{1,h}*, y_{1,h} - y_{1,b}.one> + <y_{1,h}*, y_{1,h} - y_{1,b}.one>*)
  poly_qiss_mul(tmp_poly, sum_mu_gamma[3], y1s_y1);
  poly_qiss_add(e0, e0, tmp_poly); 
  poly_qiss_conjugate(tmp_poly, y1s_s1);
  poly_qiss_add(y1s_s1, y1s_s1, tmp_poly); // (<y_{1,h}*, s_{1,h} - s_{1,b}.one> + <y_{1,h}*, s_{1,h} - s_{1,b}.one>* + <s_{1,h}*, y_{1,h} - y_{1,b}.one> + <s_{1,h}*, y_{1,h} - y_{1,b}.one>*)
  poly_qiss_mul(tmp_poly, sum_mu_gamma[3], y1s_s1);
  poly_qiss_add(e1, e1, tmp_poly);
  // {1,b}
  poly_qiss_mul(tmp_poly, chal_3_1, y1->entries[IDX_B_ISS]);
  poly_qiss_mul(y1i_star, tmp_poly, y1->entries[IDX_B_ISS]); // y1i_star used as tmp variable
  poly_qiss_add(e0, e0, y1i_star);
  poly_qiss_mul_scalar(y1i_star, s1->entries[IDX_B_ISS], 2); // y1i_star used as tmp variable
  poly_qiss_mul(y1i_star, tmp_poly, y1i_star);
  poly_qiss_add(e1, e1, y1i_star);

  /**********************************************
  * Computing linear terms in e1:
  *  sum_{i < l} mu_i.(-[Byg.y2_1]_i - 2^{-1}.(SE_{2i} + x^{n_iss/2}.SE_{2i+1}).[Byg.y2_1*]_{:256/n_iss}
  *                    -(SE_{2i}* + x^{n_iss/2}.SE_{2i+1}*).[Byg.y2_1]_{:256/n_iss}
  *                    + (SR_{2i,1} + x^{n_iss/2}.SR_{2i+1,1}).y1* + (SR_{2i,1}* + x^{n_iss/2}.SR_{2i+1,1}*).y1
  *                    + (SR_{2i,2} + x^{n_iss/2}.SR_{2i+1,2}).enc_arp_y1* + (SR_{2i,2}* + x^{n_iss/2}.SR_{2i+1,2}*).enc_arp_y1
  *                    + c_i[i].y_{1,b}* + c_i_prime[i].y_{1,b})
  *   + sum_{i < dk} mu_{l+i}.[A_theta.y_{1,1} + B_theta.y_{1,23} + D.y_{1,h} - y_{1,b}.cmt_embed]_i
  **********************************************/
  poly_qiss_mat_256l_m2_d_mul_vec_m2_d(tmp_vec_256_l, Byg, y2_1); // B_{y,g}.y2_1
  for (i = 0; i < PARAM_L_ISS; i++) {
    // using y1s_y1 to store sum before multiplying by mu_i
    poly_qiss_mul(y1s_y1, c_i_prime->entries[i], y1->entries[IDX_B_ISS]);
    poly_qiss_conjugate(tmp_poly, y1->entries[IDX_B_ISS]);
    poly_qiss_mul(tmp_poly, tmp_poly, c_i->entries[i]);
    poly_qiss_add(y1s_y1, y1s_y1, tmp_poly);

    poly_qiss_sub(y1s_y1, y1s_y1, tmp_vec_256_l->entries[PARAM_ARP_DIV_N_ISS + i]);

    poly_qiss_zero(t0); // t0 used as temp variable
    // term in SE_i
    for (j = 0; j < PARAM_ARP_DIV_N_ISS; j++) {
      poly_qiss_mul_xj(tmp_poly, sum_gamma_e_star[2*i+1]->entries[j], PARAM_N_ISS/2); 
      poly_qiss_sub(y1s_s1, sum_gamma_e_star[2*i]->entries[j], tmp_poly); // SE_{2i}* - x^{n_iss/2}.SE_{2i+1}* // y1s_s1 temp variable
      poly_qiss_add(tmp_poly, sum_gamma_e_star[2*i]->entries[j], tmp_poly); // SE_{2i}* + x^{n_iss/2}.SE_{2i+1}*

      poly_qiss_mul(tmp_poly, tmp_poly, tmp_vec_256_l->entries[j]); // (SE_{2i}* + x^{n_iss/2}.SE_{2i+1}*).[Byg.y2_1]_{:256/n_iss}
      poly_qiss_sub(t0, t0, tmp_poly);
      poly_qiss_mul(y1s_s1, y1s_s1, tmp_vec_256_l->entries[j]); // (SE_{2i}* - x^{n_iss/2}.SE_{2i+1}*).[Byg.y2_1]_{:256/n_iss}
      poly_qiss_conjugate(tmp_poly, y1s_s1); // (SE_{2i} + x^{n_iss/2}.SE_{2i+1}).[Byg.y2_1*]_{:256/n_iss}
      poly_qiss_sub(t0, t0, tmp_poly);
    }
    // term in SR_{i,1}
    for (j = 0; j < PARAM_M1_ISS; j++) {
      poly_qiss_mul_xj(tmp_poly, sum_gamma_r1_star[2*i+1]->entries[j], PARAM_N_ISS/2); 
      poly_qiss_sub(y1s_s1, sum_gamma_r1_star[2*i]->entries[j], tmp_poly); // SR_{2i,1}* - x^{n_iss/2}.SR_{2i+1,1}* // y1s_s1 temp variable
      poly_qiss_add(tmp_poly, sum_gamma_r1_star[2*i]->entries[j], tmp_poly); // SR_{2i,1}* + x^{n_iss/2}.SR_{2i+1,1}*

      poly_qiss_mul(tmp_poly, tmp_poly, y1->entries[j]); // (SR_{2i,1}* + x^{n_iss/2}.SR_{2i+1,1}*).y1
      poly_qiss_add(t0, t0, tmp_poly);
      poly_qiss_mul(y1s_s1, y1s_s1, y1->entries[j]); // (SR_{2i,1}* - x^{n_iss/2}.SR_{2i+1,1}*).y1
      poly_qiss_conjugate(tmp_poly, y1s_s1); // (SR_{2i,1} + x^{n_iss/2}.SR_{2i+1,1}).y1*
      poly_qiss_add(t0, t0, tmp_poly);
    }
    // term in SR_{i,2}
    for (j = 0; j < (PARAM_DE+1)*PARAM_K_ISS; j++) {
      poly_qiss_mul_xj(tmp_poly, sum_gamma_r2_star[2*i+1][j / PARAM_K_ISS]->entries[j % PARAM_K_ISS], PARAM_N_ISS/2); 
      poly_qiss_sub(y1s_s1, sum_gamma_r2_star[2*i][j / PARAM_K_ISS]->entries[j % PARAM_K_ISS], tmp_poly); // SR_{2i,2}* - x^{n_iss/2}.SR_{2i+1,2}* // y1s_s1 temp variable
      poly_qiss_add(tmp_poly, sum_gamma_r2_star[2*i][j / PARAM_K_ISS]->entries[j % PARAM_K_ISS], tmp_poly); // SR_{2i,2}* + x^{n_iss/2}.SR_{2i+1,2}*

      poly_qiss_mul(tmp_poly, tmp_poly, enc_arp_y1[j / PARAM_K_ISS]->entries[j % PARAM_K_ISS]); // (SR_{2i,2}* + x^{n_iss/2}.SR_{2i+1,2}*).enc_arp_y1
      poly_qiss_add(t0, t0, tmp_poly);
      poly_qiss_mul(y1s_s1, y1s_s1, enc_arp_y1[j / PARAM_K_ISS]->entries[j % PARAM_K_ISS]); // (SR_{2i,2}* - x^{n_iss/2}.SR_{2i+1,2}*).enc_arp_y1
      poly_qiss_conjugate(tmp_poly, y1s_s1); // (SR_{2i,2} + x^{n_iss/2}.SR_{2i+1,2}).enc_arp_y1*
      poly_qiss_add(t0, t0, tmp_poly);
    }
    poly_qiss_mul_scalar(t0, t0, PARAM_TWO_INVMOD_Q_ISS);
    poly_qiss_add(y1s_y1, y1s_y1, t0);

    poly_qiss_mul(tmp_poly, y1s_y1, chal_3_l->entries[i]); // multiply by mu_i
    poly_qiss_add(e1, e1, tmp_poly); // add to e1
  }

  for (i = 0; i < PARAM_D*PARAM_K_ISS; i++) {
    // y1s_s1 temp variable
    // A_theta.y_{1,1}
    poly_qiss_mul_scalar(tmp_poly, y1->entries[i], PARAM_Q1_ISS);
    for (j = 0; j < PARAM_D*PARAM_K_ISS; j++) {
      poly_qiss_mul(y1s_s1, A_embed[i / PARAM_K_ISS][j / PARAM_K_ISS]->rows[i % PARAM_K_ISS]->entries[j % PARAM_K_ISS], y1->entries[IDX_R12_ISS + j]); 
      poly_qiss_add(tmp_poly, tmp_poly, y1s_s1);
    }
    // + B_theta.y_{1,23}
    for (j = 0; j < PARAM_K*(PARAM_D+1)*PARAM_K_ISS; j++) {
      poly_qiss_mul(y1s_s1, B_embed[i / PARAM_K_ISS][j / PARAM_K_ISS]->rows[i % PARAM_K_ISS]->entries[j % PARAM_K_ISS], y1->entries[IDX_R23_ISS + j]);
      poly_qiss_add(tmp_poly, tmp_poly, y1s_s1);
    }
    // + D_theta.y_{1,h}
    for (j = 0; j < PARAM_K_ISS; j++) {
      poly_qiss_mul(y1s_s1, D_embed[i / PARAM_K_ISS]->rows[i % PARAM_K_ISS]->entries[j], y1->entries[IDX_MSG_ISS + j]);
      poly_qiss_add(tmp_poly, tmp_poly, y1s_s1);
    }
    // - y_{1,b}.u
    poly_qiss_mul(y1s_s1, cmt_embed[i / PARAM_K_ISS]->entries[i % PARAM_K_ISS], y1->entries[IDX_B_ISS]);
    poly_qiss_sub(tmp_poly, tmp_poly, y1s_s1);
    // multiply by mu_{l+i}
    poly_qiss_mul(tmp_poly, tmp_poly, chal_3_dk[i / PARAM_K_ISS]->entries[i % PARAM_K_ISS]);
    poly_qiss_add(e1, e1, tmp_poly);
  }

  // committing to garbage terms
  poly_qiss_vec_m2_d_mul_inner(t0, b, y2_1);
  poly_qiss_add(t0, t0, e0);
  poly_qiss_vec_m2_d_mul_inner(tmp_poly, b, s2_1);
  poly_qiss_add(proof_1->t1, tmp_poly, e1);

  // computing fourth challenge
  kappa_c = 0;
  buf[0] = 4;
  poly_qiss_pack(buf + CHAL3_ISS_INPUT_BYTES, t0);
  poly_qiss_pack(buf + CHAL3_ISS_INPUT_BYTES + POLYQISS_PACKEDBYTES, proof_1->t1);
  shake256(challenge_seed, SEED_BYTES, buf, CHAL4_ISS_INPUT_BYTES);
  do {
    poly_qiss_sample_challenge(proof_1->c, challenge_seed, DOMAIN_SEPARATOR_CHAL4_ISS, kappa_c++, SEED_BYTES);
  } while (challenge_size_iss(proof_1->c) > PARAM_ETA_ISS);
  proof_1->ctr_c = kappa_c - 1;

  // clean up vectors and polynomials
  for (i = 0; i < PARAM_DE+1; i++) {
    poly_qiss_vec_k_clear(enc_arp_y1[i]);
  }
  poly_qiss_vec_256_l_clear(tmp_vec_256_l);
  poly_qiss_vec_l_clear(c_i);
  poly_qiss_vec_l_clear(c_i_prime);
  for (i = 0; i < 4; i++) {
    poly_qiss_clear(sum_mu_gamma[i]);
  }
  poly_qiss_clear(tmp_poly);
  poly_qiss_clear(e0);
  poly_qiss_clear(e1);
  poly_qiss_clear(y1s_y1);
  poly_qiss_clear(y1s_s1);
  poly_qiss_clear(y1i_star);
  poly_qiss_clear(y1b_one);
  poly_qiss_clear(t0);
}

/*************************************************
* Name:        prove_1_round5 [static]
*
* Description: Compute round 5 of zero-knowledge proof of commitment opening
*              and verifiable encryption
*
* Arguments:   - proof_1_t *proof_1: pointer to issuance proof structure
*              - const poly_qiss_vec_m1 s1: polynomial vectors, witness
*              - const poly_qiss_vec_m2_d s2_1: polynomial vector, ABDLOP commitment randomness
*              - const poly_qiss_vec_d s2_2: polynomial vector, ABDLOP commitment randomness
*              - const poly_qiss_vec_m1 y1: polynomial vectors, mask for c.s1
*              - const poly_qiss_vec_m2_d y2_1: polynomial vector, mask for c.s2_1
*              - const poly_qiss_vec_d y2_2: polynomial vector, mask for c.s2_2
*              - const poly_qiss_vec_d tA_L: polynomial vector, low bits of witness commitment w
*              - const poly_qiss_vec_d w_H: polynomial vector, high bits of mask commitment w
*              - const poly_qiss_vec_d w_L: polynomial vector, low bits of mask commitment w
*              - const uint64_t sq_norm_arp: precomputed square norm of Rx
*              - const int64_t z3_Rx_inner: precomputed inner product <z3, Rx>
* 
* Returns 1 if round 5 passes, and 0 if it rejects
**************************************************/
static int prove_1_round5(
    proof_1_t                 *proof_1,
    const poly_qiss_vec_m1    s1,
    const poly_qiss_vec_m2_d  s2_1,
    const poly_qiss_vec_d     s2_2,
    const poly_qiss_vec_m1    y1,
    const poly_qiss_vec_m2_d  y2_1,
    const poly_qiss_vec_d     y2_2,
    const poly_qiss_vec_d     tA_L,
    const poly_qiss_vec_d     w_H,
    const poly_qiss_vec_d     w_L,
    const uint64_t            sq_norm_arp,
    const int64_t             z3_Rx_inner) {
  size_t i,j;
  int pass = 1;
  poly_qiss tmp_poly;
  poly_qiss_vec_d z2_2, tmp_vec_d;
  uint128 sq_norm_cs1 = 0, sq_norm_cs2 = 0;
  int128 z1_cs1_inner = 0, z2_prime_cs2_inner = 0;

  // init vectors and polynomials
  poly_qiss_init(tmp_poly);
  poly_qiss_vec_d_init(z2_2);
  poly_qiss_vec_d_init(tmp_vec_d);

  for (i = 0; i < PARAM_M1_ISS; i++) {
    poly_qiss_mul(tmp_poly, s1->entries[i], proof_1->c);
    sq_norm_cs1 += poly_qiss_sq_norm2(tmp_poly);
    poly_qiss_add(proof_1->z1->entries[i], y1->entries[i], tmp_poly);
    for (j = 0; j < PARAM_N_ISS; j++) {
      z1_cs1_inner += ((int128)poly_qiss_get_coeff_centered(proof_1->z1->entries[i], j) * (int128)poly_qiss_get_coeff_centered(tmp_poly, j));
    }
  }

  for (i = 0; i < PARAM_M2_D_ISS; i++) {
    poly_qiss_mul(tmp_poly, s2_1->entries[i], proof_1->c);
    sq_norm_cs2 += poly_qiss_sq_norm2(tmp_poly);
    poly_qiss_add(proof_1->z2_1->entries[i], y2_1->entries[i], tmp_poly);
    for (j = 0; j < PARAM_N_ISS; j++) {
      z2_prime_cs2_inner += ((int128)poly_qiss_get_coeff_centered(proof_1->z2_1->entries[i], j) * (int128)poly_qiss_get_coeff_centered(tmp_poly, j));
    }
  }
  for (i = 0; i < PARAM_D_ISS; i++) {
    poly_qiss_mul(tmp_poly, s2_2->entries[i], proof_1->c);
    sq_norm_cs2 += poly_qiss_sq_norm2(tmp_poly);
    poly_qiss_add(z2_2->entries[i], y2_2->entries[i], tmp_poly);
    for (j = 0; j < PARAM_N_ISS; j++) {
      z2_prime_cs2_inner += ((int128)poly_qiss_get_coeff_centered(z2_2->entries[i], j) * (int128)poly_qiss_get_coeff_centered(tmp_poly, j));
    }
  }

  // rejection sampling
  // sample u1 uniform in (0,1), goto proof_1_reject if u1 > exp(pi * (sq_norm_y1 - sq_norm_z1) / PARAM_S1SQ_ISS) / PARAM_REJ1_ISS
  // sample u2 uniform in (0,1), goto proof_1_reject if u2 > exp(pi * (sq_norm_y2 - sq_norm_z2) / PARAM_S2SQ_ISS) / PARAM_REJ2_ISS
  if (_reject_exp(
    expl(M_PI * ((long double)sq_norm_cs1/(long double)PARAM_S1SQ_ISS + (long double)sq_norm_arp/(long double)PARAM_S3SQ_ISS)) / ((long double)(PARAM_REJ1_ISS*PARAM_REJ3_ISS) * coshl(2*M_PI*((long double)z1_cs1_inner/(long double)PARAM_S1SQ_ISS + (long double)z3_Rx_inner/(long double)PARAM_S3SQ_ISS)))
  ) || _reject_exp(
    expl(M_PI * (long double)sq_norm_cs2/(long double)PARAM_S2SQ_ISS) / ((long double)PARAM_REJ2_ISS * coshl(2*M_PI*(long double)z2_prime_cs2_inner/(long double)PARAM_S2SQ_ISS))
  )) {
    pass = 0;
    goto prove_1_round5_cleanup;
  }

  // hint computation
  poly_qiss_vec_d_mul_poly_qiss(tmp_vec_d, tA_L, proof_1->c);
  poly_qiss_vec_d_sub(z2_2, z2_2, tmp_vec_d);
  poly_qiss_vec_d_sub(z2_2, z2_2, w_L);
  sq_norm_cs2 = poly_qiss_vec_m2_d_norm2(proof_1->z2_1) + poly_qiss_vec_d_norm2(z2_2);
  if (sq_norm_cs2 > ((uint128)PARAM_B2SQ_ISS_LOW64 + (((uint128)PARAM_B2SQ_ISS_HIGH64) << 64))) {
    pass = 0;
    goto prove_1_round5_cleanup;
  }
  poly_qiss_vec_d_mul_scalar(tmp_vec_d, w_H, PARAM_GAMMA_ISS);
  poly_qiss_vec_d_sub(tmp_vec_d, tmp_vec_d, z2_2);
  poly_qiss_vec_d_makeGhint(proof_1->hint, z2_2, tmp_vec_d);

  // clean up vectors and polynomials
prove_1_round5_cleanup:
  poly_qiss_clear(tmp_poly);
  poly_qiss_vec_d_clear(z2_2);
  poly_qiss_vec_d_clear(tmp_vec_d);
  return pass;
}

/*************************************************
* Name:        prove_1
*
* Description: Compute zero-knowledge proof of commitment opening
*              and verifiable encryption
*
* Arguments:   - proof_1_t *proof_1: pointer to issuance proof structure (initialized)
*              - const poly_qiss_mat_k_k *A_embed: array of polynomial matrices, subring embedding of q1.A'
*              - const poly_qiss_mat_k_k *B_embed: array of polynomial matrices to host subring embedding of q1.[tG - B|A3]
*              - const poly_qiss_mat_k_k *D_embed: array of polynomial matrices to host subring embedding of q1.d
*              - const poly_qiss_mat_k_k *Ae_be_embed: array of polynomial matrices to host subring embedding of [A_e | b_e]^T
*              - const poly_qiss_vec_k *cmt_embed: array of polynomial vectors to host subring embedding of q1.cmt
*              - const poly_qiss_vec_k *ct_embed: array of polynomial vectors to host subring embedding of ct
*              - const poly_qiss_vec_m1 s1: polynomial vectors, subring embedding of witness
*              - const uint8_t *crs_seed: pointer to byte array containing the CRS seed (allocated SEED_BYTES bytes)
*              - const uint8_t *seed: pointer to byte array containing the seed for public parameters (allocated SEED_BYTES bytes)
**************************************************/
void prove_1(
    proof_1_t               *proof_1, 
    const poly_qiss_mat_k_k A_embed[PARAM_D][PARAM_D], 
    const poly_qiss_mat_k_k B_embed[PARAM_D][PARAM_K*(PARAM_D+1)], 
    const poly_qiss_mat_k_k D_embed[PARAM_D], 
    const poly_qiss_mat_k_k Ae_be_embed[PARAM_DE+1][PARAM_ME],
    const poly_qiss_vec_k   cmt_embed[PARAM_D], 
    const poly_qiss_vec_k   ct_embed[PARAM_DE+1], 
    const poly_qiss_vec_m1  s1, 
    const uint8_t           crs_seed[CRS_SEED_BYTES],
    const uint8_t           seed[SEED_BYTES]) {
  size_t i,j;
  uint8_t randomness_seed[SEED_BYTES];
  uint8_t buf[CHAL4_ISS_INPUT_BYTES] = {0};
  uint32_t kappa;
  uint64_t sq_norm_arp;
  int64_t z3_Rx_inner;
  coeff_qiss tmp_coeff;
  poly_qiss tmp_poly, tmp_poly_2, one, quadratic_precomp[4];
  poly_qiss_vec_m1 y1, s1_star, sum_gamma_r1_star[2*PARAM_L_ISS];
  poly_qiss_vec_256 sum_gamma_e_star[2*PARAM_L_ISS];
  poly_qiss_vec_k sum_gamma_r2_star[2*PARAM_L_ISS][PARAM_DE+1], enc_arp_embed[PARAM_DE+1];
  poly_qiss_vec_m2_d s2_1, y2_1, b;
  poly_qiss_vec_d s2_2, y2_2, tA_L, w_H, w_L;
  poly_qiss_vec_256_l y3_g;
  poly_qiss_mat_d_m1 A1;
  poly_qiss_mat_d_m2_d A2;
  poly_qiss_mat_256l_m2_d Byg;

  // challenges
  poly_qiss_vec_m1 chal_11[PARAM_ARP_ISS];
  poly_qiss_vec_k chal_12[PARAM_ARP_ISS][PARAM_DE+1];
  coeff_qiss chal_2[2*PARAM_L_ISS][PARAM_ARP_ISS + 4 + PARAM_N_ISS - 1];
  poly_qiss_vec_l chal_3_l;
  poly_qiss_vec_k chal_3_dk[PARAM_D];
  poly_qiss chal_3_1;

  // init
  // init polynomials
  poly_qiss_init(tmp_poly);
  poly_qiss_init(tmp_poly_2);
  poly_qiss_init(one);
  for (i = 0; i < 4; i++) {
    poly_qiss_init(quadratic_precomp[i]);
    poly_qiss_zero(quadratic_precomp[i]);
  }
  // init vectors and matrices
  poly_qiss_vec_m1_init(y1);
  poly_qiss_vec_m1_init(s1_star);
  for (i = 0; i < 2*PARAM_L_ISS; i++) {
    poly_qiss_vec_m1_init(sum_gamma_r1_star[i]);    
    poly_qiss_vec_256_init(sum_gamma_e_star[i]);    
    for (j = 0; j < PARAM_DE+1; j++) {
      poly_qiss_vec_k_init(sum_gamma_r2_star[i][j]);
    }
  }
  for (i = 0; i < PARAM_DE+1; i++) {
    poly_qiss_vec_k_init(enc_arp_embed[i]);
  }
  poly_qiss_vec_m2_d_init(s2_1);
  poly_qiss_vec_m2_d_init(y2_1);
  poly_qiss_vec_m2_d_init(b);
  poly_qiss_vec_d_init(s2_2);
  poly_qiss_vec_d_init(y2_2);
  poly_qiss_vec_d_init(tA_L);
  poly_qiss_vec_d_init(w_H);
  poly_qiss_vec_d_init(w_L);
  poly_qiss_vec_256_l_init(y3_g);
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

  // generate random secret seed
  randombytes(randomness_seed, SEED_BYTES);

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

  // precomputations
  // s1*
  poly_qiss_vec_m1_conjugate(s1_star, s1);
  // p^{-1}.(b.ct - [A_e | b_e]^T.r_e - [0 | round(p/2).H(m)]^T)
  for (i = 0; i < (PARAM_DE+1)*PARAM_K_ISS; i++) {
    if (i < PARAM_DE*PARAM_K_ISS) {
      poly_qiss_zero(tmp_poly);
    } else {
      poly_qiss_mul_scalar(tmp_poly, s1->entries[IDX_MSG_ISS + (i % PARAM_K_ISS)], PARAM_HALF_P); // round(p/2)*H(m)
    }
    for (j = 0; j < PARAM_ME*PARAM_K_ISS; j++) {
      poly_qiss_mul(tmp_poly_2, Ae_be_embed[i / PARAM_K_ISS][j / PARAM_K_ISS]->rows[i % PARAM_K_ISS]->entries[j % PARAM_K_ISS], s1->entries[IDX_RE_ISS + j]);
      poly_qiss_add(tmp_poly, tmp_poly, tmp_poly_2);
    }
    // - s_{1,b}.ct'
    poly_qiss_mul(tmp_poly_2, ct_embed[i / PARAM_K_ISS]->entries[i % PARAM_K_ISS], s1->entries[IDX_B_ISS]);
    poly_qiss_sub(tmp_poly, tmp_poly, tmp_poly_2);

    // divide by -p in Z
    for (j = 0; j < PARAM_N_ISS; j++) {
      tmp_coeff = poly_qiss_get_coeff_centered(tmp_poly, j);
      poly_qiss_set_coeff(tmp_poly, j, -tmp_coeff/(coeff_qiss)PARAM_P);
    }
    poly_qiss_set(enc_arp_embed[i / PARAM_K_ISS]->entries[i % PARAM_K_ISS], tmp_poly);
  }
  // <s_{1,1}'*,s_{1,1}'> - B_{r,1}^2
  for (i = IDX_R1_ISS; i < IDX_R23_ISS; i++) {
    poly_qiss_mul(tmp_poly, s1_star->entries[i], s1->entries[i]);
    poly_qiss_add(quadratic_precomp[0], quadratic_precomp[0], tmp_poly);
  }
  poly_qiss_muladd_constant(quadratic_precomp[0], -PARAM_B_R1_SQ, 1); 
  // <s_{1,23}'*,s_{1,23}'> - B_{r,2}^2
  for (i = IDX_R23_ISS; i < IDX_RE_ISS; i++) {
    poly_qiss_mul(tmp_poly, s1_star->entries[i], s1->entries[i]);
    poly_qiss_add(quadratic_precomp[1], quadratic_precomp[1], tmp_poly);
  }
  poly_qiss_muladd_constant(quadratic_precomp[1], -PARAM_B_R2_SQ, 1); 
  // <s_{1,e}'*,s_{1,e}'> - B_{r,e}^2
  for (i = IDX_RE_ISS; i < IDX_MSG_ISS; i++) {
    poly_qiss_mul(tmp_poly, s1_star->entries[i], s1->entries[i]);
    poly_qiss_add(quadratic_precomp[2], quadratic_precomp[2], tmp_poly);
  }
  poly_qiss_muladd_constant(quadratic_precomp[2], -PARAM_B_RE_SQ, 1);
  // <s_{1,h}*, s_{1,h} - s_{1,b}.one>
  tmp_coeff = poly_qiss_get_coeff_centered(s1->entries[IDX_B_ISS], 0); // s_{1,b}
  for (i = 0; i < PARAM_N_ISS; i++) {
    poly_qiss_set_coeff(one, i, 1);
    poly_qiss_set_coeff(tmp_poly_2, i, tmp_coeff); // tmp_poly_2 contains s_{1,b}.one
  }
  for (i = IDX_MSG_ISS; i < IDX_B_ISS; i++) {
    poly_qiss_sub(tmp_poly, s1->entries[i], tmp_poly_2);
    poly_qiss_mul(tmp_poly, s1_star->entries[i], tmp_poly);
    poly_qiss_add(quadratic_precomp[3], quadratic_precomp[3], tmp_poly);
  }

  kappa = 0;
proof_1_reject:
  /****** first round ******/
  prove_1_round1(proof_1, chal_11, chal_12, buf, s2_1, s2_2, y1, y2_1, y2_2, y3_g, tA_L, w_H, w_L, s1, A1, A2, Byg, randomness_seed, kappa);
  kappa += PARAM_ARP_DIV_N_L_ISS - PARAM_ARP_DIV_N_ISS + 1;

  /****** second round ******/
  sq_norm_arp = prove_1_round2(proof_1, chal_2, &z3_Rx_inner, buf, s1, enc_arp_embed, chal_11, chal_12, y3_g);

  /****** third round ******/
  prove_1_round3(proof_1, chal_3_l, chal_3_dk, chal_3_1, buf, sum_gamma_e_star, sum_gamma_r1_star, sum_gamma_r2_star, s1, chal_11, chal_12, chal_2, y3_g, enc_arp_embed, quadratic_precomp);
  
  /****** fourth round ******/
  prove_1_round4(proof_1, buf, s1, s1_star, s2_1, Byg, b, A_embed, B_embed, D_embed, Ae_be_embed, cmt_embed, ct_embed, sum_gamma_e_star, sum_gamma_r1_star, sum_gamma_r2_star, y1, y2_1, chal_2, chal_3_l, chal_3_dk, chal_3_1, tmp_poly_2, one);

  /****** fifth round ******/
  if (!prove_1_round5(proof_1, s1, s2_1, s2_2, y1, y2_1, y2_2, tA_L, w_H, w_L, sq_norm_arp, z3_Rx_inner)) {
    goto proof_1_reject;
  }

  // clean up
  // clear polynomials
  poly_qiss_clear(tmp_poly);
  poly_qiss_clear(tmp_poly_2);
  poly_qiss_clear(one);
  for (i = 0; i < 4; i++) {
    poly_qiss_clear(quadratic_precomp[i]);
  }
  // clean up vectors and matrices
  poly_qiss_vec_m1_clear(y1);
  poly_qiss_vec_m1_clear(s1_star);
  for (i = 0; i < 2*PARAM_L_ISS; i++) {
    poly_qiss_vec_m1_clear(sum_gamma_r1_star[i]);    
    poly_qiss_vec_256_clear(sum_gamma_e_star[i]);    
    for (j = 0; j < PARAM_DE+1; j++) {
      poly_qiss_vec_k_clear(sum_gamma_r2_star[i][j]);
    }
  }
  for (i = 0; i < PARAM_DE+1; i++) {
    poly_qiss_vec_k_clear(enc_arp_embed[i]);
  }
  poly_qiss_vec_m2_d_clear(s2_1);
  poly_qiss_vec_m2_d_clear(y2_1);
  poly_qiss_vec_m2_d_clear(b);
  poly_qiss_vec_d_clear(s2_2);
  poly_qiss_vec_d_clear(y2_2);
  poly_qiss_vec_d_clear(tA_L);
  poly_qiss_vec_d_clear(w_H);
  poly_qiss_vec_d_clear(w_L);
  poly_qiss_vec_256_l_clear(y3_g);
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
}
