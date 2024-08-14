#include "arith.h"
#include "randombytes.h"
#include "poly_qshow_sampling.h"
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
* Name:        prove_2_round1 [static]
*
* Description: Compute round 1 of zero-knowledge proof included in the blind signature
*
* Arguments:   - proof_2_t *proof_2: pointer to show proof structure
*              - poly_qshow_vec_m1 *chal_1: array of polynomial vectors to host first challenge (R = R0 - R1) (initialized)
*              - uint8_t *buf: array for XOF input (allocated CHAL1_SHOW_INPUT_BYTES bytes)
*              - poly_qshow_vec_m2_d s2_1: polynomial vectors to host ABDLOP commitment randomness (initialized)
*              - poly_qshow_vec_d s2_2: polynomial vectors to host ABDLOP commitment randomness (initialized)
*              - poly_qshow_vec_m1 y1: polynomial vectors to host mask for c.s1 (initialized)
*              - poly_qshow_vec_m2_d y2_1: polynomial vectors to host mask for c.s2 (initialized)
*              - poly_qshow_vec_d y2_2: polynomial vectors to host mask for c.s2 (initialized)
*              - poly_qshow_vec_256_l y3_g: polynomial vectors to host mask for R.s1 and automorphism eq. (initialized)
*              - poly_qshow_vec_d tA_L: polynomial vector to host the lower-order bits of witness commitment tA (initialized)
*              - poly_qshow_vec_d w_L: polynomial vector to host the lower-order bits of mask commitment w (initialized)
*              - poly_qshow_vec_d w_H: polynomial vector to host the high-order bits of mask commitment w (initialized)
*              - const poly_qshow_vec_m1 s1: polynomial vector, witness
*              - const poly_qshow_mat_d_m1 A1: polynomial matrix for A1 (CRS)
*              - const poly_qshow_mat_d_m2 A2: polynomial matrix for A2 (CRS)
*              - const poly_qshow_mat_256l_m2 Byg: polynomial matrix for B_{y,g} (CRS)
*              - const uint8_t *randomness_seed: seed to extract randomness from (allocated SEED_BYTES bytes)
*              - const uint32_t kappa: XOF domain separator due to rejections 
**************************************************/
static void prove_2_round1(
    proof_2_t                    *proof_2,
    poly_qshow_vec_m1            chal_1[PARAM_ARP_SHOW],
    uint8_t                      buf[CHAL1_SHOW_INPUT_BYTES], 
    poly_qshow_vec_m2_d          s2_1, 
    poly_qshow_vec_d             s2_2, 
    poly_qshow_vec_m1            y1, 
    poly_qshow_vec_m2_d          y2_1, 
    poly_qshow_vec_d             y2_2, 
    poly_qshow_vec_256_l         y3_g,
    poly_qshow_vec_d             tA_L,
    poly_qshow_vec_d             w_H,
    poly_qshow_vec_d             w_L,
    const poly_qshow_vec_m1      s1,
    const poly_qshow_mat_d_m1    A1,
    const poly_qshow_mat_d_m2_d  A2,
    const poly_qshow_mat_256l_m2_d Byg,
    const uint8_t                randomness_seed[SEED_BYTES],
    const uint32_t               kappa) {
  size_t i, j;
  uint32_t kpp = kappa;
  poly_qshow_vec_d tmp_vec_d;
  uint8_t binomial_seed[SEED_BYTES];

  // init vectors
  poly_qshow_vec_d_init(tmp_vec_d);

  // sampling ABDLOP commitment randomness
  poly_qshow_vec_m2_d_binomial(s2_1, s2_2, randomness_seed, kpp++, DOMAIN_SEPARATOR_RAND_S2_SHOW);

  // sampling Gaussian masks for c.s1, c.s2, R.s1
  poly_qshow_vec_m1_sample_gaussian_s1(y1);
  poly_qshow_vec_m2_d_sample_gaussian_s2(y2_1, y2_2);
  for (i = 0; i < PARAM_ARP_DIV_N_SHOW; i++) {
    for (j = 0; j < PARAM_N_SHOW; j++) {
      poly_qshow_set_coeff(y3_g->entries[i], j, SampleZ(0, PARAM_S3_SHOW));
    }
  }

  // sampling uniform mask for equations with automorphisms
  for (i = PARAM_ARP_DIV_N_SHOW; i < PARAM_ARP_DIV_N_L_SHOW; i++) {
    poly_qshow_uniform_but_zero_half(y3_g->entries[i], randomness_seed, kpp++, DOMAIN_SEPARATOR_RAND_G_SHOW);
  }

  // computing commitments tA = A1.s1 + A2.s2_1 + s2_2, w = A1.y1 + A2.y2_1 + y2_2, tB = Byg.s2_2 + y3_g
  // tA
  poly_qshow_mat_d_m1_mul_vec_m1(proof_2->tA_H, A1, s1);
  poly_qshow_mat_d_m2_d_mul_vec_m2_d(tmp_vec_d, A2, s2_1); 
  poly_qshow_vec_d_add(proof_2->tA_H, proof_2->tA_H, tmp_vec_d);
  poly_qshow_vec_d_add(proof_2->tA_H, proof_2->tA_H, s2_2);
  // w
  poly_qshow_mat_d_m1_mul_vec_m1(w_H, A1, y1);
  poly_qshow_mat_d_m2_d_mul_vec_m2_d(tmp_vec_d, A2, y2_1); 
  poly_qshow_vec_d_add(w_H, w_H, tmp_vec_d);
  poly_qshow_vec_d_add(w_H, w_H, y2_2);
  // tB
  poly_qshow_mat_256l_m2_d_mul_vec_m2_d(proof_2->tB, Byg, s2_1);
  poly_qshow_vec_256_l_add(proof_2->tB, proof_2->tB, y3_g);

  // rounding
  poly_qshow_vec_d_power2round(proof_2->tA_H, tA_L, proof_2->tA_H);
  poly_qshow_vec_d_decompose(w_H, w_L, w_H);

  // computing first challenge
  buf[0] = 1;
  poly_qshow_vec_d_pack(buf + SHOW_CHALLENGE_BASE_BYTES, proof_2->tA_H);
  poly_qshow_vec_d_pack(buf + SHOW_CHALLENGE_BASE_BYTES + POLYQSHOW_VECD_PACKEDBYTES, w_H);
  poly_qshow_vec_256_l_pack(buf + SHOW_CHALLENGE_BASE_BYTES + 2*POLYQSHOW_VECD_PACKEDBYTES, proof_2->tB);
  shake256(binomial_seed, SEED_BYTES, buf, CHAL1_SHOW_INPUT_BYTES);
  for (i = 0; i < PARAM_ARP_SHOW; i++) {
    poly_qshow_vec_m1_binomial(chal_1[i], binomial_seed, DOMAIN_SEPARATOR_CHAL1_SHOW, i, SEED_BYTES);
  }

  // clean up vectors
  poly_qshow_vec_d_clear(tmp_vec_d);
}

/*************************************************
* Name:        prove_2_round2 [static]
*
* Description: Compute round 2 of zero-knowledge proof included in the blind signature
*
* Arguments:   - proof_2_t *proof_2: pointer to show proof structure
*              - coeff_qshow *chal_2: array of coeff_qshow to host second challenge (gamma_{i,j}) (allocated)
*              - uint8_t *buf: array for XOF input (allocated CHAL2_SHOW_INPUT_BYTES bytes)
*              - const poly_qshow_vec_m1 s1: polynomial vector, witness
*              - const poly_qshow_vec_m1 *chal_1: array of polynomial vectors, first challenge R = R0 - R1
*              - const poly_qshow_vec_256_l y3_g: polynomial vector, mask for R.s1 (ARP)
* 
* Returns the 64-bit unsigned integer corresponding to ||z_3 - y_3||^2
**************************************************/
static uint64_t prove_2_round2(
    proof_2_t                  *proof_2,
    coeff_qshow                chal_2[2*PARAM_L_SHOW][PARAM_ARP_SHOW + 4 + PARAM_N_SHOW - 1],
    int64_t                    *z3_Rx_inner,
    uint8_t                    buf[CHAL2_SHOW_INPUT_BYTES], 
    const poly_qshow_vec_m1    s1,
    const poly_qshow_vec_m1    chal_1[PARAM_ARP_SHOW],
    const poly_qshow_vec_256_l y3_g) {
  size_t i,j,k;
  uint64_t sq_norm_arp = 0;
  coeff_qshow tmp_coeff;
  uint8_t challenge_seed[SEED_BYTES];

  // computing z3 (in Z not in ring) and square norms
  *z3_Rx_inner = 0;
  for (i = 0; i < PARAM_ARP_SHOW; i++) {
    tmp_coeff = 0;
    for (j = 0; j < PARAM_M1_SHOW; j++) {
      for (k = 0; k < PARAM_N_SHOW; k++) {
          tmp_coeff += poly_qshow_get_coeff_centered(chal_1[i]->entries[j], k) * poly_qshow_get_coeff_centered(s1->entries[j], k);
      }
    }
    CHK_UI_OVF_ADDITION(sq_norm_arp, (uint64_t)(tmp_coeff * tmp_coeff));
    proof_2->z3[i] = poly_qshow_get_coeff_centered(y3_g->entries[i / PARAM_N_SHOW], i % PARAM_N_SHOW) + tmp_coeff;
    *z3_Rx_inner += (int64_t)(proof_2->z3[i] * tmp_coeff); // <z_3, Rx> should not overflow 63 bits
  }

  // computing second challenge
  buf[0] = 2;
  for (i = 0; i < PARAM_ARP_SHOW; i++) {
    coeff_qshow_pack(buf + CHAL1_SHOW_INPUT_BYTES + i*COEFFQSHOW_PACKEDBYTES, proof_2->z3[i]);
  }
  shake256(challenge_seed, SEED_BYTES, buf, CHAL2_SHOW_INPUT_BYTES);
  for (i = 0; i < 2*PARAM_L_SHOW; i++) {
    vec_qshow_uniform(chal_2[i], challenge_seed, DOMAIN_SEPARATOR_CHAL2_SHOW, i, SEED_BYTES); // writes PARAM_ARP_SHOW + 4 + PARAM_N_SHOW - 1 uniform numbers to the first argument
  }
  return sq_norm_arp; // holds ||z_3 - y_3||^2
}

/*************************************************
* Name:        prove_2_round3 [static]
*
* Description: Compute round 3 of zero-knowledge proof included in the blind signature
*
* Arguments:   - proof_2_t *proof_2: pointer to show proof structure
*              - poly_qshow_vec_l chal_3_l: polynomial vector to host third challenge (mu_i)_{i < l} (initialized)
*              - poly_qshow_vec_k *chal_3_dk: array of polynomial vectors to host third challenge (mu_i)_{l <= i < dk+l} (initialized)
*              - poly_qshow chal_3_1: polynomial to host third challenge mu_{dk+l+1} (initialized)
*              - uint8_t *buf: array for XOF input (allocated CHAL3_SHOW_INPUT_BYTES bytes)
*              - poly_qshow_vec_256 *sum_gamma_e_star: array polynomial vectors to host Sum_j gamma_{ij}e_j* (initialized)
*              - poly_qshow_vec_m1 *sum_gamma_r_star: array polynomial vectors to host Sum_j gamma_{ij}r_j* (initialized)
*              - const poly_qshow_vec_m1 s1: polynomial vector, witness
*              - const poly_qshow_vec_m1 *chal_1: array of polynomial vectors, first challenge (R = R0 - R1) 
*              - const coeff_qshow *chal_2: array of coeff_qshow, second challenge (gamma_{i,j})
*              - const poly_qshow_vec_256_l y3_g: polynomial vectors, mask for R.s1 and automorphism eq.
*              - const poly_qshow *quadratic_precomp: array of polynomials containing precomputations
**************************************************/
static void prove_2_round3(
    proof_2_t                  *proof_2,
    poly_qshow_vec_l           chal_3_l,
    poly_qshow_vec_k           chal_3_dk[PARAM_D],
    poly_qshow                 chal_3_1,
    uint8_t                    buf[CHAL3_SHOW_INPUT_BYTES], 
    poly_qshow_vec_256         sum_gamma_e_star[2*PARAM_L_SHOW],
    poly_qshow_vec_m1          sum_gamma_r_star[2*PARAM_L_SHOW],
    const poly_qshow_vec_m1    s1,
    const poly_qshow_vec_m1    chal_1[PARAM_ARP_SHOW],
    const coeff_qshow          chal_2[2*PARAM_L_SHOW][PARAM_ARP_SHOW + 4 + PARAM_N_SHOW - 1],
    const poly_qshow_vec_256_l y3_g,
    const poly_qshow           quadratic_precomp[4]) {
  size_t i,j,k;
  poly_qshow tmp_poly, tmp_i[2*PARAM_L_SHOW];
  poly_qshow_vec_m1 tmp_vec_m1;
  uint8_t challenge_seed[SEED_BYTES];

  // init vectors and polynomials
  poly_qshow_init(tmp_poly);
  poly_qshow_vec_m1_init(tmp_vec_m1);

  for (i = 0; i < 2*PARAM_L_SHOW; i++) {
    poly_qshow_init(tmp_i[i]);
    poly_qshow_zero(tmp_i[i]);
    
    // sum of -gamma_{ij}z3_j
    for (j = 0; j < PARAM_ARP_SHOW; j++) {
      poly_qshow_muladd_constant(tmp_i[i], chal_2[i][j], -proof_2->z3[j]);
    }

    // sum of gamma_{ij}r_j*
    poly_qshow_vec_m1_conjugate(sum_gamma_r_star[i], chal_1[0]);
    poly_qshow_vec_m1_mul_scalar(sum_gamma_r_star[i], sum_gamma_r_star[i], chal_2[i][0]);
    for (j = 1; j < PARAM_ARP_SHOW; j++) {
      poly_qshow_vec_m1_conjugate(tmp_vec_m1, chal_1[j]);
      poly_qshow_vec_m1_mul_scalar(tmp_vec_m1, tmp_vec_m1, chal_2[i][j]);
      poly_qshow_vec_m1_add(sum_gamma_r_star[i], sum_gamma_r_star[i], tmp_vec_m1);
    }

    // sum of gamma_{ij}e_j* = conjugate(tau^-1([gamma_{i,0} | ... | gamma_{i,256}]))
    for (j = 0; j < PARAM_ARP_DIV_N_SHOW; j++) {
      poly_qshow_set_coeff(sum_gamma_e_star[i]->entries[j], 0, chal_2[i][j * PARAM_N_SHOW]);
      for (k = 1; k < PARAM_N_SHOW; k++) {
        poly_qshow_set_coeff(sum_gamma_e_star[i]->entries[j], k, - chal_2[i][(j + 1) * PARAM_N_SHOW - k]); // set conjugate directly
      }
    }

    /**********************************************
    * Computing
    *   tmp_i = - sum_j gamma_{ij}.z3_j
    *         + sum_j gamma_{ij}.e_j*.y3 
    *         + sum_j gamma_{ij}.r_j*.s1 
    *         + gamma_{i,256}.(<s_{1,1}'*, s_{1,1}'> - B_1'^2)
    *         + gamma_{i,257}.(<s_{1,23}'*, s_{1,23}'> - B_2'^2)
    *         + gamma_{i,258}.(<s_{1,t}*, s_{1,t}> - w)
    *         + gamma_{i,259}.<s_{1,t}*, s_{1,t} - s_{1,b}.one>
    *         - sum_{0 < j < n_show} gamma_{i,259+j}.s_{1,b}.x^{n_show - j}
    **********************************************/
    for (j = 0; j < PARAM_ARP_DIV_N_SHOW; j++) { // + <sum_e_gamma_star,y3>
      poly_qshow_mul(tmp_poly, sum_gamma_e_star[i]->entries[j], y3_g->entries[j]);
      poly_qshow_add(tmp_i[i], tmp_i[i], tmp_poly);
    }
    poly_qshow_vec_m1_mul_inner(tmp_poly, sum_gamma_r_star[i], s1); // + <sum_r1_gamma_star,s1>
    poly_qshow_add(tmp_i[i], tmp_i[i], tmp_poly);
    for (j = 0; j < 4; j++){
      poly_qshow_mul_scalar(tmp_poly, quadratic_precomp[j], chal_2[i][PARAM_ARP_SHOW + j]); // + gamma_{i,256+j}.quadratic_precomp[j]
      poly_qshow_add(tmp_i[i], tmp_i[i], tmp_poly);
    }
    for (j = 1; j < PARAM_N_SHOW; j++) { // - sum_{0 < j < n_show} gamma_{i,259+j}.s_{1,b}.x^{n_show - j}
      poly_qshow_mul_xj(tmp_poly, s1->entries[IDX_B_SHOW], PARAM_N_SHOW - j);
      poly_qshow_mul_scalar(tmp_poly, tmp_poly, chal_2[i][PARAM_ARP_SHOW+3+j]);
      poly_qshow_sub(tmp_i[i], tmp_i[i], tmp_poly);
    }
  }

  // computing f_i = g_i + 2^{-1}(tmp_{2i-1} + tmp_{2i-1}*) + 2^{-1}.x^{n_show/2}.(tmp_{2i} + tmp_{2i}*) 
  for (i = 0; i < PARAM_L_SHOW; i++) {
    poly_qshow_conjugate(tmp_poly, tmp_i[2*i+1]);
    poly_qshow_add(tmp_poly, tmp_poly, tmp_i[2*i+1]);
    poly_qshow_mul_xj(proof_2->f->entries[i], tmp_poly, PARAM_N_SHOW/2);
    poly_qshow_conjugate(tmp_poly, tmp_i[2*i]);
    poly_qshow_add(tmp_poly, tmp_poly, tmp_i[2*i]);
    poly_qshow_add(proof_2->f->entries[i], proof_2->f->entries[i], tmp_poly);
    poly_qshow_mul_scalar(proof_2->f->entries[i], proof_2->f->entries[i], PARAM_TWO_INVMOD_Q_SHOW);
    poly_qshow_add(proof_2->f->entries[i], proof_2->f->entries[i], y3_g->entries[PARAM_ARP_DIV_N_SHOW + i]);
  }

  // computing third challenge
  buf[0] = 3;
  poly_qshow_vec_l_pack(buf + CHAL2_SHOW_INPUT_BYTES, proof_2->f);
  shake256(challenge_seed, SEED_BYTES, buf, CHAL3_SHOW_INPUT_BYTES);
  poly_qshow_vec_l_1_uniform(chal_3_l, chal_3_1, buf, DOMAIN_SEPARATOR_CHAL3_SHOW, SEED_BYTES);
  for (i = 0; i < PARAM_D; i++) {
    poly_qshow_vec_k_uniform(chal_3_dk[i], challenge_seed, DOMAIN_SEPARATOR_CHAL3_SHOW, i+PARAM_L_SHOW+1, SEED_BYTES);
  }

  // clean up vectors and polynomials
  poly_qshow_clear(tmp_poly);
  poly_qshow_vec_m1_clear(tmp_vec_m1);
  for (i = 0; i < 2*PARAM_L_SHOW; i++) {
    poly_qshow_clear(tmp_i[i]);
  }
}

/*************************************************
* Name:        prove_2_round4 [static]
*
* Description: Compute round 4 of zero-knowledge proof included in the blind signature
*
* Arguments:   - proof_2_t *proof_2: pointer to show proof structure
*              - uint8_t *buf: array for XOF input (allocated CHAL4_SHOW_INPUT_BYTES bytes)
*              - const poly_qshow_vec_m1 s1: polynomial vector, witness
*              - const poly_qshow_vec_m1 s1_star: polynomial vector, witness conjugate
*              - const poly_qshow_vec_m2_d s2_1: polynomial vector, ABDLOP commitment randomness
*              - const poly_qshow_mat_256l_m2_d Byg: polynomial matrix for B_{y,g} (CRS)
*              - const poly_qshow_vec_m2_d b: polynomial vector for b (CRS)
*              - const poly_qshow_mat_k_k *A_embed: array of polynomial matrices, subring embedding of q1.b1.A'
*              - const poly_qshow_mat_k_k *B_embed: array of polynomial matrices to host subring embedding of q1.b2.B
*              - const poly_qshow_mat_k_k *A3_embed: array of polynomial matrices to host subring embedding of q1.b2.A3
*              - const poly_qshow_mat_k_k *G_embed: array of polynomial matrices to host subring embedding of q1.G.w_{2,L}
*              - const poly_qshow_vec_256 *sum_gamma_e_star: array polynomial vectors, Sum_j gamma_{ij}e_j*
*              - const poly_qshow_vec_m1 *sum_gamma_r_star: array polynomial vectors, Sum_j gamma_{ij}r_j*
*              - const poly_qshow_vec_m1 y1: polynomial vector, mask for c.s1
*              - const poly_qshow_vec_m2_d y2_1: polynomial vector, mask for c.s2
*              - const coeff_qshow *chal_2: array of coeff_qshow, second challenge (gamma_{i,j})
*              - const poly_qshow_vec_l chal_3_l: polynomial vector, third challenge (mu_i)_{i < l} 
*              - const poly_qshow_vec_k *chal_3_dk: array of polynomial vectors, third challenge (mu_i)_{l <= i < dk+l}
*              - const poly_qshow chal_3_1: array of polynomial vectors, third challenge mu_{dk+l+1}
*              - const poly_qshow s1b_one: polynomial s_{1,b}.one
*              - const poly_qshow one: polynomial with all ones
**************************************************/
static void prove_2_round4(
    proof_2_t                      *proof_2,
    uint8_t                        buf[CHAL4_SHOW_INPUT_BYTES],
    const poly_qshow_vec_m1        s1,
    const poly_qshow_vec_m1        s1_star,
    const poly_qshow_vec_m2_d      s2_1,
    const poly_qshow_mat_256l_m2_d Byg, 
    const poly_qshow_vec_m2_d      b,
    const poly_qshow_mat_k_k       A_embed[PARAM_D][PARAM_D], 
    const poly_qshow_mat_k_k       B_embed[PARAM_D][PARAM_K*PARAM_D], 
    const poly_qshow_mat_k_k       A3_embed[PARAM_D][PARAM_K], 
    const poly_qshow_mat_k_k       G_embed[PARAM_D], 
    const poly_qshow_vec_256       sum_gamma_e_star[2*PARAM_L_SHOW],
    const poly_qshow_vec_m1        sum_gamma_r_star[2*PARAM_L_SHOW],
    const poly_qshow_vec_m1        y1,
    const poly_qshow_vec_m2_d      y2_1,
    const coeff_qshow              chal_2[2*PARAM_L_SHOW][PARAM_ARP_SHOW + 4 + PARAM_N_SHOW - 1],
    const poly_qshow_vec_l         chal_3_l,
    const poly_qshow_vec_k         chal_3_dk[PARAM_D],
    const poly_qshow               chal_3_1,
    const poly_qshow               s1b_one,
    const poly_qshow               one) {
  size_t i,j,k;
  uint32_t kappa_c;
  int64_t bexpi;
  uint8_t challenge_seed[SEED_BYTES];
  poly_qshow sum_mu_gamma[4], tmp_poly, e0, e1, y1i_star, y1s_y1, y1s_s1, y1b_one, t0;
  poly_qshow_vec_256_l tmp_vec_256_l;
  poly_qshow_vec_l c_i, c_i_prime;
  poly_qshow_vec_k tmp_vec_k, Gy1_w2, Gs1_w2, sum_vec_k_y, sum_vec_k_s;
  poly_qshow_mat_k_k chal_3_quad_matrix;
  
  // init matrices, vectors and polynomials
  poly_qshow_init(tmp_poly);
  poly_qshow_init(e0);
  poly_qshow_init(e1);
  poly_qshow_init(y1i_star);
  poly_qshow_init(y1s_y1);
  poly_qshow_init(y1s_s1);
  poly_qshow_init(y1b_one);
  poly_qshow_init(t0);
  poly_qshow_vec_256_l_init(tmp_vec_256_l);
  poly_qshow_vec_l_init(c_i);
  poly_qshow_vec_l_init(c_i_prime);
  poly_qshow_vec_k_init(tmp_vec_k);
  poly_qshow_vec_k_init(sum_vec_k_y);
  poly_qshow_vec_k_init(sum_vec_k_s);
  poly_qshow_mat_k_k_init(chal_3_quad_matrix);
  poly_qshow_vec_k_init(Gy1_w2);
  poly_qshow_vec_k_init(Gs1_w2);
  for (i = 0; i < 4; i++) {
    poly_qshow_init(sum_mu_gamma[i]);
    poly_qshow_zero(sum_mu_gamma[i]);
  }

  // Computing sum_i mu_i (gamma_{2i,256+j} + x^(n_show/2).gamma_{2i+1,256+j}) and c_i, c_i'
  for (i = 0; i < PARAM_L_SHOW; i++) {
    poly_qshow_mul_xj(e0, chal_3_l->entries[i], PARAM_N_SHOW/2); // using e0 as tmp variable
    for (j = 0; j < 4; j++) {
      poly_qshow_mul_scalar(tmp_poly, chal_3_l->entries[i], chal_2[2*i][PARAM_ARP_SHOW + j]); // mu_i gamma_{2i,256+j}
      poly_qshow_add(sum_mu_gamma[j], sum_mu_gamma[j], tmp_poly);
      poly_qshow_mul_scalar(tmp_poly, e0, chal_2[2*i+1][PARAM_ARP_SHOW + j]); // mu_i.x^(n_show/2).gamma_{2i+1,256+j}
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
  }
  poly_qshow_add(sum_mu_gamma[2], sum_mu_gamma[2], sum_mu_gamma[3]); // sum_mu_gamma[2] = sum_i mu_i(gamma_{2i,258}+gamma_{2i,259} + x^{n_show/2}.(gamma_{2i+1,258}+gamma_{2i+1,259}))
  poly_qshow_mul_scalar(sum_mu_gamma[3], sum_mu_gamma[3], PARAM_TWO_INVMOD_Q_SHOW);

  /**********************************************
  * Computing quadratic terms in e0 and e1
  *  For e0: sum_mu_gamma[0].<y_{1,1}'*, y_{1,1}'> + sum_mu_gamma[1].<y_{1,23}'*, y_{1,23}'> + sum_mu_gamma[2].<y_{1,t}*, y_{1,t}>
  *          - sum_mu_gamma[3].(<y_{1,t}*, y_{1,b}.one> + <y_{1,t}*, y_{1,b}.one>*)
  *          + y_{1,t}.chal_3_quad_matrix.y_{1,2} + sum_{i < dk} mu_{l+i}.y_{1,b}.([A_theta.y_{1,1} - B_theta.y_{1,2} + A3_theta.y_{1,3} + G_theta.y_{1,t}]_i)
  *          + mu_{l+dk+1}.y_{1,b}^2
  *   
  *  For e1: sum_mu_gamma[0].(<y_{1,1}'*, s_{1,1}'> + <y_{1,1}'*, s_{1,1}'>*)
  *          + sum_mu_gamma[1].(<y_{1,23}'*, s_{1,23}'> + <y_{1,23}'*, s_{1,23}'>*)
  *          + sum_mu_gamma[2].(<y_{1,t}*, s_{1,t}> + <y_{1,t}*, s_{1,t}>*)
  *          - sum_mu_gamma[3].(<y_{1,t}*, s_{1,b}.one> + <y_{1,t}*, s_{1,b}.one>* + <s_{1,t}*, y_{1,b}.one> + <s_{1,t}*, y_{1,b}.one>*)
  *          + y_{1,t}.chal_3_quad_matrix.s_{1,2} + s_{1,t}.chal_3_quad_matrix.y_{1,2} 
  *          + sum_{i < dk} mu_{l+i}.y_{1,b}.([A_theta.s_{1,1} - B_theta.s_{1,2} + A3_theta.s_{1,3} + G_theta.s_{1,t}]_i)
  *          + sum_{i < dk} mu_{l+i}.s_{1,b}.([A_theta.y_{1,1} - B_theta.y_{1,2} + A3_theta.y_{1,3} + G_theta.y_{1,t}]_i)
  *          + 2.mu_{l+dk+1}.s_{1,b}.y_{1.b}
  **********************************************/
  // {1,1}
  poly_qshow_zero(y1s_y1);
  poly_qshow_zero(y1s_s1);
  for (i = IDX_W1_SHOW; i < IDX_W23_SHOW; i++) {
    poly_qshow_conjugate(y1i_star, y1->entries[i]);
    poly_qshow_mul(tmp_poly, y1i_star, y1->entries[i]);
    poly_qshow_add(y1s_y1, y1s_y1, tmp_poly);
    poly_qshow_mul(tmp_poly, y1i_star, s1->entries[i]);
    poly_qshow_add(y1s_s1, y1s_s1, tmp_poly);
  }
  poly_qshow_mul(e0, sum_mu_gamma[0], y1s_y1);
  poly_qshow_conjugate(tmp_poly, y1s_s1);
  poly_qshow_add(y1s_s1, y1s_s1, tmp_poly); // <y_{1,1}'*,s_{1,1}'> + <y_{1,1}'*,s_{1,1}'>*
  poly_qshow_mul(e1, sum_mu_gamma[0], y1s_s1);
  // {1,23}
  poly_qshow_zero(y1s_y1);
  poly_qshow_zero(y1s_s1);
  for (i = IDX_W23_SHOW; i < IDX_TAG_SHOW; i++) {
    poly_qshow_conjugate(y1i_star, y1->entries[i]);
    poly_qshow_mul(tmp_poly, y1i_star, y1->entries[i]);
    poly_qshow_add(y1s_y1, y1s_y1, tmp_poly);
    poly_qshow_mul(tmp_poly, y1i_star, s1->entries[i]);
    poly_qshow_add(y1s_s1, y1s_s1, tmp_poly);
  }
  poly_qshow_mul(tmp_poly, sum_mu_gamma[1], y1s_y1);
  poly_qshow_add(e0, e0, tmp_poly); 
  poly_qshow_conjugate(tmp_poly, y1s_s1);
  poly_qshow_add(y1s_s1, y1s_s1, tmp_poly); // <y_{1,23}'*,s_{1,23}'> + <y_{1,23}'*,s_{1,23}'>*
  poly_qshow_mul(tmp_poly, sum_mu_gamma[1], y1s_s1);
  poly_qshow_add(e1, e1, tmp_poly);
  // {1,t}
  poly_qshow_zero(y1s_y1);
  poly_qshow_zero(y1s_s1);
  for (i = IDX_TAG_SHOW; i < IDX_B_SHOW; i++) {
    poly_qshow_conjugate(y1i_star, y1->entries[i]);
    poly_qshow_mul(tmp_poly, y1i_star, y1->entries[i]);
    poly_qshow_add(y1s_y1, y1s_y1, tmp_poly);
    poly_qshow_mul(tmp_poly, y1i_star, s1->entries[i]);
    poly_qshow_add(y1s_s1, y1s_s1, tmp_poly);
  }
  poly_qshow_mul(tmp_poly, sum_mu_gamma[2], y1s_y1);
  poly_qshow_add(e0, e0, tmp_poly); 
  poly_qshow_conjugate(tmp_poly, y1s_s1);
  poly_qshow_add(y1s_s1, y1s_s1, tmp_poly); // <y_{1,t}*,s_{1,t}> + <y_{1,t}*,s_{1,t}>*
  poly_qshow_mul(tmp_poly, sum_mu_gamma[2], y1s_s1);
  poly_qshow_add(e1, e1, tmp_poly);
  // one term
  poly_qshow_zero(y1s_y1);
  poly_qshow_zero(y1s_s1);
  poly_qshow_mul(y1b_one, y1->entries[IDX_B_SHOW], one);
  for (i = IDX_TAG_SHOW; i < IDX_B_SHOW; i++) {
    poly_qshow_conjugate(y1i_star, y1->entries[i]);

    poly_qshow_mul(tmp_poly, y1i_star, s1b_one);
    poly_qshow_add(y1s_s1, y1s_s1, tmp_poly); // <y_{1,t}*, s_{1,b}.one>
    poly_qshow_mul(tmp_poly, s1_star->entries[i], y1b_one);
    poly_qshow_add(y1s_s1, y1s_s1, tmp_poly); // + <s_{1,t}*, y_{1,b}.one>

    poly_qshow_mul(tmp_poly, y1i_star, y1b_one);
    poly_qshow_add(y1s_y1, y1s_y1, tmp_poly); // <y_{1,t}*, y_{1,b}.one>
  }
  poly_qshow_conjugate(tmp_poly, y1s_y1);
  poly_qshow_add(y1s_y1, y1s_y1, tmp_poly); // (<y_{1,t}*, y_{1,b}.one> + <y_{1,t}*, y_{1,b}.one>*)
  poly_qshow_mul(tmp_poly, sum_mu_gamma[3], y1s_y1);
  poly_qshow_sub(e0, e0, tmp_poly); 
  poly_qshow_conjugate(tmp_poly, y1s_s1);
  poly_qshow_add(y1s_s1, y1s_s1, tmp_poly); // (<y_{1,t}*, s_{1,b}.one> + <s_{1,t}*, y_{1,b}.one> + <y_{1,t}*, s_{1,b}.one>* + <s_{1,t}*, y_{1,b}.one>*)
  poly_qshow_mul(tmp_poly, sum_mu_gamma[3], y1s_s1);
  poly_qshow_sub(e1, e1, tmp_poly);
  // {1,b}
  poly_qshow_mul(tmp_poly, chal_3_1, y1->entries[IDX_B_SHOW]);
  poly_qshow_mul(y1i_star, tmp_poly, y1->entries[IDX_B_SHOW]); // y1i_star used as tmp variable
  poly_qshow_add(e0, e0, y1i_star);
  poly_qshow_mul_scalar(y1i_star, s1->entries[IDX_B_SHOW], 2); // y1i_star used as tmp variable
  poly_qshow_mul(y1i_star, tmp_poly, y1i_star);
  poly_qshow_add(e1, e1, y1i_star);

  // quadratic part depending on chal_3_quad_matrix
  poly_qshow_vec_k_zero(sum_vec_k_y);
  poly_qshow_vec_k_zero(sum_vec_k_s);
  for (i = 0; i < PARAM_D; i++) {
    poly_qshow_mat_k_k_chal_3_embed(chal_3_quad_matrix, chal_3_dk[i]);
    for (j = 0; j < PARAM_K_SHOW; j++) {
      poly_qshow_zero(Gy1_w2->entries[j]);
      poly_qshow_zero(Gs1_w2->entries[j]);
      bexpi = PARAM_B2*PARAM_Q1_SHOW;
      for (k = 0; k < PARAM_K; k++) {
        poly_qshow_mul_scalar(tmp_poly, y1->entries[IDX_W23_SHOW + k*PARAM_D*PARAM_K_SHOW + i*PARAM_K_SHOW + j], bexpi);
        poly_qshow_add(Gy1_w2->entries[j], Gy1_w2->entries[j], tmp_poly);
        poly_qshow_mul_scalar(tmp_poly, s1->entries[IDX_W23_SHOW + k*PARAM_D*PARAM_K_SHOW + i*PARAM_K_SHOW + j], bexpi);
        poly_qshow_add(Gs1_w2->entries[j], Gs1_w2->entries[j], tmp_poly);
        bexpi *= PARAM_B;
      }
    }
    poly_qshow_mat_k_k_mul_vec_k(tmp_vec_k, chal_3_quad_matrix, Gy1_w2);
    poly_qshow_vec_k_add(sum_vec_k_y, sum_vec_k_y, tmp_vec_k);
    poly_qshow_mat_k_k_mul_vec_k(tmp_vec_k, chal_3_quad_matrix, Gs1_w2);
    poly_qshow_vec_k_add(sum_vec_k_s, sum_vec_k_s, tmp_vec_k);

  }
  // <y_{1,t}, (sum_i mu_{l+i}G_i").y_{1,2}>, <s_{1,t}, (sum_i mu_{l+i}G_i").y_{1,2}>, and <y_{1,t}, (sum_i mu_{l+i}G_i").s_{1,2}> and adding to e0/e1
  for (i = 0; i < PARAM_K_SHOW; i++) {
    poly_qshow_mul(tmp_poly, y1->entries[IDX_TAG_SHOW + i], sum_vec_k_y->entries[i]);
    poly_qshow_add(e0, e0, tmp_poly);
    poly_qshow_mul(tmp_poly, y1->entries[IDX_TAG_SHOW + i], sum_vec_k_s->entries[i]);
    poly_qshow_add(e1, e1, tmp_poly);
    poly_qshow_mul(tmp_poly, s1->entries[IDX_TAG_SHOW + i], sum_vec_k_y->entries[i]);
    poly_qshow_add(e1, e1, tmp_poly);
  }

  // part depending on embedded matrices A_embed, B_embed, A3_embed, G_embed
  // storing sum_{i < dk} mu_{l+i}.[A_theta.y_{1,1} - B_theta.y_{1,2} + A3_theta.y_{1,3} + G_theta.y_{1,t}]_i in y1s_y1 before multiplying by y_{1,b} / s_{1,b}
  // storing sum_{i < dk} mu_{l+i}.[A_theta.s_{1,1} - B_theta.s_{1,2} + A3_theta.s_{1,3} + G_theta.s_{1,t}]_i in y1s_s1 before multiplying by y_{1,b}
  poly_qshow_zero(y1s_y1);
  poly_qshow_zero(y1s_s1);
  for (i = 0; i < PARAM_D*PARAM_K_SHOW; i++) {
    // y1i_star used as temp variable to host [A_theta.y_{1,1} - B_theta.y_{1,2} + A3_theta.y_{1,3} + G_theta.y_{1,t}]_i
    // y1b_one used as temp variable to host [A_theta.s_{1,1} - B_theta.s_{1,2} + A3_theta.s_{1,3} + G_theta.s_{1,t}]_i
    poly_qshow_mul_scalar(y1i_star, y1->entries[i], PARAM_B1*PARAM_Q1_SHOW);
    poly_qshow_mul_scalar(y1b_one, s1->entries[i], PARAM_B1*PARAM_Q1_SHOW);
    for (j = 0; j < PARAM_D*PARAM_K_SHOW; j++) {
      poly_qshow_mul(tmp_poly, A_embed[i / PARAM_K_SHOW][j / PARAM_K_SHOW]->rows[i % PARAM_K_SHOW]->entries[j % PARAM_K_SHOW], y1->entries[IDX_W12_SHOW + j]);
      poly_qshow_add(y1i_star, y1i_star, tmp_poly);
      poly_qshow_mul(tmp_poly, A_embed[i / PARAM_K_SHOW][j / PARAM_K_SHOW]->rows[i % PARAM_K_SHOW]->entries[j % PARAM_K_SHOW], s1->entries[IDX_W12_SHOW + j]);
      poly_qshow_add(y1b_one, y1b_one, tmp_poly);
    }
    for (j = 0; j < PARAM_D*PARAM_K*PARAM_K_SHOW; j++) {
      poly_qshow_mul(tmp_poly, B_embed[i / PARAM_K_SHOW][j / PARAM_K_SHOW]->rows[i % PARAM_K_SHOW]->entries[j % PARAM_K_SHOW], y1->entries[IDX_W23_SHOW + j]);
      poly_qshow_sub(y1i_star, y1i_star, tmp_poly); // substraction
      poly_qshow_mul(tmp_poly, B_embed[i / PARAM_K_SHOW][j / PARAM_K_SHOW]->rows[i % PARAM_K_SHOW]->entries[j % PARAM_K_SHOW], s1->entries[IDX_W23_SHOW + j]);
      poly_qshow_sub(y1b_one, y1b_one, tmp_poly); // substraction
    }
    for (j = 0; j < PARAM_K*PARAM_K_SHOW; j++) {
      poly_qshow_mul(tmp_poly, A3_embed[i / PARAM_K_SHOW][j / PARAM_K_SHOW]->rows[i % PARAM_K_SHOW]->entries[j % PARAM_K_SHOW], y1->entries[IDX_W3_SHOW + j]);
      poly_qshow_add(y1i_star, y1i_star, tmp_poly);
      poly_qshow_mul(tmp_poly, A3_embed[i / PARAM_K_SHOW][j / PARAM_K_SHOW]->rows[i % PARAM_K_SHOW]->entries[j % PARAM_K_SHOW], s1->entries[IDX_W3_SHOW + j]);
      poly_qshow_add(y1b_one, y1b_one, tmp_poly);
    }
    for (j = 0; j < PARAM_K_SHOW; j++) {
      poly_qshow_mul(tmp_poly, G_embed[i / PARAM_K_SHOW]->rows[i % PARAM_K_SHOW]->entries[j], y1->entries[IDX_TAG_SHOW + j]);
      poly_qshow_add(y1i_star, y1i_star, tmp_poly); 
      poly_qshow_mul(tmp_poly, G_embed[i / PARAM_K_SHOW]->rows[i % PARAM_K_SHOW]->entries[j], s1->entries[IDX_TAG_SHOW + j]);
      poly_qshow_add(y1b_one, y1b_one, tmp_poly); 
    }
    poly_qshow_mul(y1i_star, y1i_star, chal_3_dk[i / PARAM_K_SHOW]->entries[i % PARAM_K_SHOW]);
    poly_qshow_mul(y1b_one, y1b_one, chal_3_dk[i / PARAM_K_SHOW]->entries[i % PARAM_K_SHOW]);

    poly_qshow_add(y1s_y1, y1s_y1, y1i_star);
    poly_qshow_add(y1s_s1, y1s_s1, y1b_one);
  }
  poly_qshow_mul(tmp_poly, y1s_y1, y1->entries[IDX_B_SHOW]);
  poly_qshow_add(e0, e0, tmp_poly);
  poly_qshow_mul(tmp_poly, y1s_y1, s1->entries[IDX_B_SHOW]);
  poly_qshow_add(e1, e1, tmp_poly);
  poly_qshow_mul(tmp_poly, y1s_s1, y1->entries[IDX_B_SHOW]);
  poly_qshow_add(e1, e1, tmp_poly);

  /**********************************************
  * Computing linear terms in e1:
  *  sum_{i < l} mu_i.(-[Byg.y2_1]_i - 2^{-1}.(SE_{2i} + x^{n_show/2}.SE_{2i+1}).[Byg.y2_1*]_{:256/n_show}
  *                    -(SE_{2i}* + x^{n_show/2}.SE_{2i+1}*).[Byg.y2_1]_{:256/n_show}
  *                    + (SR_{2i} + x^{n_show/2}.SR_{2i+1}).y1* + (SR_{2i}* + x^{n_show/2}.SR_{2i+1}*).y1
  *                    + c_i[i].y_{1,b}* + c_i_prime[i].y_{1,b})
  **********************************************/
  poly_qshow_mat_256l_m2_d_mul_vec_m2_d(tmp_vec_256_l, Byg, y2_1); // B_{y,g}.y2_1
  for (i = 0; i < PARAM_L_SHOW; i++) {
    // using y1s_y1 to store sum before multiplying by mu_i
    poly_qshow_mul(y1s_y1, c_i_prime->entries[i], y1->entries[IDX_B_SHOW]);
    poly_qshow_conjugate(tmp_poly, y1->entries[IDX_B_SHOW]);
    poly_qshow_mul(tmp_poly, tmp_poly, c_i->entries[i]);
    poly_qshow_add(y1s_y1, y1s_y1, tmp_poly);

    poly_qshow_sub(y1s_y1, y1s_y1, tmp_vec_256_l->entries[PARAM_ARP_DIV_N_SHOW + i]);

    poly_qshow_zero(t0); // t0 used as temp variable
    // term in SE_i
    for (j = 0; j < PARAM_ARP_DIV_N_SHOW; j++) {
      poly_qshow_mul_xj(tmp_poly, sum_gamma_e_star[2*i+1]->entries[j], PARAM_N_SHOW/2); 
      poly_qshow_sub(y1s_s1, sum_gamma_e_star[2*i]->entries[j], tmp_poly); // SE_{2i}* - x^{n_show/2}.SE_{2i+1}* // y1s_s1 temp variable
      poly_qshow_add(tmp_poly, sum_gamma_e_star[2*i]->entries[j], tmp_poly); // SE_{2i}* + x^{n_show/2}.SE_{2i+1}*

      poly_qshow_mul(tmp_poly, tmp_poly, tmp_vec_256_l->entries[j]); // (SE_{2i}* + x^{n_show/2}.SE_{2i+1}*).[Byg.y2_1]_{:256/n_show}
      poly_qshow_sub(t0, t0, tmp_poly);
      poly_qshow_mul(y1s_s1, y1s_s1, tmp_vec_256_l->entries[j]); // (SE_{2i}* - x^{n_show/2}.SE_{2i+1}*).[Byg.y2_1]_{:256/n_show}
      poly_qshow_conjugate(tmp_poly, y1s_s1); // (SE_{2i} + x^{n_show/2}.SE_{2i+1}).[Byg.y2_1*]_{:256/n_show}
      poly_qshow_sub(t0, t0, tmp_poly);
    }
    // term in SR_i
    for (j = 0; j < PARAM_M1_SHOW; j++) {
      poly_qshow_mul_xj(tmp_poly, sum_gamma_r_star[2*i+1]->entries[j], PARAM_N_SHOW/2); 
      poly_qshow_sub(y1s_s1, sum_gamma_r_star[2*i]->entries[j], tmp_poly); // SR_{2i,1}* - x^{n_show/2}.SR_{2i+1,1}* // y1s_s1 temp variable
      poly_qshow_add(tmp_poly, sum_gamma_r_star[2*i]->entries[j], tmp_poly); // SR_{2i,1}* + x^{n_show/2}.SR_{2i+1,1}*

      poly_qshow_mul(tmp_poly, tmp_poly, y1->entries[j]); // (SR_{2i,1}* + x^{n_show/2}.SR_{2i+1,1}*).y1
      poly_qshow_add(t0, t0, tmp_poly);
      poly_qshow_mul(y1s_s1, y1s_s1, y1->entries[j]); // (SR_{2i,1}* - x^{n_show/2}.SR_{2i+1,1}*).y1
      poly_qshow_conjugate(tmp_poly, y1s_s1); // (SR_{2i,1} + x^{n_show/2}.SR_{2i+1,1}).y1*
      poly_qshow_add(t0, t0, tmp_poly);
    }
    poly_qshow_mul_scalar(t0, t0, PARAM_TWO_INVMOD_Q_SHOW);
    poly_qshow_add(y1s_y1, y1s_y1, t0);

    poly_qshow_mul(tmp_poly, y1s_y1, chal_3_l->entries[i]); // multiply by mu_i
    poly_qshow_add(e1, e1, tmp_poly); // add to e1
  }

  // committing to garbage terms
  poly_qshow_vec_m2_d_mul_inner(t0, b, y2_1);
  poly_qshow_add(t0, t0, e0);
  poly_qshow_vec_m2_d_mul_inner(tmp_poly, b, s2_1);
  poly_qshow_add(proof_2->t1, tmp_poly, e1);

  // computing fourth challenge
  kappa_c = 0;
  buf[0] = 4;
  poly_qshow_pack(buf + CHAL3_SHOW_INPUT_BYTES, t0);
  poly_qshow_pack(buf + CHAL3_SHOW_INPUT_BYTES + POLYQSHOW_PACKEDBYTES, proof_2->t1);
  shake256(challenge_seed, SEED_BYTES, buf, CHAL4_SHOW_INPUT_BYTES);
  do {
    poly_qshow_sample_challenge(proof_2->c, challenge_seed, DOMAIN_SEPARATOR_CHAL4_SHOW, kappa_c++, SEED_BYTES);
  } while (challenge_size_show(proof_2->c) > PARAM_ETA_SHOW);
  proof_2->ctr_c = kappa_c - 1;

  // clean up matrices, vectors and polynomials
  poly_qshow_clear(tmp_poly);
  poly_qshow_clear(e0);
  poly_qshow_clear(e1);
  poly_qshow_clear(y1i_star);
  poly_qshow_clear(y1s_y1);
  poly_qshow_clear(y1s_s1);
  poly_qshow_clear(y1b_one);
  poly_qshow_clear(t0);
  poly_qshow_vec_256_l_clear(tmp_vec_256_l);
  poly_qshow_vec_l_clear(c_i);
  poly_qshow_vec_l_clear(c_i_prime);
  poly_qshow_vec_k_clear(tmp_vec_k);
  poly_qshow_vec_k_clear(sum_vec_k_y);
  poly_qshow_vec_k_clear(sum_vec_k_s);
  poly_qshow_mat_k_k_clear(chal_3_quad_matrix);
  poly_qshow_vec_k_clear(Gy1_w2);
  poly_qshow_vec_k_clear(Gs1_w2);
  for (i = 0; i < 4; i++) {
    poly_qshow_clear(sum_mu_gamma[i]);
  }
}

/*************************************************
* Name:        prove_2_round5 [static]
*
* Description: Compute round 5 of zero-knowledge proof included in the blind signature
*
* Arguments:   - proof_2_t *proof_2: pointer to showuance proof structure
*              - const poly_qshow_vec_m1 s1: polynomial vectors, witness
*              - const poly_qshow_vec_m2_d s2_1: polynomial vector, ABDLOP commitment randomness
*              - const poly_qshow_vec_d s2_2: polynomial vector, ABDLOP commitment randomness
*              - const poly_qshow_vec_m1 y1: polynomial vectors, mask for c.s1
*              - const poly_qshow_vec_m2_d y2_1: polynomial vector, mask for c.s2_1
*              - const poly_qshow_vec_d y2_2: polynomial vector, mask for c.s2_2
*              - const poly_qshow_vec_d tA_L: polynomial vector, low bits of witness commitment w
*              - const poly_qshow_vec_d w_H: polynomial vector, high bits of mask commitment w
*              - const poly_qshow_vec_d w_L: polynomial vector, low bits of mask commitment w
*              - const uint64_t sq_norm_arp: precomputed square norm of Rx
*              - const int64_t z3_Rx_inner: precomputed inner product <z3, Rx>
* 
* Returns 1 if round 5 passes, and 0 if it rejects
**************************************************/
static int prove_2_round5(
    proof_2_t                  *proof_2,
    const poly_qshow_vec_m1    s1,
    const poly_qshow_vec_m2_d  s2_1,
    const poly_qshow_vec_d     s2_2,
    const poly_qshow_vec_m1    y1,
    const poly_qshow_vec_m2_d  y2_1,
    const poly_qshow_vec_d     y2_2,
    const poly_qshow_vec_d     tA_L,
    const poly_qshow_vec_d     w_H,
    const poly_qshow_vec_d     w_L,
    const uint64_t             sq_norm_arp,
    const int64_t              z3_Rx_inner) {
  size_t i,j;
  int pass = 1;
  poly_qshow tmp_poly;
  poly_qshow_vec_d z2_2, tmp_vec_d;
  uint128 sq_norm_cs2 = 0;
  uint64_t sq_norm_cs1 = 0;
  int128 z2_prime_cs2_inner = 0;
  int64_t z1_cs1_inner = 0;

  // init vectors and polynomials
  poly_qshow_init(tmp_poly);
  poly_qshow_vec_d_init(z2_2);
  poly_qshow_vec_d_init(tmp_vec_d);

  for (i = 0; i < PARAM_M1_SHOW; i++) {
    poly_qshow_mul(tmp_poly, s1->entries[i], proof_2->c);
    CHK_UI_OVF_ADDITION(sq_norm_cs1, (uint64_t)poly_qshow_sq_norm2(tmp_poly));
    poly_qshow_add(proof_2->z1->entries[i], y1->entries[i], tmp_poly);
    for (j = 0; j < PARAM_N_SHOW; j++) {
      z1_cs1_inner += (poly_qshow_get_coeff_centered(proof_2->z1->entries[i], j) * poly_qshow_get_coeff_centered(tmp_poly, j)); // should not overflow 63 bits
    }
  }

  for (i = 0; i < PARAM_M2_D_SHOW; i++) {
    poly_qshow_mul(tmp_poly, s2_1->entries[i], proof_2->c);
    sq_norm_cs2 += poly_qshow_sq_norm2(tmp_poly);
    poly_qshow_add(proof_2->z2_1->entries[i], y2_1->entries[i], tmp_poly);
    for (j = 0; j < PARAM_N_SHOW; j++) {
      z2_prime_cs2_inner += ((int128)poly_qshow_get_coeff_centered(proof_2->z2_1->entries[i], j) * (int128)poly_qshow_get_coeff_centered(tmp_poly, j));
    }
  }

  for (i = 0; i < PARAM_D_SHOW; i++) {
    poly_qshow_mul(tmp_poly, s2_2->entries[i], proof_2->c);
    sq_norm_cs2 += poly_qshow_sq_norm2(tmp_poly);
    poly_qshow_add(z2_2->entries[i], y2_2->entries[i], tmp_poly);
    for (j = 0; j < PARAM_N_SHOW; j++) {
      z2_prime_cs2_inner += ((int128)poly_qshow_get_coeff_centered(z2_2->entries[i], j) * (int128)poly_qshow_get_coeff_centered(tmp_poly, j));
    }
  }

  // rejection sampling
  // sample u1 uniform in (0,1), goto proof_2_reject if u1 > exp(pi * (sq_norm_y1 - sq_norm_z1) / PARAM_S1SQ_SHOW) / PARAM_REJ1_SHOW
  // sample u2 uniform in (0,1), goto proof_2_reject if u2 > exp(pi * (sq_norm_y2 - sq_norm_z2) / PARAM_S2SQ_SHOW) / PARAM_REJ2_SHOW
  if (_reject_exp(
    expl(M_PI * ((long double)sq_norm_cs1/(long double)PARAM_S1SQ_SHOW + (long double)sq_norm_arp/(long double)PARAM_S3SQ_SHOW)) / ((long double)(PARAM_REJ1_SHOW*PARAM_REJ3_SHOW) * coshl(2*M_PI*((long double)z1_cs1_inner/(long double)PARAM_S1SQ_SHOW + (long double)z3_Rx_inner/(long double)PARAM_S3SQ_SHOW)))
  ) || _reject_exp(
    expl(M_PI * (long double)sq_norm_cs2/(long double)PARAM_S2SQ_SHOW) / ((long double)PARAM_REJ2_SHOW * coshl(2*M_PI*(long double)z2_prime_cs2_inner/(long double)PARAM_S2SQ_SHOW))
  )) {
    pass = 0;
    goto prove_2_round5_cleanup;
  }

  // hint computation
  poly_qshow_vec_d_mul_poly_qshow(tmp_vec_d, tA_L, proof_2->c);
  poly_qshow_vec_d_sub(z2_2, z2_2, tmp_vec_d);
  poly_qshow_vec_d_sub(z2_2, z2_2, w_L);

  sq_norm_cs2 = poly_qshow_vec_m2_d_norm2(proof_2->z2_1) + poly_qshow_vec_d_norm2(z2_2);
  if (sq_norm_cs2 > ((uint128)PARAM_B2SQ_SHOW_LOW64 + (((uint128)PARAM_B2SQ_SHOW_HIGH64) << 64))) {
    pass = 0;
    goto prove_2_round5_cleanup;
  }
  poly_qshow_vec_d_mul_scalar(tmp_vec_d, w_H, PARAM_GAMMA_SHOW);
  poly_qshow_vec_d_sub(tmp_vec_d, tmp_vec_d, z2_2);
  poly_qshow_vec_d_makeGhint(proof_2->hint, z2_2, tmp_vec_d);
  
  // clean up vectors and polynomials
prove_2_round5_cleanup:
  poly_qshow_clear(tmp_poly);
  poly_qshow_vec_d_clear(z2_2);
  poly_qshow_vec_d_clear(tmp_vec_d);
  return pass;
}

/*************************************************
* Name:        prove_2
*
* Description: Compute zero-knowledge proof included in the blind signature
*
* Arguments:   - proof_2_t *proof_2: pointer to show proof structure
*              - const poly_qshow_mat_k_k *A_embed: array of polynomial matrices to host subring embedding of q1.b1.A'
*              - const poly_qshow_mat_k_k *B_embed: array of polynomial matrices to host subring embedding of q1.b2.B
*              - const poly_qshow_mat_k_k *A3_embed: array of polynomial matrices to host subring embedding of q1.b2.A3
*              - const poly_qshow_mat_k_k *G_embed: array of polynomial matrices to host subring embedding of q1.G.w_{2,L}
*              - const poly_qshow_vec_k *u_embed: array of polynomial vectors to host subring embedding of q1.(u + d.H(m) - [I|A'].w_{1,L} + B.w_{2,L} - A3.w_{3,L})
*              - const poly_qshow_vec_m1 s1: polynomial vector to host subring embedding of witness
*              - const uint8_t *crs_seed: pointer to byte array containing the CRS seed (allocated SEED_BYTES bytes)
*              - const uint8_t *seed: pointer to byte array containing the seed for public parameters (allocated SEED_BYTES bytes)
**************************************************/
void prove_2(
    proof_2_t                *proof_2, 
    const poly_qshow_mat_k_k A_embed[PARAM_D][PARAM_D], 
    const poly_qshow_mat_k_k B_embed[PARAM_D][PARAM_K*PARAM_D], 
    const poly_qshow_mat_k_k A3_embed[PARAM_D][PARAM_K], 
    const poly_qshow_mat_k_k G_embed[PARAM_D], 
    const poly_qshow_vec_k   u_embed[PARAM_D], 
    const poly_qshow_vec_m1  s1, 
    const uint8_t            crs_seed[CRS_SEED_BYTES], 
    const uint8_t            seed[SEED_BYTES]) {
  size_t i;
  uint8_t randomness_seed[SEED_BYTES];
  uint8_t buf[CHAL4_SHOW_INPUT_BYTES] = {0};
  uint32_t kappa;
  uint64_t sq_norm_arp;
  int64_t z3_Rx_inner;
  coeff_qshow tmp_coeff;
  poly_qshow tmp_poly, tmp_poly_2, one, quadratic_precomp[4];
  poly_qshow_vec_m1 y1, s1_star, sum_gamma_r_star[2*PARAM_L_SHOW];
  poly_qshow_vec_256 sum_gamma_e_star[2*PARAM_L_SHOW];
  poly_qshow_vec_m2_d s2_1, y2_1, b;
  poly_qshow_vec_d s2_2, y2_2, tA_L, w_H, w_L;
  poly_qshow_vec_256_l y3_g;
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
  poly_qshow_init(tmp_poly_2);
  poly_qshow_init(one);
  for (i = 0; i < 4; i++) {
    poly_qshow_init(quadratic_precomp[i]);
    poly_qshow_zero(quadratic_precomp[i]);
  }
  // init vectors and matrices
  poly_qshow_vec_m1_init(y1);
  poly_qshow_vec_m1_init(s1_star);
  for (i = 0; i < 2*PARAM_L_SHOW; i++) {
    poly_qshow_vec_m1_init(sum_gamma_r_star[i]);    
    poly_qshow_vec_256_init(sum_gamma_e_star[i]);    
  }
  poly_qshow_vec_m2_d_init(s2_1);
  poly_qshow_vec_m2_d_init(y2_1);
  poly_qshow_vec_m2_d_init(b);
  poly_qshow_vec_d_init(s2_2);
  poly_qshow_vec_d_init(y2_2);
  poly_qshow_vec_d_init(tA_L);
  poly_qshow_vec_d_init(w_H);
  poly_qshow_vec_d_init(w_L);
  poly_qshow_vec_256_l_init(y3_g);
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

  // generate random secret seed
  randombytes(randomness_seed, SEED_BYTES);

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

  // precomputations
  // s1*
  poly_qshow_vec_m1_conjugate(s1_star, s1);
  // <s_{1,1}'*,s_{1,1}'> - B_1'^2
  for (i = IDX_W1_SHOW; i < IDX_W23_SHOW; i++) {
    poly_qshow_mul(tmp_poly, s1_star->entries[i], s1->entries[i]);
    poly_qshow_add(quadratic_precomp[0], quadratic_precomp[0], tmp_poly);
  }
  poly_qshow_muladd_constant(quadratic_precomp[0], -PARAM_B1_PRIME_SQ, 1); 
  // <s_{1,23}'*,s_{1,23}'> - B_2'^2
  for (i = IDX_W23_SHOW; i < IDX_TAG_SHOW; i++) {
    poly_qshow_mul(tmp_poly, s1_star->entries[i], s1->entries[i]);
    poly_qshow_add(quadratic_precomp[1], quadratic_precomp[1], tmp_poly);
  }
  poly_qshow_muladd_constant(quadratic_precomp[1], -PARAM_B2_PRIME_SQ, 1); 
  // <s_{1,t}*,s_{1,t}> - w
  for (i = IDX_TAG_SHOW; i < IDX_B_SHOW; i++) {
    poly_qshow_mul(tmp_poly, s1_star->entries[i], s1->entries[i]);
    poly_qshow_add(quadratic_precomp[2], quadratic_precomp[2], tmp_poly);
  }
  poly_qshow_muladd_constant(quadratic_precomp[2], -PARAM_W, 1);
  // <s_{1,t}*, s_{1,t} - s_{1,b}.one>
  tmp_coeff = poly_qshow_get_coeff_centered(s1->entries[IDX_B_SHOW], 0); // s_{1,b}
  for (i = 0; i < PARAM_N_SHOW; i++) {
    poly_qshow_set_coeff(one, i, 1);
    poly_qshow_set_coeff(tmp_poly_2, i, tmp_coeff); // tmp_poly_2 contains s_{1,b}.one
  }
  for (i = IDX_TAG_SHOW; i < IDX_B_SHOW; i++) {
    poly_qshow_sub(tmp_poly, s1->entries[i], tmp_poly_2);
    poly_qshow_mul(tmp_poly, s1_star->entries[i], tmp_poly);
    poly_qshow_add(quadratic_precomp[3], quadratic_precomp[3], tmp_poly);
  }

  kappa = 0;
proof_2_reject:
  /****** first round ******/
  prove_2_round1(proof_2, chal_1, buf, s2_1, s2_2, y1, y2_1, y2_2, y3_g, tA_L, w_H, w_L, s1, A1, A2, Byg, randomness_seed, kappa);
  kappa += PARAM_ARP_DIV_N_L_SHOW - PARAM_ARP_DIV_N_SHOW + 1;

  /****** second round ******/
  sq_norm_arp = prove_2_round2(proof_2, chal_2, &z3_Rx_inner, buf, s1, chal_1, y3_g);

  /****** third round ******/
  prove_2_round3(proof_2, chal_3_l, chal_3_dk, chal_3_1, buf, sum_gamma_e_star, sum_gamma_r_star, s1, chal_1, chal_2, y3_g, quadratic_precomp);
  
  /****** fourth round ******/
  prove_2_round4(proof_2, buf, s1, s1_star, s2_1, Byg, b, A_embed, B_embed, A3_embed, G_embed, sum_gamma_e_star, sum_gamma_r_star, y1, y2_1, chal_2, chal_3_l, chal_3_dk, chal_3_1, tmp_poly_2, one);

  /****** fifth round ******/
  if (!prove_2_round5(proof_2, s1, s2_1, s2_2, y1, y2_1, y2_2, tA_L, w_H, w_L, sq_norm_arp, z3_Rx_inner)) {
    goto proof_2_reject;
  }

  // clean up
  // clean up polynomials
  poly_qshow_clear(tmp_poly);
  poly_qshow_clear(tmp_poly_2);
  poly_qshow_clear(one);
  for (i = 0; i < 4; i++) {
    poly_qshow_clear(quadratic_precomp[i]);
  }
  // clean up vectors and matrices
  poly_qshow_vec_m1_clear(y1);
  poly_qshow_vec_m1_clear(s1_star);
  for (i = 0; i < 2*PARAM_L_SHOW; i++) {
    poly_qshow_vec_m1_clear(sum_gamma_r_star[i]);    
    poly_qshow_vec_256_clear(sum_gamma_e_star[i]);    
  }
  poly_qshow_vec_m2_d_clear(s2_1);
  poly_qshow_vec_m2_d_clear(y2_1);
  poly_qshow_vec_m2_d_clear(b);
  poly_qshow_vec_d_clear(s2_2);
  poly_qshow_vec_d_clear(y2_2);
  poly_qshow_vec_d_clear(tA_L);
  poly_qshow_vec_d_clear(w_H);
  poly_qshow_vec_d_clear(w_L);
  poly_qshow_vec_256_l_clear(y3_g);
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
}