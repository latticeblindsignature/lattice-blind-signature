#ifndef BSIG_USER_H
#define BSIG_USER_H

#include <stdint.h>
#include "params.h"
#include "bsig_signer.h"
#include "arith.h"

// Starting index for r1 in s1
#define IDX_R1_ISS 0
// Starting index for r12 in s1
#define IDX_R12_ISS (IDX_R1_ISS + PARAM_D*PARAM_K_ISS)
// Starting index for r23 in s1
#define IDX_R23_ISS (IDX_R1_ISS + 2*PARAM_D*PARAM_K_ISS + 1) // +1 for four squares polynomial
// Starting index for r3 in s1
#define IDX_R3_ISS (IDX_R23_ISS + PARAM_D*PARAM_K*PARAM_K_ISS) 
// Starting index for re in s1
#define IDX_RE_ISS (IDX_R23_ISS + (PARAM_D+1)*PARAM_K*PARAM_K_ISS + 1) // +1 for four squares polynomial
// Starting index for H(m) in s1 
#define IDX_MSG_ISS (IDX_RE_ISS + PARAM_ME*PARAM_K_ISS + 1) // +1 for four squares polynomial
// Starting index for b in s1 
#define IDX_B_ISS (IDX_MSG_ISS + PARAM_K_ISS) 

// Byte-length of language packing 
#define ISS_CHALLENGE_BASE_BYTES (1 + CRS_SEED_BYTES + SEED_BYTES + (PARAM_D+PARAM_DE+1)*POLYQISS_VECK_PACKEDBYTES)
// Byte-length of first challenge input
#define CHAL1_ISS_INPUT_BYTES (ISS_CHALLENGE_BASE_BYTES + 2*POLYQISS_VECD_PACKEDBYTES + POLYQISS_VEC256L_PACKEDBYTES)
// Byte-length of second challenge input
#define CHAL2_ISS_INPUT_BYTES (CHAL1_ISS_INPUT_BYTES + PARAM_ARP_ISS*COEFFQISS_PACKEDBYTES)
// Byte-length of third challenge input 
#define CHAL3_ISS_INPUT_BYTES (CHAL2_ISS_INPUT_BYTES + POLYQISS_VECL_PACKEDBYTES)
// Byte-length of fourth challenge input 
#define CHAL4_ISS_INPUT_BYTES (CHAL3_ISS_INPUT_BYTES + 2*POLYQISS_PACKEDBYTES)

// Starting index for w1H in s1
#define IDX_W1_SHOW 0
// Starting index for w1H2 in s1
#define IDX_W12_SHOW (IDX_W1_SHOW + PARAM_D*PARAM_K_SHOW)
// Starting index for w2H in s1
#define IDX_W23_SHOW (IDX_W12_SHOW + PARAM_D*PARAM_K_SHOW + 1) // +1 for four squares polynomial
// Starting index for w3H in s1 
#define IDX_W3_SHOW (IDX_W23_SHOW + PARAM_D*PARAM_K*PARAM_K_SHOW) 
// Starting index for tag in s1 
#define IDX_TAG_SHOW (IDX_W23_SHOW + (PARAM_D+1)*PARAM_K*PARAM_K_SHOW + 1) // +1 for four squares polynomial
// Starting index for usk in s1 
#define IDX_B_SHOW (IDX_TAG_SHOW + PARAM_K_SHOW)

// Byte-length of language packing 
#define SHOW_CHALLENGE_BASE_BYTES (1 + CRS_SEED_BYTES + SEED_BYTES + PARAM_D*POLYQSHOW_VECK_PACKEDBYTES) 
// Byte-length of first challenge input
#define CHAL1_SHOW_INPUT_BYTES (SHOW_CHALLENGE_BASE_BYTES + 2*POLYQSHOW_VECD_PACKEDBYTES + POLYQSHOW_VEC256L_PACKEDBYTES)
// Byte-length of second challenge input 
#define CHAL2_SHOW_INPUT_BYTES (CHAL1_SHOW_INPUT_BYTES + PARAM_ARP_SHOW*COEFFQSHOW_PACKEDBYTES)
// Byte-length of third challenge input
#define CHAL3_SHOW_INPUT_BYTES (CHAL2_SHOW_INPUT_BYTES + POLYQSHOW_VECL_PACKEDBYTES)
// Byte-length of fourth challenge input
#define CHAL4_SHOW_INPUT_BYTES (CHAL3_SHOW_INPUT_BYTES + 2*POLYQSHOW_PACKEDBYTES)

typedef struct {
  poly_q_vec_d r1[2];
  poly_q_vec_d r2[PARAM_K];
  poly_q_vec_k r3;
} rand_t;

typedef struct {
  poly_p_vec_d ct0;
  poly_p ct1;
} ct_t;

typedef struct {
  poly_qiss_vec_d tA_H;
  poly_qiss_vec_256_l tB;
  coeff_qiss z3[PARAM_ARP_ISS];
  poly_qiss_vec_l f;
  poly_qiss t1;
  poly_qiss c;
  uint32_t ctr_c;
  poly_qiss_vec_m1 z1;
  poly_qiss_vec_m2_d z2_1;
  poly_qiss_vec_d hint;
} proof_1_t;

typedef struct {
  poly_qshow_vec_d tA_H;
  poly_qshow_vec_256_l tB;
  coeff_qshow z3[PARAM_ARP_SHOW];
  poly_qshow_vec_l f;
  poly_qshow t1;
  poly_qshow c;
  uint32_t ctr_c;
  poly_qshow_vec_m1 z1;
  poly_qshow_vec_m2_d z2_1;
  poly_qshow_vec_d hint;
} proof_2_t;

typedef struct {
  uint8_t crs_seed[CRS_SEED_BYTES];
  poly_q_vec_d w1L[2];
  poly_q_vec_d w2L[PARAM_K];
  poly_q_vec_k w3L;
  proof_2_t proof_2;
} bsig_t;

void rand_init(rand_t *rand);
void rand_clear(rand_t *rand);
void ct_init(ct_t *ct);
void ct_clear(ct_t *ct);
void proof_1_init(proof_1_t *proof_1);
void proof_1_clear(proof_1_t *proof_1);
void proof_2_init(proof_2_t *proof_2);
void proof_2_clear(proof_2_t *proof_2);
void bsig_init(bsig_t *bsig);
void bsig_clear(bsig_t *bsig);

void commit(rand_t *rand, poly_q_vec_d cmt, const poly_q tag, const uint8_t msg[PARAM_N/8], const pk_t *pk);
void encrypt(ct_t *ct, poly_p_vec_m r_e, const uint8_t msg[PARAM_N/8], const pk_t *pk);

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
    const uint8_t msg[PARAM_N/8]);
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
    const poly_q tag);
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
    const uint8_t           seed[SEED_BYTES]);
int verify_1(
    const proof_1_t         *proof_1, 
    const poly_qiss_mat_k_k A_embed[PARAM_D][PARAM_D], 
    const poly_qiss_mat_k_k B_embed[PARAM_D][PARAM_K*(PARAM_D+1)], 
    const poly_qiss_mat_k_k D_embed[PARAM_D], 
    const poly_qiss_mat_k_k Ae_be_embed[PARAM_DE+1][PARAM_ME],
    const poly_qiss_vec_k   cmt_embed[PARAM_D], 
    const poly_qiss_vec_k   ct_embed[PARAM_DE+1],  
    const uint8_t           crs_seed[CRS_SEED_BYTES], 
    const uint8_t           seed[SEED_BYTES]);

int pre_sig_verify(poly_q_vec_d v11, const poly_q tag, const pre_sig_t *pre_sig, const poly_q_vec_d cmt, const pk_t *pk);
void complete_decompose(
    bsig_t *bsig,
    poly_q_vec_d w1H[2],
    poly_q_vec_d w2H[PARAM_K],
    poly_q_vec_k w3H,
    const poly_q_vec_d v11,
    const pre_sig_t *pre_sig, 
    const rand_t *rand);

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
    const uint8_t msg[PARAM_N/8]);
void embed_2_verifier(
    poly_qshow_mat_k_k A_embed[PARAM_D][PARAM_D], 
    poly_qshow_mat_k_k B_embed[PARAM_D][PARAM_K*PARAM_D], 
    poly_qshow_mat_k_k A3_embed[PARAM_D][PARAM_K], 
    poly_qshow_mat_k_k G_embed[PARAM_D], 
    poly_qshow_vec_k u_embed[PARAM_D], 
    const pk_t *pk,
    const bsig_t *bsig, 
    const uint8_t msg[PARAM_N/8]);
void prove_2(
    proof_2_t                *proof_2, 
    const poly_qshow_mat_k_k A_embed[PARAM_D][PARAM_D], 
    const poly_qshow_mat_k_k B_embed[PARAM_D][PARAM_K*PARAM_D], 
    const poly_qshow_mat_k_k A3_embed[PARAM_D][PARAM_K], 
    const poly_qshow_mat_k_k G_embed[PARAM_D], 
    const poly_qshow_vec_k   u_embed[PARAM_D], 
    const poly_qshow_vec_m1  s1, 
    const uint8_t            crs_seed[CRS_SEED_BYTES], 
    const uint8_t            seed[SEED_BYTES]);
int verify_2(
    const proof_2_t          *proof_2, 
    const poly_qshow_mat_k_k A_embed[PARAM_D][PARAM_D], 
    const poly_qshow_mat_k_k B_embed[PARAM_D][PARAM_K*PARAM_D], 
    const poly_qshow_mat_k_k A3_embed[PARAM_D][PARAM_K], 
    const poly_qshow_mat_k_k G_embed[PARAM_D], 
    const poly_qshow_vec_k   u_embed[PARAM_D], 
    const uint8_t            crs_seed[CRS_SEED_BYTES], 
    const uint8_t            seed[SEED_BYTES]);

#endif /* BSIG_USER_H */
