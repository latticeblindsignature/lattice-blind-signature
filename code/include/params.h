#ifndef PARAMS_H
#define PARAMS_H

/*************************************************
* Domain separators for XOF expansion
**************************************************/
// sep
#define DOMAIN_SEPARATOR_A 0
#define DOMAIN_SEPARATOR_R 1
#define DOMAIN_SEPARATOR_A3 2
#define DOMAIN_SEPARATOR_U 3
#define DOMAIN_SEPARATOR_D 4
// cmt
#define DOMAIN_SEPARATOR_R1 5
#define DOMAIN_SEPARATOR_R2 6
#define DOMAIN_SEPARATOR_R3 7
// pke
#define DOMAIN_SEPARATOR_AE 8
#define DOMAIN_SEPARATOR_BE 9
#define DOMAIN_SEPARATOR_RE 10
// proof_1
#define DOMAIN_SEPARATOR_A1_ISS 11
#define DOMAIN_SEPARATOR_A2_ISS 12
#define DOMAIN_SEPARATOR_BYG_ISS 13
#define DOMAIN_SEPARATOR_B_ISS 14
#define DOMAIN_SEPARATOR_CHAL1_ISS 15
#define DOMAIN_SEPARATOR_CHAL2_ISS 16
#define DOMAIN_SEPARATOR_CHAL3_ISS 17
#define DOMAIN_SEPARATOR_CHAL4_ISS 18
#define DOMAIN_SEPARATOR_RAND_S2_ISS 19
#define DOMAIN_SEPARATOR_RAND_G_ISS 20
// proof_2
#define DOMAIN_SEPARATOR_A1_SHOW 21
#define DOMAIN_SEPARATOR_A2_SHOW 22
#define DOMAIN_SEPARATOR_BYG_SHOW 23
#define DOMAIN_SEPARATOR_B_SHOW 24
#define DOMAIN_SEPARATOR_CHAL1_SHOW 25
#define DOMAIN_SEPARATOR_CHAL2_SHOW 26
#define DOMAIN_SEPARATOR_CHAL3_SHOW 27
#define DOMAIN_SEPARATOR_CHAL4_SHOW 28
#define DOMAIN_SEPARATOR_RAND_S2_SHOW 29
#define DOMAIN_SEPARATOR_RAND_G_SHOW 30

/*************************************************
* Signature parameters
**************************************************/
// Ring degree for the signature
#define PARAM_N 256
// Modulus for the signature
#define PARAM_Q 8388581L
// Modulus bit-length upper bound for uniform sampling
#define PARAM_Q_BITLEN 23
// Module rank for the signature
#define PARAM_D 5
// Gadget dimension
#define PARAM_K 3
// Gadget base
#define PARAM_B 204
// First decomposition base
#define PARAM_B1 512
// First decomposition base (log2)
#define PARAM_B1_LOG 9
// Second decomposition base
#define PARAM_B2 8
// Second decomposition base (log2)
#define PARAM_B2_LOG 3
// Number of iterations for the spectral norm estimation
#define PARAM_IT_SPEC_NORM 5
// Hamming weight of the tags
#define PARAM_W 5
// Bound on the square spectral norm of R
#define PARAM_R_MAX_SQ_SPECTRAL_NORM 6888.15728882705207070103
// Gaussian parameter s_2 for v_2 and v_3
#define PARAM_S2 1156.13467820429173116281
// Squared Gaussian parameter s_1^2
#define PARAM_S1SQ 12436790137.75049209594726562500
// Gaussian width for p_2 (sqrt(s_2^2 - s_G^2))
#define PARAM_SQRT_S2SQ_SGSQ 922.41323638654216665600
// Negated ratio -1/(1/s_G^2 - 1/s_2^2)
#define PARAM_SGINVSQ_S2INVSQ -763175.46583292330615222454
// Negated ratio -s_G^2/(s_2^2 - s_G^2)
#define PARAM_NEGSGSQ_DIV_S2SQ_SGSQ -0.57096244617318558934
// Squared verification bound on v_1
#define PARAM_B1SQ 7222652870284UL
// Squared verification bound on [v_2 | v_3]
#define PARAM_B2SQ 1281829227UL
// Squared proof bound for w_{1,H}
#define PARAM_B1_PRIME_SQ 29168766L
// Squared proof bound for [w_{2,H} | w_{3,H}]
#define PARAM_B2_PRIME_SQ 21262196L

// Length of the public and secret seeds
#define SEED_BYTES 32
// Length of the public seed for CRS expansion
#define CRS_SEED_BYTES 32
// Length of the state
#define STATE_BYTES 64

/*************************************************
* Commitment parameters
**************************************************/
// Squared verification bound on r_1
#define PARAM_B_R1_SQ 2684354560L
// Squared verification bound on [r_2|r_3]
#define PARAM_B_R2_SQ 294912L

/*************************************************
* Encryption to the sky parameters
**************************************************/
// Modulus for the encryption to the sky
#define PARAM_P 4993L
// Modulus bit-length upper bound for uniform sampling
#define PARAM_P_BITLEN 13
// Round(p/2)
#define PARAM_HALF_P 2497L
// Module rank for the encryption to the sky
#define PARAM_DE 3
// Number of samples for the encryption to the sky
#define PARAM_ME 7
// Squared binomial tail bound 
#define PARAM_B_RE_SQ 1185L

/*************************************************
* [PROOF_1] Zero-Knowledge proof parameters
**************************************************/
// Ring degree for the issuance proof
#define PARAM_N_ISS 64
// Ring degree gap between the issuance proof and the signature (subring embedding)
#define PARAM_K_ISS 4
// Modulus for the issuance proof
#define PARAM_Q_ISS 144114722315180017L
// Modulus bit-length upper bound for uniform sampling
#define PARAM_Q_ISS_BITLEN 57
// Modulus factor for the issuance proof
#define PARAM_Q1_ISS 17179868957L
// (Negated) Inverse of p modulo q_1q_2
#define PARAM_P_INVMOD_Q_ISS_NEG 34838067261049926L
// Inverse of 2 modulo q_1q_2
#define PARAM_TWO_INVMOD_Q_ISS 72057361157590009L
// Module rank for the issuance proof
#define PARAM_D_ISS 22
// Witness dimension
#define PARAM_M1_ISS 148
// ABDLOP commitment randomness dimension m2 - d (s_{2,1})
#define PARAM_M2_D_ISS 47
// Soundness amplification dimension
#define PARAM_L_ISS 3
// Dimension for Approximate Range Proof
#define PARAM_ARP_ISS 256
// Rank for Approximate Range Proof (256 / n)
#define PARAM_ARP_DIV_N_ISS 4
// 256 / n + l
#define PARAM_ARP_DIV_N_L_ISS 7
// Gaussian mask width for cs_1
#define PARAM_S1_ISS 10258622.89452818408608436584
// Squared Gaussian mask width for cs_1
#define PARAM_S1SQ_ISS 105239343692137.81250000000000000000
// Gaussian mask width for cs_2
#define PARAM_S2_ISS 13157.08519160629293764941
// Squared Gaussian mask width for cs_2
#define PARAM_S2SQ_ISS 173108890.73918560147285461426
// Gaussian mask width for Rs_1 (ARP)
#define PARAM_S3_ISS 2220943.58923175372183322906
// Squared Gaussian mask width for Rs_1 (ARP)
#define PARAM_S3SQ_ISS 4932590426549.62500000000000000000
// Rejection sampling rate for y_1
#define PARAM_REJ1_ISS 2
// Rejection sampling rate for y_2
#define PARAM_REJ2_ISS 2
// Rejection sampling rate for y_3
#define PARAM_REJ3_ISS 2
// Squared verification bound for z_1
#define PARAM_B1SQ_ISS 767107606508928768UL
// Squared verification bound for z_2 (high bits)
#define PARAM_B2SQ_ISS_LOW64 13530852096220594176UL
// Squared verification bound for z_2 (low bits)
#define PARAM_B2SQ_ISS_HIGH64 48UL
// Squared verification bound for z_3
#define PARAM_B3SQ_ISS 543467640738730UL
// Infinity norm of challenges
#define PARAM_RHO_ISS 8
// Manhattan-like norm of challenges
#define PARAM_ETA_ISS 93
// Compression parameter for t_A (D)
#define PARAM_D_ROUND_ISS 21
// Compression parameter for w (gamma)
#define PARAM_GAMMA_ISS 603990638L
// Highbits bound
#define PARAM_Q_GAMMA_ISS 238604232L

/*************************************************
* [PROOF_2] Zero-Knowledge proof parameters
**************************************************/
// Ring degree for the show proof
#define PARAM_N_SHOW 64
// Ring degree gap between the show proof and the signature (subring embedding)
#define PARAM_K_SHOW 4
// Modulus for the show proof
#define PARAM_Q_SHOW 2251790057742217L
// Modulus bit-length upper bound for uniform sampling
#define PARAM_Q_SHOW_BITLEN 51
// Modulus factor for the show proof
#define PARAM_Q1_SHOW 268435157L
// (Negated) Inverse of p modulo q_1q_2
#define PARAM_P_INVMOD_Q_SHOW_NEG 887547132713135L
// Inverse of 2 modulo q_1q_2
#define PARAM_TWO_INVMOD_Q_SHOW 1125895028871109L
// Module rank for the show proof
#define PARAM_D_SHOW 22
// Witness dimension
#define PARAM_M1_SHOW 119
// ABDLOP commitment randomness dimension m2 - d (s_{2,1})
#define PARAM_M2_D_SHOW 43
// Soundness amplification dimension
#define PARAM_L_SHOW 3
// Dimension for Approximate Range Proof
#define PARAM_ARP_SHOW 256
// Rank for Approximate Range Proof (256 / n)
#define PARAM_ARP_DIV_N_SHOW 4
// 256 / n + l
#define PARAM_ARP_DIV_N_L_SHOW 7
// Gaussian mask width for cs_1
#define PARAM_S1_SHOW 1406027.47268197103403508663
// Squared Gaussian mask width for cs_1
#define PARAM_S1SQ_SHOW 1976913253936.45068359375000000000
// Gaussian mask width for cs_2
#define PARAM_S2_SHOW 12770.02712312388757709414
// Squared Gaussian mask width for cs_2
#define PARAM_S2SQ_SHOW 163073592.72531974315643310547
// Gaussian mask width for Rs_1 (ARP)
#define PARAM_S3_SHOW 277540.14345016190782189369
// Squared Gaussian mask width for Rs_1 (ARP)
#define PARAM_S3SQ_SHOW 77028531226.33645629882812500000
// Rejection sampling rate for y_1
#define PARAM_REJ1_SHOW 2
// Rejection sampling rate for y_2
#define PARAM_REJ2_SHOW 2
// Rejection sampling rate for y_3
#define PARAM_REJ3_SHOW 2
// Squared verification bound for z_1
#define PARAM_B1SQ_SHOW 11833315492058180UL
// Squared verification bound for z_2 (high bits)
#define PARAM_B2SQ_SHOW_LOW64 16830841499513323520UL
// Squared verification bound for z_2 (low bits)
#define PARAM_B2SQ_SHOW_HIGH64 2UL
// Squared verification bound for z_3
#define PARAM_B3SQ_SHOW 8486922796148UL
// Infinity norm of challenges
#define PARAM_RHO_SHOW 8
// Manhattan-like norm of challenges
#define PARAM_ETA_SHOW 93
// Compression parameter for t_A (D)
#define PARAM_D_ROUND_SHOW 19
// Compression parameter for w (gamma)
#define PARAM_GAMMA_SHOW 146557902L
// Highbits bound
#define PARAM_Q_GAMMA_SHOW 15364508L

#endif /* PARAMS_H */
