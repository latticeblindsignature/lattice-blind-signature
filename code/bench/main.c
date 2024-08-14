#include <stdio.h>
#include <time.h>

#include "bsig_user.h"
#include "bsig_signer.h"
#include "randombytes.h"
#include "random.h"

#include "bench_bsig_signer.h"
#include "bench_bsig_user.h"
#include "bench_bsig_verify.h"
#include "bench_proof_1.h"
#include "bench_proof_2.h"

#define BENCH_ITERATIONS 10

int main(void) {
  arith_setup();
  random_init();

  printf("[+] Running ./build/bench with N = %d iterations\n\n", BENCH_ITERATIONS);
  printf("_____________________________________________________________________________________________________________________________________________________________________________________\n");

  printf("%-50s\t%26s%4s%27s\t%27s%6s%28s\n", "", "", "Time", "", "", "Cycles", "");
  printf("%-50s\t%57s\t%61s\n", "", "---------------------------------------------------------", "-------------------------------------------------------------");
  printf("%-50s\t%9s\t%9s\t%9s\t%9s\t%13s\t%13s\t%13s\t%13s\n", "Benchmarked Functionality", "mean (ms)", "med (ms)", "min (ms)", "max (ms)", "mean (cycles)", "med (cycles)", "min (cycles)", "max (cycles)");
  printf("_____________________________________________________________________________________________________________________________________________________________________________________\n");

  printf("\n");

  benchmark("[SIGNER] KEYGEN", BENCH_ITERATIONS, keygen_bench);
  benchmark("[SIGNER] TAG_GEN", BENCH_ITERATIONS, tag_gen_bench);
  benchmark("[ USER ] TAG_VERIFY", BENCH_ITERATIONS, tag_verify_bench);
  benchmark("[ USER ] COMMIT", BENCH_ITERATIONS, commit_bench);
  benchmark("[ USER ] ENCRYPT", BENCH_ITERATIONS, encrypt_bench);
  benchmark("[ USER ] EMBED_1 (user side)", BENCH_ITERATIONS, embed_1_bench);
  benchmark("[SIGNER] EMBED_1 (verifier side)", BENCH_ITERATIONS, embed_1_verifier_bench);
  benchmark("[ USER ] PROVE_1", BENCH_ITERATIONS, prove_1_bench);
  benchmark("[ USER ] VERIFY_1 (valid)", BENCH_ITERATIONS, verify_1_valid_bench);
  benchmark("[ USER ] VERIFY_1 (invalid c)", BENCH_ITERATIONS, verify_1_invalid_c_bench);
  benchmark("[ USER ] VERIFY_1 (invalid f)", BENCH_ITERATIONS, verify_1_invalid_f_bench);
  benchmark("[SIGNER] PRE_SIGN_COMMITMENT", BENCH_ITERATIONS, pre_sign_commitment_bench);
  benchmark("[ USER ] PRE_SIG_VERIFY_FROM_COMMITMENT (valid)", BENCH_ITERATIONS, pre_sig_verify_from_commitment_valid_bench);
  benchmark("[ USER ] PRE_SIG_VERIFY_FROM_COMMITMENT (invalid)", BENCH_ITERATIONS, pre_sig_verify_from_commitment_invalid_bench);
  benchmark("[ USER ] COMPLETE_DECOMPOSE", BENCH_ITERATIONS, complete_decompose_bench);
  benchmark("[ USER ] EMBED_2 (user side)", BENCH_ITERATIONS, embed_2_bench);
  benchmark("[ USER ] PROVE_2", BENCH_ITERATIONS, prove_2_bench);
  benchmark("[ USER ] VERIFY_2 (valid)", BENCH_ITERATIONS, verify_2_valid_bench);
  benchmark("[ USER ] VERIFY_2 (invalid c)", BENCH_ITERATIONS, verify_2_invalid_c_bench);
  benchmark("[ USER ] VERIFY_2 (invalid f)", BENCH_ITERATIONS, verify_2_invalid_f_bench);
  benchmark("[VERIF ] EMBED_2 (verifier side)", BENCH_ITERATIONS, embed_2_verifier_bench);
  benchmark("[VERIF ] BSIG_VERIFY (valid)(includes embed_2)", BENCH_ITERATIONS, bsig_verify_valid_bench);
  benchmark("[VERIF ] BSIG_VERIFY (invalid)(includes embed_2)", BENCH_ITERATIONS, bsig_verify_invalid_bench);

  printf("_____________________________________________________________________________________________________________________________________________________________________________________\n");


  arith_teardown();
  return 0;
}


