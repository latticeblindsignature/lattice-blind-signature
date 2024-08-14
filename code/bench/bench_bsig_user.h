#ifndef BENCH_BSIG_USER_H
#define BENCH_BSIG_USER_H

#include "benchmark.h"

double tag_verify_bench(timer* t);
double commit_bench(timer* t);
double encrypt_bench(timer* t);
double embed_1_bench(timer *t);
double embed_1_verifier_bench(timer *t);
double complete_decompose_bench(timer* t);
double embed_2_bench(timer *t);
double embed_2_verifier_bench(timer *t);

/*double osig_user_embed_bench(timer* t);
double prove_1_bench(timer* t);
double prove_1_valid_bench(timer* t);
double prove_1_invalid_bench(timer* t);*/

#endif /* BENCH_BSIG_USER_H */

