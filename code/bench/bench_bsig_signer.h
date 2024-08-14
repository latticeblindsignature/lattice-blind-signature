#ifndef BENCH_BSIG_SIGNER_H
#define BENCH_BSIG_SIGNER_H

#include "benchmark.h"

double tag_gen_bench(timer* t);
double keygen_bench(timer* t);
double pre_sign_commitment_bench(timer* t);
double pre_sig_verify_from_commitment_valid_bench(timer* t);
double pre_sig_verify_from_commitment_invalid_bench(timer* t);

#endif /* BENCH_BSIG_SIGNER_H */

