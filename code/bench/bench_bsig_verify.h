#ifndef BENCH_BSIG_VERIFY_H
#define BENCH_BSIG_VERIFY_H

#include "benchmark.h"

double bsig_verify_valid_bench(timer* t);
double bsig_verify_invalid_bench(timer* t);

#endif /* BENCH_BSIG_VERIFY_H */

