#ifndef BSIG_VERIFY_H
#define BSIG_VERIFY_H

#include <stdint.h>
#include "params.h"
#include "bsig_signer.h"
#include "bsig_user.h"
#include "arith.h"

int bsig_verify(const bsig_t *bsig, const pk_t *pk, const uint8_t msg[PARAM_N/8]);

#endif /* BSIG_VERIFY_H */