/**
 *******************************************************************************
 * @version 0.1
 * @date 2021-12-16
 * @copyright Copyright © 2021 by University of Luxembourg
 * @author Hao Cheng
 *******************************************************************************
 */

#ifndef _SIKE_H
#define _SIKE_H

#include "sidh.h"
#include "fips202.h"

void crypto_kem_keypair(uint8_t *pk, uint8_t *sk);
void crypto_kem_enc(uint8_t *ct, uint8_t *ss, const uint8_t *pk);
void crypto_kem_dec(uint8_t *ss, const uint8_t *ct, const uint8_t *sk);

#endif