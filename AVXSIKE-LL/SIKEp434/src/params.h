/**
 *******************************************************************************
 * @version 0.1
 * @date 2021-12-16
 * @copyright Copyright © 2021 by University of Luxembourg
 * @author Hao Cheng
 *******************************************************************************
 */

#ifndef _PARAMS_H
#define _PARAMS_H

#define NBITS_FIELD             434 
#define ALICE                   0
#define BOB                     1 
#define OALICE_BITS             216  
#define OBOB_BITS               218
#define MASK_ALICE              0xFF 
#define MASK_BOB                0x01 
#define GFP_BYTES               55
#define NWORDS_ORDER            4
#define FP2_ENCODED_BYTES       2*((NBITS_FIELD + 7) / 8)

// fixed parameters for isogeny tree computation
#define MAX_INT_POINTS_ALICE    7        
#define MAX_INT_POINTS_BOB      8  
#define MAX_Alice               108
#define MAX_Bob                 137  

// batch parameters
#define INSTANCES               8
#define MSG_BYTES               16
#define MSG_VECTS               2
#define SECRETKEY_A_BYTES       27
#define SECRETKEY_B_BYTES       28
#define SK_A_VECTS              4
#define SK_B_VECTS              4

#define CRYPTO_SECRETKEYBYTES   374    // MSG_BYTES + SECRETKEY_B_BYTES + CRYPTO_PUBLICKEYBYTES bytes
#define CRYPTO_PUBLICKEYBYTES   330
#define CRYPTO_BYTES             16
#define CRYPTO_CIPHERTEXTBYTES  346    // CRYPTO_PUBLICKEYBYTES + MSG_BYTES bytes  


#endif 
