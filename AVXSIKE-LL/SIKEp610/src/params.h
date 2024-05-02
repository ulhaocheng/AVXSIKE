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

#define NBITS_FIELD             610 
#define ALICE                   0
#define BOB                     1 
#define OALICE_BITS             305  
#define OBOB_BITS               305
#define MASK_ALICE              0x01 
#define MASK_BOB                0xFF 
#define GFP_BYTES               77
#define NWORDS_ORDER            5
#define FP2_ENCODED_BYTES       2*((NBITS_FIELD + 7) / 8)


// fixed parameters for isogeny tree computation
#define MAX_INT_POINTS_ALICE    8        
#define MAX_INT_POINTS_BOB      10  
#define MAX_Alice               152
#define MAX_Bob                 192  

// batch parameters
#define INSTANCES               8
#define MSG_BYTES               24
#define MSG_VECTS               3
#define SECRETKEY_A_BYTES       39
#define SECRETKEY_B_BYTES       38
#define SK_A_VECTS              5
#define SK_B_VECTS              5

#define CRYPTO_SECRETKEYBYTES   524    // MSG_BYTES + SECRETKEY_B_BYTES + CRYPTO_PUBLICKEYBYTES bytes
#define CRYPTO_PUBLICKEYBYTES   462
#define CRYPTO_BYTES             24
#define CRYPTO_CIPHERTEXTBYTES  486    // CRYPTO_PUBLICKEYBYTES + MSG_BYTES bytes  


#endif 
