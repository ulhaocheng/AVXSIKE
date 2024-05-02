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

#define NBITS_FIELD             751 
#define ALICE                   0
#define BOB                     1 
#define OALICE_BITS             372  
#define OBOB_BITS               379
#define MASK_ALICE              0x0F 
#define MASK_BOB                0x03 
#define GFP_BYTES               94

// fixed parameters for isogeny tree computation
#define MAX_INT_POINTS_ALICE    8        
#define MAX_INT_POINTS_BOB      10  
#define MAX_Alice               186
#define MAX_Bob                 239  

// batch parameters
#define INSTANCES               8
#define MSG_BYTES               32
#define MSG_VECTS               4
#define SECRETKEY_A_BYTES       47
#define SECRETKEY_B_BYTES       48
#define SK_A_VECTS              6
#define SK_B_VECTS              6

#define CRYPTO_SECRETKEYBYTES   644    // MSG_BYTES + SECRETKEY_B_BYTES + CRYPTO_PUBLICKEYBYTES bytes
#define CRYPTO_PUBLICKEYBYTES   564
#define CRYPTO_BYTES             32
#define CRYPTO_CIPHERTEXTBYTES  596    // CRYPTO_PUBLICKEYBYTES + MSG_BYTES bytes  


#endif 
