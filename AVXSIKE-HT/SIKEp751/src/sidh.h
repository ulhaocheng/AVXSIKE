/**
 *******************************************************************************
 * @version 0.1
 * @date 2021-12-16
 * @copyright Copyright © 2021 by University of Luxembourg
 * @author Hao Cheng
 *******************************************************************************
 */

#ifndef _SIDH_H
#define _SIDH_H

#include "curve.h"
#include "utils.h"
#include "random.h"

static const uint64_t A_gen[6*NWORDS] = {
  0x1E8D6E8001755, 0x43E4E77B08221, 0x38A5026DD2931, 0x33DBC16FB97BA,
  0x4104A6CE8B635, 0x5C55FD11BE7A0, 0x310D00C248088, 0x5544BCDD2F27D,
  0x3928D64956A47, 0x64C1F5AA8BAA1, 0x04B394A947907, 0x7868F3255CE1A,
  0x37FF14B12C88D, 0x70289007665AC, 0x0000077C3C352, // XPA0
  0x6F16772DEB272, 0x078A7242DD2BA, 0x6E25678D59334, 0x111BE520AC30F,
  0x70E1555C74970, 0x56F6DCE3593BF, 0x3F66BF6B5E7C3, 0x6EE399F0AD09A,
  0x6949B8F144BC4, 0x541B8DF5285C0, 0x78A915A68C3F2, 0x761623CE66BDA,
  0x60E0587BD515D, 0x748EE049B24EE, 0x0000F0737DE0A, // XPA1
  0x279870D015195, 0x4E08B9FFB5C3E, 0x6056DB3BB1739, 0x24EDE549F005A,
  0x3294E3FBA5887, 0x19170AC11C871, 0x7F26A6526E6A4, 0x472F99131DDA1,
  0x2935410BF7C9E, 0x23571C07343D4, 0x0C63A3CEC3A29, 0x374C6928B4F6E,
  0x68AA53400D3EC, 0x5A21054D143CE, 0x0000489E8F53B, // XQA0
  0x7E6EB7E60F3E9, 0x2C10FA03BCFB3, 0x62E3D68FDFF3A, 0x5BF217376CBE3,
  0x180C64F4FC539, 0x4592710E2FCF5, 0x4130DD174B908, 0x1B1B4667CB929,
  0x3E7B5EA48DFB2, 0x1AFE81EE0644E, 0x399B9A1ECEF3D, 0x63BCAE910E35B,
  0x4CF924FE81F7B, 0x2A35524BA0B16, 0x000110A1FAE00, // XQA1 
  0x018A4BA820B3E, 0x66828760FDC58, 0x1591F9AD2CB44, 0x39106BFDEED6D,
  0x2B8FCE0E4E184, 0x328693A0DA504, 0x0D301492D013B, 0x69A68F3A04365,
  0x5A4810537448C, 0x0BB4D87D165C0, 0x214AD7E020330, 0x5437330ADDD5F,
  0x17177C4B21C5F, 0x0E4DEAA75C946, 0x0000528809CF2, // XRA0
  0x729F3FA3987DD, 0x02A6CBF48FA6D, 0x17215B3ED9972, 0x3CE726192A254,
  0x3431D29F9408B, 0x706F0E39F368B, 0x5BA9243EF14B6, 0x7436DA39A1D91,
  0x2A3B98C6EC730, 0x04C59C5EAB331, 0x1B45231957F39, 0x5C9CCD5710C4C,
  0x58F793496DD50, 0x26FDE84B2B560, 0x00000E3F1200A, // XRA1
};

static const uint64_t B_gen[6*NWORDS] = {
  0x52355E802BF11, 0x6E30DBA58C615, 0x726EB713A3C62, 0x1DC77EF2A185D,
  0x528777DC8D7D7, 0x7B189E61E4B48, 0x5842598EB48D1, 0x4CBFA6833DBC4,
  0x1EBD271C80CEA, 0x1BB2C1C3AB516, 0x68BB1DD51F967, 0x1600FFF19D5E8,
  0x2AE9F824199DB, 0x380F7466CBBA9, 0x0000F9713A575, // XPB0
  0x0000000000000, 0x0000000000000, 0x0000000000000, 0x0000000000000, 
  0x0000000000000, 0x0000000000000, 0x0000000000000, 0x0000000000000,
  0x0000000000000, 0x0000000000000, 0x0000000000000, 0x0000000000000, 
  0x0000000000000, 0x0000000000000, 0x0000000000000, // XPB1
  0x7967F4A642268, 0x53DE817BB8239, 0x7237E1C6B3E72, 0x5F1EFDDC2C8BC,
  0x610B0FF5C6221, 0x6091BE88AC162, 0x66482D2EDEEB5, 0x5CED3E84DAC3A,
  0x084AA3A17E5F5, 0x3968901B7D921, 0x78BBAC27B38DC, 0x7FCCED051C907,
  0x6FC90A68E7145, 0x34F4AC019ED05, 0x00001403564A4, // XQB0
  0x0000000000000, 0x0000000000000, 0x0000000000000, 0x0000000000000, 
  0x0000000000000, 0x0000000000000, 0x0000000000000, 0x0000000000000,
  0x0000000000000, 0x0000000000000, 0x0000000000000, 0x0000000000000, 
  0x0000000000000, 0x0000000000000, 0x0000000000000, // XQB1
  0x48AE02DA3ADA3, 0x2CE3CE71246D5, 0x52ACA2901EF66, 0x4883DFB910AC8,
  0x153A1E2FA3FE4, 0x0AEEBD4570923, 0x36ED8E5B3EC2B, 0x32865FFDF5B2A,
  0x6B7981C813584, 0x52D631287517F, 0x0222BE7424E46, 0x18D2780F000BF,
  0x4C1DD2199526F, 0x560CC20A3F35C, 0x0000FADB0F5C3, // XRB0
  0x12CF133DA2DDE, 0x5DB7188A92D63, 0x0589340C0CCC1, 0x17303D1C1DCFE,
  0x7268FBD772A7C, 0x6CD134EA0B0B5, 0x694102A9594D5, 0x3BAE50A9F481B,
  0x2F011FD2FD1C9, 0x2621F0D5804DC, 0x25423560BE2E9, 0x4BEE487A1D19A,
  0x4BCD02B4AE04C, 0x70B8CA1D08878, 0x0001146291C6E, // XRB1  
};

// Fixed parameters for isogeny tree computation
static const unsigned int strat_Alice[MAX_Alice-1] = { 
  80, 48, 27, 15, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 7, 4, 2, 1, 1, 2, 1, 
  1, 3, 2, 1, 1, 1, 1, 12, 7, 4, 2, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 5, 3, 2, 1, 1, 
  1, 1, 2, 1, 1, 1, 21, 12, 7, 4, 2, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 5, 3, 2, 1, 
  1, 1, 1, 2, 1, 1, 1, 9, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1, 
  33, 20, 12, 7, 4, 2, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 2, 1, 
  1, 1, 8, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 2, 1, 1, 16, 8, 4, 2, 1, 1, 
  1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1 };

static const unsigned int strat_Bob[MAX_Bob-1] = { 
  112, 63, 32, 16, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 
  1, 4, 2, 1, 1, 2, 1, 1, 16, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 
  1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 31, 16, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 
  1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 15, 8, 4, 2, 1, 1, 2, 1, 1, 4, 
  2, 1, 1, 2, 1, 1, 7, 4, 2, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 49, 31, 16, 8, 4, 2, 
  1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 
  15, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 7, 4, 2, 1, 1, 2, 1, 1, 3, 2, 1, 
  1, 1, 1, 21, 12, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 5, 3, 2, 1, 1, 1, 1, 
  2, 1, 1, 1, 9, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1 };

  

void random_mod_order_B(__m512i *random_digits, uint8_t *sk);

void EphemeralKeyGeneration_A(const __m512i *PrivateKeyA, __m512i *PublicKeyA);
void EphemeralKeyGeneration_B(const __m512i *PrivateKeyB, __m512i *PublicKeyB);
void EphemeralSecretAgreement_A(const __m512i *PrivateKeyA, const __m512i *PublicKeyB, __m512i *SharedSecretA);
void EphemeralSecretAgreement_B(const __m512i *PrivateKeyB, const __m512i *PublicKeyA, __m512i *SharedSecretB);

#endif

// p := 0x6FE5D541F71C0E12909F97BADC668562B5045CB25748084E9867D6EBE876DA959B1A13F7CC76E3EC968549F878A8EEAFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF;

// XPA0 := 0x4514F8CC94B140F24874F8B87281FA6004CA5B3637C68AC0C0BDB29838051F385FBBCC300BBB24BFBBF6710D7DC8B29ACB81E429BD1BD5629AD0ECAD7C90622F6BB801D0337EE6BC78A7F12FDCB09DECFAE8BFD643C89C3BAC1D87F8B6FA;
// XPA1 := 0x158ABF500B5914B3A96CED5FDB37D6DD925F2D6E4F7FEA3CC16E1085754077737EA6F8CC74938D971DA289DCF2435BCAC1897D2627693F9BB167DC01BE34AC494C60B8A0F65A28D7A31EA0D54640653A8099CE5A84E4F0168D818AF02041;
// XQA0 := 0x1723D2BFA01A78BF4E39E3A333F8A7E0B415A17F208D3419E7591D59D8ABDB7EE6D2B2DFCB21AC29A40F837983C0F057FD041AD93237704F1597D87F074F682961A38B5489D1019924F8A0EF5E4F1B2E64A7BA536E219F5090F76276290E;
// XQA1 := 0x2569D7EAFB6C60B244EF49E05B5E23F73C4F44169A7E02405E90CEB680CB0756054AC0E3DCE95E2950334262CC973235C2F87D89500BCD465B078BD0DEBDF322A2F86AEDFDCFEE65C09377EFBA0C5384DD837BEDB710209FBC8DDB8C35C7;
// XRA0 := 0x6066E07F3C0D964E8BC963519FAC8397DF477AEA9A067F3BE343BC53C883AF29CCF008E5A30719A29357A8C33EB3600CD078AF1C40ED5792763A4D213EBDE44CC623195C387E0201E7231C529A15AF5AB743EE9E7C9C37AF3051167525BB;
// XRA1 := 0x50E30C2C06494249BC4A144EB5F31212BD05A2AF0CB3064C322FC3604FC5F5FE3A08FB3A02B05A48557E15C992254FFC8910B72B8E1328B4893CDCFBFC003878881CE390D909E39F83C5006E0AE979587775443483D13C65B107FADA5165;

// XPB0 := 0x605D4697A245C394B98024A5554746DC12FF56D0C6F15D2F48123B6D9C498EEE98E8F7CD6E216E2F1FF7CE0C969CCA29CAA2FAA57174EF985AC0A504260018760E9FDF67467E20C13982FF5B49B8BEAB05F6023AF873F827400E453432FE;
// XQB0 := 0x5BF9544781803CBD7E0EA8B96D934C5CBCA970F9CC327A0A7E4DAD931EC29BAA8A854B8A9FDE5409AF96C5426FA375D99C68E9AE714172D7F04502D45307FA4839F39A28338BBAFD54A461A535408367D5132E6AA0D3DA6973360F8CD0F1;
// XRB0 := 0x55E5124A05D4809585F67FE9EA1F02A06CD411F38588BB631BF789C3F98D1C3325843BB53D9B011D8BD1F682C0E4D8A5E723364364E40DAD1B7A476716AC7D1BA705CCDD680BFD4FE4739CC21A9A59ED544B82566BF633E8950186A79FE3;
// XRB1 := 0x5AC57EAFD6CC7569E8B53A148721953262C5B404C143380ADCC184B6C21F0CAFE095B7E9C79CA88791F9A72F1B2F3121829B2622515B694A16875ED637F421B539E66F2FEF1CE8DCEFC8AEA608055E9C44077266AB64611BF851BA06C821;