/* Correctly-rounded logarithm function for binary32 value.

Copyright (c) 2023 Alexei Sibidanov.
Copyright (c) 2026 Maher El Hmissi (fast path without rounding test).

This file is part of the CORE-MATH project
(https://core-math.gitlabpages.inria.fr/).

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#include <stdint.h>
#include <errno.h>

// Warning: clang also defines __GNUC__
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic ignored "-Wunknown-pragmas"
#endif

#pragma STDC FENV_ACCESS ON

typedef union {float f; uint32_t u;} b32u32_u;
typedef union {double f; uint64_t u;} b64u64_u;
static __attribute__((noinline)) float as_special(float x){
  b32u32_u t = {.f = x};
  uint32_t ux = t.u;
  if(ux == 0u){// +0.0
    errno = ERANGE;
    return -1.0f/0.0f; // to raise FE_DIVBYZERO
  }
  if(ux == 0x7f800000u) return x; // +inf
  uint32_t ax = ux<<1;
  if(ax == 0u) { // -0.0
    errno = ERANGE;
    return -1.0f/0.0f; // to raise FE_DIVBYZERO
  }
  if(ax > 0xff000000u) return x; // nan
  errno = EDOM;
  return 0.0f/0.0f; // to raise FE_INVALID and return nan
}

/* previous CORE-MATH code, used for special inputs and as a fallback */
static float logf_accurate(float x){
  static const double tr[] = {
    0x1p+0, 0x1.f81f82p-1, 0x1.f07c1fp-1, 0x1.e9131acp-1,
    0x1.e1e1e1ep-1, 0x1.dae6077p-1, 0x1.d41d41dp-1, 0x1.cd85689p-1,
    0x1.c71c71cp-1, 0x1.c0e0704p-1, 0x1.bacf915p-1, 0x1.b4e81b5p-1,
    0x1.af286bdp-1, 0x1.a98ef6p-1, 0x1.a41a41ap-1, 0x1.9ec8e95p-1,
    0x1.999999ap-1, 0x1.948b0fdp-1, 0x1.8f9c19p-1, 0x1.8acb90fp-1,
    0x1.8618618p-1, 0x1.8181818p-1, 0x1.7d05f41p-1, 0x1.78a4c81p-1,
    0x1.745d174p-1, 0x1.702e05cp-1, 0x1.6c16c17p-1, 0x1.6816817p-1,
    0x1.642c859p-1, 0x1.605816p-1, 0x1.5c9882cp-1, 0x1.58ed231p-1,
    0x1.5555555p-1, 0x1.51d07ebp-1, 0x1.4e5e0a7p-1, 0x1.4afd6ap-1,
    0x1.47ae148p-1, 0x1.446f865p-1, 0x1.4141414p-1, 0x1.3e22cbdp-1,
    0x1.3b13b14p-1, 0x1.3813814p-1, 0x1.3521cfbp-1, 0x1.323e34ap-1,
    0x1.2f684bep-1, 0x1.2c9fb4ep-1, 0x1.29e412ap-1, 0x1.27350b9p-1,
    0x1.2492492p-1, 0x1.21fb781p-1, 0x1.1f7047ep-1, 0x1.1cf06aep-1,
    0x1.1a7b961p-1, 0x1.1811812p-1, 0x1.15b1e5fp-1, 0x1.135c811p-1,
    0x1.1111111p-1, 0x1.0ecf56cp-1, 0x1.0c9715p-1, 0x1.0a6810ap-1,
    0x1.0842108p-1, 0x1.0624dd3p-1, 0x1.041041p-1, 0x1.0204081p-1, 0.5};
  static const double tl[] = {
    -0x1.3b40815cd0628p-45, 0x1.fc0a890fbb514p-7, 0x1.f829b1e780b98p-6, 0x1.77458f532c948p-5,
    0x1.f0a30c2114ef2p-5, 0x1.341d793bbc7f7p-4, 0x1.6f0d28d256172p-4, 0x1.a926d3a6acb89p-4,
    0x1.e2707722ae90cp-4, 0x1.0d77e7a90896cp-3, 0x1.29552f6fff036p-3, 0x1.44d2b6c5b7831p-3,
    0x1.5ff306ee78ee7p-3, 0x1.7ab890410d41cp-3, 0x1.9525a9e3451c7p-3, 0x1.af3c94ed0bb06p-3,
    0x1.c8ff7c59a9535p-3, 0x1.e27076d5aedf9p-3, 0x1.fb9186b5e393ep-3, 0x1.0a324e38b8e6dp-2,
    0x1.1675cacaba398p-2, 0x1.22941fc0f76efp-2, 0x1.2e8e2bc311abap-2, 0x1.3a64c56b14373p-2,
    0x1.4618bc31c5c4cp-2, 0x1.51aad874df5b7p-2, 0x1.5d1bdbea80754p-2, 0x1.686c81d331238p-2,
    0x1.739d7f6dbcd9p-2, 0x1.7eaf83c82ad4dp-2, 0x1.89a3385813fe4p-2, 0x1.947941aa91484p-2,
    0x1.9f323edbf95d5p-2, 0x1.a9cec9a4205d3p-2, 0x1.b44f77c5c8cecp-2, 0x1.beb4d9ea71905p-2,
    0x1.c8ff7c69a97abp-2, 0x1.d32fe7f38e95fp-2, 0x1.dd46a0501c22ap-2, 0x1.e7442617e8511p-2,
    0x1.f128f5eaf0476p-2, 0x1.faf588dd8f0a8p-2, 0x1.02552a5edcfc4p-1, 0x1.0723e5c64de05p-1,
    0x1.0be72e3852947p-1, 0x1.109f39d554b5cp-1, 0x1.154c3d2c4d4aep-1, 0x1.19ee6b38bc834p-1,
    0x1.1e85f5ef03f95p-1, 0x1.23130d7fabe07p-1, 0x1.2795e1219afep-1, 0x1.2c0e9ec9c8d5p-1,
    0x1.307d7337f0f83p-1, 0x1.34e289cb4e098p-1, 0x1.393e0d42e28dep-1, 0x1.3d9026ad555bfp-1,
    0x1.41d8fe8667173p-1, 0x1.4618bc1ec5d87p-1, 0x1.4a4f85d303d8p-1, 0x1.4e7d8127f5a75p-1,
    0x1.52a2d26dbc47p-1, 0x1.56bf9d597f25ep-1, 0x1.5ad404cb59df2p-1, 0x1.5ee02a928153ap-1,
    0x1.62e42fefa38b4p-1};
  static const double b[] = {0x1.00000006342eap+0, -0x1.0001f7fdc3977p-1, 0x1.554a4e5cae9cfp-2};
  static const double c[] =
    {-0x1p-1, 0x1.55555555571cap-2, -0x1.0000000002d85p-2, 0x1.9999987d0c963p-3,
     -0x1.555554059a8bbp-3, 0x1.24aebcf71a38fp-3, -0x1.001c73915d758p-3};
  b32u32_u t = {.f = x};
  uint32_t ux = t.u;
  if(__builtin_expect(ux<(1<<23) || ux >= 0x7f800000u, 0)){
    if(ux==0 || ux >= 0x7f800000u) return as_special(x); // <=0, nan, inf
    // subnormal
    int n = __builtin_clz(ux) - 8;
    ux <<= n;
    ux -= n<<23;
  }
  uint32_t m = ux&((1<<23)-1), j = (m + (1<<(23-7)))>>(23-6);
  int32_t e = ((int32_t)ux>>23)-127;
  b64u64_u tz = {.u = ((uint64_t)m|(1023l<<23))<<(52-23)};
  double z = tz.f*tr[j] - 1, z2 = z*z;
  double r = ((e*0x1.62e42fefa39efp-1 + tl[j]) + z*b[0]) + z2*(b[1] + z*b[2]);
  float ub = r, lb = r + 0x1.f06p-33;
  if(__builtin_expect(ub != lb, 0)){
    double f = z2*((c[0] + z*c[1]) + z2*((c[2] + z*c[3]) + z2*(c[4] + z*c[5] + z2*c[6])));
    if(__builtin_expect(__builtin_fabsf(x-1.0f)<0x1p-10f, 0)) {
      if(x == 1.0f) return 0.0f;
      return z + f;
    }
    f -= 0x1.0ca86c3898dp-49*e;
    f += z;
    f += tl[j]-tl[0];
    double el = e*0x1.62e42fefa3ap-1, r = el + f;
    ub = r;
    tz.f = r;
    if(__builtin_expect(!(tz.u&((1u<<28)-1u)), 0) ){
      double dr = (el - r) + f;
      r += dr*64.0;
      ub = r;
    }
  }
  return ub;
}

/* Fast path without rounding test: z = x/2^e*ix[j] - 1 is exact (|z| < 2^-8),
   log(x) = e*log(2) - lix[j] + log1p(z), with lix[j] = log(ix[j]) (or log(2*ix[j])
   when e was incremented), and log1p(z) ~ z*(c[0] + ... + c[5]*z^5). The tables
   ix and the polynomial are those of CORE-MATH log2f, lix[j] are correctly rounded.
   For a few inputs (5 among the 2^31 positive floats), the double result r lies
   too close to a binary32 rounding boundary for any double computation to decide;
   the cheap check below catches them (and x = 1, where the sign of zero matters):
   if r is within one ulp of a binary32 number or of a midpoint, the accurate code
   is used. This happens for about 2^-27 of the inputs. The result is correct for
   all binary32 inputs in all four rounding modes (checked exhaustively). */
float cr_logf(float x){
  static const double ix[] = {
    0x1.0000000000000p+0, 0x1.fc07f01fc0000p-1, 0x1.f81f81f820000p-1, 0x1.f44659e4a0000p-1,
    0x1.f07c1f07c0000p-1, 0x1.ecc07b3020000p-1, 0x1.e9131abf00000p-1, 0x1.e573ac9020000p-1,
    0x1.e1e1e1e1e0000p-1, 0x1.de5d6e3f80000p-1, 0x1.dae6076ba0000p-1, 0x1.d77b654b80000p-1,
    0x1.d41d41d420000p-1, 0x1.d0cb58f6e0000p-1, 0x1.cd85689040000p-1, 0x1.ca4b3055e0000p-1,
    0x1.c71c71c720000p-1, 0x1.c3f8f01c40000p-1, 0x1.c0e0703820000p-1, 0x1.bdd2b89940000p-1,
    0x1.bacf914c20000p-1, 0x1.b7d6c3dda0000p-1, 0x1.b4e81b4e80000p-1, 0x1.b2036406c0000p-1,
    0x1.af286bca20000p-1, 0x1.ac5701ac60000p-1, 0x1.a98ef606a0000p-1, 0x1.a6d01a6d00000p-1,
    0x1.a41a41a420000p-1, 0x1.a16d3f97a0000p-1, 0x1.9ec8e95100000p-1, 0x1.9c2d14ee40000p-1,
    0x1.99999999a0000p-1, 0x1.970e4f80c0000p-1, 0x1.948b0fcd60000p-1, 0x1.920fb49d00000p-1,
    0x1.8f9c18f9c0000p-1, 0x1.8d3018d300000p-1, 0x1.8acb90f6c0000p-1, 0x1.886e5f0ac0000p-1,
    0x1.8618618620000p-1, 0x1.83c977ab20000p-1, 0x1.8181818180000p-1, 0x1.7f405fd020000p-1,
    0x1.7d05f417e0000p-1, 0x1.7ad2208e00000p-1, 0x1.78a4c81780000p-1, 0x1.767dce4340000p-1,
    0x1.745d1745e0000p-1, 0x1.724287f460000p-1, 0x1.702e05c0c0000p-1, 0x1.6e1f76b440000p-1,
    0x1.6c16c16c20000p-1, 0x1.6a13cd1540000p-1, 0x1.6816816820000p-1, 0x1.661ec6a520000p-1,
    0x1.642c8590c0000p-1, 0x1.623fa77020000p-1, 0x1.6058160580000p-1, 0x1.5e75bb8d00000p-1,
    0x1.5c9882b940000p-1, 0x1.5ac056b020000p-1, 0x1.58ed230820000p-1, 0x1.571ed3c500000p-1,
    0x1.5555555560000p-1, 0x1.5390948f40000p-1, 0x1.51d07eae20000p-1, 0x1.5015015020000p-1,
    0x1.4e5e0a7300000p-1, 0x1.4cab887260000p-1, 0x1.4afd6a0520000p-1, 0x1.49539e3b20000p-1,
    0x1.47ae147ae0000p-1, 0x1.460cbc7f60000p-1, 0x1.446f865620000p-1, 0x1.42d6625d60000p-1,
    0x1.4141414140000p-1, 0x1.3fb013fb00000p-1, 0x1.3e22cbce40000p-1, 0x1.3c995a47c0000p-1,
    0x1.3b13b13b20000p-1, 0x1.3991c2c180000p-1, 0x1.3813813820000p-1, 0x1.3698df3de0000p-1,
    0x1.3521cfb2c0000p-1, 0x1.33ae45b580000p-1, 0x1.323e34a2c0000p-1, 0x1.30d1901300000p-1,
    0x1.2f684bda20000p-1, 0x1.2e025c04c0000p-1, 0x1.2c9fb4d820000p-1, 0x1.2b404ad020000p-1,
    0x1.29e4129e40000p-1, 0x1.288b012880000p-1, 0x1.27350b8820000p-1, 0x1.25e2270800000p-1,
    0x1.24924924a0000p-1, 0x1.23456789a0000p-1, 0x1.21fb781220000p-1, 0x1.20b470c680000p-1,
    0x1.1f7047dc20000p-1, 0x1.1e2ef3b400000p-1, 0x1.1cf06ada20000p-1, 0x1.1bb4a40460000p-1,
    0x1.1a7b9611a0000p-1, 0x1.19453808c0000p-1, 0x1.1811811820000p-1, 0x1.16e0689420000p-1,
    0x1.15b1e5f760000p-1, 0x1.1485f0e0a0000p-1, 0x1.135c811360000p-1, 0x1.12358e75e0000p-1,
    0x1.1111111120000p-1, 0x1.0fef010fe0000p-1, 0x1.0ecf56be60000p-1, 0x1.0db20a8900000p-1,
    0x1.0c9714fbc0000p-1, 0x1.0b7e6ec260000p-1, 0x1.0a6810a680000p-1, 0x1.0953f39020000p-1,
    0x1.0842108420000p-1, 0x1.073260a480000p-1, 0x1.0624dd2f20000p-1, 0x1.05197f7d80000p-1,
    0x1.0410410420000p-1, 0x1.03091b5200000p-1, 0x1.0204081020000p-1, 0x1.0101010100000p-1,
    0x1.0000000000000p-1};
  static const double lix[] = {
    0x0.0p+0, -0x1.fe02a6b146789p-8, -0x1.fc0a8b0fa03e4p-7, -0x1.7b91b07de311bp-6,
    -0x1.f829b0e7c3300p-6, -0x1.39e87b9fd7d60p-5, -0x1.77458f63edcfcp-5, -0x1.b42dd7117b1bfp-5,
    -0x1.f0a30c01362a6p-5, -0x1.16536eea7fae1p-4, -0x1.341d7961791d1p-4, -0x1.51b073f07983fp-4,
    -0x1.6f0d28ae3eb4cp-4, -0x1.8c345d6383b21p-4, -0x1.a926d3a475563p-4, -0x1.c5e548f63a743p-4,
    -0x1.e27076e28f2e6p-4, -0x1.fec9131dbaabbp-4, -0x1.0d77e7ccf6e59p-3, -0x1.1b72ad52f87a0p-3,
    -0x1.29552f81eb523p-3, -0x1.371fc201f7f74p-3, -0x1.44d2b6ccbfd1ep-3, -0x1.526e5e3a41438p-3,
    -0x1.5ff3070a613d4p-3, -0x1.6d60fe717221dp-3, -0x1.7ab890212b909p-3, -0x1.87fa065214911p-3,
    -0x1.9525a9cf296b4p-3, -0x1.a23bc1fe42563p-3, -0x1.af3c94e81bff3p-3, -0x1.bc2867430acd6p-3,
    -0x1.c8ff7c7989a22p-3, -0x1.d5c216b535b91p-3, -0x1.e27076e2f92e6p-3, -0x1.ef0adcbe0d936p-3,
    -0x1.fb9186d5ebe2bp-3, -0x1.0402594b51041p-2, -0x1.0a324e27370e3p-2, -0x1.1058bf9ad7ad5p-2,
    -0x1.1675cabaa660ep-2, -0x1.1c898c16b91fbp-2, -0x1.22941fbcfb966p-2, -0x1.2895a13dd2ea3p-2,
    -0x1.2e8e2bade7d31p-2, -0x1.347dd9a9afd55p-2, -0x1.3a64c556b05eap-2, -0x1.40430868877e4p-2,
    -0x1.4618bc219dec2p-2, -0x1.4be5f9579e0a1p-2, -0x1.51aad872c982dp-2, -0x1.5767717432a6cp-2,
    -0x1.5d1bdbf5669cap-2, -0x1.62c82f2b83795p-2, 0x1.5d5bddf5b0f30p-2, 0x1.57bf753cb49fbp-2,
    0x1.522ae073b23d8p-2, 0x1.4c9e09e18f43cp-2, 0x1.4718dc271841bp-2, 0x1.419b423d5a8c7p-2,
    0x1.3c2527735f184p-2, 0x1.36b6776bff917p-2, 0x1.314f1e1d54ce4p-2, 0x1.2bef07cdb5354p-2,
    0x1.269621136db92p-2, 0x1.214456d0e88d4p-2, 0x1.1bf9963577b95p-2, 0x1.16b5ccbaf1373p-2,
    0x1.1178e822ae47cp-2, 0x1.0c42d67625ae3p-2, 0x1.07138604b0862p-2, 0x1.01eae56243e91p-2,
    0x1.f991c6cb33379p-3, 0x1.ef5ade4de2fe6p-3, 0x1.e530effe1b012p-3, 0x1.db13db0da1940p-3,
    0x1.d1037f264de7bp-3, 0x1.c6ffbc6ef8f71p-3, 0x1.bd087383798adp-3, 0x1.b31d8575dee3dp-3,
    0x1.a93ed3c8fd9e3p-3, 0x1.9f6c407055664p-3, 0x1.95a5adcfc217fp-3, 0x1.8beafeb38ce8cp-3,
    0x1.823c1655523c2p-3, 0x1.7898d85460c73p-3, 0x1.6f0128b7baabcp-3, 0x1.6574ebe86933ap-3,
    0x1.5bf406b59bdb2p-3, 0x1.527e5e4a5158dp-3, 0x1.4913d83395561p-3, 0x1.3fb45a59ed8ccp-3,
    0x1.365fcb0151016p-3, 0x1.2d1610c81c13ap-3, 0x1.23d712a4fa202p-3, 0x1.1aa2b7e1ff72ap-3,
    0x1.1178e822de47cp-3, 0x1.08598b5990a07p-3, 0x1.fe89139dc1566p-4, 0x1.ec739830d9120p-4,
    0x1.da7276390c6a2p-4, 0x1.c885801c04b23p-4, 0x1.b6ac88da61b1cp-4, 0x1.a4e7640a45c38p-4,
    0x1.9335e5d524989p-4, 0x1.8197e2f37a3f0p-4, 0x1.700d30af800e1p-4, 0x1.5e95a4d90f1cbp-4,
    0x1.4d3115d2cfeacp-4, 0x1.3bdf5a7c60e64p-4, 0x1.2aa04a44a57a5p-4, 0x1.1973bd1527567p-4,
    0x1.08598b5ac3a07p-4, 0x1.eea31bfea787cp-5, 0x1.ccb73cdcb32ccp-5, 0x1.aaef2d11110fcp-5,
    0x1.894aa1485b343p-5, 0x1.67c94f2e07b58p-5, 0x1.466aed42be3eap-5, 0x1.252f32faad83fp-5,
    0x1.0415d89e54444p-5, 0x1.c63d2ec16aaf2p-6, 0x1.8492528ddcabfp-6, 0x1.432a925ca0cc1p-6,
    0x1.0205658d15847p-6, 0x1.82448a3d8a2aap-7, 0x1.010157586de71p-7, 0x1.0080559488b35p-8,
    0x0.0p+0};
  static const double c[] =
    {0x1.0000000000000p+0, -0x1.0000000000000p-1, 0x1.55555555030bcp-2, -0x1.ffffffff2b4e5p-3, 0x1.999b5076a42f2p-3, -0x1.55570c45a647dp-3};
  b32u32_u t = {.f = x};
  uint32_t ux = t.u;
  if (__builtin_expect(ux < 1u<<23 || ux >= 0xffu<<23, 0)) return logf_accurate(x); // subnormal, <= 0, inf, nan
  uint64_t m = ux&(~0u>>9); m <<= 52-23;
  int e = (ux>>23) - 0x7f;
  int j = (m + (1l<<(52-8)))>>(52-7), k = j>53;
  e += k;
  b64u64_u xd = {.u = m | 0x3ffull<<52};
  double z = __builtin_fma(xd.f, ix[j], -1.0); // exact
  double z2 = z*z;
  double c0 = __builtin_fma(z, c[1], c[0]);
  double c2 = __builtin_fma(z, c[3], c[2]);
  double c4 = __builtin_fma(z, c[5], c[4]);
  c2 = __builtin_fma(z2, c4, c2);
  c0 = __builtin_fma(z2, c2, c0);
  const double ln2 = 0x1.62e42fefa39efp-1;
  double y0 = __builtin_fma((double)e, ln2, -lix[j]);
  b64u64_u r = {.f = __builtin_fma(z, c0, y0)};
  if (__builtin_expect(((r.u + 1) & 0x0fffffffu) < 2u, 0)) return logf_accurate(x);
  return r.f; // implicit conversion to float
}
