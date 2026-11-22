/* Correctly-rounded 2^x function for binary32 value.

Copyright (c) 2023 Alexei Sibidanov.
Copyright (c) 2026 Maher El Hmissi (rounding-test-free fast path).

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

static float as_special(float x){
  b32u32_u t = {.f = x};
  uint32_t ux = t.u<<1;
  if(ux >= 0xffu<<24) { // x is inf or nan
    if(ux > 0xffu<<24) return x; // x = nan
    static const float ir[] = {__builtin_inff(), 0.0f};
    return ir[t.u>>31]; // x = +-inf
  }
  if(t.u>=0xc3150000u){
    double z = x, y = 0x1p-149 + (z + 149)*0x1p-150;
    y = __builtin_fmax(y, 0x1p-151);
    float r = y;
    if(r==0.0f) errno = ERANGE;
    return r;
  } 
  float r = 0x1p127f * 0x1p127f;
  if(r>0x1.fffffep127f) errno = ERANGE;
  return r;
}

/* No rounding test and no accurate path: the final rounding of r is correct
   for all binary32 inputs and all four rounding modes (checked exhaustively).
   tb[i] approximates 2^(i/64), c[] approximates (2^h - 1)/h on |h| <= 2^-6
   (in directed rounding modes the reduced argument h can reach 2^-6).
   Starting from the correctly rounded tb[i] and a near-minimax polynomial of
   degree 5, two table entries (i=13: -2 ulps, i=62: +2 ulps) and two
   coefficients (c[0]: +111 ulps, c[1]: +55 ulps) were adjusted by a search
   over the hard-to-round inputs (distance < 2^-22 ulp to a rounding boundary,
   computed with MPFR), so that every hard case rounds correctly.
   All operations are explicit, so the result does not depend on FMA contraction. */
float cr_exp2f(float x){
  static const double tb[] =
    {0x1.0000000000000p+0, 0x1.02c9a3e778061p+0, 0x1.059b0d3158574p+0, 0x1.0874518759bc8p+0,
     0x1.0b5586cf9890fp+0, 0x1.0e3ec32d3d1a2p+0, 0x1.11301d0125b51p+0, 0x1.1429aaea92de0p+0,
     0x1.172b83c7d517bp+0, 0x1.1a35beb6fcb75p+0, 0x1.1d4873168b9aap+0, 0x1.2063b88628cd6p+0,
     0x1.2387a6e756238p+0, 0x1.26b4565e27cdbp+0, 0x1.29e9df51fdee1p+0, 0x1.2d285a6e4030bp+0,
     0x1.306fe0a31b715p+0, 0x1.33c08b26416ffp+0, 0x1.371a7373aa9cbp+0, 0x1.3a7db34e59ff7p+0,
     0x1.3dea64c123422p+0, 0x1.4160a21f72e2ap+0, 0x1.44e086061892dp+0, 0x1.486a2b5c13cd0p+0,
     0x1.4bfdad5362a27p+0, 0x1.4f9b2769d2ca7p+0, 0x1.5342b569d4f82p+0, 0x1.56f4736b527dap+0,
     0x1.5ab07dd485429p+0, 0x1.5e76f15ad2148p+0, 0x1.6247eb03a5585p+0, 0x1.6623882552225p+0,
     0x1.6a09e667f3bcdp+0, 0x1.6dfb23c651a2fp+0, 0x1.71f75e8ec5f74p+0, 0x1.75feb564267c9p+0,
     0x1.7a11473eb0187p+0, 0x1.7e2f336cf4e62p+0, 0x1.82589994cce13p+0, 0x1.868d99b4492edp+0,
     0x1.8ace5422aa0dbp+0, 0x1.8f1ae99157736p+0, 0x1.93737b0cdc5e5p+0, 0x1.97d829fde4e50p+0,
     0x1.9c49182a3f090p+0, 0x1.a0c667b5de565p+0, 0x1.a5503b23e255dp+0, 0x1.a9e6b5579fdbfp+0,
     0x1.ae89f995ad3adp+0, 0x1.b33a2b84f15fbp+0, 0x1.b7f76f2fb5e47p+0, 0x1.bcc1e904bc1d2p+0,
     0x1.c199bdd85529cp+0, 0x1.c67f12e57d14bp+0, 0x1.cb720dcef9069p+0, 0x1.d072d4a07897cp+0,
     0x1.d5818dcfba487p+0, 0x1.da9e603db3285p+0, 0x1.dfc97337b9b5fp+0, 0x1.e502ee78b3ff6p+0,
     0x1.ea4afa2a490dap+0, 0x1.efa1bee615a27p+0, 0x1.f50765b6e4542p+0, 0x1.fa7c1819e90d8p+0};
  static const double c[] =
    {0x1.62e42fefa3a5ep-1, 0x1.ebfbdff8131fap-3, 0x1.c6b08d70400d0p-5, 0x1.3b2b1bee88b09p-7, 0x1.5d884e708e6aap-10};
  b32u32_u t = {.f = x};
  uint32_t ux = t.u<<1;
  if (__builtin_expect(ux>=0x86000000u || ux<0x65000000u, 0)){
    if(__builtin_expect(ux<0x65000000u, 1)) return 1.0f + x; // |x| < 2^-26
    if(!(t.u>=0xc3000000 && t.u<0xc3150000u)) return as_special(x);
    // -150 < x <= -128: the result is subnormal, the generic path below is correct
  }
  const double shift = 0x1.8p46; // 1.5*2^52/64: rounds x*64 to an integer k
  double xd = x;
  b64u64_u kd = {.f = xd + shift};
  double h = xd - (kd.f - shift); // exact, |h| <= 2^-6 in all rounding modes
  b64u64_u sv = {.f = tb[kd.u & 63]};
  sv.u += (uint64_t)((int64_t)kd.u >> 6) << 52; // multiply by 2^floor(k/64)
  double h2 = h*h;
  double a = __builtin_fma(h, c[1], c[0]);
  double b = __builtin_fma(h, c[3], c[2]);
  b = __builtin_fma(h2, c[4], b);
  double p = __builtin_fma(h2, b, a);
  double sh = sv.f*h;
  return __builtin_fma(sh, p, sv.f); // implicit conversion to float: the only rounding that matters
}
