// cr2.h : exp2f à arrondi correct (arrondi au plus proche), calcul en double, utilisable CPU et GPU
#pragma once
#include <stdint.h>
#include <string.h>
#include <math.h>
#ifdef __CUDACC__
#define HD __host__ __device__ __forceinline__
#else
#define HD static inline
#endif
HD uint64_t d2u(double d){
#ifdef __CUDA_ARCH__
  return (uint64_t)__double_as_longlong(d);
#else
  uint64_t u; memcpy(&u,&d,8); return u;
#endif
}
HD double u2d(uint64_t u){
#ifdef __CUDA_ARCH__
  return __longlong_as_double((long long)u);
#else
  double d; memcpy(&d,&u,8); return d;
#endif
}
HD uint32_t f2u(float f){
#ifdef __CUDA_ARCH__
  return __float_as_uint(f);
#else
  uint32_t u; memcpy(&u,&f,4); return u;
#endif
}
HD float u2f(uint32_t u){
#ifdef __CUDA_ARCH__
  return __uint_as_float(u);
#else
  float f; memcpy(&f,&u,4); return f;
#endif
}
static const double CR2_TB[64] = {
  0x1.0000000000000p+0, 0x1.02c9a3e778061p+0, 0x1.059b0d3158574p+0, 0x1.0874518759bc8p+0,
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
HD float cr2_exp2f(float x, const double* tb){
  const double c0 = 0x1.62e42fefa3a5ep-1, c1 = 0x1.ebfbdff8131fap-3, c2 = 0x1.c6b08d70400d0p-5,
               c3 = 0x1.3b2b1bee88b09p-7, c4 = 0x1.5d884e708e6aap-10;
  uint32_t u = f2u(x), ux = u << 1;
  if (ux >= 0x86000000u || ux < 0x65000000u){
    if (ux < 0x65000000u) return 1.0f + x;                    // |x| < 2^-26
    if (!(u >= 0xc3000000u && u < 0xc3150000u)){               // cas spéciaux (CORE-MATH as_special)
      if (ux >= (0xffu << 24)){ if (ux > (0xffu << 24)) return x + x; return (u >> 31) ? 0.0f : x; }
      if (u >= 0xc3150000u){ double z = x, y = 0x1p-149 + (z + 149.0) * 0x1p-150; y = fmax(y, 0x1p-151); return (float)y; }
      return 0x1p127f * 0x1p127f;
    }
  }
  double xd = x; const double shift = 0x1.8p46;
  double kdd = xd + shift; uint64_t kd = d2u(kdd);
  double h = xd - (kdd - shift);
  double s = u2d(d2u(tb[kd & 63]) + ((uint64_t)((int64_t)kd >> 6) << 52));
  double h2 = h * h;
  double a = fma(h, c1, c0), b = fma(h, c3, c2);
  b = fma(h2, c4, b);
  double p = fma(h2, b, a);
  double sh = s * h;
  return (float) fma(sh, p, s);
}
