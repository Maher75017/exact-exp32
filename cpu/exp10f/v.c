// exp10f sans test d'arrondi : variante à paramètres modifiables (compiler avec -ffp-contract=off)
// VSCHEME : 1 = réduction par fma (3 opérations), 2 = fma + correction de log2(10) (4 opérations)
#include <stdint.h>
typedef union {float f; uint32_t u;} b32u32_u;
typedef union {double f; uint64_t u;} b64u64_u;
int VB = 6, VDEG = 5, VSCHEME = 1;
double VT[256], VC[10];
double VL = 0x1.a934f0979a371p+1, VLLO = 0x1.7f2495fb7fa6dp-53;
float cr_exp10f(float);
float v_exp10f(float x){
  b32u32_u t = {.f = x};
  uint32_t ux = t.u<<1;
  if (__builtin_expect(ux>0x8466d3e8u || ux<0x72adf1c6u, 0)) return cr_exp10f(x);
  if (__builtin_expect(!(t.u<<12), 0)){ // cas exacts (x entier de 1 à 10), même test que CORE-MATH
    unsigned k = (t.u >> 20) - 1016;
    if(k<=26){ unsigned bt = 1<<k, msk = 0x7551101; if(bt&msk) return cr_exp10f(x); }
  }
  const double big = 0x1.8p46;
  double z = x;
  b64u64_u u = {.f = __builtin_fma(z, VL, big)};
  double kf = u.f - big;
  double h = __builtin_fma(z, VL, -kf);
  if (VSCHEME == 2) h = __builtin_fma(z, VLLO, h);
  b64u64_u sv = {.f = VT[u.u & 63]};
  sv.u += (uint64_t)((int64_t)u.u >> 6) << 52;
  double h2 = h*h;
  double a = __builtin_fma(h, VC[1], VC[0]);
  double b = __builtin_fma(h, VC[3], VC[2]);
  b = __builtin_fma(h2, VC[4], b);
  double p = __builtin_fma(h2, b, a);
  double sh = sv.f*h;
  return (float) __builtin_fma(sh, p, sv.f);
}
