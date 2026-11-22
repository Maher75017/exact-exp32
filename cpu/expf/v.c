// expf sans test d'arrondi : variante à paramètres modifiables (compiler avec -ffp-contract=off)
// VSCHEME : 0 = réduction façon CORE-MATH (a = iln2*z arrondi),
//           1 = réduction par fma (3 opérations), 2 = fma + correction iln2lo (4 opérations)
#include <stdint.h>
typedef union {float f; uint32_t u;} b32u32_u;
typedef union {double f; uint64_t u;} b64u64_u;
int VB = 6, VDEG = 5, VSCHEME = 1;
double VT[256], VC[10];
double VILN2 = 0x1.71547652b82fep+0, VILN2LO = 0x1.777d0ffda0d24p-56;
float cr_expf(float);
float v_expf(float x){
  b32u32_u t = {.f = x};
  uint32_t ux = t.u<<1;
  if (__builtin_expect(ux>0x859d1d80u || ux<0x6f93813eu, 0)) return cr_expf(x);
  const double big = 0x1.8p46; // ulp = 2^-6 : arrondit à un multiple de 1/64
  double z = x, h;
  b64u64_u u;
  if (VSCHEME == 0){
    double a = VILN2*z; u.f = a + big; double ia = big - u.f; h = a + ia;
  } else {
    u.f = __builtin_fma(z, VILN2, big);
    double kf = u.f - big;
    h = __builtin_fma(z, VILN2, -kf);
    if (VSCHEME == 2) h = __builtin_fma(z, VILN2LO, h);
  }
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
