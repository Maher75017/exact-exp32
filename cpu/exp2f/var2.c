// Variante "structure glibc" sans test, opérations toutes explicites (compiler avec -ffp-contract=off)
#include <stdint.h>
#include <string.h>
typedef union {float f; uint32_t u;} b32u32_u;
typedef union {double f; uint64_t u;} b64u64_u;
int VB = 6, VDEG = 5, VSCHEME = 0;
double VT[256], VC[10];
float cr_exp2f(float);
float v_exp2f(float x){
  b32u32_u t = {.f = x};
  if (__builtin_expect((t.u<<1) - 0x65000000u >= 0x86000000u - 0x65000000u, 0)) return cr_exp2f(x);
  const int B = VB, N = 1 << B;
  double xd = x, SH = __builtin_ldexp(1.5, 52 - B);
  b64u64_u kd = {.f = xd + SH};
  double h = xd - (kd.f - SH);
  b64u64_u sv = {.f = VT[kd.u & (N - 1)]};
  sv.u += (uint64_t)((int64_t)kd.u >> B) << 52;
  double h2 = h*h;
  double a = __builtin_fma(h, VC[1], VC[0]);
  double b = __builtin_fma(h, VC[3], VC[2]);
  if (VDEG == 5) b = __builtin_fma(h2, VC[4], b);
  double p = __builtin_fma(h2, b, a);
  double sh = sv.f*h;
  double r = __builtin_fma(sh, p, sv.f);
  return (float) r;
}
