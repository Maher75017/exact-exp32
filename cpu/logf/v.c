// logf sans test d'arrondi, structure du log2f de CORE-MATH (compiler avec -ffp-contract=off)
#include <stdint.h>
typedef union {float f; uint32_t u;} b32u32_u;
typedef union {double f; uint64_t u;} b64u64_u;
double VC[6], VIX[129], VLIX[129];
int VONE = 1, VSTRUCT = 0, VGRID = 0; long VNFALL = 0; // 1 : renvoie +0 pour x = 1
float cr_logf(float);
float v_logf(float x){
  b32u32_u t = {.f = x};
  uint32_t ux = t.u;
  if (__builtin_expect(ux < 1u<<23 || ux >= 0xffu<<23, 0)) return cr_logf(x);
  uint64_t m = ux&(~0u>>9); m <<= 52-23;
  int e = (ux>>23) - 0x7f;
  if (VONE && __builtin_expect(!m && !e, 0)) return 0.0f;
  int j = (m + (1l<<(52-8)))>>(52-7), k = j>53;
  e += k;
  b64u64_u xd = {.u = m | 0x3ffull<<52};
  double z = __builtin_fma(xd.f, VIX[j], -1.0); // exact
  double z2 = z*z;
  double c0 = __builtin_fma(z, VC[1], VC[0]);
  double c2 = __builtin_fma(z, VC[3], VC[2]);
  double c4 = __builtin_fma(z, VC[5], VC[4]);
  c2 = __builtin_fma(z2, c4, c2);
  c0 = __builtin_fma(z2, c2, c0);
  const double ln2 = 0x1.62e42fefa39efp-1, ln2l = 0x1.abc9e3b39803fp-56;
  b64u64_u r;
  if (VSTRUCT == 0){
    double y0 = __builtin_fma((double)e, ln2, -VLIX[j]);
    r.f = __builtin_fma(z, c0, y0);
  } else {
    double s = __builtin_fma(z, c0, -VLIX[j]);   // |s| < 0.36
    if (VSTRUCT == 2) s = __builtin_fma((double)e, ln2l, s);
    r.f = __builtin_fma((double)e, ln2, s);      // un seul arrondi final pour les grands e
  }
  // contrôle léger : r à moins de VGRID ulp double d'un float ou d'un milieu de floats -> calcul sûr
  if (VGRID && __builtin_expect(((r.u + VGRID) & 0x0fffffffu) < 2u*VGRID, 0)){ VNFALL++; return cr_logf(x); }
  return r.f;
}
