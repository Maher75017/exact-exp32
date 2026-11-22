# Écrit exp10f.c (format CORE-MATH) à partir d'un fichier de paramètres réglés
import sys
L = open(sys.argv[1]).read().split(); s = int(L[2]); d = 5
C = L[3:3+d]; T = L[3+d:3+d+64]
src = open("cm_exp10f_ref.c").read()
head = src[:src.index("#include <stdlib.h>")]
head = head.replace("Copyright (c) 2023 Alexei Sibidanov.", "Copyright (c) 2023 Alexei Sibidanov.\nCopyright (c) 2026 Maher El Hmissi (rounding-test-free fast path).")
inc = """#include <stdint.h>
#include <errno.h>

// Warning: clang also defines __GNUC__
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic ignored "-Wunknown-pragmas"
#endif

#pragma STDC FENV_ACCESS ON

typedef union {float f; uint32_t u;} b32u32_u;
typedef union {double f; uint64_t u;} b64u64_u;

"""
tb = ",\n     ".join(", ".join(T[i+j] for j in range(4)) for i in range(0,64,4))
cc = ", ".join(C)
fixnote = {1: "one table entry (i=40, +12 ulps) was adjusted by a search over\n   the hard-to-round inputs", 2: "no adjustment was needed for\n   the hard-to-round inputs"}[s]
red2 = "\n  h = __builtin_fma(z, l2l, h);                       // add the low part of log2(10)" if s == 2 else ""
l2l = ", l2l = 0x1.7f2495fb7fa6dp-53" if s == 2 else ""
body = f'''/* No rounding test and no accurate path: the final rounding of the result is
   correct for all binary32 inputs and all four rounding modes (checked
   exhaustively against the previous CORE-MATH code).
   x*log2(10) = k/64 + h with |h| <= 2^-6 (h is computed with fma, the product
   z*l2 is not rounded), tb[i] approximates 2^(i/64) and c[] approximates
   (2^h - 1)/h. Starting from the correctly rounded tb[i] and a near-minimax
   polynomial of degree 5, {fixnote} (distance < 2^-22 ulp to a rounding
   boundary, computed with MPFR). The exact cases x = 1, ..., 10 are handled
   separately. All operations are explicit, so the result does not depend on
   FMA contraction. */
float cr_exp10f(float x){{
  static const double tb[] =
    {{{tb}}};
  static const double c[] =
    {{{cc}}};
  static const float ex[] =
    {{10, 100, 1000, 10000, 100000, 1000000, 10000000, 100000000, 1000000000, 10000000000}};
  const double l2 = 0x1.a934f0979a371p+1{l2l}, big = 0x1.8p46;
  b32u32_u t = {{.f = x}};
  double z = x;
  uint32_t ux = t.u<<1;
  if (__builtin_expect(ux>0x8466d3e8u || ux<0x72adf1c6u, 0)){{
    if(__builtin_expect(ux<0x72adf1c6u, 1))
      return 1.0 + z*(0x1.26bb1bbb55516p+1 + z*(0x1.53524c73cea69p+1 + z*0x1.0470591de2ca4p+1));
    if(ux >= 0xffu<<24) {{ // x is inf or nan
      if(ux > 0xffu<<24) return x; // x = nan
      static const float ir[] = {{__builtin_inff(), 0.0f}};
      return ir[t.u>>31]; // x = +-inf
    }}
    if(t.u>0xc23369f4u){{
      double y = 0x1p-149 + (z + 0x1.66d3e7bd9a403p+5)*0x1.a934f0979a37p-149;
      y = __builtin_fmax(y, 0x1p-151);
      float r = y;
      if(r==0.0f) errno = ERANGE;
      return r;
    }}
    if(t.u>0x421a209au){{
      float r = 0x1p127f * 0x1p127f;
      if(r>0x1.fffffep127f) errno = ERANGE;
      return r;
    }}
  }}
  if(__builtin_expect(!(t.u<<12), 0)){{
    unsigned k = (t.u >> 20) - 1016;
    if(k<=26){{
      unsigned bt = 1<<k, msk = 0x7551101;
      if(bt&msk) return ex[__builtin_popcount(msk&(bt-1))];
    }}
  }}
  b64u64_u u = {{.f = __builtin_fma(z, l2, big)}};    // rounds x*log2(10) to a multiple of 1/64
  double kf = u.f - big;                              // exact
  double h = __builtin_fma(z, l2, -kf);               // |h| <= 2^-6 in all rounding modes{red2}
  b64u64_u sv = {{.f = tb[u.u&0x3f]}};
  sv.u += (uint64_t)((int64_t)u.u >> 6) << 52;        // multiply by 2^floor(k/64)
  double h2 = h*h;
  double a = __builtin_fma(h, c[1], c[0]);
  double b = __builtin_fma(h, c[3], c[2]);
  b = __builtin_fma(h2, c[4], b);
  double p = __builtin_fma(h2, b, a);
  double sh = sv.f*h;
  return __builtin_fma(sh, p, sv.f); // implicit conversion to float: the only rounding that matters
}}
'''
open(sys.argv[2],"w").write(head + inc + body)
