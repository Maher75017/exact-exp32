# Écrit expf.c (format CORE-MATH) à partir d'un fichier de paramètres réglés
import sys
B,d,s = 6,5,None
L = open(sys.argv[1]).read().split(); s = int(L[2])
C = L[3:3+d]; T = L[3+d:3+d+64]
src = open("cm_expf_ref.c").read()
head = src[:src.index("float cr_expf(float x){")]
head = head.replace("Copyright (c) 2023 Alexei Sibidanov.", "Copyright (c) 2023 Alexei Sibidanov.\nCopyright (c) 2026 Maher El Hmissi (rounding-test-free fast path).")
def fmt(v): return float.fromhex(v).hex() if False else "0x%sp%s" % (v.split("0x")[1].split("p")[0].ljust(15,"0"), v.split("p")[1] if v.split("p")[1].startswith("-") else "+"+v.split("p")[1].lstrip("+"))
tb = ",\n     ".join(", ".join(fmt(T[i+j]) for j in range(4)) for i in range(0,64,4))
cc = ", ".join(fmt(c) for c in C)
red = {1: """  b64u64_u u = {.f = __builtin_fma(z, iln2, big)}; // rounds x/log(2) to a multiple of 1/64
  double kf = u.f - big;                              // exact
  double h = __builtin_fma(z, iln2, -kf);             // |h| <= 2^-6 in all rounding modes""",
       2: """  b64u64_u u = {.f = __builtin_fma(z, iln2, big)}; // rounds x/log(2) to a multiple of 1/64
  double kf = u.f - big;                              // exact
  double h = __builtin_fma(z, iln2, -kf);             // |h| <= 2^-6 in all rounding modes
  h = __builtin_fma(z, iln2l, h);                     // add the low part of 1/log(2)"""}[s]
note = {1: "two table entries (i=5 and i=6, -2 ulps each) were adjusted",
        2: "no adjustment was needed"}[s]
body = f'''/* No rounding test and no accurate path: the final rounding of the result is
   correct for all binary32 inputs and all four rounding modes (checked
   exhaustively against the previous CORE-MATH code).
   x/log(2) = k/64 + h with |h| <= 2^-6 (h is computed with fma, the product
   z*iln2 is not rounded), tb[i] approximates 2^(i/64) and c[] approximates
   (2^h - 1)/h. Starting from the correctly rounded tb[i] and a near-minimax
   polynomial of degree 5, {note} by a search over the hard-to-round inputs
   (distance < 2^-22 ulp to a rounding boundary, computed with MPFR).
   All operations are explicit, so the result does not depend on FMA contraction. */
float cr_expf(float x){{
  static const double tb[] =
    {{{tb}}};
  static const double c[] =
    {{{cc}}};
  const double iln2 = 0x1.71547652b82fep+0, iln2l = 0x1.777d0ffda0d24p-56, big = 0x1.8p46;
  b32u32_u t = {{.f = x}};
  double z = x;
  uint32_t ux = t.u<<1;
  if (__builtin_expect(ux>0x859d1d80u || ux<0x6f93813eu, 0)){{
    if(__builtin_expect(ux<0x6f93813eu, 1)) return 1.0 + z*(1 + z*0.5);
    if(ux >= 0xffu<<24) {{ // x is inf or nan
      if(ux > 0xffu<<24) return x; // x = nan
      static const float ir[] = {{__builtin_inff(), 0.0f}};
      return ir[t.u>>31]; // x = +-inf
    }}
    if(t.u>0xc2ce8ec0u){{
      double y = 0x1p-149 + (z + 0x1.9d1d9fccf477p+6)*0x1.71547652b82edp-150;
      y = __builtin_fmax(y, 0x1p-151);
      float r = y;
      if(r==0.0f) errno = ERANGE;
      return r;
    }}
    if(t.u>0x42b17217u){{
      float r = 0x1p127f * 0x1p127f;
      if(r>0x1.fffffep127f) errno = ERANGE;
      return r;
    }}
  }}
{red}
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
if s == 1: body = body.replace(", iln2l = 0x1.777d0ffda0d24p-56", "")
open(sys.argv[2],"w").write(head + body)
