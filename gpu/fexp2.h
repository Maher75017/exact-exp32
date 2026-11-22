// fexp2.h : exp2f à arrondi correct (au plus proche) calculé EN FLOAT (double-float),
// avec test d'arrondi ; les cas rares partent vers cr2_exp2f (calcul en double, exact).
#pragma once
#include "cr2.h"
#include "tabf.h"
#ifdef __CUDA_ARCH__
#define FADD(a,b) __fadd_rn(a,b)
#define FSUB(a,b) __fsub_rn(a,b)
#define FMUL(a,b) __fmul_rn(a,b)
#define FFMA(a,b,c) __fmaf_rn(a,b,c)
#else
#define FADD(a,b) ((a)+(b))
#define FSUB(a,b) ((a)-(b))
#define FMUL(a,b) ((a)*(b))
#define FFMA(a,b,c) fmaf(a,b,c)
#endif
#ifndef FEXP2_E
#define FEXP2_E 0x1p-40f      // marge du test d'arrondi (erreur relative garantie à vérifier)
#endif
// renvoie 1 si le chemin rapide a conclu (résultat dans *res), 0 sinon
HD int fexp2_fast(float x, const float* tf, float* res){   // tf : table entrelacée (hi, lo)
  uint32_t u = f2u(x), ux = u << 1;
  // domaine du chemin rapide : 2^-26 <= |x| < 128 et x > -125 (résultat normal)
  if (ux < 0x65000000u || ux >= 0x85fffc00u || u > 0xc2fa0000u) return 0;   // x >= 127.996 : 2^e déborderait
  const float MAGIC = 0x1.8p23f;
  float t  = FMUL(x, 128.0f);                 // exact
  float tk = FADD(t, MAGIC);
  int32_t k = (int32_t)(f2u(tk) - f2u(MAGIC)); // k = arrondi(128 x)
  float kf = FSUB(tk, MAGIC);
  float h  = FMUL(FSUB(t, kf), 0x1p-7f);      // exact, |h| <= 2^-8
  int i = k & 127, e = k >> 7;
  float th = tf[2*i], tl = tf[2*i+1];
  // q = 2^h - 1 = (ph + pe) + qlo
  float ph = FMUL(h, C1H);
  float pe = FFMA(h, C1H, -ph);               // erreur exacte du produit
  float h2 = FMUL(h, h);
  float w  = FFMA(h, C4, C3);
  w = FFMA(h, w, C2);
  float qlo = FFMA(h, C1L, FFMA(h2, w, pe));
  // y/2^e = (th + tl)(1 + ph + qlo)
  float mh = FMUL(th, ph);
  float ml = FFMA(th, ph, -mh);               // th*ph = mh + ml exactement
  float sh = FADD(th, mh);
  float sl = FSUB(mh, FSUB(sh, th));          // Fast2Sum (|th| >= 1 > |mh|)
  float L  = FADD(sl, ml);
  L = FADD(L, tl);
  L = FFMA(tl, ph, L);
  L = FFMA(th, qlo, L);
  float eps = FMUL(sh, FEXP2_E);
  float r1 = FADD(sh, FADD(L, eps)), r2 = FADD(sh, FSUB(L, eps));
  if (r1 != r2) return 0;                     // trop près d'une frontière d'arrondi
  *res = FMUL(r1, u2f((uint32_t)(e + 127) << 23)); // exact (résultat normal)
  return 1;
}
HD float fexp2f(float x, const float* tf, const double* tb){
  float r;
  if (fexp2_fast(x, tf, &r)) return r;
  return cr2_exp2f(x, tb);
}
