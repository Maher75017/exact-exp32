Activated conda/uv virtual environment at /venv/main
(main) root@C.61801330:/workspace$ cat best_progs.cuh
cat: best_progs.cuh: No such file or directory
(main) root@C.61801330:/workspace$ cat /workspace/c3c_bundle/best_progs.cuh
// Généré par c3b : chemins rapides exacts trouvés par programmation génétique (vérifiés exhaustivement)
#pragma once
#include "fexp2.h"
__device__ __forceinline__ int gp_fast_0(float x, const float* tf, float* res){
  uint32_t u = f2u(x), ux = u << 1;
  if (ux < 0x65000000u || ux >= 0x85fffc00u || u > 0xc2fa0000u) return 0;
  const float MAGIC = 0x1.8p23f;
  float a6 = 0.f, a7 = 0.f, a8 = 0.f, a9 = 0.f, a10 = 0.f, a11 = 0.f;
  a7 = FFMA(x, 128.0f, MAGIC);
  a8 = FSUB(a7, MAGIC);
  a10 = FFMA(-a8, 0x1p-7f, x);
  float tk = a7, h = a10;
  int32_t k = (int32_t)(f2u(tk) - f2u(MAGIC));
  if (k < -16384 || k > 16384) return 0;
  int i = k & 127, e = k >> 7;
  const float K2 = 0x1.ebfbdep-3f, K3 = 0x1.c6b2cap-5f, K4 = 0x1.3a2298p-7f;
  float r9 = 0.f, r10 = 0.f, r11 = 0.f, r12 = 0.f, r13 = 0.f, r14 = 0.f, r15 = 0.f;
  r12 = FFMA(h, K4, K3);
  r12 = FFMA(h, r12, K2);
  r9 = FMUL(h, C1H);
  r10 = FFMA(h, C1H, -r9);
  r14 = FFMA(h, r12, C1L);
  r11 = FFMA(tf[2*i], r9, tf[2*i]);
  r12 = FSUB(r11, tf[2*i]);
  r13 = FFMA(h, r14, r10);
  r12 = FFMA(tf[2*i], r9, -r12);
  r12 = FADD(r12, tf[2*i+1]);
  r12 = FFMA(tf[2*i], r13, r12);
  r12 = FFMA(tf[2*i+1], r9, r12);
  const float eps = 0x1p-40f;
  float r1 = FADD(r11, FADD(r12, eps)), r2 = FADD(r11, FSUB(r12, eps));
  if (!(r1 == r2)) return 0;
  *res = FMUL(r1, u2f((uint32_t)(e + 127) << 23));
  return 1;
}
__device__ __forceinline__ int gp_fast_1(float x, const float* tf, float* res){
  uint32_t u = f2u(x), ux = u << 1;
  if (ux < 0x65000000u || ux >= 0x85fffc00u || u > 0xc2fa0000u) return 0;
  const float MAGIC = 0x1.8p23f;
  float a6 = 0.f, a7 = 0.f, a8 = 0.f, a9 = 0.f, a10 = 0.f, a11 = 0.f;
  a7 = FFMA(x, 128.0f, MAGIC);
  a8 = FSUB(a7, MAGIC);
  a9 = FFMA(x, 128.0f, -a8);
  a10 = FMUL(a9, 0x1p-7f);
  float tk = a7, h = a10;
  int32_t k = (int32_t)(f2u(tk) - f2u(MAGIC));
  if (k < -16384 || k > 16384) return 0;
  int i = k & 127, e = k >> 7;
  const float K2 = 0x1.ebfbdep-3f, K3 = 0x1.c6b104p-5f, K4 = 0x1.3a306p-7f;
  float r9 = 0.f, r10 = 0.f, r11 = 0.f, r12 = 0.f, r13 = 0.f, r14 = 0.f, r15 = 0.f;
  r9 = FMUL(h, C1H);
  r10 = FFMA(h, C1H, -r9);
  r12 = FFMA(h, K4, K3);
  r12 = FFMA(h, r12, K2);
  r14 = FFMA(h, r12, C1L);
  r13 = FFMA(h, r14, r10);
  r11 = FFMA(tf[2*i], r9, tf[2*i]);
  r12 = FSUB(r11, tf[2*i]);
  r12 = FFMA(tf[2*i], r9, -r12);
  r12 = FADD(r12, tf[2*i+1]);
  r12 = FFMA(tf[2*i], r13, r12);
  r12 = FFMA(tf[2*i+1], r9, r12);
  const float eps = 0x1p-40f;
  float r1 = FADD(r11, FADD(r12, eps)), r2 = FADD(r11, FSUB(r12, eps));
  if (!(r1 == r2)) return 0;
  *res = FMUL(r1, u2f((uint32_t)(e + 127) << 23));
  return 1;
}
__device__ __forceinline__ int gp_fast_2(float x, const float* tf, float* res){
  uint32_t u = f2u(x), ux = u << 1;
  if (ux < 0x65000000u || ux >= 0x85fffc00u || u > 0xc2fa0000u) return 0;
  const float MAGIC = 0x1.8p23f;
  float a6 = 0.f, a7 = 0.f, a8 = 0.f, a9 = 0.f, a10 = 0.f, a11 = 0.f;
  a7 = FFMA(x, 128.0f, MAGIC);
  a8 = FSUB(a7, MAGIC);
  a9 = FFMA(x, 128.0f, -a8);
  a10 = FMUL(a9, 0x1p-7f);
  float tk = a7, h = a10;
  int32_t k = (int32_t)(f2u(tk) - f2u(MAGIC));
  if (k < -16384 || k > 16384) return 0;
  int i = k & 127, e = k >> 7;
  const float K2 = 0x1.ebfbdep-3f, K3 = 0x1.c6b086p-5f, K4 = 0x1.3a3044p-7f;
  float r9 = 0.f, r10 = 0.f, r11 = 0.f, r12 = 0.f, r13 = 0.f, r14 = 0.f, r15 = 0.f;
  r9 = FMUL(h, C1H);
  r10 = FFMA(h, C1H, -r9);
  r12 = FFMA(h, K4, K3);
  r12 = FFMA(h, r12, K2);
  r13 = FMUL(h, r12);
  r13 = FFMA(h, r13, r10);
  r13 = FFMA(h, C1L, r13);
  r11 = FFMA(tf[2*i], r9, tf[2*i]);
  r12 = FSUB(r11, tf[2*i]);
  r12 = FFMA(tf[2*i], r9, -r12);
  r12 = FADD(r12, tf[2*i+1]);
  r12 = FFMA(tf[2*i+1], r9, r12);
  r12 = FFMA(tf[2*i], r13, r12);
  const float eps = 0x1p-40f;
  float r1 = FADD(r11, FADD(r12, eps)), r2 = FADD(r11, FSUB(r12, eps));
  if (!(r1 == r2)) return 0;
  *res = FMUL(r1, u2f((uint32_t)(e + 127) << 23));
  return 1;
}
__device__ __forceinline__ int gp_fast_3(float x, const float* tf, float* res){
  uint32_t u = f2u(x), ux = u << 1;
  if (ux < 0x65000000u || ux >= 0x85fffc00u || u > 0xc2fa0000u) return 0;
  const float MAGIC = 0x1.8p23f;
  float a6 = 0.f, a7 = 0.f, a8 = 0.f, a9 = 0.f, a10 = 0.f, a11 = 0.f;
  a7 = FFMA(x, 128.0f, MAGIC);
  a8 = FSUB(a7, MAGIC);
  a9 = FFMA(x, 128.0f, -a8);
  a10 = FMUL(a9, 0x1p-7f);
  float tk = a7, h = a10;
  int32_t k = (int32_t)(f2u(tk) - f2u(MAGIC));
  if (k < -16384 || k > 16384) return 0;
  int i = k & 127, e = k >> 7;
  const float K2 = 0x1.ebfbep-3f, K3 = 0x1.c6b092p-5f, K4 = 0x1.3a3044p-7f;
  float r9 = 0.f, r10 = 0.f, r11 = 0.f, r12 = 0.f, r13 = 0.f, r14 = 0.f, r15 = 0.f;
  r9 = FMUL(h, C1H);
  r10 = FFMA(h, C1H, -r9);
  r11 = FMUL(h, h);
  r12 = FFMA(h, K4, K3);
  r12 = FFMA(h, r12, K2);
  r13 = FFMA(r11, r12, r10);
  r13 = FFMA(h, C1L, r13);
  r14 = FMUL(tf[2*i], r9);
  r15 = FFMA(tf[2*i], r9, -r14);
  r11 = FADD(tf[2*i], r14);
  r12 = FSUB(r11, tf[2*i]);
  r12 = FSUB(r14, r12);
  r12 = FADD(r12, r15);
  r12 = FADD(r12, tf[2*i+1]);
  r12 = FFMA(tf[2*i+1], r9, r12);
  r12 = FFMA(tf[2*i], r13, r12);
  const float eps = 0x1p-39f;
  float r1 = FADD(r11, FADD(r12, eps)), r2 = FADD(r11, FSUB(r12, eps));
  if (!(r1 == r2)) return 0;
  *res = FMUL(r1, u2f((uint32_t)(e + 127) << 23));
  return 1;
}
__device__ __forceinline__ int gp_fast_4(float x, const float* tf, float* res){
  uint32_t u = f2u(x), ux = u << 1;
  if (ux < 0x65000000u || ux >= 0x85fffc00u || u > 0xc2fa0000u) return 0;
  const float MAGIC = 0x1.8p23f;
  float a6 = 0.f, a7 = 0.f, a8 = 0.f, a9 = 0.f, a10 = 0.f, a11 = 0.f;
  a6 = FMUL(x, 128.0f);
  a7 = FADD(a6, MAGIC);
  a8 = FSUB(a7, MAGIC);
  a9 = FSUB(a6, a8);
  a10 = FMUL(a9, 0x1p-7f);
  float tk = a7, h = a10;
  int32_t k = (int32_t)(f2u(tk) - f2u(MAGIC));
  if (k < -16384 || k > 16384) return 0;
  int i = k & 127, e = k >> 7;
  const float K2 = 0x1.ebfbep-3f, K3 = 0x1.c6b092p-5f, K4 = 0x1.3a3044p-7f;
  float r9 = 0.f, r10 = 0.f, r11 = 0.f, r12 = 0.f, r13 = 0.f, r14 = 0.f, r15 = 0.f;
  r9 = FMUL(h, C1H);
  r10 = FFMA(h, C1H, -r9);
  r11 = FMUL(h, h);
  r12 = FFMA(h, K4, K3);
  r12 = FFMA(h, r12, K2);
  r13 = FFMA(r11, r12, r10);
  r13 = FFMA(h, C1L, r13);
  r14 = FMUL(tf[2*i], r9);
  r15 = FFMA(tf[2*i], r9, -r14);
  r11 = FADD(tf[2*i], r14);
  r12 = FSUB(r11, tf[2*i]);
  r12 = FSUB(r14, r12);
  r12 = FADD(r12, r15);
  r12 = FADD(r12, tf[2*i+1]);
  r12 = FFMA(tf[2*i+1], r9, r12);
  r12 = FFMA(tf[2*i], r13, r12);
  const float eps = 0x1p-40f;
  float r1 = FADD(r11, FADD(r12, eps)), r2 = FADD(r11, FSUB(r12, eps));
  if (!(r1 == r2)) return 0;
  *res = FMUL(r1, u2f((uint32_t)(e + 127) << 23));
  return 1;
}
#define NGP 5
(main) root@C.61801330:/workspace$ 