// ESSAI C3a : (1) coût en conditions réelles (tableau en mémoire) ; (2) diagnostic : où part le temps ?
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cmath>
#include <omp.h>
#include <cuda_runtime.h>
#include "fexp2.h"
extern "C" float cr_exp2f(float);
#define CK(c) do { cudaError_t e_ = (c); if (e_ != cudaSuccess) { printf("ERREUR CUDA %s ligne %d\n", cudaGetErrorString(e_), __LINE__); exit(1); } } while (0)

// chemin lent séparé (non intégré) : allège les registres du chemin rapide
__device__ __noinline__ float slow_path(float x, const double* tb){ return cr2_exp2f(x, tb); }

// V : 0 = complet (actuel) ; 1 = chemin lent non intégré (exact) ; 2 = sans chemin lent (diagnostic, FAUX)
//     3 = sans chemin lent ni test d'arrondi (diagnostic) ; 4 = idem + table remplacée par une constante (diagnostic)
template<int V>
__device__ __forceinline__ float fvar(float x, const float* tf, const double* tb){
  if (V == 0) return fexp2f(x, tf, tb);
  float r;
  if (V == 1){ if (fexp2_fast(x, tf, &r)) return r; return slow_path(x, tb); }
  uint32_t u = f2u(x), ux = u << 1;
  if (ux < 0x65000000u || ux >= 0x85fffc00u || u > 0xc2fa0000u) return __int_as_float(0x7fc00000);
  const float MAGIC = 0x1.8p23f;
  float t = FMUL(x, 128.0f), tk = FADD(t, MAGIC);
  int32_t k = (int32_t)(f2u(tk) - f2u(MAGIC));
  float kf = FSUB(tk, MAGIC), h = FMUL(FSUB(t, kf), 0x1p-7f);
  int i = k & 127, e = k >> 7;
  float th, tl;
  if (V == 4){ th = FFMA((float)i, 0x1p-8f, 1.0f); tl = FMUL((float)i, 0x1p-30f); } else { th = tf[2*i]; tl = tf[2*i+1]; }
  float ph = FMUL(h, C1H), pe = FFMA(h, C1H, -ph), h2 = FMUL(h, h);
  float w = FFMA(h, C4, C3); w = FFMA(h, w, C2);
  float qlo = FFMA(h, C1L, FFMA(h2, w, pe));
  float mh = FMUL(th, ph), ml = FFMA(th, ph, -mh), sh = FADD(th, mh), sl = FSUB(mh, FSUB(sh, th));
  float L = FADD(sl, ml); L = FADD(L, tl); L = FFMA(tl, ph, L); L = FFMA(th, qlo, L);
  float res;
  if (V == 2){ float eps = FMUL(sh, FEXP2_E); float r1 = FADD(sh, FADD(L, eps)), r2 = FADD(sh, FSUB(L, eps));
               res = (r1 == r2) ? r1 : __int_as_float(0x7fc00000); }
  else res = FADD(sh, L);
  return FMUL(res, __uint_as_float((uint32_t)(e + 127) << 23));
}

template<int F>
__global__ void kthr(const double* __restrict__ tbg, const float* __restrict__ tfg, uint32_t* out, int iters, float lo, float span){
  __shared__ double tb[64];
  __shared__ __align__(16) float tf[256];
  if (threadIdx.x < 64) tb[threadIdx.x] = tbg[threadIdx.x];
  if (threadIdx.x < 256) tf[threadIdx.x] = tfg[threadIdx.x];
  __syncthreads();
  uint32_t idx = blockIdx.x * blockDim.x + threadIdx.x;
  uint32_t s = idx * 2654435761u + 12345u, acc = 0;
  for (int k = 0; k < iters; k++){
    s = s * 1664525u + 1013904223u;
    float x = fmaf((float)(s >> 8), span * 0x1p-24f, lo);
    float y;
    if (F == 9) y = exp2f(x); else if (F == 8) y = x; else y = fvar<F>(x, tf, tb);
    acc ^= f2u(y);
  }
  out[idx] = acc;
}

// conditions réelles : y[i] = somme_{j<D} f(x[i] + j/1024), tableau de 2^28 float en mémoire
template<int F, int D>
__global__ void kmem(const double* __restrict__ tbg, const float* __restrict__ tfg, const float* __restrict__ X, float* Y, uint32_t n){
  __shared__ double tb[64];
  __shared__ __align__(16) float tf[256];
  if (threadIdx.x < 64) tb[threadIdx.x] = tbg[threadIdx.x];
  if (threadIdx.x < 256) tf[threadIdx.x] = tfg[threadIdx.x];
  __syncthreads();
  for (uint32_t i = blockIdx.x * blockDim.x + threadIdx.x; i < n; i += gridDim.x * blockDim.x){
    float x = X[i], acc = 0.f;
    #pragma unroll
    for (int j = 0; j < D; j++){ float xj = x + j * 0x1p-10f; acc += (F == 9) ? exp2f(xj) : (F == 3 ? cr2_exp2f(xj, tb) : fvar<F>(xj, tf, tb)); }
    Y[i] = acc;
  }
}
__global__ void kinit(float* X, uint32_t n){
  for (uint32_t i = blockIdx.x * blockDim.x + threadIdx.x; i < n; i += gridDim.x * blockDim.x){
    uint32_t s = i * 2654435761u + 777u; s ^= s >> 13; s *= 1664525u; X[i] = -10.f + 20.f * (float)(s >> 8) * 0x1p-24f; }
}
__global__ void kcheck1(const double* __restrict__ tbg, const float* __restrict__ tfg, uint32_t base, uint32_t n, uint32_t* o){
  __shared__ double tb[64];
  __shared__ __align__(16) float tf[256];
  if (threadIdx.x < 64) tb[threadIdx.x] = tbg[threadIdx.x];
  if (threadIdx.x < 256) tf[threadIdx.x] = tfg[threadIdx.x];
  __syncthreads();
  uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < n) o[i] = f2u(fvar<1>(u2f(base + i), tf, tb));
}

static cudaEvent_t e0, e1;
template<typename K> static double timeit(K launch, int reps = 5){ double best = 1e9; float ms;
  for (int r = 0; r < reps; r++){ cudaEventRecord(e0); launch(); cudaEventRecord(e1); CK(cudaEventSynchronize(e1)); cudaEventElapsedTime(&ms, e0, e1); if (ms < best) best = ms; } return best; }
template<typename T> static int regs(T k){ cudaFuncAttributes a; cudaFuncGetAttributes(&a, k); return a.numRegs; }

int main(){
  cudaDeviceProp p; CK(cudaGetDeviceProperties(&p, 0)); int nsm = p.multiProcessorCount;
  printf("Carte : %s, %d SM\n", p.name, nsm);
  cudaEventCreate(&e0); cudaEventCreate(&e1);
  double* dtb; float* dtf;
  CK(cudaMalloc(&dtb, 64 * 8)); CK(cudaMemcpy(dtb, CR2_TB, 64 * 8, cudaMemcpyHostToDevice));
  CK(cudaMalloc(&dtf, 256 * 4)); CK(cudaMemcpy(dtf, TF, 256 * 4, cudaMemcpyHostToDevice));

  // ---------- (1) conditions réelles
  { const uint32_t n = 1u << 28; float *X, *Y; CK(cudaMalloc(&X, n * 4ull)); CK(cudaMalloc(&Y, n * 4ull));
    kinit<<<nsm * 8, 256>>>(X, n); CK(cudaDeviceSynchronize());
    int G = nsm * 8, T = 256;
    printf("\n(1) CONDITIONS RÉELLES : y[i] = somme de D appels exp2f, tableau de 268 millions de float (1 Go lu + 1 Go écrit)\n");
    printf("   %-34s %14s %14s %14s\n", "", "D = 1", "D = 4", "D = 16");
    const char* nm[4] = {"exp2f de CUDA (faux sur 4,4 %)", "EXACT float (actuel)", "EXACT float, chemin lent séparé", "EXACT double"};
    double ref[3] = {0, 0, 0};
    for (int f = 0; f < 4; f++){
      double t1, t4, t16;
      if (f == 0){ t1 = timeit([&]{ kmem<9,1><<<G,T>>>(dtb, dtf, X, Y, n); }); t4 = timeit([&]{ kmem<9,4><<<G,T>>>(dtb, dtf, X, Y, n); }); t16 = timeit([&]{ kmem<9,16><<<G,T>>>(dtb, dtf, X, Y, n); }); ref[0] = t1; ref[1] = t4; ref[2] = t16; }
      if (f == 1){ t1 = timeit([&]{ kmem<0,1><<<G,T>>>(dtb, dtf, X, Y, n); }); t4 = timeit([&]{ kmem<0,4><<<G,T>>>(dtb, dtf, X, Y, n); }); t16 = timeit([&]{ kmem<0,16><<<G,T>>>(dtb, dtf, X, Y, n); }); }
      if (f == 2){ t1 = timeit([&]{ kmem<1,1><<<G,T>>>(dtb, dtf, X, Y, n); }); t4 = timeit([&]{ kmem<1,4><<<G,T>>>(dtb, dtf, X, Y, n); }); t16 = timeit([&]{ kmem<1,16><<<G,T>>>(dtb, dtf, X, Y, n); }); }
      if (f == 3){ t1 = timeit([&]{ kmem<3,1><<<G,T>>>(dtb, dtf, X, Y, n); }); t4 = timeit([&]{ kmem<3,4><<<G,T>>>(dtb, dtf, X, Y, n); }); t16 = timeit([&]{ kmem<3,16><<<G,T>>>(dtb, dtf, X, Y, n); }); }
      printf("   %-34s %7.2f ms x%4.2f %7.2f ms x%4.2f %7.2f ms x%4.2f\n", nm[f], t1, t1 / ref[0], t4, t4 / ref[1], t16, t16 / ref[2]);
    }
    printf("   (x = temps relatif à CUDA ; D = nombre d'exp2f par nombre lu en mémoire)\n");
    cudaFree(X); cudaFree(Y); }

  // ---------- (2) diagnostic en calcul pur
  { int blocks = nsm * 16, th = 256, it = 2048; uint32_t* o; CK(cudaMalloc(&o, blocks * th * 4));
    double c = (double)blocks * th * it;
    double tc = timeit([&]{ kthr<9><<<blocks, th>>>(dtb, dtf, o, it, -10.f, 20.f); });
    double t0 = timeit([&]{ kthr<0><<<blocks, th>>>(dtb, dtf, o, it, -10.f, 20.f); });
    double t1 = timeit([&]{ kthr<1><<<blocks, th>>>(dtb, dtf, o, it, -10.f, 20.f); });
    double t2 = timeit([&]{ kthr<2><<<blocks, th>>>(dtb, dtf, o, it, -10.f, 20.f); });
    double t3 = timeit([&]{ kthr<3><<<blocks, th>>>(dtb, dtf, o, it, -10.f, 20.f); });
    double t4 = timeit([&]{ kthr<4><<<blocks, th>>>(dtb, dtf, o, it, -10.f, 20.f); });
    double tv = timeit([&]{ kthr<8><<<blocks, th>>>(dtb, dtf, o, it, -10.f, 20.f); });
    printf("\n(2) DIAGNOSTIC en calcul pur (entrées dans [-10,10]) : débit, temps relatif à CUDA, registres par fil\n");
    printf("   %-52s %8.1f G/s  x%5.2f  %3d reg\n", "exp2f de CUDA", c / tc / 1e6, 1.0, regs(kthr<9>));
    printf("   %-52s %8.1f G/s  x%5.2f  %3d reg\n", "EXACT float actuel", c / t0 / 1e6, t0 / tc, regs(kthr<0>));
    printf("   %-52s %8.1f G/s  x%5.2f  %3d reg\n", "EXACT float, chemin lent séparé (exact aussi)", c / t1 / 1e6, t1 / tc, regs(kthr<1>));
    printf("   %-52s %8.1f G/s  x%5.2f  %3d reg\n", "  diagnostic : sans chemin lent du tout", c / t2 / 1e6, t2 / tc, regs(kthr<2>));
    printf("   %-52s %8.1f G/s  x%5.2f  %3d reg\n", "  diagnostic : sans chemin lent ni test d'arrondi", c / t3 / 1e6, t3 / tc, regs(kthr<3>));
    printf("   %-52s %8.1f G/s  x%5.2f  %3d reg\n", "  diagnostic : idem, sans lecture de table", c / t4 / 1e6, t4 / tc, regs(kthr<4>));
    printf("   %-52s %8.1f G/s  x%5.2f  %3d reg\n", "  référence vide (génération des x)", c / tv / 1e6, tv / tc, regs(kthr<8>));
    cudaFree(o); }

  // ---------- (3) exactitude de la variante « chemin lent séparé »
  { const uint32_t n = 1u << 28; uint32_t* d; CK(cudaMalloc(&d, n * 4ull)); uint32_t* h = (uint32_t*)malloc(n * 4ull); long bad = 0;
    for (uint32_t cc = 0; cc < 16; cc++){ uint32_t base = cc << 28;
      kcheck1<<<n / 256, 256>>>(dtb, dtf, base, n, d); CK(cudaGetLastError()); CK(cudaMemcpy(h, d, n * 4ull, cudaMemcpyDeviceToHost));
      #pragma omp parallel for reduction(+:bad) schedule(static, 1 << 16)
      for (long i = 0; i < (long)n; i++){ float x = u2f(base + (uint32_t)i), rf = cr_exp2f(x), g = u2f(h[i]); if (!(rf != rf && g != g) && h[i] != f2u(rf)) bad++; } }
    printf("\n(3) Exactitude « chemin lent séparé » sur les 4 294 967 296 entrées : %ld erreurs (doit être 0)\n", bad);
    free(h); cudaFree(d); }
  printf("\nFIN\n");
  return 0;
}
