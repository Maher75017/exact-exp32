// ESSAI C1 (GPU) : quel écart de vitesse entre l'exp2f de CUDA (rapide, pas toujours exact)
// et un exp2f exact (calcul en double) sur cette carte ? Et combien de résultats CUDA sont faux ?
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <cmath>
#include <omp.h>
#include <cuda_runtime.h>
#include "cr2.h"
extern "C" float cr_exp2f(float);          // CORE-MATH (référence CPU, prouvée exacte)

#define CK(c) do { cudaError_t e = (c); if (e != cudaSuccess) { printf("ERREUR CUDA %s ligne %d\n", cudaGetErrorString(e), __LINE__); exit(1); } } while (0)

__device__ __forceinline__ float ex2raw(float x){ float r; asm("ex2.approx.f32 %0, %1;" : "=f"(r) : "f"(x)); return r; }

template<int F>
__global__ void kthr(const double* __restrict__ tbg, uint32_t* out, int iters, float lo, float span){
  __shared__ double tb[64];
  if (threadIdx.x < 64) tb[threadIdx.x] = tbg[threadIdx.x];
  __syncthreads();
  uint32_t idx = blockIdx.x * blockDim.x + threadIdx.x;
  uint32_t s = idx * 2654435761u + 12345u, acc = 0;
  for (int k = 0; k < iters; k++){
    s = s * 1664525u + 1013904223u;
    float x = fmaf((float)(s >> 8), span * 0x1p-24f, lo);
    float y;
    if (F == 0) y = exp2f(x);
    else if (F == 1) y = ex2raw(x);
    else if (F == 2) y = cr2_exp2f(x, tb);
    else y = x;
    acc ^= f2u(y);
  }
  out[idx] = acc;
}

template<typename T>
__global__ void kfma(T* out, int iters){
  T a = (T)threadIdx.x * (T)1e-3, b = a + (T)1, c = a + (T)2, d = a + (T)3, m = (T)0.999, q = (T)1e-4;
  for (int k = 0; k < iters; k++){ a = a*m + q; b = b*m + q; c = c*m + q; d = d*m + q; }
  out[blockIdx.x * blockDim.x + threadIdx.x] = a + b + c + d;
}

__global__ void kcheck(const double* __restrict__ tbg, uint32_t base, uint32_t n, uint32_t* oc, uint32_t* orr){
  __shared__ double tb[64];
  if (threadIdx.x < 64) tb[threadIdx.x] = tbg[threadIdx.x];
  __syncthreads();
  uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= n) return;
  float x = u2f(base + i);
  oc[i] = f2u(exp2f(x));
  orr[i] = f2u(cr2_exp2f(x, tb));
}

int main(){
  cudaDeviceProp p; CK(cudaGetDeviceProperties(&p, 0));
  int nsm = p.multiProcessorCount;
  printf("Carte : %s, %d SM\n", p.name, nsm);
  double* dtb; CK(cudaMalloc(&dtb, 64 * sizeof(double))); CK(cudaMemcpy(dtb, CR2_TB, 64 * sizeof(double), cudaMemcpyHostToDevice));
  cudaEvent_t e0, e1; cudaEventCreate(&e0); cudaEventCreate(&e1);
  float ms;

  // --- 0) rapport float / double sur cette carte
  { int blocks = nsm * 16, th = 256, it = 4096; float* of; double* od;
    CK(cudaMalloc(&of, blocks * th * sizeof(float))); CK(cudaMalloc(&od, blocks * th * sizeof(double)));
    double bf = 1e9, bd = 1e9;
    for (int r = 0; r < 5; r++){
      cudaEventRecord(e0); kfma<float><<<blocks, th>>>(of, it); cudaEventRecord(e1); CK(cudaEventSynchronize(e1)); cudaEventElapsedTime(&ms, e0, e1); if (ms < bf) bf = ms;
      cudaEventRecord(e0); kfma<double><<<blocks, th>>>(od, it / 16); cudaEventRecord(e1); CK(cudaEventSynchronize(e1)); cudaEventElapsedTime(&ms, e0, e1); if (ms < bd) bd = ms;
    }
    double nf = 4.0 * 2 * blocks * th * (double)it, nd = 4.0 * 2 * blocks * th * (double)(it / 16);
    printf("\n0) Puissance brute : float %.1f TFLOPS, double %.2f TFLOPS (rapport %.0f)\n", nf / bf / 1e9, nd / bd / 1e9, (nf / bf) / (nd / bd));
    cudaFree(of); cudaFree(od); }

  // --- 1) débit des fonctions (entrées générées dans les registres, pas de mémoire)
  { int blocks = nsm * 16, th = 256, it = 2048; uint32_t* o; CK(cudaMalloc(&o, blocks * th * 4));
    const char* nm[4] = {"exp2f de CUDA (pas toujours exact)", "instruction brute ex2.approx", "exp2f EXACT (calcul en double)", "référence vide (génération des x)"};
    float lo[2] = {-10.f, -1.f}, sp[2] = {20.f, 2.f};
    for (int dset = 0; dset < 2; dset++){
      double best[4];
      for (int f = 0; f < 4; f++){ best[f] = 1e9;
        for (int r = 0; r < 5; r++){
          cudaEventRecord(e0);
          if (f == 0) kthr<0><<<blocks, th>>>(dtb, o, it, lo[dset], sp[dset]);
          if (f == 1) kthr<1><<<blocks, th>>>(dtb, o, it, lo[dset], sp[dset]);
          if (f == 2) kthr<2><<<blocks, th>>>(dtb, o, it, lo[dset], sp[dset]);
          if (f == 3) kthr<3><<<blocks, th>>>(dtb, o, it, lo[dset], sp[dset]);
          cudaEventRecord(e1); CK(cudaEventSynchronize(e1)); cudaEventElapsedTime(&ms, e0, e1); if (ms < best[f]) best[f] = ms; } }
      double calls = (double)blocks * th * it;
      printf("\n1) Débit, entrées dans [%g, %g] (milliards d'appels par seconde ; « net » = sans le coût de génération des x)\n", lo[dset], lo[dset] + sp[dset]);
      for (int f = 0; f < 4; f++){
        double net = best[f] - best[3];
        printf("   %-40s %8.1f G/s", nm[f], calls / best[f] / 1e6);
        if (f < 3 && net > 0) printf("   net %8.1f G/s", calls / net / 1e6);
        printf("\n"); }
      double r = (best[2] - best[3]) / (best[0] - best[3]);
      printf("   => la version exacte en double est %.1f fois plus lente que l'exp2f de CUDA (net)\n", r);
    }
    cudaFree(o); }

  // --- 2) exactitude sur les 2^32 entrées (arrondi au plus proche)
  { const uint32_t n = 1u << 28; uint32_t *dc, *dr; CK(cudaMalloc(&dc, n * 4ull)); CK(cudaMalloc(&dr, n * 4ull));
    uint32_t* hc = (uint32_t*)malloc(n * 4ull); uint32_t* hr = (uint32_t*)malloc(n * 4ull);
    long bad_port = 0, wrong_cuda = 0; uint32_t maxd = 0;
    for (uint32_t c = 0; c < 16; c++){
      uint32_t base = c << 28;
      kcheck<<<n / 256, 256>>>(dtb, base, n, dc, dr); CK(cudaGetLastError());
      CK(cudaMemcpy(hc, dc, n * 4ull, cudaMemcpyDeviceToHost)); CK(cudaMemcpy(hr, dr, n * 4ull, cudaMemcpyDeviceToHost));
      #pragma omp parallel for reduction(+:bad_port,wrong_cuda) reduction(max:maxd) schedule(static, 1 << 16)
      for (long i = 0; i < (long)n; i++){
        uint32_t u = base + (uint32_t)i; float x = u2f(u); float rf = cr_exp2f(x); uint32_t rb = f2u(rf);
        bool rn = (rf != rf);
        float gf = u2f(hr[i]), cf = u2f(hc[i]);
        if (!(rn && gf != gf) && hr[i] != rb) bad_port++;
        if (!(rn && cf != cf) && hc[i] != rb){ wrong_cuda++;
          if (!rn && ((hc[i] ^ rb) >> 31) == 0){ uint32_t dd = hc[i] > rb ? hc[i] - rb : rb - hc[i]; if (dd > maxd) maxd = dd; } }
      }
    }
    printf("\n2) Exactitude sur les 4 294 967 296 entrées (arrondi au plus proche), référence CORE-MATH sur CPU\n");
    printf("   exp2f EXACT sur GPU : %ld différences (doit être 0)\n", bad_port);
    printf("   exp2f de CUDA       : %ld résultats faux (%.2f %%), écart maximal %u ulp\n", wrong_cuda, 100.0 * wrong_cuda / 4294967296.0, maxd);
    free(hc); free(hr); cudaFree(dc); cudaFree(dr); }
  printf("\nFIN\n");
  return 0;
}
