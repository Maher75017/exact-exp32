// ESSAI C2 (GPU) : exp2f EXACT calculé en float (double-float) contre l'exp2f de CUDA
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cmath>
#include <omp.h>
#include <cuda_runtime.h>
#include "fexp2.h"
extern "C" float cr_exp2f(float);
#define CK(c) do { cudaError_t e_ = (c); if (e_ != cudaSuccess) { printf("ERREUR CUDA %s ligne %d\n", cudaGetErrorString(e_), __LINE__); exit(1); } } while (0)

__device__ unsigned long long g_slow;

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
    if (F == 0) y = exp2f(x);                     // CUDA
    else if (F == 1) y = fexp2f(x, tf, tb);       // exact, float, table en mémoire partagée
    else if (F == 2) y = fexp2f(x, tfg, tb);      // exact, float, table en mémoire globale (cache)
    else if (F == 3) y = cr2_exp2f(x, tb);        // exact, double
    else y = x;                                   // référence vide
    acc ^= f2u(y);
  }
  out[idx] = acc;
}

__global__ void kcheck(const double* __restrict__ tbg, const float* __restrict__ tfg, uint32_t base, uint32_t n, uint32_t* o){
  __shared__ double tb[64];
  __shared__ __align__(16) float tf[256];
  if (threadIdx.x < 64) tb[threadIdx.x] = tbg[threadIdx.x];
  if (threadIdx.x < 256) tf[threadIdx.x] = tfg[threadIdx.x];
  __syncthreads();
  uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= n) return;
  float x = u2f(base + i), r;
  if (!fexp2_fast(x, tf, &r)){ r = cr2_exp2f(x, tb); atomicAdd(&g_slow, 1ull); }
  o[i] = f2u(r);
}

int main(){
  cudaDeviceProp p; CK(cudaGetDeviceProperties(&p, 0));
  int nsm = p.multiProcessorCount;
  printf("Carte : %s, %d SM\n", p.name, nsm);
  double* dtb; float* dtf;
  CK(cudaMalloc(&dtb, 64 * 8)); CK(cudaMemcpy(dtb, CR2_TB, 64 * 8, cudaMemcpyHostToDevice));
  CK(cudaMalloc(&dtf, 256 * 4)); CK(cudaMemcpy(dtf, TF, 256 * 4, cudaMemcpyHostToDevice));
  cudaEvent_t e0, e1; cudaEventCreate(&e0); cudaEventCreate(&e1); float ms;

  { int blocks = nsm * 16, th = 256, it = 2048; uint32_t* o; CK(cudaMalloc(&o, blocks * th * 4));
    const char* nm[5] = {"exp2f de CUDA (faux sur 4,4 % des entrées)", "EXACT en float, table partagée", "EXACT en float, table globale", "EXACT en double (essai C1)", "référence vide"};
    float lo[2] = {-10.f, -1.f}, sp[2] = {20.f, 2.f};
    for (int dset = 0; dset < 2; dset++){
      double best[5];
      for (int f = 0; f < 5; f++){ best[f] = 1e9;
        for (int r = 0; r < 5; r++){
          cudaEventRecord(e0);
          if (f == 0) kthr<0><<<blocks, th>>>(dtb, dtf, o, it, lo[dset], sp[dset]);
          if (f == 1) kthr<1><<<blocks, th>>>(dtb, dtf, o, it, lo[dset], sp[dset]);
          if (f == 2) kthr<2><<<blocks, th>>>(dtb, dtf, o, it, lo[dset], sp[dset]);
          if (f == 3) kthr<3><<<blocks, th>>>(dtb, dtf, o, it, lo[dset], sp[dset]);
          if (f == 4) kthr<4><<<blocks, th>>>(dtb, dtf, o, it, lo[dset], sp[dset]);
          cudaEventRecord(e1); CK(cudaEventSynchronize(e1)); cudaEventElapsedTime(&ms, e0, e1); if (ms < best[f]) best[f] = ms; } }
      double calls = (double)blocks * th * it;
      printf("\nDébit, entrées dans [%g, %g] (milliards d'appels/s)\n", lo[dset], lo[dset] + sp[dset]);
      for (int f = 0; f < 5; f++){
        printf("   %-44s %8.1f G/s", nm[f], calls / best[f] / 1e6);
        if (f > 0 && f < 4) printf("   %5.2f fois le temps de CUDA (brut), %5.2f (net)", best[f] / best[0], (best[f] - best[4]) / (best[0] - best[4]));
        printf("\n"); }
    }
    cudaFree(o); }

  { const uint32_t n = 1u << 28; uint32_t* d; CK(cudaMalloc(&d, n * 4ull)); uint32_t* h = (uint32_t*)malloc(n * 4ull);
    unsigned long long zero = 0; CK(cudaMemcpyToSymbol(g_slow, &zero, sizeof(zero)));
    long bad = 0;
    for (uint32_t c = 0; c < 16; c++){
      uint32_t base = c << 28;
      kcheck<<<n / 256, 256>>>(dtb, dtf, base, n, d); CK(cudaGetLastError());
      CK(cudaMemcpy(h, d, n * 4ull, cudaMemcpyDeviceToHost));
      #pragma omp parallel for reduction(+:bad) schedule(static, 1 << 16)
      for (long i = 0; i < (long)n; i++){
        float x = u2f(base + (uint32_t)i); float rf = cr_exp2f(x); float g = u2f(h[i]);
        if (!(rf != rf && g != g) && h[i] != f2u(rf)) bad++; }
    }
    unsigned long long slow; CK(cudaMemcpyFromSymbol(&slow, g_slow, sizeof(slow)));
    printf("\nExactitude sur GPU, 4 294 967 296 entrées (arrondi au plus proche) : %ld erreurs (doit être 0)\n", bad);
    printf("Appels passés par le chemin double (cas limites + cas spéciaux) : %llu\n", slow);
    free(h); cudaFree(d); }
  printf("\nFIN\n");
  return 0;
}
