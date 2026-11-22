// Chronométrage et re-vérification (code compilé) des chemins rapides trouvés par c3b
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cuda_runtime.h>
#include <algorithm>
#include "best_progs.cuh"
#include "hand_progs.cuh"
#define CK(c) do { cudaError_t e_ = (c); if (e_ != cudaSuccess) { printf("ERREUR CUDA %s ligne %d\n", cudaGetErrorString(e_), __LINE__); exit(1); } } while (0)
__device__ __noinline__ float slow_path(float x, const double* tb){ return cr2_exp2f(x, tb); }
template<int K> __device__ __forceinline__ int gpf(float x, const float* tf, float* r);
#define SPEC(k) template<> __device__ __forceinline__ int gpf<k>(float x, const float* tf, float* r){ return gp_fast_##k(x, tf, r); }
#if NGP > 0
SPEC(0)
#endif
#if NGP > 1
SPEC(1)
#endif
#if NGP > 2
SPEC(2)
#endif
#if NGP > 3
SPEC(3)
#endif
#if NGP > 4
SPEC(4)
#endif
#if NGP > 5
SPEC(5)
#endif
// F = -2 : CUDA ; F = -1 : graine ; F = -3 : main (cœur à 15) ; F = -4 : main (toutes astuces) ; F >= 0 : candidat génétique F
template<int F> __device__ __forceinline__ float f_(float x, const float* tf, const double* tb){
  if (F == -2) return exp2f(x);
  float r; int ok;
  if (F == -1) ok = fexp2_fast(x, tf, &r); else if (F == -3) ok = hand15_fast(x, tf, &r); else if (F == -4) ok = handall_fast(x, tf, &r); else ok = gpf<(F < 0 ? 0 : F)>(x, tf, &r);
  return ok ? r : slow_path(x, tb);
}
template<int F> __global__ void kthr(const double* tbg, const float* tfg, uint32_t* out, int iters){
  __shared__ double tb[64]; __shared__ __align__(16) float tf[256];
  if (threadIdx.x < 64) tb[threadIdx.x] = tbg[threadIdx.x]; if (threadIdx.x < 256) tf[threadIdx.x] = tfg[threadIdx.x]; __syncthreads();
  uint32_t idx = blockIdx.x * blockDim.x + threadIdx.x, s = idx * 2654435761u + 12345u, acc = 0;
  for (int k = 0; k < iters; k++){ s = s * 1664525u + 1013904223u; float x = fmaf((float)(s >> 8), 20.f * 0x1p-24f, -10.f); acc ^= f2u(f_<F>(x, tf, tb)); }
  out[idx] = acc;
}
template<int F> __global__ void kmem(const double* tbg, const float* tfg, const float* X, float* Y, uint32_t n){
  __shared__ double tb[64]; __shared__ __align__(16) float tf[256];
  if (threadIdx.x < 64) tb[threadIdx.x] = tbg[threadIdx.x]; if (threadIdx.x < 256) tf[threadIdx.x] = tfg[threadIdx.x]; __syncthreads();
  for (uint32_t i = blockIdx.x * blockDim.x + threadIdx.x; i < n; i += gridDim.x * blockDim.x){
    float x = X[i], acc = 0.f;
    #pragma unroll
    for (int j = 0; j < 16; j++) acc += f_<F>(x + j * 0x1p-10f, tf, tb);
    Y[i] = acc; }
}
template<int F> __global__ void kver(const double* tbg, const float* tfg, uint64_t n, unsigned long long* bad){
  __shared__ double tb[64]; __shared__ __align__(16) float tf[256];
  if (threadIdx.x < 64) tb[threadIdx.x] = tbg[threadIdx.x]; if (threadIdx.x < 256) tf[threadIdx.x] = tfg[threadIdx.x]; __syncthreads();
  unsigned long long b = 0;
  for (uint64_t i = blockIdx.x * (uint64_t)blockDim.x + threadIdx.x; i < n; i += (uint64_t)gridDim.x * blockDim.x){
    float x = u2f((uint32_t)i), a = f_<F>(x, tf, tb), r = cr2_exp2f(x, tb); if (f2u(a) != f2u(r) && !(a != a && r != r)) b++; }
  if (b) atomicAdd(bad, b);
}
__global__ void kinit(float* X, uint32_t n){ for (uint32_t i = blockIdx.x * blockDim.x + threadIdx.x; i < n; i += gridDim.x * blockDim.x){ uint32_t s = i * 2654435761u + 777u; s ^= s >> 13; s *= 1664525u; X[i] = -10.f + 20.f * (float)(s >> 8) * 0x1p-24f; } }
static cudaEvent_t e0, e1;
template<typename L> static double tm(L f){ double b = 1e9; float ms; for (int r = 0; r < 7; r++){ cudaEventRecord(e0); f(); cudaEventRecord(e1); CK(cudaEventSynchronize(e1)); cudaEventElapsedTime(&ms, e0, e1); if (ms < b) b = ms; } return b; }
double *dtb; float *dtf, *X, *Y; uint32_t* o; unsigned long long* dbad; int nsm; const uint32_t NM = 1u << 26;
template<int F> static double row(const char* name, double tc, double mc){
  int blocks = nsm * 16, th = 256, it = 2048;
  double t = tm([&]{ kthr<F><<<blocks, th>>>(dtb, dtf, o, it); });
  double m = tm([&]{ kmem<F><<<nsm * 8, 256>>>(dtb, dtf, X, Y, NM); });
  unsigned long long bad = 0; CK(cudaMemset(dbad, 0, 8)); kver<F><<<nsm * 32, 256>>>(dtb, dtf, 1ull << 32, dbad); CK(cudaMemcpy(&bad, dbad, 8, cudaMemcpyDeviceToHost));
  printf("  %-40s calcul pur x%5.2f   16 appels/nombre x%5.2f   erreurs sur 2^32 : %llu\n", name, tc > 0 ? t / tc : 1.0, mc > 0 ? m / mc : 1.0, bad);
  return bad ? 1e9 : t;
}
int main(){
  cudaDeviceProp p; CK(cudaGetDeviceProperties(&p, 0)); nsm = p.multiProcessorCount; cudaEventCreate(&e0); cudaEventCreate(&e1);
  CK(cudaMalloc(&dtb, 512)); CK(cudaMemcpy(dtb, CR2_TB, 512, cudaMemcpyHostToDevice)); CK(cudaMalloc(&dtf, 1024)); CK(cudaMemcpy(dtf, TF, 1024, cudaMemcpyHostToDevice));
  CK(cudaMalloc(&X, NM * 4ull)); CK(cudaMalloc(&Y, NM * 4ull)); CK(cudaMalloc(&o, nsm * 16 * 256 * 4)); CK(cudaMalloc(&dbad, 8));
  kinit<<<nsm * 8, 256>>>(X, NM); CK(cudaDeviceSynchronize());
  int blocks = nsm * 16, th = 256, it = 2048;
  double tc = tm([&]{ kthr<-2><<<blocks, th>>>(dtb, dtf, o, it); });
  double mc = tm([&]{ kmem<-2><<<nsm * 8, 256>>>(dtb, dtf, X, Y, NM); });
  printf("\nTemps relatifs à l'exp2f de CUDA (x1 = même vitesse) ; erreurs comptées contre la version exacte en double, sur la carte\n");
  printf("  %-36s calcul pur x%5.2f   16 appels/nombre x%5.2f   (CUDA : faux sur 4,4 %% des entrées)\n", "exp2f de CUDA", 1.0, 1.0);
  row<-1>("graine (version de départ)", tc, mc);
  double thand = row<-3>("main : cœur à 15 opérations", tc, mc);
  double ta = row<-4>("main : toutes astuces connues", tc, mc);
  double best = 1e9;
#if NGP > 0
  best = std::min(best, row<0>("Koza : candidat 0 (le plus court)", tc, mc));
#endif
#if NGP > 1
  best = std::min(best, row<1>("Koza : candidat 1", tc, mc));
#endif
#if NGP > 2
  best = std::min(best, row<2>("Koza : candidat 2", tc, mc));
#endif
#if NGP > 3
  best = std::min(best, row<3>("Koza : candidat 3", tc, mc));
#endif
  printf("\nCRITÈRE fixé d'avance : succès si un candidat Koza exact est au moins 10 %% plus rapide que la meilleure version\n");
  printf("écrite à la main (cœur à 15 opérations). Bonus : x4,00 ou moins par rapport à CUDA.\n");
  printf("  meilleur Koza : %.1f %% plus rapide que la version main à 15 ; %.1f %% par rapport à « toutes astuces » ; x%.2f CUDA\n",
         100 * (1 - best / thand), 100 * (1 - best / ta), best / tc);
  printf("  VERDICT : %s%s\n", best <= 0.90 * thand ? "SUCCÈS" : "NON ATTEINT", best / tc <= 4.0 ? " + BONUS (x4 atteint)" : "");
  printf("FIN\n");
  return 0;
}
