// Chronos réels : débit (appels indépendants) et latence (appels enchaînés) de plusieurs expf
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
float cr_expf(float); float new_expf(float); float s2_expf(float); float glibc_expf(float);
static float sys_expf(float x){ return expf(x); }
typedef float (*fn)(float);
static double now(void){ struct timespec t; clock_gettime(CLOCK_MONOTONIC, &t); return t.tv_sec + 1e-9*t.tv_nsec; }
#define N 4096
static float X[N], Y[N];
static double thr(fn f, int R){ double t0 = now(); for (int r = 0; r < R; r++) for (int i = 0; i < N; i++) Y[i] = f(X[i]); return (now() - t0) / ((double)N*R) * 1e9; }
static double lat(fn f, int R){ float acc = 0; double t0 = now(); for (int r = 0; r < R; r++) for (int i = 0; i < N; i++) acc = f(X[i] + acc*0.0f); Y[0] = acc; return (now() - t0) / ((double)N*R) * 1e9; }
int main(void){
  const char* names[] = {"glibc système (non exacte)", "clone glibc (non exact)", "CORE-MATH actuel (exact)", "nouveau, réduction 3 op. (exact)", "nouveau, réduction 4 op. (exact)"};
  fn fs[] = {sys_expf, glibc_expf, cr_expf, new_expf, s2_expf};
  int nf = 5;
  double lo[3] = {-10, -1, -80}, hi[3] = {10, 1, 80};
  for (int dset = 0; dset < 3; dset++){
    srand(42); for (int i = 0; i < N; i++) X[i] = lo[dset] + (hi[dset]-lo[dset]) * (rand() / (double)RAND_MAX);
    double bt[5], bl[5]; for (int k = 0; k < nf; k++){ bt[k] = bl[k] = 1e9; }
    for (int rep = 0; rep < 25; rep++) for (int k = 0; k < nf; k++){
      double a = thr(fs[k], 200), b = lat(fs[k], 100); if (a < bt[k]) bt[k] = a; if (b < bl[k]) bl[k] = b; }
    printf("\nEntrées uniformes dans [%g, %g] (meilleur de 25 essais)\n", lo[dset], hi[dset]);
    printf("  %-38s %9s %11s %14s %14s\n", "fonction", "débit ns", "latence ns", "débit vs CM", "latence vs CM");
    for (int k = 0; k < nf; k++) printf("  %-38s %9.3f %11.3f %+13.1f%% %+13.1f%%\n", names[k], bt[k], bl[k], 100*(bt[k]/bt[2]-1), 100*(bl[k]/bl[2]-1));
  }
  return 0;
}
