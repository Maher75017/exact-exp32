// Chronos réels de logf : débit (appels indépendants) et latence (appels enchaînés)
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
float cr_logf(float); float new_logf(float);
static float sys_logf(float x){ return logf(x); }
typedef float (*fn)(float);
static double now(void){ struct timespec t; clock_gettime(CLOCK_MONOTONIC, &t); return t.tv_sec + 1e-9*t.tv_nsec; }
#define N 4096
static float X[N], Y[N];
static double thr(fn f, int R){ double t0 = now(); for (int r = 0; r < R; r++) for (int i = 0; i < N; i++) Y[i] = f(X[i]); return (now() - t0) / ((double)N*R) * 1e9; }
static double lat(fn f, int R){ float acc = 0; double t0 = now(); for (int r = 0; r < R; r++) for (int i = 0; i < N; i++) acc = f(X[i] + acc*0.0f); Y[0] = acc; return (now() - t0) / ((double)N*R) * 1e9; }
int main(void){
  const char* names[] = {"glibc système (non exacte)", "CORE-MATH actuel (exact)", "nouveau (exact)"};
  fn fs[] = {sys_logf, cr_logf, new_logf}; int nf = 3;
  for (int dset = 0; dset < 3; dset++){
    srand(42);
    for (int i = 0; i < N; i++){ float u = rand()/(float)RAND_MAX;
      X[i] = dset == 0 ? ldexpf(1 + u, rand()%40 - 20) : dset == 1 ? 0.5f + u : ldexpf(1 + u, rand()%250 - 125); }
    double bt[3], bl[3]; for (int k = 0; k < nf; k++) bt[k] = bl[k] = 1e9;
    for (int rep = 0; rep < 25; rep++) for (int k = 0; k < nf; k++){ double a = thr(fs[k], 200), b = lat(fs[k], 100); if (a < bt[k]) bt[k] = a; if (b < bl[k]) bl[k] = b; }
    const char* dn[] = {"2^-20 à 2^20", "0,5 à 1,5", "2^-125 à 2^125"};
    printf("\nEntrées de %s (meilleur de 25 essais)\n", dn[dset]);
    printf("  %-30s %9s %11s %14s %14s\n", "fonction", "débit ns", "latence ns", "débit vs CM", "latence vs CM");
    for (int k = 0; k < nf; k++) printf("  %-30s %9.3f %11.3f %+13.1f%% %+13.1f%%\n", names[k], bt[k], bl[k], 100*(bt[k]/bt[1]-1), 100*(bl[k]/bl[1]-1));
  }
}
