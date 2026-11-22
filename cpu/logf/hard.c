// Entrées "difficiles" de logf sur le domaine du chemin rapide : log(x) très proche de la grille demi-ulp binary32.
// Filtre rapide en double (exp glibc), confirmation MPFR 160 bits. Résultats de référence : CORE-MATH, 4 modes.
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <string.h>
#include <fenv.h>
#include <mpfr.h>
#include <omp.h>
float cr_logf(float);
static inline int in_fast(uint32_t u){ return !(u < 1u<<23 || u >= 0xffu<<23) && u != 0x3f800000u; }
static double grid_dist_ulp_d(double y){
  int e; frexp(y,&e); int ue = e-24; if(ue < -149) ue = -149;
  double q = ldexp(y, -(ue-1)); double f = q - floor(q); double d = f<0.5? f : 1-f; return d/2; }
int main(int argc, char** argv){
  const double TH_HARD = argc>1 ? atof(argv[1]) : 0x1p-22, TH_FILTER = 0x1p-16;
  long cnt=0, nd=0, nm=0;
  FILE *out = fopen("hard.txt","w");
  #pragma omp parallel reduction(+:cnt,nd,nm)
  {
    mpfr_t X, Y, Q, F; mpfr_inits2(160, X, Y, Q, F, (mpfr_ptr)0);
    #pragma omp for schedule(dynamic,65536)
    for (long i=0;i<(1L<<32);i++){
      uint32_t u=(uint32_t)i; if(!in_fast(u)) continue; nd++;
      float x; memcpy(&x,&u,4);
      if (grid_dist_ulp_d(fabs(log((double)x))) > TH_FILTER) continue;
      nm++;
      mpfr_set_flt(X, x, MPFR_RNDN); mpfr_log(Y, X, MPFR_RNDN); mpfr_abs(Y, Y, MPFR_RNDN);
      long e = mpfr_get_exp(Y); long ue = e-24; if (ue < -149) ue = -149;
      mpfr_mul_2si(Q, Y, -(ue-1), MPFR_RNDN); mpfr_frac(F, Q, MPFR_RNDN);
      double f = mpfr_get_d(F, MPFR_RNDN); double d = (f<0.5? f:1-f)/2;
      if (d <= TH_HARD){
        float o[4]; int modes[4]={FE_TONEAREST,FE_TOWARDZERO,FE_UPWARD,FE_DOWNWARD};
        for(int m=0;m<4;m++){ fesetround(modes[m]); o[m]=cr_logf(x); } fesetround(FE_TONEAREST);
        uint32_t ob[4]; memcpy(ob,o,16);
        #pragma omp critical
        fprintf(out, "%08x %08x %08x %08x %08x %.3e\n", u, ob[0],ob[1],ob[2],ob[3], d);
        cnt++;
      }
    }
    mpfr_clears(X,Y,Q,F,(mpfr_ptr)0);
  }
  fclose(out);
  printf("entrées du chemin rapide : %ld ; passées au filtre : %ld ; difficiles (<= %.3g ulp) : %ld\n", nd, nm, TH_HARD, cnt);
  return 0;
}
