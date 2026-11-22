// Trouve les entrées "difficiles" de exp2f : 2^x très proche d'un point de la grille demi-ulp binary32.
// Filtre rapide en double (exp2 glibc, ~1 ulp double), confirmation MPFR 128 bits.
#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <string.h>
#include <fenv.h>
#include <mpfr.h>
#include <omp.h>
float cr_exp2f(float);
static inline int in_domain(uint32_t u){ uint32_t ux=u<<1;
  if (ux>=0x86000000u || ux<0x65000000u){ if(!(u>=0xc3000000u && u<0xc3150000u)) return 0; }
  return 1; }
// distance (en ulp binary32) de y à la grille demi-ulp ; y>0
static double grid_dist_ulp_d(double y){
  int e; frexp(y,&e);             // y = m*2^e, m in [0.5,1)
  int ue = e-24; if(ue < -149) ue = -149;   // ulp = 2^ue (sous-normaux inclus)
  double q = ldexp(y, -(ue-1));   // y / (ulp/2)
  double f = q - floor(q); double d = f<0.5? f : 1-f; return d/2; }
int main(void){
  const double TH_FILTER = 0x1p-16, TH_HARD = 0x1p-22;
  long cnt=0, nd=0, nm=0;
  FILE *out = fopen("hard.txt","w");
  #pragma omp parallel reduction(+:cnt,nd,nm)
  {
    mpfr_t X, Y, Q, F; mpfr_inits2(160, X, Y, Q, F, (mpfr_ptr)0);
    #pragma omp for schedule(dynamic,65536)
    for (long i=0;i<(1L<<32);i++){
      uint32_t u=(uint32_t)i; if(!in_domain(u)) continue; nd++;
      float x; memcpy(&x,&u,4);
      double y = exp2((double)x);
      if (grid_dist_ulp_d(y) > TH_FILTER) continue;
      nm++;
      mpfr_set_flt(X, x, MPFR_RNDN); mpfr_exp2(Y, X, MPFR_RNDN);
      long e = mpfr_get_exp(Y);          // Y = m*2^e, m in [0.5,1)
      long ue = e-24; if (ue < -149) ue = -149;
      mpfr_mul_2si(Q, Y, -(ue-1), MPFR_RNDN);
      mpfr_frac(F, Q, MPFR_RNDN);
      double f = mpfr_get_d(F, MPFR_RNDN); double d = (f<0.5? f:1-f)/2;
      if (d <= TH_HARD){
        float o[4]; int modes[4]={FE_TONEAREST,FE_TOWARDZERO,FE_UPWARD,FE_DOWNWARD};
        for(int m=0;m<4;m++){ fesetround(modes[m]); o[m]=cr_exp2f(x); } fesetround(FE_TONEAREST);
        uint32_t ob[4]; memcpy(ob,o,16);
        #pragma omp critical
        fprintf(out, "%08x %08x %08x %08x %08x %.3e\n", u, ob[0],ob[1],ob[2],ob[3], d);
        cnt++;
      }
    }
    mpfr_clears(X,Y,Q,F,(mpfr_ptr)0);
  }
  fclose(out);
  printf("entrées du domaine : %ld ; passées au filtre : %ld ; difficiles (< 2^-22 ulp) : %ld\n", nd, nm, cnt);
  return 0;
}
