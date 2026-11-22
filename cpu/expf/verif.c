// Comparaison exhaustive : new_expf (fichier propre) contre cr_expf (CORE-MATH actuel), 4 modes, 2^32 entrées
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <fenv.h>
float cr_expf(float); float new_expf(float);
int main(void){ int modes[4]={FE_TONEAREST,FE_TOWARDZERO,FE_UPWARD,FE_DOWNWARD}; const char* mn[4]={"RN","RZ","RU","RD"}; long tot=0;
  for(int m=0;m<4;m++){ long f=0;
    #pragma omp parallel reduction(+:f)
    { fesetround(modes[m]);
      #pragma omp for schedule(static,1<<16)
      for(long i=0;i<(1L<<32);i++){ uint32_t u=(uint32_t)i; float x; memcpy(&x,&u,4);
        float a=new_expf(x), b=cr_expf(x); uint32_t ab,bb; memcpy(&ab,&a,4); memcpy(&bb,&b,4);
        if(ab!=bb && !(a!=a && b!=b)) f++; }
      fesetround(FE_TONEAREST); }
    printf("%s:%ld ", mn[m], f); tot+=f; }
  printf("=> total %ld échecs\n", tot); return tot!=0; }
