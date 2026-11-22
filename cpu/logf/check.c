// Vérification de v_logf : liste difficile (hard.txt) ou exhaustive (--full) contre cr_logf de CORE-MATH
// usage : check params.txt SCHEME [--full]
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <fenv.h>
#include <omp.h>
float cr_logf(float); float v_logf(float);
#include "io.h"
int main(int argc, char** argv){
  load(argv[1]); { extern int VSTRUCT; char* s=getenv("VSTRUCT"); if(s) VSTRUCT=atoi(s); extern int VGRID; s=getenv("VGRID"); if(s) VGRID=atoi(s); }
  int full = argc>2 && !strcmp(argv[2],"--full");
  int modes[4]={FE_TONEAREST,FE_TOWARDZERO,FE_UPWARD,FE_DOWNWARD};
  const char* mn[4]={"RN","RZ","RU","RD"};
  if(!full){
    FILE*f=fopen("hard.txt","r"); unsigned u,o[4]; double d; long fails=0, n=0;
    while(fscanf(f,"%x %x %x %x %x %le",&u,&o[0],&o[1],&o[2],&o[3],&d)==6){ n++;
      float x; memcpy(&x,&u,4);
      for(int m=0;m<4;m++){ fesetround(modes[m]); float y=v_logf(x); fesetround(FE_TONEAREST);
        uint32_t yb; memcpy(&yb,&y,4); if(yb!=o[m]){ fails++; if(fails<=20) printf("ECHEC %s x=%08x (%a) obtenu %08x attendu %08x dist=%.2e\n",mn[m],u,x,yb,o[m],d);} } }
    { extern long VNFALL; printf("liste difficile : %ld entrées, %ld échecs, %ld renvois\n", n, fails, VNFALL); } return fails!=0;
  }
  long tot_fail=0; FILE* fo = fopen("fails_full.txt","w");
  for(int m=0;m<4;m++){ long fails=0;
    #pragma omp parallel reduction(+:fails)
    { fesetround(modes[m]);
      #pragma omp for schedule(static,1<<16)
      for(long i=0;i<(1L<<32);i++){ uint32_t u=(uint32_t)i; float x; memcpy(&x,&u,4);
        float a=v_logf(x), b=cr_logf(x); uint32_t ab,bb; memcpy(&ab,&a,4); memcpy(&bb,&b,4);
        if(ab!=bb && !(a!=a && b!=b)){ fails++;
          #pragma omp critical
          { if(fails<=3) printf("  ECHEC %s x=%08x (%a) obtenu %08x attendu %08x\n",mn[m],u,x,ab,bb);
            fprintf(fo,"%08x %d\n",u,m); } } }
      fesetround(FE_TONEAREST); }
    printf("exhaustif %s : %ld échecs sur 2^32 entrées\n", mn[m], fails); fflush(stdout); tot_fail+=fails; }
  fclose(fo);
  return tot_fail!=0;
}
