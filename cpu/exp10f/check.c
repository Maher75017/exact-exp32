// Vérification de v_exp10f : liste difficile (hard.txt) ou exhaustive (--full) contre cr_exp10f de CORE-MATH
// usage : check params.txt SCHEME [--full]
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <fenv.h>
#include <omp.h>
float cr_exp10f(float); float v_exp10f(float);
extern int VB, VDEG, VSCHEME; extern double VT[256], VC[10];
static void load(const char* fn){ FILE*f=fopen(fn,"r"); if(!f){perror(fn);exit(1);} int s;
  if(fscanf(f,"%d %d %d",&VB,&VDEG,&s)!=3) exit(2);
  for(int j=0;j<VDEG;j++) if(fscanf(f,"%la",&VC[j])!=1) exit(3);
  for(int i=0;i<(1<<VB);i++) if(fscanf(f,"%la",&VT[i])!=1) exit(4);
  fclose(f); }
int main(int argc, char** argv){
  load(argv[1]); VSCHEME = atoi(argv[2]);
  int full = argc>3 && !strcmp(argv[3],"--full");
  int modes[4]={FE_TONEAREST,FE_TOWARDZERO,FE_UPWARD,FE_DOWNWARD};
  const char* mn[4]={"RN","RZ","RU","RD"};
  if(!full){
    FILE*f=fopen("hard.txt","r"); unsigned u,o[4]; double d; long fails=0, n=0;
    while(fscanf(f,"%x %x %x %x %x %le",&u,&o[0],&o[1],&o[2],&o[3],&d)==6){ n++;
      float x; memcpy(&x,&u,4);
      for(int m=0;m<4;m++){ fesetround(modes[m]); float y=v_exp10f(x); fesetround(FE_TONEAREST);
        uint32_t yb; memcpy(&yb,&y,4); if(yb!=o[m]){ fails++; if(fails<=20) printf("ECHEC %s x=%08x (%a) obtenu %08x attendu %08x dist=%.2e\n",mn[m],u,x,yb,o[m],d);} } }
    printf("liste difficile : %ld entrées, %ld échecs\n", n, fails); return fails!=0;
  }
  long tot_fail=0; FILE* fo = fopen("fails_full.txt","w");
  for(int m=0;m<4;m++){ long fails=0;
    #pragma omp parallel reduction(+:fails)
    { fesetround(modes[m]);
      #pragma omp for schedule(static,1<<16)
      for(long i=0;i<(1L<<32);i++){ uint32_t u=(uint32_t)i; float x; memcpy(&x,&u,4);
        float a=v_exp10f(x), b=cr_exp10f(x); uint32_t ab,bb; memcpy(&ab,&a,4); memcpy(&bb,&b,4);
        if(ab!=bb && !(a!=a && b!=b)){ fails++;
          #pragma omp critical
          { if(fails<=3) printf("  ECHEC %s x=%08x (%a) obtenu %08x attendu %08x\n",mn[m],u,x,ab,bb);
            fprintf(fo,"%08x %d\n",u,m); } } }
      fesetround(FE_TONEAREST); }
    printf("exhaustif %s : %ld échecs sur 2^32 entrées\n", mn[m], fails); fflush(stdout); tot_fail+=fails; }
  fclose(fo);
  return tot_fail!=0;
}
