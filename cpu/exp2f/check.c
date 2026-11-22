#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <fenv.h>
#include <omp.h>
float cr_exp2f(float); float v_exp2f(float);
extern int VB, VDEG, VSCHEME; extern double VT[256], VC[10];
// charge "params.txt" : B DEG SCHEME, puis d coefficients (hex), puis N entrées de table (hex)
static void load(const char* fn){ FILE*f=fopen(fn,"r"); if(!f){perror(fn);exit(1);}
  if(fscanf(f,"%d %d %d",&VB,&VDEG,&VSCHEME)!=3) exit(2);
  for(int j=0;j<VDEG;j++) if(fscanf(f,"%la",&VC[j])!=1) exit(3);
  for(int i=0;i<(1<<VB);i++) if(fscanf(f,"%la",&VT[i])!=1) exit(4);
  fclose(f); }
int main(int argc, char** argv){
  load(argv[1]);
  int full = argc>2 && !strcmp(argv[2],"--full");
  int modes[4]={FE_TONEAREST,FE_TOWARDZERO,FE_UPWARD,FE_DOWNWARD};
  const char* mn[4]={"RN","RZ","RU","RD"};
  if(!full){ // sur la liste des entrées difficiles (hard.txt : x + 4 résultats de référence)
    FILE*f=fopen("hard.txt","r"); unsigned u,o[4]; double d; long fails=0, n=0;
    while(fscanf(f,"%x %x %x %x %x %le",&u,&o[0],&o[1],&o[2],&o[3],&d)==6){ n++;
      float x; memcpy(&x,&u,4);
      for(int m=0;m<4;m++){ fesetround(modes[m]); float y=v_exp2f(x); fesetround(FE_TONEAREST);
        uint32_t yb; memcpy(&yb,&y,4); if(yb!=o[m]){ fails++; printf("ECHEC %s x=%08x obtenu %08x attendu %08x dist=%.2e\n",mn[m],u,yb,o[m],d);} } }
    printf("liste difficile : %ld entrées, %ld échecs\n", n, fails); return fails!=0;
  }
  long tot_fail=0;
  for(int m=0;m<4;m++){ long fails=0;
    #pragma omp parallel reduction(+:fails)
    { fesetround(modes[m]);
      #pragma omp for schedule(static,1<<16)
      for(long i=0;i<(1L<<32);i++){ uint32_t u=(uint32_t)i; float x; memcpy(&x,&u,4);
        float a=v_exp2f(x), b=cr_exp2f(x); uint32_t ab,bb; memcpy(&ab,&a,4); memcpy(&bb,&b,4);
        if(ab!=bb && !(a!=a && b!=b)){ fails++; if(fails<=5){
          #pragma omp critical
          printf("  ECHEC %s x=%08x obtenu %08x attendu %08x\n",mn[m],u,ab,bb);} } }
      fesetround(FE_TONEAREST); }
    printf("exhaustif %s : %ld échecs sur 2^32 entrées\n", mn[m], fails); tot_fail+=fails; }
  return tot_fail!=0;
}
