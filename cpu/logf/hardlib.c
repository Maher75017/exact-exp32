// Bibliothèque pour la recherche logf : charge la liste difficile, compte les échecs pour des paramètres donnés
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <fenv.h>
float v_logf(float);
extern double VC[6], VIX[129], VLIX[129];
static uint32_t HX[8192], HO[8192][4]; static int NH=0;
int load_hard(const char* fn){ FILE*f=fopen(fn,"r"); unsigned u,o0,o1,o2,o3; double d; NH=0;
  while(fscanf(f,"%x %x %x %x %x %le",&u,&o0,&o1,&o2,&o3,&d)==6){ HX[NH]=u; HO[NH][0]=o0;HO[NH][1]=o1;HO[NH][2]=o2;HO[NH][3]=o3; NH++; }
  fclose(f); return NH; }
void set_params(const double* C, const double* IX, const double* L){ memcpy(VC,C,48); memcpy(VIX,IX,129*8); memcpy(VLIX,L,129*8); }
int count_fails(int* idx_out){ int modes[4]={FE_TONEAREST,FE_TOWARDZERO,FE_UPWARD,FE_DOWNWARD}; int nf=0;
  for(int m=0;m<4;m++){ fesetround(modes[m]);
    for(int i=0;i<NH;i++){ float x; memcpy(&x,&HX[i],4); float y=v_logf(x); uint32_t yb; memcpy(&yb,&y,4);
      if(yb!=HO[i][m]){ if(idx_out) idx_out[nf]=i*4+m; nf++; } } }
  fesetround(FE_TONEAREST); return nf; }
