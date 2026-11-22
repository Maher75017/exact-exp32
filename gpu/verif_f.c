#include <stdio.h>
#include "fexp2.h"
float cr_exp2f(float);
int main(){ long bad=0, slow=0, dom=0;
  #pragma omp parallel for reduction(+:bad,slow,dom) schedule(static,1<<16)
  for (long i=0;i<(1L<<32);i++){ uint32_t u=(uint32_t)i; float x=u2f(u); float r;
    int ok = fexp2_fast(x, TF, &r);
    uint32_t uxx = u<<1; int indom = !(uxx < 0x65000000u || uxx >= 0x85fffc00u || u > 0xc2fa0000u);
    dom += indom; if (indom && !ok) slow++;
    float y = ok ? r : cr2_exp2f(x, CR2_TB); float z = cr_exp2f(x);
    if (f2u(y)!=f2u(z) && !(y!=y && z!=z)) { bad++; } }
  printf("E=%a : erreurs %ld ; chemin rapide refusé %ld fois sur %ld entrées du domaine (1 sur %.0f)\n", (double)FEXP2_E, bad, slow, dom, (double)dom/(slow?slow:1));
  return bad!=0; }
