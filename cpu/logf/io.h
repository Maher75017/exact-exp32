#include <stdio.h>
#include <stdlib.h>
extern double VC[6], VIX[129], VLIX[129];
static void load(const char* fn){ FILE*f=fopen(fn,"r"); if(!f){perror(fn);exit(1);} int n,d;
  if(fscanf(f,"%d %d",&n,&d)!=2) exit(2);
  for(int j=0;j<6;j++) if(fscanf(f,"%la",&VC[j])!=1) exit(3);
  for(int i=0;i<129;i++) if(fscanf(f,"%la",&VIX[i])!=1) exit(4);
  for(int i=0;i<129;i++) if(fscanf(f,"%la",&VLIX[i])!=1) exit(5);
  fclose(f); }
