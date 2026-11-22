#!/bin/bash
cd "$(dirname "$0")"
exec > >(tee resultats_c1.txt) 2>&1
which nvcc >/dev/null || export PATH=/usr/local/cuda/bin:$PATH
gcc -O2 -march=native -frounding-math -c exp2f_orig.c -o exp2f_orig.o || exit 1
nvcc -O3 -arch=native -Xcompiler -fopenmp c1_double_baseline.cu exp2f_orig.o -o c1 -lgomp 2>&1 | grep -v "warning" \
  || true
[ -x c1 ] || nvcc -O3 -arch=sm_120 -Xcompiler -fopenmp c1_double_baseline.cu exp2f_orig.o -o c1 -lgomp || exit 1
./c1
