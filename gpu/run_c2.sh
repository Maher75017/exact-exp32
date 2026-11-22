#!/bin/bash
cd "$(dirname "$0")"
exec > >(tee resultats_c2.txt) 2>&1
which nvcc >/dev/null || export PATH=/usr/local/cuda/bin:$PATH
gcc -O2 -march=native -frounding-math -c exp2f_orig.c -o exp2f_orig.o || exit 1
rm -f c2
nvcc -O3 -arch=native -Xcompiler -fopenmp c2.cu exp2f_orig.o -o c2 -lgomp 2>&1 | grep -iv warning
[ -x c2 ] || nvcc -O3 -arch=sm_120 -Xcompiler -fopenmp c2.cu exp2f_orig.o -o c2 -lgomp || exit 1
./c2
