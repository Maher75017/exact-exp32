#!/bin/bash
cd "$(dirname "$0")"
exec > >(tee resultats_c3a.txt) 2>&1
which nvcc >/dev/null || export PATH=/usr/local/cuda/bin:$PATH
gcc -O2 -march=native -frounding-math -c exp2f_orig.c -o exp2f_orig.o || exit 1
rm -f c3a
nvcc -O3 -std=c++17 -arch=native -Xcompiler -fopenmp c3a.cu exp2f_orig.o -o c3a -lgomp 2>&1 | grep -iv warning
[ -x c3a ] || nvcc -O3 -std=c++17 -arch=sm_120 -Xcompiler -fopenmp c3a.cu exp2f_orig.o -o c3a -lgomp || exit 1
./c3a
