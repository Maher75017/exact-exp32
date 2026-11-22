#!/bin/bash
# Chronométrage et re-vérification des programmes trouvés par la recherche génétique (moteur non publié)
cd "$(dirname "$0")"
which nvcc >/dev/null || export PATH=/usr/local/cuda/bin:$PATH
for p in c3b_time c3c_time; do
  nvcc -O3 -std=c++17 -I../gpu -arch=native $p.cu -o $p 2>&1 | grep -iv warning
  [ -x $p ] || nvcc -O3 -std=c++17 -I../gpu -arch=sm_120 $p.cu -o $p || exit 1
  ./$p
done
