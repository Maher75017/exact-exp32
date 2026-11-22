#!/bin/bash
# Official CORE-MATH benchmark (perf.sh): current CORE-MATH vs this repository, for exp2f, expf, exp10f.
# Usage (from cpu/bench): bash run_officiel.sh            -> timings (about 30-40 min)
#                         bash run_officiel.sh --check    -> also exhaustive MPFR check of the 3 new functions (hours)
set -e
ICI="$(cd "$(dirname "$0")" && pwd)"
cd "$ICI"
( cd .. && sha256sum exp2f.c expf.c exp10f.c )
command -v git >/dev/null || { apt-get update -qq && apt-get install -y -qq git; }
dpkg -s libmpfr-dev >/dev/null 2>&1 || { apt-get update -qq && apt-get install -y -qq libmpfr-dev libgmp-dev; }
command -v clang >/dev/null || apt-get install -y -qq clang || true
if [ ! -d core-math ]; then
  git clone -q https://gitlab.inria.fr/core-math/core-math.git || git clone -q https://github.com/mockingbirdnest/core-math.git
fi
rm -rf old new; cp -r core-math old; cp -r core-math new
cp ../exp2f.c new/src/binary32/exp2/exp2f.c
cp ../expf.c new/src/binary32/exp/expf.c
cp ../exp10f.c new/src/binary32/exp10/exp10f.c
echo "CORE-MATH : $(cd core-math && git log -1 --format='%h %ci')" | tee resultats.txt
lscpu | grep -m1 "Model name" | tee -a resultats.txt
export CORE_MATH_PERF_MODE=rdtsc CORE_MATH_QUIET=1 CORE_MATH_LAUNCHER="taskset -c 2"
for cc in gcc clang; do
  command -v $cc >/dev/null || continue
  if [ $cc = clang ]; then export CC=clang CFLAGS="-O3 -march=native -ffp-contract=on"; else unset CC CFLAGS; fi
  for f in exp2f expf exp10f; do
    for essai in 1 2 3; do
      for v in old new; do
        t=$(cd $v && ./perf.sh $f 2>/dev/null | tail -1)
        l=$(cd $v && PERF_ARGS=--latency ./perf.sh $f 2>/dev/null | tail -1)
        echo "$cc $f essai$essai $v debit=$t latence=$l" | tee -a resultats.txt
      done
    done
  done
done
echo "--- meilleur de 3 (cycles par appel ; plus petit = plus rapide) ---" | tee -a resultats.txt
awk '{split($5,a,"="); split($6,b,"="); a[2]+=0; b[2]+=0; k=$1" "$2" "$4; if(!(k in T)||a[2]<T[k])T[k]=a[2]; if(!(k in L)||b[2]<L[k])L[k]=b[2]}
     END{for(k in T){split(k,p," "); if(p[3]=="old"){n=p[1]" "p[2]" new"; printf "%-6s %-7s debit %6.3f -> %6.3f (%+5.1f %%)   latence %6.3f -> %6.3f (%+5.1f %%)\n", p[1], p[2], T[k], T[n], 100*(T[n]/T[k]-1), L[k], L[n], 100*(L[n]/L[k]-1)}}}' \
     <(grep essai resultats.txt) | sort | tee -a resultats.txt
if [ "$1" = "--check" ]; then
  for f in exp2f expf exp10f; do (cd new && ./check.sh --exhaustive $f) 2>&1 | tail -5 | sed "s/^/$f : /" | tee -a resultats.txt; done
fi
echo "Fini. Résultats dans $ICI/resultats.txt"
