# exact-exp32

Correctly rounded single-precision exponentials (`exp2f`, `expf`, `exp10f`) without a rounding test on the CPU, a float-only correctly rounded `exp2f` on an NVIDIA GPU, and a genetic-programming (Koza-style) search for faster exact code, every candidate checked on all 2^32 inputs.

Author: Maher El Hmissi, November 2026. Français : [README.fr.md](README.fr.md).

## Results at a glance

| Part | What it is | Verified | Speed (honest) |
| --- | --- | --- | --- |
| `cpu/exp2f.c`, `cpu/expf.c`, `cpu/exp10f.c` | Drop-in CORE-MATH files: table of 64, degree-5 polynomial, **no rounding test, no slow path** | 0 errors on all 2^32 inputs × 4 rounding modes, against CORE-MATH, with gcc and clang | No net gain with CORE-MATH's official `perf.sh` on a dedicated Ryzen 7 7800X3D: −9 % (gcc), +9 to +10 % (clang), same latency |
| `cpu/logf.c` | Same idea for `logf`, plus a cheap midpoint check for 5 inputs no double computation can decide | 0 errors, 2^32 inputs × 4 modes | No clear gain |
| `gpu/fexp2.h` | **Float-only** (double-float) correctly rounded `exp2f`, round-to-nearest, for GPUs with weak FP64 | 0 errors on all 2^32 inputs, on an RTX 5090 | ×1.02 the time of CUDA's `exp2f` at 1 call per element read, ×1.17 at 4, ×2.88 at 16, ×5.7 in pure compute |
| `koza/` | Linear genetic programming on the GPU over the fast path; promising candidates verified exhaustively, counterexamples fed back | Winner exact on 2^32 inputs | 15 operations instead of 21; 7.2 % faster than the best hand-written version (goal was 10 %) |

CUDA's own `exp2f` returns a wrong result for 189,192,458 of the 4,294,967,296 inputs (4.40 %), with errors up to 2 ulp.

## Prior art: the CPU part is an independent reproduction, not a first

A search done after the work found that:

- **Arm optimized-routines** ships correctly rounded AdvSIMD `expf`, `exp2f`, `exp10f` built with the same recipe (x/ln2 rounded to 1/64, table of 64 values of 2^(i/64), degree-5 polynomial, computed in double).
- **LLVM libc** contains a float-only `expf` correctly rounded in round-to-nearest and designed for GPUs (`float_eval_cr_rn`). As of November 2026 its float `exp2f` is not correctly rounded; the one in `gpu/fexp2.h` is.
- **Intel** (arXiv 2605.15547, 2026) reports correctly rounded `expf` on its GPUs, computed mainly in double precision, and names 32-bit internal computation as future work.
- CORE-MATH's `log2f` is already free of a rounding test; RLIBM generates correctly rounded functions by constrained polynomial fitting.
- An integer-only, statically rounded `expf` was proposed for LLVM libc in August 2026 (PR #209406): another way to avoid double precision, for `expf` rather than `exp2f`.
- Evolutionary search for faster code already exists, for example GEVO (evolutionary optimization of GPU kernels) and STOKE (stochastic superoptimization of floating-point code with tunable precision). Neither requires every candidate to be exact on all inputs.

To our knowledge (search done in November 2026), `gpu/fexp2.h` is the first correctly rounded float-only `exp2f` for CUDA, and `koza/` the first GPU genetic-programming search in which every candidate is verified on all 2^32 inputs.

What this repository adds: the data (hard-case lists, exhaustive checks, official benchmarks on a dedicated machine), a float-only correctly rounded `exp2f` for CUDA, and the GPU genetic-programming search, including a 3-operation compensated step it found: `sh = fma(th, ph, th); d = sh - th; err ≈ fma(th, ph, -d)`.

## Method (CPU)

1. Hard cases: every input whose exact result lies within 2^-22 ulp of a rounding boundary, computed with MPFR (`hard.c`, a few seconds on two cores).
2. Start from correctly rounded table entries and a near-minimax polynomial (`coeffs.py`).
3. A (1+λ) evolutionary search shifts table entries and coefficients by a few ulps until every hard case rounds correctly in all four modes (`search.py`).
4. Exhaustive check of all 2^32 inputs in the 4 modes against CORE-MATH (`check.c --full`, `verif.c`).

Adjustments found: `exp2f` table 13 −2 ulp, table 62 +2 ulp, c0 +111 ulp, c1 +55 ulp; `expf` table 5 and 6 −2 ulp; `exp10f` table 40 +12 ulp.

## Reproduce

Requirements: gcc or clang, libmpfr-dev, Python 3 with numpy and mpmath; for the GPU part, CUDA 12+ and an NVIDIA GPU (tested on an RTX 5090, sm_120).

```sh
# CPU, e.g. expf: hard cases, then exhaustive check of the tuned parameters
cd cpu/expf
gcc -O2 -fopenmp -o hard hard.c cm_expf_ref.c -lmpfr -lm && ./hard
gcc -O2 -fopenmp -ffp-contract=off -o check check.c v.c cm_expf_ref.c -lm
./check s1_fix.txt 1 --full          # 0 errors expected in RN, RZ, RU, RD

# Clean file against current CORE-MATH (all 2^32 inputs, 4 modes)
gcc -O2 -fopenmp -Dcr_expf=new_expf -c ../expf.c -o n.o
gcc -O2 -fopenmp -c cm_expf_ref.c -o r.o && gcc -O2 -fopenmp verif.c n.o r.o -o verif_x -lm && ./verif_x

# logf: the cheap midpoint check is switched on with VGRID=1
cd ../logf && gcc -O2 -fopenmp -ffp-contract=off -o check check.c v.c cm_logf_ref.c -lm
VGRID=1 ./check p_base.txt --full

# Official CORE-MATH benchmark, old vs new (Linux, about 30-40 minutes)
cd ../bench && bash run_officiel.sh

# GPU: exhaustive check and timings of the float-only exp2f against CUDA
cd ../../gpu && bash run_c2.sh

# Timing and re-verification of the programs found by the Koza search
cd ../koza && bash run_verif.sh
```

Logs of the RTX 5090 runs are in `koza/results/`, and the official Ryzen benchmark is in `cpu/bench/results_ryzen7800x3d_official.txt`.

## A lesson about benchmarks

A home-made harness on a shared cloud machine first suggested 17 to 26 % faster CPU code. CORE-MATH's official tool on a dedicated machine showed no net gain. Always measure with the maintainers' tool, on dedicated hardware, before claiming a speed-up.

## Layout

```
cpu/          drop-in CORE-MATH files: exp2f.c expf.c exp10f.c logf.c
cpu/<fn>/     hard cases, search, exhaustive checkers, reference copy of CORE-MATH
cpu/bench/    official perf.sh comparison script and Ryzen results
gpu/          float-only exact exp2f (fexp2.h, tabf.h), double fallback (cr2.h), benchmarks
koza/         programs found by the GPU genetic-programming search, timing, re-verification, RTX 5090 logs
SHA256SUMS    fingerprints of every file
```

## See also

[tt-erf-fp32](https://github.com/Maher75017/tt-erf-fp32): FP32 `erf` for the Tenstorrent SFPU (Blackhole, Wormhole), at most 1 ULP on every input, checked with the same exhaustive method on a bit-exact model of the hardware arithmetic. Proposed upstream in [tenstorrent/tt-metal#66809](https://github.com/tenstorrent/tt-metal/pull/66809).

## Note on the search engine

The GPU genetic-programming engine (Koza-style, linear GP with exhaustive verification and counterexample feedback) is the author's private research tool and is not published here. The programs it found, their exhaustive re-verification and their timings are published.

## License

MIT (see `LICENSE`). Files derived from CORE-MATH keep their original MIT notice.
