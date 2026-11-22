# Règle coefficients (continus) + décalages entiers de table (en ulp) pour que l'arrondi final soit correct
# sur toutes les entrées difficiles et les 4 modes : programme linéaire en nombres entiers (HiGHS via scipy).
# Modèle : v/y - 1 ≈ e_tab(i) + δ_i + e_poly(h) + h·Δp(h)/2^h   (erreurs d'évaluation ~2^-60, négligées)
import sys, struct, numpy as np, mpmath as mp
from scipy.optimize import milp, LinearConstraint, Bounds
mp.mp.prec = 200
B, d = int(sys.argv[1]), int(sys.argv[2]); N = 1 << B
extra = sys.argv[4] if len(sys.argv) > 4 else None
L = open(sys.argv[3]).read().split(); C0 = [float.fromhex(v) for v in L[3:3+d]]; T0 = [float.fromhex(v) for v in L[3+d:3+d+N]]
def f2u(f): return struct.unpack('<I', struct.pack('<f', f))[0]
def u2f(u): return struct.unpack('<f', struct.pack('<I', u))[0]
def ulp_d(v): return mp.mpf(2)**(mp.floor(mp.log(abs(v), 2)) - 52)
rows = []  # (i, h, Lrel, Urel)
lines = [l.split() for l in open("hard.txt")]
if extra: lines += [l.split() for l in open(extra)]
U52 = mp.mpf(2)**-52
for parts in lines:
    u = int(parts[0], 16); x = mp.mpf(u2f(u)); outs = [int(p, 16) for p in parts[1:5]]
    y = mp.power(2, x)
    for m in range(4):
        # indice k selon le mode (somme x*N + SH positive : RZ = vers le bas)
        xn = x*N
        k = int(mp.nint(xn)) if m == 0 else (int(mp.floor(xn)) if m in (1, 3) else int(mp.ceil(xn)))
        if m == 0 and abs(xn - mp.floor(xn) - mp.mpf(0.5)) == 0: k = int(mp.floor(xn)) + (int(mp.floor(xn)) % 2)
        i = k % N; h = x - mp.mpf(k)/N
        o = mp.mpf(u2f(outs[m]))
        if o == y: continue                      # cas exacts : h = 0, i = 0, rien à régler
        nxt = mp.mpf(u2f(outs[m] + 1)); prv = mp.mpf(u2f(outs[m] - 1))
        if m == 0:   lo, hi = (o + prv)/2, (o + nxt)/2; mg = ulp_d(o) if (outs[m] & 1) else 0          # RN : marge d'un demi-ulp double de chaque côté (≥)
        elif m == 2: lo, hi = prv, o; mg = 0                                      # RU : v ∈ (prv, o]
        else:        lo, hi = o, nxt; mg = 0                                      # RD/RZ : v ∈ [o, nxt)
        eps = mp.mpf(2)**-58                                                     # erreurs d'évaluation négligées + bords stricts
        Lr = ((lo + mg/2) / y - 1) + eps; Ur = ((hi - mg/2) / y - 1) - eps
        rows.append((i, h, Lr, Ur))
print("contraintes :", len(rows))
# variables : u_0..u_{d-1} (coefficients, échelle 2^-52 * 2^(7(j+1))), k_1..k_{N-1} (entiers, ulps de table), t (marge)
nv = d + (N - 1) + 1
A = []; lb = []; ub = []
scale_c = [mp.mpf(2)**(7*(j+1)) for j in range(d)]
for (i, h, Lr, Ur) in rows:
    yh = mp.power(2, h)
    e_tab = mp.mpf(T0[i]) / mp.power(2, mp.mpf(i)/N) - 1
    p0 = sum(mp.mpf(C0[j]) * h**j for j in range(d))
    e_poly = (1 + h*p0) / yh - 1
    base = (e_tab + e_poly) / U52
    row = np.zeros(nv)
    for j in range(d): row[j] = float(h**(j+1) / yh * scale_c[j])
    if i > 0: row[d + i - 1] = float(ulp_d(T0[i]) / T0[i] / U52)
    # Lr/U52 - base <= row·x - t  et  row·x + t <= Ur/U52 - base
    r1 = row.copy(); r1[-1] = -1; A.append(r1); lb.append(float(Lr/U52 - base)); ub.append(np.inf)
    r2 = row.copy(); r2[-1] = +1; A.append(r2); lb.append(-np.inf); ub.append(float(Ur/U52 - base))
A = np.array(A)
z = np.zeros(nv); viol = sum(1 for k in range(0, len(A), 2) if (A[k]@z < lb[k]) or (A[k+1]@z > ub[k+1]))
print('contraintes violées par les paramètres de départ :', viol)
w = np.array([ub[k+1]-lb[k] for k in range(0,len(A),2)]); print('largeur min des intervalles (unités 2^-52) : %.3f' % w.min(), '; nb < 1 :', (w<1).sum())
integ = np.zeros(nv); integ[d:d+N-1] = 1
lo_b = np.full(nv, -1e6); hi_b = np.full(nv, 1e6); lo_b[d:d+N-1] = -64; hi_b[d:d+N-1] = 64; lo_b[-1] = -50; hi_b[-1] = 50
c = np.zeros(nv); c[-1] = -1                       # maximiser la marge t
res = milp(c, constraints=LinearConstraint(A, lb, ub), integrality=integ, bounds=Bounds(lo_b, hi_b), options={"time_limit": 300})
print("statut :", res.status, res.message)
if res.x is None: sys.exit(1)
x = res.x; t = x[-1]
print("marge minimale (en unités 2^-52) : %.4f" % t)
C = [float(mp.mpf(C0[j]) + mp.mpf(x[j]) * U52 * scale_c[j]) for j in range(d)]
def shift(v, k): q = struct.unpack('<q', struct.pack('<d', v))[0]; return struct.unpack('<d', struct.pack('<q', q + int(k)))[0]
ks = [0] + [int(round(x[d + i - 1])) for i in range(1, N)]
T = [shift(T0[i], ks[i]) for i in range(N)]
print("décalages de table non nuls :", sum(1 for k in ks if k), "; max |k| =", max(abs(k) for k in ks))
out = sys.argv[3].replace(".txt", "_milp.txt")
open(out, "w").write(f"{B} {d} 0\n" + "\n".join(c.hex() for c in C) + "\n" + "\n".join(t.hex() for t in T) + "\n")
print("écrit", out)
# --- diagnostic : contraintes les plus violées, relaxation continue
sl = []
for k in range(0, len(A), 2):
    v = A[k][:-1] @ x[:-1]
    sl.append((min(v - lb[k], ub[k+1] - v), k//2))
sl.sort()
from collections import Counter
cnt = Counter(r[0] for r in rows)
for s_, r in sl[:8]:
    i, h, Lr, Ur = rows[r]; print("  marge %.2f  tranche %3d (contraintes dans la tranche : %d)  h = %+.6f" % (s_, i, cnt[i], float(h)))
res2 = milp(c, constraints=LinearConstraint(A, lb, ub), integrality=np.zeros(nv), bounds=Bounds(lo_b, hi_b))
print("relaxation continue : marge %.4f" % res2.x[-1])
