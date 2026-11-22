# Recherche évolutive (1+λ) : décalages entiers (ulp) des coefficients c1..c5 et des entrées lix[j]
# pour éliminer les échecs d'arrondi correct sur la liste difficile de logf. c0 = 1 et lix[0] = 0 restent exacts.
import ctypes, numpy as np, struct, sys
lib = ctypes.CDLL("./libhard.so"); print("entrées difficiles :", lib.load_hard(b"hard.txt"))
dp = ctypes.POINTER(ctypes.c_double)
L = open(sys.argv[1]).read().split()
C0 = np.array([float.fromhex(v) for v in L[2:8]]); IX = np.array([float.fromhex(v) for v in L[8:137]]); LX0 = np.array([float.fromhex(v) for v in L[137:266]])
def sh(v, k): u = struct.unpack('<q', struct.pack('<d', v))[0]; return struct.unpack('<d', struct.pack('<q', u + int(k)))[0]
def build(dc, dl):
    C = np.array([C0[0]] + [sh(C0[j], dc[j]) for j in range(1, 6)])
    LX = np.array([LX0[i] if i == 0 else sh(LX0[i], dl[i]) for i in range(129)])
    return C, LX
def fails(C, LX, want=False):
    lib.set_params(C.ctypes.data_as(dp), IX.ctypes.data_as(dp), LX.ctypes.data_as(dp))
    buf = (ctypes.c_int * 32768)(); n = lib.count_fails(buf if want else None)
    return (n, [buf[i] for i in range(n)]) if want else n
HX = [int(l.split()[0], 16) for l in open("hard.txt")]
def jidx(u):  # indice de table utilisé par l'entrée u (même calcul que le code)
    m = (u & 0x7fffff) << 29; return (m + (1 << 44)) >> 45
rng = np.random.default_rng(int(sys.argv[2]) if len(sys.argv) > 2 else 1)
dc = np.zeros(6, int); dl = np.zeros(129, int)
C, LX = build(dc, dl); n0, idx = fails(C, LX, True); best = (n0, dc, dl, idx)
print("départ :", n0, "échecs")
scale = [0, 2**2, 2**6, 2**12, 2**18, 2**24]
for it in range(int(sys.argv[3]) if len(sys.argv) > 3 else 5000):
    dc = best[1].copy(); dl = best[2].copy(); r = rng.random()
    if r < 0.6 and best[3]:   # viser l'entrée de table d'un échec
        f = best[3][rng.integers(len(best[3]))] // 4; j = jidx(HX[f])
        if j > 0: dl[j] += int(rng.choice([-1, 1]) * rng.integers(1, 4))
    elif r < 0.9:
        j = rng.integers(1, 6); dc[j] += int(rng.normal(0, scale[j]) or rng.choice([-1, 1]))
    else:
        j = rng.integers(1, 129); dl[j] += int(rng.integers(-2, 3))
    C, LX = build(dc, dl); n, idx = fails(C, LX, True)
    cost = lambda a, b: np.sum(np.abs(a[1:]) / np.array(scale[1:])) + 0.01 * np.sum(np.abs(b))
    if n < best[0] or (n == best[0] and cost(dc, dl) < cost(best[1], best[2])):
        if n < best[0]: print(f"  it {it}: {n} échecs ; coefs {list(map(int, dc[1:]))} ; table {[(i, int(v)) for i, v in enumerate(dl) if v]}", flush=True)
        best = (n, dc, dl, idx)
    if best[0] == 0 and it > 500: break
n, dc, dl, _ = best; C, LX = build(dc, dl)
print("final :", fails(C, LX), "échecs ; coefs (ulp) :", list(map(int, dc[1:])), "; table :", [(i, int(v)) for i, v in enumerate(dl) if v])
out = sys.argv[1].replace(".txt", "_fix.txt")
open(out, "w").write("129 6\n" + "\n".join(c.hex() for c in C) + "\n" + "\n".join(v.hex() for v in IX) + "\n" + "\n".join(v.hex() for v in LX) + "\n")
print("écrit", out)
