# Recherche locale (évolutive, 1+λ) : décalages entiers en ulp des coefficients et des entrées de table
# pour éliminer les échecs d'arrondi correct sur la liste des entrées difficiles. T[0] reste exactement 1.
import ctypes, numpy as np, struct, sys
lib = ctypes.CDLL("./libhard.so")
lib.load_hard.restype = ctypes.c_int
print("entrées difficiles :", lib.load_hard(b"hard.txt"))
dp = ctypes.POINTER(ctypes.c_double)
def rd(fn):
    L = open(fn).read().split(); B, d, s = int(L[0]), int(L[1]), int(L[2])
    C = [float.fromhex(v) for v in L[3:3+d]]; T = [float.fromhex(v) for v in L[3+d:3+d+(1<<B)]]
    return B, d, s, np.array(C), np.array(T)
def shift_ulp(v, k):  # décale un double de k ulps
    u = struct.unpack('<q', struct.pack('<d', v))[0]; return struct.unpack('<d', struct.pack('<q', u + int(k)))[0]
def fails(B, d, s, C, T, want_idx=False):
    lib.set_params(B, d, s, C.ctypes.data_as(dp), T.ctypes.data_as(dp))
    buf = (ctypes.c_int * 8192)()
    n = lib.count_fails(buf if want_idx else None)
    return (n, [buf[i] for i in range(n)]) if want_idx else n
fn = sys.argv[1]; B, d, s, C0, T0 = rd(fn)
n0, idx = fails(B, d, s, C0, T0, True)
print("départ :", n0, "échecs")
rng = np.random.default_rng(int(sys.argv[2]) if len(sys.argv) > 2 else 1)
best = (n0, np.zeros(d, int), np.zeros(1 << B, int))
def build(dc, dt):
    C = np.array([shift_ulp(C0[j], dc[j]) for j in range(d)])
    T = np.array([T0[i] if i == 0 else shift_ulp(T0[i], dt[i]) for i in range(1 << B)])
    return C, T
# échelles de mutation par coefficient (en ulp) : c1 plus fin, les suivants plus grossiers
scale = [64, 2**10, 2**16, 2**22, 2**28, 2**34][:d]
for it in range(int(sys.argv[3]) if len(sys.argv) > 3 else 4000):
    dc = best[1].copy(); dt = best[2].copy()
    r = rng.random()
    if r < 0.5:
        j = rng.integers(d); dc[j] += int(rng.normal(0, scale[j]))
    elif r < 0.8:
        i = rng.integers(1, 1 << B); dt[i] += int(rng.integers(-3, 4))
    else:
        j = rng.integers(d); dc[j] += int(rng.normal(0, scale[j])); i = rng.integers(1, 1 << B); dt[i] += int(rng.integers(-2, 3))
    C, T = build(dc, dt); n = fails(B, d, s, C, T)
    # à nombre d'échecs égal, préférer les plus petits décalages
    cost = lambda a, b: np.sum(np.abs(a) / np.array(scale)) + 0.01 * np.sum(np.abs(b))
    if n < best[0] or (n == best[0] and cost(dc, dt) < cost(best[1], best[2])):
        if n < best[0]: print(f"  it {it}: {n} échecs  décalages coefs {list(dc)} table {[(i, int(v)) for i, v in enumerate(dt) if v]}")
        best = (n, dc, dt)
    if best[0] == 0 and it > 300: break
n, dc, dt = best; C, T = build(dc, dt)
print("final :", fails(B, d, s, C, T), "échecs ; décalages coefs (ulp) :", list(dc), "; table :", [(i, int(v)) for i, v in enumerate(dt) if v])
out = fn.replace(".txt", "_fix.txt")
with open(out, "w") as f:
    f.write(f"{B} {d} {s}\n" + "\n".join(float(c).hex() for c in C) + "\n" + "\n".join(float(t).hex() for t in T) + "\n")
print("écrit", out)
