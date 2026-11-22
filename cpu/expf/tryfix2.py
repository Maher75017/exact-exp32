import ctypes, numpy as np, struct
lib = ctypes.CDLL("./libhard.so"); lib.load_hard(b"hard.txt")
dp = ctypes.POINTER(ctypes.c_double)
L = open("s1_6_5.txt").read().split(); B,d,s = int(L[0]),int(L[1]),int(L[2])
C = np.array([float.fromhex(v) for v in L[3:3+d]]); T0 = np.array([float.fromhex(v) for v in L[3+d:3+d+64]])
def sh(v,k): u=struct.unpack('<q',struct.pack('<d',v))[0]; return struct.unpack('<d',struct.pack('<q',u+k))[0]
best=None
for a in range(-4,1):
  for b in range(-4,1):
    T=T0.copy(); T[5]=sh(T0[5],a); T[6]=sh(T0[6],b); lib.set_params(B,d,s,C.ctypes.data_as(dp),T.ctypes.data_as(dp))
    n=lib.count_fails(None)
    if n==0: print(f"T[5] {a:+d}, T[6] {b:+d} : 0 échec")
    if n==0 and (best is None or abs(a)+abs(b) < abs(best[0])+abs(best[1])): best=(a,b,T.copy())
a,b,T=best; print("retenu", a, b)
open("s1_fix.txt","w").write(f"{B} {d} {s}\n"+"\n".join(float(c).hex() for c in C)+"\n"+"\n".join(float(t).hex() for t in T)+"\n")
