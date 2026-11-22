# Essaie des décalages plus larges des entrées de table utilisées par l'entrée en échec (dans les 4 modes)
import ctypes, numpy as np, struct, math
lib = ctypes.CDLL("./libhard.so"); lib.load_hard(b"hard.txt")
dp = ctypes.POINTER(ctypes.c_double)
L = open("s1_6_5.txt").read().split(); B,d,s = int(L[0]),int(L[1]),int(L[2])
C = np.array([float.fromhex(v) for v in L[3:3+d]]); T0 = np.array([float.fromhex(v) for v in L[3+d:3+d+64]])
def sh(v,k): u=struct.unpack('<q',struct.pack('<d',v))[0]; return struct.unpack('<d',struct.pack('<q',u+k))[0]
x = float.fromhex("0x1.fafecp+3"); k = round(x*float.fromhex("0x1.a934f0979a371p+1")*64); i = k % 64
print("k =", k, "indice", i)
res=[]
for j in [i-1, i, i+1]:
  for kk in range(-40, 41):
    T=T0.copy(); T[j%64]=sh(T0[j%64],kk); lib.set_params(B,d,s,C.ctypes.data_as(dp),T.ctypes.data_as(dp))
    n=lib.count_fails(None)
    if n==0: res.append((abs(kk), j%64, kk))
print(sorted(res)[:10])
if res:
  _, j, kk = sorted(res)[0]; T=T0.copy(); T[j]=sh(T0[j],kk)
  open("s1_fix.txt","w").write(f"{B} {d} {s}\n"+"\n".join(float(c).hex() for c in C)+"\n"+"\n".join(float(t).hex() for t in T)+"\n")
  print("retenu : T[%d] %+d ulp" % (j, kk))
