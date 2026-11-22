# Table 2^(i/128) en double-float (hi, lo), coefficients de q = 2^h - 1 sur |h| <= 2^-8
import mpmath as mp, struct
mp.mp.prec = 200
def f32(v): return struct.unpack('<f', struct.pack('<f', float(v)))[0]
N = 128
hi = []; lo = []
for i in range(N):
    v = mp.power(2, mp.mpf(i)/N); t = f32(v); hi.append(t); lo.append(f32(v - mp.mpf(t)))
ln2 = mp.log(2); c1h = f32(ln2); c1l = f32(ln2 - mp.mpf(c1h))
A_ = mp.mpf(2)**-8
C1 = mp.mpf(c1h) + mp.mpf(c1l)
def g2(h):
    if abs(h) < mp.mpf(10)**-20: return ln2**2/2 + ln2**3*h/6 + ln2**4*h**2/24
    return (mp.expm1(h*ln2) - ln2*h)/h**2
nodes = [A_*mp.cos(mp.pi*(2*k+1)/6) for k in range(3)]
M = mp.matrix([[x**j for j in range(3)] for x in nodes]); rhs = mp.matrix([g2(x) for x in nodes])
c = mp.lu_solve(M, rhs); print("c =", c.T)
c2 = f32(c[0])
nodes2 = [A_*mp.cos(mp.pi*(2*k+1)/4) for k in range(2)]
M2 = mp.matrix([[x, x**2] for x in nodes2]); rhs2 = mp.matrix([g2(x) - mp.mpf(c2) for x in nodes2])
cc = mp.lu_solve(M2, rhs2); c3 = f32(cc[0]); c4 = f32(cc[1])
m = 0
for k in range(-2000, 2001):
    h = A_*k/2000
    if h == 0: continue
    q = C1*h + h**2*(mp.mpf(c2) + h*(mp.mpf(c3) + h*mp.mpf(c4)))
    m = max(m, abs((1+q)/mp.power(2,h) - 1))
print("erreur relative du modèle (coefficients float) : 2^%.2f" % float(mp.log(m,2)))
def fx(v): return float(v).hex() + "f"
with open("tabf.h","w") as f:
    f.write("static const float TFH[128] = {" + ", ".join(fx(v) for v in hi) + "};\n")
    f.write("static const float TFL[128] = {" + ", ".join(fx(v) for v in lo) + "};\n")
    f.write(f"#define C1H {fx(c1h)}\n#define C1L {fx(c1l)}\n#define C2 {fx(c2)}\n#define C3 {fx(c3)}\n#define C4 {fx(c4)}\n")
print(open("tabf.h").read()[-160:])
