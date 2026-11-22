# Coefficients de 2^h ≈ 1 + h*q(h), q de degré d-1, par interpolation de Tchebychev de g(h)=(2^h-1)/h
# sur |h|<=a (quasi-minimax), puis erreur relative maximale vérifiée en multiprécision.
import mpmath as mp, sys
mp.mp.prec = 300
def g(h): return mp.expm1(h*mp.log(2))/h if h != 0 else mp.log(2)
def cheb(d, a):
    xs = [a*mp.cos(mp.pi*(2*k+1)/(2*d)) for k in range(d)]
    A = mp.matrix(d, d); b = mp.matrix(d, 1)
    for i,x in enumerate(xs):
        for j in range(d): A[i,j] = x**j
        b[i] = g(x)
    s = mp.lu_solve(A, b); return [s[j] for j in range(d)]
def relerr(c, a, n=4000, dbl=False):
    cc = [mp.mpf(float(x)) for x in c] if dbl else c
    m = 0
    for k in range(-n, n+1):
        h = a*k/n
        v = 1 + h*sum(cc[j]*h**j for j in range(len(cc)))
        m = max(m, abs(v/mp.power(2,h) - 1))
    return m
if __name__ == "__main__":
    d, B = int(sys.argv[1]), int(sys.argv[2]); a = mp.mpf(2)**(-(int(sys.argv[3]) if len(sys.argv)>3 else B+1))
    c = cheb(d, a)
    e = relerr(c, a, 1000, dbl=True)
    print("// deg=%d B=%d  erreur relative (coeffs arrondis en double) 2^%.2f" % (d, B, float(mp.log(e,2))))
    print("static const double C[] = {" + ", ".join(float(x).hex() for x in c) + "};")
