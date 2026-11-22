import sys, subprocess, mpmath as mp
mp.mp.prec = 200
B, d, scheme = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
out = subprocess.run(["python3","coeffs.py",str(d),str(B),str(B)],capture_output=True,text=True).stdout
line = [l for l in out.splitlines() if l.startswith("static")][0]
C = line.split("{")[1].split("}")[0].split(", ")
N = 1<<B
T = [float(mp.power(2, mp.mpf(i)/N)).hex() for i in range(N)]
with open(sys.argv[4],"w") as f:
    f.write(f"{B} {d} {scheme}\n" + "\n".join(C) + "\n" + "\n".join(T) + "\n")
print(out.splitlines()[0])
