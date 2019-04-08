import numpy as np
import matplotlib.pyplot as plt
def Sup(n):
    total = 0
    vals = np.float64(np.linspace(1, n, n))
    for num in vals:
        total += np.float64(num**-1)
    return total
def Sdown(n):
    total = 0
    vals = np.float64(np.linspace(n, 1, n))
    for num in vals:
        total += num**-1
    return total
def logcomp(up, down):
    up = np.float64(up)
    down = np.float64(down)
    if up > down: return (up - down)/(np.abs(up) + np.abs(down))
    return (down - up)/(np.abs(up) + np.abs(down))
comps = []
nvals = []
print("Summing up: \t \t Summing down:")
for n in range(2, 8):
    N = 10**n
    nvals.append(N)
    sup = Sup(N)
    sdown = Sdown(N)
    comps.append(logcomp(sup, sdown))
    print(sup, "\t", sdown)
comps = np.array(comps)
nvals = np.array(nvals)
plt.loglog(nvals, comps)
plt.title("Plot of Percent Error of Summations vs Number of Terms in Summation")
plt.xlabel("sum limit $N = 10^n$")
plt.ylabel("Sum of $n^{-1}$")
plt.savefig("3b.pdf")
