import numpy as np
import scipy as sp
import math as m
import matplotlib.pyplot as plt
def mac_sin(x, N):
    # Calculate sin(x) from the maclaurin series up to N terms
    term = []
    for n in range(1, N+1):
        num = (x)**(2*n-1) * (-1)**(n-1)
        den = sp.math.factorial(2*n - 1)
        term.append(num / den)
    return np.array(term)
def ntherm(x,N,n):
    term = mac_sin(x, N)
    return term[n]
yvals = []
xvals = []
for x in [1, 10, 50]:
    for N in range(2, 100):
        xvals.append(N)
        yvals.append(ntherm(x,N,N-1))
    plt.figure()
    plt.plot(xvals, yvals, linestyle="none", marker=".")
    plt.title("Approximation of sin(x) for x=" + str(x))
    plt.xlabel("term n")
    plt.ylabel("nth term from maclaurin series")
    plt.savefig(str(x)+"_p4a.pdf")
