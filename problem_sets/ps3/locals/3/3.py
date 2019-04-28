import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as scint
E = np.linspace(0, 180, 10000)
M = 90
G = 10
def BW(E, M, G):
    g = np.sqrt(M**2*(M**2 + G**2))
    k = (2*np.sqrt(2)*M*G*g)
    k /= (np.pi * np.sqrt(M**2 + g))
    den = (E**2 - M**2)**2 + (M*G)**2
    return k/den
def max_x_y(x, y):
    target = max(y)
    for i in range(len(x)):
        if y[i] == target:
            print("Target: ", target, "\nLocated at index: ",
                  i, "\nOf X-value: ", x[i])
            return (i, x[i], y[i])
fE = BW(E, M, G)
plt.plot(E, fE)
plt.title("Breit Wigner Distribution")
plt.xlabel("Energy in GeV")
plt.ylabel("Width of Resonance")
plt.savefig("BW.pdf")
# Part B
# Find domain of integration is +/- 3G away from xmax
index, xmax, ymax = max_x_y(E, fE)
A = xmax - G*3
B = xmax + G*3
N = 1000
h = (B-A)/(N-1)
w = np.ones(N, np.float128)*(h/3)
for i in range(1, N-1):
    w[i] *= 2
    if i % 2 == 1: w[i] *= 2
sums = 0
wsums = 0
for i in range(1, N+1):
    x = A + (i-1)*h
    sums += w[i-1]*BW(x, M, G)
    wsums += w[i-1]
print("Area under curve = ", sums)
print(wsums)
