import numpy as np
import matplotlib.pyplot as plt
# Starter code taken from Textbook (Landau, Paez, Bordenianu)
xmax = 40
xmin = 0.25
step = 0.1
order = 10
start = 100
def down(x, n, m):
    j = np.zeros(start+2, np.float128)
    j[m+1] = j[m] = 1
    for k in range(m, 0, -1):
        j[k-1] = np.float(((2*k+1) /x)*j[k] - j[k+1])
    scale = np.float128(np.sin(x)/x/j[0])
    return j[int(n)]*scale

def up(x, n, m):
    j = np.zeros(start+2, np.float128)
    j[m+1] = j[m] = 1
    for k in range(m, start+1, 1):
        j[k+1] = np.float128(((2.*k+1.) /x)*j[k] - j[k-1])
    scale = np.float128(np.sin(x)/x/j[0])
    return j[int(n)]*scale

def diff(u, d):
    return np.abs(u - d) / (np.abs(u) + np.abs(d))
print("l \t \t x \t \t j_l(x) down \t \t \t j_l(x) up \t \t \t % diff")
for l in np.linspace(1, 25, 25):
    for x in [0.1, 1, 10]:
        d = down(x, l, start)
        u = up(x, l, 0)
        percdiff = diff(u, d)
        print(l, "\t", x, "\t", d, "\t", u, "\t", percdiff)
# Generate Graphs of first 100 bessel functions
xdata = np.linspace(-95, 95, 1000)
ysum = np.zeros(len(xdata))
plt.figure()
for l in np.linspace(0, 99, 100):
    ydata = []
    for x in xdata:
        ydata.append(down(x, l, start))
    ysum += ydata
    plt.plot(xdata, ydata)
plt.title("First 100th order bessel functions")
plt.xlabel("x")
plt.ylabel("$J_l(x)$")
plt.savefig("bessel.pdf")
plt.figure()
plt.plot(xdata, ysum)
plt.title("Sum of the first 100 order bessel functions")
plt.xlabel("x")
plt.ylabel("$\Sigma J_l(x)$")
plt.savefig("Sum.pdf")
