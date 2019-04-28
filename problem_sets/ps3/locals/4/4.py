import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as scint

def x_i(r_i):
    # \int_a^{x_i} w(x) dx = r_i
    # For w(x) = Ax, Ax^2/2 = r_i
    # Therefore x = (2r_i/A)^0.5
    return (r_i/.01)**0.5

def drawfour(x):
    # Draw four numbers from distribution, add them together,
    # Make new histogram from this
    sums = []
    for i in range(int((len(x)-4)/ 4)):
        a = x[4*i+0]
        b = x[4*i+1]
        c = x[4*i+2]
        d = x[4*i+3] 
        sums.append(a+b+c+d)
    return np.array(sums)

def gauss(x, sig, mu):
    a = np.sqrt(1/(2*np.pi*sig**2))
    b = np.exp(-((x-mu)**2)/(2*sig**2))
    return a*b
x = []

for i in range(10000000):
    x.append(x_i(np.random.rand()))
x = np.array(x)
plt.hist(x, bins=100)
plt.savefig("Ax.pdf")

sums = drawfour(x)
plt.figure()
h1 = plt.hist(sums, bins=1000)
plt.savefig("sums.pdf")


stdv = np.std(sums)
avg = np.average(sums)
t = np.linspace(0, 40, 100000)
y = gauss(t, stdv, avg)
y = np.array(y)
scale = max(h1[0]) / max(y)
y *= scale
plt.plot(t, y)
plt.savefig("Guass.pdf")
print("The shape of the summed distribution is a gaussian with",
        "Standard Deviation ", stdv, " and a mean of ", avg, "The ",
        "Histogram has a a peak frequency of ",
        h1[0][646], " for the sum equivalent to",
        h1[1][646])
