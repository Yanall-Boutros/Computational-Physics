import numpy as np
import matplotlib.pyplot as plt
answers = dict() # For dynamics programming purposes, stores already calculated answers
def totsin(x):
    # Returns the approximation of sin(x)
    n = 1
    tot = 0
    while answers.get((x,n)) is not None:
        tot += answers[(x,n)]
        n += 1
    return tot
def sin(x, n):
    # function calculates the first n terms for sin(x)
    if n == 1: # Base Case
        answers[(x,n)] = x # Update the book of answers
        return x # sin(x) \approx x for n=1
    elif answers.get((x,n)) is not None: return answers[(x, n)] # Don't calculate anything we've aready calculated...
    else: # Recursive Case
        entry = ((-x**2)/((2*n - 1)*(2*n - 2))) * sin(x, n-1) # recursive call
        answers[(x,n)] = entry
        return entry

xdata = np.linspace(-10,10,500)
for N in np.linspace(1, 15, 15):
    ydata = []
    yexact = []
    for x in xdata:
        sin(x, N)
        ydata.append(totsin(x))
        yexact.append(np.sin(x))
    ydata = np.array(ydata)
    yexact = np.array(yexact)
    yflip = -1*ydata
    yexactflip = -1*yexact
    plt.plot(xdata, ydata, alpha = 0.5, label="Maclaurin Series Approximation to N = " + str(N))
    plt.plot(xdata, yflip, alpha = 0.5, label="Maclaurin Series Approximation to N = " + str(N))
    plt.plot(xdata, yexact, alpha = 0.5, label="Numpy Answer", color='red')
    # plt.legend(loc=1)
    plt.grid(True)
    plt.xlim(-11, 11)
    plt.ylim(-11, 11)
plt.savefig("4c.pdf")
