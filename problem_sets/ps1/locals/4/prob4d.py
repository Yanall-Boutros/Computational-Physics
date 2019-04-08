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
def n_needed_terms(x):
    n = 1
    while np.abs(sin(x, n)) > 10**-8:
        n+=1
    return n
xdata = np.linspace(-10,10,500)
ydata = []
for x in xdata:
    ydata.append(n_needed_terms(x))
    print(ydata[-1])
plt.plot(xdata,ydata, marker=".", linestyle="none")
plt.savefig("4d.pdf")
