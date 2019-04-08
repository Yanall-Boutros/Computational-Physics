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
print("sin(0) is roughly: ",sin(0,1)) # Test
print("sin(10) is roughly: ",sin(10,10))
xaxis = []
yaxis = []
for x in [1, 10, 50, 100]:
    for n in range(1, 200):
        yaxis.append(sin(x,n))
        xaxis.append(n)
    plt.figure()
    plt.plot(xaxis, yaxis)
    plt.title("x="+str(x))
    plt.xlabel("term $n$")
    plt.ylabel("nth term in series solution for $sin(x)$")
print("sin(100) is (not roughly): ", totsin(100))
