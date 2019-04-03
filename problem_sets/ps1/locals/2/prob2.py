import numpy as np
import matplotlib.pyplot as plt
def x(a, b, c):
    # Find the roots of ax^2 + bx + c = 0
    sol1 = (-b + (b**2 - 4*a*c)**0.5)/(2*a)
    sol2 = (-b - (b**2 - 4*a*c)**0.5)/(2*a)
    return (sol1, sol2)

def xprime(a, b ,c):
    # Find the roots of ax^2 + bx + c = 0 via the second equation
    sol1 = -2*c/(b + (b**2 - 4*a*c)**0.5)
    sol2 = -2*c/(b - (b**2 - 4*a*c)**0.5)
    return (sol1, sol2)

def xact(a, b, c): # get it? xact? Exact?
    a = np.float128(a)
    b = np.float128(b)
    c = np.float128(c)
    sol1 = -2*c/(b + (b**2 - 4*a*c)**0.5)
    sol2 = -2*c/(b - (b**2 - 4*a*c)**0.5)
    return (sol1, sol2)

def err(xact, notexact):
    return np.abs((xact - notexact))/np.abs(xact)

def errplot(errset, i):
    n = np.linspace(1, len(errset), len(errset))
    plt.plot(n, np.log10(errset))
    if i == 0:
        plt.title("Percent Error for Each Root determined from function x as $a = b = 1$, $c = 10^{-n}$")
        plt.xlabel("n")
        plt.ylabel("Percent Error of the Roots")
    elif i == 1:
        plt.title("Percent Error for Each Root determined from function xprime as $a = b = 1$, $c = 10^{-n}$")
        plt.xlabel("n")
        plt.ylabel("Percent Error of the Roots")
    elif i == 2:
        plt.title("Percent Error for Each Root determined from function x as $a = 10^{-3}$, $b = 1$, $c = 10^{-n}$")
        plt.xlabel("n")
        plt.ylabel("Percent Error of the Roots")
    else:
        plt.title("Percent Error for Each Root determined from function xprime as $a = 10^{-3}$, $b = 1$, $c = 10^{-n}")
        plt.xlabel("n")
        plt.ylabel("Percent Error of the Roots")
    plt.savefig(str(i)+'.pdf')
# First run with a = b = 1
a = 1
b = 1
xset = []
xprimeset = []
xactset = []
for n in range(1, 14): # counting starts at 0
    c = 10**-n
    xset.append(x(a, b, c))
    xprimeset.append(xprime(a, b, c))
    xactset.append(xact(a, b, c))
    
# Second run with a = 10^-3, b = 1
a = 10**-3
b = 1
xsetrun2 = []
xprimesetrun2 = []
xactsetrun2 = []
for n in range(1, 14): # counting starts at 0
    c = 10**-n
    xsetrun2.append(x(a, b, c))
    xprimesetrun2.append(xprime(a, b, c))
    xactsetrun2.append(xact(a, b, c))

# Convert lists to np arrays
xset = np.array(xset)
xprimeset = np.array(xprimeset)
xactset = np.array(xactset)
xsetrun2 = np.array(xsetrun2)
xprimesetrun2 = np.array(xprimesetrun2)
xactsetrun2 = np.array(xactsetrun2)

# Calculate Error
errset = err(xactset, xset)
errprimeset = err(xactset, xprimeset)
errset2 = err(xactsetrun2, xsetrun2)
errprimeset2 = err(xactsetrun2, xprimesetrun2)

allsets = [errset, errprimeset, errset2, errprimeset2]
# Plot the log10 of the error with respect to n
for i, error_set in enumerate(allsets):
    plt.figure()
    errplot(error_set, i)
