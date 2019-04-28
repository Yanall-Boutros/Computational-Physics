# ---------------------------------------------------------------------
# Import Statements
# ---------------------------------------------------------------------
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as scint
# ---------------------------------------------------------------------
# Constant Definitions
# ---------------------------------------------------------------------
k = 1
m = 1
t = np.linspace(-3*np.pi*np.sqrt(k/m), 3*np.pi*np.sqrt(k/m), 10000)
F_0 = 30
w = 8*np.pi
# ---------------------------------------------------------------------
# Function Definitions
# ---------------------------------------------------------------------
def V(x, p, k):
    return (k/p)*(x**p)

# Describe our function as y'' is a function of its derivatives
def ddotx(X, t, k, m, p):
    # m * x'' = F_ext - (deriv V)
    # with f_ext = 0, deriv V = kx**p-1...
    # mx'' = -kx**p-1
    # x'' = -k/m * x**p-1
    # Now convert into two first order diffs
    # Define y_2 = x'
    # Define y_1 = x
    # Define y_1' = x' = y_2
    # Define y_2' = x''
    # solve for y_1' and y_2'
    # The equation is now  y_2' = -k/m * y_1**p-1
    # In summary, the nth element is the nth derivative of X
    return [X[1], (-k/m)*X[0]**(p-1)]

def fext_ddotx(X, t, F_0, w, k, m, p):
    # m * x'' = F_ext - (deriv V)
    # with f_ext != 0, deriv V = kx**p-1...
    # mx'' = F_0sin(\omega t) - kx**p-1
    # x'' = F_0sin(omega t))/m -k/m * x**p-1
    # Now convert into two first order diffs
    # Define y_2 = x'
    # Define y_1 = x
    # Define y_1' = x' = y_2
    # Define y_2' = x''
    # solve for y_1' and y_2'
    # The equation is now  y_2' = -k/m * y_1**p-1 + F_0sin(omega t)/m)
    # In summary, the nth element is the nth derivative of X
    return [X[1], (-k/m)*X[0]**(p-1) + (F_0/m)*np.sin(w*t)]
# ---------------------------------------------------------------------
# Main Loop
# ---------------------------------------------------------------------
for p in [2, 4, 8]: 
    y = scint.odeint(ddotx, [0, 1], t, (k, m, p))
    v = V(t, p, 1)
    # Graph of phasespace
    plt.figure()
    plt.plot(y[:,0], y[:,1])
    plt.title("x vs $\dot{x}$ for p = " + str(p))
    plt.xlabel("$\dot{x}$")
    plt.ylabel("x")
    plt.savefig(str(p)+"phasespace.pdf")
    # Graph of potential
    plt.figure()
    plt.plot(t, v)
    plt.title("V(x) vs x for p = "+str(p))
    plt.xlabel("x")
    plt.ylabel("V(x)")
    plt.savefig(str(p)+"pot.pdf")
    # Graph of xdot
    plt.figure()
    plt.plot(t, y[:,0])
    plt.xlabel("Time")
    plt.ylabel("$\dot{x}$")
    plt.title("$\dot{x}$ vs time for p = " + str(p))
    plt.savefig(str(p)+"y0.pdf")
    # Graph of x
    plt.figure()
    plt.plot(t, y[:,1])
    plt.title("x vs time for p = " + str(p))
    plt.xlabel("Time")
    plt.ylabel("x")
    plt.savefig(str(p)+"y1.pdf")
    # Graph of energy 
    plt.figure()
    Edata = y[:,1]**2 + y[:,0]**2
    plt.plot(t, Edata)
    plt.title("Graph of 2E for p = " + str(p))
    plt.xlabel("Time")
    plt.ylabel("2E")
    plt.savefig(str(p)+"e.pdf")
# Main loop but with F_ext != 0
for F_0 in [30, 1]:
    for p in [2, 4]:
        y = scint.odeint(fext_ddotx, [0, 1], t, (F_0, w, k, m, p))
        # Graph of phase space
        plt.figure()
        plt.plot(y[:,0], y[:,1])
        plt.title("x vs $\dot{x}$ for p = " + str(p) + " and F_0 = " + str(F_0))
        plt.xlabel("$\dot{x}$")
        plt.ylabel("x")
        plt.savefig(str(F_0)+str(p)+"phasespacef0.pdf")
        # Graph of x
        plt.figure()
        plt.plot(t, y[:,1])
        plt.title("x vs time for p = " + str(p)+ " and F_0 = " + str(F_0))
        plt.xlabel("Time")
        plt.ylabel("x")
        plt.savefig(str(F_0)+str(p)+"f0.pdf")
# Plot the max amp as a function of the driving frequency
omega = np.linspace(0, 20, 2000)
for p in [2, 4]:
    maxamp = []
    for w in omega:
        y = scint.odeint(fext_ddotx, [0, 1], t, (F_0, w, k, m, p))
        maxamp.append(max(y[:, 1]))
    plt.figure()
    plt.plot(omega, maxamp)
    plt.title("omega max amp for p = " + str(p))
    plt.xlabel("$\omega$")
    plt.ylabel("Max Amplitidue")
    plt.savefig(str(p)+"omegavsmaxamp.pdf")
    
