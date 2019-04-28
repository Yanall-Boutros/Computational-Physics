import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as scint
# Generalized rk2 ODE solver
h = 0.0001
def rk2n(func_pntr, func_args=[0, 0], init_cond=[0, 0]):
    # y_n+1 = y_n + k_2
    # k_2 = hf(t_n + h/2, y_n + k_1/2)
    # k_1 = hf(t_n, y_n)
    k_1 = h*func_pntr(*func_args)
    func_args[0] += h/2
    k_2 = h*func_pntr(*func_args)
    return init_cond[1] + (k_1 + 2*k_2)/3

def sin(w, x):
    # Assume \omega to equal 1
    return np.sin(w*x)

def scisin(y, t):
    # Returns the derivative of y = cos at a given y, t
    return np.sin(t)
w = 1
t_0 = 0
args = [w, t_0]
in_cond = [t_0, 0]
approx = [rk2n(sin, args, in_cond)]
t = 0
tvals = [t]

while t < 40*np.pi:
    t += h
    tvals.append(t)
    args = [w, tvals[-1]]
    in_cond = [tvals[-1], approx[-1]]
    approx.append(rk2n(sin, args, in_cond))

plt.plot(tvals, approx, label="rk2", alpha=0.5)

# Compare to scipy odeint
init_cond = 0
y = scint.odeint(scisin, init_cond, tvals)
plt.plot(tvals, y, label="Scipy", alpha=0.5)
plt.legend(loc=1)
plt.xlabel("x")
plt.ylabel("integration of sin(w,x)")
plt.savefig("tvals.png")
plt.figure()

# Part b, solve with f(t, y) = 1 + y^2 + t^3
def fty(t, y):
    return 1 + y**2 + t**3

def scifty(y, t):
    return 1 + y**2 + t**3
# Init variables
t = 1
btvals = [t]
bargs = [1, -4]
init_cond = [1, -4]
# Init conditions same as init args...
bapprox = [-4]

while t < 2:
    t += h
    btvals.append(t)
    bargs = [btvals[-1], bapprox[-1]]
    init_cond = bargs
    bapprox.append(rk2n(fty, bargs, init_cond))

plt.plot(btvals, bapprox, label="rk2n", alpha=0.5)
# Compare to scipy odeint
y = scint.odeint(scifty, y0=-4, t=btvals)
plt.plot(btvals, y, label="Scipy", alpha=0.5)
plt.legend(loc=1)
plt.xlabel("t")
plt.ylabel("integration of f(t,y)")
plt.savefig("b.png")
