import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as scint
step = 0.0001
def trap_int(func_pntr, a, b, step):
    total = 0
    domain = np.linspace(a, b, (b-a)/step) 
    # Area of a trapezoid is func(a) + func(b) / 2) * step
    for val in range(len(domain)-1): # Stop 1 step early for final b val
        left_endpoint = func_pntr(domain[val])
        right_endpoint = func_pntr(domain[val+1])
        trap_area = ((left_endpoint + right_endpoint)/2)
        trap_area *= (domain[val+1]-domain[val])
        if val == 0 or val == len(domain)-2:
            trap_area *= 0.5
        total += trap_area
    return total

def Cprime(v):
    return np.cos(np.pi * (v**2)/2)

def Sprime(v):
    return np.sin(np.pi * (v**2)/2)

def I(a, v):
    first_term = trap_int(Cprime, a, v, step) + 0.5
    first_term *= first_term

    second_term = trap_int(Sprime, a, v, step) + 0.5
    second_term *= second_term

    return 0.5*(first_term + second_term)

vdomain = np.linspace(0, 75, 10000)
I_output = []
I_output.append(I(0, vdomain[1]))
a = 0
for v in range(1, len(vdomain)):
    I_output.append(I(vdomain[v-1], vdomain[v]))
    if v == 1070:
        plt.plot(vdomain[:1071], I_output, linewidth=0.5)
        plt.xlabel("v")
        plt.ylabel("I(v)")
        plt.title("Plot of I(v) vs v up to v=8")
        plt.savefig("v=8I.pdf")
plt.figure()
plt.xlabel("v")
plt.ylabel("I(v)")
plt.title("Plot of I(v) vs v up to v=75")
plt.plot(vdomain, I_output, linewidth=0.25)
plt.savefig("I.pdf")
