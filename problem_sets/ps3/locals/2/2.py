import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as scint
step = 0.0001
cprimevals = [0] # Maintain a list of vals where the most recent 
sprimevals = [0] # item in the list is area under the curve from
# zero up until the most recently updated boundary value
def trap_int(func_pntr, a, b, step, func_flag=""):
    total = 0
    domain = np.linspace(a, b, (b-a)/step) 
    # Area of a trapezoid is func(a) + func(b) / 2) * step
    if func_flag != "":
        if func_flag == "Cprime":
            total = cprimevals[-1]
        else:
            total = sprimevals[-1]
    for val in range(len(domain)-1): # Stop 1 step early for final b val
        left_endpoint = func_pntr(domain[val])
        right_endpoint = func_pntr(domain[val+1])
        trap_area = ((left_endpoint + right_endpoint)/2)
        trap_area *= (domain[val+1]-domain[val])
        if val == 0 or val == len(domain)-2:
            trap_area *= 0.5
        total += trap_area
    if func_flag == "Cprime":
        cprimevals.append(total)
    if func_flag == "Sprime":
        sprimevals.append(total)
    return total

def Cprime(v):
    return np.cos(np.pi * (v**2)/2)

def Sprime(v):
    return np.sin(np.pi * (v**2)/2)

def I(a, v):
    first_term = trap_int(Cprime, a, v, step, "Cprime") + 0.5
    first_term *= first_term

    second_term = trap_int(Sprime, a, v, step,"Sprime") + 0.5
    second_term *= second_term

    return 0.5*(first_term + second_term)

vdomain = np.linspace(0, 75, 10000)
I_output = []
a = 0
for v in range(1, len(vdomain)-1):
    I_output.append(I(vdomain[v-1], vdomain[v]))
    if v == 1070:
        plt.plot(vdomain[:1070], I_output, linewidth=0.5)
        plt.xlabel("v")
        plt.ylabel("I(v)")
        plt.title("Plot of I(v) vs v up to v=8")
        plt.savefig("v=8I.pdf")
plt.figure()
plt.xlabel("v")
plt.ylabel("I(v)")
plt.title("Plot of I(v) vs v up to v=75")
plt.plot(vdomain[1:-1], I_output, linewidth=0.25)
plt.savefig("I.pdf")
