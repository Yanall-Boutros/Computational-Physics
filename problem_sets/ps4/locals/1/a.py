# ---------------------------------------------------------------------
# Import Statements
# ---------------------------------------------------------------------
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as scint
# ---------------------------------------------------------------------
# Function Definitions
# ---------------------------------------------------------------------
def TenDPoly(x1, x2, x3, x4, x5, x6, x7, x8, x9, x10):
    return (x1+x2+x3+x4+x5+x6+x7+x8+x9+x10)**2
# I \approx Q_N \equiv \frac{V}{N}\Sum_i^N \bar{f(x_i)}
def MCIntegration(func_pntr, dim_bounds, N_points):
    # Boudaries should be immutable. Therefore dim_bouds should be
    # an array of tuples
    if type(dim_bounds) != type(list()):
        print("Boundaries must be passed as list of tuples")
    if type(dim_bounds[0]) != type(tuple()):
        print("Boundary components must be tuples")
    if len(dim_bounds[0]) != 2:
        print("Tuples must be length(2)")
    # The dimensionality of a function is the size of its arguments.
    # For example, if f depends on x, y, it's dimensionality is 2.
    # And has area ~ R^2
    # The volume is then the boundary lengths multiplied by each other
    V = 1
    for int_bounds in dim_bounds:
        V *= (int_bounds[1] - int_bounds[0])
    # Now, a function f which takes D = len(dim_bounds)) independent
    # arguments is summed together over N points
    D = len(dim_bounds)
    total = 0
    for i in range(N_points):
        # Must calculate func_pntr with D random arguments.
        total += func_pntr(*np.random.rand(D))
    return V*total/N_points 
# ---------------------------------------------------------------------
# Initialization
# ---------------------------------------------------------------------
bounds = [(0,1), (0,1), (0,1), (0,1), (0,1),
          (0,1), (0,1), (0,1), (0,1), (0,1),]
results = []
nvals = 2**np.arange(5, 19)
answer = np.float128(155/6)
# ---------------------------------------------------------------------
# Main Loop
# ---------------------------------------------------------------------
for N in nvals:
    results.append(MCIntegration(TenDPoly, bounds, N))
# Plot the error, the correct answer is 155/6
results = np.array(results)
error = (np.abs(results-answer)/answer)*100
plt.loglog(np.sqrt(nvals),error)
plt.title("Results of Monte Carlo Integration on log-log scale")
plt.xlabel("$\sqrt{N}$")
plt.ylabel("$\% Error$")
plt.savefig("mc.pdf")
