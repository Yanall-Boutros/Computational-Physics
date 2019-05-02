import numpy as np
import matplotlib.pyplot as plt
def weights(dim, nweights):
    ws = []
    if dim == 1:
        for i in range(nweights):
            if i == 0 or i == nweights-1:
                ws.append(1)
            else:
                ws.append(2)
        return np.array(ws)
    else:
        for i in range(nweights):
            if i == 0 or i == nweights-1:
                ws.append(weights(dim-1, nweights))
            else:
                ws.append(2*(weights(dim-1, nweights)))
    return(np.array(ws))
def SevenDTrapInt(func_pntr, dim_bounds, step):
    # We always evaluate the func_pntr at the starting boudnaries
    # We iterate forward in steps of step
    # Count in binary to find all points to be summed
    w = weights(7, 5)
    endpoints = np.ones(2**len(dim_bounds))
    domains = []
    for boundary in dim_bounds:
        b = boundary[1]
        a = boundary[0]
        domains.append(np.linspace(a, b, (b-a)/step))
    domains = np.array(domains)
    total = 0
    # SevenD implies 7 for loops
    for a_i, a in enumerate(domains[0]):
        for b_i, b in enumerate(domains[1]):
            for c_i, c in enumerate(domains[2]):
                for d_i, d in enumerate(domains[3]):
                    for e_i, e in enumerate(domains[4]):
                        for f_i, f in enumerate(domains[5]):
                            for g_i, g in enumerate(domains[6]):
                                val = func_pntr(a, b, c, d, e, f, g)
                                val *= (step**6/2**7) 
                                val *= w[a_i][b_i][c_i][d_i][e_i][f_i][g_i]
                                total += val
    return total
def SevenDPoly(x1, x2, x3, x4, x5, x6, x7):
    return (x1+x2+x3+x4+x5+x6+x7)**2
bounds = [(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1)]
plsdeargodwork = SevenDTrapInt(SevenDPoly, bounds, 0.2)
print(plsdeargodwork)
# Compare to montecarlo integration
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
