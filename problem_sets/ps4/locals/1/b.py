import numpy as np
import matplotlib.pyplot as plt
def weights(dim, nweights):
    # Originally I tried to accomplish this program by writing a
    # function which could solve any N dimensional trapezoidal 
    # integration. Unfortunately I ran out of time and had to
    # switch to solving only this specific 7D case.
    # To generalize trapezoidal integration, recursion is necessary
    # which is used in the weights function (the recursive call is
    # in the nested if else statements inside the outer most 
    # else statement
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
