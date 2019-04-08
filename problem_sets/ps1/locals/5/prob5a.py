import numpy as np
import matplotlib.pyplot as plt
# Starter code taken from Textbook (Landau, Paez, Bordenianu)
xmax = 40
xmin = 0.25
step = 0.1
order = 10
start = 50
def down(x, n, m):
    j = np.zeros(start+2, np.float128)
    j[m+1] = j[m] = 1
    for k in range(m, 0, -1):
        j[k-1] = ((2*k+1) /x)*j[k] - j[k+1]
    scale = np.sin(x)/x/j[0]
    return j[int(n)]*scale

def up(x, n, m):
    j = np.zeros(start+2, np.float128)
    j[m+1] = j[m] = 1
    for k in range(m, start+1, 1):
        j[k+1] = ((2.*k+1.) /x)*j[k] - j[k-1]
    scale = np.sin(x)/x/j[0]
    return j[int(n)]*scale

print("l \t \t x \t \t j_l(x) down \t \t j_l(x) up")
for l in np.linspace(1, 25, 25):
    for x in [0.1, 1, 10]:
        print(l, "\t \t", x, "\t \t", down(x, l, start), "\t"
                ,up(x, l, 0))
