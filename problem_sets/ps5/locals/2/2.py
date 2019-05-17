import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt
from mpl_toolkits.mplot3d import Axes3D 
def functz(V):  
    # Function returns V(x, y)
    z = V[X,Y]                        
    return z


# improve (or rewrite) the code so that the relaxation iteration stops when
# the solution has converged to some tolerance. This means that the grid points
# all change by less than some amount between iterations. Plot the number of
# iterations needed for a tolerance of 10−4, 10−3, 10−2, and 10−1.
Nmax = 100
x = range(0, Nmax)
y = range(0, Nmax)
X, Y = np.meshgrid(x,y)          
def c_f_sol(x, y):
    L = 100
    total = 0
    n = 1
    while n < 200 or np.abs(term) > 0.00001:
        term = np.sin(n*np.pi*x/L) * 400/(n*np.pi) 
        term *= np.sinh(n*np.pi*y/L)
        term /= np.sinh(n*np.pi)
        total += term
        n += 2
    return total
# Analytically
x = np.linspace(0, Nmax-1, Nmax)
y = x.copy()
vprime = np.zeros((Nmax, Nmax))
for x_a in x:
    for y_a in y:
        vprime[int(y_a), int(x_a)] = c_f_sol(y_a, x_a)
Z = functz(vprime)                          
fig = plt.figure()                              # Create figure
ax = Axes3D(fig)                                # plot axes
ax.plot_wireframe(X, Y, Z, color = 'r', rstride=5, cstride=1)  # red wireframe, setting the stride values explicitly
ax.set_xlabel('X')                              # label axes
ax.set_ylabel('Y')
ax.set_zlabel('Potential U (analytic)')
plt.savefig("pot_analytic.pdf")
# Not analytically
max_delta = 1 # Assume max delta is large so that we keep iterating
needed_iterations = []
V = np.zeros((Nmax, Nmax))   
V[0,:] = 100.0 
iterations = 0 # count the iterations. Initalize outside for loop so it stays
# updated on tolerance change.
trange = [1, 0.1, 0.01, 0.001]
for tolerance in trange:
    while max_delta > tolerance: # Keep looping until acceptable delta reached
        max_delta = 0 # reinit to 0 so we can see max change per iter
        for i in range(1, Nmax-1):                                                
            for j in range(1,Nmax-1):
                old = V[i,j]
                V[i,j] = 0.25*(V[i+1,j]+V[i-1,j]+V[i,j+1]+V[i,j-1])  
                # Using np.concactenate and np.sum you can vectorize
                # the above code, I'm too lazy to do that as of now however
                # and choose to let my code run for 2 minutes instead
                delta = np.abs(V[i, j] - old)
                if delta > max_delta:
                    max_delta = delta
        iterations += 1
    needed_iterations.append(iterations)

Z = functz(V)                          
fig = plt.figure()                              # Create figure
ax = Axes3D(fig)                                # plot axes
ax.plot_wireframe(X, Y, Z, color = 'r', rstride=5, cstride=1)  # red wireframe, setting the stride values explicitly
ax.set_xlabel('X')                              # label axes
ax.set_ylabel('Y')
ax.set_zlabel('Potential U')
plt.savefig("pot.pdf")

plt.figure()
plt.plot(trange, needed_iterations)
plt.xlabel("Tolerance")
plt.ylabel("Numer of iterations")
plt.title("Number of iterations for potential change to be under tolerance vs tolerance")
plt.savefig("n_it.pdf")

plt.figure()
plt.title("Heatmap of potential in $X, Y$")
plt.xlabel("X coordinate")
plt.ylabel("Y coordinate")
plt.contour(Z, cmap="jet", levels=998)
plt.savefig("cont.pdf")

