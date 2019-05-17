import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt
from mpl_toolkits.mplot3d import Axes3D 
# improve (or rewrite) the code so that the relaxation iteration stops when
# the solution has converged to some tolerance. This means that the grid points
# all change by less than some amount between iterations. Plot the number of
# iterations needed for a tolerance of 10−4, 10−3, 10−2, and 10−1.
def functz(V):  
    # Function returns V(x, y)
    z = V[X,Y]                        
    return z
Nmax = 150
x = range(0, Nmax)
y = range(0, Nmax)
X, Y = np.meshgrid(x,y)          
max_delta = 1 # Assume max delta is large so that we keep iterating
needed_iterations = []
V = np.ones((Nmax, Nmax))   
V[0,:] = 100.0 
V[-1,:] = 100.0
V[:,0] = 100.0
V[:,-1] = 100.0
# The length / width is 150x150. The square whole is a third this value
# Or at 50x50 located in the center. That means at a point , 75x75
# our mesh spans , +- 25 and is initated to zero
iterrange = np.ones((Nmax, Nmax))
iterrange[:,:] = True
iterrange[0,:] = False
iterrange[-1,:] = False
iterrange[:,0] = False
iterrange[:,-1] = False
iterrange[49:99, 49:99] = False
iterations = 0 # count the iterations. Initalize outside for loop so it stays
# updated on tolerance change.
trange = [0.01] # this is fine
for tolerance in trange:
    while max_delta > tolerance and iterations < 1000: # Keep looping until
        # acceptable delta reached or timeout.
        max_delta = 0 # reinit to 0 so we can see max change per iter
        for i in range(1, len(iterrange)-1):                                                
            for j in range(1, len(iterrange)-1):
                #if iterrange[i][j] == True:
                old = V[i,j]
                V[i,j] = 0.25*(V[i+1,j]+V[i-1,j]+V[i,j+1]+V[i,j-1])  
                    # Using np.concactenate and np.sum you can vectorize
                    # the above code, I'm too lazy to do that as of now however
                    # and choose to let my code run for 6 minutes instead
                delta = np.abs(V[i, j] - old)
                if delta > max_delta:
                    max_delta = delta
        iterations += 1
        if iterations % 100 == 0:
            print (iterations/10, "%.........................................")
        V[49:99, 49:99] = 0 

Z = functz(V)                          
fig = plt.figure()                              # Create figure
ax = Axes3D(fig)                                # plot axes
ax.plot_wireframe(X, Y, Z, color = 'r', rstride=5, cstride=3)  # red wireframe, setting the stride values explicitly
ax.set_xlabel('X')                              # label axes
ax.set_ylabel('Y')
ax.set_zlabel('Potential U')
plt.savefig("pot.pdf")

plt.figure()
plt.title("Heatmap of potential in $X, Y$")
plt.xlabel("X coordinate")
plt.ylabel("Y coordinate")
plt.contour(Z, cmap="jet", levels=998)
plt.savefig("cont.pdf")

plt.figure()
plt.title("Gradient of Z")
gradZ = -1*np.array(np.gradient(Z))
plt.streamplot(X, Y, gradZ[1], gradZ[0])
plt.xlabel("X")
plt.ylabel("Y")
plt.savefig("grad.pdf")
