import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
# Parameters (in problem set 5, you'll add L and derive the x and t steps from
# N and CFL)
L = float(input("Enter a value for L:\n"))
N = int(input("Enter a value for N:\n"))
mu = float(input("Enter a value for mu:\n"))
ten = float(input("Enter a value for tension:\n"))
stability = float(input("Enter a value for stability:\n"))
c = np.sqrt(ten/mu)                                       # Propagation speed
deltaX = L/(N-1)
deltaT = deltaX/(stability*c)
c1 = deltaX/deltaT  # CFL criterion. You'll make better use of this in Problem set 5.
# c is deltax / deltaT
ratio =  c*c/(c1*c1)
# Initialization
# N x's & 3 t's to save memory (only need 3 times for algorithm)

xi = np.zeros( (N, 3), float)
k = np.linspace(0, N-1, N)
ipluck=0.5*L
def init():
    for i in range(0, int(ipluck*N/L)):
        xi[i, 0] = 0.1*i/int(ipluck*N)          # Initial condition: string plucked,shape
    for i in range (int(ipluck*N/L), N):                           # first part of string
        xi[i, 0] = 0.1 - 0.1/(N*(1-ipluck))*(i - int(ipluck*N))                 # second part of string

init()                                     # plot string initial position   
fig=plt.figure()                           # figure to plot (a changing line)
# select axis; 111: only one plot, x,y, scales given
ax = fig.add_subplot(111, autoscale_on=False, xlim=(0, N), ylim=(-0.15, 0.15))
ax.grid()                                                       # plot a grid
plt.title("Vibrating String with Asymmetric Pluck")
line, = ax.plot(k, xi[:,0], lw=2)             # x axis, y values, linewidth=2     

# Later time steps
for i in range(1, N-1):                                      # use  algorithm
    xi[i, 1] = xi[i, 0] + 0.5*ratio*(xi[i + 1, 0] + xi[i - 1, 0] - 2*xi[i, 0])   

def animate(num):               #num: dummy,  algorithm, will plot (x, xi)            
    for i in range(1, N-1):              
        xi[i,2] = 2.*xi[i,1]-xi[i,0]+ratio*(xi[i+1,1]+xi[i-1,1]-2*xi[i,1])
    line.set_data(k,xi[:,2])                              # data to plot ,x,y           
    for m in range (0,N):                               # part of algorithm
        xi[m, 0] = xi[m, 1]                               # recycle array 
        xi[m, 1] = xi[m, 2]
    return line,
# next: animation(figure, function,dummy argument: 1      
ani = animation.FuncAnimation(fig, animate,1)           
plt.show()             
print("finished")
