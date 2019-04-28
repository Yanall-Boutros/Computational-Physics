import matplotlib.pyplot as plt
import numpy as np
import scipy.integrate as scint
import random
# Start code grabbed from 
# https://en.wikipedia.org/wiki/Lorenz_system#Python_simulation
from mpl_toolkits.mplot3d import Axes3D

rho = 28.0
sigma = 10.0
beta = 8.0 / 3.0
delta = 0.000000000000001

def f(state, t):
  x, y, z = state  # unpack the state vector
  return sigma * (y - x), x * (rho - z) - y, x * y - beta * z  # derivatives

def modify_state(state):
    # Increment just one of the state params
    state[0] += delta
    # Then shuffle
    random.shuffle(state)

state0 = [1.0, 1.0, 1.0]
t = np.arange(0.0, 40.0, 0.01)

states = scint.odeint(f, state0, t)

fig = plt.figure()
ax = fig.gca(projection='3d')
ax.plot(states[:,0], states[:,1], states[:,2],
        label="Function from Initial State")
ax.scatter(states[-1,0], states[-1,1], states[-1,2],
            c='r', marker='o', s=1, label="End Point")
ax.scatter(states[0,0], states[0,1], states[0,2],
            c='g', marker='x', s=3, label="Starting Point")
for i in range(1000):
    ax.scatter(states[-1,0], states[-1,1], states[-1,2],
               c='r', marker='o', s=1)
    # Now modify state a little bit and replot
    modify_state(state0)
    states = scint.odeint(f, state0, t)
plt.legend(loc=3)
plt.title("Lorentz Attraction")
print("Final Initial State")
print(state0)
plt.savefig("Wiki.pdf")
