import numpy as np
import matplotlib.pyplot as plt
# As seen in lecture:
# Ising model implementation with Metropolis algorithm
# Based on IsingViz.py from Landau, et al.

N     = (2000, 2000) # number of spin sites                                   
num_steps = 20*N[0]*N[1]  # number of iterations
B     = 0.05         # magnetic field                               
mu    = .33          # g mu (not needed if B=0)
J     = 1.           # exchange energy                              
k     = 1.
t     = 1.
np.random.seed() 
state0 = -1*np.ones(N) # Start with an arbitrary spin configuration
Evals = []
def energy(State, J, mu, B):
    # Energy will call row_energy on every row in State, and on 
    # every row of its transpose.
    total_energy = 0
    for row in State:
        total_energy += row_energy(row, J, mu, B)
    for col in State.transpose():
        total_energy += row_energy(col, J, mu, B)
    return total_energy
def row_energy(S, J, mu, B):
    first_set = np.concatenate([np.array([S[-1]]), S[:-1]])
    FirstTerm = np.sum(-J*first_set[:-1]*first_set[1:])
    SecondTerm = np.sum(-mu*S*B)
    return (FirstTerm + SecondTerm)

def energy_change(S, coor):
    """Determine the change in energy if spin `i` is flipped
    
    `(i+1) % len(S)` implements the periodic boundary condition.
    """
    x = coor[0]
    y = coor[1]
    # Multiply spin site by all adjacent elements
    S_left  = S[x-1,y]
    S_right = S[(x+1) % len(S),y]
    S_up    = S[x,(y+1) % len(S)]
    S_down  = S[x,y-1]
    return 2*J*S[x,y]*(S_left + S_right + S_up + S_down % len(S)) + 2*B*mu*S[x,y]

def TwoDIsing(state0, num_steps, J, mu, B, kT):
    ES = energy(state0, J, mu, B)
    energy_values = []
    energy_values.append(ES)
    # Contains a copy of the state configuration so we don't have to store
    # 2**N**2 * # of time step elements
    state_configs = np.array([state0, state0])
    deltas = [] # A lighter way of keeping track of how the state changes.
    randx = np.random.randint(2000, size=num_steps)
    randy = np.random.randint(2000, size=num_steps)
    count = 1
    for x, y in np.stack(randx, randy):
        #test_state = state_configs[-1]
        # Trial step: flip spin at one random site
        #test_state[x,y] *= -1.
        """Our key modification is here: instead of calculating the energy
        afresh, we only calculate the change therein.
        """
        state_configs[-1][x,y] *= -1
        ET = ES + energy_change(state_configs[-1], (x,y))
        R = np.exp((ES-ET)/(kT))           # Boltzmann test
        if R > np.random.random():
            #state_configs[-1] = test_state      # replace the state, or
            ES = ET
            deltas.append((x,y))
        else:
            # advance the previous state forward
            state_configs[-1]=state_configs[-2]
            deltas.append(())
        energy_values.append(ES)
        count += 1
        if count % 1000 == 0: print((count / num_steps)*100," %.....................")
    return state_configs, energy_values, deltas
