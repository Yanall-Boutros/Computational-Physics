# Author: Yanall Boutros
# Ising model implementation with Metropolis algorithm
# Based on IsingViz.py from Landau, et al.
# Starter code taken from Professor Ritz's pset4 solutions
# ==============================================================================
# Import Statements
# ==============================================================================
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
# ==============================================================================
# Function Definitions
# ==============================================================================
def calc_cv(E, kT, N):
    # specific heat is given by 1/N^2
    # times average (E^2) - (average E)^2
    # Divided by kT^2
    A = 1/(N**2)
    B = np.average(E**2)
    C = np.average(E)**2
    D = 1/(kT**2)
    return A*(B-C)*D

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

# Different configurations for what spin sites are identified as 
# nearest neighbors
def all_adjacents(S, coor):
    x = coor[0]
    y = coor[1]
    length = len(S)
    # Multiply spin site by all adjacent elements
    total_adjacent_spins = 0
    total_adjacent_spins += S[x-1,y]
    total_adjacent_spins += S[(x+1) % length,y]
    total_adjacent_spins += S[x,(y+1) % length]
    total_adjacent_spins += S[x,y-1]
    total_adjacent_spins += S[(x+1) % length, (y+1) % length]
    total_adjacent_spins += S[(x-1),(y+1) % length]
    total_adjacent_spins += S[(x-1),(y-1)]  
    total_adjacent_spins += S[(x+1) % length,(y-1)]
    return total_adjacent_spins

def grid_adjacents(S, coor):
    x = coor[0]
    y = coor[1]
    length = len(S)
    # Multiply spin site by all adjacent elements
    total_adjacent_spins = 0
    total_adjacent_spins += S[x-1,y]
    total_adjacent_spins += S[(x+1) % length,y]
    total_adjacent_spins += S[x,(y+1) % length]
    total_adjacent_spins += S[x,y-1]
    return total_adjacent_spins

def diag_adjacents_1(S, coor):
    x = coor[0]
    y = coor[1]
    length = len(S)
    # Multiply spin site by all adjacent elements
    total_adjacent_spins = 0
    total_adjacent_spins += S[x-1,y-1]
    total_adjacent_spins += S[(x+1) % length,(y+1)%length]
    total_adjacent_spins += S[x-1,(y+1) % length]
    total_adjacent_spins += S[(x+1)%length,y-1]
    return total_adjacent_spins

def diag_adjacents_2(S, coor):
    x = coor[0]
    y = coor[1]
    length = len(S)
    # Multiply spin site by all adjacent elements
    total_adjacent_spins = 0
    total_adjacent_spins += S[x-1,y-1]
    total_adjacent_spins += S[(x+1) % length,(y+1)%length]
    total_adjacent_spins += S[x-1,(y+1) % length]
    total_adjacent_spins += S[(x+1)%length,y-1]
    total_adjacent_spins += S[x-2,y-2]
    total_adjacent_spins += S[(x+2) % length,(y+2)%length]
    total_adjacent_spins += S[x-2,(y+2) % length]
    total_adjacent_spins += S[(x+2)%length,y-2]
    return total_adjacent_spins
    
def energy_change(S, coor):
    # Multiply spin site by all adjacent elements
    total_adjacent_spins = grid_adjacents(S, coor)
    return 2*S[coor[0],coor[1]]*(J*(total_adjacent_spins) + B*mu)

def make_ani(deltas):
    # Makes an animation from inital state, to delta values
    init_state = -1*np.ones(N)
    all_states = []
    all_states.append(init_state.copy())
    # Average the state over each 100 time steps
    for i in range(len(deltas)):
        if deltas[i] != ():
            init_state[deltas[i]] *= -1
        if (i+1) % 10000 == 0:
            all_states.append(init_state.copy())
    return all_states

def save_ani(all_states):
    ims = []
    fig = plt.figure()
    for i, img in enumerate(all_states):
        ims.append([plt.imshow(img, cmap=binary,animated=True)])
    im_ani = animation.ArtistAnimation(fig, ims,
                                       interval=35,
                                       repeat_delay=1000,
                                       blit=True)
    im_ani.save("StateChange.mp4")

def TwoDIsing(state0, num_steps, J, mu, B, kT):
    ES = energy(state0, J, mu, B)
    energy_values = [ES]
    # Contains a copy of the state configuration so we don't have to store
    # 2**N**2 * # of time step elements
    #state_configs = np.array([state0, state0])
    deltas = [] # A lighter way of keeping track of how the state changes.
    rands = np.random.randint(N[0], size=(num_steps,2))
    count = 1
    magsi = []
    evals = []
    for x, y in rands:
        if count/num_steps > 0.5:
            magsi.append(np.abs(np.sum(state0)))
        # Trial step: flip spin at one random site
        state0[x][y] *= -1
        ET = ES + energy_change(state0, (x,y))
        if np.exp((ES-ET)/(kT)) > np.random.random():
            #state_configs[-1] = test_state      # replace the state, or
            ES = ET
            deltas.append((x,y))
        else:
            # advance the previous state forward
            state0[x,y] *= -1
            deltas.append(())
        energy_values.append(ES)
        count += 1
        if count % 100000 == 0: print((count / num_steps)*100,
                                      " %..................")
    magsi = np.array(magsi)
    return (state0, energy_values, deltas, np.average(magsi),
            np.average(energy_values[int(num_steps/2):]))
# ==============================================================================
# Constant Definitions
# ==============================================================================
N           = (200, 200)                # number of spin sites                                   
num_steps   = 20*N[0]*N[1]              # number of iterations
B           = 0.05                      # magnetic field                               
mu          = .33                       # g mu (not needed if B=0)
J           = 1.                        # exchange energy                              
state0      = -1*np.ones(N)             # Initalize the model
ktvals      = np.arange(0.05, 10, 0.05) # Domain of kT to iterate over
mags        = []                        # Average Magnetization
avgevals    = []                        # Average Energy
cv_vals     = []                        # Specific Heat
deltas      = []                        # Record the changes with each time step
big_deltas  = []                        # Then concatenate (for animation purposes)
np.random.seed()                        # random seed 
# ==============================================================================
# Main Loop
# ==============================================================================
for i, kt in enumerate(ktvals):
    rvals = TwoDIsing(state0, num_steps, J, mu, B, kt)
    E_set = np.array(rvals[1])
    plt.figure()
    plt.plot(E_set)
    plt.xlabel("Time Step")
    plt.ylabel("Energy")
    plt.savefig("Energy"+str(i)+".pdf")
    plt.close()
    cv_vals.append(calc_cv(E_set[int(num_steps/2):], kt, N[0]*N[1]))
    deltas.append(rvals[2])
    mags.append(rvals[3])
    avgevals.append(rvals[4])
    print(100*(i+1)/len(ktvals), "%..........................................")
# -----------------------------------------------------------------------------
# Plotting and Animating
# -----------------------------------------------------------------------------
# Magnetization vs kT
plt.figure()
plt.plot(ktvals[1:], mags[1:])
plt.xlabel("$kT$")
plt.ylabel("Average Magnetization")
plt.title("Average Magnetization vs $kT$")
plt.savefig("mags.pdf")
plt.close()
# Average Energy vs kT
plt.figure()
plt.plot(ktvals[1:], avgevals[1:])
plt.xlabel("$kT$")
plt.ylabel("Average Energy")
plt.title("Average Energy vs $kT$")
plt.savefig("avgEnergy.pdf")
plt.close()
# Specific Heat vs kT
plt.figure()
plt.plot(ktvals. cv_vals)
plt.xlabel("$kT$")
plt.ylabel("CV")
plt.title("Specific Heat vs Temperature")
plt.savefig("cv.pdf")
plt.close()
# Animate the change from beginning to end
print("Animating", 72*'=')
for delta in deltas:
    big_deltas.extend(delta)
frames = make_ani(big_deltas)
save_ani(frames)
plt.close()
# ==============================================================================
# Induce a magnetic pulse, and watch how the simulation changes
# ==============================================================================
#state0[89:109,89:109] = -1
#init_state = state0.copy()
#rvals = TwoDIsing(state0, num_steps, J, mu, B, kt)
#E_set = np.array(rvals[1])
#plt.plot(E_set)
#plt.xlabel("Time Step")
#plt.ylabel("Energy")
#plt.savefig("Energy_Pulse.pdf")
#plt.close()
#cv_vals.append(calc_cv(E_set[int(num_steps/2):], kt, N[0]*N[1]))
#deltas.append(rvals[2])
#mags.append(rvals[3])
#avgevals.append(rvals[4])
#
#all_states = []
#all_states.append(init_state.copy())
## Average the state over each 100 time steps
#for i in range(len(deltas[-1])):
#    if deltas[-1][i] != ():
#        init_state[deltas[-1][i]] *= -1
#    if (i+1) % 100 == 0:
#        all_states.append(init_state.copy())
#
#ims = []
#fig = plt.figure()
#for img in all_states:
#    ims.append([plt.imshow(img, animated=True)])
#im_ani = animation.ArtistAnimation(fig, ims,
#                                   interval=35,
#                                   repeat_delay=1000,
#                                   blit=True)
#im_ani.save("PulseChange.mp4")
