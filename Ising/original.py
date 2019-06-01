# As seen in lecture:
# Ising model implementation with Metropolis algorithm
# Based on IsingViz.py from Landau, et al.

N     = 2000                 # number of spin sites                                   
num_steps = 20*N           # number of iterations
B     = 0.05                # magnetic field                               
mu    = .33                # g mu (not needed if B=0)
J     = 1.                 # exchange energy                              
np.random.seed() 
        
def energy(S, J, mu, B):
    FirstTerm = 0.
    SecondTerm = 0.
    for i in range(-1,N-1):
        FirstTerm += -J * S[i]*S[i+1]
    for i in range(0,N):   
        SecondTerm += -mu*S[i]*B
    return (FirstTerm + SecondTerm);

def energy_change(S, i):
    """Determine the change in energy if spin `i` is flipped
    
    `(i+1) % len(S)` implements the periodic boundary condition.
    """
    return 2*J*(S[i-1]*S[i] + S[i]*S[(i+1) % len(S)]) + 2*B*mu*S[i]

def ising(state0, num_steps, J, mu, B, kT):
    ES = energy(state0, J, mu, B)
    energy_values = []
    energy_values.append(ES)
    state = np.zeros([num_steps,N])        # spin states: up(1), down (-1)
    state[0,:]=state0
    for istep in range(1,num_steps):
        test_state = list(state[istep-1,:])
        random_site = int(N*np.random.random())
        # Trial step: flip spin at one random site
        test_state[random_site] *= -1.
        """Our key modification is here: instead of calculating the energy
        afresh, we only calculate the change therein.
        """
        ET = ES + energy_change(state[istep-1], random_site)
        R = np.exp((ES-ET)/(kT))           # Boltzmann test
        if R > np.random.random():
            state[istep,:] = test_state      # replace the state, or
            ES = ET
        else:
            # advance the previous state forward
            state[istep,:]=state[istep-1,:]
        energy_values.append(ES)
    return state, energy_values
