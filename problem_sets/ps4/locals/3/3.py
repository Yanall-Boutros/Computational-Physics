import numpy as np
import matplotlib.pyplot as plt
def theoretical_M(N,J, k, T, B):
    a = N*np.exp(J/(k*T))
    b = np.sinh(B/(k*T))
    num = a*b
    c = np.exp(2*J/(k*T))
    d = np.sinh(B/(k*T))**2
    e = np.exp(-2*J/(k*T))
    den = np.sqrt(c*d + e)
    return num/den
def M(mu, spin_states):
    return mu*np.sum(spin_states)
# Ising model implementation with Metropolis algorithm
# Based on IsingViz.py from Landau, et al.
# As well as S. Ritz starter code
N     = 2000                 # number of spin sites                                     
B     = 0.05               # magnetic field                               
mu    = .33                # g mu (not needed if B=0)
J     = 1.                 # exchange energy                              
k     = 1.                 # Boltzmann constant
T     = 1.                 # temperature                                 
state = np.ones(N)        # spin states: up(1), down (-1)
np.random.seed() 
        
def energy(S):                                  
    first_set = np.concatenate([np.array([S[-1]]), S[:-1]])
    FirstTerm = np.sum(-J*first_set[:-1]*first_set[1:])
    SecondTerm = np.sum(-mu*S*B)                                          
#    for i in range(-1,N-1):         # by starting with index -1, we can implement periodic boundary conditions
#        FirstTerm += -J * S[i]*S[i+1]
    return (FirstTerm + SecondTerm); 

state = -1*state
ES = energy(state)
num_steps = 40*N
energy_values = []
step_index_values = []
m = []
for istep in range(num_steps):
    test_state = state.copy()
    random_site = int(N*np.random.random())
    test_state[random_site] *= -1.   # Trial step: flip spin at one random site
    ET = energy(test_state)
    R = np.exp((ES-ET)/(k*T))           # Boltzmann test
    #print(ES, ET, k*T, R)
    if R > np.random.random():
        state = test_state
        ES = ET
    # else stay as is
    
    energy_values.append(ES)
    step_index_values.append(istep)
    m.append(M(mu, state))
    #row = ''
    #for site in state:
    #    if site==1.0:
    #        row += 'X'
    #    else:
    #        # fancy codes for red O's (thanks, J. Nielsen). You could also make a scatter plot with matplotlib.
    #        row += '\x1b[31mO\x1b[0m'
    #print(row)
plt.plot(step_index_values, energy_values)
plt.xlabel('time step')
plt.ylabel('energy')
plt.title('%d spins: kT = %.1f' % (N, k*T))
plt.savefig("spins.pdf")
print("Theoretical M value:")
print(theoretical_M(N,J, k, T, B))
plt.figure()
plt.plot(step_index_values, m)
plt.title("Magnetization vs Time for kT = 1")
plt.xlabel("time_step")
plt.ylabel("$M = \mu \sum_i^N \overline{s_i}$")
plt.savefig("m.pdf")
