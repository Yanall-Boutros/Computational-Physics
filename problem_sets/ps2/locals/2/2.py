import numpy as np
import scipy.integrate as scint
import matplotlib.pyplot as plt
w = np.sqrt(3)*np.pi
u = 0.334
I = 1.567
def van_der_pol(Y,t, w, u, I):
    # y'' = -w^2*y - u(y^2 - I^2)y'
    # Define x_2 = y'
    # Define x_1 = y
    # Define x_1' = x_2
    # Define x_2' = y''
    # Our equation is now
    # x_2' = -w^2*x_1 - u(x_1^2 - I^2)*x_2
    # Solve for x_1' and x_2'
    return [Y[1], -w**2*Y[0]-u*(Y[0]**2 - I**2)*Y[1]]
t = np.linspace(-100, 100, 100000)
y = scint.odeint(van_der_pol, [0.3,5], t, (w, u, I))
plt.plot(y[:,0], y[:,1])
plt.savefig("vdp.png")
plt.figure()
plt.plot(t, y[:,0])
plt.savefig("y0.png")
plt.figure()
plt.plot(y, y[:,1])
plt.savefig("y1.png")
