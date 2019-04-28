import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as scint
def ddp(phi, t, A, l, w_0, w, g):
    # ddotphi + g*dotphi + w0^2sin(phi) = (A/l)w0^2cos(wt)
    # Solve for y''
    # ddot phi = (A/l)w0^2cos(wt) - w0^2sin(phi) - g*dotphi
    # Define y_2 = phidot
    # Define y_1 = phi
    # Define y_1' = y_2 = phidot
    # Define y_2' = ddotphi
    # The equation is now y_2' = (A/l)w0^2cos(wt) - w0^2sin(y_1) - g*y_1'
    a = (A/l)*(w_0**2)*np.cos(w*t)
    b = (w_0**2) * np.sin(phi[0])
    c = g * phi[1]
    return (phi[1], a - b - c)
# No driving force (A = 0)
# No damping (g = 0)
t = np.linspace(0, 300, 30000)
y = scint.odeint(ddp, [10,-1*np.pi/2], t, args=(0, 1, 1, 0, .1))
# Plot for Phase Space
plt.plot(y[:,0], y[:,1])
plt.title("Phase Space of Damped Pendulum")
plt.xlabel("$\dot{\phi}$")
plt.ylabel("$\phi$")
plt.savefig("phasespace.pdf")
# Plot for phi vs time
plt.figure()
plt.plot(t, y[:,1])
plt.xlabel("Time")
plt.ylabel("$\phi$")
plt.title("$\phi$ vs Time")
plt.savefig("phivtime.pdf")

# Plot for phidot vs time
plt.figure()
plt.plot(t, y[:,0])
plt.xlabel("Time")
plt.ylabel("$\dot{\phi}$")
plt.title("$\dot{\phi}$ vs time")
plt.savefig("phidotvtime.pdf")

# Part b, adjust errything
y = scint.odeint(ddp, [-20, -1*np.pi/8], t, args=(1, 4, 1.5, 2, .1))
plt.figure()
plt.plot(t, y[:,1])
plt.xlabel("Time")
plt.ylabel("$\phi$")
plt.title("$\phi$ for a Driven Damped Pendulum")
plt.savefig("bphi.pdf")

plt.figure()
plt.plot(y[:,0], y[:,1])
plt.xlabel("$\dot{\phi}$")
plt.ylabel("$\phi$")
plt.savefig("bphasespace.pdf")

# Plot for phidot vs time
plt.figure()
plt.plot(t, y[:,0])
plt.xlabel("Time")
plt.title("$\dot{\phi}$ vs time")
plt.ylabel("$\dot{\phi}$")
plt.savefig("bphidotvtime.pdf")
