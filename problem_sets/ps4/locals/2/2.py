import numpy as np
import matplotlib.pyplot as plt
AllCollisions = np.arange(0, 100000)
FREEDOMHAHAHAHA = []
FREEDOMESCAPEPOINT = []
FREEDOMENTRANCEPOINT = []
for neutron in AllCollisions:
    Nsteps = 15
    x=0.
    y=0.
    xvalues = [x]
    yvalues = [y]
    sum_of_R = 0
    theta = 2.*np.pi*np.random.random(Nsteps)
    for i in range(1, Nsteps):
        x += np.cos(theta[i])
        y += np.sin(theta[i])
        xvalues.append(x)
        yvalues.append(y)
        if x > 5 or x < 0:
            if x > 5:
                FREEDOMHAHAHAHA.append(neutron)
                FREEDOMESCAPEPOINT.append(y)
            break
print(str(100*len(FREEDOMHAHAHAHA) / len(AllCollisions))+"% Escape",
      "through to the other side")
plt.hist(FREEDOMESCAPEPOINT, bins=20)
plt.title("Counts for y-position escape values")
plt.xlabel("y-position")
plt.ylabel("Counts")
plt.savefig("hist.pdf")
