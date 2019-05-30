import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate
import matplotlib.animation as animation
############## Schrod and V definition. Psi is vectorized as (psi,psi') ##########
def f(y, x, k2fac,E):
    return [y[1],-k2fac*(E-V(x))*y[0]]

def V(x):
    if (abs(x) < 10.):  return (-16.0) # Depth of well
    else:               return (0.)
########## boundary condition match error ###############        
def diff(E, k2fac, xs, imatch, icount, Lwf, Rwf):
    # start on the left.
    Lwf[icount,0,0] = 1.e-5                                      
    Lwf[icount,0,1]= Lwf[icount,0,0]*np.sqrt(k2fac*abs(E))
    state = ([Lwf[icount,0,0], Lwf[icount,0,1]])
    Lwf[icount,:,:] = integrate.odeint(f, state, xs[0:len(Lwf[0,:,0])], args=(k2fac,E))
    # and now the right. odeint doesn't integrate backwards
    # apparently so flip by hand!:
    yb=np.zeros((2,len(Rwf[0,:,0])))
    yb[0,0] = 1.e-5                                      
    yb[0,1]=  yb[0,0]*np.sqrt(k2fac*abs(E)) # no minus sign for the slope,
                                            # as it's already going backward
    state = ([yb[0,0],yb[0,1]])
    xb=-xs[:(len(Lwf[0,:,0])-2):-1]
    yb = integrate.odeint(f, state,xb, args=(k2fac,E)) # step backward to imatch
    Rwf[icount,:,:]=yb[::-1]
    # get the relative normalization right. In effect, the search is
    # really just matching curvatures.
    norma=Lwf[icount,-1,0]/Rwf[icount,0,0]
    Rwf[icount,:,:]=norma*Rwf[icount,:,:]
    # calculate logarithmic derivative
    left = Lwf[icount,-1,1]/Lwf[icount,-1,0]        
    right = Rwf[icount,0,1]/Rwf[icount,0,0]                     
    ldif=(left - right)/(left + right)
    # print("left and right:",left,right)
    return ldif, Lwf[icount,:,:], Rwf[icount,:,:]

### parameters  ############
count_max = 100  # maximun number of search iterations
eps       = 1e-6 # Precision
n_steps   = 501  # number of (spatial) integration steps
xmax=20.         # set the size of the full +/- x range 
xmatch = 9.      # the x value at which the left and right sides should match
# guess search range
Eguesses = [(-15.5, -15.7), (-14.3, -14.4), (-12.0, -13.3)]
k2fac=0.04829
for i, minmax in enumerate(Eguesses):
    Emin = minmax[0]
    Emax = minmax[1]
    ############
    h=2.*xmax/(n_steps-1)
    xs=np.arange(-xmax,xmax,h)
    imatch=np.nonzero((xs>xmatch))[0][0]
    # arrays for the left and right wf(index 0) and wf'(index 1)
    Lwf = np.zeros((count_max,imatch+1,2))
    # we need these to overlap at one point, so add 1 more Lwf. 
    Rwf = np.zeros((count_max,n_steps-imatch-1,2))
    icount=0
    Diffmax, Lwf[icount,:,:], Rwf[icount,:,:] = diff(Emax, k2fac, xs, imatch,
                                                     icount, Lwf, Rwf)
    Diffmin, Lwf[icount,:,:], Rwf[icount,:,:] = diff(Emin, k2fac, xs, imatch,
                                                     icount, Lwf, Rwf)
    if ((Diffmax*Diffmin)>0) : 
        print("\n\n\n WHOOOOPS looks as if there is not one zero in this range\n\n\n",
                Diffmax, Diffmin)
        sys.exit()
    print("There is one zero between %.2f and %.2f. Let's proceed.\n"%(Emax,Emin))
    print("Diffmax and Diffmin:",Diffmax,Diffmin)
    for icount in range(0,count_max+1):
        # Iteration loop
        E = (Emax + Emin)/2.                                 # Divide E range
        Diffmax, Lwf[icount,:,:], Rwf[icount,:,:] = diff(Emax, k2fac, xs, imatch,
                                                         icount, Lwf, Rwf)
        Diff, Lwf[icount,:,:], Rwf[icount,:,:] = diff(E, k2fac, xs, imatch, icount,
                                                      Lwf, Rwf) 
        if (Diffmax*Diff > 0):  Emax = E  # Bisection algorithm
        else:                   Emin = E     
        if (abs((Emax-Emin)/E)  <  eps ): break
        print("Count: %d   Energy: %.5f    Diff %.2f "%(icount, E, Diff))
    
    print("\n Final eigenvalue E = %.5f "%E)
    print("iterations = ",icount)
    
    # let's plot the final wave functions
    fig = plt.figure()
    fig.suptitle('final left and right functions', fontsize=16)
    plt.scatter(xs[0:imatch+1],(Lwf[icount-1,:,0]))
    plt.scatter(xs[imatch:],(Rwf[icount-1,:,0]))
    plt.savefig("plt"+str(i)+".png")
    plt.show()
