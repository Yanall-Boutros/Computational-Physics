# ==============================================================================
# Import Statements
# ==============================================================================
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt
# ==============================================================================
# Function Definitions
# ==============================================================================
# Code from chi_squared testing grabbed from
# old lab code
def W(n, xdata, ydata, yerror):
    w_of_n = np.sum(((xdata**n)*ydata)/(yerror**2))
    return w_of_n

def U(n, xdata, yerror):
    u_of_n = np.sum((xdata**n)/(yerror**2))
    return u_of_n

def chi_squared(ydata, y_bestfit, sigma):
    cs = np.sum(((ydata - y_bestfit)**2)/(sigma**2))
    csr = cs / (len(ydata)-1)
    return (cs, csr)
    
def bretWig(e, a_1, a_2, a_3):
    # Credit to Jeff for the function signature
    return a_1 / ((e-a_2)**2 + a_3**2 / 4)

def LeastSquaresFit(xdata, ydata, y_sigma, func_pntr):
   # Least Squares Fit was originally a python program written by Prof.
   # David Smith, I have adapted and altered it to be generalized as an
   # individual function.
   # This is some old code I used in Phys-134
   if type(func_pntr) is not type(LeastSquaresFit):
       print("Function Pointer (func_pntr) not provided.")
       return
   xsmooth = np.linspace(np.min(xdata),np.max(xdata), 1000)
   popt, pcov = opt.curve_fit(func_pntr, xdata, ydata,
           sigma=y_sigma, absolute_sigma=1)
   fsmooth_next = func_pntr(xsmooth, *popt)
   plt.plot(xsmooth, fsmooth_next, color='green',
            label='Smoothed Line of Best Fit', alpha=0.5)
   plt.savefig("LSF.pdf")
   return xsmooth, fsmooth_next, popt, pcov
# ==============================================================================
# Data Initalization and Plotting
# ==============================================================================
E = np.array([0., 25., 50., 75., 100., 125., 150., 175., 200.])
y = np.array([10.6, 16.0, 45.0, 83.5, 52.8, 19.9, 10.8, 8.25, 4.7])
err=np.array([.934, 1.79, 4.15, 8.55, 5.15, 2.15, 1.08, 6.29, 4.14 ])

plt.errorbar(E, y, yerr=err, label="Raw data", linestyle="--", alpha=0.5)
plt.title("Fitting of bretWig function")
plt.xlabel("Energy (MeV)")
plt.ylabel("Cross-section (mb)")

xbest, ybest, popt, popc = LeastSquaresFit(E, y, err, bretWig)
ygood = bretWig(E, *popt)
plt.plot(E, ygood, label="Line of Best Fit (Shape = size of dataset)")
plt.legend(loc=1)
plt.savefig("LSF.pdf")
print("E_r = ", popt[0], ", f_r = ", popt[1], ", Gamma = ", popt[2])
print("Covariance Matrix: ")
print(popc)

cs, csr = chi_squared(y, ygood, err)
print("Chi squared = ", cs, "Reduced chi squared = ", csr)
