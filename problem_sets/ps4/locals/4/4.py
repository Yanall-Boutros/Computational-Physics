import numpy as np
import matplotlib.pyplot as plt

def bisect(func, left, right, Nmax, eps):
    for i in range(Nmax):
        mid = (left+right)/2
        if func(mid)*func(right) > eps:
            right = mid
        else:
            left = mid
        if np.abs(func(mid)) < eps:
            print("Root Found at ",mid)
            return mid
        if i == Nmax - 1:
            print("Search Failed")
            return mid

def fw(z):
    return np.tan(z) - np.sqrt((8/z)**2 - 1)

z = np.linspace(-8, 8, 10000)
plt.plot(z, fw(z))
plt.ylim(-30, 30)
# Should return -5.463
a = bisect(fw, -5.5, -5.3, 1000000, 0.0000000001)
print("x = ", str(a), " f(x) = ", str(fw(a)))
# Should return 6.830
b = bisect(fw, 6.7, 6.9, 1000000, 0.0000000001)
print("x = ", b, " f(x) = ", fw(b))
c = bisect(fw, 1, 1.5, 1000000, 0.0000000001)
print("x = ", c, " f(x) = ", fw(c))
# Should return 1.395
plt.savefig("plt.pdf")
