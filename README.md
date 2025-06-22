# cubicspline

this is a code for natural cubic sline in numerical methods that finds a spline that connects all the given (x,y) points

 import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline

def natural_cubic_spline(x, y):
    n = len(x) - 1
    h = [x[i+1] - x[i] for i in range(n)]
    alpha = [0] * (n+1)
    for i in range(1, n):
        alpha[i] = (3/h[i]) * (y[i+1] - y[i]) - (3/h[i-1]) * (y[i] - y[i-1])
    l = [1] + [0]*n
    mu = [0]*(n+1)
    z = [0]*(n+1)
    for i in range(1, n):
        l[i] = 2*(x[i+1] - x[i-1]) - h[i-1]*mu[i-1]
        mu[i] = h[i]/l[i]
        z[i] = (alpha[i] - h[i-1]*z[i-1])/l[i]
    l[n] = 1
    z[n] = 0
    c = [0]*(n+1)
    b = [0]*n
    d = [0]*n
    a = y[:]
    for j in reversed(range(n)):
        c[j] = z[j] - mu[j]*c[j+1]
        b[j] = (a[j+1] - a[j])/h[j] - h[j]*(c[j+1] + 2*c[j])/3
        d[j] = (c[j+1] - c[j])/(3*h[j])
    return a, b, c, d, x

def spline_eval(x_val, a, b, c, d, x):
    for i in range(len(x) - 1):
        if x[i] <= x_val <= x[i+1]:
            dx = x_val - x[i]
            return a[i] + b[i]*dx + c[i]*dx**2 + d[i]*dx**3
    return None

# թեստային օրինակ
x = [0,1,2,3,4]
y = [1,2,0,2,1]
def f(x):
    return 1 / (1 + 25 * x**2)
n = 10 
x = np.linspace(-1, 1, n+1)
y = f(x)
xx = np.linspace(x[0], x[-1], 500)
cs_natural = CubicSpline(x, y, bc_type='natural')
#սպլայների հաշվարկում
a_nat, b_nat, c_nat, d_nat, _ = natural_cubic_spline(x, y)

# գնահատում
yy_nat = [spline_eval(v, a_nat, b_nat, c_nat, d_nat, x) for v in xx]
cs = CubicSpline(x, y, bc_type='natural')  
yy_scipy_natural = cs_natural(xx)

# պատկերում
plt.figure(figsize=(10, 8))
plt.subplot(1, 1, 1)
plt.plot(x, y, 'o', label='Given Points')
plt.plot(xx, yy_nat, label='Natural Spline')
plt.title("natural Cubic Splines (No SciPy)")
plt.plot(xx, yy_scipy_natural, '--', label='SciPy Natural Spline')
plt.title("SciPy CubicSpline Built-in")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()
