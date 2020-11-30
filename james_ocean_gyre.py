import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# define the simulation constants
r = 1.0 
beta = 1.0
lambd = - r / beta
n = 1.0
t0 = 1.0
rho = 1.0
H = 1.0
# and the boundary conditions 
xmin = 0
xmax = 10
ymin = 0
ymax = 1 # setting this to 1, as per the suggestion that r/beta*ymax = 1
# and now the simulation parameters
xsteps = 100
ysteps = xsteps
xstep = (xmax - xmin) / xsteps
ystep = (ymax - ymin) / ysteps
# define the surface stressors
def tx(x, y):
    return t0 * np.cos(2 * np.pi * y / ymax)
# don't define ty as it's 0 

# define the domain
X = np.linspace(xmin, xmax, xsteps)
Y = np.linspace(ymin, ymax, ysteps)
# define the array in which we store the psi(x,y) values
simulation = np.zeros((xsteps, ysteps), dtype=float)
# define a function which computes the coefficients on the LHS of the scheme
def left_coeff(lambd, xstep, ystep):
    return 1 - lambd / xstep
l1 = left_coeff(lambd, xsteps, ysteps)
# now solve for each value 
for x in range(1, 90):
    for y in range(2, ysteps-1):
        # calculate terms in the scheme
        frac1 = (-2 * simulation[x, y] + simulation[x-1, y]) / xstep**2
        frac2 = (simulation[x, y+1] - 2 * simulation[x, y] + simulation[x, y-1]) / ystep**2 
        stressors = (tx(X[x], Y[y]) - tx(X[x], Y[y-1])) / ystep
        stressors = - stressors * xstep / (rho * H * beta)
        # now calculate the new value
        new_val = simulation[x, y] + stressors + xstep * lambd * (frac1 + frac2)
        new_value = new_val / l1
        # now append 
        simulation[x+1, y] = new_value
        print("\nNew value is ", new_value, " at y=", y)
        if np.isnan(new_value):
            print("\nNaN found at x=", x, ", y=", y)
            break  

# now visualise
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(X=X, Y=Y, Z=simulation)
plt.savefig("surface.png")






