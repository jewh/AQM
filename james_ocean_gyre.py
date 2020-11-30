import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# define the simulation constants
r = 1.0 
beta = 1.0
lambd = - r / beta
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
# define the array in which we store the psi(x,y) values
simulation = np.zeros((xsteps, ysteps), dtype=float)
# define a function which computes the coefficients on the LHS of the scheme
def left_coeff(lambd, xstep, ystep):
    return 1 - lambd / xstep
l1 = left_coeff(lambd, xsteps, ysteps)
# now solve for each value 
for x in range(1, xsteps-1):
    for y in range(2, ysteps-1):
        # calculate terms in the scheme
        frac1 = (-2 * simulation[x, y] + simulation[x-1, y]) / xstep**2
        frac2 = (simulation[x, y+1] - 2 * simulation[x, y] + simulation[x, y-1]) / ystep**2 
        # now calculate the new value
        new_val = simulation[x, y] + xstep * lambd * (frac1 + frac2)
        new_value = new_val / l1
        # now append 
        simulation[x+1, y] = new_value

print(simulation)
# now visualise
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
X = np.linspace(xmin, xmax, xsteps)
Y = np.linspace(ymin, ymax, ysteps)
ax.plot_surface(X=X, Y=Y, Z=simulation)
plt.savefig("test.png")






