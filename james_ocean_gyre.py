import numpy as np
from scipy.linalg import solve
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# define the simulation constants
r = 1.0 
beta = 2.0
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
xsteps = 1000
ysteps = xsteps
xstep = (xmax - xmin) / xsteps
ystep = (ymax - ymin) / ysteps
# define the surface stressors
def dtxdy(y):
    return t0 * 2 * xstep * n * np.pi * np.sin(2 * np.pi * n * y / ymax) / (rho * H * ymax * beta)
# don't define ty as it's 0 

# define the domain
X = np.linspace(xmin, xmax, xsteps)
Y = np.linspace(ymin, ymax, ysteps)
# define the array in which we store the psi(x,y) values
simulation = np.zeros((xsteps, ysteps), dtype=float)
# Compute the coefficients on the LHS of the scheme
l1 = r * xstep / (beta * ystep**2) 
l2 = 1 + r / (beta * xstep) - 2 * r * xstep / (beta * ystep**2)
# now the RHS
r1 = 1 + 2 * r / (beta * xstep)
r2 = - r / (beta * xstep)
# now define a function that computes the problem matrix 
def problem(ysteps):
    matrix = np.zeros((ysteps, ysteps), dtype=float)
    y_position = 0
    for i in range(0, ysteps):
        counter = 0
        # put in the initial values in the solution vector b, so set these to 1
        if i in [0, 1, ysteps-2, ysteps-1]:
            matrix[i,i] = 1
        else:
            for j in range(y_position, y_position+3):
                if counter % 2 == 0: # if it's 0 or 2
                    matrix[i, j] = l1
                else:
                    matrix[i, j] = l2
                counter += 1
            y_position += 1 # now shift along one
    return matrix
# and now a function that creates a solution (b) vector 
def solution_vect(ysteps, psiy1, psiyn_1, previous_simulation_slices):
    # previous_simulation_slices is a 2 x ysteps array which contains the psi^y_{x} and psi^y_{x-1} values
    # make sure previous_simulation_slices[1] = x, and previous_simulation_slices[0] = x-1
    sol_vect = np.zeros(ysteps, dtype=float)
    sol_vect[1] = psiy1
    sol_vect[ysteps-2] = psiyn_1 # set the initial values 
    # Notice we don't set the initial values for y[0], y[ysteps-1] as end values are 0 by the boundary value
    for i in range(2, ysteps-2): 
        sol_vect[i] = dtxdy(Y[i]) + r1 * previous_simulation_slices[1, i] + r2 * previous_simulation_slices[0, i]
    return sol_vect
        

# now solve for each value 
for x in range(1, xsteps-1):
    # first need to calculate the values of psi at the edges
    # just do so with a forward euler scheme where we assume the laplacian of psi is 0 at the boundary
    if x == 1:
        for y in range(1, ysteps-1):
            simulation[x, y] = simulation[x-1, y] + xstep * dtxdy(Y[y])
    else:
        # want to generate a matrix problem for each x step in the solution
        # this matrix problem will solve for the values at x+1
        A = problem(ysteps)
        # get the edge values
        psi_1 = simulation[x-1, 1] + xstep * dtxdy(Y[1])
        psi_y_1 = simulation[x-1, ysteps-1] + xstep * dtxdy(Y[ysteps-1])
        b = solution_vect(ysteps, psi_1, psi_y_1, simulation[x-1:x+1])
        # now solve the problem 
        solution = solve(A, b)
        simulation[x+1] = solution

# now visualise
plt.contourf(X, Y, simulation)
plt.savefig("contour.png")
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(X=X, Y=Y, Z=simulation)
plt.savefig("surface.png")







