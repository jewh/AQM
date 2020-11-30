import numpy as np

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
# define the array in which we store the psi(x,y) values
simulation = np.zeros((xsteps, ysteps), dtype=float)
# define a series of functions which compute the coefficients in the problem matrix
def first_coeff(lambd, xstep, ystep):
    return -lambd * xtsep**2
def second_coeff(lambd, xstep, ystep):
    return (xstep - lambd) * ystep**2 + 2 * lambd * xstep**2
def third_coeff(lambd, xstep, ystep):
    return -lambd * xtsep**2
# and likewise for the right hand side (rhs)
def first_rhs(lambd, xstep, ystep):
    return (2*lambd - xstep) * ystep**2
def second_rhs(lambd, xstep, ystep):
    return (2*lambd - xstep) * ystep**2






