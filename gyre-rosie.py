import numpy as np 
import matplotlib.pyplot as plt   
import scipy.linalg as sp
from mpl_toolkits.mplot3d import Axes3D    

def get_directions(sol_matrix, x, y):
    # this is to get the horizontal and vertical components for visualising with streamplot()
    # approximate dpsi/dx and dpsi/dy with forward euler
    u = np.zeros((len(x), len(y)), dtype=float)
    v = np.zeros((len(x), len(y)), dtype=float)
    for j in range(0, len(y)-1):
        for i in range(0, len(x)-1):
            v[i, j] = -(sol_matrix[i+1, j] - sol_matrix[i, j]) / (x[i+1]-x[i])
            u[i, j] = (sol_matrix[i, j+1] - sol_matrix[i, j]) / (y[i+1]-y[i])
    return u, v

class Gyre:
        
    def __init__(self, beta, rho_0, H, curl_tau, tau_0, r, W, L, n, m):
        self.beta = beta
        self.rho_0 = rho_0
        self.H = H
        self.curl_tau = curl_tau
        self.tau_0 = tau_0
        self.r = r
        self.W = W
        self.L = L
        self.n = n
        self.m = m
        self.rhs = []
        self.lhs = []
        self.dx = self.W/self.m
        self.dy = self.L/self.n
        self.sol = []
        self.sol_matrix = []
        
    def get_lhs(self): # without wind stress for now 
        
        A = self.beta/self.dx + (self.r)/(self.dx**2)
        B = -self.beta/self.dx - (2*self.r)/(self.dx**2) - (2*self.r)/(self.dy**2)
        C = self.r/(self.dx**2)
        D = self.r/(self.dy**2)
        
        matrix = np.zeros(((self.n+1)*(self.m+1),(self.n+1)*(self.m+1)))

        for x in range(len(matrix)): # run through rows of matrix
            i = x%(self.n+1) # y
            j = (x-i)/(self.n+1) # x
            
            if i==0 or i==self.n or j==0 or j==self.m:
                matrix[x][x] = 1
                
            else:
                matrix[x][x] = B
                matrix[x][x-1] = D
                matrix[x][x+1] = D
                matrix[x][x-(self.n+1)] = C
                matrix[x][x+(self.n+1)] = A
                
        self.lhs = matrix
        
    def get_lhs_periodic(self): # without wind stress for now 
        
        A = self.beta/self.dx + (self.r)/(self.dx**2)
        B = -self.beta/self.dx - (2*self.r)/(self.dx**2) - (2*self.r)/(self.dy**2)
        C = self.r/(self.dx**2)
        D = self.r/(self.dy**2)
        
        size = (self.n+1)*(self.m+1)
        matrix = np.zeros((size,size))

        for x in range(len(matrix)): # run through rows of matrix
            i = x%(self.n+1) # y
            j = (x-i)/(self.n+1) # x
            
            if i==0 or i==self.n:
                matrix[x][x] = 1
                
            else:
                matrix[x][x] = B
                matrix[x][(x-1)%size] = D
                matrix[x][(x+1)%size] = D
                matrix[x][(x-(self.n+1))%size] = C
                matrix[x][(x+(self.n+1))%size] = A
                
        self.lhs = matrix
        
    def get_rhs(self):
        
        vector = np.zeros((self.m+1)*(self.n+1))
        
        for a in range(len(vector)):
            
            i = a%(self.n+1) # y co-ord
            y = self.dy*i # y-value
            
            j = (a-i)/(self.n+1) # x co-ord
            x = self.dx*j # x-value
            
            if i!=0 and i!=self.n and j!=0 and j!=self.m:
                
                vector[a] = (1/(self.rho_0*self.H))*self.curl_tau(x, y, self.W, self.L, self.tau_0)
            
            
        self.rhs = vector
        
                   
    def solve(self, bc):
        
        if bc=='BC':
            self.get_lhs()
        elif bc=='periodic':
            self.get_lhs_periodic()
        
        self.get_rhs()
        
        sol = sp.solve(self.lhs,self.rhs)             
        sol_matrix = np.transpose(sol.reshape((self.n+1,self.m+1)))

        x = np.arange(0, self.W+self.dx, self.dx)/1000
        y = np.arange(0, self.L+self.dy, self.dy)/1000

        # adding in 3D surface plot
        '''X,Y = np.meshgrid(x, y) # don't know why it needs meshgrid axes but for some reason it does
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.plot_surface(X,Y,sol_matrix, rstride=1, cstride=1,
                cmap='viridis', edgecolor='none')
        ax.set_xlabel("x (km)")
        ax.set_ylabel("y (km)")
        ax.set_zlabel("$\psi (x,y)$")
        ax.set_title(r"Ocean Gyre for $\tau^{x} (x, y) = cos(2 \pi n y / L)$," + f"\n r = {self.r}")
        plt.show()'''


        # contour plot
        fig,ax = plt.subplots()
        plt.contour(x,y,sol_matrix, cmap='RdBu', levels=15)
        # Want to make a streamplot in the background
        # So calculate the direction vectors u, v
        u, v = get_directions(sol_matrix, x, y)
        # now streamplot, behind the contours
        ax.streamplot(x, y, u, v, color='0.8', density=2)
        ax.quiver(18,733,0,-1)
        ax.quiver(18,281,0,1)
        plt.xlabel('x (km)')
        plt.ylabel('y (km)')
        ax.set_aspect(1)
        cbar = plt.colorbar()
        cbar.set_label('$\psi$')
        #plt.show()
        
    def plot_curl(self):
        
        curl = np.zeros((self.n+1,self.m+1))
        
        for i in range(self.n+1):
            for j in range(self.m+1):
                    x = self.dx*j
                    y = self.dy*i
                    curl[i][j] = self.curl_tau(x, y, self.W, self.L, self.tau_0)
                    
        x = np.arange(0, self.W+self.dx, self.dx)/1000
        y = np.arange(0, self.L+self.dy, self.dy)/1000
        # adding in 3D surface plot
        X, Y = np.meshgrid(x, y)
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.plot_surface(X,Y,curl, rstride=1, cstride=1,
                cmap='viridis', edgecolor='none')

        
        # # and now plot contours
        # fig,ax = plt.subplots()
        # plt.contourf(x,y,curl)
        # plt.xlabel('x (km)')
        # plt.ylabel('y (km)')
        # ax.set_aspect(1)
        # cbar = plt.colorbar()
        # cbar.set_label('curl $\\tau$')
        # plt.show()
            
 
def curl_tau(x,y,W,L,tau_0):
    return ((tau_0*2*np.pi)/L)*np.sin((2*np.pi*y)/L) 
         

test_gyre = Gyre(beta=2e-11, rho_0=1000, H=1000, curl_tau=curl_tau, tau_0=1, r=2e-7, W=10**6, L=10**6, n=50, m=50)
test_gyre.solve('BC')
# plt.savefig("tau=default BC=periodic.png")
plt.show()





            
            
            
            
            
            
        
        
        
        
        
        

        
    
    
    
        