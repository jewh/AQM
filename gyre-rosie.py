import numpy as np 
import matplotlib.pyplot as plt   
import scipy.linalg as sp
from mpl_toolkits.mplot3d import Axes3D    
        
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
        
    def get_lhs(self): # without wind stress for now 
        
        A = self.beta/self.dx - (self.r)/(self.dx**2)
        B = -self.beta/self.dx + (2*self.r)/(self.dx**2) + (2*self.r)/(self.dy**2)
        C = -self.r/(self.dx**2)
        D = -self.r/(self.dy**2)
        
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
                matrix[x][x-(self.m+1)] = C
                matrix[x][x+(self.m+1)] = A
                
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
        
                   
    def solve(self):
        
        self.get_lhs()
        self.get_rhs()
        
        sol = sp.solve(self.lhs,self.rhs)              
        sol_matrix = np.transpose(sol.reshape((self.n+1,self.m+1)))
        
        x = np.arange(0, self.W+self.dx, self.dx)/1000
        y = np.arange(0, self.L+self.dy, self.dy)/1000
        # adding in 3D surface plot
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.plot_surface(X=x, Y=y, Z=sol_matrix)
        plt.savefig("surface.png")
        plt.show()

        # contour plot
        fig,ax = plt.subplots()
        plt.contourf(x,y,sol_matrix)
        plt.xlabel('x (km)')
        plt.ylabel('y (km)')
        ax.set_aspect(1)
        cbar = plt.colorbar()
        cbar.set_label('$\psi$')
        plt.show()
        plt.savefig("contour.png") # for the WSL users
        
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
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.plot_surface(x,y,curl)
        plt.show()
        
        # and now plot contours
        fig,ax = plt.subplots()
        plt.contourf(x,y,curl)
        plt.xlabel('x (km)')
        plt.ylabel('y (km)')
        ax.set_aspect(1)
        cbar = plt.colorbar()
        cbar.set_label('curl $\\tau$')
        plt.show()
            
 
def curl_tau(x,y,W,L,tau_0):
    return ((tau_0*2*np.pi)/L)*np.sin((2*np.pi*y)/L)
         
test_gyre = Gyre(2e-11, 1000, 1000, curl_tau, 1, 2e-5, 10**6, 10**6, 50, 50)
test_gyre.plot_curl()
test_gyre.solve()


            
            
            
            
            
            
        
        
        
        
        
        

        
    
    
    
        