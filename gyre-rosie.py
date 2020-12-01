import numpy as np 
import matplotlib.pyplot as plt   
import scipy.linalg as sp    
        
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
            self.dx = self.W/self.n
            self.dy = self.L/self.m
            
        def get_lhs(self): # without wind stress for now 
            
            A = self.beta/self.dx + (2*self.r)/(self.dx**2)
            B = self.beta + (4*self.r)/(self.dx**2) + (4*self.r)/(self.dy**2)
            C = (2*self.r)/(self.dy**2)
            D = (2*self.r)/(self.dx**2)
            
            matrix = np.zeros(((self.n+1)*(self.m+1),(self.n+1)*(self.m+1)))

            for x in range(len(matrix)): # run through rows of matrix
                i = x%self.m # y
                j = (x-i)/self.m # x
                
                if i==0 or i==self.m+1 or j==0 or j==self.n+1:
                    matrix[x][x] = 1
                    
                else:
                    matrix[x][x] = B
                    matrix[x][x-1] = C
                    matrix[x][x+1] = D
                    matrix[x][x-self.n] = D
                    matrix[x][x+self.n] = A
                    
            self.rhs = matrix
            
        def get_rhs(self):
            
            vector = np.zeros((self.n+1)*(self.m+1))
            
            for a in range(len(vector)):
                
                i = a%self.m # y co-ord
                y = self.dy*i # y-value
                
                j = (a-i)/self.m # x co-ord
                x = self.dx*j # x-value
                
                vector[a] = (1/(self.rho_0*self.H))*self.curl_tau(x, y, self.W, self.L, self.n, self.m, self.tau_0)
                
            self.lhs = vector
            
            
        def solve(self):
            
            self.get_lhs()
            self.get_rhs()
            
            sol = sp.solve(self.rhs,self.lhs)              
            sol_matrix = sol.reshape((self.m+1,self.n+1))
            
            x = np.arange(0, self.W+self.dx, self.dx)
            y = np.arange(0, self.L+self.dy, self.dy)
            
            fig,ax = plt.subplots()
            plt.contourf(x,y,sol_matrix)
            plt.xlabel('x')
            plt.ylabel('y')
            ax.set_aspect(1)
            cbar = plt.colorbar()
            cbar.set_label('$\psi$')
            plt.show()
 
def curl_tau(x,y,W,L,n,m,tau_0):
    return -((tau_0*2*np.pi*n)/L)*np.sin((2*np.pi*n*y)/L)

         
test_gyre = Gyre(2e-11, 1000, 1000, curl_tau, 1, 2e-8, 1000, 1000, 100, 100)
test_gyre.solve()

            
            
            
            
            
            
        
        
        
        
        
        

        
    
    
    
        