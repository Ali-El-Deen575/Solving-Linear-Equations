import numpy as np
from Method import Method
import math

class Jacobi(Method):
    def __init__(self,coff,sol,guess,iter,tol,sig):
        super().__init__(coff, sol, sig) 
        if guess is not None:
            self.guess=np.array(guess, dtype=float)
        else:
            self.guess=np.zeros(len(sol))
        self.iter=iter
        self.tol=tol
        self.n=len(sol)

    def apply(self):
        if self.guess is not None:
            x=np.array(self.guess,dtype=float)
        else:
            x=np.zeros(self.n)
        if (self.iter !=None and self.tol !=None):
            for z in range(self.iter):               
                y=np.zeros_like(x)
                for i in range(self.n):
                    s=sum(self.coff[i][j] * x[j] for j in range(self.n) if j != i)
                    s = self.sign(s)
                    y[i]=(self.sol[i]-s)/self.coff[i,i]
                    y[i] = self.sign(y[i])
                if np.linalg.norm(y-x,np.inf) < (self.tol*max(1.0,np.linalg.norm(y,np.inf))):
                    return y
                x = y
            return x
        elif self.iter !=None:
            for z in range(self.iter):
                y=np.zeros_like(x)
                for i in range(self.n):
                    s=sum(self.coff[i][j] * x[j] for j in range(self.n) if j != i)
                    s = self.sign(s)
                    y[i]=(self.sol[i]-s)/self.coff[i,i]
                    y[i] = self.sign(y[i])
                x = y        
            return x 
        else :
            y = np.zeros_like(x)
            while True:
                y = np.zeros_like(x)
                for i in range(self.n):
                    s = sum(self.coff[i][j] * x[j] for j in range(self.n) if j != i)
                    s = self.sign(s)
                    y[i] = (self.sol[i] - s) / self.coff[i, i]
                    y[i] = self.sign(y[i])
                if np.linalg.norm(y - x, np.inf) < (self.tol*max(1.0,np.linalg.norm(y,np.inf))):
                    break
                x = y
            return x
    def sign (self,value):
        if value == 0:
            return 0
        magnitude = math.floor(math.log10(abs(value)))
        scale =  10 ** (self.sig -1 - magnitude)
        rounded = round(value * scale) / scale
        return rounded     

    
       
sol = np.array([7, 5])
coff = np.array([[3, 2],[1, 3]])
guess = np.array([0,0])
jacobi =Jacobi(coff,sol,guess,iter=3,tol=None,sig=5)   
print(jacobi.apply()) 