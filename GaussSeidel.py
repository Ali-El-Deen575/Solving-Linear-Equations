import numpy as np
from Method import Method
import math
class GaussSeidel(Method):
    def __init__(self,coff,sol,guess,iter,tol,sig=5):
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
        if self.iter !=None and self.tol !=None :
            for z in range(self.iter):
                for i in range(self.n):
                    s1=sum(self.coff[i][j]*x[j] for j in range(i))
                    s2=sum(self.coff[i][j]*x[j] for j in range(i + 1, self.n))
                    s1 = self.sign(s1)
                    s2 = self.sign(s2)
                    x[i]=(self.sol[i]-s1-s2)/self.coff[i,i]
                    x[i] =self.sign(x[i])
                if np.linalg.norm(np.dot(self.coff,x)-self.sol,np.inf) < (self.tol*max(1.0,np.linalg.norm(x,np.inf))):
                    return x , z+1
            return x , self.iter
        elif self.iter !=None:
            for z in range(self.iter):
                for i in range(self.n):
                    s1=sum(self.coff[i][j]*x[j] for j in range(i))
                    s2=sum(self.coff[i][j]*x[j] for j in range(i + 1, self.n))
                    s1 = self.sign(s1)
                    s2 = self.sign(s2)
                    x[i]=(self.sol[i]-s1-s2)/self.coff[i,i]
                    x[i] =self.sign(x[i])
            return x ,self.iter
        else:
            iteration = 0
            while True:
                iteration += 1
                for i in range(self.n):
                    s1=sum(self.coff[i][j]*x[j] for j in range(i))
                    s2=sum(self.coff[i][j]*x[j] for j in range(i + 1, self.n))
                    s1 = self.sign(s1)
                    s2 = self.sign(s2)
                    x[i]=(self.sol[i]-s1-s2)/self.coff[i,i]
                    x[i] =self.sign(x[i])
                if np.linalg.norm(np.dot(self.coff,x)-self.sol,np.inf) < (self.tol*max(1.0,np.linalg.norm(x,np.inf))):
                    return x , iteration
                if iteration > 100:
                    break
            return x , iteration

