import numpy as np
from Method import Method
import math
class GaussSeidel(Method):
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
                    return x
            return x
        elif self.iter !=None:
            for z in range(self.iter):
                for i in range(self.n):
                    s1=sum(self.coff[i][j]*x[j] for j in range(i))
                    s2=sum(self.coff[i][j]*x[j] for j in range(i + 1, self.n))
                    s1 = self.sign(s1)
                    s2 = self.sign(s2)
                    x[i]=(self.sol[i]-s1-s2)/self.coff[i,i]
                    x[i] =self.sign(x[i])
            return x
        else:
            while True:
                for i in range(self.n):
                    s1=sum(self.coff[i][j]*x[j] for j in range(i))
                    s2=sum(self.coff[i][j]*x[j] for j in range(i + 1, self.n))
                    s1 = self.sign(s1)
                    s2 = self.sign(s2)
                    x[i]=(self.sol[i]-s1-s2)/self.coff[i,i]
                    x[i] =self.sign(x[i])
                if np.linalg.norm(np.dot(self.coff,x)-self.sol,np.inf) < (self.tol*max(1.0,np.linalg.norm(x,np.inf))):
                    break
            return x
    def sign (self,value):
        if value == 0:
            return 0
        magnitude = math.floor(math.log10(abs(value)))
        scale =  10 ** (self.sig -1 - magnitude)
        rounded = round(value * scale) / scale
        return rounded         
sol = np.array([7, 12, 13])
coff = np.array([[3, 2, -1],[1, 3, 2],[2, -1, 4]])
guess = np.array([0,0,0])
seidel =GaussSeidel(coff,sol,guess,iter=None,tol=5e-2,sig =3)   
print(seidel.apply()) 

        