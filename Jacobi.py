import numpy as np
from Method import Method
import math

class Jacobi(Method):
    def __init__(self,coff,sol,guess,iter,tol,sig=5):
        super().__init__(coff, sol, sig)
        if guess is not None:
            self.guess=np.array(guess, dtype=float)
        else:
            self.guess=np.zeros(len(sol))
        self.iter=iter
        self.tol=tol
        self.n=len(sol)
        for i in range(self.n):
            row_sum = sum(abs(self.coff[i][j]) for j in range(self.n) if j != i)
            if abs(self.coff[i][i]) < row_sum:
                raise ValueError(f"Matrix is not diagonally dominant at row {i}")

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
                    return y , z+1
                x = y
            return x , self.iter
        elif self.iter !=None:
            for z in range(self.iter):
                y=np.zeros_like(x)
                for i in range(self.n):
                    s=sum(self.coff[i][j] * x[j] for j in range(self.n) if j != i)
                    s = self.sign(s)
                    y[i]=(self.sol[i]-s)/self.coff[i,i]
                    y[i] = self.sign(y[i])
                x = y
            return x , self.iter
        else :
            y = np.zeros_like(x)
            iteration = 0
            while True:
                iteration += 1
                y = np.zeros_like(x)
                for i in range(self.n):
                    s = sum(self.coff[i][j] * x[j] for j in range(self.n) if j != i)
                    s = self.sign(s)
                    y[i] = (self.sol[i] - s) / self.coff[i, i]
                    y[i] = self.sign(y[i])
                if np.linalg.norm(y - x, np.inf) < (self.tol*max(1.0,np.linalg.norm(y,np.inf))):
                    break
                x = y
            return x , iteration

if __name__ == "__main__":
    sol = np.array([7, 5])
    coff = np.array([[3, 2],[1, 3]])
    guess = np.array([0,0])
    jacobi =Jacobi(coff,sol,guess,iter=3,tol=None,sig=5)
    print(jacobi.apply())