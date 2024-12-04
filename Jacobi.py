import numpy as np
from Method import Method
import math

class Jacobi(Method):
    def __init__(self,coff,sol,guess,iter,tol,sig=5, step_by_step=False):
        super().__init__(coff, sol, sig, step_by_step)
        if guess is not None:
            self.guess=np.array(guess, dtype=float)
        else:
            self.guess=np.zeros(len(sol))
        self.iter=iter
        self.tol=tol
        self.n=len(sol)
        for i in range(self.n):
            row_sum = sum(abs(self.coff[i][j]) for j in range(self.n) if j != i)
            # if abs(self.coff[i][i]) < row_sum:
                # raise ValueError(f"Matrix is not diagonally dominant at row {i}")

    def apply(self):
        if self.step_by_step:
            print("**** Jacobi start ****")
            print("a = ")
            print(self.coff)
            print("b = ")
            print(self.sol)

        if self.guess is not None:
            x=np.array(self.guess,dtype=float)
        else:
            x=np.zeros(self.n)
        
        if self.step_by_step:
            print("Initial guess = ")
            print(x)
        
        if (self.iter !=None and self.tol !=None):
            for z in range(self.iter):
                if self.step_by_step:
                    print(f"* Iteration {z+1} *")

                y=np.zeros_like(x)
                for i in range(self.n):
                    s=sum(self.coff[i][j] * x[j] for j in range(self.n) if j != i)
                    s = self.sign(s)
                    y[i]=(self.sol[i]-s)/self.coff[i,i]
                    y[i] = self.sign(y[i])

                    # if self.step_by_step:
                        # print(f"y{i+1} = {y[i]}")
                if self.step_by_step:
                    print(f"y = {y}")

                if np.linalg.norm(y-x,np.inf) < (self.tol*max(1.0,np.linalg.norm(y,np.inf))):
                    if self.step_by_step:
                        print(f"Solution found at iteration {z+1}")
                        print("**** Jacobi end ****")

                    return y , z+1
                
                if self.step_by_step:
                    print(f"didn't converge enough, we have to carry on.")
                
                x = y
            
            if self.step_by_step:
                print(f"Solution not found after {self.iter} iterations")
                print("**** Jacobi end ****")

            return x , self.iter
        elif self.iter !=None:
            for z in range(self.iter):
                if self.step_by_step:
                    print(f"* Iteration {z+1} *")

                y=np.zeros_like(x)
                for i in range(self.n):
                    s=sum(self.coff[i][j] * x[j] for j in range(self.n) if j != i)
                    s = self.sign(s)
                    y[i]=(self.sol[i]-s)/self.coff[i,i]
                    y[i] = self.sign(y[i])

                    if self.step_by_step:
                        print(f"y{i+1} = {y[i]}")
                if self.step_by_step:
                    print(f"y = {y}")
                    
                x = y        
            
            if self.step_by_step:
                print("**** Jacobi end ****")

            return x , self.iter
        
        else :
            y = np.zeros_like(x)
            iteration = 0
            while True:
                if self.step_by_step:
                    print(f"* Iteration {iteration+1} *")
                    
                iteration += 1
                y = np.zeros_like(x)
                for i in range(self.n):
                    s = sum(self.coff[i][j] * x[j] for j in range(self.n) if j != i)
                    s = self.sign(s)
                    y[i] = (self.sol[i] - s) / self.coff[i, i]
                    y[i] = self.sign(y[i])

                    if self.step_by_step:
                        print(f"y{i+1} = {y[i]}")
                if self.step_by_step:
                    print(f"y = {y}")

                if np.linalg.norm(y - x, np.inf) < (self.tol*max(1.0,np.linalg.norm(y,np.inf))):
                    if self.step_by_step:
                        print(f"Solution found at iteration {iteration}")
                        print("**** Jacobi end ****")

                    break

                if self.step_by_step:
                    print(f"didn't converge enough, we have to carry on.")

                x = y

            return x , iteration

if __name__ == "__main__":
    sol = np.array([7, 5])
    coff = np.array([[3, 2],[1, 3]])
    guess = np.array([0,0])
    jacobi =Jacobi(coff,sol,guess,iter=3,tol=None,sig=5)
    print(jacobi.apply())