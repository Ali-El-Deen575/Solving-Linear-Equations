import numpy as np
from GaussElemination import GaussElemination
from Method import Method,Equations
TOL = 1e-7

class LU(Method): 
    def __init__(self, coff, sol, method, sig, step_by_step):
        super().__init__(coff, sol, sig, step_by_step)

        self.method = method

    def apply(self):
        #self.coff = self.sign_array(self.coff)
        if self.method == "Doolittle":
            
            return self.dooLittle()
        elif self.method == "Crout":
            return self.crout()
        elif self.method == "Cholesky":
            return self.cholesky()
        else:
            raise ValueError("Invalid method")
        

    def dooLittle(self):
        if self.step_by_step:
            print("**** Doolittle start ****")
            print("a = ")
            print(self.coff)
            print("b = ")
            print(self.sol)


        if TOL <= 0:
            raise ValueError("Tolerance must be a positive number")
        
        n = len(self.coff)
        lower = np.zeros((n, n))
        upper = np.zeros((n, n))
        if self.step_by_step:
            print("Lower = ")
            print(lower)
            print("Upper = ")
            print(upper)

        for i in range(n):
            # Upper Triangular
            for j in range(i, n):
                sum = 0
                for k in range(i):
                    sum += self.sign(lower[i, k] * upper[k, j])
                upper[i, j] = self.sign(self.coff[i, j] - sum)
                if self.step_by_step:
                    print(f"upper[{i}][{j}] = {upper[i][j]}")
                    print(f"upper = ")
                    print(upper)

            # Lower Triangular, i and j are swapped
            for j in range(i, n):
                if j == i:
                    lower[j, i] = 1
                    if self.step_by_step:
                        print(f"lower[{j}][{i}] =   1  (diagonal)")
                        print(f"lower = ")
                        print(lower)
                else:
                    sum = 0
                    for k in range(i):
                        sum += self.sign(lower[j, k] * upper[k, i])
                    if abs(upper[i, i]) < TOL:
                        raise ValueError("Singular Matrix")
                    lower[j, i] = self.sign((self.coff[j, i] - sum) / upper[i, i])

                    if self.step_by_step:
                        print(f"lower[{j}][{i}] = {lower[j][i]}")
                        print(f"lower = ")
                        print(lower)
        
        # Solving using upper
        if self.step_by_step:
            print("Solving using back substitution on upper") 


        outer = GaussElemination(lower,self.sol,self.sig,self.step_by_step)
        Y = outer.forwardSub()
        inner = GaussElemination(upper,Y,self.sig,self.step_by_step)
        sol = inner.backSub()
        return sol
    
    def crout(self):
        if self.step_by_step:
            print("**** Crout start ****")
            print("a = ")
            print(self.coff)
            print("b = ")
            print(self.sol)

        if TOL <= 0:
            raise ValueError("Tolerance must be a positive number")
        
        n = len(self.coff)
        lower = np.zeros((n, n))
        upper = np.eye(n)

        if self.step_by_step:
            print("Lower = ")
            print(lower)
            print("Upper = ")
            print(upper)

        for j in range(n):
            for i in range(j, n):
                sum = 0
                for k in range(j):
                    sum += self.sign(lower[i, k] * upper[k, j])
                lower[i, j] = self.sign(self.coff[i, j] - sum)

                if self.step_by_step:
                    print(f"lower[{i}][{j}] = {lower[i][j]}")
                    print(f"lower = ")
                    print(lower)

            for i in range(j + 1, n):
                sum = 0
                for k in range(j):
                    sum += self.sign(lower[j, k] * upper[k, i])
                
                if self.step_by_step:
                    print(f"product = {sum}")
                if abs(lower[j, j]) < TOL:


                if self.step_by_step:
                        print(f"product is almost 0 which means it's a singular matrix")
                    raise ZeroDivisionError("Singular Matrix")
                upper[j, i] = self.sign((self.coff[j, i] - sum) / lower[j, j])

                if self.step_by_step:
                    print(f"upper[{j}][{i}] = {upper[j][i]}")
                    print(f"upper = ")
                    print(upper)

        if self.step_by_step:
            print("Solving using back substitution on upper")
            
        outer = GaussElemination(lower,self.sol,self.sig)
        Y = outer.forwardSub()
        inner = GaussElemination(upper,Y,self.sig)
        sol = inner.backSub()
        return sol
    
    def cholesky(self):
        if self.step_by_step:
            print("**** Cholesky start ****")
            print("a = ")
            print(self.coff)
            print("b = ")
            print(self.sol)

        if TOL <= 0:
            raise ValueError("Tolerance must be a positive number")
        n = len(self.coff)

        if not np.allclose(self.coff, self.coff.T, atol=TOL):

            if self.step_by_step:
                print("Matrix is not symmetric")

            raise ValueError("Input matrix must be symmetric")
        
        lower = np.zeros((n, n))

        if self.step_by_step:
            print("Lower = ")
            print(lower)

        for i in range(n):
            for j in range(i + 1):
                sum = 0
                if j == i:  # Diagonal elements
                    for k in range(j):
                        sum += self.sign(lower[j, k] ** 2)
                    diff = self.sign(self.coff[j, j] - sum)

                    if self.step_by_step:
                        print(f"for [{i}][{j}] , diff = {diff} (diagonal)")

                    if diff < 0:
                        if self.step_by_step:
                            print("As diff < 0 , Matrix is not positive definite")

                        raise ValueError("Matrix is not positive definite")
                    lower[j, j] = self.sign(np.sqrt(diff))
                    
                    if self.step_by_step:
                        print(f"lower[{j}][{j}] = sqrt(diff) = {lower[j][j]}")
                        print(f"lower = ")
                        print(lower)
                else:
                    for k in range(j):
                        sum += self.sign(lower[i, k] * lower[j, k])
                    if abs(lower[j, j]) < TOL:

                        if self.step_by_step:
                            print(f"lower[{j}][{j}] is almost 0 which means it's a singular matrix")
                        raise ZeroDivisionError("Singular Matrix")

                    lower[i, j] = self.sign((self.coff[i, j] - sum) / lower[j, j])

                    if self.step_by_step:
                        print(f"lower[{i}][{j}] = {lower[i][j]}")
                        print(f"lower = ")
                        print(lower)
        
        if self.step_by_step:
            print("Upper = lower.T = ")
            print(lower.T)
            print("Solving using back substitution on upper")



        # Solving using upper
        outer = GaussElemination(lower,self.sol,self.sig)
        Y = outer.forwardSub()
        inner = GaussElemination(lower.T,Y,self.sig)
        sol = inner.backSub()
        return sol
sol = np.array([7, 12, 13])
#sol = sol.astype(float)
coff = np.array([[6, 15, 55],[15, 55, 225],[55, 225, 979]])
#coff = coff.astype(float)
jr =LU(coff,sol,"Crout",5)
print(jr.apply())
print(np.linalg.solve(coff,sol))

