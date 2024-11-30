import numpy as np
from Method import Method,Equations
TOL = 1e-7
class LU(Method): 
    def __init__(self, coff, sol, method, sig):
        super().__init__(coff, sol, sig)
        self.method = method

    def apply(self):
        self.coff = self.sign_array(self.coff)
        if self.method == "Doolittle":
            return self.dooLittle()
        elif self.method == "Crout":
            return self.crout()
        elif self.method == "Cholesky":
            return self.cholesky()
        else:
            raise ValueError("Invalid method")
        

    def dooLittle(self):
        if TOL <= 0:
            raise ValueError("Tolerance must be a positive number")
        
        n = len(self.coff)
        lower = np.zeros((n, n))
        upper = np.zeros((n, n))

        for i in range(n):
            # Upper Triangular
            for j in range(i, n):
                sum = 0
                for k in range(i):
                    sum += self.sign(lower[i, k] * upper[k, j])
                upper[i, j] = self.sign(self.coff[i, j] - sum)

            # Lower Triangular, i and j are swapped
            for j in range(i, n):
                if j == i:
                    lower[j, i] = 1
                else:
                    sum = 0
                    for k in range(i):
                        sum += self.sign(lower[j, k] * upper[k, i])
                    if abs(upper[i, i]) < TOL:
                        raise ZeroDivisionError("Singular Matrix")
                    lower[j, i] = self.sign((self.coff[j, i] - sum) / upper[i, i])
        
        # Solving using lower
        temp = self.coff 
        self.coff = lower
        self.sol = self.backSub()
        self.coff = temp
        return self.sol
    
    def crout(self):
        if TOL <= 0:
            raise ValueError("Tolerance must be a positive number")
        
        n = len(self.coff)
        lower = np.zeros((n, n))
        upper = np.eye(n)

        for j in range(n):
            for i in range(j, n):
                sum = 0
                for k in range(j):
                    sum += self.sign(lower[i, k] * upper[k, j])
                lower[i, j] = self.sign(self.coff[i, j] - sum)

            for i in range(j + 1, n):
                sum = 0
                for k in range(j):
                    sum += self.sign(lower[j, k] * upper[k, i])
                if abs(lower[j, j]) < TOL:
                    raise ZeroDivisionError("Singular Matrix")
                upper[j, i] = self.sign((self.coff[j, i] - sum) / lower[j, j])

        # Solving using lower
        temp = self.coff 
        self.coff = lower
        self.sol = self.backSub()
        self.coff = temp
        return self.sol
    
    def cholesky(self):
        if TOL <= 0:
            raise ValueError("Tolerance must be a positive number")
        n = len(self.coff)

        if not np.allclose(self.coff, self.coff.T, atol=TOL):
            raise ValueError("Input matrix must be symmetric")
        
        lower = np.zeros((n, n))

        for i in range(n):
            for j in range(i + 1):
                sum = 0
                if j == i:  # Diagonal elements
                    for k in range(j):
                        sum += self.sign(lower[j, k] ** 2)
                    diff = self.sign(self.coff[j, j] - sum)
                    if diff < 0:
                        raise ValueError("Matrix is not positive definite")
                    lower[j, j] = self.sign(np.sqrt(diff))
                else:
                    for k in range(j):
                        sum += self.sign(lower[i, k] * lower[j, k])
                    if abs(lower[j, j]) < TOL:
                        raise ZeroDivisionError("Singular Matrix")
                    lower[i, j] = self.sign((self.coff[i, j] - sum) / lower[j, j])

        # return lower, lower.T

        # Solving using lower
        temp = self.coff 
        self.coff = lower
        self.sol = self.backSub()
        self.coff = temp
        return self.sol

