import numpy as np
from Method import Method,Equations

class LU(Method): 
    def __init__(self,coff,sol):
        super().__init__(coff,sol)

    
    def dooLittle(self , a, tolerance=1e-7):
        n = len(a)
        lower = np.zeros((n, n))
        upper = np.zeros((n, n))

        for i in range(n):
            # Upper Triangular
            for j in range(i, n):
                sum = 0
                for k in range(i):
                    sum += lower[i, k] * upper[k, j]
                upper[i, j] = a[i, j] - sum

            # Lower Triangular, i and j are swapped
            for j in range(i, n):
                if j == i:
                    lower[j, i] = 1
                else:
                    sum = 0
                    for k in range(i):
                        sum += lower[j, k] * upper[k, i]
                    if abs(upper[i, i]) < tolerance:
                        raise ZeroDivisionError("Singular Matrix")
                    lower[j, i] = (a[j, i] - sum) / upper[i, i]
        return lower, upper
    
    def crout(self, a, tolerance=1e-7):
        n = len(a)
        lower = np.zeros((n, n))
        upper = np.eye(n)

        for j in range(n):
            for i in range(j, n):
                sum = 0
                for k in range(j):
                    sum += lower[i, k] * upper[k, j]
                lower[i, j] = a[i, j] - sum

            for i in range(j + 1, n):
                sum = 0
                for k in range(j):
                    sum += lower[j, k] * upper[k, i]
                if abs(lower[j, j]) < tolerance:
                    raise ZeroDivisionError("Singular Matrix")
                upper[j, i] = (a[j, i] - sum) / lower[j, j]

        return lower, upper
    
    def cholesky(self, a, tolerance=1e-7):
        n = len(a)
        lower = np.zeros((n, n))

        for i in range(n):
            for j in range(i + 1):
                sum = 0
                if j == i:  # Diagonal elements
                    for k in range(j):
                        sum += lower[j, k] ** 2
                    diff = a[j, j] - sum
                    if diff < 0:
                        raise ValueError("Matrix is not positive definite")
                    lower[j, j] = np.sqrt(diff)
                else:
                    for k in range(j):
                        sum += lower[i, k] * lower[j, k]
                    if abs(lower[j, j]) < tolerance:
                        raise ZeroDivisionError("Singular Matrix")
                    lower[i, j] = (a[i, j] - sum) / lower[j, j]
        return lower, lower.T

