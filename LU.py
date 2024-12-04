import numpy as np
from GaussElemination import GaussElemination
from Method import Method,Equations
TOL = 1e-7
class LU(Method):
    def __init__(self, coff, sol, method, sig=5):
        super().__init__(coff, sol, sig)
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
        if TOL <= 0:
            raise ValueError("Tolerance must be a positive number")
          
        n = len(self.coff)
        lower = np.zeros((n, n))
        upper = np.zeros((n, n))
        P = np.eye(n)

        for i in range(n):
            # Partial Pivoting
            max_index = np.argmax(abs(self.coff[i:, i])) + i
            if i != max_index:
                self.coff[[i, max_index]] = self.coff[[max_index, i]]
                P[[i, max_index]] = P[[max_index, i]]
                self.sol[[i, max_index]] = self.sol[[max_index, i]]

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
                        raise ValueError("Singular Matrix")
                    lower[j, i] = self.sign((self.coff[j, i] - sum) / upper[i, i])
        
        # Solving using lower

        outer = GaussElemination(lower,self.sol,self.sig)
        Y = outer.forwardSub()
        inner = GaussElemination(upper,Y,self.sig)
        sol = inner.backSub()
        return sol
    
    def crout(self):
        if TOL <= 0:
            raise ValueError("Tolerance must be a positive number")
        
        n = len(self.coff)
        lower = np.zeros((n, n))
        upper = np.eye(n)
        P = np.eye(n)

        for j in range(n):
            # Partial Pivoting
            max_index = np.argmax(abs(self.coff[j:, j])) + j
            if j != max_index:
                self.coff[[j, max_index]] = self.coff[[max_index, j]]
                P[[j, max_index]] = P[[max_index, j]]
                self.sol[[j, max_index]] = self.sol[[max_index, j]]

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
                    raise ValueError("Singular Matrix")
                upper[j, i] = self.sign((self.coff[j, i] - sum) / lower[j, j])

        # Solving using lower
        outer = GaussElemination(lower,self.sol,self.sig)
        Y = outer.forwardSub()
        inner = GaussElemination(upper,Y,self.sig)
        sol = inner.backSub()
        return sol
    
    def cholesky(self):
        if TOL <= 0:
            raise ValueError("Tolerance must be a positive number")
        n = len(self.coff)

        if not np.allclose(self.coff, self.coff.T, atol=TOL):
            raise ValueError("Input matrix must be symmetric")
        
        lower = np.zeros((n, n))
        P = np.eye(n)

        for i in range(n):
            # Partial Pivoting
            max_index = np.argmax(abs(self.coff[i:, i])) + i
            if i != max_index:
                self.coff[[i, max_index]] = self.coff[[max_index, i]]
                P[[i, max_index]] = P[[max_index, i]]
                self.sol[[i, max_index]] = self.sol[[max_index, i]]

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
                        raise ValueError("Singular Matrix")
                    lower[i, j] = self.sign((self.coff[i, j] - sum) / lower[j, j])

        # return lower, lower.T

        # Solving using lower
        outer = GaussElemination(lower,self.sol,self.sig)
        Y = outer.forwardSub()
        inner = GaussElemination(lower.T,Y,self.sig)
        sol = inner.backSub()
        return sol

if __name__ == "__main__":
    sol = np.array([7, 12, 13])
    #sol = sol.astype(float)
    coff = np.array([[6, 15, 55],[15, 55, 225],[55, 225, 979]])
    #coff = coff.astype(float)
    jr =LU(coff,sol,"Crout")
    print(jr.apply())
    print(np.linalg.solve(coff,sol))

