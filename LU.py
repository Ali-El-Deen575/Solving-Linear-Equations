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
          
        n = len(self.coff)
        o = np.zeros(n)
        for i in range (n):
            o[i] =i
        P    = np.identity(n)
        L    = np.identity(n)
        U    = self.coff.copy()
        PF   = np.identity(n)
        LF   = np.zeros((n,n))
        for k in range(0, n - 1):
            index = np.argmax(abs(U[k:, k]))
            index = index + k 
            if index != k:
                P = np.identity(n)
                o[index], o[k] = o[k], o[index]
                P[[index, k], k:n] = P[[k, index], k:n]
                U[[index, k], k:n] = U[[k, index], k:n] 
                PF = np.dot(P, PF)
                LF = np.dot(P, LF)
            L = np.identity(n)
            for j in range(k+1,n):
                L[j, k]  = -(U[j, k] / U[k, k])
                L[j, k]  = self.sign(L[j, k])
                LF[j, k] =  (U[j, k] / U[k, k])
                LF[j, k] =  self.sign(LF[j, k])
            U = np.dot(L,U)
        np.fill_diagonal(LF, 1)
        
        # Solving using lower

        order = np.array(o, dtype=int)
        sort = self.sol[order]
        outer = GaussElemination(LF,sort,self.sig)
        Y = outer.forwardSub()
        inner = GaussElemination(U,Y,self.sig)
        sol = inner.backSub()
        return sol
    
    def crout(self):
        
        n = len(self.coff)
        lower = np.zeros((n, n))
        upper = np.eye(n)
        P = np.eye(n)

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
                    raise ValueError("Singular Matrix")
                upper[j, i] = self.sign((self.coff[j, i] - sum) / lower[j, j])

        # Solving using lower
        print(upper)
        print(lower)
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



