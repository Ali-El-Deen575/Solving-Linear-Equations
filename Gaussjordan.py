import numpy as np
from Method import Method, Equations

import numpy as np

class GaussJordan(Method):
    def __init__(self, coff, sol, sig):
        super().__init__(coff, sol, sig)

    def apply(self):
        self.reducedEchelon()
        return self.sol  # Round final results for presentation

    def forwardElimination(self):
        n = len(self.coff)
        for i in range(n):
            self.pivoting(i, i)  # Ensure pivoting for numerical stability
            pivot = self.coff[i, i]
            if pivot == 0:
                return False

            # Normalize pivot row
            self.coff[i] = self.sign_array(self.coff[i] / pivot)  # Normalize and round
            self.sol[i] = self.sign(self.sol[i] / pivot)          # Round solution

            # Eliminate all rows below
            for j in range(i + 1, n):
                factor = self.coff[j, i]
                self.coff[j] = self.sign_array(self.coff[j] - factor * self.coff[i])  # Eliminate and round
                self.sol[j] = self.sign(self.sol[j] - factor * self.sol[i])           # Round solution
        return True        

    def reducedEchelon(self):
        n = len(self.coff)
        for i in range(n - 1, -1, -1):  # Start from the last row
            # Eliminate all rows above
            for j in range(i - 1, -1, -1):
                factor = self.coff[j, i]
                self.coff[j] = self.sign_array(self.coff[j] - factor * self.coff[i])
                self.sol[j] = self.sign(self.sol[j] - factor * self.sol[i])
    
               
sol = np.array([7, 12, 13])
sol = sol.astype(float)
coff = np.array([[3, 2, -1],[1, 3, 2],[2, -1, 4]])
coff = coff.astype(float)
jr =GaussJordan(coff,sol,7)   
print(jr.apply())                 