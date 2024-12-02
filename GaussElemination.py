import numpy as np
from Method import Method,Equations

class GaussElemination(Method):
    def __init__(self,coff,sol,sig=5):
        super().__init__(coff,sol,sig)
    def apply(self):
        return self.backSub()

if __name__ == "__main__":
    sol = np.array([7, 12, 13])
    sol = sol.astype(float)
    coff = np.array([[3, 2, -1],[1, 3, 2],[2, -1, 4]])
    coff = coff.astype(float)
    jr =GaussElemination(coff,sol,4)
    print(jr.apply())
