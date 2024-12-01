import numpy as np
from Method import Method,Equations

class GaussElemination(Method):
    def __init__(self,coff,sol,sig , step_by_step):
        super().__init__(coff,sol,sig, step_by_step)
    def apply(self):
        if self.step_by_step:
            print("**** Gauss Elemination start ****")
            print("a = ")
            print(self.coff)
            print("b = ")
            print(self.sol)
        x = self.backSub()
        if self.step_by_step:
            print("**** Gauss Elemination end ****")
        return x

sol = np.array([7, 12, 13])
sol = sol.astype(float)
coff = np.array([[3, 2, -1],[1, 3, 2],[2, -1, 4]])
coff = coff.astype(float)
jr =GaussElemination(coff,sol,4)   
print(jr.apply())                 
       
