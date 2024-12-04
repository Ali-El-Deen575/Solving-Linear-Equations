import numpy as np
from Method import Method,Equations

class GaussElemination(Method):
    def __init__(self,coff,sol,sig=5):
        super().__init__(coff,sol,sig)
    def apply(self):
        return self.backSub()


