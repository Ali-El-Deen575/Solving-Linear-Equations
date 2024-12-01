import numpy as np
import math
class Equations():

    def __init__(self,num):
        self.i=0
        self.j=0
        self.sig= 5
        self.coff = np.zeros((num,num))
        self.sol = np.zeros(num)
        self.num =num
    def setCoff(self,b):
        if(self.j == self.num):
            self.sol[self.i] = b
            self.i +=1
            self.j = 0
        else:
            self.coff[self.i,self.j] = b
            self.j += 1

class Method():
    def __init__(self,coff,solu,sig):
        self.coff = coff
        self.sol =  solu
        self.sig = sig

    def __str__(self):
        return f"{self.coff} = {self.sol}"
    
    def pivoting(self,i,j):
        row = i
        max = self.coff[i,j]
        for k in range(i,len(self.coff)):
            if (abs(self.coff[k,j]) > max):
                max =self.coff[k,j]
                row = k
        self.coff[[i, row], :] = self.coff[[row, i], :]

        self.sol[i],self.sol[row] = self.sol[row],self.sol[i]
    
    
    def forwardElemination(self,i,j):
        self.pivoting(i,j)
        pivot = self.coff[i, i]
        if pivot == 0 and self.sol[i] == 0:
            raise ValueError("Infinite Number of Solutions")
        if pivot == 0:

          if self.step_by_step:
                print(f"pivot = 0 , No Solution (Singular matrix)")
                print(f"**** Forward elemination on a{i}{j} end ****")
            raise ValueError("No Solution (Singular matrix)")
        for k in range (i+1,len(self.coff)):
            m = self.coff[k,j]/self.coff[i,j]
            self.coff[k,j] =0
            self.sol[k] -= m*self.sol[i]
            self.sol[k]= self.sign(self.sol[k])
            for l in range(j+1,len(self.coff)):
                self.coff[k,l] -= m*self.coff[i,l]
                self.coff[k,l]= self.sign(self.coff[k,l])
        if(i == len(self.coff)-1 and j == len(self.coff)-1):
            return
        else:
            return self.forwardElemination(i+1,j+1)

    def backSub(self):

        x = np.zeros(len(self.sol))

        for k in range(len(self.coff) - 1, -1, -1):

            sum = self.sol[k]
        

            for l in range(k + 1, len(self.coff)):
                sum -= self.coff[k, l] * x[l]
        

            x[k] =self.sign( sum / self.coff[k, k])

        return x
    def forwardSub(self):

        x = np.zeros(len(self.sol))
        for k in range(len(self.coff)):
            sum = self.sol[k]
            for l in range(k):
                sum -= self.coff[k, l] * x[l]
        
            x[k] =self.sign( sum / self.coff[k, k])

        return x

    def sign (self,value):
        if value == 0:
            return 0
        magnitude = math.floor(math.log10(abs(value)))
        scale =  10 ** (self.sig -1 - magnitude)
        rounded = round(value * scale) / scale
        return rounded
    
    def sign_array(self, array):
        return np.array([self.sign(val) for val in array])