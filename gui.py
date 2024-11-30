
import numpy as np
import time
from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
from PyQt5.QtGui import *
from GaussElemination import GaussElemination
from Gaussjordan import GaussJordan
from Jacobi import Jacobi
from GaussSeidel import GaussSeidel
from Method import Equations  
from PyQt5.uic import loadUiType
from os import path 
import sys 
        
mainWindowFileName = "test.ui"                
FORM_CLASS, _ = loadUiType(path.join(path.dirname(__file__), mainWindowFileName))
    
app = QApplication(sys.argv)
         

class Ui_MainWindow(QMainWindow,FORM_CLASS):
    def __init__(self, parent=None):
        super(Ui_MainWindow, self).__init__(parent)
        QMainWindow.__init__(self)
        self.setupUi(self)
        self.system=Equations(0)
        self.EnterNo.clicked.connect(self.setNoEqn)
        self.EnterEqn.clicked.connect(self.setEqn)
        self.FindSol.clicked.connect(self.findSol)
        self.Clear.clicked.connect(self.clear)
        self.EnterSigFigures.clicked.connect(self.setSigNo)
        self.method.currentTextChanged.connect(self.parameters)
        self.parameters()

   
    def parameters (self):
        if(self.method.currentText()=="Gauss Elimination" or self.method.currentText()=="Gauss Jordan") :
            self.ParametersLabel.setEnabled(False)
            self.Parameters.setEnabled(False)
            self.InitialGuess.setEnabled(False)
            self.InitiaGuessLabel.setEnabled(False)
            self.IterationNumber.setEnabled(False)
            self.StoppingCondition.setEnabled(False)
            self.IterationLabel.setEnabled(False)
            self.StoppingConditionLabel.setEnabled(False)

        elif(self.method.currentText()=="LU decompostion") :
            self.ParametersLabel.setEnabled(True)
            self.Parameters.setEnabled(True)
            self.InitialGuess.setEnabled(False)
            self.InitiaGuessLabel.setEnabled(False)
            self.IterationNumber.setEnabled(False)
            self.StoppingCondition.setEnabled(False)
            self.IterationLabel.setEnabled(False)
            self.StoppingConditionLabel.setEnabled(False)

  
        elif(self.method.currentText()=="Jacobi" or self.method.currentText()=="Gauss sidel") :
            self.ParametersLabel.setEnabled(False)
            self.Parameters.setEnabled(False)
            self.InitialGuess.setEnabled(True)
            self.InitiaGuessLabel.setEnabled(True)
            self.IterationNumber.setEnabled(True)
            self.StoppingCondition.setEnabled(True)
            self.IterationLabel.setEnabled(True)
            self.StoppingConditionLabel.setEnabled(True)
    
    

    def is_number(self,s):
        try:
          float(s)
          return True
        except ValueError:
          return False
    def setNoEqn (self): 
      if(self.isWholeNumber(self.NoEqn.text())) :      
        eqnNo = int(self.NoEqn.text())
        self.No.setText("n= "+str(eqnNo)) 
        self.var.setText("a11 = ")
        self.system =Equations(eqnNo)
    def setSigNo (self): 
      if(self.isWholeNumber(self.SigFigures.text())) :      
        SigNo = int(self.SigFigures.text())
        self.system.sig = SigNo    
    
    def setEqn(self):
        if(self.is_number(self.Eqn.text())) :
            if(self.system.i<self.system.num):
                self.system.setCoff(float(self.Eqn.text()))
                self.Eqn.setText("")
                if(self.system.j < self.system.num and self.system.i < self.system.num):
                    self.var.setText(("a"+str(self.system.i+1)+str(self.system.j+1)+" ="))
                    
                elif(self.system.i < self.system.num) :               
                    self.var.setText("b"+str(self.system.i+1)+" =")
                    
                else:
                    self.var.setText("")

    def findSol(self):
        if(self.method.currentText()=="Gauss Elimination"):
            gaussElem = GaussElemination(self.system.coff,self.system.sol,self.system.sig)
            startTime = time.time()
            if(gaussElem.forwardElemination(0,0)):
                res = gaussElem.apply()
                EndTime = time.time()
                for i in range(len(res)):
                    self.result.setText(self.result.text()+f"X{i+1} = {res[i]}\n")
                self.time.setText(f"{EndTime - startTime}")    
            else:
                self.result.setText("System has No Solution OR infinite Number of Solutions")  
        
        elif(self.method.currentText()=="Gauss Jordan"):
            gaussJor = GaussJordan(self.system.coff,self.system.sol,self.system.sig)
            startTime = time.time()
            if(gaussJor.forwardElimination()):
                res = gaussJor.apply()
                EndTime = time.time()
                for i in range(len(res)):
                    self.result.setText(self.result.text()+f"X{i+1} = {res[i]}\n")
                self.time.setText(f"{EndTime - startTime}")    
            else:
                self.result.setText("System has No Solution OR infinite Number of Solutions")          
        
        elif(self.method.currentText()=="LU Decompostion"):
            gauss = GaussElemination(self.system.coff,self.system.sol,self.system.sig)
            res = gauss.apply()
            for i in range(len(res)):
                self.result.setText(self.result.text()+f"X{i+1} = {res[i]}\n")
        
        elif(self.method.currentText()=="Jacobi"):
            guess = self.InitialGuess.text()
            Guess = None if guess == "" else self.extractNumbers(guess) 
            iterations =self.IterationNumber.text()
            error = self.StoppingCondition.text()           
            if not(iterations == "" and error == ""):
                if(self.isWholeNumber(iterations) or iterations == "") and (self.is_number(error) or error == ""):
                    if(Guess is None or (all(isinstance(x, (int, float)) for x in Guess) and len(Guess) == len(self.system.sol))):
                        Iteration = None if iterations =="" else int(iterations)
                        Error = None if error =="" else float(error)
                        startTime = time.time()
                        jacobi = Jacobi(self.system.coff,self.system.sol,Guess,Iteration,Error,self.system.sig)
                        res = jacobi.apply()
                        EndTime = time.time()
                        for i in range(len(res)):
                            self.result.setText(self.result.text()+f"X{i+1} = {res[i]}\n")
                        self.time.setText(f"{EndTime - startTime}")     
        
        elif(self.method.currentText()=="Gauss sidel"):
            guess = self.InitialGuess.text()
            Guess = None if guess == "" else self.extractNumbers(guess) 
            iterations =self.IterationNumber.text()
            error = self.StoppingCondition.text()           
            if not(iterations == "" and error == ""):
                if(self.isWholeNumber(iterations) or iterations == "") and (self.is_number(error) or error == ""):
                    if(Guess is None or (all(isinstance(x, (int, float)) for x in Guess) and len(Guess) == len(self.system.sol))):
                        Iteration = None if iterations =="" else int(iterations)
                        Error = None if error =="" else float(error)
                        startTime = time.time()
                        sidel = GaussSeidel(self.system.coff,self.system.sol,Guess,Iteration,Error,self.system.sig)
                        res = sidel.apply()
                        EndTime = time.time()
                        for i in range(len(res)):
                            self.result.setText(self.result.text()+f"X{i+1} = {res[i]}\n")
                        self.time.setText(f"{EndTime - startTime}")                               

    def clear(self):
        self.system=Equations(0)
        self.No.setText("")
        self.var.setText("")
        self.result.setText("")
        self.Eqn.setText("")
        self.NoEqn.setText("")
    
    def isWholeNumber(self, s):
        try:
            number = float(s)  
            if number % 1 == 0: 
                return True
            else:
                return False
        except ValueError:
            return False
    def extractNumbers(self, s):
        try:
            numbers = s.split(',')
            numbers = [float(num) for num in numbers]           
            return numbers
        except ValueError:
            return [] 
    
                                
            

def main():   
    window = Ui_MainWindow() 
    window.show()
    app.exec_()

if __name__ == '__main__':
    main()
