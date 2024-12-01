
import numpy as np
import time
from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
from PyQt5.QtGui import *
from GaussElemination import GaussElemination
from Gaussjordan import GaussJordan
from LU import LU
from Jacobi import Jacobi
from GaussSeidel import GaussSeidel
from Method import Equations  
from PyQt5.uic import loadUiType
from os import path 
import sys 
        
mainWindowFileName = "test.ui"                
FORM_CLASS, _ = loadUiType(path.join(path.dirname(__file__), mainWindowFileName))
    
#app = QApplication(sys.argv)
         

class Ui_MainWindow(QMainWindow,FORM_CLASS):
    def __init__(self, parent=None):
        super(Ui_MainWindow, self).__init__(parent)
        QMainWindow.__init__(self)
        self.setupUi(self)
        self.system=Equations(0)
        self.scaling = False
        self.EnterNo.clicked.connect(self.setNoEqn)
        self.EnterEqn.clicked.connect(self.setEqn)
        self.FindSol.clicked.connect(self.findSol)
        self.Clear.clicked.connect(self.clear)
        self.EnterSigFigures.clicked.connect(self.setSigNo)
        self.method.currentTextChanged.connect(self.parameters)
        self.ScalingCheckBox.stateChanged.connect(self.toggleScaling)
        self.parameters()

    def toggleScaling(self, state):
        if state == QtCore.Qt.Checked:
            self.scaling = True

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
        if self.scaling:
            self.system.coff, self.system.sol = self.scale_matrix(self.system.coff,self.system.sol)

        print(self.system.coff)
        print(self.system.sol)
        if(self.method.currentText()=="Gauss Elimination"):
            gaussElem = GaussElemination(self.system.coff,self.system.sol,self.system.sig)
            startTime = time.time()
            try:
                gaussElem.forwardElemination(0,0)
            except ValueError as e:
                self.result.setText(f"{e}")
            else:
                res = gaussElem.apply()
                EndTime = time.time()
                for i in range(len(res)):
                    self.result.setText(self.result.toPlainText()+f"X{i+1} = {res[i]}\n")
                self.time.setText(f"{EndTime - startTime}")

        
        elif(self.method.currentText()=="Gauss Jordan"):
            gaussJor = GaussJordan(self.system.coff,self.system.sol,self.system.sig)
            startTime = time.time()
            try:
                gaussJor.forwardElimination()
            except ValueError as e:
                self.result.setText(f"{e}")
            else:
                res = gaussJor.apply()
                EndTime = time.time()
                for i in range(len(res)):
                    self.result.setText(self.result.toPlainText()+f"X{i+1} = {res[i]}\n")
                self.time.setText(f"{EndTime - startTime}")

        elif(self.method.currentText()=="LU decompostion"):
            method = str(self.Parameters.currentText())
            lu = LU(self.system.coff,self.system.sol,method,self.system.sig)
            startTime = time.time()
            try:
                res = lu.apply()
            except ValueError as e:
                self.result.setText(f"{e}")
            else:
                EndTime = time.time()
                for i in range(len(res)):
                    self.result.setText(self.result.toPlainText()+f"X{i+1} = {res[i]}\n")
                self.time.setText(f"{EndTime - startTime}")
        
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
                        res,it = jacobi.apply()
                        EndTime = time.time()
                        for i in range(len(res)):
                            self.result.setText(self.result.toPlainText()+f"X{i+1} = {res[i]}\n")
                        self.time.setText(f"{EndTime - startTime}")
                        self.Iterations.setText(f"{it}")     
        
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
                        res,it = sidel.apply()
                        EndTime = time.time()
                        for i in range(len(res)):
                            self.result.setText(self.result.toPlainText()+f"X{i+1} = {res[i]}\n")
                        self.time.setText(f"{EndTime - startTime}")
                        self.Iterations.setText(f"{it}")                               

    def scale_matrix(self, A, B):
        n = A.shape[0]
        for i in range(n):
            row_max = max(abs(A[i, j]) for j in range(n))
            if row_max != 0:
                A[i, :] = A[i, :] / row_max
                B[i] = B[i] / row_max
        return A, B

    def clear(self):
        self.system=Equations(0)
        self.No.setText("n=")
        self.NoEqn.setText("")
        self.var.setText("")
        self.result.setText("")
        self.Eqn.setText("")
        self.InitialGuess.setText("")
        self.IterationNumber.setText("")
        self.StoppingCondition.setText("")
        self.time.setText("")
        self.Iterations.setText("")
        self.SigFigures.setText("")
        self.result.setText("")
    
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
    
def loadStylesheet():
    try:
        with open("StyleSheet.qss", "r") as file:
            stylesheet = file.read()
            QApplication.instance().setStyleSheet(stylesheet)
    except FileNotFoundError:
        print("Stylesheet file not found.")


def main():
    app = QApplication(sys.argv)
    loadStylesheet()
    window = Ui_MainWindow()
    window.show()
    app.exec_()

if __name__ == '__main__':
    main()
