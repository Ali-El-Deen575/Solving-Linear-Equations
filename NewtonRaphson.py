from sympy import diff, symbols
from sympy import sympify as sp
import numpy as np
from math import factorial  
from math import log10, floor
# from Method import sign

class NewtonRaphson:
    def __init__(self, function, initial_guess , m=1, tolerance=1e-5 , sign=None, max_iterations=1000 , step_by_step=False):
        self.function = function
        self.initial_guess = initial_guess
        self.m = m
        self.tolerance = tolerance
        self.sign = sign
        self.max_iterations = max_iterations
        self.step_by_step = step_by_step
    
    def sign_func(self,value):
        if self.sign == None:
            return value
        if value == 0:
            return 0
        magnitude = floor(log10(abs(value)))
        scale =  10 ** (self.sign -1 - magnitude)
        rounded = round(value * scale) / scale
        return rounded

    def solve(self ):
        if not isinstance(self.m, int) or self.m < 1:
            raise ValueError(f"The multiplicity 'm' must be an integer greater than or equal to 1 , entered {self.m}")
        
        

        # Initialize the variables
        x = symbols('x')
        current_guess = self.initial_guess
        steps_string = ""

        # Derivatives of the function
        derivatives = [self.function]      
        for i in range(self.m):
            derivatives.append(derivatives[-1].diff())

        
        if self.step_by_step:
            steps_string +=(f"simplified equation: {self.function}\n")
            steps_string +=(f"1st derivative of the equation: {derivatives[1]}\n")
            steps_string +=(f"Initial guess: {current_guess}\n")

        for i in range(self.max_iterations):
            derivatives_values = [self.sign_func(float(derivative.subs(x, current_guess))) for derivative in derivatives]

            if self.step_by_step:
                steps_string +=(f"\nIteration {i+1}:\n\t f({current_guess}) = {derivatives_values[0]} , f'({current_guess}) = {derivatives_values[1]}\n")

            if derivatives_values[1] == 0:
                raise ValueError(f"Derivative is zero at iteration {i+1}. No solution found.")
            
            
            correction = self.sign_func(derivatives_values[0] / derivatives_values[1])
            for k in range(2, self.m + 1):
                term = self.sign_func(derivatives_values[0] ** (k - 1) * derivatives_values[k] / (factorial(k) * derivatives_values[1] ** k))
                correction *= self.sign_func((1 - term))

            next_guess = self.sign_func(current_guess - correction)


            error = self.sign_func(abs(next_guess - current_guess)/abs(next_guess))

            if self.step_by_step:
                steps_string +=(f"\t Next guess: {next_guess}\n")
                steps_string +=(f"\t Error: {error}\n")

            if error < self.tolerance:
                if self.step_by_step:
                    steps_string +=(f"\nSolution found at iteration {i+1}: {next_guess}\n")

                return next_guess , i+1 , steps_string
            current_guess = next_guess
        
        print(steps_string)
        raise ValueError("Exceeded maximum iterations. No solution found.")


if __name__ == "__main__":
    x = symbols('x')
    expression = "x**3 - 2*x"
    parsed_expr = sp(expression)
    newton = NewtonRaphson(parsed_expr, 16,  m=3 , max_iterations=25 , step_by_step=True , sign=4)
    root , iterations_num , steps = newton.solve()
    print(f"Solution: {root} , Iterations: {iterations_num}")
    print(steps)