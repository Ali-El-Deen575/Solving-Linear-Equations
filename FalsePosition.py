import time
from sympy import symbols, sympify, lambdify

class FalsePosition:
    def __init__(self, equation, lower, upper, epsilon=1e-5, max_iterations=50, step_by_step=False):
        self.equation = sympify(equation)
        self.lower = lower
        self.upper = upper
        self.epsilon = epsilon
        self.max_iterations = max_iterations
        self.step_by_step = step_by_step
        self.x = symbols('x')
        self.function = lambdify(self.x, self.equation)

    def solve(self):
        a, b = self.lower, self.upper
        if self.function(a) * self.function(b) > 0:
            raise ValueError("No root found. The function must have opposite signs at the interval boundaries.")

        c_old = 0
        for iteration in range(1, self.max_iterations + 1):
            fa = self.function(a)
            fb = self.function(b)
            c = (a * fb - b * fa) / (fb - fa)
            fc = self.function(c)
            if self.step_by_step:
                print(f"Iteration {iteration}: a = {a}, b = {b}, c = {c}, f(c) = {fc}")

            if abs((c-c_old)/c) < self.epsilon:
                
                if self.step_by_step:
                    print(f"Converged to {c} after {iteration} iterations")
                return c, iteration, abs((c-c_old)/c)*100
            c_old = c
            if self.function(a) * fc < 0:
                b = c
            else:
                a = c

        raise ValueError(f"Iteration did not converge within {self.max_iterations} iterations. Last interval: [{a}, {b}]")