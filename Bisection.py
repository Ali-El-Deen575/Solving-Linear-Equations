import time
from sympy import symbols, sympify, lambdify

class Bisection:
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

        start_time = time.time()
        for iteration in range(1, self.max_iterations + 1):
            c = (a + b) / 2
            fc = self.function(c)
            if self.step_by_step:
                print(f"Iteration {iteration}: a = {a}, b = {b}, c = {c}, f(c) = {fc}")

            if abs(fc) < self.epsilon or abs(b - a) < self.epsilon:
                execution_time = time.time() - start_time
                if self.step_by_step:
                    print(f"Converged to {c} after {iteration} iterations")
                return c, iteration, abs(b - a), execution_time

            if self.function(a) * fc < 0:
                b = c
            else:
                a = c

        raise ValueError(f"Iteration did not converge within {self.max_iterations} iterations. Last interval: [{a}, {b}]")