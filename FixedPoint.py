import math

from sympy import lambdify, symbols, sympify

class FixedPoint:
    def __init__(self, g, x0, tol=1e-5, max_iter=100, step_by_step=False):
        self.g = sympify(g)
        self.x0 = x0
        self.tol = tol
        self.max_iter = max_iter
        self.step_by_step = step_by_step
        self.x = symbols('x')
        self.g = lambdify(self.x, self.g)

    def apply(self):
        x = self.x0
        iter_count = 0
        ea = 100.0

        if self.step_by_step:
            print("**** Fixed Point Iteration start ****")
            print(f"Initial guess: x0 = {x}")

        while ea > self.tol and iter_count < self.max_iter:
            x_old = x
            try:
                x = self.g(x_old)
            except (OverflowError, ZeroDivisionError) as e:
                raise ValueError(f"Error during evaluation of g(x): {e}")
            except Exception as e:
                raise ValueError(f"Unexpected error during evaluation of g(x): {e}")

            if x != 0:
                ea = abs((x - x_old) / x) 

            iter_count += 1

            if self.step_by_step:
                print(f"Iteration {iter_count}: x = {x}, relative error = {ea}")

            if math.isinf(x) or math.isnan(x):
                raise ValueError(f"Iteration diverged at step {iter_count}: x = {x}")

        if ea <= self.tol:
            if self.step_by_step:
                print(f"Converged to {x} after {iter_count} iterations")
                print("**** Fixed Point Iteration end ****")
            return x, iter_count ,ea
        else:
            raise ValueError(f"Iteration did not converge within {self.max_iter} iterations. Last value: {x}")