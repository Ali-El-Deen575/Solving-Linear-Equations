import math

class FixedPoint:
    def __init__(self, g, x0, tol=1e-5, max_iter=100, step_by_step=False):
        self.g = g
        self.x0 = x0
        self.tol = tol
        self.max_iter = max_iter
        self.step_by_step = step_by_step

    def apply(self):
        x = self.x0
        if self.step_by_step:
            print("**** Fixed Point Iteration start ****")
            print(f"Initial guess: x0 = {x}")

        for i in range(self.max_iter):
            try:
                x_new = self.g(x)
            except (OverflowError, ZeroDivisionError) as e:
                raise ValueError(f"Error during evaluation of g(x): {e}")
            except Exception as e:
                raise ValueError(f"Unexpected error during evaluation of g(x): {e}")

            error = abs(x_new - x) / max(abs(x_new), 1.0)
            if self.step_by_step:
                print(f"Iteration {i+1}: x_new = {x_new}, relative error = {error}")

            if error < self.tol:
                if self.step_by_step:
                    print(f"Converged to {x_new} after {i+1} iterations")
                    print("**** Fixed Point Iteration end ****")
                return x_new, i + 1

            if math.isinf(x_new) or math.isnan(x_new):
                raise ValueError(f"Iteration diverged at step {i+1}: x_new = {x_new}")

            x = x_new

        raise ValueError(f"Iteration did not converge within {self.max_iter} iterations. Last value: {x}")