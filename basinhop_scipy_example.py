"""Scipy's basinhopping example."""

import numpy as np
from scipy.optimize import basinhopping

def func2d(x):
    f = np.cos(14.5 * x[0] - 0.3) + (x[1] + 0.2) * x[1] + (x[0] +
                                                           0.2) * x[0]
    df = np.zeros(2)
    df[0] = -14.5 * np.sin(14.5 * x[0] - 0.3) + 2. * x[0] + 0.2
    df[1] = 2. * x[1] + 0.2
    return f, df


minimizer_kwargs = {"method":"L-BFGS-B", "jac":True}
x0 = [1.0, 1.0]
ret = basinhopping(func2d, x0, minimizer_kwargs=minimizer_kwargs,
                   niter=200)
print("global minimum: x = [%.4f, %.4f], f(x) = %.4f" % (ret.x[0],
                                                          ret.x[1],
                                                          ret.fun))


class MyTakeStep:
   def __init__(self, stepsize=0.5):
       self.stepsize = stepsize
       self.rng = np.random.default_rng()
   def __call__(self, x):
       s = self.stepsize
       x[0] += self.rng.uniform(-2.*s, 2.*s)
       x[1:] += self.rng.uniform(-s, s, x[1:].shape)
       return x
   

mytakestep = MyTakeStep()
ret = basinhopping(func2d, x0, minimizer_kwargs=minimizer_kwargs,
                   niter=200, take_step=mytakestep)
print("global minimum: x = [%.4f, %.4f], f(x) = %.4f" % (ret.x[0],
                                                          ret.x[1],
                                                          ret.fun))


def print_fun(x, f, accepted):
        print("at minimum %.4f accepted %d" % (f, int(accepted)))


rng = np.random.default_rng()
ret = basinhopping(func2d, x0, minimizer_kwargs=minimizer_kwargs,
                   niter=10, callback=print_fun, seed=rng)



import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Define the function to plot
def func2d_plot(x, y):
    return np.cos(14.5 * x - 0.3) + (y + 0.2) * y + (x + 0.2) * x

# Create a 3D plot of the function
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
x = np.linspace(-2, 2, 100)
y = np.linspace(-2, 2, 100)
X, Y = np.meshgrid(x, y)
Z = func2d_plot(X, Y)
ax.plot_surface(X, Y, Z, cmap='viridis')

# Plot the global minimum found by basinhopping
ax.scatter(ret.x[0], ret.x[1], ret.fun, color='r', s=100)

# Set the labels and title
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('f(x, y)')
ax.set_title('Function and global minimum')

plt.show()
 