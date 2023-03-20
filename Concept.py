import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import matplotlib.pyplot as plt

# Set up the grid
N = 50  # number of grid points in each direction
L = 10  # length of room
h = L / (N-1)  # grid spacing
x = np.linspace(0, L, N)
y = np.linspace(0, L, N)
X, Y = np.meshgrid(x, y)

# Define the source function
f = np.zeros((N, N))
f[1:15, 1:15] = 1  # point source in corner
f[1:15, 1:15] = 1  # point source in corner
f[1:15, 35:49] = 1  # point source in other corner
# f[20:49, 23:26] = 1

# Set up the finite difference matrix
A = sp.diags([-1, 2, -1], [-1, 0, 1], shape=(N-2, N-2))
I = sp.eye(N-2)
Lap = sp.kron(A, I) + sp.kron(I, A)
Lap = Lap / h**2

# Impose Dirichlet boundary conditions
u = np.zeros((N, N))
u[:, 0] = 0.0  # left boundary
u[:, -1] = 0.0  # right boundary
u[0, :] = 0.0  # top boundary
u[-1, :] = 0.0  # bottom boundary

# Solve the system of equations
b = f[1:-1, 1:-1].flatten()  # rhs of the system
u[1:-1, 1:-1] = spla.spsolve(Lap, b).reshape(N-2, N-2)

# Plot the solution
plt.contourf(X, Y, u)
plt.colorbar()
plt.show()