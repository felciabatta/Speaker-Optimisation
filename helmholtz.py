"""Trying to solve Helmholtz, using fenicsx"""

from numpy import pi
from utils import get_path
# meshing stuff
import gmsh
from dolfinx.io import gmshio
from mpi4py import MPI
# fem stuff
from dolfinx.fem import FunctionSpace, Function
from dolfinx.fem.petsc import LinearProblem
# problem definition
import ufl
# plotting
import pyvista
import dolfinx.plot

meshdir = "geomeshes"
geo_filename = get_path("geomeshes/square room/square room.geo")
# geo_filename = get_path("geomeshes/experimental/enclosed speaker.geo")

gmsh.initialize()
gmsh.open(geo_filename)
gmsh.model.occ.synchronize()
gdim = 2
gmsh.option.setNumber("Mesh.MeshSizeMax", 0.1)
gmsh.model.mesh.generate(gdim)

# convert to dolfinx mesh
gmsh_model_rank = 0
mesh_comm = MPI.COMM_WORLD
omega, cell_markers, facet_markers = gmshio.model_to_mesh(
    gmsh.model, mesh_comm, gmsh_model_rank, gdim=gdim)

gmsh.finalize()

# define function space (field)
P = FunctionSpace(omega, ("Lagrange", 1))
P_source = FunctionSpace(omega, ("DG", 0))

# define unknown and test variables
p = ufl.TrialFunction(P)
q = ufl.TestFunction(P)

# constants
c = 343  # speed of sound
rho = 1.293  # density of air
f = 20  # frequency
w = 2*pi*f  # radial frequency
k = w/c  # wave number
Q_strength = 1 # monopole source strength
Q = Function(P_source)
source_cells = cell_markers.find(2)
Q.x.array[source_cells] = Q_strength

# define weak form
a = ufl.inner(ufl.grad(p), ufl.grad(q))*ufl.dx - k**2 * ufl.inner(p, q) * ufl.dx
L = rho * ufl.inner(Q, q) * ufl.dx

# solve
problem = LinearProblem(a, L, petsc_options={"ksp_type": "preonly", "pc_type": "lu"})
uh = problem.solve()

# SET UP SOLUTION PLOT DATA
# make grid
u_topology, u_cell_types, u_geometry = dolfinx.plot.create_vtk_mesh(P)
u_grid = pyvista.UnstructuredGrid(u_topology, u_cell_types, u_geometry)
# add point data
u_grid.point_data["u"] = uh.x.array.real
u_grid.set_active_scalars("u")
warped = u_grid.warp_by_scalar()

# PLOT SOLUTION
# initialise plotter
u_plotter = pyvista.Plotter()
# add data
u_plotter.add_mesh(u_grid, show_edges=0, cmap='seismic')
# view
u_plotter.view_xy()
u_plotter.show()
