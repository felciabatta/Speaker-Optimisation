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

# CONSTANTS
C = 343  # speed of sound
RHO = 1.293  # density of air
MESHDIR = "geomeshes/"

def geo2mesh(geofile="square room/square room", max_size=0.1):

    geo_path = get_path(MESHDIR+geofile+".geo")
    gdim = 2

    gmsh.initialize()
    gmsh.open(geo_path)
    gmsh.model.occ.synchronize()
    gmsh.option.setNumber("Mesh.MeshSizeMax", max_size)
    gmsh.model.mesh.generate(gdim)

    # convert to dolfinx mesh
    gmsh_model_rank = 0
    mesh_comm = MPI.COMM_WORLD
    omega_mesh, cell_markers, facet_markers = gmshio.model_to_mesh(
        gmsh.model, mesh_comm, gmsh_model_rank, gdim=gdim)

    gmsh.finalize()
    
    return omega_mesh, cell_markers, facet_markers

def solve_helmholtz(omega_mesh, cell_markers, facet_markers=None, f=20):
    # define function space (field)
    P = FunctionSpace(omega_mesh, ("Lagrange", 1))
    P_source = FunctionSpace(omega_mesh, ("DG", 0))

    # define unknown and test variables
    p = ufl.TrialFunction(P)
    q = ufl.TestFunction(P)

    # constants
    w = 2*pi*f  # radial frequency
    k = w/C  # wave number
    Q_strength = 1 # monopole source strength
    Q = Function(P_source)
    source_cells = cell_markers.find(2)
    Q.x.array[source_cells] = Q_strength

    # define weak form
    a = ufl.inner(ufl.grad(p), ufl.grad(q))*ufl.dx - k**2 * ufl.inner(p, q) * ufl.dx
    L = RHO * ufl.inner(Q, q) * ufl.dx

    # solve
    problem = LinearProblem(a, L, petsc_options={"ksp_type": "preonly", "pc_type": "lu"})
    return problem.solve()

def plot_solution(ph, warp=False):
    P = ph.ufl_function_space()
    
    # SET UP SOLUTION PLOT DATA
    # make grid
    u_topology, u_cell_types, u_geometry = dolfinx.plot.create_vtk_mesh(P)
    u_grid = pyvista.UnstructuredGrid(u_topology, u_cell_types, u_geometry)
    # add point data
    u_grid.point_data["p"] = ph.x.array.real
    u_grid.set_active_scalars("p")
    if warp:
        u_grid = u_grid.warp_by_scalar(factor=warp)

    # PLOT SOLUTION
    # initialise plotter
    u_plotter = pyvista.Plotter()
    # add data
    u_plotter.add_mesh(u_grid, show_edges=0, cmap='seismic')
    # view
    u_plotter.view_xy()
    u_plotter.show()
