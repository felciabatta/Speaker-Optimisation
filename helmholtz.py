"""Trying to solve Helmholtz, using fenicsx"""

from numpy import pi
from utils import get_path
# meshing stuff
import gmsh
from gmsh import model as mdl
from dolfinx.io import gmshio
from mpi4py import MPI
# fem stuff
from dolfinx.fem import FunctionSpace, Function
from dolfinx.fem.petsc import LinearProblem
# problem definition
from ufl import TrialFunction, TestFunction, inner, grad, dx
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
    mdl.occ.synchronize()
    gmsh.option.setNumber("Mesh.MeshSizeMax", max_size)
    mdl.mesh.generate(gdim)
    mdl.mesh.refine()

    # convert to dolfinx mesh
    gmsh_model_rank = 0
    mesh_comm = MPI.COMM_WORLD
    omega_mesh, cell_markers, facet_markers = gmshio.model_to_mesh(
        mdl, mesh_comm, gmsh_model_rank, gdim=gdim)

    gmsh.finalize()
    
    return omega_mesh, cell_markers, facet_markers


def helmholtz_problem(omega_mesh, cell_markers, facet_markers, 
                      source_tag=2, f=20):
    # define function space (field)
    P = FunctionSpace(omega_mesh, ("Lagrange", 1))
    P_source = FunctionSpace(omega_mesh, ("DG", 0))

    # define unknown and test variables
    p = TrialFunction(P)
    q = TestFunction(P)

    # constants
    w = 2*pi*f  # radial frequency
    k = w/C  # wave number
    Q_strength = 1 # monopole source strength
    Q = Function(P_source)
    source_cells = cell_markers.find(source_tag)
    Q.x.array[source_cells] = Q_strength

    # define weak form
    a = inner(grad(p), grad(q))*dx - k**2 * inner(p, q) * dx
    L = RHO * inner(Q, q) * dx

    # formulate problem
    return LinearProblem(
        a, L, petsc_options={"ksp_type": "preonly", "pc_type": "lu"})


def plot_solution(ph, warp=False):
    P = ph.ufl_function_space()
    
    # SET UP SOLUTION PLOT DATA
    # make grid
    p_topology, p_cell_types, p_geometry = dolfinx.plot.create_vtk_mesh(P)
    p_grid = pyvista.UnstructuredGrid(p_topology, p_cell_types, p_geometry)
    # add point data
    p_grid.point_data["p"] = ph.x.array.real
    p_max = abs(p_grid.point_data["p"]).max()
    p_grid.set_active_scalars("p")
    if warp:
        scale = warp/p_max
        p_grid = p_grid.warp_by_scalar(factor=scale)

    # PLOT SOLUTION
    # initialise plotter
    p_plotter = pyvista.Plotter()
    # add data
    p_plotter.add_mesh(
        p_grid, show_edges=0, cmap='seismic', clim=[-p_max, p_max])
    # view
    p_plotter.view_xy()
    p_plotter.show()
    
    return p_grid
