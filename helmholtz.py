"""Trying to solve Helmholtz, using fenicsx"""

import numpy as np
from numpy import pi
from scipy.stats import wasserstein_distance
from utils import get_path, normalize_rows, distance_mat
import ot
# meshing stuff
import gmsh
from gmsh import model as mdl
from dolfinx.io import gmshio
from mpi4py import MPI
from geometry import basic_room
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
FIGUREDIR = "figures/"

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


def helmproblem(omega_mesh, cell_markers, facet_markers, 
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


def plot_solution(ph, warp=False, save_as=False):
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
    if save_as:
        p_plotter.save_graphic(FIGUREDIR+save_as+".pdf", title=save_as)
    return p_grid


def helmsolve(geofunc=basic_room, pos=(5,5), angle=0*pi, f=[20], max_size=0.05, 
              warp=0, refineit=3, plot=False, save_as=False):
    
    # number of frequencies
    Nf = len(f)
    
    # make mesh
    (_, source), mesh, cells, facets = geofunc(
        pos, angle, max_size=max_size, to_dolfin=True, refineit=refineit)
    
    # find solutions
    ph = np.empty(Nf, dtype=object)
    for i, fi in enumerate(f):
        ph[i] = helmproblem(mesh, cells, facets, source, fi).solve()

    if plot:
        p_grid = np.empty(Nf, dtype=object)
        for i, pi in enumerate(ph):
            p_grid[i] = plot_solution(pi, warp=warp, save_as=save_as+f" {f[i]}Hz")
        return ph, p_grid
    return ph


def get_spectra(ph):
    return abs(np.column_stack([p.x.array for p in ph]))


def wasserstein_rms(array):
    array = normalize_rows(array)
    N = len(array)
    W = np.empty((N, N))
    for i, ai in enumerate(array):
        for j, aj, in enumerate(array):
            W[i, j] = wasserstein_distance(ai, aj)
    # extract above non-duplicates (upper triangle)
    tri_idx = np.triu_indices(len(W), 1)
    return np.sqrt((W[tri_idx]**2).mean())


def emd_rms(array, distance):
    array = normalize_rows(array)
    # number of samples
    N = len(array)
    # initialise emd array
    E = np.empty((N, N))
    for i, ai in enumerate(array):
        E[i, :] = ot.sinkhorn2(ai, array.T, distance, 1e-1)
    # extract above non-duplicates (upper triangle)
    tri_idx = np.triu_indices(len(E), 1)
    return np.sqrt((E[tri_idx]**2).mean())
    
