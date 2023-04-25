"""Run Helmholtz simulations."""

from helmholtz import geo2mesh, helmholtz_problem, plot_solution
from geometry import basic_room
import numpy as np
from numpy import pi

def helmsolve(geofunc=basic_room, pos=(5,5), angle=0*pi, f=[20], max_size=0.05, 
              warp=0, refineit=3, plot=False):
    
    # number of frequencies
    Nf = len(f)
    
    # make mesh
    (_, source), mesh, cells, facets = geofunc(
        pos, angle, max_size=max_size, to_dolfin=True, refineit=refineit)
    
    # find solutions
    ph = np.empty(Nf, dtype=object)
    for i, fi in enumerate(f):
        ph[i] = helmholtz_problem(mesh, cells, facets, source, fi).solve()

    if plot:
        p_grid = np.empty(Nf, dtype=object)
        for i, pi in enumerate(ph):
            p_grid[i] = plot_solution(pi, warp=warp)
        return ph, p_grid
    return ph    
