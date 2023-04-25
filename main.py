"""Run Helmholtz simulations."""

from helmholtz import geo2mesh, helmholtz_problem, plot_solution
from geometry import basic_room
from numpy import pi

def one_frequency_example(geofunc=basic_room, pos=(5,5), angle=0*pi, f=20, 
                          max_size=0.1, warp=0):
    (_, source), mesh, cells, facets = geofunc(pos, angle, max_size=max_size, 
                                             to_dolfin=True)
    ph = helmholtz_problem(mesh, cells, facets, source, f).solve()
    p_grid = plot_solution(ph, warp=warp)
    return ph, p_grid
