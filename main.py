"""Run Helmholtz simulations."""

from helmholtz import geo2mesh, solve_helmholtz, plot_solution

def one_frequency_example(geofile="square room/square room", f=20, max_size=0.1,
                          warp=0):
    mesh, cells, facets = geo2mesh(geofile, max_size)
    ph = solve_helmholtz(mesh, cells, facets, f)
    plot_solution(ph, warp=warp)
    return ph
