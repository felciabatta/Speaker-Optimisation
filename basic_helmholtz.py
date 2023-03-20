"""Trying to solve Helmholtz"""

from numpy import pi
from utils import get_path

filename_mesh = get_path("geomeshes/square room/square room 2D.mesh")

c = 1
f = 1/(2*pi)
w = 2*pi*f
k = 1

regions = {
    'Omega': 'all',
    'Gamma': ('all', 'facet', 'Omega')
}

fields = {
    'accoustic_pressure': ('complex', 1, 'Omega', 1),
}

variables = {
    'p': ('unknown field', 'accoustic_pressure', 0),
    'q': ('test field',    'accoustic_pressure', 'p'),
}

# options = {
#     'nls' : 'newton',
#     'ls' : 'ls',
# }

materials = {
    'air': ({'c': c},),
    'frequecy': ({'w': w},),
    'k': ({'k_squared': k**2},)
}

ebcs = {
}

integrals = {
    'i': 2,
}

equations = {
    'Acoustic pressure':
    """dw_laplace.i.Omega(q, p)
    -  dw_dot.i.Omega(k.k_squared, q, p) 
    =  0"""
}

solvers = {
    'ls': ('ls.scipy_direct', {}),
    'newton': ('nls.newton', {
        'i_max': 1,
        'eps_a': 1e-1,
        'eps_r': 1.0,
        'macheps': 1e-16,
        'lin_red': 1e-1,  # Linear system error < (eps_a * lin_red).
        'ls_red': 0.1,
        'ls_red_warp': 0.001,
        'ls_on': 1.1,
        'ls_min': 1e-5,
        'check': 0,
        'delta': 1e-6,
    })
}
