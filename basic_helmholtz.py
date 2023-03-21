"""Trying to solve Helmholtz"""

from numpy import pi
from utils import get_path

meshdir = "geomeshes"
# filename_mesh = get_path("geomeshes/square room/square room 2D.mesh")
filename_mesh = get_path(meshdir+"/square hole/square hole.mesh")

c = 343  # speed of sound
rho = 1.293  # density of air
f = 800  # frequency
w = 2*pi*f  # radial frequency
k = w/c  # wave number
Q = 1000  # monopole source strength

regions = {
    'Omega': 'all',
    'Gamma': ('vertices of group 1 +v vertices of group 2 +v vertices of group 3', 'facet'),
    'Gamma_d': ('vertices of group 1', 'facet'),
    # BUG: dont think the vertex group work properly
    'LR': ('vertices in (x > 9.99) +v vertices in (x<0.01)', 'facet'),
    'Source': ('cells of group 2')
}

fields = {
    'acoustic_pressure': ('real', 1, 'Omega', 1),
}

variables = {
    'p': ('unknown field', 'acoustic_pressure', 0),
    'q': ('test field',    'acoustic_pressure', 'p'),
}

materials = {
    'air': ({'c': c,
             'rho': rho},),

    'freq': ({'f': f,
              'w': w},),

    'wave_num': ({'k_square': k**2},),

    'source': ({'mono': rho*Q},)
}

ebcs = {}  # {'p0': ('TB', {'p.all': 0.0})}

integrals = {
    'i': 2,
}

equations = {
    'Acoustic pressure':
    """
    - dw_laplace.i.Omega(q, p)
    + dw_dot.i.Omega(wave_num.k_square, q, p)
    + dw_integrate.i.Source(source.mono, q)
    = 0
    """
}

# solvers = {
#     'ls': ('ls.scipy_direct', {}),
#     'newton': ('nls.newton', {
#         'i_max': 1,
#         'eps_a': 1e-1,
#         'eps_r': 1.0,
#         'macheps': 1e-16,
#         'lin_red': 1e-1,  # Linear system error < (eps_a * lin_red).
#         'ls_red': 0.1,
#         'ls_red_warp': 0.001,
#         'ls_on': 1.1,
#         'ls_min': 1e-5,
#         'check': 0,
#         'delta': 1e-6,
#     })
# }

# options = {
#     'nls' : 'newton',
#     'ls' : 'ls',
# }

solvers = {
    'ls': ('ls.auto_direct', {}),
    'newton': ('nls.newton', {
        'i_max': 1,
    }),
}

options = {
    'refinement_level': 0,  # refines mesh even more
}
