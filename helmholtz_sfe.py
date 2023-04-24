"""Trying to solve Helmholtz with sfepy"""

from numpy import pi
from utils import get_path

meshdir = "geomeshes"
filename_mesh = get_path("geomeshes/square room/square room.mesh")
# filename_mesh = get_path(meshdir+"/square hole/square hole.mesh")

c = 343  # speed of sound
rho = 1.293  # density of air
f = 20  # frequency
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

solvers = {
    'ls': ('ls.auto_direct', {}),
    'newton': ('nls.newton', {
        'i_max': 1,
    }),
}

options = {
    'refinement_level': 0,  # refines mesh even more
}
