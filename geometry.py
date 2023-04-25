"""Create geometries."""

import gmsh
from gmsh import model as mdl
from dolfinx.io import gmshio
from mpi4py import MPI
# import numpy as np
from functools import wraps

def add_rectangle(dim=(10, 10)):
    # boundary coordinates 
    left = 0
    right = dim[0]
    bottom = 0
    top = dim[1]
    
    # define the corners
    lb = gmsh.model.occ.addPoint(left, bottom, 0)
    rb = gmsh.model.occ.addPoint(right, bottom, 0)
    rt = gmsh.model.occ.addPoint(right, top, 0)
    lt = gmsh.model.occ.addPoint(left, top, 0)

    # define the walls
    wall_b = gmsh.model.occ.addLine(lb, rb)
    wall_r = gmsh.model.occ.addLine(rb, rt)
    wall_t = gmsh.model.occ.addLine(rt, lt)
    wall_l = gmsh.model.occ.addLine(lt, lb)

    # connect walls
    return gmsh.model.occ.addCurveLoop([wall_b, wall_r, wall_t, wall_l])


def add_speaker(pos=(0, 0), angle=0):
    
    r = 0.2
    thickness = 0.1
    gap = 0.1
    inner = r+gap
    outer = inner+thickness
    
    # speaker box points
    lb = mdl.occ.add_point(-outer, -outer, 0)
    rb = mdl.occ.add_point(outer, -outer, 0)
    rt = mdl.occ.add_point(outer, inner, 0)
    rt_in = mdl.occ.add_point(inner, inner, 0)
    rb_in = mdl.occ.add_point(inner, -inner, 0)
    lb_in = mdl.occ.add_point(-inner, -inner, 0)
    lt_in = mdl.occ.add_point(-inner, inner, 0)
    lt = mdl.occ.add_point(-outer, inner, 0)

    # speaker box walls
    b_wall = mdl.occ.add_line(lb, rb)
    r_wall = mdl.occ.add_line(rb, rt)
    rt_wall = mdl.occ.add_line(rt, rt_in)
    r_inwall = mdl.occ.add_line(rt_in, rb_in)
    b_inwall = mdl.occ.add_line(rb_in, lb_in)
    l_inwall = mdl.occ.add_line(lb_in, lt_in)
    lt_wall = mdl.occ.add_line(lt_in, lt)
    l_wall = mdl.occ.add_line(lt, lb)
    
    wall_vec = [b_wall, r_wall, rt_wall, r_inwall, b_inwall, 
             l_inwall, lt_wall, l_wall]
    dimTags = [(1, wall) for wall in wall_vec]
    
    mdl.occ.rotate(dimTags, *(0, 0, 0), *(0, 0, 1), angle)
    mdl.occ.translate(dimTags, *pos, 0)
    circle = mdl.occ.add_circle(*pos, 0, r)
    
    # connect walls
    speaker_walls = mdl.occ.add_curve_loop(wall_vec)
    source_boundary = mdl.occ.add_curve_loop([circle])
    
    return speaker_walls, source_boundary


def refine(N):
    # NOTE: using multi refine doesn't always appear quicker, N=3 appears best
    for n in range(N):
        mdl.mesh.refine()


def init_finalize(func):
    """Add initialization and finalization processes to a geometry `func`.
    """
    @wraps(func)
    def wrapper(*args, save_as=False, **kwargs):
        gmsh.initialize()
        
        mdl.add("room")
        
        tags = func(*args, **kwargs)
        
        mdl.occ.synchronize()
        
        if save_as:
            gmsh.write(save_as+".brep")
        
        gmsh.finalize()
        return tags
    return wrapper


def add_mesh(func):
    """Add a mesh to a geometry `func`."""
    @wraps(func)
    def wrapper(*args, max_size=0.1, gdim=2, to_dolfin=False, refineit=3, 
                **kwargs):
        
        tags = func(*args, **kwargs)
        
        mdl.occ.synchronize()
        gmsh.option.set_number("Mesh.MeshSizeMax", max_size*2**(refineit+1))
        gmsh.option.set_number("Mesh.Algorithm", 5) # delaunay
        gmsh.option.set_number("General.Verbosity", 3)
        mdl.mesh.generate(gdim)
        refine(refineit)
        
        if to_dolfin:
            # define rank & comm parallelization
            gmsh_model_rank = 0
            mesh_comm = MPI.COMM_WORLD
            # convert to dolfinx mesh
            omega_mesh, cell_markers, facet_markers = gmshio.model_to_mesh(
                mdl, mesh_comm, gmsh_model_rank, gdim=gdim)
            return tags, omega_mesh, cell_markers, facet_markers
        
        return tags
    return wrapper


@init_finalize
@add_mesh
def basic_room(pos=(5,5), angle=0):
    walls = add_rectangle()
    speaker_walls, source_boundary = add_speaker(pos, angle)
    
    air_surf = mdl.occ.add_plane_surface([walls, speaker_walls, source_boundary])
    source_surf = mdl.occ.add_plane_surface([source_boundary])
    
    mdl.occ.synchronize()
    
    air_tag = mdl.add_physical_group(2, [air_surf])
    source_tag = mdl.add_physical_group(2, [source_surf])
    
    return air_tag, source_tag
