"""Create geometries."""

import gmsh
from gmsh import model as mdl
import numpy as np

def add_rectangle(dim=(10, 10)):
    return mdl.occ.add_rectangle(0,0,0,dim[0],dim[1])


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


def basic_room(pos=(5,5), angle=0):
    gmsh.initialize()
    
    mdl.add("basic_room")
    
    add_rectangle()
    add_speaker((5,5), angle)
    
    mdl.occ.synchronize()
    gmsh.write("basic_room.brep")
    
    gmsh.finalize()
