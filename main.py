## <one line to give the program's name and a brief idea of what it does.>
##     Copyright (C) 2018-2020  Alexis Boris Jawtuschenko.
## 
## This file is part of ?????.
## 
## ????? is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
## 
## ????? is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU General Public License
## along with ?????.  If not, see <https://www.gnu.org/licenses/>.

import numpy as np
import mesh
import mesh_write
import mesh_connectivity
import time
import mesh_graphics


"""
The dictionary of parsers works like a switch with three default cases:

default 0: macroel_hybrid
default 1: macroel_tetrahedra
default 2: macroel_prisms

"""

parsers                         = { 3: mesh.split_parallel_brick_into_tetrahedra,  # Only for bricks parallel to the coordinates!
                                    4: mesh.split_preordered_hexahedron_into_tetrahedra,
                                    5: mesh.split_brick_into_prisms,
                                    6: mesh.split_lshape_into_prisms,
                                    7: mesh.reorder_prism_vertices,
                                    8: mesh.morph_flexible_polyhedron_with_control_points }

# CONTINUE HERE: make a handy way of typing the initial partition,
# perhaps a class with functions like "tau_zero.translate(x=x0, y=y0, z=z0)"


## TODO: split_parallel_brick_into_tetrahedra works fine if the edges of the
## brick are paralell to the x,y,z axes. Research for an algorithm that allows
## to put the vertices in any order even for oblique hexahedra.
## Then, continue with parsers 5 and 6.
## ## after finishing all that, make cute examples to put
## in .md file and document or go to mesh.macroel_tetrahedra

elements_by_vertices_writers    = { 0: mesh_connectivity.elements_by_vertices_hybrid,
                                    1: mesh_connectivity.elements_by_vertices_tetra,
                                    2: mesh_connectivity.elements_by_vertices_prisms }

local_meshers                   = { 0: mesh.macroel_hybrid, 
                                    1: mesh.macroel_tetrahedra, 
                                    2: mesh.macroel_prisms }

physical_vertices_writers       = { 0: mesh_write.vertices_macro_hybrid,
                                    1: mesh_write.vertices_macro_tetra,
                                    2: mesh_write.vertices_macro_prism } 

def load_partition (in_file, levels):
    """ in_file is a csv with:
    macroelement_type, macro_vertices, local_mu 

    The returned macro_elements list is built with parsers
    or with a default lambda function for the case we put
    any of our original macro--elements directly.
    
    colors = [ "green", "red", "blue", "purple"]
    """
    with open(in_file,'r') as infile:
        inlist = infile.readlines()
    pre_list = [line.strip(' \n').split(',') for line in inlist\
                                               if (line not in ['','\n','\t'])]
    pre_list = [[int(float(st[0]))] + [float(st[k]) for k in range(1,len(st))] for st in pre_list]
    macro_elements = []
    for macro in pre_list:
        macro_elements += parsers.get(macro[0], lambda l: [{\
                                  0 : l[0],
                                  1 : np.array(l[1:-1]).reshape(3,(len(l)-2)//3,order='F'),
                                  2 : l[-1] }])(macro)
    ## TODO: learn how to do a singleton for this:
    # If our initial partition contains an hybrid element, write the graph topology.
    if (0 in [e[0] for e in macro_elements]):
        mesh_write.write_element_indices(in_file+".elem", levels)
    return macro_elements

def filter_repeated_vertices (in_file,n_vert_prism = 6):
    """	look for repetitions of vertices.
    at this point the program is already general:
    this 'inverts' the table of elements_by_vertices.txt
    writes on disc: vertices_by_elements.txt
    mesh_connectivity.vertices_by_elements('elements_by_vertices.txt', 'Octave')
    """
    with open(in_file+'.ebvr','r') as inp:
        things = inp.readlines()
    
    n_elem = len(things)
    
    elem_vert_repeated = np.zeros((n_elem,n_vert_prism+1),dtype=int)

    for k in range(len(things)):
        ele                     = np.fromstring(things[k],dtype=int,sep=' ')
        elem_vert_repeated[k]   = np.concatenate((ele,np.zeros((n_vert_prism+1-len(ele)),dtype=int)))
    
    print ('mesh_connectivity.kill_repeated()')
    t0 = time.time()
    replace_verts = mesh_connectivity.kill_repeated(in_file+'.ver')  ## repetitions in vertices.txt leave a dictionary
    print (time.time() - t0)
    print ('\n')
    
    print ('vertices replacement loop version 2')
    t0      = time.time()
    counter = 1
    elem_vert_dscnt_indcs = np.copy(elem_vert_repeated[:,1:]).reshape(n_elem*6)
    for key in replace_verts:
        print('progress: {}/{}\r'.format(counter,len(replace_verts)),sep=' ',\
                end='',flush=True)
        counter += 1
        for r in replace_verts[key]:
            elem_vert_dscnt_indcs[elem_vert_dscnt_indcs == r+1] = (key +1)
    print ('\r')
    print (time.time() - t0)
    print ('\n')
    elem_vert_dscnt_indcs = elem_vert_dscnt_indcs.reshape((n_elem,6))
    col = elem_vert_repeated[:,0].reshape(elem_vert_repeated.shape[0],1)
    elem_vert_dscnt_indcs = np.hstack((col,elem_vert_dscnt_indcs))
    del elem_vert_repeated
    ## <--- whithout repetitions
    with open(in_file+'.ebv','wb') as out:
        for e in elem_vert_dscnt_indcs:
            np.savetxt(out, e[0:e[0]+1].reshape((1,e[0]+1)), fmt='%d')
    return n_elem

def filter_repeated_faces (in_file,n_elem):
    """ Takes the output of mesh_connectivity.vertices_by_elements().
    writes on disc: in_file.nf          ----> number of faces
                    in_file.faces_rep 	----> faces indices with repetitions
                    in_file.fltg        ----> fltg stands for a
                                              "faces_local_to_global" correspondence
                    in_file.ebf         ----> elements_by_faces correspondence """	
    print ('Face Enumeration')
    mesh_connectivity.face_enumeration(in_file)
    print ('\r')
    with open(in_file+'.fltg', 'r') as ent:
        stuff = ent.readlines()
    elem_faces_repeated = np.zeros((n_elem, 6),dtype=int)
    for elem in range(len(stuff)):
        el                          = np.fromstring(stuff[elem],dtype=int,sep=' ')
        elem_faces_repeated[elem]   = np.concatenate((el,np.zeros((6-len(el)),dtype=int)))

    ## reads repetitions in faces_repeated.txt and leaves a dictionary
    ## writes file: faces.txt -----> global unique enumeration of faces
    print('Killing repeated faces')
    replace_faces, indices, num_faces = mesh_connectivity.kill_repeated_faces(in_file)  
    print( '\r')
    with open (in_file+'.nf','w') as n_of_faces:
        n_of_faces.write(str(num_faces))
    #### uniquifying faces:
    print('Face replacements loop version 2.')
    counter = 1
    # discontinuous indices
    elem_faces = np.copy(elem_faces_repeated[:,1:]).reshape(n_elem*5)
    first_col = elem_faces_repeated[:,0].reshape(elem_faces_repeated.shape[0],1)
    for key in replace_faces:
        print('progress: {0}/{1}\r'.format(counter,len(replace_faces)), sep = ' ', end = '', flush = True)
        counter += 1
        for r in replace_faces[key]:
            elem_faces[elem_faces==r] = key
    print ('\r')
    del elem_faces_repeated
    ## bijection indices ----> [1,...,n_faces]
    ## unique, continuous indices with filling zeros
    print ('loop for indexing faces with {1 ... n_faces}')    
    for i in range(len(indices)):
        elem_faces[elem_faces==indices[i]] = i+1
    print ('\r')
    elem_faces = np.hstack((first_col,elem_faces.reshape(n_elem,5)))
    with open(in_file+'.ebf','ab') as ex:
        for elem in elem_faces:
            np.savetxt(ex, elem.reshape((1,6)),fmt='%d')
    return 

def mesh_domain(in_file, levels=1):
    """
    1st: set levels (if levels == 1, then only macro--elements).
    elements_by_vertices_writers:    write elements_by_vertices_repeated.txt, GLOBAL INDICES per element
    physical_vertices_writers:         write vertices.txt, global list of vertices
    """
    print('<program>  Copyright (C) 2018-2022  Alexis Boris Jawtuschenko\n' +
          'This program comes with ABSOLUTELY NO WARRANTY; for details type ????.\n' +
          'This is free software, and you are welcome to redistribute it\n' +
          'under certain conditions; type show c for details.')

    initial_partition = load_partition(in_file, levels)
    init = 0
    for E in initial_partition:
        # for case E[0] == 1: the following writes contiguous indices with repetitions on 
        # shared faces.
        elements_by_vertices_writers[E[0]](in_file, levels, "Octave", init)

        # for case E[0] == 1: the following calculates coordinates with repetitions on shared 
        # faces, with the backtracking for tetrahedra. 
        points = local_meshers[E[0]](E[1],E[2],levels)

        # for case E[0] = 1: the following writes coordinates with repetitions on shared 
        # faces, with the backtracking for tetrahedra. 
        # Maybe we can overlap nicely the elements_by_vertices_writers and 
        # the local_meshers for this case.
        init += physical_vertices_writers[E[0]](points, in_file+".ver")
    filter_repeated_faces(in_file, filter_repeated_vertices(in_file))
    mesh_graphics.plot_lines(in_file + ".ver", in_file + ".ebv")
    return
