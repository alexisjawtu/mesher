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

from typing import List


# {0,..,7} ---> {columns of reflections}
# Meaning: we take the integer whos binary representation
# equals the R3 coordinates of the singular vertex
# viewed in the unitary cube and returns the octant number
# which has the origin in that orientation:
# for example: if the input I_0 is
# (1,1,1)  (1,2,1)  (0,1,1)  (0,2,1)  
# (1,1,2)  (1,2,2)  (0,1,2)  (0,2,2)
# with singular vertex (1,2,2), then (1,2,2) is oriented
# as the vertex (1,1,1) in the unitary cube, so it is number 7.
# Now, octant_encoding[7] = 6,
# then at the end we have: I_0 ---> octant number 6.

# group = np.array([[6,5,7,4,2,1,3,0],  ## these are the succesive output arrays in
#                   [2,1,3,0,6,5,7,4],  ## reference_permutations()
#                   [3,0,2,1,7,4,6,5],  ## we leave them just because
#                   [7,4,6,5,3,0,2,1],
#                   [5,6,4,7,1,2,0,3],
#                   [1,2,0,3,5,6,4,7],
#                   [0,3,1,2,4,7,5,6], 
#                   [4,7,5,6,0,3,1,2]])

#    """ This builds and returns the points in an hexahedron divided into five
#    tetrahedral macroelements by calling the proper macroelement functions.
#    Four of these tetrahedral macroelements are hybrid and include pyramids and
#    prisms, while the fifth, 'inner', macroelement is built refined 
#    tetrahedra. If it is a cube, then the 'inner' is a regular tetrahedron. 
#    
#    vertices: a 3 x 8 matrix with the vertices of the hexahedron. The first
#    of these is the singular vertex (if any singular vertex is present in
#    this part of the mesh), and we perform the graduation towards that vertex. 
#
#    mu:       the graduation parameter
#
#    levels:   levels of refinement, that is, the number of subintervals in any
#    edge of the whole macroelement.
#
#    positions_encoding ->  
#    dato en lugar                  [7 3 4 5 1 6 0 2] 
#                                    | | | | | | | |
#    "en base 2 vale"                | | | | | | | |   
#                                    | | | | | | | | 
#                                    v v v v v v v v
#                                    0 1 2 3 4 5 6 7 
#                                    | | | | | | | |   
#                                    v v v v v v v v
#    que en perm. ref oc=2 es:      [3,0,2,1,7,4,6,5] == group[2] == reference_permutations (2)
#    """

octant_encoding = { 0 : 0, 1 : 4, 2 : 3, 3 : 7, 4 : 1, 5 : 5, 6 : 2, 7 : 6 }
p_              = np.array([[ 1,  1,  1,  1,  0,  0,  0,  0],   
                            [-1,  0,  0, -1, -1,  0,  0, -1],   
                            [-1, -1,  0,  0, -1, -1,  0,  0]])
# 0. hybrid --never graded--;
# 1. hybrid;
# 2. hybrid;
# 3. hybrid;
# 4. tetra;
std_macro_elems = np.array([[0,1,3,4],[2,3,1,6],[7,3,4,6],[5,4,1,6],[6,1,3,4]])
reflections     = np.array([[ 1, -1, -1,  1,  1, -1, -1, 1],
                            [-1, -1,  1,  1, -1, -1,  1, 1],
                            [-1, -1, -1, -1,  1,  1,  1, 1]])

def bit_arrays (nodes):
    """
    This is to play with permutations to determine the orientation of the brick with respect
    to some singular vertex.

    Warning: this packbits only work for bricks with faces parallel to the coordinate planes!
    """

    bit_arr = np.zeros(8,dtype=int)
    for i in range(8): # TODO transform to one line
        bits = ( nodes[0,i]==np.max(nodes[0,:]),\
                 nodes[1,i]==np.max(nodes[1,:]),\
                 nodes[2,i]==np.max(nodes[2,:]) )
        bit_arr[int(np.right_shift(np.packbits(bits),5))] = i
    return bit_arr 

def reference_permutations(orientation):
    """ the vertices permutations: this are the transformed
    enumerations of the vertices, one for each singular_vertex_position,
    so that we keep the numbers in std_macro_elems. Comparisons with 
    the maximum x, y and z coordinates and the transforming from base 2 
    to base 10. """
    nodes       = p_*reflections[:,orientation].reshape((3,1))
    permutation = [0]*8
    for i in range(8):
        permutation[bit_arrays(nodes)[i]] = i
    return permutation

def convex_coef (i, ijk, n, mu):
    """ ijk is the list [i,j,k] in the standard notation for this
    grading technique 
    TODO: in python 3 I can remove the float() calls
    """
    return ijk[i-1]/n * (sum(ijk)/n)**(1/mu-1)
    
def macroel_tetrahedra (vertices, mu, n):
    """ 
        n == number of levels, i.e. n + 1 == number of nodes per edge of 
        the macro--element.

        OBS: the graduation is, by default, towards 'P0',
        so we have to put the corresponding octant_encoding of vertices when 
        calling the function.
        
        TODO:  calculations of the form lambda_[0,i,j,k]*P0 could be done
        a priori and then call in the following manner:

            initialize an (4 x I x J x K x 3)-D array Lambda_ with

                    Lambda_[l,i,j,k,:] <--- lambda_[l,i,j,k]*P_l

            and then call it. Make sure the dimension is the most comfortable.
    """
    vertices = vertices.transpose()

    P0 = vertices[0,:].reshape(1,3)
    P1 = vertices[1,:].reshape(1,3)
    P2 = vertices[2,:].reshape(1,3)
    P3 = vertices[3,:].reshape(1,3)

    lambda_ = np.zeros((4, n+1, n+1, n+1))

    points  = np.zeros((1,3))

    ## todo: fix the order of the indices

    for k in range(n+1):
        for j in range(n+1):
            for i in range(n+1):
                lambda_[1,i,j,k] = convex_coef(1, [i,j,k], n, mu)
                lambda_[2,i,j,k] = convex_coef(2, [i,j,k], n, mu)
                lambda_[3,i,j,k] = convex_coef(3, [i,j,k], n, mu)
                lambda_[0,i,j,k] = 1 - lambda_[1,i,j,k] \
                                   - lambda_[2,i,j,k] - lambda_[3,i,j,k]

    # now: element vertices
    for k in range(n):
        for j in range(n-k):
            for i in range(n-j-k):
                q0 = lambda_[0,i,j,k]  *P0 + lambda_[1,i,j,k]  *P1 \
                     + lambda_[2,i,j,k]*P2 + lambda_[3,i,j,k]*P3
                q1 = lambda_[0,i+1,j,k]*P0 + lambda_[1,i+1,j,k]*P1 \
                     + lambda_[2,i+1,j,k]*P2 + lambda_[3,i+1,j,k]*P3
                q2 = lambda_[0,i,j+1,k]*P0 + lambda_[1,i,j+1,k]*P1 \
                     + lambda_[2,i,j+1,k]*P2 + lambda_[3,i,j+1,k]*P3
                q3 = lambda_[0,i,j,k+1]*P0 + lambda_[1,i,j,k+1]*P1 \
                     + lambda_[2,i,j,k+1]*P2 + lambda_[3,i,j,k+1]*P3

                points = np.vstack((points,q0,q1,q2,q3))
                del(q0,q1,q2,q3)

    for k in range(n-1):
        for j in range(n-1-k):
            for i in range(n-1-j-k):
                q0 = lambda_[0,i+1,j,k]*P0 + lambda_[1,i+1,j,k]*P1 \
                     + lambda_[2,i+1,j,k]*P2 + lambda_[3,i+1,j,k]*P3
                q1 = lambda_[0,i,j+1,k]*P0 + lambda_[1,i,j+1,k]*P1 \
                     + lambda_[2,i,j+1,k]*P2 + lambda_[3,i,j+1,k]*P3
                q2 = lambda_[0,i,j,k+1]*P0 + lambda_[1,i,j,k+1]*P1 \
                     + lambda_[2,i,j,k+1]*P2 + lambda_[3,i,j,k+1]*P3
                q3 = lambda_[0,i+1,j,k+1]*P0 + lambda_[1,i+1,j,k+1]*P1 \
                     + lambda_[2,i+1,j,k+1]*P2 + lambda_[3,i+1,j,k+1]*P3

                points = np.vstack((points,q0,q1,q2,q3))

                q0 = lambda_[0,i,j+1,k]*P0 + lambda_[1,i,j+1,k]*P1 \
                     + lambda_[2,i,j+1,k]*P2 + lambda_[3,i,j+1,k]*P3
                q1 = lambda_[0,i,j,k+1]*P0 + lambda_[1,i,j,k+1]*P1 \
                     + lambda_[2,i,j,k+1]*P2 + lambda_[3,i,j,k+1]*P3
                q2 = lambda_[0,i+1,j,k+1]*P0 + lambda_[1,i+1,j,k+1]*P1 \
                     + lambda_[2,i+1,j,k+1]*P2 + lambda_[3,i+1,j,k+1]*P3
                q3 = lambda_[0,i,j+1,k+1]*P0 + lambda_[1,i,j+1,k+1]*P1 \
                     + lambda_[2,i,j+1,k+1]*P2 + lambda_[3,i,j+1,k+1]*P3

                points = np.vstack((points,q0,q1,q2,q3))
                
                q0 = lambda_[0,i+1,j,k]*P0 + lambda_[1,i+1,j,k]*P1 \
                     + lambda_[2,i+1,j,k]*P2 + lambda_[3,i+1,j,k]*P3
                q1 = lambda_[0,i,j+1,k]*P0 + lambda_[1,i,j+1,k]*P1 \
                     + lambda_[2,i,j+1,k]*P2 + lambda_[3,i,j+1,k]*P3
                q2 = lambda_[0,i+1,j+1,k]*P0 + lambda_[1,i+1,j+1,k]*P1 \
                     + lambda_[2,i+1,j+1,k]*P2 + lambda_[3,i+1,j+1,k]*P3
                q3 = lambda_[0,i+1,j,k+1]*P0 + lambda_[1,i+1,j,k+1]*P1 \
                     + lambda_[2,i+1,j,k+1]*P2 + lambda_[3,i+1,j,k+1]*P3

                points = np.vstack((points,q0,q1,q2,q3))
                
                q0 = lambda_[0,i,j+1,k]*P0 + lambda_[1,i,j+1,k]*P1 \
                     + lambda_[2,i,j+1,k]*P2 + lambda_[3,i,j+1,k]*P3
                q1 = lambda_[0,i+1,j+1,k]*P0 + lambda_[1,i+1,j+1,k]*P1 \
                     + lambda_[2,i+1,j+1,k]*P2 + lambda_[3,i+1,j+1,k]*P3
                q2 = lambda_[0,i+1,j,k+1]*P0 + lambda_[1,i+1,j,k+1]*P1 \
                     + lambda_[2,i+1,j,k+1]*P2 + lambda_[3,i+1,j,k+1]*P3
                q3 = lambda_[0,i,j+1,k+1]*P0 + lambda_[1,i,j+1,k+1]*P1 \
                     + lambda_[2,i,j+1,k+1]*P2 + lambda_[3,i,j+1,k+1]*P3
                
                points = np.vstack((points,q0,q1,q2,q3))

    for k in range(n-2):
        for j in range(n-2-k):
            for i in range(n-2-j-k):
                q0 = lambda_[0,i+1,j+1,k]*P0 + lambda_[1,i+1,j+1,k]*P1 \
                     + lambda_[2,i+1,j+1,k]*P2 + lambda_[3,i+1,j+1,k]*P3
                q1 = lambda_[0,i+1,j,k+1]*P0 + lambda_[1,i+1,j,k+1]*P1 \
                     + lambda_[2,i+1,j,k+1]*P2 + lambda_[3,i+1,j,k+1]*P3
                q2 = lambda_[0,i,j+1,k+1]*P0 + lambda_[1,i,j+1,k+1]*P1 \
                     + lambda_[2,i,j+1,k+1]*P2 + lambda_[3,i,j+1,k+1]*P3
                q3 = lambda_[0,i+1,j+1,k+1]*P0 + lambda_[1,i+1,j+1,k+1]*P1 \
                     + lambda_[2,i+1,j+1,k+1]*P2 + lambda_[3,i+1,j+1,k+1]*P3

                points = np.vstack((points,q0,q1,q2,q3))

    points = np.delete(points, 0, 0)
    return points

def macroel_hybrid (macroel_vertices, mu, n):
    """ 
        vertices = ( local_origin | base_vrtx_1 | base_vrtx_2 | sing_vrtx )
        singular edge the one which is parallel to v = (sing_vrtx - local_origin)
        n  == levels
        mu == grading param
    """
    points = np.zeros((n+1, n+1, 3, n+1))  # level, i, coordinate, j
    for k in range(n+1):
        for i in range(n-k+1):
            for j in range(n-k-i+1):
                # it can be done with much lesser calls to these lambda
                coef = (convex_coef (1, [i,j,0], n, mu), convex_coef (2, [i,j,0], n, mu))
                # the sub-mesh is just the following two lines
                points[k,i,:,j] = coef[0]*(macroel_vertices[:,1]-macroel_vertices[:,0]) \
                                  + coef[1]*(macroel_vertices[:,2]-macroel_vertices[:,0])
                points[k,i,:,j] += (1-(float(n-k)/n)**(1/mu))*(macroel_vertices[:,3]-macroel_vertices[:,0]) \
                                  + macroel_vertices[:,0]
    return points

def macroel_prisms (macroel_vertices, mu, levels):
    """ 
    levels := number of subintervals betweeen nodes of each horizontal edge of the
    macroelement
    n_vertical := number of subintervals betweeen nodes of each vertical edge of the
    macroelement. Warning: my Thesis states that n_vertical == levels. Otherwise it's not proved.

    TODO: In a future version we will have n_vertical independent of levels

    order of the columns in M := macroel_vertices:
        (M[:,0], M[:,1], M[:,2]) == a triangle
        (M[:,3], M[:,4], M[:,5]) == the other triangle
    
    edges perpendicular to the triangles:
        singular edge (if it is singular) == [M[0],M[3]] 
        the others:             [M[1],M[4]], [M[2],M[5]]

    OBS: to sum faster, in the array of points, the levels
    greater than zero end up filled as a rectangle of points. 
    Be careful not to use that coordinates. 

    At the end we translate level 0 to the levels above and make a stack
    with the "upper triangle" of the points array in each level.
    """
    n_vertical  = levels 
    top_layer   = []
    for y in range(levels+1):
        for z in range(levels+1-y):
            lambda_1, lambda_2 = convex_coef (1, [y,z,0], levels, mu), convex_coef (2, [y,z,0], levels, mu)
            convex_sum         = lambda_1*(macroel_vertices[:,1] - macroel_vertices[:,0]) \
                                 + lambda_2*(macroel_vertices[:,2] - macroel_vertices[:,0])
            top_layer         += [convex_sum + macroel_vertices[:,0]]
    top_layer   = np.vstack(tuple(top_layer))
    points      = np.vstack(tuple([top_layer+(x/n_vertical)*(macroel_vertices[:,3] \
                    - macroel_vertices[:,0]).T for x in range(n_vertical+1)]))
    return points

def split_preordered_hexahedron_into_tetrahedra(eight_hexahedron_vertices: List[float]) -> List:
    """
    Requieres the eight points in eight_hexahedron_vertices to be ordered according to the
    standard macro--elements defined in main.py. 
    If this ordering is correct, then the hexaedron doesn't need to have faces parallel to
    the coordinates.
    """
    nodes = np.array(eight_hexahedron_vertices[1:-1]).reshape(3,(len(eight_hexahedron_vertices)-2)//3,order='F')
    tetrahedra    = []
    types_in_cube = [0,0,0,0,1]
    grads_in_cube = [1] + [eight_hexahedron_vertices[-1]]*4
    
    for t in range (5):
        tetrahedra += [ { 0 : types_in_cube[t],
         1 : np.array([nodes[:,std_macro_elems[t][j]] for j in range(4)]).T,
         2 : grads_in_cube[t] } ]

    return tetrahedra


def split_brick_into_prisms (macroel_raw):
    pass

def split_lshape_into_prisms (macroel_raw):
    pass

def reorder_prism_vertices (macroel_raw):
    pass

def split_parallel_brick_into_tetrahedra (macroel_raw):
    """ nodes: the array of eight vertices of an hexahedron which
    has faces parallel to the axes. The columns of this input
    array may be in any order, as long as the first one remains as
    _the singular vertex_ (if any singular vertex is present in this
    part of the mesh), and the program performs the local graduation towards
    that vertex. Compare with max{x}, max{y}, max{z} to decide where is the
    singular vertex pointing to """
    nodes = np.array(macroel_raw[1:-1]).reshape(3,(len(macroel_raw)-2)//3,order='F')
    positions_encoding       = bit_arrays(nodes)
    singular_vertex_position = octant_encoding[np.where(positions_encoding==0)[0][0]]
    pi                       = reference_permutations (singular_vertex_position)
    tetrahedra    = []
    types_in_cube = [0,0,0,0,1]
    grads_in_cube = [1] + [macroel_raw[-1]]*4
    for t in range (5):
        tetrahedra += [ { 0 : types_in_cube[t],
         1 : np.array([nodes[:,positions_encoding[pi[std_macro_elems[t][j]]]]\
                       for j in range(4)]).T,
         2 :  grads_in_cube[t] } ]
    return tetrahedra

def macroel_more_flexible_tending_more_generality ():
    pass
