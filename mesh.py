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

# {0,..,7} ---> {columns of trans}
permutation1     = { 0 : 0, 1 : 4, 2 : 3, 3 : 7, 4 : 1, 5 : 5, 6 : 2, 7 : 6 }

#permutation2     = { 6 : 0, 7 : 1, 5 : 2, 0 : 3, 2 : 4, 1 : 5, 4 : 6, 3 : 7 }

std_macro_elems = np.array([[0,1,3,4],[2,3,1,6],[7,3,4,6],[5,4,1,6],[6,1,3,4]])

reflections     = np.array([[ 1, -1, -1,  1,  1, -1, -1, 1],
                            [-1, -1,  1,  1, -1, -1,  1, 1],
                            [-1, -1, -1, -1,  1,  1,  1, 1]])

def lambda1 (i, j, k, n, mu):
    return float(i)/n * (float(i+j+k)/n)**((1/float(mu))-1)

def lambda2 (i, j, k, n, mu):
    return float(j)/n * (float(i+j+k)/n)**((1/float(mu))-1)

def lambda3 (i, j, k, n, mu):
    return float(k)/n * (float(i+j+k)/n)**((1/float(mu))-1)

def macroel_tetrahedra (vertices, mu, n):
    """ 
        n == number of levels, i.e. n + 1 == number of nodes per edge of 
        the macro--element.

        OBS: the graduation is, by default, towards 'P0',
        so we have to put the corresponding permutation1 of vertices when 
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
                lambda_[1,i,j,k] = lambda1(i,j,k,n,mu)
                lambda_[2,i,j,k] = lambda2(i,j,k,n,mu)
                lambda_[3,i,j,k] = lambda3(i,j,k,n,mu)
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
                ## TODO sefuir aca continuo 

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
                coef = (lambda1 (i,j,0,n,mu), lambda2 (i,j,0,n,mu))
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
            lambda_1, lambda_2 = lambda1 (y,z,0,levels,mu), lambda2 (y,z,0,levels,mu)
            convex_sum         = lambda_1*(macroel_vertices[:,1] - macroel_vertices[:,0]) \
                                 + lambda_2*(macroel_vertices[:,2] - macroel_vertices[:,0])
            top_layer         += [convex_sum + macroel_vertices[:,0]]
    top_layer   = np.vstack(tuple(top_layer))
    points      = np.vstack(tuple([top_layer+(x/n_vertical)*(macroel_vertices[:,3] \
                    - macroel_vertices[:,0]).T for x in range(n_vertical+1)]))
    return points

def split_cube_into_tetrahedra (nodes):
    """ nodes: the array of eight vertices of an hexahedron which
    has faces parallel to the axes. The order of the column of this input
    array is arbitrary, as long as the first one remains _the singular vertex_
    (if any singular vertex is present in this part of the mesh), and the 
    program performs the local graduation towards that vertex. """
    
    #singular_p = nodes[:,0]
    
    # compare with max{x}, max{y}, max{z} to decide
    # where is the singular vertex pointing




    bit_index = np.zeros(8,dtype=int)
    for i in range(8): # TODO transform to one line
        bits = ( nodes[0,i]==np.max(nodes[0,:]),\
                 nodes[1,i]==np.max(nodes[1,:]),\
                 nodes[2,i]==np.max(nodes[2,:]) ) 
        bit_index[i] = int(np.right_shift(np.packbits(bits),5))
    singular_vertex_orientation = permutation1[bit_index[0]]
    nodes = nodes*reflections[:,singular_vertex_orientation].reshape((3,1))

    # CONTINUE HERE: make one beatiful permutation,
    # once we have the hexahedron translated to the octant number 8
    # ONLY THEN we will be able to use the std_macro_elems

    ## REMOVE these 5 lines after tests
    # t0 = nodes[:,std_macro_elems[0]]  
    # t1 = nodes[:,std_macro_elems[1]]  
    # t2 = nodes[:,std_macro_elems[2]]  
    # t3 = nodes[:,std_macro_elems[3]]  
    # t4 = nodes[:,std_macro_elems[4]] 
    # 0. hybrid --never graded--; 1. hybrid; 2. hybrid; 3. hybrid; 4. tetra
    tetrahedra = [nodes[:,std_macro_elems[t]] for t in range(5)]
    return tetrahedra

def macroel_hybridhexa (vertices, mu, levels):
    """ This builds and returns the points in an hexahedron divided into five
    tetrahedral macroelements by calling the proper macroelement functions.
    Four of these tetrahedral macroelements are hybrid and include pyramids and
    prisms, while the fifth, 'inner', macroelement is built refined 
    tetrahedra. If it is a cube, then the 'inner' is a regular tetrahedron. 
    
    vertices: a 3 x 8 matrix with the vertices of the hexahedron. The first
    of these is the singular vertex (if any singular vertex is present in
    this part of the mesh), and we perform the graduation towards that vertex. 

    mu:       the graduation parameter

    levels:   levels of refinement, that is, the number of subintervals in any
    edge of the whole macroelement.
    """
    split_verts = split_cube_into_tetrahedra(vertices)
    print(split_verts)
    # with open ('split.txt','w') as out:
    #     out.write(split_verts)
    points      = { i : macroel_hybrid(split_verts[i], mu, levels) for i in range(4) }
    points[4]   = macroel_tetrahedra(split_verts[4], mu, levels)
    return points

def macroel_more_flexible_tending_more_generality ():
    ## TODO
    pass