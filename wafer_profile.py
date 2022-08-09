"""
To run:

In [i]: run wafer_profile.py n file_in


----------------------------------------------------------------------------------------------------
Ex. of the previous experiments run:

initial file, one line per node layer

#1              0,  -1,-.85,  0,   0,  -1,    0,   0,   0,  0,-1,  -.15

#2              0,-1.3,-.67,  0,  -1,-.85,    0,  -1,-.15,  0,-1.3,-.33

#3              0,-1.4,-.53,  0,-1.3,-.67,    0,-1.3,-.33,  0,-1.4,-.47
----------------------------------------------------------------------------------------------------


Para la primera version, inputs requeridos:
==========================================

Para toda layer: #1 z_upper,
                 #2 z_lower, 
                 #3 delta_x (radial).

Pasos:
=====

#1 agrupar de a cuatro sobre la 1ra layer

    Ejemplo para que quede en sentido antihorario:

    (Obs. LOS OVERLAPS ya fueron resueltos con el (j + 1) % n :))

               [-1.43, -1.36, -1.05, -0.61,    0,  0.61,  1.05,  1.36,  1.44, 1.5, 1.44, 1.36, -1.5]
               [-0.43, -0.62, -1.06, -1.36, -1.5, -1.36, -1.06, -0.62, -0.40,   0, 0.40, 0.62,    0]
               [  0.5,   0.5,   0.5,   0.5,  0.5,   0.5,   0.5,   0.5,   0.5, 0.5,  0.5,  0.5,  0.5]
                
                     ^      ^      ^
                     |      |      |
                                   |            
First square  >      1      4      |
                                   | 
                            ^      |
                            |      |
            
Second square >             1      4


                           ...

                     ^                                                                            ^
Last square   >      4                                                                            1


               [-1.43, -1.36, -1.05, -0.61,  0. ,  0.61,  1.05,  1.36,  1.44, 1.5, 1.44, 1.36, -1.5]
               [-0.43, -0.62, -1.06, -1.36, -1.5, -1.36, -1.06, -0.62, -0.4 ,  0., 0.4 , 0.62,  0. ]
               [ 0.  ,  0.  ,  0.  ,  0.  ,  0. ,  0.  ,  0.  ,  0.  ,  0.  ,  0., 0.  , 0.  ,  0. ]
                
                     ^      ^      ^
                     |      |      |
                                   | 
First square  >      2      3      |
                                   |
                            ^      |
                            |      |
                
Second square >             2      3


                           ...

                     ^                                                                            ^ 
Last square   >      3                                                                            2      



#2 Para cada lista de 4, hacer los cuatro de adelante según z_up, z_lo, dx
   (así, enumerativo, porque LIU pidió simpleza)

   Hay que tomar la direcci'on radial, y en Point.h ya hay una norm() y operator/ para esto :)

#3 Para cada uno, hacemos "split_PREORDERED_hexahedron_into_tetrahedra"

TODO:
====

Make it possible to have this as an input:

                                            1, 1, 1
                                            1, 0, 1
                                            1, 0, 0
                                            1, 1, 0
                                            0, 1, 1
                                            0, 0, 1
                                            0, 0, 0
                                            0, 1, 0

Perhaps some small dictionary to select the parsers, for example:

                    input: List = [{
                                      "type": 4,
                                      "grading": 1,
                                      "outer_vertices": [1, 1, 1,
                                                         1, 0, 1,
                                                         1, 0, 0,
                                                         1, 1, 0,
                                                         0, 1, 1,
                                                         0, 0, 1,
                                                         0, 0, 0,
                                                         0, 1, 0]
                                  }]

wafer cases:
                "edge_a"
                "edge_bc"
                "orientation_flat_a"
                "orientation_flat_b"
                "notch"

"""

import sys
import numpy as np

from typing import List, Dict


def wafer_profile_generated_by_regular_polygon(n: int, o: np.array) -> None:
    """
    o is some "istream" with just four points in R3.
    
    This is just the use of the case of a regular polygon: 
    TODO: fix this to be the last part of this method, to use it with just one call.
    =====================================================

    o               = np.loadtxt(file_in, delimiter=",")

    flip = -1
    paired_squares_of_points: Dict = None  # as we see the eight nodes of a hexahedron

    out_hexahedra = np.vstack([np.vstack([
        make_brick(paired_squares_of_points[l][j % n], 
                   paired_squares_of_points[l][(j + 1) % n], 
                   flip**(j + l + 1))
        for j in range(n)
    ]) for l in range(surround_layers)])

    np.savetxt(file_out, out_hexahedra, delimiter=",", fmt="%5.2f")
    """

    def rotation(alpha: float) -> np.array:
        """
        This returns the matrix for the rotation of alpha radians around the Z axis.
        """
        return np.array([[np.cos(alpha), -np.sin(alpha), 0],
                         [np.sin(alpha), np.cos(alpha), 0],
                         [0, 0, 1]])
    
    def make_brick(left_points: np.array, right_points: np.array, flip: int) -> np.array:
        if flip > 0:
            left_points, right_points = right_points, left_points

        right_indices           = [4, 5, 6, 7]
        left_indices            = [0, 1, 2, 3]
        brick                   = np.zeros((3, 8))
        brick[:, right_indices] = right_points
        brick[:, left_indices]  = left_points
        return np.hstack((4, brick.flatten(order="F"), 1)).reshape((1, 26))

    p               = o.reshape((surround_layers, 3, 4), order="F") \
                        + 3 * np.array([0, -1, 0]).reshape((3, 1))
    alpha           = 2 * np.pi/n
    rotation_matrix = rotation(alpha)
    rotated_points  = {l: {0: p[l].copy()} for l in range(surround_layers)}

    for l in range(surround_layers):
        j = 1
        while j < n:
            p[l] = np.dot(rotation_matrix, p[l])
            rotated_points[l][j] = p[l].copy()
            j = j + 1

################################################################################################
################################################################################################
################################################################################################
################################################################################################

"""
surround_layers = 3  # sys.argv[*]
wafer_case      = 1  #sys.argv[1]
file_in         = sys.argv[2] 
file_out        = "bricks_for_" + file_in
"""

std_unitary_cube = np.array([[1, 1, 1, 1, 0, 0, 0, 0],
                             [1, 0, 0, 1, 1, 0, 0, 1],
                             [1, 1, 0, 0, 1, 1, 0, 0]])

std_tetrahedra_inside_a_cube = np.array([[0, 1, 3, 4], 
                                         [2, 3, 1, 6], 
                                         [7, 3, 4, 6], 
                                         [5, 4, 1, 6], 
                                         [6, 1, 3, 4]])

def layered_mesh(
        sorted_surrounding_nodes: str, 
        upper_z: float = .4,
        lower_z: float = .1,
        delta_over_plane_xy: float = .8
    ) -> None:

    """
    Documenting example:

    Take the following

    >>> points_in_first_layer = np.array([
                                    [-1.43645, -0.431997, .5],
                                    [-1.43645, -0.431997, 0],
                                    [-1.36396, -0.624191, 0],
                                    [-1.36396, -0.624191, .5]
                                ])

    and then do the following
    
    >>> p[[0, 3], 2] = upper_z                                                                      
    >>> p[[1, 2], 2] = lower_z   
    >>> p[0, 0:2]   += (delta_over_plane_xy/np.linalg.norm(p[0, 0:2]))*p[0, 0:2]                      
    >>> p[1, 0:2]   += (delta_over_plane_xy/np.linalg.norm(p[1, 0:2]))*p[1, 0:2]                      
    >>> p[2, 0:2]   += (delta_over_plane_xy/np.linalg.norm(p[2, 0:2]))*p[2, 0:2]                      
    >>> p[3, 0:2]   += (delta_over_plane_xy/np.linalg.norm(p[3, 0:2]))*p[3, 0:2]                      

    and the points in the second layer should be

    >>> points_in_second_layer = np.array([
                                    [-2.20255496, -0.66239489, 0.4],
                                    [-2.20255496, -0.66239489, 0.1],
                                    [-2.09140513, -0.95709277, 0.1],
                                    [-2.09140513, -0.95709277, 0.4]
                                ])
    """
    
    sorted_surrounding_nodes: np.array = np.loadtxt(sorted_surrounding_nodes, delimiter=",")

    quantity_of_squares_per_layer: int = sorted_surrounding_nodes.shape[0]//2

    def construct_second_layer(fst: np.array) -> np.array:
        second_layer: np.array  = fst.copy()
        second_layer[[0, 3], 2] = upper_z                                                                      
        second_layer[[1, 2], 2] = lower_z   
        second_layer[0, 0:2]   += (delta_over_plane_xy/np.linalg.norm(second_layer[0, 0:2])) \
                                    * second_layer[0, 0:2]
        second_layer[1, 0:2]   += (delta_over_plane_xy/np.linalg.norm(second_layer[1, 0:2])) \
                                    * second_layer[1, 0:2]
        second_layer[2, 0:2]   += (delta_over_plane_xy/np.linalg.norm(second_layer[2, 0:2])) \
                                    * second_layer[2, 0:2]
        second_layer[3, 0:2]   += (delta_over_plane_xy/np.linalg.norm(second_layer[3, 0:2])) \
                                    * second_layer[3, 0:2] 
        return second_layer

    stack: List = []
    first_layer: np.array
    second_layer: np.array
    this_brick: np.array

    for s in range(quantity_of_squares_per_layer - 1):
        # TODO: still the last square missing, wrapping around.
        first_layer = sorted_surrounding_nodes[[
                                    s + 0, 
                                    s + quantity_of_squares_per_layer,
                                    s + quantity_of_squares_per_layer + 1,
                                    s + 1
                                ]]

        # The flip is just the horizontal stacking here:
        second_layer = construct_second_layer(first_layer) 
        this_brick = np.hstack((4, first_layer.flatten(), second_layer.flatten(), 1))
        stack.append(this_brick)

    del this_brick

    # The last brick, which uses the first nodes again.
    first_layer = sorted_surrounding_nodes[[
                      quantity_of_squares_per_layer - 1,
                      -1,
                      quantity_of_squares_per_layer,
                      0
                  ]]

    second_layer = construct_second_layer(first_layer)
    last_brick: np.array = np.hstack((4, first_layer.flatten(), second_layer.flatten(), 1))
    stack.append(last_brick)

    np.savetxt("bricks.txt", np.array(stack), delimiter=",", fmt="%5.2f")

"""
Now the classic mesher enters the game.
"""

CONTINUE HERE:
                1. eliminar las escrituras intermedias de archivos y poner pasajes
                   de contenedores.
                
                2. correr tests de esta escritura para ver si no pierdo nada.
                    2.1 resumir e integrar las implementaciones de 
                            *elements_by_vertices_hybrid
                            *macroel_hybrid
                            *vertices_macro_hybrid
                        y lo mismo para las referentes a tetra

                3. eliminar el "4" y el "1" en las puntas de cada fila de los bricks


def kill_repeated(vertices_file_name: str) -> Dict:
    """ 
        searches:       repeated vertices
                        It is merely an iterative cuadratic exploration.

        return value:   dictionary of the vertices that have to be re-numbered in file
                        elements_by_vertices.txt

        obs:            now we don't delete the entries on vertices.txt. We simply
                        don't request them.
    """ 

    counter: int       = 1
    vertices: np.array = np.loadtxt(vertices_file_name)
    d_out: Dict        = {}
    Nv: int            = vertices.shape[0]

    for v in range(Nv):
        counter += 1
        for w in range(v + 1, Nv):
            if np.all(np.equal(vertices[v], vertices[w])):
                d_out[v] = d_out.get(v, [])
                d_out[v].append(w)

    return d_out

def filter_repeated_vertices(in_file: str, max_number_of_vertices: int = 6) -> int:

    with open(in_file + ".ebvr", 'r') as inp:
        things: List = inp.readlines()
    
    n_elem: int = len(things)  # number of elements with repetition of vertices.
    
    elem_vert_repeated: np.array = np.zeros((n_elem, max_number_of_vertices + 1), dtype=int)

    for k in range(n_elem):
        ele: np.array            = np.fromstring(things[k], dtype=int, sep=" ")
        trailing_zeros: np.array = np.zeros(max_number_of_vertices + 1 - len(ele), dtype=int)
        elem_vert_repeated[k]    = np.concatenate((ele, trailing_zeros))
        # TODO:                    ^ change this to elem_vert_repeated[k, 0:len(ele)],
        #                            eliminate each trailing_zeros array and test :)

    replace_verts = kill_repeated(in_file + ".ver")
    
    counter: int = 1
    elem_vert_dscnt_indcs = np.copy(elem_vert_repeated[:, 1:]).reshape(n_elem * 6)
    for key in replace_verts:

        counter += 1
        for r in replace_verts[key]:
            # Make_by_hand
            elem_vert_dscnt_indcs[elem_vert_dscnt_indcs == r + 1] = (key + 1)

    elem_vert_dscnt_indcs = elem_vert_dscnt_indcs.reshape((n_elem,6))
    col = elem_vert_repeated[:, 0].reshape(elem_vert_repeated.shape[0], 1)
    elem_vert_dscnt_indcs = np.hstack((col, elem_vert_dscnt_indcs))

    ## <--- whithout repetitions
    # Make_by_hand: check if it is correct to np.savetxt something for each e
    with open(in_file + ".ebv", "wb") as out:
        for e in elem_vert_dscnt_indcs:
            np.savetxt(out, e[0:e[0] + 1].reshape((1, e[0] + 1)), fmt="%d")

    return n_elem


initial_partition: np.array = np.loadtxt("bricks.txt", levels=1)

tetrahedra_inside_bricks: List = []

for brick in initial_partition:
    split: List = mesh.split_preordered_hexahedron_into_tetrahedra(brick)
    tetrahedra_inside_bricks.append(split)

init: int = 0
for s in tetrahedra_inside_bricks:
    mesh_connectivity.elements_by_vertices_hybrid("wafer_profile", 1, "Octave", init) # type 0
    points = mesh.macroel_hybrid(s[0][1], 1, 1)
    init += mesh_write.vertices_macro_hybrid(points, "wafer_profile_vertices.txt")

    mesh_connectivity.elements_by_vertices_hybrid("wafer_profile", 1, "Octave", init) # type 0
    points = mesh.macroel_hybrid(s[1][1], 1, 1)
    init += mesh_write.vertices_macro_hybrid(points, "wafer_profile_vertices.txt")

    mesh_connectivity.elements_by_vertices_hybrid("wafer_profile", 1, "Octave", init) # type 0
    points = mesh.macroel_hybrid(s[2][1], 1, 1)
    init += mesh_write.vertices_macro_hybrid(points, "wafer_profile_vertices.txt")
    
    mesh_connectivity.elements_by_vertices_hybrid("wafer_profile", 1, "Octave", init) # type 0
    points = mesh.macroel_hybrid(s[3][1], 1, 1)
    init += mesh_write.vertices_macro_hybrid(points, "wafer_profile_vertices.txt")
    
    mesh_connectivity.elements_by_vertices_tetra("wafer_profile", 1, "Octave", init)  # type 1
    points = mesh.macroel_tetrahedra(s[4][1], 1, 1)
    init += mesh_write.vertices_macro_tetra(points, "wafer_profile_vertices.txt")

number_of_elements: int = filter_repeated_vertices("wafer_profile")