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

from typing import List, Dict, Tuple


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

def kill_repeated(physical_vertices: str) -> Dict:
    """ 
        searches:       repeated vertices
                        It is merely an iterative cuadratic exploration.

        return value:   dictionary of the vertices that have to be re-numbered in file
                        elements_by_vertices.txt

        obs:            now we don't delete the entries on vertices.txt. We simply
                        don't request them.
    """ 

    counter: int       = 1
    vertices: np.array = np.loadtxt(physical_vertices)
    d_out: Dict        = {}
    Nv: int            = vertices.shape[0]

    for v in range(Nv):
        counter += 1
        for w in range(v + 1, Nv):
            if np.all(np.equal(vertices[v], vertices[w])):
                d_out[v] = d_out.get(v, [])
                d_out[v].append(w)

    return d_out


def split_preordered_hexahedron_into_tetrahedra(eight_hexahedron_vertices: List[float]) -> List:
    """ Here we pick subsequences of the vertices of the hexahedron to determine,
        one by one, the five tetrahedra.
    """
    number_of_r3_points = 8
    nodes               = np.array(eight_hexahedron_vertices[1:-1]).reshape(3, number_of_r3_points, order='F')
    tetrahedra          = []

    for t in range(5):
        # rows (4 x 3)
        tetrahedra += [{1: np.array([nodes[:, std_tetrahedra_inside_a_cube[t][j]] for j in range(4)])}]

    return tetrahedra

def elements_by_vertices_tetra(elements: List, init: int) -> None:
    """ TODO: we should write an algorithm that passes just one time per vertex
        directly from the four 'corner tetrahedra'
    """ 
    n_vert_repeated = 4
    arr_out = np.array(range(init + 1, init + n_vert_repeated + 1)).reshape((1, 4))

    elements.append(arr_out)
    

def elements_by_vertices_hybrid(corners: List, initial: int) -> None:
    """ 
    initial: tracks the number of vertices of the previous element.

    corners is passed in by reference, we mutate it with the calls.

    Writes GLOBAL INDICES per element appending in the list corners.
    """
    
    def current_initial(k: int) -> int:
        return 4 * k + 1
        
    line: List = [0, 0, 0, 0]

    for i in range(4):  # OJO QUE AHORA no escribo mas el 4.
        line[i] = initial + const_calc[i] + 1
        #                                   ^ OCTAVE_LANG_INI: int = 1
    corners.append(line)

def filter_repeated_vertices(elem_vert_repeated: np.array, in_file: str) -> int:

    CONTINUE HERE:

                1. Revisar todo el flujo y eliminar las escrituras intermedias de archivos y poner pasajes
                   de contenedores.

                    1.1. eliminar el "4" y el "1" en las puntas de cada fila de los bricks
                
                        1.1.1. controlar si comienzo indices en 0 o 1, de acuerdo a Liu

                1.2. el uso de kill_repeated() puede ser ciclando desde el min index de sorted_surrounding_nodes
                     o bien trayendo la indexation de inner_mesh

                2. correr tests de esta escritura para ver si no pierdo nada.
            
                3. Hacer el ejemplo para enviar
                
                4. resumir e integrar las implementaciones de 
                        *elements_by_vertices_hybrid

                        *elements_by_vertices_tetra
                        *macroel_tetra
                        *vertices_macro_tetra
                    y lo mismo para las referentes a tetra



    n_elem: int = len(elem_vert_repeated)  # number of elements with repetition of vertices indices.

    replace_verts: Dict = kill_repeated(in_file + ".ver")
    
    counter: int = 1
    
    # elem_vert_dscnt_indcs = np.copy(elem_vert_repeated[:, 1:]).reshape(n_elem * 6)

    elem_vert_repeated = elem_vert_repeated.reshape(4 * n_elem)

    for key in replace_verts:
        counter += 1
        for r in replace_verts[key]:
            # elem_vert_dscnt_indcs[elem_vert_dscnt_indcs == r + 1] = (key + 1)
            elem_vert_repeated[elem_vert_repeated == r + 1] = (key + 1)

    # whithout repetitions
    np.savetxt("wafer_profile_elements.dat", elem_vert_repeated.reshape((n_elem, 4)), fmt="%d")

    return n_elem



def macroel_hybrid(macroel_vertices):
    """ vertices = ( local_origin | base_vrtx_1 | base_vrtx_2 | sing_vrtx )
        TODO: el input de macroel_hybrid, macroel_vertices, deberian ser exactamente los points a retornar! :D """
    points = np.zeros((2, 2, 3, 2))  # level, i, coordinate, j
    for k in range(2):
        for i in range(2 - k):
            for j in range(2 - k - i):

                points[k,i,:,j] = i * (macroel_vertices[:, 1] - macroel_vertices[:, 0]) \
                                  + j * (macroel_vertices[:, 2] - macroel_vertices[:, 0])

                points[k,i,:,j] += k * (macroel_vertices[:, 3] - macroel_vertices[:, 0]) + macroel_vertices[:, 0]
    # This always returns (V0, V2, V1, V3)
    return points


def vertices_macro_hybrid(points, f_write):
    """ The file being constructed has the repeated physical vertices. This is just an enumeration of the points to write """
    with open (f_write, 'ab') as out:
        for l in range(2):
            for i in range(2 - l):
                for j in range(2 - l - i):
                    np.savetxt(out, points[l, i, :, j].reshape((1, 3)), fmt = "%8.5f")
    return 4


def macroel_tetrahedra(vertices) -> np.array:
    # TODO: don't let this as function
    return vertices.T


def vertices_macro_tetra(points, f_write):
    with open(f_write, 'ab') as out:
        np.savetxt(out, points, fmt = "%8.5f")

    return len(points)

physical_vertices_repeated: List        = []
all_elements_by_vertices_repeated: List = []
const_calc: Tuple                       = (0, 2, 1, 3)

initial_partition: np.array = np.loadtxt("bricks.txt")

partitioned_bricks: List = []

for brick in initial_partition:
    split: List = split_preordered_hexahedron_into_tetrahedra(brick)
    partitioned_bricks.append(split)

# Next loop is for building individual tetrahedra (with repeated nodes)
init: int = 0  # tracks the number of vertices
               # init should be number_of_vertices_inner_mesh
for s in partitioned_bricks:
    # Each s is a list with five dictionaries
    
    for t in range(4):
        # corner tetrahedra
        elements_by_vertices_hybrid(all_elements_by_vertices_repeated, init)
        v = s[t][1]
        for j in const_calc:
            # vertices_macro_hybrid(macroel_hybrid(v)) === hstack((v[0],v[2],v[1],v[3]))
            physical_vertices_repeated.append(v[j,:])
        init += 4
        
    # inner tetrahedron
    elements_by_vertices_tetra(all_elements_by_vertices_repeated, init) 
    v = s[4][1]
    for j in range(4):
        physical_vertices_repeated.append(v[j,:])
    init += 4


    """
    elements_by_vertices_hybrid(all_elements_by_vertices_repeated, init)
    init += vertices_macro_hybrid(macroel_hybrid(s[1][1]), "wafer_profile_vertices.txt")


    
    elements_by_vertices_hybrid(all_elements_by_vertices_repeated, init)
    init += vertices_macro_hybrid(macroel_hybrid(s[2][1]), "wafer_profile_vertices.txt")

    
    
    elements_by_vertices_hybrid(all_elements_by_vertices_repeated, init)
    init += vertices_macro_hybrid(macroel_hybrid(s[3][1]), "wafer_profile_vertices.txt")
    """    

# Here we filter repetitions
number_of_elements: int = filter_repeated_vertices(all_elements_by_vertices_repeated)
