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

__format__ = '%8.5f'

def write_element_indices (file_name, levels):
    """ Depends only on n. Not on which macroelement or octant is. 

    takes n and writes on disc:

    file_name.elem

    ---------------------------------------------

    6|(l,i,j)_{prism}|(l,i,j)_{prism}|(l,i,j)_{prism}|(l,i,j)_{prism}|(l,i,j)_{prism}|(l,i,j)_{prism}|

    4|(l,i,j)_{tetra}|(l,i,j)_{tetra}|(l,i,j)_{tetra}|(l,i,j)_{tetra}|

    5|(l,i,j)_{pyra} |(l,i,j)_{pyra} |(l,i,j)_{pyra} |(l,i,j)_{pyra} |(l,i,j)_{pyra} |

    etc...

    This funtion is initial in the program, for the case
    the mesh contains a hybrid macroelement. For a general mesh, the algorithm
    starts directly in the next step (with element_by_vertices.txt given 
    somehow).

    """
    string 		= ''
    n 			= levels
    for level in range(n-1):
        for i in range(n-level-1): # prisms
            for j in range(n-level-i-1):  
                nodos = [6,level,i,j,level,i+1,j,level,i,j+1,level+1,i,j,
                        level+1,i+1,j,level+1,i,j+1]
                string += ' '.join([str(r) for r in nodos]) +'\n'
        for i in range(n-level-2):
            for j in range(n-level-i-2):
                nodos = [6,level,i+1,j,level,i+1,j+1,level,i,j+1,
                        level+1,i+1,j,level+1,i+1,j+1,level+1,i,j+1]
                string += ' '.join([str(r) for r in nodos]) +'\n'
        for i in range(n - level): # tetra
            nodos = [4,
                    level, i, n-level-i-1,
                    level, i+1, n-level-i-1,
                    level, i, n-level-i,
                    level+1, i, n-level-i-1]
            string += ' '.join([str(r) for r in nodos]) + '\n'
        for i in range(n - level - 1): # pyramids: invariant: level+i+j == n-1
            nodos = [5, level,   i,   n-i-1-level,
                        level,   i+1, n-i-2-level,
                        level+1, i+1, n-i-2-level,
                        level+1, i,   n-i-1-level,
                        level,   i+1, n-i-1-level]
            string += ' '.join([str(r) for r in nodos]) + '\n'
    nodos = [4,n-1,0,0,n-1,1,0,n-1,0,1,n,0,0] # singular tetra 
    string += ' '.join([str(r) for r in nodos])
    with open(file_name, 'w') as elements:
        elements.write(string)
    return

def vertices_macro_hybrid (points, f_write):
    """ Writes the coordinates of the vertices of the local mesh in f_write.
        The file being constructed has the repeated physical vertices.
        Return value: number of vertices in this macroelement. """
    nVertices = 0
    L = points.shape[0]
    with open (f_write,'ab') as out:
        for l in range(L):
            for i in range(L-l):
                for j in range(L-l-i):
                    np.savetxt(out, points[l,i,:,j].reshape((1,3)), fmt = __format__)
                    nVertices += 1
        #out.write('\n')
    return nVertices

def vertices_macro_tetra (points, f_write):
    """ points: np.array with the points of the 
    tetrahedral non--hybrid macro--element """
    with open (f_write, 'ab') as out:
        np.savetxt(out, points, fmt = __format__)
    return len(points)

def vertices_macro_prism (points, f_write):
    """ points: np.array with the points of the 
    nested--prismatic macro--element """
    with open (f_write, 'ab') as out:
        np.savetxt(out, points, fmt = __format__)
    return len(points)

def vertices_macro_hybridhexa (points, f_write):
    """ points is a dictionary of length five of point arrays. 
    points == { 0 : p_0, ...} """
    return sum ([vertices_macro_hybrid (points[i], f_write) for i in range(4)])\
     + vertices_macro_tetra (points[4], f_write)
