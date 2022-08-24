## <one line to give the program's name and a brief idea of what it does.>
##     Copyright (C) 2020  Ignacio Ojea.
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

"""
plot_all_tetrahedra("experiments/wafer30jul22/nodes.dat",
                    "experiments/wafer30jul22/elements.dat",
                    "experiments/wafer30jul22/bricks.txt.ver",
                    "experiments/wafer30jul22/bricks.txt.ebv",
                    vert_delim=",")
"""

import sys

import numpy as np

from typing import Tuple
from mayavi import mlab


def plot_lines(
        vertices_file: str, 
        connectivity_file: str, 
        isolated_points: str = None,
        vert_delim: str = None, 
        colors: Tuple = (.8,.8,.8)) -> None:

    # Reads .ver and .ebv files. col is a color definition. 
    p = np.loadtxt(vertices_file, delimiter=vert_delim)
    with open(connectivity_file,'r') as infile:
        inlist = infile.readlines()
    con_list = [line.strip(' \n').split(' ') for line in inlist]
    con_list = [[int(st[k]) for k in range(len(st))] for st in con_list]
    x = p[:,0]
    y = p[:,1]
    z = p[:,2]
    # s could be modified in order to give different colors to different elements. 
    s = np.ones(len(x), dtype="float64")
    cant_edges = {6:9, 5:8, 4:6}  # dictionary {n_nodes: n_edges}
    n_con = 0
    for i in range(len(con_list)):
        n_con = n_con + cant_edges[con_list[i][0]]
    connections = np.zeros((n_con,2))
    
    # Connections for each type of element
    connections_prism = np.array([
        [0,1,2,0,3,4,5,1,2], 
        [1,2,0,3,4,5,3,4,5]
    ]).T
    
    connections_tetra = np.array([
        [0,1,2,3,0,1],
        [1,2,3,0,2,3]
    ]).T
    
    connections_pyrad = np.array([
        [0,1,2,3,0,1,2,3],
        [1,2,3,0,4,4,4,4]
    ]).T

    # dictionary
    new_connections = {6: connections_prism, 5: connections_pyrad, 4: connections_tetra}
    last = 0
    for i in range(len(con_list)):
        row = np.array(con_list[i])
        #add connections depending on the type of elelement

        # TODO: this [[]+1]-1 at the end of the line should be depending of lang = C/Octave
        connections[last:last+cant_edges[row[0]],:] = row[new_connections[row[0]]+1]-1
        last = last + cant_edges[row[0]]
    #plot:
    fig = mlab.figure(1, size=(400, 400), bgcolor=(1, 1, 1))
    src = mlab.pipeline.scalar_scatter(x, y, z)
    src.mlab_source.dataset.lines = connections 
    src.update()
    mlab.pipeline.surface(src, color=colors)

    # Put some balls to visualize the axes
    mlab.points3d([0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1], scale_factor=.06, color=(.5, 1, .5))

    if isolated_points:
        isolated_points = np.loadtxt(isolated_points, delimiter=vert_delim)
        
        mlab.points3d(
            isolated_points[:,0],
            isolated_points[:,1],
            isolated_points[:,2], 
            np.linspace(.001, .3, len(isolated_points)),
            scale_factor=1, 
            color=(1, 0, 0)
        )

    mlab.show()


def plot_all_tetrahedra(
        vertices_file: str,
        connectivity_file: str,
        nodes: str = None,
        elements_file: str = None,
        isolated_points: str = None,
        vert_delim: str = None,
        nodes_delim: str = None,
        colors: Tuple = (.8,.8,.8)) -> None:


    # CONTINUE HERE:
    # ** 1- poner el brick faltante
    # 2- hacer una layer más
    # ** 3- enviar a Liu la propuesta de input con upper_z, delta_xy, etc.
    # 4- pedir una parte del dinero.

    pair: Tuple = (
        np.loadtxt(vertices_file, delimiter=vert_delim),
        np.loadtxt(nodes, delimiter=nodes_delim)
    )

    vertices: np.array = np.vstack(pair)

    x = vertices[:, 0]
    y = vertices[:, 1]
    z = vertices[:, 2]

    del vertices
    
    connectivity: np.array = np.loadtxt(connectivity_file)
    elements: np.array = np.loadtxt(elements_file)

    # Translate indices to have all in the same picture.
    elements[:, 1:] += len(pair[0])

    all_elements: np.array = np.vstack((connectivity, elements))

    # Here we have all tetrahedra, which have 6 edges.
    cant_edges = {4: 6}
    n_con = 0
    for e in all_elements:
        n_con = n_con + cant_edges[e[0]]

    connections = np.zeros((n_con,2))

    connections_tetra = np.array([
        [0,1,2,3,0,1],
        [1,2,3,0,2,3]
    ]).T

    new_connections = {4: connections_tetra}

    last = 0

    for i in range(len(all_elements)):
        row = np.array(all_elements[i])
        #add connections depending on the type of elelement

        # TODO: this [[]+1]-1 at the end of the line should be depending of lang = C/Octave
        connections[last: last + cant_edges[row[0]], :] = row[new_connections[row[0]] + 1] - 1
        last = last + cant_edges[row[0]]

    #plot:
    fig = mlab.figure(1, size=(400, 400), bgcolor=(1, 1, 1))

    src = mlab.pipeline.scalar_scatter(x, y, z)
    
    src.mlab_source.dataset.lines = connections 

    mlab.pipeline.surface(src, color=colors)

    src.update()

    # Put some balls to visualize the axes
    mlab.points3d(
        [0, 1, 0, 0],
        [0, 0, 1, 0],
        [0, 0, 0, 1],
        scale_factor=.06,
        color=(.5, 1, .5)
    )

    if isolated_points:
        isolated_points = np.loadtxt(isolated_points, delimiter=vert_delim)
        
        mlab.points3d(
            isolated_points[:,0],
            isolated_points[:,1],
            isolated_points[:,2], 
            np.linspace(.001, .3, len(isolated_points)),
            scale_factor=1, 
            color=(1, 0, 0)
        )

    mlab.show()
