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


correction = {"C" : 0, "Octave" : 1}

def plot_lines(
        vertices_file: str, 
        connectivity_file: str, 
        isolated_points: str = None,
        vert_delim: str = None, 
        colors: Tuple = (.8,.8,.8)
    ) -> None:

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
        connections[last:last+cant_edges[row[0]],:] = row[new_connections[row[0]] + 1] - 1
        last = last + cant_edges[row[0]]
    #plot:
    fig = mlab.figure(1, size=(400, 400), bgcolor=(1, 1, 1))
    src = mlab.pipeline.scalar_scatter(x, y, z)
    src.mlab_source.dataset.lines = connections
    print("connections:")
    print(connections)
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
        isolated_points: str = None,
        number_of_isolated_points_to_draw: int = 3,
        vert_delim: str = ",",
        elem_delim: str = ",",
        colors: Tuple = (.2,.3,.4)
    ) -> None:

    vertices: np.array = np.loadtxt(vertices_file, delimiter=vert_delim)

    x = vertices[:, 0]
    y = vertices[:, 1]
    z = vertices[:, 2]

    del vertices
    
    all_elements: np.array = np.loadtxt(connectivity_file, delimiter=elem_delim)

    # Here we have all tetrahedra, which have 6 edges.
    n_con = 0
    for e in all_elements:
        n_con = n_con + 6

    connections = np.zeros((n_con,2))

    links: np.array = np.array([
                            [0,1,2,3,0,1],
                            [1,2,3,0,2,3]
                        ]).T

    last: int = 0

    try:
        for i in range(len(all_elements)):
            row = np.array(all_elements[i])
            connections[last: last + 6, :] = row[links]
            last = last + 6

    except IndexError:
        print(i, row)
        input()

    except KeyError:
        print(i, row)
        input()

    # plot:
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
        isolated_points = np.loadtxt(isolated_points, delimiter=vert_delim)[0:number_of_isolated_points_to_draw,:]
        
        mlab.points3d(
            isolated_points[:,0],
            isolated_points[:,1],
            isolated_points[:,2], 
            np.linspace(.01, .03, len(isolated_points)),
            scale_factor=1, 
            color=(1, 0, 0)
        )

    mlab.show()


def plot_separate(
        vertices_file_inner: str,
        elements_file_inner: str,
        vertices_file_surrounding: str = None,
        elements_file_surrounding: str = None,
        isolated_points: str = None,
        vert_delim: str = None,
        nodes_delim: str = None,
        colors_inner: Tuple = (.2,.7,.2),
        colors_surrounding: Tuple = (.7,.2,.2)
    ) -> None:
    
    """
    plot_separate(
        "experiments/wafer24ago22/nodes_inner.dat",
        "experiments/wafer24ago22/elements_inner.dat",
        "experiments/wafer24ago22/physical_vertices_16613627703281198.dat",
        "experiments/wafer24ago22/wafer_profile_elements_16613627703281198.dat"
    ) 
    """

    vertices_inner = np.loadtxt(vertices_file_inner, delimiter=",")

    x = vertices_inner[:, 0]
    y = vertices_inner[:, 1]
    z = vertices_inner[:, 2]

    del vertices_inner
    
    elements_inner: np.array = np.loadtxt(elements_file_inner, delimiter=",")

    # all_elements: np.array = np.vstack((elements_inner, elements))

    # Here we have all tetrahedra, which have 6 edges.
    cant_edges = {4: 6}
    n_con = 0
    for e in elements_inner:
        n_con = n_con + 6

    connections: np.array = np.zeros((n_con, 2))

    connections_tetra: np.array = np.array([
        [0,1,2,3,0,1],
        [1,2,3,0,2,3]
    ]).T

    new_connections = {4: connections_tetra}

    last = 0

    for i in range(len(elements_inner)):
        row = np.array(elements_inner[i])
        #add connections depending on the type of elelement

        # TODO: this [[]+1]-1 at the end of the line should be depending of lang = C/Octave
        connections[last: last + cant_edges[row[0]], :] = row[new_connections[row[0]] + 1] - 1
        last = last + cant_edges[row[0]]

    #plot:
    fig = mlab.figure(1, size=(400, 400), bgcolor=(.92, .92, .92))

    src = mlab.pipeline.scalar_scatter(x, y, z)
    
    src.mlab_source.dataset.lines = connections 

    mlab.pipeline.surface(src, color=colors_inner)

    src.update()

    # Put some balls to visualize the axes
    mlab.points3d(
        [0, 1, 0, 0],
        [0, 0, 1, 0],
        [0, 0, 0, 1],
        scale_factor=.06,
        color=(.5, 1, .5)
    )

    if vertices_file_surrounding:
        vertices_surrounding: np.array = np.loadtxt(vertices_file_surrounding, delimiter=",")
        elements_surrounding: np.array = np.loadtxt(elements_file_surrounding)
        # prepend a legacy 4
        rows: int = elements_surrounding.shape[0]
        elements_surrounding = np.hstack((4 * np.ones((rows, 1)), elements_surrounding))

        x_surr = vertices_surrounding[:, 0]
        y_surr = vertices_surrounding[:, 1]
        z_surr = vertices_surrounding[:, 2]

        del vertices_surrounding
        
        n_con_surr: int = 6 * elements_surrounding.shape[0]

        connections_surr = np.zeros((n_con_surr, 2))

        new_connections_surr = {4: connections_tetra}

        last_surr: int = 0

        for i in range(len(elements_surrounding)):
            row = np.array(elements_surrounding[i])

            # TODO: this [[]+1]-1 at the end of the line should be depending of lang = C/Octave

            try:
                connections_surr[last_surr: last_surr + 6, :] = row[connections_tetra + 1] - 1
                last_surr = last_surr + 6
            except:
                print(i, row, row.shape)
                exit()

        #plot:
        fig = mlab.figure(1, size=(400, 400), bgcolor=(.92, .92, .92))

        src = mlab.pipeline.scalar_scatter(x_surr, y_surr, z_surr)
        
        src.mlab_source.dataset.lines = connections_surr 

        mlab.pipeline.surface(src, color=colors_surrounding)

        src.update()

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
