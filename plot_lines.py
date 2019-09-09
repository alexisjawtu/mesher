# Author: Gael Varoquaux <gael dot varoquaux at normalesup dot org>
# Copyright (c) 2010, Enthought
# License: BSD style

import numpy as np
from mayavi import mlab


def plot_macroel(vertices_file,connectivity_file,col=(1,1,1)):
    p = np.loadtxt(vertices_file)
    with open(connectivity_file,'r') as infile:
        inlist = infile.readlines()
    con_list = [line.strip(' \n').split(' ') for line in inlist]
    con_list = [[int(st[k]) for k in range(len(st))] for st in con_list]
    x = p[:,0]
    y = p[:,1]
    z = p[:,2]
    s = np.ones(len(x),dtype="float64")
    cant_edges = {6:9,5:8,4:6}
    n_con = 0
    for i in range(len(con_list)):
        n_con = n_con + cant_edges[con_list[i][0]]
    connections = np.zeros((n_con,2))
    
    connections_prism = np.array([[0,1,2,0,3,4,5,1,2],[1,2,0,3,4,5,3,4,5]]).T
    connections_tetra = np.array([[0,1,2,3,0,1],[1,2,3,0,2,3]]).T
    connections_pyrad = np.array([[0,1,2,3,0,1,2,3],[1,2,3,0,4,4,4,4]]).T
    new_connections = {6:connections_prism,5:connections_pyrad,4:connections_tetra}
    last = 0
    for i in range(len(con_list)):
        row = np.array(con_list[i])
        connections[last:last+cant_edges[row[0]],:] = row[new_connections[row[0]]+1]-1
        last = last + cant_edges[row[0]]
    #fig = mlab.figure(1, size=(400, 400), bgcolor=(1, 1, 1))
    #mlab.clf()
    src = mlab.pipeline.scalar_scatter(x, y, z)
    src.mlab_source.dataset.lines = connections 
    src.update()
    mlab.pipeline.surface(src,color=col)
    mlab.show()

def pildora():
    pass