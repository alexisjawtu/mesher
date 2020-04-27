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
from mayavi import mlab
import main
#import ellipse

def plot(vertices_file,connectivity_file,col=(0.5,0,0.5)):
    # Reads .ver and .ebv files. col is a color definition. 
    p = np.loadtxt(vertices_file)
    with open(connectivity_file,'r') as infile:
        inlist = infile.readlines()
    con_list = [line.strip(' \n').split(' ') for line in inlist]
    con_list = [[int(st[k]) for k in range(len(st))] for st in con_list]
    x = p[:,0]
    y = p[:,1]
    z = p[:,2]
    #s could be modified in order to give different colors to different elements. 
    s = np.ones(len(x),dtype="float64")
    cant_edges = {6:9,5:8,4:6} #dictionary: #nodes:#edges
    n_con = 0
    for i in range(len(con_list)):
        n_con = n_con + cant_edges[con_list[i][0]]
    connections = np.zeros((n_con,2))
    # List of connections for each type of element
    connections_prism = np.array([[0,1,2,0,3,4,5,1,2],[1,2,0,3,4,5,3,4,5]]).T
    connections_tetra = np.array([[0,1,2,3,0,1],[1,2,3,0,2,3]]).T
    connections_pyrad = np.array([[0,1,2,3,0,1,2,3],[1,2,3,0,4,4,4,4]]).T
    # dictionary
    new_connections = {6:connections_prism,5:connections_pyrad,4:connections_tetra}
    last = 0
    for i in range(len(con_list)):
        row = np.array(con_list[i])
        #add connections depending on the type of elelement
        connections[last:last+cant_edges[row[0]],:] = row[new_connections[row[0]]+1]-1
        last = last + cant_edges[row[0]]
    #plot:
    fig = mlab.figure(1, size=(400, 400), bgcolor=(1, 1, 1))
    src = mlab.pipeline.scalar_scatter(x, y, z)
    src.mlab_source.dataset.lines = connections 
    src.update()
    mlab.pipeline.surface(src,color=col)
    mlab.show()

"""
#######
# These funtions should be moved to another file. examples.py?
def ball(file_name,levels=3,mu=1,plotear=True):
    # Creates a ball (centered at the origin, and with radious 1), from scratch. It is graded toward the center.
    with open(file_name,'w') as inp:
        # Create the partition file, with 8 tetrahedra.
        inp.write('1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, '+str(mu)+'\n')
        inp.write('1, 0, 0, 0, 0, 1, 0,-1, 0, 0, 0, 0, 1, '+str(mu)+'\n')
        inp.write('1, 0, 0, 0,-1, 0, 0, 0,-1, 0, 0, 0, 1, '+str(mu)+'\n')
        inp.write('1, 0, 0, 0, 1, 0, 0, 0,-1, 0, 0, 0, 1, '+str(mu)+'\n')
        inp.write('1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0,-1, '+str(mu)+'\n')
        inp.write('1, 0, 0, 0, 0, 1, 0,-1, 0, 0, 0, 0,-1, '+str(mu)+'\n')
        inp.write('1, 0, 0, 0,-1, 0, 0, 0,-1, 0, 0, 0,-1, '+str(mu)+'\n')
        inp.write('1, 0, 0, 0, 1, 0, 0, 0,-1, 0, 0, 0,-1, '+str(mu)+'\n')
    # Compute mesh:
    main.omega(file_name,levels)
    # Push nodes to the the edge: 
    p = np.loadtxt(file_name+'.ver')
    nz = np.where(np.any(p!=0,axis=1))[0] # nodes different from (0,0,0)
    # First 4 octants ant its oposites are marked with 1 and -1, respectively. 
    oct1 = 1*((p[nz,0]>=0)*(p[nz,1]>=0)*(p[nz,2]>=0))-1*((p[nz,0]<=0)*(p[nz,1]<=0)*(p[nz,2]<0))
    oct2 = 1*((p[nz,0]<0)*(p[nz,1]>=0)*(p[nz,2]>=0))-1*((p[nz,0]>0)*(p[nz,1]<=0)*(p[nz,2]<0))
    oct3 = 1*((p[nz,0]<=0)*(p[nz,1]<0)*(p[nz,2]>=0))-1*((p[nz,0]>=0)*(p[nz,1]>0)*(p[nz,2]<0))
    oct4 = 1*((p[nz,0]>0)*(p[nz,1]<0)*(p[nz,2]>=0))-1*((p[nz,0]<0)*(p[nz,1]>0)*(p[nz,2]<0))
    norma = np.sqrt(np.sum(p[nz,:]**2,1))
    # Value of the constants that defines the planes where each node lie. 
    s1 = np.sum(p[nz,:],1)
    s2 = -p[nz,0]+p[nz,1]+p[nz,2]
    s3 = -p[nz,0]-p[nz,1]+p[nz,2]
    s4 = p[nz,0]-p[nz,1]+p[nz,2]
    # Correct nodes possitions:
    p[nz,:] = p[nz,:]*((oct1*s1+oct2*s2+oct3*s3+oct4*s4)/norma).reshape((p[nz,:].shape[0],1))
    np.savetxt(file_name+'.ver',p)
    if plotear:
        plot(file_name+'.ver',file_name+'.ebv')

def many_balls(_levels,_mu,folder=''):
    for levels in _levels:
        for mu in _mu:
            name = folder+'/ball_'+str(levels)+'_'+str(mu)
            ball(name,levels,mu,False)


def ellipsoid(file_name,levels,mu,plotear=True):
    with open(file_name,'w') as inp:
        # Polo Sur
        inp.write('1, 0, 0, -1, 1.5, 0, -1, 0, 1.5, -1, 0, 0, -2, '+str(mu)+'\n')
        inp.write('1, 0, 0, -1, 0, 1.5, -1,-1.5, 0, -1, 0, 0, -2, '+str(mu)+'\n')
        inp.write('1, 0, 0, -1,-1.5, 0, -1, 0,-1.5, -1, 0, 0, -2, '+str(mu)+'\n')
        inp.write('1, 0, 0, -1, 0,-1.5, -1, 1.5, 0, -1, 0, 0, -2, '+str(mu)+'\n')
        # Polo Norte
        inp.write('1, 0, 0, 1, 1.5, 0, 1, 0, 1.5, 1, 0, 0, 2, '+str(mu)+'\n')
        inp.write('1, 0, 0, 1, 0, 1.5, 1,-1.5, 0, 1, 0, 0, 2, '+str(mu)+'\n')
        inp.write('1, 0, 0, 1,-1.5, 0, 1, 0,-1.5, 1, 0, 0, 2, '+str(mu)+'\n')
        inp.write('1, 0, 0, 1, 0,-1.5, 1, 1.5, 0, 1, 0, 0, 2, '+str(mu)+'\n') 
        # Cilindro inferior:
        inp.write('0, 0, 0, 0, 1.5,    0,  0,    0, 1.5,  0,   0,    0, -1, '+str(mu)+'\n')
        inp.write('1, 0, 0,-1, 1.5,    0, -1,    0, 1.5, -1, 1.5,    0,  0, '+str(mu)+'\n')
        inp.write('1, 0, 0,-1,   0,  1.5, -1,    0, 1.5,  0, 1.5,    0,  0, '+str(mu)+'\n')
        inp.write('0, 0, 0, 0,   0,  1.5,  0, -1.5,   0,  0,   0,    0, -1, '+str(mu)+'\n')
        inp.write('1, 0, 0,-1,   0,  1.5, -1, -1.5,   0, -1,   0,  1.5,  0, '+str(mu)+'\n')
        inp.write('1, 0, 0,-1,   0,  1.5,  0, -1.5,   0,  0,-1.5,    0, -1, '+str(mu)+'\n')
        inp.write('0, 0, 0, 0,   0, -1.5,  0, -1.5,   0,  0,   0,    0, -1, '+str(mu)+'\n')
        inp.write('1, 0, 0,-1,   0, -1.5, -1, -1.5,   0, -1,   0, -1.5,  0, '+str(mu)+'\n')
        inp.write('1, 0, 0,-1,   0, -1.5,  0, -1.5,   0,  0,   0,    0,  0, '+str(mu)+'\n')
        inp.write('0, 0, 0, 0,   0, -1.5,  0,  1.5,   0,  0,   0,    0, -1, '+str(mu)+'\n')
        inp.write('1, 0, 0,-1,   0, -1.5, -1,  1.5,   0, -1,   0, -1.5,  0, '+str(mu)+'\n')
        inp.write('1, 0, 0,-1,   0, -1.5,  0,  1.5,   0,  0,   0,    0,  0, '+str(mu)+'\n')
        # Cilindro superior
        ## inp.write('0, 0, 0, 0, 1.5,    0,  0,    0, 1.5,  0,   0,    0,  1, '+str(mu)+'\n')
        ## inp.write('1, 0, 0, 1, 1.5,    0,  1,    0, 1.5,  1, 1.5,    0,  0, '+str(mu)+'\n')
        ## inp.write('1, 0, 0, 1,   0,  1.5,  1,    0, 1.5,  0, 1.5,    0,  0, '+str(mu)+'\n')
        ## inp.write('0, 0, 0, 0,   0,  1.5,  0, -1.5,   0,  0,   0,    0,  1, '+str(mu)+'\n')
        ## inp.write('1, 0, 0, 1,   0,  1.5,  1, -1.5,   0,  1,   0,  1.5,  0, '+str(mu)+'\n')
        ## inp.write('1, 0, 0, 1,   0,  1.5,  0, -1.5,   0,  0,-1.5,    0,  1, '+str(mu)+'\n')
        ## inp.write('0, 0, 0, 0,   0, -1.5,  0, -1.5,   0,  0,   0,    0,  1, '+str(mu)+'\n')
        ## inp.write('1, 0, 0, 1,   0, -1.5,  1, -1.5,   0,  1,   0, -1.5,  0, '+str(mu)+'\n')
        ## inp.write('1, 0, 0, 1,   0, -1.5,  0, -1.5,   0,  0,   0,    0,  0, '+str(mu)+'\n')
        ## inp.write('0, 0, 0, 0,   0, -1.5,  0,  1.5,   0,  0,   0,    0,  1, '+str(mu)+'\n')
        ## inp.write('1, 0, 0, 1,   0, -1.5,  1,  1.5,   0,  1,   0, -1.5,  0, '+str(mu)+'\n')
        ## inp.write('1, 0, 0, 1,   0, -1.5,  0,  1.5,   0,  0,   0,    0,  0, '+str(mu)+'\n')
        

    main.omega(file_name,levels)
    p = np.loadtxt(file_name+'.ver') 
    p = p.T  
    q,coeff,height,classification = ellipse.project_and_scale(p)
    q = ellipse.stereo_projection(q)
    q = ellipse.scale_ellipse(q,coeff,height,classification)
    #q = np.vstack((q,np.zeros(q.shape[1])))
    np.savetxt(file_name+'.ver',q.T)
    if plotear:
        plot(file_name+'.ver',file_name+'.ebv')
"""