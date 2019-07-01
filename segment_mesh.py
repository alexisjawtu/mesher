#!/usr/bin/python3

import numpy as np
import mesh
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt


def build_pill(mu,levels):
	vert_tetra = np.array([[0,0,0.5],[1,0,0.5],[0,1,0.5],[0,0,1.5]]).transpose()
	points = mesh.macroel_tetrahedra(vert_tetra,mu,levels)
	
	vert_prism = np.array([[0,0,0.5],[0,0,0],[1,0,0],[1,0,0.5]],[0,1,0],[0,1,0.5]).transpose()
	nodos = macro_prism(vert_prism,mu,levels)
	nodos.update({'t0':points})

	pbc_tetra(plt,nodos['t0'],'blue')
	pbc_tetra(plt,nodos['t1'],'green')
	pbc_tetra(plt,nodos['t2'],'green')
	pbc_hybrid(plt,nodos['tp'],'black')

	pbc_tetra(plt,points,_color)
	plt.show()
	#p_new = []
	#for k in range(points.shape[0]):
	#	punto = points[k,:]
	#	if np.linalg.norm(punto)>0:
	#		scale = np.sum(punto)/(np.linalg.norm(punto))
	#	else:
	#		scale = 1
	#	p_new.append(punto*scale)
	#points_r = np.array(p_new)
	#pbc_tetra(plt,points_r)
	#r3 = np.sqrt(3)
	#points_e = np.concatenate([r3*points_r[:,0,None],r3*points_r[:,1,None],2*points_r[:,2,None]],1)
	#print(points_e)

	#pbc_tetra(plt,points_e)
	

def macro_prism(vertices,mu,levels):
	dict_out = {}
	t1_tetra = vertices[[0,2,3,4],:]
	t2_tetra = vertices[[0,3,4,5],:]
	t_prism = vertices[[1,2,4,0],:]
	points_t1 = mesh.macroel_tetrahedra(t1_tetra.transpose(),mu,levels)
	points_t2 = mesh.macroel_tetrahedra(t2_tetra.transpose(),mu,levels)
	points_t_p = mesh.macro_hybrid(t_prism.transpose(),mu,levels)
	dict_out['tp'] = points_t_p
	dict_out['t1'] = points_t1
	dict_out['t2'] = points_t2
	return dict_out

def pbc_hybrid(ax,points,_color):
	drawing = []
	n = points.shape[0] - 1
    for k in range (n+1):
        for i in range (n-k+1):
            for j in range (n-k-i+1):
                x, y, z = points[k,0:n-k-j+1,0,j], points[k,0:n-k-j+1,1,j], points[k,0:n-k-j+1,2,j]
                drawing.append([x,y,z])
            x, y, z = points[k,i,0,0:n-k-i+1], points[k,i,1,0:n-k-i+1], points[k,i,2,0:n-k-i+1]
            drawing.append([x,y,z])
        for c in range (n-k+1):
            x, y, z = np.zeros(c+1), np.zeros(c+1), np.zeros(c+1)
            for i in range (c+1):
                x[i] = points[k,i,0,c-i]
                y[i] = points[k,i,1,c-i]
                z[i] = points[k,i,2,c-i]
            drawing.append([x,y,z])  #  transversals
    for i in range (n):
        for j in range (n-i):
            stop = n - (i + j) + 1
            x = points[0:stop,i,0,j]
            y = points[0:stop,i,1,j]
            z = points[0:stop,i,2,j]
            drawing.append([x,y,z])  #  verticals
        x, y, z = np.zeros(n-i+1), np.zeros(n-i+1), np.zeros(n-i+1)
        for k in range (n-i+1):
            x[k] = points[k,i,0,n-k-i]
            y[k] = points[k,i,1,n-k-i]
            z[k] = points[k,i,2,n-k-i]
        drawing.append([x,y,z])  #  pyramidals
    # simply interchange the roles of i and j
        x, y, z = np.zeros(n-i+1), np.zeros(n-i+1), np.zeros(n-i+1)
        for k in range (n-i+1):
            x[k] = points[k,n-k-i,0,i]
            y[k] = points[k,n-k-i,1,i]
            z[k] = points[k,n-k-i,2,i]
        drawing.append([x,y,z])  #  pyramidals
    for dr in range(len(drawing)):
        ## TODO: this can be done putting (array_of_X, array_of_Y, array_of_Z, ...)
        ## and not one by one as is now
        ax.plot(drawing[dr][0],drawing[dr][1],drawing[dr][2], color = _color)
    return


def pbc_tetra(ax,points,_color):
	for i in [4*nn for nn in range(points.shape[0]//4)]:
		a = np.array([0,1,2,3]+[0,2,3,1])+i*np.ones(8,dtype=int)
		z = points[a,np.array([2]*8)]
		y = points[a,np.array([1]*8)]
		x = points[a,np.array([0]*8)]
		ax.plot(x, y, z,color=_color)
	


def pbc_all(data)
	fig = plt.figure()
	ax  = fig.add_subplot(1,1,1,projection='3d')
	ax.plot([],[],[],label = " ")
	legend = ax.legend()
	ax.set_xlabel(' X ')
	ax.set_ylabel(' Y ')
	ax.set_zlabel(' Z ')

