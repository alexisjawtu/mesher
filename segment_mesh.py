#!/usr/bin/python3

import numpy as np
import mesh
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

def plot_oct(mu,levels):
	vert = np.array([[0,0,0],[1,0,0],[0,1,0],[0,0,1]]).transpose()
	points = mesh.macroel_tetrahedra(vert,mu,levels)
	pbc_tetra(plt,points)
	p_new = []
	for k in range(points.shape[0]):
		punto = points[k,:]
		if np.linalg.norm(punto)>0:
			scale = np.sum(punto)/(np.linalg.norm(punto))
		else:
			scale = 1
		p_new.append(punto*scale)
	points_r = np.array(p_new)
	pbc_tetra(plt,points_r)
	r3 = np.sqrt(3)
	points_e = np.concatenate([r3*points_r[:,0,None],r3*points_r[:,1,None],2*points_r[:,2,None]],1)
	print(points_e)

	pbc_tetra(plt,points_e)
	


def pbc_tetra(plt_axes,points):
	fig = plt_axes.figure()
	ax  = fig.add_subplot(1,1,1,projection='3d')
	ax.plot([],[],[],label = " ")
	legend = ax.legend()
	ax.set_xlabel(' X ')
	ax.set_ylabel(' Y ')
	ax.set_zlabel(' Z ')
	color_name = "blue"
	for i in [4*nn for nn in range(points.shape[0]//4)]:
		a = np.array([0,1,2,3]+[0,2,3,1])+i*np.ones(8,dtype=int)
		z = points[a,np.array([2]*8)]
		y = points[a,np.array([1]*8)]
		x = points[a,np.array([0]*8)]
		ax.plot(x, y, z,color=color_name)
	plt.show()
