import numpy as np
import mesh
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt


plot_function = {
	4: 'pbc_tetra',
	5: 'pbc_pyramid',
	6: 'pbc_prism'
}

def build_pill(mu,levels):
	vert_tetra = np.array([[0,0,0.5],[1,0,0.5],[0,1,0.5],[0,0,1.5]]).transpose()
	points = mesh.macroel_tetrahedra(vert_tetra,mu,levels)
	
	vert_prism = np.array([[0,0,0.5],[0,0,0],[1,0,0],[1,0,0.5],[0,1,0],[0,1,0.5]])
	nodos = macro_prism(vert_prism,mu,levels)
	nodos.update({'t0':points})

	return nodos


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
	print(vertices.shape)
	t1_tetra = vertices[[0,2,3,4],:]
	t2_tetra = vertices[[0,3,4,5],:]
	t_prism = vertices[[1,2,4,0],:]
	points_t1 = mesh.macroel_tetrahedra(t1_tetra.transpose(),mu,levels)
	points_t2 = mesh.macroel_tetrahedra(t2_tetra.transpose(),mu,levels)
	points_t_p = mesh.macroel_hybrid(t_prism.transpose(),mu,levels)
	dict_out['tp'] = points_t_p
	dict_out['t1'] = points_t1
	dict_out['t2'] = points_t2
	return dict_out


def pbc_tetra(ax,points,tetra):
	a = np.array([0,1,2,3]+[0,2,3,1])
	z = points[a,np.array([2]*8)]
	y = points[a,np.array([1]*8)]
	x = points[a,np.array([0]*8)]
	ax.plot(x, y, z,color='blue')


def pbc_pyramid(ax,points,tetra):
	a = np.array([0,1,2,0,3,4,0,2,4,3,1])
	z = points[a,np.array([2]*8)]
	y = points[a,np.array([1]*8)]
	x = points[a,np.array([0]*8)]
	ax.plot(x, y, z,color='black')

def pbc_prism(ax,points,prism):
	a = np.array([0,1,2,5,4,3,0,1,4,3,5])
	x = points[a,np.array([0]*11)]
	y = points[a,np.array([1]*11)]
	z = points[a,np.array([2]+11)]
	ax.plot(x,y,z,color='green')

def pbc_all(file_points,file_conectivity):
	points = np.fromfile(file_points)
	with open(file_conectivity,'r') as elements:
		elems = elements.readlines()
	n_elems = len(elems)
	
	fig = plt.figure()
	ax  = fig.add_subplot(1,1,1,projection='3d')
	ax.plot([],[],[],label = " ")
	legend = ax.legend()
	ax.set_xlabel(' X ')
	ax.set_ylabel(' Y ')
	ax.set_zlabel(' Z ')

	for r in range(n_elems):
		 elem = [int(x) for x in elems[r].rstrip().split(' ')]
		 plot_function[elem[0]](ax,points,elem[1:])

	plt.show()


