import numpy as np                              
import scipy.io as sio
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


mu_ 		= .65
p_ 			= np.array([[ 1,  1,  1,  1,  0,  0,  0,  0], 	# octante 8
       					[-1,  0,  0, -1, -1,  0,  0, -1],   # x > 0, y < 0, z < 0
       					[-1, -1,  0,  0, -1, -1,  0,  0]])

## organization of the 5 tetrahedra resolving a cube: 
## the vertices of the macro--cube are: 0 ... 7.
macro_el 	= np.array([[0,1,3,4],[2,3,1,6],[7,3,4,6],[5,4,1,6],[6,1,3,4]])

def n_elem_macrotetra (lev):
	""" sum([sum([(2*j+1) for j in range(l)]) for l in range(1,a+1)]) 

	this works just for the macro_element with both singularities

	"""
	dic = {}
	dic['number_of_elements'] 	= lev*(lev+1)*(2*lev+1)/6
	dic['number_of_prisms']		= lev*(lev*(2*lev-3)+1)/6
	dic['number_of_tetr']		= lev*(lev+1)/2 
	dic['number_of_pyrs']		= lev*(lev-1)/2
	dic['number_of_vertices']	= sum([r*(r+1)/2 for r in range(lev+2)])
	return dic

def n_faces_macrotetra (lev):
	Nel = n_elem_macrotetra(lev)['number_of_elements']
	return 2*lev**2+lev*(lev+1)+Nel-lev**2+2*(lev-1)**2+(lev-1)*lev*(lev+1)/6

def octant (o, points): # mas corta esta funcion ??
	## takes a fixed octant and affine--transforms it to the other six.
	if   o == 2: q = points*np.array([-1,-1,-1]).reshape((3,1))
	elif o == 3: q = points*np.array([-1, 1,-1]).reshape((3,1))
	elif o == 4: q = points*np.array([ 1, 1,-1]).reshape((3,1))
	elif o == 5: q = points*np.array([ 1,-1, 1]).reshape((3,1))
	elif o == 6: q = points*np.array([-1,-1, 1]).reshape((3,1))
	elif o == 7: q = points*np.array([-1, 1, 1]).reshape((3,1))
	else:		 q = points
	return q

def lambda1 (i, j, k, n, mu):
	return float(i)/n * (float(i+j+k)/n)**((1/float(mu))-1)

def lambda2 (i, j, k, n, mu):
	return float(j)/n * (float(i+j+k)/n)**((1/float(mu))-1)

def lambda3 (i, j, k, n, mu):
	return float(k)/n * (float(i+j+k)/n)**((1/float(mu))-1)

def macroel_sing_vrtx (P0, P1, P2, P3, mu, n):
	"""	
		return value: Nel == nmbr of elmts

		TODO:  las cuentas tipo lambda_[0,i,j,k]*P0 tambien se pueden hacer 
		juntas al ppio y llamar.

	"""
	P0 = np.array(P0).reshape((1,3))
	P1 = np.array(P1).reshape((1,3))
	P2 = np.array(P2).reshape((1,3))
	P3 = np.array(P3).reshape((1,3))

	lambda_ = np.zeros((4, n+1, n+1, n+1))

	points  = np.zeros((1,3))



## hacer esto:falta acomodar los indices

	Nel 	= 0
	for k in range(n+1):
		for j in range(n+1):
			for i in range(n+1):
				lambda_[1,i,j,k] = lambda1(i,j,k,n,mu)
				lambda_[2,i,j,k] = lambda2(i,j,k,n,mu)
				lambda_[3,i,j,k] = lambda3(i,j,k,n,mu)
				lambda_[0,i,j,k] = 1 - lambda_[1,i,j,k] - lambda_[2,i,j,k] - lambda_[3,i,j,k]

	# now: element vertices
	for k in range(n):
		for j in range(n-k):
			for i in range(n-j-k):
				q0 = lambda_[0,i,j,k]  *P0 + lambda_[1,i,j,k]  *P1 + lambda_[2,i,j,k]*P2 + lambda_[3,i,j,k]*P3
				q1 = lambda_[0,i+1,j,k]*P0 + lambda_[1,i+1,j,k]*P1 + lambda_[2,i+1,j,k]*P2 + lambda_[3,i+1,j,k]*P3
				q2 = lambda_[0,i,j+1,k]*P0 + lambda_[1,i,j+1,k]*P1 + lambda_[2,i,j+1,k]*P2 + lambda_[3,i,j+1,k]*P3
				q3 = lambda_[0,i,j,k+1]*P0 + lambda_[1,i,j,k+1]*P1 + lambda_[2,i,j,k+1]*P2 + lambda_[3,i,j,k+1]*P3

				points = np.concatenate((points,q0,q1,q2,q3))
				Nel += 1

	for k in range(n-1):
		for j in range(n-1-k):
			for i in range(n-1-j-k):
				q0 = lambda_[0,i+1,j,k]*P0 + lambda_[1,i+1,j,k]*P1 + lambda_[2,i+1,j,k]*P2 + lambda_[3,i+1,j,k]*P3
				q1 = lambda_[0,i,j+1,k]*P0 + lambda_[1,i,j+1,k]*P1 + lambda_[2,i,j+1,k]*P2 + lambda_[3,i,j+1,k]*P3
				q2 = lambda_[0,i,j,k+1]*P0 + lambda_[1,i,j,k+1]*P1 + lambda_[2,i,j,k+1]*P2 + lambda_[3,i,j,k+1]*P3
				q3 = lambda_[0,i+1,j,k+1]*P0 + lambda_[1,i+1,j,k+1]*P1 + lambda_[2,i+1,j,k+1]*P2 + lambda_[3,i+1,j,k+1]*P3

				points = np.concatenate((points,q0,q1,q2,q3))
				Nel += 1

				q0 = lambda_[0,i,j+1,k]*P0 + lambda_[1,i,j+1,k]*P1 + lambda_[2,i,j+1,k]*P2 + lambda_[3,i,j+1,k]*P3
				q1 = lambda_[0,i,j,k+1]*P0 + lambda_[1,i,j,k+1]*P1 + lambda_[2,i,j,k+1]*P2 + lambda_[3,i,j,k+1]*P3
				q2 = lambda_[0,i+1,j,k+1]*P0 + lambda_[1,i+1,j,k+1]*P1 + lambda_[2,i+1,j,k+1]*P2 + lambda_[3,i+1,j,k+1]*P3
				q3 = lambda_[0,i,j+1,k+1]*P0 + lambda_[1,i,j+1,k+1]*P1 + lambda_[2,i,j+1,k+1]*P2 + lambda_[3,i,j+1,k+1]*P3

				points = np.concatenate((points,q0,q1,q2,q3))
				Nel += 1
				
				q0 = lambda_[0,i+1,j,k]*P0 + lambda_[1,i+1,j,k]*P1 + lambda_[2,i+1,j,k]*P2 + lambda_[3,i+1,j,k]*P3
				q1 = lambda_[0,i,j+1,k]*P0 + lambda_[1,i,j+1,k]*P1 + lambda_[2,i,j+1,k]*P2 + lambda_[3,i,j+1,k]*P3
				q2 = lambda_[0,i+1,j+1,k]*P0 + lambda_[1,i+1,j+1,k]*P1 + lambda_[2,i+1,j+1,k]*P2 + lambda_[3,i+1,j+1,k]*P3
				q3 = lambda_[0,i+1,j,k+1]*P0 + lambda_[1,i+1,j,k+1]*P1 + lambda_[2,i+1,j,k+1]*P2 + lambda_[3,i+1,j,k+1]*P3

				points = np.concatenate((points,q0,q1,q2,q3))
				Nel += 1
				
				q0 = lambda_[0,i,j+1,k]*P0 + lambda_[1,i,j+1,k]*P1 + lambda_[2,i,j+1,k]*P2 + lambda_[3,i,j+1,k]*P3
				q1 = lambda_[0,i+1,j+1,k]*P0 + lambda_[1,i+1,j+1,k]*P1 + lambda_[2,i+1,j+1,k]*P2 + lambda_[3,i+1,j+1,k]*P3
				q2 = lambda_[0,i+1,j,k+1]*P0 + lambda_[1,i+1,j,k+1]*P1 + lambda_[2,i+1,j,k+1]*P2 + lambda_[3,i+1,j,k+1]*P3
				q3 = lambda_[0,i,j+1,k+1]*P0 + lambda_[1,i,j+1,k+1]*P1 + lambda_[2,i,j+1,k+1]*P2 + lambda_[3,i,j+1,k+1]*P3
				
				points = np.concatenate((points,q0,q1,q2,q3))
				Nel += 1

	for k in range(n-2):
		for j in range(n-2-k):
			for i in range(n-2-j-k):
				q0 = lambda_[0,i+1,j+1,k]*P0 + lambda_[1,i+1,j+1,k]*P1 + lambda_[2,i+1,j+1,k]*P2 + lambda_[3,i+1,j+1,k]*P3
				q1 = lambda_[0,i+1,j,k+1]*P0 + lambda_[1,i+1,j,k+1]*P1 + lambda_[2,i+1,j,k+1]*P2 + lambda_[3,i+1,j,k+1]*P3
				q2 = lambda_[0,i,j+1,k+1]*P0 + lambda_[1,i,j+1,k+1]*P1 + lambda_[2,i,j+1,k+1]*P2 + lambda_[3,i,j+1,k+1]*P3
				q3 = lambda_[0,i+1,j+1,k+1]*P0 + lambda_[1,i+1,j+1,k+1]*P1 + lambda_[2,i+1,j+1,k+1]*P2 + lambda_[3,i+1,j+1,k+1]*P3

				points = np.concatenate((points,q0,q1,q2,q3))
				Nel += 1

	points = np.delete(points, 0, 0)
	return points

def line (x, y, z):
	""" (x[0], y[0], z[0]) -- ... -- (x[n-1], y[n-1], z[n-1]) """
	s = '\n\t\\draw '
	l = len(x)
	for i in range (l - 1):
		s += '('+str(x[i])[0:6]+','+str(y[i])[0:6]+','+str(z[i])[0:6]+')'+' -- '
	s += '('+str(x[l-1])[0:6]+','+str(y[l-1])[0:6]+','+str(z[l-1])[0:6]+');'
	return s

def cube2tex (obj, name = 'cube.tex'):
	s = ''
	for tetra in obj:
		for dr in tetra:
			s += line(dr[0],dr[1],dr[2])
	with open (name,'w') as d:
		d.write('\\documentclass{article}\n\\usepackage{tikz}\n\\begin{document}\n')
		d.write('\\begin{tikzpicture}[scale=3]\n')
		d.write(s)
		d.write('\n\n\\end{tikzpicture}\n\\end{document}')
	return

def cube2mat (obj, file_name = 'data.mat'):
	## to draw in Octave with cubo.m
	d = {}	
	i = -1
	for tetra in obj:
		for dr in tetra:
			i += 1
			d['dr'+str(i)] = [dr[0],dr[1],dr[2]]
	sio.savemat(file_name, d)
	return i

colours = ['lightgreen']*4

def plot_mesh (obj, elev = 30, azim = 45, colors = ['darkgreen']+['black']+['fuchsia']+['blue']):
	fig = plt.figure()
	ax 	= fig.add_subplot(1,1,1, projection='3d')
	ax.axis('equal')
	ax.view_init(elev, azim)
	#corregir esto!!
	## plt.title(str(n) + ' niveles.')
	###
	colors = iter(colors)
	for tetra in obj:
		col = next(colors)#for te in tetra:
		for dr in tetra:
			ax.plot(dr[0],dr[1],dr[2],color = col)
	plt.show()
	return fig

def cube_mesh_2 (n, mu, p, tetrahedra):
	""" here we calculate the mesh of the whole fichera	
		n == levels
		mu == grading param
		p == vertices of the cubes; the octants
		tetrahedra ==  
	"""
	mu_vec = [1, mu, mu, mu]            # first one is not graded!
	dict_save = {}
	for o in xrange(2,9):
		q = octant(o, p)
		for t in xrange(4):
			point0 = q[:,tetrahedra[t,0]]
			point1 = q[:,tetrahedra[t,1]]
			point2 = q[:,tetrahedra[t,2]]
			e3 	   = q[:,tetrahedra[t,3]]
			points = np.zeros((n+1, n+1, 3, n+1))  # nivel, i, coordenada, j
			for k in xrange(n+1):
				for i in xrange(n-k+1):
					for j in xrange(n-k-i+1):
						lambda_1 = lambda1 (i,j,0,n,mu_vec[t]); # se puede con muchas menos llamadas a estos lambda
						lambda_2 = lambda2 (i,j,0,n,mu_vec[t]);
						# the mesh are just the following two lines
						points[k,i,:,j] = lambda_1*(point1-point0) + lambda_2*(point2-point0)
						points[k,i,:,j] += (1-(float(n-k)/n)**(1/mu_vec[t]))*(e3-point0) + point0
			dict_save['points_T'+str(t)+'_C'+str(o)] = points.copy()

		# CALCULATION FOR t = 4 (T5)
		P0 	  = q[:,tetrahedra[4,0]]
		P1 	  = q[:,tetrahedra[4,1]]
		P2 	  = q[:,tetrahedra[4,2]]
		P3 	  = q[:,tetrahedra[4,3]]
		points = macroel_sing_vrtx(P0, P1, P2, P3, mu, n)
		dict_save['points_T4_C' + str(o)] = points.copy()

	return dict_save

def cube_drawing (coord, oct_range = range(2,9)):
	drawing = [[],[],[],[]]
	for o in oct_range:
		for t in [0,1,2,3]:
			points = coord['points_T'+str(t)+'_C'+str(o)] # np.zeros((n+1,n+1,3,n+1))  # nivel, i, coordenada, j
			n = points.shape[0] - 1
			for k in range (n+1):
				for i in range (n-k+1):
					for j in range (n-k-i+1):
						x, y, z = points[k,0:n-k-j+1,0,j], points[k,0:n-k-j+1,1,j], points[k,0:n-k-j+1,2,j]
						drawing[t].append([x,y,z])
			# ax.scatter3D(points[k][i][0,0:n-k-i+1],points[k][i][1,0:n-k-i+1],points[k][i][2,0:n-k-i+1]) # opcional
					x, y, z = points[k,i,0,0:n-k-i+1], points[k,i,1,0:n-k-i+1], points[k,i,2,0:n-k-i+1]
					drawing[t].append([x,y,z])
				for c in range (n-k+1):
					x, y, z = np.zeros(c+1), np.zeros(c+1), np.zeros(c+1)
					for i in range (c+1):
						x[i] = points[k,i,0,c-i]
						y[i] = points[k,i,1,c-i]
						z[i] = points[k,i,2,c-i]
					drawing[t].append([x,y,z])  #  transversales
			for i in range (n):
				for j in range (n-i):
					stop = n - (i + j) + 1
					x = points[0:stop,i,0,j]
					y = points[0:stop,i,1,j]
					z = points[0:stop,i,2,j]
					drawing[t].append([x,y,z])  #  verticales
				x, y, z = np.zeros(n-i+1), np.zeros(n-i+1), np.zeros(n-i+1)
				for k in range (n-i+1):
					x[k] = points[k,i,0,n-k-i]
					y[k] = points[k,i,1,n-k-i]
					z[k] = points[k,i,2,n-k-i]
				drawing[t].append([x,y,z])  #  piramidales
			# ahora simplemente intercambio los papeles de i y j
				x, y, z = np.zeros(n-i+1), np.zeros(n-i+1), np.zeros(n-i+1)
				for k in range (n-i+1):
					x[k] = points[k,n-k-i,0,i]
					y[k] = points[k,n-k-i,1,i]
					z[k] = points[k,n-k-i,2,i]
				drawing[t].append([x,y,z])	#  piramidales
	return np.array(drawing)

#########################################################################################################
## Drawing test
##
## n = 3
## coords 	= cube_mesh_2(n,mu_,p_,macro_el)
## drawing = cube_drawing(coords)
## fig 	= plot_mesh(drawing, elev = 30, azim = 45)
## fig.savefig('refine' + str(n) + '.png')