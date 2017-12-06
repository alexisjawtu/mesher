import numpy as np
from mesh import n_elem_macrotetra

__format__ = '%8.5f'

def write_element_indices (file_name, levels):
	""" Depende solo de n. No depende de cual macrotetraedro o cubo es. 

	takes n and writes on disc:

	file_name.txt

	---------------------------------------------

	6|(l,i,j)_{prism}|(l,i,j)_{prism}|(l,i,j)_{prism}|(l,i,j)_{prism}|(l,i,j)_{prism}|(l,i,j)_{prism}|

	4|(l,i,j)_{tetra}|(l,i,j)_{tetra}|(l,i,j)_{tetra}|(l,i,j)_{tetra}|

	5|(l,i,j)_{pyra} |(l,i,j)_{pyra} |(l,i,j)_{pyra} |(l,i,j)_{pyra} |(l,i,j)_{pyra} |

	etc...

	This funtion is initial in the program, for the case
	we start with the mesh we proposed. For a general mesh, the algorithm
	starts directly in the next step (with element_by_vertices.txt given 
	somehow).

	"""
	string 		= ''
	n 			= levels
	for level in range(n-1):
		for i in range(n-level-1): # prismas
			for j in range(n-level-i-1):  
				nodos = [6,level,i,j,level,i+1,j,level,i,j+1,level+1,i,j,
				 		level+1,i+1,j,level+1,i,j+1]
				string += ' '.join([str(r) for r in nodos]) +'\n'
		for i in range(n-level-2):
			for j in range(n-level-i-2):
				nodos = [6,level,i+1,j,level,i+1,j+1,level,i,j+1,
						level+1,i+1,j,level+1,i+1,j+1,level+1,i,j+1]
				string += ' '.join([str(r) for r in nodos]) +'\n'
		for i in range(n - level): # tetraedros
			nodos = [4,
					level, i, n-level-i-1,
					level, i+1, n-level-i-1,
					level, i, n-level-i,
					level+1, i, n-level-i-1]
			string += ' '.join([str(r) for r in nodos]) + '\n'
		for i in range(n - level - 1): # piramides: invariante: level+i+j == n-1
			nodos = [5, level,   i,   n-i-1-level,
						level,   i+1, n-i-2-level,
						level+1, i+1, n-i-2-level,
						level+1, i,   n-i-1-level,
						level,   i+1, n-i-1-level]
			string += ' '.join([str(r) for r in nodos]) + '\n'
	nodos = [4,n-1,0,0,n-1,1,0,n-1,0,1,n,0,0] # tetraedro singular
	string += ' '.join([str(r) for r in nodos])
	with open(file_name, 'w') as elements:
		elements.write(string)
	return

def write_face_indices (n, file_name):
	""" Aparentemente est'a deprecada. i,j son los enteros que definen las
	coordenadas baricentricas. """
	face_type 	= 0 # 0 == triangular
	Nfaces 		= 0
	Nel 		= 0
	string = ''
	for level in range(n):		# caras horizontales
		for i in range(n-level):
			for j in range(n-level-i):
				Nfaces += 1
				Nel    += 1
				dat = [Nel,face_type,level,i,j,level,i+1,j,level,i,j+1]
				string += ' '.join([str(d) for d in dat]) + '\n'
		for i in range(n-level-1):
			for j in range(n-level-i-1):
				dat = [Nel,face_type,level,i+1,j,level,i+1,j+1,level,i,j+1]
				string += ' '.join([str(d) for d in dat]) + '\n'
				Nfaces += 1
				Nel    += 1
		# fin caras horizontales
		# caras triangulares no horizontales
	# fin caras triangulares no horizontales
	with open(file_name, 'w') as fn:
		fn.write(string)
	return Nfaces

def vertices (points):
	""" 
		
		This only works for type I macro--element

		writes the coordinates of the vertices of the mesh.
		the rows in the file output 'vertices.txt' define
		the global enumerations.

		open ( '...', mode = 'a') for appending

		return value: nmbr of vrtcs at the moment
	"""
	nVertices = 0
	L = points.shape[0]
	with open ('vertices.txt','a') as out:
		for l in range(L):
			for i in range(L-l):
				for j in range(L-l-i):
					np.savetxt(out, points[l,i,:,j].reshape((1,3)), fmt = __format__)
					nVertices += 1
		#out.write('\n')
	return nVertices

def verticesT4 (points, f_write):
	"""
		points: dictionary with the points of the type II macro--element
	"""
	L = len(points)
	with open (f_write, 'a') as out:
		np.savetxt(out, points, fmt = __format__)
	return L		

def write_element_coordinates (index_file, tetr):
	""" no la entremezclo con write_element_indices porque los indices
	son los mismos para todos lados """
	indices = [] 
	with open(index_file, 'r') as data:
		lines = data.readlines()
	for line in lines:
		indices.append(np.fromstring(line, dtype=int, sep=' '))
	# print(np.array	indices)  <<<---- para que queria saber esto ????
	with open('opc_2_elem_coordinates.txt','w') as out:
		out.write('# Esto se lee:\n#        x           y           z\n')
		for elem in indices:   # nivel, i, coordenada, j elem[0] dice que poliedro es
			for j in range(elem[0]):
				np.savetxt(out, tetr[elem[3*j+1],elem[3*j+2],:,elem[3*j+3]].reshape((1,3)), fmt='%.8f')
			out.write('# paso de un elemento a otro__\n')
	with open('elem_coordinates.txt','w') as out: # opc 2 es una matriz seguida de una linea en bco, etc.
		for elem in indices:   # nivel, i, coordenada, j
			aux = np.zeros((3,elem[0]))
			for j in range(elem[0]):
				aux[:,j] = tetr[elem[3*j+1],elem[3*j+2],:,elem[3*j+3]]
			np.savetxt(out, aux, fmt=__format__)
	return indices

def write_face_coordinates (indices, tetr):
	""" 
	indices es el valor devuelto por write_element_coordinates (index_file, tetr).
	De su primer argumento es importante la primera columna  
						
							[6,4,4,5,4]

	para saber de que poliedro se trata. tambien podria calcular una
	longitud y listo. n == tetr.shape[0] - 1 
		
	Para 0 <= k <= n(n+1)(2n+1)/6 el elemento k esta en las lineas

							3k, 3k+1, 3k+2

	del archivo elem_coordinates.txt.
	"""

	n_dic 	  = n_elem_macrotetra(tetr.shape[0] - 1)

	with open ('face_coordinates.txt','w') as fc:

		ini  = 0
		stop = n_dic['number_of_prisms']
		print(range (ini,stop))
		for pr in range (ini,stop):
			face = np.zeros((3,4))
			face[:,3]  	= pr*np.ones(3)
			## Poner aca algo que dependa de p0, p1, ... , p6 para
			## mayor legibilidad.
			##					  p0,p1,p2	
			z = tetr[indices[pr][[1,4,7]],indices[pr][[2,5,8]],:,indices[pr][[3,6,9]]]
			face[:,0:3] = np.transpose(z)
			np.savetxt(fc,face,fmt = __format__)
			##					  p3,p4,p5	
			z = tetr[indices[pr][[10,13,16]],indices[pr][[11,14,17]],:,indices[pr][[12,15,18]]]
			face[:,0:3] = np.transpose(z)
			np.savetxt(fc,face,fmt = __format__)

			face = np.zeros((3,5))
			face[:,4]  	= pr*np.ones(3)
			##					 p0,p1,p3,p4	
			z = tetr[indices[pr][[1,4,10,13]],indices[pr][[2,5,11,14]],:,indices[pr][[3,6,12,15]]]
			face[:,0:4] = np.transpose(z)
			np.savetxt(fc,face,fmt = __format__)
			##					 p0,p2,p3,p5	
			z = tetr[indices[pr][[1,7,10,16]],indices[pr][[2,8,11,17]],:,indices[pr][[3,9,12,18]]]
			face[:,0:4] = np.transpose(z)
			np.savetxt(fc,face,fmt = __format__)
			##					 p1,p2,p4,p5	
			z = tetr[indices[pr][[4,7,13,16]],indices[pr][[5,8,14,17]],:,indices[pr][[6,9,15,18]]]
			face[:,0:4] = np.transpose(z)
			np.savetxt(fc,face,fmt = __format__)

		face 		= np.zeros((3,4))
		
		ini 	= ini + stop
		stop 	= stop + n_dic['number_of_tetr'] - 1

		#print(range(ini,stop) + [n_dic['number_of_elements'] - 1])
		for te in range(ini,stop) + [n_dic['number_of_elements'] - 1]: ## the top tetrahedron appears last
		##	lo dejo asi por legibilidad. Hay una manera de hacer esto con algoritmos
		##  en permutaciones???

		##
		## tal vez convenga pasar por alto el indices <--> elements.txt
		##

		##	estoy contando cada cara una vez por cada elemento que la tiene en su borde 

			face[:,3] 	= te*np.ones(3)			
			##					p0,p1,p2
			z = tetr[indices[te][[1,4,7]],indices[te][[2,5,8]],:,indices[te][[3,6,9]]]
			face[:,0:3] = np.transpose(z)
			np.savetxt(fc,face,fmt = __format__)
			##					 p0,p1,p3
			z = tetr[indices[te][[1,4,10]],indices[te][[2,5,11]],:,indices[te][[3,6,12]]]
			face[:,0:3] = np.transpose(z)
			np.savetxt(fc,face,fmt = __format__)
			##					 p0,p2,p3
			z = tetr[indices[te][[1,7,10]],indices[te][[2,8,11]],:,indices[te][[3,9,12]]]
			face[:,0:3] = np.transpose(z)
			np.savetxt(fc,face,fmt = __format__)
			##					 p1,p2,p3
			z = tetr[indices[te][[4,7,10]],indices[te][[5,8,11]],:,indices[te][[6,9,12]]]
			face[:,0:3] = np.transpose(z)
			np.savetxt(fc,face,fmt = __format__)

		ini  = ini + stop - 1
		stop = n_dic['number_of_elements'] - 1
		print(range(ini,stop))		
		for py in range(ini,stop):
			face[:,3]  	= py*np.ones(3)
			##					 p0,p1,p4
			z = tetr[indices[py][[1,4,13]],indices[py][[2,5,14]],:,indices[py][[3,5,15]]]
			face[:,0:3] = np.transpose(z)
			np.savetxt(fc,face,fmt = __format__)
			##					 p0,p3,p4
			z = tetr[indices[py][[1,10,13]],indices[py][[2,11,14]],:,indices[py][[3,12,15]]]
			face[:,0:3] = np.transpose(z)
			np.savetxt(fc,face,fmt = __format__)
			##					 p2,p3,p4
			z = tetr[indices[py][[7,10,13]],indices[py][[8,11,14]],:,indices[py][[9,12,15]]]
			face[:,0:3] = np.transpose(z)
			np.savetxt(fc,face,fmt = __format__)
			##					 p1,p2,p4
			z = tetr[indices[py][[4,7,13]],indices[py][[5,8,14]],:,indices[py][[6,9,15]]]
			face[:,0:3] = np.transpose(z)
			np.savetxt(fc,face,fmt = __format__)
			face 		= np.zeros((3,5)) # square face
			face[:,4] 	= py*np.ones(3)
			##					 p0,p1,p2,p3
			z = tetr[indices[py][[1,4,7,10]],indices[py][[2,5,8,11]],:,indices[py][[3,6,9,12]]]
			np.savetxt(fc,face,fmt = __format__)

		return