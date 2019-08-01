import numpy as np
#from mesh import n_elem_macrotetra

__format__ = '%8.5f'

def write_element_indices (file_name, levels):
	""" Depends only on n. Not on which macroelement or octant is. 

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
	with open (f_write,'ab') as out:
		for l in range(L):
			for i in range(L-l):
				for j in range(L-l-i):
					np.savetxt(out, points[l,i,:,j].reshape((1,3)), fmt = __format__)
					nVertices += 1
		#out.write('\n')
	return nVertices

def vertices_macro_tetra (points, f_write):
	""" points: dictionary with the points of the 
	tetrahedral non--hybrid macro--element """
	L = len(points)
	with open (f_write, 'ab') as out:
		np.savetxt(out, points, fmt = __format__)
	return L		

def vertices_macro_prism (points, f_write):
    return 0 # <------ local_n_vertices


"""
TODO: deprecate this
def write_element_coordinates (index_file, tetr):
	# I don't merge this with write_element_indices because indices are the same everywhere
	indices = [] 
	with open(index_file, 'r') as data:
		lines = data.readlines()
	for line in lines:
		indices.append(np.fromstring(line, dtype=int, sep=' '))
	with open('opc_2_elem_coordinates.txt','w') as out:
		out.write('# Esto se lee:\n#        x           y           z\n')
		for elem in indices:   # level, i, coord, j elem[0] tells what polyhedron it is
			for j in range(elem[0]):
				np.savetxt(out, tetr[elem[3*j+1],elem[3*j+2],:,elem[3*j+3]].reshape((1,3)), fmt='%.8f')
			out.write('# paso de un elemento a otro__\n')
	with open('elem_coordinates.txt','w') as out: # option 2 is to make a matrix followed by blank line, etc.
		for elem in indices:   # level, i, coord, j
			aux = np.zeros((3,elem[0]))
			for j in range(elem[0]):
				aux[:,j] = tetr[elem[3*j+1],elem[3*j+2],:,elem[3*j+3]]
			np.savetxt(out, aux, fmt=__format__)
	return indices
"""
