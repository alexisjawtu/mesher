# hacer unas formulas
import numpy as np
#from boltons.iterutils import remap

def vertex_global_index (n, l, i, j):
	"""
	*indices in the macro--element with sing v and sing e into naturals

	*a one by one function that maps indices of the tensor
	into naturals 
	l = level
	"""
	#return ((n+1)**2)*level + (n+1)*i + j 

	## TAREA: esto tiene una formula cerrada corta. Hacerla.
	return sum([(n-r+1)*(n-r+2)/2 for r in range(l)]) + sum([n-l+1-k for k in range(i)]) + j

def write_elements_by_vertices (f_name, n, lang, initial):
	""" 
	This funtion is at the beggining in the program, for the case
	we start with the mesh we proposed. For a general mesh, the algorithm
	starts directly in the next step (with element_by_vertices.txt given 
	somehow).

	initial: tracks the number of vertices of the previous macro_element

	esto es solo para el macro elemento
	con un vertice singular y una arista singular en el np.array:

	Takes the file f_name written by
	mesh_write.write_element_indices() 

	n == current quantity of levels
	
	Writes GLOBAL INDICES per element appending in the file

	elements_by_vertices_repeated.txt:
	--------------------------------------------------------

	6|v1_{prism}|v2_{prism}|v3_{prism}|v4_{prism}|v5_{prism}|v6_{prism}|

	4|v1_{tetra}|v2_{tetra}|v3_{tetra}|v4_{tetra}|

	5|v1_{pyra} |v2_{pyra} |v3_{pyra} |v4_{pyra} |v5_{pyra} |

	etc...
	"""
	language  = {'C' : 0, 'Octave' : 1}
	with open (f_name, 'r') as inp:
		indices = inp.readlines()
	with open ('elements_by_vertices_repeated.txt', 'a') as out:
		for j in range(len(indices)):
			l = [int(c) for c in indices[j].rstrip().split(' ')]
			line = np.zeros((1,l[0]+1))
			line[0] = l[0]
			for i in range (1,l[0]+1):
				calc = vertex_global_index(n,l[3*(i-1)+1],l[3*(i-1)+2],l[3*(i-1)+3])
				line[0,i] =  initial + calc + language[lang]
			np.savetxt(out,line,fmt='%d')
	return len(indices)

def write_elements_by_vertices_T4 (n_vert_T4, init, f_name_write):
	"""
		here we assume that in the (Npts x 3) array of points
		the tetrahedra appear in order taking the rows
		four at a time . This is how it si done in
		mesh.macroel_sing_vrtx()
	"""
	arr_out = np.array(range(init + 1, init + n_vert_T4 + 1)).reshape((n_vert_T4/4, 4))
	arr_out	= np.concatenate((4*np.ones((n_vert_T4/4, 1)), arr_out), axis=1)
	with open (f_name_write, 'a') as tgt:
		np.savetxt(tgt, arr_out, fmt='%d')
	return

def vertices_by_elements (f_name, lang):
	""" 
	This function intends to be GENERAL. For any mesh given, not
	only our mesh. 

	input: text file like elements_by_vertices.txt
	
	return value: 

	output: vertices_by_elements.txt

	returns THE ELEMENTS touching each VERTEX IN GLOBAL ENUMERATION
	in the file

	vertices_by_elements.txt:
	--------------------------------------------------------
	
	+++ Each line is a vertex +++
	+++ contents of each line l: elements with vertex l +++

	e1|e2|e3|e4
	e1|e2|
	e1|e2|e3|
	e1|e2|e3|e4
	e1|

	etc

	In [2]: elem
	Out[2]: 
	['6 1 4 2 7 9 8\n',
	 '4 2 5 3 8\n',
	 '4 4 6 5 9\n',
	 '5 2 4 9 8 5\n',
	 '4 7 9 8 10\n']
	"""
	language  = {'C' : 0, 'Octave' : 1}
	elem_vert = {}
	with open (f_name, 'r') as elements:
		elem = elements.readlines()
	for j in range(len(elem)):
		e = [int(x) for x in elem[j].rstrip().split(' ')]
		for v_index in e[1:len(e)]:  ## con .get()
			if v_index in elem_vert: elem_vert[v_index].append(j+language[lang])
			else 				   : elem_vert[v_index] = [j+language[lang]] 
	with open ('vertices_by_elements.txt','w') as out:
		for vertex in elem_vert:
			out.write(str(vertex) + ' ' + ' '.join([str(r) for r in elem_vert[vertex]]) + '\n')
			#out.write(' '.join([str(r) for r in elem_vert[vertex]]) + '\n')
	return len(elem_vert)

def faces (f_name, n_elem, lang):
	""" 
	f_name == vertices_by_elements.txt

	n_elem == sum of write_elements_by_vertices (f,n) in this way it won't depend on the topology of the mesh.

	language[lang]: this is for the case starting with ones or zeroes

	First we make a dictionary with the table:  v | el1 el2 ... elk_v,
	that is, from the input f_name == vertices_by_elements.txt

	por ahora testeado en un macro--elemento tetra con \bv y \be.
	anda bien registrando las caras interelementales internas a Lambda_ell

	returns faces.txt as:
	--------------------------------------------------
	el | el | v | v | v | (v)
	1    4    2   4   8    9
	2    4    2   5   8
	3    4    4   5   9
	1    5    7   8   9
	"""
	language  = {'C' : 0, 'Octave' : 1}
	d = {}
	with open (f_name,'r') as inp:	
		lines = inp.readlines()
	n_lines = len(lines)
	
	for r in range(n_lines):
		chars = lines[r].rstrip().split(' ')
		d[chars[0]]  = dict(zip(range(len(chars)-1), [int(c) for c in chars[1:]]))
	## ver si se puede hacer mas lindo esto sin usar la lista values()
	
	s = ''
	for n in range(language[lang], n_elem + language[lang]):
		for m in range(language[lang], n):
			count = 0
			face_vertices = []
			for key, val in d.iteritems():
				if (m in val.values()) and (n in val.values()):
					count += 1
					face_vertices.append(str(int(key)))
			if count == 3:
				face_type = ' share triangle'
				s += str(m) + ' ' + str(n) + ' ' + ' '.join(face_vertices) + '\n'
			elif count == 4: 
				face_type = ' share rect'
				s += str(m) + ' ' + str(n) + ' ' + ' '.join(face_vertices) + '\n'
			if n==13 and m==2: 
				print('count')
				print(count)
	with open ('shared_faces.txt','w') as f:
		f.write(s)
	return d

def traverse (obj):
	""" This is something like a deep copy. """
	if isinstance(obj, dict):
		return {k: traverse(v) for k,v in obj.iteritems()}
	elif isinstance(obj, list):
		return [traverse(r) for r in obj]
	else:
		return obj

def face_enumeration (file_name):
	"""
	input:  elements_by_vertices.txt
	writes: faces_global.txt   
			faces_local_to_global.txt  
	"""
	face_index = 1
	global_string = ''
	local_to_global_string = ''
	with open (file_name, 'r') as data:
		elem_by_vert = data.readlines()
	n_el = len(elem_by_vert)
	for el in range(n_el):
		vertices = elem_by_vert[el].rstrip().split(' ')
		if vertices[0] == '6':
			global_string += '3 ' + vertices[1] + ' ' + vertices[2] + ' ' + vertices[3] + '\n'
			local_to_global_string += '2 ' + str(face_index) + ' ' 
			face_index += 1
			global_string += '3 ' + vertices[4] + ' ' + vertices[5] + ' ' + vertices[6] + '\n'
			local_to_global_string += str(face_index) + ' '
			face_index += 1
			global_string += '4 ' + vertices[2] + ' ' + vertices[1] + ' ' + vertices[5] + ' ' + vertices[4] + '\n'
			local_to_global_string += str(face_index) + ' '
			face_index += 1
			global_string += '4 ' + vertices[1] + ' ' + vertices[3] + ' ' + vertices[4] + ' ' + vertices[6] + '\n'
			local_to_global_string += str(face_index) + ' '
			face_index += 1
			global_string += '4 ' + vertices[3] + ' ' + vertices[2] + ' ' + vertices[6] + ' ' + vertices[5] + '\n'
			local_to_global_string += str(face_index) + '\n'
			face_index += 1
		elif vertices[0] == '5':
			global_string += '3 ' + vertices[1] + ' ' + vertices[2] + ' ' + vertices[5] + ' ' + '\n'
			local_to_global_string += '1 ' + str(face_index) + ' ' 
			face_index += 1
			global_string += '3 ' + vertices[4] + ' ' + vertices[1] + ' ' + vertices[5] + ' ' + '\n'
			local_to_global_string += str(face_index) + ' ' 
			face_index += 1
			global_string += '3 ' + vertices[3] + ' ' + vertices[4] + ' ' + vertices[5] + ' ' + '\n'
			local_to_global_string += str(face_index) + ' ' 
			face_index += 1
			global_string += '3 ' + vertices[2] + ' ' + vertices[3] + ' ' + vertices[5] + ' ' + '\n'
			local_to_global_string += str(face_index) + ' ' 
			face_index += 1
			global_string += '4 ' + vertices[1] + ' ' + vertices[2] + ' ' + vertices[4] + ' ' + vertices[3] + '\n'
			local_to_global_string += str(face_index) + '\n'
			face_index += 1
		elif vertices[0] == '4':
			global_string += '3 ' + vertices[1] + ' ' + vertices[2] + ' ' + vertices[3] + ' ' + '\n'
			local_to_global_string += '0 ' + str(face_index) + ' ' 
			face_index += 1
			global_string += '3 ' + vertices[1] + ' ' + vertices[2] + ' ' + vertices[4] + ' ' + '\n'
			local_to_global_string += str(face_index) + ' ' 
			face_index += 1
			global_string += '3 ' + vertices[1] + ' ' + vertices[3] + ' ' + vertices[4] + ' ' + '\n'
			local_to_global_string += str(face_index) + ' ' 
			face_index += 1
			global_string += '3 ' + vertices[2] + ' ' + vertices[3] + ' ' + vertices[4] + ' ' + '\n'
			local_to_global_string += str(face_index) + '\n'
			face_index += 1
		else: print('something went wrong')
	with open ('faces_global.txt', 'w') as out_global:
		out_global.write(global_string)
	with open ('faces_local_to_global.txt', 'w') as out_loc_to_global:
		out_loc_to_global.write(local_to_global_string)
	return

def kill_repeated (vertices_file_name):
	""" 
		searches: 		repeated vertices

		return value: 	dictionary of the vertices that have to be re-numbered in file
						elements_by_vertices.txt

		obs:            now we don't delete the entries on vertices.txt. We simply
						don't request them.
	""" 
	vertices = np.loadtxt(vertices_file_name)
	d_out 	 = {}
	Nv 		 = vertices.shape[0]

	for v in range(Nv):
		for w in range(v + 1, Nv):
			if np.all(np.equal(vertices[v], vertices[w])):
				d_out[v] = d_out.get(v,[])
				d_out[v].append(w)
	# remove_indices = []
	# for key in d_out:
	# 	remove_indices += d_out[key]
	# vertices_unique = np.delete(vertices, remove_indices, axis = 0)
	# with open ('vertices_unique.txt', 'w') as out_file:
	#     np.savetxt(out_file, vertices_unique, fmt = '%8.5f')
	return d_out

def kill_repeated_faces (faces_file_name):
	with open(faces_file_name,'r') as inp:
		faces = inp.readlines()

	face_dict = {}
	n_faces = len(faces)

	for f in range(n_faces):
		face_dict[f+1] = np.array([c for c in faces[f].rstrip().split(' ')],dtype=int)	

	d_out = {}
	for f in range(1, n_faces+1):
		for g in range(f + 1, n_faces+1):
			if (face_dict[f][0] == face_dict[g][0]):
				cmp1 = np.sort(face_dict[f][1:])
				cmp2 = np.sort(face_dict[g][1:])
				if np.all(np.equal(cmp1, cmp2)):
					d_out[f] = d_out.get(f, [])
					d_out[f].append(g)
	return d_out