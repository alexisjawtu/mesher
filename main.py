"""
1st: set levels (if levels == 1, then only macro--elements)

searching repeated vertices:
take from second column in file elem_by_vert.txt because 6,5 y 4 are not 
things to replace. the analogue for faces repetitions
another option is to put -3 -2 y -1 etc ...

Obs: write_elements_by_vertices_tetra() has 4*ones 
has to be changed to: -np.ones( whatever , dtype=int) etc
## macro_elements type 1: prisms, pyrs and tetra
## macro_elements type 2: prisms
## macro_elements type 3: tetra
"""

import numpy as np
import mesh
import mesh_write
import mesh_conectivity
import time

## TODO: como hacer el flujo de ejecucion. 
## Dejar para que ande con $ python module.py?
## Poner un txt con un diccionario con los refinamientos y el parametro de graduacion y otros,
## tipo los diccionarios para configurar el sublime?

#    Despues deberia andar tambien desde un csv hasta
#    dejar la malla en los txt como de fichera

local_meshers = { 0 : mesh.macroel_hybrid, 
				  1 : mesh.macroel_tetrahedra, 
				  2 : mesh.macroel_prisms }

def load_partition (in_file):
    with open(in_file,'r') as infile:
        inlist = infile.readlines()
    pre_list = [line.strip(' \n').split(',') for line in inlist]
    pre_list = [[int(st[0])] + [float(st[k]) for k in range(1,len(st)-1)]+[float(st[-1])] for st in pre_list]
    colors = [ "green", "red", "blue"]
	## TODO : reorder indices in this dictionary to be the same of file partition
    macro_elements = { key : 
                        { 0 : np.array(pre_list[key][1:-1]).reshape((len(pre_list[key])-2)//3,3), 
                          1 : pre_list[key][-1], 
                          2 : colors[pre_list[key][0]], 
                          3 : pre_list[key][0] } 
                       for key in range(len(pre_list)) }
    return macro_elements

def filter_repeated_vertices (n_vert_prism = 6):
	"""	look for repetitions of vertices
	###########################################################################
	at this point the program is already general:
	this 'inverts' the table of elements_by_vertices.txt
	writes on disc: vertices_by_elements.txt
	mesh_conectivity.vertices_by_elements('elements_by_vertices.txt', 'Octave')
	########################################################################"""
	with open('elements_by_vertices_repeated.txt','r') as inp:
		things = inp.readlines()
	
	n_elem = len(things)
	
	elem_vert_repeated = np.zeros((n_elem,n_vert_prism+1),dtype=int)
	for k in range(len(things)):
		ele 					= np.fromstring(things[k],dtype=int,sep=' ')
		elem_vert_repeated[k] 	= np.concatenate((ele,np.zeros((n_vert_prism+1-len(ele)),dtype=int)))
	
	print ('mesh_conectivity.kill_repeated()')
	t0 = time.time()
	replace_verts = mesh_conectivity.kill_repeated('vertices.txt')  ## repetitions in vertices.txt leave a dictionary
	print (time.time() - t0)
	print ('\n')
	
	print ('vertices replacement loop version 2')
	t0 		= time.time()
	counter = 1
	elem_vert_dscnt_indcs = np.copy(elem_vert_repeated).reshape(n_elem*7)
	for key in replace_verts:
		print ('progress: {0}/{1}\r'.format(counter,len(replace_verts))),
		counter += 1
		for r in replace_verts[key]:
			elem_vert_dscnt_indcs[elem_vert_dscnt_indcs == r+1] = (key +1)
	print ('\r')
	print (time.time() - t0)
	print ('\n')
	
	elem_vert_dscnt_indcs = elem_vert_dscnt_indcs.reshape((n_elem,7))
	
	del elem_vert_repeated
	
	## <--- whithout repetitions
	with open('elements_by_vertices.txt','w') as out:
		for e in elem_vert_dscnt_indcs:
			np.savetxt(out, e[0:e[0]+1].reshape((1,e[0]+1)), fmt='%d')

	return n_elem

def filter_repeated_faces (n_elem):
	""" general function. takes output of 
	mesh_conectivity.vertices_by_elements()	writes on disc: shared_faces.txt
	writes: faces_repeated.txt 			----> with repetitions
	        faces_local_to_global.txt   ----> for macro--element of type 1 """	
	print ('mesh_conectivity.face_enumeration')
	t0 = time.time()
	mesh_conectivity.face_enumeration('elements_by_vertices.txt')
	print (time.time() - t0)
	print ('\r')

	with open('faces_local_to_global.txt', 'r') as ent:
		stuff = ent.readlines()

	elem_faces_repeated = np.zeros((n_elem, 6),dtype=int)
	for elem in range(len(stuff)):
		el 							= np.fromstring(stuff[elem],dtype=int,sep=' ')
		elem_faces_repeated[elem] 	= np.concatenate((el,np.zeros((6-len(el)),dtype=int)))

	## reads repetitions in faces_repeated.txt and leaves a dictionary
	## writes file: faces.txt -----> global unique enumeration of faces
	print( 'mesh_conectivity.kill_repeated_faces')
	t0 = time.time()
	replace_faces, indices, num_faces = mesh_conectivity.kill_repeated_faces('faces_repeated.txt')  
	print( time.time() - t0)
	print( '\r')

	with open ('num_faces.txt','w') as n_of_faces:
		n_of_faces.write(str(num_faces))

	#### uniquifying faces:
	print( 'faces replacement loop version 2')
	t0 = time.time()
	elem_faces_discnt_index = np.copy(elem_faces_repeated).reshape(n_elem*6)
	counter = 1
	for key in replace_faces:
		print( 'progress: {0}/{1}\r'.format(counter,len(replace_faces))),
		counter += 1
		for r in replace_faces[key]:
			elem_faces_discnt_index[elem_faces_discnt_index==r] = key
	print ('\r')
	print (time.time() - t0)
	print ('\r')

	#elem_faces_discnt_index = elem_faces_discnt_index.reshape((n_elem,6))

	del elem_faces_repeated

	## bijection indices ----> [1,...,n_faces]
	## unique, continuous indices with filling zeros
	print ('loop for indexing faces with {1 ... n_faces}')
	t0 = time.time()
	elem_faces = np.copy(elem_faces_discnt_index)
	for i in range(len(indices)):
		elem_faces[elem_faces_discnt_index==indices[i]] = i+1
	print (time.time() - t0)
	print ('\r')

	elem_faces = np.copy(elem_faces).reshape((n_elem,6))
	del elem_faces_discnt_index

	with open('elements_by_faces.txt','a') as ex:
		for elem in elem_faces:
			np.savetxt(ex, elem.reshape((1,6)),fmt='%d')
	
	return 

def fichera (levels = 3, mu_ = .35, n_vert_prism = 6):
	print('levels: %s\r' % levels)
	## dictionary with the fichera mesh.
	print ('mesh.cube_mesh_2()')
	t0 = time.time()
	fichera_coords_ = mesh.cube_mesh_2(levels, mu_, mesh.p_, mesh.macro_el) 
	print (time.time() - t0)
	print ('\r')
	
	# indices only for type I macro--el. writes file: elements.txt 
	mesh_write.write_element_indices('elements.txt', levels)
	
	init 	= 0
	
	for oc in range(2,9):    # integers 2, 3, 4, 5, 6, 7, 8
		for t in range(4):
			coords 	= fichera_coords_['points_T' + str(t) + '_C' + str(oc)] 
			mesh_conectivity.write_elements_by_vertices_hybrid('elements.txt', levels, 'Octave', init) # writes elements_by_vertices_repeated.txt: GLOBAL INDICES per element
			init   += mesh_write.vertices(coords)	# writes 'vertices.txt' global list of vertices
		## Type II macro--element
		n_vertT4 = np.shape(fichera_coords_['points_T4_C' + str(oc)])[0]
	
		mesh_conectivity.write_elements_by_vertices_tetra(n_vertT4, init ,'elements_by_vertices_repeated.txt')
		init += mesh_write.vertices_macro_tetra(fichera_coords_['points_T4_C' + str(oc)], 'vertices.txt')

	filter_repeated_faces(filter_repeated_vertices())

def bBrick (levels = 3, n_vert_prism = 6, mu_ = .35):
	tau_zero = load_partition ("partition")

	for i, E in iter(tau_zero.items()):
		pass
		##local_meshers[E[3]]( CONTINUE HERE )


	return




####### SLOW (not used for now)
# print('mesh_conectivity.faces()')
# print(time.time())
# d = mesh_conectivity.faces('vertices_by_elements.txt', n_elem, 'Octave')
# print(time.time())
# print('\n')
##############################################################################

