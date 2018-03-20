"""
searching repeated vertices:
 take from second column in file elem_by_vert.txt
 (because 6,5 y 4 are not things to repĺace)

the analogue for faces repetitions

another option is to put -3 -2 y -1 etc ...

Obs: write_elements_by_vertices_T4() has 4*ones etc...
has to be changed to: "-np.ones( ? ,dtype=int)" etc
"""

import mesh
import mesh_write
import mesh_conectivity

## macro_elements type 1: prisms, pyrs and tetra
## macro_elements type 2: prisms
## macro_elements type 3: tetra

levels = 2

## dictionary with the fichera mesh.
coords_ = mesh.cube_mesh_2(levels, mesh.mu_, mesh.p_, mesh.macro_el) 

# indices only for type I macro--el. writes file: elements.txt 
mesh_write.write_element_indices('elements.txt', levels)

init 	= 0

for oc in xrange(2,9):    # integers 2, 3, 4 ,5 ,6, 7, 8 represent octants
	for t in xrange(4):
		coords 	= coords_['points_T' + str(t) + '_C' + str(oc)] 
		mesh_conectivity.write_elements_by_vertices('elements.txt',2,'Octave',init) # writes elements_by_vertices_repeated.txt: GLOBAL INDICES per element
		init   += mesh_write.vertices(coords)	# writes 'vertices.txt' global list of vertices
	## Type II macro--element
	n_vertT4 = mesh.np.shape(coords_['points_T4_C' + str(oc)])[0]

	if not ((n_vertT4 % 4) == 0):
		print('incorrect number of points in T4')
		print('oc:' + str(oc) + ' t:' + str(t))
		exit()

	mesh_conectivity.write_elements_by_vertices_T4(n_vertT4, init ,'elements_by_vertices_repeated.txt')
	init += mesh_write.verticesT4(coords_['points_T4_C' + str(oc)], 'vertices.txt')

## look for repetitions of vertices
with open('elements_by_vertices_repeated.txt','r') as inp:
	things = inp.readlines()

n_elem = len(things)

elem_vert_repeated = mesh.np.zeros((n_elem,7),dtype=int)
for k in range(len(things)):
	ele 					= mesh.np.fromstring(things[k],dtype=int,sep=' ')
	elem_vert_repeated[k] 	= mesh.np.concatenate((ele,mesh.np.zeros((7-len(ele)),dtype=int)))

replace_verts = mesh_conectivity.kill_repeated('vertices.txt')  ## reads repetitions in vertices.txt and leaves a dictionary

source = mesh.np.copy(elem_vert_repeated).reshape(n_elem*7)
elem_vert_dscnt_indcs = mesh.np.copy(source)
for key in replace_verts:
	for r in replace_verts[key]:
		for l in range(len(elem_vert_dscnt_indcs)):
			if (elem_vert_dscnt_indcs[l] == r+1): elem_vert_dscnt_indcs[l] = key +1

elem_vert_dscnt_indcs = elem_vert_dscnt_indcs.reshape((n_elem,7))

del source
del elem_vert_repeated

## <--- whithout repetitions
with open('elements_by_vertices.txt','w') as out:
	for e in elem_vert_dscnt_indcs:
		mesh.np.savetxt(out, e[0:e[0]+1].reshape((1,e[0]+1)), fmt='%d')

# at this point the program is already general:
# this 'inverts' the table of elements_by_vertices.txt
# writes on disc: vertices_by_elements.txt
mesh_conectivity.vertices_by_elements('elements_by_vertices.txt', 'Octave')

# general function. 
# takes output of mesh_conectivity.vertices_by_elements()
# writes on disc: shared_faces.txt
d = mesh_conectivity.faces('vertices_by_elements.txt', n_elem, 'Octave')
# writes files: faces_repeated.txt 			----> with repetitions
#               faces_local_to_global.txt   ----> for macro--element of type 1
mesh_conectivity.face_enumeration('elements_by_vertices.txt')

with open('faces_local_to_global.txt', 'r') as ent:
	stuff = ent.readlines()

elem_faces_repeated = mesh.np.zeros((n_elem, 6),dtype=int)
for elem in range(len(stuff)):
	el 							= mesh.np.fromstring(stuff[elem],dtype=int,sep=' ')
	elem_faces_repeated[elem] 	= mesh.np.concatenate((el,mesh.np.zeros((6-len(el)),dtype=int)))

## reads repetitions in faces_repeated.txt and leaves a dictionary
## writes file: faces.txt                  -----> global unique enumeration of faces
replace_faces, indices, num_faces = mesh_conectivity.kill_repeated_faces('faces_repeated.txt')  

with open ('num_faces.txt','w') as n_of_faces:
	n_of_faces.write(str(num_faces))

## uniquifying faces:
source = mesh.np.copy(elem_faces_repeated).reshape(n_elem*6)
elem_faces_discnt_index = mesh.np.copy(source)
for key in replace_faces:
	for r in replace_faces[key]:
		for l in range(len(elem_faces_discnt_index)):
			if (elem_faces_discnt_index[l] == r):  elem_faces_discnt_index[l] = key

elem_faces_discnt_index = elem_faces_discnt_index.reshape((n_elem,6))

del elem_faces_repeated
del source

## bijection indices ----> [1,...,n_faces]
## unique, continuous indices with filling zeros
source = mesh.np.copy(elem_faces_discnt_index).reshape(n_elem*6)
elem_faces = mesh.np.copy(source)

for i in range(len(indices)):
	elem_faces[source==indices[i]] = i+1

elem_faces = mesh.np.copy(elem_faces).reshape((n_elem,6))
del elem_faces_discnt_index
del source

with open('elements_by_faces.txt','a') as ex:
	for elem in elem_faces:
		mesh.np.savetxt(ex, elem.reshape((1,6)),fmt='%d')
