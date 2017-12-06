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

for oc in xrange(2,9):    # integers 2, 3, 4 ,5 ,6, 7, 8
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

## look for repetitions
with open('elements_by_vertices_repeated.txt','r') as inp:
	things = inp.readlines()
a = mesh.np.array(things)
a = mesh.np.core.defchararray.replace(a,'\n','')

replace_verts = mesh_conectivity.kill_repeated('vertices.txt')  ## reads repetitions in vertices.txt and leaves a dictionary

for key in replace_verts:
	for r in replace_verts[key]:
		a = mesh.np.core.defchararray.replace(a, str(r + 1), str(key + 1))
with open('elements_by_vertices.txt','w') as out:						## <--- whithout repetitions
	mesh.np.savetxt(out, a, delimiter=" ",fmt='%s')

print('a')
print(a)

# at this point the program is already general:
# this 'inverts' the table of elements_by_vertices.txt
# writes on disc: vertices_by_elements.txt
mesh_conectivity.vertices_by_elements('elements_by_vertices.txt', 'Octave')


n_elem = mesh.np.shape(a)[0]
print('n_elem:')
print(n_elem)

# general function. 
# takes output of mesh_conectivity.vertices_by_elements()
# writes on disc: shared_faces.txt
d = mesh_conectivity.faces('vertices_by_elements.txt', n_elem, 'Octave')
# writes files: faces_global.txt 			----> with repetitions
#               faces_local_to_global.txt   ----> for macro--element of type 1
mesh_conectivity.face_enumeration('elements_by_vertices.txt')

with open('faces_local_to_global.txt', 'r') as ent:
	stuff = ent.readlines()
b = mesh.np.array(stuff)
b = mesh.np.core.defchararray.replace(b,'\n','')

replace_faces = mesh_conectivity.kill_repeated_faces('faces_global.txt')  ## reads repetitions in faces_global.txt and leaves a dictionary

for key in replace_faces:
	for r in replace_faces[key]:
		b = mesh.np.core.defchararray.replace(b, str(r), str(key))
with open('elements_by_faces.txt','w') as ex:                             ## <--- whithout repetitions
	mesh.np.savetxt(ex, b, delimiter=" ",fmt='%s')