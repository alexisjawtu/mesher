## <one line to give the program's name and a brief idea of what it does.>
##     Copyright (C) 2018-2020  Alexis Boris Jawtuschenko.
## 
## This file is part of ?????.
## 
## ????? is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
## 
## ????? is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU General Public License
## along with ?????.  If not, see <https://www.gnu.org/licenses/>.

import numpy as np

def elements_by_vertices_hybrid (f_name_out, n, lang, initial):
    """ 
    initial: tracks the number of vertices of the previous macro_element

    Requires the file f_name_out.elem, written by
    mesh_write.write_element_indices(), to be in the same directory.

    n == current quantity of levels
    
    Writes GLOBAL INDICES per element appending in the file

    f_name_out in the following way:
    --------------------------------------------------------

    6|v1_{prism}|v2_{prism}|v3_{prism}|v4_{prism}|v5_{prism}|v6_{prism}|

    4|v1_{tetra}|v2_{tetra}|v3_{tetra}|v4_{tetra}|

    5|v1_{pyra} |v2_{pyra} |v3_{pyra} |v4_{pyra} |v5_{pyra} |

    etc...
    """
    language  = {'C' : 0, 'Octave' : 1}
    with open (f_name_out+".elem", 'r') as inp:
        indices = inp.readlines()
    with open (f_name_out+'.ebvr', 'ab') as out:
        for j in range(len(indices)):
            l = [int(c) for c in indices[j].rstrip().split(' ')]
            line = np.zeros((1,l[0]+1),dtype=int)
            line[0] = l[0]
            for i in range (1,l[0]+1):
                # (positions in the graph of the macro--element)
                h = l[3*(i-1)+1] ## height  
                d = l[3*(i-1)+2] ## depth   
                w = l[3*(i-1)+3] ## with    
                calc = int((h*(n**2+3*n+2) + h*(h-1)*(2*h-1)/6 - (2*n+3)*h*(h-1)/2) /2 + w + d*(n-h-d/2+3/2))
                line[0,i] =  initial + calc + language[lang]
            np.savetxt(out,line,fmt='%d')
    return len(indices)

def elements_by_vertices_tetra (f_name_out, levels, lang, init):
    """ TODO: we should write an algorithm that passes just one time per vertex
        recording every element that vertex belongs to in the correct order,
        then mark wich are the vertices in the frontier of macroelements and
        so reduce the complexity of the kill_repeated from cubic to quadratic.
 
    Writes rows of repeated indices, one row per element.
    levels**3   == number of elements. 
    levels**3*4 == number of vertices WITH REPETITIONS.
    We assume that in the (Npts x 3) array of points the tetrahedra appear in order taking the rows
    four at a time. This is how it is done in mesh.macroel_tetrahedra().
    """
    n_vert_repeated = levels**3*4
    arr_out = np.array(range(init + 1, init + n_vert_repeated + 1)).reshape((n_vert_repeated//4, 4))
    arr_out = np.concatenate((4*np.ones((n_vert_repeated//4, 1),dtype=int), arr_out), axis=1)
    with open (f_name_out+".ebvr", 'ab') as tgt:
        np.savetxt(tgt, arr_out, fmt='%d')
    return len(arr_out)

def elements_by_vertices_prisms (f_name_out, levels, lang, init):
    """ f_name_out: we append the elements represented
    by lists of six vertex indices. 
    
    The sum of the first k odd numbers is k**2, then
    levels**3 equals the number of elements of the present macro--element.

    nodes_per_layer: number of nodes in each layer 
    
    The algorithm implements six Node --> Element affine transforms
    according to this hexagon:
              
                 ---------
                . .  5  . .
               . 4 .   . 6 .
               ----- . -----
               . 1 .   . 3 .
                . .  2  . .
                 --------- 
    At the end of the function we add init to the whole np.array
    to have an inyective increasing enumeration.
    The return value is the number of elements in the prismatic macro--element
    """
    if levels==1:
        local_elmnts_by_vertices = np.array([6]+list(range(1,7))).reshape(1,6)
    else:    
        elems_per_level  = levels**2
        nodes_per_layer  = (levels+1)*(levels+2)//2    
        # Dictionary of 3-lists to append the first
        # nodes_per_layer nodes
        local_3_lists = {}
        for k in range(1,elems_per_level+1):
            local_3_lists[k] = []
        def assign (node, elements):
            for element in elements:
                local_3_lists[element]+=[node]
        # LOWEST ROW
        #  head base step
        assign(node=1, elements=[1])
        #  inductive middle steps:
        current_node = 2
        while current_node < levels+1:
            # 3 affine transforms
            left_above  = 2*current_node-3
            above       = left_above + 1
            right_above = above + 1
            assign(current_node, elements=[left_above,above,right_above]) 
            current_node += 1
        #  tail base step
        left_above = 2*levels-1
        assign(current_node, elements=[left_above])
        # INDUCTIVE MIDDLE ROWS
        row = levels #it also works as the odd sum limit
        while row > 2:
            extra_odd = 2*row-1
            # ROW HEAD BASE CASE
            current_node += 1
            below 	    = sum([2*k-1 for k in range(levels,row,-1)]) + 1
            right_below = below + 1
            right_above = below + extra_odd
            assign(current_node, [below,right_below,right_above])
            # INDUCTIVE MIDDLE STEPS. START HEXAGONS
            step = 1
            while step < row - 1:  # row 'row' has 'row' nodes
                current_node += 1
                left_below  = sum([2*k-1 for k in range(levels,row,-1)]) + step*2
                below       = left_below + 1
                right_below = below + 1
                left_above  = sum([2*k-1 for k in range(levels,row,-1)]) + extra_odd + (2*step-1)
                above       = left_above + 1
                right_above = above + 1
                assign(current_node, [left_below,below,right_below,left_above,above,right_above])
                step += 1
            # ROW TAIL BASE CASE
            current_node += 1
            below      = sum([2*k-1 for k in range(levels,row,-1)]) + extra_odd
            left_below = below - 1
            left_above = left_below + extra_odd - 1
            assign(current_node, [left_below,below,left_above])
            row -= 1
        # TWO UPPER ROWS
        # antepenultimate node
        assign(nodes_per_layer-2, [elems_per_level-3,elems_per_level-2,elems_per_level])
        # penultimate node
        assign(nodes_per_layer-1, [elems_per_level-2,elems_per_level-1,elems_per_level])
        # last node
        assign(nodes_per_layer, [elems_per_level])
        # INDEX REFLECTIONS
        local_elmnts_by_vertices = np.array(list(local_3_lists.values()))
        local_elmnts_by_vertices = np.hstack((np.array(list(local_3_lists.values())),
                                     local_elmnts_by_vertices + nodes_per_layer))
        for l in range(1,levels):
            local_elmnts_by_vertices = np.vstack((local_elmnts_by_vertices,
                                         local_elmnts_by_vertices[0:elems_per_level]
                                         +l*nodes_per_layer))
        local_elmnts_by_vertices = np.hstack((6*np.ones((levels**3,1)),
                                     local_elmnts_by_vertices+init))
    with open (f_name_out+".ebvr",'ab') as target:
        np.savetxt(target,local_elmnts_by_vertices.astype(int),fmt='%d')
    return levels**3 

def elements_by_vertices_hybridhexa (f_name_out, levels, lang, init):
    for i in range(4):
        init += elements_by_vertices_hybrid (f_name_out, levels, lang, init)
    elements_by_vertices_tetra (f_name_out, levels, lang, init)
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
	"""
	language  = {'C' : 0, 'Octave' : 1}
	elem_vert = {}
	with open (f_name+".elem", 'r') as elements:
		elem = elements.readlines()
	for j in range(len(elem)):
		e = [int(x) for x in elem[j].rstrip().split(' ')]
		for v_index in e[1:len(e)]:  ## con .get()
			if v_index in elem_vert: elem_vert[v_index].append(j+language[lang])
			else 				   : elem_vert[v_index] = [j+language[lang]] 
	with open (f_name+'.vbe','w') as out:
		for vertex in elem_vert:
			out.write(str(vertex) + ' ' + ' '.join([str(r) for r in elem_vert[vertex]]) + '\n')
			#out.write(' '.join([str(r) for r in elem_vert[vertex]]) + '\n')
	return len(elem_vert)

def faces (f_name, n_elem, lang):
	""" 
	returns shared_faces.txt as:
	--------------------------------------------------
	el | el | v | v | v | (v)
	1    4    2   4   8    9
	2    4    2   5   8
	3    4    4   5   9
	1    5    7   8   9

	f_name == vertices_by_elements.txt

	n_elem == sum of elements_by_vertices_hybrid (f,n) in this way it won't depend on the topology of the mesh.

	language[lang]: this is for the case starting with ones or zeroes

	First we make a dictionary with the table:  v | el1 el2 ... elk_v,
	that is, from the input f_name == vertices_by_elements.txt

	For now: tested against macroelement with both singularities
	"""
	language  = {'C' : 0, 'Octave' : 1}
	d = {}
	with open (f_name+".elem",'r') as inp:	
		lines = inp.readlines()
	n_lines = len(lines)
	
	for r in range(n_lines):
		chars = lines[r].rstrip().split(' ')
		d[chars[0]]  = dict(zip(range(len(chars)-1), [int(c) for c in chars[1:]]))
	## TODO: see if this can be done without the list values()
	
	s = ''
	for n in range(language[lang], n_elem + language[lang]):
		for m in range(language[lang], n):
			count = 0
			face_vertices = []
			for key, val in iter(d.items()):
				if (m in val.values()) and (n in val.values()):
					count += 1
					face_vertices.append(str(int(key)))
			if count == 3:
				s += str(m) + ' ' + str(n) + ' ' + ' '.join(face_vertices) + '\n'
			elif count == 4: 
				s += str(m) + ' ' + str(n) + ' ' + ' '.join(face_vertices) + '\n'
	with open ('shared_faces.txt','w') as f:
		f.write(s)
	return d

def face_enumeration (file_name):
	"""
	input:  elements_by_vertices.txt
	writes: faces_repeated.txt   
			faces_local_to_global.txt  
	"""
	def str_prism (vertices):
		ret  = '3 ' + vertices[1] + ' ' + vertices[2] + ' ' + vertices[3] + '\n'
		ret += '3 ' + vertices[4] + ' ' + vertices[5] + ' ' + vertices[6] + '\n'
		ret += '4 ' + vertices[2] + ' ' + vertices[1] + ' ' + vertices[5] + ' ' + vertices[4] + '\n'
		ret += '4 ' + vertices[1] + ' ' + vertices[3] + ' ' + vertices[4] + ' ' + vertices[6] + '\n'
		ret += '4 ' + vertices[3] + ' ' + vertices[2] + ' ' + vertices[6] + ' ' + vertices[5] + '\n'
		return ret

	def str_pyram (vertices):
		ret  = '3 ' + vertices[1] + ' ' + vertices[2] + ' ' + vertices[5] + ' ' + '\n'
		ret += '3 ' + vertices[4] + ' ' + vertices[1] + ' ' + vertices[5] + ' ' + '\n'
		ret += '3 ' + vertices[3] + ' ' + vertices[4] + ' ' + vertices[5] + ' ' + '\n'
		ret += '3 ' + vertices[2] + ' ' + vertices[3] + ' ' + vertices[5] + ' ' + '\n'
		ret += '4 ' + vertices[1] + ' ' + vertices[2] + ' ' + vertices[4] + ' ' + vertices[3] + '\n'
		return ret

	def str_tetra (vertices):
		ret  = '3 ' + vertices[1] + ' ' + vertices[2] + ' ' + vertices[3] + ' ' + '\n'
		ret += '3 ' + vertices[1] + ' ' + vertices[2] + ' ' + vertices[4] + ' ' + '\n'
		ret += '3 ' + vertices[1] + ' ' + vertices[3] + ' ' + vertices[4] + ' ' + '\n'
		ret += '3 ' + vertices[2] + ' ' + vertices[3] + ' ' + vertices[4] + ' ' + '\n'
		return ret

	writers 	= {'6' : str_prism, '5' : str_pyram, '4' : str_tetra}
	nr_faces 	= {6 : 5, 5 : 5, 4 : 4}
	face_index 	= 1
	
	global_string 			= ''
	local_to_global_string 	= ''

	with open (file_name+".ebv", 'r') as data:
		elem_by_vert = data.readlines()
	n_el = len(elem_by_vert)
	
	for el in range(n_el):
		vertices = elem_by_vert[el].rstrip().split(' ')
		nr_vert  = int(vertices[0])
		nr_face  = nr_faces[nr_vert]
		local_to_global_string	+= str(nr_vert - 4) + ' ' + ' '.join([str(face_index + m) for m in range(nr_face)]) + '\n'
		global_string 			+= writers[vertices[0]](vertices)
		face_index 				+= nr_face

	with open (file_name+'.faces_rep', 'w') as out_global:
		out_global.write(global_string)
	with open (file_name+'.fltg', 'w') as out_loc_to_global:
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
	counter = 1
	vertices = np.loadtxt(vertices_file_name)
	d_out 	 = {}
	Nv 		 = vertices.shape[0]
	for v in range(Nv):
		print ('progress: {0}/{1}\r'.format(counter,Nv), sep = ' ', end = '', flush=True)
		counter += 1
		for w in range(v + 1, Nv):
			if np.all(np.equal(vertices[v], vertices[w])):
				d_out[v] = d_out.get(v,[])
				d_out[v].append(w)
	print ('\r')
	return d_out

def kill_repeated_faces (faces_file_name):
	with open(faces_file_name+".faces_rep",'r') as inp:
		faces = inp.readlines()

	face_dict = {}
	n_faces = len(faces)

	for f in range(n_faces):
		face_dict[f+1] = np.array([c for c in faces[f].rstrip().split(' ')],dtype=int)	

	d_out = {}
	
	indices = list(range(1, n_faces+1))

	counter = 1
	for f in indices:
		print ('progress: {0}/{1}\r'.format(counter,n_faces), sep = ' ', end = '', flush=True)
		counter += 1
		for g in range(f + 1, n_faces+1):
			if (face_dict[f][0] == face_dict[g][0]):	
				if set(face_dict[f][1:]) == set(face_dict[g][1:]):
					d_out[f] = d_out.get(f, [])
					d_out[f].append(g)
	print ('\r')

	duplicates = []
	for ff in d_out:
		duplicates += d_out[ff]

	#duplicates = sorted(duplicates, reverse=True)
	duplicates = np.unique(duplicates)
	duplicates = np.fliplr([duplicates])[0]

	for i in duplicates:
		del faces[i-1]
		if i in indices:
			indices.remove(i)

	with open (faces_file_name+'.faces', 'w') as face_list:
		for x in range(len(faces)):
			face_list.write(faces[x])

	with open(faces_file_name+'.aidx', 'w') as idx_file:
		idx_file.write(str(indices))

	num_faces = len(faces)

	return d_out, indices, num_faces

def kill_repeated_edges():
    # TODO
    pass
