For plotting, ipython run plot_mesh.py

**instrucciones:**
------------------


macroelements.py:
----------------

For i in [0,1,2], D[i] is the dictionary of all the i-type macroelements:

Type 0 : hybrid
Type 1 : isotropic
Type 2 : prismatic

D[i][j] is a dictionary of lists with the coordinates of the vertices
of macro--element j of type i

main.py:
-------

* levels:  habr'a (levels + 1) ptos por arista de macro--el
* mu est'a en la variable mesh.mu_, o sea en el archivo mesh.py, linea 7, bajo 
el nombre mu_.

* si modificas un archivo que no es el que ejecutas con
'run' (untitled.py) entonces tenes que cerrar ipython y
volver abrirlo porque memoriza las importaciones.

* si modificas el archivo que ejecutas con 'run' (untitled.py), entonces anda 
dinamico, sin cerrar nada.

* C2 es cuadrante 2 (octante 2)
  En cada cuadrante vale que:
	  -T0 ... T3 son los macro-el con los tres poliedros y T4 es el regular que 
	  se subdivide solo con tetraedros.
	  -T0 es el macro-el uniforme opuesto al (0,0,0).
	  -T1 es el macro-el que toca al (0,0,0) con el tetraedro chiquito 
	  "de arriba de todo".

** en los archivos de la malla, para revisar junto con las fotos. **

* elements.txt: solo sirve por el dise~no del macro-el con los tres poliedros. 
Es algo ad-hoc y no hay que controlarlo.

* vertices.txt: cada linea es el numero (con repeticion) del vertice fisico que 
se lee hacia la derecha. De aca el programa lee solo las filas que necesita 
(o sea saltea repeticiones)

* elements_by_vertices_repeated.txt: 
	-col 1 == tipo de elem por cant de vert.
	-despues tiene los elementos enumeracion repetitiva de vertices

* elements_by_vertices.txt: el resultado de pasarle al anterior la funcion 
kill_repeated(). Son los elementos de toda la malla con la malla con la 
enumeracion global final con unicidad. El mayor indice de
vertice es mayor que la cantidad de vertices. 

* vertices_by_elements.txt: la inversa de la tabla anterior.

* shared_faces.txt: el archivo para el cual hacemo el anterior.
	-col1 y col2 == elem (a) y elem(b)
	-cols >= 3: la cara compartida por (a) y (b)

* faces_global.txt: todas las caras como van apareciendo, con repeticion, 
pero los vertices ya sin repeticion!
	-col 1 == tipo de poligono
	-cols > 1: vertices en algun recorrido.

* faces_local_to_global.txt: definicion de elementos por caras, caras con 
repeticion.
	-col1: tipo de elem que representa la presente linea.
		2== prisma, 1==piramde, 0== tetra.

* elements_by_faces.txt: el resultado de aplicarle a la tabla
anterior la funcion kill_repeated().
Para revisar elements_by_faces.txt sugiero tenerlo abierto al mismo tiempo 
que:
	-faces_local_to_global.txt
	-faces_global.txt

** Todo **
----------

The orientation of the non graded hybrid macrotetrahedron is irrelevant.
It remains to CHECK WITH EXPERIMENTS whether it is irrelevant to write the txt files