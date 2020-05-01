#### TODOs
* elements_by_vertices_tetra() has 4\*ones 
has to be changed to: -np.ones( whatever , dtype=int) etc
macro_elements type 1: prisms, pyrs and tetra
macro_elements type 2: prisms
macro_elements type 3: tetra

##### Prioritary Items
0. ~~cubic macroelement 0: four tetrahedral and one hybrid macroel.~~
1. TEST the new cubic macroelement
2. document what is the input
2. Execution flow? $ python module.py?
3. Documentation about order and positions 
   of the vertices in each element? For example, in a pyramid is it the top always
   the _last_?
4. start overall documentation
5. choose a name for the project
6. joss?
8. cubic macroelement 1: tetrahedra within five tetrahedra 
9. cubic macroelement 2: two prisms
10. make tests ready to include in the installation
11. release under a license
  1. ~~Put a license notice in each file.~~
  2. ~~startup notice.~~
12. sourceforge  


about the elements_by_vertices_[type]() functions:
This funtion is at the beginning in the program, for the case
    we start with the mesh we proposed. For a general mesh, the algorithm
    starts directly in the next step (with the elements-by-vertices file 
    partition.ebv given somehow).