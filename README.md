## The Project
This software implements the meshing procedure proposed in chapter five of this work: http://cms.dm.uba.ar/academico/carreras/doctorado/thesisJawtuschenko.pdf.

### Basic tutorial
#### To build the anisotropic hybrid mesh for _the Fichera domain_
1. In an ipython3 shell:
   ```
   run main
   omega("experiments/fichera",levels=3)
   ```
2. That's it.

#### To _see_ the plotted mesh:
The present paragraph is for the case you comment 
the plotting line in main.omega() and want to know how to
plot after. Observe that now you have the **ten mesh files** in 
you working directory _experiments/_. This is the mesh, 
built as in the upper paragraph.
1. To visualize the mesh, in an ipython3 shell:
   ```
   run plot_lines.py
   plot("experiments/fichera.ver","experiments/fichera.ebv")`
   ```
2. That's it.

Repeat with any of the example domain files in _experiments/_, or with any initial 
partition you write for the polihedral domain you want, as long as you write your csv 
using the same convention as in the included examples, and with the levels of 
refinement you want.

#### The input file _fichera_
1. The rows are just the lists of vertices of the seven 
cubes whose union is the meshed domain. 
2. The leading number _3_
indicates each row is a graded hexahedral macroelement, one for
each octant in the euclidean three--dimensional space. 
3. The trailing number _0.65_ is the grading parameter 
defined within the theory of meshing procedures.
4. The same applies to the rest of input files.
