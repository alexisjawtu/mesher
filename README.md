## The Project
This software implements the meshing procedure proposed in chapter five of this work: http://cms.dm.uba.ar/academico/carreras/doctorado/thesisJawtuschenko.pdf.

### Basic tutorial
#### To build the anisotropic hybrid mesh
1. In an ipython3 shell:
   ```
   run main
   omega("experiments/fichera",levels=3)
   ```
2. That's it.

#### To _see_ the plotted mesh:
Now that have the **ten mesh files** in you working directory _experiments_ built in the 
upper paragraph. This is the mesh. To visualize the mesh, in an ipython3 shell:
   ```
   run plot_lines.py
   plot("experiments/fichera.ver","experiments/fichera.ebv")`
   ```
That's it.

Repeat with any of the example domain files in experiments, or with any initial partition you write for the polihedral domain you want, and with the levels of refinement you want.
