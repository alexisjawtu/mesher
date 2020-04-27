## The Project
This software implements the meshing procedure proposed in chapter five of this work: http://cms.dm.uba.ar/academico/carreras/doctorado/thesisJawtuschenko.pdf.

### Basic tutorial
1. In an ipython3 shell:
   ```
   run main
   omega("experiments/fichera",levels=3)
   ```
2. Now you have the **ten mesh files** in you working directory _experiments_. This is the mesh.
3. To _see_ the plotted mesh:
   ```
   run plot_lines.py
   plot("experiments/fichera.ver","experiments/fichera.ebv")`
   ```
4. Repeat with any of the example domain files in experiments, or with any initial partition you write for the polihedral domain you want, and with the levels of refinement you want.
