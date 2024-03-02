import sys

import numpy as np


local 	= "jul22"

# to read
# adding 1 at the end because the connectivity was defined for Octave

elems 	= np.loadtxt("/home/alexis/work/CFEst_Mesh/CFEstMesh_10/bin/elements.dat", delimiter=",") + 1

# to save
to_save = "/home/alexis/work/mesher/experiments/wafer" + local

_elems_ = np.hstack((4 * np.ones((elems.shape[0], 1)), elems)).astype(int)                                                  

np.savetxt(to_save + "/elements.dat", _elems_, fmt='%d')
