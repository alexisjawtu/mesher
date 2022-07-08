import sys

import numpy as np

TODO:

sys.argv[1]
sys.argv[1] etc ... con el nombre del dir

# local = "6jul22_1"
local0 = "5jul22"
local = "7jul22"

to_save = "/home/alexis/work/mesher/experiments/wafer" + local
to_save0 = "/home/alexis/work/mesher/experiments/wafer" + local0
# to_save1 = "/home/alexis/work/mesher/experiments/wafer" + local1



elems = np.loadtxt("/home/alexis/work/CFEst_Mesh/CFEstMesh_1.0/bin/wafer" 
		+ constants.local 
		+ "/elements.dat", delimiter=",") + 1

nodes = np.loadtxt("/home/alexis/work/CFEst_Mesh/CFEstMesh_1.0/bin/wafer" 
		+ constants.local 
		+ "/nodes.dat", delimiter=",") + 1

_elems_ = np.hstack((4 * np.ones((elems.shape[0], 1)), elems)).astype(int)                                                  

np.savetxt(constants.to_save + "/elements.dat", _elems_, fmt='%d')
np.savetxt(constants.to_save + "/nodes.dat", nodes, delimiter=",")
