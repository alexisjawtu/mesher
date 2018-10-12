from mesh import *
import random

mu          = .65
macro_elems = [0,1,2,3]  # 0,1,2 or 3 in each cube
angle_steps = range(8,12)
refinements = [1,2,3]
octants     = range(5,6) # any sublist in range(2,9)

for n in refinements:
    coords  = cube_mesh_2(n,mu,p_,macro_el,octants,macro_elems)
    # uncomment for just the second octant
    drawing = cube_drawing(coords,octants,macro_elems)
    # drawing = cube_drawing(coords)
    
    fig = plt.figure()
    
    elev = 30
    
    #for azim in [49, 79,109,139]:
    colors  = ['brown','darkgreen','red','black','fuchsia','blue']*7
    #random.shuffle(colors)
    for azim in angle_steps:
        c   = 0
        ax  = fig.add_subplot(1,1,1, projection='3d')
        # ax.axis('equal')
    	# ax.set_xlim3d(0.1,-1.1)
    	# ax.set_ylim3d(-1.0,1.0)
    	# ax.set_zlim3d(-0.2,1.2)
        angle = 49 + 15*(azim-1)
        ax.view_init(elev, angle)
        print(drawing.shape)
        for tetra in drawing:
            col_interval = int(len(tetra)/len(octants))
            for o in range(len(octants)):
                col = colors[c]
                c   = c + 1
                for dr in range(col_interval*o, col_interval*(1+o)):
                    ax.plot(tetra[dr][0],tetra[dr][1],tetra[dr][2],color = col)

        ax.plot([],[],[],label = "mu = " + str(mu) + str(macro_elems) + str(octants))
        legend = ax.legend()
        fig.savefig('mesh' + str(azim) + '-' + str(n) + '.png')

