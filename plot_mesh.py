###	TODO: refactor everything economizing code
from mesh import *
from main import load_partition
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import random

def plot_hybrid_macroel(plt_axes, vertices, n, local_mu = 1, color_name = "green"):
    points = macroel_sing_vrtx_and_edge (vertices[0], vertices[1], vertices[2],vertices[3], local_mu, n)
    drawing = []
    n = points.shape[0] - 1
    for k in range (n+1):
        for i in range (n-k+1):
            for j in range (n-k-i+1):
                x, y, z = points[k,0:n-k-j+1,0,j], points[k,0:n-k-j+1,1,j], points[k,0:n-k-j+1,2,j]
                drawing.append([x,y,z])
            x, y, z = points[k,i,0,0:n-k-i+1], points[k,i,1,0:n-k-i+1], points[k,i,2,0:n-k-i+1]
            drawing.append([x,y,z])
        for c in range (n-k+1):
            x, y, z = np.zeros(c+1), np.zeros(c+1), np.zeros(c+1)
            for i in range (c+1):
                x[i] = points[k,i,0,c-i]
                y[i] = points[k,i,1,c-i]
                z[i] = points[k,i,2,c-i]
            drawing.append([x,y,z])  #  transversals
    for i in range (n):
        for j in range (n-i):
            stop = n - (i + j) + 1
            x = points[0:stop,i,0,j]
            y = points[0:stop,i,1,j]
            z = points[0:stop,i,2,j]
            drawing.append([x,y,z])  #  verticals
        x, y, z = np.zeros(n-i+1), np.zeros(n-i+1), np.zeros(n-i+1)
        for k in range (n-i+1):
            x[k] = points[k,i,0,n-k-i]
            y[k] = points[k,i,1,n-k-i]
            z[k] = points[k,i,2,n-k-i]
        drawing.append([x,y,z])  #  pyramidals
    # simply interchange the roles of i and j
        x, y, z = np.zeros(n-i+1), np.zeros(n-i+1), np.zeros(n-i+1)
        for k in range (n-i+1):
            x[k] = points[k,n-k-i,0,i]
            y[k] = points[k,n-k-i,1,i]
            z[k] = points[k,n-k-i,2,i]
        drawing.append([x,y,z])  #  pyramidals
    for dr in range(len(drawing)):
        ## TODO: this can be done putting (array_of_X, array_of_Y, array_of_Z, ...)
        ## and not one by one as is now
        plt_axes.plot(drawing[dr][0],drawing[dr][1],drawing[dr][2], color = color_name)
    return

def plot_prism_macroel(plt_axes, vertices, n, local_mu = 1, color_name = "blue"):
    n_vertical = n  ## TODO: FIX THIS
    local_grid_points = macroel_sing_edge(vertices, local_mu, n, n_vertical)
    for j in range(n_vertical+1):
        for k in range(n+1):
            plt_axes.plot(local_grid_points[j,k,0,0:n+1-k], local_grid_points[j,k,1,0:n+1-k], local_grid_points[j,k,2,0:n+1-k], color=color_name)
            plt_axes.plot(local_grid_points[j,0:n+1-k,0,k], local_grid_points[j,0:n+1-k,1,k], local_grid_points[j,0:n+1-k,2,k], color=color_name)
            x = np.array([local_grid_points[j,l,:,n-l-k] for l in range(n-k,-1,-1)])
            plt_axes.plot(x[:,0],x[:,1],x[:,2],color=color_name)
    
    for j in range(n+1):  # <<
        for i in range(n+1-j):  # verticals: TODO FIX PLOTTING BUG when setting independent horix and vertic. refinements
            plt_axes.plot(local_grid_points[:,j,0,i],local_grid_points[:,j,1,i],local_grid_points[:,j,2,i], color=color_name)
    return local_grid_points

def plot_tetra_macroel(plt_axes, vertices, n, local_mu = 1, color_name = "red"):
    points_T5 = macroel_sing_vrtx(vertices[0], vertices[1], vertices[2], vertices[3], local_mu, n)
    for i in [4*nn for nn in range(points_T5.shape[0]//4)]:
        a = np.array([0,1,2,3]+[0,2,3,1])+i*np.ones(8,dtype=int)
        z = points_T5[a,np.array([2]*8)]
        y = points_T5[a,np.array([1]*8)]
        x = points_T5[a,np.array([0]*8)]
        plt_axes.plot(x, y, z,color=color_name)
    return

# dictionary with the three previous plotting functs
plot_functions = { 0 : plot_hybrid_macroel,
                   1 : plot_tetra_macroel,
                   2 : plot_prism_macroel}

def plot(initial_partition = "partition.txt", angle_steps = [9], refinements = [3], vertical_prism_refinement = 1):
    """ initial_partition is a dictionary with the macroelements, that is, the
    first of the sequence of meshes. A record in initial_partition has to be:

    { 0 : np.array([P0,..,PN]), 1 : mu, 2 : color_string, 3 : type_of_macro_element}

    mu == [1,.4,.4,.4,.4] example for the graded case """
    macro_elements = load_partition (initial_partition)
    
    macro_elements[0][2] = "black"

    elev    = 90

    trans = {}

    for n in refinements:
        fig = plt.figure()
        ax  = fig.add_subplot(1,1,1,projection='3d')
        for azim in angle_steps:
            #ax.view_init(elev,49+15*(azim-1))
            ax.view_init(elev,0)
       

####### CONTINUE HERE: REVISAR ESTO!!! ojo que piso m por referencia
            for k, m in iter(macro_elements.items()):
                #if m[3] < 2:
                #    trans [k] = m
                #    trans [k][0] = m[0]*np.array([1,1,-1]) + np.array([0,0,-6])
                #if k in range(41,52):
                plot_functions[m[3]](ax, m[0], n, m[1], m[2])

            with open ('translations.txt','w') as tr__:
                tr__.write(trans.__str__())      

            m = macro_elements[0]
            plot_functions[m[3]](ax, m[0], n, m[1], m[2])

            ax.plot([],[],[],label = " ")
            legend = ax.legend()
            ax.set_xlabel(' X ')
            ax.set_ylabel(' Y ')
            ax.set_zlabel(' Z ')
#           fig.savefig('fichera-' + str(azim) + '-' + str(n) + str(mu) + '.png')
            plt.show()
    return

def plot_fichera():
    mu          = .3
    macro_elems = [0,1,2,3]  # 0,1,2 or 3 in each cube
    angle_steps = range(1,2)
    refinements = range(1,4)
    octants     = [2]#range(2,9) # any sublist in range(2,9)
    
    for n in refinements:
        coords  = cube_mesh_2(n,mu,p_,macro_el,octants,macro_elems)
        # uncomment for just the second octant
        drawing = cube_drawing(coords,octants,macro_elems)
        del(coords)
        fig = plt.figure()
        elev = 30
        
        #for azim in [49, 79,109,139]:
        colors  = ['brown','darkgreen','red','black','fuchsia','blue']*7
        #colors  = ['brown','darkgreen','red','black','fuchsia','blue']*7
        #random.shuffle(colors)
        for azim in angle_steps:
            c   = 0
            ax  = fig.add_subplot(1,1,1, projection='3d')
            plt.tight_layout()
            # ax.axis('equal')
            #ax.set_xlim3d(0,-3)
            # ax.set_ylim3d(-1.0,1.0)
            # ax.set_zlim3d(-0.2,1.2)
            ax.set_xlabel(' X ')
            ax.set_ylabel(' Y ')
            ax.set_zlabel(' Z ')
            angle = 49 + 15*(azim-1)
            ax.view_init(elev, angle)
            for tetra in drawing:
                col_interval = int(len(tetra)//len(octants))
                for o in range(len(octants)):
                    col = colors[c]
                    c   = c + 1
                    for dr in range(col_interval*o, col_interval*(1+o)):
                        ax.plot(tetra[dr][0],tetra[dr][1],tetra[dr][2],color = col)
    
            ax.plot([],[],[],label = "mu = " + str(mu) + str(macro_elems) + str(octants))
            legend = ax.legend()
            fig.savefig('fichera-' + str(azim) + '-' + str(n) + str(mu) + '.png')
            plt.show()
        plt.close(fig)