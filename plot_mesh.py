    ###	TODO: refactor everything economizing code
from mesh import *
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import random

def plot_hybrid_macroel(plt_axes, vertices, n, local_mu = 1):
    pass

def plot_prism_macroel(plt_axes, vertices, n, local_mu = 1):
    local_grid_points = macroel_sing_edge(vertices, local_mu, n)
    for j in range(n+1):
        for k in xrange(n+1):
            plt_axes.plot(local_grid_points[j,k,0,0:n+1-k], local_grid_points[j,k,1,0:n+1-k], local_grid_points[j,k,2,0:n+1-k], color="red")
            plt_axes.plot(local_grid_points[j,0:n+1-k,0,k], local_grid_points[j,0:n+1-k,1,k], local_grid_points[j,0:n+1-k,2,k], color="green")
            x = np.array([local_grid_points[j,l,:,n-l-k] for l in range(n-k,-1,-1)])
            plt_axes.plot(x[:,0],x[:,1],x[:,2], color="black")
        for i in xrange(n+1-j):
            plt_axes.plot(local_grid_points[:,j,0,i],local_grid_points[:,j,1,i],local_grid_points[:,j,2,i], color="blue")
    return

def plot_tetra_macroel(plt_axes, vertices, n, local_mu = 1):
    points_T5 = macroel_sing_vrtx(vertices[:,0], vertices[:,1], vertices[:,2], vertices[:,3], local_mu, n)
    for i in [4*nn for nn in range(points_T5.shape[0]/4)]:
        a = np.array([0,1,2,3]+[0,2,3,1])+i*np.ones(8,dtype=int)
        z = points_T5[a,np.array([2]*8)]
        y = points_T5[a,np.array([1]*8)]
        x = points_T5[a,np.array([0]*8)]
        plt_axes.plot(x, y, z)
    return

def plot_bbrick(mu = [1,1,1,1,.65], angle_steps = [9], refinements = [3]):
    """
    mu == [1,.4,.4,.4,.4] example for the graded case
    """
    macro_elems = [0,1,2,3]  # 0,1,2 or 3 in each cube
    octants     = [6] # range(6,9) # any sublist in range(2,9)
    permutation_of_vertices = np.array([[0,1,2,3],[3,1,2,0],[3,1,2,0],[0,1,2,3],[0,1,2,3]])
    for n in refinements:
        # here macro_elems is 3, just one hybrid
        coords  = cube_mesh_2(n,mu[3],p_,macro_el,octants,macro_elems)
        drawing = cube_drawing(coords,octants,macro_elems)
        del(coords)
        fig = plt.figure()
        ax  = fig.add_subplot(1,1,1, projection='3d')
        
        elev = 30
        
        #for azim in [49, 79,109,139]:
        colors  = ['brown','darkgreen','red','black','fuchsia','blue']*7
        #random.shuffle(colors)
        for azim in angle_steps:
            c   = 0
            # ax.axis('equal')
            # ax.set_xlim3d(0.1,-1.1)
            # ax.set_ylim3d(-1.0,1.0)
            # ax.set_zlim3d(-0.2,1.2)
            ax.view_init(elev, 49 + 15*(azim-1))
            for tetra in drawing:
                col_interval = int(len(tetra)/len(octants))
                for o in range(len(octants)):
                    col = colors[c]
                    c   = c + 1
                    for dr in range(col_interval*o, col_interval*(1+o)):
                        pass
                        ## TODO: this can be done putting (array_of_X, array_of_Y, array_of_Z, ...)
                        ## and not one by one as is now
                        ax.plot(tetra[dr][0],tetra[dr][1],tetra[dr][2],color = "green")

            # now cubic macro-els nr 0,1,2 and 4 with tetrahedra
            for oc in octants:
                for m in [4]:
                    q       = octant(oc, p_)
                    vertices = q[:,macro_el[m,permutation_of_vertices[m,:]]]
                    plot_tetra_macroel(ax,vertices,n,mu[m])
            
            # CONTINUE HERE:
            # figure out how to put universally the points in 'prism' for macroel_sing_edge()
            
            # now prismatic macroels
            trans           = np.array([0,0,3])
            Q0              = np.array([0,0,-1])
            Q1              = np.array([0,1,-1])
            Q2              = np.array([-1,0,-1])
            points_prisms   = np.array([Q0,Q1,Q2])
            points_prisms   = np.concatenate((points_prisms,points_prisms - trans)).transpose()
            plot_prism_macroel(ax, points_prisms, n, mu[m])  
            R0    = Q1 - Q0 + Q2
            points_prisms = np.array([R0,Q1,Q2])
            points_prisms = np.concatenate((points_prisms,points_prisms - trans)).transpose()
            plot_prism_macroel(ax, points_prisms, n, 1)  

            ax.plot([],[],[],label = "mu[3] = " + str(mu[3]) + " " + str(macro_elems) + " " + str(octants))
            legend = ax.legend()
            ax.set_xlabel(' X ')
            ax.set_ylabel(' Y ')
            ax.set_zlabel(' Z ')
            plt.show()
            fig.savefig('test-bbrick' + str(azim) + '-' + str(n) + '.png')
    return

def plot_fichera():
    mu          = .3
    macro_elems = [0,1,2,3]  # 0,1,2 or 3 in each cube
    angle_steps = range(1,2)
    refinements = range(1,2)
    octants     = range(2,9) # any sublist in range(2,9)
    
    for n in refinements:
        coords  = cube_mesh_2(n,mu,p_,macro_el,octants,macro_elems)
        # uncomment for just the second octant
        drawing = cube_drawing(coords,octants,macro_elems)
        # drawing = cube_drawing(coords)
        
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
            # ax.set_xlim3d(0.1,-1.1)
            # ax.set_ylim3d(-1.0,1.0)
            # ax.set_zlim3d(-0.2,1.2)
            angle = 49 + 15*(azim-1)
            ax.view_init(elev, angle)
            for tetra in drawing:
                col_interval = int(len(tetra)/len(octants))
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