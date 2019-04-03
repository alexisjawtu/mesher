    ###	TODO: refactor everything economizing code
from mesh import *
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import random

def plot_hybrid_macroel(plt_axes, vertices, n, local_mu = 1):
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
        plt_axes.plot(drawing[dr][0],drawing[dr][1],drawing[dr][2],color = "green")
    return

def plot_prism_macroel(plt_axes, vertices, n, n_vertical, local_mu = 1):
    local_grid_points = macroel_sing_edge(vertices, local_mu, n, n_vertical)
    for j in range(n_vertical+1):
        for k in xrange(n+1):
            plt_axes.plot(local_grid_points[j,k,0,0:n+1-k], local_grid_points[j,k,1,0:n+1-k], local_grid_points[j,k,2,0:n+1-k], color="red")
            plt_axes.plot(local_grid_points[j,0:n+1-k,0,k], local_grid_points[j,0:n+1-k,1,k], local_grid_points[j,0:n+1-k,2,k], color="green")
            x = np.array([local_grid_points[j,l,:,n-l-k] for l in range(n-k,-1,-1)])
            plt_axes.plot(x[:,0],x[:,1],x[:,2],color="black")
    
    for j in range(n+1):  # <<
        for i in xrange(n+1-j):  # verticals: TODO FIX PLOTTING BUG when setting independent horix and vertic. refinements
            plt_axes.plot(local_grid_points[:,j,0,i],local_grid_points[:,j,1,i],local_grid_points[:,j,2,i], color="blue")
    return local_grid_points

def plot_tetra_macroel(plt_axes, vertices, n, local_mu = 1):
    points_T5 = macroel_sing_vrtx(vertices[:,0], vertices[:,1], vertices[:,2], vertices[:,3], local_mu, n)
    for i in [4*nn for nn in range(points_T5.shape[0]/4)]:
        a = np.array([0,1,2,3]+[0,2,3,1])+i*np.ones(8,dtype=int)
        z = points_T5[a,np.array([2]*8)]
        y = points_T5[a,np.array([1]*8)]
        x = points_T5[a,np.array([0]*8)]
        plt_axes.plot(x, y, z,color="red")
    return

def plot_bbrick(mu = .65, angle_steps = [9], refinements = [3], vert_prims_refinements = 1):
    """ mu == [1,.4,.4,.4,.4] example for the graded case """

    ############ CONTINUE HERE: put variables for the repersentative lengths,
    ############ for example "prism_h" or the "2" between corners
    ############ and distances between singularities.
    permutation_of_vertices = np.array([[0,1,2,3],[3,1,2,0],[3,1,2,0],[0,1,2,3],[0,1,2,3]])
    elev    = 30
    #colors  = ['brown','darkgreen','red','black','fuchsia','blue']*7
    prism_h   = np.array([0,0,4])
    horiz1  = 2    
    y_max   = 1
    y_min   = -3
    y_int_min   = -2
    A0 = np.zeros(3)
    Q0 = np.array([0,0,-1])
    Q1 = np.array([0,1,-1])
    Q2 = np.array([-1,0,-1])
    R0 = Q1 - Q0 + Q2
    P1_hybrid_4 = Q2+Q1-2*Q0+A0

    vertices_hybrid_1   = np.array([Q0,Q2,Q1,A0])
    vertices_hybrid_2   = np.array([Q2-Q0+A0,P1_hybrid_4,Q2,A0]) 
    vertices_hybrid_3   = np.array([Q1+A0-Q0,Q1,P1_hybrid_4,A0])
    vertices_hybrid_4   = np.array([Q2+Q1-Q0,Q2,P1_hybrid_4,Q1]) # -----> oposite to a singular vertex
    vertices_tetra_1    = np.array([[0,-1,-1,0],[0,0,1,1],[0,-1,0,-1]])
    points_prisms   = np.array([Q0,Q1,Q2])
    points_prisms_1 = np.concatenate((points_prisms,points_prisms - prism_h)).transpose()
    points_prisms   = np.array([R0,Q1,Q2])
    points_prisms_2 = np.concatenate((points_prisms,points_prisms - prism_h)).transpose()
###########
    vertices_hybrid_11   = np.array([[0,-2,-1],[-1,-2,-1],[0,-3,-1],[0,-2,0]])
    vertices_tetra_2     = np.array([[0,-1,0,-1],[-2,-2,-3,-3],[0,-1,-1,0]])
    vertices_hybrid_12   = np.array([[0,-3,0],[0,-3,-1],[-1,-3,0],[0,-2,0]]) 
    vertices_hybrid_13   = np.array([[-1,-2,0],[-1,-3,0],[-1,-2,-1],[0,-2,0]])
    vertices_hybrid_14   = np.array([[-1,-3,-1],[-1,-3,0],[0,-3,-1],[-1,-2,-1]])
    points_prisms_3      = np.array([[-1,-3,-1],[0,-3,-1],[-1,-2,-1],[-1,-3,-1]-prism_h,[0,-3,-1]-prism_h,[-1,-2,-1]-prism_h]).transpose()
    points_prisms_4      = np.array([[0,-2,-1],[-1,-2,-1],[0,-3,-1],[0,-2,-1]-prism_h,[-1,-2,-1]-prism_h,[0,-3,-1]-prism_h]).transpose()
###########
    points_prisms_5      = np.array([[0,-2,-1],[0,-1,-1], [-1,-2,-1], [0,-2,-1]-prism_h,[0,-1,-1]-prism_h,[0,-2,-1]-prism_h]).transpose()
    points_prisms_6      = np.array([[0,-2,-1],[0,-3,-1],[1,-2,-1],[0,-2,-1]-prism_h,[0,-3,-1]-prism_h,[1,-2,-1]-prism_h]).transpose()

    points_prisms_7      = np.array([[-1,-1,-1],[-1,-2,-1],[0,-1,-1],[-1,-1,-1]-prism_h,[-1,-2,-1]-prism_h,[0,-1,-1]-prism_h]).transpose()
    points_prisms_8      = np.array([[1,y_min,-1],[1,y_int_min,-1],[0,y_min,-1],
                                     [1,y_min,-1]-prism_h,[1,y_int_min,-1]-prism_h,[0,y_min,-1]-prism_h]).transpose()



    for n in refinements:
        vert_prims_refinements = n
        fig = plt.figure()
        ax  = fig.add_subplot(1,1,1, projection='3d')
        for azim in angle_steps:
            ax.view_init(elev, 49 + 15*(azim-1))
            
            plot_hybrid_macroel(ax, vertices_hybrid_1, n, mu)
            plot_hybrid_macroel(ax, vertices_hybrid_2, n, mu)
            plot_hybrid_macroel(ax, vertices_hybrid_3, n, mu)
            plot_hybrid_macroel(ax, vertices_hybrid_4, n, 1)
            plot_tetra_macroel (ax, vertices_tetra_1, n, mu)
            
            plot_prism_macroel(ax, points_prisms_1,n,vert_prims_refinements, mu)
            plot_prism_macroel(ax, points_prisms_2,n,vert_prims_refinements, 1)  

            plot_hybrid_macroel(ax, vertices_hybrid_11, n, mu)
            plot_hybrid_macroel(ax, vertices_hybrid_12, n, mu)
            plot_hybrid_macroel(ax, vertices_hybrid_13, n, mu)
            plot_hybrid_macroel(ax, vertices_hybrid_14, n, 1)
            plot_tetra_macroel (ax, vertices_tetra_2, n, mu)

            plot_prism_macroel(ax,points_prisms_3,n,vert_prims_refinements, 1)
            plot_prism_macroel(ax,points_prisms_4,n,vert_prims_refinements, mu)

#########   non corner part
            plot_prism_macroel(ax, points_prisms_5,n,vert_prims_refinements, mu)
            plot_prism_macroel(ax, points_prisms_6,n,vert_prims_refinements, mu)
            plot_prism_macroel(ax, points_prisms_7,n,vert_prims_refinements, 1)
            plot_prism_macroel(ax, points_prisms_8,n,vert_prims_refinements,1)
#           ax.scatter(A0[0],A0[1],A0[2],color="black")
#           ax.scatter(0,-2,0,color="red")
#           ax.scatter(0,-2,-1,color="green")
#           ax.scatter(-1,-2,-1,color="blue")
#           ax.scatter(0,-3,-1,color="blue")
            #ax.scatter(P1_hybrid_4[0],P1_hybrid_4[1],P1_hybrid_4[2],color="blue")
            #ax.scatter((Q1+A0-Q0)[0],(Q1+A0-Q0)[1],(Q1+A0-Q0)[2],color="violet")

            ax.plot([],[],[],label = " ")
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
            ax.set_xlim3d(0,-3)
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