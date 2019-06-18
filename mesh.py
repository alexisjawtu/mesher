import numpy as np                              
import scipy.io as sio
# import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d import Axes3D

#mu_        = .35
p_          = np.array([[ 1,  1,  1,  1,  0,  0,  0,  0],   # octant 8
                        [-1,  0,  0, -1, -1,  0,  0, -1],   # x > 0, y < 0, z < 0
                        [-1, -1,  0,  0, -1, -1,  0,  0]])

### organization of the 5 tetrahedra resolving a cube: 
### the vertices of the macro--cube are: 0 ... 7.
macro_el    = np.array([[0,1,3,4],[2,3,1,6],[7,3,4,6],[5,4,1,6],[6,1,3,4]])

def n_elem_macrotetra (lev):
    """ sum([sum([(2*j+1) for j in range(l)]) for l in range(1,a+1)]) 
    this works just for the macro_element with both singularities """
    dic = {}
    dic['number_of_elements']   = lev*(lev+1)*(2*lev+1)//6
    dic['number_of_prisms']     = lev*(lev*(2*lev-3)+1)//6
    dic['number_of_tetr']       = lev*(lev+1)//2 
    dic['number_of_pyrs']       = lev*(lev-1)//2
    dic['number_of_vertices']   = sum([r*(r+1)//2 for r in range(lev+2)])
    return dic

def n_faces_macrotetra (lev):
    Nel = n_elem_macrotetra(lev)['number_of_elements']
    return 2*lev**2+lev*(lev+1)+Nel-lev**2+2*(lev-1)**2+(lev-1)*lev*(lev+1)//6

def octant (o, points):
    """ takes a fixed octant [points] and affine--transforms it to the other six.
        the last one is ones(3,1) to leave it unchanged """
    trans = np.array([[ 0,  0, -1, -1,  1,  1, -1, -1,  1],
                      [ 0,  0, -1,  1,  1, -1, -1,  1,  1],
                      [ 0,  0, -1, -1, -1,  1,  1,  1,  1]])
    return points*trans[:,o].reshape((3,1))

def lambda1 (i, j, k, n, mu):
    return float(i)/n * (float(i+j+k)/n)**((1/float(mu))-1)

def lambda2 (i, j, k, n, mu):
    return float(j)/n * (float(i+j+k)/n)**((1/float(mu))-1)

def lambda3 (i, j, k, n, mu):
    return float(k)/n * (float(i+j+k)/n)**((1/float(mu))-1)

#def macroel_tetrahedra (P0, P1, P2, P3, mu, n):
def macroel_tetrahedra (vertices, mu, n):
    """ 
        return value: Nel == nmbr of elmts
        
        OBS: the graduation is, by default, towards 'P0',
        so we have to put the corresponding permutation of vertices when 
        calling the function.
        
        TODO:  calculations of the form lambda_[0,i,j,k]*P0 could be done
        a priori and then call in the following manner:

            initialize an (4 x I x J x K x 3)-D array Lambda_ with

                    Lambda_[l,i,j,k,:] <--- lambda_[l,i,j,k]*P_l

            and then call it. Make sure the dimension is the most comfortable.

    P0 = np.array(P0).reshape((1,3))
    P1 = np.array(P1).reshape((1,3))
    P2 = np.array(P2).reshape((1,3))
    P3 = np.array(P3).reshape((1,3))

    """


    ## TODO: remove this reshape() thing
    P0 = np.array(vertices[:,0]).reshape((1,3))
    P1 = np.array(vertices[:,1]).reshape((1,3))
    P2 = np.array(vertices[:,2]).reshape((1,3))
    P3 = np.array(vertices[:,3]).reshape((1,3))


    lambda_ = np.zeros((4, n+1, n+1, n+1))

    points  = np.zeros((1,3))

## todo: fix the order of the indices

    Nel     = 0
    for k in range(n+1):
        for j in range(n+1):
            for i in range(n+1):
                lambda_[1,i,j,k] = lambda1(i,j,k,n,mu)
                lambda_[2,i,j,k] = lambda2(i,j,k,n,mu)
                lambda_[3,i,j,k] = lambda3(i,j,k,n,mu)
                lambda_[0,i,j,k] = 1 - lambda_[1,i,j,k] - lambda_[2,i,j,k] - lambda_[3,i,j,k]

    # now: element vertices
    for k in range(n):
        for j in range(n-k):
            for i in range(n-j-k):
                q0 = lambda_[0,i,j,k]  *P0 + lambda_[1,i,j,k]  *P1 + lambda_[2,i,j,k]*P2 + lambda_[3,i,j,k]*P3
                q1 = lambda_[0,i+1,j,k]*P0 + lambda_[1,i+1,j,k]*P1 + lambda_[2,i+1,j,k]*P2 + lambda_[3,i+1,j,k]*P3
                q2 = lambda_[0,i,j+1,k]*P0 + lambda_[1,i,j+1,k]*P1 + lambda_[2,i,j+1,k]*P2 + lambda_[3,i,j+1,k]*P3
                q3 = lambda_[0,i,j,k+1]*P0 + lambda_[1,i,j,k+1]*P1 + lambda_[2,i,j,k+1]*P2 + lambda_[3,i,j,k+1]*P3

                points = np.concatenate((points,q0,q1,q2,q3))
                del(q0,q1,q2,q3)
                Nel += 1

    for k in range(n-1):
        for j in range(n-1-k):
            for i in range(n-1-j-k):
                q0 = lambda_[0,i+1,j,k]*P0 + lambda_[1,i+1,j,k]*P1 + lambda_[2,i+1,j,k]*P2 + lambda_[3,i+1,j,k]*P3
                q1 = lambda_[0,i,j+1,k]*P0 + lambda_[1,i,j+1,k]*P1 + lambda_[2,i,j+1,k]*P2 + lambda_[3,i,j+1,k]*P3
                q2 = lambda_[0,i,j,k+1]*P0 + lambda_[1,i,j,k+1]*P1 + lambda_[2,i,j,k+1]*P2 + lambda_[3,i,j,k+1]*P3
                q3 = lambda_[0,i+1,j,k+1]*P0 + lambda_[1,i+1,j,k+1]*P1 + lambda_[2,i+1,j,k+1]*P2 + lambda_[3,i+1,j,k+1]*P3

                points = np.concatenate((points,q0,q1,q2,q3))
                Nel += 1

                q0 = lambda_[0,i,j+1,k]*P0 + lambda_[1,i,j+1,k]*P1 + lambda_[2,i,j+1,k]*P2 + lambda_[3,i,j+1,k]*P3
                q1 = lambda_[0,i,j,k+1]*P0 + lambda_[1,i,j,k+1]*P1 + lambda_[2,i,j,k+1]*P2 + lambda_[3,i,j,k+1]*P3
                q2 = lambda_[0,i+1,j,k+1]*P0 + lambda_[1,i+1,j,k+1]*P1 + lambda_[2,i+1,j,k+1]*P2 + lambda_[3,i+1,j,k+1]*P3
                q3 = lambda_[0,i,j+1,k+1]*P0 + lambda_[1,i,j+1,k+1]*P1 + lambda_[2,i,j+1,k+1]*P2 + lambda_[3,i,j+1,k+1]*P3

                points = np.concatenate((points,q0,q1,q2,q3))
                Nel += 1
                
                q0 = lambda_[0,i+1,j,k]*P0 + lambda_[1,i+1,j,k]*P1 + lambda_[2,i+1,j,k]*P2 + lambda_[3,i+1,j,k]*P3
                q1 = lambda_[0,i,j+1,k]*P0 + lambda_[1,i,j+1,k]*P1 + lambda_[2,i,j+1,k]*P2 + lambda_[3,i,j+1,k]*P3
                q2 = lambda_[0,i+1,j+1,k]*P0 + lambda_[1,i+1,j+1,k]*P1 + lambda_[2,i+1,j+1,k]*P2 + lambda_[3,i+1,j+1,k]*P3
                q3 = lambda_[0,i+1,j,k+1]*P0 + lambda_[1,i+1,j,k+1]*P1 + lambda_[2,i+1,j,k+1]*P2 + lambda_[3,i+1,j,k+1]*P3

                points = np.concatenate((points,q0,q1,q2,q3))
                Nel += 1
                
                q0 = lambda_[0,i,j+1,k]*P0 + lambda_[1,i,j+1,k]*P1 + lambda_[2,i,j+1,k]*P2 + lambda_[3,i,j+1,k]*P3
                q1 = lambda_[0,i+1,j+1,k]*P0 + lambda_[1,i+1,j+1,k]*P1 + lambda_[2,i+1,j+1,k]*P2 + lambda_[3,i+1,j+1,k]*P3
                q2 = lambda_[0,i+1,j,k+1]*P0 + lambda_[1,i+1,j,k+1]*P1 + lambda_[2,i+1,j,k+1]*P2 + lambda_[3,i+1,j,k+1]*P3
                q3 = lambda_[0,i,j+1,k+1]*P0 + lambda_[1,i,j+1,k+1]*P1 + lambda_[2,i,j+1,k+1]*P2 + lambda_[3,i,j+1,k+1]*P3
                
                points = np.concatenate((points,q0,q1,q2,q3))
                Nel += 1

    for k in range(n-2):
        for j in range(n-2-k):
            for i in range(n-2-j-k):
                q0 = lambda_[0,i+1,j+1,k]*P0 + lambda_[1,i+1,j+1,k]*P1 + lambda_[2,i+1,j+1,k]*P2 + lambda_[3,i+1,j+1,k]*P3
                q1 = lambda_[0,i+1,j,k+1]*P0 + lambda_[1,i+1,j,k+1]*P1 + lambda_[2,i+1,j,k+1]*P2 + lambda_[3,i+1,j,k+1]*P3
                q2 = lambda_[0,i,j+1,k+1]*P0 + lambda_[1,i,j+1,k+1]*P1 + lambda_[2,i,j+1,k+1]*P2 + lambda_[3,i,j+1,k+1]*P3
                q3 = lambda_[0,i+1,j+1,k+1]*P0 + lambda_[1,i+1,j+1,k+1]*P1 + lambda_[2,i+1,j+1,k+1]*P2 + lambda_[3,i+1,j+1,k+1]*P3

                points = np.concatenate((points,q0,q1,q2,q3))
                Nel += 1

    points = np.delete(points, 0, 0)
    return points

#def macroel_hybrid (local_origin, base_vrtx_1, base_vrtx_2, sing_vrtx, mu, n):
def macroel_hybrid (macroel_vertices, mu, n):
    """ 
        vertices = ( local_origin | base_vrtx_1 | base_vrtx_2 | sing_vrtx )
        singular edge the one which is parallel to v = (sing_vrtx - local_origin)
        n  == levels
        mu == grading param
    """
    points = np.zeros((n+1, n+1, 3, n+1))  # level, i, coordinate, j
    for k in range(n+1):
        for i in range(n-k+1):
            for j in range(n-k-i+1):
                # it can be done with much lesser calls to these lambda
                coef = (lambda1 (i,j,0,n,mu), lambda2 (i,j,0,n,mu))
                # the sub-mesh is just the following two lines
                points[k,i,:,j] = coef[0]*(macroel_vertices[:,1]-macroel_vertices[:,0]) + coef[1]*(macroel_vertices[:,2]-macroel_vertices[:,0])
                points[k,i,:,j] += (1-(float(n-k)/n)**(1/mu))*(macroel_vertices[:,3]-macroel_vertices[:,0]) + macroel_vertices[:,0]
    return points

def macroel_prisms (macroel_vertices, mu, n):
    """ 
    n := number of subintervals betweeen nodes of each horizontal edge of the
    macroelement
    n_vertical := number of subintervals betweeen nodes of each vertical edge of the
    macroelement. Warning: my Thesis states that n_vertical == n. Otherwise it's not proved.

    TODO: In a future version we will have n_vertical independent of n

    order of the columns in M := macroel_vertices:
        (M[0], M[1], M[2]) == a triangle
        (M[3], M[4], M[5]) == the other triangle
    
        edges perpendicular to the triangles:
            edge in the singular direction == [M[0],M[3]]
            the others:                       [M[1],M[4]], [M[2],M[5]]

    obs: to sum faster, the levels in the points array which are 
    greater than zero end up filled as a rectangle of points. 
    Be careful not to use that coordinates. 
    """
    n_vertical = n 
    points = np.zeros((n_vertical+1,n+1,3,n+1))
    for y in range(n+1):
        for z in range(n+1-y):
            lambda_1, lambda_2 = lambda1 (y,z,0,n,mu), lambda2 (y,z,0,n,mu)
            temp = lambda_1*(macroel_vertices[:,1] - macroel_vertices[:,0]) + lambda_2*(macroel_vertices[:,2] - macroel_vertices[:,0])
            #points[0,y,:,z] = temp + macroel_vertices[:,0]
            points[0,y,:,z] = temp + macroel_vertices[:,0]

    for x in range(1,n_vertical+1): # translating level 0 to the levels above
    	points[x,:,:,:] = points[0,:,:,:] + (float(x)/n_vertical)*(macroel_vertices[:,3] - macroel_vertices[:,0]).reshape((3,1))
    return points

def line (x, y, z):
    """ (x[0], y[0], z[0]) -- ... -- (x[n-1], y[n-1], z[n-1]) """
    s = '\n\t\\draw '
    l = len(x)
    for i in range (l - 1):
        s += '('+str(x[i])[0:6]+','+str(y[i])[0:6]+','+str(z[i])[0:6]+')'+' -- '
    s += '('+str(x[l-1])[0:6]+','+str(y[l-1])[0:6]+','+str(z[l-1])[0:6]+');'
    return s

def cube2tex (obj, name = 'cube.tex'):
    s = ''
    for tetra in obj:
        for dr in tetra:
            s += line(dr[0],dr[1],dr[2])
    with open (name,'w') as d:
        d.write('\\documentclass{article}\n\\usepackage{tikz}\n\\begin{document}\n')
        d.write('\\begin{tikzpicture}[scale=3]\n')
        d.write(s)
        d.write('\n\n\\end{tikzpicture}\n\\end{document}')
    return

def cube2mat (obj, file_name = 'data.mat'):
    ## to draw in Octave with cubo.m
    d = {}  
    i = -1
    for tetra in obj:
        for dr in tetra:
            i += 1
            d['dr'+str(i)] = [dr[0],dr[1],dr[2]]
    sio.savemat(file_name, d)
    return i

def cube_mesh_2 (n, mu, p, tetrahedra, octants = range(2,9), macro_elems = [0,1,2,3,4]):
    """ TODO: 
        
        this has to be an independent "hexaedral" macroelement. see point 6)
        in the notebook

        n == levels
        mu == grading param
        p == vertices of the cubes; the octants
        tetrahedra ==  
    """
    mu_vec = [1, mu, mu, mu]            # first one is not graded!
    dict_save = {}
    for o in octants:
        q = octant(o, p)
        ## TODO: fix these ugly two FORs
        for m in [z for z in macro_elems if z < 4]:
            dict_save['points_T'+str(m)+'_C'+str(o)] = macroel_hybrid (q[:,tetrahedra[m,0:4]],mu_vec[m],n)
        for m in [z for z in macro_elems if z == 4]: # CALCULATION FOR t = 4 (T5)
            dict_save['points_tetra_C' + str(o)] = macroel_tetrahedra(q[:,tetrahedra[4,0:4]], mu, n)
    return dict_save

################################################################################
################################################################################
### not likely to be used again:

#def plot_mesh (obj, elev = 30, azim = 45, colors = ['darkgreen']+['black']+['fuchsia']+['blue']):
#   fig = plt.figure()
#   ax  = fig.add_subplot(1,1,1, projection='3d')
#   ax.axis('equal')
#   ax.view_init(elev, azim)
#   #  fix this!!
#   ## plt.title(str(n) + ' levels.')
#   colors = iter(colors)
#   for tetra in obj:
#       col = next(colors)#for te in tetra:
#       for dr in tetra:
#           ax.plot(dr[0],dr[1],dr[2],color = col)
#   plt.show()
#   return fig
################################################################################
################################################################################
