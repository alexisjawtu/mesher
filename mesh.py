import numpy as np                              

def lambda1 (i, j, k, n, mu):
    return float(i)/n * (float(i+j+k)/n)**((1/float(mu))-1)

def lambda2 (i, j, k, n, mu):
    return float(j)/n * (float(i+j+k)/n)**((1/float(mu))-1)

def lambda3 (i, j, k, n, mu):
    return float(k)/n * (float(i+j+k)/n)**((1/float(mu))-1)

def macroel_tetrahedra (vertices, mu, n):
    """ 
        n == number of levels, i.e. n + 1 == number of nodes per edge of 
        the macro--element.

        OBS: the graduation is, by default, towards 'P0',
        so we have to put the corresponding permutation of vertices when 
        calling the function.
        
        TODO:  calculations of the form lambda_[0,i,j,k]*P0 could be done
        a priori and then call in the following manner:

            initialize an (4 x I x J x K x 3)-D array Lambda_ with

                    Lambda_[l,i,j,k,:] <--- lambda_[l,i,j,k]*P_l

            and then call it. Make sure the dimension is the most comfortable.
    """
    vertices = vertices.transpose()

    P0 = vertices[0,:].reshape(1,3)
    P1 = vertices[1,:].reshape(1,3)
    P2 = vertices[2,:].reshape(1,3)
    P3 = vertices[3,:].reshape(1,3)

    lambda_ = np.zeros((4, n+1, n+1, n+1))

    points  = np.zeros((1,3))

## todo: fix the order of the indices

    for k in range(n+1):
        for j in range(n+1):
            for i in range(n+1):
                lambda_[1,i,j,k] = lambda1(i,j,k,n,mu)
                lambda_[2,i,j,k] = lambda2(i,j,k,n,mu)
                lambda_[3,i,j,k] = lambda3(i,j,k,n,mu)
                lambda_[0,i,j,k] = 1 - lambda_[1,i,j,k] \
                                   - lambda_[2,i,j,k] - lambda_[3,i,j,k]

    # now: element vertices
    for k in range(n):
        for j in range(n-k):
            for i in range(n-j-k):
                q0 = lambda_[0,i,j,k]  *P0 + lambda_[1,i,j,k]  *P1 \
                     + lambda_[2,i,j,k]*P2 + lambda_[3,i,j,k]*P3
                q1 = lambda_[0,i+1,j,k]*P0 + lambda_[1,i+1,j,k]*P1 \
                     + lambda_[2,i+1,j,k]*P2 + lambda_[3,i+1,j,k]*P3
                q2 = lambda_[0,i,j+1,k]*P0 + lambda_[1,i,j+1,k]*P1 \
                     + lambda_[2,i,j+1,k]*P2 + lambda_[3,i,j+1,k]*P3
                q3 = lambda_[0,i,j,k+1]*P0 + lambda_[1,i,j,k+1]*P1 \
                     + lambda_[2,i,j,k+1]*P2 + lambda_[3,i,j,k+1]*P3

                points = np.concatenate((points,q0,q1,q2,q3))
                del(q0,q1,q2,q3)

    for k in range(n-1):
        for j in range(n-1-k):
            for i in range(n-1-j-k):
                q0 = lambda_[0,i+1,j,k]*P0 + lambda_[1,i+1,j,k]*P1 \
                     + lambda_[2,i+1,j,k]*P2 + lambda_[3,i+1,j,k]*P3
                q1 = lambda_[0,i,j+1,k]*P0 + lambda_[1,i,j+1,k]*P1 \
                     + lambda_[2,i,j+1,k]*P2 + lambda_[3,i,j+1,k]*P3
                q2 = lambda_[0,i,j,k+1]*P0 + lambda_[1,i,j,k+1]*P1 \
                     + lambda_[2,i,j,k+1]*P2 + lambda_[3,i,j,k+1]*P3
                q3 = lambda_[0,i+1,j,k+1]*P0 + lambda_[1,i+1,j,k+1]*P1 \
                     + lambda_[2,i+1,j,k+1]*P2 + lambda_[3,i+1,j,k+1]*P3

                points = np.concatenate((points,q0,q1,q2,q3))

                q0 = lambda_[0,i,j+1,k]*P0 + lambda_[1,i,j+1,k]*P1 \
                     + lambda_[2,i,j+1,k]*P2 + lambda_[3,i,j+1,k]*P3
                q1 = lambda_[0,i,j,k+1]*P0 + lambda_[1,i,j,k+1]*P1 \
                     + lambda_[2,i,j,k+1]*P2 + lambda_[3,i,j,k+1]*P3
                q2 = lambda_[0,i+1,j,k+1]*P0 + lambda_[1,i+1,j,k+1]*P1 \
                     + lambda_[2,i+1,j,k+1]*P2 + lambda_[3,i+1,j,k+1]*P3
                q3 = lambda_[0,i,j+1,k+1]*P0 + lambda_[1,i,j+1,k+1]*P1 \
                     + lambda_[2,i,j+1,k+1]*P2 + lambda_[3,i,j+1,k+1]*P3

                points = np.concatenate((points,q0,q1,q2,q3))
                
                q0 = lambda_[0,i+1,j,k]*P0 + lambda_[1,i+1,j,k]*P1 \
                     + lambda_[2,i+1,j,k]*P2 + lambda_[3,i+1,j,k]*P3
                q1 = lambda_[0,i,j+1,k]*P0 + lambda_[1,i,j+1,k]*P1 \
                     + lambda_[2,i,j+1,k]*P2 + lambda_[3,i,j+1,k]*P3
                q2 = lambda_[0,i+1,j+1,k]*P0 + lambda_[1,i+1,j+1,k]*P1 \
                     + lambda_[2,i+1,j+1,k]*P2 + lambda_[3,i+1,j+1,k]*P3
                q3 = lambda_[0,i+1,j,k+1]*P0 + lambda_[1,i+1,j,k+1]*P1 \
                     + lambda_[2,i+1,j,k+1]*P2 + lambda_[3,i+1,j,k+1]*P3

                points = np.concatenate((points,q0,q1,q2,q3))
                
                q0 = lambda_[0,i,j+1,k]*P0 + lambda_[1,i,j+1,k]*P1 \
                     + lambda_[2,i,j+1,k]*P2 + lambda_[3,i,j+1,k]*P3
                q1 = lambda_[0,i+1,j+1,k]*P0 + lambda_[1,i+1,j+1,k]*P1 \
                     + lambda_[2,i+1,j+1,k]*P2 + lambda_[3,i+1,j+1,k]*P3
                q2 = lambda_[0,i+1,j,k+1]*P0 + lambda_[1,i+1,j,k+1]*P1 \
                     + lambda_[2,i+1,j,k+1]*P2 + lambda_[3,i+1,j,k+1]*P3
                q3 = lambda_[0,i,j+1,k+1]*P0 + lambda_[1,i,j+1,k+1]*P1 \
                     + lambda_[2,i,j+1,k+1]*P2 + lambda_[3,i,j+1,k+1]*P3
                
                points = np.concatenate((points,q0,q1,q2,q3))

    for k in range(n-2):
        for j in range(n-2-k):
            for i in range(n-2-j-k):
                q0 = lambda_[0,i+1,j+1,k]*P0 + lambda_[1,i+1,j+1,k]*P1 \
                     + lambda_[2,i+1,j+1,k]*P2 + lambda_[3,i+1,j+1,k]*P3
                q1 = lambda_[0,i+1,j,k+1]*P0 + lambda_[1,i+1,j,k+1]*P1 \
                     + lambda_[2,i+1,j,k+1]*P2 + lambda_[3,i+1,j,k+1]*P3
                q2 = lambda_[0,i,j+1,k+1]*P0 + lambda_[1,i,j+1,k+1]*P1 \
                     + lambda_[2,i,j+1,k+1]*P2 + lambda_[3,i,j+1,k+1]*P3
                q3 = lambda_[0,i+1,j+1,k+1]*P0 + lambda_[1,i+1,j+1,k+1]*P1 \
                     + lambda_[2,i+1,j+1,k+1]*P2 + lambda_[3,i+1,j+1,k+1]*P3

                points = np.concatenate((points,q0,q1,q2,q3))

    points = np.delete(points, 0, 0)
    return points

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
                points[k,i,:,j] = coef[0]*(macroel_vertices[:,1]-macroel_vertices[:,0]) \
                                  + coef[1]*(macroel_vertices[:,2]-macroel_vertices[:,0])
                points[k,i,:,j] += (1-(float(n-k)/n)**(1/mu))*(macroel_vertices[:,3]-macroel_vertices[:,0]) \
                                  + macroel_vertices[:,0]
    return points

def macroel_prisms (macroel_vertices, mu, levels):
    """ 
    levels := number of subintervals betweeen nodes of each horizontal edge of the
    macroelement
    n_vertical := number of subintervals betweeen nodes of each vertical edge of the
    macroelement. Warning: my Thesis states that n_vertical == levels. Otherwise it's not proved.

    TODO: In a future version we will have n_vertical independent of levels

    order of the columns in M := macroel_vertices:
        (M[:,0], M[:,1], M[:,2]) == a triangle
        (M[:,3], M[:,4], M[:,5]) == the other triangle
    
    edges perpendicular to the triangles:
        singular edge (if it is singular) == [M[0],M[3]] 
        the others:             [M[1],M[4]], [M[2],M[5]]

    OBS: to sum faster, in the array of points, the levels
    greater than zero end up filled as a rectangle of points. 
    Be careful not to use that coordinates. 

    At the end we translate level 0 to the levels above.
    """
    n_vertical = levels 
    points = np.zeros((n_vertical+1,levels+1,3,levels+1))

    for y in range(levels+1):
        for z in range(levels+1-y):
            lambda_1, lambda_2 = lambda1 (y,z,0,levels,mu), lambda2 (y,z,0,levels,mu)
            temp = lambda_1*(macroel_vertices[:,1] - macroel_vertices[:,0]) + lambda_2*(macroel_vertices[:,2] - macroel_vertices[:,0])
            points[0,y,:,y+z] = temp + macroel_vertices[:,0]
    for x in range(1,n_vertical+1):
        points[x,:,:,:] = points[0,:,:,:] + (float(x)/n_vertical)*(macroel_vertices[:,3] - macroel_vertices[:,0]).reshape((3,1))
    
    ## CONTINUE HERE: 
    ## *use:    vstack((first_level[np.triu_indices(levels+1)], scnd_level[np.triu_indices(levels+1)], ...))
    ## to make the (local_nvert x 3) shaped output
    ## return points

def macroel_more_flexible_tending_more_generality ():
    ## TODO
    pass