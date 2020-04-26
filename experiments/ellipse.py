import numpy as np
import matplotlib.pyplot as plt 
from mpl_toolkits.mplot3d import Axes3D

# Equations of the planes, according to quadrant.
plane_equation = { 0: lambda p: 2*(p[0,:]+p[1,:])/3+p[2,:],
                   1: lambda p: 2*(-p[0,:]+p[1,:])/3+p[2,:],
                   2: lambda p: -2*(p[0,:]+p[1,:])/3+p[2,:],
                   3: lambda p: 2*(p[0,:]-p[1,:])/3+p[2,:] }
# Equations of the lines, after projecting. 
proj_equation = { 0: lambda q: q[0,:]+q[1,:],
                  1: lambda q: -q[0,:]+q[1,:],
                  2: lambda q: -(q[0,:]+q[1,:]),
                  3: lambda q: q[0,:]-q[1,:] }

def classify(p,dec=2):
    # classify(p) classifies every point in the grid. 
    # output: classification, dictionary with three vectors: 
    #     pole = -1,0,1 if the node is in the south pole, cylinder or north pole
    #     quadrant = 0,1,2,3 if the node is in the first to forth quadrat (seen from the north)
    #     b = value of the right hand side of the equation of the plane where the node belongs. 
    corner = np.zeros(p.shape[1])
    polo = np.zeros(p.shape[1])
    corner[np.where((p[0,:]<=0)*(p[1,:]>0))] = 1
    corner[np.where((p[0,:]<0)*(p[1,:]<=0))] = 2
    corner[np.where((p[0,:]>=0)*(p[1,:]<0))] = 3
    polo[np.where(p[2,:]>=1)] = 1
    polo[np.where(p[2,:]<=-1)] = -1
    b = np.zeros(p.shape[1])
    for i in range(4):
        # North pole:
        tramo = (corner==i)*(polo==1)
        b[tramo] = plane_equation[i](p[:,tramo])
        # South pole: equation from the oposite quadrant.
        tramo = (corner==i)*(polo==-1)
        b[tramo] = plane_equation[np.mod(i-2,4)](p[:,tramo])
    b = np.round(b,decimals = dec)
    classification = { 'pole' : polo,
                       'quadrant' : corner,
                       'b' : b }
    return classification

def project_and_scale(p):
    N = p.shape[1]
    classification = classify(p)
    b = classification['b']
    b_set = np.unique(b)
    line = np.zeros(N)
    polo = classification['pole']
    for i in range(4):
        # Equation of the line
        tramo = (classification['quadrant']==i)
        line[tramo] = proj_equation[i](p[:2,tramo])
    coeff = np.zeros(N)
    height = np.zeros(N)
    for b_elem in b_set:
        # coeff: scale coefficient, it tracks the projected triangle of the plane. 
        coeff[b==b_elem] = np.max(np.abs(line[b==b_elem]))
        height[b==b_elem] = polo[b==b_elem]*np.max(np.abs(p[2,b==b_elem]))-polo[b==b_elem]
    coeff[coeff==0]=1 # to avoid divide by 0. 
    norm = np.linalg.norm(p[:2,:],axis=0)
    norm[norm==0]=1
    # project (delete z coordinate), divide by norm (circle of radious 1) and multiple by 
    #   (line/coeff) which push the points to arcs of circles proportional to line and scaled
    #   by coeff
    t = (line/coeff)*p[:2,:]/norm
    return t,coeff,height,classification

def stereo_projection(t):
    return np.array([2*t[0,:]/(1+t[0,:]**2+t[1,:]**2),\
                     2*t[1,:]/(1+t[0,:]**2+t[1,:]**2),\
                     (1-t[0,:]**2-t[1,:]**2)/(1+t[0,:]**2+t[1,:]**2)])

def scale_ellipse(p,coeff,height,classification):
    polo = classification['pole']
    # coeff push the x and y coordinates to the ellipsoid. height do the same with z.  
    p[0,:] = coeff*p[0,:]
    p[1,:] = coeff*p[1,:]
    p[2,:] = height*p[2,:]
    # move to the singular vertex (1 or -1)
    correction = np.vstack((np.zeros((2,p.shape[1])),polo))
    return p+correction



#v =np.array(p[2,:]>=1)
#p = p[:,v]
#fig = plt.figure()
#ax  = fig.add_subplot(1,1,1,projection='3d')
#ax.scatter(p[0,:],p[1,:],p[2,:],'b')

#q,coeff,height,classification = project_and_scale(p)
#q = stereo_projection(q)
#q = scale_ellipse(q,coeff,height,classification)
#ax.scatter(q[0,:],q[1,:],q[2,:],'r',marker='v')
#plt.show()



