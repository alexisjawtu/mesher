""" is this good as input
0
x y z
x y z
x y z
x y z
mu
'color'? 
"""

mu = .65
hybrid_color = "green"
isotropic_color = "red"
prism_h   = np.array([0,0,4])
horiz1  = 2    
x_max = 3
x_int_max = 2
x_min   = -1
x_int_min   = 0
y_max   = 1
y_min   = -3
y_int_min   = -2
y_int_max   = 0
z_max   = 0
A0 = np.zeros(3)
Q0 = np.array([x_int_min,y_int_max,-1])
Q1 = np.array([0,1,-1])
Q2 = np.array([-1,0,-1])
R0 = Q1 - Q0 + Q2
P1_hybrid_4 = Q2+Q1-2*Q0+A0
D = { 
0 : { 0 : np.array([[x_int_min,y_int_max,-1],[-1,0,-1],[0,1,-1],[0,0,0]]), 1 : mu, 2 : hybrid_color,    3 : 0},
1 : { 0 : np.array([[-1,0,0],[-1,1,0],[-1,0,-1],[0,0,0]]),                 1 : mu, 2 : hybrid_color,    3 : 0}, 
2 : { 0 : np.array([[0,1,0],[0,1,-1],[-1,1,0],[0,0,0]]),                   1 : mu, 2 : hybrid_color,    3 : 0},
3 : { 0 : np.array([[-1,1,-1],[-1,0,-1],[-1,1,0],[0,1,-1]]),               1 : 1,  2 : hybrid_color,    3 : 0},
4 : { 0 : np.array([[0,0,0],[-1,0,-1],[-1,1,0],[0,1,-1]]),                 1 : mu, 2 : isotropic_color, 3 : 1}

5 : { 0 : np.array([-1,1,1])*np.array([[x_int_min,y_int_max,-1],[-1,0,-1],[0,1,-1],[0,0,0]])  + np.array([2,0,0]), 1 : mu,2 : hybrid_color,    3 : 0},
6 : { 0 : np.array([-1,1,1])*np.array([[-1,0,0],[-1,1,0],[-1,0,-1],[0,0,0]])                  + np.array([2,0,0]), 1 : mu,2 : hybrid_color,    3 : 0},
7 : { 0 : np.array([-1,1,1])*np.array([[0,1,0],[0,1,-1],[-1,1,0],[0,0,0]])                    + np.array([2,0,0]), 1 : mu,2 : hybrid_color,    3 : 0},
8 : { 0 : np.array([-1,1,1])*np.array([[-1,1,-1],[-1,0,-1],[-1,1,0],[0,1,-1]])                + np.array([2,0,0]), 1 : 1 ,2 : hybrid_color,    3 : 0},
9 : { 0 : np.array([-1,1,1])*np.array([[0,0,0],[-1,0,-1],[-1,1,0],[0,1,-1]])                  + np.array([2,0,0]), 1 : mu,2 : isotropic_color, 3 : 1}

10 : {0: np.array([[0,-2,-1],[-1,-2,-1],[0,-3,-1],[0,-2,z_max]])        ,  1 : mu,2 : hybrid_color,    3 : 0},
11 : {0: np.array([[0,-3,0],[0,-3,-1],[-1,-3,0],[0,-2,z_max]])          ,  1 : mu,2 : hybrid_color,    3 : 0},
12 : {0: np.array([[-1,-2,z_max],[-1,-3,z_max],[-1,-2,-1],[0,-2,z_max]]),  1 : mu,2 : hybrid_color,    3 : 0},
13 : {0: np.array([[-1,-3,-1],[-1,-3,z_max],[0,-3,-1],[-1,-2,-1]])      ,  1 : 1 ,2 : hybrid_color,    3 : 0},
14 : {0: np.array([[ 0, -2,  0],[-1, -2, -1],[ 0, -3, -1],[-1, -3,  0]]),  1 : mu,2 : isotropic_color, 3 : 1}

}
"""

vertices_tetra_1_reflected   = np.array([-1,1,1])*vertices_tetra_1   + np.array([2,0,0])
vertices_hybrid_11_reflected = np.array([-1,1,1])*vertices_hybrid_11 + np.array([2,0,0])
vertices_hybrid_12_reflected = np.array([-1,1,1])*vertices_hybrid_12 + np.array([2,0,0])
vertices_hybrid_13_reflected = np.array([-1,1,1])*vertices_hybrid_13 + np.array([2,0,0])
vertices_hybrid_14_reflected = np.array([-1,1,1])*vertices_hybrid_14 + np.array([2,0,0])


vertices_hybrid_21   = np.array([[x_int_min,-1,z_max],[x_int_min,-1,-1],[x_min,-1,z_max],[x_int_min, y_int_max,z_max]])
vertices_hybrid_22   = np.array([[x_int_min,y_int_max,-1],[x_min,y_int_max,-1],[x_int_min,-1,-1],[x_int_min,y_int_max,0]])
vertices_hybrid_23   = np.array([[-1,0,0],[-1,-1,0],[-1,0,-1],[0,0,0]])
vertices_hybrid_24   = np.array([[-1,-1,-1],[-1,-1,0],[0,-1,-1],[-1,0,-1]])
vertices_tetra_2     = np.array([[x_int_min,y_int_max,z_max],[-1,-1,z_max],[x_min,y_int_max,-1],[x_int_min,-1,-1]])


vertices_tetra_2_reflected   = np.array([-1,1,1])*vertices_tetra_2+np.array([2,0,0])
vertices_hybrid_21_reflected   = np.array([-1,1,1])*vertices_hybrid_21+np.array([2,0,0])
vertices_hybrid_22_reflected   = np.array([-1,1,1])*vertices_hybrid_22+np.array([2,0,0])
vertices_hybrid_23_reflected   = np.array([-1,1,1])*vertices_hybrid_23+np.array([2,0,0])
vertices_hybrid_24_reflected   = np.array([-1,1,1])*vertices_hybrid_24+np.array([2,0,0])

vertices_tetra_3    = np.array([[x_int_min,y_int_min,z_max],[-1,-1,z_max],[x_int_min,-1,-1],[x_min,y_int_min,-1]])
vertices_hybrid_31   = np.array([[x_int_min,y_int_min,-1],[x_int_min,-1,-1],[x_min,y_int_min,-1],[x_int_min,y_int_min,0]])
vertices_hybrid_32   = np.array([[x_int_min, -1, z_max],[x_min, -1, z_max],[x_int_min, -1, -1],[x_int_min, y_int_min, z_max]])
vertices_hybrid_33   = np.array([[-1,-2,0],[-1,-2,-1],[-1,-1,0],[0,-2,0]]) 
vertices_hybrid_34  = np.array([[-1,-1,-1],[0,-1,-1],[-1,-1,0],[-1,-2,-1]])


vertices_tetra_3_reflected = np.array([-1,1,1])*vertices_tetra_3+np.array([2,0,0])
vertices_hybrid_31_reflected = np.array([-1,1,1])*vertices_hybrid_31+np.array([2,0,0])
vertices_hybrid_32_reflected = np.array([-1,1,1])*vertices_hybrid_32+np.array([2,0,0])
vertices_hybrid_33_reflected = np.array([-1,1,1])*vertices_hybrid_33+np.array([2,0,0])
vertices_hybrid_34_reflected = np.array([-1,1,1])*vertices_hybrid_34+np.array([2,0,0])

vertices_prisms = []
points_prisms   = np.array([Q0,Q1,Q2])
vertices_prisms = vertices_prisms + [np.concatenate((points_prisms,points_prisms - prism_h))]
points_prisms   = np.array([R0,Q1,Q2])
vertices_prisms = vertices_prisms + [np.concatenate((points_prisms,points_prisms - prism_h))]
vertices_prisms      = vertices_prisms + [np.array([[-1,-3,-1],[0,-3,-1],[-1,-2,-1],[-1,-3,-1]-prism_h,[0,-3,-1]-prism_h,[-1,-2,-1]-prism_h])]
vertices_prisms      = vertices_prisms + [np.array([[0,-2,-1],[-1,-2,-1],[0,-3,-1],[0,-2,-1]-prism_h,[-1,-2,-1]-prism_h,[0,-3,-1]-prism_h])]
vertices_prisms      = vertices_prisms + [np.array([[0,-2,-1],[0,-1,-1], [-1,-2,-1], [0,-2,-1]-prism_h,[0,-1,-1]-prism_h,[0,-2,-1]-prism_h])]
vertices_prisms      = vertices_prisms + [np.array([[0,-2,-1],[0,-3,-1],[1,-2,-1],[0,-2,-1]-prism_h,[0,-3,-1]-prism_h,[1,-2,-1]-prism_h])]

vertices_prisms      = vertices_prisms + [np.array([[-1,-1,-1],[-1,-2,-1],[0,-1,-1],[-1,-1,-1]-prism_h,[-1,-2,-1]-prism_h,[0,-1,-1]-prism_h])]
vertices_prisms      = vertices_prisms + [np.array([[1,y_min,-1],[1,y_int_min,-1],[0,y_min,-1],
                                         [1,y_min,-1]-prism_h,[1,y_int_min,-1]-prism_h,[0,y_min,-1]-prism_h])]
vertices_prisms      = vertices_prisms +  [np.array([[x_int_min,y_int_max,-1],[x_min,y_int_max,-1],[x_int_min,-1,-1],
                                         [x_int_min,y_int_max,-1]-prism_h,[x_min,y_int_max,-1]-prism_h,[x_int_min,-1,-1]-prism_h])]
vertices_prisms      = vertices_prisms + [np.array([[x_min,-1,-1],[x_int_min,-1,-1],[x_min,y_int_max,-1],
                                        [x_min,-1,-1]-prism_h,[x_int_min,-1,-1]-prism_h,[x_min,y_int_max,-1]-prism_h])]
vertices_prisms      = vertices_prisms + [np.array([[1,y_min,-1],[2,y_min,-1],[1,y_int_min,-1],
                                            [1,y_min,-1]-prism_h,[2,y_min,-1]-prism_h,[1,y_int_min,-1]-prism_h])]
vertices_prisms      = vertices_prisms + [np.array([[2,y_int_min,-1],[1,y_int_min,-1],[2,y_min,-1],
                                            [2,y_int_min,-1]-prism_h,[1,y_int_min,-1]-prism_h,[2,y_min,-1]-prism_h])]
vertices_prisms      = vertices_prisms + [np.array([[2,y_int_min,-1],[2,y_min,-1],[x_max,y_int_min,-1],
                                         [2,y_int_min,-1]-prism_h,[2,y_min,-1]-prism_h,[x_max,y_int_min,-1]-prism_h])]
vertices_prisms      = vertices_prisms + [np.array([[x_max,y_min,-1],[x_max,y_int_min,-1],[x_int_max,y_min,-1],
                                         [x_max,y_min,-1]-prism_h,[x_max,y_int_min,-1]-prism_h,[x_int_max,y_min,-1]-prism_h])]
vertices_prisms      = vertices_prisms + [np.array([[2,y_int_min,-1],[x_max,y_int_min,-1],[x_int_max,-1,-1],
                                         [2,y_int_min,-1]-prism_h,[x_max,y_int_min,-1]-prism_h,[x_int_max,-1,-1]-prism_h])]
vertices_prisms      = vertices_prisms + [np.array([[x_max,-1,-1],[x_int_max,-1,-1],[x_max,y_int_min,-1],
                                         [x_max,-1,-1]-prism_h,[x_int_max,-1,-1]-prism_h,[x_max,y_int_min,-1]-prism_h])]
vertices_prisms      = vertices_prisms + [np.array([[x_max,-1,-1],[x_max,y_int_max,-1],[x_int_max,-1,-1],
                                         [x_max,-1,-1]-prism_h,[x_max,y_int_max,-1]-prism_h,[x_int_max,-1,-1]-prism_h])]

vertices_prisms      = vertices_prisms + [np.array([[x_int_max,y_int_max,-1],[x_int_max,-1,-1],[x_max,y_int_max,-1],
                                         [x_int_max,y_int_max,-1]-prism_h,[x_int_max,-1,-1]-prism_h,[x_max,y_int_max,-1]-prism_h])]

vertices_prisms      = vertices_prisms + [np.array([[x_int_max,y_int_max,-1],           [x_max,y_int_max,-1],[x_int_max,y_max,-1],
                                         [x_int_max,y_int_max,-1]-prism_h,   [x_max,y_int_max,-1]-prism_h,[x_int_max,y_max,-1]-prism_h])]

vertices_prisms      = vertices_prisms + [np.array([[x_max,y_max,-1],[x_int_max,y_max,-1],[x_max,y_int_max,-1],
                                        [x_max,y_max,-1]-prism_h,[x_int_max,y_max,-1]-prism_h,[x_max,y_int_max,-1]-prism_h])]

vertices_prisms      = vertices_prisms + [np.array([[x_int_max,y_int_max,-1],[x_int_max,y_max,-1],[1,y_int_max,-1],
                            [x_int_max,y_int_max,-1]-prism_h,[x_int_max,y_max,-1]-prism_h,[1,y_int_max,-1]-prism_h])]

vertices_prisms      = vertices_prisms + [np.array([[1,y_max,-1],[1,y_int_max,-1],[x_int_max,y_max,-1],
                            [1,y_max,-1]-prism_h,[1,y_int_max,-1]-prism_h,[x_int_max,y_max,-1]-prism_h])]
    
vertices_prisms      = vertices_prisms + [np.array([[1,y_max,-1],[x_int_min,y_max,-1],[1,y_int_max,-1],
                            [1,y_max,-1]-prism_h,[x_int_min,y_max,-1]-prism_h,[1,y_int_max,-1]-prism_h])]

vertices_prisms      = vertices_prisms + [np.array([[x_int_min,y_int_max,-1],[1,y_int_max,-1],[x_int_min,y_max,-1],
                            [x_int_min,y_int_max,-1]-prism_h,[1,y_int_max,-1]-prism_h,[x_int_min,y_max,-1]-prism_h])]
"""