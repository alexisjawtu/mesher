import numpy as np
mu = .65
hybrid_color = "green"
isotropic_color = "red"
prism_h   = np.array([0,0,4])

#horiz1  = 2    

x_max = 3
x_int_max = 2
x_min   = -1
x_int_min   = 0
y_max   = 1
y_min   = -3
y_int_min   = -2
y_int_max   = 0
z_max   = 0

#A0 = np.zeros(3)
#Q0 = np.array([x_int_min,y_int_max,-1])
#Q1 = np.array([0,1,-1])
#Q2 = np.array([-1,0,-1])
#R0 = Q1 - Q0 + Q2
#P1_hybrid_4 = Q2+Q1-2*Q0+A0

macro_elements = { 
0 : { 0 : np.array([[x_int_min,y_int_max,-1],[-1,0,-1],[0,1,-1],[0,0,0]]), 1 : mu, 2 : hybrid_color,    3 : 0},
1 : { 0 : np.array([[-1,0,0],[-1,1,0],[-1,0,-1],[0,0,0]]),                 1 : mu, 2 : hybrid_color,    3 : 0}, 
2 : { 0 : np.array([[0,1,0],[0,1,-1],[-1,1,0],[0,0,0]]),                   1 : mu, 2 : hybrid_color,    3 : 0},
3 : { 0 : np.array([[-1,1,-1],[-1,0,-1],[-1,1,0],[0,1,-1]]),               1 : 1,  2 : hybrid_color,    3 : 0},
4 : { 0 : np.array([[0,0,0],[-1,0,-1],[-1,1,0],[0,1,-1]]),                 1 : mu, 2 : isotropic_color, 3 : 1},

5 : { 0 : np.array([-1,1,1])*np.array([[x_int_min,y_int_max,-1],[-1,0,-1],[0,1,-1],[0,0,0]])  + np.array([2,0,0]), 1 : mu,2 : hybrid_color,    3 : 0},
6 : { 0 : np.array([-1,1,1])*np.array([[-1,0,0],[-1,1,0],[-1,0,-1],[0,0,0]])                  + np.array([2,0,0]), 1 : mu,2 : hybrid_color,    3 : 0},
7 : { 0 : np.array([-1,1,1])*np.array([[0,1,0],[0,1,-1],[-1,1,0],[0,0,0]])                    + np.array([2,0,0]), 1 : mu,2 : hybrid_color,    3 : 0},
8 : { 0 : np.array([-1,1,1])*np.array([[-1,1,-1],[-1,0,-1],[-1,1,0],[0,1,-1]])                + np.array([2,0,0]), 1 : 1 ,2 : hybrid_color,    3 : 0},
9 : { 0 : np.array([-1,1,1])*np.array([[0,0,0],[-1,0,-1],[-1,1,0],[0,1,-1]])                  + np.array([2,0,0]), 1 : mu,2 : isotropic_color, 3 : 1},

10 : {0: np.array([[0,-2,-1],[-1,-2,-1],[0,-3,-1],[0,-2,z_max]])        ,  1 : mu,2 : hybrid_color,    3 : 0},
11 : {0: np.array([[0,-3,0],[0,-3,-1],[-1,-3,0],[0,-2,z_max]])          ,  1 : mu,2 : hybrid_color,    3 : 0},
12 : {0: np.array([[-1,-2,z_max],[-1,-3,z_max],[-1,-2,-1],[0,-2,z_max]]),  1 : mu,2 : hybrid_color,    3 : 0},
13 : {0: np.array([[-1,-3,-1],[-1,-3,z_max],[0,-3,-1],[-1,-2,-1]])      ,  1 : 1 ,2 : hybrid_color,    3 : 0},
14 : {0: np.array([[ 0, -2,  0],[-1, -2, -1],[ 0, -3, -1],[-1, -3,  0]]),  1 : mu,2 : isotropic_color, 3 : 1},

15 : {0 : np.array([[x_int_min,-1,z_max],[x_int_min,-1,-1],[x_min,-1,z_max],[x_int_min, y_int_max,z_max]])   ,  1 : mu,2 : hybrid_color,    3 : 0},
16 : {0 : np.array([[x_int_min,y_int_max,-1],[x_min,y_int_max,-1],[x_int_min,-1,-1],[x_int_min,y_int_max,0]]),  1 : mu,2 : hybrid_color,    3 : 0},
17 : {0 : np.array([[-1,0,0],[-1,-1,0],[-1,0,-1],[0,0,0]])                                                   ,  1 : mu,2 : hybrid_color,    3 : 0},
18 : {0 : np.array([[-1,-1,-1],[-1,-1,0],[0,-1,-1],[-1,0,-1]])                                               ,  1 : 1 ,2 : hybrid_color,    3 : 0},
19 : {0 : np.array([[x_int_min,y_int_max,z_max],[-1,-1,z_max],[x_min,y_int_max,-1],[x_int_min,-1,-1]])       ,  1 : mu,2 : isotropic_color, 3 : 1},

20 : {0 : np.array([[x_int_min,y_int_min,-1],[x_int_min,-1,-1],[x_min,y_int_min,-1],[x_int_min,y_int_min,0]]),    1 : mu,2 : hybrid_color,    3 : 0},
21 : {0 : np.array([[x_int_min, -1, z_max],[x_min, -1, z_max],[x_int_min, -1, -1],[x_int_min, y_int_min, z_max]]),1 : mu,2 : hybrid_color,    3 : 0},
22 : {0 : np.array([[-1,-2,0],[-1,-2,-1],[-1,-1,0],[0,-2,0]]) ,                                                   1 : mu,2 : hybrid_color,    3 : 0},
23 : {0 : np.array([[-1,-1,-1],[0,-1,-1],[-1,-1,0],[-1,-2,-1]]),                                                  1 : 1 ,2 : hybrid_color,    3 : 0},
24 : {0 : np.array([[x_int_min,y_int_min,z_max],[-1,-1,z_max],[x_int_min,-1,-1],[x_min,y_int_min,-1]]),           1 : mu,2 : isotropic_color, 3 : 1},

25 : {0 : np.array([[ 2,-2,-1],[ 3,-2,-1],[ 2,-3,-1],[ 2,-2, 0]]),1 : mu,2 : hybrid_color,    3 : 0},
26 : {0 : np.array([[ 2,-3, 0],[ 2,-3,-1],[ 3,-3, 0],[ 2,-2, 0]]),1 : mu,2 : hybrid_color,    3 : 0},
27 : {0 : np.array([[ 3,-2, 0],[ 3,-3, 0],[ 3,-2,-1],[ 2,-2, 0]]),1 : mu,2 : hybrid_color,    3 : 0},
28 : {0 : np.array([[ 3,-3,-1],[ 3,-3, 0],[ 2,-3,-1],[ 3,-2,-1]]),1 : 1 ,2 : hybrid_color,    3 : 0},
29 : {0 : np.array([[ 2,-2, 0],[ 3,-2,-1],[ 2,-3,-1],[ 3,-3, 0]]),1 : mu,2 : isotropic_color, 3 : 1},


30 : {0: np.array([[2,-1, 0],[2,-1,-1],[3,-1, 0],[2,0,0]]),1 : mu,2 : hybrid_color,    3 : 0},
31 : {0: np.array([[2,0,-1],[3,0,-1],[2,-1,-1],[2,0,0]]),1 : mu,2 : hybrid_color,    3 : 0},
32: {0: np.array([[3,0,0],[3,-1,0],[3,0,-1],[2,0,0]]),1 : mu,2 : hybrid_color,    3 : 0},
33 : {0: np.array([[3,-1,-1],[3,-1,0],[2,-1,-1],[3,0,-1]]),1 : 1 ,2 : hybrid_color,    3 : 0},
34 : {0: np.array([[2,0,0],[3,-1,0],[3,0,-1],[2,-1,-1]]),1 : mu,2 : isotropic_color, 3 : 1},


35 : {0 : np.array([[2,-2,-1],[2,-1,-1],[3,-2,-1],[2,-2, 0]]),1 : mu,2 : hybrid_color,    3 : 0},
36 : {0 : np.array([[2,-1, 0],[3,-1, 0],[2,-1,-1],[2,-2, 0]]),1 : mu,2 : hybrid_color,    3 : 0},
37 : {0 : np.array([[3,-2, 0],[3,-2,-1],[3,-1, 0],[2,-2, 0]]),1 : mu,2 : hybrid_color,    3 : 0},
38 : {0 : np.array([[3,-1,-1],[2,-1,-1],[3,-1, 0],[3,-2,-1]]),1 : 1 ,2 : hybrid_color,    3 : 0},
39 : {0 : np.array([[2,-2, 0],[3,-1, 0],[2,-1,-1],[3,-2,-1]]),1 : mu,2 : isotropic_color, 3 : 1},

}
