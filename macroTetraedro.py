from numpy import *
from mpl_toolkits.mplot3d import Axes3D

import matplotlib.pyplot as plt

color0 	= '#A2B5CD' # lightsteelblue 3
red 	= '#ff0000' # red
green 	= '#006400' # dark green
black  	= '#000000' # black
orange 	= '#EE9A00' # orange 2
blue 	= '#4169E1' # royalblue
red2 	= '#B0171F' # indian red
gold	= '#FFD700'

origin 	= array([0,0,1])

""" Idea: poner funciones f(x,y) que dibujen planos """

def coeficientekj (k,j):
	c3 = ((j+k)/3.)**(.5)*k/3.
	c2 = ((j+k)/3.)**(.5)*j/3.
	c1 = 1 - c2 - c3
	return [c1, c2, c3]

def planoSup (x,y):
	""" coordenada z tq (xyz) est'e en el plano
	antiguo planoSup: return 2 - x
	"""
	X0 		= array([0,0,2])
	X1 		= array([1,0,1])
	X2 		= array([sqrt(2)/2,sqrt(2)/2,1])
	normal 	= cross(X1-X0, X2-X0)
	d 		= inner(normal, X0)
	return (d - x*normal[0] - y*normal[1])/normal[2]

def planoPi0 (x,y):
	""" coordenada z tq (xyz) est'e en el plano
	antiguo Pi0: return 1 - y*(1 - sqrt(2)) """
	return 1

N 			= 8
pisos 		= 4
mallaZ 		= linspace(-1,1,pisos + 1)
complexNum	= exp(1j*array(range(N))*2*pi/N)
reales 		= real(complexNum)
imaginarias = imag(complexNum)

fig 		= plt.figure()
ax 			= fig.add_subplot(1,1,1, projection='3d')
#ax1 		= fig.add_subplot(2,2,2, projection='3d')
#ax2 		= fig.add_subplot(2,2,3, projection='3d')
#ax3 		= fig.add_subplot(2,2,4, projection='3d')

#plt.title("Mesh with singularities")

del(complexNum)

# 1 tri'angulo macro
t = 0

v1 = [0,0]
v2 = [reales[t],imaginarias[t]]
v3 = [reales[t + 1],imaginarias[t + 1]]
	
# primera franja de tri'angulos
listaPares1 = [(0,2),(1,2),(1,1),(2,1),(2,0)]
triangs1 = [[],[]]
for a, b in listaPares1:
	c = coeficientekj(a,b)
	triangs1[0].append(c[0]*v1[0] + c[1]*v2[0] + c[2]*v3[0])
	triangs1[1].append(c[0]*v1[1] + c[1]*v2[1] + c[2]*v3[1])

# segunda franja de tri'angulos
listaPares2 = [(0,1),(1,1),(1,0)]
triangs2 = [[],[]]
for a, b in listaPares2:
	c = coeficientekj(a,b)
	triangs2[0].append(c[0]*v1[0] + c[1]*v2[0] + c[2]*v3[0])
	triangs2[1].append(c[0]*v1[1] + c[1]*v2[1] + c[2]*v3[1])

if t == 0:

	p10 = [triangs1[0][0], triangs1[1][0]]   # (x2, y2)
	p11 = [triangs1[0][1], triangs1[1][1]]   
	p12 = [triangs1[0][2], triangs1[1][2]]   
	p13 = [triangs1[0][3], triangs1[1][3]]   
	p14 = [triangs1[0][4], triangs1[1][4]]   

	q20 = [triangs2[0][0], triangs2[1][0]]
	q21 = [triangs2[0][1], triangs2[1][1]]
	q22 = [triangs2[0][2], triangs2[1][2]]

	ax.plot([p10[1],triangs1[0][4]],
			[triangs1[1][0],triangs1[1][4]],
			[planoSup(p10[1],triangs1[1][0]),
				planoSup(triangs1[0][4],triangs1[1][4])],
			color = color0)
	ax.plot([q20[0], triangs2[0][2]],
			[q20[1], q22[1]],
			[planoSup(q20[0],q20[1]),
				planoSup(triangs2[0][2],q22[1])],
			color = color0)
	# verticales (z=1) ---> (z = 2-x)
	# INTERCALAR EL PLANO Pi4
	ax.plot([q20[0], q20[0]],
			[q20[1], q20[1]],
			[1, planoSup(q20[0],q20[1])],
			color = color0)
	# subdivido primera parte en 2 tetra grandes
	ax.plot([reales[0],reales[1],0,reales[0]],
			[imaginarias[0],imaginarias[1],0,imaginarias[0]],
			[1,planoSup(reales[1],imaginarias[1]),1,1],
			color = color0)
	# plano: pi2  (sqrt(2)/2 - 1)*Y + (sqrt(2)/2)*Z = 1
	# TODO poner este plano en funci'on de tres puntos, etc.
	ax.plot([q20[0],triangs2[0][2]],
			[q20[1],q22[1]],
			[planoPi0(q20[0],q20[1]),
				planoPi0(triangs2[0][2],q22[1])],
			color = color0)
	ax.plot([p10[1],triangs1[0][4]],
			[triangs1[1][0],triangs1[1][4]],
			[planoPi0(p10[1],triangs1[1][0]),
				planoPi0(triangs1[0][4],triangs1[1][4])],
			color = color0)
	
	"""esto va como un plano Pi2
	def planoPi2(x,y):
		return z
	"""
	ax.plot([q20[0],origin[1],q22[0]],
			[q20[1],q20[1],q22[1]],
			[planoSup(q20[0],q20[1]),planoSup(q20[0],q20[1]),planoSup(q22[0],q22[1])],
			color = black)

	xCoords = [p10[0],origin[1],p14[0]]
	yCoords = [p10[1],p10[1],p14[1]]
	zCoords = [planoSup(p10[0],p10[1]),planoSup(p10[0],p10[1]),planoSup(p14[0],p14[1])]
	ax.plot(xCoords,
			yCoords,
			zCoords,
			color = color0)

	def planoPi1(x,y):
		""" coordenada z tq (xyz) est'e en el plano """
		P0 		= array([xCoords[1],yCoords[1],zCoords[1]])
		P1 		= array([xCoords[0],yCoords[0],zCoords[0]])
		P2 		= array([xCoords[2],yCoords[2],zCoords[2]])
		normal 	= cross(P1-P0, P2-P0)
		d 		= inner(normal, P0)
		return (d - x*normal[0] - y*normal[1])/normal[2]
	
	ax.plot([q22[0],q22[0]],
			[q22[1],q22[1]],
			[planoPi0(q22[0],q22[1]), planoPi1(q22[0],q22[1])],
			color = color0)

	# verticales
	ax.plot([p10[0],p10[0]],[p10[1],p10[1]],[planoPi0(p10[0],p10[1]),planoSup(p10[0],p10[1])],color = red)
	ax.plot([p12[0],p12[0]],
			[p12[1],p12[1]],
			[planoPi0(p12[0],p12[1]),planoSup(p12[0],p12[1])],
			color = red)
	ax.plot([p14[0],p14[0]],[p14[1],p14[1]],[planoPi0(p14[0],p14[1]),planoSup(p14[0],p14[1])],color = red)

	# ------------ pyramid !
	ax.plot([q20[0],p12[0],q22[0],q20[0]],
			[q20[1],p12[1],q22[1],q20[1]],
			[planoPi1(q20[0],q20[1]),planoSup(p12[0],p12[1]),planoPi1(q22[0],q22[1]),planoPi1(q20[0],q20[1])],
			color = red)

	ax.plot([q20[0],p12[0],q22[0],q20[0]],
			[q20[1],p12[1],q22[1],q20[1]],
			[planoSup(q20[0],q20[1]),planoSup(p12[0],p12[1]),planoSup(q22[0],q22[1]),planoSup(q20[0],q20[1])],
			color = red)

	ax.plot([q20[0],q20[0]],[q20[1],q20[1]],[planoSup(q20[0],q20[1]),planoPi1(q20[0],q20[1])],color=red)
	ax.plot([q22[0],q22[0]],[q22[1],q22[1]],[planoSup(q22[0],q22[1]),planoPi1(q22[0],q22[1])],color=red)
	# ------------
	
	# "V" y "W" en el plano Pi2
	ax.plot([q22[0],p12[0],q20[0]],
			[q22[1],p12[1],q20[1]],
			[planoPi0(q22[0],q22[1]),planoPi0(p12[0],p12[1]),1],
			color = color0)

	# ------------ pyramid !
	ax.plot([p10[0],p11[0],p12[0],p13[0],p14[0]],
			[p10[1],p11[1],p12[1],p13[1],p14[1]],
			[planoPi0(p10[0],p10[1]),
				planoPi0(p11[0],p11[1]),
				planoPi0(p12[0],p12[1]),
				planoPi0(p13[0],p13[1]),
				planoPi0(p14[0],p14[1])],
			color = red)

	ax.plot([p14[0],p13[0],p12[0],p11[0],p10[0]],
			[p14[1],p13[1],p12[1],p11[1],p10[1]],
			[planoSup(p14[0],p14[1]),
				planoSup(p13[0],p13[1]),
				planoSup(p12[0],p12[1]),
				planoSup(p11[0],p11[1]),
				planoSup(p10[0],p10[1])],
			color = red)

	ax.plot([p10[0], p12[0], p14[0]],
			[p10[1], p12[1], p14[1]],
			[planoPi0(p10[0],p10[1]),
				planoPi0(p12[0],p12[1]),
				planoPi0(p14[0], p14[1])],
			color = red)

	ax.plot([p10[0], p12[0], p14[0]],
			[p10[1], p12[1], p14[1]],
			[planoSup(p10[0],p10[1]),
				planoSup(p12[0],p12[1]),
				planoSup(p14[0], p14[1])],
			color = red)

 # s'i
ax.plot([1,0,v3[0]],[0,0,v3[1]],[1,2,planoPi0(v3[0],v3[1])],color = color0)

ax.plot([0,0],
		[0,0],
		[mallaZ[pisos],2],
		color = orange)

plt.show()
