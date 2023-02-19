from sympy.physics.mechanics import ReferenceFrame
from sympy import symbols
from math import asin,cos,sin,pi,radians
from kinmod4 import ourpoint, oureqn

#Step 1
N=ReferenceFrame('N')
N1=N.orientnew('N1','Axis',[pi/4,N.z])
N2=N.orientnew('N2','Axis',[2*pi/3,N.z])

#Step 2
O=ourpoint('O',N,[0,0])
A1=ourpoint('A1',N1,[0.250,0])
A2=ourpoint('A2',N1,[0.250,0])
BA=0.250*sin(radians(45))/sin(radians(60))
OB=0.250*cos(radians(45))+BA*cos(radians(60))
B=ourpoint('B',N,[OB,0])

#Step 3
O.setvel_direct(N,[0,0])
B.setvel_direct(N,[0,0])

#Step 4
omBC,vrel=symbols('omBC,vrel')

#Step 5
A1.calvel_rigid(2*N.z,O)
A2.calvel_sliding(omBC*N.z,B,N2,[vrel,0])

#Step 6
e=oureqn()
e.get_eqnv(A1,A2,N)
e.solvev('omBC','vrel')

#Step 7
O.setaccl_direct(N,[0,0])
B.setaccl_direct(N,[0,0])

#Step 8
alBC,arel=symbols('alBC,arel')

#Step 9
A1.calaccl_rigid(2*N.z,0.0*N.z,O)
A2.calaccl_sliding(omBC*N.z,alBC*N.z,B,N2,[vrel,0],N2,[arel,0])

#Step 10
e.get_eqna(A1,A2,N)
e.subsa(omBC=e.vsol[omBC],vrel=e.vsol[vrel])
e.solvea(alBC,arel)

#Step 11
e.draw(N,'a')