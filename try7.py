from sympy.physics.mechanics import ReferenceFrame
from kinmod4 import *
from sympy import symbols
'''
Step 1: Select Reference Frames
'''
N=ReferenceFrame('N')
N1=N.orientnew('N1','Axis',[pi/6,N.z])
'''
Step 2: Use symbolic variables
'''
Cx,Cy,Dx,Dy,Ax,Ay,Bx,By,M=symbols('Cx,Cy,Dx,Dy,Ax,Ay,Bx,By,M')
'''
Step 3: Define points where forces are applied
'''
C=ourpoint('C',N,[0,0])
A1=ourpoint('A1',N1,[1.5,0])
A2=ourpoint('A2',N1,[1.5,0])
D=ourpoint('D',N,[0,1.8])
B1=ourpoint('B1',N,[0,1.8],N1,[1.5,0])
B2=ourpoint('B2',N,[0,1.8],N1,[1.5,0])
G=ourpoint('G',N1,[1.5,0],N,[0,1.2])
C.setfor_xy(N1,[Cx,Cy])
C.setmom_direct(N, M)
D.setfor_xy(N1,[Dx,0])
A1.setfor_xy(N1,[Ax,Ay])
A2.setfor_xy(N1,[-Ax,-Ay])
B1.setfor_xy(N1,[Bx,0])
B2.setfor_xy(N1,[-Bx,0])
G.setfor_xy(N,[0,-150*9.8])
'''
Step 4: Define bodies
'''
BodyCA=ourbody('CA',[C,A1],C)
BodyDB=ourbody('DB',[D,B1],D)
BodyAB=ourbody('AB',[A2,G,B2],G)
m=150
Ig=0
BodyAB.set_for_dynamics(G,m,G,Ig,N1,[-13.455,13.725],0)
'''
Step 5: Generate equation for each body
'''
BodyCA.generate_eqn(N1)
BodyDB.generate_eqn(N1)
BodyAB.generate_eqn(N1)
'''
Step 6: Combine constituent bodies in a group
'''
group1 = group_of_bodies('group1',[BodyCA,BodyDB,BodyAB])
'''
Step 7: Solution for the group
'''
group1.gr_solve([Cx,Cy,Dx,Ax,Ay,Bx,M])