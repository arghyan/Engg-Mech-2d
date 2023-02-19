from sympy.physics.mechanics import ReferenceFrame
from kinmod4 import *
from sympy import symbols
'''
Step 1: Select Reference Frames
'''
N=ReferenceFrame('N')
'''
Step 2: Use symbolic variables
'''
Ax,Ay=symbols('Ax,Ay')
By=symbols('By')
'''
Step 3: Define points where forces are applied
'''
A=ourpoint('A',N,[0,0])
B=ourpoint('B',N,[10,0])
A.setfor_xy(N,[Ax,Ay])
A.setmom_direct(N,-80)
B.setfor_xy(N,[0,By])
G=ourpoint('G',N,[5,0])
G.setfor_xy(N,[0,-15])
'''
Step 4: Define bodies
'''
Bd=ourbody('Bd',[A,G,B],A)
Bd.input_dist_load([G,B,-5,-5])
'''
Step 5: Generate equation for each body
'''
Bd.generate_eqn(N)
'''
Step 6: Solution
'''
Bd.solve_eqn([Ax,Ay,By])
'''
Step 7: SF BM diagram
'''
Bd.draw_sfbm(N)