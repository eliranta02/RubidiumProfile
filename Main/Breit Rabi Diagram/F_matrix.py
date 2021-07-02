from __future__ import division
from numpy import *
from J_matrix import getJ
from I_matrix import *



I=3/2
J=1/2
'''
Fx=kron(Jx,eye(2*I+1))+kron(eye(2*J+1),Ix)
Fy=kron(Jy,eye(2*I+1))+kron(eye(2*J+1),Iy)
Fz=kron(Jz,eye(2*I+1))+kron(eye(2*J+1),Iz)
'''


Jvals=getJ(0)
Jx=Jvals[0]
Jy=Jvals[1]
Jz=Jvals[2]

Fx=kron(Jx,eye(2*I+1))+kron(eye(2*J+1),Ix)
Fy=kron(Jy,eye(2*I+1))+kron(eye(2*J+1),Iy)
Fz=kron(Jz,eye(2*I+1))+kron(eye(2*J+1),Iz)

F_2=Fx.dot(Fx)+Fy.dot(Fy)+Fz.dot(Fz)
F_2=F_2/(hbar**2)

print(F_2/hbar)