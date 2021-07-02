from __future__ import division
from numpy import *
from L_matrix import *
from S_matrix import *

S=0.5
L=1



Jx=kron(Lx,eye(2*S+1))+kron(eye(2*L+1),Sx)
Jy=kron(Ly,eye(2*S+1))+kron(eye(2*L+1),Sy)
Jz=kron(Lz,eye(2*S+1))+kron(eye(2*L+1),Sz)

J_2=Jx.dot(Jx)+Jy.dot(Jy)+Jz.dot(Jz)




def getJ(L):
    if L==0:
        Jx=kron(eye(2*L+1),Sx)
        Jy=kron(eye(2*L+1),Sy)
        Jz=kron(eye(2*L+1),Sz)
    elif L==1:
        Jx=kron(Lx,eye(2*S+1))+kron(eye(2*L+1),Sx)
        Jy=kron(Ly,eye(2*S+1))+kron(eye(2*L+1),Sy)
        Jz=kron(Lz,eye(2*S+1))+kron(eye(2*L+1),Sz)

    return Jx,Jy,Jz






