from __future__ import division
from scipy.constants import hbar
from numpy import *


Lx=(sqrt(0.5)*hbar)*array([[0,1,0],[1,0,1],[0,1,0]])
Ly=(sqrt(0.5)*hbar)*array([[0,-1j,0],[1j,0,-1j],[0,1j,0]])
Lz=hbar*array([[1,0,0],[0,0,0],[0,0,-1]])


L_2=Lx.dot(Lx)+Ly.dot(Ly)+Lz.dot(Lz)

