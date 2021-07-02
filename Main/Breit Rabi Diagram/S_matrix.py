from __future__ import division
from scipy.constants import hbar
from numpy import *


Sx=(hbar/2)*array([[0,1],[1,0]])
Sy=(hbar/2)*array([[0,-1j],[1j,0]])
Sz=(hbar/2)*array([[1,0],[0,-1]])

S_2=Sx.dot(Sx)+Sy.dot(Sy)+Sz.dot(Sz)

