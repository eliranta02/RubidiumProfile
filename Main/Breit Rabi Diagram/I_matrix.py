from __future__ import division
from scipy.constants import hbar
from numpy import *


Ix=(0.5*hbar)*array([[0,sqrt(3),0,0],[sqrt(3),0,2,0],[0,2,0,sqrt(3)],[0,0,sqrt(3),0]])
Iy=(0.5*hbar)*array([[0,-sqrt(3)*1j,0,0],[sqrt(3)*1j,0,-2j,0],[0,2j,0,-sqrt(3)*1j],[0,0,sqrt(3)*1j,0]])
Iz=(0.5*hbar)*array([[3,0,0,0],[0,1,0,0],[0,0,-1,0],[0,0,0,-3]])

I_2=Ix.dot(Ix)+Iy.dot(Iy)+Iz.dot(Iz)


