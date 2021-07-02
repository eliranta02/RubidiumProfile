from __future__ import division
from numpy import linalg as LA
from J_matrix import *
from S_matrix import *
from L_matrix import *

S=0.5
L=1

LS=0.5*(J_2-kron(L_2,eye(2*S+1))-kron(S_2,eye(2*L+1)))

eigVal,eigVec=LA.eig(LS)

