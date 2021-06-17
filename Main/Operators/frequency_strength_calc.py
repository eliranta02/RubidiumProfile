from __future__ import division, print_function
from Wigner_Calc import *
from numpy import sqrt, pi
from Main.Constants.Rb_constants import *
from scipy.constants import c, hbar

#calculate the rabi frequency of define transition. the calculation is based on aritcle "On the consistency of Rabi frequency calculations"
def rabi_freq(I,lhs,rhs):
    # TODO maybe define it as a calss of state
    f1 = lhs[0]
    m_f1 = lhs[1]
    j_f1= lhs[2]
    l_f1 = lhs[2]

    f2 = rhs[0]
    mf2= rhs[1]
    j_f1 = rhs[2]

    ret_val = sqrt(3 * wavelen**3 * gamma * I / (4 * pi**2 * c * hbar)) * sqrt()

    return ret_val


