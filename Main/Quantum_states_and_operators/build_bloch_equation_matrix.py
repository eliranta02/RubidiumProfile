from __future__ import division, print_function
from numpy import kron, eye, transpose

def buildRhoMatrix(H, N, L):
    '''
        dp/dt = i/h [H, p] + Lp
    :param H: hamiltonian
    :param N: number of levels
    :param L: decay matrix
    :return:
    '''
    Hrho = kron(H, eye(N))
    rhoH = kron(eye(N), transpose(H))
    ret_val = -1j * (Hrho - rhoH) + L
    return ret_val

