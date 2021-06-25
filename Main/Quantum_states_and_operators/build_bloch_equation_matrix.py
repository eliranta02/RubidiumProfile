from __future__ import division, print_function
from numpy import kron, eye, transpose


def buildGammaMatrix(g_matrix, N):
    '''
    Lp = {L, p}
    :param g_matrix:
    :param N:
    :return:
    '''
    LR = kron(g_matrix, eye(N))
    RL = kron(eye(N), transpose(g_matrix))
    ret_val = -1/2 * (LR + RL)
    return ret_val

def buildRhoMatrix(H, N):
    '''
        [H, p]
    :param H: hamiltonian
    :param N: number of levels
    :param L: decay matrix
    :return:
    '''
    Hrho = kron(H, eye(N))
    rhoH = kron(eye(N), transpose(H))
    ret_val = -1j * (Hrho - rhoH)
    return ret_val

#def build