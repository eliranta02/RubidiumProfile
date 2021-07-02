from __future__ import division, print_function
from Main.Quantum_states_and_operators import q_state, qunatum_states_dictionary
from Main.Quantum_states_and_operators.build_bloch_equation_matrix import *
from Main.Quantum_states_and_operators.Linblad_master_equation_solver import  *

from Main.Constants.Rb_constants import *

import pylab as plt

VEE_TYPE = 1
LAMBDA_TYPE = 2
LADDER_TYPE = 3

DIAGRAM_LEVEL_TYPE = VEE_TYPE
states = None
N = 0


def init_states():
    state1 = q_state.State(3,0,1)
    state2 = q_state.State(3,1,2)
    state3 = q_state.State(3,2,3)
    return state1, state2, state3

def init():
    global states, N
    states = init_states()
    N = 3
    return 1
#-----------------------------------------------------------------------#

def H_vee(delta_pr, delta_pu, omega_pr, omega_pu):
    '''
            |2> -----       |3> -----
                 ||          ||
                   ||       ||
                |1> -----

    :param delta_pr:
    :param delta_pu:
    :param omega_pr:
    :param omega_pu:
    :return:
    '''
    global states
    a1, a2, a3 = states

    rho22 = a2 *a2
    rho12 = a1 * a2
    rho21 = a2 * a1

    rho33 = a3 * a3
    rho13 = a1 * a3
    rho31 = a3 * a1

    H = (delta_pr) * rho22 + (delta_pu) * rho33 + 0.5 * omega_pr * rho12 + 0.5 * omega_pr * rho21+ \
        0.5 * omega_pu * rho13 + 0.5 * omega_pu * rho31
    return H

def decay_martrix_vee(gamma2, gamma3):

    a1, a2, a3 = states
    rho22 = a2 * a2
    rho33 = a3 * a3
    L = gamma2 * rho22 + gamma3 * rho33
    return L

def repopulation_decay_matrix_vee(gamma2, gamma3):
    a1, a2, a3 = states
    rho11 = a1 * a1
    rho22 = a2 * a2
    rho33 = a3 * a3
    ret_val = gamma3 * outer(rho11, rho33) + gamma2 * outer(rho11, rho22)
    return ret_val

#-----------------------------------------------------------------------#
def H_Lambda(delta_pr, delta_pu, omega_pr, omega_pu):
    '''
                   |3> -----
                    ||    ||
                   ||       ||
                |1> -----   |2> -----

    :param delta_pr:
    :param delta_pu:
    :param omega_pr:
    :param omega_pu:
    :return:
    '''
    global states

    a1, a2, a3 = states

    rho22 = a2 * a2
    rho33 = a3 * a3
    rho13 = a1 * a3
    rho31 = a3 * a1

    rho23 = a2 * a3
    rho32 = a3 * a2

    H = (delta_pr) * rho33 + (delta_pr - delta_pu) * rho22 + 0.5 * omega_pr * rho13 + 0.5 * omega_pr * rho31 + \
        0.5 * omega_pu * rho23 + 0.5 * omega_pu * rho32
    return H

def decay_martrix_Lambda(gamma3):

    a1, a2, a3 = states
    rho33 = a3 * a3
    L = gamma3 * rho33
    return L

def repopulation_decay_matrix_Lambda(gamma3):
    a1, a2, a3 = states
    rho11 = a1 * a1
    rho22 = a2 * a2
    rho33 = a3 * a3
    ret_val = gamma3 * outer(rho11, rho33) + gamma3 * outer(rho22, rho33)
    return ret_val

#-----------------------------------------------------------------------#
def H_Ladder(delta_pr, delta_pu, omega_pr, omega_pu):
    global states

    a1, a2, a3 = states

    rho22 = a2 * a2
    rho33 = a3 * a3

    rho12 = a1 * a2
    rho21 = a2 * a1

    rho23 = a2 * a3
    rho32 = a3 * a2

    H = (delta_pu) * rho22 + (delta_pr + delta_pu) * rho33 + 0.5 * omega_pr * rho23 + 0.5 * omega_pr * rho32 + \
        0.5 * omega_pu * rho12 + 0.5 * omega_pu * rho21

    return H

def decay_martrix_Ladder(gamma2, gamma3):

    a1, a2, a3 = states
    rho33 = a3 * a3
    rho22 = a2 * a2
    L = gamma3 * rho33 + gamma2 * rho22
    return L

def repopulation_decay_matrix_Ladder(gamma2, gamma3):
    a1, a2, a3 = states
    rho11 = a1 * a1
    rho22 = a2 * a2
    rho33 = a3 * a3
    ret_val = gamma2 * outer(rho11, rho22) + gamma3 * outer(rho22, rho33)
    return ret_val

#-----------------------------------------------------------------------#

def callback(param):
    k_pump = 795e-9
    k_probe = 780e-9
    omegaProbe = 2 * pi * 10e5
    omegaPump = 2 * pi * 1e6
    delta_pu = 0
    (del_val, velocity) = param

    if DIAGRAM_LEVEL_TYPE == VEE_TYPE:
        gamma2, gamma3 = 2 * pi * 6.06e6, 2 * pi * 5.75e6
        ret_val = buildRhoMatrix(H_vee(del_val-k_probe * velocity, delta_pu + k_pump * velocity, omegaProbe, omegaPump), N) + buildGammaMatrix(decay_martrix_vee(gamma2, gamma3), N)
        ret_val += repopulation_decay_matrix_vee(gamma2, gamma3)

    if DIAGRAM_LEVEL_TYPE == LAMBDA_TYPE:
        gamma3 = 2 * pi * 6.06e6
        ret_val = buildRhoMatrix(H_Lambda(param, delta_pu, omegaProbe, omegaPump), N) + buildGammaMatrix(
            decay_martrix_Lambda(gamma3), N)
        ret_val += repopulation_decay_matrix_Lambda(gamma3)

    if DIAGRAM_LEVEL_TYPE == LADDER_TYPE:
        gamma2 = 2 * pi * 6.06e6
        gamma3 = 2 * pi * 1.8e6
        ret_val = buildRhoMatrix(H_Ladder(param, delta_pu, omegaProbe, omegaPump), N) + buildGammaMatrix(
            decay_martrix_Ladder(gamma2, gamma3), N)
        ret_val += repopulation_decay_matrix_Ladder(gamma2, gamma3)

    return ret_val


if __name__ == "__main__":
    # init all states
    init()

    states_name = qunatum_states_dictionary.rhoMatrixNames(N)
    rho11 = states_name.getLocationByName('rho11')
    rho12 = states_name.getLocationByName('rho12')
    rho22 = states_name.getLocationByName('rho22')
    print(states_name)

    if DIAGRAM_LEVEL_TYPE == VEE_TYPE:

        y0 = zeros((N * N, 1))
        y0[rho11] = 0.5
        y0[rho22] = 0.5

        lmes = Linblad_master_equation_solver(False)

        returnDic = [rho12, rho22]

        running_param = linspace(-2 * pi * 0.1e9, 2 * pi * 0.1e9, 500)  # (frequency scaning) detuning array
        v_param = linspace(-600, 600, 500)  # atomic velocities array
        time_val = 0.1
        Tc = 50
        results = lmes.solve_master_equation_with_Doppler_effect(callback, running_param, v_param, y0, time_val,
                                                                 Tc, returnDic)
        solution = [res.real for res in results[rho22]]
        plt.plot(running_param, solution)
        plt.show()


    if DIAGRAM_LEVEL_TYPE == LAMBDA_TYPE:
        pass

    if DIAGRAM_LEVEL_TYPE == LADDER_TYPE:
        pass
