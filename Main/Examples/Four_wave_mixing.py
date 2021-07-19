from __future__ import division, print_function
from __future__ import division, print_function
from Main.Quantum_states_and_operators import q_state, qunatum_states_dictionary
from Main.Quantum_states_and_operators.build_bloch_equation_matrix import *
from Main.Quantum_states_and_operators.frequency_strength_calc import *
from Main.Quantum_states_and_operators.Linblad_master_equation_solver import  *
from Main.Constants.Rb_constants import *
import pylab as plt

states = None
N = 4

def init_states():
    state1 = q_state.State(N,0,2, is_Ground = True)
    state2 = q_state.State(N,1,3, is_Ground = True)
    state3 = q_state.State(N,2,2)
    state4 = q_state.State(N,3,3)
    return state1, state2, state3, state4

def init():
    global states
    states = init_states()
    return True

#--------------------------------------------------------------------------------------------------------#

def H(delta, Delta_1, Delta_2, omega_c, omega_p, omega_1, omega_2):
    global states
    a1, a2, a3, a4 = states

    rho11 = a1 * a1
    rho22 = a2 * a2
    rho33 = a3 * a3
    rho44 = a4 * a4

    rho13 = a1 * a3
    rho14 = a1 * a4

    rho23 = a2 * a3
    rho24 = a2 * a4


    H0 = (delta) * rho22 + (Delta_1) * rho33 + (Delta_2) * rho44

    V = (omega_1)  * rho13 + (omega_c) * rho14 + (omega_p) * rho23 + (omega_2) * rho24
    V += transpose(V)
    return H0 + V

def decay_martrix(gamma1, gamma2):

    a1, a2, a3, a4 = states

    rho33 = a3 * a3
    rho44 = a4 * a4

    L = gamma1 * rho33 + gamma2 * rho44

    return L

def repopulation_decay_matrix(gamma1, gamma2):
    a1, a2, a3, a4 = states
    rho11 = a1 * a1
    rho22 = a2 * a2
    rho33 = a3 * a3
    rho44 = a4 * a4
    ret_val = 0.5 * gamma1 * outer(rho11, rho33) + 0.5 * gamma1 * outer(rho22, rho44) + \
              0.5 * gamma2 * outer(rho22, rho33) + 0.5 * gamma2 * outer(rho22, rho44)
    return ret_val

#--------------------------------------------------------------------------------------------------------#

def callback(param):
    ret_val  = 1
    return ret_val


if __name__ == "__main__":
    init()
    delta = 1
    Delta_1 = 1
    Delta_2 = 1
    omega_c = 1
    omega_c = 1
    omega_1 = 1
    omega_2 = 1
    a = H(delta, Delta_1, Delta_2, omega_c, omega_c, omega_1, omega_2)
    print(a)