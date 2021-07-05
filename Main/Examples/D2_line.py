from __future__ import division, print_function
from __future__ import division, print_function
from Main.Quantum_states_and_operators import q_state, qunatum_states_dictionary
from Main.Quantum_states_and_operators.build_bloch_equation_matrix import *
from Main.Quantum_states_and_operators.frequency_strength_calc import *
from Main.Quantum_states_and_operators.Linblad_master_equation_solver import  *
from Main.Constants.Rb_constants import *
import pylab as plt

#Constants
DELTA_1 = -2 * pi * 1.770e9
DELTA_2 = 2 * pi * 1.264e9
DELTA_3 = -2 * pi * 113.208e6
DELTA_4 = -2 * pi * 83.835e6
DELTA_5 = -2 * pi * 20.435e6
DELTA_6 = 2 * pi * 100.205e6



states = None
N = 6

def init_states():
    state1 = q_state.State(N,0,2, is_Ground = True)
    state1.set_S(0.5)
    state2 = q_state.State(N,1,3, is_Ground = True)
    state2.set_S(0.5)
    state3 = q_state.State(N,2,1)
    state3.set_L(1)
    state3.set_S(0.5)
    state4 = q_state.State(N,3,2)
    state4.set_L(1)
    state4.set_S(0.5)
    state5 = q_state.State(N,4,3)
    state5.set_L(1)
    state5.set_S(0.5)
    state6 = q_state.State(N,5,4)
    state6.set_L(1)
    state6.set_S(0.5)
    return state1, state2, state3, state4, state5, state6

def init():
    global states
    states = init_states()
    return True

#-----------------------------------------------------------------------#
def H(delta_pr, omega_pr):
    global states
    a1, a2, a3, a4, a5, a6 = states

    rho11 = a1 * a1
    rho22 = a2 * a2
    rho33 = a3 * a3
    rho44 = a4 * a4
    rho55 = a5 * a5
    rho66 = a6 * a6

    rho13 = a1 * a3
    rho14 = a1 * a4
    rho15 = a1 * a5

    rho24 = a2 * a4
    rho25 = a2 * a5
    rho26 = a2 * a6

    H0 = (DELTA_1) * rho11 + (DELTA_2) * rho22 + (delta_pr + DELTA_3) * rho33 + (delta_pr + DELTA_4) * rho44\
        + (delta_pr + DELTA_5) * rho55 + (delta_pr + DELTA_5) * rho66

    V = (omega_pr)  * rho13 + (omega_pr) * rho14 + (omega_pr) * rho15 + (omega_pr)  * rho24 + (omega_pr) * rho25 + (omega_pr) * rho26
    V += transpose(V)
    return H0 + V

init()

state1, state2, state3, state4, state5, state6 = states


print(transition_strength(state1,state3))
