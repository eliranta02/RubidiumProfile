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
N = 11

def init_states():
    state1 = q_state.State(N,0,2, is_Ground = True)
    state1.set_S(0.5)
    state2 = q_state.State(N,1,3, is_Ground = True)

    #probe
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

    # pump

    state7 = q_state.State(N,2,1)
    state7.set_L(1)
    state7.set_S(0.5)
    state8 = q_state.State(N,3,2)
    state8.set_L(1)
    state8.set_S(0.5)
    state9 = q_state.State(N,4,3)
    state9.set_L(1)
    state9.set_S(0.5)
    state10 = q_state.State(N,5,4)
    state10.set_L(1)
    state10.set_S(0.5)

    return state1, state2, state3, state4, state5, state6, state7, state8, state9, state10

def init():
    global states
    states = init_states()
    return True

#--------------------------------------------------------------------------------------------------------#

def H(delta_pr, delta_pu, omega_pr, omega_pu):
    global states
    a1, a2, a3, a4, a5, a6, a7, a8, a9, a10 = states

    rho11 = a1 * a1
    rho22 = a2 * a2
    rho33 = a3 * a3
    rho44 = a4 * a4
    rho55 = a5 * a5
    rho66 = a6 * a6
    rho77 = a7 * a7
    rho88 = a8 * a8
    rho99 = a9 * a9
    rho1010 = a10 * a10

    #probe
    rho13 = a1 * a3
    rho14 = a1 * a4
    rho15 = a1 * a5

    rho24 = a2 * a4
    rho25 = a2 * a5
    rho26 = a2 * a6

    # pump
    rho17 = a1 * a7
    rho18 = a1 * a8
    rho19 = a1 * a9

    rho28 = a2 * a8
    rho29 = a2 * a9
    rho210 = a2 * a10

    # probe
    H0 = (DELTA_1) * rho11 + (DELTA_2) * rho22 + (delta_pr + DELTA_3) * rho33 + (delta_pr + DELTA_4) * rho44\
        + (delta_pr + DELTA_5) * rho55 + (delta_pr + DELTA_5) * rho66

    # pump
    H0 = (delta_pu + DELTA_3) * rho77 + (delta_pu + DELTA_4) * rho88 \
         + (delta_pu + DELTA_5) * rho99 + (delta_pu + DELTA_5) * rho1010

    #probe
    omega_pr_13 = transition_strength(a1,a3) * omega_pr
    omega_pr_14 = transition_strength(a1,a4) * omega_pr
    omega_pr_15 = transition_strength(a1,a5) * omega_pr

    omega_pr_24 = transition_strength(a2, a4) * omega_pr
    omega_pr_25 = transition_strength(a2, a5) * omega_pr
    omega_pr_26 = transition_strength(a2, a6) * omega_pr

    # pump
    omega_pu_17 = transition_strength(a1, a7) * omega_pu
    omega_pu_18 = transition_strength(a1, a8) * omega_pu
    omega_pu_19 = transition_strength(a1, a9) * omega_pu

    omega_pu_28 = transition_strength(a2, a8) * omega_pu
    omega_pu_29 = transition_strength(a2, a9) * omega_pu
    omega_pu_210 = transition_strength(a2, a10) * omega_pu

    #probe
    V = (omega_pr_13)  * rho13 + (omega_pr_14) * rho14 + (omega_pr_15) * rho15 + (omega_pr_24)  * rho24 + (omega_pr_25) * rho25 \
        + (omega_pr_26) * rho26
    V += transpose(V)

    #pump
    V = (omega_pu_17) * rho17 + (omega_pu_18) * rho18 + (omega_pu_19) * rho19 + (omega_pu_28) * rho28 + (
        omega_pu_29) * rho29 + (omega_pu_210) * rho210
    V += transpose(V)
    return H0 + V

def decay_martrix(gamma_val):
    a1, a2, a3, a4, a5, a6, a7, a8, a9, a10 = states

    rho33 = a3 * a3
    rho44 = a4 * a4
    rho55 = a5 * a5
    rho66 = a6 * a6
    rho77 = a7 * a7
    rho88 = a8 * a8
    rho99 = a9 * a9
    rho1010 = a10 * a10

    L = gamma_val * (rho33 + rho44 + rho55 + rho66 + rho77 + rho88 + rho99 + rho1010)

    return L

def repopulation_decay_matrix_vee(gamma_val):
    a1, a2, a3, a4, a5, a6, a7, a8, a9, a10 = states

    rho11 = a1 * a1
    rho22 = a2 * a2
    rho33 = a3 * a3
    rho44 = a4 * a4
    rho55 = a5 * a5
    rho66 = a6 * a6
    rho77 = a7 * a7
    rho88 = a8 * a8
    rho99 = a9 * a9
    rho1010 = a10 * a10

    gamma_val_13 = transition_strength(a1,a3) * gamma_val
    gamma_val_14 = transition_strength(a1,a4) * gamma_val
    gamma_val_15 = transition_strength(a1,a5) * gamma_val
    gamma_val_24 = transition_strength(a2, a4) *gamma_val
    gamma_val_25 = transition_strength(a2, a5) *gamma_val
    gamma_val_26 = transition_strength(a2, a6) *gamma_val

    gamma_val_17 = transition_strength(a1, a7) * gamma_val
    gamma_val_18 = transition_strength(a1, a8) * gamma_val
    gamma_val_19 = transition_strength(a1, a9) * gamma_val
    gamma_val_28 = transition_strength(a2, a8) * gamma_val
    gamma_val_29 = transition_strength(a2, a9) * gamma_val
    gamma_val_210 = transition_strength(a2, a10) * gamma_val

    ret_val = gamma_val_13 * outer(rho11, rho33) + gamma_val_14 * outer(rho11, rho44) + gamma_val_15 * outer(rho11, rho55) + \
              gamma_val_24 * outer(rho22, rho44) + gamma_val_25 * outer(rho22, rho55) + gamma_val_26 * outer(rho22,rho66)

    ret_val += gamma_val_17 * outer(rho11, rho77) + gamma_val_18 * outer(rho11, rho88) + gamma_val_19 * outer(rho11,
                                                                                                             rho99) + \
              gamma_val_28 * outer(rho22, rho88) + gamma_val_29 * outer(rho22, rho99) + gamma_val_210* outer(rho22,
                                                                                                             rho1010)

    return ret_val
#--------------------------------------------------------------------------------------------------------#


def callback(param):
    k_probe = 2 * pi /(780e-9)
    (del_val, velocity) = param
    omega_pr = 2 * pi * 0.5e5
    gamma_val = 2 * pi * 6.06e6
    delta_pr = del_val
    delta_pu = delta_pr
    omega_pu = 2 * pi * 5e6

    ret_val = buildRhoMatrix(H(delta_pr - k_probe * velocity, delta_pu + k_probe * velocity, omega_pr, omega_pu),N) + buildGammaMatrix(decay_martrix(gamma_val), N)
    ret_val += repopulation_decay_matrix_vee(gamma_val)

    return ret_val

if __name__ == "__main__":
    # init all states
    init()

    states_name = qunatum_states_dictionary.rhoMatrixNames(N)

    rho11 = states_name.getLocationByName('rho11')
    rho22 = states_name.getLocationByName('rho22')

    rho13 = states_name.getLocationByName('rho13')
    rho14 = states_name.getLocationByName('rho14')
    rho15 = states_name.getLocationByName('rho15')

    rho24 = states_name.getLocationByName('rho24')
    rho25 = states_name.getLocationByName('rho25')
    rho26 = states_name.getLocationByName('rho26')

    y0 = zeros((N * N, 1))
    y0[rho11] = 0.5
    y0[rho22] = 0.5

    lmes = Linblad_master_equation_solver(False)

    returnDic = [rho13, rho14, rho15, rho24, rho25, rho26]

    running_param = linspace(-2 * pi * 4e9, 2 * pi * 4e9, 500) # (frequency scaning) detuning array
    v_param = linspace(-600, 600, 200)  # atomic velocities array
    time_val = 1
    Tc = 25
    results = lmes.solve_master_equation_with_Doppler_effect(callback, running_param,v_param, y0, time_val, Tc, returnDic)


    absor = []

    for idx, _ in enumerate(running_param):
        absor.append(results[rho13][idx].imag + results[rho14][idx].imag + results[rho15][idx].imag + results[rho24][idx].imag + results[rho25][idx].imag + \
                       results[rho26][idx].imag)

    plt.plot(running_param / (2 * pi), absor)
    plt.show()

