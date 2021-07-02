from __future__ import division, print_function
from Main.Quantum_states_and_operators import q_state, qunatum_states_dictionary
from Main.Quantum_states_and_operators.build_bloch_equation_matrix import *
from Main.Quantum_states_and_operators.Linblad_master_equation_solver import  *

from Main.Constants.Rb_constants import *

import pylab as plt

states = None
N = 2

def init_states():
    state1 = q_state.State(2,0,1)
    state2 = q_state.State(2,1,2)
    return state1, state2


def H(delta,omega):
    global states
    if states == None:
        states = init_states()

    a1, a2 = states

    rho22 = a2 *a2
    rho12 = a1 * a2
    rho21 = a2 * a1
    H = (delta) * rho22 + 0.5 * omega * rho12 + 0.5 * omega * rho21
    return H

def decay_martrix(gamma):
    global states
    if states == None:
        states = init_states()

    a1, a2 = states
    rho22 = a2 * a2
    L = gamma * rho22
    return L


def callback(param):
    omegaProbe = 2 * pi * 0.5e5
    ret_val = buildRhoMatrix(H(param, omegaProbe), N) + buildGammaMatrix(decay_martrix(gamma), N)
    return ret_val

if __name__ == "__main__":
    states_name = qunatum_states_dictionary.rhoMatrixNames(N)
    rho11 = states_name.getLocationByName('rho11')
    rho12 = states_name.getLocationByName('rho12')
    rho22 = states_name.getLocationByName('rho22')

    y0 = zeros((N * N, 1))
    y0[rho11] = 1


    lmes = Linblad_master_equation_solver(False)

    returnDic = [rho12, rho22]

    running_param = linspace(-2 * pi * 2e9, 2 * pi * 2e9, 300) #(frequency scaning) detuning array
    v_param = linspace(-600, 600, 500) #atomic velocities array
    time_val = 0.1
    k_wave = k_num
    Tc = 50
    results = lmes.solve_master_equation_with_Doppler_effect(callback, running_param, v_param, y0, time_val, k_wave, Tc, returnDic)

    solution = [res.real for res in results[rho22]]
    plt.plot(running_param, solution)
    plt.show()

