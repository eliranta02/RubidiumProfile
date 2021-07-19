from __future__ import division, print_function
from __future__ import division, print_function
from Main.Quantum_states_and_operators import q_state, qunatum_states_dictionary
from Main.Quantum_states_and_operators.build_bloch_equation_matrix import *
from Main.Quantum_states_and_operators.frequency_strength_calc import *
from Main.Quantum_states_and_operators.Linblad_master_equation_solver import  *
from Main.Constants.Rb_constants import *

from Main.Operators.tools import *
import pylab as plt

states = None
N = 4
HFS = 3.035e9

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

# All the calculation is based on

def H(delta, Delta_1, omega_c, omega_p, omega_1, omega_2):
    global states
    a1, a2, a3, a4 = states

    rho22 = a2 * a2
    rho33 = a3 * a3
    rho44 = a4 * a4

    rho13 = a1 * a3
    rho14 = a1 * a4

    rho23 = a2 * a3
    rho24 = a2 * a4

    Delta_2 = Delta_1 + HFS/2 + delta

    H0 = -HFS/2 * rho11 + (HFS/2 + delta) * rho22 + (Delta_1) * rho33 + (Delta_2) * rho44

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
    ret_val = 0.5 * gamma1 * outer(rho11, rho33) + 0.5 * gamma2 * outer(rho11, rho44) + \
              0.5 * gamma2 * outer(rho22, rho33) + 0.5 * gamma2 * outer(rho22, rho44)
    return ret_val

#--------------------------------------------------------------------------------------------------------#

def callback(param):
    k_pump = 2 * pi / (795e-9)
    (del_val, velocity) = param
    delta = del_val- (k_pump * velocity)/(2 * pi)
    gamma1 = gamma2 =  5.7e6
    Delta_1 =  800e6
    omega_c = 0
    omega_p = 0
    omega_1 =  430e6
    omega_2 =  430e6

    ret_val = buildRhoMatrix(H(delta, Delta_1, omega_c, omega_p, omega_1, omega_2), N) + buildGammaMatrix(
        decay_martrix(gamma1, gamma2), N) + repopulation_decay_matrix(gamma1, gamma2)

    return ret_val


if __name__ == "__main__":
    init()

    states_name = qunatum_states_dictionary.rhoMatrixNames(N)
    rho11 = states_name.getLocationByName('rho11')
    rho22 = states_name.getLocationByName('rho22')
    rho23 = states_name.getLocationByName('rho23')
    rho14 = states_name.getLocationByName('rho14')

    y0 = zeros((N * N, 1))
    y0[rho11] = 0.5
    y0[rho22] = 0.5

    lmes = Linblad_master_equation_solver(False)

    returnDic = [rho23,rho14]

    running_param = linspace( -6e9, 6e9, 2500)  # (frequency scaning) detuning array
    v_param = linspace(-800, 800, 1000)  # atomic velocities array
    #time_val = 1
    Tc = 100
    #callback, detuning_param, velocity_param, N, Tc, keys
    results = lmes.solve_master_equation_steady_state_with_Doppler_effect(callback, running_param, v_param, N, Tc, returnDic)
    #results = lmes.solve_master_equation_steady_state_without_Doppler_effect(callback, running_param, N, returnDic)

    solution = [res.imag for res in results[rho23]]
    solution1 = [res.imag for res in results[rho14]]

    val = []
    for idx, _ in enumerate(solution):
        val.append(solution1[idx] + solution[idx])
    plt.figure(1)
    plt.plot(running_param , solution,label='probe')
    plt.plot(running_param , solution1,label='conj.')
    plt.legend(loc=0)
    plt.figure(2)
    plt.plot(running_param, val)

    plt.show()
