from __future__ import division, print_function
from Main.Quantum_states_and_operators import q_state, qunatum_states_dictionary
from Main.Quantum_states_and_operators.build_bloch_equation_matrix import *
from Main.Quantum_states_and_operators.Linblad_master_equation_solver import  *

from Main.Constants.Rb_constants import *

from Main.Tools.PlottingTemplate import *


states = None
N = 2

def init_states():
    state1 = q_state.State(2,0,1)
    state2 = q_state.State(2,1,2)
    return state1, state2


def H(delta,omega):
    #global states
    #if states == None:
    #    states = init_states()

    a1, a2 = states

    rho22 = a2 *a2
    rho12 = a1 * a2
    rho21 = a2 * a1
    H = (delta) * rho22 + 0.5 * omega * rho12 + 0.5 * omega * rho21
    return H

def decay_martrix(gamma):
    #global states
    #if states == None:
    #    states = init_states()

    a1, a2 = states
    rho22 = a2 * a2
    L = gamma * rho22
    return L


def callback(param):
    omegaProbe = 2 * pi * 0.5e5
    (del_val, velocity) = param

    global states
    if states == None:
        states = init_states()

    a1, a2 = states
    rho11 = a1*a1
    rho22 = a2*a2
    k_wave = k_num
    ret_val = buildRhoMatrix(H(del_val-k_wave * velocity, omegaProbe), N) + buildGammaMatrix(decay_martrix(gamma), N) + gamma * outer(rho11 ,rho22)
    return ret_val

def time_dependent_TLS():

    global states
    if states == None:
        states = init_states()

    lmes = Linblad_master_equation_solver(False)

    states_name = qunatum_states_dictionary.rhoMatrixNames(N)

    rho11_num = states_name.getLocationByName('rho11')
    rho12_num = states_name.getLocationByName('rho12')
    rho22_num = states_name.getLocationByName('rho22')

    y0 = zeros((N * N, 1))
    y0[rho11_num] = 1

    a1, a2 = states
    rho11 = a1 * a1
    rho22 = a2 * a2

    returnDic = [rho12_num, rho22_num,rho11_num]

    time_arr = linspace(0,0.2e-6,1000)
    omegaProbe = 2 * pi * 50e6
    g = gamma + 2 * pi * 10e6
    matrix_val = buildRhoMatrix(H(0, omegaProbe), N) + buildGammaMatrix(decay_martrix(g), N) + g * outer(rho11 ,rho22)

    ret_val = lmes.solve_density_matrix_evolution(matrix_val, y0, time_arr, returnDic)
    solution = [res.imag for res in ret_val[rho12_num]]
    solution1 = [res.real for res in ret_val[rho22_num]]
    #plt.plot(time_arr, solution)
    plt.figure(1,figsize=(10, 8))
    plt.subplot(121)
    plt.plot(time_arr*200e6, solution1,lw=2,color=d_red)
    plt.xlim([0,35])
    plt.xlabel('L ($\mu m$)')
    plt.ylabel(r'Population ($\rho_{22}$)')
    plt.subplot(122)
    plt.plot(time_arr * 200e6, solution, lw=2, color=d_blue)
    plt.xlim([0, 35])
    plt.xlabel('L ($\mu m$)')
    plt.ylabel(r'Population ($\rho_{12}$)')

    g = gamma
    matrix_val = buildRhoMatrix(H(0, omegaProbe), N) + buildGammaMatrix(decay_martrix(g), N) + g * outer(rho11, rho22)

    ret_val = lmes.solve_density_matrix_evolution(matrix_val, y0, time_arr, returnDic)
    solution = [res.real for res in ret_val[rho11_num]]
    solution1 = [res.real for res in ret_val[rho22_num]]
    #plt.plot(time_arr, solution)
    #plt.plot(time_arr*1e5, solution1,color=d_blue,lw=2)
    plt.savefig('TLS.png')
    plt.show()

def run():
    states_name = qunatum_states_dictionary.rhoMatrixNames(N)

    rho11 = states_name.getLocationByName('rho11')
    rho12 = states_name.getLocationByName('rho12')
    rho22 = states_name.getLocationByName('rho22')

    y0 = zeros((N * N, 1))
    y0[rho11] = 1

    lmes = Linblad_master_equation_solver(False)

    returnDic = [rho12, rho22]

    running_param = linspace(-2 * pi * 2.5e9, 2 * pi * 2.5e9, 3000) #(frequency scaning) detuning array
    v_param = linspace(-600, 600, 500) #atomic velocities array
    time_val = 0.1
    k_wave = k_num
    Tc = 50
    #results = lmes.solve_master_equation_with_Doppler_effect(callback, running_param, v_param, y0, time_val, Tc, returnDic)
    #results = lmes.solve_master_equation_steady_state_without_Doppler_effect(callback, running_param, N, returnDic)
    results = lmes.solve_master_equation_steady_state_with_Doppler_effect(callback, running_param, v_param, N, Tc, returnDic)
    plt.figure(1,figsize=(5,4))
    solution = [res.real for res in results[rho12]]
    plt.plot(running_param/(2 * pi*1e9), solution,color=d_blue,lw=2,label=r'$\chi_R$')

    solution = [res.imag for res in results[rho12]]
    plt.plot(running_param / (2 * pi*1e9), solution,color=d_red,lw=2,label=r'$\chi_I$')

    plt.xlabel('Detunning (GHz)')
    plt.yticks([])
    plt.legend(loc=0)
    plt.tight_layout()
    plt.savefig('abs.png')
    plt.show()

if __name__ == "__main__":
    run()

