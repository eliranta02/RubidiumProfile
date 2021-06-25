from __future__ import division, print_function
from Main.Quantum_states_and_operators import q_state, qunatum_states_dictionary
from Main.Quantum_states_and_operators.build_bloch_equation_matrix import *
from Main.Quantum_states_and_operators.Linblad_master_equation_solver import  *

states = None

def init_states():
    state1 = q_state.State(2,0,1)
    state2 = q_state.State(2,1,2)
    return state1, state2


def H(delta,omega,gamma):
    global states
    if states == None:
        states = init_states()

    a1, a2 = states

    rho22 = a2 *a2
    rho12 = a1 * a2
    rho21 = a2 * a1
    H = (delta - 1j * gamma)* rho22 + 0.5 * omega * rho12 + 0.5 * omega * rho21
    return H


def callback(param):
    ret_val = buildRhoMatrix(H(param, 1, 0.9), 2)
    return ret_val

if __name__ == "__main__":
    start = time.perf_counter()
    N = 2
    states_name = qunatum_states_dictionary.rhoMatrixNames(N)
    rho11 = states_name.getLocationByName('rho11')
    rho12 = states_name.getLocationByName('rho12')
    rho22 = states_name.getLocationByName('rho22')

    y0 = zeros((N * N,1))
    y0[rho11] = 1


    temp = Linblad_master_equation_solver(False)

    returnDic = {rho12: [], rho22: []}

    running_param = linspace(-100, 100, 10000)

    results = temp.solve_master_equation_without_Doppler_effect(callback, running_param, y0, returnDic)

    finish = time.perf_counter()

    print(f'Finished in {round(finish - start, 2)} second(s)')

    import pylab as plt
    plt.plot(results[rho22])
    plt.show()
    '''
    temp = Linblad_master_equation_solver(False)
    returnDic = {rho12: [], rho22: [], rho11: []}
    time_arr = linspace(0,100,100)
    t, dic = temp.timeDependentSolver(callback(0), y0, time_arr, returnDic)
    print(t)

    import pylab as plt

    plt.plot(dic[rho11])
    plt.plot(dic[rho22])
    plt.show()
    '''

