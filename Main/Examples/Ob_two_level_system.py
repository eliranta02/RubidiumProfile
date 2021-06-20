from __future__ import division, print_function
from Main.Quantum_states_and_operators import q_state

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

    print(a1.outer_product(a2))
    rho22 = a2 *a2
    rho12 = a1 * a2
    rho21 = a2 * a1
    H = (delta - 1j * gamma)* rho22 + 0.5 * omega * rho12 + 0.5 * omega * rho21
    return H


if __name__ == "__main__":
    H(1,3,1)

