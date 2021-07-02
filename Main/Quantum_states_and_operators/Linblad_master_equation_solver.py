from __future__ import division, print_function

from numpy import zeros, transpose, kron, sqrt, dot, array, linspace
from numpy.linalg import solve, lstsq
from pylab import eig,mat,inv,exp,diag
from scipy.integrate import ode, simps

import concurrent.futures
import functools

from Main.Constants.Velocities_distribution import *

from Main.Unit_converters.Rb_unit_converter import *
from Main.Unit_converters.Global_unit_converter import *
from tqdm import tqdm


class Ode_time_dependent_solver(object):

    def timeDependentSolver(self,matrix,y0,time_arr,keys):
        for time_val in time_arr:
            eval, evec = eig(matrix)
            solT = evec * mat(diag(exp(eval * time_val))) * inv(evec) * y0
            ret_val = {}
            for key in keys:
                ret_val[key] = solT[key].item()
        return time_arr, ret_val

    def solveSteadyState(self,matrix):
        '''
        Returns
        ------
        the solution of steady state solution
        '''

        N = len(matrix[0])
        vec = zeros((N, ))
        for i in range(int(sqrt(N))):
            psi = zeros((int(sqrt(N)), ))
            psi[i] = 1
            vec += transpose(kron(psi, psi))
        matrix[0] = vec
        res = zeros((N,))
        res[0] = 1
        A = lstsq(matrix, res)
        #TODO still needed to check
        return A

    def timeDependentSolverWithParam(self,jac, userSupply, param):
        # initate state
        t0 = 0
        points = param[0]
        t_scan = param[1]
        y0 = param[2]

        H = t_scan / points
        r = ode(userSupply, jac).set_integrator('zvode', method='bdf', with_jacobian=True, rtol=1e-6, order=5)
        tm = param
        r.set_initial_value(y0, t0).set_f_params(tm).set_jac_params(tm)
        x = []
        y = []

        loc = 0
        while r.successful() and r.t < t_scan:
            if (t_scan - r.t) < H:
                H = t_scan - r.t
            r.integrate(r.t + H)
            x.append(r.t)
            a = (r.y).reshape(1, len(r.y))
            y.append(a[0])
            loc += 1

        return x, y

    def odeSolver(self,matrix,y0,time_val,keys):

        '''

        Parameters
        ----------
        y0 : y0 = zeros((N,1)) ; y0[1] = 1
        time_val: runing integration parameters (e.g. time parameter)
        keys : all the return values

        Returns
        -------
        returnDic
        '''
        eval, evec = eig(matrix)
        solT = evec * mat(diag(exp(eval * time_val))) * inv(evec) * y0
        retval = {}
        for key in keys:
            retval[key]= solT[key].item()

        return retval

    @staticmethod
    def print_info():
        print('This class is used in order to solve differential equations.')

    def __str__(self):
        N = len(self.matrix)
        return 'The number of levels : {} \nThe '.format(int(sqrt(N)))

class Linblad_master_equation_solver(Ode_time_dependent_solver):

    def __init__(self, enable_multiprocessing):
        self.is_multi_processing_enabled = enable_multiprocessing

    '''def test_solve(self,callback, delta_array,v_array, y0):
        time_val = 3
        k_wave = 1
        rho22 = zeros((4 ,1))
        rho22[3] = 1
        rho22_del = []

        for idx, delt in enumerate(delta_array):
            returnDic = {3: []}
            rho22_vel = []
            for v in v_array:
                delta = delt - k_wave * v
                matrix = callback(delta)
                #eval, evec = eig(matrix)
                #solT = evec * mat(diag(exp(eval * time_val))) * inv(evec) * y0
                #sol_with_doppler = ((transpose(rho22) * solT)).item()
                #rho22_vel.append(sol_with_doppler)
                sol_with_doppler = self.odeSolver(matrix,y0,time_val,returnDic)

            velocity_dist = [maxwell(param, 300) for param in v_array]
            #product = [a * b for a, b in zip(rho22_vel, velocity_dist)]
            product1 = [a * b for a, b in zip(returnDic[3], velocity_dist)]

            #rho22_del.append(simps(product,v_array))
            rho22_del.append(simps(product1,v_array))
        return rho22_del'''
    def solve_master_equation_without_Doppler_effect(self, callback, detuning_param, y0, time_val, keys):
        '''
        :param callback:
        :param detuning_param:
        :param y0: initial vector  y0 = zeros((N,)) ; y0[1] = 1
        :param returnDic:
        :return:
        '''

        ret_val = {}
        for key in keys:
            ret_val[key] = []

        mat_solver = [callback((param, 0)) for param in detuning_param]

        if self.is_multi_processing_enabled == True:
            with concurrent.futures.ProcessPoolExecutor() as executor:
                results = executor.map(functools.partial(self.odeSolver, y0 = y0, time_val = time_val, keys = keys), mat_solver)
                temp_list = list(results)
        else:
            results = map(functools.partial(self.odeSolver, y0=y0, time_val = time_val, keys = keys), mat_solver)
            temp_list = list(results)

        for key in keys:
            ret_val[key] = [item[key] for item in temp_list]

        return ret_val

    def solve_master_equation_with_Doppler_effect(self, callback, detuning_param, velocity_param, y0, time_val, Tc, keys):
        ret_val = {}
        for key in keys:
            ret_val[key] = []

        velocity_dist = [maxwell(param, temp2velocity(celsius2kelvin(Tc))) for param in velocity_param]

        for del_val in tqdm(detuning_param):
            #mat_solver = [callback(param) for param in (del_val-k_wave * velocity_param)]
            mat_solver = [callback((del_val, velocity)) for velocity in velocity_param]
            if self.is_multi_processing_enabled == True:
                with concurrent.futures.ProcessPoolExecutor() as executor:
                    results = executor.map(functools.partial(self.odeSolver, y0=y0, time_val=time_val, keys=keys),
                                           mat_solver)
                    temp_list = list(results)
            else:
                results = map(functools.partial(self.odeSolver, y0=y0, time_val=time_val, keys = keys), mat_solver)
                temp_list = list(results)

            vell = {}
            for key in keys:
                rho = [item[key] for item in temp_list]
                vell[key] = [a * b for a, b in zip(rho, velocity_dist)]

            for key in keys:
                ret_val[key].append(simps(vell[key],velocity_param))

        return ret_val

'''
def jac(t,y,param):
    delta = param[0] * t
    retVal=buildRhoMatrix(delta)
    return retVal


def userSupply(t,y,param):
    delta = param[0] * t
    rhoDot = buildRhoMatrix(delta)
    N = len(rhoDot[0])
    retVal=zeros((N,), dtype=complex)
    for i in range(N):
        sumVal=0
        for j in range(N):
            sumVal+=rhoDot[i][j]*y[j]
        retVal[i]=sumVal
    return retVal
'''