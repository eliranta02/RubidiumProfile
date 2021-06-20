from __future__ import division, print_function

from numpy import zeros, transpose, kron, sqrt, dot, array
from numpy.linalg import solve
from pylab import eig,mat,inv,exp,diag
from scipy.integrate import ode


#from abc import ABCMeta, abstractmethod

class Linblad_Solver(object):


    def __init__(self):
        self._matrix = 0

    def solveSteadyState(self):
        '''

        Returns
        ------
        the solution of steady state solution
        '''
        matrix = self.matrix
        N = len(matrix[0])
        vec = zeros((N, ))
        for i in range(int(sqrt(N))):
            psi = zeros((int(sqrt(N)), ))
            psi[i] = 1
            vec += transpose(kron(psi, psi))
        matrix[0] = vec
        res = zeros((N,))
        res[0] = 1
        A = solve(matrix, res)
        return A

    def odeSolverWithParam(self,jac, userSupply, param):
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

    def odeSolver(self,y0,time_arr,returnDic):
        '''

        Parameters
        ----------
        y0 : y0 = zeros((N,)) ; y0[1] = 1
        time : runing integration parameters (e.g. time parameter)
        returnDic : all the return values

        Returns
        -------
        returnDic
        '''
        dtloop = 0
        for idx, t in enumerate(time_arr):
            #print(dtloop)
            eval,evec = eig(self.matrix)
            solT = dot((evec * (diag(exp(eval * t))) * inv(evec)), y0.transpose())
            for idx1, key in enumerate(returnDic):
                returnDic[key].append(solT[idx1])

            dtloop += 1

        return returnDic

    @property
    def matrix(self):
        return self._matrix

    @matrix.setter
    def matrix(self,temp):
        self._matrix = temp


    @staticmethod
    def print_info():
        print('This class is used in order to solve differential equations.')

    def __str__(self):
        N = len(self.matrix)
        return 'The number of levels : {} \nThe '.format(int(sqrt(N)))




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

def buildRhoMatrix(d):
    mat = array([[1,0,0,0],\
                 [0,2,3,0],\
                 [0,3,2,0],\
                 [0,0,0,d]])
    return mat


a = Linblad_Solver()

points = 100
t_scan = 1
y0 = zeros((4,))
y0[1] = 1
delta = -2

param = [points,t_scan,y0,delta]
x,y = a.odeSolverWithParam(jac,userSupply,param)

print(x)
print('*'*20)
print(y)


'''
from numpy import array, linspace

def buildRhoMatrix(d):
    mat = array([[1,0,0,0],\
                 [0,2,-3,0],\
                 [0,3,-1,0],\
                 [0,0,0,d]])
    return mat


a = Linblad_Solver()
a.matrix = buildRhoMatrix(1)
y0 = zeros((4,))
y0[1] = 1

time_arr = linspace(0, 100, 1000)


from pylab import *

myDict = {}
myDict['1,1'] = []
myDict['1,0'] = []

myDict = a.odeSolver(y0,time_arr,myDict)

plot(time_arr,real(myDict['1,1']))
plot(time_arr,real(myDict['1,0']))
show()
'''

'''
from numpy import array, dot

def buildRhoMatrix(d):
    tempmat = array([[1,0,0,0],\
                 [0,2,-3,0],\
                 [0,3,-1,0],\
                 [0,0,0,d]])
    return tempmat




a = Linblad_Solver()
a.matrix = buildRhoMatrix(2)

#print(a)

d = array([[1, 0],\
           [0,1]])

cc = d.reshape((4,))
c = transpose(cc)


print(dot(c,cc))

'''