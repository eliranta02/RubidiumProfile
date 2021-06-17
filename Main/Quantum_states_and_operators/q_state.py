from __future__ import division,print_function
from numpy import zeros,transpose,tensordot
from Main.Constants.Rb_constants import S
class State:


    def __init__(self,levels_number,key,F=0,mF=0,L=0,gF=0,is_Ground = False):
        '''
        Parameters
        ----------
        levels_number : how many levels total in your system in order to define the vector size
        key - unique value of the state (default the state number)
        F - F number of the state
        mF - qunatum number of F
        L - angular momentum number
        gF -
        is_Ground - boolean if it is the ground state or not
        '''
        self.level = zeros((levels_number, 1 ))
        self.L = L
        self.F = F
        self.mF = mF
        self.gF = gF
        self.is_Ground = is_Ground
        if 0 <= key < levels_number:
            self.level[key] = 1

    def set_S(self,S):
        self.S = S

    def set_L(self,L):
        self.L = L

    def set_gF(self,gF):
        self.gF = gF


    def set_F(self,F):
        self.F = F

    def set_mF(self, mF):
        self.mF = mF

    def set_isGround(self,is_Ground):
        self.is_Ground = is_Ground

    def set_key(self,key):
        self.key = key


    def get_L(self):
        return self.L

    def get_J(self):
        J = self.L + self.S
        return J

    def get_gF(self):
        return self.gF

    def get_F(self):
        return (self.F)

    def get_mF(self):
        return (self.mF)

    def get_isGround(self):
        return (self.is_Ground)

    def get_key(self):
        return (self.key)

    def get_level(self):
        return self.level

    def outer_product(self,state):
        lhs = state.get_level()
        ret_val = tensordot(self.get_level(),lhs,axes=0)
        ret_val = ret_val.reshape(len(self.get_level())*len(lhs),1)
        return ret_val

    def __str__(self):
        if self.is_Ground:
            ret_val = '|Fg:{},mFg:{}>'.format(self.F,self.mF)
        else:
            ret_val = '|Fe:{},mFe:{}>'.format(self.F, self.mF)
        return ret_val

    def __mul__(self, other):
        ret_val = None
        if isinstance(other, State):
            ret_val = self.level * transpose(other.level)
        else:
            print ('Multiply must be between to states')
        return ret_val

    def __rmul__(self, other):
        ret_val = None
        if isinstance(other, State):
            ret_val = other.level * transpose(self.level)
        else:
            print('Multiply must be between to states')
        return ret_val
