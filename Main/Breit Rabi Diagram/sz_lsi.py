"""
Functions to calculate the z-projection matrices.

Calculates the full size matrices for the projection on the quantization
axis for electron spin, orbital angular momentum, nuclear spin, and the
coupled angular momentum F. Essentially takes results from jz and puts
them in the full Hilbert space.

Calls jz from ang_mom.
"""

from numpy import identity
from scipy.linalg import kron
from ang_mon import jz

def sz(L,S,I):
    Sz=jz(S)
    gL=int(2*L+1)
    Li=identity(gL)
    gI=int(2*I+1)
    Ii=identity(gI)
    sz=kron(kron(Li,Sz),Ii)
    return sz

def lz(L,S,I):
    gS=int(2*S+1)
    Si=identity(gS)
    Lz=jz(L)
    gI=int(2*I+1)
    Ii=identity(gI)
    lz=kron(kron(Lz,Si),Ii)
    return lz

def Iz(L,S,I):
    gS=int(2*S+1)
    gL=int(2*L+1)
    Si=identity(gS)
    Li=identity(gL)
    Iz_num=jz(I)
    Iz=kron(kron(Li,Si),Iz_num)
    return Iz

def fz(L,S,I):
    gS=int(2*S+1)
    Sz=jz(S)
    Si=identity(gS)
    gL=int(2*L+1)
    Lz=jz(L)
    Li=identity(gL)
    gJ=gL*gS
    Jz=kron(Lz,Si)+kron(Li,Sz)
    Ji=identity(gJ)
    gI=int(2*I+1)
    Iz=jz(I)
    Ii=identity(gI)
    Fz=kron(Jz,Ii)+kron(Ji,Iz)
    return Fz
