from scipy.linalg import eig, eigh
from numpy import pi, append, transpose, identity

import AtomicConstants as AC
from FundamentalConstants import *

from sz_lsi import sz, lz, Iz
from fs_hfs import Hfs,Hhfs,Bbhfs

import pylab as plt

def groundStateManifold(A_hyp_coeff,IsotopeShift,Bfield):
    """Function to produce the ground state manifold"""
    ds = int((2*S+1)*(2*I+1))  # total dimension of matrix
    #print 'Matrix dim:', ds
    As = A_hyp_coeff
    # Add the S-term hyperfine interaction
    S_StateHamiltonian = As*Hhfs(0.0,S,I)+IsotopeShift*identity(ds)
    Ez = muB*Bfield*1.e-4/(hbar*2.0*pi*1.0e6)
    S_StateHamiltonian += Ez*(gs*sz(0.0,S,I)+gI*Iz(0.0,S,I)) # Add Zeeman
    EigenSystem = eigh(S_StateHamiltonian)
    EigenValues = EigenSystem[0].real
    EigenVectors = EigenSystem[1]
    stateManifold = append([EigenValues],EigenVectors,axis=0)
    sortedManifold = sorted(transpose(stateManifold),key=(lambda i:i[0]))
    return sortedManifold, EigenValues

def excitedStateManifold(gL,A_hyp_coeff,B_hyp_coeff,Bfield):
    """Function to produce the excited state manifold"""
    dp = int(3*(2*S+1)*(2*I+1))  # total dimension of matrix
    # The actual value of FS is unimportant.
    FS = 7.123e6       # Fine structure splitting (of Rb - careful when using other elements at high B fields)
    Ap = A_hyp_coeff
    Bp = B_hyp_coeff
    # Add P-term fine and hyperfine interactions
    if Bp==0.0:
        P_StateHamiltonian=FS*Hfs(1.0,S,I)+FS*identity(dp)+Ap*Hhfs(1.0,S,I)
    if Bp!=0.0:
        P_StateHamiltonian=FS*Hfs(1.0,S,I)-(FS/2.0)*identity(dp)+Ap*Hhfs(1.0,S,I)
        P_StateHamiltonian+=Bp*Bbhfs(1.0,S,I) # add p state quadrupole
    E=muB*(Bfield*1.0e-4)/(hbar*2.0*pi*1.0e6)
    # Add magnetic interaction
    P_StateHamiltonian+=E*(gL*lz(1.0,S,I)+gs*sz(1.0,S,I)+gI*Iz(1.0,S,I))
    ep=eigh(P_StateHamiltonian)
    EigenValues=ep[0].real
    EigenVectors=ep[1]
    stateManifold=append([EigenValues],EigenVectors,axis=0)
    sortedManifold=sorted(transpose(stateManifold),key=(lambda i:i[0]))
    return sortedManifold, EigenValues


I  = AC.I        #Nuclear spin
A_hyp_coeff = AC.As #Ground state hyperfine constant in units of MHz
gI = AC.gI #nuclear spin g-factor
IsotopeShift=21.734
#Bfield=10

from pylab import *
B =linspace(1e-7,0.1e5,1000)
val=[]
val1=[]
val2=[]
val3=[]
val4=[]
val5=[]
val6=[]
val7=[]
val8=[]
val9=[]
val10=[]
val11=[]
for Bfield in B:
    x=groundStateManifold(A_hyp_coeff,IsotopeShift,Bfield)[1]
    val.append(x[0]/1e3)
    val1.append(x[1]/1e3)
    val2.append(x[2]/1e3)
    val3.append(x[3]/1e3)
    val4.append(x[4]/1e3)
    val5.append(x[5]/1e3)
    val6.append(x[6]/1e3)
    val7.append(x[7]/1e3)
    #val8.append(x[8])
    #val9.append(x[9])
    #val10.append(x[10])
    #val11.append(x[11])

'''
#Example
plt.figure(1,facecolor='w',figsize=(5, 4))

B = 0.0001 * B
plt.plot(B,val,label='$F_{1,-1}$',lw=2)
plt.plot(B,val1,label='$F_{1,0}$',lw=2)
plt.plot(B,val2,label='$F_{1,1}$',lw=2)
plt.plot(B,val3,label='$F_{2,-2}$',lw=2)
plt.plot(B,val4,label='$F_{2,-1}$',lw=2)
plt.plot(B,val5,label='$F_{2,0}$',lw=2)
plt.plot(B,val6,label='$F_{2,1}$',lw=2)
plt.plot(B,val7,label='$F_{2,2}$',lw=2)
#plot(B,val8,label='$F_{3,0}$',lw=2)
#plot(B,val9,label='$F_{3,1}$',lw=2)
#plot(B,val10,label='$F_{3,2}$',lw=2)
#plot(B,val11,label='$F_{3,3}$',lw=2)
plt.grid()
plt.legend(loc=0)
plt.xlabel('B [Tesla]')
plt.ylabel('E/h [GHz]')
plt.tight_layout()
#plt.savefig('magnetic_sublevel.png')
plt.show()
'''