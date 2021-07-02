from numpy import identity,dot
from scipy.linalg import kron
from ang_mon import jx,jy,jz


def Hfs(L,S,I):
    """Provides the L dot S matrix (fine structure)"""
    gS=int(2*S+1) #number of mS values
    Sx=jx(S)
    Sy=jy(S)
    Sz=jz(S)
    Si=identity(gS)

    gL=int(2*L+1)
    Lx=jx(L)
    Ly=jy(L)
    Lz=jz(L)
    Li=identity(gL)

    gJ=gL*gS
    Jx=kron(Lx,Si)+kron(Li,Sx)
    Jy=kron(Ly,Si)+kron(Li,Sy)
    Jz=kron(Lz,Si)+kron(Li,Sz)
    J2=dot(Jx,Jx)+dot(Jy,Jy)+dot(Jz,Jz)

    gI=int(2*I+1)
    Ii=identity(gI)
    gF=gJ*gI
    Fi=identity(gF)
    Hfs=0.5*(kron(J2,Ii)-L*(L+1)*Fi-S*(S+1)*Fi) # fine structure in m_L,m_S,m_I basis
    return Hfs

def Hhfs(L,S,I):
    """Provides the I dot J matrix (hyperfine structure interaction)"""
    gS=int(2*S+1)
    Sx=jx(S)
    Sy=jy(S)
    Sz=jz(S)
    Si=identity(gS)

    gL=int(2*L+1)
    Lx=jx(L)
    Ly=jy(L)
    Lz=jz(L)
    Li=identity(gL)

    gJ=gL*gS
    Jx=kron(Lx,Si)+kron(Li,Sx)
    Jy=kron(Ly,Si)+kron(Li,Sy)
    Jz=kron(Lz,Si)+kron(Li,Sz)
    Ji=identity(gJ)
    J2=dot(Jx,Jx)+dot(Jy,Jy)+dot(Jz,Jz)

    gI=int(2*I+1)
    gF=gJ*gI
    Ix=jx(I)
    Iy=jy(I)
    Iz=jz(I)
    Ii=identity(gI)
    Fx=kron(Jx,Ii)+kron(Ji,Ix)
    Fy=kron(Jy,Ii)+kron(Ji,Iy)
    Fz=kron(Jz,Ii)+kron(Ji,Iz)
    Fi=identity(gF)
    F2=dot(Fx,Fx)+dot(Fy,Fy)+dot(Fz,Fz)
    Hhfs=0.5*(F2-I*(I+1)*Fi-kron(J2,Ii))
    return Hhfs

def Bbhfs(L,S,I):
    """Calculates electric quadrupole matrix.

    Calculates the part in square brakets from
    equation (8) in manual
    """
    gS=int(2*S+1)
    Sx=jx(S)
    Sy=jy(S)
    Sz=jz(S)
    Si=identity(gS)

    gL=int(2*L+1)
    Lx=jx(L)
    Ly=jy(L)
    Lz=jz(L)
    Li=identity(gL)

    gJ=gL*gS
    Jx=kron(Lx,Si)+kron(Li,Sx)
    Jy=kron(Ly,Si)+kron(Li,Sy)
    Jz=kron(Lz,Si)+kron(Li,Sz)

    gI=int(2*I+1)
    gF=gJ*gI
    Ix=jx(I)
    Iy=jy(I)
    Iz=jz(I)

    Fi=identity(gF)

    IdotJ=kron(Jx,Ix)+kron(Jy,Iy)+kron(Jz,Iz)
    IdotJ2=dot(IdotJ,IdotJ)

    if I != 0:
        Bbhfs=1./(6*I*(2*I-1))*(3*IdotJ2+3./2*IdotJ-I*(I+1)*15./4*Fi)
    else:
        Bbhfs = 0
    return Bbhfs
