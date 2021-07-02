from __future__ import division
from numpy import transpose,dot
import ang_mon_p

def jx(jj):
    jp=ang_mon_p.jp(jj)
    jm=transpose(jp)
    jx=0.5*(jp+jm)
    return jx

def jy(jj):
    jp=ang_mon_p.jp(jj)
    jm=transpose(jp)
    jy=0.5j*(jm-jp)
    return jy

def jz(jj):
    jp=ang_mon_p.jp(jj)
    jm=transpose(jp)
    jz=0.5*(dot(jp,jm)-dot(jm,jp))
    return jz
