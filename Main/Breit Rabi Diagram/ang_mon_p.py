from __future__ import division
from numpy import zeros,sqrt,arange

def jp(jj):
    b = 0
    dim = int(2*jj+1)
    jp = zeros((dim,dim))
    z = arange(dim)
    m = jj-z
    while b<dim-1.0:
        mm = m[b+1]
        jp[b,b+1] = sqrt(jj*(jj+1)-mm*(mm+1))
        b = b+1
    return jp

