from __future__ import division,print_function
from numpy import *
from Buffer_gas_constants import *
from Relaxation_Mechanisim import *
from pylab import *

N2 = N2()
Ar = Ar()

pressures = linspace(0.1,100,1000)
T2 = []
T1 = []
s_ds = []


rubidium_num = 87
T = celsius2kelvin(90)

for pressure in pressures:
    param = [0.06, 3.48/2]
    D0 = N2.get_D0()
    N2.set_pressure(pressure)
    N2.set_density(T,pressure)
    #Ar.set_pressure(30)
    #Ar.set_density(T,30)
    y1 = wall_relaxation(T, pressure, 'cylind', param, D0)
    y2 = spin_exchange_relaxation(rubidium_num,T)
    y3 = spin_destruction_relaxation(rubidium_num,T,[N2])
    T2.append(1e3/(y1+y2+y3))
    T1.append(1e3/(y1+y3))
    s_ds.append(y3)

def wall_collision(L,r,T):
    V = pi * r **2 * L
    A = 2 * pi * r **2 + 2 * pi * r * L
    m = 87 / N_A  # [gram]
    m = 0.001 * m  # convert to Kg
    M = m / 2
    v = sqrt(8 * k * T/ (pi * M))

    return 1e3 * 4 * V / (v * A)




#print (wall_collision(0.6e-3, 1.5e-3,T))
#print (wall_collision(2e-3, 1.5e-3,T))

plot(pressures,T1,label='T1')
plot(pressures,T2,label='T2')
#plot(pressures,s_ds,label='T2')
legend()
#yscale("log")
ylabel('sec [ms]')
#ylim([10**-2,10**3])

show()