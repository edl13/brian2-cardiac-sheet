from brian2 import *
import numpy as np

# Constants (Endocardium)
Dl = 0.001 * cm ** 2 / ms
Dt = 0.001 * cm ** 2 / ms
Cm = 1 * uF / (cm**2)

tauhp = 10.8 * ms
tauhm = 10.8 * ms
taufp = 355 * ms
taufm = 52.3 * ms
taurp = 7.54 * ms
taurm = 6.07 * ms
tausp = 29.1 * ms
tausm = 10.4 * ms
gfi = 1.72 / ms
vfi = 1.24
gso = 0.00891 / ms
gsi = 0.414 / ms
beta1 = 22.9
beta2 = 2.95
v1 = 0.522
v2 = 0.596
gto = 0.300 / ms
vc = 0.130
vs = 0.6


vr = 0.6
vto = 0
vrest = -85 * mv

# Membrane Equations
eqs = '''
dv/dt = Idiff - (Iion - Istim)/Cm : volt
vu = (v - vrest) / (100 * mV) : 1
Iion = Ifi + Isi + Iso + Ito : amp/cm**2

# Scaled Currents
Jfi = -gfi * h * minf * (vfi - vu) : Hz
Jsi = -gsi * dinf * f * finf2 : Hz
Jto = gto * r * s * (vu - vto) : Hz
Jso = gso * kinf : Hz

# Gate variables
dh/dt = (hinf - h)/tauh : 1
df/dt = (finf - f)/tauf : 1
dr/dt = (rinf - r)/taur : 1
ds/dt = (sinf - s)/taus : 1

# Heaviside
theta = int((vu-vc) >= 0): 1

# Constants Equations
minf = (vu-vc) * theta >= 0): 1
hinf = 1 - theta : 1
finf = 1 - theta : 1
rinf = 1 - theta : 1
sinf = 1 - theta : 1
dinf = theta * (1 + tanh(beta1 * (vu-v1))/2 : 1
finf2 = (1 - tanh(beta2 * (vu-v2))/2 : 1
kinf = vu/vc + theta * (1 - vu/vc) : 1
rinf = int((vu-vr) >= 0 : 1
tauh = tauhp - (tauhp - tauhm) * theta : Hz
tauf = taufp - (taufp - taufm) * theta : Hz
taur = taurp - (taurp - taurm) * theta : Hz
taus = tausp - (tausp - tausm) * theta : Hz
'''
