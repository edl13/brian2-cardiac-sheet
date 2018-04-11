from brian2 import *
defaultclock.dt = 1 * ms


def getCoords(ind, xmax, ymax):
    x = ind % xmax
    y = floor((ind / xmax)) % ymax
    return (x, y)


def getInd(coords, xmax, ymax):
    x = coords[0]
    y = coords[1]
    ind = xmax * y + x
    return (ind)


@network_operation(dt=defaultclock.dt)
def longDer(t):
    print(t)
    for i in range(len(neurons)):
        left = 0 * mV
        right = 0 * mV
        center = neurons.v[i]
        if (i - 1) % int(len(neurons)**0.5) >= 0:
            left = neurons.v[i - 1]

        if (i + 1) % int(len(neurons)**0.5) != 0:
            right = neurons.v[i + 1]

        neurons[i:i + 1].xDer = (left + right - 2 * center) / dx**2


@network_operation(dt=defaultclock.dt)
def transDer():
    for i in range(len(neurons)):
        left = 0 * mV
        right = 0 * mV
        center = neurons.v[i]
        if (i - numRow) >= 0:
            left = neurons.v[i - 1]
        if (i + numRow) < len(neurons):
            right = neurons.v[i + 1]

        neurons[i:i + 1].yDer = (left + right - 2 * center) / dx**2


# Constants (Endocardium)
Dl = 0.001 * cm**2 / ms
Dt = 0.001 * cm**2 / ms
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
vrest = -85 * mV

# Membrane Equations
eqs = '''
dv/dt = Idiff - (Iion - Istim)/Cm : volt
xDer: volt / meter**2
yDer: volt / meter**2
Idiff = Dl * xDer + Dt * yDer: volt/second
vu = (v - vrest) / (100 * mV) : 1
scale = Cm * 100 * mV: coulomb / meter ** 2
Iion = Jfi * scale + Jsi * scale + Jto * scale + Jso * scale: amp / meter ** 2
Istim: amp / meter ** 2

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
minf = (vu-vc) * theta : 1
hinf = 1 - theta : 1
finf = 1 - theta : 1
sinf = 1 - theta : 1
dinf = theta * (1 + tanh(beta1 * (vu-v1)))/2 : 1
finf2 = (1 - tanh(beta2 * (vu-v2)))/2 : 1
kinf = vu/vc + theta * (1 - vu/vc) : 1
rinf = int((vu-vr) >= 0): 1
tauh = tauhp - (tauhp - tauhm) * theta : second
tauf = taufp - (taufp - taufm) * theta : second
taur = taurp - (taurp - taurm) * theta : second
taus = tausp - (tausp - tausm) * theta : second
'''

sheetLength = 5 * cm
dx = 0.1 * cm
numRow = sheetLength / dx
numNeurons = (numRow)**2

neurons = NeuronGroup(numNeurons, eqs)
M = StateMonitor(
    neurons, ['vu', 'Idiff', 'Istim', 'xDer', 'yDer'], record=True, dt=2 * ms)

neurons.Istim = 0
run(10 * ms, report='text')
center = int(getInd((25, 25), numRow, numRow))
neurons[center - 1:center].Istim = 5 * mA / cm**2
run(1 * ms, report='text')
neurons.Istim = 0 * mA / cm**2
run(1000 * ms, report='text')
