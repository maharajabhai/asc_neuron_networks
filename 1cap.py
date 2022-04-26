from neuron import h
from neuron import crxd as rxd
from matplotlib import pyplot as plt
import numpy as np
from neuron.crxd import rxdmath
h.load_file('stdrun.hoc')
rxd.options.enable.extracellular = True

tottime = 1000

sc=1000
oxdiff=.054

nsg = 21
eps  = 1/nsg

nvess = 1

class Vessel:
    def __init__(self, x1, y1, z1, x2, y2, z2, q):
        self.x1 = x1
        self.y1 = y1
        self.z1 = z1
        self.x2 = x2
        self.y2 = y2
        self.z2 = z2
        self.vess = h.Section(name='vess', cell=self)
        h.pt3dadd( x1, y1, z1, 1, sec=self.vess)
        h.pt3dadd( x2, y2, z2, 1, sec=self.vess)
        for mech in ("deriv","out","cap"):
            self.vess.insert(mech)
        for k in range(nsg):
            self.vess(eps*(.5+k)).cap.qv = q*sc
            self.vess(eps*(.5+k)).deriv.delta = 1*sc
#            self.vess(eps*(.5+k)).out.kout = .001*sc*rho
            self.vess(eps*(.5+k)).out.kout = .1*sc
#            self.vess(eps*(.5+k)).out.nkout = 0

vess_params = [
    [0, 0, 0, 0, 0,  nsg, .002],
]

vessels = []
caps = []

for param in vess_params:
    vessels.append(Vessel(*param))

for (n, v) in enumerate(vessels):
    caps.append(rxd.Region([v.vess], name='cap' + str(n), nrn_region='i'))

for j in range(nvess):
    vessels[j].vess.nseg = nsg

ecs = rxd.Extracellular(-20, -20, 0, 20, 20, nsg, dx=3)

secs3d = [ecs]
secs1d = [sec for sec in h.allsec() if sec not in secs3d]

ox = rxd.Species([ecs], name='ox', d= oxdiff, charge=1, initial=20)
oxc = rxd.Species(caps, name='oxc', d=0, charge=0, initial = 70)

ak=3
m=-5

r = rxd.Rate(ox,m * ox/(ak + ox), membrane_flux=False)

# oxc is constant at nodes where blood enter network 

vessels[0].vess(eps/2).cap.qv = 0
vessels[0].vess(eps/2).deriv.delta = 0

for cap in caps:
    for j in range(nsg):
        nd = oxc[cap].nodes[j]
        seg = nd.segment
        nd.include_flux(seg.deriv._ref_oxflux)
        nd.include_flux(seg.cap._ref_ioxc)

#  oxm and oxmb must point to something at every node                                                  
for i in range(1):
    for j in range(nsg):
        vessels[i].vess(eps*(.5+j)).deriv._ref_oxm = oxc[caps[i]].nodes[j]._ref_value
        vessels[i].vess(eps*(.5+j)).deriv._ref_oxmb = oxc[caps[i]].nodes[j]._ref_value

# derivative along each vessel                                                                         
for i in range(1):
    for j in range(1,nsg):
        vessels[i].vess(eps*(.5+j)).deriv._ref_oxm = oxc[caps[i]].nodes[j-1]._ref_value 


#states_init = ox[ecs].states3d.copy()

t = h.Vector().record(h._ref_t)
y=list(range(nsg))
for j in range(nsg):
    y[j] = h.Vector().record(vessels[0].vess(eps*(.1+j))._ref_oxci)

#h.finitialize(-70)

h.tstop = tottime
h.run()

ss = h.SaveState()
ss.save()

sf = h.File('state.bin')
ss.fwrite(sf)

plt.subplot(2, 2, 1)
for j in range(nsg):
    plt.plot(t, y[j])

y1=list(range(nsg))
for i in range(nsg):
    y1[i] = oxc[caps[0]].nodes[i].concentration
a1=list(range(0,nsg))


plt.subplot(2, 2, 2)
plt.plot(a1, y1)

states_mid = ox[ecs].states3d.copy()

z=states_mid[:, int(ecs._ny/2), :]

plt.subplot(2,1,2)
plt.imshow(states_mid[:, int(ecs._ny/2), :])

plt.colorbar()

plt.show()




