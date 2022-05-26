from neuron import h
from neuron import crxd as rxd
from matplotlib import pyplot as plt
import numpy as np
from neuron.crxd import rxdmath
h.load_file('stdrun.hoc')
rxd.options.enable.extracellular = True

tottime=1000

diam=5

sc=1000

oxdiff = .054

with open("pos.csv") as pos:
   pos = np.loadtxt(pos, delimiter=",")
nvess = np.shape(pos)[0]

with open("conn.csv") as conn:
   conn = np.loadtxt(conn, delimiter=",")

with open("dim.csv") as dim:
   dim = np.loadtxt(dim, delimiter=",")

class Vessel:
    def __init__(self, x1, y1, z1, x2, y2, z2, q):
        self.x1 = x1
        self.y1 = y1
        self.z1 = z1
        self.x2 = x2
        self.y2 = y2
        self.z2 = z2
        self.vess = h.Section(name='vess', cell=self)
        h.pt3dadd( x1, y1, z1, diam, sec=self.vess)
        h.pt3dadd( x2, y2, z2, diam, sec=self.vess)
        for mech in ("deriv","out","cap"):
            self.vess.insert(mech)

vess_params = pos
vessels = []
caps = []
for param in vess_params:
    vessels.append(Vessel(*param))
for (n, v) in enumerate(vessels):
    caps.append(rxd.Region([v.vess], name='cap' + str(n), nrn_region='i'))

# The nseg of each vessel approx = vessel length
nseg = np.zeros(nvess, dtype = int)
eps  = list(range(nvess))

for j in range(nvess):
    p1 = np.array((pos[j][0], pos[j][1], pos[j][2]))
    p2 = np.array((pos[j][3], pos[j][4], pos[j][5]))
    dist = np.linalg.norm(p1-p2)
    nseg[j] = np.ceil(dist)
    eps[j] = 1/nseg[j]

for j in range(nvess):
    for i in range(nseg[j]):
            vessels[j].vess(eps[j]*(.5+i)).cap.qv = .001*diam*sc
            vessels[j].vess(eps[j]*(.5+i)).cap.qflow = pos[j][6]
            vessels[j].vess(eps[j]*(.5+i)).deriv.delta = 1*sc
            vessels[j].vess(eps[j]*(.5+i)).out.kout = .001*diam*sc
            vessels[j].vess(eps[j]*(.5+i)).out.nkout = 0

for j in range(nvess):
    vessels[j].vess.nseg = nseg[j]
    vessels[j].vess.diam = diam

ecs = rxd.Extracellular(dim[0]-10, dim[1], dim[2], dim[3]+10, dim[4], dim[5], dx=3, 
                        volume_fraction=0.2,tortuosity=1.6)
secs3d = [ecs]
secs1d = [sec for sec in h.allsec() if sec not in secs3d]

ox = rxd.Species([ecs], name='ox', d= oxdiff, charge=1, initial=
          lambda nd: 30 if nd.region == ecs else 20)
oxc = rxd.Species(caps, name='oxc', d=0, charge=0, initial = 70)


for i in range(nvess):
    for j in range(nseg[i]):
        nd = oxc[caps[i]].nodes[j]
        seg = nd.segment
        nd.include_flux(seg.deriv._ref_oxflux)
        nd.include_flux(seg.cap._ref_ioxc)

# oxc is constant at nodes where blood enter network                                      
for j in range(nvess):
       if conn[j][0] == 0:
          vessels[j].vess(eps[j]/2).cap.qv = 0
          vessels[j].vess(eps[j]/2).deriv.delta = 0

#  oxm and oxmb must point to something at every node                                                  
for i in range(nvess):
    for j in range(nseg[i]):
        vessels[i].vess(eps[i]*(.5+j)).deriv._ref_oxm = oxc[caps[i]].nodes[j]._ref_value
        vessels[i].vess(eps[i]*(.5+j)).deriv._ref_oxmb = oxc[caps[i]].nodes[j]._ref_value

# derivative along each vessel                                                                         
for i in range(nvess):
    for j in range(1,nseg[i]):
        vessels[i].vess(eps[i]*(.5+j)).deriv._ref_oxm = oxc[caps[i]].nodes[j-1]._ref_value 
        vessels[i].vess(eps[i]/2).deriv.delta = 1*sc

# connectivity                                                                            
for j in range(nvess):
      if conn[j][0] != 0:
         i=int(conn[j][0])
         vessels[j].vess(eps[j]/2).deriv._ref_oxm = oxc[caps[i-1]].nodes[nseg[i-1]-1]._ref_value
         vessels[j].vess(eps[j]/2).deriv.delta = 1*sc
         if conn[j][1] != 0:
            i=int(conn[j][1])
            n=int(conn[j][0])
            vessels[j].vess(eps[j]/2).deriv._ref_oxm = oxc[caps[n]].nodes[nseg[n]-1]._ref_value
            vessels[j].vess(eps[j]/2).deriv._ref_oxmb = oxc[caps[i]].nodes[nseg[i]-1]._ref_value
            vessels[j].vess(eps[j]/2).deriv.delta = pos[n][6]*sc
            vessels[j].vess(eps[j]/2).deriv.deltab = pos[i][6]*sc

states_init = ox[ecs].states3d.copy()

t = h.Vector().record(h._ref_t)

h.finitialize(-70)

h.continuerun(tottime)

ad=int(dim[3])+10
zox=list(range(ad))
for j in range(ad):
    zox[j] = ox[ecs].node_by_location(j-2,0, 5).value

plt.subplot(2, 1, 1)
plt.plot(zox)

states_mid = ox[ecs].states3d.copy()

z=states_mid[:, int(ecs._ny/2), :]

plt.subplot(2,1,2)
plt.imshow(states_mid[:, int(ecs._ny/2), :])

plt.colorbar()

import csv

#with open('data.csv', mode='w') as f:
#    data = csv.writer(f)
#    data.writerows(z)

plt.show()




