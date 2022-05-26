from neuron import h, gui
from matplotlib import pyplot

soma = h.Section()

soma.insert('hh')

h.pt3dadd( 0, 0, 0, 2, sec=soma)
h.pt3dadd( 0, 0, 2, 2, sec=soma)

#Insert Clamp                                                                                      
iclamp = h.IClamp(0.5, sec=soma)
iclamp.dur = 100
iclamp.delay = 20
iclamp.amp = .01

v_soma = h.Vector()
t = h.Vector()
v_soma.record(soma(0.5)._ref_v)
t.record(h._ref_t)

h.finitialize(-70)
h.continuerun(500)

pyplot.plot(t, v_soma)
pyplot.show()
