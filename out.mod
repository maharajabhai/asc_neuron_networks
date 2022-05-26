NEURON {
    SUFFIX out
    USEION ox  READ oxo, oxi WRITE iox VALENCE 0
    USEION oxc READ oxci
    RANGE  kout,  iox
}

ASSIGNED {
    iox 
    oxo
    oxi
    oxci
}

PARAMETER {
    kout = 0
}

BREAKPOINT {
    iox = kout*(oxci-oxo)
}








