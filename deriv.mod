NEURON {
    SUFFIX deriv
    USEION oxc READ oxci  VALENCE 0
    POINTER oxm, oxmb
    RANGE oxflux, delta, oxfluxb, deltab
}

ASSIGNED {
    oxci 
    oxm 
    oxflux
    oxmb 
    oxfluxb
}

PARAMETER {
    delta = 1000
    deltab = 0
}

BREAKPOINT {
   oxflux = -delta * (oxci - oxm)
   oxfluxb = -deltab * (oxci - oxmb)

}