NEURON {
    SUFFIX cap
    USEION oxc READ oxci 
    USEION ox READ oxo
    RANGE qv, qflow, ioxc, oxo
}

ASSIGNED {
    oxo
    oxci
    ioxc
}

PARAMETER {
    qv = 2
    qflow = 2
    alphab = .000031
    hd = .2
    po50=38
}

BREAKPOINT {
    ioxc = -qv*(oxci - oxo) / (qflow * (alphab + hd * 3 * pow(po50,3) * pow(oxci,2) / pow((pow(po50,3) + pow(oxci,3)),2) ))
}









