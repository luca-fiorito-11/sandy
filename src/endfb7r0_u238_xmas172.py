#!/usr/local/bin/python
from PyNjoy import *
from os import uname
endfb = PyNjoy()
endfb.evaluationName = "/tmp/xmas172_u238_endfb7r0"
endfb.execDir = "../" + uname()[0]
endfb.nstr = 22
endfb.iwt = 4
endfb.Espectra = None
endfb.autolib = (2.76792, 677.2873, 0.00125)

endfb.scatteringLaw = None
endfb.temperatures = ( 293., )
endfb.fission = None
endfb.dilutions = None

endfb.hmat = "U238"
endfb.mat = 9237
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r0/n-ENDF-VII0.endf/n-092_U_238.endf"
endfb.fission = 2 # fission with delayed neutrons
endfb.ss = (2.76792, 1.22773e5)
endfb.potential = 11.17103
endfb.dilutions = ( 1.e10, 94.5317612 )
endfb.pendf()
endfb.gendf()
endfb.draglib()
