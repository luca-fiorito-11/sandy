#!/usr/local/bin/python
from PyNjoy import *
from os import uname
jef2p2 = PyNjoy()
jef2p2.evaluationName = "Jef2.2dummy"
jef2p2.execDir = "../" + uname()[0]
jef2p2.legendre = 1
jef2p2.nstr = 22
jef2p2.iwt = 4
jef2p2.Espectra = None
jef2p2.autolib = (2.76792, 677.2873, 0.00125)

jef2p2.hmat = "H1_H2O"
jef2p2.mat = 125
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape1"
jef2p2.scatteringLaw = "$HOME/evaluations/Jef2.2/tape21"
jef2p2.scatteringMat = 1
jef2p2.temperatures = ( 293.6, 623.6 )
jef2p2.fission = None
jef2p2.dilutions = None
jef2p2.pendf()
jef2p2.gendf()
jef2p2.draglib()

jef2p2.scatteringLaw = None
jef2p2.temperatures = ( 293., 1200. )
jef2p2.fission = None
jef2p2.dilutions = None

jef2p2.hmat = "H3"
jef2p2.mat = 131
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape1"
jef2p2.pendf()
jef2p2.gendf()
jef2p2.draglib()

jef2p2.hmat = "He3"
jef2p2.mat = 225
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape1"
jef2p2.pendf()
jef2p2.gendf()
jef2p2.draglib()

# Process the burnup chain:

jef2p2.fissionFile = "$HOME/evaluations/Jef2.2/tape24"
jef2p2.decayFile = "$HOME/evaluations/Jef2.2/tape22"
jef2p2.burnup()
