#!/usr/local/bin/python
from PyNjoy import *
from os import uname
jef2p2 = PyNjoy()
jef2p2.evaluationName = "Jef2.2"
jef2p2.execDir = "../" + uname()[0]
jef2p2.Espectra = None

# Process the burnup chain:

jef2p2.fissionFile = "$HOME/evaluations/Jef2.2/tape24"
jef2p2.decayFile = "$HOME/evaluations/Jef2.2/tape22"
jef2p2.burnup()
