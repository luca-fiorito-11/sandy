#!/usr/local/bin/python
from PyNjoy import *
#from os import uname
jef2p2 = PyNjoy()
jef2p2.evaluationName = "Jef2.2dummy"
jef2p2.NjoyExec = "njoy2016"
jef2p2.legendre = 1
jef2p2.nstr = 22
jef2p2.iwt = 4
jef2p2.Espectra = None
jef2p2.autolib = (2.76792, 677.2873, 0.00125)

def read_float(x):
    try:
        return float(x[0] + x[1:].replace('+', 'E+').replace('-', 'E-'))
    except:
        return x


jef2p2.scatteringLaw = None
jef2p2.temperatures = ( 293., 1200. )
jef2p2.fission = None
jef2p2.dilutions = (1e10, 1e1, 1e0)

jef2p2.hmat = "H3"
jef2p2.mat = 131
jef2p2.evaluationFile = "1-H-3g.jeff33"
import pandas as pd
aaa=1
fwidths = [11,11,11,11,11,11,4,2,3]
columns = ["C1", "C2", "L1", "L2", "N1", "N2", "MAT", "MF", "MT"]
df = pd.read_fwf("1-H-3g.jeff33",
                 widths=fwidths,
                 names=columns,
                 converters={"C1" : read_float})
aaa=1
jef2p2.pendf()
jef2p2.dirName = "acefiles"
jef2p2.suffixes = ( ".03", ".12" )
jef2p2.acer(iprint=1)
#jef2p2.gendf()
#jef2p2.draglib()
