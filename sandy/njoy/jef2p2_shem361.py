#!/usr/local/bin/python
from PyNjoy import *
from os import uname
jef2p2 = PyNjoy()
jef2p2.evaluationName = "/tmp/shem361_Jef2.2"
jef2p2.execDir = "../" + uname()[0]
jef2p2.nstr = 26
jef2p2.iwt = 4
jef2p2.Espectra = None
jef2p2.autolib = (22.53556, 1.11377e4, 0.0005)

jef2p2.legendre = 3
jef2p2.hmat = "H1_H2O"
jef2p2.mat = 125
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape1"
jef2p2.scatteringLaw = "$HOME/evaluations/Jef2.2/tape21"
jef2p2.scatteringMat = 1
jef2p2.temperatures = ( 293.6, 373.6, 523.6, 623.6 )
jef2p2.fission = None
jef2p2.dilutions = None
jef2p2.pendf()
jef2p2.gendf()
jef2p2.draglib()

jef2p2.legendre = 1
jef2p2.hmat = "H2_D2O"
jef2p2.mat = 128
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape1"
jef2p2.scatteringLaw = "$HOME/evaluations/Jef2.2/tape21"
jef2p2.scatteringMat = 11
jef2p2.temperatures = ( 293.6, 373.6, 523.6, 673.6 )
jef2p2.fission = None
jef2p2.dilutions = None
jef2p2.pendf()
jef2p2.gendf()
jef2p2.draglib()

jef2p2.hmat = "H1_CH2"
jef2p2.mat =  125
jef2p2.temperatures = ( 293.6,  350. )
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape1"
jef2p2.scatteringLaw = "$HOME/evaluations/Jef2.2/tape21"
jef2p2.scatteringMat = 37
jef2p2.fission = None
jef2p2.dilutions = None
jef2p2.pendf()
jef2p2.gendf()
jef2p2.draglib()

jef2p2.hmat = "Be9"
jef2p2.mat = 425
jef2p2.temperatures = ( 293.6 , 400 , 500 , 600 , 700 , 800 , 1000 , 1200)
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape1"
jef2p2.scatteringLaw = "$HOME/evaluations/Jef2.2/tape21"
jef2p2.scatteringMat = 26
jef2p2.fission = None
jef2p2.dilutions = None
jef2p2.pendf()
jef2p2.gendf()
jef2p2.draglib()

jef2p2.hmat = "C0_GR"
jef2p2.mat = 600
jef2p2.temperatures = ( 293.6 , 400 , 500 , 600 , 700 , 800 , 1000 , 1200 , 1600 , 2000 )
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape1"
jef2p2.scatteringLaw = "$HOME/evaluations/Jef2.2/tape21"
jef2p2.scatteringMat = 31
jef2p2.fission = None
jef2p2.dilutions = None
jef2p2.pendf()
jef2p2.gendf()
jef2p2.draglib()

jef2p2.scatteringLaw = None
jef2p2.temperatures = ( 293., 550., 900., 1200., 2000. )
jef2p2.fission = None
jef2p2.dilutions = None

jef2p2.hmat = "H1"
jef2p2.mat = 125
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape1"
jef2p2.pendf()
jef2p2.gendf()
jef2p2.draglib()

jef2p2.hmat = "H2"
jef2p2.mat = 128
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape1"
jef2p2.pendf()
jef2p2.gendf()
jef2p2.draglib()

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

jef2p2.hmat = "He4"
jef2p2.mat = 228
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape1"
jef2p2.pendf()
jef2p2.gendf()
jef2p2.draglib()

jef2p2.hmat = "Li6"
jef2p2.mat = 325
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape1"
jef2p2.pendf()
jef2p2.gendf()
jef2p2.draglib()

jef2p2.hmat = "Li7"
jef2p2.mat = 328
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape1"
jef2p2.pendf()
jef2p2.gendf()
jef2p2.draglib()

jef2p2.hmat = "B10"
jef2p2.mat = 525
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape1"
jef2p2.pendf()
jef2p2.gendf()
jef2p2.draglib()

jef2p2.hmat = "B11"
jef2p2.mat = 528
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape1"
jef2p2.pendf()
jef2p2.gendf()
jef2p2.draglib()

jef2p2.hmat = "C0"
jef2p2.mat = 600
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape1"
jef2p2.pendf()
jef2p2.gendf()
jef2p2.draglib()

jef2p2.hmat = "N14"
jef2p2.mat = 725
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape1"
jef2p2.pendf()
jef2p2.gendf()
jef2p2.draglib()

jef2p2.hmat = "N15"
jef2p2.mat = 728
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape1"
jef2p2.pendf()
jef2p2.gendf()
jef2p2.draglib()

jef2p2.hmat = "O16"
jef2p2.mat = 825
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape1"
jef2p2.pendf()
jef2p2.gendf()
jef2p2.draglib()

jef2p2.hmat = "O17"
jef2p2.mat = 828
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape1"
jef2p2.pendf()
jef2p2.gendf()
jef2p2.draglib()

jef2p2.hmat = "F19"
jef2p2.mat = 925
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape1"
jef2p2.pendf()
jef2p2.gendf()
jef2p2.draglib()

jef2p2.hmat = "Na23"
jef2p2.mat = 1125
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape1"
jef2p2.pendf()
jef2p2.gendf()
jef2p2.draglib()

jef2p2.hmat = "Mg0"
jef2p2.mat = 1200
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape1"
jef2p2.pendf()
jef2p2.gendf()
jef2p2.draglib()

jef2p2.hmat = "Al27"
jef2p2.mat = 1325
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape1"
jef2p2.pendf()
jef2p2.gendf()
jef2p2.draglib()

jef2p2.hmat = "Si0"
jef2p2.mat = 1400
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape1"
jef2p2.pendf()
jef2p2.gendf()
jef2p2.draglib()

jef2p2.hmat = "P31"
jef2p2.mat = 1525
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape1"
jef2p2.pendf()
jef2p2.gendf()
jef2p2.draglib()

jef2p2.hmat = "S32"
jef2p2.mat = 1625
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape1"
jef2p2.pendf()
jef2p2.gendf()
jef2p2.draglib()

jef2p2.hmat = "S33"
jef2p2.mat = 1628
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape1"
jef2p2.pendf()
jef2p2.gendf()
jef2p2.draglib()

jef2p2.hmat = "S34"
jef2p2.mat = 1631
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape1"
jef2p2.pendf()
jef2p2.gendf()
jef2p2.draglib()

jef2p2.hmat = "S36"
jef2p2.mat = 1637
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape1"
jef2p2.pendf()
jef2p2.gendf()
jef2p2.draglib()

jef2p2.hmat = "Cl0"
jef2p2.mat = 1700
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape1"
jef2p2.pendf()
jef2p2.gendf()
jef2p2.draglib()

jef2p2.hmat = "K0"
jef2p2.mat = 1900
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape1"
jef2p2.pendf()
jef2p2.gendf()
jef2p2.draglib()

jef2p2.hmat = "Ca0"
jef2p2.mat = 2000
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape1"
jef2p2.pendf()
jef2p2.gendf()
jef2p2.draglib()

jef2p2.hmat = "Ti0"
jef2p2.mat = 2200
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape1"
jef2p2.pendf()
jef2p2.gendf()
jef2p2.draglib()

jef2p2.hmat = "V0"
jef2p2.mat = 2300
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape1"
jef2p2.pendf()
jef2p2.gendf()
jef2p2.draglib()

jef2p2.hmat = "Cr50"
jef2p2.mat = 2425
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape2"
jef2p2.pendf()
jef2p2.gendf()
jef2p2.draglib()

jef2p2.hmat = "Cr52"
jef2p2.mat = 2431
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape2"
jef2p2.pendf()
jef2p2.gendf()
jef2p2.draglib()

jef2p2.hmat = "Cr53"
jef2p2.mat = 2434
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape2"
jef2p2.pendf()
jef2p2.gendf()
jef2p2.draglib()

jef2p2.hmat = "Cr54"
jef2p2.mat = 2437
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape2"
jef2p2.pendf()
jef2p2.gendf()
jef2p2.draglib()

jef2p2.hmat = "Mn55"
jef2p2.mat = 2525
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape2"
jef2p2.pendf()
jef2p2.gendf()
jef2p2.draglib()

jef2p2.hmat = "Fe54"
jef2p2.mat = 2625
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape2"
jef2p2.pendf()
jef2p2.gendf()
jef2p2.draglib()

jef2p2.hmat = "Fe56"
jef2p2.mat = 2631
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape2"
jef2p2.pendf()
jef2p2.gendf()
jef2p2.draglib()

jef2p2.hmat = "Fe57"
jef2p2.mat = 2634
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape2"
jef2p2.pendf()
jef2p2.gendf()
jef2p2.draglib()

jef2p2.hmat = "Fe58"
jef2p2.mat = 2637
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape2"
jef2p2.pendf()
jef2p2.gendf()
jef2p2.draglib()

jef2p2.hmat = "Co59"
jef2p2.mat = 2725
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape2"
jef2p2.pendf()
jef2p2.gendf()
jef2p2.draglib()

jef2p2.hmat = "Ni58"
jef2p2.mat = 2825
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape2"
jef2p2.pendf()
jef2p2.gendf()
jef2p2.draglib()

jef2p2.hmat = "Ni60"
jef2p2.mat = 2831
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape2"
jef2p2.pendf()
jef2p2.gendf()
jef2p2.draglib()

jef2p2.hmat = "Ni61"
jef2p2.mat = 2834
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape2"
jef2p2.pendf()
jef2p2.gendf()
jef2p2.draglib()

jef2p2.hmat = "Ni62"
jef2p2.mat = 2837
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape2"
jef2p2.pendf()
jef2p2.gendf()
jef2p2.draglib()

jef2p2.hmat = "Ni64"
jef2p2.mat = 2843
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape2"
jef2p2.pendf()
jef2p2.gendf()
jef2p2.draglib()

jef2p2.hmat = "Cu0"
jef2p2.mat = 2900
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape2"
jef2p2.pendf()
jef2p2.gendf()
jef2p2.draglib()

# Reactions mt=22 (n,na) and 28 (n,np) were manually removed from evaluation.
jef2p2.hmat = "Zn64"
jef2p2.mat = 3025
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape2"
jef2p2.pendf()
jef2p2.gendf()
jef2p2.draglib()

jef2p2.hmat = "W182"
jef2p2.mat =  7431  
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape6"
jef2p2.pendf()
jef2p2.gendf()
jef2p2.draglib()

jef2p2.hmat = "W183"
jef2p2.mat = 7434 
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape6"
jef2p2.pendf()
jef2p2.gendf()
jef2p2.draglib()

jef2p2.hmat = "W184"
jef2p2.mat = 7437
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape6"
jef2p2.pendf()
jef2p2.gendf()
jef2p2.draglib()

jef2p2.hmat = "W186"
jef2p2.mat =  7443
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape6"
jef2p2.pendf()
jef2p2.gendf()
jef2p2.draglib()

jef2p2.hmat = "Pb0"
jef2p2.mat =  8200 
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape6"
jef2p2.pendf()
jef2p2.gendf()
jef2p2.draglib()

jef2p2.hmat = "Bi209"
jef2p2.mat =  8325
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape6"
jef2p2.pendf()
jef2p2.gendf()
jef2p2.draglib()

jef2p2.hmat = "Zr0"
jef2p2.mat = 4000
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape3"
jef2p2.fission = None
jef2p2.ss = (22.53556, 1.858471e4)
jef2p2.potential = 6.8813
jef2p2.dilutions = ( 1.e10, 10000.0, 3549.18335, 1259.67004, 447.079956, \
158.676849, 56.3173141, 19.9880447, 7.09412289, 2.51783395 )
jef2p2.pendf()
jef2p2.gendf()
jef2p2.draglib()

jef2p2.hmat = "Th230"
jef2p2.mat = 9034
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape7"
jef2p2.fission = 1 # fission without delayed neutrons
jef2p2.ss = (22.53556, 3.206464e5)
jef2p2.potential = 8.7040
jef2p2.eFiss = 1.7406E+02
jef2p2.dilutions = ( 1.e10, 10000.0, 3549.18335, 1259.67004, 447.079956, \
158.676849, 56.3173141, 19.9880447, 7.09412289, 2.51783395 )
jef2p2.pendf()
jef2p2.gendf()
jef2p2.draglib()

jef2p2.hmat = "Th232"
jef2p2.mat = 9040
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape7"
jef2p2.fission = 2 # fission with delayed neutrons
jef2p2.ss = (22.53556, 3.206464e5)
jef2p2.potential = 11.8699
jef2p2.eFiss = 1.7193E+02
jef2p2.dilutions = ( 1.e10, 94.5317612, 56.3173141, 33.5510521, 19.9880447, \
11.9078817, 7.09412289, 4.22632504, 2.51783395, 1.5 )
jef2p2.pendf()
jef2p2.gendf()
jef2p2.draglib()
jef2p2.eFiss = 1.7193E+02
jef2p2.dilutions = ( 1.e10, 10000.0, 5957.50244, 3549.18335, 2114.42676, \
1259.67004, 750.448669, 447.079956, 266.347961, 158.676849 )
jef2p2.pendf()
jef2p2.gendf()
jef2p2.draglib()

jef2p2.hmat = "Pa231"
jef2p2.mat = 9131
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape7"
jef2p2.fission = 1 # fission without delayed neutrons
jef2p2.ss = (22.53556, 3.206464e5)
jef2p2.potential = 8.7266
jef2p2.eFiss = 1.7E2
jef2p2.dilutions = ( 1.e10, 10000.0, 3549.18335, 1259.67004, 447.079956, \
158.676849, 56.3173141, 19.9880447, 7.09412289, 2.51783395 )
jef2p2.pendf()
jef2p2.gendf()
jef2p2.draglib()

jef2p2.hmat = "Pa233"
jef2p2.mat = 9137
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape7"
jef2p2.fission = 1 # fission without delayed neutrons
jef2p2.ss = (22.53556, 3.206464e5)
jef2p2.potential = 9.9950
jef2p2.eFiss = 1.7558E+02
jef2p2.dilutions = ( 1.e10, 10000.0, 3549.18335, 1259.67004, 447.079956, \
158.676849, 56.3173141, 19.9880447, 7.09412289, 2.51783395 )
jef2p2.pendf()
jef2p2.gendf()
jef2p2.draglib()

jef2p2.hmat = "U232"
jef2p2.mat = 9219
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape7"
jef2p2.fission = 1 # fission without delayed neutrons
jef2p2.ss = (22.53556, 3.206464e5)
jef2p2.potential = 12.0687
jef2p2.eFiss = 1.8E2
jef2p2.dilutions = ( 1.e10, 10000.0, 3549.18335, 1259.67004, 447.079956, \
158.676849, 56.3173141, 19.9880447, 7.09412289, 2.51783395 )
jef2p2.pendf()
jef2p2.gendf()
jef2p2.draglib()

jef2p2.hmat = "U233"
jef2p2.mat = 9222
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape7"
jef2p2.fission = 2 # fission with delayed neutrons
jef2p2.ss = (22.53556, 3.206464e5)
jef2p2.potential = 12.2989
jef2p2.eFiss = 1.8087E+02
jef2p2.dilutions = ( 1.e10, 10000.0, 3549.18335, 1259.67004, 447.079956, \
158.676849, 56.3173141, 19.9880447, 7.09412289, 2.51783395 )
jef2p2.pendf()
jef2p2.gendf()
jef2p2.draglib()

jef2p2.hmat = "U234"
jef2p2.mat = 9225
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape7"
jef2p2.fission = 1 # fission without delayed neutrons
jef2p2.ss = (22.53556, 3.206464e5)
jef2p2.potential = 10.0210
jef2p2.eFiss = 1.7946E+02
jef2p2.dilutions = ( 1.e10, 10000.0, 3549.18335, 1259.67004, 447.079956, \
158.676849, 56.3173141, 19.9880447, 7.09412289, 2.51783395 )
jef2p2.pendf()
jef2p2.gendf()
jef2p2.draglib()

jef2p2.hmat = "U235"
jef2p2.mat = 9228
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape7"
jef2p2.fission = 2 # fission with delayed neutrons
jef2p2.ss = (22.53556, 2.499908e4)
jef2p2.potential = 11.6070
jef2p2.eFiss = 1.8089E+02
jef2p2.dilutions = ( 1.e10, 94.5317612, 56.3173141, 33.5510521, 19.9880447, \
11.9078817, 7.09412289, 4.22632504, 2.51783395, 1.5 )
jef2p2.pendf()
jef2p2.gendf()
jef2p2.draglib()
jef2p2.eFiss = 1.8089E+02
jef2p2.dilutions = ( 1.e10, 10000.0, 5957.50244, 3549.18335, 2114.42676, \
1259.67004, 750.448669, 447.079956, 266.347961, 158.676849 )
jef2p2.pendf()
jef2p2.gendf()
jef2p2.draglib()

jef2p2.hmat = "U236"
jef2p2.mat = 9231
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape7"
jef2p2.fission = 1 # fission without delayed neutrons
jef2p2.ss = (22.53556, 3.206464e5)
jef2p2.potential = 10.9954
jef2p2.eFiss = 1.7951E+02
jef2p2.dilutions = ( 1.e10, 10000.0, 3549.18335, 1259.67004, 447.079956, \
158.676849, 56.3173141, 19.9880447, 7.09412289, 2.51783395 )
jef2p2.pendf()
jef2p2.gendf()
jef2p2.draglib()

jef2p2.hmat = "U237"
jef2p2.mat = 9234
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape7"
jef2p2.fission = 1 # fission without delayed neutrons
jef2p2.ss = (22.53556, 3.206464e5)
jef2p2.potential = 10.5000
jef2p2.eFiss = 1.8E2
jef2p2.dilutions = ( 1.e10, 10000.0, 3549.18335, 1259.67004, 447.079956, \
158.676849, 56.3173141, 19.9880447, 7.09412289, 2.51783395 )
jef2p2.pendf()
jef2p2.gendf()
jef2p2.draglib()

jef2p2.hmat = "U238"
jef2p2.mat = 9237
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape7"
jef2p2.fission = 2 # fission with delayed neutrons
jef2p2.ss = (22.53556, 3.206464e5)
jef2p2.potential = 11.1710
jef2p2.eFiss = 1.8133E+02
jef2p2.dilutions = ( 1.e10, 94.5317612, 56.3173141, 33.5510521, 19.9880447, \
11.9078817, 7.09412289, 4.22632504, 2.51783395, 1.5 )
jef2p2.pendf()
jef2p2.gendf()
jef2p2.draglib()
jef2p2.eFiss = 1.8133E+02
jef2p2.dilutions = ( 1.e10, 10000.0, 5957.50244, 3549.18335, 2114.42676, \
1259.67004, 750.448669, 447.079956, 266.347961, 158.676849 )
jef2p2.pendf()
jef2p2.gendf()
jef2p2.draglib()

jef2p2.hmat = "Np237"
jef2p2.mat = 9346
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape8"
jef2p2.fission = 1 # fission without delayed neutrons
jef2p2.ss = (22.53556, 3.206464e5)
jef2p2.potential = 11.4369
jef2p2.eFiss = 1.8368E+02
jef2p2.dilutions = ( 1.e10, 10000.0, 3549.18335, 1259.67004, 447.079956, \
158.676849, 56.3173141, 19.9880447, 7.09412289, 2.51783395 )
jef2p2.branchingN2N = 0.801
jef2p2.pendf()
jef2p2.gendf()
jef2p2.draglib()
jef2p2.branchingN2N = None

jef2p2.hmat = "Np238"
jef2p2.mat = 9349
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape8"
jef2p2.fission = 1 # fission without delayed neutrons
jef2p2.ss = (22.53556, 3.206464e5)
jef2p2.eFiss = 1.8E2
jef2p2.potential = 10.4000
jef2p2.dilutions = ( 1.e10, 10000.0, 3549.18335, 1259.67004, 447.079956, \
158.676849, 56.3173141 )
jef2p2.pendf()
jef2p2.gendf()
jef2p2.draglib()

jef2p2.hmat = "Np239"
jef2p2.mat = 9352
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape8"
jef2p2.fission = 1 # fission without delayed neutrons
jef2p2.ss = (22.53556, 3.206464e5)
jef2p2.potential = 10.4979
jef2p2.eFiss = 1.8E2
jef2p2.dilutions = ( 1.e10, 10000.0, 3549.18335, 1259.67004, 447.079956, \
158.676849, 56.3173141, 19.9880447, 7.09412289, 2.51783395 )
jef2p2.pendf()
jef2p2.gendf()
jef2p2.draglib()

jef2p2.hmat = "Pu236"
jef2p2.mat = 9428
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape8"
jef2p2.fission = 1 # fission without delayed neutrons
jef2p2.ss = (22.53556, 3.206464e5)
jef2p2.potential = 11.2458
jef2p2.eFiss = 1.8E2
jef2p2.dilutions = ( 1.e10, 10000.0, 3549.18335, 1259.67004, 447.079956, \
158.676849, 56.3173141, 19.9880447, 7.09412289, 2.51783395 )
jef2p2.pendf()
jef2p2.gendf()
jef2p2.draglib()

jef2p2.hmat = "Pu237"
jef2p2.mat = 9431
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape8"
jef2p2.fission = 1 # fission without delayed neutrons
jef2p2.ss = (22.53556, 3.206464e5)
jef2p2.potential = 10.5209
jef2p2.eFiss = 1.8444E+02
jef2p2.dilutions = ( 1.e10, 10000.0, 3549.18335, 1259.67004, 447.079956, \
158.676849, 56.3173141, 19.9880447, 7.09412289, 2.51783395 )
jef2p2.pendf()
jef2p2.gendf()
jef2p2.draglib()

jef2p2.hmat = "Pu238"
jef2p2.mat = 9434
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape8"
jef2p2.fission = 1 # fission without delayed neutrons
jef2p2.ss = (22.53556, 3.206464e5)
jef2p2.potential = 10.8897
jef2p2.eFiss = 1.8665E+02
jef2p2.dilutions = ( 1.e10, 94.5317612, 56.3173141, 33.5510521, 19.9880447, \
11.9078817, 7.09412289, 4.22632504, 2.51783395, 1.5 )
jef2p2.pendf()
jef2p2.gendf()
jef2p2.draglib()
jef2p2.eFiss = 1.8665E+02
jef2p2.dilutions = ( 1.e10, 10000.0, 5957.50244, 3549.18335, 2114.42676, \
1259.67004, 750.448669, 447.079956, 266.347961, 158.676849 )
jef2p2.pendf()
jef2p2.gendf()
jef2p2.draglib()

jef2p2.hmat = "Pu239"
jef2p2.mat = 9437
#--------------------------------------------
#Note: Use the CEA-modified version of Pu-239
#--------------------------------------------
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape101"
jef2p2.fission = 2 # fission with delayed neutrons
jef2p2.ss = (22.53556, 3.206464e5)
jef2p2.potential = 10.8897
jef2p2.eFiss = 1.8944E+02
jef2p2.dilutions = ( 1.e10, 158.887822, 100.279413, 63.2896882, \
39.9442369, 25.2101426, 15.9109633, 10.0419406, 6.33780423, 4.0 )
jef2p2.pendf()
jef2p2.gendf()
jef2p2.draglib()
jef2p2.eFiss = 1.8944E+02
jef2p2.dilutions = ( 1.e10, 10000.0, 6311.33494, 3983.29436, 2513.99016, \
1586.66319, 1001.39615, 632.014573, 398.885514, 251.749976 )
jef2p2.pendf()
jef2p2.gendf()
jef2p2.draglib()

jef2p2.hmat = "Pu240"
jef2p2.mat = 9440
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape8"
jef2p2.fission = 2 # fission with delayed neutrons
jef2p2.ss = (22.53556, 3.206464e5)
jef2p2.potential = 9.9091
jef2p2.eFiss = 1.8636E+02
jef2p2.dilutions = ( 1.e10, 94.5317612, 56.3173141, 33.5510521, 19.9880447, \
11.9078817, 7.09412289, 4.22632504, 2.51783395, 1.5 )
jef2p2.pendf()
jef2p2.gendf()
jef2p2.draglib()
jef2p2.eFiss = 1.8636E+02
jef2p2.dilutions = ( 1.e10, 10000.0, 5957.50244, 3549.18335, 2114.42676, \
1259.67004, 750.448669, 447.079956, 266.347961, 158.676849 )
jef2p2.pendf()
jef2p2.gendf()
jef2p2.draglib()

jef2p2.hmat = "Pu241"
jef2p2.mat = 9443
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape8"
jef2p2.fission = 2 # fission with delayed neutrons
jef2p2.ss = (22.53556, 3.206464e5)
jef2p2.potential = 11.2156
jef2p2.eFiss = 1.8900E+02
jef2p2.dilutions = ( 1.e10, 94.5317612, 56.3173141, 33.5510521, 19.9880447, \
11.9078817, 7.09412289, 4.22632504, 2.51783395, 1.5 )
jef2p2.pendf()
jef2p2.gendf()
jef2p2.draglib()
jef2p2.eFiss = 1.8900E+02
jef2p2.dilutions = ( 1.e10, 10000.0, 5957.50244, 3549.18335, 2114.42676, \
1259.67004, 750.448669, 447.079956, 266.347961, 158.676849 )
jef2p2.pendf()
jef2p2.gendf()
jef2p2.draglib()

jef2p2.hmat = "Pu242"
jef2p2.mat = 9446
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape8"
jef2p2.fission = 1 # fission without delayed neutrons
jef2p2.ss = (22.53556, 3.206464e5)
jef2p2.potential = 10.6961
jef2p2.eFiss = 1.8599E+02
jef2p2.dilutions = ( 1.e10,  469.546659, 258.807233, 142.650752, \
78.6270025, 43.3380508, 23.8872981, 13.1663284, 7.25708715, 4.0 )
jef2p2.pendf()
jef2p2.gendf()
jef2p2.draglib()
jef2p2.eFiss = 1.8599E+02
jef2p2.dilutions = ( 1.e10, 1.e5, 55118.5161, 30380.5178, 16745.2959, \
9229.76155, 5087.30922, 2804.05024, 1545.55137, 851.885253 )
jef2p2.pendf()
jef2p2.gendf()
jef2p2.draglib()

jef2p2.hmat = "Pu243"
jef2p2.mat = 9449
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape8"
jef2p2.fission = 1 # fission without delayed neutrons
jef2p2.ss = (22.53556, 3.206464e5)
jef2p2.potential = 10.2000
jef2p2.eFiss = 1.9E2
jef2p2.dilutions = ( 1.e10, 10000.0, 3549.18335, 1259.67004, 447.079956, \
158.676849, 56.3173141, 19.9880447, 7.09412289, 2.51783395 )
jef2p2.pendf()
jef2p2.gendf()
jef2p2.draglib()

jef2p2.hmat = "Am241"
jef2p2.mat = 9543
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape9"
jef2p2.fission = 1 # fission without delayed neutrons
jef2p2.ss = (22.53556, 3.206464e5)
jef2p2.potential = 11.0329
jef2p2.eFiss = 1.9083E+02
jef2p2.branchingNG = 0.112
jef2p2.dilutions = ( 1.e10, 94.5317612, 56.3173141, 33.5510521, 19.9880447, \
11.9078817, 7.09412289, 4.22632504, 2.51783395, 1.5 )
jef2p2.pendf()
jef2p2.gendf()
jef2p2.draglib()
jef2p2.eFiss = 1.9083E+02
jef2p2.branchingNG = 0.112
jef2p2.dilutions = ( 1.e10, 10000.0, 5957.50244, 3549.18335, 2114.42676, \
1259.67004, 750.448669, 447.079956, 266.347961, 158.676849 )
jef2p2.pendf()
jef2p2.gendf()
jef2p2.draglib()
jef2p2.branchingNG = None

jef2p2.hmat = "Am242"
jef2p2.mat = 9546
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape9"
jef2p2.fission = 1 # fission without delayed neutrons
jef2p2.ss = (22.53556, 3.206464e5)
jef2p2.potential = 10.2000
jef2p2.eFiss = 1.9E2
jef2p2.dilutions = ( 1.e10, 10000.0, 3549.18335, 1259.67004, 447.079956, \
158.676849, 56.3173141, 19.9880447, 7.09412289, 2.51783395 )
jef2p2.pendf()
jef2p2.gendf()
jef2p2.draglib()

jef2p2.hmat = "Am242m"
jef2p2.mat = 9547
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape9"
jef2p2.fission = 1 # fission without delayed neutrons
jef2p2.ss = (22.53556, 3.206464e5)
jef2p2.potential = 10.2000
jef2p2.eFiss = 1.8759E+02
jef2p2.dilutions = ( 1.e10, 10000.0, 3549.18335, 1259.67004, 447.079956, \
158.676849, 56.3173141 )
jef2p2.pendf()
jef2p2.gendf()
jef2p2.draglib()

jef2p2.hmat = "Am243"
jef2p2.mat = 9549
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape9"
jef2p2.fission = 1 # fission without delayed neutrons
jef2p2.ss = (22.53556, 3.206464e5)
jef2p2.potential = 11.8237
jef2p2.eFiss = 1.9025E+02
jef2p2.dilutions = ( 1.e10, 10000.0, 3549.18335, 1259.67004, 447.079956, \
158.676849, 56.3173141, 19.9880447, 7.09412289, 2.51783395 )
jef2p2.pendf()
jef2p2.gendf()
jef2p2.draglib()

jef2p2.hmat = "Cm241"
jef2p2.mat = 9628
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape9"
jef2p2.fission = 1 # fission without delayed neutrons
jef2p2.ss = (22.53556, 3.206464e5)
jef2p2.potential = 10.1788
jef2p2.eFiss = 1.9245E+02
jef2p2.dilutions = ( 1.e10, 10000.0, 3549.18335, 1259.67004, 447.079956, \
158.676849, 56.3173141 )
jef2p2.pendf()
jef2p2.gendf()
jef2p2.draglib()

jef2p2.hmat = "Cm242"
jef2p2.mat = 9631
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape9"
jef2p2.fission = 1 # fission without delayed neutrons
jef2p2.ss = (22.53556, 3.206464e5)
jef2p2.potential = 10.2000
jef2p2.eFiss = 1.9142E+02
jef2p2.dilutions = ( 1.e10, 10000.0, 3549.18335, 1259.67004, 447.079956, \
158.676849, 56.3173141 )
jef2p2.pendf()
jef2p2.gendf()
jef2p2.draglib()

jef2p2.hmat = "Cm243"
jef2p2.mat = 9634
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape9"
jef2p2.fission = 1 # fission without delayed neutrons
jef2p2.ss = (22.53556, 3.206464e5)
jef2p2.potential = 11.2832
jef2p2.eFiss = 1.9E2
jef2p2.dilutions = ( 1.e10, 10000.0, 3549.18335, 1259.67004, 447.079956, \
158.676849, 56.3173141 )
jef2p2.pendf()
jef2p2.gendf()
jef2p2.draglib()

jef2p2.hmat = "Cm244"
jef2p2.mat = 9637
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape9"
jef2p2.fission = 1 # fission without delayed neutrons
jef2p2.ss = (22.53556, 3.206464e5)
jef2p2.potential = 10.3200
jef2p2.eFiss = 1.9049E+02
jef2p2.dilutions = ( 1.e10, 10000.0, 3549.18335, 1259.67004, 447.079956, \
158.676849, 56.3173141 )
jef2p2.pendf()
jef2p2.gendf()
jef2p2.draglib()

jef2p2.hmat = "Cm245"
jef2p2.mat = 9640
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape9"
jef2p2.fission = 1 # fission without delayed neutrons
jef2p2.ss = (22.53556, 3.206464e5)
jef2p2.potential = 11.3900
jef2p2.eFiss = 1.9E2
jef2p2.dilutions = ( 1.e10, 10000.0, 3549.18335, 1259.67004, 447.079956, \
158.676849, 56.3173141 )
jef2p2.pendf()
jef2p2.gendf()
jef2p2.draglib()

jef2p2.hmat = "Cm246"
jef2p2.mat = 9643
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape9"
jef2p2.fission = 1 # fission without delayed neutrons
jef2p2.ss = (22.53556, 3.206464e5)
jef2p2.potential = 10.2758
jef2p2.eFiss = 1.9E2
jef2p2.dilutions = ( 1.e10, 10000.0, 3549.18335, 1259.67004, 447.079956, \
158.676849, 56.3173141 )
jef2p2.pendf()
jef2p2.gendf()
jef2p2.draglib()

jef2p2.hmat = "Cm247"
jef2p2.mat = 9646
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape9"
jef2p2.fission = 1 # fission without delayed neutrons
jef2p2.ss = (22.53556, 3.206464e5)
jef2p2.potential = 10.0000
jef2p2.eFiss = 1.9E2
jef2p2.dilutions = ( 1.e10, 10000.0, 3549.18335, 1259.67004, 447.079956, \
158.676849, 56.3173141 )
jef2p2.pendf()
jef2p2.gendf()
jef2p2.draglib()

jef2p2.hmat = "Cm248"
jef2p2.mat = 9649
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape9"
jef2p2.fission = 1 # fission without delayed neutrons
jef2p2.ss = (22.53556, 3.206464e5)
jef2p2.potential = 10.3971
jef2p2.eFiss = 1.8885E+02
jef2p2.dilutions = ( 1.e10, 10000.0, 3549.18335, 1259.67004, 447.079956, \
158.676849, 56.3173141 )
jef2p2.pendf()
jef2p2.gendf()
jef2p2.draglib()

jef2p2.hmat = "Bk249"
jef2p2.mat = 9752
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape9"
jef2p2.fission = 1 # fission without delayed neutrons
jef2p2.ss = (22.53556, 3.206464e5)
jef2p2.potential = 11.0983
jef2p2.eFiss = 1.9E2
jef2p2.dilutions = ( 1.e10, 10000.0, 3549.18335, 1259.67004, 447.079956, \
158.676849, 56.3173141 )
jef2p2.pendf()
jef2p2.gendf()
jef2p2.draglib()

jef2p2.hmat = "Cf249"
jef2p2.mat = 9852
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape9"
jef2p2.fission = 1 # fission without delayed neutrons
jef2p2.ss = (22.53556, 3.206464e5)
jef2p2.potential = 11.1510
jef2p2.eFiss = 1.9E2
jef2p2.dilutions = ( 1.e10, 10000.0, 3549.18335, 1259.67004, 447.079956, \
158.676849, 56.3173141 )
jef2p2.pendf()
jef2p2.gendf()
jef2p2.draglib()

jef2p2.hmat = "Cf250"
jef2p2.mat = 9855
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape9"
jef2p2.fission = 1 # fission without delayed neutrons
jef2p2.ss = (22.53556, 3.206464e5)
jef2p2.potential = 9.8800
jef2p2.eFiss = 1.9E2
jef2p2.dilutions = ( 1.e10, 10000.0, 3549.18335, 1259.67004, 447.079956, \
158.676849, 56.3173141 )
jef2p2.pendf()
jef2p2.gendf()
jef2p2.draglib()

jef2p2.hmat = "Cf251"
jef2p2.mat = 9858
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape9"
jef2p2.fission = 1 # fission without delayed neutrons
jef2p2.ss = (22.53556, 3.206464e5)
jef2p2.potential = 9.8800
jef2p2.eFiss = 1.9E2
jef2p2.dilutions = ( 1.e10, 10000.0, 3549.18335, 1259.67004, 447.079956, \
158.676849, 56.3173141 )
jef2p2.pendf()
jef2p2.gendf()
jef2p2.draglib()

jef2p2.hmat = "Cf252"
jef2p2.mat = 9861
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape9"
jef2p2.fission = 1 # fission without delayed neutrons
jef2p2.ss = (22.53556, 3.206464e5)
jef2p2.potential = 9.8000
jef2p2.eFiss = 1.9E2
jef2p2.dilutions = ( 1.e10, 10000.0, 3549.18335, 1259.67004, 447.079956, \
158.676849, 56.3173141 )
jef2p2.pendf()
jef2p2.gendf()
jef2p2.draglib()

jef2p2.hmat = "Cf253"
jef2p2.mat = 9864
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape9"
jef2p2.fission = 1 # fission without delayed neutrons
jef2p2.ss = (22.53556, 3.206464e5)
jef2p2.potential = 9.7600
jef2p2.eFiss = 1.9E2
jef2p2.dilutions = ( 1.e10, 10000.0, 3549.18335, 1259.67004, 447.079956, \
158.676849, 56.3173141 )
jef2p2.pendf()
jef2p2.gendf()
jef2p2.draglib()

# Process the fission products:

jef2p2.scatteringLaw = None
jef2p2.legendre = 0
jef2p2.fission = None
jef2p2.dilutions = None
jef2p2.eFiss = None

jef2p2.hmat = "Ge72"
jef2p2.mat = 3231
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape3"
jef2p2.makeFp()

jef2p2.hmat = "Ge73"
jef2p2.mat = 3234
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape3"
jef2p2.makeFp()

jef2p2.hmat = "Ge74"
jef2p2.mat = 3237
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape3"
jef2p2.makeFp()

jef2p2.hmat = "Ge76"
jef2p2.mat = 3243
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape3"
jef2p2.makeFp()

jef2p2.hmat = "As75"
jef2p2.mat = 3325
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape3"
jef2p2.makeFp()

jef2p2.hmat = "Se76"
jef2p2.mat = 3431
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape3"
jef2p2.makeFp()

jef2p2.hmat = "Se77"
jef2p2.mat = 3434
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape3"
jef2p2.makeFp()

jef2p2.hmat = "Se78"
jef2p2.mat = 3437
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape3"
jef2p2.makeFp()

jef2p2.hmat = "Se80"
jef2p2.mat = 3443
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape3"
jef2p2.makeFp()

jef2p2.hmat = "Se82"
jef2p2.mat = 3449
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape3"
jef2p2.makeFp()

jef2p2.hmat = "Br79"
jef2p2.mat = 3525
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape3"
jef2p2.makeFp()

jef2p2.hmat = "Br81"
jef2p2.mat = 3531
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape3"
jef2p2.makeFp()

jef2p2.hmat = "Kr80"
jef2p2.mat = 3631
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape3"
jef2p2.makeFp()

jef2p2.hmat = "Kr82"
jef2p2.mat = 3637
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape3"
jef2p2.makeFp()

jef2p2.hmat = "Kr83"
jef2p2.mat = 3640
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape3"
jef2p2.makeFp()

jef2p2.hmat = "Kr84"
jef2p2.mat = 3643
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape3"
jef2p2.makeFp()

jef2p2.hmat = "Kr85"
jef2p2.mat = 3646
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape3"
jef2p2.makeFp()

jef2p2.hmat = "Kr86"
jef2p2.mat = 3649
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape3"
jef2p2.makeFp()

jef2p2.hmat = "Rb85"
jef2p2.mat = 3725
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape3"
jef2p2.makeFp()

jef2p2.hmat = "Rb87"
jef2p2.mat = 3731
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape3"
jef2p2.makeFp()

jef2p2.hmat = "Sr86"
jef2p2.mat = 3831
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape3"
jef2p2.makeFp()

jef2p2.hmat = "Sr87"
jef2p2.mat = 3834
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape3"
jef2p2.makeFp()

jef2p2.hmat = "Sr88"
jef2p2.mat = 3837
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape3"
jef2p2.makeFp()

jef2p2.hmat = "Sr89"
jef2p2.mat = 3840
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape3"
jef2p2.makeFp()

jef2p2.hmat = "Sr90"
jef2p2.mat = 3843
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape3"
jef2p2.makeFp()

jef2p2.hmat = "Y89"
jef2p2.mat = 3925
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape3"
jef2p2.makeFp()

jef2p2.hmat = "Y90"
jef2p2.mat = 3928
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape3"
jef2p2.makeFp()

jef2p2.hmat = "Y91"
jef2p2.mat = 3931
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape3"
jef2p2.makeFp()

jef2p2.hmat = "Zr90"
jef2p2.mat = 4025
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape3"
jef2p2.makeFp()

jef2p2.hmat = "Zr91"
jef2p2.mat = 4028
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape3"
jef2p2.makeFp()

jef2p2.hmat = "Zr92"
jef2p2.mat = 4031
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape3"
jef2p2.makeFp()

jef2p2.hmat = "Zr93"
jef2p2.mat = 4034
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape3"
jef2p2.makeFp()

jef2p2.hmat = "Zr94"
jef2p2.mat = 4037
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape3"
jef2p2.makeFp()

jef2p2.hmat = "Zr95"
jef2p2.mat = 4040
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape3"
jef2p2.makeFp()

jef2p2.hmat = "Zr96"
jef2p2.mat = 4043
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape3"
jef2p2.makeFp()

jef2p2.hmat = "Nb93"
jef2p2.mat = 4125
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape3"
jef2p2.makeFp()

jef2p2.hmat = "Nb94"
jef2p2.mat = 4128
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape3"
jef2p2.makeFp()

jef2p2.hmat = "Nb95"
jef2p2.mat = 4131
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape3"
jef2p2.makeFp()

jef2p2.hmat = "Mo0"
jef2p2.mat = 4200
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape3"
jef2p2.ss = (22.53556, 3.206464e5)
jef2p2.potential = 6.1140
jef2p2.dilutions = ( 1.e10, 10000.0, 3549.18335, 1259.67004, 447.079956, \
158.676849, 56.3173141, 19.9880447, 7.09412289, 2.51783395 )
jef2p2.pendf()
jef2p2.gendf()
jef2p2.draglib()

jef2p2.hmat = "Mo95"
jef2p2.mat = 4234
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape3"
jef2p2.makeFp()

jef2p2.hmat = "Mo96"
jef2p2.mat = 4237
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape3"
jef2p2.makeFp()

jef2p2.hmat = "Mo97"
jef2p2.mat = 4240
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape3"
jef2p2.makeFp()

jef2p2.hmat = "Mo98"
jef2p2.mat = 4243
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape3"
jef2p2.makeFp()

jef2p2.hmat = "Mo99"
jef2p2.mat = 4246
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape3"
jef2p2.makeFp()

jef2p2.hmat = "Mo100"
jef2p2.mat = 4249
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape3"
jef2p2.makeFp()

jef2p2.hmat = "Tc99"
jef2p2.mat = 4331
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape3"
jef2p2.makeFp()

jef2p2.hmat = "Ru99"
jef2p2.mat = 4434
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape3"
jef2p2.makeFp()

jef2p2.hmat = "Ru100"
jef2p2.mat = 4437
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape3"
jef2p2.makeFp()

jef2p2.hmat = "Ru101"
jef2p2.mat = 4440
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape3"
jef2p2.makeFp()

jef2p2.hmat = "Ru102"
jef2p2.mat = 4443
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape3"
jef2p2.makeFp()

jef2p2.hmat = "Ru103"
jef2p2.mat = 4446
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape3"
jef2p2.makeFp()

jef2p2.hmat = "Ru104"
jef2p2.mat = 4449
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape3"
jef2p2.makeFp()

jef2p2.hmat = "Ru105"
jef2p2.mat = 4452
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape3"
jef2p2.makeFp()

jef2p2.hmat = "Ru106"
jef2p2.mat = 4455
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape3"
jef2p2.makeFp()

jef2p2.hmat = "Rh103"
jef2p2.mat = 4525
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape3"
jef2p2.makeFp()

jef2p2.hmat = "Rh105"
jef2p2.mat = 4531
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape3"
jef2p2.makeFp()

jef2p2.hmat = "Pd104"
jef2p2.mat = 4631
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape3"
jef2p2.makeFp()

jef2p2.hmat = "Pd105"
jef2p2.mat = 4634
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape3"
jef2p2.makeFp()

jef2p2.hmat = "Pd106"
jef2p2.mat = 4637
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape3"
jef2p2.makeFp()

jef2p2.hmat = "Pd107"
jef2p2.mat = 4640
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape3"
jef2p2.makeFp()

jef2p2.hmat = "Pd108"
jef2p2.mat = 4643
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape3"
jef2p2.makeFp()

jef2p2.hmat = "Pd110"
jef2p2.mat = 4649
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape3"
jef2p2.makeFp()

jef2p2.hmat = "Ag107"
jef2p2.mat = 4725
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape3"
jef2p2.ss = (22.53556, 1.858471e4)
jef2p2.potential = 5.4739
jef2p2.dilutions = ( 1.e10, 10000.0, 4216.96552, 1778.27959, 749.894278, \
316.227791, 133.352152,  56.2341357, 23.7137381, 10.0 )
jef2p2.pendf()
jef2p2.gendf()
jef2p2.draglib()

jef2p2.hmat = "Ag109"
jef2p2.mat = 4731
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape3"
jef2p2.ss = (22.53556, 1.858471e4)
jef2p2.potential = 5.4739
jef2p2.dilutions = ( 1.e10, 10000.0, 4216.96552, 1778.27959, 749.894278, \
316.227791, 133.352152,  56.2341357, 23.7137381, 10.0 )
jef2p2.pendf()
jef2p2.gendf()
jef2p2.draglib()

jef2p2.hmat = "Ag111"
jef2p2.mat = 4737
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape3"
jef2p2.makeFp()

jef2p2.hmat = "Cd110"
jef2p2.mat = 4837
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape4"
jef2p2.makeFp()

jef2p2.hmat = "Cd111"
jef2p2.mat = 4840
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape4"
jef2p2.makeFp()

jef2p2.hmat = "Cd112"
jef2p2.mat = 4843
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape4"
jef2p2.makeFp()

jef2p2.hmat = "Cd113"
jef2p2.mat = 4846
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape4"
jef2p2.makeFp()

jef2p2.hmat = "Cd114"
jef2p2.mat = 4849
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape4"
jef2p2.makeFp()

jef2p2.hmat = "Cd115m"
jef2p2.mat = 4852
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape4"
jef2p2.makeFp()

jef2p2.hmat = "Cd116"
jef2p2.mat = 4855
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape4"
jef2p2.makeFp()

jef2p2.hmat = "In113"
jef2p2.mat = 4925
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape4"
jef2p2.makeFp()

jef2p2.hmat = "In115"
jef2p2.mat = 4931
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape4"
jef2p2.ss = (22.53556, 1.858471e4)
jef2p2.potential = 5.0439
jef2p2.dilutions = ( 1.e10, 10000.0, 4216.96552, 1778.27959, 749.894278, \
316.227791, 133.352152,  56.2341357, 23.7137381, 10.0 )
jef2p2.pendf()
jef2p2.gendf()
jef2p2.draglib()

jef2p2.hmat = "Sn115"
jef2p2.mat = 5034
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape4"
jef2p2.makeFp()

jef2p2.hmat = "Sn116"
jef2p2.mat = 5037
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape4"
jef2p2.makeFp()

jef2p2.hmat = "Sn117"
jef2p2.mat = 5040
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape4"
jef2p2.makeFp()

jef2p2.hmat = "Sn118"
jef2p2.mat = 5043
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape4"
jef2p2.makeFp()

jef2p2.hmat = "Sn119"
jef2p2.mat = 5046
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape4"
jef2p2.makeFp()

jef2p2.hmat = "Sn120"
jef2p2.mat = 5049
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape4"
jef2p2.makeFp()

jef2p2.hmat = "Sn122"
jef2p2.mat = 5055
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape4"
jef2p2.makeFp()

jef2p2.hmat = "Sn123"
jef2p2.mat = 5058
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape4"
jef2p2.makeFp()

jef2p2.hmat = "Sn124"
jef2p2.mat = 5061
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape4"
jef2p2.makeFp()

jef2p2.hmat = "Sn125"
jef2p2.mat = 5064
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape4"
jef2p2.makeFp()

jef2p2.hmat = "Sn126"
jef2p2.mat = 5067
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape4"
jef2p2.makeFp()

jef2p2.hmat = "Sb121"
jef2p2.mat = 5125
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape4"
jef2p2.makeFp()

jef2p2.hmat = "Sb123"
jef2p2.mat = 5131
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape4"
jef2p2.makeFp()

jef2p2.hmat = "Sb124"
jef2p2.mat = 5134
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape4"
jef2p2.makeFp()

jef2p2.hmat = "Sb125"
jef2p2.mat = 5137
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape4"
jef2p2.makeFp()

jef2p2.hmat = "Sb126"
jef2p2.mat = 5140
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape4"
jef2p2.makeFp()

jef2p2.hmat = "Te122"
jef2p2.mat = 5231
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape4"
jef2p2.makeFp()

jef2p2.hmat = "Te123"
jef2p2.mat = 5234
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape4"
jef2p2.makeFp()

jef2p2.hmat = "Te124"
jef2p2.mat = 5237
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape4"
jef2p2.makeFp()

jef2p2.hmat = "Te125"
jef2p2.mat = 5240
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape4"
jef2p2.makeFp()

jef2p2.hmat = "Te126"
jef2p2.mat = 5243
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape4"
jef2p2.makeFp()

jef2p2.hmat = "Te127m"
jef2p2.mat = 5246
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape4"
jef2p2.makeFp()

jef2p2.hmat = "Te128"
jef2p2.mat = 5249
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape4"
jef2p2.makeFp()

jef2p2.hmat = "Te129m"
jef2p2.mat = 5252
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape4"
jef2p2.makeFp()

jef2p2.hmat = "Te130"
jef2p2.mat = 5255
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape4"
jef2p2.makeFp()

jef2p2.hmat = "Te132"
jef2p2.mat = 5261
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape4"
jef2p2.makeFp()

jef2p2.hmat = "I127"
jef2p2.mat = 5325
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape4"
jef2p2.makeFp()

jef2p2.hmat = "I129"
jef2p2.mat = 5331
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape4"
jef2p2.makeFp()

jef2p2.hmat = "I130"
jef2p2.mat = 5334
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape4"
jef2p2.makeFp()

jef2p2.hmat = "I131"
jef2p2.mat = 5337
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape4"
jef2p2.makeFp()

jef2p2.hmat = "I135"
jef2p2.mat = 5349
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape4"
jef2p2.makeFp()

jef2p2.hmat = "Xe128"
jef2p2.mat = 5437
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape4"
jef2p2.makeFp()

jef2p2.hmat = "Xe129"
jef2p2.mat = 5440
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape4"
jef2p2.makeFp()

jef2p2.hmat = "Xe130"
jef2p2.mat = 5443
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape4"
jef2p2.makeFp()

jef2p2.hmat = "Xe131"
jef2p2.mat = 5446
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape4"
jef2p2.makeFp()

jef2p2.hmat = "Xe132"
jef2p2.mat = 5449
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape4"
jef2p2.makeFp()

jef2p2.hmat = "Xe133"
jef2p2.mat = 5452
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape4"
jef2p2.makeFp()

jef2p2.hmat = "Xe134"
jef2p2.mat = 5455
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape4"
jef2p2.makeFp()

jef2p2.hmat = "Xe135"
jef2p2.mat = 5458
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape4"
jef2p2.makeFp()

jef2p2.hmat = "Xe136"
jef2p2.mat = 5461
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape4"
jef2p2.makeFp()

jef2p2.hmat = "Cs133"
jef2p2.mat = 5525
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape4"
jef2p2.makeFp()

jef2p2.hmat = "Cs134"
jef2p2.mat = 5528
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape4"
jef2p2.makeFp()

jef2p2.hmat = "Cs135"
jef2p2.mat = 5531
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape4"
jef2p2.makeFp()

jef2p2.hmat = "Cs136"
jef2p2.mat = 5534
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape4"
jef2p2.makeFp()

jef2p2.hmat = "Cs137"
jef2p2.mat = 5537
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape4"
jef2p2.makeFp()

jef2p2.hmat = "Ba134"
jef2p2.mat = 5637
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape4"
jef2p2.makeFp()

jef2p2.hmat = "Ba135"
jef2p2.mat = 5640
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape4"
jef2p2.makeFp()

jef2p2.hmat = "Ba136"
jef2p2.mat = 5643
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape4"
jef2p2.makeFp()

jef2p2.hmat = "Ba137"
jef2p2.mat = 5646
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape4"
jef2p2.makeFp()

jef2p2.hmat = "Ba138"
jef2p2.mat = 5649
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape4"
jef2p2.makeFp()

jef2p2.hmat = "Ba140"
jef2p2.mat = 5655
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape4"
jef2p2.makeFp()

jef2p2.hmat = "La139"
jef2p2.mat = 5728
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape4"
jef2p2.makeFp()

jef2p2.hmat = "La140"
jef2p2.mat = 5731
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape4"
jef2p2.makeFp()

jef2p2.hmat = "Ce140"
jef2p2.mat = 5837
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape4"
jef2p2.makeFp()

jef2p2.hmat = "Ce141"
jef2p2.mat = 5840
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape4"
jef2p2.makeFp()

jef2p2.hmat = "Ce142"
jef2p2.mat = 5843
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape4"
jef2p2.makeFp()

jef2p2.hmat = "Ce143"
jef2p2.mat = 5846
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape4"
jef2p2.makeFp()

jef2p2.hmat = "Ce144"
jef2p2.mat = 5849
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape4"
jef2p2.makeFp()

jef2p2.hmat = "Pr141"
jef2p2.mat = 5925
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape4"
jef2p2.makeFp()

jef2p2.hmat = "Pr142"
jef2p2.mat = 5928
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape4"
jef2p2.makeFp()

jef2p2.hmat = "Pr143"
jef2p2.mat = 5931
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape4"
jef2p2.makeFp()

jef2p2.hmat = "Nd142"
jef2p2.mat = 6025
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape5"
jef2p2.makeFp()

jef2p2.hmat = "Nd143"
jef2p2.mat = 6028
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape5"
jef2p2.makeFp()

jef2p2.hmat = "Nd144"
jef2p2.mat = 6031
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape5"
jef2p2.makeFp()

jef2p2.hmat = "Nd145"
jef2p2.mat = 6034
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape5"
jef2p2.makeFp()

jef2p2.hmat = "Nd146"
jef2p2.mat = 6037
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape5"
jef2p2.makeFp()

jef2p2.hmat = "Nd147"
jef2p2.mat = 6040
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape5"
jef2p2.makeFp()

jef2p2.hmat = "Nd148"
jef2p2.mat = 6043
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape5"
jef2p2.makeFp()

jef2p2.hmat = "Nd150"
jef2p2.mat = 6049
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape5"
jef2p2.makeFp()

jef2p2.hmat = "Pm147"
jef2p2.mat = 6149
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape5"
jef2p2.branchingNG = 0.47
jef2p2.makeFp()
jef2p2.branchingNG = None

jef2p2.hmat = "Pm148"
jef2p2.mat = 6152
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape5"
jef2p2.makeFp()

jef2p2.hmat = "Pm148m"
jef2p2.mat = 6153
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape5"
jef2p2.makeFp()

jef2p2.hmat = "Pm149"
jef2p2.mat = 6155
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape5"
jef2p2.makeFp()

jef2p2.hmat = "Pm151"
jef2p2.mat = 6161
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape5"
jef2p2.makeFp()

jef2p2.hmat = "Sm147"
jef2p2.mat = 6234
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape5"
jef2p2.makeFp()

jef2p2.hmat = "Sm148"
jef2p2.mat = 6237
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape5"
jef2p2.makeFp()

jef2p2.hmat = "Sm149"
jef2p2.mat = 6240
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape5"
jef2p2.makeFp()

jef2p2.hmat = "Sm150"
jef2p2.mat = 6243
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape5"
jef2p2.makeFp()

jef2p2.hmat = "Sm151"
jef2p2.mat = 6246
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape5"
jef2p2.makeFp()

jef2p2.hmat = "Sm152"
jef2p2.mat = 6249
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape5"
jef2p2.makeFp()

jef2p2.hmat = "Sm153"
jef2p2.mat = 6252
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape5"
jef2p2.makeFp()

jef2p2.hmat = "Sm154"
jef2p2.mat = 6255
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape5"
jef2p2.makeFp()

jef2p2.hmat = "Eu151"
jef2p2.mat = 6325
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape5"
jef2p2.makeFp()

jef2p2.hmat = "Eu152"
jef2p2.mat = 6328
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape5"
jef2p2.makeFp()

jef2p2.hmat = "Eu153"
jef2p2.mat = 6331
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape5"
jef2p2.makeFp()

jef2p2.hmat = "Eu154"
jef2p2.mat = 6334
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape5"
jef2p2.makeFp()

jef2p2.hmat = "Eu155"
jef2p2.mat = 6337
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape5"
jef2p2.makeFp()

jef2p2.hmat = "Eu156"
jef2p2.mat = 6340
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape5"
jef2p2.makeFp()

jef2p2.hmat = "Eu157"
jef2p2.mat = 6343
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape5"
jef2p2.makeFp()

jef2p2.hmat = "Gd154"
jef2p2.mat = 6431
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape5"
jef2p2.ss = (22.53556, 1.858471e4)
jef2p2.potential = 7.6587
jef2p2.dilutions = ( 1.e10, 10000.0, 4216.96552, 1778.27959, 749.894278, \
316.227791, 133.352152,  56.2341357, 23.7137381, 10.0 )
jef2p2.pendf()
jef2p2.gendf()
jef2p2.draglib()

jef2p2.hmat = "Gd155"
jef2p2.mat = 6434
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape5"
jef2p2.ss = (22.53556, 1.858471e4)
jef2p2.potential = 8.0425
jef2p2.dilutions = ( 1.e10, 10000.0, 4216.96552, 1778.27959, 749.894278, \
316.227791, 133.352152,  56.2341357, 23.7137381, 10.0 )
jef2p2.pendf()
jef2p2.gendf()
jef2p2.draglib()

jef2p2.hmat = "Gd156"
jef2p2.mat = 6437
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape5"
jef2p2.ss = (22.53556, 1.858471e4)
jef2p2.potential = 8.2448
jef2p2.dilutions = ( 1.e10, 10000.0, 4216.96552, 1778.27959, 749.894278, \
316.227791, 133.352152,  56.2341357, 23.7137381, 10.0 )
jef2p2.pendf()
jef2p2.gendf()
jef2p2.draglib()

jef2p2.hmat = "Gd157"
jef2p2.mat = 6440
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape5"
jef2p2.ss = (22.53556, 1.858471e4)
jef2p2.potential = 7.6454
jef2p2.dilutions = ( 1.e10, 10000.0, 4216.96552, 1778.27959, 749.894278, \
316.227791, 133.352152,  56.2341357, 23.7137381, 10.0 )
jef2p2.pendf()
jef2p2.gendf()
jef2p2.draglib()

jef2p2.hmat = "Gd158"
jef2p2.mat = 6443
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape5"
jef2p2.ss = (22.53556, 1.858471e4)
jef2p2.potential = 5.3093
jef2p2.dilutions = ( 1.e10, 10000.0, 4216.96552, 1778.27959, 749.894278, \
316.227791, 133.352152,  56.2341357, 23.7137381, 10.0 )
jef2p2.pendf()
jef2p2.gendf()
jef2p2.draglib()

jef2p2.hmat = "Gd160"
jef2p2.mat = 6449
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape5"
jef2p2.ss = (22.53556, 1.858471e4)
jef2p2.potential = 5.8107
jef2p2.dilutions = ( 1.e10, 10000.0, 4216.96552, 1778.27959, 749.894278, \
316.227791, 133.352152,  56.2341357, 23.7137381, 10.0 )
jef2p2.pendf()
jef2p2.gendf()
jef2p2.draglib()

jef2p2.dilutions = None

jef2p2.hmat = "Tb159"
jef2p2.mat = 6525
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape5"
jef2p2.makeFp()

jef2p2.hmat = "Tb160"
jef2p2.mat = 6528
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape5"
jef2p2.makeFp()

jef2p2.hmat = "Dy160"
jef2p2.mat = 6637
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape6"
jef2p2.ss = (22.53556, 1.858471e4)
jef2p2.potential = 7.0686
jef2p2.dilutions = ( 1.e10, 10000.0, 4216.96552, 1778.27959, 749.894278, \
316.227791, 133.352152,  56.2341357, 23.7137381, 10.0 )
jef2p2.pendf()
jef2p2.gendf()
jef2p2.draglib()

jef2p2.hmat = "Dy161"
jef2p2.mat = 6640
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape6"
jef2p2.ss = (22.53556, 1.858471e4)
jef2p2.potential = 7.0686
jef2p2.dilutions = ( 1.e10, 10000.0, 4216.96552, 1778.27959, 749.894278, \
316.227791, 133.352152,  56.2341357, 23.7137381, 10.0 )
jef2p2.pendf()
jef2p2.gendf()
jef2p2.draglib()

jef2p2.hmat = "Dy162"
jef2p2.mat = 6643
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape6"
jef2p2.ss = (22.53556, 1.858471e4)
jef2p2.potential = 7.6846
jef2p2.dilutions = ( 1.e10, 10000.0, 4216.96552, 1778.27959, 749.894278, \
316.227791, 133.352152,  56.2341357, 23.7137381, 10.0 )
jef2p2.pendf()
jef2p2.gendf()
jef2p2.draglib()

jef2p2.hmat = "Dy163"
jef2p2.mat = 6646
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape6"
jef2p2.ss = (22.53556, 1.858471e4)
jef2p2.potential = 7.0686
jef2p2.dilutions = ( 1.e10, 10000.0, 4216.96552, 1778.27959, 749.894278, \
316.227791, 133.352152,  56.2341357, 23.7137381, 10.0 )
jef2p2.pendf()
jef2p2.gendf()
jef2p2.draglib()

jef2p2.hmat = "Dy164"
jef2p2.mat = 6649
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape6"
jef2p2.ss = (22.53556, 1.858471e4)
jef2p2.potential = 7.6846
jef2p2.dilutions = ( 1.e10, 10000.0, 4216.96552, 1778.27959, 749.894278, \
316.227791, 133.352152,  56.2341357, 23.7137381, 10.0 )
jef2p2.pendf()
jef2p2.gendf()
jef2p2.draglib()

jef2p2.hmat = "Hf174"
jef2p2.mat = 7225
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape6"
jef2p2.ss = (22.53556, 1.858471e4)
jef2p2.potential = 7.1385
jef2p2.dilutions = ( 1.e10, 10000.0, 4216.96552, 1778.27959, 749.894278, \
316.227791, 133.352152,  56.2341357, 23.7137381, 10.0 )
jef2p2.pendf()
jef2p2.gendf()
jef2p2.draglib()

jef2p2.hmat = "Hf176"
jef2p2.mat = 7231
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape6"
jef2p2.ss = (22.53556, 1.858471e4)
jef2p2.potential = 7.1935
jef2p2.dilutions = ( 1.e10, 10000.0, 4216.96552, 1778.27959, 749.894278, \
316.227791, 133.352152,  56.2341357, 23.7137381, 10.0 )
jef2p2.pendf()
jef2p2.gendf()
jef2p2.draglib()

jef2p2.hmat = "Hf177"
jef2p2.mat = 7234
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape6"
jef2p2.ss = (22.53556, 1.858471e4)
jef2p2.potential = 7.2202
jef2p2.dilutions = ( 1.e10, 10000.0, 4216.96552, 1778.27959, 749.894278, \
316.227791, 133.352152,  56.2341357, 23.7137381, 10.0 )
jef2p2.pendf()
jef2p2.gendf()
jef2p2.draglib()

jef2p2.hmat = "Hf178"
jef2p2.mat = 7237
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape6"
jef2p2.ss = (22.53556, 1.858471e4)
jef2p2.potential = 7.2469
jef2p2.dilutions = ( 1.e10, 10000.0, 4216.96552, 1778.27959, 749.894278, \
316.227791, 133.352152,  56.2341357, 23.7137381, 10.0 )
jef2p2.pendf()
jef2p2.gendf()
jef2p2.draglib()

jef2p2.hmat = "Hf179"
jef2p2.mat = 7240
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape6"
jef2p2.ss = (22.53556, 1.858471e4)
jef2p2.potential = 7.2736
jef2p2.dilutions = ( 1.e10, 10000.0, 4216.96552, 1778.27959, 749.894278, \
316.227791, 133.352152,  56.2341357, 23.7137381, 10.0 )
jef2p2.pendf()
jef2p2.gendf()
jef2p2.draglib()

jef2p2.hmat = "Hf180"
jef2p2.mat = 7243
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape6"
jef2p2.ss = (22.53556, 1.858471e4)
jef2p2.potential = 7.3004
jef2p2.dilutions = ( 1.e10, 10000.0, 4216.96552, 1778.27959, 749.894278, \
316.227791, 133.352152,  56.2341357, 23.7137381, 10.0 )
jef2p2.pendf()
jef2p2.gendf()
jef2p2.draglib()

jef2p2.dilutions = None

jef2p2.hmat = "Ho165"
jef2p2.mat = 6725
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape6"
jef2p2.makeFp()

jef2p2.hmat = "Er166"
jef2p2.mat = 6837
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape6"
jef2p2.makeFp()

jef2p2.hmat = "Er167"
jef2p2.mat = 6840
jef2p2.evaluationFile = "$HOME/evaluations/Jef2.2/tape6"
jef2p2.makeFp()

# Process the burnup chain:

jef2p2.fissionFile = "$HOME/evaluations/Jef2.2/tape24"
jef2p2.decayFile = "$HOME/evaluations/Jef2.2/tape22"
jef2p2.burnup()
