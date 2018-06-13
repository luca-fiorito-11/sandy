#!/usr/local/bin/python
from PyNjoy import *
from os import uname
endfb = PyNjoy()
endfb.evaluationName = "/tmp/shem281_endfb6r8"
endfb.execDir = "../" + uname()[0]
endfb.nstr = 24
endfb.iwt = 4
endfb.Espectra = None
endfb.autolib = (22.53556, 748.5173, 0.001)

endfb.legendre = 3
endfb.hmat = "H1_H2O"
endfb.mat = 125
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.scatteringLaw = "$HOME/evaluations/ENDFB6r8/TAPE.101"
endfb.scatteringMat = 1
endfb.temperatures = ( 296.0, 350.0, 400.0, 450.0, 500.0, 600.0 )
endfb.fission = None
endfb.dilutions = None
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.legendre = 1
endfb.hmat = "H2_D2O"
endfb.mat = 128
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.scatteringLaw = "$HOME/evaluations/ENDFB6r8/TAPE.101"
endfb.scatteringMat = 11
endfb.temperatures = ( 296.0, 350.0, 400.0, 450.0, 500.0, 600.0 )
endfb.fission = None
endfb.dilutions = None
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.scatteringLaw = None
endfb.temperatures = ( 293., 550., 900., 1200., 2000. )
endfb.fission = None
endfb.dilutions = None

endfb.hmat = "H1"
endfb.mat = 125
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "H2"
endfb.mat = 128
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "H3"
endfb.mat = 131
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "He3"
endfb.mat = 225
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "He4"
endfb.mat = 228
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Li6"
endfb.mat = 325
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Li7"
endfb.mat = 328
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Be9"
endfb.mat = 425
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "B10"
endfb.mat = 525
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "B11"
endfb.mat = 528
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "C0"
endfb.mat = 600
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "N14"
endfb.mat = 725
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "N15"
endfb.mat = 728
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "O16"
endfb.mat = 825
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "O17"
endfb.mat = 828
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "F19"
endfb.mat = 925
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Na23"
endfb.mat = 1125
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Mg0"
endfb.mat = 1200
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Al27"
endfb.mat = 1325
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Si0"
endfb.mat = 1400
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "P31"
endfb.mat = 1525
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "S0"
endfb.mat = 1600
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "S32"
endfb.mat = 1625
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Cl0"
endfb.mat = 1700
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "K0"
endfb.mat = 1900
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Ca0"
endfb.mat = 2000
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Ti0"
endfb.mat = 2200
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "V0"
endfb.mat = 2300
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Cr50"
endfb.mat = 2425
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Cr52"
endfb.mat = 2431
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Cr53"
endfb.mat = 2434
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Cr54"
endfb.mat = 2437
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Mn55"
endfb.mat = 2525
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Fe54"
endfb.mat = 2625
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Fe56"
endfb.mat = 2631
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Fe57"
endfb.mat = 2634
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Fe58"
endfb.mat = 2637
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Co59"
endfb.mat = 2725
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Ni58"
endfb.mat = 2825
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Ni60"
endfb.mat = 2831
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Ni61"
endfb.mat = 2834
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Ni62"
endfb.mat = 2837
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Ni64"
endfb.mat = 2843
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Cu63"
endfb.mat = 2925
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Cu65"
endfb.mat = 2931
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Pb206"
endfb.mat = 8231
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Pb207"
endfb.mat = 8234
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Pb208"
endfb.mat = 8237
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Bi209"
endfb.mat = 8325
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Zr0"
endfb.mat = 4000
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.fission = None
endfb.ss = (22.53556, 1.858471e4)
endfb.potential = 6.8813
endfb.dilutions = ( 1.e10, 10000.0, 3549.18335, 1259.67004, 447.079956, \
158.676849, 56.3173141, 19.9880447, 7.09412289, 2.51783395 )
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Th230"
endfb.mat = 9034
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.fission = 1 # fission without delayed neutrons
endfb.ss = (22.53556, 1.227732e5)
endfb.potential = 8.7040
endfb.eFiss = 1.7406E+02
endfb.dilutions = ( 1.e10, 10000.0, 3549.18335, 1259.67004, 447.079956, \
158.676849, 56.3173141, 19.9880447, 7.09412289, 2.51783395 )
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Th232"
endfb.mat = 9040
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.fission = 2 # fission with delayed neutrons
endfb.ss = (22.53556, 1.227732e5)
endfb.potential = 11.8699
endfb.eFiss = 1.7193E+02
endfb.dilutions = ( 1.e10, 94.5317612, 56.3173141, 33.5510521, 19.9880447, \
11.9078817, 7.09412289, 4.22632504, 2.51783395, 1.5 )
endfb.pendf()
endfb.gendf()
endfb.draglib()
endfb.eFiss = 1.7193E+02
endfb.dilutions = ( 1.e10, 10000.0, 5957.50244, 3549.18335, 2114.42676, \
1259.67004, 750.448669, 447.079956, 266.347961, 158.676849 )
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Pa231"
endfb.mat = 9131
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.fission = 1 # fission without delayed neutrons
endfb.ss = (22.53556, 1.227732e5)
endfb.potential = 8.7266
endfb.eFiss = 1.7E2
endfb.dilutions = ( 1.e10, 10000.0, 3549.18335, 1259.67004, 447.079956, \
158.676849, 56.3173141, 19.9880447, 7.09412289, 2.51783395 )
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Pa233"
endfb.mat = 9137
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.fission = 1 # fission without delayed neutrons
endfb.ss = (22.53556, 1.227732e5)
endfb.potential = 9.9950
endfb.eFiss = 1.7558E+02
endfb.dilutions = ( 1.e10, 10000.0, 3549.18335, 1259.67004, 447.079956, \
158.676849, 56.3173141, 19.9880447, 7.09412289, 2.51783395 )
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "U232"
endfb.mat = 9219
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.fission = 1 # fission without delayed neutrons
endfb.ss = (22.53556, 1.227732e5)
endfb.potential = 12.0687
endfb.eFiss = 1.8E2
endfb.dilutions = ( 1.e10, 10000.0, 3549.18335, 1259.67004, 447.079956, \
158.676849, 56.3173141, 19.9880447, 7.09412289, 2.51783395 )
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "U233"
endfb.mat = 9222
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.fission = 2 # fission with delayed neutrons
endfb.ss = (22.53556, 1.227732e5)
endfb.potential = 12.2989
endfb.eFiss = 1.8087E+02
endfb.dilutions = ( 1.e10, 10000.0, 3549.18335, 1259.67004, 447.079956, \
158.676849, 56.3173141, 19.9880447, 7.09412289, 2.51783395 )
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "U234"
endfb.mat = 9225
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.fission = 1 # fission without delayed neutrons
endfb.ss = (22.53556, 1.227732e5)
endfb.potential = 10.0210
endfb.dilutions = ( 1.e10, 10000.0, 3549.18335, 1259.67004, 447.079956, \
158.676849, 56.3173141, 19.9880447, 7.09412289, 2.51783395 )
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "U235"
endfb.mat = 9228
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.fission = 2 # fission with delayed neutrons
endfb.ss = (22.53556, 1.227732e5)
endfb.potential = 11.6070
endfb.dilutions = ( 1.e10, 94.5317612, 56.3173141, 33.5510521, 19.9880447, \
11.9078817, 7.09412289, 4.22632504, 2.51783395, 1.5 )
endfb.pendf()
endfb.gendf()
endfb.draglib()
endfb.dilutions = ( 1.e10, 10000.0, 5957.50244, 3549.18335, 2114.42676, \
1259.67004, 750.448669, 447.079956, 266.347961, 158.676849 )
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "U236"
endfb.mat = 9231
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.fission = 1 # fission without delayed neutrons
endfb.ss = (22.53556, 1.227732e5)
endfb.potential = 10.9954
endfb.dilutions = ( 1.e10, 10000.0, 3549.18335, 1259.67004, 447.079956, \
158.676849, 56.3173141, 19.9880447, 7.09412289, 2.51783395 )
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "U237"
endfb.mat = 9234
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.fission = 1 # fission without delayed neutrons
endfb.ss = (22.53556, 1.227732e5)
endfb.potential = 10.5000
endfb.dilutions = ( 1.e10, 10000.0, 3549.18335, 1259.67004, 447.079956, \
158.676849, 56.3173141, 19.9880447, 7.09412289, 2.51783395 )
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "U238"
endfb.mat = 9237
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.fission = 2 # fission with delayed neutrons
endfb.ss = (22.53556, 1.227732e5)
endfb.potential = 11.17103
endfb.dilutions = ( 1.e10, 94.5317612, 56.3173141, 33.5510521, 19.9880447, \
11.9078817, 7.09412289, 4.22632504, 2.51783395, 1.5 )
endfb.pendf()
endfb.gendf()
endfb.draglib()
endfb.dilutions = ( 1.e10, 10000.0, 5957.50244, 3549.18335, 2114.42676, \
1259.67004, 750.448669, 447.079956, 266.347961, 158.676849 )
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Np236"
endfb.mat = 9343
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.fission = 1 # fission without delayed neutrons
endfb.ss = (22.53556, 1.227732e5)
endfb.eFiss = 1.8E2
endfb.potential = 12.1427
endfb.dilutions = ( 1.e10, 10000.0, 3549.18335, 1259.67004, 447.079956, \
158.676849, 56.3173141, 19.9880447, 7.09412289, 2.51783395 )
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Np237"
endfb.mat = 9346
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.fission = 1 # fission without delayed neutrons
endfb.ss = (22.53556, 1.227732e5)
endfb.potential = 11.4369
endfb.branchingN2N = 0.801
endfb.dilutions = ( 1.e10, 10000.0, 3549.18335, 1259.67004, 447.079956, \
158.676849, 56.3173141, 19.9880447, 7.09412289, 2.51783395 )
endfb.pendf()
endfb.gendf()
endfb.draglib()
endfb.branchingN2N = None

endfb.hmat = "Np238"
endfb.mat = 9349
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.fission = 1 # fission without delayed neutrons
endfb.ss = (22.53556, 1.227732e5)
endfb.eFiss = 1.8E2
endfb.potential = 10.4000
endfb.dilutions = ( 1.e10, 10000.0, 3549.18335, 1259.67004, 447.079956, \
158.676849, 56.3173141 )
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Np239"
endfb.mat = 9352
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.fission = 1 # fission without delayed neutrons
endfb.ss = (22.53556, 1.227732e5)
endfb.eFiss = 1.8E2
endfb.potential = 10.4979
endfb.dilutions = ( 1.e10, 10000.0, 3549.18335, 1259.67004, 447.079956, \
158.676849, 56.3173141, 19.9880447, 7.09412289, 2.51783395 )
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Pu236"
endfb.mat = 9428
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.fission = 1 # fission without delayed neutrons
endfb.ss = (22.53556, 1.227732e5)
endfb.eFiss = 1.8E2
endfb.potential = 11.2458
endfb.dilutions = ( 1.e10, 10000.0, 3549.18335, 1259.67004, 447.079956, \
158.676849, 56.3173141, 19.9880447, 7.09412289, 2.51783395 )
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Pu237"
endfb.mat = 9431
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.fission = 1 # fission without delayed neutrons
endfb.ss = (22.53556, 1.227732e5)
endfb.potential = 10.5209
endfb.dilutions = ( 1.e10, 10000.0, 3549.18335, 1259.67004, 447.079956, \
158.676849, 56.3173141, 19.9880447, 7.09412289, 2.51783395 )
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Pu238"
endfb.mat = 9434
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.fission = 1 # fission without delayed neutrons
endfb.ss = (22.53556, 1.227732e5)
endfb.potential = 10.8897
endfb.dilutions = ( 1.e10, 94.5317612, 56.3173141, 33.5510521, 19.9880447, \
11.9078817, 7.09412289, 4.22632504, 2.51783395, 1.5 )
endfb.pendf()
endfb.gendf()
endfb.draglib()
endfb.dilutions = ( 1.e10, 10000.0, 5957.50244, 3549.18335, 2114.42676, \
1259.67004, 750.448669, 447.079956, 266.347961, 158.676849 )
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Pu239"
endfb.mat = 9437
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.fission = 2 # fission with delayed neutrons
endfb.ss = (22.53556, 1.227732e5)
endfb.potential = 10.8897
endfb.dilutions = ( 1.e10, 158.887822, 100.279413, 63.2896882, \
39.9442369, 25.2101426, 15.9109633, 10.0419406, 6.33780423, 4.0 )
endfb.pendf()
endfb.gendf()
endfb.draglib()
endfb.dilutions = ( 1.e10, 10000.0, 6311.33494, 3983.29436, 2513.99016, \
1586.66319, 1001.39615, 632.014573, 398.885514, 251.749976 )
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Pu240"
endfb.mat = 9440
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.fission = 2 # fission with delayed neutrons
endfb.ss = (22.53556, 1.227732e5)
endfb.potential = 9.9091
endfb.dilutions = ( 1.e10, 94.5317612, 56.3173141, 33.5510521, 19.9880447, \
11.9078817, 7.09412289, 4.22632504, 2.51783395, 1.5 )
endfb.pendf()
endfb.gendf()
endfb.draglib()
endfb.dilutions = ( 1.e10, 10000.0, 5957.50244, 3549.18335, 2114.42676, \
1259.67004, 750.448669, 447.079956, 266.347961, 158.676849 )
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Pu241"
endfb.mat = 9443
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.fission = 2 # fission with delayed neutrons
endfb.ss = (22.53556, 1.227732e5)
endfb.potential = 11.2156
endfb.dilutions = ( 1.e10, 94.5317612, 56.3173141, 33.5510521, 19.9880447, \
11.9078817, 7.09412289, 4.22632504, 2.51783395, 1.5 )
endfb.pendf()
endfb.gendf()
endfb.draglib()
endfb.dilutions = ( 1.e10, 10000.0, 5957.50244, 3549.18335, 2114.42676, \
1259.67004, 750.448669, 447.079956, 266.347961, 158.676849 )
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Pu242"
endfb.mat = 9446
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.fission = 1 # fission without delayed neutrons
endfb.ss = (22.53556, 1.227732e5)
endfb.potential = 10.6961
endfb.dilutions = ( 1.e10,  469.546659, 258.807233, 142.650752, \
78.6270025, 43.3380508, 23.8872981, 13.1663284, 7.25708715, 4.0 )
endfb.pendf()
endfb.gendf()
endfb.draglib()
endfb.dilutions = ( 1.e10, 1.e5, 55118.5161, 30380.5178, 16745.2959, \
9229.76155, 5087.30922, 2804.05024, 1545.55137, 851.885253 )
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Pu243"
endfb.mat = 9449
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.fission = 1 # fission without delayed neutrons
endfb.ss = (22.53556, 1.227732e5)
endfb.potential = 10.2000
endfb.eFiss = 1.9E2
endfb.dilutions = ( 1.e10, 10000.0, 3549.18335, 1259.67004, 447.079956, \
158.676849, 56.3173141, 19.9880447, 7.09412289, 2.51783395 )
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Pu244"
endfb.mat = 9452
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.fission = 1 # fission without delayed neutrons
endfb.ss = (22.53556, 1.227732e5)
endfb.potential = 10.1000
endfb.dilutions = ( 1.e10, 10000.0, 3549.18335, 1259.67004, 447.079956, \
158.676849, 56.3173141, 19.9880447, 7.09412289, 2.51783395 )
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Am241"
endfb.mat = 9543
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.fission = 1 # fission without delayed neutrons
endfb.ss = (22.53556, 1.227732e5)
endfb.potential = 11.0329
endfb.branchingNG = 0.112
endfb.dilutions = ( 1.e10, 94.5317612, 56.3173141, 33.5510521, 19.9880447, \
11.9078817, 7.09412289, 4.22632504, 2.51783395, 1.5 )
endfb.pendf()
endfb.gendf()
endfb.draglib()
endfb.branchingNG = 0.112
endfb.dilutions = ( 1.e10, 10000.0, 5957.50244, 3549.18335, 2114.42676, \
1259.67004, 750.448669, 447.079956, 266.347961, 158.676849 )
endfb.pendf()
endfb.gendf()
endfb.draglib()
endfb.branchingNG = None

endfb.hmat = "Am242"
endfb.mat = 9546
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.fission = 1 # fission without delayed neutrons
endfb.ss = (22.53556, 1.227732e5)
endfb.potential = 10.2000
endfb.eFiss = 1.9E2
endfb.dilutions = ( 1.e10, 10000.0, 3549.18335, 1259.67004, 447.079956, \
158.676849, 56.3173141, 19.9880447, 7.09412289, 2.51783395 )
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Am242m"
endfb.mat = 9547
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.fission = 1 # fission without delayed neutrons
endfb.ss = (22.53556, 1.227732e5)
endfb.potential = 10.2000
endfb.eFiss = 1.9E2
endfb.dilutions = ( 1.e10, 10000.0, 3549.18335, 1259.67004, 447.079956, \
158.676849, 56.3173141 )
endfb.purr = 1
endfb.pendf()
endfb.gendf()
endfb.draglib()
endfb.purr = None

endfb.hmat = "Am243"
endfb.mat = 9549
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.fission = 1 # fission without delayed neutrons
endfb.ss = (22.53556, 1.227732e5)
endfb.potential = 11.8237
endfb.dilutions = ( 1.e10, 10000.0, 3549.18335, 1259.67004, 447.079956, \
158.676849, 56.3173141, 19.9880447, 7.09412289, 2.51783395 )
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Cm241"
endfb.mat = 9628
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.fission = 1 # fission without delayed neutrons
endfb.ss = (22.53556, 1.227732e5)
endfb.potential = 10.1788
endfb.eFiss = 1.9245E+02
endfb.dilutions = ( 1.e10, 10000.0, 3549.18335, 1259.67004, 447.079956, \
158.676849, 56.3173141 )
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Cm242"
endfb.mat = 9631
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.fission = 1 # fission without delayed neutrons
endfb.ss = (22.53556, 1.227732e5)
endfb.potential = 10.2000
endfb.eFiss = 1.9142E+02
endfb.dilutions = ( 1.e10, 10000.0, 3549.18335, 1259.67004, 447.079956, \
158.676849, 56.3173141 )
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Cm243"
endfb.mat = 9634
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.fission = 1 # fission without delayed neutrons
endfb.ss = (22.53556, 1.227732e5)
endfb.potential = 11.2832
endfb.eFiss = 1.9E2
endfb.dilutions = ( 1.e10, 10000.0, 3549.18335, 1259.67004, 447.079956, \
158.676849, 56.3173141 )
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Cm244"
endfb.mat = 9637
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.fission = 1 # fission without delayed neutrons
endfb.ss = (22.53556, 1.227732e5)
endfb.potential = 10.3200
endfb.eFiss = 1.9049E+02
endfb.dilutions = ( 1.e10, 10000.0, 3549.18335, 1259.67004, 447.079956, \
158.676849, 56.3173141 )
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Cm245"
endfb.mat = 9640
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.fission = 1 # fission without delayed neutrons
endfb.ss = (22.53556, 1.227732e5)
endfb.potential = 11.3900
endfb.eFiss = 1.9E2
endfb.dilutions = ( 1.e10, 10000.0, 3549.18335, 1259.67004, 447.079956, \
158.676849, 56.3173141 )
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Cm246"
endfb.mat = 9643
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.fission = 1 # fission without delayed neutrons
endfb.ss = (22.53556, 1.227732e5)
endfb.potential = 10.2758
endfb.eFiss = 1.9E2
endfb.dilutions = ( 1.e10, 10000.0, 3549.18335, 1259.67004, 447.079956, \
158.676849, 56.3173141 )
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Cm247"
endfb.mat = 9646
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.fission = 1 # fission without delayed neutrons
endfb.ss = (22.53556, 1.227732e5)
endfb.potential = 10.0000
endfb.eFiss = 1.9E2
endfb.dilutions = ( 1.e10, 10000.0, 3549.18335, 1259.67004, 447.079956, \
158.676849, 56.3173141 )
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Cm248"
endfb.mat = 9649
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.fission = 1 # fission without delayed neutrons
endfb.ss = (22.53556, 1.227732e5)
endfb.potential = 10.3971
endfb.eFiss = 1.8885E+02
endfb.dilutions = ( 1.e10, 10000.0, 3549.18335, 1259.67004, 447.079956, \
158.676849, 56.3173141 )
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Bk249"
endfb.mat = 9752
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.fission = 1 # fission without delayed neutrons
endfb.ss = (22.53556, 1.227732e5)
endfb.potential = 11.0983
endfb.eFiss = 1.9E2
endfb.dilutions = ( 1.e10, 10000.0, 3549.18335, 1259.67004, 447.079956, \
158.676849, 56.3173141 )
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Cf249"
endfb.mat = 9852
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.fission = 1 # fission without delayed neutrons
endfb.ss = (22.53556, 1.227732e5)
endfb.potential = 11.1510
endfb.eFiss = 1.9E2
endfb.dilutions = ( 1.e10, 10000.0, 3549.18335, 1259.67004, 447.079956, \
158.676849, 56.3173141 )
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Cf250"
endfb.mat = 9855
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.fission = 1 # fission without delayed neutrons
endfb.ss = (22.53556, 1.227732e5)
endfb.potential = 9.8800
endfb.eFiss = 1.9E2
endfb.dilutions = ( 1.e10, 10000.0, 3549.18335, 1259.67004, 447.079956, \
158.676849, 56.3173141 )
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Cf251"
endfb.mat = 9858
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.fission = 1 # fission without delayed neutrons
endfb.ss = (22.53556, 1.227732e5)
endfb.potential = 9.8800
endfb.eFiss = 1.9E2
endfb.dilutions = ( 1.e10, 10000.0, 3549.18335, 1259.67004, 447.079956, \
158.676849, 56.3173141 )
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Cf252"
endfb.mat = 9861
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.fission = 1 # fission without delayed neutrons
endfb.ss = (22.53556, 1.227732e5)
endfb.potential = 9.8000
endfb.eFiss = 1.9E2
endfb.dilutions = ( 1.e10, 10000.0, 3549.18335, 1259.67004, 447.079956, \
158.676849, 56.3173141 )
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Cf253"
endfb.mat = 9864
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.fission = 1 # fission without delayed neutrons
endfb.ss = (22.53556, 1.227732e5)
endfb.potential = 9.7600
endfb.eFiss = 1.9E2
endfb.dilutions = ( 1.e10, 10000.0, 3549.18335, 1259.67004, 447.079956, \
158.676849, 56.3173141 )
endfb.pendf()
endfb.gendf()
endfb.draglib()

# Process the fission products:

endfb.scatteringLaw = None
endfb.legendre = 0
endfb.fission = None
endfb.dilutions = None
endfb.eFiss = None

endfb.hmat = "Ge72"
endfb.mat = 3231
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Ge73"
endfb.mat = 3234
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Ge74"
endfb.mat = 3237
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Ge76"
endfb.mat = 3243
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "As75"
endfb.mat = 3325
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Se76"
endfb.mat = 3431
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Se77"
endfb.mat = 3434
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Se78"
endfb.mat = 3437
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

# CAUTION: Use Se79 from Jeff3.0
endfb.hmat = "Se79"
endfb.mat = 3440
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/JEFF30N3440_1.ASC"
endfb.makeFp()

endfb.hmat = "Se80"
endfb.mat = 3443
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Se82"
endfb.mat = 3449
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Br79"
endfb.mat = 3525
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Br81"
endfb.mat = 3531
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Kr80"
endfb.mat = 3631
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Kr82"
endfb.mat = 3637
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Kr83"
endfb.mat = 3640
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Kr84"
endfb.mat = 3643
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Kr85"
endfb.mat = 3646
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Kr86"
endfb.mat = 3649
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Rb85"
endfb.mat = 3725
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Rb87"
endfb.mat = 3731
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Sr86"
endfb.mat = 3831
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Sr87"
endfb.mat = 3834
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Sr88"
endfb.mat = 3837
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Sr89"
endfb.mat = 3840
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Sr90"
endfb.mat = 3843
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Y89"
endfb.mat = 3925
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Y90"
endfb.mat = 3928
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Y91"
endfb.mat = 3931
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Zr90"
endfb.mat = 4025
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Zr91"
endfb.mat = 4028
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Zr92"
endfb.mat = 4031
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Zr93"
endfb.mat = 4034
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Zr94"
endfb.mat = 4037
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Zr95"
endfb.mat = 4040
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Zr96"
endfb.mat = 4043
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Nb93"
endfb.mat = 4125
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Nb94"
endfb.mat = 4128
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Nb95"
endfb.mat = 4131
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Mo0"
endfb.mat = 4200
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.ss = (22.53556, 1.227732e5)
endfb.potential = 6.1140
endfb.dilutions = ( 1.e10, 10000.0, 3549.18335, 1259.67004, 447.079956, \
158.676849, 56.3173141, 19.9880447, 7.09412289, 2.51783395 )
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Mo95"
endfb.mat = 4234
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Mo96"
endfb.mat = 4237
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Mo97"
endfb.mat = 4240
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Mo99"
endfb.mat = 4246
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Tc99"
endfb.mat = 4325
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Ru99"
endfb.mat = 4434
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Ru100"
endfb.mat = 4437
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Ru101"
endfb.mat = 4440
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Ru102"
endfb.mat = 4443
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Ru103"
endfb.mat = 4446
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Ru104"
endfb.mat = 4449
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Ru105"
endfb.mat = 4452
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Ru106"
endfb.mat = 4455
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Rh103"
endfb.mat = 4525
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Rh105"
endfb.mat = 4531
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Pd104"
endfb.mat = 4631
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Pd105"
endfb.mat = 4634
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Pd106"
endfb.mat = 4637
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Pd107"
endfb.mat = 4640
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Pd108"
endfb.mat = 4643
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Pd110"
endfb.mat = 4649
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Ag107"
endfb.mat = 4725
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.ss = (22.53556, 1.858471e4)
endfb.potential = 5.4739
endfb.dilutions = ( 1.e10, 10000.0, 4216.96552, 1778.27959, 749.894278, \
316.227791, 133.352152,  56.2341357, 23.7137381, 10.0 )
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Ag109"
endfb.mat = 4731
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.ss = (22.53556, 1.858471e4)
endfb.potential = 5.4739
endfb.dilutions = ( 1.e10, 10000.0, 4216.96552, 1778.27959, 749.894278, \
316.227791, 133.352152,  56.2341357, 23.7137381, 10.0 )
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Ag111"
endfb.mat = 4737
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Cd110"
endfb.mat = 4837
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Cd111"
endfb.mat = 4840
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Cd112"
endfb.mat = 4843
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Cd113"
endfb.mat = 4846
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Cd114"
endfb.mat = 4849
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Cd115m"
endfb.mat = 4853
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Cd116"
endfb.mat = 4855
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "In113"
endfb.mat = 4925
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "In0"
endfb.mat = 4900
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.ss = (22.53556, 1.858471e4)
endfb.potential = 5.0439
endfb.dilutions = ( 1.e10, 10000.0, 4216.96552, 1778.27959, 749.894278, \
316.227791, 133.352152,  56.2341357, 23.7137381, 10.0 )
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Sn115"
endfb.mat = 5034
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Sn116"
endfb.mat = 5037
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Sn117"
endfb.mat = 5040
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Sn118"
endfb.mat = 5043
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Sn119"
endfb.mat = 5046
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Sn120"
endfb.mat = 5049
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Sn122"
endfb.mat = 5055
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Sn123"
endfb.mat = 5058
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Sn124"
endfb.mat = 5061
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Sn125"
endfb.mat = 5064
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Sn126"
endfb.mat = 5067
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Sb121"
endfb.mat = 5125
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Sb123"
endfb.mat = 5131
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Sb124"
endfb.mat = 5134
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Sb125"
endfb.mat = 5137
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Sb126"
endfb.mat = 5140
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Te122"
endfb.mat = 5231
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Te123"
endfb.mat = 5234
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Te124"
endfb.mat = 5237
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Te125"
endfb.mat = 5240
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Te126"
endfb.mat = 5243
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Te127m"
endfb.mat = 5247
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Te128"
endfb.mat = 5249
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Te129m"
endfb.mat = 5253
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Te130"
endfb.mat = 5255
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Te132"
endfb.mat = 5261
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "I127"
endfb.mat = 5325
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "I129"
endfb.mat = 5331
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "I130"
endfb.mat = 5334
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "I131"
endfb.mat = 5337
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "I135"
endfb.mat = 5349
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Xe128"
endfb.mat = 5437
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Xe129"
endfb.mat = 5440
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Xe130"
endfb.mat = 5443
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Xe131"
endfb.mat = 5446
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Xe132"
endfb.mat = 5449
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Xe133"
endfb.mat = 5452
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Xe134"
endfb.mat = 5455
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Xe135"
endfb.mat = 5458
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Xe136"
endfb.mat = 5461
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Cs133"
endfb.mat = 5525
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Cs134"
endfb.mat = 5528
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Cs135"
endfb.mat = 5531
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Cs136"
endfb.mat = 5534
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Cs137"
endfb.mat = 5537
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Ba134"
endfb.mat = 5637
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Ba135"
endfb.mat = 5640
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Ba136"
endfb.mat = 5643
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Ba137"
endfb.mat = 5646
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Ba138"
endfb.mat = 5649
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Ba140"
endfb.mat = 5655
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

# CAUTION: Use La138 from Jeff3.0
endfb.hmat = "La138"
endfb.mat = 5725
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/JEFF30N5725_1.ASC"
endfb.makeFp()

endfb.hmat = "La139"
endfb.mat = 5728
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "La140"
endfb.mat = 5731
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Ce140"
endfb.mat = 5837
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Ce141"
endfb.mat = 5840
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Ce142"
endfb.mat = 5843
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Ce143"
endfb.mat = 5846
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Ce144"
endfb.mat = 5849
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Pr141"
endfb.mat = 5925
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Pr142"
endfb.mat = 5928
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Pr143"
endfb.mat = 5931
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Nd142"
endfb.mat = 6025
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Nd143"
endfb.mat = 6028
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Nd144"
endfb.mat = 6031
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Nd145"
endfb.mat = 6034
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Nd146"
endfb.mat = 6037
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Nd147"
endfb.mat = 6040
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Nd148"
endfb.mat = 6043
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Nd150"
endfb.mat = 6049
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Pm147"
endfb.mat = 6149
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.branchingNG = 0.47
endfb.makeFp()
endfb.branchingNG = None

endfb.hmat = "Pm148"
endfb.mat = 6152
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Pm148m"
endfb.mat = 6153
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Pm149"
endfb.mat = 6155
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Pm151"
endfb.mat = 6161
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Sm147"
endfb.mat = 6234
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Sm148"
endfb.mat = 6237
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Sm149"
endfb.mat = 6240
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Sm150"
endfb.mat = 6243
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Sm151"
endfb.mat = 6246
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Sm152"
endfb.mat = 6249
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Sm153"
endfb.mat = 6252
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Sm154"
endfb.mat = 6255
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Eu151"
endfb.mat = 6325
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Eu152"
endfb.mat = 6328
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Eu153"
endfb.mat = 6331
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Eu154"
endfb.mat = 6334
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Eu155"
endfb.mat = 6337
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Eu156"
endfb.mat = 6340
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Eu157"
endfb.mat = 6343
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Gd154"
endfb.mat = 6431
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.ss = (22.53556, 1.858471e4)
endfb.potential = 7.6587
endfb.dilutions = ( 1.e10, 10000.0, 4216.96552, 1778.27959, 749.894278, \
316.227791, 133.352152,  56.2341357, 23.7137381, 10.0 )
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Gd155"
endfb.mat = 6434
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.ss = (22.53556, 1.858471e4)
endfb.potential = 8.0425
endfb.dilutions = ( 1.e10, 10000.0, 4216.96552, 1778.27959, 749.894278, \
316.227791, 133.352152,  56.2341357, 23.7137381, 10.0 )
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Gd156"
endfb.mat = 6437
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.ss = (22.53556, 1.858471e4)
endfb.potential = 8.2448
endfb.dilutions = ( 1.e10, 10000.0, 4216.96552, 1778.27959, 749.894278, \
316.227791, 133.352152,  56.2341357, 23.7137381, 10.0 )
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Gd157"
endfb.mat = 6440
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.ss = (22.53556, 1.858471e4)
endfb.potential = 7.6454
endfb.dilutions = ( 1.e10, 10000.0, 4216.96552, 1778.27959, 749.894278, \
316.227791, 133.352152,  56.2341357, 23.7137381, 10.0 )
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Gd158"
endfb.mat = 6443
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.ss = (22.53556, 1.858471e4)
endfb.potential = 5.3093
endfb.dilutions = ( 1.e10, 10000.0, 4216.96552, 1778.27959, 749.894278, \
316.227791, 133.352152,  56.2341357, 23.7137381, 10.0 )
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Gd160"
endfb.mat = 6449
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.ss = (22.53556, 1.858471e4)
endfb.potential = 5.8107
endfb.dilutions = ( 1.e10, 10000.0, 4216.96552, 1778.27959, 749.894278, \
316.227791, 133.352152,  56.2341357, 23.7137381, 10.0 )
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.dilutions = None

endfb.hmat = "Tb159"
endfb.mat = 6525
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Tb160"
endfb.mat = 6528
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Dy160"
endfb.mat = 6637
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.ss = (22.53556, 1.858471e4)
endfb.potential = 7.0686
endfb.dilutions = ( 1.e10, 10000.0, 4216.96552, 1778.27959, 749.894278, \
316.227791, 133.352152,  56.2341357, 23.7137381, 10.0 )
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Dy161"
endfb.mat = 6640
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.ss = (22.53556, 1.858471e4)
endfb.potential = 7.0686
endfb.dilutions = ( 1.e10, 10000.0, 4216.96552, 1778.27959, 749.894278, \
316.227791, 133.352152,  56.2341357, 23.7137381, 10.0 )
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Dy162"
endfb.mat = 6643
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.ss = (22.53556, 1.858471e4)
endfb.potential = 7.6846
endfb.dilutions = ( 1.e10, 10000.0, 4216.96552, 1778.27959, 749.894278, \
316.227791, 133.352152,  56.2341357, 23.7137381, 10.0 )
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Dy164"
endfb.mat = 6649
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.ss = (22.53556, 1.858471e4)
endfb.potential = 7.6846
endfb.dilutions = ( 1.e10, 10000.0, 4216.96552, 1778.27959, 749.894278, \
316.227791, 133.352152,  56.2341357, 23.7137381, 10.0 )
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Hf174"
endfb.mat = 7225
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.ss = (22.53556, 1.858471e4)
endfb.potential = 7.1385
endfb.dilutions = ( 1.e10, 10000.0, 4216.96552, 1778.27959, 749.894278, \
316.227791, 133.352152,  56.2341357, 23.7137381, 10.0 )
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Hf176"
endfb.mat = 7231
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.ss = (22.53556, 1.858471e4)
endfb.potential = 7.1935
endfb.dilutions = ( 1.e10, 10000.0, 4216.96552, 1778.27959, 749.894278, \
316.227791, 133.352152,  56.2341357, 23.7137381, 10.0 )
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Hf177"
endfb.mat = 7234
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.ss = (22.53556, 1.858471e4)
endfb.potential = 7.2202
endfb.dilutions = ( 1.e10, 10000.0, 4216.96552, 1778.27959, 749.894278, \
316.227791, 133.352152,  56.2341357, 23.7137381, 10.0 )
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Hf178"
endfb.mat = 7237
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.ss = (22.53556, 1.858471e4)
endfb.potential = 7.2469
endfb.dilutions = ( 1.e10, 10000.0, 4216.96552, 1778.27959, 749.894278, \
316.227791, 133.352152,  56.2341357, 23.7137381, 10.0 )
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Hf179"
endfb.mat = 7240
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.ss = (22.53556, 1.858471e4)
endfb.potential = 7.2736
endfb.dilutions = ( 1.e10, 10000.0, 4216.96552, 1778.27959, 749.894278, \
316.227791, 133.352152,  56.2341357, 23.7137381, 10.0 )
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Hf180"
endfb.mat = 7243
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.ss = (22.53556, 1.858471e4)
endfb.potential = 7.3004
endfb.dilutions = ( 1.e10, 10000.0, 4216.96552, 1778.27959, 749.894278, \
316.227791, 133.352152,  56.2341357, 23.7137381, 10.0 )
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.dilutions = None

endfb.hmat = "Ho165"
endfb.mat = 6725
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Er166"
endfb.mat = 6837
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

endfb.hmat = "Er167"
endfb.mat = 6840
endfb.evaluationFile = "$HOME/evaluations/ENDFB6r8/%d.OUT"%endfb.mat
endfb.makeFp()

# Process the burnup chain:

endfb.fissionFile = "$HOME/evaluations/ENDFB6r8/TAPE.107"
endfb.decayFile = "$HOME/evaluations/ENDFB6r8/TAPE.106"
endfb.burnup()
