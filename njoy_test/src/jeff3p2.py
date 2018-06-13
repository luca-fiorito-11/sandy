#!/usr/local/bin/python
from PyNjoy import *
from os import uname
jeff3p2 = PyNjoy()
jeff3p2.evaluationName = "/tmp/Jeff3.2"
jeff3p2.execDir = "../" + uname()[0]
jeff3p2.nstr = 22
jeff3p2.iwt = 4
jeff3p2.Espectra = None
jeff3p2.autolib = (2.76792, 677.2873, 0.00125)

jeff3p2.legendre = 3
jeff3p2.hmat = "H1_H2O"
jeff3p2.mat = 125
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-1-H-001.jeff32"
jeff3p2.scatteringLaw =  "$HOME/evaluations/Jeff3.2/JEFF31TS_INDIV/JEFF31TS0001_1.ASC"
jeff3p2.scatteringMat = 1
jeff3p2.temperatures = ( 293.6, 323.6, 373.6, 423.6, 473.6, 523.6, 573.6, 647.2, 800.0 )
jeff3p2.fission = None
jeff3p2.dilutions = None
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.legendre = 1
jeff3p2.hmat = "H2_D2O"
jeff3p2.mat = 128
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-1-H-002.jeff32"
jeff3p2.scatteringLaw =  "$HOME/evaluations/Jeff3.2/JEFF31TS_INDIV/JEFF31TS0011_1.ASC"
jeff3p2.scatteringMat = 11
jeff3p2.temperatures = ( 293.6, 323.6, 373.6, 423.6, 473.6, 523.6, 573.6, 643.9 )
jeff3p2.fission = None
jeff3p2.dilutions = None
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "H1_CH2"
jeff3p2.mat = 125
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-1-H-001.jeff32"
jeff3p2.scatteringLaw =  "$HOME/evaluations/Jeff3.2/JEFF31TS_INDIV/JEFF31TS0037_1.ASC"
jeff3p2.scatteringMat = 37
jeff3p2.temperatures = ( 293.6,  350. )
jeff3p2.fission = None
jeff3p2.dilutions = None
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "H1_ZRH"
jeff3p2.mat = 125
jeff3p2.temperatures = ( 293.6, 400., 500., 600., 700., 800., 1000., 1200. )
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-1-H-001.jeff32"
jeff3p2.scatteringLaw =  "$HOME/evaluations/Jeff3.2/JEFF31TS_INDIV/JEFF31TS0007_1.ASC"
jeff3p2.scatteringMat = 7
jeff3p2.fission = None
jeff3p2.dilutions = None
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "Be9"
jeff3p2.mat = 425
jeff3p2.temperatures = ( 293.6 , 400 , 500 , 600 , 700 , 800 , 1000 , 1200)
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-4-Be-009.jeff32"
jeff3p2.scatteringLaw =  "$HOME/evaluations/Jeff3.2/JEFF31TS_INDIV/JEFF31TS0026_1.ASC"
jeff3p2.scatteringMat = 26
jeff3p2.fission = None
jeff3p2.dilutions = None
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "C0_GR"
jeff3p2.mat = 600
jeff3p2.temperatures = ( 293.6 , 400 , 500 , 600 , 700 , 800 , 1000 , 1200 , 1600 , 2000 )
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-6-C-000.jeff32"
jeff3p2.scatteringLaw =  "$HOME/evaluations/Jeff3.2/JEFF31TS_INDIV/JEFF31TS0031_1.ASC"
jeff3p2.scatteringMat = 31
jeff3p2.fission = None
jeff3p2.dilutions = None
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.scatteringLaw = None
jeff3p2.temperatures = ( 293., 550., 900., 1200., 2000. )
jeff3p2.fission = None
jeff3p2.dilutions = None

jeff3p2.hmat = "H1"
jeff3p2.mat = 125
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-1-H-001.jeff32"
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "H2"
jeff3p2.mat = 128
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-1-H-002.jeff32"
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "H3"
jeff3p2.mat = 131
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-1-H-003.jeff32"
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "He3"
jeff3p2.mat = 225
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-2-He-003.jeff32"
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "He4"
jeff3p2.mat = 228
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-2-He-004.jeff32"
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "Li6"
jeff3p2.mat = 325
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-3-Li-006.jeff32"
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "Li7"
jeff3p2.mat = 328
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-3-Li-007.jeff32"
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "B10"
jeff3p2.mat = 525
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-5-B-010.jeff32"
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "B11"
jeff3p2.mat = 528
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-5-B-011.jeff32"
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "C0"
jeff3p2.mat = 600
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-6-C-000.jeff32"
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "N14"
jeff3p2.mat = 725
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-7-N-014.jeff32"
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "N15"
jeff3p2.mat = 728
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-7-N-015.jeff32"
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "O16"
jeff3p2.mat = 825
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-8-O-016.jeff32"
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "O17"
jeff3p2.mat = 828
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-8-O-017.jeff32"
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "F19"
jeff3p2.mat = 925
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-9-F-019.jeff32"
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "Na23"
jeff3p2.mat = 1125
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-11-Na-023.jeff32"
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "Mg24"
jeff3p2.mat = 1225
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-12-Mg-024.jeff32"
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "Mg25"
jeff3p2.mat = 1228
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-12-Mg-025.jeff32"
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "Mg26"
jeff3p2.mat = 1231
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-12-Mg-026.jeff32"
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "Al27"
jeff3p2.mat = 1325
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-13-Al-027.jeff32"
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "Si28"
jeff3p2.mat = 1425
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-14-Si-028.jeff32"
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "Si29"
jeff3p2.mat = 1428
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-14-Si-029.jeff32"
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "Si30"
jeff3p2.mat = 1431
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-14-Si-030.jeff32"
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "P31"
jeff3p2.mat = 1525
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-15-P-031.jeff32"
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "S32"
jeff3p2.mat = 1625
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-16-S-032.jeff32"
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "S33"
jeff3p2.mat = 1628
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-16-S-033.jeff32"
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "S34"
jeff3p2.mat = 1631
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-16-S-034.jeff32"
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "S36"
jeff3p2.mat = 1637
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-16-S-036.jeff32"
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "Cl35"
jeff3p2.mat = 1725
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-17-Cl-035.jeff32"
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "Cl37"
jeff3p2.mat = 1731
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-17-Cl-037.jeff32"
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "K39"
jeff3p2.mat = 1925
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-19-K-039.jeff32"
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "K40"
jeff3p2.mat = 1928
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-19-K-040.jeff32"
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "K41"
jeff3p2.mat = 1931
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-19-K-041.jeff32"
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "Ca40"
jeff3p2.mat = 2025
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-20-Ca-040.jeff32"
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "Ca42"
jeff3p2.mat = 2031
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-20-Ca-042.jeff32"
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "Ca43"
jeff3p2.mat = 2034
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-20-Ca-043.jeff32"
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "Ca44"
jeff3p2.mat = 2037
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-20-Ca-044.jeff32"
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "Ca46"
jeff3p2.mat = 2043
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-20-Ca-046.jeff32"
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "Ca48"
jeff3p2.mat = 2049
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-20-Ca-048.jeff32"
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "Ti46"
jeff3p2.mat = 2225
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-22-Ti-046.jeff32"
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "Ti47"
jeff3p2.mat = 2228
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-22-Ti-047.jeff32"
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "Ti48"
jeff3p2.mat = 2231
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-22-Ti-048.jeff32"
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "Ti49"
jeff3p2.mat = 2234
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-22-Ti-049.jeff32"
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "Ti50"
jeff3p2.mat = 2237
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-22-Ti-050.jeff32"
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "V00"
jeff3p2.mat = 2300
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-23-V-000.jeff32"
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "Cr50"
jeff3p2.mat = 2425
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-24-Cr-050.jeff32"
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "Cr51"
jeff3p2.mat = 2428
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-24-Cr-051.jeff32"
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "Cr52"
jeff3p2.mat = 2431
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-24-Cr-052.jeff32"
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "Cr53"
jeff3p2.mat = 2434
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-24-Cr-053.jeff32"
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "Cr54"
jeff3p2.mat = 2437
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-24-Cr-054.jeff32"
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "Mn55"
jeff3p2.mat = 2525
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-25-Mn-055.jeff32"
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "Fe54"
jeff3p2.mat = 2625
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-26-Fe-054.jeff32"
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "Fe56"
jeff3p2.mat = 2631
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-26-Fe-056.jeff32"
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "Fe57"
jeff3p2.mat = 2634
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-26-Fe-057.jeff32"
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "Fe58"
jeff3p2.mat = 2637
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-26-Fe-058.jeff32"
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "Co59"
jeff3p2.mat = 2725
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-27-Co-059.jeff32"
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "Ni58"
jeff3p2.mat = 2825
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-28-Ni-058.jeff32"
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "Ni60"
jeff3p2.mat = 2831
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-28-Ni-060.jeff32"
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "Ni61"
jeff3p2.mat = 2834
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-28-Ni-061.jeff32"
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "Ni62"
jeff3p2.mat = 2837
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-28-Ni-062.jeff32"
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "Ni64"
jeff3p2.mat = 2843
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-28-Ni-064.jeff32"
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "Cu63"
jeff3p2.mat = 2925
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-29-Cu-063.jeff32"
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "Cu65"
jeff3p2.mat = 2931
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-29-Cu-065.jeff32"
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "Zn66"
jeff3p2.mat = 3031
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-30-Zn-066.jeff32"
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "Zn67"
jeff3p2.mat = 3034
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-30-Zn-067.jeff32"
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "Zn68"
jeff3p2.mat = 3037
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-30-Zn-068.jeff32"
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "Zn70"
jeff3p2.mat = 3043
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-30-Zn-070.jeff32"
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "Ga69"
jeff3p2.mat = 3125
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-31-Ga-069.jeff32"
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "Ga71"
jeff3p2.mat = 3131
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-31-Ga-071.jeff32"
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "W182"
jeff3p2.mat =  7431  
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-74-W-182.jeff32"
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "W183"
jeff3p2.mat = 7434 
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-74-W-183.jeff32"
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "W184"
jeff3p2.mat = 7437
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-74-W-184.jeff32"
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "W186"
jeff3p2.mat =  7443
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-74-W-186.jeff32"
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "Re187"
jeff3p2.mat =  7531
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-75-Re-187.jeff32"
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "Au197"
jeff3p2.mat =  7925
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-79-Au-197.jeff32"
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "Pb204"
jeff3p2.mat =  8225
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-82-Pb-204.jeff32"
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "Pb206"
jeff3p2.mat =  8231 
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-82-Pb-206.jeff32"
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "Pb207"
jeff3p2.mat =  8234
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-82-Pb-207.jeff32"
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "Pb208"
jeff3p2.mat =  8237
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-82-Pb-208.jeff32"
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "Bi209"
jeff3p2.mat =  8325
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-83-Bi-209.jeff32"
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "Th230"
jeff3p2.mat = 9034
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-90-Th-230.jeff32"
jeff3p2.fission = 1 # fission without delayed neutrons
jeff3p2.ss = (2.76792, 1.22773e5)
jeff3p2.potential = 8.7040
jeff3p2.eFiss = 1.7406E+02
jeff3p2.dilutions = ( 1.e10, 10000.0, 3546.31, 1257.43, 445.8898, 158.1139, 56.0677, 19.8818, 7.0501, 2.5 )
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

#Use Jeff3.1.2 Th232 to have 8 delayed neutron groups
jeff3p2.hmat = "Th232"
jeff3p2.mat = 9040
#jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-90-Th-232.jeff32"
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.1.2/JEFF312N/JEFF312N9040_0.ASC"
jeff3p2.fission = 2 # fission with delayed neutrons
jeff3p2.ss = (2.76792, 1.22773e5)
jeff3p2.potential = 11.8699
jeff3p2.eFiss = 1.7193E+02
jeff3p2.dilutions = ( 1.e10, 94.5317612, 56.3173141, 33.5510521, 19.9880447, \
11.9078817, 7.09412289, 4.22632504, 2.51783395, 1.5 )
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()
jeff3p2.eFiss = 1.7193E+02
jeff3p2.dilutions = ( 1.e10, 10000.0, 5957.50244, 3549.18335, 2114.42676, \
1259.67004, 750.448669, 447.079956, 266.347961, 158.676849 )
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "Th233"
jeff3p2.mat = 9043
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-90-Th-233.jeff32"
jeff3p2.fission = 1 # fission without delayed neutrons
jeff3p2.ss = (22.53556, 1.227732e5)
jeff3p2.potential = 11.8699
jeff3p2.eFiss = 1.7193E+02
jeff3p2.dilutions = ( 1.e10, 10000.0, 3546.31, 1257.43, 445.8898, 158.1139, 56.0677, 19.8818, 7.0501, 2.5 )
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "Pa231"
jeff3p2.mat = 9131
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-91-Pa-231.jeff32"
jeff3p2.fission = 1 # fission without delayed neutrons
jeff3p2.ss = (2.76792, 1.22773e5)
jeff3p2.potential = 8.7266
jeff3p2.eFiss = 1.7E2
jeff3p2.dilutions = ( 1.e10, 10000.0, 3546.31, 1257.43, 445.8898, 158.1139, 56.0677, 19.8818, 7.0501, 2.5 )
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "Pa233"
jeff3p2.mat = 9137
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-91-Pa-233.jeff32"
jeff3p2.fission = 1 # fission without delayed neutrons
jeff3p2.ss = (2.76792, 1.22773e5)
jeff3p2.potential = 9.9950
jeff3p2.eFiss = 1.7558E+02
jeff3p2.dilutions = ( 1.e10, 10000.0, 3546.31, 1257.43, 445.8898, 158.1139, 56.0677, 19.8818, 7.0501, 2.5 )
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "U232"
jeff3p2.mat = 9219
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-92-U-232.jeff32"
jeff3p2.fission = 1 # fission without delayed neutrons
jeff3p2.ss = (2.76792, 1.22773e5)
jeff3p2.potential = 12.0687
jeff3p2.eFiss = 1.8E2
jeff3p2.dilutions = ( 1.e10, 10000.0, 3546.31, 1257.43, 445.8898, 158.1139, 56.0677, 19.8818, 7.0501, 2.5 )
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "U233"
jeff3p2.mat = 9222
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-92-U-233.jeff32"
jeff3p2.fission = 2 # fission with delayed neutrons
jeff3p2.ss = (2.76792, 1.22773e5)
jeff3p2.potential = 12.2989
jeff3p2.eFiss = 1.8087E+02
jeff3p2.dilutions = ( 1.e10, 10000.0, 3546.31, 1257.43, 445.8898, 158.1139, 56.0677, 19.8818, 7.0501, 2.5 )
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "U234"
jeff3p2.mat = 9225
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-92-U-234.jeff32"
jeff3p2.fission = 1 # fission without delayed neutrons
jeff3p2.ss = (2.76792, 1.22773e5)
jeff3p2.potential = 10.0210
jeff3p2.dilutions = ( 1.e10, 10000.0, 3546.31, 1257.43, 445.8898, 158.1139, 56.0677, 19.8818, 7.0501, 2.5 )
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "U235"
jeff3p2.mat = 9228
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-92-U-235.jeff32"
jeff3p2.fission = 2 # fission with delayed neutrons
jeff3p2.ss = (2.76792, 1.22773e5)
jeff3p2.potential = 11.6070
jeff3p2.dilutions = ( 1.e10, 94.5317612, 56.3173141, 33.5510521, 19.9880447, \
11.9078817, 7.09412289, 4.22632504, 2.51783395, 1.5 )
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()
jeff3p2.dilutions = ( 1.e10, 10000.0, 5957.50244, 3549.18335, 2114.42676, \
1259.67004, 750.448669, 447.079956, 266.347961, 158.676849 )
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "U236"
jeff3p2.mat = 9231
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-92-U-236.jeff32"
jeff3p2.fission = 1 # fission without delayed neutrons
jeff3p2.ss = (2.76792, 1.22773e5)
jeff3p2.potential = 10.9954
jeff3p2.dilutions = ( 1.e10, 10000.0, 3546.31, 1257.43, 445.8898, 158.1139, 56.0677, 19.8818, 7.0501, 2.5 )
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "U237"
jeff3p2.mat = 9234
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-92-U-237.jeff32"
jeff3p2.fission = 1 # fission without delayed neutrons
jeff3p2.ss = (2.76792, 1.22773e5)
jeff3p2.potential = 10.5000
jeff3p2.dilutions = ( 1.e10, 10000.0, 3546.31, 1257.43, 445.8898, 158.1139, 56.0677, 19.8818, 7.0501, 2.5 )
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "U238"
jeff3p2.mat = 9237
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-92-U-238.jeff32"
jeff3p2.fission = 2 # fission with delayed neutrons
jeff3p2.ss = (2.76792, 1.22773e5)
jeff3p2.potential = 11.17103
jeff3p2.dilutions = ( 1.e10, 94.5317612, 56.3173141, 33.5510521, 19.9880447, \
11.9078817, 7.09412289, 4.22632504, 2.51783395, 1.5 )
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()
jeff3p2.dilutions = ( 1.e10, 10000.0, 5957.50244, 3549.18335, 2114.42676, \
1259.67004, 750.448669, 447.079956, 266.347961, 158.676849 )
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "Np236"
jeff3p2.mat = 9343
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-93-Np-236.jeff32"
jeff3p2.fission = 1 # fission without delayed neutrons
jeff3p2.ss = (2.76792, 1.22773e5)
jeff3p2.eFiss = 1.8E2
jeff3p2.potential = 12.1427
jeff3p2.dilutions = ( 1.e10, 10000.0, 3546.31, 1257.43, 445.8898, 158.1139, 56.0677, 19.8818, 7.0501, 2.5 )
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "Np237"
jeff3p2.mat = 9346
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-93-Np-237.jeff32"
jeff3p2.fission = 1 # fission without delayed neutrons
jeff3p2.ss = (2.76792, 1.22773e5)
jeff3p2.potential = 11.4369
jeff3p2.branchingN2N = 0.801
jeff3p2.dilutions = ( 1.e10, 10000.0, 3546.31, 1257.43, 445.8898, 158.1139, 56.0677, 19.8818, 7.0501, 2.5 )
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()
jeff3p2.branchingN2N = None

jeff3p2.hmat = "Np238"
jeff3p2.mat = 9349
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-93-Np-238.jeff32"
jeff3p2.fission = 1 # fission without delayed neutrons
jeff3p2.ss = (2.76792, 1.22773e5)
jeff3p2.eFiss = 1.8E2
jeff3p2.potential = 10.4000
jeff3p2.dilutions = ( 1.e10, 10000.0, 3549.18335, 1259.67004, 447.079956, \
158.676849, 56.3173141 )
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "Np239"
jeff3p2.mat = 9352
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-93-Np-239.jeff32"
jeff3p2.fission = 1 # fission without delayed neutrons
jeff3p2.ss = (2.76792, 1.22773e5)
jeff3p2.eFiss = 1.8E2
jeff3p2.potential = 10.4979
jeff3p2.dilutions = ( 1.e10, 10000.0, 3546.31, 1257.43, 445.8898, 158.1139, 56.0677, 19.8818, 7.0501, 2.5 )
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "Pu236"
jeff3p2.mat = 9428
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-94-Pu-236.jeff32"
jeff3p2.fission = 1 # fission without delayed neutrons
jeff3p2.ss = (2.76792, 1.22773e5)
jeff3p2.eFiss = 1.8E2
jeff3p2.potential = 11.2458
jeff3p2.dilutions = ( 1.e10, 10000.0, 3546.31, 1257.43, 445.8898, 158.1139, 56.0677, 19.8818, 7.0501, 2.5 )
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "Pu237"
jeff3p2.mat = 9431
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-94-Pu-237.jeff32"
jeff3p2.fission = 1 # fission without delayed neutrons
jeff3p2.ss = (2.76792, 1.22773e5)
jeff3p2.potential = 10.5209
jeff3p2.dilutions = ( 1.e10, 10000.0, 3546.31, 1257.43, 445.8898, 158.1139, 56.0677, 19.8818, 7.0501, 2.5 )
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "Pu238"
jeff3p2.mat = 9434
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-94-Pu-238.jeff32"
jeff3p2.fission = 1 # fission without delayed neutrons
jeff3p2.ss = (2.76792, 1.22773e5)
jeff3p2.potential = 10.8897
jeff3p2.dilutions = ( 1.e10, 94.5317612, 56.3173141, 33.5510521, 19.9880447, \
11.9078817, 7.09412289, 4.22632504, 2.51783395, 1.5 )
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()
jeff3p2.dilutions = ( 1.e10, 10000.0, 5957.50244, 3549.18335, 2114.42676, \
1259.67004, 750.448669, 447.079956, 266.347961, 158.676849 )
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "Pu239"
jeff3p2.mat = 9437
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-94-Pu-239.jeff32"
jeff3p2.fission = 2 # fission with delayed neutrons
jeff3p2.ss = (2.76792, 1.22773e5)
jeff3p2.potential = 10.8897
jeff3p2.dilutions = ( 1.e10, 158.887822, 100.279413, 63.2896882, \
39.9442369, 25.2101426, 15.9109633, 10.0419406, 6.33780423, 4.0 )
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()
jeff3p2.dilutions = ( 1.e10, 10000.0, 6311.33494, 3983.29436, 2513.99016, \
1586.66319, 1001.39615, 632.014573, 398.885514, 251.749976 )
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "Pu240"
jeff3p2.mat = 9440
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-94-Pu-240.jeff32"
jeff3p2.fission = 2 # fission with delayed neutrons
jeff3p2.ss = (2.76792, 1.22773e5)
jeff3p2.potential = 9.9091
jeff3p2.dilutions = ( 1.e10, 94.5317612, 56.3173141, 33.5510521, 19.9880447, \
11.9078817, 7.09412289, 4.22632504, 2.51783395, 1.5 )
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()
jeff3p2.dilutions = ( 1.e10, 10000.0, 5957.50244, 3549.18335, 2114.42676, \
1259.67004, 750.448669, 447.079956, 266.347961, 158.676849 )
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "Pu241"
jeff3p2.mat = 9443
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-94-Pu-241.jeff32"
jeff3p2.fission = 2 # fission with delayed neutrons
jeff3p2.ss = (2.76792, 1.22773e5)
jeff3p2.potential = 11.2156
jeff3p2.dilutions = ( 1.e10, 94.5317612, 56.3173141, 33.5510521, 19.9880447, \
11.9078817, 7.09412289, 4.22632504, 2.51783395, 1.5 )
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()
jeff3p2.dilutions = ( 1.e10, 10000.0, 5957.50244, 3549.18335, 2114.42676, \
1259.67004, 750.448669, 447.079956, 266.347961, 158.676849 )
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "Pu242"
jeff3p2.mat = 9446
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-94-Pu-242.jeff32"
jeff3p2.fission = 1 # fission without delayed neutrons
jeff3p2.ss = (2.76792, 1.22773e5)
jeff3p2.potential = 10.6961
jeff3p2.dilutions = ( 1.e10,  469.546659, 258.807233, 142.650752, \
78.6270025, 43.3380508, 23.8872981, 13.1663284, 7.25708715, 4.0 )
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()
jeff3p2.dilutions = ( 1.e10, 1.e5, 55118.5161, 30380.5178, 16745.2959, \
9229.76155, 5087.30922, 2804.05024, 1545.55137, 851.885253 )
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "Pu243"
jeff3p2.mat = 9449
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-94-Pu-243.jeff32"
jeff3p2.fission = 1 # fission without delayed neutrons
jeff3p2.ss = (2.76792, 1.22773e5)
jeff3p2.potential = 10.2000
jeff3p2.eFiss = 1.9E2
jeff3p2.dilutions = ( 1.e10, 10000.0, 3546.31, 1257.43, 445.8898, 158.1139, 56.0677, 19.8818, 7.0501, 2.5 )
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "Pu244"
jeff3p2.mat = 9452
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-94-Pu-244.jeff32"
jeff3p2.fission = 1 # fission without delayed neutrons
jeff3p2.ss = (2.76792, 1.22773e5)
jeff3p2.potential = 10.1000
jeff3p2.dilutions = ( 1.e10, 10000.0, 3546.31, 1257.43, 445.8898, 158.1139, 56.0677, 19.8818, 7.0501, 2.5 )
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "Am241"
jeff3p2.mat = 9543
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-95-Am-241.jeff32"
jeff3p2.fission = 1 # fission without delayed neutrons
jeff3p2.ss = (2.76792, 1.22773e5)
jeff3p2.potential = 11.0329
jeff3p2.branchingNG = 0.115
jeff3p2.dilutions = ( 1.e10, 94.5317612, 56.3173141, 33.5510521, 19.9880447, \
11.9078817, 7.09412289, 4.22632504, 2.51783395, 1.5 )
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()
jeff3p2.branchingNG = 0.115
jeff3p2.dilutions = ( 1.e10, 10000.0, 5957.50244, 3549.18335, 2114.42676, \
1259.67004, 750.448669, 447.079956, 266.347961, 158.676849 )
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()
jeff3p2.branchingNG = None

jeff3p2.hmat = "Am242"
jeff3p2.mat = 9546
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-95-Am-242.jeff32"
jeff3p2.fission = 0 # no fission matrix!!!
jeff3p2.ss = (2.76792, 1.22773e5)
jeff3p2.potential = 10.2000
jeff3p2.eFiss = 1.9E2
jeff3p2.dilutions = ( 1.e10, 10000.0, 3546.31, 1257.43, 445.8898, 158.1139, 56.0677, 19.8818, 7.0501, 2.5 )
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "Am242m"
jeff3p2.mat = 9547
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-95-Am-242M.jeff32"
jeff3p2.fission = 1 # fission without delayed neutrons
jeff3p2.ss = (2.76792, 1.22773e5)
jeff3p2.potential = 10.2000
jeff3p2.eFiss = 1.9E2
jeff3p2.dilutions = ( 1.e10, 10000.0, 3549.18335, 1259.67004, 447.079956, \
158.676849, 56.3173141 )
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "Am243"
jeff3p2.mat = 9549
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-95-Am-243.jeff32"
jeff3p2.fission = 1 # fission without delayed neutrons
jeff3p2.ss = (2.76792, 1.22773e5)
jeff3p2.potential = 11.8237
jeff3p2.dilutions = ( 1.e10, 10000.0, 3546.31, 1257.43, 445.8898, 158.1139, 56.0677, 19.8818, 7.0501, 2.5 )
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "Cm241"
jeff3p2.mat = 9628
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-96-Cm-241.jeff32"
jeff3p2.fission = 1 # fission without delayed neutrons
jeff3p2.ss = (2.76792, 1.22773e5)
jeff3p2.potential = 10.1788
jeff3p2.eFiss = 1.9245E+02
jeff3p2.dilutions = ( 1.e10, 10000.0, 3549.18335, 1259.67004, 447.079956, \
158.676849, 56.3173141 )
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "Cm242"
jeff3p2.mat = 9631
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-96-Cm-242.jeff32"
jeff3p2.fission = 1 # fission without delayed neutrons
jeff3p2.ss = (2.76792, 1.22773e5)
jeff3p2.potential = 10.2000
jeff3p2.eFiss = 1.9142E+02
jeff3p2.dilutions = ( 1.e10, 10000.0, 3549.18335, 1259.67004, 447.079956, \
158.676849, 56.3173141 )
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "Cm243"
jeff3p2.mat = 9634
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-96-Cm-243.jeff32"
jeff3p2.fission = 1 # fission without delayed neutrons
jeff3p2.ss = (2.76792, 1.22773e5)
jeff3p2.potential = 11.2832
jeff3p2.eFiss = 1.9E2
jeff3p2.dilutions = ( 1.e10, 10000.0, 3549.18335, 1259.67004, 447.079956, \
158.676849, 56.3173141 )
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "Cm244"
jeff3p2.mat = 9637
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-96-Cm-244.jeff32"
jeff3p2.fission = 1 # fission without delayed neutrons
jeff3p2.ss = (2.76792, 1.22773e5)
jeff3p2.potential = 10.3200
jeff3p2.eFiss = 1.9049E+02
jeff3p2.dilutions = ( 1.e10, 10000.0, 3549.18335, 1259.67004, 447.079956, \
158.676849, 56.3173141 )
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "Cm245"
jeff3p2.mat = 9640
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-96-Cm-245.jeff32"
jeff3p2.fission = 1 # fission without delayed neutrons
jeff3p2.ss = (2.76792, 1.22773e5)
jeff3p2.potential = 11.3900
jeff3p2.eFiss = 1.9E2
jeff3p2.dilutions = ( 1.e10, 10000.0, 3549.18335, 1259.67004, 447.079956, \
158.676849, 56.3173141 )
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "Cm246"
jeff3p2.mat = 9643
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-96-Cm-246.jeff32"
jeff3p2.fission = 1 # fission without delayed neutrons
jeff3p2.ss = (2.76792, 1.22773e5)
jeff3p2.potential = 10.2758
jeff3p2.eFiss = 1.9E2
jeff3p2.dilutions = ( 1.e10, 10000.0, 3549.18335, 1259.67004, 447.079956, \
158.676849, 56.3173141 )
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "Cm247"
jeff3p2.mat = 9646
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-96-Cm-247.jeff32"
jeff3p2.fission = 1 # fission without delayed neutrons
jeff3p2.ss = (2.76792, 1.22773e5)
jeff3p2.potential = 10.0000
jeff3p2.eFiss = 1.9E2
jeff3p2.dilutions = ( 1.e10, 10000.0, 3549.18335, 1259.67004, 447.079956, \
158.676849, 56.3173141 )
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "Cm248"
jeff3p2.mat = 9649
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-96-Cm-248.jeff32"
jeff3p2.fission = 1 # fission without delayed neutrons
jeff3p2.ss = (2.76792, 1.22773e5)
jeff3p2.potential = 10.3971
jeff3p2.eFiss = 1.8885E+02
jeff3p2.dilutions = ( 1.e10, 10000.0, 3549.18335, 1259.67004, 447.079956, \
158.676849, 56.3173141 )
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "Bk249"
jeff3p2.mat = 9752
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-97-Bk-249.jeff32"
jeff3p2.fission = 1 # fission without delayed neutrons
jeff3p2.ss = (2.76792, 1.22773e5)
jeff3p2.potential = 11.0983
jeff3p2.eFiss = 1.9E2
jeff3p2.dilutions = ( 1.e10, 10000.0, 3549.18335, 1259.67004, 447.079956, \
158.676849, 56.3173141 )
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "Cf249"
jeff3p2.mat = 9852
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-98-Cf-249.jeff32"
jeff3p2.fission = 1 # fission without delayed neutrons
jeff3p2.ss = (2.76792, 1.22773e5)
jeff3p2.potential = 11.1510
jeff3p2.eFiss = 1.9E2
jeff3p2.dilutions = ( 1.e10, 10000.0, 3549.18335, 1259.67004, 447.079956, \
158.676849, 56.3173141 )
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "Cf250"
jeff3p2.mat = 9855
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-98-Cf-250.jeff32"
jeff3p2.fission = 1 # fission without delayed neutrons
jeff3p2.ss = (2.76792, 1.22773e5)
jeff3p2.potential = 9.8800
jeff3p2.eFiss = 1.9E2
jeff3p2.dilutions = ( 1.e10, 10000.0, 3549.18335, 1259.67004, 447.079956, \
158.676849, 56.3173141 )
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "Cf251"
jeff3p2.mat = 9858
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-98-Cf-251.jeff32"
jeff3p2.fission = 1 # fission without delayed neutrons
jeff3p2.ss = (2.76792, 1.22773e5)
jeff3p2.potential = 9.8800
jeff3p2.eFiss = 1.9E2
jeff3p2.dilutions = ( 1.e10, 10000.0, 3549.18335, 1259.67004, 447.079956, \
158.676849, 56.3173141 )
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "Cf252"
jeff3p2.mat = 9861
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-98-Cf-252.jeff32"
jeff3p2.fission = 1 # fission without delayed neutrons
jeff3p2.ss = (2.76792, 1.22773e5)
jeff3p2.potential = 9.8000
jeff3p2.eFiss = 1.9E2
jeff3p2.dilutions = ( 1.e10, 10000.0, 3549.18335, 1259.67004, 447.079956, \
158.676849, 56.3173141 )
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "Cf253"
jeff3p2.mat = 9864
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-98-Cf-253.jeff32"
jeff3p2.fission = 1 # fission without delayed neutrons
jeff3p2.ss = (2.76792, 1.22773e5)
jeff3p2.potential = 9.7600
jeff3p2.eFiss = 1.9E2
jeff3p2.dilutions = ( 1.e10, 10000.0, 3549.18335, 1259.67004, 447.079956, \
158.676849, 56.3173141 )
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

# Process the fission products:

jeff3p2.scatteringLaw = None
jeff3p2.legendre = 0
jeff3p2.fission = None
jeff3p2.dilutions = None
jeff3p2.eFiss = None

jeff3p2.hmat = "Ge72"
jeff3p2.mat = 3231
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-32-Ge-072.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Ge73"
jeff3p2.mat = 3234
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-32-Ge-073.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Ge74"
jeff3p2.mat = 3237
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-32-Ge-074.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Ge76"
jeff3p2.mat = 3243
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-32-Ge-076.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "As75"
jeff3p2.mat = 3325
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-33-As-075.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Se76"
jeff3p2.mat = 3431
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-34-Se-076.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Se77"
jeff3p2.mat = 3434
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-34-Se-077.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Se78"
jeff3p2.mat = 3437
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-34-Se-078.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Se79"
jeff3p2.mat = 3440
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-34-Se-079.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Se80"
jeff3p2.mat = 3443
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-34-Se-080.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Se82"
jeff3p2.mat = 3449
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-34-Se-082.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Br79"
jeff3p2.mat = 3525
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-35-Br-079.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Br81"
jeff3p2.mat = 3531
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-35-Br-081.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Kr80"
jeff3p2.mat = 3631
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-36-Kr-080.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Kr82"
jeff3p2.mat = 3637
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-36-Kr-082.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Kr83"
jeff3p2.mat = 3640
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-36-Kr-083.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Kr84"
jeff3p2.mat = 3643
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-36-Kr-084.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Kr85"
jeff3p2.mat = 3646
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-36-Kr-085.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Kr86"
jeff3p2.mat = 3649
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-36-Kr-086.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Rb85"
jeff3p2.mat = 3725
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-37-Rb-085.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Rb87"
jeff3p2.mat = 3731
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-37-Rb-087.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Sr86"
jeff3p2.mat = 3831
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-38-Sr-086.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Sr87"
jeff3p2.mat = 3834
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-38-Sr-087.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Sr88"
jeff3p2.mat = 3837
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-38-Sr-088.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Sr89"
jeff3p2.mat = 3840
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-38-Sr-089.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Sr90"
jeff3p2.mat = 3843
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-38-Sr-090.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Y89"
jeff3p2.mat = 3925
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-39-Y-089.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Y90"
jeff3p2.mat = 3928
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-39-Y-090.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Y91"
jeff3p2.mat = 3931
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-39-Y-091.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Zr90"
jeff3p2.mat = 4025
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-40-Zr-090.jeff32"
jeff3p2.fission = None
jeff3p2.ss = (2.76792, 1.66156e4)
jeff3p2.potential = 6.5144
jeff3p2.dilutions = ( 1.e10, 10000.0,  3866.97, 1495.35, 578.2475, 223.6068, 86.4682, 33.4370, 12.9300, 5.0 )
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "Zr91"
jeff3p2.mat = 4028
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-40-Zr-091.jeff32"
jeff3p2.fission = None
jeff3p2.ss = (2.76792, 1.66156e4)
jeff3p2.potential = 6.5144
jeff3p2.dilutions = ( 1.e10, 10000.0, 3546.31, 1257.43, 445.8898, 158.1139, 56.0677, 19.8818, 7.0501, 2.5 )
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "Zr92"
jeff3p2.mat = 4031
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-40-Zr-092.jeff32"
jeff3p2.fission = None
jeff3p2.ss = (2.76792, 1.66156e4)
jeff3p2.potential = 6.5144
jeff3p2.dilutions = ( 1.e10, 10000.0, 3546.31, 1257.43, 445.8898, 158.1139, 56.0677, 19.8818, 7.0501, 2.5 )
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "Zr93"
jeff3p2.mat = 4034
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-40-Zr-093.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Zr94"
jeff3p2.mat = 4037
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-40-Zr-094.jeff32"
jeff3p2.fission = None
jeff3p2.ss = (2.76792, 1.66156e4)
jeff3p2.potential = 6.5144
jeff3p2.dilutions = ( 1.e10, 10000.0, 3546.31, 1257.43, 445.8898, 158.1139, 56.0677, 19.8818, 7.0501, 2.5 )
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "Zr95"
jeff3p2.mat = 4040
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-40-Zr-095.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Zr96"
jeff3p2.mat = 4043
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-40-Zr-096.jeff32"
jeff3p2.fission = None
jeff3p2.ss = (2.76792, 1.66156e4)
jeff3p2.potential = 6.5144
jeff3p2.dilutions = ( 1.e10, 10000.0, 3546.31, 1257.43, 445.8898, 158.1139, 56.0677, 19.8818, 7.0501, 2.5 )
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "Nb93"
jeff3p2.mat = 4125
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-41-Nb-093.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Nb94"
jeff3p2.mat = 4128
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-41-Nb-094.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Nb95"
jeff3p2.mat = 4131
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-41-Nb-095.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Mo92"
jeff3p2.mat = 4225
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-42-Mo-092.jeff32"
jeff3p2.fission = None
jeff3p2.ss = (2.76792, 1.22773e5)
jeff3p2.potential = 6.1140
jeff3p2.dilutions = ( 1.e10, 10000.0, 3546.31, 1257.43, 445.8898, 158.1139, 56.0677, 19.8818, 7.0501, 2.5 )
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "Mo94"
jeff3p2.mat = 4231
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-42-Mo-094.jeff32"
jeff3p2.fission = None
jeff3p2.ss = (2.76792, 1.22773e5)
jeff3p2.potential = 6.1140
jeff3p2.dilutions = ( 1.e10, 10000.0, 3546.31, 1257.43, 445.8898, 158.1139, 56.0677, 19.8818, 7.0501, 2.5 )
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "Mo95"
jeff3p2.mat = 4234
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-42-Mo-095.jeff32"
jeff3p2.fission = None
jeff3p2.ss = (2.76792, 1.22773e5)
jeff3p2.potential = 6.1140
jeff3p2.dilutions = ( 1.e10, 10000.0, 3546.31, 1257.43, 445.8898, 158.1139, 56.0677, 19.8818, 7.0501, 2.5 )
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "Mo96"
jeff3p2.mat = 4237
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-42-Mo-096.jeff32"
jeff3p2.fission = None
jeff3p2.ss = (2.76792, 1.22773e5)
jeff3p2.potential = 6.1140
jeff3p2.dilutions = ( 1.e10, 10000.0, 3546.31, 1257.43, 445.8898, 158.1139, 56.0677, 19.8818, 7.0501, 2.5 )
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "Mo97"
jeff3p2.mat = 4240
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-42-Mo-097.jeff32"
jeff3p2.fission = None
jeff3p2.ss = (2.76792, 1.22773e5)
jeff3p2.potential = 6.1140
jeff3p2.dilutions = ( 1.e10, 10000.0, 3546.31, 1257.43, 445.8898, 158.1139, 56.0677, 19.8818, 7.0501, 2.5 )
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "Mo98"
jeff3p2.mat = 4243
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-42-Mo-098.jeff32"
jeff3p2.fission = None
jeff3p2.ss = (2.76792, 1.22773e5)
jeff3p2.potential = 6.1140
jeff3p2.dilutions = ( 1.e10, 10000.0, 3546.31, 1257.43, 445.8898, 158.1139, 56.0677, 19.8818, 7.0501, 2.5 )
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "Mo99"
jeff3p2.mat = 4246
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-42-Mo-099.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Mo100"
jeff3p2.mat = 4249
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-42-Mo-100.jeff32"
jeff3p2.fission = None
jeff3p2.ss = (2.76792, 1.22773e5)
jeff3p2.potential = 6.1140
jeff3p2.dilutions = ( 1.e10, 10000.0, 3546.31, 1257.43, 445.8898, 158.1139, 56.0677, 19.8818, 7.0501, 2.5 )
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "Tc99"
jeff3p2.mat = 4331
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-43-Tc-099.jeff32"
jeff3p2.ss = (4.632489, 1.66156e4)
jeff3p2.potential = 6.4675
jeff3p2.dilutions = ( 1.e10, 10000.0, 4216.96552, 1778.27959, 749.894278, \
316.227791, 133.352152,  56.2341357, 23.7137381, 10.0 )
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "Ru99"
jeff3p2.mat = 4434
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-44-Ru-099.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Ru100"
jeff3p2.mat = 4437
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-44-Ru-100.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Ru101"
jeff3p2.mat = 4440
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-44-Ru-101.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Ru102"
jeff3p2.mat = 4443
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-44-Ru-102.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Ru103"
jeff3p2.mat = 4446
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-44-Ru-103.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Ru104"
jeff3p2.mat = 4449
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-44-Ru-104.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Ru105"
jeff3p2.mat = 4452
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-44-Ru-105.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Ru106"
jeff3p2.mat = 4455
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-44-Ru-106.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Rh103"
jeff3p2.mat = 4525
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-45-Rh-103.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Rh105"
jeff3p2.mat = 4531
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-45-Rh-105.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Pd104"
jeff3p2.mat = 4631
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-46-Pd-104.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Pd105"
jeff3p2.mat = 4634
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-46-Pd-105.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Pd106"
jeff3p2.mat = 4637
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-46-Pd-106.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Pd107"
jeff3p2.mat = 4640
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-46-Pd-107.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Pd108"
jeff3p2.mat = 4643
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-46-Pd-108.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Pd110"
jeff3p2.mat = 4649
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-46-Pd-110.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Ag107"
jeff3p2.mat = 4725
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-47-Ag-107.jeff32"
jeff3p2.ss = (2.76792, 1.66156e4)
jeff3p2.potential = 5.4739
jeff3p2.dilutions = ( 1.e10, 10000.0, 3546.31, 1257.43, 445.8898, 158.1139, 56.0677, 19.8818, 7.0501, 2.5 )
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "Ag109"
jeff3p2.mat = 4731
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-47-Ag-109.jeff32"
jeff3p2.ss = (2.76792, 1.66156e4)
jeff3p2.potential = 5.3316
jeff3p2.dilutions = ( 1.e10, 10000.0, 3546.31, 1257.43, 445.8898, 158.1139, 56.0677, 19.8818, 7.0501, 2.5 )
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "Ag111"
jeff3p2.mat = 4737
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-47-Ag-111.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Cd106"
jeff3p2.mat = 4825
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-48-Cd-106.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Cd108"
jeff3p2.mat = 4831
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-48-Cd-108.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Cd110"
jeff3p2.mat = 4837
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-48-Cd-110.jeff32"
jeff3p2.ss = (2.76792, 1.66156e4)
jeff3p2.potential = 5.1762
jeff3p2.dilutions = ( 1.e10, 10000.0, 3546.31, 1257.43, 445.8898, 158.1139, 56.0677, 19.8818, 7.0501, 2.5 )
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "Cd111"
jeff3p2.mat = 4840
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-48-Cd-111.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Cd112"
jeff3p2.mat = 4843
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-48-Cd-112.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Cd113"
jeff3p2.mat = 4846
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-48-Cd-113.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Cd114"
jeff3p2.mat = 4849
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-48-Cd-114.jeff32"
jeff3p2.branchingNG = 0.079383
jeff3p2.makeFp()
jeff3p2.branchingNG = None

jeff3p2.hmat = "Cd115m"
jeff3p2.mat = 4853
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-48-Cd-115M.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Cd116"
jeff3p2.mat = 4855
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-48-Cd-116.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "In113"
jeff3p2.mat = 4925
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-49-In-113.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "In115"
jeff3p2.mat = 4931
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-49-In-115.jeff32"
jeff3p2.ss = (2.76792, 1.66156e4)
jeff3p2.potential = 5.0695
jeff3p2.dilutions = ( 1.e10, 10000.0, 3546.31, 1257.43, 445.8898, 158.1139, 56.0677, 19.8818, 7.0501, 2.5 )
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "Sn112"
jeff3p2.mat = 5025
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-50-Sn-112.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Sn114"
jeff3p2.mat = 5031
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-50-Sn-114.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Sn115"
jeff3p2.mat = 5034
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-50-Sn-115.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Sn116"
jeff3p2.mat = 5037
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-50-Sn-116.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Sn117"
jeff3p2.mat = 5040
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-50-Sn-117.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Sn118"
jeff3p2.mat = 5043
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-50-Sn-118.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Sn119"
jeff3p2.mat = 5046
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-50-Sn-119.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Sn120"
jeff3p2.mat = 5049
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-50-Sn-120.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Sn122"
jeff3p2.mat = 5055
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-50-Sn-122.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Sn123"
jeff3p2.mat = 5058
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-50-Sn-123.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Sn124"
jeff3p2.mat = 5061
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-50-Sn-124.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Sn125"
jeff3p2.mat = 5064
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-50-Sn-125.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Sn126"
jeff3p2.mat = 5067
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-50-Sn-126.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Sb121"
jeff3p2.mat = 5125
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-51-Sb-121.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Sb123"
jeff3p2.mat = 5131
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-51-Sb-123.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Sb124"
jeff3p2.mat = 5134
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-51-Sb-124.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Sb125"
jeff3p2.mat = 5137
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-51-Sb-125.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Sb126"
jeff3p2.mat = 5140
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-51-Sb-126.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Te122"
jeff3p2.mat = 5231
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-52-Te-122.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Te123"
jeff3p2.mat = 5234
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-52-Te-123.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Te124"
jeff3p2.mat = 5237
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-52-Te-124.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Te125"
jeff3p2.mat = 5240
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-52-Te-125.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Te126"
jeff3p2.mat = 5243
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-52-Te-126.jeff32"
jeff3p2.branchingNG = 0.091528
jeff3p2.makeFp()
jeff3p2.branchingNG = None

jeff3p2.hmat = "Te127m"
jeff3p2.mat = 5247
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-52-Te-127M.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Te128"
jeff3p2.mat = 5249
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-52-Te-128.jeff32"
jeff3p2.branchingNG = 0.031894
jeff3p2.makeFp()
jeff3p2.branchingNG = None

jeff3p2.hmat = "Te129m"
jeff3p2.mat = 5253
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-52-Te-129M.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Te130"
jeff3p2.mat = 5255
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-52-Te-130.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Te131m"
jeff3p2.mat = 5259
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-52-Te-131M.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Te132"
jeff3p2.mat = 5261
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-52-Te-132.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "I127"
jeff3p2.mat = 5325
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-53-I-127.jeff32"
jeff3p2.ss = (2.76792, 1.66156e4)
jeff3p2.potential = 4.5239
jeff3p2.dilutions = ( 1.e10, 10000.0, 4216.96552, 1778.27959, 749.894278, \
316.227791, 133.352152,  56.2341357, 23.7137381, 10.0 )
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "I129"
jeff3p2.mat = 5331
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-53-I-129.jeff32"
jeff3p2.ss = (2.76792, 1.66156e4)
jeff3p2.potential = 5.8221
jeff3p2.dilutions = ( 1.e10, 10000.0, 4216.96552, 1778.27959, 749.894278, \
316.227791, 133.352152,  56.2341357, 23.7137381, 10.0 )
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "I130"
jeff3p2.mat = 5334
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-53-I-130.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "I131"
jeff3p2.mat = 5337
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-53-I-131.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "I135"
jeff3p2.mat = 5349
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-53-I-135.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Xe128"
jeff3p2.mat = 5437
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-54-Xe-128.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Xe129"
jeff3p2.mat = 5440
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-54-Xe-129.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Xe130"
jeff3p2.mat = 5443
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-54-Xe-130.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Xe131"
jeff3p2.mat = 5446
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-54-Xe-131.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Xe132"
jeff3p2.mat = 5449
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-54-Xe-132.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Xe133"
jeff3p2.mat = 5452
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-54-Xe-133.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Xe134"
jeff3p2.mat = 5455
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-54-Xe-134.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Xe135"
jeff3p2.mat = 5458
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-54-Xe-135.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Xe136"
jeff3p2.mat = 5461
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-54-Xe-136.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Cs133"
jeff3p2.mat = 5525
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-55-Cs-133.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Cs134"
jeff3p2.mat = 5528
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-55-Cs-134.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Cs135"
jeff3p2.mat = 5531
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-55-Cs-135.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Cs136"
jeff3p2.mat = 5534
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-55-Cs-136.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Cs137"
jeff3p2.mat = 5537
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-55-Cs-137.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Ba134"
jeff3p2.mat = 5637
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-56-Ba-134.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Ba135"
jeff3p2.mat = 5640
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-56-Ba-135.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Ba136"
jeff3p2.mat = 5643
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-56-Ba-136.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Ba137"
jeff3p2.mat = 5646
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-56-Ba-137.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Ba138"
jeff3p2.mat = 5649
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-56-Ba-138.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Ba140"
jeff3p2.mat = 5655
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-56-Ba-140.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "La138"
jeff3p2.mat = 5725
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-57-La-138.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "La139"
jeff3p2.mat = 5728
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-57-La-139.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "La140"
jeff3p2.mat = 5731
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-57-La-140.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Ce140"
jeff3p2.mat = 5837
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-58-Ce-140.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Ce141"
jeff3p2.mat = 5840
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-58-Ce-141.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Ce142"
jeff3p2.mat = 5843
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-58-Ce-142.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Ce143"
jeff3p2.mat = 5846
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-58-Ce-143.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Ce144"
jeff3p2.mat = 5849
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-58-Ce-144.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Pr141"
jeff3p2.mat = 5925
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-59-Pr-141.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Pr142"
jeff3p2.mat = 5928
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-59-Pr-142.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Pr143"
jeff3p2.mat = 5931
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-59-Pr-143.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Nd142"
jeff3p2.mat = 6025
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-60-Nd-142.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Nd143"
jeff3p2.mat = 6028
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-60-Nd-143.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Nd144"
jeff3p2.mat = 6031
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-60-Nd-144.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Nd145"
jeff3p2.mat = 6034
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-60-Nd-145.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Nd146"
jeff3p2.mat = 6037
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-60-Nd-146.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Nd147"
jeff3p2.mat = 6040
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-60-Nd-147.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Nd148"
jeff3p2.mat = 6043
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-60-Nd-148.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Nd150"
jeff3p2.mat = 6049
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-60-Nd-150.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Pm147"
jeff3p2.mat = 6149
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-61-Pm-147.jeff32"
jeff3p2.branchingNG = 0.470
jeff3p2.makeFp()
jeff3p2.branchingNG = None

jeff3p2.hmat = "Pm148"
jeff3p2.mat = 6152
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-61-Pm-148.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Pm148m"
jeff3p2.mat = 6153
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-61-Pm-148M.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Pm149"
jeff3p2.mat = 6155
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-61-Pm-149.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Pm151"
jeff3p2.mat = 6161
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-61-Pm-151.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Sm147"
jeff3p2.mat = 6234
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-62-Sm-147.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Sm148"
jeff3p2.mat = 6237
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-62-Sm-148.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Sm149"
jeff3p2.mat = 6240
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-62-Sm-149.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Sm150"
jeff3p2.mat = 6243
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-62-Sm-150.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Sm151"
jeff3p2.mat = 6246
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-62-Sm-151.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Sm152"
jeff3p2.mat = 6249
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-62-Sm-152.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Sm153"
jeff3p2.mat = 6252
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-62-Sm-153.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Sm154"
jeff3p2.mat = 6255
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-62-Sm-154.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Eu151"
jeff3p2.mat = 6325
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-63-Eu-151.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Eu152"
jeff3p2.mat = 6328
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-63-Eu-152.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Eu153"
jeff3p2.mat = 6331
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-63-Eu-153.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Eu154"
jeff3p2.mat = 6334
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-63-Eu-154.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Eu155"
jeff3p2.mat = 6337
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-63-Eu-155.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Eu156"
jeff3p2.mat = 6340
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-63-Eu-156.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Eu157"
jeff3p2.mat = 6343
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-63-Eu-157.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Gd152"
jeff3p2.mat = 6425
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-64-Gd-152.jeff32"
jeff3p2.ss = (2.76792, 1.66156e4)
jeff3p2.potential = 8.0425
jeff3p2.dilutions = ( 1.e10, 10000.0, 4216.96552, 1778.27959, 749.894278, \
316.227791, 133.352152,  56.2341357, 23.7137381, 10.0 )
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "Gd154"
jeff3p2.mat = 6431
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-64-Gd-154.jeff32"
jeff3p2.ss = (2.76792, 1.66156e4)
jeff3p2.potential = 6.6723
jeff3p2.dilutions = ( 1.e10, 10000.0, 4216.96552, 1778.27959, 749.894278, \
316.227791, 133.352152,  56.2341357, 23.7137381, 10.0 )
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "Gd155"
jeff3p2.mat = 6434
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-64-Gd-155.jeff32"
jeff3p2.ss = (2.76792, 1.66156e4)
jeff3p2.potential = 6.3376
jeff3p2.dilutions = ( 1.e10, 10000.0, 4216.96552, 1778.27959, 749.894278, \
316.227791, 133.352152,  56.2341357, 23.7137381, 10.0 )
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "Gd156"
jeff3p2.mat = 6437
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-64-Gd-156.jeff32"
jeff3p2.ss = (2.76792, 1.66156e4)
jeff3p2.potential = 7.3792
jeff3p2.dilutions = ( 1.e10, 10000.0, 4216.96552, 1778.27959, 749.894278, \
316.227791, 133.352152,  56.2341357, 23.7137381, 10.0 )
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "Gd157"
jeff3p2.mat = 6440
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-64-Gd-157.jeff32"
jeff3p2.ss = (2.76792, 1.66156e4)
jeff3p2.potential = 6.6063
jeff3p2.dilutions = ( 1.e10, 10000.0, 4216.96552, 1778.27959, 749.894278, \
316.227791, 133.352152,  56.2341357, 23.7137381, 10.0 )
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "Gd158"
jeff3p2.mat = 6443
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-64-Gd-158.jeff32"
jeff3p2.ss = (2.76792, 1.66156e4)
jeff3p2.potential = 7.6454
jeff3p2.dilutions = ( 1.e10, 10000.0, 4216.96552, 1778.27959, 749.894278, \
316.227791, 133.352152,  56.2341357, 23.7137381, 10.0 )
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "Gd160"
jeff3p2.mat = 6449
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-64-Gd-160.jeff32"
jeff3p2.ss = (2.76792, 1.66156e4)
jeff3p2.potential = 7.0241
jeff3p2.dilutions = ( 1.e10, 10000.0, 4216.96552, 1778.27959, 749.894278, \
316.227791, 133.352152,  56.2341357, 23.7137381, 10.0 )
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.dilutions = None

jeff3p2.hmat = "Tb159"
jeff3p2.mat = 6525
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-65-Tb-159.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Tb160"
jeff3p2.mat = 6528
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-65-Tb-160.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Dy160"
jeff3p2.mat = 6637
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-66-Dy-160.jeff32"
jeff3p2.ss = (2.76792, 1.66156e4)
jeff3p2.potential = 6.9861
jeff3p2.dilutions = ( 1.e10, 10000.0, 4216.96552, 1778.27959, 749.894278, \
316.227791, 133.352152,  56.2341357, 23.7137381, 10.0 )
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "Dy161"
jeff3p2.mat = 6640
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-66-Dy-161.jeff32"
jeff3p2.ss = (2.76792, 1.66156e4)
jeff3p2.potential = 7.0121
jeff3p2.dilutions = ( 1.e10, 10000.0, 4216.96552, 1778.27959, 749.894278, \
316.227791, 133.352152,  56.2341357, 23.7137381, 10.0 )
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "Dy162"
jeff3p2.mat = 6643
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-66-Dy-162.jeff32"
jeff3p2.ss = (2.76792, 1.66156e4)
jeff3p2.potential = 4.5681
jeff3p2.dilutions = ( 1.e10, 10000.0, 4216.96552, 1778.27959, 749.894278, \
316.227791, 133.352152,  56.2341357, 23.7137381, 10.0 )
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "Dy163"
jeff3p2.mat = 6646
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-66-Dy-163.jeff32"
jeff3p2.ss = (2.76792, 1.66156e4)
jeff3p2.potential = 7.0639
jeff3p2.dilutions = ( 1.e10, 10000.0, 4216.96552, 1778.27959, 749.894278, \
316.227791, 133.352152,  56.2341357, 23.7137381, 10.0 )
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "Dy164"
jeff3p2.mat = 6649
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-66-Dy-164.jeff32"
jeff3p2.ss = (2.76792, 1.66156e4)
jeff3p2.potential = 7.0897
jeff3p2.dilutions = ( 1.e10, 10000.0, 4216.96552, 1778.27959, 749.894278, \
316.227791, 133.352152,  56.2341357, 23.7137381, 10.0 )
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "Dy165"
jeff3p2.mat = 6652
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-66-Dy-165.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Hf174"
jeff3p2.mat = 7225
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-72-Hf-174.jeff32"
jeff3p2.ss = (2.76792, 1.66156e4)
jeff3p2.potential = 7.1385
jeff3p2.dilutions = ( 1.e10, 10000.0, 4216.96552, 1778.27959, 749.894278, \
316.227791, 133.352152,  56.2341357, 23.7137381, 10.0 )
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "Hf176"
jeff3p2.mat = 7231
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-72-Hf-176.jeff32"
jeff3p2.ss = (2.76792, 1.66156e4)
jeff3p2.potential = 7.1935
jeff3p2.dilutions = ( 1.e10, 10000.0, 4216.96552, 1778.27959, 749.894278, \
316.227791, 133.352152,  56.2341357, 23.7137381, 10.0 )
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "Hf177"
jeff3p2.mat = 7234
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-72-Hf-177.jeff32"
jeff3p2.ss = (2.76792, 1.66156e4)
jeff3p2.potential = 7.2202
jeff3p2.dilutions = ( 1.e10, 10000.0, 4216.96552, 1778.27959, 749.894278, \
316.227791, 133.352152,  56.2341357, 23.7137381, 10.0 )
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "Hf178"
jeff3p2.mat = 7237
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-72-Hf-178.jeff32"
jeff3p2.ss = (2.76792, 1.66156e4)
jeff3p2.potential = 7.2469
jeff3p2.dilutions = ( 1.e10, 10000.0, 4216.96552, 1778.27959, 749.894278, \
316.227791, 133.352152,  56.2341357, 23.7137381, 10.0 )
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "Hf179"
jeff3p2.mat = 7240
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-72-Hf-179.jeff32"
jeff3p2.ss = (2.76792, 1.66156e4)
jeff3p2.potential = 7.2736
jeff3p2.dilutions = ( 1.e10, 10000.0, 4216.96552, 1778.27959, 749.894278, \
316.227791, 133.352152,  56.2341357, 23.7137381, 10.0 )
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.hmat = "Hf180"
jeff3p2.mat = 7243
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-72-Hf-180.jeff32"
jeff3p2.ss = (2.76792, 1.66156e4)
jeff3p2.potential = 7.3004
jeff3p2.dilutions = ( 1.e10, 10000.0, 4216.96552, 1778.27959, 749.894278, \
316.227791, 133.352152,  56.2341357, 23.7137381, 10.0 )
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

jeff3p2.dilutions = None

jeff3p2.hmat = "Ho165"
jeff3p2.mat = 6725
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-67-Ho-165.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Er166"
jeff3p2.mat = 6837
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-68-Er-166.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Er167"
jeff3p2.mat = 6840
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-68-Er-167.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Er168"
jeff3p2.mat = 6843
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-68-Er-168.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Er170"
jeff3p2.mat = 6849
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-68-Er-170.jeff32"
jeff3p2.makeFp()

jeff3p2.hmat = "Lu176"
jeff3p2.mat = 7128
jeff3p2.evaluationFile = "$HOME/evaluations/Jeff3.2/JEFF32N/n-71-Lu-176.jeff32"
jeff3p2.makeFp()

# Process the burnup chain (Jeff3.1.2):

jeff3p2.fissionFile = "$HOME/evaluations/Jeff3.1.2/JEFF311NFY"
jeff3p2.decayFile = "$HOME/evaluations/Jeff3.1.2/JEFF311RDD_ALL"
jeff3p2.burnup()
