#!/usr/local/bin/python
from PyNjoy import *
from os import uname
endfb = PyNjoy()
endfb.evaluationName = "/tmp/shem295_endfb7r1"
endfb.execDir = "../" + uname()[0]
endfb.nstr = 25
endfb.iwt = 4
endfb.Espectra = None
endfb.autolib = (4.632489, 1.11377e4, 0.0005)

endfb.legendre = 3
endfb.hmat = "H1_H2O"
endfb.mat = 125
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-001_H_001.endf"
endfb.scatteringLaw = "$HOME/evaluations/ENDFB7r1/thermal_scatt/tsl-HinH2O.endf"
endfb.scatteringMat = 1
endfb.temperatures = ( 293.6, 350.0, 400.0, 450.0, 500.0, 600.0 )
endfb.fission = None
endfb.dilutions = None
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.legendre = 1
endfb.hmat = "H2_D2O"
endfb.mat = 128
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-001_H_002.endf"
endfb.scatteringLaw = "$HOME/evaluations/ENDFB7r1/thermal_scatt/tsl-DinD2O.endf"
endfb.scatteringMat = 11
endfb.temperatures = ( 293.6, 350.0, 400.0, 450.0, 500.0, 600.0 )
endfb.fission = None
endfb.dilutions = None
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "H1_CH2"
endfb.mat = 125
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-001_H_001.endf"
endfb.scatteringLaw = "$HOME/evaluations/ENDFB7r1/thermal_scatt/tsl-HinCH2.endf"
endfb.scatteringMat = 37
endfb.temperatures = ( 296.,  350. )
endfb.fission = None
endfb.dilutions = None
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "H1_ZRH"
endfb.mat = 125
endfb.temperatures = ( 296., 400., 500., 600., 700., 800., 1000., 1200.)
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-001_H_001.endf"
endfb.scatteringLaw = "$HOME/evaluations/ENDFB7r1/thermal_scatt/tsl-HinZrH.endf"
endfb.scatteringMat = 7
endfb.fission = None
endfb.dilutions = None
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Zr90_ZRH"
endfb.mat =  4025
endfb.temperatures = ( 296., 400., 500., 600., 700., 800., 1000., 1200.)
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-040_Zr_090.endf"
endfb.scatteringLaw =  "$HOME/evaluations/ENDFB7r1/thermal_scatt/tsl-ZrinZrH.endf"
endfb.scatteringMat = 58
endfb.fission = None
endfb.dilutions = None
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "C0_GR"
endfb.mat = 600
endfb.temperatures = ( 296., 400., 500., 600., 700., 800., 1000., 1200., 1600., 2000.)
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-006_C_000.endf"
endfb.scatteringLaw = "$HOME/evaluations/ENDFB7r1/thermal_scatt/tsl-graphite.endf"
endfb.scatteringMat = 31
endfb.fission = None
endfb.dilutions = None
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Be9"
endfb.mat = 425  
endfb.temperatures = ( 296., 400., 500., 600., 700., 800., 1000., 1200.)
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-004_Be_009.endf"
endfb.scatteringLaw =  "$HOME/evaluations/ENDFB7r1/thermal_scatt/tsl-Be-metal.endf"
endfb.scatteringMat = 26
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
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-001_H_001.endf"
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "H2"
endfb.mat = 128
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-001_H_002.endf"
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "H3"
endfb.mat = 131
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-001_H_003.endf"
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "He3"
endfb.mat = 225
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-002_He_003.endf"
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "He4"
endfb.mat = 228
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-002_He_004.endf"
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Li6"
endfb.mat = 325
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-003_Li_006.endf"
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Li7"
endfb.mat = 328
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-003_Li_007.endf"
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Be7"
endfb.mat = 419
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-004_Be_007.endf"
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "B10"
endfb.mat = 525
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-005_B_010.endf"
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "B11"
endfb.mat = 528
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-005_B_011.endf"
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "C0"
endfb.mat = 600
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-006_C_000.endf"
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "N14"
endfb.mat = 725
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-007_N_014.endf"
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "N15"
endfb.mat = 728
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-007_N_015.endf"
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "O16"
endfb.mat = 825
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-008_O_016.endf"
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "O17"
endfb.mat = 828
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-008_O_017.endf"
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "F19"
endfb.mat = 925
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-009_F_019.endf"
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Na23"
endfb.mat = 1125
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-011_Na_023.endf"
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Mg24"
endfb.mat = 1225
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-012_Mg_024.endf"
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Mg25"
endfb.mat = 1228
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-012_Mg_025.endf"
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Mg26"
endfb.mat = 1231
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-012_Mg_026.endf"
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Al27"
endfb.mat = 1325
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-013_Al_027.endf"
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Si28"
endfb.mat = 1425
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-014_Si_028.endf"
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Si29"
endfb.mat = 1428
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-014_Si_029.endf"
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Si30"
endfb.mat = 1431
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-014_Si_030.endf"
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "P31"
endfb.mat = 1525
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-015_P_031.endf"
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "S32"
endfb.mat = 1625
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-016_S_032.endf"
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "S33"
endfb.mat = 1628
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-016_S_033.endf"
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "S34"
endfb.mat = 1631
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-016_S_034.endf"
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "S36"
endfb.mat = 1637
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-016_S_036.endf"
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Cl35"
endfb.mat = 1725
#endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-017_Cl_035.endf"
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r0/n-ENDF-VII0.endf/n-017_Cl_035.endf"
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Cl37"
endfb.mat = 1731
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-017_Cl_037.endf"
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "K39"
endfb.mat = 1925
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-019_K_039.endf"
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "K40"
endfb.mat = 1928
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-019_K_040.endf"
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "K41"
endfb.mat = 1931
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-019_K_041.endf"
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Ca40"
endfb.mat = 2025
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-020_Ca_040.endf"
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Ca42"
endfb.mat = 2031
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-020_Ca_042.endf"
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Ca43"
endfb.mat = 2034
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-020_Ca_043.endf"
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Ca44"
endfb.mat = 2037
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-020_Ca_044.endf"
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Ca46"
endfb.mat = 2043
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-020_Ca_046.endf"
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Ca48"
endfb.mat = 2049
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-020_Ca_048.endf"
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Ti46"
endfb.mat = 2225
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-022_Ti_046.endf"
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Ti47"
endfb.mat = 2228
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-022_Ti_047.endf"
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Ti48"
endfb.mat = 2231
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-022_Ti_048.endf"
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Ti49"
endfb.mat = 2234
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-022_Ti_049.endf"
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Ti50"
endfb.mat = 2237
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-022_Ti_050.endf"
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "V50"
endfb.mat = 2325
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-023_V_050.endf"
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "V51"
endfb.mat = 2328
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-023_V_051.endf"
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Cr50"
endfb.mat = 2425
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-024_Cr_050.endf"
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Cr52"
endfb.mat = 2431
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-024_Cr_052.endf"
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Cr53"
endfb.mat = 2434
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-024_Cr_053.endf"
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Cr54"
endfb.mat = 2437
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-024_Cr_054.endf"
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Mn55"
endfb.mat = 2525
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-025_Mn_055.endf"
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Fe54"
endfb.mat = 2625
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-026_Fe_054.endf"
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Fe56"
endfb.mat = 2631
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-026_Fe_056.endf"
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Fe57"
endfb.mat = 2634
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-026_Fe_057.endf"
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Fe58"
endfb.mat = 2637
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-026_Fe_058.endf"
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Co59"
endfb.mat = 2725
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-027_Co_059.endf"
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Ni58"
endfb.mat = 2825
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-028_Ni_058.endf"
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Ni60"
endfb.mat = 2831
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-028_Ni_060.endf"
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Ni61"
endfb.mat = 2834
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-028_Ni_061.endf"
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Ni62"
endfb.mat = 2837
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-028_Ni_062.endf"
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Ni64"
endfb.mat = 2843
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-028_Ni_064.endf"
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Cu63"
endfb.mat = 2925
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-029_Cu_063.endf"
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Cu65"
endfb.mat = 2931
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-029_Cu_065.endf"
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Zn66"
endfb.mat = 3031
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-030_Zn_066.endf"
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Zn67"
endfb.mat = 3034
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-030_Zn_067.endf"
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Zn68"
endfb.mat = 3037
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-030_Zn_068.endf"
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Zn70"
endfb.mat = 3043
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-030_Zn_070.endf"
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Ga69"
endfb.mat = 3125
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-031_Ga_069.endf"
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Ga71"
endfb.mat = 3131
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-031_Ga_071.endf"
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "W182"
endfb.mat =  7431  
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-074_W_182.endf"
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "W183"
endfb.mat = 7434 
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-074_W_183.endf"
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "W184"
endfb.mat = 7437
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-074_W_184.endf"
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "W186"
endfb.mat =  7443
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-074_W_186.endf"
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Pb204"
endfb.mat =  8225
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-082_Pb_204.endf"
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Pb206"
endfb.mat =  8231 
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-082_Pb_206.endf"
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Pb207"
endfb.mat =  8234
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-082_Pb_207.endf"
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Pb208"
endfb.mat =  8237
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-082_Pb_208.endf"
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Bi209"
endfb.mat =  8325
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-083_Bi_209.endf"
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Th230"
endfb.mat = 9034
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-090_Th_230.endf"
endfb.fission = 1 # fission without delayed neutrons
endfb.ss = (4.632489, 3.206464e5)
endfb.potential = 8.7040
endfb.eFiss = 1.7406E+02
endfb.dilutions = ( 1.e10, 10000.0, 3549.18335, 1259.67004, 447.079956, \
158.676849, 56.3173141, 19.9880447, 7.09412289, 2.51783395 )
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Th232"
endfb.mat = 9040
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-090_Th_232.endf"
endfb.fission = 2 # fission with delayed neutrons
endfb.ss = (4.632489, 3.206464e5)
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
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-091_Pa_231.endf"
endfb.fission = 1 # fission without delayed neutrons
endfb.ss = (4.632489, 3.206464e5)
endfb.potential = 8.7266
endfb.eFiss = 1.7E2
endfb.dilutions = ( 1.e10, 10000.0, 3549.18335, 1259.67004, 447.079956, \
158.676849, 56.3173141, 19.9880447, 7.09412289, 2.51783395 )
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Pa233"
endfb.mat = 9137
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-091_Pa_233.endf"
endfb.fission = 1 # fission without delayed neutrons
endfb.ss = (4.632489, 3.206464e5)
endfb.potential = 9.9950
endfb.eFiss = 1.7558E+02
endfb.dilutions = ( 1.e10, 10000.0, 3549.18335, 1259.67004, 447.079956, \
158.676849, 56.3173141, 19.9880447, 7.09412289, 2.51783395 )
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "U232"
endfb.mat = 9219
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-092_U_232.endf"
endfb.fission = 1 # fission without delayed neutrons
endfb.ss = (4.632489, 3.206464e5)
endfb.potential = 12.0687
endfb.eFiss = 1.8E2
endfb.dilutions = ( 1.e10, 10000.0, 3549.18335, 1259.67004, 447.079956, \
158.676849, 56.3173141, 19.9880447, 7.09412289, 2.51783395 )
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "U233"
endfb.mat = 9222
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-092_U_233.endf"
endfb.fission = 2 # fission with delayed neutrons
endfb.ss = (4.632489, 3.206464e5)
endfb.potential = 12.2989
endfb.eFiss = 1.8087E+02
endfb.dilutions = ( 1.e10, 10000.0, 3549.18335, 1259.67004, 447.079956, \
158.676849, 56.3173141, 19.9880447, 7.09412289, 2.51783395 )
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "U234"
endfb.mat = 9225
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-092_U_234.endf"
endfb.fission = 1 # fission without delayed neutrons
endfb.ss = (4.632489, 3.206464e5)
endfb.potential = 10.0210
endfb.dilutions = ( 1.e10, 10000.0, 3549.18335, 1259.67004, 447.079956, \
158.676849, 56.3173141, 19.9880447, 7.09412289, 2.51783395 )
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "U235"
endfb.mat = 9228
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-092_U_235.endf"
endfb.fission = 2 # fission with delayed neutrons
endfb.ss = (4.632489, 3.206464e5)
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
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-092_U_236.endf"
endfb.fission = 1 # fission without delayed neutrons
endfb.ss = (4.632489, 3.206464e5)
endfb.potential = 10.9954
endfb.dilutions = ( 1.e10, 10000.0, 3549.18335, 1259.67004, 447.079956, \
158.676849, 56.3173141, 19.9880447, 7.09412289, 2.51783395 )
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "U237"
endfb.mat = 9234
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-092_U_237.endf"
endfb.fission = 1 # fission without delayed neutrons
endfb.ss = (4.632489, 3.206464e5)
endfb.potential = 10.5000
endfb.dilutions = ( 1.e10, 10000.0, 3549.18335, 1259.67004, 447.079956, \
158.676849, 56.3173141, 19.9880447, 7.09412289, 2.51783395 )
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "U238"
endfb.mat = 9237
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-092_U_238.endf"
endfb.fission = 2 # fission with delayed neutrons
endfb.ss = (4.632489, 3.206464e5)
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
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-093_Np_236.endf"
endfb.fission = 1 # fission without delayed neutrons
endfb.ss = (4.632489, 3.206464e5)
endfb.eFiss = 1.8E2
endfb.potential = 12.1427
endfb.dilutions = ( 1.e10, 10000.0, 3549.18335, 1259.67004, 447.079956, \
158.676849, 56.3173141, 19.9880447, 7.09412289, 2.51783395 )
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Np237"
endfb.mat = 9346
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-093_Np_237.endf"
endfb.fission = 1 # fission without delayed neutrons
endfb.ss = (4.632489, 3.206464e5)
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
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-093_Np_238.endf"
endfb.fission = 1 # fission without delayed neutrons
endfb.ss = (4.632489, 3.206464e5)
endfb.eFiss = 1.8E2
endfb.potential = 10.4000
endfb.dilutions = ( 1.e10, 10000.0, 3549.18335, 1259.67004, 447.079956, \
158.676849, 56.3173141 )
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Np239"
endfb.mat = 9352
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-093_Np_239.endf"
endfb.fission = 1 # fission without delayed neutrons
endfb.ss = (4.632489, 3.206464e5)
endfb.eFiss = 1.8E2
endfb.potential = 10.4979
endfb.dilutions = ( 1.e10, 10000.0, 3549.18335, 1259.67004, 447.079956, \
158.676849, 56.3173141, 19.9880447, 7.09412289, 2.51783395 )
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Pu236"
endfb.mat = 9428
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-094_Pu_236.endf"
endfb.fission = 1 # fission without delayed neutrons
endfb.ss = (4.632489, 3.206464e5)
endfb.eFiss = 1.8E2
endfb.potential = 11.2458
endfb.dilutions = ( 1.e10, 10000.0, 3549.18335, 1259.67004, 447.079956, \
158.676849, 56.3173141, 19.9880447, 7.09412289, 2.51783395 )
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Pu237"
endfb.mat = 9431
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-094_Pu_237.endf"
endfb.fission = 1 # fission without delayed neutrons
endfb.ss = (4.632489, 3.206464e5)
endfb.potential = 10.5209
endfb.dilutions = ( 1.e10, 10000.0, 3549.18335, 1259.67004, 447.079956, \
158.676849, 56.3173141, 19.9880447, 7.09412289, 2.51783395 )
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Pu238"
endfb.mat = 9434
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-094_Pu_238.endf"
endfb.fission = 1 # fission without delayed neutrons
endfb.ss = (4.632489, 3.206464e5)
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
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-094_Pu_239.endf"
endfb.fission = 2 # fission with delayed neutrons
endfb.ss = (4.632489, 3.206464e5)
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
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-094_Pu_240.endf"
endfb.fission = 2 # fission with delayed neutrons
endfb.ss = (4.632489, 3.206464e5)
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
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-094_Pu_241.endf"
endfb.fission = 2 # fission with delayed neutrons
endfb.ss = (4.632489, 3.206464e5)
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
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-094_Pu_242.endf"
endfb.fission = 1 # fission without delayed neutrons
endfb.ss = (4.632489, 3.206464e5)
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
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-094_Pu_243.endf"
endfb.fission = 1 # fission without delayed neutrons
endfb.ss = (4.632489, 3.206464e5)
endfb.potential = 10.2000
endfb.eFiss = 1.9E2
endfb.dilutions = ( 1.e10, 10000.0, 3549.18335, 1259.67004, 447.079956, \
158.676849, 56.3173141, 19.9880447, 7.09412289, 2.51783395 )
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Pu244"
endfb.mat = 9452
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-094_Pu_244.endf"
endfb.fission = 1 # fission without delayed neutrons
endfb.ss = (4.632489, 3.206464e5)
endfb.potential = 10.1000
endfb.dilutions = ( 1.e10, 10000.0, 3549.18335, 1259.67004, 447.079956, \
158.676849, 56.3173141, 19.9880447, 7.09412289, 2.51783395 )
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Am241"
endfb.mat = 9543
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-095_Am_241.endf"
endfb.fission = 1 # fission without delayed neutrons
endfb.ss = (4.632489, 3.206464e5)
endfb.potential = 11.0329
endfb.branchingNG = 0.115
endfb.dilutions = ( 1.e10, 94.5317612, 56.3173141, 33.5510521, 19.9880447, \
11.9078817, 7.09412289, 4.22632504, 2.51783395, 1.5 )
endfb.pendf()
endfb.gendf()
endfb.draglib()
endfb.branchingNG = 0.115
endfb.dilutions = ( 1.e10, 10000.0, 5957.50244, 3549.18335, 2114.42676, \
1259.67004, 750.448669, 447.079956, 266.347961, 158.676849 )
endfb.pendf()
endfb.gendf()
endfb.draglib()
endfb.branchingNG = None

endfb.hmat = "Am242"
endfb.mat = 9546
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-095_Am_242.endf"
endfb.fission = 0 # no fission matrix!!!
endfb.ss = (4.632489, 3.206464e5)
endfb.potential = 10.2000
endfb.eFiss = 1.9E2
endfb.dilutions = ( 1.e10, 10000.0, 3549.18335, 1259.67004, 447.079956, \
158.676849, 56.3173141, 19.9880447, 7.09412289, 2.51783395 )
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Am242m"
endfb.mat = 9547
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-095_Am_242m1.endf"
endfb.fission = 1 # fission without delayed neutrons
endfb.ss = (4.632489, 3.206464e5)
endfb.potential = 10.2000
endfb.eFiss = 1.9E2
endfb.dilutions = ( 1.e10, 10000.0, 3549.18335, 1259.67004, 447.079956, \
158.676849, 56.3173141 )
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Am243"
endfb.mat = 9549
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-095_Am_243.endf"
endfb.fission = 1 # fission without delayed neutrons
endfb.ss = (4.632489, 3.206464e5)
endfb.potential = 11.8237
endfb.dilutions = ( 1.e10, 10000.0, 3549.18335, 1259.67004, 447.079956, \
158.676849, 56.3173141, 19.9880447, 7.09412289, 2.51783395 )
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Cm241"
endfb.mat = 9628
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-096_Cm_241.endf"
endfb.fission = 1 # fission without delayed neutrons
endfb.ss = (4.632489, 3.206464e5)
endfb.potential = 10.1788
endfb.eFiss = 1.9245E+02
endfb.dilutions = ( 1.e10, 10000.0, 3549.18335, 1259.67004, 447.079956, \
158.676849, 56.3173141 )
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Cm242"
endfb.mat = 9631
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-096_Cm_242.endf"
endfb.fission = 1 # fission without delayed neutrons
endfb.ss = (4.632489, 3.206464e5)
endfb.potential = 10.2000
endfb.eFiss = 1.9142E+02
endfb.dilutions = ( 1.e10, 10000.0, 3549.18335, 1259.67004, 447.079956, \
158.676849, 56.3173141 )
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Cm243"
endfb.mat = 9634
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-096_Cm_243.endf"
endfb.fission = 1 # fission without delayed neutrons
endfb.ss = (4.632489, 3.206464e5)
endfb.potential = 11.2832
endfb.eFiss = 1.9E2
endfb.dilutions = ( 1.e10, 10000.0, 3549.18335, 1259.67004, 447.079956, \
158.676849, 56.3173141 )
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Cm244"
endfb.mat = 9637
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-096_Cm_244.endf"
endfb.fission = 1 # fission without delayed neutrons
endfb.ss = (4.632489, 3.206464e5)
endfb.potential = 10.3200
endfb.eFiss = 1.9049E+02
endfb.dilutions = ( 1.e10, 10000.0, 3549.18335, 1259.67004, 447.079956, \
158.676849, 56.3173141 )
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Cm245"
endfb.mat = 9640
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-096_Cm_245.endf"
endfb.fission = 1 # fission without delayed neutrons
endfb.ss = (4.632489, 3.206464e5)
endfb.potential = 11.3900
endfb.eFiss = 1.9E2
endfb.dilutions = ( 1.e10, 10000.0, 3549.18335, 1259.67004, 447.079956, \
158.676849, 56.3173141 )
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Cm246"
endfb.mat = 9643
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-096_Cm_246.endf"
endfb.fission = 1 # fission without delayed neutrons
endfb.ss = (4.632489, 3.206464e5)
endfb.potential = 10.2758
endfb.eFiss = 1.9E2
endfb.dilutions = ( 1.e10, 10000.0, 3549.18335, 1259.67004, 447.079956, \
158.676849, 56.3173141 )
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Cm247"
endfb.mat = 9646
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-096_Cm_247.endf"
endfb.fission = 1 # fission without delayed neutrons
endfb.ss = (4.632489, 3.206464e5)
endfb.potential = 10.0000
endfb.eFiss = 1.9E2
endfb.dilutions = ( 1.e10, 10000.0, 3549.18335, 1259.67004, 447.079956, \
158.676849, 56.3173141 )
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Cm248"
endfb.mat = 9649
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-096_Cm_248.endf"
endfb.fission = 1 # fission without delayed neutrons
endfb.ss = (4.632489, 3.206464e5)
endfb.potential = 10.3971
endfb.eFiss = 1.8885E+02
endfb.dilutions = ( 1.e10, 10000.0, 3549.18335, 1259.67004, 447.079956, \
158.676849, 56.3173141 )
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Bk249"
endfb.mat = 9752
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-097_Bk_249.endf"
endfb.fission = 1 # fission without delayed neutrons
endfb.ss = (4.632489, 3.206464e5)
endfb.potential = 11.0983
endfb.eFiss = 1.9E2
endfb.dilutions = ( 1.e10, 10000.0, 3549.18335, 1259.67004, 447.079956, \
158.676849, 56.3173141 )
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Cf249"
endfb.mat = 9852
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-098_Cf_249.endf"
endfb.fission = 1 # fission without delayed neutrons
endfb.ss = (4.632489, 3.206464e5)
endfb.potential = 11.1510
endfb.eFiss = 1.9E2
endfb.dilutions = ( 1.e10, 10000.0, 3549.18335, 1259.67004, 447.079956, \
158.676849, 56.3173141 )
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Cf250"
endfb.mat = 9855
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-098_Cf_250.endf"
endfb.fission = 1 # fission without delayed neutrons
endfb.ss = (4.632489, 3.206464e5)
endfb.potential = 9.8800
endfb.eFiss = 1.9E2
endfb.dilutions = ( 1.e10, 10000.0, 3549.18335, 1259.67004, 447.079956, \
158.676849, 56.3173141 )
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Cf251"
endfb.mat = 9858
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-098_Cf_251.endf"
endfb.fission = 1 # fission without delayed neutrons
endfb.ss = (4.632489, 3.206464e5)
endfb.potential = 9.8800
endfb.eFiss = 1.9E2
endfb.dilutions = ( 1.e10, 10000.0, 3549.18335, 1259.67004, 447.079956, \
158.676849, 56.3173141 )
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Cf252"
endfb.mat = 9861
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-098_Cf_252.endf"
endfb.fission = 1 # fission without delayed neutrons
endfb.ss = (4.632489, 3.206464e5)
endfb.potential = 9.8000
endfb.eFiss = 1.9E2
endfb.dilutions = ( 1.e10, 10000.0, 3549.18335, 1259.67004, 447.079956, \
158.676849, 56.3173141 )
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Cf253"
endfb.mat = 9864
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-098_Cf_253.endf"
endfb.fission = 1 # fission without delayed neutrons
endfb.ss = (4.632489, 3.206464e5)
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
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-032_Ge_072.endf"
endfb.makeFp()

endfb.hmat = "Ge73"
endfb.mat = 3234
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-032_Ge_073.endf"
endfb.makeFp()

endfb.hmat = "Ge74"
endfb.mat = 3237
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-032_Ge_074.endf"
endfb.makeFp()

endfb.hmat = "Ge76"
endfb.mat = 3243
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-032_Ge_076.endf"
endfb.makeFp()

endfb.hmat = "As75"
endfb.mat = 3325
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-033_As_075.endf"
endfb.makeFp()

endfb.hmat = "Se76"
endfb.mat = 3431
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-034_Se_076.endf"
endfb.makeFp()

endfb.hmat = "Se77"
endfb.mat = 3434
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-034_Se_077.endf"
endfb.makeFp()

endfb.hmat = "Se78"
endfb.mat = 3437
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-034_Se_078.endf"
endfb.makeFp()

endfb.hmat = "Se79"
endfb.mat = 3440
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-034_Se_079.endf"
endfb.makeFp()

endfb.hmat = "Se80"
endfb.mat = 3443
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-034_Se_080.endf"
endfb.makeFp()

endfb.hmat = "Se82"
endfb.mat = 3449
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-034_Se_082.endf"
endfb.makeFp()

endfb.hmat = "Br79"
endfb.mat = 3525
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-035_Br_079.endf"
endfb.makeFp()

endfb.hmat = "Br81"
endfb.mat = 3531
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-035_Br_081.endf"
endfb.makeFp()

endfb.hmat = "Kr80"
endfb.mat = 3631
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-036_Kr_080.endf"
endfb.makeFp()

endfb.hmat = "Kr82"
endfb.mat = 3637
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-036_Kr_082.endf"
endfb.makeFp()

endfb.hmat = "Kr83"
endfb.mat = 3640
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-036_Kr_083.endf"
endfb.makeFp()

endfb.hmat = "Kr84"
endfb.mat = 3643
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-036_Kr_084.endf"
endfb.makeFp()

endfb.hmat = "Kr85"
endfb.mat = 3646
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-036_Kr_085.endf"
endfb.makeFp()

endfb.hmat = "Kr86"
endfb.mat = 3649
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-036_Kr_086.endf"
endfb.makeFp()

endfb.hmat = "Rb85"
endfb.mat = 3725
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-037_Rb_085.endf"
endfb.makeFp()

endfb.hmat = "Rb87"
endfb.mat = 3731
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-037_Rb_087.endf"
endfb.makeFp()

endfb.hmat = "Sr86"
endfb.mat = 3831
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-038_Sr_086.endf"
endfb.makeFp()

endfb.hmat = "Sr87"
endfb.mat = 3834
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-038_Sr_087.endf"
endfb.makeFp()

endfb.hmat = "Sr88"
endfb.mat = 3837
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-038_Sr_088.endf"
endfb.makeFp()

endfb.hmat = "Sr89"
endfb.mat = 3840
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-038_Sr_089.endf"
endfb.makeFp()

endfb.hmat = "Sr90"
endfb.mat = 3843
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-038_Sr_090.endf"
endfb.makeFp()

endfb.hmat = "Y89"
endfb.mat = 3925
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-039_Y_089.endf"
endfb.makeFp()

endfb.hmat = "Y90"
endfb.mat = 3928
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-039_Y_090.endf"
endfb.makeFp()

endfb.hmat = "Y91"
endfb.mat = 3931
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-039_Y_091.endf"
endfb.makeFp()

endfb.hmat = "Zr90"
endfb.mat = 4025
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-040_Zr_090.endf"
endfb.fission = None
endfb.ss = (4.632489, 1.858471e4)
endfb.potential = 6.8813
endfb.dilutions = ( 1.e10, 10000.0,  3866.97, 1495.35, 578.2475, 223.6068, 86.4682, 33.4370, 12.9300, 5.0 )
endfb.autolib = (4.632489, 3.481068e3, 0.0005)
endfb.pendf()
endfb.gendf()
endfb.draglib()
endfb.autolib = (4.632489, 1.11377e4, 0.0005)

endfb.hmat = "Zr91"
endfb.mat = 4028
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-040_Zr_091.endf"
endfb.fission = None
endfb.ss = (4.632489, 1.858471e4)
endfb.potential = 6.8813
endfb.dilutions = ( 1.e10, 10000.0, 3549.18335, 1259.67004, 447.079956, \
158.676849, 56.3173141, 19.9880447, 7.09412289, 2.51783395 )
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Zr92"
endfb.mat = 4031
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-040_Zr_092.endf"
endfb.fission = None
endfb.ss = (4.632489, 1.858471e4)
endfb.potential = 6.8813
endfb.dilutions = ( 1.e10, 10000.0, 3549.18335, 1259.67004, 447.079956, \
158.676849, 56.3173141, 19.9880447, 7.09412289, 2.51783395 )
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Zr93"
endfb.mat = 4034
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-040_Zr_093.endf"
endfb.makeFp()

endfb.hmat = "Zr94"
endfb.mat = 4037
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-040_Zr_094.endf"
endfb.fission = None
endfb.ss = (4.632489, 1.858471e4)
endfb.potential = 6.8813
endfb.dilutions = ( 1.e10, 10000.0, 3549.18335, 1259.67004, 447.079956, \
158.676849, 56.3173141, 19.9880447, 7.09412289, 2.51783395 )
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Zr95"
endfb.mat = 4040
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-040_Zr_095.endf"
endfb.makeFp()

endfb.hmat = "Zr96"
endfb.mat = 4043
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-040_Zr_096.endf"
endfb.fission = None
endfb.ss = (4.632489, 1.858471e4)
endfb.potential = 6.8813
endfb.dilutions = ( 1.e10, 10000.0, 3549.18335, 1259.67004, 447.079956, \
158.676849, 56.3173141, 19.9880447, 7.09412289, 2.51783395 )
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Nb93"
endfb.mat = 4125
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-041_Nb_093.endf"
endfb.makeFp()

endfb.hmat = "Nb94"
endfb.mat = 4128
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-041_Nb_094.endf"
endfb.makeFp()

endfb.hmat = "Nb95"
endfb.mat = 4131
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-041_Nb_095.endf"
endfb.makeFp()

endfb.hmat = "Mo92"
endfb.mat = 4225
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-042_Mo_092.endf"
endfb.makeFp()

endfb.hmat = "Mo94"
endfb.mat = 4231
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-042_Mo_094.endf"
endfb.makeFp()

endfb.hmat = "Mo95"
endfb.mat = 4234
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-042_Mo_095.endf"
endfb.makeFp()

endfb.hmat = "Mo96"
endfb.mat = 4237
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-042_Mo_096.endf"
endfb.makeFp()

endfb.hmat = "Mo97"
endfb.mat = 4240
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-042_Mo_097.endf"
endfb.makeFp()

endfb.hmat = "Mo98"
endfb.mat = 4243
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-042_Mo_098.endf"
endfb.makeFp()

endfb.hmat = "Mo99"
endfb.mat = 4246
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-042_Mo_099.endf"
endfb.makeFp()

endfb.hmat = "Mo100"
endfb.mat = 4249
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-042_Mo_100.endf"
endfb.makeFp()

endfb.hmat = "Tc99"
endfb.mat = 4325
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-043_Tc_099.endf"
endfb.ss = (4.632489, 1.858471e4)
endfb.potential = 6.4675
endfb.dilutions = ( 1.e10, 10000.0, 4216.96552, 1778.27959, 749.894278, \
316.227791, 133.352152,  56.2341357, 23.7137381, 10.0 )
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Ru99"
endfb.mat = 4434
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-044_Ru_099.endf"
endfb.makeFp()

endfb.hmat = "Ru100"
endfb.mat = 4437
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-044_Ru_100.endf"
endfb.makeFp()

endfb.hmat = "Ru101"
endfb.mat = 4440
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-044_Ru_101.endf"
endfb.makeFp()

endfb.hmat = "Ru102"
endfb.mat = 4443
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-044_Ru_102.endf"
endfb.makeFp()

endfb.hmat = "Ru103"
endfb.mat = 4446
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-044_Ru_103.endf"
endfb.makeFp()

endfb.hmat = "Ru104"
endfb.mat = 4449
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-044_Ru_104.endf"
endfb.makeFp()

endfb.hmat = "Ru105"
endfb.mat = 4452
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-044_Ru_105.endf"
endfb.makeFp()

endfb.hmat = "Ru106"
endfb.mat = 4455
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-044_Ru_106.endf"
endfb.makeFp()

endfb.hmat = "Rh103"
endfb.mat = 4525
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-045_Rh_103.endf"
endfb.makeFp()

endfb.hmat = "Rh105"
endfb.mat = 4531
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-045_Rh_105.endf"
endfb.makeFp()

endfb.hmat = "Pd104"
endfb.mat = 4631
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-046_Pd_104.endf"
endfb.makeFp()

endfb.hmat = "Pd105"
endfb.mat = 4634
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-046_Pd_105.endf"
endfb.makeFp()

endfb.hmat = "Pd106"
endfb.mat = 4637
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-046_Pd_106.endf"
endfb.makeFp()

endfb.hmat = "Pd107"
endfb.mat = 4640
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-046_Pd_107.endf"
endfb.makeFp()

endfb.hmat = "Pd108"
endfb.mat = 4643
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-046_Pd_108.endf"
endfb.makeFp()

endfb.hmat = "Pd110"
endfb.mat = 4649
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-046_Pd_110.endf"
endfb.makeFp()

endfb.hmat = "Ag107"
endfb.mat = 4725
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-047_Ag_107.endf"
endfb.ss = (4.632489, 1.858471e4)
endfb.potential = 5.4739
endfb.dilutions = ( 1.e10, 10000.0, 3549.18335, 1259.67004, 447.079956, \
158.676849, 56.3173141, 19.9880447, 7.09412289, 2.51783395 )
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Ag109"
endfb.mat = 4731
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-047_Ag_109.endf"
endfb.ss = (4.632489, 1.858471e4)
endfb.potential = 5.3316
endfb.dilutions = ( 1.e10, 10000.0, 3549.18335, 1259.67004, 447.079956, \
158.676849, 56.3173141, 19.9880447, 7.09412289, 2.51783395 )
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Ag111"
endfb.mat = 4737
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-047_Ag_111.endf"
endfb.makeFp()

endfb.hmat = "Cd106"
endfb.mat = 4825
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-048_Cd_106.endf"
endfb.makeFp()

endfb.hmat = "Cd108"
endfb.mat = 4831
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-048_Cd_108.endf"
endfb.makeFp()

endfb.hmat = "Cd110"
endfb.mat = 4837
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-048_Cd_110.endf"
jeff3p2.ss = (4.632489, 1.858471e4)
jeff3p2.potential = 5.1762
endfb.dilutions = ( 1.e10, 10000.0, 3549.18335, 1259.67004, 447.079956, \
158.676849, 56.3173141, 19.9880447, 7.09412289, 2.51783395 )
jeff3p2.pendf()
jeff3p2.gendf()
jeff3p2.draglib()

endfb.hmat = "Cd111"
endfb.mat = 4840
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-048_Cd_111.endf"
endfb.makeFp()

endfb.hmat = "Cd112"
endfb.mat = 4843
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-048_Cd_112.endf"
endfb.makeFp()

endfb.hmat = "Cd113"
endfb.mat = 4846
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-048_Cd_113.endf"
endfb.makeFp()

endfb.hmat = "Cd114"
endfb.mat = 4849
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-048_Cd_114.endf"
endfb.branchingNG = 0.079383
endfb.makeFp()
endfb.branchingNG = None

endfb.hmat = "Cd115m"
endfb.mat = 4853
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-048_Cd_115m1.endf"
endfb.makeFp()

endfb.hmat = "Cd116"
endfb.mat = 4855
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-048_Cd_116.endf"
endfb.makeFp()

endfb.hmat = "In113"
endfb.mat = 4925
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-049_In_113.endf"
endfb.makeFp()

endfb.hmat = "In115"
endfb.mat = 4931
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-049_In_115.endf"
endfb.ss = (4.632489, 1.858471e4)
endfb.potential = 5.0439
endfb.dilutions = ( 1.e10, 10000.0, 3549.18335, 1259.67004, 447.079956, \
158.676849, 56.3173141, 19.9880447, 7.09412289, 2.51783395 )
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Sn112"
endfb.mat = 5025
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-050_Sn_112.endf"
endfb.makeFp()

endfb.hmat = "Sn114"
endfb.mat = 5031
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-050_Sn_114.endf"
endfb.makeFp()

endfb.hmat = "Sn115"
endfb.mat = 5034
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-050_Sn_115.endf"
endfb.makeFp()

endfb.hmat = "Sn116"
endfb.mat = 5037
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-050_Sn_116.endf"
endfb.makeFp()

endfb.hmat = "Sn117"
endfb.mat = 5040
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-050_Sn_117.endf"
endfb.makeFp()

endfb.hmat = "Sn118"
endfb.mat = 5043
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-050_Sn_118.endf"
endfb.makeFp()

endfb.hmat = "Sn119"
endfb.mat = 5046
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-050_Sn_119.endf"
endfb.makeFp()

endfb.hmat = "Sn120"
endfb.mat = 5049
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-050_Sn_120.endf"
endfb.makeFp()

endfb.hmat = "Sn122"
endfb.mat = 5055
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-050_Sn_122.endf"
endfb.makeFp()

endfb.hmat = "Sn123"
endfb.mat = 5058
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-050_Sn_123.endf"
endfb.makeFp()

endfb.hmat = "Sn124"
endfb.mat = 5061
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-050_Sn_124.endf"
endfb.makeFp()

endfb.hmat = "Sn125"
endfb.mat = 5064
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-050_Sn_125.endf"
endfb.makeFp()

endfb.hmat = "Sn126"
endfb.mat = 5067
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-050_Sn_126.endf"
endfb.makeFp()

endfb.hmat = "Sb121"
endfb.mat = 5125
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-051_Sb_121.endf"
endfb.makeFp()

endfb.hmat = "Sb123"
endfb.mat = 5131
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-051_Sb_123.endf"
endfb.makeFp()

endfb.hmat = "Sb124"
endfb.mat = 5134
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-051_Sb_124.endf"
endfb.makeFp()

endfb.hmat = "Sb125"
endfb.mat = 5137
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-051_Sb_125.endf"
endfb.makeFp()

endfb.hmat = "Sb126"
endfb.mat = 5140
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-051_Sb_126.endf"
endfb.makeFp()

endfb.hmat = "Te122"
endfb.mat = 5231
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-052_Te_122.endf"
endfb.makeFp()

endfb.hmat = "Te123"
endfb.mat = 5234
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-052_Te_123.endf"
endfb.makeFp()

endfb.hmat = "Te124"
endfb.mat = 5237
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-052_Te_124.endf"
endfb.makeFp()

endfb.hmat = "Te125"
endfb.mat = 5240
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-052_Te_125.endf"
endfb.makeFp()

endfb.hmat = "Te126"
endfb.mat = 5243
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-052_Te_126.endf"
endfb.branchingNG = 0.091528
endfb.makeFp()
endfb.branchingNG = None

endfb.hmat = "Te127m"
endfb.mat = 5247
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-052_Te_127m1.endf"
endfb.makeFp()

endfb.hmat = "Te128"
endfb.mat = 5249
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-052_Te_128.endf"
endfb.branchingNG = 0.031894
endfb.makeFp()
endfb.branchingNG = None

endfb.hmat = "Te129m"
endfb.mat = 5253
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-052_Te_129m1.endf"
endfb.makeFp()

endfb.hmat = "Te130"
endfb.mat = 5255
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-052_Te_130.endf"
endfb.makeFp()

endfb.hmat = "Te132"
endfb.mat = 5261
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-052_Te_132.endf"
endfb.makeFp()

endfb.hmat = "I127"
endfb.mat = 5325
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-053_I_127.endf"
endfb.ss = (4.632489, 1.858471e4)
endfb.potential = 4.5239
endfb.dilutions = ( 1.e10, 10000.0, 4216.96552, 1778.27959, 749.894278, \
316.227791, 133.352152,  56.2341357, 23.7137381, 10.0 )
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "I129"
endfb.mat = 5331
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-053_I_129.endf"
endfb.ss = (4.632489, 1.858471e4)
endfb.potential = 5.8221
endfb.dilutions = ( 1.e10, 10000.0, 4216.96552, 1778.27959, 749.894278, \
316.227791, 133.352152,  56.2341357, 23.7137381, 10.0 )
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "I130"
endfb.mat = 5334
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-053_I_130.endf"
endfb.makeFp()

endfb.hmat = "I131"
endfb.mat = 5337
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-053_I_131.endf"
endfb.makeFp()

endfb.hmat = "I135"
endfb.mat = 5349
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-053_I_135.endf"
endfb.makeFp()

endfb.hmat = "Xe128"
endfb.mat = 5437
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-054_Xe_128.endf"
endfb.makeFp()

endfb.hmat = "Xe129"
endfb.mat = 5440
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-054_Xe_129.endf"
endfb.makeFp()

endfb.hmat = "Xe130"
endfb.mat = 5443
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-054_Xe_130.endf"
endfb.makeFp()

endfb.hmat = "Xe131"
endfb.mat = 5446
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-054_Xe_131.endf"
endfb.makeFp()

endfb.hmat = "Xe132"
endfb.mat = 5449
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-054_Xe_132.endf"
endfb.makeFp()

endfb.hmat = "Xe133"
endfb.mat = 5452
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-054_Xe_133.endf"
endfb.makeFp()

endfb.hmat = "Xe134"
endfb.mat = 5455
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-054_Xe_134.endf"
endfb.makeFp()

endfb.hmat = "Xe135"
endfb.mat = 5458
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-054_Xe_135.endf"
endfb.makeFp()

endfb.hmat = "Xe136"
endfb.mat = 5461
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-054_Xe_136.endf"
endfb.makeFp()

endfb.hmat = "Cs133"
endfb.mat = 5525
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-055_Cs_133.endf"
endfb.makeFp()

endfb.hmat = "Cs134"
endfb.mat = 5528
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-055_Cs_134.endf"
endfb.makeFp()

endfb.hmat = "Cs135"
endfb.mat = 5531
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-055_Cs_135.endf"
endfb.makeFp()

endfb.hmat = "Cs136"
endfb.mat = 5534
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-055_Cs_136.endf"
endfb.makeFp()

endfb.hmat = "Cs137"
endfb.mat = 5537
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-055_Cs_137.endf"
endfb.makeFp()

endfb.hmat = "Ba134"
endfb.mat = 5637
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-056_Ba_134.endf"
endfb.makeFp()

endfb.hmat = "Ba135"
endfb.mat = 5640
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-056_Ba_135.endf"
endfb.makeFp()

endfb.hmat = "Ba136"
endfb.mat = 5643
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-056_Ba_136.endf"
endfb.makeFp()

endfb.hmat = "Ba137"
endfb.mat = 5646
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-056_Ba_137.endf"
endfb.makeFp()

endfb.hmat = "Ba138"
endfb.mat = 5649
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-056_Ba_138.endf"
endfb.makeFp()

endfb.hmat = "Ba140"
endfb.mat = 5655
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-056_Ba_140.endf"
endfb.makeFp()

endfb.hmat = "La138"
endfb.mat = 5725
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-057_La_138.endf"
endfb.makeFp()

endfb.hmat = "La139"
endfb.mat = 5728
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-057_La_139.endf"
endfb.makeFp()

endfb.hmat = "La140"
endfb.mat = 5731
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-057_La_140.endf"
endfb.makeFp()

endfb.hmat = "Ce140"
endfb.mat = 5837
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-058_Ce_140.endf"
endfb.makeFp()

endfb.hmat = "Ce141"
endfb.mat = 5840
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-058_Ce_141.endf"
endfb.makeFp()

endfb.hmat = "Ce142"
endfb.mat = 5843
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-058_Ce_142.endf"
endfb.makeFp()

endfb.hmat = "Ce143"
endfb.mat = 5846
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-058_Ce_143.endf"
endfb.makeFp()

endfb.hmat = "Ce144"
endfb.mat = 5849
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-058_Ce_144.endf"
endfb.makeFp()

endfb.hmat = "Pr141"
endfb.mat = 5925
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-059_Pr_141.endf"
endfb.makeFp()

endfb.hmat = "Pr142"
endfb.mat = 5928
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-059_Pr_142.endf"
endfb.makeFp()

endfb.hmat = "Pr143"
endfb.mat = 5931
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-059_Pr_143.endf"
endfb.makeFp()

endfb.hmat = "Nd142"
endfb.mat = 6025
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-060_Nd_142.endf"
endfb.makeFp()

endfb.hmat = "Nd143"
endfb.mat = 6028
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-060_Nd_143.endf"
endfb.makeFp()

endfb.hmat = "Nd144"
endfb.mat = 6031
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-060_Nd_144.endf"
endfb.makeFp()

endfb.hmat = "Nd145"
endfb.mat = 6034
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-060_Nd_145.endf"
endfb.makeFp()

endfb.hmat = "Nd146"
endfb.mat = 6037
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-060_Nd_146.endf"
endfb.makeFp()

endfb.hmat = "Nd147"
endfb.mat = 6040
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-060_Nd_147.endf"
endfb.makeFp()

endfb.hmat = "Nd148"
endfb.mat = 6043
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-060_Nd_148.endf"
endfb.makeFp()

endfb.hmat = "Nd150"
endfb.mat = 6049
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-060_Nd_150.endf"
endfb.makeFp()

endfb.hmat = "Pm147"
endfb.mat = 6149
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-061_Pm_147.endf"
endfb.branchingNG = 0.470
endfb.makeFp()
endfb.branchingNG = None

endfb.hmat = "Pm148"
endfb.mat = 6152
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-061_Pm_148.endf"
endfb.makeFp()

endfb.hmat = "Pm148m"
endfb.mat = 6153
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-061_Pm_148m1.endf"
endfb.makeFp()

endfb.hmat = "Pm149"
endfb.mat = 6155
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-061_Pm_149.endf"
endfb.makeFp()

endfb.hmat = "Pm151"
endfb.mat = 6161
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-061_Pm_151.endf"
endfb.makeFp()

endfb.hmat = "Sm147"
endfb.mat = 6234
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-062_Sm_147.endf"
endfb.makeFp()

endfb.hmat = "Sm148"
endfb.mat = 6237
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-062_Sm_148.endf"
endfb.makeFp()

endfb.hmat = "Sm149"
endfb.mat = 6240
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-062_Sm_149.endf"
endfb.makeFp()

endfb.hmat = "Sm150"
endfb.mat = 6243
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-062_Sm_150.endf"
endfb.makeFp()

endfb.hmat = "Sm151"
endfb.mat = 6246
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-062_Sm_151.endf"
endfb.makeFp()

endfb.hmat = "Sm152"
endfb.mat = 6249
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-062_Sm_152.endf"
endfb.makeFp()

endfb.hmat = "Sm153"
endfb.mat = 6252
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-062_Sm_153.endf"
endfb.makeFp()

endfb.hmat = "Sm154"
endfb.mat = 6255
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-062_Sm_154.endf"
endfb.makeFp()

endfb.hmat = "Eu151"
endfb.mat = 6325
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-063_Eu_151.endf"
endfb.makeFp()

endfb.hmat = "Eu152"
endfb.mat = 6328
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-063_Eu_152.endf"
endfb.makeFp()

endfb.hmat = "Eu153"
endfb.mat = 6331
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-063_Eu_153.endf"
endfb.makeFp()

endfb.hmat = "Eu154"
endfb.mat = 6334
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-063_Eu_154.endf"
endfb.makeFp()

endfb.hmat = "Eu155"
endfb.mat = 6337
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-063_Eu_155.endf"
endfb.makeFp()

endfb.hmat = "Eu156"
endfb.mat = 6340
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-063_Eu_156.endf"
endfb.makeFp()

endfb.hmat = "Eu157"
endfb.mat = 6343
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-063_Eu_157.endf"
endfb.makeFp()

endfb.hmat = "Gd152"
endfb.mat = 6425
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-064_Gd_152.endf"
endfb.ss = (4.632489, 1.858471e4)
endfb.potential = 8.0425
endfb.dilutions = ( 1.e10, 10000.0, 4216.96552, 1778.27959, 749.894278, \
316.227791, 133.352152,  56.2341357, 23.7137381, 10.0 )
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Gd154"
endfb.mat = 6431
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-064_Gd_154.endf"
endfb.ss = (4.632489, 1.858471e4)
endfb.potential = 6.6723
endfb.dilutions = ( 1.e10, 10000.0, 4216.96552, 1778.27959, 749.894278, \
316.227791, 133.352152,  56.2341357, 23.7137381, 10.0 )
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Gd155"
endfb.mat = 6434
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-064_Gd_155.endf"
endfb.ss = (4.632489, 1.858471e4)
endfb.potential = 6.3376
endfb.dilutions = ( 1.e10, 10000.0, 4216.96552, 1778.27959, 749.894278, \
316.227791, 133.352152,  56.2341357, 23.7137381, 10.0 )
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Gd156"
endfb.mat = 6437
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-064_Gd_156.endf"
endfb.ss = (4.632489, 1.858471e4)
endfb.potential = 7.3792
endfb.dilutions = ( 1.e10, 10000.0, 4216.96552, 1778.27959, 749.894278, \
316.227791, 133.352152,  56.2341357, 23.7137381, 10.0 )
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Gd157"
endfb.mat = 6440
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-064_Gd_157.endf"
endfb.ss = (4.632489, 1.858471e4)
endfb.potential = 6.6063
endfb.dilutions = ( 1.e10, 10000.0, 4216.96552, 1778.27959, 749.894278, \
316.227791, 133.352152,  56.2341357, 23.7137381, 10.0 )
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Gd158"
endfb.mat = 6443
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-064_Gd_158.endf"
endfb.ss = (4.632489, 1.858471e4)
endfb.potential = 7.6454
endfb.dilutions = ( 1.e10, 10000.0, 4216.96552, 1778.27959, 749.894278, \
316.227791, 133.352152,  56.2341357, 23.7137381, 10.0 )
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Gd160"
endfb.mat = 6449
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-064_Gd_160.endf"
endfb.ss = (4.632489, 1.858471e4)
endfb.potential = 7.0241
endfb.dilutions = ( 1.e10, 10000.0, 4216.96552, 1778.27959, 749.894278, \
316.227791, 133.352152,  56.2341357, 23.7137381, 10.0 )
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.dilutions = None

endfb.hmat = "Tb159"
endfb.mat = 6525
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-065_Tb_159.endf"
endfb.makeFp()

endfb.hmat = "Tb160"
endfb.mat = 6528
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-065_Tb_160.endf"
endfb.makeFp()

endfb.hmat = "Dy160"
endfb.mat = 6637
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-066_Dy_160.endf"
endfb.ss = (4.632489, 1.858471e4)
endfb.potential = 6.9861
endfb.dilutions = ( 1.e10, 10000.0, 4216.96552, 1778.27959, 749.894278, \
316.227791, 133.352152,  56.2341357, 23.7137381, 10.0 )
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Dy161"
endfb.mat = 6640
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-066_Dy_161.endf"
endfb.ss = (4.632489, 1.858471e4)
endfb.potential = 7.0121
endfb.dilutions = ( 1.e10, 10000.0, 4216.96552, 1778.27959, 749.894278, \
316.227791, 133.352152,  56.2341357, 23.7137381, 10.0 )
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Dy162"
endfb.mat = 6643
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-066_Dy_162.endf"
endfb.ss = (4.632489, 1.858471e4)
endfb.potential = 4.5681
endfb.dilutions = ( 1.e10, 10000.0, 4216.96552, 1778.27959, 749.894278, \
316.227791, 133.352152,  56.2341357, 23.7137381, 10.0 )
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Dy163"
endfb.mat = 6646
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-066_Dy_163.endf"
endfb.ss = (4.632489, 1.858471e4)
endfb.potential = 7.0639
endfb.dilutions = ( 1.e10, 10000.0, 4216.96552, 1778.27959, 749.894278, \
316.227791, 133.352152,  56.2341357, 23.7137381, 10.0 )
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Dy164"
endfb.mat = 6649
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-066_Dy_164.endf"
endfb.ss = (4.632489, 1.858471e4)
endfb.potential = 7.0897
endfb.dilutions = ( 1.e10, 10000.0, 4216.96552, 1778.27959, 749.894278, \
316.227791, 133.352152,  56.2341357, 23.7137381, 10.0 )
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Hf174"
endfb.mat = 7225
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-072_Hf_174.endf"
endfb.ss = (4.632489, 1.858471e4)
endfb.potential = 7.1385
endfb.dilutions = ( 1.e10, 10000.0, 4216.96552, 1778.27959, 749.894278, \
316.227791, 133.352152,  56.2341357, 23.7137381, 10.0 )
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Hf176"
endfb.mat = 7231
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-072_Hf_176.endf"
endfb.ss = (4.632489, 1.858471e4)
endfb.potential = 7.1935
endfb.dilutions = ( 1.e10, 10000.0, 4216.96552, 1778.27959, 749.894278, \
316.227791, 133.352152,  56.2341357, 23.7137381, 10.0 )
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Hf177"
endfb.mat = 7234
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-072_Hf_177.endf"
endfb.ss = (4.632489, 1.858471e4)
endfb.potential = 7.2202
endfb.dilutions = ( 1.e10, 10000.0, 4216.96552, 1778.27959, 749.894278, \
316.227791, 133.352152,  56.2341357, 23.7137381, 10.0 )
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Hf178"
endfb.mat = 7237
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-072_Hf_178.endf"
endfb.ss = (4.632489, 1.858471e4)
endfb.potential = 7.2469
endfb.dilutions = ( 1.e10, 10000.0, 4216.96552, 1778.27959, 749.894278, \
316.227791, 133.352152,  56.2341357, 23.7137381, 10.0 )
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Hf179"
endfb.mat = 7240
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-072_Hf_179.endf"
endfb.ss = (4.632489, 1.858471e4)
endfb.potential = 7.2736
endfb.dilutions = ( 1.e10, 10000.0, 4216.96552, 1778.27959, 749.894278, \
316.227791, 133.352152,  56.2341357, 23.7137381, 10.0 )
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.hmat = "Hf180"
endfb.mat = 7243
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-072_Hf_180.endf"
endfb.ss = (4.632489, 1.858471e4)
endfb.potential = 7.3004
endfb.dilutions = ( 1.e10, 10000.0, 4216.96552, 1778.27959, 749.894278, \
316.227791, 133.352152,  56.2341357, 23.7137381, 10.0 )
endfb.pendf()
endfb.gendf()
endfb.draglib()

endfb.dilutions = None

endfb.hmat = "Ho165"
endfb.mat = 6725
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-067_Ho_165.endf"
endfb.makeFp()

endfb.hmat = "Er166"
endfb.mat = 6837
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-068_Er_166.endf"
endfb.makeFp()

endfb.hmat = "Er167"
endfb.mat = 6840
endfb.evaluationFile = "$HOME/evaluations/ENDFB7r1/neutrons/n-068_Er_167.endf"
endfb.makeFp()

# Process the burnup chain:

endfb.fissionFile = "$HOME/evaluations/ENDFB7r1/nfy/"
endfb.decayFile = "$HOME/evaluations/ENDFB7r1/decay/"
endfb.burnup()
