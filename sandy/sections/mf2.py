"""
This module contains only two public functions:

    * `read_mf2`
    * `write_mf2`

Function `read` reads a MF2/MT section from a string and produces a content
object with a dictionary-like structure.
The content object can be accessed using most of the keywords specified in
the ENDF6 manual for this specific MF section.

Function `write_mf2` writes a content object for a MF2/MT section into a
string.
MAT, MF, MT and line numbers are also added (each line ends with a `\n`).
"""

import sandy


__author__ = "Rayan HADDAD"
__all__ = [
        "read_mf2",
        "write_mf2"
        ]


mf = 2
mt = 151

def read_mf2(tape, mat):
    """
    Write MT section for MF2

    Parameters
    ----------
    tape : `sandy.Endf6`
        endf6 object containing requested section
    mat : `int`
        MAT number

    Returns
    -------
    `dict`
        Content of the ENDF-6 tape structured as nested `dict`.
    Examples
    --------
    Resonance parameters of the Thorium 233 
    LRU = 0 
    >>> tape = sandy.get_endf6_file("jeff_33", "xs", 902330)
    >>> sandy.read_mf2(tape, 9043) 
    {'MAT': 9043,
     'MF': 2,
     'MT': 151,
     'INTRO': {'ZA': 90233.0, 'AWR': 231.04, 'NIS': 1},
     'NIS': {'ISO': {90233: {'ZAI': 90233.0,
        'ABN': 1.0,
        'LFW': 0,
        'NER': {(1e-05, 1.9): {'EL': 1e-05,
          'EH': 1.9,
          'LRU': 0,
          'LRF': 0,
          'NRO': 0,
          'NAPS': 0,
          'SPI': 0.5,
          'AP': 0.9765,
          'NLS': 0}}}}}}
    
    
    Resonance parameters of the Thorium 230 
    LRU = 1 LRF = 2
    >>> tape = sandy.get_endf6_file("jeff_33", "xs", 902300)
    >>> dic = sandy.read_mf2(tape, 9034)
    >>> print (dic["NIS"]["ISO"][90230]["NER"][(1e-05, 251.0)]["L"][0]["ER"])
         [-1.427, 1.427, 17.27, 23.84, 32.2, 39.8, 48.1, 64.5, 75.6, 83.3, 103.0, 116.1, 133.8, 
          139.0, 148.2, 171.4, 184.2, 195.0, 209.0, 226.0, 241.0, 248.0]
        
    Resonance parameters of the iron 58
    LRU = 1 LRF = 3 
    >>> tape = sandy.get_endf6_file("endfb_80", "xs", 260580)
    >>> dic = sandy.read_mf2(tape, 2637)
    >>> print (dic["NIS"]["ISO"][26058]["NER"][(1e-05, 350000.0)]["L"][0]["AJ"])
     [0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 
      0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5]
     
    
    Resonance parameters of the iron 58
    LRU = 2 LFW = 0 LRF = 1
    >>> tape = sandy.get_endf6_file("endfb_80", "xs", 260580)
    >>> dic = sandy.read_mf2(tape, 2637)
    >>> print (dic["NIS"]["ISO"][26058]["NER"][(350000.0, 3000000.0)]["L"][3])
    {'AWRI': 57.43561, 'L': 3, 'NJS': 2, 'D': [8033.33, 6025.0], 'AJ': [2.5, 3.5],
     'AMUN': [1.0, 1.0], 'GNO': [0.482, 0.3615], 'GG': [0.7, 0.7]}
    
    
    Resonance parameters of the iron 54
    LRU = 1 LRF = 7
    >>> tape = sandy.get_endf6_file("endfb_80", "xs", 260540)
    >>> dic = sandy.read_mf2(tape, 2625)
    >>> print (dic["NIS"]["ISO"][26054]["NER"][(1e-05, 1036000.0)]["J"][2]["ER"])
    [3099.0, 0.0, 13581.0, 0.0, 19278.0, 0.0, 23030.0, 0.0, 28220.0, 0.0, 35260.0, 0.0, 38440.0, 0.0, 41200.0, 0.0, 
     50140.0, 0.0, 59200.0, 0.0, 68760.0, 0.0, 75800.0, 0.0, 77200.0, 0.0, 81280.0, 0.0, 83155.0, 0.0, 87330.0, 0.0, 
     91355.0, 0.0, 97770.0, 0.0, 99840.0, 0.0, 101710.0, 0.0, 112580.0, 0.0, 112960.0, 0.0, 115850.0, 0.0, 119780.0, 0.0, 
     126530.0, 0.0, 137930.0, 0.0, 141010.0, 0.0, 145510.0, 0.0, 153190.0, 0.0, 157200.0, 0.0, 164510.0, 0.0, 189090.0, 0.0, 
     194480.0, 0.0, 203970.0, 0.0, 207600.0, 0.0, 210010.0, 0.0, 225710.0, 0.0, 231220.0, 0.0, 233710.0, 0.0, 242060.0, 0.0, 
     245290.0, 0.0, 254890.0, 0.0, 263020.0, 0.0, 282940.0, 0.0, 309140.0, 0.0, 322020.0, 0.0, 332690.0, 0.0, 345270.0, 0.0, 
     357090.0, 0.0, 361090.0, 0.0, 384840.0, 0.0, 389720.0, 0.0, 389910.0, 0.0, 397910.0, 0.0, 398800.0, 0.0, 412390.0, 0.0, 
     427530.0, 0.0, 427700.0, 0.0, 446970.0, 0.0, 458970.0, 0.0, 461460.0, 0.0, 465450.0, 0.0, 470330.0, 0.0, 478690.0, 0.0, 
     490150.0, 0.0, 498060.0, 0.0, 511379.9, 0.0, 524220.0, 0.0, 530528.3, 0.0, 552571.6, 0.0, 559890.8, 0.0, 565199.3, 0.0, 
     568736.5, 0.0, 580930.9, 0.0, 595230.0, 0.0, 604249.2, 0.0, 604484.3, 0.0, 621213.9, 0.0, 638550.0, 0.0, 647068.3, 0.0, 
     656720.0, 0.0, 669905.7, 0.0, 676842.2, 0.0, 680048.7, 0.0, 680575.4, 0.0, 681811.0, 0.0, 683546.4, 0.0, 684087.7, 0.0, 
     687628.4, 0.0, 689409.5, 0.0, 697832.6, 0.0, 708779.0, 0.0, 722606.7, 0.0, 726551.9, 0.0, 738710.0, 0.0, 749476.6, 0.0, 
     773557.5, 0.0, 790978.6, 0.0, 795409.0, 0.0, 807088.2, 0.0, 812907.7, 0.0, 818600.5, 0.0, 825406.2, 0.0, 829051.9, 0.0, 
     835768.2, 0.0, 841278.0, 0.0, 843465.7, 0.0, 847479.5, 0.0, 852237.1, 0.0, 859129.1, 0.0, 876376.3, 0.0, 888111.6, 0.0, 
     896010.1, 0.0, 897719.0, 0.0, 909989.2, 0.0, 914450.8, 0.0, 916498.1, 0.0, 918795.3, 0.0, 922752.7, 0.0, 927389.0, 0.0, 
     944893.5, 0.0, 963621.3, 0.0, 968002.0, 0.0, 978423.1, 0.0, 983843.8, 0.0, 990014.5, 0.0, 1002368.0, 0.0, 1005913.0, 0.0, 
     1010327.0, 0.0, 1012920.0, 0.0, 1022815.0, 0.0, 1026431.0, 0.0, 1028429.0, 0.0, 1031001.0, 0.0, 1031465.0, 0.0, 1032329.0, 
     0.0, 1051001.0, 0.0, 1056050.0, 0.0, 1058038.0, 0.0, 1062925.0, 0.0, 1072620.0, 0.0, 1083488.0, 0.0, 1085559.0, 0.0, 
     1091677.0, 0.0, 1095462.0, 0.0, 1099070.0, 0.0, 1101233.0, 0.0, 1104403.0, 0.0, 1106613.0, 0.0, 1107490.0, 0.0, 1109121.0, 
     0.0, 1110144.0, 0.0, 1117903.0, 0.0, 1118632.0, 0.0, 1119608.0, 0.0, 1125225.0, 0.0, 1128490.0, 0.0, 1135245.0, 0.0, 
     1141235.0, 0.0, 1144014.0, 0.0, 1160850.0, 0.0, 1166727.0, 0.0, 1172094.0, 0.0, 1179700.0, 0.0, 1190647.0, 0.0, 1191495.0, 
     0.0, 1195039.0, 0.0, 1197252.0, 0.0, 1200500.0, 0.0, 1214900.0, 0.0, 1220700.0, 0.0, 1274800.0, 0.0, 1288100.0, 0.0, 
     1297300.0, 0.0, 1300387.0, 0.0, 1306103.0, 0.0, 1307309.0, 0.0, 1312399.0, 0.0, 1318355.0, 0.0, 1324614.0, 0.0, 1362008.0, 
     0.0, 1370898.0, 0.0, 1382499.0, 0.0, 1396057.0, 0.0, 1403500.0, 0.0, 1407800.0, 0.0, 1421700.0, 0.0, 1430000.0, 0.0, 
     1432400.0, 0.0, 1439300.0, 0.0, 1441300.0, 0.0, 1465800.0, 0.0, 1467900.0, 0.0, 1468900.0, 0.0, 1475500.0, 0.0, 1496700.0, 0.0]
    
    
    Resonance parameters of the Uranium 235
    LRU = 2 LRF = 2
    >>> tape = sandy.get_endf6_file("jeff_33", "xs",922350)
    >>> dic = sandy.read_mf2(tape,9228)
    >>> print (dic["NIS"]["ISO"][92235]["NER"][(2250.0, 46200.0)]['L'][0]["J"][0])
    {'AJ': 3.0, 'INT': 2, 'NE': 14, 'AMUX': 1.0, 'AMUN': 1.0, 'AMUG': 0.0, 'AMUF': 1.0, 
     'ES': [2250.0, 3000.0, 4000.0, 5000.0, 6000.0, 7000.0, 8000.001, 9000.0, 10000.0, 13000.0, 15000.0, 20000.0, 25000.0, 46200.0], 
     'D': [1.058, 1.0558, 1.0536, 1.0515, 1.0493, 1.0471, 1.0449, 1.0428, 1.0406, 1.0341, 1.0298, 1.0192, 1.0087, 0.96532], 
     'GX': [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 6.3286e-05, 0.0004222, 0.00093359, 0.0038198], 
     'GNO': [0.000107789, 0.0001072693, 0.0001067929, 0.0001063697, 0.0001059583, 0.0001055686, 0.0001051901, 
             0.0001048431, 0.0001044866, 0.0001035031, 0.000102875, 0.0001014063, 0.0001000348, 9.486007e-05], 
     'GG': [0.038513, 0.038529, 0.038545, 0.038561, 0.038577, 0.038593, 0.038609, 0.038626, 0.038641, 0.03869, 
            0.038722, 0.038803, 0.038884, 0.03921], 
     'GF': [0.40102, 0.40099, 0.40096, 0.40093, 0.40091, 0.40088, 0.40085, 0.40082, 0.4008, 0.40072, 0.40067, 
            0.40055, 0.40045, 0.40018]}
    """
    df = tape._get_section_df(mat, mf, mt)
    out = {
            "MAT": mat,         
            "MF": mf,
            "MT": mt,
            }
    i = 0
    I = {}
    C, i = sandy.read_cont(df, i)
    add = {
            "ZA": C.C1,             #designation for an isotope
            "AWR": C.C2,            #AWR is defines as the ratio of the mass of the material to that of the neutron
            "NIS": C.N1,            #NIS is the number of isotopes in the material
            }
    NIS = int(C.N1)
    I.update(add)
    NER = int(C.N1)
    LFW = int(C.L2)
    out.update({"INTRO": I})
    P = {}
    for l in range(NIS): 
        M = {}
        NER1 = {}
        ISO = {}
        C, i = sandy.read_cont(df, i)
        header1 = {
                "ZAI": C.C1,            #designation for an isotope
                "ABN": C.C2,            #Abundance of an isotope in the material 
                "LFW": C.L2,            #indication whether average fission wifths are given in the unresolbed resonance region
                } 
        NER = int(C.N1)
        LFW = int(C.L2)
        ZAI = int(C.C1)
        M.update(header1)
        dico = {}
        for j in range(NER):
            info = {}
            C, i = sandy.read_cont(df, i)
            header2 = {
                    "EL": C.C1,             #Lower limit for an energy range
                    "EH": C.C2,             #Upper limit for an energy range
                    "LRU": C.L1,            #Flag indicating whether this energy range contains data for resolved or unresolved resonance parameters:
                    "LRF": C.L2,            #Flag indicating which representation has been used for the energy range. 
                    "NRO": C.N1,            #Flag designating possible energy dependence of the scattering radiu
                    "NAPS":C.N2,            #Flag controlling the use of the two radii
                    }
            EL = C.C1
            EH = C.C2
            LRU = int(C.L1)
            LRF = int(C.L2)
            NRO = int(C.N1)
            if LRU == 0:
                LRU0 = {}
                LRU0.update(header2)
                C, i = sandy.read_cont(df, i)
                add = {
                        "SPI": C.C1,        #Spin, I, of the target nucleus.
                        "AP": C.C2,         #Scattering radius in units of 10e-12cm.
                        "NLS": C.N1,        #Number of l-values (neutron orbital angular momentum) in this energy region.
                        }
                LRU0.update(add)
                dico[(EL, EH)] = LRU0
            elif LRU == 1:
               if LRF == 1 or LRF == 2 :
                    if NRO == 0 :
                        LRU1_LRF1_2_NRO0 = {}
                        LRU1_LRF1_2_NRO0.update(header2)
                        C, i = sandy.read_cont(df, i)
                        add = {
                                "SPI": C.C1,
                                "AP": C.C2,
                                "NLS": C.N1,
                                }
                        NLS = int(C.N1)
                        LRU1_LRF1_2_NRO0.update(add)
                        LRU1_LRF1_2_NRO0_NLS = {}
                        for k in range(NLS) :
                            L, i = sandy.read_list(df, i)
                            add = {
                                    "AWRI": L.C1,       #Ratio of the mass of a particular isotope to that of a neutron
                                    "QX": L.C2,         #Q-value to be added to the incident particle’s center-of-mass energy to determine the channel energy for use in the penetrability factor.
                                    "L": L.L1,          #Value of l.
                                    "LRX": L.L2,        #Flag indicating whether this energy range contains a competitive width
                                    "NRS": L.N2,        #Number of resolved resonances for a given l-value.
                                    }
                            
                            LIST ={}
                            add_2 =  {
                                       "ER": L.B[0::6],     #Resonance energy
                                       "AJ": L.B[1::6],     #The absolute value of AJ is the floating-point value of J
                                       "GT": L.B[2::6],     #Resonance total width
                                       "GN": L.B[3::6],     #Neutron width 
                                       "GG": L.B[4::6],     #Radiation width
                                       "GF": L.B[5::6],     #Fission width
                                      }
                            LIST.update(add_2)
                            add.update(LIST)
                            LRU1_LRF1_2_NRO0_NLS.update({k:add})
                            LRU1_LRF1_2_NRO0.update({"L":LRU1_LRF1_2_NRO0_NLS})
                        dico[(EL, EH)] = LRU1_LRF1_2_NRO0  
                    else :   
                        LRU1_LRF1_2_NRO1 = {}
                        LRU1_LRF1_2_NRO1.update(header2)
                        T, i = sandy.read_tab1(df, i)
                        add = {
                                "NR": T.NBT,    
                                "NP": T.INT,
                                "E_int": T.x,
                                "AP": T.y,
                                }
                        LRU1_LRF1_2_NRO1.update(add)
                        C, i = sandy.read_cont(df, i)
                        add = {
                                "SPI": C.C1,
                                "AP": C.C2,
                                "NLS": C.N1,
                                }
                        NLS = int(C.N1)
                        LRU1_LRF1_2_NRO1.update(add)
                        for k in range(NLS) :
                            L, i = sandy.read_list(df, i)
                            add = {
                                    "AWRI": L.C1,       #Ratio of the mass of a particular isotope to that of a neutron
                                    "QX": L.C2,         #Q-value to be added to the incident particle’s center-of-mass energy to determine the channel energy for use in the penetrability factor.
                                    "L": L.L1,          #Value of l.
                                    "LRX": L.L2,        #Flag indicating whether this energy range contains a competitive width
                                    "NRS": L.N2,        #Number of resolved resonances for a given l-value.
                                    }
                            
                            LIST ={}
                            add_2 =  {
                                       "ER": L.B[0::6],
                                       "AJ": L.B[1::6],
                                       "GT": L.B[2::6],
                                       "GN": L.B[3::6],
                                       "GG": L.B[4::6],
                                       "GF": L.B[5::6],
                                      }
                            add.update(add_2)
                            LIST.update(add)
                            add.update({k:LIST})
                        LRU1_LRF1_2_NRO1.update({"L":add})
                        dico[(EL, EH)] = LRU1_LRF1_2_NRO1
               elif LRF == 3 :
                    if NRO == 0 :
                        LRU1_LRF3_NRO0 = {}
                        LRU1_LRF3_NRO0.update(header2)
                        C, i = sandy.read_cont(df, i)
                        add = {
                                "SPI": C.C1,
                                "AP": C.C2,
                                "LAD": C.L1,        #Flag indicating whether these parameters can be used to compute angular distributions.
                                "NLS": C.N1,        
                                "NLSC": C.N2,       #Number of l-values which must be used to converge the calculation with respect to the incident l-value in order to obtain accurate elastic angular distributions.
                                }
                        NLS = int(C.N1)
                        LRU1_LRF3_NRO0.update(add)
                        LRU1_LRF3_NRO0_NLS = {}
                        for k in range(NLS) :
                            L, i = sandy.read_list(df, i)
                            add = {
                                    "AWRI": L.C1,       #Ratio of the mass of a particular isotope to that of a neutron
                                    "APL": L.C2,         #Q-value to be added to the incident particle’s center-of-mass energy to determine the channel energy for use in the penetrability factor.
                                    "L": L.L1,          #Value of l.   
                                    "NRS": L.N2,        #Number of resolved resonances for a given l-value.
                                    }
                            LIST ={}
                            add_2 =  {
                                       "ER": L.B[0::6],
                                       "AJ": L.B[1::6],
                                       "GN": L.B[2::6],
                                       "GG": L.B[3::6],
                                       "GFA": L.B[4::6],        #First partial fission width, a constant
                                       "GFB": L.B[5::6],        #Second partial fission width, a constant.
                                      }
                            LIST.update(add_2)
                            add.update(LIST)
                            LRU1_LRF3_NRO0_NLS.update({k:add})
                        LRU1_LRF3_NRO0.update({"L":LRU1_LRF3_NRO0_NLS})
                        dico[(EL, EH)] = LRU1_LRF3_NRO0
                    else :
                        LRU1_LRF3_NRO1={}
                        T, i = sandy.read_tab1(df, i)
                        add = { 
                               "NR": T.NBT,
                               "NP": T.INT, 
                               "E_int": T.x, 
                               "AP": T.y, 
                               }
                        LRU1_LRF3_NRO1.update(add)
                        C, i = sandy.read_cont(df, i)
                        add = {
                                "SPI": C.C1,
                                "AP": C.C2,
                                "LAD": C.L1,
                                "NLS": C.N1,
                                "NLSC": C.N2
                                }
                        NLS = int(C.N1)
                        LRU1_LRF3_NRO1.update(add)
                        LRU1_LRF3_NRO1_NLS = {}
                        for k in range(NLS) :
                            L, i = sandy.read_list(df, i)
                            add = {
                                    "AWRI": L.C1,       #Ratio of the mass of a particular isotope to that of a neutron
                                    "APL": L.C2,         #Q-value to be added to the incident particle’s center-of-mass energy to determine the channel energy for use in the penetrability factor.
                                    "L": L.L1,          #Value of l.     
                                    "NRS": L.N2,        #Number of resolved resonances for a given l-value.
                                    }
                            LIST ={}
                            add_2 =  {
                                       "ER": L.B[0::6],
                                       "AJ": L.B[1::6],
                                       "GN": L.B[2::6],
                                       "GG": L.B[3::6],
                                       "GFA": L.B[4::6],
                                       "GFB": L.B[5::6],
                                      }
                            LIST.update(add_2)
                            add.update(LIST)
                            LRU1_LRF3_NRO1_NLS.update({k:add})
                        LRU1_LRF3_NRO1.update({"L":LRU1_LRF3_NRO1_NLS})
                        dico[(EL, EH)] = LRU1_LRF3_NRO1
               elif LRF == 4 :
                     raise ValueError("LRF = 4 is not supported in SANDY")
               elif LRF == 7 :
                    LRU1_LRF7={}
                    LRU1_LRF7.update(header2)
                    C, i = sandy.read_cont(df, i)
                    add = {
                                "IFG": C.L1,
                                "KRM": C.L2,        #Flag to specify which formulae for the R-matrix are to be used
                                "NJS": C.N1,        #Number of values of J
                                "KRL": C.N2,        #Flag is zero for non-relativistic kinematics
                                }
                    NJS = int(C.N1)
                    LRU1_LRF7.update(add)
                    L, i = sandy.read_list(df, i)
                    add = {"NPP": L.L1,                #Total number of particle-pairs.
                           }
                    LIST ={}
                    add_2 = {"MA": L.B[0::12],          #Mass of first particle in the pair
                             "MB": L.B[1::12],          #Mass of second particle
                             "ZA": L.B[2::12],          #Charge of first particle.
                             "ZB": L.B[3::12],          #Charge of second particle
                             "IA": L.B[4::12],          #Spin (and parity, if non-zero) of one particle in the pair
                             "IB": L.B[5::12],          #Spin of the other particle in the pair
                             "Q": L.B[6::12],           
                             "PNT": L.B[7::12],         #Flag if penetrability is to be calculated;
                             "SHF": L.B[8::12],         #Flag if shift factor is to be calculated
                             "MT": L.B[9::12],          #Reaction type associated with this particle-pair
                             "PA": L.B[10::12],         #Parity for first particle in the pair
                             "PB": L.B[11::12],         #Parity for second particle
                             }
                    LIST.update(add_2)
                    add.update(LIST)
                    LRU1_LRF7.update(add)
                    LRU1_LRF7_NJS = {}
                    for k in range(NJS):
                            LISTS = {}
                            L, i = sandy.read_list(df, i)
                            add1 = {
                                    "AJ": L.C1,             #Floating point value of J (spin); sign indicates parity
                                    "PJ": L.C2,             #Parity (used only if AJ = 0.0).
                                    "KBK": L.L1,            #Non-zero if background R-matrix exists
                                    "KPS": L.L2,            #Non-zero if non-hard-sphere phase shift are to be specified.
                                    "NCH": L.N2,            #Number of channels for the given J pi. 
                                    }
                            NCH = int (L.N2)
                            LIST1 = {}
                            add_2 =  {
                                       "PPI": L.B[0::6],    #Particle-pair number for this channel
                                       "L": L.B[1::6],      #Orbital angular momentum
                                       "SCH": L.B[2::6],    #Channel spin
                                       "BND": L.B[3::6],    #Boundary condition for this channel
                                       "APE": L.B[4::6],    #Effective channel radius
                                       "APT": L.B[5::6],    #True channel radius
                                      }
                            LIST1.update(add_2)
                            add1.update(LIST1)
                            L, i = sandy.read_list(df, i)
                            add2 = {
                                    "NRS": L.L2,            #Number of resonances for the given J pi 
                                    "NX": L.N2,             #Number of lines required for all resonances for the given J
                                    }
                            LIST2 ={}
                            Ep = L.B[::NCH+1]
                            b = L.B[::1]
                            del b[::NCH+1]
                            add_3 =  {
                                       "ER": Ep,            #Resonance energy in eV
                                       "GAM": b,            #Channel width in eV
                                      }
                            LIST2.update(add_3)
                            add2.update(LIST2)
                            LISTS.update(add1)
                            LISTS.update(add2)
                            LRU1_LRF7_NJS.update({k:LISTS})
                    LRU1_LRF7.update({"J":LRU1_LRF7_NJS})
                    dico[(EL, EH)] = LRU1_LRF7
            elif LRU == 2 :
                if LFW == 0 and LRF == 1 : 
                    LRU2_LFW0_LRF1 ={}
                    LRU2_LFW0_LRF1.update(header2)
                    C, i = sandy.read_cont(df, i)
                    add = {
                            "SPI": C.C1,
                            "AP": C.C2,
                            "LSSF": C.L1,                   #Flag governing the interpretation of the File 3 cross sections
                            "NLS": C.N1,                    #Number of l-values.
                            }
                    NLS = int(C.N1)
                    LRU2_LFW0_LRF1.update(add)
                    LRU2_LFW0_LRF1_NLS = {}
                    for k in range(NLS):
                        L, i = sandy.read_list(df, i)
                        add = {
                                "AWRI": L.C1,
                                "L": L.L1,
                                "NJS": L.N2,                    #Number of J-states for a particular l-state
                                    }
                        LIST ={}
                        add_2 =  {
                                   "D": L.B[0::6],
                                   "AJ": L.B[1::6],
                                   "AMUN": L.B[2::6],
                                   "GNO": L.B[3::6],
                                   "GG": L.B[4::6],
                                  }
                        LIST.update(add_2)
                        add.update(LIST)
                        LRU2_LFW0_LRF1_NLS.update({k : add})
                    LRU2_LFW0_LRF1.update({"L": LRU2_LFW0_LRF1_NLS})
                    dico[(EL, EH)] = LRU2_LFW0_LRF1
                elif LFW == 1 and LRF == 1 :
                    raise ValueError("LFW = 1  LRF = 1 is not supported in SANDY")
                elif LRF == 2 : 
                    LRU2_LRF2 = {}
                    LRU2_LRF2.update(header2)
                    C, i = sandy.read_cont(df, i)
                    add = {
                            "SPI": C.C1,
                            "AP": C.C2,
                            "LSSF": C.L1,
                            "NLS": C.N1,
                            }
                    NLS = int (C.N1)
                    LRU2_LRF2.update(add)
                    LRU2_LRF2_NLS = {}
                    for m in range(NLS) :
                        C, i = sandy.read_cont(df, i)
                        add_1 = {
                                "AWRI": C.C1,
                                "L": C.L1,
                                "NJS": C.N1,
                                }
                        NJS = int(C.N1)
                        LRU2_LRF2_NjS = {}
                        LIST ={}
                        for k in range(NJS):
                            L, i = sandy.read_list(df, i)
                            add_2 = {
                                    "AJ": L.C1,
                                    "INT": L.L1,                #Interpolation scheme to be used for interpolating between the cross sections obtained from average resonance parameters.
                                    "NE": L.N2,
                                    }
                            add_3 ={
                                    "AMUX": L.B[2],             #Number of degrees of freedom used in the competitive width distribution.
                                    "AMUN": L.B[3],             #Number of degrees of freedom in the neutron width distribution.
                                    "AMUG": L.B[4],             #Number of degrees of freedom in the radiation width distribution.
                                    "AMUF": L.B[5],             #Integer value of the number of degrees of freedom for fission widths.
                                    "ES": L.B[6::6],
                                    "D": L.B[7::6],
                                    "GX": L.B[8::6],
                                    "GNO": L.B[9::6],
                                    "GG": L.B[10::6],
                                    "GF": L.B[11::6],
                                    }
                            add_2.update(add_3)
                            LIST.update({k:add_2})
                        LRU2_LRF2_NjS.update(add_1)    
                        LRU2_LRF2_NjS.update({"J":LIST})
                        LRU2_LRF2_NLS.update({ m :LRU2_LRF2_NjS })
                    LRU2_LRF2.update({"L":LRU2_LRF2_NLS})
                    dico[(EL, EH)] = LRU2_LRF2
            NER1.update(dico)
            info.update({"NER":NER1})
            M.update(info)
            ISO.update({ZAI : M})
        P.update({"ISO": ISO})
    out.update({"NIS":P})
    return out
        
def write_mf2(sec): 
    """
    Given the content of MF2 write it to string.

    Parameters
    ----------
    sec : 'dict'
        Multiline string reproducing the content of a ENDF-6 section.

    Returns
    -------
    `str`
        Multiline string reproducing the content of a ENDF-6 section.
    
    Examples
    --------
    resonance parameters of Thorium 233
    LRU = 0
    >>> tape = sandy.get_endf6_file("jeff_33", "xs", 902330)
    >>> dic = sandy.read_mf2(tape, 9043)
    >>> text = sandy.write_mf2(dic)
    >>> print(text[:1000])
    90233.0000 231.040000          0          0          1          09043 2151    1
    90233.0000 1.00000000          0          0          1          09043 2151    2
    1.000000-5 1.90000000          0          0          0          09043 2151    3
    5.000000-1 9.765000-1          0          0          0          09043 2151    4
    
    resonance parameters of Thorium 230
    LRU = 1 LRF = 2
    >>> tape = sandy.get_endf6_file("jeff_33", "xs", 902300)
    >>> dic = sandy.read_mf2(tape, 9034)
    >>> text = sandy.write_mf2(dic)
    >>> print(text[:1054])
    90230.0000 228.060000          0          0          1          09034 2151    1
    90230.0000 1.00000000          0          0          1          09034 2151    2
    1.000000-5 251.000000          1          2          0          09034 2151    3
    0.00000000 8.322500-1          0          0          1          09034 2151    4
    228.060000 0.00000000          0          0        132         229034 2151    5
    -1.42700000 5.000000-1 2.510000-2 2.000000-4 2.490000-2 0.000000009034 2151    6
    1.42700000 5.000000-1 2.520000-2 3.000000-4 2.490000-2 0.000000009034 2151    7
    17.2700000 5.000000-1 3.510000-2 1.310000-2 2.200000-2 0.000000009034 2151    8
    23.8400000 5.000000-1 3.710000-2 1.110000-2 2.600000-2 0.000000009034 2151    9
    32.2000000 5.000000-1 3.230000-2 3.300000-3 2.900000-2 0.000000009034 2151   10
    39.8000000 5.000000-1 3.750000-2 8.500000-3 2.900000-2 0.000000009034 2151   11
    48.1000000 5.000000-1 3.900000-2 1.000000-2 2.900000-2 0.000000009034 2151   12
    64.5000000 5.000000-1 2.870000-2 3.100000-3 2.560000-2 0.000000009034 2151   13
              
     
    resonance parameters of Uranium 235
    LRU = 1 LRF = 3
    >>> tape = sandy.get_endf6_file("jeff_33", "xs", 922350)
    >>> dic = sandy.read_mf2(tape, 9228)
    >>> text = sandy.write_mf2(dic)
    >>> print(text[:1054])
    92235.0000 233.025000          0          0          1          09228 2151    1
    92235.0000 1.00000000          0          1          2          09228 2151    2
    1.000000-5 2250.00000          1          3          0          19228 2151    3
    3.50000000 9.602000-1          1          0          1          49228 2151    4
    233.025000 9.602000-1          0          0      19080       31809228 2151    5
    -75.4054200-3.00000000 5.072748-1 4.778164-2-4.870903-1-4.433456-19228 2151    6
    -5.25305700 4.00000000 1.217027-2 3.679714-2 1.956806-1-1.600387-19228 2151    7
    -4.808332-1-3.00000000 8.876626-5 3.922852-2 1.296617-1-8.053507-29228 2151    8
    -4.320587-1 4.00000000 3.325503-5 3.802396-2 1.670726-1-8.283233-39228 2151    9
    -3.657000-5 4.00000000 6.46080-11 3.998846-2-5.091200-4 9.353600-49228 2151   10
    2.638505-1-3.00000000 4.252450-6 4.573547-2 1.232188-1 6.114063-59228 2151   11
    1.13613600 4.00000000 1.526560-5 3.855000-2 6.673820-5 1.192506-19228 2151   12
    1.29868300 4.00000000 3.816210-7 3.860000-2-1.671125-4 1.735410-29228 2151   13
     
     
    resonance parameters of Iron 54
    LRU = 1 LRF = 7
    >>> tape = sandy.get_endf6_file("endfb_80", "xs", 260540)
    >>> dic = sandy.read_mf2(tape, 2625)
    >>> text = sandy.write_mf2(dic)
    >>> print(text[:1053])
    26054.0000 53.4762400          0          0          1          02625 2151    1
    26054.0000 1.00000000          0          0          1          02625 2151    2
    1.000000-5 1036000.00          1          7          0          12625 2151    3
    0.00000000 0.00000000          0          3          5          02625 2151    4
    0.00000000 0.00000000          2          0         24          42625 2151    5
    0.00000000 54.4663500 0.00000000 26.0000000 1.00000000 0.000000002625 2151    6
    0.00000000 0.00000000 0.00000000 102.000000 0.00000000 0.000000002625 2151    7
    1.00000000 53.4762400 0.00000000 26.0000000 5.000000-1 0.000000002625 2151    8
    0.00000000 1.00000000 0.00000000 2.00000000 0.00000000 1.000000002625 2151    9
    5.000000-1 0.00000000          0          0         12          22625 2151   10
    1.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.000000002625 2151   11
    2.00000000 0.00000000 5.000000-1 0.00000000 5.437300-1 5.437300-12625 2151   12
    0.00000000 0.00000000          0        148        888        1482625 2151   13

     
     
     resonance parameters of Iron 58
     LRU = 2 LFW = 0 LRF = 1
     >>> tape = sandy.get_endf6_file("endfb_80", "xs", 260580)
     >>> dic = sandy.read_mf2(tape, 2637)
     >>> text = sandy.write_mf2(dic)
     >>> print(text[:1053])
     26058.0000 57.4356000          0          0          1          02637 2151    1
     26058.0000 1.00000000          0          0          2          02637 2151    2
     1.000000-5 350000.000          1          3          0          02637 2151    3
     0.00000000 5.226000-1          1          0          5          52637 2151    4
     57.4356100 5.226000-1          0          0        204         342637 2151    5
     -552450.000 5.000000-1 6986.00000 2.446000-1 0.00000000 0.000000002637 2151    6
     -524307.500 5.000000-1 6805.70000 2.446000-1 0.00000000 0.000000002637 2151    7
     -496165.000 5.000000-1 6620.50000 2.446000-1 0.00000000 0.000000002637 2151    8
     -468022.500 5.000000-1 6430.00000 2.446000-1 0.00000000 0.000000002637 2151    9
     -439880.000 5.000000-1 6233.70000 2.446000-1 0.00000000 0.000000002637 2151   10
     -411737.500 5.000000-1 6031.00000 2.446000-1 0.00000000 0.000000002637 2151   11
     -383595.000 5.000000-1 5821.20000 2.446000-1 0.00000000 0.000000002637 2151   12
     -355452.500 5.000000-1 5603.60000 2.446000-1 0.00000000 0.000000002637 2151   13
    
    
    
     resonance parameters of Uranium 235
     LRU = 2 LRF = 2
     >>> tape = sandy.get_endf6_file("jeff_33", "xs", 922350)
     >>> dic = sandy.read_mf2(tape, 9228)
     >>> text = sandy.write_mf2(dic)
     >>> print(text[10046:12068])
     3.4229400-3.00000000 5.878595-4 3.672970-2-6.820955-2-1.940129-29228 2151  125
     54.1870000 4.00000000 1.004536-4 4.111675-2-5.317350-2-7.705571-39228 2151  126
     54.9713800-3.00000000 9.482413-4 4.426821-2-4.932297-2 2.162867-29228 2151  127
     55.1077500 4.00000000 2.049826-3 3.423696-2 1.238615-4-3.312716-29228 2151  128
     55.8415900 4.00000000 2.400311-3 3.885363-2-2.229854-1-2.704497-29228 2151  129
     56.0333400-3.00000000 6.820638-4 4.058219-2-1.184018-1 1.145587-29228 2151  130
     56.4703900 4.00000000 4.075371-3 4.062175-2 2.073727-2-8.683868-29228 2151  131
     57.3011700 4.00000000 1.692560-7 3.488284-2 1.545313-1-6.059518-29228 2151  132
     57.6945900-3.00000000 6.997812-4 3.866154-2-1.248931-1 6.105844-29228 2151  133
     58.0726200-3.00000000 1.766641-3 4.201754-2 5.672139-2 2.816238-39228 2151  134
     58.6701500 4.00000000 1.251797-3 3.985913-2-1.282329-1-1.093225-39228 2151  135
     59.7017500 4.00000000 7.955139-5 3.753891-2-6.452942-2-4.530963-39228 2151  136
     60.1662700-3.00000000 1.631662-3 4.269916-2 6.939703-2 2.568978-19228 2151  137
     60.7998400 4.00000000 5.306126-4 3.771821-2-1.692058-1 4.601213-59228 2151  138
     61.0760900-3.00000000 4.019216-4 4.062546-2-8.796944-2 4.825410-49228 2151  139
     61.3131700 4.00000000 1.820134-4 4.249260-2-2.965849-1-4.730334-19228 2151  140
     61.7583400-3.00000000 9.078133-6 4.451335-2 9.764961-3 2.378250-19228 2151  141
     62.5350600 4.00000000 5.510943-5 3.875537-2 5.874718-2 7.565506-29228 2151  142
     62.9074700-3.00000000 2.433310-6 4.539639-2-4.230665-1 1.806913-19228 2151  143
     63.6294000 4.00000000 7.504646-5 3.657835-2 1.533564-3 1.998567-19228 2151  144
     63.6562800-3.00000000 1.445787-3 4.294969-2 2.159020-1 7.735164-19228 2151  145
     64.2980700 4.00000000 1.013067-3 3.407947-2 4.462792-3-1.769018-59228 2151  146
     65.1483700 4.00000000 2.463481-8 4.170096-2 1.453683-2-2.180552-19228 2151  147
     65.7741300-3.00000000 4.208330-4 3.762744-2-1.429333-2-1.591547-29228 2151  148
     66.2503700-3.00000000 7.161637-5 3.759448-2 2.357995-2 1.894606-19228 2151  149
    """
    lines = sandy.write_cont(
        sec["INTRO"]["ZA"],
        sec["INTRO"]["AWR"],
        0,
        0,
        sec["INTRO"]["NIS"],
        0,
        )
    for ZAI, sec2 in sec["NIS"]["ISO"].items():
        lines += sandy.write_cont(
            ZAI,
            sec2["ABN"],
            0,
            sec2["LFW"],
            len(sec2["NER"]),
            0,
            )
        for (EL, EH), sec3 in sec2["NER"].items():
            lines += sandy.write_cont(
                sec3["EL"],
                sec3["EH"],
                sec3["LRU"],
                sec3["LRF"],
                sec3["NRO"],
                sec3["NAPS"],
                )
            if sec3["LRU"] == 0:
                lines += sandy.write_cont(
                    sec3["SPI"],
                    sec3["AP"],
                    0,
                    0,
                    sec3["NLS"],
                    0,
                    )
            elif sec3["LRU"] == 1:
                if sec3["LRF"] == 1 or sec3["LRF"] == 2 :
                    if sec3["NRO"] == 0 :
                        lines += sandy.write_cont(
                            sec3["SPI"],
                            sec3["AP"],
                            0,
                            0,
                            sec3["NLS"],
                            0,
                            )
                        NLS = int(sec3["NLS"])
                        for k in range(NLS): 
                            S_NLS = len(sec3["L"][k]["ER"])
                            add = [0] * (6 * S_NLS)
                            add[0::6] = sec3["L"][k]["ER"]
                            add[1::6] = sec3["L"][k]["AJ"]
                            add[2::6] = sec3["L"][k]["GT"]
                            add[3::6] = sec3["L"][k]["GN"]
                            add[4::6] = sec3["L"][k]["GG"]
                            add[5::6] = sec3["L"][k]["GF"]
                            lines +=sandy.write_list(
                                sec3["L"][k]["AWRI"],
                                sec3["L"][k]["QX"],
                                sec3["L"][k]["L"],
                                sec3["L"][k]["LRX"],
                                sec3["L"][k]["NRS"],
                                add,
                                )
                    else :
                        lines += sandy.write_tab1(
                                0,
                                0,
                                0,
                                0,
                                0,
                                sec3["NR"],
                                sec3["NP"],
                                sec3["E_int"],
                                sec3["AP"],
                                )
                        lines += sandy.write_cont(
                            sec3["SPI"],
                            sec3["AP"],
                            0,
                            0,
                            sec3["NLS"],
                            0,
                            )
                        NLS = int(sec3["NLS"])
                        for k in range(NLS): 
                            S_NLS = len(sec3["L"][k]["ER"])
                            add = [0] * (6 * S_NLS)
                            add[0::6] = sec3["L"][k]["ER"]
                            add[1::6] = sec3["L"][k]["AJ"]
                            add[2::6] = sec3["L"][k]["GT"]
                            add[3::6] = sec3["L"][k]["GN"]
                            add[4::6] = sec3["L"][k]["GG"]
                            add[5::6] = sec3["L"][k]["GF"]
                            lines +=sandy.write_list(
                                sec3["L"][k]["AWRI"],
                                sec3["L"][k]["QX"],
                                sec3["L"][k]["L"],
                                sec3["L"][k]["LRX"],
                                sec3["L"][k]["NRS"],
                                add,
                                )
                elif sec3["LRF"] == 3:
                    if sec3["NRO"] == 0: 
                        lines += sandy.write_cont(
                            sec3["SPI"],
                            sec3["AP"],
                            sec3["LAD"],
                            0,
                            sec3["NLS"],
                            sec3["NLSC"],
                            )
                        NLS = int(sec3["NLS"])
                        for k in range(NLS): 
                            S_ER = len (sec3["L"][k]["ER"])
                            add = [0]*(6*S_ER)
                            add[0::6] = sec3["L"][k]["ER"]
                            add[1::6] = sec3["L"][k]["AJ"]
                            add[2::6] = sec3["L"][k]["GN"]
                            add[3::6] = sec3["L"][k]["GG"]
                            add[4::6] = sec3["L"][k]["GFA"]
                            add[5::6] = sec3["L"][k]["GFB"]
                            lines +=sandy.write_list(
                                sec3["L"][k]["AWRI"],
                                sec3["L"][k]["APL"],
                                sec3["L"][k]["L"],
                                0,
                                sec3["L"][k]["NRS"],
                                add,
                                )
                    else :
                        lines += sandy.write_tab1(
                                0,
                                0,
                                0,
                                0,
                                0,
                                sec3["NR"],
                                sec3["NP"],
                                sec3["E_int"],
                                sec3["AP"],
                                )
                        lines += sandy.write_cont(
                            sec3["SPI"],
                            sec3["AP"],
                            sec3["LAD"],
                            0,
                            sec3["NLS"],
                            sec3["NLSC"],
                            )
                        NLS = int(sec3["NLS"])
                        for k in range(NLS): 
                            S_ER = len (sec3["L"][k]["ER"])
                            add = [0]*(6*S_ER)
                            add[0::6] = sec3["L"][k]["ER"]
                            add[1::6] = sec3["L"][k]["AJ"]
                            add[2::6] = sec3["L"][k]["GN"]
                            add[3::6] = sec3["L"][k]["GG"]
                            add[4::6] = sec3["L"][k]["GFA"]
                            add[5::6] = sec3["L"][k]["GFB"]
                            lines +=sandy.write_list(
                                sec3["L"][k]["AWRI"],
                                sec3["L"][k]["APL"],
                                sec3["L"][k]["L"],
                                0,
                                sec3["L"][k]["NRS"],
                                add,
                                )
                elif sec3["LRF"] == 7:
                    lines += sandy.write_cont(
                        0,
                        0,
                        sec3["IFG"],
                        sec3["KRM"],
                        sec3["NJS"],
                        sec3["KRL"],
                        )
                    NJS = int(sec3["NJS"]) 
                    S_NPP = len (sec3["MA"])
                    add = [0]*(12*S_NPP)
                    add[0::12] = sec3["MA"]
                    add[1::12] = sec3["MB"]
                    add[2::12] = sec3["ZA"]
                    add[3::12] = sec3["ZB"]
                    add[4::12] = sec3["IA"]
                    add[5::12] = sec3["IB"]
                    add[6::12] = sec3["Q"]
                    add[7::12] = sec3["PNT"]
                    add[8::12] = sec3["SHF"]
                    add[9::12] = sec3["MT"]
                    add[10::12] = sec3["PA"]
                    add[11::12] = sec3["PB"]
                    lines += sandy.write_list(
                        0,
                        0,
                        sec3["NPP"],
                        0,
                        2 * sec3["NPP"],
                        add,
                        )
                    for k in range(NJS):
                        S_NJS = len (sec3["J"][k]["PPI"])
                        add = [0]*(6*S_NJS)
                        add[0::6] = sec3["J"][k]["PPI"]
                        add[1::6] = sec3["J"][k]["L"]
                        add[2::6] = sec3["J"][k]["SCH"]
                        add[3::6] = sec3["J"][k]["BND"]
                        add[4::6] = sec3["J"][k]["APE"]
                        add[5::6] = sec3["J"][k]["APT"]
                        lines += sandy.write_list(
                            sec3["J"][k]["AJ"],
                            sec3["J"][k]["PJ"],
                            sec3["J"][k]["KBK"],
                            sec3["J"][k]["KPS"],
                            sec3["J"][k]["NCH"],
                            add,
                            )
                        NCH = sec3["J"][k]["NCH"]
                        size_E = len(sec3["J"][k]["ER"])
                        size_b = len(sec3["J"][k]["GAM"])
                        add2 = []
                        n = NCH 
                        x = [sec3["J"][k]["GAM"][i:i + n] for i in range(0, size_b, n)]
                        for i in range(size_E):
                            add2.append(sec3["J"][k]["ER"][i])
                            for m in range(len(x[i])):
                                add2.append(x[i][m])
                        lines += sandy.write_list(
                            0,
                            0,
                            0,
                            sec3["J"][k]["NRS"],
                            sec3["J"][k]["NX"],
                            add2,
                            )
            elif sec3["LRU"] == 2:
                if sec2["LFW"] == 0 and sec3["LRF"] == 1:
                    lines += sandy.write_cont(
                        sec3["SPI"],
                        sec3["AP"],
                        sec3["LSSF"],
                        0,
                        sec3["NLS"],
                        0,
                        )
                    NLS = int (sec3["NLS"])
                    for k in range(NLS):
                        S_NLS = len (sec3["L"][k]["AJ"])
                        add = [0]*(6*S_NLS)
                        add[0::6] = sec3["L"][k]["D"]
                        add[1::6] = sec3["L"][k]["AJ"]
                        add[2::6] = sec3["L"][k]["AMUN"]
                        add[3::6] = sec3["L"][k]["GNO"]
                        add[4::6] = sec3["L"][k]["GG"]
                        lines += sandy.write_list(
                            sec3["L"][k]["AWRI"],
                            0,
                            sec3["L"][k]["L"],
                            0,
                            sec3["L"][k]["NJS"],
                            add,
                            )
                elif sec3["LRF"] == 2:
                    lines += sandy.write_cont(
                        sec3["SPI"],
                        sec3["AP"],
                        sec3["LSSF"],
                        0,
                        sec3["NLS"],
                        0,
                        )
                    for m, sec4 in sec3["L"].items():
                        lines += sandy.write_cont(
                            sec4["AWRI"],
                            0,
                            sec4["L"],
                            0,
                            sec4["NJS"],
                            0,
                            )
                        NJS = int(sec4["NJS"])
                        for k in range(NJS):
                            S_NJS = len (sec4["J"][k]["ES"])
                            add1 = [0]*6
                            add2 = [0]*(6*S_NJS)
                            add1[2] = sec4["J"][k]["AMUX"]
                            add1[3] = sec4["J"][k]["AMUN"]
                            add1[4] = sec4["J"][k]["AMUG"]
                            add1[5] = sec4["J"][k]["AMUF"]
                            add2[0::6] = sec4["J"][k]["ES"]
                            add2[1::6] = sec4["J"][k]["D"]
                            add2[2::6] = sec4["J"][k]["GX"]
                            add2[3::6] = sec4["J"][k]["GNO"]
                            add2[4::6] = sec4["J"][k]["GG"]
                            add2[5::6] = sec4["J"][k]["GF"]
                            add = add1 + add2
                            lines += sandy.write_list(
                                sec4["J"][k]["AJ"],
                                0,
                                sec4["J"][k]["INT"],
                                0,
                                sec4["J"][k]["NE"],
                                add,
                                )
    return "\n".join(sandy.write_eol(lines, sec["MAT"], 2, 151))
