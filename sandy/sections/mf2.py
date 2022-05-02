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

import numpy as np 

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
    >>> tape = sandy.get_endf6_file("jeff_33", "xs",902330)
    >>> df = sandy.read_mf2(tape,9043)
    >>> df = {'MAT': 9043,
             'MF': 2,
             'MT': 151,
             'Intro': {'ZA': 90233.0,
              'AWR': 231.04,
              'NIS': 1,
              'ZAI': 90233.0,
              'ABN': 1.0,
              'LFW': 0,
              'NER': 1},
             'NIS 0': {'NER 0': {'header 0': {'EL': 1e-05,
                'EH': 1.9,
                'LRU': 0,
                'LRF': 0,
                'NRO': 0,
                'NAPS': 0},
               'LRU = 0': {'SPI': 0.5, 'AP': 0.9765, 'NLS': 0}}}}
    
    
    
    Resonance parameters of the Thorium 230 
    LRU = 1 LRF = 2
    >>> tape = sandy.get_endf6_file("jeff_33", "xs",902300)
    >>> df = sandy.read_mf2(tape,9034)
    >>> df[f"NIS {0}"][f"NER {0}"]["LRU = 1 LRF = 1 or 2 NRO = 0"]
        {'SPI': 0.0,
         'AP': 0.83225,
         'NLS': 1,
         'List 0': {'AWRI': 228.06,
          'QX': 0.0,
          'L': 0,
          'LRX': 0,
          '6*NRS': 132,
          'NRS': 22,
          'ER': array([ -1.427,   1.427,  17.27 ,  23.84 ,  32.2  ,  39.8  ,  48.1  ,
                  64.5  ,  75.6  ,  83.3  , 103.   , 116.1  , 133.8  , 139.   ,
                 148.2  , 171.4  , 184.2  , 195.   , 209.   , 226.   , 241.   ,
                 248.   ]),
          'AJ': array([0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5,
                 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5]),
          'GT': array([0.0251, 0.0252, 0.0351, 0.0371, 0.0323, 0.0375, 0.039 , 0.0287,
                 0.0283, 0.0504, 0.0308, 0.0659, 0.0331, 0.028 , 0.0312, 0.0509,
                 0.0489, 0.0709, 0.1009, 0.0484, 0.0311, 0.0809]),
          'GN': array([0.0002, 0.0003, 0.0131, 0.0111, 0.0033, 0.0085, 0.01  , 0.0031,
                 0.0027, 0.0248, 0.0052, 0.0403, 0.0075, 0.0024, 0.0056, 0.0253,
                 0.0233, 0.0453, 0.0753, 0.0228, 0.0055, 0.0553]),
          'GG': array([0.0249, 0.0249, 0.022 , 0.026 , 0.029 , 0.029 , 0.029 , 0.0256,
                 0.0256, 0.0256, 0.0256, 0.0256, 0.0256, 0.0256, 0.0256, 0.0256,
                 0.0256, 0.0256, 0.0256, 0.0256, 0.0256, 0.0256]),
          'GF': array([0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                 0., 0., 0., 0., 0.])}}
        
        
        
    Resonance parameters of the iron 58
    LRU = 1 LRF = 3 
    >>> tape = sandy.get_endf6_file("endfb_80", "xs",260580)
    >>> df = sandy.read_mf2(tape,2637)
    >>> df[f"NIS {0}"][f"NER {0}"]["LRU = 1 LRF = 3 NRO = 0"][f"List {0}"]["ER"]
    [-552450.0,
     -524307.5,
     -496165.0,
     -468022.5,
     -439880.0,
     -411737.5,
     -383595.0,
     -355452.5,
     -327310.0,
     -299167.5,
     -296610.0,
     -265909.0,
     -235208.0,
     -204507.0,
     -173806.0,
     -143105.0,
     -112404.0,
     -81703.0,
     -56180.0,
     -30660.0,
     -3612.8,
     10426.21,
     37665.85,
     43818.0,
     67056.94,
     93642.29,
     121186.2,
     160244.7,
     178828.5,
     239987.5,
     265800.0,
     310193.7,
     320294.1,
     348105.0]
    
    
    Resonance parameters of the iron 58
    LRU = 2 LFW = 0 LRF = 1
    >>> tape = sandy.get_endf6_file("endfb_80", "xs",260580)
    >>> df = sandy.read_mf2(tape,2637)
    >>> df[f"NIS {0}"][f"NER {1}"]["LRU = 2 LFW = 0 LRF = 1"][f"List {0}"]
    {'AWRI': 57.43561,
     'L': 0,
     '6*NJS': 6,
     'NJS': 1,
     'D': [24100.0],
     'AJ': [0.5],
     'AMUN': [1.0],
     'GNO': [7.23],
     'GG': [0.9]}
    
    
    Resonance parameters of the iron 54
    LRU = 1 LRF = 7
    >>> tape = sandy.get_endf6_file("endfb_80", "xs",260540)
    >>> df = sandy.read_mf2(tape,2625)
    >>> df[f"NIS {0}"][f"NER {0}"]["LRU = 1 LRF = 7"][f"NJS {0}"]["PPI"]
    [1.0, 2.0]
    
    
    Resonance parameters of the Uranium 235
    LRU = 2 LRF = 2
    >>> tape = sandy.get_endf6_file("jeff_33", "xs",922350)
    >>> df = sandy.read_mf2(tape,9228)
    >>> df[f"NIS {0}"][f"NER {1}"]["LRU = 2 LRF = 2"][f"NJS {0, 0}"][f"List {0}"]["AMUX"]
    1.0
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
            "ZA": C.C1,             
            "AWR": C.C2,            #AWR is defines as the ratio of the mass of the material to that of the neutron
            "NIS": C.N1,            #NIS is the number of isotopes in the material
            }
    NIS = int(C.N1)
    I.update(add)
    C, i = sandy.read_cont(df, i)
    add = {
            "ZAI": C.C1,
            "ABN": C.C2,           #Abundance of an isotope in the material 
            "LFW": C.L2,           #indication whether average fission wifths are given in the unresolbed resonance region
            "NER": C.N1,           #Number of resonance energy ranges for this isotope. 
            }   
    NER = int(C.N1)
    LFW = int(C.L2)
    I.update(add)
    out.update({"Intro": I})
    P = {}
    for l in range(NIS): 
        M = {}
        for j in range(NER):
            info = {}
            C, i = sandy.read_cont(df, i)
            add = {
                    "EL": C.C1,             #Lower limit for an energy range
                    "EH": C.C2,             #Upper limit for an energy range
                    "LRU": C.L1,            #Flag indicating whether this energy range contains data for resolved or unresolved resonance parameters:
                    "LRF": C.L2,            #Flag indicating which representation has been used for the energy range. 
                    "NRO": C.N1,            #Flag designating possible energy dependence of the scattering radiu
                    "NAPS":C.N2,            #Flag controlling the use of the two radii
                    }
            LRU = int(C.L1)
            LRF = int(C.L2)
            NRO = int(C.N1)
            info.update({f"header {j}": add})
            M.update({f"NER {j}": info})
            if LRU == 0:
                LRU0 = {}
                C, i = sandy.read_cont(df, i)
                add = {
                        "SPI": C.C1,        #Spin, I, of the target nucleus.
                        "AP": C.C2,         #Scattering radius in units of 10e-12cm.
                        "NLS": C.N1,        #Number of l-values (neutron orbital angular momentum) in this energy region.
                        }
                LRU0.update(add)
                info.update({"LRU = 0":LRU0})
                M.update({f"NER {j}": info})
            if LRU == 1:
               if LRF == 1 or LRF == 2 :
                    if NRO == 0 :
                        LRU1_LRF1_2_NRO0 = {}
                        C, i = sandy.read_cont(df, i)
                        add = {
                                "SPI": C.C1,
                                "AP": C.C2,
                                "NLS": C.N1,
                                }
                        NLS = int(C.N1)
                        LRU1_LRF1_2_NRO0.update(add)
                        for k in range(NLS) :
                            L, i = sandy.read_list(df, i)
                            add = {
                                    "AWRI": L.C1,       #Ratio of the mass of a particular isotope to that of a neutron
                                    "QX": L.C2,         #Q-value to be added to the incident particle’s center-of-mass energy to determine the channel energy for use in the penetrability factor.
                                    "L": L.L1,          #Value of l.
                                    "LRX": L.L2,        #Flag indicating whether this energy range contains a competitive width
                                    "6*NRS": L.NPL,     
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
                            LIST.update(add_2)
                            add.update(LIST)
                            LRU1_LRF1_2_NRO0.update({f"List {k}":add})
                        info.update({"LRU = 1 LRF = 1 or 2 NRO = 0" : LRU1_LRF1_2_NRO0})
                        M.update({f"NER {j}": info})
                    else :   
                        LRU1_LRF1_2_NRO1 = {}
                        T, i = sandy.read_tab1(df, i)
                        add = {
                                "NR": T.NBT,    
                                "NP": T.INT,
                                "E_int": T.x,
                                "AP": T.y,
                                }
                        LRU1_LRF1_2_NRO1.update(add)
                        out.update({"LRU = 1 LRF = 1 or 2": LRU1_LRF1_2_NRO1})
                        M.update({f"NER {j}": info})
               if LRF == 3 :
                    if NRO == 0 :
                        
                        LRU1_LRF3_NRO0 = {}
                        C, i = sandy.read_cont(df, i)
                        add = {
                                "SPI": C.C1,
                                "AP": C.C2,
                                "LAD": C.L1,
                                "NLS": C.N1,
                                "NLSC": C.N2
                                }
                        NLS = int(C.N1)
                        LRU1_LRF3_NRO0.update(add)
                        for k in range(NLS) :
                            L, i = sandy.read_list(df, i)
                            add = {
                                    "AWRI": L.C1,       #Ratio of the mass of a particular isotope to that of a neutron
                                    "APL": L.C2,         #Q-value to be added to the incident particle’s center-of-mass energy to determine the channel energy for use in the penetrability factor.
                                    "L": L.L1,          #Value of l.
                                    "6*NRS": L.NPL,     
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
                            add.update( LIST)
                            LRU1_LRF3_NRO0.update({f"List {k}":add})
                        info.update({"LRU = 1 LRF = 3 NRO = 0" : LRU1_LRF3_NRO0})
                        M.update({f"NER {j}": info})
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
                        out.update({"LRU = 1 LRF = 3": LRU1_LRF3_NRO1})
                        M.update({f"NER {j}": info})
               if LRF == 4 :
                     raise ValueError("LRF = 4 is not supported in SANDY")
               if LRF == 7 :
                    LRU1_LRF7={}
                    C, i = sandy.read_cont(df, i)
                    add = {
                                "IFG": C.L1,
                                "KRM": C.L2,
                                "NJS": C.N1,
                                "KRL": C.N2,
                                }
                    NJS = int(C.N1)
                    LRU1_LRF7.update(add)
                    L, i = sandy.read_list(df, i)
                    add = {"NPP": L.L1,
                           "12*NPP": L.NPL,
                           "2*NPP": L.N2,                  #Total number of particle-pairs.
                           }
                    LIST ={}
                    add_2 = {"MA": L.B[0::12],
                             "MB": L.B[1::12],
                             "ZA": L.B[2::12],
                             "ZB": L.B[3::12],
                             "IA": L.B[4::12],
                             "IB": L.B[5::12],
                             "Q": L.B[6::12],
                             "PNT": L.B[7::12],
                             "SHF": L.B[8::12],
                             "MT": L.B[9::12],
                             "PA": L.B[10::12],
                             "PB": L.B[11::12],
                             }
                    LIST.update(add_2)
                    add.update(LIST)
                    LRU1_LRF7.update(add)
                    for k in range(NJS):
                            L, i = sandy.read_list(df, i)
                            add = {
                                    "AJ": L.C1,             #Floating point value of J (spin); sign indicates parity
                                    "PJ": L.C2,             #Parity (used only if AJ = 0.0).
                                    "KBK": L.L1,            #Non-zero if background R-matrix exists
                                    "KPS": L.L2,            #Non-zero if non-hard-sphere phase shift are to be specified.
                                    "6*NCH": L.NPL,
                                    "NCH": L.N2,            #Number of channels for the given J pi. 
                                    }
                            NCH = int (L.N2)
                            LIST = {}
                            add_2 =  {
                                       "PPI": L.B[0::6],
                                       "L": L.B[1::6],
                                       "SCH": L.B[2::6],
                                       "BND": L.B[3::6],
                                       "APE": L.B[4::6],
                                       "APT": L.B[5::6],
                                      }
                            LIST.update(add_2)
                            add.update(LIST)
                            L, i = sandy.read_list(df, i)
                            add2 = {
                                    "NRS": L.L2,            #Number of resonances for the given J pi 
                                    "6*NX": L.NPL,
                                    "NX": L.N2,
                                    }
                            LIST ={}
                            Ep = np.array(L.B[::NCH+1])
                            b = L.B[::1]
                            del b[::NCH+1]
                            add_2 =  {
                                       "ER": Ep,
                                       "GAM": b,
                                      }
                            LIST.update(add_2)
                            add2.update({f"List {k}" :LIST})
                            add.update(add2) 
                            LRU1_LRF7.update({f"NJS {k}" :add})  
                    info.update({"LRU = 1 LRF = 7": LRU1_LRF7}) 
                    M.update({f"NER {j}": info})
            if LRU == 2 :
                if LFW == 0 and LRF == 1 : 
                    LRU2_LFW0_LRF1 ={}
                    C, i = sandy.read_cont(df, i)
                    add = {
                            "SPI": C.C1,
                            "AP": C.C2,
                            "LSSF": C.L1,                   #Flag governing the interpretation of the File 3 cross sections
                            "NLS": C.N1,                    #Number of l-values.
                            }
                    NLS = int(C.N1)
                    LRU2_LFW0_LRF1.update(add)
                    for k in range(NLS):
                        L, i = sandy.read_list(df, i)
                        add = {
                                "AWRI": L.C1,
                                "L": L.L1,
                                "6*NJS": L.NPL,
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
                        LRU2_LFW0_LRF1.update({f"List {k}" : add})
                    info.update({"LRU = 2 LFW = 0 LRF = 1" : LRU2_LFW0_LRF1})
                    M.update({f"NER {j}": info})
                if LFW == 1 and LRF == 1 :
                    raise ValueError("LRF = 4 is not supported in SANDY")
                if LRF == 2 : 
                    LRU2_LRF2 = {}
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
                        add1 = {
                                "AWRI": C.C1,
                                "L": C.L1,
                                "NJS": C.N1,
                                }
                        NJS = int(C.N1)
                        LRU2_LRF2_NLS.update({f"NLS {m}" : add1})
                        LRU2_LRF2_NJS = {}
                        for k in range(NJS):
                            L, i = sandy.read_list(df, i)
                            add2 = {
                                    "AJ": L.C1,
                                    "INT": L.L1,                #Interpolation scheme to be used for interpolating between the cross sections obtained from average resonance parameters.
                                    "6*NE+6": L.NPL,
                                    "NE": L.N2,
                                    }
                            LIST ={}
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
                            LIST.update(add_3)       
                            add2.update( {f"List {k}" :LIST})
                            LRU2_LRF2_NJS.update({f"NJS {m,k}": add2})
                            LRU2_LRF2_NLS.update(LRU2_LRF2_NJS)
                        LRU2_LRF2.update(LRU2_LRF2_NLS)
                    info.update({"LRU = 2 LRF = 2": LRU2_LRF2})
                    M.update({f"NER {j}": info})
    P.update({f"NIS {l}": M})
    out.update(P)
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
    >>> tape = sandy.get_endf6_file("jeff_33", "xs",902330)
    >>> df = sandy.read_mf2(tape, 9043)
    >>> text = sandy.write_mf2(df)
    >>> print(text[:1000])
    90233.0000 231.040000          0          0          1          09043 2151    1
    90233.0000 1.00000000          0          0          1          09043 2151    2
    1.000000-5 1.90000000          0          0          0          09043 2151    3
    5.000000-1 9.765000-1          0          0          0          09043 2151    4
    
    resonance parameters of Thorium 230
    LRU = 1 LRF = 2
    >>> tape = sandy.get_endf6_file("jeff_33", "xs",902300)
    >>> df = sandy.read_mf2(tape, 9034)
    >>> text = sandy.write_mf2(df)
    >>> print(text[:1000])
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
     64.5000000 5.000000-1 2.870
              
     
    resonance parameters of Uranium 235
    LRU = 1 LRF = 3
    >>> tape = sandy.get_endf6_file("jeff_33", "xs",922350)
    >>> df = sandy.read_mf2(tape, 9228)
    >>> text = sandy.write_mf2(df)
    >>> print(text[:1000])
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
    1.29868300 4.00000000 3.816
     
     
     resonance parameters of Iron 54
     LRU = 1 LRF = 7
     >>> tape = sandy.get_endf6_file("endfb_80", "xs",260540)
     >>> df = sandy.read_mf2(tape, 2625)
     >>> text = sandy.write_mf2(df)
     >>> print(text[:1000])
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
     0.00000000 0.00000000
     
     
     resonance parameters of Iron 58
     LRU = 2 LFW = 0 LRF = 1
     >>> tape = sandy.get_endf6_file("endfb_80", "xs",260580)
     >>> df = sandy.read_mf2(tape, 2637)
     >>> text = sandy.write_mf2(df)
     >>> print(text[:1000])
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
    -355452.500 5.000000-1 5603.
    
    
    
     resonance parameters of Uranium 235
     LRU = 2 LRF = 2
     >>> tape = sandy.get_endf6_file("jeff_33", "xs",922350)
     >>> df = sandy.read_mf2(tape, 9228)
     >>> text = sandy.write_mf2(df)
     >>> print(text[:1000])
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
     1.29868300 4.00000000 3.816
    """
    lines = sandy.write_cont(
        sec["Intro"]["ZA"],
        sec["Intro"]["AWR"],
        0,
        0,
        sec["Intro"]["NIS"],
        0,
        )
    NIS = int(sec["Intro"]["NIS"])
    lines += sandy.write_cont(
        sec["Intro"]["ZAI"],
        sec["Intro"]["ABN"],
        0,
        sec["Intro"]["LFW"],
        sec["Intro"]["NER"],
        0,
        )
    NER = int(sec["Intro"]["NER"])
    LFW = int(sec["Intro"]["LFW"])
    for l in range(NIS): 
        for j in range(NER):
            lines += sandy.write_cont(
                sec[f"NIS {l}"][f"NER {j}"][f"header {j}"]["EL"],
                sec[f"NIS {l}"][f"NER {j}"][f"header {j}"]["EH"],
                sec[f"NIS {l}"][f"NER {j}"][f"header {j}"]["LRU"],
                sec[f"NIS {l}"][f"NER {j}"][f"header {j}"]["LRF"],
                sec[f"NIS {l}"][f"NER {j}"][f"header {j}"]["NRO"],
                sec[f"NIS {l}"][f"NER {j}"][f"header {j}"]["NAPS"],
                )
            LRU = int(sec[f"NIS {l}"][f"NER {j}"][f"header {j}"]["LRU"])
            LRF = int(sec[f"NIS {l}"][f"NER {j}"][f"header {j}"]["LRF"])
            NRO = int(sec[f"NIS {l}"][f"NER {j}"][f"header {j}"]["NRO"])
            if LRU == 0:
                lines += sandy.write_cont(
                    sec[f"NIS {l}"][f"NER {j}"]["LRU = 0"]["SPI"],
                    sec[f"NIS {l}"][f"NER {j}"]["LRU = 0"]["AP"],
                    0,
                    0,
                    sec[f"NIS {l}"][f"NER {j}"]["LRU = 0"]["NLS"],
                    0,
                    )
            if LRU == 1:
                if LRF == 1 or LRF == 2 :
                    if NRO == 0 :
                        lines += sandy.write_cont(
                            sec[f"NIS {l}"][f"NER {j}"]["LRU = 1 LRF = 1 or 2 NRO = 0"]["SPI"],
                            sec[f"NIS {l}"][f"NER {j}"]["LRU = 1 LRF = 1 or 2 NRO = 0"]["AP"],
                            0,
                            0,
                            sec[f"NIS {l}"][f"NER {j}"]["LRU = 1 LRF = 1 or 2 NRO = 0"]["NLS"],
                            0,
                            )
                        NLS = int(sec[f"NIS {l}"][f"NER {j}"]["LRU = 1 LRF = 1 or 2 NRO = 0"]["NLS"])
                        for k in range(NLS): 
                            S_NLS = len (sec[f"NIS {l}"][f"NER {j}"]["LRU = 1 LRF = 1 or 2 NRO = 0"][f"List {k}"]["ER"])
                            add = [0]*(6*S_NLS)
                            add[0::6] = sec[f"NIS {l}"][f"NER {j}"]["LRU = 1 LRF = 1 or 2 NRO = 0"][f"List {k}"]["ER"]
                            add[1::6] = sec[f"NIS {l}"][f"NER {j}"]["LRU = 1 LRF = 1 or 2 NRO = 0"][f"List {k}"]["AJ"]
                            add[2::6] = sec[f"NIS {l}"][f"NER {j}"]["LRU = 1 LRF = 1 or 2 NRO = 0"][f"List {k}"]["GT"]
                            add[3::6] = sec[f"NIS {l}"][f"NER {j}"]["LRU = 1 LRF = 1 or 2 NRO = 0"][f"List {k}"]["GN"]
                            add[4::6] = sec[f"NIS {l}"][f"NER {j}"]["LRU = 1 LRF = 1 or 2 NRO = 0"][f"List {k}"]["GG"]
                            add[5::6] = sec[f"NIS {l}"][f"NER {j}"]["LRU = 1 LRF = 1 or 2 NRO = 0"][f"List {k}"]["GF"]
                            lines +=sandy.write_list(
                                sec[f"NIS {l}"][f"NER {j}"]["LRU = 1 LRF = 1 or 2 NRO = 0"][f"List {k}"]["AWRI"],
                                sec[f"NIS {l}"][f"NER {j}"]["LRU = 1 LRF = 1 or 2 NRO = 0"][f"List {k}"]["QX"],
                                sec[f"NIS {l}"][f"NER {j}"]["LRU = 1 LRF = 1 or 2 NRO = 0"][f"List {k}"]["L"],
                                sec[f"NIS {l}"][f"NER {j}"]["LRU = 1 LRF = 1 or 2 NRO = 0"][f"List {k}"]["LRX"],
                                sec[f"NIS {l}"][f"NER {j}"]["LRU = 1 LRF = 1 or 2 NRO = 0"][f"List {k}"]["NRS"],
                                add,
                                )
                    else :
                        lines += sandy.write_tab1(
                                0,
                                0,
                                0,
                                0,
                                0,
                                sec[f"NIS {l}"][f"NER {j}"]["LRU = 1 LRF = 1 or 2"]["NR"],
                                sec[f"NIS {l}"][f"NER {j}"]["LRU = 1 LRF = 1 or 2"]["NP"],
                                sec[f"NIS {l}"][f"NER {j}"]["LRU = 1 LRF = 1 or 2"]["E_int"],
                                sec[f"NIS {l}"][f"NER {j}"]["LRU = 1 LRF = 1 or 2"]["AP"],
                                )
                if LRF == 3:
                    if NRO == 0: 
                        lines += sandy.write_cont(
                            sec[f"NIS {l}"][f"NER {j}"]["LRU = 1 LRF = 3 NRO = 0"]["SPI"],
                            sec[f"NIS {l}"][f"NER {j}"]["LRU = 1 LRF = 3 NRO = 0"]["AP"],
                            sec[f"NIS {l}"][f"NER {j}"]["LRU = 1 LRF = 3 NRO = 0"]["LAD"],
                            0,
                            sec[f"NIS {l}"][f"NER {j}"]["LRU = 1 LRF = 3 NRO = 0"]["NLS"],
                            sec[f"NIS {l}"][f"NER {j}"]["LRU = 1 LRF = 3 NRO = 0"]["NLSC"],
                            )
                        NLS = int(sec[f"NIS {l}"][f"NER {j}"]["LRU = 1 LRF = 3 NRO = 0"]["NLS"])
                        for k in range(NLS): 
                            S_ER = len (sec[f"NIS {l}"][f"NER {j}"]["LRU = 1 LRF = 3 NRO = 0"][f"List {k}"]["ER"])
                            add = [0]*(6*S_ER)
                            add[0::6] = sec[f"NIS {l}"][f"NER {j}"]["LRU = 1 LRF = 3 NRO = 0"][f"List {k}"]["ER"]
                            add[1::6] = sec[f"NIS {l}"][f"NER {j}"]["LRU = 1 LRF = 3 NRO = 0"][f"List {k}"]["AJ"]
                            add[2::6] = sec[f"NIS {l}"][f"NER {j}"]["LRU = 1 LRF = 3 NRO = 0"][f"List {k}"]["GN"]
                            add[3::6] = sec[f"NIS {l}"][f"NER {j}"]["LRU = 1 LRF = 3 NRO = 0"][f"List {k}"]["GG"]
                            add[4::6] = sec[f"NIS {l}"][f"NER {j}"]["LRU = 1 LRF = 3 NRO = 0"][f"List {k}"]["GFA"]
                            add[5::6] = sec[f"NIS {l}"][f"NER {j}"]["LRU = 1 LRF = 3 NRO = 0"][f"List {k}"]["GFB"]
                            lines +=sandy.write_list(
                                sec[f"NIS {l}"][f"NER {j}"]["LRU = 1 LRF = 3 NRO = 0"][f"List {k}"]["AWRI"],
                                sec[f"NIS {l}"][f"NER {j}"]["LRU = 1 LRF = 3 NRO = 0"][f"List {k}"]["APL"],
                                sec[f"NIS {l}"][f"NER {j}"]["LRU = 1 LRF = 3 NRO = 0"][f"List {k}"]["L"],
                                0,
                                sec[f"NIS {l}"][f"NER {j}"]["LRU = 1 LRF = 3 NRO = 0"][f"List {k}"]["NRS"],
                                add,
                                )
                    else :
                        lines += sandy.write_tab1(
                                0,
                                0,
                                0,
                                0,
                                0,
                                sec[f"NIS {l}"][f"NER {j}"]["LRU = 1 LRF = 3"]["NR"],
                                sec[f"NIS {l}"][f"NER {j}"]["LRU = 1 LRF = 3"]["NP"],
                                sec[f"NIS {l}"][f"NER {j}"]["LRU = 1 LRF = 3"]["E_int"],
                                sec[f"NIS {l}"][f"NER {j}"]["LRU = 1 LRF = 3"]["AP"],
                                )
                if LRF == 7:
                    lines += sandy.write_cont(
                        0,
                        0,
                        sec[f"NIS {l}"][f"NER {j}"]["LRU = 1 LRF = 7"]["IFG"],
                        sec[f"NIS {l}"][f"NER {j}"]["LRU = 1 LRF = 7"]["KRM"],
                        sec[f"NIS {l}"][f"NER {j}"]["LRU = 1 LRF = 7"]["NJS"],
                        sec[f"NIS {l}"][f"NER {j}"]["LRU = 1 LRF = 7"]["KRL"],
                        )
                    NJS = int(sec[f"NIS {l}"][f"NER {j}"]["LRU = 1 LRF = 7"]["NJS"]) 
                    S_NPP = len (sec[f"NIS {l}"][f"NER {j}"]["LRU = 1 LRF = 7"]["MA"])
                    add = [0]*(12*S_NPP)
                    add[0::12] = sec[f"NIS {l}"][f"NER {j}"]["LRU = 1 LRF = 7"]["MA"]
                    add[1::12] = sec[f"NIS {l}"][f"NER {j}"]["LRU = 1 LRF = 7"]["MB"]
                    add[2::12] = sec[f"NIS {l}"][f"NER {j}"]["LRU = 1 LRF = 7"]["ZA"]
                    add[3::12] = sec[f"NIS {l}"][f"NER {j}"]["LRU = 1 LRF = 7"]["ZB"]
                    add[4::12] = sec[f"NIS {l}"][f"NER {j}"]["LRU = 1 LRF = 7"]["IA"]
                    add[5::12] = sec[f"NIS {l}"][f"NER {j}"]["LRU = 1 LRF = 7"]["IB"]
                    add[6::12] = sec[f"NIS {l}"][f"NER {j}"]["LRU = 1 LRF = 7"]["Q"]
                    add[7::12] = sec[f"NIS {l}"][f"NER {j}"]["LRU = 1 LRF = 7"]["PNT"]
                    add[8::12] = sec[f"NIS {l}"][f"NER {j}"]["LRU = 1 LRF = 7"]["SHF"]
                    add[9::12] = sec[f"NIS {l}"][f"NER {j}"]["LRU = 1 LRF = 7"]["MT"]
                    add[10::12] = sec[f"NIS {l}"][f"NER {j}"]["LRU = 1 LRF = 7"]["PA"]
                    add[11::12] = sec[f"NIS {l}"][f"NER {j}"]["LRU = 1 LRF = 7"]["PB"]
                    lines += sandy.write_list(
                        0,
                        0,
                        sec[f"NIS {l}"][f"NER {j}"]["LRU = 1 LRF = 7"]["NPP"],
                        0,
                        sec[f"NIS {l}"][f"NER {j}"]["LRU = 1 LRF = 7"]["2*NPP"],
                        add,
                        )
                    for k in range(NJS):
                        S_NJS = len (sec[f"NIS {l}"][f"NER {j}"]["LRU = 1 LRF = 7"][f"NJS {k}"]["PPI"])
                        add = [0]*(6*S_NJS)
                        add[0::6] = sec[f"NIS {l}"][f"NER {j}"]["LRU = 1 LRF = 7"][f"NJS {k}"]["PPI"]
                        add[1::6] = sec[f"NIS {l}"][f"NER {j}"]["LRU = 1 LRF = 7"][f"NJS {k}"]["L"]
                        add[2::6] = sec[f"NIS {l}"][f"NER {j}"]["LRU = 1 LRF = 7"][f"NJS {k}"]["SCH"]
                        add[3::6] = sec[f"NIS {l}"][f"NER {j}"]["LRU = 1 LRF = 7"][f"NJS {k}"]["BND"]
                        add[4::6] = sec[f"NIS {l}"][f"NER {j}"]["LRU = 1 LRF = 7"][f"NJS {k}"]["APE"]
                        add[5::6] = sec[f"NIS {l}"][f"NER {j}"]["LRU = 1 LRF = 7"][f"NJS {k}"]["APT"]
                        lines += sandy.write_list(
                            sec[f"NIS {l}"][f"NER {j}"]["LRU = 1 LRF = 7"][f"NJS {k}"]["AJ"],
                            sec[f"NIS {l}"][f"NER {j}"]["LRU = 1 LRF = 7"][f"NJS {k}"]["PJ"],
                            sec[f"NIS {l}"][f"NER {j}"]["LRU = 1 LRF = 7"][f"NJS {k}"]["KBK"],
                            sec[f"NIS {l}"][f"NER {j}"]["LRU = 1 LRF = 7"][f"NJS {k}"]["KPS"],
                            sec[f"NIS {l}"][f"NER {j}"]["LRU = 1 LRF = 7"][f"NJS {k}"]["NCH"],
                            add,
                            )
                        NCH = sec[f"NIS {l}"][f"NER {j}"]["LRU = 1 LRF = 7"][f"NJS {k}"]["NCH"]
                        size_E = len(sec[f"NIS {l}"][f"NER {j}"]["LRU = 1 LRF = 7"][f"NJS {k}"][f"List {k}"]["ER"])
                        size_b = len(sec[f"NIS {l}"][f"NER {j}"]["LRU = 1 LRF = 7"][f"NJS {k}"][f"List {k}"]["GAM"])
                        add2 = []
                        n = NCH 
                        x = [sec[f"NIS {l}"][f"NER {j}"]["LRU = 1 LRF = 7"][f"NJS {k}"][f"List {k}"]["GAM"][i:i + n] for i in range(0, size_b, n)]
                        for i in range(size_E):
                            add2.append(sec[f"NIS {l}"][f"NER {j}"]["LRU = 1 LRF = 7"][f"NJS {k}"][f"List {k}"]["ER"][i])
                            for m in range(len(x[i])):
                                add2.append(x[i][m])
                        lines += sandy.write_list(
                            0,
                            0,
                            0,
                            sec[f"NIS {l}"][f"NER {j}"]["LRU = 1 LRF = 7"][f"NJS {k}"]["NRS"],
                            sec[f"NIS {l}"][f"NER {j}"]["LRU = 1 LRF = 7"][f"NJS {k}"]["NX"],
                            add2,
                            )
            if LRU == 2:
                if LFW == 0 and LRF == 1:
                    lines += sandy.write_cont(
                        sec[f"NIS {l}"][f"NER {j}"]["LRU = 2 LFW = 0 LRF = 1"]["SPI"],
                        sec[f"NIS {l}"][f"NER {j}"]["LRU = 2 LFW = 0 LRF = 1"]["AP"],
                        sec[f"NIS {l}"][f"NER {j}"]["LRU = 2 LFW = 0 LRF = 1"]["LSSF"],
                        0,
                        sec[f"NIS {l}"][f"NER {j}"]["LRU = 2 LFW = 0 LRF = 1"]["NLS"],
                        0,
                        )
                    NLS = int (sec[f"NIS {l}"][f"NER {j}"]["LRU = 2 LFW = 0 LRF = 1"]["NLS"])
                    for k in range(NLS):
                        S_NLS = len (sec[f"NIS {l}"][f"NER {j}"]["LRU = 2 LFW = 0 LRF = 1"][f"List {k}"]["AJ"])
                        add = [0]*(6*S_NLS)
                        add[0::6] = sec[f"NIS {l}"][f"NER {j}"]["LRU = 2 LFW = 0 LRF = 1"][f"List {k}"]["D"]
                        add[1::6] = sec[f"NIS {l}"][f"NER {j}"]["LRU = 2 LFW = 0 LRF = 1"][f"List {k}"]["AJ"]
                        add[2::6] = sec[f"NIS {l}"][f"NER {j}"]["LRU = 2 LFW = 0 LRF = 1"][f"List {k}"]["AMUN"]
                        add[3::6] = sec[f"NIS {l}"][f"NER {j}"]["LRU = 2 LFW = 0 LRF = 1"][f"List {k}"]["GNO"]
                        add[4::6] = sec[f"NIS {l}"][f"NER {j}"]["LRU = 2 LFW = 0 LRF = 1"][f"List {k}"]["GG"]
                        lines += sandy.write_list(
                            sec[f"NIS {l}"][f"NER {j}"]["LRU = 2 LFW = 0 LRF = 1"][f"List {k}"]["AWRI"],
                            0,
                            sec[f"NIS {l}"][f"NER {j}"]["LRU = 2 LFW = 0 LRF = 1"][f"List {k}"]["L"],
                            0,
                            sec[f"NIS {l}"][f"NER {j}"]["LRU = 2 LFW = 0 LRF = 1"][f"List {k}"]["NJS"],
                            add,
                            )
                if LRF == 2:
                    lines += sandy.write_cont(
                        sec[f"NIS {l}"][f"NER {j}"]["LRU = 2 LRF = 2"]["SPI"],
                        sec[f"NIS {l}"][f"NER {j}"]["LRU = 2 LRF = 2"]["AP"],
                        sec[f"NIS {l}"][f"NER {j}"]["LRU = 2 LRF = 2"]["LSSF"],
                        0,
                        sec[f"NIS {l}"][f"NER {j}"]["LRU = 2 LRF = 2"]["NLS"],
                        0,
                        )
                    NLS = int(sec[f"NIS {l}"][f"NER {j}"]["LRU = 2 LRF = 2"]["NLS"])
                    for m in range(NLS) :
                        lines += sandy.write_cont(
                            sec[f"NIS {l}"][f"NER {j}"]["LRU = 2 LRF = 2"][f"NLS {m}"]["AWRI"],
                            0,
                            sec[f"NIS {l}"][f"NER {j}"]["LRU = 2 LRF = 2"][f"NLS {m}"]["L"],
                            0,
                            sec[f"NIS {l}"][f"NER {j}"]["LRU = 2 LRF = 2"][f"NLS {m}"]["NJS"],
                            0,
                            )
                        NJS = int(sec[f"NIS {l}"][f"NER {j}"]["LRU = 2 LRF = 2"][f"NLS {m}"]["NJS"])
                        for k in range(NJS):
                            S_NJS = len (sec[f"NIS {l}"][f"NER {j}"]["LRU = 2 LRF = 2"][f"NJS {m, k}"][f"List {k}"]["ES"])
                            add1 = [0]*6
                            add2 = [0]*(6*S_NJS)
                            add1[2] = sec[f"NIS {l}"][f"NER {j}"]["LRU = 2 LRF = 2"][f"NJS {m, k}"][f"List {k}"]["AMUX"]
                            add1[3] = sec[f"NIS {l}"][f"NER {j}"]["LRU = 2 LRF = 2"][f"NJS {m, k}"][f"List {k}"]["AMUN"]
                            add1[4] = sec[f"NIS {l}"][f"NER {j}"]["LRU = 2 LRF = 2"][f"NJS {m, k}"][f"List {k}"]["AMUG"]
                            add1[5] = sec[f"NIS {l}"][f"NER {j}"]["LRU = 2 LRF = 2"][f"NJS {m, k}"][f"List {k}"]["AMUF"]
                            add2[0::6] = sec[f"NIS {l}"][f"NER {j}"]["LRU = 2 LRF = 2"][f"NJS {m, k}"][f"List {k}"]["ES"]
                            add2[1::6] = sec[f"NIS {l}"][f"NER {j}"]["LRU = 2 LRF = 2"][f"NJS {m, k}"][f"List {k}"]["D"]
                            add2[2::6] = sec[f"NIS {l}"][f"NER {j}"]["LRU = 2 LRF = 2"][f"NJS {m, k}"][f"List {k}"]["GX"]
                            add2[3::6] = sec[f"NIS {l}"][f"NER {j}"]["LRU = 2 LRF = 2"][f"NJS {m, k}"][f"List {k}"]["GNO"]
                            add2[4::6] = sec[f"NIS {l}"][f"NER {j}"]["LRU = 2 LRF = 2"][f"NJS {m, k}"][f"List {k}"]["GG"]
                            add2[5::6] = sec[f"NIS {l}"][f"NER {j}"]["LRU = 2 LRF = 2"][f"NJS {m, k}"][f"List {k}"]["GF"]
                            add = add1 + add2
                            lines += sandy.write_list(
                                sec[f"NIS {l}"][f"NER {j}"]["LRU = 2 LRF = 2"][f"NJS {m, k}"]["AJ"],
                                0,
                                sec[f"NIS {l}"][f"NER {j}"]["LRU = 2 LRF = 2"][f"NJS {m, k}"]["INT"],
                                0,
                                sec[f"NIS {l}"][f"NER {j}"]["LRU = 2 LRF = 2"][f"NJS {m, k}"]["NE"],
                                add,
                                )
    return "\n".join(sandy.write_eol(lines, sec["MAT"], 2, 151))
