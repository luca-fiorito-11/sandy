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
import pandas as pd

__author__ = "Rayan HADDAD"
__all__ = [
    "read_mf2",
    "write_mf2"
]

mf = 2
mt = 151

pd.options.display.float_format = '{:.5e}'.format


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
     'ZA': 90233.0,
     'AWR': 231.04,
     'NIS': {90233: {'ABN': 1.0,
       'LFW': 0,
       'NER': {(1e-05, 1.9): {'LRU': 0,
         'LRF': 0,
         'NRO': 0,
         'NAPS': 0,
         'SPI': 0.5,
         'AP': 0.9765,
         'NLS': 0}}}}}
   


    Resonance parameters of the Thorium 230
    LRU = 1 LRF = 2
    >>> tape = sandy.get_endf6_file("jeff_33", "xs", 902300)
    >>> dic = sandy.read_mf2(tape, 9034)
    >>> print (dic["NIS"][90230]["NER"][(1e-05, 251.0)]["L"][0])
    {'AWRI': 228.06, 'QX': 0.0, 'L': 0, 'LRX': 0, 'NRS': 22,
     'RES_PAR': [{'ER': -1.427, 'AJ': 0.5, 'GT': 0.0251, 'GN': 0.0002, 'GG': 0.0249, 'GF': 0.0},
                 {'ER': 1.427, 'AJ': 0.5, 'GT': 0.0252, 'GN': 0.0003, 'GG': 0.0249, 'GF': 0.0},
                 {'ER': 17.27, 'AJ': 0.5, 'GT': 0.0351, 'GN': 0.0131, 'GG': 0.022, 'GF': 0.0},
                 {'ER': 23.84, 'AJ': 0.5, 'GT': 0.0371, 'GN': 0.0111, 'GG': 0.026, 'GF': 0.0},
                 {'ER': 32.2, 'AJ': 0.5, 'GT': 0.0323, 'GN': 0.0033, 'GG': 0.029, 'GF': 0.0},
                 {'ER': 39.8, 'AJ': 0.5, 'GT': 0.0375, 'GN': 0.0085, 'GG': 0.029, 'GF': 0.0},
                 {'ER': 48.1, 'AJ': 0.5, 'GT': 0.039, 'GN': 0.01, 'GG': 0.029, 'GF': 0.0},
                 {'ER': 64.5, 'AJ': 0.5, 'GT': 0.0287, 'GN': 0.0031, 'GG': 0.0256, 'GF': 0.0},
                 {'ER': 75.6, 'AJ': 0.5, 'GT': 0.0283, 'GN': 0.0027, 'GG': 0.0256, 'GF': 0.0},
                 {'ER': 83.3, 'AJ': 0.5, 'GT': 0.0504, 'GN': 0.0248, 'GG': 0.0256, 'GF': 0.0},
                 {'ER': 103.0, 'AJ': 0.5, 'GT': 0.0308, 'GN': 0.0052, 'GG': 0.0256, 'GF': 0.0},
                 {'ER': 116.1, 'AJ': 0.5, 'GT': 0.0659, 'GN': 0.0403, 'GG': 0.0256, 'GF': 0.0},
                 {'ER': 133.8, 'AJ': 0.5, 'GT': 0.0331, 'GN': 0.0075, 'GG': 0.0256, 'GF': 0.0},
                 {'ER': 139.0, 'AJ': 0.5, 'GT': 0.028, 'GN': 0.0024, 'GG': 0.0256, 'GF': 0.0},
                 {'ER': 148.2, 'AJ': 0.5, 'GT': 0.0312, 'GN': 0.0056, 'GG': 0.0256, 'GF': 0.0},
                 {'ER': 171.4, 'AJ': 0.5, 'GT': 0.0509, 'GN': 0.0253, 'GG': 0.0256, 'GF': 0.0},
                 {'ER': 184.2, 'AJ': 0.5, 'GT': 0.0489, 'GN': 0.0233, 'GG': 0.0256, 'GF': 0.0},
                 {'ER': 195.0, 'AJ': 0.5, 'GT': 0.0709, 'GN': 0.0453, 'GG': 0.0256, 'GF': 0.0},
                 {'ER': 209.0, 'AJ': 0.5, 'GT': 0.1009, 'GN': 0.0753, 'GG': 0.0256, 'GF': 0.0},
                 {'ER': 226.0, 'AJ': 0.5, 'GT': 0.0484, 'GN': 0.0228, 'GG': 0.0256, 'GF': 0.0},
                 {'ER': 241.0, 'AJ': 0.5, 'GT': 0.0311, 'GN': 0.0055, 'GG': 0.0256, 'GF': 0.0},
                 {'ER': 248.0, 'AJ': 0.5, 'GT': 0.0809, 'GN': 0.0553, 'GG': 0.0256, 'GF': 0.0}]}


    Resonance parameters of the iron 58
    LRU = 1 LRF = 3
    >>> tape = sandy.get_endf6_file("endfb_80", "xs", 260580)
    >>> dic = sandy.read_mf2(tape, 2637)
    >>> print (dic["NIS"][26058]["NER"][(1e-05, 350000.0)]["L"][0]["RES_PAR"][20])
    {'ER': -3612.8, 'AJ': 0.5, 'GN': 271.55, 'GG': 0.74169, 'GFA': 0.0, 'GFB': 0.0}



    Resonance parameters of the iron 58
    LRU = 2 LFW = 0 LRF = 1
    >>> tape = sandy.get_endf6_file("endfb_80", "xs", 260580)
    >>> dic = sandy.read_mf2(tape, 2637)
    >>> print (dic["NIS"][26058]["NER"][(350000.0, 3000000.0)]["L"][3])
    {'AWRI': 57.43561, 'L': 3, 'NJS': 2,
     'RES_PAR': [{'D': 8033.33, 'AJ': 2.5, 'AMUN': 1.0, 'GN0': 0.482, 'GG': 0.7},
                 {'D': 6025.0, 'AJ': 3.5, 'AMUN': 1.0, 'GN0': 0.3615, 'GG': 0.7}]}


    Resonance parameters of the iron 54
    LRU = 1 LRF = 7
    >>> tape = sandy.get_endf6_file("endfb_80", "xs", 260540)
    >>> dic = sandy.read_mf2(tape, 2625)
    >>> print (dic["NIS"][26054]["NER"][(1e-05, 1036000.0)]["J"][2]["SPIN_GROUP"])
    [{'PPI': 1.0, 'L': 0.0, 'SCH': 0.0, 'BND': 0.0, 'APE': 0.0, 'APT': 0.0},
     {'PPI': 2.0, 'L': 1.0, 'SCH': 0.5, 'BND': 0.0, 'APE': 0.54373, 'APT': 0.54373}]


    Resonance parameters of the Uranium 235
    LRU = 2 LRF = 2
    >>> tape = sandy.get_endf6_file("jeff_33", "xs", 922350)
    >>> dic = sandy.read_mf2(tape, 9228)
    >>> print (dic["NIS"][92235]["NER"][(2250.0, 46200.0)]['L'][0]["J"][0]["RES_PAR"][0])
    {'ES': 2250.0, 'D': 1.058, 'GX': 0.0, 'GN0': 0.000107789, 'GG': 0.038513, 'GF': 0.40102}

    """
    df = tape._get_section_df(mat, mf, mt)
    out = {
        "MAT": mat,
        "MF": mf,
        "MT": mt,
    }
    i = 0
    C, i = sandy.read_cont(df, i)
    add = {
        "ZA": C.C1,  # designation for an isotope
        "AWR": C.C2,  # AWR is defines as the ratio of the mass of the material to that of the neutron
    }
    NIS = int(C.N1)

    out.update(add)
    P = {}
    for l in range(NIS):
        M = {}
        NER1 = {}
        ISO = {}
        C, i = sandy.read_cont(df, i)
        header1 = {
            "ABN": C.C2,  # Abundance of an isotope in the material
            "LFW": C.L2,  # indication whether average fission wifths are given in the unresolbed resonance region
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
                # Flag indicating whether this energy range contains data for
                # resolved or unresolved resonance parameters:
                "LRU": C.L1,
                # Flag indicating which representation has been used for the
                # energy range.
                "LRF": C.L2,
                "NRO": C.N1,  # Flag designating possible energy dependence of the scattering radiu
                "NAPS": C.N2,  # Flag controlling the use of the two radii
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
                    "SPI": C.C1,  # Spin, I, of the target nucleus.
                    "AP": C.C2,  # Scattering radius in units of 10e-12cm.
                    # Number of l-values (neutron orbital angular momentum) in
                    # this energy region.
                    "NLS": C.N1,
                }
                LRU0.update(add)
                dico[(EL, EH)] = LRU0
            elif LRU == 1:
                if LRF == 1 or LRF == 2:
                    if NRO == 0:
                        LRU1_LRF1_2_NRO0 = {}
                        LRU1_LRF1_2_NRO0.update(header2)
                        C, i = sandy.read_cont(df, i)
                        add = {
                            "SPI": C.C1,
                            "AP": C.C2,
                        }
                        NLS = int(C.N1)
                        LRU1_LRF1_2_NRO0.update(add)
                        LRU1_LRF1_2_NRO0_NLS = {}
                        for k in range(NLS):
                            L, i = sandy.read_list(df, i)
                            add = {
                                "AWRI": L.C1,  # Ratio of the mass of a particular isotope to that of a neutron
                                # Q-value to be added to the incident
                                # particle’s center-of-mass energy to determine
                                # the channel energy for use in the
                                # penetrability factor.
                                "QX": L.C2,
                                "L": L.L1,  # Value of l.
                                "LRX": L.L2,  # Flag indicating whether this energy range contains a competitive width
                                # Number of resolved resonances for a given
                                # l-value.
                                "NRS": L.N2,
                            }
                            keys = ["ER", "AJ", "GT", "GN", "GG", "GF"]
                            RES_PAR = [dict(zip(keys, items))
                                       for items in sandy.utils.grouper(L.B, 6)]
                            add.update({"RES_PAR": RES_PAR})
                            LRU1_LRF1_2_NRO0_NLS.update({k: add})
                            LRU1_LRF1_2_NRO0.update(
                                {"L": LRU1_LRF1_2_NRO0_NLS})
                        dico[(EL, EH)] = LRU1_LRF1_2_NRO0
                    else:
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
                        }
                        NLS = int(C.N1)
                        LRU1_LRF1_2_NRO1.update(add)
                        for k in range(NLS):
                            L, i = sandy.read_list(df, i)
                            add = {
                                "AWRI": L.C1,  # Ratio of the mass of a particular isotope to that of a neutron
                                # Q-value to be added to the incident
                                # particle’s center-of-mass energy to determine
                                # the channel energy for use in the
                                # penetrability factor.
                                "QX": L.C2,
                                "L": L.L1,  # Value of l.
                                "LRX": L.L2,  # Flag indicating whether this energy range contains a competitive width
                                "NRS": L.N2,
                            }
                            keys = ["ER", "AJ", "GT", "GN", "GG", "GF"]
                            RES_PAR = [dict(zip(keys, items))
                                       for items in sandy.utils.grouper(L.B, 6)]
                            add.update({"RES_PAR": RES_PAR})
                        LRU1_LRF1_2_NRO1.update({"L": add})
                        dico[(EL, EH)] = LRU1_LRF1_2_NRO1
                elif LRF == 3:
                    if NRO == 0:
                        LRU1_LRF3_NRO0 = {}
                        LRU1_LRF3_NRO0.update(header2)
                        C, i = sandy.read_cont(df, i)
                        add = {
                            "SPI": C.C1,
                            "AP": C.C2,
                            # Flag indicating whether these parameters can be
                            # used to compute angular distributions.
                            "LAD": C.L1,
                            # Number of l-values which must be used to converge
                            # the calculation with respect to the incident
                            # l-value in order to obtain accurate elastic
                            # angular distributions.
                            "NLSC": C.N2,
                        }
                        NLS = int(C.N1)
                        LRU1_LRF3_NRO0.update(add)
                        LRU1_LRF3_NRO0_NLS = {}
                        for k in range(NLS):
                            L, i = sandy.read_list(df, i)
                            add = {
                                "AWRI": L.C1,  
                                "APL": L.C2,
                                "L": L.L1,  # Value of l.
                                # Number of resolved resonances for a given
                                # l-value.
                                "NRS": L.N2,
                            }

                            keys = ["ER", "AJ", "GN", "GG", "GFA", "GFB"]
                            RES_PAR = [dict(zip(keys, items))
                                       for items in sandy.utils.grouper(L.B, 6)]
                            add.update({"RES_PAR": RES_PAR})
                            LRU1_LRF3_NRO0_NLS.update({k: add})
                        LRU1_LRF3_NRO0.update({"L": LRU1_LRF3_NRO0_NLS})
                        dico[(EL, EH)] = LRU1_LRF3_NRO0
                    else:
                        LRU1_LRF3_NRO1 = {}
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
                            "NLSC": C.N2
                        }
                        NLS = int(C.N1)
                        LRU1_LRF3_NRO1.update(add)
                        LRU1_LRF3_NRO1_NLS = {}
                        for k in range(NLS):
                            L, i = sandy.read_list(df, i)
                            add = {
                                "AWRI": L.C1, 
                                "APL": L.C2,
                                "L": L.L1,  # Value of l.
                                # Number of resolved resonances for a given
                                # l-value.
                                "NRS": L.N2,
                            }

                            keys = ["ER", "AJ", "GN", "GG", "GFA", "GFB"]
                            RES_PAR = [dict(zip(keys, items))
                                       for items in sandy.utils.grouper(L.B, 6)]
                            add.update({"RES_PAR": RES_PAR})
                            LRU1_LRF3_NRO1_NLS.update({k: add})
                        LRU1_LRF3_NRO1.update({"L": LRU1_LRF3_NRO1_NLS})
                        dico[(EL, EH)] = LRU1_LRF3_NRO1
                elif LRF == 4:
                    raise ValueError("LRF = 4 is not supported in SANDY")
                elif LRF == 7:
                    LRU1_LRF7 = {}
                    LRU1_LRF7.update(header2)
                    C, i = sandy.read_cont(df, i)
                    add = {
                        "IFG": C.L1,
                        "KRM": C.L2,  # Flag to specify which formulae for the R-matrix are to be used
                        "KRL": C.N2,  # Flag is zero for non-relativistic kinematics
                    }
                    NJS = int(C.N1)
                    LRU1_LRF7.update(add)
                    L, i = sandy.read_list(df, i)
                    add = {
                        "NPP": L.L1,  # Total number of particle-pairs.
                    }
                    keys = [
                        "MA",
                        "MB",
                        "ZA",
                        "ZB",
                        "IA",
                        "IB",
                        "Q",
                        "PNT",
                        "SHF",
                        "MT",
                        "PA",
                        "PB"]
                    PAIR_PART = [dict(zip(keys, items))
                                 for items in sandy.utils.grouper(L.B, 12)]
                    add.update({"PAIR_PART": PAIR_PART})
                    LRU1_LRF7.update(add)
                    LRU1_LRF7_NJS = {}
                    for k in range(NJS):
                        LISTS = {}
                        L, i = sandy.read_list(df, i)
                        add1 = {
                            # Floating point value of J (spin); sign indicates
                            # parity
                            "AJ": L.C1,
                            "PJ": L.C2,  # Parity (used only if AJ = 0.0).
                            "KBK": L.L1,  # Non-zero if background R-matrix exists
                            # Non-zero if non-hard-sphere phase shift are to be
                            # specified.
                            "KPS": L.L2,
                            # Number of channels for the given J pi.
                            "NCH": L.N2,
                        }
                        NCH = int(L.N2)
                        keys = ["PPI", "L", "SCH", "BND", "APE", "APT"]
                        SPIN_GROUP = [dict(zip(keys, items))
                                      for items in sandy.utils.grouper(L.B, 6)]
                        add1.update({"SPIN_GROUP": SPIN_GROUP})

                        L, i = sandy.read_list(df, i)
                        add2 = {
                            "NRS": L.L2,  # Number of resonances for the given J pi
                            "NX": L.N2,  # Number of lines required for all resonances for the given J
                        }

                        LIST2 = {}
                        Ep = L.B[::NCH + 1]
                        b = L.B[::1]
                        del b[::NCH + 1]
                        add_3 = {
                            "ER": Ep,  # Resonance energy in eV
                            "GAM": b,  # Channel width in eV
                        }
                        LIST2.update(add_3)
                        add2.update(LIST2)
                        LISTS.update(add1)
                        LISTS.update(add2)
                        LRU1_LRF7_NJS.update({k: LISTS})
                    LRU1_LRF7.update({"J": LRU1_LRF7_NJS})
                    dico[(EL, EH)] = LRU1_LRF7
            elif LRU == 2:
                if LFW == 0 and LRF == 1:
                    LRU2_LFW0_LRF1 = {}
                    LRU2_LFW0_LRF1.update(header2)
                    C, i = sandy.read_cont(df, i)
                    add = {
                        "SPI": C.C1,
                        "AP": C.C2,
                        "LSSF": C.L1,  # Flag governing the interpretation of the File 3 cross sections
                    }
                    NLS = int(C.N1)
                    LRU2_LFW0_LRF1.update(add)
                    LRU2_LFW0_LRF1_NLS = {}
                    for k in range(NLS):
                        L, i = sandy.read_list(df, i)
                        add = {
                            "AWRI": L.C1,
                            "L": L.L1,
                            "NJS": L.N2,  # Number of J-states for a particular l-state
                        }
                        keys = ["D", "AJ", "AMUN", "GN0", "GG", ]
                        RES_PAR = [dict(zip(keys, items))
                                   for items in sandy.utils.grouper(L.B, 6)]
                        add.update({"RES_PAR": RES_PAR})
                        LRU2_LFW0_LRF1_NLS.update({k: add})
                    LRU2_LFW0_LRF1.update({"L": LRU2_LFW0_LRF1_NLS})
                    dico[(EL, EH)] = LRU2_LFW0_LRF1
                elif LFW == 1 and LRF == 1:
                    raise ValueError(
                        "LFW = 1  LRF = 1 is not supported in SANDY")
                elif LRF == 2:
                    LRU2_LRF2 = {}
                    LRU2_LRF2.update(header2)
                    C, i = sandy.read_cont(df, i)
                    add = {
                        "SPI": C.C1,
                        "AP": C.C2,
                        "LSSF": C.L1,
                    }
                    NLS = int(C.N1)
                    LRU2_LRF2.update(add)
                    LRU2_LRF2_NLS = {}
                    for m in range(NLS):
                        C, i = sandy.read_cont(df, i)
                        add_1 = {
                            "AWRI": C.C1,
                            "L": C.L1,
                        }
                        NJS = int(C.N1)
                        LRU2_LRF2_NjS = {}
                        LIST = {}
                        for k in range(NJS):
                            L, i = sandy.read_list(df, i)
                            add_2 = {
                                "AJ": L.C1,
                                # Interpolation scheme to be used for
                                # interpolating between the cross sections
                                # obtained from average resonance parameters.
                                "INT": L.L1,
                                "NE": L.N2,
                                # Number of degrees of freedom used in the
                                # competitive width distribution.
                                "AMUX": L.B[2],
                                # Number of degrees of freedom in the neutron
                                # width distribution.
                                "AMUN": L.B[3],
                                # Number of degrees of freedom in the radiation
                                # width distribution.
                                "AMUG": L.B[4],
                                # Integer value of the number of degrees of
                                # freedom for fission widths.
                                "AMUF": L.B[5],
                            }
                            keys = ["ES", "D", "GX", "GN0", "GG", "GF"]
                            RES_PAR = [dict(zip(keys, items))
                                       for items in sandy.utils.grouper(L.B[6::], 6)]
                            add_2.update({"RES_PAR": RES_PAR})
                            LIST.update({k: add_2})
                        LRU2_LRF2_NjS.update(add_1)
                        LRU2_LRF2_NjS.update({"J": LIST})
                        LRU2_LRF2_NLS.update({m: LRU2_LRF2_NjS})
                    LRU2_LRF2.update({"L": LRU2_LRF2_NLS})
                    dico[(EL, EH)] = LRU2_LRF2
            NER1.update(dico)
            info.update({"NER": NER1})
            M.update(info)
            ISO.update({ZAI: M})
        P.update(ISO)
    out.update({"NIS": P})
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
        sec["ZA"],
        sec["AWR"],
        0,
        0,
        len(sec["NIS"]),
        0,
    )
    for ZAI, sec2 in sec["NIS"].items():
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
                EL,
                EH,
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
                if sec3["LRF"] == 1 or sec3["LRF"] == 2:
                    if sec3["NRO"] == 0:
                        lines += sandy.write_cont(
                            sec3["SPI"],
                            sec3["AP"],
                            0,
                            0,
                            len(sec3["L"]),
                            0,
                        )
                        for l, sec4 in sec3["L"].items():
                            keys = ["ER", "AJ", "GT", "GN", "GG", "GF"]
                            tab = [res[k] for res in sec4["RES_PAR"]
                                   for k in keys]
                            lines += sandy.write_list(
                                sec4["AWRI"],
                                sec4["QX"],
                                sec4["L"],
                                sec4["LRX"],
                                sec4["NRS"],
                                tab,
                            )
                    else:
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
                            len(sec3["L"]),
                            0,
                        )
                        for l, sec4 in sec3["L"].items():
                            keys = ["ER", "AJ", "GT", "GN", "GG", "GF"]
                            tab = [res[k] for res in sec4["RES_PAR"]
                                   for k in keys]
                            lines += sandy.write_list(
                                sec4["AWRI"],
                                sec4["QX"],
                                sec4["L"],
                                sec4["LRX"],
                                sec4["NRS"],
                                tab,
                            )
                elif sec3["LRF"] == 3:
                    if sec3["NRO"] == 0:
                        lines += sandy.write_cont(
                            sec3["SPI"],
                            sec3["AP"],
                            sec3["LAD"],
                            0,
                            len(sec3["L"]),
                            sec3["NLSC"],
                        )
                        for l, sec4 in sec3["L"].items():
                            keys = ["ER", "AJ", "GN", "GG", "GFA", "GFB"]
                            tab = [res[k] for res in sec4["RES_PAR"]
                                   for k in keys]
                            lines += sandy.write_list(
                                sec4["AWRI"],
                                sec4["APL"],
                                sec4["L"],
                                0,
                                sec4["NRS"],
                                tab,
                            )
                    else:
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
                            len(sec3["L"]),
                            sec3["NLSC"],
                        )
                        for l, sec4 in sec3["L"].items():
                            keys = ["ER", "AJ", "GN", "GG", "GFA", "GFB"]
                            tab = [res[k] for res in sec4["RES_PAR"]
                                   for k in keys]
                            lines += sandy.write_list(
                                sec4["AWRI"],
                                sec4["APL"],
                                sec4["L"],
                                0,
                                sec4["NRS"],
                                tab,
                            )
                elif sec3["LRF"] == 7:
                    lines += sandy.write_cont(
                        0,
                        0,
                        sec3["IFG"],
                        sec3["KRM"],
                        len(sec3["J"]),
                        sec3["KRL"],
                    )
                    keys = [
                        "MA",
                        "MB",
                        "ZA",
                        "ZB",
                        "IA",
                        "IB",
                        "Q",
                        "PNT",
                        "SHF",
                        "MT",
                        "PA",
                        "PB"]
                    tab = [res[k] for res in sec3["PAIR_PART"] for k in keys]
                    lines += sandy.write_list(
                        0,
                        0,
                        sec3["NPP"],
                        0,
                        2 * sec3["NPP"],
                        tab,
                    )
                    for l, sec4 in sec3["J"].items():
                        keys = ["PPI", "L", "SCH", "BND", "APE", "APT"]
                        tab = [res[k] for res in sec4["SPIN_GROUP"]
                               for k in keys]
                        lines += sandy.write_list(
                            sec4["AJ"],
                            sec4["PJ"],
                            sec4["KBK"],
                            sec4["KPS"],
                            sec4["NCH"],
                            tab,
                        )
                        NCH = sec4["NCH"]
                        size_E = len(sec4["ER"])
                        size_b = len(sec4["GAM"])
                        add2 = []
                        n = NCH
                        x = [sec4["GAM"][i:i + n] for i in range(0, size_b, n)]
                        for i in range(size_E):
                            add2.append(sec4["ER"][i])
                            for m in range(len(x[i])):
                                add2.append(x[i][m])
                        lines += sandy.write_list(
                            0,
                            0,
                            0,
                            sec4["NRS"],
                            sec4["NX"],
                            add2,
                        )
            elif sec3["LRU"] == 2:
                if sec2["LFW"] == 0 and sec3["LRF"] == 1:
                    lines += sandy.write_cont(
                        sec3["SPI"],
                        sec3["AP"],
                        sec3["LSSF"],
                        0,
                        len(sec3["L"]),
                        0,
                    )
                    for l, sec4 in sec3["L"].items():
                        keys = ["D", "AJ", "AMUN", "GN0", "GG", ]
                        tab = [res[k] for res in sec4["RES_PAR"] for k in keys]
                        lines += sandy.write_list(
                            sec4["AWRI"],
                            0,
                            sec4["L"],
                            0,
                            sec4["NJS"],
                            tab,
                        )
                elif sec3["LRF"] == 2:
                    lines += sandy.write_cont(
                        sec3["SPI"],
                        sec3["AP"],
                        sec3["LSSF"],
                        0,
                        len(sec3["L"]),
                        0,
                    )
                    for m, sec4 in sec3["L"].items():
                        lines += sandy.write_cont(
                            sec4["AWRI"],
                            0,
                            sec4["L"],
                            0,
                            len(sec4["J"]),
                            0,
                        )
                        for k, sec5 in sec4["J"].items():
                            keys = ["ES", "D", "GX", "GN0", "GG", "GF"]
                            tab = [res[k] for res in sec5["RES_PAR"]
                                   for k in keys]
                            add1 = [0] * 6
                            add1[2] = sec5["AMUX"]
                            add1[3] = sec5["AMUN"]
                            add1[4] = sec5["AMUG"]
                            add1[5] = sec5["AMUF"]
                            add = add1 + tab
                            lines += sandy.write_list(
                                sec5["AJ"],
                                0,
                                sec5["INT"],
                                0,
                                sec5["NE"],
                                add,
                            )
    return "\n".join(sandy.write_eol(lines, sec["MAT"], 2, 151))
