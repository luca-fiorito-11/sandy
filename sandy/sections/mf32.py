"""
This module contains only two public functions:

    * `read_mf32`

Function `read` reads a MF32/MT section from a string and produces a content
object with a dictionary-like structure.
The content object can be accessed using most of the keywords specified in
the ENDF6 manual for this specific MF section.

Function `write_mf32` writes a content object for a MF32/MT section into a
string.
MAT, MF, MT and line numbers are also added (each line ends with a `\n`).
"""

import sandy
from collections import namedtuple
import numpy as np
import pandas as pd
import io

__author__ = "Rayan HADDAD"
__all__ = [
    "read_mf32",
]

mf = 32
mt = 151

pd.options.display.float_format = '{:.5e}'.format

def read_intg(tape, mat, NDIGIT,ipos):
    """
    Read ``ENDF-6`` ``INTG`` record in formatted fortran.

    Returns
    -------
    `CONT` : `collections.namedtuple`
        - `II` : `int`
            first five elements of string
        - `JJ` : `int`
            second five elements of string
        - `KIJ` : `array`
            depending of the format imposed by NDIGIT
            for NDIGIT = 2 KIJ is an array of 18 str 
            for NDIGIT = 3 KIJ is an array of 13 str 
            for NDIGIT = 3 KIJ is an array of 11 str 
            for NDIGIT = 4 KIJ is an array of 9 str 
            for NDIGIT = 5 KIJ is an array of 8 str 
            
    """
    INTG = namedtuple('INTG', 'II JJ KIJ')
    
    if NDIGIT == 2:
        
        nparam = 18
        widths = [5] * 2 + [1] + [3] * nparam
        names = range (21)
        line = tape.data[(mat, 32, 151)].splitlines()[ipos]
        df = pd.read_fwf(io.StringIO(line), widths=widths, dtype=str, na_filter=True, names=names).fillna(0)
        II = df.iloc[:, 0].astype(int).squeeze()
        JJ = df.iloc[:, 1].astype(int).squeeze()
        KIJ = df.iloc[:, 3:].astype(int).squeeze().values
        
    elif NDIGIT == 3:
        nparam = 13
        widths = [5] * 2 + [1] + [4] * nparam
        names = range (nparam + 3)
        line = tape.data[(mat, 32, 151)].splitlines()[ipos]
        df = pd.read_fwf(io.StringIO(line), widths=widths, dtype=str, na_filter=True, names=names).fillna(0)
        II = df.iloc[:, 0].astype(int).squeeze()
        JJ = df.iloc[:, 1].astype(int).squeeze()
        KIJ = df.iloc[:, 3:].astype(int).squeeze().values
        
    elif NDIGIT == 4:
        nparam = 11
        widths = [5] * 2 + [1] + [5] * nparam
        names = range (nparam + 3)
        line = tape.data[(mat, 32, 151)].splitlines()[ipos]
        df = pd.read_fwf(io.StringIO(line), widths=widths, dtype=str, na_filter=True, names=names).fillna(0)
        II = df.iloc[:, 0].astype(int).squeeze()
        JJ = df.iloc[:, 1].astype(int).squeeze()
        KIJ = df.iloc[:, 3:].astype(int).squeeze().values
        
    elif NDIGIT == 5:
        nparam = 9
        widths = [5] * 2 + [1] + [6] * nparam
        names = range (nparam + 3)
        line = tape.data[(mat, 32, 151)].splitlines()[ipos]
        df = pd.read_fwf(io.StringIO(line), widths=widths, dtype=str, na_filter=True, names=names).fillna(0)
        II = df.iloc[:, 0].astype(int).squeeze()
        JJ = df.iloc[:, 1].astype(int).squeeze()
        KIJ = df.iloc[:, 3:].astype(int).squeeze().values
        
    elif NDIGIT == 6:
        nparam = 8
        widths = [5] * 2 + [7] * nparam
        names = range (nparam + 2)
        line = tape.data[(mat, 32, 151)].splitlines()[ipos]
        df = pd.read_fwf(io.StringIO(line), widths=widths, dtype=str, na_filter=True, names=names).fillna(0)
        II = df.iloc[:, 0].astype(int).squeeze()
        JJ = df.iloc[:, 1].astype(int).squeeze()
        KIJ = df.iloc[:, 3:].astype(int).squeeze().values
    ipos += 1
    
    return INTG(II, JJ, KIJ), ipos

    


def read_mf32(tape, mat):
    """
    Write MT section for MF32

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
    Covariances of resonance parameters of the Curium 2245
    LCOMP = 0
    >>> tape = sandy.get_endf6_file("jeff_33", "xs", 962450)
    >>> dic = sandy.read_mf32(tape, 9640)
    >>> print( dic["NIS"][96245]['NER'][(1e-05, 100.0)]["L"][0]['COVAR_PAR'][0:2])
    [{'ER': -0.1, 'AJ': 3.0, 'GT': 0.230946, 'GN': 4.61e-05, 'GG': 0.0359, 'GF': 0.195,
      'DE²': 1e-08, 'DN²': 2.12521e-11, 'DNDG': 0.0, 'DG²': 5.15524e-05, 'DNDF': 0.0, 'DGDF': 0.0, 
      'DF²': 0.00038025, 'DJDN': 0.0, 'DJDG': 0.0, 'DJDF': 0.0, 'DJ²': 0.0}, 
     {'ER': 0.85, 'AJ': 4.0, 'GT': 0.84409, 'GN': 9e-05, 'GG': 0.044, 'GF': 0.8,
      'DE²': 0.0009, 'DN²': 6.30623e-11, 'DNDG': 0.0, 'DG²': 7.744e-05, 'DNDF': 0.0, 'DGDF': 0.0,
      'DF²': 0.0025, 'DJDN': 0.0, 'DJDG': 0.0, 'DJDF': 0.0, 'DJ²': 0.0}]
    
    Covariances of resonance parameters of the Americium 241
    LCOMP = 1 LRF = 2
    >>> tape = sandy.get_endf6_file("jeff_33", "xs", 952410)
    >>> dic = sandy.read_mf32(tape, 9543)
    >>> print(dic["NIS"][95241]['NER'][(1e-05, 150.0)]["COVAR_PAR"][0][10:50])
    [ 2.494710e-10  7.847580e-11 -9.012000e-10 -9.512310e-12 -2.029530e-10
      1.497910e-11  7.083150e-09 -1.567680e-11  1.527096e-09  7.121410e-11
      6.664252e-09 -1.968360e-11 -2.291857e-09  3.975690e-11  9.398297e-09
     -2.101620e-11 -4.383682e-09 -6.580010e-11  1.123283e-08 -2.688670e-11
     -2.484173e-09 -3.560260e-11 -3.078174e-08 -6.652760e-11 -1.722510e-08
      5.561780e-11 -4.403603e-08 -3.319870e-11 -1.722510e-08  1.024590e-10
      4.308922e-09  4.159470e-10 -1.722510e-08 -1.178230e-10 -8.631673e-08
     -1.642100e-11 -1.722510e-08  7.081580e-11 -9.513253e-08 -3.665940e-12]
    
    Covariances of resonance parameters of the Caesium 137
    LCOMP = 2 LRF = 2
    >>> tape = sandy.get_endf6_file("jeff_33", "xs", 551370)
    >>> dic = sandy.read_mf32(tape, 5537)
    >>> dic['NIS'][55137]["NER"][(1e-05, 56666.38)]['RES_PAR'][20::30]
    [{'ER': -1342.1, 'AJ': 4.0, 'GT': 10.75884, 'GN': 10.09489, 'GG': 0.66395, 'GF': 0.0, 
      'DER': 0.0005, 'DGN': 0.5, 'DGG': 0.1, 'DGF': 0.0},
     {'ER': 20098.0, 'AJ': 4.0, 'GT': 36.19939, 'GN': 36.163, 'GG': 0.036395, 'GF': 0.0, 
      'DER': 10.049, 'DGN': 0.00036163, 'DGG': 3.6395e-07, 'DGF': 0.0},
     {'ER': 36830.02, 'AJ': 4.0, 'GT': 3.388613, 'GN': 3.298613, 'GG': 0.09000001, 'GF': 0.0,
      'DER': 36.83, 'DGN': 0.0659723, 'DGG': 0.0018, 'DGF': 0.0},
     {'ER': 51650.52, 'AJ': 3.0, 'GT': 17.43718, 'GN': 17.34718, 'GG': 0.09000001, 'GF': 0.0,
      'DER': 51.6505, 'DGN': 0.346944, 'DGG': 0.0018, 'DGF': 0.0}]
   
    Covariances of resonance parameters of the Caesium 137
    LCOMP = 2 LRF = 3
    >>> tape = sandy.get_endf6_file("jeff_33", "xs", 902320)
    >>> dic = sandy.read_mf32(tape, 9040)
    >>> dic['NIS'][90232]["NER"][(1e-05, 4000.0)]["INTG"][350]
    {'II': 795,
     'JJ': 793,
     'KIJ': array([44, 69,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0])}
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
        NER = C.N1
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
            NRO = C.N1
            LRF = C.L2
            if NRO != 0:
                NRO1 = {}
                C, i = sandy.read_cont(df, i)
                add = {
                    "NI": C.N2,
                }
                NRO1.update(add)
                dico[(EL, EH)] = NRO1
            else:
                C, i = sandy.read_cont(df, i)
                header3 = {
                    "SPI": C.C1,  # Flag controlling the use of the two radii
                    "AP": C.C2,
                    "LCOMP": C.L2,
                    "NLS": C.N1,
                    "ISR": C.N2,
                }
                LCOMP = int(C.L2)
                ISR = int(C.N2)
                NLS = int(C.N1)
                if LCOMP == 0:
                    LCOMP0 = {}
                    LCOMP0.update(header2)
                    LCOMP0.update(header3)
                    if ISR > 0:
                        C, i = sandy.read_cont(df, i)
                        add = {
                            "DAP": C.C2,
                        }
                        LCOMP0.update(add)
                        LCOMP0_NLS = {}
                        for k in range(NLS):
                            L, i = sandy.read_list(df, i)
                            add = {
                                "AWRI": L.C1,  # Ratio of the mass of a particular isotope to that of a neutron
                                # Number of resolved resonances for a given
                                # l-value.
                                "NRS": L.N2,
                            }
                            keys = [
                                "ER",
                                "AJ",
                                "GT",
                                "GN",
                                "GG",
                                "GF",
                                "DE²",
                                "DN²",
                                "DNDG",
                                "DG²",
                                "DNDF",
                                "DGDF",
                                "DF²",
                                "DJDN",
                                "DJDG",
                                "DJDF",
                                "DJ²",
                            ]
                            COVAR_PAR = [dict(zip(keys, items))
                                         for items in sandy.utils.grouper(L.B, 18)]
                            add.update({"COVAR_PAR": COVAR_PAR})
                            LCOMP0_NLS.update({L.L1: add})
                            LCOMP0.update({"L": LCOMP0_NLS})
                        dico[(EL, EH)] = LCOMP0
                    else:
                        LCOMP0_NLS = {}
                        for k in range(NLS):
                            L, i = sandy.read_list(df, i)
                            add = {
                                "AWRI": L.C1,  # Ratio of the mass of a particular isotope to that of a neutron
                                # Number of resolved resonances for a given
                                # l-value.
                                "NRS": L.N2,
                            }
                            keys = [
                                "ER",
                                "AJ",
                                "GT",
                                "GN",
                                "GG",
                                "GF",
                                "DE²",
                                "DN²",
                                "DNDG",
                                "DG²",
                                "DNDF",
                                "DGDF",
                                "DF²",
                                "DJDN",
                                "DJDG",
                                "DJDF",
                                "DJ²",
                            ]
                            COVAR_PAR = [dict(zip(keys, items))
                                         for items in sandy.utils.grouper(L.B, 18)]
                            add.update({"COVAR_PAR": COVAR_PAR})
                            LCOMP0_NLS.update({L.L1: add})
                            LCOMP0.update({"L": LCOMP0_NLS})
                        dico[(EL, EH)] = LCOMP0
                elif LCOMP == 1:
                    LCOMP1 = {}
                    LCOMP1.update(header2)
                    LCOMP1.update(header3)
                    if LRF == 1 or LRF == 2:  # Breit-Wigner
                        if ISR > 0:
                            C, i = sandy.read_cont(df, i)
                            add = {
                                "DAP": C.C2,
                            }
                            LCOMP1.update(add)
                            C, i = sandy.read_cont(df, i)
                            add = {
                                "AWRI": C.C1,
                                "NSRS": C.N1,
                                "NLRS": C.N2,
                            }
                            LCOMP1.update(add)
                            L, i = sandy.read_list(df, i)
                            add = {
                                "MPAR": L.L1,
                                "NRB": L.N2,
                            }
                            NRB = int(L.N2)
                            MPAR = int(L.L1)
                            keys = ["ER", "AJ", "GT", "GN", "GG", "GF"]
                            RES_PAR = [dict(zip(keys, items)) for items in sandy.utils.grouper(
                                L.B[:6 * NRB], 6)]
                            add.update({"RES_PAR": RES_PAR})
                            tri = np.zeros((MPAR * NRB, MPAR * NRB))
                            tri[np.triu_indices(
                                MPAR * NRB, 0)] = np.array(L.B[6 * NRB:])
                            COVAR_PAR = {
                                "COVAR_PAR": tri
                            }
                            add.update(COVAR_PAR)
                            LCOMP1.update(add)
                            dico[(EL, EH)] = LCOMP1
                        if ISR == 0:
                            C, i = sandy.read_cont(df, i)
                            add = {
                                "AWRI": C.C1,
                                "NSRS": C.N1,
                                "NLRS": C.N2,
                            }
                            LCOMP1.update(add)
                            L, i = sandy.read_list(df, i)
                            add = {
                                "MPAR": L.L1,
                                "NRB": L.N2,
                            }
                            NRB = int(L.N2)
                            MPAR = int(L.L1)
                            keys = ["ER", "AJ", "GT", "GN", "GG", "GF"]
                            RES_PAR = [dict(zip(keys, items)) for items in sandy.utils.grouper(
                                L.B[:6 * NRB], 6)]
                            add.update({"RES_PAR": RES_PAR})
                            tri = np.zeros((MPAR * NRB, MPAR * NRB))
                            tri[np.triu_indices(
                                MPAR * NRB, 0)] = np.array(L.B[6 * NRB:])
                            COVAR_PAR = {
                                "COVAR_PAR": tri
                            }
                            add.update(COVAR_PAR)
                            LCOMP1.update(add)
                            dico[(EL, EH)] = LCOMP1
                    elif LRF == 3:  # Reich-Moore
                        add = {}
                        if ISR > 0:
                            L, i = sandy.read_list(df, i)
                            # DAP is uncertainty on scattering radius
                            DAP = [dict(zip(keys, items)) for items in sandy.utils.grouper(L.B, 1)]
                            add.update({
                                "MLS": L.NPL,
                                "DAP": DAP,
                            })
                        C, i = sandy.read_cont(df, i)
                        add.update({
                            "AWRI": C.C1,
                            "NSRS": C.N1,
                            "NLRS": C.N2,
                        })
                        L, i = sandy.read_list(df, i)
                        NRB = int(L.N2)
                        MPAR = int(L.L1)
                        keys = ["ER", "AJ", "GN", "GG", "GFA", "GFB"]
                        
                        # resonance parameters
                        RES_PAR = [dict(zip(keys, items)) for items in sandy.utils.grouper(L.B[:6 * NRB], 6)]
                        # resonance parameters covariance data
                        tri = np.zeros((MPAR * NRB, MPAR * NRB))
                        tri[np.triu_indices(MPAR * NRB, 0)] = np.array(L.B[6 * NRB:])

                        add.update({
                            "MPAR": L.L1,
                            "NRB": L.N2,
                            "RES_PAR": RES_PAR,
                            "COVAR_PAR": tri,
                            })

                        LCOMP1.update(add)

                        dico[(EL, EH)] = LCOMP1
                    elif LRF == 4:
                        L, i = sandy.read_list(df, i)
                        add = {
                            "MPAR": L.L1,
                            "NRB": L.N2,
                        }
                        NRB = int(L.N2)
                        MPAR = int(L.L1)
                        keys = [
                            "DET",
                            "DWT",
                            "GRT",
                            "GIT",
                            "DEF",
                            "DWF",
                            "GRF",
                            "GIF",
                            "DEC",
                            "DWC",
                            "GRC",
                            "GIC"]
                        RES_PAR = [dict(zip(keys, items))
                                   for items in sandy.utils.grouper(L.B[:6 * NRB], 12)]
                        add.update({"RES_PAR": RES_PAR})
                        tri = np.zeros((MPAR * NRB, MPAR * NRB))
                        tri[np.triu_indices(MPAR * NRB, 0)
                            ] = np.array(L.B[6 * NRB:])
                        COVAR_PAR = {
                            "COVAR_PAR": tri
                        }
                        add.update({"COVAR_PAR": COVAR_PAR})
                        LCOMP1.update(add)
                        dico[(EL, EH)] = LCOMP1
                    elif LRF == 7:
                        if ISR > 0:
                            L, i = sandy.read_list(df, i)
                            add = {
                                "JCH": L.NPL,
                                "(1+(NJCH-1)/6)": L.N2,
                                "DAP": L.B,
                            }
                            C, i = sandy.read_cont(df, i)
                            add = {
                                "NJSX": C.L1,
                            }
                            add.update(add)
                            NJSX = int(C.L1)
                            for k in range(NJSX):
                                L, i = sandy.read_list(df, i)
                                add = {
                                    "NCH": L.L1,
                                    "NRB": L.L2,
                                    "NX": L.N2,
                                }
                                NCH = int(L.L1)
                                LIST2 = {}
                                Ep = L.B[::NCH + 1]
                                b = L.B[::1]
                                del b[::NCH + 1]
                                add_3 = {
                                    "ER": Ep,  # Resonance energy in eV
                                    "GAM": b,  # Channel width in eV
                                }
                                LIST2.update({k: add_3})
                            add.update({"J": LIST2})
                            L, i = sandy.read_list(df, i)
                            add = {
                                "N": L.NPL,
                                "NPARB": L.N2,
                            }
                            add.update(add)
                            NPARB = int(L.N2)
                            tri = np.zeros((NPARB * NPARB, NPARB * NPARB))
                            tri[np.triu_indices(
                                NPARB * NPARB, 0)] = np.array(L.B[6 * NRB:])
                            COVAR_PAR = {
                                "COVAR_PAR": tri
                            }
                            add.update({"COVAR_PAR": COVAR_PAR})
                            LCOMP1.update(add)
                            dico[(EL, EH)] = LCOMP1
                        else:
                            C, i = sandy.read_cont(df, i)
                            add = {
                                "NJSX": C.L1,
                            }
                            add.update(add)
                            NJSX = int(C.L1)
                            for k in range(NJSX):
                                L, i = sandy.read_list(df, i)
                                add = {
                                    "NCH": L.L1,
                                    "NRB": L.L2,
                                    "NX": L.N2,
                                }
                                NCH = int(L.L1)
                                LIST2 = {}
                                Ep = L.B[::NCH + 1]
                                b = L.B[::1]
                                del b[::NCH + 1]
                                add_3 = {
                                    "ER": Ep,  # Resonance energy in eV
                                    "GAM": b,  # Channel width in eV
                                }
                                LIST2.update({k: add_3})
                            add.update({"J": LIST2})
                            L, i = sandy.read_list(df, i)
                            add = {
                                "N": L.NPL,
                                "NPARB": L.N2,
                            }
                            add.update(add)
                            NPARB = int(L.N2)
                            tri = np.zeros((NPARB * NPARB, NPARB * NPARB))
                            tri[np.triu_indices(
                                NPARB * NPARB, 0)] = np.array(L.B[6 * NRB:])
                            COVAR_PAR = {
                                "COVAR_PAR": tri
                            }
                            add.update({"COVAR_PAR": COVAR_PAR})
                            LCOMP1.update(add)
                            dico[(EL, EH)] = LCOMP1
                elif LCOMP == 2:
                    LCOMP2 = {}
                    LCOMP2.update(header2)
                    LCOMP2.update(header3)
                    if LRF == 1 or LRF == 2:
                        if ISR > 0:
                            C, i = sandy.read_cont(df, i)
                            add = {
                                "DAP": C.C2,
                            }
                            LCOMP2.update(add)
                            L, i = sandy.read_list(df, i)
                            add = {
                                "AWRI": C.C1,
                                "QX": C.C2,
                                "LRX": L.L2,
                            }
                            NRSA = L.N2
                            keys = ["ER", "AJ", "GT", "GN", "GG", "GF",
                                    "DER", "01", "02", "DGN", "DGG", "DGF"]
                            RES_PAR = [dict(zip(keys, items))
                                       for items in sandy.utils.grouper(L.B, 12)]
                            RES_PAR1 = RES_PAR
                            for h in range(NRSA):
                                del RES_PAR1[h]["01"]
                                del RES_PAR1[h]["02"]
                            RES_PAR = RES_PAR1
                            add.update({"RES_PAR": RES_PAR})
                            LCOMP2.update(add)
                            C, i = sandy.read_cont(df, i)
                            add = {
                                "NDIGIT": C.L1,
                                "NNN": C.L2,
                                "NM": C.N1,
                            }
                            NDIGIT = C.L1
                            NM = C.N1
                            add.update(add)
                            LCOMP2.update(add)
                            INTG = {}
                            for j in range(NM):
                                I, i = read_intg(tape, mat, NDIGIT, i)
                                add_intg = {
                                    "II": int(I.II),
                                    "JJ": int(I.JJ),
                                    "KIJ": I.KIJ,
                                }
                                INTG.update({j: add_intg})
                            add.update({"INTG": INTG})
                            LCOMP2.update(add)
                            dico[(EL, EH)] = LCOMP2
                        elif ISR == 0:
                            L, i = sandy.read_list(df, i)
                            add = {
                                "AWRI": C.C1,
                                "QX": C.C2,
                                "LRX": L.L2,
                            }
                            keys = ["ER", "AJ", "GT", "GN", "GG", "GF",
                                    "DER", "01", "02", "DGN", "DGG", "DGF"]
                            RES_PAR = [dict(zip(keys, items))
                                       for items in sandy.utils.grouper(L.B, 12)]
                            RES_PAR1 = RES_PAR
                            for h in range(NRSA):
                                del RES_PAR1[h]["01"]
                                del RES_PAR1[h]["02"]
                            RES_PAR = RES_PAR1
                            add.update({"RES_PAR": RES_PAR})
                            C, i = sandy.read_cont(df, i)
                            add = {
                                "NDIGIT": C.L1,
                                "NNN": C.L2,
                                "NM": C.N1,
                            }
                            NDIGIT = C.L1
                            NM = C.N1
                            add.update(add)
                            LCOMP2.update(add)
                            INTG = {}
                            for j in range(NM):
                                I, i = read_intg(tape, mat, NDIGIT, i)
                                add_intg = {
                                    "II": int(I.II),
                                    "JJ": int(I.JJ),
                                    "KIJ": I.KIJ,
                                }
                                INTG.update({j: add_intg})
                            add.update({"INTG": INTG})
                            LCOMP2.update(add)
                            dico[(EL, EH)] = LCOMP2
                    elif LRF == 3:
                        if ISR == 1:
                            L, i = sandy.read_list(df, i)
                            add = {
                                "MLS": L.NPL,
                                "DAP": L.B
                            }
                            L, i = sandy.read_list(df, i)
                            add = {
                                "AWRI": C.C1,
                                "APL": C.C2,
                                "NRSA": L.N2,
                            }
                            keys = [
                                "ER",
                                "AJ",
                                "GN",
                                "GG",
                                "GFA",
                                "GFB",
                                "DER",
                                "0",
                                "DGN",
                                "DGG",
                                "DGFA",
                                "DGFB"]
                            RES_PAR = [dict(zip(keys, items))
                                       for items in sandy.utils.grouper(L.B, 12)]
                            RES_PAR1 = RES_PAR
                            for h in range(NRSA):
                                del RES_PAR1[h]["0"]
                            RES_PAR = RES_PAR1
                            add.update({"RES_PAR": RES_PAR})
                            C, i = sandy.read_cont(df, i)
                            add = {
                                "NDIGIT": C.L1,
                                "NNN": C.L2,
                                "NM": C.N1,
                            }
                            NDIGIT = C.L1
                            NM = C.N1
                            add.update(add)
                            LCOMP2.update(add)
                            INTG = {}
                            for j in range(NM):
                                I, i = read_intg(tape, mat, NDIGIT, i)
                                add_intg = {
                                    "II": int(I.II),
                                    "JJ": int(I.JJ),
                                    "KIJ": I.KIJ,
                                }
                                INTG.update({j: add_intg})
                            add.update({"INTG": INTG})
                            LCOMP2.update(add)
                            dico[(EL, EH)] = LCOMP2
                        else:
                            L, i = sandy.read_list(df, i)
                            add = {
                                "AWRI": C.C1,
                                "APL": C.C2,
                                "NRSA": L.N2,
                            }
                            NRSA = L.N2
                            keys = [
                                "ER",
                                "AJ",
                                "GN",
                                "GG",
                                "GFA",
                                "GFB",
                                "DER",
                                "0",
                                "DGN",
                                "DGG",
                                "DGFA",
                                "DGFB"]
                            RES_PAR = [dict(zip(keys, items))
                                       for items in sandy.utils.grouper(L.B, 12)]
                            RES_PAR1 = RES_PAR
                            for h in range(NRSA):
                                del RES_PAR1[h]["0"]
                            RES_PAR = RES_PAR1
                            add.update({"RES_PAR": RES_PAR})
                            C, i = sandy.read_cont(df, i)
                            add = {
                                "NDIGIT": C.L1,
                                "NNN": C.L2,
                                "NM": C.N1,
                            }
                            NDIGIT = C.L1
                            NM = C.N1
                            add.update(add)
                            LCOMP2.update(add)
                            INTG = {}
                            for j in range(NM):
                                I, i = read_intg(tape, mat, NDIGIT, i)
                                add_intg = {
                                    "II": int(I.II),
                                    "JJ": int(I.JJ),
                                    "KIJ": I.KIJ,
                                }
                                INTG.update({j: add_intg})
                            add.update({"INTG": INTG})
                            LCOMP2.update(add)
                            dico[(EL, EH)] = LCOMP2
                    elif LRF == 7:
                        if ISR > 0:
                            L, i = sandy.read_list(df, i)
                            add = {

                                "JCH": L.NPL,
                                "(1+(NJCH-1)/6)": L.N2,
                                "DAP": L.B,
                            }
                            LCOMP2.update(add)
                            L, i = sandy.read_list(df, i)
                            add = {
                                "NPP": L.L1,
                                "NJSX": L.L2,
                            }
                            keys = ["MA", "MB", "ZA", "ZB", "IA", "IB",
                                    "Q", "PNT", "SHF", "MT", "PA", "PB"]
                            RES_PAR = [dict(zip(keys, items))
                                       for items in sandy.utils.grouper(L.B, 12)]
                            add.update({"RES_PAR": RES_PAR})
                            LCOMP2.update(add)
                            for k in range(NLS):
                                LIST = {}
                                L, i = sandy.read_list(df, i)
                                add = {
                                    "AJ": L.C1,
                                    "PJ": L.C2,
                                    "NCH": L.N2,
                                }
                                NCH = int(L.N1)
                                keys = ["PPI", "L", "SCH", "BND", "APE", "APT"]
                                RES_PAR = [dict(zip(keys, items))
                                           for items in sandy.utils.grouper(L.B, 6)]
                                add.update({"RES_PAR": RES_PAR})
                                LCOMP2.update(add)
                                L, i = sandy.read_list(df, i)
                                add = {
                                    "NCH": L.L1,
                                    "NRB": L.L2,
                                    "NX": L.N2,
                                }

                                Ep = L.B[::NCH + 1]
                                b = L.B[::1]
                                del b[::NCH + 1]
                                add_2 = {
                                    "ER": Ep,  # Resonance energy in eV
                                    "GAM": b,  # Channel width in eV
                                }
                                add.update(add_2)
                                LIST.update({k: add})
                            LCOMP2.update({"J": LIST})
                            C, i = sandy.read_cont(df, i)
                            add = {
                                "NDIGIT": C.L1,
                                "NNN": C.L2,
                                "NM": C.N1,
                            }
                            NDIGIT = C.L1
                            NM = C.N1
                            LCOMP2.update(add)
                            INTG = {}
                            for j in range(NM):
                                I, i = read_intg(tape, mat, NDIGIT, i)
                                add_intg = {
                                    "II": int(I.II),
                                    "JJ": int(I.JJ),
                                    "KIJ": I.KIJ,
                                }
                                INTG.update({j: add_intg})
                            add.update({"INTG": INTG})
                            LCOMP2.update(add)
                            dico[(EL, EH)] = LCOMP2
                        else:
                            L, i = sandy.read_list(df, i)
                            add = {
                                "NPP": L.L1,
                                "NJSX": L.L2,
                            }
                            keys = ["MA", "MB", "ZA", "ZB", "IA", "IB",
                                    "Q", "PNT", "SHF", "MT", "PA", "PB"]
                            RES_PAR = [dict(zip(keys, items))
                                       for items in sandy.utils.grouper(L.B, 12)]
                            add.update({"RES_PAR": RES_PAR})
                            LCOMP2.update(add)
                            for k in range(NLS):
                                LIST = {}
                                L, i = sandy.read_list(df, i)
                                add = {
                                    "AJ": L.C1,
                                    "PJ": L.C2,
                                    "NCH": L.N2,
                                }
                                NCH = int(L.N1)
                                keys = ["PPI", "L", "SCH", "BND", "APE", "APT"]
                                RES_PAR = [dict(zip(keys, items))
                                           for items in sandy.utils.grouper(L.B, 6)]
                                add.update({"RES_PAR": RES_PAR})
                                LCOMP2.update(add)
                                L, i = sandy.read_list(df, i)
                                add = {
                                    "NCH": L.L1,
                                    "NRB": L.L2,
                                    "NX": L.N2,
                                }
                                Ep = L.B[::NCH + 1]
                                b = L.B[::1]
                                del b[::NCH + 1]
                                add_2 = {
                                    "ER": Ep,  # Resonance energy in eV
                                    "GAM": b,  # Channel width in eV
                                }
                                add.update(add_2)
                                LIST.update({k: add})
                            LCOMP2.update({"J": LIST})
                            C, i = sandy.read_cont(df, i)
                            add = {
                                "NDIGIT": C.L1,
                                "NNN": C.L2,
                                "NM": C.N1,
                            }
                            NDIGIT = C.L1
                            NM = C.N1
                            LCOMP2.update(add)
                            INTG = {}
                            for j in range(NM):
                                I, i = read_intg(tape, mat, NDIGIT, i)
                                add_intg = {
                                    "II": int(I.II),
                                    "JJ": int(I.JJ),
                                    "KIJ": I.KIJ,
                                }
                                INTG.update({j: add_intg})
                            add.update({"INTG": INTG})
                            LCOMP2.update(add)
                            dico[(EL, EH)] = LCOMP2
        NER1.update(dico)
        info.update({"NER": NER1})
        M.update(info)
        ISO.update({ZAI: M})
    P.update(ISO)
    out.update({"NIS": P})

    return out
