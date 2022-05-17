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


def read_intg(df, NDIGIT, ipos):
    INTG = namedtuple('INTG', 'II JJ KIJ')
    series = df.iloc[ipos]
    line = "".join(series.values.astype(str))
    if NDIGIT == 2:
        df = pd.read_fwf(
            io.StringIO(line),
            widths=[
                5,
                5,
                1,
                3,
                3,
                3,
                3,
                3,
                3,
                3,
                3,
                3,
                3,
                3,
                3,
                3,
                3,
                3,
                3,
                3,
                3],
            names=[
                "II",
                "JJ",
                "/",
                "KIJ1",
                "KIJ2",
                "KIJ3",
                "KIJ4",
                "KIJ5",
                "KIJ6",
                "KIJ7",
                "KIJ8",
                "KIJ9",
                "KIJ10",
                "KIJ11",
                "KIJ12",
                "KIJ13",
                "KIJ14",
                "KIJ15",
                "KIJ16",
                "KIJ17",
                "KIJ18"],
            dtype={
                "II": int,
                "JJ": int,
                "/": str,
                "KIJ1": str,
                "KIJ2": str,
                "KIJ3": str,
                "KIJ4": str,
                "KIJ5": str,
                "KIJ6": str,
                "KIJ7": str,
                "KIJ8": str,
                "KIJ9": str,
                "KIJ10": str,
                "KIJ11": str,
                "KIJ12": str,
                "KIJ13": str,
                "KIJ14": str,
                "KIJ15": str,
                "KIJ16": str,
                "KIJ17": str,
                "KIJ18": str},
            na_filter=True,
        )
        KIJ1 = []
        for i in range(1, 19):
            KIJ1.append(list(df[f"KIJ{i}"]))
        KIJ = np.array(KIJ1).squeeze()
        II = df.II
        JJ = df.JJ

    elif NDIGIT == 3:
        df = pd.read_fwf(
            io.StringIO(line),
            widths=[
                5,
                5,
                1,
                4,
                4,
                4,
                4,
                4,
                4,
                4,
                4,
                4,
                4,
                4,
                4,
                4,
                4],
            names=[
                "II",
                "JJ",
                "/",
                "KIJ1",
                "KIJ2",
                "KIJ3",
                "KIJ4",
                "KIJ5",
                "KIJ6",
                "KIJ7",
                "KIJ8",
                "KIJ9",
                "KIJ10",
                "KIJ11",
                "KIJ12",
                "KIJ13",
                "KIJ14"],
            dtype={
                "II": int,
                "JJ": int,
                "/": str,
                "KIJ1": str,
                "KIJ2": str,
                "KIJ3": str,
                "KIJ4": str,
                "KIJ5": str,
                "KIJ6": str,
                "KIJ7": str,
                "KIJ8": str,
                "KIJ9": str,
                "KIJ10": str,
                "KIJ11": str,
                "KIJ12": str,
                "KIJ13": str,
                "KIJ14": str},
            na_filter=False,
        )
        KIJ1 = []
        for i in range(1, 15):
            KIJ1.append(list(df[f"KIJ{i}"]))
        KIJ = np.array(KIJ1).squeeze()
        II = df.II
        JJ = df.JJ

    elif NDIGIT == 4:
        df = pd.read_fwf(
            io.StringIO(line),
            widths=[
                5,
                5,
                1,
                5,
                5,
                5,
                5,
                5,
                5,
                5,
                5,
                5,
                5,
                5],
            names=[
                "II",
                "JJ",
                "/",
                "KIJ1",
                "KIJ2",
                "KIJ3",
                "KIJ4",
                "KIJ5",
                "KIJ6",
                "KIJ7",
                "KIJ8",
                "KIJ9",
                "KIJ10",
                "KIJ11"],
            dtype={
                "II": int,
                "JJ": int,
                "/": str,
                "KIJ1": str,
                "KIJ2": str,
                "KIJ3": str,
                "KIJ4": str,
                "KIJ5": str,
                "KIJ6": str,
                "KIJ7": str,
                "KIJ8": str,
                "KIJ9": str,
                "KIJ10": str,
                "KIJ11": str},
            na_filter=False,
        )
        KIJ1 = []
        for i in range(1, 12):
            KIJ1.append(list(df[f"KIJ{i}"]))
        KIJ = np.array(KIJ1).squeeze()
        II = df.II
        JJ = df.JJ

    elif NDIGIT == 5:
        df = pd.read_fwf(
            io.StringIO(line),
            widths=[
                5,
                5,
                1,
                6,
                6,
                6,
                6,
                6,
                6,
                6,
                6,
                6],
            names=[
                "II",
                "JJ",
                "/",
                "KIJ1",
                "KIJ2",
                "KIJ3",
                "KIJ4",
                "KIJ5",
                "KIJ6",
                "KIJ7",
                "KIJ8",
                "KIJ9"],
            dtype={
                "II": int,
                "JJ": int,
                "/": str,
                "KIJ1": str,
                "KIJ2": str,
                "KIJ3": str,
                "KIJ4": str,
                "KIJ5": str,
                "KIJ6": str,
                "KIJ7": str,
                "KIJ8": str,
                "KIJ9": str},
            na_filter=False,
        )
        KIJ1 = []
        for i in range(1, 10):
            KIJ1.append(list(df[f"KIJ{i}"]))
        KIJ = np.array(KIJ1).squeeze()
        II = df.II
        JJ = df.JJ
    elif NDIGIT == 6:
        df = pd.read_fwf(
            io.StringIO(line),
            widths=[
                5,
                5,
                7,
                7,
                7,
                7,
                7,
                7,
                7,
                7],
            names=[
                "II",
                "JJ",
                "KIJ1",
                "KIJ2",
                "KIJ3",
                "KIJ4",
                "KIJ5",
                "KIJ6",
                "KIJ7"],
            dtype={
                "II": int,
                "JJ": int,
                "KIJ1": str,
                "KIJ2": str,
                "KIJ3": str,
                "KIJ4": str,
                "KIJ5": str,
                "KIJ6": str,
                "KIJ7": str},
            na_filter=False,
        )
        KIJ1 = []
        for i in range(1, 9):
            KIJ1.append(list(df[f"KIJ{i}"]))
        KIJ = np.array(KIJ1).squeeze()
        II = df.II
        JJ = df.JJ
    ipos += 1
    return INTG(II, JJ, KIJ), ipos


def read_mf32(tape, mat):
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
                                "DF²",
                                "DJDN",
                                "DJDG",
                                "DJDF",
                                "DJ²"]
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
                    if LRF == 1 or LRF == 2:
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
                    elif LRF == 3:
                        if ISR > 0:
                            L, i = sandy.read_list(df, i)
                            add = {
                                "MLS": L.NPL
                            }
                            keys = ["DAP"]
                            DAP = [dict(zip(keys, items))
                                   for items in sandy.utils.grouper(L.B, 1)]
                            add.update({"DAP": DAP})
                            C, i = sandy.read_cont(df, i)
                            add = {
                                "AWRI": C.C1,
                                "NSRS": C.N1,
                                "NLRS": C.N2,
                            }
                            add.update(add)
                            LCOMP1.update(add)
                            dico[(EL, EH)] = LCOMP1
                        if ISR == 0:
                            C, i = sandy.read_cont(df, i)
                            add = {
                                "AWRI": C.C1,
                                "NSRS": C.N1,
                                "NLRS": C.N2,
                            }
                            add.update(add)
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
                                I, i = read_intg(df, NDIGIT, i)
                                add_intg = {
                                    "II": I.II,
                                    "JJ": I.JJ,
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
                                I, i = read_intg(df, NDIGIT, i)
                                add_intg = {
                                    "II": I.II,
                                    "JJ": I.JJ,
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
                                I, i = read_intg(df, NDIGIT, i)
                                add_intg = {
                                    "II": I.II,
                                    "JJ": I.JJ,
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
                                I, i = read_intg(df, NDIGIT, i)
                                add_intg = {
                                    "II": I.II,
                                    "JJ": I.JJ,
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
                                I, i = read_intg(df, NDIGIT, i)
                                add_intg = {
                                    "II": I.II,
                                    "JJ": I.JJ,
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
                                I, i = read_intg(df, NDIGIT, i)
                                add_intg = {
                                    "II": I.II,
                                    "JJ": I.JJ,
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
