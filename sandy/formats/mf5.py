# -*- coding: utf-8 -*-
"""
This module contains only two public functions:

    * `read`
    * `write`

Function `read` reads a MF5/MT section from a string and produces a content object with a dictionary-like 
structure.
The content object can be accessed using most of the keywords specified in the ENDF6 manual for this specific 
MF section.

Function `write` writes a content object for a MF5/MT section into a string.
MAT, MF, MT and line numbers are also added (each line ends with a `\n`).
"""

import pdb

from sandy.formats.records import *
from sandy.formats.utils import Section

__author = "Luca Fiorito"
__all__ = ["read", "write"]

def read(text):
    """Read MT section for MF5
    
    Parameters
    ----------
    text : `str`
        one string containing the whole section
    
    Returns
    -------
    `sandy.formats.utils.Section`
        MF5 content sructured as a collection of dictionaries 
    """
    str_list = text.splitlines()
    MAT, MF, MT = read_control(str_list[0])[:3]
    out = {"MAT" : MAT, "MF" : MF, "MT" : MT}
    i = 0
    C, i = read_cont(str_list, i)
    out["ZA"] = C.C1    # ZA = Z*1000 + A
    out["AWR"] = C.C2   # Ratio of the nuclide mass to that of neutron
    nk = C.N1           # Number of partial energy distributions. There will be one subsection for each partial distribution.
    pdistr = {}
    for j in range(nk):
        Tp, i = read_tab1(str_list, i)
        sub = {"NBT_P" : Tp.NBT, "INT_P" : Tp.INT, "E_P" : Tp.x, "P" : Tp.y} # Fractional part of the particular cross section which can be described by the kth partial energy distribution at the i-th incident energy point
        sub["LF"] = Tp.L2  # Flag specifying the energy distribution law used for a particular subsection (partial energy distribution)
        if sub["LF"] == 5: # General Evaporation Spectrum (LF=5)
            """
            Found in:
                100-Fm-255g.jeff33 (x6)
                88-Ra-226g.jeff33 (x6)
                91-Pa-233g.jeff33 (x6)
                92-U-239g.jeff33
                92-U-240g.jeff33
            """
            sub.update({'U' : Tp.C1})
            T, i = read_tab1(str_list, i)
            sub.update({"NBT_THETA" : T.NBT, "INT_THETA" : T.INT, "E_THETA" : T.x, "THETA" : T.y})
            T, i = read_tab1(str_list, i)
            sub.update({"NBT_G" : T.NBT, "INT_G" : T.INT, "E_G" : T.x, "G" : T.y})
        elif sub["LF"] in (7,9): # Simple Maxwellian Fission Spectrum (LF=7) / Evaporation Spectrum (LF=9)
            """
            Found in:
                27-Co-59g.jeff33
            """
            sub.update({'U' : Tp.C1})
            T, i = read_tab1(str_list, i)
            sub.update({"NBT_THETA" : T.NBT, "INT_THETA" : T.INT, "E_THETA" : T.x, "THETA" : T.y})
        elif sub["LF"] == 11: # Energy-Dependent Watt Spectrum (LF=11)
            sub.update({'U' : Tp.C1})
            T, i = read_tab1(str_list, i)
            sub.update({"NBT_A" : T.NBT, "INT_A" : T.INT, "E_A" : T.x, "A" : T.y})
            T, i = read_tab1(str_list, i)
            sub.update({"NBT_B" : T.NBT, "INT_B" : T.INT, "E_B" : T.x, "B" : T.y})
        elif sub["LF"] == 12: # Energy-Dependent Fission Neutron Spectrum (Madland and Nix) (LF=12)
            TM, i = read_tab1(str_list, i)
            sub.update({"EFL" : T.C1, "EHL" : T.C2, "NBT_TM" : T.NBT, "INT_TM" : T.INT, "E_TM" : T.x, "TM" : T.y})
        elif sub["LF"] == 1: # Arbitrary Tabulated Function (LF=1)
            T2, i = read_tab2(str_list, i)
            sub.update({ "NBT_EIN" : T2.NBT, "INT_EIN" : T2.INT, "EIN" : {} })
            for k in range(T2.NZ):
                T1, i = read_tab1(str_list, i)
                sub["EIN"].update({ T1.C2 : {"EOUT" : T1.x, "EDISTR" : T1.y, "NBT" : T1.NBT, "INT" : T1.INT}})
        pdistr[j] = sub
    if pdistr:
        out["PDISTR"] = pdistr
    return Section(out)



def write(sec):
    """Write MT section for MF5
    
    Parameters
    ----------
    sec : `sandy.utils.Section`
        dictionary with MT section for MF5
    
    Returns
    -------
    `str`
        section content in a single string
    """
    text = write_cont(sec["ZA"], sec["AWR"], 0, 0, len(sec["PDISTR"]), 0)
    for k, sub in sorted(sec["PDISTR"].items()):
        U = sub['U'] if 'U' in sub else 0
        text += write_tab1(U, 0, 0, sub["LF"], sub["NBT_P"], sub["INT_P"], sub["E_P"], sub["P"])
        if sub["LF"] == 1:
            text += write_tab2(0, 0, 0, 0, len(sub['EIN']), sub["NBT_EIN"], sub["INT_EIN"])
            for ein, distr in sorted(sub['EIN'].items()):
                text += write_tab1(0, ein, 0, 0, distr["NBT"], distr["INT"], distr["EOUT"], distr["EDISTR"])
        elif sub["LF"] == 5:
            text += write_tab1(0, 0, 0, 0, sub["NBT_THETA"], sub["INT_THETA"], sub["E_THETA"], sub["THETA"])
            text += write_tab1(0, 0, 0, 0, sub["NBT_G"], sub["INT_G"], sub["E_G"], sub["G"])
        elif sub["LF"] in (7,9):
            text += write_tab1(0, 0, 0, 0, sub["NBT_THETA"], sub["INT_THETA"], sub["E_THETA"], sub["THETA"])
        elif sub["LF"] == 11:
            text += write_tab1(0, 0, 0, 0, sub["NBT_A"], sub["INT_A"], sub["E_A"], sub["A"])
            text += write_tab1(0, 0, 0, 0, sub["NBT_B"], sub["INT_B"], sub["E_B"], sub["B"])
        elif sub["LF"] == 12:
            text += write_tab1(0, 0, 0, 0, sub["NBT_TM"], sub["INT_TM"], sub["E_TM"], sub["TM"])
    textout = []
    iline = 1
    for line in text:
        if iline > 99999:
            iline = 1
        textout.append("{:<66}{:4}{:2}{:3}{:5}\n".format(line, sec["MAT"], sec["MF"], sec["MT"], iline))
        iline += 1
    return "".join(textout)