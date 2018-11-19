# -*- coding: utf-8 -*-
"""
Created on Thu Jun 14 09:23:33 2018

@author: fiorito_l
"""

import numpy as np

from ..records import read_cont, read_list, read_control, write_cont, write_list
from ..utils import Section

__author__ = "Luca Fiorito"
__all__ = ["read", "write"]

def read(text):
    """Read MT section for MF8
    
    Parameters
    ----------
    text: `str`
        one string containing the whole section
    
    Returns
    -------
    `sandy.utils.Section`
    """
    str_list = text.splitlines()
    MAT, MF, MT = read_control(str_list[0])[:3]
    if MT in (454, 459):
        return _read_fy(text)
    elif MT == 457:
        return _read_rdd(text)

def write(sec):
    """Write MT section for MF8
    
    Parameters
    ----------
    sec: `sandy.utils.Section`
        dictionary with MT section for MF8
    
    Returns
    -------
    `str`
    """
    if sec["MT"] in (454, 459):
        return _write_fy(sec)
    elif sec["MT"] == 457:
        return _write_rdd(sec)

def _read_fy(text):
    str_list = text.splitlines()
    MAT, MF, MT = read_control(str_list[0])[:3]
    out = {"MAT" : MAT, "MF" : MF, "MT" : MT}
    i = 0
    C, i = read_cont(str_list, i)
    out.update({"ZA" : C.C1, "AWR" : C.C2, "E" : {}})
    for j in range(C.L1):
        L, i = read_list(str_list, i)
        FY = { int(ZAFP*10+FPS) : {"ZAFP" : int(ZAFP), "FPS" : int(FPS), "YI" : YI, "DYI" : DYI} for ZAFP,FPS,YI,DYI in  zip(*[iter(L.B)]*4)}
        out["E"].update({ L.C1 : { "FY" : FY } })
        if j > 0:
            out["E"][L.C1].update({ "I" : L.L1 })
    return Section(out)

def _write_fy(sec):
    LE = len(sec["E"])
    text = write_cont(sec["ZA"], sec["AWR"], LE, 0, 0, 0)
    for i,(e,esec) in enumerate(sorted(sec["E"].items())):
        tab = [ksec[j] for k,ksec in sorted(esec["FY"].items()) for j in ("ZAFP","FPS","YI","DYI")]
        NFP = len(esec["FY"])
        I = LE-1 if i == 0 else esec["I"]
        text += write_list(e, 0, I, 0, NFP, tab)
    TextOut = []; iline = 1
    for line in text:
        if iline > 99999:
            iline = 1
        TextOut.append("{:<66}{:4}{:2}{:3}{:5}\n".format(line, sec["MAT"], sec["MF"], sec["MT"], iline))
        iline += 1
    return "".join(TextOut)

#def read_rdd(text):
#    str_list = split_endf(text)
#    i = 0
#    out = {"MAT" : str_list["MAT"].iloc[0],
#           "MF" : str_list["MF"].iloc[0],
#           "MT" : str_list["MT"].iloc[0]}
#    C, i = read_cont(str_list, i)
#    out.update({"ZA" : C.C1, "AWR" : C.C2, "LIS" : C.L1, "LISO" : C.L2, "NST" :C.N1, "NSP" : C.N2})
#    L, i = read_list(str_list, i)
#    out.update({"HL" : L.C1, "DHL" : L.C2, "E" : L.B[::2], "DE" : L.B[1::2]})
#    L, i = read_list(str_list, i)
#    out.update({"SPI" : L.C1, "PAR" : L.C2, "DK" : []})
#    if out["NST"] == 0:
#        # Update list of decay modes when nuclide is radioactive
#        out.update({ "DK" : [ {"RTYP" : RTYP, "RFS" : RFS, "Q" : Q, "DQ" : DQ, "BR" : BR, "DBR" : DBR } for RTYP, RFS, Q, DQ, BR, DBR in  zip(*[iter(L.B)]*6) ] })
#    return out