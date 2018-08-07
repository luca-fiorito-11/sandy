# -*- coding: utf-8 -*-
"""
Created on Thu Jun 14 09:23:33 2018

@author: fiorito_l
"""
from ..records import read_cont, read_list, read_control
from ..utils import Section

__author__ = "Luca Fiorito"
__all__ = ["read"]

def read(text):
    str_list = text.splitlines()
    MAT, MF, MT = read_control(str_list[0])[:3]
    if MT in (454, 459):
        return read_fy(text)
    elif MT == 457:
        return read_rdd(text)

def read_fy(text):
    str_list = text.splitlines()
    MAT, MF, MT = read_control(str_list[0])[:3]
    out = {"MAT" : MAT, "MF" : MF, "MT" : MT}
    i = 0
    C, i = read_cont(str_list, i)
    out.update({"ZA" : C.C1, "AWR" : C.C2, "E" : {}})
    for j in range(C.L1):
        L, i = read_list(str_list, i)
        FY = { int(ZAFP*10+FPS) : {"ZAFP" : ZAFP, "FPS" : FPS, "YI" : YI, "DYI" : DYI} for ZAFP,FPS,YI,DYI in  zip(*[iter(L.B)]*4)}
        out["E"].update({ L.C1 : { "FY" : FY } })
        if j > 0:
            out["E"][L.C1].update({ "I" : L.L1 })
    return Section(out)

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