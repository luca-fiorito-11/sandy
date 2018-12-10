# -*- coding: utf-8 -*-
"""
Created on Thu Jun 14 09:23:33 2018

@author: fiorito_l
"""
import pdb

from .records import *
from .utils import Section

__author__ = "Luca Fiorito"
__all__ = ["read", "write", "read_errorr", "read_groupr"]

def read(text):
    str_list = text.splitlines()
    MAT, MF, MT = read_control(str_list[0])[:3]
    out = {"MAT" : MAT, "MF" : MF, "MT" : MT}
    i = 0
    C, i = read_cont(str_list, i)
    out.update({"ZA" : C.C1, "AWR" : C.C2, "PFLAG" : C.L2})
    T, i = read_tab1(str_list, i)
    out.update({"QM" : T.C1, "QI" : T.C2, "LR" : T.L2, "NBT" : T.NBT, "INT" : T.INT, "E" : T.x, "XS" : T.y})
    return Section(out)

def write(sec):
    text = write_cont(sec["ZA"], sec["AWR"], 0, sec["PFLAG"], 0, 0)
    text += write_tab1(sec["QM"], sec["QI"], 0, sec["LR"], sec["NBT"], sec["INT"], sec["E"], sec["XS"])
    TextOut = []; iline = 1
    for line in text:
        if iline > 99999:
            iline = 1
        TextOut.append("{:<66}{:4}{:2}{:3}{:5}\n".format(line, sec["MAT"], sec["MF"], sec["MT"], iline))
        iline += 1
    return "".join(TextOut)

def read_errorr(text):
    str_list = text.splitlines()
    MAT, MF, MT = read_control(str_list[0])[:3]
    out = {"MAT" : MAT, "MF" : MF, "MT" : MT}
    i = 0
    L, i = read_list(str_list, i)
    out.update({"XS" : L.B})
    return out

def read_groupr(text):
    str_list = text.splitlines()
    MAT, MF, MT = read_control(str_list[0])[:3]
    out = {"MAT" : MAT, "MF" : MF, "MT" : MT}
    i = 0
    C, i = read_cont(str_list, i)
    out.update({"ZA" : C.C1, "AWR" : C.C2, "NL" : C.L1, "NZ" : C.L2, "LRFLAG" : C.N1, "NGN" : C.N2})
    groups = {}
    for ig in range(out["NGN"]):
        L, i = read_list(str_list, i)
        group = {"TEMPIN" : L.C1, "NG2" : L.L1, "IG2LO" : L.L2, "IG" : L.N2, "DATA" : L.B}
        groups[L.N2] = group
    out["GROUPS"] = groups
    return out