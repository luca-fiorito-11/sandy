# -*- coding: utf-8 -*-
"""
Created on Thu Jun 14 09:23:33 2018

@author: fiorito_l
"""
from .records import *
from .utils import Section

def read(text):
    str_list = text.splitlines()
    MAT, MF, MT = read_control(str_list[0])[:3]
    out = {"MAT" : MAT, "MF" : MF, "MT" : MT}
    i = 0
    C, i = read_cont(str_list, i)
    out.update({ "ZA" : C.C1, "AWR" : C.C2, "NK" : C.N1, "SUB" : {}})
    for k in range(out["NK"]):
        L, i = read_list(str_list, i)
        D = {"ELO" : L.C1, "EHI" : L.C2, "LS" : L.L1, "LB" : L.L2, "NE" : L.N2, "EK" : L.B[:L.N2], "FKK" : L.B[L.N2:]}
        out["SUB"].update({k+1 : D})
    return Section(out)