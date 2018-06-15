# -*- coding: utf-8 -*-
"""
Created on Thu Jun 14 09:23:33 2018

@author: fiorito_l
"""
from ..records2 import read_cont, read_tab1, read_control
from ..utils import Section

def read(text):
    str_list = text.splitlines()
    MAT, MF, MT = read_control(str_list[0])[:3]
    out = {"MAT" : MAT, "MF" : MF, "MT" : MT}
    i = 0
    C, i = read_cont(str_list, i)
    out.update({"ZA" : C.C1, "AWR" : C.C2})
    T, i = read_tab1(str_list, i)
    out.update({"QM" : T.C1, "QI" : T.C2, "LR" : T.L2, "NBT" : T.NBT, "INT" : T.INT, "E" : T.x, "XS" : T.y})
    return Section(out)

def read_mf35_mt(str_list):
#    str_list = text.splitlines()
    i = 0
    out = {"MAT" : str_list["MAT"].iloc[0],
           "MF" : str_list["MF"].iloc[0],
           "MT" : str_list["MT"].iloc[0]}
#    out = {"MAT" : int(str_list[i][66:70]),
#           "MF" : int(str_list[i][70:72]),
#           "MT" : int(str_list[i][72:75])}
    C, i = read_cont(str_list, i)
    out.update({ "ZA" : C.C1, "AWR" : C.C2, "NK" : C.N1, "SUB" : []})
    for k in range(out["NK"]):
        L, i = read_list(str_list, i)
        out["SUB"].append({ "Elo" : L.C1, "Ehi" : L.C2, "NE" : L.N2, "Ek" : L.B[:L.N2],
           "Fkk" : L.B[L.N2:] })
    return out