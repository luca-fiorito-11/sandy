# -*- coding: utf-8 -*-
"""
Created on Thu Jun 14 09:23:33 2018

@author: fiorito_l
"""
from ..records import read_cont, read_tab1, read_tab2, read_list, read_control
from ..utils import Section
import pdb

def read(text):
    str_list = text.splitlines()
    MAT, MF, MT = read_control(str_list[0])[:3]
    out = {"MAT" : MAT, "MF" : MF, "MT" : MT}
    i = 0
    pdb.set_trace()
    C, i = read_cont(str_list, i)
    out.update({"ZA" : C.C1, "AWR" : C.C2, "LTT" : C.L2})
    C, i = read_cont(str_list, i)
    out.update({"LI" : C.L1, "LCT" : C.L2})
    if out["LTT"] in (1,3): # LEGENDRE
        T2, i = read_tab2(str_list, i)
        sub = {"NE" : T2.NZ, "NBT" : T2.NBT, "INT" : T2.INT, "E" : {} }
        for j in range(sub["NE"]):
            L, i = read_list(str_list, i)
            sub["E"].update({ L.C2 : {"COEFF" : L.B, "T" : L.C1, "LT" : L.L1}})
        out.update({"LPC" : sub})
    if out["LTT"] in (2,3): # TABULATED
        T2, i = read_tab2(str_list, i)
        sub = {"NE" : T2.NZ, "NBT" : T2.NBT, "INT" : T2.INT, "E" : {}}
        for i in range(sub["NE"]):
            T, i = read_tab1(str_list, i)
            distr = {"T" : T.C1, "LT" : T.L1, "NBT" : T.NBT, "INT" : T.INT, "E" : T.x, "ADISTR" : T.y}
            sub["E"].update({T.C2 : distr})
        out.update({"TAB" : sub})
    return Section(out)
