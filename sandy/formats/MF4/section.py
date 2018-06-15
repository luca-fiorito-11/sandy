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
def read_mf4_mt(text):
    str_list = split_endf(text)
    i = 0
    out = {"MAT" : str_list["MAT"].iloc[0],
           "MF" : str_list["MF"].iloc[0],
           "MT" : str_list["MT"].iloc[0]}
    C, i = read_cont(str_list, i)
    out.update({"ZA" : C.C1, "AWR" : C.C2, "LTT" : C.L2})
    if out["LTT"] in (1,3):
        C, i = read_cont(str_list, i)
        lpc = {"LI" : C.L1, "LCT" : C.L2}
        T2, i = read_tab2(str_list, i)
        lpc.update({"NE" : T2.NZ, "NBT" : T2.NBT, "INT" : T2.INT, "E" : {} })
        for j in range(lpc["NE"]):
            L, i = read_list(str_list, i)
            lpc["E"].update({ L.C2 : {"P" : L.B, "T" : L.C1, "LT" : L.L1}})
        out.update({"LPC" : lpc})
    if out["LTT"] in (2,3):
        C, i = read_cont(str_list, i)
        sub = {"LI" : C.L1, "LCT" : C.L2}
        T2, i = read_tab2(str_list, i)
        sub.update({"NE" : T2.NZ, "NBT" : T2.NBT, "INT" : T2.INT, "E" : {}})
        for i in range(sub["NE"]):
            T1, i = read_tab1(str_list, i)
            distr = pd.Series(T1.y, index = T1.x, name=T1.C2).rename_axis("mu")
            sub["E"].update({ T1.C2 : {"T" : T1.C1, "LT" : T1.L1, "PDF" : distr, "NBT" : T1.NBT, "INT" : T1.INT}})
        out.update({"TAB" : sub})
    if out["LTT"] == 0:
        C, i = read_cont(str_list, i)
        out.update({"ISO" : {"LI" : C.L1, "LCT" : C.L2}})
    return out