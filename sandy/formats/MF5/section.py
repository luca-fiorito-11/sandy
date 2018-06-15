# -*- coding: utf-8 -*-
"""
Created on Thu Jun 14 09:23:33 2018

@author: fiorito_l
"""
from ..records2 import read_cont, read_tab1, read_control, read_tab2
from ..utils import Section

def read(text):
    str_list = text.splitlines()
    MAT, MF, MT = read_control(str_list[0])[:3]
    out = {"MAT" : MAT, "MF" : MF, "MT" : MT}
    i = 0
    C, i = read_cont(str_list, i)
    # subsections for partial energy distributions are given in a list
    out.update({"ZA" : C.C1, "AWR" : C.C2, "NK" : C.N1, "PDISTR" : {} })
    for j in range(out["NK"]):
        Tp, i = read_tab1(str_list, i)
        sub = {"LF" : Tp.L2, "NBT_P" : Tp.NBT, "INT_P" : Tp.INT, "E_P" : Tp.x, "P" : Tp.y}
        if sub["LF"] == 5:
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
        elif sub["LF"] in (7,9):
            """
            Found in:
                27-Co-59g.jeff33
            """
            sub.update({'U' : Tp.C1})
            T, i = read_tab1(str_list, i)
            sub.update({"NBT_THETA" : T.NBT, "INT_THETA" : T.INT, "E_THETA" : T.x, "THETA" : T.y})
        elif sub["LF"] == 11:
            sub.update({'U' : Tp.C1})
            T, i = read_tab1(str_list, i)
            sub.update({"NBT_A" : T.NBT, "INT_A" : T.INT, "E_A" : T.x, "A" : T.y})
            T, i = read_tab1(str_list, i)
            sub.update({"NBT_B" : T.NBT, "INT_B" : T.INT, "E_B" : T.x, "B" : T.y})
        elif sub["LF"] == 12:
            TM, i = read_tab1(str_list, i)
            sub.update({"EFL" : T.C1, "EHL" : T.C2, "NBT_TM" : T.NBT, "INT_TM" : T.INT, "E_TM" : T.x, "TM" : T.y})
        elif sub["LF"] == 1:
            T2, i = read_tab2(str_list, i)
            sub.update({ "NBT_EIN" : T2.NBT, "INT_EIN" : T2.INT, "EIN" : {} })
            for k in range(T2.NZ):
                T1, i = read_tab1(str_list, i)
                sub["EIN"].update({ T1.C2 : {"EOUT" : T1.x, "PDF" : T1.y, "NBT" : T1.NBT, "INT" : T1.INT}})
        out["SUB"].update({j : sub})
    return Section(out)