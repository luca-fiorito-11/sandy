# -*- coding: utf-8 -*-
"""
Created on Thu Jun 14 09:23:33 2018

@author: fiorito_l
"""
from .records import *
from .utils import Section

__author = "Luca Fiorito"
__all__ = ["read", "write"]

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
                sub["EIN"].update({ T1.C2 : {"EOUT" : T1.x, "EDISTR" : T1.y, "NBT" : T1.NBT, "INT" : T1.INT}})
        out["PDISTR"].update({j+1 : sub})
    return Section(out)

def write(sec):
    text = write_cont(sec["ZA"], sec["AWR"], 0, 0, len(sec["PDISTR"]), 0)
    for k,sub in sorted(sec["PDISTR"].items()):
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
    textOut = []; iline = 1
    for line in text:
        if iline > 99999:
            iline = 1
        textOut.append("{:<66}{:4}{:2}{:3}{:5}\n".format(line, sec["MAT"], sec["MF"], sec["MT"], iline))
        iline += 1
    return "".join(textOut)