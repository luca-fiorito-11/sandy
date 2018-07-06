# -*- coding: utf-8 -*-
"""
Created on Thu Jun 14 09:23:33 2018

@author: fiorito_l
"""
from ..records import read_cont, read_tab1, read_tab2, read_list, read_control, write_cont, write_tab1, write_tab2, write_list
from ..utils import Section

def read(text):
    str_list = text.splitlines()
    MAT, MF, MT = read_control(str_list[0])[:3]
    out = {"MAT" : MAT, "MF" : MF, "MT" : MT}
    i = 0
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
        for j in range(sub["NE"]):
            T, i = read_tab1(str_list, i)
            distr = {"T" : T.C1, "LT" : T.L1, "NBT" : T.NBT, "INT" : T.INT, "MU" : T.x, "ADISTR" : T.y}
            sub["E"].update({T.C2 : distr})
        out.update({"TAB" : sub})
    return Section(out)

def write(sec):
    if "LCP" in sec and "TAB" in sec:
        sec["LTT"] == 3
        sec["LI"] == 0
    elif "LCP" in sec:
        sec["LTT"] == 1
        sec["LI"] == 0
    elif "TAB" in sec:
        sec["LTT"] == 2
        sec["LI"] == 0
    else:
        sec["LTT"] == 0
        sec["LI"] == 1
    text = write_cont(sec["ZA"], sec["AWR"], 0, sec["LTT"], 0, 0)
    text += write_cont(0, sec["AWR"], sec["LI"], sec["LCT"], 0, 0)
    if sec["LTT"] in (1,3):
        lpc = sec["LPC"]
        text += write_tab2(0, 0, 0, 0, len(lpc["E"]), lpc["NBT"], lpc["INT"])
        for e,sub in lpc["E"].items():
            text += write_list(sub["T"], e, sub["LT"], 0, 0, sub["COEFF"])
    if sec["LTT"] in (2,3):
        tab = sec["TAB"]
        text += write_tab2(0, 0, 0, 0, len(tab["E"]), tab["NBT"], tab["INT"])
        for e,sub in tab["E"].items():
            text += write_tab1(sub["T"], e, sub["LT"], 0, sub["NBT"], sub["INT"], sub["MU"], sub["ADISTR"])
    TextOut = []; iline = 1
    for line in text:
        if iline > 99999:
            iline = 1
        TextOut.append("{:<66}{:4}{:2}{:3}{:5}\n".format(line, sec["MAT"], sec["MF"], sec["MT"], iline))
        iline += 1
    return "".join(TextOut)