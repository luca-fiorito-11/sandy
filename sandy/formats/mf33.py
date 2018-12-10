# -*- coding: utf-8 -*-
"""
Created on Thu Jun 14 09:23:33 2018

@author: fiorito_l
"""
from .records import *
from .utils import Section
import pdb
import numpy as np

def read(text):
    str_list = text.splitlines()
    MAT, MF, MT = read_control(str_list[0])[:3]
    out = {"MAT" : MAT, "MF" : MF, "MT" : MT}
    i = 0
    C, i = read_cont(str_list, i)
    # Subsections are given as dictionary values.
    # Keys are MAT1*100+MT1
    out.update({"ZA" : C.C1, "AWR" : C.C2, "MTL" : C.L2, "SUB" : {}})
    for j in range(C.N2):
        C, i = read_cont(str_list, i)
        sub = {"XMF1" : C.C1, "XLFS1" : C.C2, "MAT1" : C.L1, "MT1" : C.L2}
        NC = C.N1
        NI = C.N2
        NCDICT = {}
        for k in range(NC):
            C, i = read_cont(str_list, i)
            subsub = {"LTY" : C.L2}
            if subsub["LTY"] == 0:
                L, i = read_list(str_list, i)
                subsub.update({"E1" : L.C1, "E2" : L.C2,
                               "CI" : L.B[:L.N2], "XMTI" : L.B[L.N2:]})
            elif subsub["LTY"] in (1,2,3):
                L, i = read_list(str_list, i)
                subsub.update({"E1" : L.C1, "E2" : L.C2, "MATS" : L.L1, "MTS" : L.L2,
                               "XMFS" : L.B[0], "XLFSS" : L.B[1],
                               "EI" : L.B[2:2+L.N2], "WEI" : L.B[2+L.N2:]})
            NCDICT.update({k : subsub})
        sub.update({"NC" : NCDICT})
        NIDICT = {}
        for k in range(NI):
            L, i = read_list(str_list, i)
            subsub = {"LB" : L.L2}
            if subsub["LB"] in range(5):
                subsub.update({"LT" : L.L1, "NT" : L.NPL, "NP" : L.N2})
                if subsub["LT"] == 0:
                    subsub["EK"] = L.B[::2]
                    subsub["FK"] = L.B[1::2]
                else: # found in 26-Fe-54g.jeff33, MAT2625/MF33/MT103, LB=4 LT=1
                    nk = subsub["NP"] - subsub["LT"]
                    arrk = L.B[:nk*2]
                    subsub["EK"] = arrk[::2]
                    subsub["FK"] = arrk[1::2]
                    arrl = L.B[nk*2:]
                    subsub["EL"] = arrl[::2]
                    subsub["FL"] = arrl[1::2]
            elif subsub["LB"] == 5:
                subsub.update({"LS" : L.L1, "NT" : L.NPL, "NE" : L.N2,
                               "EK" : L.B[:L.N2], "FKK" : L.B[L.N2:]})
            elif subsub["LB"] == 6:
                subsub.update({"NT" : L.NPL, "NER" : L.N2, "NEC" : (L.NPL-1)//L.N2})
                subsub.update({"EK" : L.B[:subsub["NER"]]})
                subsub.update({"EL" : L.B[subsub["NER"]:subsub["NER"]+subsub["NEC"]]})
                subsub.update({"FKL" : L.B[subsub["NER"]+subsub["NEC"]:]})
            elif subsub["LB"] in (8,9):
                subsub.update({"LT" : L.L1, "NT" : L.NPL, "NP" : L.N2})
                subsub.update({"EK" : L.B[:subsub["NP"]], "FK" : L.B[subsub["NP"]:]})
            else:
                pdb.set_trace()
            NIDICT.update({k : subsub})
        sub.update({"NI" : NIDICT})
        out["SUB"].update({sub["MAT1"]*1000+sub["MT1"] : sub})
    return Section(out)

def read_errorr(text):
    str_list = text.splitlines()
    MAT, MF, MT = read_control(str_list[0])[:3]
    out = {"MAT" : MAT, "MF" : MF, "MT" : MT}
    i = 0
    C, i = read_cont(str_list, i)
    out.update({"ZA" : C.C1, "AWR" : C.C2, "RP" : {}})
    for rp in range(C.N2): # number of reaction pairs
        C, i = read_cont(str_list, i)
        MT1 = C.L2
        NG = C.N2
        M = np.zeros((NG,NG))
        while True:
            L, i = read_list(str_list, i)
            NGCOL = L.L1
            GROW = L.N2
            GCOL = L.L2
            M[GROW-1, GCOL-1:GCOL+NGCOL-1] = L.B
            if GCOL+NGCOL >= NG and GROW >= NG: break
        out["RP"].update({MT1 : M})
    return out