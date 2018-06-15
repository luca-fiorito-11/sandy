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

def read_mf34_mt(str_list):
#    str_list = text.splitlines()
    i = 0
    out = {"MAT" : str_list["MAT"].iloc[0],
           "MF" : str_list["MF"].iloc[0],
           "MT" : str_list["MT"].iloc[0]}
#    out = {"MAT" : int(str_list[i][66:70]),
#           "MF" : int(str_list[i][70:72]),
#           "MT" : int(str_list[i][72:75])}
    C, i = read_cont(str_list, i)
    # Subsections are given as dictionary values.
    # Keys are MAT1*100+MT1
    out.update({"ZA" : C.C1, "AWR" : C.C2, "MTL" : C.L2, "SUB" : {}})
    for j in range(C.N2):
        sub = {}
        C, i = read_cont(str_list, i)
        sub.update({"XMF1" : C.C1, "XLFS1" : C.C2, "MAT1" : C.L1, "MT1" : C.L2})
        NC = C.N1
        NI = C.N2
        NCLIST = []
        for k in range(NC):
            C, i = read_cont(str_list, i)
            subsub = {"LTY" : C.L2}
            if subsub["LTY"] == 0:
                L, i = read_list(str_list, i)
                subsub.update({"E1" : L.C1, "E2" : L.C2,
                               "CI" : L.B[:L.N2], "XMTI" : L.B[L.N2:]})
                NCLIST.append(subsub)
            elif subsub["LTY"] in (1,2,3):
                L, i = read_list(str_list, i)
                subsub.update({"E1" : L.C1, "E2" : L.C2, "MATS" : L.L1, "MTS" : L.L2,
                               "XMFS" : L.B[0], "XLFSS" : L.B[1],
                               "EI" : L.B[2:2+L.N2], "WEI" : L.B[2+L.N2:]})
                NCLIST.append(subsub)
            NCLIST.append(subsub)
        sub.update({"NC" : NCLIST})
        NILIST = []
        for k in range(NI):
            L, i = read_list(str_list, i)
            subsub = {"LB" : L.L2}
            if subsub["LB"] in range(5):
                subsub.update({"LT" : L.L1, "NT" : L.NPL, "NP" : L.N2})
                if subsub["LT"] == 0:
                    subsub.update({"Ek" : L.B[::2], "Fk" : L.B[1::2]})
                else:
                    pdb.set_trace()
                    Nk = subsub["NP"] - subsub["LT"]
                    ARRk = L.B[:Nk]
                    ARRl = L.B[Nk:]
                    subsub.update({"Ek" : ARRk[:Nk/2], "Fk" : ARRk[Nk/2:],
                                   "El" : ARRl[:subsub["LT"]], "Fl" : ARRl[subsub["LT"]:]})
            elif subsub["LB"] == 5:
                subsub.update({"LS" : L.L1, "NT" : L.NPL, "NE" : L.N2,
                               "Ek" : L.B[:L.N2], "Fkk" : L.B[L.N2:]})
            elif subsub["LB"] == 6:
                subsub.update({"NT" : L.NPL, "NER" : L.N2, "NEC" : (L.NPL-1)//L.N2})
                subsub.update({"Ek" : L.B[:subsub["NER"]]})
                subsub.update({"El" : L.B[subsub["NER"]:subsub["NER"]+subsub["NEC"]]})
                subsub.update({"Fkl" : L.B[subsub["NER"]+subsub["NEC"]:]})
            elif subsub["LB"] in (8,9):
                subsub.update({"LT" : L.L1, "NT" : L.NPL, "NP" : L.N2})
                subsub.update({"Ek" : L.B[:subsub["NP"]], "Fk" : L.B[subsub["NP"]:]})
            else:
                pdb.set_trace()
            NILIST.append(subsub)
        sub.update({"NI" : NILIST})
        out["SUB"].update({sub["MAT1"]*1000+sub["MT1"] : sub})
    return out