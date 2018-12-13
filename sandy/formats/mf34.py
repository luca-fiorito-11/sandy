# -*- coding: utf-8 -*-
"""
Created on Thu Jun 14 09:23:33 2018

@author: fiorito_l
"""
from .records import *
from .utils import Section
import pdb

def read(text):
    str_list = text.splitlines()
    MAT, MF, MT = read_control(str_list[0])[:3]
    out = {"MAT" : MAT, "MF" : MF, "MT" : MT}
    i = 0
    C, i = read_cont(str_list, i)
    out.update({"ZA" : C.C1, "AWR" : C.C2, "LTT" : C.L2, "NMT1" : C.N2, "REAC" : {}})
    for j in range(out["NMT1"]):
        C, i = read_cont(str_list, i)
        mat1 = C.L1
        mt1 = C.L2
        sub = {"NL" : C.N1, "NL1" : C.N2, "P" : {}}
        nss = C.N1*(C.N1+1)//2 if C.L2 == out["MT"] else C.N1*C.N2
        for k in range(nss):
            C, i = read_cont(str_list, i)
            l = C.L1
            l1 = C.L2
            ni = C.N2
            ssub = {"LCT" : C.N1, "NI" : {}}
            for m in range(ni):
                L, i = read_list(str_list, i)
                sssub = {"LS" : L.L1, "LB" : L.L2, "NT" : L.NPL, "NE" : L.N2}
                if sssub["LB"] in range(5):
                    if sssub["LS"] == 0:
                        sssub.update({"EK" : L.B[::2], "FK" : L.B[1::2]})
                    else:
                        pdb.set_trace()
                        Nk = sssub["NE"] - sssub["LS"]
                        ARRk = L.B[:Nk]
                        ARRl = L.B[Nk:]
                        sssub.update({"EK" : ARRk[:Nk/2], "FK" : ARRk[Nk/2:], "EL" : ARRl[:sssub["LS"]], "FL" : ARRl[sssub["LS"]:]})
                elif sssub["LB"] == 5:
                    sssub.update({"EK" : L.B[:L.N2], "FKK" : L.B[L.N2:]})
                elif sssub["LB"] == 6:
                    sssub.update({"NT" : L.NPL, "NER" : L.N2, "NEC" : (L.NPL-1)//L.N2})
                    sssub.update({"EK" : L.B[:sssub["NER"]]})
                    sssub.update({"EL" : L.B[sssub["NER"]:sssub["NER"]+sssub["NEC"]]})
                    sssub.update({"FKL" : L.B[sssub["NER"]+sssub["NEC"]:]})
                ssub["NI"].update({m : sssub})
            sub["P"].update({(l,l1) : ssub})
        out["REAC"].update({(mat1,mt1) : sub})
    return Section(out)