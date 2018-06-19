# -*- coding: utf-8 -*-
"""
Created on Thu Jun 14 09:23:33 2018

@author: fiorito_l
"""
from ..records import read_cont, read_tab1, read_control
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
def read_mf2_mt151(text):
    str_list = split_endf(text)
    i = 0
    out = {"MAT" : str_list["MAT"].iloc[0],
           "MF" : str_list["MF"].iloc[0],
           "MT" : str_list["MT"].iloc[0]}
    C, i = read_cont(str_list, i)
    out.update({ "ZA" : C.C1, "AWR" : C.C2, "NIS" : C.N1, "ZAI" : {} })
    for i_iso in range(out["NIS"]): # LOOP ISOTOPES
        C, i = read_cont(str_list, i)
        zai = { "ZAI" : C.C1, "ABN" : C.C2, "LFW" : C.L2, "NER" : C.N1, "ERANGE" : {} }
        for i_erange in range(zai["NER"]): # LOOP ENERGY RANGES
            C, i = read_cont(str_list, i)
            sub = {"EL" : C.C1, "EH" : C.C2, "LRU" : C.L1, "LRF" : C.L2, "NRO" : C.N1, "NAPS" : C.N2}
            if sub["NRO"] != 0: # Tabulated scattering radius
                T, i = read_tab1(str_list, i)
                sub.update({"NBT" : T.NBT, "INT" : T.INT})
                sub["AP"] = pd.Series(T.y, index = T.x).rename_axis("E")
            if sub["LRU"] == 0: # ONLY SCATTERING RADIUS
                C, i = read_cont(str_list, i)
                sub.update({"SPI" : C.C1, "SR" : C.C2, "NLS" : C.N1})
            if sub["LRU"] == 1: # RESOLVED RESONANCES
                if sub["LRF"] in (1,2): # BREIT-WIGNER
                    C, i = read_cont(str_list, i)
                    sub.update({"SPI" : C.C1, "SR" : C.C2, "NLS" : C.N1})
                    L, i = read_list(str_list, i)
                elif sub["LRF"] == 3: # REICH-MOORE
                    C, i = read_cont(str_list, i)
                    sub.update({"SPI" : C.C1, "SR" : C.C2, "LAD" : C.L1, "NLS" : C.N1, "NLSC" : C.N2})
                    L, i = read_list(str_list, i)
                elif sub["LRF"] == 4: # ADLER-ADLER
                    sys.exit("ERROR: SANDY cannot read resonance parameters in Adler-Adler formalism")
                elif sub["LRF"] == 5: # GENERAL R-MATRIX
                    sys.exit("ERROR: General R-matrix formalism no longer available in ENDF-6")
                elif sub["LRF"] == 6: # HYBRID R-FUNCTION
                    sys.exit("ERROR: Hybrid R-function formalism no longer available in ENDF-6")
                elif sub["LRF"] == 7: # LIMITED R-MATRIX
                    C, i = read_cont(str_list, i)
                    sub.update({"IFG" : C.L1, "KRM" : C.L2, "NJS" : C.N1, "KRL" : C.N2})
                    for j in range(sub["NJS"]):
                        L1, i = read_list(str_list, i)
                        L2, i = read_list(str_list, i)
                        L3, i = read_list(str_list, i)
            elif sub["LRU"] == 2: # UNRESOLVED RESONANCES
                if sub["LRF"] == 1 and out["LFW"] == 0: # CASE A
                    C, i = read_cont(str_list, i)
                    sub.update({"SPI" : C.C1, "SR" : C.C2, "LSSF" : C.L1, "NLS" : C.N1})
                    for k in range(sub["NLS"]):
                        L, i = read_list(str_list, i)
                elif sub["LRF"] == 1 and out["LFW"] == 1: # CASE B
                    C, i = read_cont(str_list, i)
                    sub.update({"SPI" : C.C1, "SR" : C.C2, "LSSF" : C.L1, "NE" : C.N1, "NLS" : C.N2})
                    L, i = read_list(str_list, i)
                    for k in range(sub["NLS"]):
                        C, i = read_cont(str_list, i)
                        for l in range(C.N1):
                            L, i = read_list(str_list, i)
                elif sub["LRF"] == 2: # CASE C
                    C, i = read_cont(str_list, i)
                    sub.update({"SPI" : C.C1, "SR" : C.C2, "LSSF" : C.L1, "NLS" : C.N1})
                    for k in range(sub["NLS"]):
                        C, i = read_cont(str_list, i)
                        for l in range(C.N1):
                            L, i = read_list(str_list, i)
            zai["ERANGE"].update({ (sub["EL"],sub["EH"]) : sub })
        out["ZAI"].update({ zai["ZAI"] : zai })
    return out