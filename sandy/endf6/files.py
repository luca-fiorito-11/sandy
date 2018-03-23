# -*- coding: utf-8 -*-
"""
Created on Mon Jan 16 18:03:13 2017

@author: lfiorito
"""
import sys
import logging
import numpy as np
from sandy.endf6.records import read_cont, read_tab1, read_tab2, read_list, read_text, write_cont, write_tab1, write_list#, add_records
import matplotlib.pyplot as plt
import pandas as pd
import pdb


def split(file):
    """
    Split ``ENDF-6`` file  into MFMT sections.
    """
    import re
    pattern = ".{74}0.{5}\n?"
    text = open(file).read()
    U = re.split(pattern, text)
    return list(filter(None, U)) # remove empty lines


"""
        Found error in:
            - n-17-Cl-035.jeff32
            - n-3-Li-007.jeff32
            - n-63-Eu-152.jeff32
            - n-63-Eu-153.jeff32
            - n-64-Gd-155.jeff32
            - n-77-Ir-193.jeff32
            - n-90-Th-229.jeff32
            - n-94-Pu-238.jeff32
            - n-94-Pu-241.jeff32
            - n-94-Pu-242.jeff32
            - n-97-Bk-250.jeff32
            - n-98-Cf-254.jeff32
"""

def process_endf_section(text):
    mf = int(text[70:72])
    mt = int(text[72:75])
#    if mf not in (5,35):
#        return None
    if mf ==1 and mt == 451:
        return read_mf1_mt451(text)
    elif mf ==1 and mt in (452, 455, 456):
        return read_mf1_nubar(text)
    elif mf == 3:
        return read_mf3_mt(text)
    elif mf == 4:
        return read_mf4_mt(text)
    elif mf == 5:
        return read_mf5_mt(text)
    elif mf == 8 and mt == 457:
        return read_mf8_mt457(text)
    elif mf == 31 or mf == 33:
        return read_mf33_mt(text)
    elif mf == 35:
        return read_mf35_mt(text)
    else:
        return None

def endf2df(file):
    tape = pd.DataFrame([[int(x[66:70]), int(x[70:72]), int(x[72:75]), x] for x in split(file)],
            columns=('MAT', 'MF', 'MT','TEXT'))
    tape = tape.set_index(['MAT','MF','MT']).sort_index() # Multi-indexing
    tape['DATA'] = tape['TEXT'].apply(process_endf_section)
#    pool = mp.Pool(processes=settings.args.processes)
#    AAA = pool.map( e6.process_endf_section, tape['TEXT'].tolist())
#    data = [ pool.apply( e6.process_endf_section, args=(x)) for x in tape['TEXT'].tolist() ]
    return tape


def read_mf1_mt451(text):
    str_list = text.splitlines()
    i = 0
    out = {"MAT" : int(str_list[i][66:70]),
           "MF" : int(str_list[i][70:72]),
           "MT" : int(str_list[i][72:75])}
    C, i = read_cont(str_list, i)
    out.update({"ZA" : C.C1, "AWR" : C.C2, "LRP" : C.L1, "LFI" : C.L2, "NLIB" :C.N1, "NMOD" : C.N2})
    C, i = read_cont(str_list, i)
    out.update({"ELIS" : C.C1, "STA" : C.C2, "LIS" : C.L1, "LISO" : C.L2, "NFOR" : C.N2})
    C, i = read_cont(str_list, i)
    out.update({"AWI" : C.C1, "EMAX" : C.C2, "LREL" : C.L1, "NSUB" : C.N1, "NVER" : C.N2})
    C, i = read_cont(str_list, i)
    out.update({"TEMP" : C.C1, "LDRV" : C.L1, "NWD" : C.N1, "NXC" : C.N2})
    TEXT = []
    for j in range(out["NWD"]):
        T, i = read_text(str_list, i)
        TEXT.append(T)
    out.update({ "TEXT" : TEXT })
    out["Z"] = int(TEXT[0][:3])
    out["SYM"] = TEXT[0][4:6].rstrip()
    out["A"] = int(TEXT[0][7:10])
    out["M"] =  'g' if TEXT[0][10:11] is ' ' else TEXT[0][10:11].lower()
    out['ALAB'] = TEXT[0][11:22]
    out['EDATE'] = TEXT[0][22:32]
    out['AUTH'] = TEXT[0][33:66]
    out['REF'] = TEXT[1][1:22]
    out['DDATE'] = TEXT[1][22:32]
    out['RDATE'] = TEXT[1][33:43]
    out['ENDATE'] = TEXT[1][55:63]
    out['LIBVER'] = TEXT[2][:22].strip('- ')
    out['SUB'] = TEXT[3].strip('- ')
    out['FOR'] = TEXT[4].strip('- ')
    return out

def read_mf1_nubar(text):
    str_list = text.splitlines()
    i = 0
    out = {"MAT" : int(str_list[i][66:70]),
           "MF" : int(str_list[i][70:72]),
           "MT" : int(str_list[i][72:75])}
    C, i = read_cont(str_list, i)
    out.update({"ZA" : C.C1, "AWR" : C.C2, "LDG" : C.L1, "LNU" : C.L2})
    if out["MT"] == 455:
        if out["LDG"] == 0:
            L, i = read_list(str_list, i)
            out.update({ "NNF" : L.NPL, "LAMBDAS" : L.B })
        elif out["LDG"] == 1:
            # None found in JEFF33 and ENDFB8, hence not implemented
            pass
    if out["LNU"] == 1:
        # None found in JEFF33 and ENDFB8 neither for MT455 nor for MT456
        L, i = read_list(str_list, i)
        out.update({ "NC" : L.NPL, "C" : L.B})
    else:
        # JEFF33 and ENDFB8 only have lin-lin interpolation schemes
        T, i = read_tab1(str_list, i)
        out.update({"NBR" : T.NBT, "INT" : T.INT})
        out["NUBAR"] = pd.Series(T.y, index = T.x, name = out["MT"]).rename_axis("E")
    return out

def write_mf1_nubar(tape):
    for (mat,mf,mt),df in tape.loc[(slice(None),1),:].iterrows():
        TEXT = write_cont(df.DATA["ZA"], df.DATA["AWR"], df.DATA["LDG"], df.DATA["LNU"], 0, 0)
        if df.DATA["MT"] == 455:
            if df.DATA["LDG"] == 0:
                TEXT += write_list(0, 0, 0, 0, 0, df.DATA["LAMBDAS"])
            elif df.DATA["LDG"] == 1:
                # Not found in JEFF33 and ENDFB8, hence not implemented
                pass
        if df.DATA["LNU"] == 1:
            TEXT += write_list(0, 0, 0, 0, 0, df.DATA["C"])
        else:
            TEXT += write_tab1(0, 0, 0, 0, df.DATA["NBR"], df.DATA["INT"],
                           df.DATA["NUBAR"].index, df.DATA["NUBAR"])
        TEXT = [ "{:<66}{:4}{:2}{:3}{:5}".format(l, mat, mf, mt, i+1) for i,l in enumerate(TEXT) ]
        tape.at[(mat,mf,mt),'TEXT'] = "\n".join(TEXT) + '\n'
    return tape

def read_mf3_mt(text):
    str_list = text.splitlines()
    i = 0
    out = {"MAT" : int(str_list[i][66:70]),
           "MF" : int(str_list[i][70:72]),
           "MT" : int(str_list[i][72:75])}
    C, i = read_cont(str_list, i)
    out.update({"ZA" : C.C1, "AWR" : C.C2})
    T, i = read_tab1(str_list, i)
    out.update({"QM" : T.C1, "QI" : T.C2, "LR" : T.L2, "NBR" : T.NBT, "INT" : T.INT})
    XS = pd.Series(T.y, index = T.x, name = out["MT"]).rename_axis("E")
    out.update({"XS" : XS})
    return out

def write_mf3_mt(tape):
    for (mat,mf,mt),df in tape.loc[(slice(None),3),:].iterrows():
        TEXT = write_cont(df.DATA["ZA"], df.DATA["AWR"], 0, 0, 0, 0)
        TEXT += write_tab1(df.DATA["QM"], df.DATA["QI"], 0, df.DATA["LR"],
                           df.DATA["NBR"], df.DATA["INT"],
                           df.DATA["XS"].index, df.DATA["XS"])
        TEXT = [ "{:<66}{:4}{:2}{:3}{:5}".format(l, mat, mf, mt, i+1) for i,l in enumerate(TEXT) ]
        tape.at[(mat,mf,mt),'TEXT'] = "\n".join(TEXT) + '\n'
    return tape

def read_mf4_mt(text):
    str_list = text.splitlines()
    i = 0
    out = {"MAT" : int(str_list[i][66:70]),
           "MF" : int(str_list[i][70:72]),
           "MT" : int(str_list[i][72:75])}
    C, i = read_cont(str_list, i)
    out.update({"ZA" : C.C1, "AWR" : C.C2, "LTT" : C.L2})
    if out["LTT"] in (1,3):
        C, i = read_cont(str_list, i)
        sub = {"LI" : C.L1, "LCT" : C.L2}
        T2, i = read_tab2(str_list, i)
        sub.update({"NE" : T2.NZ, "NBR" : T2.NBT, "INT" : T2.INT})
        NLMAX = 0
        E = []
        P = np.zeros((sub["NE"], 64)) # 64 is the largest NL allowed
        for j in range(sub["NE"]):
            L, i = read_list(str_list, i)
            E.append(L.C2)
            P[j,:L.NPL] = L.B
            NLMAX = max(NLMAX, L.NPL)
        columns = ["P{}".format(j+1) for j in range(NLMAX)]
        sub.update({"P" : pd.DataFrame(P[:,:NLMAX], index=E, columns=columns)})
        out.update({"LPC" : sub})
    if out["LTT"] in (2,3):
        C, i = read_cont(str_list, i)
        sub = {"LI" : C.L1, "LCT" : C.L2}
        T2, i = read_tab2(str_list, i)
        sub.update({"NE" : T2.NZ, "NBR" : T2.NBT, "INT" : T2.INT})
        TL = []
        for i in range(sub["NE"]):
            T1, i = read_tab1(str_list, i)
            TL.append(T1)
            sub.update({"TL" : TL})
        out.update({"TPD" : sub})
    if out["LTT"] == 0:
        C, i = read_cont(str_list, i)
        out.update({"ISO" : {"LI" : C.L1, "LCT" : C.L2}})
    return out

def read_mf5_mt(text):
    from sandy.functions import union_grid
    str_list = text.splitlines()
    i = 0
    out = {"MAT" : int(str_list[i][66:70]),
           "MF" : int(str_list[i][70:72]),
           "MT" : int(str_list[i][72:75])}
    C, i = read_cont(str_list, i)
    # subsections for partial energy distributions are given in a list
    out.update({"ZA" : C.C1, "AWR" : C.C2, "NK" : C.N1, "SUB" : [] })
    for j in range(out["NK"]):
        Tp, i = read_tab1(str_list, i)
        P = pd.Series(Tp.y, index = Tp.x, name = "p").rename_axis("E")
        sub = { "LF" : Tp.L2, "NBR_p" : Tp.NBT, "INT_p" : Tp.INT, "P" : P }
        if sub["LF"] == 5:
            Ttheta, i = read_tab1(str_list, i)
            Tg, i = read_tab1(str_list, i)
            sub.update({ "Ttheta" : Ttheta, "Tg" : Tg })
        elif sub["LF"] in (7,9):
            Ttheta, i = read_tab1(str_list, i)
            sub.update({ "Ttheta" : Ttheta})
        elif sub["LF"] == 11:
            Ta, i = read_tab1(str_list, i)
            Tb, i = read_tab1(str_list, i)
            sub.update({ "Ta" : Ta, "Tb" : Tb })
        elif sub["LF"] == 12:
            TTm, i = read_tab1(str_list, i)
            sub.update({ "TTm" : TTm })
        elif sub["LF"] == 1:
            T2, i = read_tab2(str_list, i)
            sub.update({ "NBT_Ein" : T2.NBT, "INT_Ein" : T2.INT, "Ein" : {} })
            for k in range(T2.NZ):
                T1, i = read_tab1(str_list, i)
                distr = pd.Series(T1.y, index = T1.x, name=T1.C2).rename_axis("Eout")
                sub["Ein"].update({ T1.C2 : {"PDF" : distr, "NBT" : T1.NBT, "INT" : T1.INT}})
        out["SUB"].append(sub)
    return out

def read_mf8_mt457(text):
    str_list = text.splitlines()
    i = 0
    out = {"MAT" : int(str_list[i][66:70]),
           "MF" : int(str_list[i][70:72]),
           "MT" : int(str_list[i][72:75])}
    C, i = read_cont(str_list, i)
    out.update({"ZA" : C.C1, "AWR" : C.C2, "LIS" : C.L1, "LISO" : C.L2, "NST" :C.N1, "NSP" : C.N2})
    L, i = read_list(str_list, i)
    out.update({"HL" : L.C1, "DHL" : L.C2, "E" : L.B[::2], "DE" : L.B[1::2]})
    return out

def read_mf33_mt(text):
    str_list = text.splitlines()
    i = 0
    out = {"MAT" : int(str_list[i][66:70]),
           "MF" : int(str_list[i][70:72]),
           "MT" : int(str_list[i][72:75])}
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

def read_mf34_mt(text):
    str_list = text.splitlines()
    i = 0
    out = {"MAT" : int(str_list[i][66:70]),
           "MF" : int(str_list[i][70:72]),
           "MT" : int(str_list[i][72:75])}
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

def read_mf35_mt(text):
    str_list = text.splitlines()
    i = 0
    out = {"MAT" : int(str_list[i][66:70]),
           "MF" : int(str_list[i][70:72]),
           "MT" : int(str_list[i][72:75])}
    C, i = read_cont(str_list, i)
    out.update({ "ZA" : C.C1, "AWR" : C.C2, "NK" : C.N1, "SUB" : []})
    for k in range(out["NK"]):
        L, i = read_list(str_list, i)
        out["SUB"].append({ "Elo" : L.C1, "Ehi" : L.C2, "NE" : L.N2, "Ek" : L.B[:L.N2],
           "Fkk" : L.B[L.N2:] })
    return out

def pandas_interpolate(df, interp_column, method='zero', axis='both'):
    # interp_column is a list
    dfout = df.copy()
    if axis in ['rows', 'both']:
        ug = np.unique(list(dfout.index) + interp_column)
        dfout = dfout.reindex(ug)
        dfout = dfout.interpolate(method=method)
        dfout = dfout.reindex(interp_column)
    # Now transpose columns to rows and reinterpolate
    if axis in ['cols', 'both']:
        dfout = dfout.transpose()
        ug = np.unique(list(df.index) + interp_column)
        dfout = dfout.reindex(ug)
        dfout = dfout.interpolate(method=method)
        dfout = dfout.reindex(interp_column)
        dfout = dfout.transpose()
    dfout.fillna(0, inplace=True)
    return dfout


def extract_xs(tape):
    XS = {} # DIctionary by MAT
#    XS = pd.DataFrame() # dataframe with MAT as index
    for mat in np.unique(tape.index.get_level_values("MAT")):
        xsdf = pd.DataFrame(tape.loc[mat,3,1].DATA["XS"])
        if len(np.unique(xsdf.index)) != len(xsdf.index):
            raise NotImplementedError()
#        xsdf.columns = xsdf.columns.droplevel() # keep only MT, drop MAT
#        xsdf.reset_index(inplace=True)
#        xsdf["MAT"] = mat
#        xsdf = xsdf.set_index(['MAT','E']).sort_index()
        # No interpolation is done because xs are on unionized grid
        for chunk in tape.query('MF==3 & MT!=1').DATA:
            df = pd.DataFrame(chunk["XS"])
#            df.columns = df.columns.droplevel() # keep only MT, drop MAT
#            df.reset_index(inplace=True)
#            df["MAT"] = mat
#            df = df.set_index(['MAT','E']).sort_index()
            # PROBLEM! "add" does not work for duplicated energy points
            if len(np.unique(df.index)) != len(df.index):
                raise NotImplementedError()
            xsdf = xsdf.add(df, fill_value=0)
        xsdf.fillna(0, inplace=True)
        xsdf.columns = pd.Index(xsdf.columns, name="MT")
#        XS = pd.concat((XS, xsdf)) # To be tested
        XS.update({ mat : xsdf })
    return XS

def extract_chi(tape):
    """
    Extract chi cov for all MAT,MT,SUB found in tape.
    Return a df with MAT,MT,SUB as index and COV as value
    Each COV is a df with Ein on rows and Eout on columns.
    """
    from sandy.functions import union_grid
    DictDf = {}
    for chunk in tape.query('MF==5').DATA:
        for k,sub in enumerate(chunk["SUB"]):
            if sub["LF"] != 1:
                continue
            if list(filter(lambda x:x["INT"] != [2], sub["Ein"].values())):
#                print("WARNING: found non-linlin interpolation, skip energy distr. for MAT {}, MT {}, subsec {}".format(chunk["MAT"],chunk["MT"],k))
                continue
            chi_dict =  { ein : ssub["PDF"] for ein,ssub in sorted(sub["Ein"].items()) }
            # merge chi_E(E') distributions on df[E,E'] with unique E' grid
            chi = pd.DataFrame.from_dict(chi_dict).interpolate(method="slinear").fillna(0).transpose()
            # include default points in E grid and interpolate
            eg_new = list(union_grid(chi.index, np.logspace(-5, 7, 13)))
            chi = pandas_interpolate(chi, eg_new, method='slinear', axis='rows')
            chi.index.name = "Ein"
            chi.columns.name = "Eout"
            DictDf.update({ (chunk["MAT"],chunk["MT"],k) : chi })
    DfChi = pd.DataFrame.from_dict(DictDf, orient='index')
    DfChi.index.names = ["MAT", "MT", "K"]
    DfChi.columns.names = ["COV"]
    return DfChi


def extract_cov33(tape, mt=[102]):
    from sandy.cov import triu_matrix
    from sandy.functions import union_grid
    columns = ('MAT', 'MT', 'MAT1', 'MT1', 'COV')
    covdf = pd.DataFrame(columns=columns)
    for chunk in tape.query('MF==33 | MF==31').DATA:
        for sub in chunk["SUB"].values():
            covs = []
            if len(sub["NI"]) == 0:
                continue
            for nisec in sub["NI"]:
                if nisec["LB"] == 5:
                    Fkk = np.array(nisec["Fkk"])
                    if nisec["LS"] == 0: # to be tested
                        cov = Fkk.reshape(nisec["NE"]-1, nisec["NE"]-1)
                    else:
                        cov = triu_matrix(Fkk, nisec["NE"]-1)
                    # add zero row and column at the end of the matrix
                    cov = np.insert(cov, cov.shape[0], [0]*cov.shape[1], axis=0)
                    cov = np.insert(cov, cov.shape[1], [0]*cov.shape[0], axis=1)
                    covs.append(pd.DataFrame(cov, index=nisec["Ek"], columns=nisec["Ek"]))
                elif nisec["LB"] == 1:
                    cov = np.diag(nisec["Fk"])
                    covs.append(pd.DataFrame(cov, index=nisec["Ek"], columns=nisec["Ek"]))
                elif nisec["LB"] == 2:
                    f = np.array(nisec["Fk"])
                    cov = f*f.reshape(-1,1)
                    covs.append(pd.DataFrame(cov, index=nisec["Ek"], columns=nisec["Ek"]))
                elif nisec["LB"] == 6:
                    cov = np.array(nisec["Fkl"]).reshape(nisec["NER"]-1, nisec["NEC"]-1)
                    # add zero row and column at the end of the matrix
                    cov = np.insert(cov, cov.shape[0], [0]*cov.shape[1], axis=0)
                    cov = np.insert(cov, cov.shape[1], [0]*cov.shape[0], axis=1)
                    covs.append(pd.DataFrame(cov, index=nisec["Ek"], columns=nisec["El"]))
            if len(covs) == 0:
                continue
            # Al union covariance matrices have the same grid (uxx) on both axis
            uxx = union_grid(*[list(cov.index) + list(cov.columns) for cov in covs])
            cov = np.sum([ pandas_interpolate(cov, list(uxx)).as_matrix() for cov in covs ], 0)
            cov = pd.DataFrame(cov, index=uxx, columns=uxx)
            objects = chunk['MAT'], chunk['MT'], sub['MAT1'], sub['MT1'], cov
            covdf = covdf.append(dict(zip(columns, objects)), ignore_index=True)
    covdf = covdf.set_index(['MAT','MT','MAT1','MT1']).sort_index() # Multi-indexing
    return None if covdf.empty else covdf




def extract_cov35(tape):
    from sandy.cov import triu_matrix, corr2cov
    # Get covariances (df) from each MAT, MT and Erange (these are the keys) in a dictionary.
    DictCov = {}
    for chunk in tape.query('MF==35').DATA:
        for sub in chunk["SUB"]:
            # Ek grid is one unit longer than covariance.
            Ek = np.array(sub["Ek"])
            Fkk = np.array(sub["Fkk"])
            NE = sub["NE"]
            cov = triu_matrix(Fkk, NE-1)
            # Normalize covariance matrix dividing by the energy bin.
            dE = 1./(Ek[1:]-Ek[:-1])
            cov = corr2cov(cov, dE)
            # Add zero row and column at the end of the matrix
            cov = np.insert(cov, cov.shape[0], [0]*cov.shape[1], axis=0)
            cov = np.insert(cov, cov.shape[1], [0]*cov.shape[0], axis=1)
            cov = pd.DataFrame(cov, index=Ek, columns=Ek)
            cov.index.name = "Ek"; cov.columns.name = "El"
            DictCov.update({ (chunk["MAT"], chunk["MT"], sub["Elo"], sub["Ehi"]) :
                cov })
    # Collect covs in a list to pass to block_diag (endf6 allows only covs in diag)
    # Recreate indexes with lists of mat, mt, elo, ehi and e.
    from scipy.linalg import block_diag
    covs = []; mat = []; mt = []; elo = []; ehi = []; e= []
    for k,c in sorted(DictCov.items()):
        covs.append(c)
        mat.extend([k[0]]*len(c.index))
        mt.extend([k[1]]*len(c.index))
        elo.extend([k[2]]*len(c.index))
        ehi.extend([k[3]]*len(c.index))
        e.extend(c.index)
    DfCov = block_diag(*covs)
    DfCov = pd.DataFrame(DfCov)
    DfCov['MAT'] = pd.Series(mat)
    DfCov['MT'] = pd.Series(mt)
    DfCov['Elo'] = pd.Series(elo)
    DfCov['Ehi'] = pd.Series(ehi)
    DfCov['E'] = pd.Series(e)
    DfCov.set_index(['MAT', 'MT', 'Elo', 'Ehi', 'E'], inplace=True)
    DfCov.columns = DfCov.index
    return DfCov

def merge_covs(covdf):
    query_diag = 'MT==MT1 & (MAT1==0 | MAT==MAT1)'
    query_offdiag = 'MT!=MT1 | MAT1!=0'
    diags = covdf.query(query_diag)
    # reset indices to be consistent with enumerate in the next loop
    diags.reset_index(inplace=True)
    # Process diagonal blocks
    # This part is extracted from scipy.linalg.block_diag
    shapes = np.array([a.COV.shape for i,a in diags.iterrows()])
    ndim = np.sum(shapes, axis=0)[0]
    C = np.zeros((ndim,ndim)); E = np.zeros(ndim)
    MATS = np.zeros(ndim, dtype=int); MTS = np.zeros(ndim, dtype=int)
    r, c = 0, 0
    beg, end = [], []
    for i, (rr, cc) in enumerate(shapes):
        d = diags.iloc[i]
        C[r:r + rr, c:c + cc] = d.COV
        E[r:r + rr] = d.COV.index
        MATS[r:r + rr] = d.MAT
        MTS[r:r + rr] = d.MT
        beg.append(r)
        end.append(r + rr)
        r += rr
        c += cc
    diags = diags.assign(BEG = beg)
    diags = diags.assign(END = end)
    # reset multindex to use loc method
    diags = diags.set_index(['MAT','MT','MAT1','MT1']).sort_index()
    # Process off-diagonal blocks
    for (mat,mt,mat1,mt1),row in covdf.query(query_offdiag).iterrows():
        cov = row.COV
        # interpolate x axis (rows)
        covk = diags.loc[mat,mt,0,mt]
        Ek = list(covk.COV.index)
        cov = pandas_interpolate(cov, Ek, method='zero', axis='rows')
        # interpolate y axis (cols)
        if mat1 == 0:
            mat1 = mat
        covl = diags.loc[mat1,mt1,0,mt1]
        El = list(covl.COV.index)
        cov = pandas_interpolate(cov, El, method='zero', axis='cols')
        C[covk.BEG:covk.END,covl.BEG:covl.END] = cov
        C[covl.BEG:covl.END,covk.BEG:covk.END,] = cov.T
    C = pd.DataFrame(C, index=[MATS,MTS,E], columns=[MATS,MTS,E])
    C.index.names = ["MAT", "MT", "E"]
    return C

def update_xs(tape, xsdf):
    for mat,df in xsdf.items():
        for mt in df:
            if (mat, 3, mt) not in tape.index:
                continue
            XS = df[mt]
            D = tape.DATA.loc[mat,3,mt]
            # Assume all xs have only 1 interpolation region and it is linear
#            if len(D["INT"]) != 1:
#                raise NotImplementedError("Cannot update xs with more than 1 interp. region")
#            if D["INT"][0] != 1:
#                raise NotImplementedError("Cannot update xs with non-linear interpolation")
            D["XS"] = XS
            D["NBR"] = [len(XS)]
            D["INT"] = [2]
    return tape

def write_tape(tape, file, title=" "*66):
    string = "{:<66}{:4}{:2}{:3}{:5}\n".format(title, 1, 0, 0, 0)
    for mat,dfmat in tape.groupby('MAT', sort=True):
        for mf,dfmf in dfmat.groupby('MF', sort=True):
            for mt,dfmt in dfmf.groupby('MT', sort=True):
                for text in dfmt.TEXT:
                    string += text.encode('ascii', 'replace').decode('ascii')
                string += "{:<66}{:4}{:2}{:3}{:5}\n".format(*write_cont(*[0]*6), mat, mf, 0, 99999)
            string += "{:<66}{:4}{:2}{:3}{:5}\n".format(*write_cont(*[0]*6), mat, 0, 0, 0)
        string += "{:<66}{:4}{:2}{:3}{:5}\n".format(*write_cont(*[0]*6), 0, 0, 0, 0)
    string += "{:<66}{:4}{:2}{:3}{:5}".format(*write_cont(*[0]*6), -1, 0, 0, 0)
    with open(file, 'w', encoding="ascii") as f:
        f.write(string)

def write_decay_data_csv(tape, filename):
    df = tape.query("MF==8 & MT==457")
    df['Z'] = df.DATA.apply(lambda x : int(x['ZA']//1000))
    df['A'] = df.DATA.apply(lambda x : int(x['ZA'] - x['ZA']//1000*1000))
    df['M'] = df.DATA.apply(lambda x : 'g' if x['LISO'] == 0 else 'm' if x['LISO'] == 1 else 'n')
    df['HL'] = df.DATA.apply(lambda x : x['HL'])
    df['DHL'] = df.DATA.apply(lambda x : x['DHL'])
    df['ELP'] = df.DATA.apply(lambda x : x['E'][0])
    df['DELP'] = df.DATA.apply(lambda x : x['DE'][0])
    df['EEM'] = df.DATA.apply(lambda x : x['E'][1])
    df['DEEM'] = df.DATA.apply(lambda x : x['DE'][2])
    df['EHP'] = df.DATA.apply(lambda x : x['E'][2])
    df['DEHP'] = df.DATA.apply(lambda x : x['DE'][2])
    df[['Z','A','M','HL','DHL','ELP','DELP','EEM','DEEM','EHP','DEHP']].to_csv(filename, index=False)