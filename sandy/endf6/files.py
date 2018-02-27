# -*- coding: utf-8 -*-
"""
Created on Mon Jan 16 18:03:13 2017

@author: lfiorito
"""
import sys
import logging
import numpy as np
from records import read_cont, read_tab1, read_list, read_text, write_cont, write_tab1#, add_records
import matplotlib.pyplot as plt
import pandas as pd
import pdb
import copy


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

def process_section(text):
    mf = int(text[70:72])
    mt = int(text[72:75])
    if mf == 3:
        return read_mf3_mt(text)
    elif mf == 33:
        return read_mf33_mt(text)
    else:
        return None

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
#    XS = pd.Series(T.y, index = T.x, name = (out["MAT"],out["MT"])).rename_axis("E")
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
    for mat in np.unique(tape.index.get_level_values(0)):
        xsdf = pd.DataFrame(tape.loc[mat,3,1].DATA["XS"])
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
            xsdf = xsdf.add(df, fill_value=0)
        xsdf.fillna(0, inplace=True)
#        XS = pd.concat((XS, xsdf)) # To be tested
        XS.update({ mat : xsdf })
    return XS

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


def merge_covs(covdf):
    query_diag = 'MT==MT1 & (MAT1==0 | MAT==MAT1)'
    query_offdiag = 'MT!=MT1 | MAT1!=0'
    diags = covdf.query(query_diag)
    # reset indices to be consistent with enumerate in the next loop
    diags.reset_index(inplace=True)
    # Process diagonal blocks
    # This part is extracted from scipy.linalg.diagblock
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
    diags["BEG"] = beg
    diags["END"] = end
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
    C.sort_index(inplace=True)
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
    string += "{:<66}{:4}{:2}{:3}{:5}".format(*write_cont(*[0]*6), mat, 0, 0, -1)
    with open(file, 'w', encoding="ascii") as f:
        f.write(string)



file = "H1.txt"
#file = "26-Fe-56g.jeff33"
tape = pd.DataFrame([[int(x[66:70]), int(x[70:72]), int(x[72:75]), x, int(x[66:70])*100000+int(x[70:72])*1000+int(x[72:75])] for x in split(file)],
        columns=('MAT', 'MF', 'MT','TEXT', 'ID'))
tape = tape.set_index(['MAT','MF','MT']).sort_index() # Multi-indexing
tape['DATA'] = tape['TEXT'].apply(process_section)

df_cov_xs = merge_covs( extract_cov33(tape) )
df_xs = extract_xs(tape)

from sandy.cov import Cov
NSMP = 100
# Must add names to indexes
df_samples_xs = pd.DataFrame( Cov(df_cov_xs.as_matrix()).sampling(NSMP) + 1 , index = df_cov_xs.index )
unique_ids = list(set(zip(df_cov_xs.index.get_level_values(0), df_cov_xs.index.get_level_values(1))))
#idx = pd.IndexSlice
#df_samples.loc[idx[:, :, ['C1', 'C3']], idx[:, 'foo']]
for ismp in range(NSMP):
#    ixs = df_xs.copy()
    ixs = copy.deepcopy(df_xs)
    for mat, df_samples_bymat in df_samples_xs.groupby(level=0, sort=True):
        for mt, P in df_samples_bymat.groupby(level=1, sort=True):
            if mat not in ixs:
                continue
            if mt not in ixs[mat]:
                continue
            XS = ixs[mat][mt]
            P = pandas_interpolate(P.loc[mat,mt][ismp],
                                   list(XS.index.get_level_values("E")),
                                   axis="rows")
            XS = XS.multiply(P)
            # Negative values are set to mean
            ixs[mat][mt][XS > 0] = XS[XS > 0]
        # Redundant XS
        columns = ixs[mat].columns
        daughters = [ x for x in range(800,850) if x in ixs[mat]]
        if daughters:
            ixs[mat][107] = ixs[mat][daughters].sum(axis=1)
        daughters = [ x for x in range(750,800) if x in ixs[mat]]
        if daughters:
            ixs[mat][106] = ixs[mat][daughters].sum(axis=1)
        daughters = [ x for x in range(700,750) if x in ixs[mat]]
        if daughters:
            ixs[mat][105] = ixs[mat][daughters].sum(axis=1)
        daughters = [ x for x in range(650,700) if x in ixs[mat]]
        if daughters:
            ixs[mat][104] = ixs[mat][daughters].sum(axis=1)
        daughters = [ x for x in range(600,650) if x in ixs[mat]]
        if daughters:
            ixs[mat][103] = ixs[mat][daughters].sum(axis=1)
        daughters = [ x for x in range(102,118) if x in ixs[mat]]
        if daughters:
            ixs[mat][101] = ixs[mat][daughters].sum(axis=1)
        daughters = [ x for x in (19,20,21,38) if x in ixs[mat]]
        if daughters:
            ixs[mat][18] = ixs[mat][daughters].sum(axis=1)
        daughters = [ x for x in (18,101) if x in ixs[mat]]
        if daughters:
            ixs[mat][27] = ixs[mat][daughters].sum(axis=1)
        daughters = [ x for x in range(50,92) if x in ixs[mat]]
        if daughters:
            ixs[mat][4] = ixs[mat][daughters].sum(axis=1)
        daughters = [ x for x in (4,5,11,*range(16,19),*range(22,38),41,42,44,45) if x in ixs[mat]]
        if daughters:
            ixs[mat][3] = ixs[mat][daughters].sum(axis=1)
        daughters = [ x for x in (2,3) if x in ixs[mat]]
        if daughters:
            ixs[mat][1] = ixs[mat][daughters].sum(axis=1)
        ixs[mat] = ixs[mat][columns]
    tape = update_xs(tape, ixs)
    tape = write_mf3_mt(tape)
    write_tape(tape, "AAA-{}".format(ismp))