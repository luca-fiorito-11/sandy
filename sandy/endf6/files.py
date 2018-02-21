# -*- coding: utf-8 -*-
"""
Created on Mon Jan 16 18:03:13 2017

@author: lfiorito
"""
import sys
import logging
import numpy as np
from records import read_cont, read_tab1, read_list, read_text
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
    out.update({"XS" : pd.Series(T.y, index = T.x,
                                 name = (out["MAT"],out["MT"]))})
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
    xs_dict = {}
    for mat in np.unique(tape.index.get_level_values(0)):
        xsdf = pd.DataFrame(tape.loc[mat,3,1].DATA["XS"])
        xsdf.columns = xsdf.columns.droplevel() # keep only MT
        # No interpolation is done because xs are on unionized grid
        for chunk in tape.query('MF==3 & MT!=1').DATA:
            df = pd.DataFrame(chunk["XS"])
            df.columns = df.columns.droplevel() # keep only MT
            xsdf = xsdf.add(df, fill_value=0)
        xsdf.fillna(0, inplace=True)
        xs_dict.update( { mat : xsdf })
    return xs_dict

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



#columns = ('MAT', 'MT', 'MAT1', 'MT1', 'COV')
#tape = pd.DataFrame(columns=('MAT', 'MF', 'MT','SEC'))
#for chunk in split("H1.txt"):
#    tape = tape.append({
#        "MAT" : int(chunk[66:70]),
#        "MF" : int(chunk[70:72]),
#        "MT" : int(chunk[72:75]),
#        }, ignore_index=True)
file = "H1.txt"
#file = "26-Fe-56g.jeff33"
tape = pd.DataFrame([[int(x[66:70]), int(x[70:72]), int(x[72:75]), x, int(x[66:70])*100000+int(x[70:72])*1000+int(x[72:75])] for x in split(file)],
        columns=('MAT', 'MF', 'MT','TEXT', 'ID'))
tape = tape.set_index(['MAT','MF','MT']).sort_index() # Multi-indexing
tape['DATA'] = tape['TEXT'].apply(process_section)

df_cov_xs = merge_covs( extract_cov33(tape) )
dict_xs = extract_xs(tape)

from sandy.cov import Cov
NSMP = 100
df_samples_xs = pd.DataFrame( Cov(df_cov_xs.as_matrix()).sampling(NSMP), index = df_cov_xs.index )
unique_ids = list(set(zip(df_cov_xs.index.get_level_values(0), df_cov_xs.index.get_level_values(1))))
#idx = pd.IndexSlice
#df_samples.loc[idx[:, :, ['C1', 'C3']], idx[:, 'foo']]
for ismp in range(NSMP):
    ixs = dict_xs.copy()
    for (mat,mt), P in df_samples_xs.groupby(level=(0,1)):
        if mat not in ixs:
            continue
        if mt not in ixs[mat]:
            continue
        XS = ixs[mat][mt]
        P = pandas_interpolate(P.loc[mat,mt][ismp], list(XS.index), axis="rows")
        XS = XS.multiply(P)
        # Negative values are set to mean
        XS[XS < 0] = ixs[mat][mt]
        ixs[mat][mt] = XS
    pass

write_mf33_mt(P)