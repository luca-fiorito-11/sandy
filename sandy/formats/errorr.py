# -*- coding: utf-8 -*-
"""
Created on Fri May 11 15:08:25 2018

@author: fiorito_l
"""
from sandy.formats.records import read_cont, read_list, read_float
from os.path import join, dirname, realpath
import pandas as pd
from sandy.formats.endf6 import split
import numpy as np

def process_errorr_section(text, keep_mf=None, keep_mt=None):
    mf = int(text[70:72])
    mt = int(text[72:75])
    if mf == 1 and mt == 451: # read always
        return read_mf1_mt451(text)
    if keep_mf:
        if mf not in keep_mf:
            return None
    if keep_mt:
        if mt not in keep_mt:
            return None
    elif mf == 3:
        return read_mf3_mt(text)
    elif mf == 33:
        return read_mf33_mt(text)
    else:
        return None

class Errorr(pd.DataFrame):

    @classmethod
    def from_file(cls, file):
        """
        Read ENDF-6 formatted file and split it into MAT/MF/MT sections.
        Produce Endf6 instance (pandas.DataFrame) with columns
            MAT MF MT TEXT
        """
        columns = ('MAT', 'MF', 'MT','TEXT')
        rows = []
        for x in split(file):
            mat = int(x[66:70])
            mf = int(x[70:72])
            mt = int(x[72:75])
            text = "\n".join([ y for y in x.split("\n") ])
            rows.append([mat, mf, mt, text])
        frame = pd.DataFrame(rows, columns=columns)
        frame["DATA"] = None
        return cls(frame)

    def process(self, keep_mf=None, keep_mt=None):
        """
        Parse TEXT column.
        """
        tape = self.copy()
        tape['DATA'] = tape['TEXT'].apply(process_errorr_section, keep_mf=keep_mf, keep_mt=keep_mt)
        return Errorr(tape)



def read_mf1_mt451(text):
    str_list = text.splitlines()
    i = 0
    out = {"MAT" : int(str_list[i][66:70]),
           "MF" : int(str_list[i][70:72]),
           "MT" : int(str_list[i][72:75])}
    C, i = read_cont(str_list, i)
    out.update({"ZA" : C.C1, "AWR" : C.C2, "ERRFLAG" :C.N1})
    L, i = read_list(str_list, i)
    out.update({"EG" : L.B})
    return out

def read_mf3_mt(text):
    str_list = text.splitlines()
    i = 0
    out = {"MAT" : int(str_list[i][66:70]),
           "MF" : int(str_list[i][70:72]),
           "MT" : int(str_list[i][72:75])}
    L, i = read_list(str_list, i)
    out.update({"XS" : L.B})
    return out

def read_mf33_mt(text):
    str_list = text.splitlines()
    i = 0
    out = {"MAT" : int(str_list[i][66:70]),
           "MF" : int(str_list[i][70:72]),
           "MT" : int(str_list[i][72:75])}
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



class ErrorrCov(pd.DataFrame):

    def from_tape(tape):
        from functools import reduce
        mat = tape.index.get_level_values("MAT")[0]
        eg = tape[mat,1,451]["EG"]
        columns = ('MAT', 'MT', 'MAT1', 'MT1', 'COV')
        DfCov = pd.DataFrame(columns=columns)
        frame = pd.DataFrame([[x["MT"], mt1, y["RP"][mt1]] for x in tape.query('MF==33 | MF==31').DATA for mt1,y in x["RP"].items()])
            if not chunk: continue
            mt = chunk["MT"]
            for mt1, sub in chunk["RP"].items():
                cov = sub["RP"][mt1]
                # add zero row and column at the end of the matrix
                cov = np.insert(cov, cov.shape[0], [0]*cov.shape[1], axis=0)
                cov = np.insert(cov, cov.shape[1], [0]*cov.shape[0], axis=1)
                e1 = nisec["Ek"]
                e2 = nisec["El"]
                    else:
                        warn("skipped NI-type covariance with flag LB={} (MAT{}/MF{}/MT{})".fomat(nisec["LB"], chunk['MAT'], chunk['MF'], chunk['MT']), category=Warning)
                        continue
                    cov = pd.DataFrame(cov, index=e1, columns=e2)
                    covs.append(cov)
                if len(covs) == 0:
                    continue
                if len(covs) > 1:
                    import pdb
                    pdb.set_trace()
                # All union covariance matrices have the same grid (uxx) on both axis
                uxx = sorted(set().union(*[ [*list(x.index), *list(x.columns)] for x in covs ]))
                covs = list(map( lambda x: pandas_interpolate(x, uxx), covs ))
                cov = reduce(lambda x, y: x.add(y, fill_value=0), covs)
    #            index = pd.MultiIndex.from_product([[chunk["MAT"]],
    #                                                [chunk["MT"]],
    #                                                [sub["MAT1"]],
    #                                                [sub["MT1"]],
    #                                                cov.index], names=["MAT", "MT", "MAT1", "MT1", "E"])
                objects = chunk['MAT'], chunk['MT'], sub['MAT1'], sub['MT1'], cov
                DfCov = DfCov.append(dict(zip(columns, objects)), ignore_index=True)
        DfCov = DfCov.set_index(['MAT','MT','MAT1','MT1']).sort_index() # Multi-indexing
        if DfCov.empty:
            warn("no MF[31,33] covariances found", category=Warning)
            return pd.DataFrame()
        # Create big cov
        query_diag = 'MT==MT1 & (MAT1==0 | MAT==MAT1)'
        query_offdiag = 'MT!=MT1 | MAT1!=0'
    #    import scipy as sp
    #    C = sp.linalg.block_diag(*map(lambda x:x.as_matrix(), DfCov.query(query_diag).COV))
    #    idxs = []
    #    for (mat,mt,mat1,mt1), row in DfCov.query(query_diag).iterrows():
    #        idxs.append( pd.MultiIndex.from_product([[mat],[mt],[mat1],[mt1], row.COV.index]) )
    #    index = reduce(lambda x, y: x.append(y), idxs)
    #    C = pd.DataFrame(C, index=index, columns=index)
        diags = DfCov.query(query_diag)
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
        for (mat,mt,mat1,mt1),row in DfCov.query(query_offdiag).iterrows():
            cov = row.COV
            # interpolate x axis (rows)
            try:
                covk = diags.loc[mat,mt,0,mt]
            except:
                warn("cannot find covariance for MAT{}/MT{}".format(mat,mt), category=Warning)
                continue
            Ek = list(covk.COV.index)
            cov = pandas_interpolate(cov, Ek, method='zero', axis='rows')
            # interpolate y axis (cols)
            if mat1 == 0:
                mat1 = mat
            try:
                covl = diags.loc[mat1,mt1,0,mt1]
            except:
                warn("cannot find covariance for MAT{}/MT{}".format(mat1,mt1), category=Warning)
                continue
            El = list(covl.COV.index)
            cov = pandas_interpolate(cov, El, method='zero', axis='cols')
            C[covk.BEG:covk.END,covl.BEG:covl.END] = cov
            C[covl.BEG:covl.END,covk.BEG:covk.END,] = cov.T
        C = pd.DataFrame(C, index=[MATS,MTS,E], columns=[MATS,MTS,E])
        C.index.names = ["MAT", "MT", "E"]
        C.columns.names = ["MAT", "MT", "E"]
        return C

##############
# UNIT TESTS #
##############

def test_read_errorr():
    from sandy.data_test import __file__ as td
    td = dirname(realpath(td))
    Errorr.from_file(join(td, r"fe56.errorr")).process()