# -*- coding: utf-8 -*-
"""
Created on Fri May 11 15:08:25 2018

@author: fiorito_l
"""
from sandy.formats.records import read_cont, read_list, read_float
from os.path import join, dirname, realpath
import pandas as pd
from sandy.formats.endf6 import split_endf, XsCov, Xs
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
        from sandy.formats.endf6 import Endf6
        return cls(Endf6.from_file(file))

    @classmethod
    def from_text(cls, text):
        from sandy.formats.endf6 import Endf6
        return cls(Endf6.from_text(text))

    def process(self, keep_mf=None, keep_mt=None):
        """
        Parse TEXT column.
        """
        tape = self.copy()
        tape['DATA'] = tape['TEXT'].apply(process_errorr_section, keep_mf=keep_mf, keep_mt=keep_mt)
        return Errorr(tape)

    def get_cov(self):
        """
        Extract xs covariances from errorr file into XsCov instance.
        """
        mat = self.index.get_level_values("MAT")[0]
        eg = self.loc[mat,1,451].DATA["EG"]
        List = []
        for x in self.query('MF==33 | MF==31').DATA:
            for mt1,y in x["RP"].items():
                List.append([mat, x["MT"], mat, mt1, y])
        frame = pd.DataFrame(List, columns=('MAT', 'MT','MAT1', 'MT1', 'COV'))
        MI = [(mat,mt,e) for mat,mt in sorted(set(zip(frame.MAT, frame.MT))) for e in eg]
        index = pd.MultiIndex.from_tuples(MI, names=("MAT", "MT", "E"))
        # initialize union matrix
        matrix = np.zeros((len(index),len(index)))
        for i,row in frame.iterrows():
            ix = index.get_loc((row.MAT,row.MT))
            ix1 = index.get_loc((row.MAT1,row.MT1))
            matrix[ix.start:ix.stop-1,ix1.start:ix1.stop-1] = row.COV
        i_lower = np.tril_indices(len(index), -1)
        matrix[i_lower] = matrix.T[i_lower]  # make the matrix symmetric
        return XsCov(matrix, index=index, columns=index)

    def get_xs(self):
        """
        Extract xs from errorr file into Xs instance.
        """
        mat = self.index.get_level_values("MAT")[0]
        eg = self.loc[mat,1,451].DATA["EG"]
        XsDict = dict(map(lambda x: ((x["MAT"],x["MT"]), x["XS"]), self.query("MF==3").DATA))
        frame = pd.DataFrame.from_dict(XsDict)
        frame.index = eg[:-1]
        frame = frame.reindex(eg, method='ffill')
        return Xs(frame)

    def get_std(self):
        """
        Extract xs and std from errorr file into dataframe.
        """
        xs = self.get_xs()
        cov = self.get_cov()
        stdvals = np.sqrt(np.diag(cov.values))
        xsvals =  xs.values.T.flatten()
        frame = pd.DataFrame.from_dict({"XS" : xsvals, "STD" : stdvals})
        frame.columns.name = "DATA"
        frame.index = cov.index
        frame = frame.unstack(level=["MAT","MT"])
        frame.columns = frame.columns.swaplevel(i=0, j=2).swaplevel(i=0, j=1)
        return frame



def read_mf1_mt451(text):
    str_list = split_endf(text)
    i = 0
    out = {"MAT" : str_list["MAT"].iloc[0],
           "MF" : str_list["MF"].iloc[0],
           "MT" : str_list["MT"].iloc[0]}
    C, i = read_cont(str_list, i)
    out.update({"ZA" : C.C1, "AWR" : C.C2, "ERRFLAG" :C.N1})
    L, i = read_list(str_list, i)
    out.update({"EG" : L.B})
    return out

def read_mf3_mt(text):
    str_list = split_endf(text)
    i = 0
    out = {"MAT" : str_list["MAT"].iloc[0],
           "MF" : str_list["MF"].iloc[0],
           "MT" : str_list["MT"].iloc[0]}
    L, i = read_list(str_list, i)
    out.update({"XS" : L.B})
    return out

def read_mf33_mt(text):
    str_list = split_endf(text)
    i = 0
    out = {"MAT" : str_list["MAT"].iloc[0],
           "MF" : str_list["MF"].iloc[0],
           "MT" : str_list["MT"].iloc[0]}
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





##############
# UNIT TESTS #
##############

#from sandy.data_test import __path__ as td
#A = Errorr.from_file(join(td[0], r"fe56.errorr")).process()
#xs = A.get_std()
#aaa=1

