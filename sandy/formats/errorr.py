# -*- coding: utf-8 -*-
"""
Created on Fri May 11 15:08:25 2018

@author: fiorito_l
"""
import logging
from functools import reduce

import pandas as pd
import numpy as np

from .utils import BaseFile, Xs, XsCov
from ..settings import SandyError

__author__ = "Luca Fiorito"
__all__ = ["Errorr"]

class Errorr(BaseFile):

    Format = "errorr"

    def read_section(self, mat, mf, mt):
        """
        Parse MAT/MF/MT section
        """
        if mf == 1:
            from .MF1 import read_errorr as read
        elif mf == 3:
            from .MF3 import read_errorr as read
        elif mf == 33 or mf == 31 or mf == 35:
            from .MF33 import read_errorr as read
        else:
            raise SandyError("SANDY cannot parse section MAT{}/MF{}/MT{}".format(mat,mf,mt))
        if (mat,mf,mt) not in self.index:
            raise SandyError("section MAT{}/MF{}/MT{} is not in tape".format(mat,mf,mt))
        return read(self.loc[mat,mf,mt].TEXT)

    def get_xs(self, listmat=None, listmt=None):
        """
        Extract xs from errorr file into Xs instance.
        """
        condition = self.index.get_level_values("MF") == 3
        tape = self[condition]
        if listmat is not None:
            conditions = [tape.index.get_level_values("MAT") == x for x in listmat]
            condition = reduce(lambda x,y: np.logical_or(x, y), conditions)
            tape = tape[condition]
        if listmt is not None:
            conditions = [tape.index.get_level_values("MT") == x for x in listmt]
            condition = reduce(lambda x,y: np.logical_or(x, y), conditions)
            tape = tape[condition]
#        query = "MF==3"
#        if listmat is not None:
#            query_mats = " | ".join(["MAT=={}".format(x) for x in listmat])
#            query += " & ({})".format(query_mats)
#        if listmt is not None:
#            query_mts = " | ".join(["MT=={}".format(x) for x in listmt])
#            query += " & ({})".format(query_mts)
#        tape = self.query(query)
        mat = self.index.get_level_values("MAT")[0]
        eg = self.read_section(mat,1,451)["EG"]
        ListXs = []
        for ix,text in tape.TEXT.iteritems():
            mat,mf,mt = ix
            X = self.read_section(*ix)
            xs = pd.Series(X["XS"], index=eg[:-1], name=(X["MAT"],X["MT"])).rename_axis("E").to_frame()
            ListXs.append(xs)
        if not ListXs:
            logging.warn("requested cross sections were not found")
            return pd.DataFrame()
        # Use concat instead of merge because indexes are the same
        frame = pd.concat(ListXs, axis=1).reindex(eg, method="ffill")
        return Xs(frame)

    def get_xs_cov(self, listmat=None, listmt=None):
        """
        Extract xs covariances from errorr file into XsCov instance.
        """
        conditions = [self.index.get_level_values("MF") == x for x in [31, 33, 35]]
        condition = reduce(lambda x,y: np.logical_or(x, y), conditions)
        tape = self[condition]
        if listmat is not None:
            conditions = [tape.index.get_level_values("MAT") == x for x in listmat]
            condition = reduce(lambda x,y: np.logical_or(x, y), conditions)
            tape = tape[condition]
        if listmt is not None:
            conditions = [tape.index.get_level_values("MT") == x for x in listmt]
            condition = reduce(lambda x,y: np.logical_or(x, y), conditions)
            tape = tape[condition]        
#        query = "(MF==33 | MF==31)"
#        if listmat is not None:
#            query_mats = " | ".join(["MAT=={}".format(x) for x in listmat])
#            query += " & ({})".format(query_mats)
#        if listmt is not None:
#            query_mts = " | ".join(["MT=={}".format(x) for x in listmt])
#            query += " & ({})".format(query_mts)
#        tape = self.query(query)
        mat = self.index.get_level_values("MAT")[0]
        eg = self.read_section(mat,1,451)["EG"]
        List = []
        for ix,text in tape.TEXT.iteritems():
            mat,mf,mt = ix
            X = self.read_section(*ix)
            for mt1,y in X["RP"].items():
                List.append([mat, X["MT"], mat, mt1, y])
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

    def get_std(self):
        """
        Extract xs and std from errorr file into dataframe:
            index = energy
            columns = (MAT, MT, DATA)
        """
        xs = self.get_xs()
        cov = self.get_cov()
        stdvals = np.sqrt(np.diag(cov.values))
        xsvals =  xs.values.T.flatten()
        frame = pd.DataFrame.from_dict({"XS" : xsvals, "STD" : stdvals})
        frame.columns.name = "DATA"
        frame.index = cov.index
        frame = frame.unstack(level=["MAT","MT"])
        frame.columns = frame.columns.reorder_levels(["MAT","MT","DATA"])
        return frame

    def get_edistr_cov(self):
        return pd.DataFrame()

    def get_lpc_cov(self):
        return pd.DataFrame()