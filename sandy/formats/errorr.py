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
        frame = pd.DataFrame(rows, columns=columns).set_index(['MAT','MF','MT']).sort_index()
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

    @classmethod
    def from_tape(cls, tape):
        mat = tape.index.get_level_values("MAT")[0]
        eg = tape.loc[mat,1,451].DATA["EG"]
        List = []
        for x in tape.query('MF==33 | MF==31').DATA:
            for mt1,y in x["RP"].items():
                List.append([x["MT"], mt1, y[mt1]])
        frame = pd.DataFrame(List, columns=('MT', 'MT1', 'COV'))
        return

##############
# UNIT TESTS #
##############

def test_read_errorr():
    from sandy.data_test import __file__ as td
    td = dirname(realpath(td))
    tape = Errorr.from_file(join(td, r"fe56.errorr")).process()
    ErrorrCov.from_tape(tape)

test_read_errorr()