# -*- coding: utf-8 -*-
"""
Created on Wed Apr  4 21:43:21 2018

@author: lucaf
"""
import pandas as pd
from sandy.endf6 import files as e6
from sandy import settings
from sandy.sampling.cov import Cov
import numpy as np
import sys
import os
import multiprocessing as mp
from copy import deepcopy
import shutil
import time
from os.path import join, dirname, realpath
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

def run(path, mat, mf, mt, orig, filecov=None):
    key = "NUBAR" if mt in (452,455,456) else "XS"
    A = []
    for inp in os.listdir(path):
        Xs = e6.endf2df(join(path,inp), keep_mf=[mf], keep_mt=[mt]).DATA.loc[mat,mf,mt][key]
        print ("read file ", inp)
        Xs.name = inp
        A.append(Xs)
    Xs = e6.endf2df(orig, keep_mf=[mf], keep_mt=[mt]).DATA.loc[mat,mf,mt][key]
    Xs.name = orig
    A.append(Xs)
    A = pd.DataFrame(A).T
    for i in range(1,6):
        A[i] /= A[0]
    ListPerts = []
    rowold = pd.Series()
    for e,row in A.iterrows():
        row.drop(0, inplace=True)
        if rowold.empty:
            e0 = deepcopy(e)
        elif not np.allclose(row, rowold):
            ListPerts.append(( e0, e, *rowold.tolist() ))
            e0 = deepcopy(e)
        rowold = deepcopy(row)
    B = pd.DataFrame.from_records(ListPerts)
    B.set_index([0, 1], inplace=True)
    C = B.as_matrix()
    C = np.insert(C, C.shape[0], C[-1]*C.shape[1], axis=0)
    C = np.insert(C, C.shape[1], C[:,-1]*C.shape[0], axis=1)
    C = np.corrcoef(C)
    x = B.index.get_level_values(0).tolist() + [B.index.get_level_values(1).tolist()[-1]]
    e6.plot_heatmap(x, x, C,
                    xscale="log", yscale="log",
                    vmin=-1, vmax=1,
                    cmap="bwr",
                    xlabel=None, ylabel=None, title=None)

if __name__ == '__main__':
    from sandy import __file__ as sd
    from sandy.data_test import __file__ as td
    sd = dirname(realpath(sd))
    td = dirname(realpath(td))
    run(join(sd, r"..\U5t"), 9228, 1, 456, join(td, r"92-U-235g.jeff33"))