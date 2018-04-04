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

def run():
    t0 = time.time()
    if len(sys.argv) == 1:
        from sandy.data_test import __file__ as td
        from sandy import __file__ as sd
        from os.path import join
        sd = os.path.dirname(os.path.realpath(sd))
        td = os.path.dirname(os.path.realpath(td))
        sys.argv.extend([join(td, r"H1.txt"),
                         "--pendf", join(td, r"H1.txt.pendf"),
                         "--outdir", r"tmp-dir",
                         "--njoy", join(sd, r"njoy2012_50.exe"),
                         "--eig", "10",
                         "--samples", "5",
                         "--processes", "1",
                         "-e", "1e-5",
                         "-e", "5e-5",
                         "-e", "1e-4",
                         "-e", "5e-4",
                         "-e", "1e-3",
                         "-e", "5e-3",
                         "-e", "1e-2",
                         "-e", "5e-2",
                         "-e", "1e-1",
                         "-e", "5e-1",
                         "-e", "1e0",
                         "-e", "5e0",
                         "-e", "1e1",
                         "-e", "5e1",])
    settings.init_sampling()

    A = []
    for i in range(1,6):
        Xs = e6.endf2df(join(settings.args.outdir, r"H1.txt.pendf-{}".format(i))).DATA.loc[125,3,102]["XS"]
        Xs.name = i
        A.append(Xs)
    Xs = e6.endf2df(join(td, r"H1.txt.pendf")).DATA.loc[125,3,102]["XS"]
    Xs.name = 0
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
    run()