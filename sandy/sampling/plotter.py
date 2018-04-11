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
import matplotlib.pyplot as plt

def run(path, mat, mf, mt, orig=None, cov=None):
    if orig is not None:
        print ("read best estimates from", orig)
        BE = e6.extract_xs(e6.endf2df(orig, keep_mf=[mf], keep_mt=[mt]))[mat,mt]
    if cov is not None:
        print ("read covariances from", cov)
        DfCov = e6.extract_cov33(e6.endf2df(cov, keep_mf=[mf+30], keep_mt=[456]))
        DfCov.index = DfCov.index.droplevel("MAT").droplevel("MT")
        DfCov.columns = DfCov.columns.droplevel("MAT").droplevel("MT")
    Std = pd.Series(np.sqrt(np.diag(DfCov))*100., index=DfCov.index)
    key = "NUBAR" if mt in (452,455,456) else "XS"
    ListXs = []
    for inp in os.listdir(path):
        print ("read file ", inp)
        D = e6.endf2df(join(path,inp), keep_mf=[mf], keep_mt=[mt]).DATA.loc[mat,mf,mt][key]
        D.name = inp
        ListXs.append(D)
    A = pd.DataFrame(ListXs)
    SmpMean = A.mean()
    SmpStd = A.std()/SmpMean*100.
    StdDiff = (SmpMean/BE-1)*100.

    fig, ax = plt.subplots()
    Std.plot(logx=True, drawstyle="steps", label="data", legend=True, ax=ax)
    SmpStd.plot(logx=True, linestyle="None", marker='o', markerfacecolor='None', label="samples", legend=True, ax=ax)
    ax.set_xlabel("E (eV)")
    ax.set_ylabel("stdev (%)")
    ax.legend(loc='best', fancybox=True, shadow=True)
    plt.show()
#    (Std/SS-1).plot(logx=True, drawstyle="steps")
#    plt.show()

    fig, ax = plt.subplots()
    BE.plot(logx=True, label="data", legend=True, ax=ax)
    SmpMean.plot(logx=True, label="samples", legend=True, ax=ax)
    ax.set_xlabel("E (eV)")
    ax.set_ylabel("NUBAR (-)")
    ax.legend(loc='best', fancybox=True, shadow=True)
#    plt.show()
    ax2 = ax.twinx()
    (SmpMean/BE-1).plot(logx=True, linestyle="None", marker='o', markerfacecolor='None', ax=ax2, alpha=0.8)
    ax2.set_ylabel('DIFF (%)')
    plt.show()
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
    run(join(sd, r"sampling//pu9-tmpdir"), 9437, 1, 452, orig=join(td,r"pu239.endf"), cov=join(td,r"pu239.endf"))