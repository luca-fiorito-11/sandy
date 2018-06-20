# -*- coding: utf-8 -*-
"""
Created on Wed Jun 20 13:10:47 2018

@author: fiorito_l
"""
import pdb, os
import numpy as np
import pandas as pd

def numerator_integral(x0, x1, s1, x2, s2):
    """
    input for http://www.wolframalpha.com
    ```
    integrate (s1 + (x-x1)(s2-s1)/(x2-x1)) (x/x0^2) exp(-(x/x0)) dx
    ```
    """
    def foo(x):
        from math import exp
        C = exp(-x/x0) / x0 / (x1-x2)
        term1 = s1 * (-(x**2) + x*(x2-2*x0) + x0*(x2-2*x0))
        term2 = s2 * (x**2 + 2*x*x0 - x*x1 + 2*x0**2 - x0*x1)
        return C*(term1 + term2)
    return foo(x2) - foo(x1)

def denominator_integral(x0, x1, x2):
    """
    input for http://www.wolframalpha.com
    ```
    integrate (x/x0^2) exp(-(x/x0)) dx
    ```
    """
    def foo(x):
        from math import exp
        return - exp(-x/x0) * (x + x0) / x0
    return foo(x2) - foo(x1)

def group_integral(x0, x1, s1, x2, s2):
    from math import sqrt, pi
    return (2/sqrt(pi)) * numerator_integral(x0, x1, s1, x2, s2) / denominator_integral(x0, x1, x2)

def calculate_integral(xs, E0=0.0253, Elo=1E-5, Ehi=1E1, listmat=None, listmt=None):
    from math import sqrt, pi
    index = set(xs.index.values)
    index.update([Elo, Ehi])
    index = np.array(sorted(index))
    index = index[(index >= Elo) & (index <= Ehi)]
    xs = xs.reindex(index).interpolate(method='slinear', axis=0).fillna(0)
    D = np.sum([denominator_integral(E0, xs.index[i], xs.index[i+1]) for i in range(len(xs)-1)])
    records = []
    for (mat,mt),x in xs.items():
        N = np.sum([numerator_integral(E0, x.index[i], x.iloc[i], x.index[i+1], x.iloc[i+1]) for i in range(len(x)-1)])
        I = N/D*(2/sqrt(pi))
        records.append([mat, mt, I, D, Elo, Ehi, E0])
    return pd.DataFrame.from_records(records, columns=["MAT", "MT", "MACS", "FLUX", "Elo (eV)", "Ehi (eV)", "E0 (eV)"])

def calculate_uncertainty(cov, E0=0.0253, Elo=1E-5, Ehi=1E1, listmat=None, listmt=None):
    records = []
    for (mat,mt),sec in cov.groupby(["MAT","MT"]):
        C = sec[mat,mt].loc[mat,mt]
        E = set(C.index.values)
        E.update([Elo, Ehi])
        E = np.array(sorted(E))
        E = E[(E >= Elo) & (E <= Ehi)]
        C = C.reindex(E).ffill().fillna(0).T.reindex(E).ffill().fillna(0)
        D = np.array([denominator_integral(E0, E[i], E[i+1]) for i in range(len(E)-1)])
        S = D/np.sum(D)
        rvar = S.dot(C.values[:-1,:-1].dot(S))
        rstd = np.sqrt(rvar)
        records.append([mat, mt, rvar, rstd])
    return pd.DataFrame.from_records(records, columns=["MAT", "MT", "VAR", "STD"])

def run(iargs=None):
    from ..formats import Endf6, Errorr
    from .. import settings
    init = settings.init_ri(iargs)
    if init.pendf is not None:
        xs = Endf6.from_file(init.pendf).get_xs()
        frame = calculate_integral(xs)
        print(frame)
    elif init.errorr is not None:
        cov = Errorr.from_file(init.errorr).get_xs_cov()
        frame = calculate_uncertainty(cov)
        print(frame)

if __name__ == "__main__":
    run()
