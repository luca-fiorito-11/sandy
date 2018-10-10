# -*- coding: utf-8 -*-
"""
Created on Wed Jun 20 13:10:47 2018

@author: fiorito_l
"""
import pdb, os, pytest
import pandas as pd
import numpy as np

def maxw_xs_int(x0, x1, s1, x2, s2):
    """
    input for http://www.wolframalpha.com
    ```
    integrate (s1 + (x-x1)(s2-s1)/(x2-x1)) (x/x0^2) exp(-(x/x0)) dx
    ```
    """
    from math import exp
    c2 = (s2-s1)/(x2-x1)
    c1 = (s1-x1*c2)
    E1 = exp(-x1/x0)
    E2 = exp(-x2/x0)
    RE = exp((x1-x2)/x0)
    definit_x_exp = ((x1/x0*E1 - x2/x0*E2) + (E1 - E2))
    definit_x2_exp = E1*((x1**2 - x2**2*RE)/x0 + 2*(x1 - x2*RE) + 2*x0*(1 - RE))
    I = c1*definit_x_exp + c2*definit_x2_exp
    return I

def x2_exp(x, x0):
    """
    input for http://www.wolframalpha.com
    ```
    integrate x^2 exp(-(x/x0)) dx
    ```
    """
    from math import exp
    return x0 * (- exp(-x/x0)) * (x**2 + 2*x*x0 + 2*x0**2)

def x_exp(x, x0):
    """
    input for http://www.wolframalpha.com
    ```
    integrate x exp(-(x/x0)) dx
    ```
    x + x0 = x0 * (x/x0 + 1) ~ x0 * exp(x/x0)
    """
    from math import exp
    return - exp(-x/x0) * (x + x0) * x0

def maxw_int(x0, x1, x2):
    """
    input for http://www.wolframalpha.com
    ```
    integrate (x/x0^2) exp(-(x/x0)) dx
    ```
    """
    return (x_exp(x2, x0) - x_exp(x1, x0))/x0**2

def macs(iargs=None):
    from ..formats import Endf6, Errorr
    from .. import settings
    MACS = None; MACS_UNC = None
    init = settings.init_macs(iargs)
    if init.pendf is not None:
        xs = Endf6.from_file(init.pendf).get_xs(listmt=init.listmt)
        MACS = pd.concat([xs.macs(E0=kT) for kT in init.kT]).reset_index(drop=True)
    if init.errorr is not None:
        cov = Errorr.from_file(init.errorr).get_xs_cov(listmt=init.listmt)
        MACS_UNC = pd.concat([cov.macs(E0=kT) for kT in init.kT]).reset_index(drop=True)
    return MACS, MACS_UNC

if __name__ == "__main__":
    m, mu = macs()
    if m is not None: print(m)
    if mu is not None: print(mu)


@pytest.mark.macs
@pytest.mark.xs
def test_H1(tmpdir):
    from ..data import H1
    iargs = ["--pendf", os.path.join(H1.__path__[0], r"h1.pendf"),
             "--errorr", os.path.join(H1.__path__[0], r"h1.errorr"),
             "--kT", "2.53E-2", "30000",
             "-mt", "102"]
    mframe, muframe = macs(iargs)
    assert np.isclose(mframe.query("E0 == 30000").MACS[1], 0.02505175)
    assert np.isclose(muframe.STD, 0.0255297).all()