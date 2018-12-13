# -*- coding: utf-8 -*-
"""
Created on Thu Jul 12 09:47:01 2018

@author: Fiorito_L
"""

import pytest
import os
from random import randint

import numpy as np

from sandy.formats.endf6 import Endf6
from sandy.formats.errorr import Errorr
from sandy.formats.utils import XsCov
from sandy.sampling import sampling
from sandy.data import H1, Cm242, U5, U8, Fe56, FY

@pytest.mark.sampling
@pytest.mark.xs
def test_sample_xs():
    errtape = Errorr.from_text("\n".join(H1.errorr))
    nsmp = 1000
    perts = XsCov.from_errorr(errtape.filter_by(listmt=[102,451])).get_samples(nsmp, eig=10)
    pendftape = Endf6.from_text("\n".join(H1.pendf))
    xs = pendftape.get_xs()
    ismp = randint(1, nsmp);
    pert = perts[ismp]
    pxs = xs.perturb(pert)
    newtape = pendftape.update_xs(pxs)
    assert perts.shape[1] == nsmp
    mat = 125; mt = 102
    ratio = pxs/xs.values
    mask1 = np.in1d(ratio[mat,mt].index, pert[mat,mt].index)
    mask2 = np.in1d(pert[mat,mt].index, ratio[mat,mt].index)
    assert np.isclose(ratio[mat,mt].values[mask1], pert[mat,mt].values[mask2]).all()
    assert newtape.loc[125,3,102].TEXT != pendftape.loc[125,3,102].TEXT
#    assert newtape.loc[125,3,2].TEXT == pendftape.loc[125,3,2].TEXT

@pytest.mark.sampling
@pytest.mark.fy
def test_sample_fy():
    tape = Endf6.from_text("\n".join(FY.endf6))
    nsmp = 1000
    ismp = randint(1, nsmp);
    fyu5 = tape.get_fy(listenergy=[4e5]).filter_by("MAT", 9228)
    fyall = tape.get_fy()
    cov = fyu5.get_cov(mat=9228, mt=454, energy=4e5)
    perts = cov.get_samples(nsmp)
    pert = perts[ismp]
    fynew = fyall.perturb(pert)
    assert (fyall.query("MAT!=9228 & MT!=454 & E!=4e5") == fynew.query("MAT!=9228 & MT!=454 & E!=4e5")).all().all()
    fy = fynew.query("MAT==9228 & MT==454 & E==4e5")
    assert not (fyall.query("MAT==9228 & MT==454 & E==4e5") == fy).all().all()
    assert (fy.YI >= 0).all()
    assert (fy.YI <= fy.YI*2).all()

@pytest.mark.sampling
@pytest.mark.errorr
@pytest.mark.xs
def test_H1(tmpdir):
    iargs = [os.path.join(H1.__path__[0], r"h1.pendf"),
             "--cov", os.path.join(H1.__path__[0], r"h1.errorr"),
             "--outdir", str(tmpdir),
             "--processes", str(os.cpu_count()),
             "--eig", "10",
             "--samples", "10",]
    sampling(iargs)

@pytest.mark.sampling
@pytest.mark.endf6
@pytest.mark.nubar
def test_Cm242(tmpdir):
    iargs = [os.path.join(Cm242.__path__[0], r"cm242.endf"),
             "--cov", os.path.join(Cm242.__path__[0], r"cm242.endf"),
             "--outdir", str(tmpdir),
             "--processes", str(os.cpu_count()),
             "--eig", "10",
             "--samples", "10",]
    sampling(iargs)

@pytest.mark.sampling
@pytest.mark.errorr
@pytest.mark.xs
def test_Fe56_errorr(tmpdir):
    iargs = [os.path.join(Fe56.__path__[0], r"fe56.pendf"),
             "--cov", os.path.join(Fe56.__path__[0], r"fe56.errorr"),
             "--outdir", str(tmpdir),
             "--processes", str(os.cpu_count()),
             "--eig", "10",
             "--samples", "10",]
#             "--mt", "2", "102"]
    sampling(iargs)

@pytest.mark.sampling
@pytest.mark.errorr
@pytest.mark.xs
@pytest.mark.slow
def test_U5_errorr(tmpdir):
    iargs = [os.path.join(U5.__path__[0], r"u235.pendf"),
             "--cov", os.path.join(U5.__path__[0], r"u235.errorr"),
             "--outdir", str(tmpdir),
             "--processes", str(os.cpu_count()) if os.cpu_count() < 10 else str(10),
             "--eig", "10",
             "--samples", "10",]
    sampling(iargs)

@pytest.mark.sampling
@pytest.mark.chi
def test_U5_chi(tmpdir):
    iargs = [os.path.join(U5.__path__[0], r"u235.endf"),
             "--cov", os.path.join(U5.__path__[0], r"u235.endf"),
             "--outdir", str(tmpdir),
             "--processes", str(os.cpu_count()),
             "--samples", "10",
             "--mf", "35"]
    sampling(iargs)

@pytest.mark.sampling
@pytest.mark.lpc
def test_Fe56_lpc(tmpdir):
    iargs = [os.path.join(Fe56.__path__[0], r"fe56.endf"),
             "--cov", os.path.join(Fe56.__path__[0], r"fe56.endf"),
             "--outdir", str(tmpdir),
             "--processes", str(os.cpu_count()),
             "--samples", "10",
             "--mf", "34"]
    sampling(iargs)

@pytest.mark.sampling
@pytest.mark.lpc
def test_U238_lpc(tmpdir):
    iargs = [os.path.join(U8.__path__[0], r"u238.endf"),
             "--cov", os.path.join(U8.__path__[0], r"u238.endf"),
             "--outdir", str(tmpdir),
             "--verbose",
             "--processes", str(os.cpu_count()),
             "--eig", "15",
             "--samples", "10",
             "--mf", "34"]
    sampling(iargs)

@pytest.mark.sampling
@pytest.mark.fy
def test_jeff33_fy(tmpdir):
    iargs = [os.path.join(FY.__path__[0], r"FY.jeff33"),
             "--outdir", str(tmpdir),
             "--mat", "9228",
             "--processes", str(os.cpu_count()),
             "--samples", "10",
             "--fission-yields"]
    sampling(iargs)