# -*- coding: utf-8 -*-
"""
Created on Thu Jul 12 09:47:01 2018

@author: Fiorito_L
"""

import pytest
import os
from random import randint

import numpy as np

from ...formats import Errorr
from ...formats import Endf6
from .. import sampling
from ...data import H1
from ...data import Cm242
from ...data import U5
from ...data import U8
from ...data import Fe56

@pytest.mark.sampling
@pytest.mark.xs
def test_sample_xs():
    errtape = Errorr.from_text("\n".join(H1.errorr))
    nsmp = 1000
    perts = errtape.get_xs_cov(listmt=[102]).get_samples(nsmp, eig=10)
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
    assert newtape.loc[125,3,2].TEXT == pendftape.loc[125,3,2].TEXT

@pytest.mark.sampling
@pytest.mark.errorr
@pytest.mark.xs
def test_H1(tmpdir):
    iargs = [os.path.join(H1.__path__[0], r"h1.pendf"),
             "--errorr-cov", os.path.join(H1.__path__[0], r"h1.errorr"),
             "--outdir", str(tmpdir),
             "--processes", str(os.cpu_count()),
             "--eig", "10",
             "--samples", "100",]
    sampling(iargs)

@pytest.mark.sampling
@pytest.mark.endf6
@pytest.mark.nubar
def test_Cm242(tmpdir):
    iargs = [os.path.join(Cm242.__path__[0], r"cm242.endf"),
             "--endf6-cov", os.path.join(Cm242.__path__[0], r"cm242.endf"),
             "--outdir", str(tmpdir),
             "--processes", str(os.cpu_count()),
             "--eig", "10",
             "--samples", "100",]
    sampling(iargs)

@pytest.mark.sampling
@pytest.mark.errorr
@pytest.mark.xs
def test_Fe56_errorr(tmpdir):
    iargs = [os.path.join(Fe56.__path__[0], r"fe56.pendf"),
             "--errorr-cov", os.path.join(Fe56.__path__[0], r"fe56.errorr"),
             "--outdir", str(tmpdir),
             "--processes", str(os.cpu_count()),
             "--eig", "10",
             "--samples", "10",]
    sampling(iargs)

@pytest.mark.sampling
@pytest.mark.errorr
@pytest.mark.xs
@pytest.mark.slow
def test_U5_errorr(tmpdir):
    iargs = [os.path.join(U5.__path__[0], r"u235.pendf"),
             "--errorr-cov", os.path.join(U5.__path__[0], r"u235.errorr"),
             "--outdir", str(tmpdir),
             "--processes", str(os.cpu_count()) if os.cpu_count() < 10 else str(10),
             "--eig", "10",
             "--samples", "100",]
    sampling(iargs)

@pytest.mark.sampling
@pytest.mark.chi
def test_U5_chi(tmpdir):
    iargs = [os.path.join(U5.__path__[0], r"u235.endf"),
             "--endf6-cov", os.path.join(U5.__path__[0], r"u235.endf"),
             "--outdir", str(tmpdir),
             "--processes", str(os.cpu_count()),
             "--eig", "10",
             "--samples", "100",]
    sampling(iargs)

@pytest.mark.sampling
@pytest.mark.lpc
def test_Fe56_lpc(tmpdir):
    iargs = [os.path.join(Fe56.__path__[0], r"fe56.endf"),
             "--endf6-cov", os.path.join(Fe56.__path__[0], r"fe56.endf"),
             "--outdir", str(tmpdir),
             "--processes", str(os.cpu_count()),
             "--eig", "10",
             "--samples", "100",]
    sampling(iargs)

@pytest.mark.sampling
@pytest.mark.lpc
def test_U238_lpc(tmpdir):
    iargs = [os.path.join(U8.__path__[0], r"u238.endf"),
             "--endf6-cov", os.path.join(U8.__path__[0], r"u238.endf"),
             "--outdir", str(tmpdir),
             "-e", "1e-5", "1e-4", "1e-3", "1e-2", "1e-1", "1e0", "1e1", "1e2", "1e3", "1e4", "1e5",
             "--verbose",
             "--processes", "1",
             "--eig", "10",
             "--mf", "34",
             "--samples", "10",]
    sampling(iargs)