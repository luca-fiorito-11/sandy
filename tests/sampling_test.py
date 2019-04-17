# -*- coding: utf-8 -*-
"""
Created on Thu Jul 12 09:47:01 2018

@author: Fiorito_L
"""

import pytest
import os
from random import randint
from io import StringIO
import numpy as np
import pandas as pd

import sandy
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
    assert np.allclose(ratio[mat,mt].values[mask1], pert[mat,mt].values[mask2])
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

@pytest.mark.sampling
def test_h1(tmpdir):
    name = "1-H-1g.jeff33"
    file = os.path.join(os.path.dirname(__file__), "data", name)
    endftape = sandy.read_formatted_file(file)
    covname = "h1cov.3x3"
    filecov = os.path.join(os.path.dirname(__file__), "data", covname)
    iargs = [file,
             "--cov", filecov,
             "--outdir", str(tmpdir),
             "--debug", 
             "--seed33", "555", 
             "--samples", "1",]
    sampling(iargs)
    # Check that perturbed files exist
    assert "{}-1".format(name) in os.listdir(str(tmpdir))
    tape = sandy.read_formatted_file(os.path.join(str(tmpdir),"{}-1".format(name)))
    # Check that covariances were removed
    assert tape.index.equals(pd.MultiIndex(levels=[[125], [1, 2, 3, 4, 6], [1, 2, 102, 151, 451]],
           codes=[[0, 0, 0, 0, 0, 0, 0], [0, 1, 2, 2, 2, 3, 4], [4, 3, 0, 1, 2, 1, 2]],
           names=['MAT', 'MF', 'MT']))
    infosec = tape.read_section(125, 1, 451)
    # Check that PENDF flag was changed
    assert infosec["TEMP"] == 0
    assert infosec["LRP"] == 2
    # Check that non perturbed data did not change
    assert endftape.TEXT.loc[125,2,151] == tape.TEXT.loc[125,2,151]
    assert endftape.TEXT.loc[125,4,2] == tape.TEXT.loc[125,4,2]
    assert endftape.TEXT.loc[125,6,102] == tape.TEXT.loc[125,6,102]
    # Check that ratio between perturbed files and original pendf are equal to perturbations
    pendf = sandy.read_formatted_file(os.path.join(str(tmpdir), "tape30"))
    xs1 = tape.get_xs()
    xs = pendf.get_xs()
    perts = pd.read_csv(os.path.join(str(tmpdir), 'perts_mf33.csv'))["1"].values
    ratio = xs1[(125,102)]/xs[(125,102)]
    egrid = ratio.index.get_level_values("E")
    assert np.allclose(ratio[(egrid < 1e1) & (egrid >=1e-5)], perts[0], rtol=1e-6)
    assert np.allclose(ratio[(egrid < 1e7) & (egrid >=1e1)], perts[1], rtol=1e-6)
    assert np.allclose(ratio[(egrid >= 1e7)], perts[2], rtol=1e-6)
    # Check that ratio between perturbed files and original pendf are equal to perturbations
    # also for redundant cross sections, although with a lower tolerance
    assert xs1[(125,2)].equals(xs[(125,2)])
    ratio = (xs1[(125,1)]-xs[(125,2)])/xs[(125,102)]
    assert np.allclose(ratio[(egrid < 1e1) & (egrid >=1e-5)], perts[0], rtol=1e-3)
    assert np.allclose(ratio[(egrid < 1e7) & (egrid >=1e1)], perts[1], rtol=1e-3)
    assert np.allclose(ratio[(egrid >= 1e7)], perts[2], rtol=1e-3)

@pytest.mark.sampling
def test_h1_errorr(tmpdir):
    name = "1-H-1g.jeff33"
    file = os.path.join(os.path.dirname(__file__), "data", name)
    endftape = sandy.read_formatted_file(file)
    iargs = [file,
             "--outdir", str(tmpdir),
             "--error",
             "--debug", 
             "--seed33", "555", 
             "--samples", "2",]
    sampling(iargs)
    # Check that perturbed files exist
    assert "{}-1".format(name) in os.listdir(str(tmpdir))
    assert "{}-2".format(name) in os.listdir(str(tmpdir))
    tape = sandy.read_formatted_file(os.path.join(str(tmpdir),"{}-1".format(name)))
    # Check that ratio between perturbed files and original pendf are equal to perturbations
    pendf = sandy.read_formatted_file(os.path.join(str(tmpdir), "tape30"))
    xs = pendf.get_xs()
    perts = pd.read_csv(os.path.join(str(tmpdir), 'perts_mf33.csv'))
    for i in range(2):
        tape = sandy.read_formatted_file(os.path.join(str(tmpdir),"{}-{}".format(name, i+1)))
        xs1 = tape.get_xs()
        for mt in (2, 102):
            pert = perts[perts.MT == mt][str(i+1)]
            ratio = xs1[(125,mt)]/xs[(125,mt)]
            egrid = ratio.index.get_level_values("E")#perts[perts.MT == mt].E
            for j in range(len(pert)-1):
                emin = perts[perts.MT == mt].E.values[j]
                emax = perts[perts.MT == mt].E.values[j+1]
                assert np.allclose(ratio[(egrid < emax) & (egrid >= emin)], pert.values[j], rtol=1e-6)
            j = len(pert) - 1
            emin = perts[perts.MT == mt].E.values[j]
            assert np.allclose(ratio[(egrid >= emin)], pert.values[j], rtol=1e-6)
        assert np.allclose(xs1[(125,2)]+xs1[(125,102)], xs1[(125,1)], rtol=1e-6)