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


def check_xs(xs, xspert, perts, mat, mt, ismp):
    pert = perts[perts.MT == mt][str(ismp)]
    ratio = xspert[(mat,mt)]/xs[(mat,mt)]
    egrid = ratio.index.get_level_values("E")
    for j in range(len(pert)-1):
        emin = perts[perts.MT == mt].E.values[j]
        emax = perts[perts.MT == mt].E.values[j+1]
        assert np.allclose(ratio[(egrid < emax) & (egrid >= emin)], pert.values[j], rtol=1e-6)
    j = len(pert) - 1
    emin = perts[perts.MT == mt].E.values[j]
    assert np.allclose(ratio[(egrid >= emin)], pert.values[j], rtol=1e-6)

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
@pytest.mark.njoy_exe
def test_h1(tmpdir):
    """Sampling test to check the following:
    
    - sampling runs correctly with keywords:
        
        * `outdir`
        * `cov`
        * `debug`
        * `seed33`
        * `samples`
    
    - `outdir` directory was created
    - only 1 sample was produced for MAT 125
    - pertubed files exist and are named correctly
    - perturbed files do not contain covariance MF sections
    - NJOY ran RECONR
    - NJOY did not run RECONR
    - perturbed files are actual PENDF 
    - perturbation do not change (seed33 option)
    - only cross sections were perturbed, other data remain unchanged
    - perturbations are dump to file (debug option)
    - ratio between perturbed files and original pendf are equal to perturbations
    - redundant cross sections are correclty summed up
    """
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
    ftape, covtape, outs = sampling(iargs)
    # check that only 1 sample was produced for MAT 125
    assert sorted(outs.keys()) == [125]
    assert sorted(outs[125].keys()) == [1]
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
    assert tape.get_file_format() == "pendf"
    # Check that covariance file is not errorr
    assert covtape.get_file_format() != "errorr"
    # Check that non perturbed data did not change
    assert endftape.TEXT.loc[125,2,151] == tape.TEXT.loc[125,2,151]
    assert endftape.TEXT.loc[125,4,2] == tape.TEXT.loc[125,4,2]
    assert endftape.TEXT.loc[125,6,102] == tape.TEXT.loc[125,6,102]
    # Check that perts do not change (seed33 option)
    perts = pd.read_csv(os.path.join(str(tmpdir), 'perts_mf33.csv'))["1"].values
    assert perts.tolist() == [1.0148227202065168, 1.0353486552229474, 0.8912750106863625]
    # Check that ratio between perturbed files and original pendf are equal to perturbations
    xs1 = tape.get_xs()
    xs = ftape.get_xs()
    ratio = xs1[(125,102)]/xs[(125,102)]
    egrid = ratio.index.get_level_values("E")
    assert np.allclose(ratio[(egrid < 1e1) & (egrid >=1e-5)], perts[0], rtol=1e-6)
    assert np.allclose(ratio[(egrid < 1e7) & (egrid >=1e1)], perts[1], rtol=1e-6)
    assert np.allclose(ratio[(egrid >= 1e7)], perts[2], rtol=1e-6)
    # Check that ratio between perturbed files and original pendf are equal to perturbations
    # also for redundant cross sections, although with a lower tolerance
    assert np.allclose(xs1.reconstruct_sums(), xs1)
    assert xs1[(125,2)].equals(xs[(125,2)])
    ratio = (xs1[(125,1)]-xs[(125,2)])/xs[(125,102)]
    assert np.allclose(ratio[(egrid < 1e1) & (egrid >=1e-5)], perts[0], rtol=1e-3)
    assert np.allclose(ratio[(egrid < 1e7) & (egrid >=1e1)], perts[1], rtol=1e-3)
    assert np.allclose(ratio[(egrid >= 1e7)], perts[2], rtol=1e-3)

@pytest.mark.sampling
@pytest.mark.njoy_exe
def test_h1_run_errorr(tmpdir):
    """Sampling test to check the following:
    
    - sampling runs correctly with keywords:
        
        * `outdir`
        * errorr
        * `debug`
        * `seed33`
        * `samples`
    
    - `outdir` directory was created
    - only 2 sample were produced for MAT 125
    - pertubed files exist and are named correctly
    - NJOY ran RECONR and ERRORR
    - perturbations are dump to file (debug option)
    - ratio between perturbed files and original pendf are equal to perturbations
    - redundant cross sections are correclty summed up
    """
    name = "1-H-1g.jeff33"
    file = os.path.join(os.path.dirname(__file__), "data", name)
    endftape = sandy.read_formatted_file(file)
    iargs = [file,
             "--outdir", str(tmpdir),
             "--error",
             "--debug", 
             "--seed33", "555", 
             "--samples", "2",]
    ftape, covtape, outs = sampling(iargs)
    # Check that covariance file is  errorr
    assert covtape.get_file_format() == "errorr"
    # Check that PENDF flag was changed
    tape = sandy.read_formatted_file(os.path.join(str(tmpdir),"{}-1".format(name)))
    infosec = tape.read_section(125, 1, 451)
    assert infosec["TEMP"] == 0
    assert tape.get_file_format() == "pendf"
    # Check that ratio between perturbed files and original pendf are equal to perturbations
    xs = ftape.get_xs()
    perts = pd.read_csv(os.path.join(str(tmpdir), 'perts_mf33.csv'))
    for i in range(2):
        tape = sandy.read_formatted_file(os.path.join(str(tmpdir),"{}-{}".format(name, i+1)))
        xs1 = tape.get_xs()
        for mt in (2, 102):
            check_xs(xs, xs1, perts, 125, mt, i+1)
        assert np.allclose(xs1[(125,2)]+xs1[(125,102)], xs1[(125,1)], rtol=1e-6)

@pytest.mark.sampling
@pytest.mark.njoy_exe
def test_h1_cov_errorr(tmpdir):
    """Sampling test to check the following:
    
    - sampling runs correctly with keywords:
        
        * `outdir`
        * `cov`
        * `debug`
        * `seed33`
        * `samples`
    
    - SANDY reads covariances in ERRORR format
    - NJOY ran RECONR
    - perturbations are dump to file (debug option)
    - ratio between perturbed files and original pendf are equal to perturbations
    - redundant cross sections are correclty summed up
    """
    name = "1-H-1g.jeff33"
    file = os.path.join(os.path.dirname(__file__), "data", name)
    endftape = sandy.read_formatted_file(file)
    covname = "1-H-1g.jeff33.errorr"
    filecov = os.path.join(os.path.dirname(__file__), "data", covname)
    iargs = [file,
             "--outdir", str(tmpdir),
             "--cov", filecov, 
             "--debug", 
             "--seed33", "555", 
             "--samples", "2",]
    ftape, covtape, outs = sampling(iargs)
    # Check that covariance file is  errorr
    assert covtape.get_file_format() == "errorr"
    # Check that NJOY ran RECONR
    assert ftape.get_file_format() == "pendf"
    # Check that ratio between perturbed files and original pendf are equal to perturbations
    xs = ftape.get_xs()
    perts = pd.read_csv(os.path.join(str(tmpdir), 'perts_mf33.csv'))
    for i in range(2):
        tape = sandy.read_formatted_file(os.path.join(str(tmpdir),"{}-{}".format(name, i+1)))
        xs1 = tape.get_xs()
        for mt in (2, 102):
            check_xs(xs, xs1, perts, 125, mt, i+1)
        assert np.allclose(xs1[(125,2)]+xs1[(125,102)], xs1[(125,1)], rtol=1e-6)

@pytest.mark.sampling
def test_h1_pendf_cov_errorr(tmpdir):
    """Sampling test to check the following:
    
    - sampling runs correctly with keywords:
        
        * `outdir`
        * `cov`
        * `samples`
    
    - SANDY reads covariances in ERRORR format
    - SANDY reads input file in PENDF format
    - NJOY does not run
    - perturbations are dump to file (debug option)
    - ratio between perturbed files and original pendf are equal to perturbations
    - redundant cross sections are correclty summed up
    """
    name = "1-H-1g.jeff33.pendf"
    file = os.path.join(os.path.dirname(__file__), "data", name)
    endftape = sandy.read_formatted_file(file)
    covname = "1-H-1g.jeff33.errorr"
    filecov = os.path.join(os.path.dirname(__file__), "data", covname)
    iargs = [file,
             "--outdir", str(tmpdir),
             "--debug",
             "--cov", filecov, 
             "--samples", "2",]
    ftape, covtape, outs = sampling(iargs)
    # Check that covariance file is  errorr
    assert covtape.get_file_format() == "errorr"
    # Check that NJOY ran RECONR
    assert ftape.get_file_format() == "pendf"
    # Check that ratio between perturbed files and original pendf are equal to perturbations
    xs = ftape.get_xs()
    perts = pd.read_csv(os.path.join(str(tmpdir), 'perts_mf33.csv'))
    for i in range(2):
        tape = sandy.read_formatted_file(os.path.join(str(tmpdir),"{}-{}".format(name, i+1)))
        xs1 = tape.get_xs()
        for mt in (2, 102):
            check_xs(xs, xs1, perts, 125, mt, i+1)
        assert np.allclose(xs1[(125,2)]+xs1[(125,102)], xs1[(125,1)], rtol=1e-6)

@pytest.mark.sampling
def test_h1_pendf_cov_errorr_only_mt102(tmpdir):
    """Sampling test to check the following:
    
    - sampling runs correctly with keywords:
        
        * `outdir`
        * `cov`
        * `samples`
    
    - SANDY reads covariances in ERRORR format
    - SANDY reads input file in PENDF format
    - NJOY does not run
    - perturbations are dump to file (debug option)
    - ratio between perturbed files and original pendf are equal to perturbations
    - redundant cross sections are correclty summed up
    """
    name = "1-H-1g.jeff33.pendf"
    file = os.path.join(os.path.dirname(__file__), "data", name)
    endftape = sandy.read_formatted_file(file)
    covname = "1-H-1g.jeff33.errorr"
    filecov = os.path.join(os.path.dirname(__file__), "data", covname)
    iargs = [file,
             "--outdir", str(tmpdir),
             "--debug",
             "--cov", filecov, 
             "--mt", "102",
             "--samples", "2",]
    ftape, covtape, outs = sampling(iargs)
    # Check that ratio between perturbed files and original pendf are equal to perturbations
    xs = ftape.get_xs()
    perts = pd.read_csv(os.path.join(str(tmpdir), 'perts_mf33.csv'))
    for i in range(2):
        tape = sandy.read_formatted_file(os.path.join(str(tmpdir),"{}-{}".format(name, i+1)))
        # MT102 is perturbed
        xs1 = tape.get_xs()
        check_xs(xs, xs1, perts, 125, 102, i+1)
        # MT2 is not perturbed
        assert np.allclose((xs1/xs)[125,2], 1, rtol=1e-3)
        # MT1 is the sum
        assert np.allclose(xs1[(125,2)]+xs1[(125,102)], xs1[(125,1)], rtol=1e-6)