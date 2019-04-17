# -*- coding: utf-8 -*-
"""
Collection of functions to test parsing options for the sampling module.
"""

import pytest
import os

import sandy
from sandy.sampling import parse

__author__ = "Luca Fiorito"


@pytest.mark.sampling
def test_sampling_parse_version():
    """Check that parser correctly handles argument `version`."""
    with pytest.raises(SystemExit):
        parse(["--version"])
    with pytest.raises(SystemExit):
        parse(["-v"])

@pytest.mark.sampling
def test_sampling_parse_help():
    """Check that parser correctly handles argument `help`."""
    with pytest.raises(SystemExit):
        parse(["--help"])
    with pytest.raises(SystemExit):
        parse(["-h"])

@pytest.fixture(scope="module")
@pytest.mark.sampling
def file():
    return os.path.join(os.path.dirname(__file__), "data", "n-002_He_003.endf")

@pytest.mark.sampling
def test_sampling_parse_file(file):
    """Check that parser correctly handles argument `file`.
    Raise error if file does not exist."""
    init = parse([file])
    assert init.file == file
    with pytest.raises(SystemExit):
        init = parse(["abc"])

@pytest.mark.sampling
def test_sampling_parse_cov(file):
    """Check that parser correctly handles argument `cov`.
    Raise error if file does not exist."""
    init = parse([file, "--cov", file])
    assert init.cov == file
    init = parse([file, "-C", file])
    assert init.cov == file
    init = parse([file])
    assert init.cov == None
    with pytest.raises(SystemExit):
        init = parse([file, "--cov", "abc"])

@pytest.mark.sampling
def test_sampling_samples(file):
    """Check that parser correctly handles argument `samples`.
    Raise error if file does not exist."""
    init = parse([file, "--samples", "3"])
    assert init.samples == 3
    init = parse([file, "-S", "5"])
    assert init.samples == 5
    init = parse([file])
    assert init.samples == 200

@pytest.mark.sampling
def test_sampling_outdir(tmpdir, file):
    """Check that parser correctly handles argument `outdir`."""
    outdir = os.path.join(str(tmpdir), "AAA")
    assert not os.path.isdir(outdir)
    init = parse([file, "--outdir", outdir])
    assert init.outdir == outdir
    assert os.path.isdir(outdir)
    init = parse([file, "-D", outdir])
    assert init.outdir == outdir
    init = parse([file])
    assert init.outdir == os.getcwd()

@pytest.mark.sampling
def test_sampling_processes(file):
    """Check that parser correctly handles argument `processes`."""
    init = parse([file, "--processes", "4"])
    assert init.processes == 4
    init = parse([file, "-N", "4"])
    assert init.processes == 4
    init = parse([file])
    assert init.processes == 1

@pytest.mark.sampling
def test_sampling_seed31(file):
    """Check that parser correctly handles argument `seed31`."""
    init = parse([file, "--seed31", "4"])
    assert init.seed31 == 4
    init = parse([file])
    assert init.seed31 is None

@pytest.mark.sampling
def test_sampling_seed33(file):
    """Check that parser correctly handles argument `seed33`."""
    init = parse([file, "--seed33", "4"])
    assert init.seed33 == 4
    init = parse([file])
    assert init.seed33 is None

@pytest.mark.sampling
def test_sampling_seed34(file):
    """Check that parser correctly handles argument `seed34`."""
    init = parse([file, "--seed34", "4"])
    assert init.seed34 == 4
    init = parse([file])
    assert init.seed34 is None

@pytest.mark.sampling
def test_sampling_seed35(file):
    """Check that parser correctly handles argument `seed35`."""
    init = parse([file, "--seed35", "4"])
    assert init.seed35 == 4
    init = parse([file])
    assert init.seed35 is None

@pytest.mark.sampling
def test_sampling_mt(file):
    """Check that parser correctly handles argument `mt`."""
    init = parse([file, "--mt", "18"])
    assert init.mt == [18]
    init = parse([file, "--mt", "18", "20"])
    assert init.mt == [18, 20]
    init = parse([file])
    assert init.mt == range(1,1000)

@pytest.mark.sampling
def test_sampling_mat(file):
    """Check that parser correctly handles argument `mat`."""
    init = parse([file, "--mat", "18"])
    assert init.mat == [18]
    init = parse([file, "--mat", "18", "20"])
    assert init.mat == [18, 20]
    init = parse([file])
    assert init.mat == range(1,10000)

@pytest.mark.sampling
def test_sampling_mf(file):
    """Check that parser correctly handles argument `mf`."""
    init = parse([file, "--mf", "18"])
    assert init.mf == [18]
    init = parse([file, "--mf", "31", "32"])
    assert init.mf == [31, 32]
    init = parse([file])
    assert init.mf == [31, 33, 34, 35]
