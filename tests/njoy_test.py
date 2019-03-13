# -*- coding: utf-8 -*-
"""
Created on Tue Mar 12 09:33:25 2019

@author: Luca Fiorito
"""

import pytest
import os

import sandy

__author__ = "Luca Fiorito"


@pytest.mark.njoy
def test_njoy_process_dryrun():
    """Test default options for njoy.process"""
    endftape = os.path.join(os.path.dirname(__file__), "data", "n-002_He_003.endf")
    input = sandy.njoy.process(endftape, dryrun=True)
    text = """moder
20 -21 /
reconr
-21 -22 /
'sandy runs njoy'/
225 0 0 /
0.001 0. /
0/
broadr
-21 -22 -23 /
225 1 0 0 0. /
0.001 /
293.6 /
0 /
thermr
0 -23 -24 /
0 225 20 1 1 0 0 1 221 0 /
293.6 /
0.001 10 /
heatr
-21 -24 -25 0 /
225 7 0 0 0 0 /
302 303 304 318 402 442 443 /
heatr
-21 -25 -26 0 /
225 4 0 0 0 0 /
444 445 446 447 /
gaspr
-21 -26 -27 /
purr
-21 -27 -28 /
225 1 1 20 32 /
293.6 /
1.00E+10 /
0 /
moder
-28 29 /
acer
-21 -28 0 50 70 /
1 0 1 .00 0 /
'sandy runs acer'/
225 293.6 /
1 1 /
/
stop"""
    assert input == text

@pytest.mark.njoy
def test_njoy_process_no_broadr():
    """Test njoy.process without broadr"""
    endftape = os.path.join(os.path.dirname(__file__), "data", "n-002_He_003.endf")
    input = sandy.njoy.process(endftape, dryrun=True, broadr=False)
    text = """moder
20 -21 /
reconr
-21 -22 /
'sandy runs njoy'/
225 0 0 /
0.001 0. /
0/
thermr
0 -22 -23 /
0 225 20 1 1 0 0 1 221 0 /
293.6 /
0.001 10 /
heatr
-21 -23 -24 0 /
225 7 0 0 0 0 /
302 303 304 318 402 442 443 /
heatr
-21 -24 -25 0 /
225 4 0 0 0 0 /
444 445 446 447 /
gaspr
-21 -25 -26 /
purr
-21 -26 -27 /
225 1 1 20 32 /
293.6 /
1.00E+10 /
0 /
moder
-27 28 /
acer
-21 -27 0 50 70 /
1 0 1 .00 0 /
'sandy runs acer'/
225 293.6 /
1 1 /
/
stop"""
    assert input == text

@pytest.mark.njoy
def test_njoy_process_no_gaspr():
    """Test njoy.process without gaspr"""
    endftape = os.path.join(os.path.dirname(__file__), "data", "n-002_He_003.endf")
    input = sandy.njoy.process(endftape, dryrun=True, broadr=False, gaspr=False)
    text = """moder
20 -21 /
reconr
-21 -22 /
'sandy runs njoy'/
225 0 0 /
0.001 0. /
0/
thermr
0 -22 -23 /
0 225 20 1 1 0 0 1 221 0 /
293.6 /
0.001 10 /
heatr
-21 -23 -24 0 /
225 7 0 0 0 0 /
302 303 304 318 402 442 443 /
heatr
-21 -24 -25 0 /
225 4 0 0 0 0 /
444 445 446 447 /
purr
-21 -25 -26 /
225 1 1 20 32 /
293.6 /
1.00E+10 /
0 /
moder
-26 27 /
acer
-21 -26 0 50 70 /
1 0 1 .00 0 /
'sandy runs acer'/
225 293.6 /
1 1 /
/
stop"""
    assert input == text

@pytest.mark.njoy
def test_njoy_process_no_thermr():
    """Test njoy.process without thermr"""
    endftape = os.path.join(os.path.dirname(__file__), "data", "n-002_He_003.endf")
    input = sandy.njoy.process(endftape, dryrun=True, broadr=False, gaspr=False,
                               thermr=False)
    text = """moder
20 -21 /
reconr
-21 -22 /
'sandy runs njoy'/
225 0 0 /
0.001 0. /
0/
heatr
-21 -22 -23 0 /
225 7 0 0 0 0 /
302 303 304 318 402 442 443 /
heatr
-21 -23 -24 0 /
225 4 0 0 0 0 /
444 445 446 447 /
purr
-21 -24 -25 /
225 1 1 20 32 /
293.6 /
1.00E+10 /
0 /
moder
-25 26 /
acer
-21 -25 0 50 70 /
1 0 1 .00 0 /
'sandy runs acer'/
225 293.6 /
1 1 /
/
stop"""
    assert input == text

@pytest.mark.njoy
def test_njoy_process_no_acer():
    """Test njoy.process without acer"""
    endftape = os.path.join(os.path.dirname(__file__), "data", "n-002_He_003.endf")
    input = sandy.njoy.process(endftape, dryrun=True, broadr=False, gaspr=False,
                               thermr=False, acer=False)
    text = """moder
20 -21 /
reconr
-21 -22 /
'sandy runs njoy'/
225 0 0 /
0.001 0. /
0/
heatr
-21 -22 -23 0 /
225 7 0 0 0 0 /
302 303 304 318 402 442 443 /
heatr
-21 -23 -24 0 /
225 4 0 0 0 0 /
444 445 446 447 /
purr
-21 -24 -25 /
225 1 1 20 32 /
293.6 /
1.00E+10 /
0 /
moder
-25 26 /
stop"""
    assert input == text

@pytest.mark.njoy
def test_njoy_process_no_purr():
    """Test njoy.process without acer"""
    endftape = os.path.join(os.path.dirname(__file__), "data", "n-002_He_003.endf")
    input = sandy.njoy.process(endftape, dryrun=True, broadr=False, gaspr=False, 
                               thermr=False, acer=False, purr=False)
    text = """moder
20 -21 /
reconr
-21 -22 /
'sandy runs njoy'/
225 0 0 /
0.001 0. /
0/
heatr
-21 -22 -23 0 /
225 7 0 0 0 0 /
302 303 304 318 402 442 443 /
heatr
-21 -23 -24 0 /
225 4 0 0 0 0 /
444 445 446 447 /
moder
-24 25 /
stop"""
    assert input == text

@pytest.mark.njoy
def test_njoy_process_no_heatr():
    """Test njoy.process without heatr"""
    endftape = os.path.join(os.path.dirname(__file__), "data", "n-002_He_003.endf")
    input = sandy.njoy.process(endftape, dryrun=True, broadr=False, gaspr=False,
                               thermr=False, acer=False, purr=False, heatr=False)
    text = """moder
20 -21 /
reconr
-21 -22 /
'sandy runs njoy'/
225 0 0 /
0.001 0. /
0/
moder
-22 23 /
stop"""
    assert input == text

@pytest.mark.njoy
def test_njoy_process_no_keep_pendf():
    """Test njoy.process and do not keep pendf"""
    endftape = os.path.join(os.path.dirname(__file__), "data", "n-002_He_003.endf")
    input = sandy.njoy.process(endftape, dryrun=True, broadr=False, gaspr=False,
                               thermr=False, acer=False, purr=False, heatr=False, keep_pendf=False)
    text = """moder
20 -21 /
reconr
-21 -22 /
'sandy runs njoy'/
225 0 0 /
0.001 0. /
0/
stop"""
    assert input == text

@pytest.mark.njoy
def test_njoy_process_pendftape():
    """Test njoy.process using argument pendftape (skip reconr)"""
    endftape = os.path.join(os.path.dirname(__file__), "data", "n-002_He_003.endf")
    input = sandy.njoy.process(endftape, pendftape="pendf", dryrun=True, broadr=False, gaspr=False,
                               thermr=False, acer=False, purr=False, heatr=False, keep_pendf=False)
    text = """moder
20 -21 /
moder
99 -22 /
stop"""

@pytest.mark.njoy
def test_njoy_process_temperatures():
    """Test njoy.process for different temperatures"""
    endftape = os.path.join(os.path.dirname(__file__), "data", "n-002_He_003.endf")
    input = sandy.njoy.process(endftape, pendftape="pendf", dryrun=True, gaspr=False,
                               thermr=False, acer=False, purr=False, heatr=False, keep_pendf=False,
                               temperatures=[300, 600.0000, 900.001])
    text = """moder
20 -21 /
moder
99 -22 /
broadr
-21 -22 -23 /
225 3 0 0 0. /
0.001 /
300.0 600.0 900.0 /
0 /
stop"""
    assert input == text

@pytest.mark.njoy
def test_njoy_process_acer():
    """Test njoy.process for acer at different temperatures"""
    endftape = os.path.join(os.path.dirname(__file__), "data", "n-002_He_003.endf")
    input = sandy.njoy.process(endftape, pendftape="pendf", dryrun=True, broadr=False, gaspr=False,
                               thermr=False, acer=True, purr=False, heatr=False, keep_pendf=False,
                               temperatures=[300, 600.0000, 900.001])
    text = """moder
20 -21 /
moder
99 -22 /
acer
-21 -22 0 50 70 /
1 0 1 .00 0 /
'sandy runs acer'/
225 300.0 /
1 1 /
/
acer
-21 -22 0 51 71 /
1 0 1 .01 0 /
'sandy runs acer'/
225 600.0 /
1 1 /
/
acer
-21 -22 0 52 72 /
1 0 1 .02 0 /
'sandy runs acer'/
225 900.0 /
1 1 /
/
stop"""
    assert input == text

@pytest.mark.njoy
def test_njoy_process_suffixes():
    """Test njoy.process for acer at different temperatures"""
    endftape = os.path.join(os.path.dirname(__file__), "data", "n-002_He_003.endf")
    input = sandy.njoy.process(endftape, pendftape="pendf", dryrun=True, broadr=False, gaspr=False,
                               thermr=False, acer=True, purr=False, heatr=False, keep_pendf=False,
                               temperatures=[300, 600.0000, 900.001], suffixes=[3, 6, 9])
    text = """moder
20 -21 /
moder
99 -22 /
acer
-21 -22 0 50 70 /
1 0 1 .03 0 /
'sandy runs acer'/
225 300.0 /
1 1 /
/
acer
-21 -22 0 51 71 /
1 0 1 .06 0 /
'sandy runs acer'/
225 600.0 /
1 1 /
/
acer
-21 -22 0 52 72 /
1 0 1 .09 0 /
'sandy runs acer'/
225 900.0 /
1 1 /
/
stop"""
    assert input == text

@pytest.mark.njoy
def test_njoy_process_sig0():
    """Test njoy.process for different sig0"""
    endftape = os.path.join(os.path.dirname(__file__), "data", "n-002_He_003.endf")
    input = sandy.njoy.process(endftape, pendftape="pendf", dryrun=True, broadr=False, gaspr=False,
                               thermr=False, acer=False, purr=True, heatr=False, keep_pendf=False,
                               sig0=[1e10, 1E9, 100000000])
    text = """moder
20 -21 /
moder
99 -22 /
purr
-21 -22 -23 /
225 1 3 20 32 /
293.6 /
1.00E+10 1.00E+09 1.00E+08 /
0 /
stop"""
    assert input == text

@pytest.mark.njoy
@pytest.mark.njoy_exe
def test_njoy_process(tmpdir):
    """Test njoy.process for ENDF/B-VIII.0 He-3.
    Check that desired outputs are produced and that xsdir files are correctly updated.
    """
    endftape = os.path.join(os.path.dirname(__file__), "data", "n-002_He_003.endf")
    input = sandy.njoy.process(endftape, temperatures=[300, 600, 900], suffixes=[3, 6, 15], tag="_b71", wdir=str(tmpdir),
                               thermr=False)
    outputs = os.listdir(str(tmpdir))
    assert '2003_b71.pendf' in outputs
    for ace in ['2003_b71.03c', '2003_b71.06c', '2003_b71.15c']:
        assert ace in outputs
        assert "{}.xsd".format(ace) in outputs
        xsdargs = open(os.path.join(str(tmpdir), ace) + ".xsd").read().split()
        assert len(xsdargs) == 10
        assert xsdargs[0] == "2003{}".format(os.path.splitext(ace)[1])
        assert xsdargs[2] == os.path.join(str(tmpdir), ace)
        assert xsdargs[3] == "0"

@pytest.mark.njoy
@pytest.mark.njoy_exe
def test_njoy_process_1(tmpdir):
    """Test njoy.process for ENDF/B-VIII.0 He-3.
    Check that desired outputs are produced and that xsdir files are correctly updated.
    """
    endftape = os.path.join(os.path.dirname(__file__), "data", "n-002_He_003.endf")
    input = sandy.njoy.process(endftape, wdir=str(tmpdir), thermr=False, acer=False, keep_pendf=False)
    assert not os.listdir(str(tmpdir))