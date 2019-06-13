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
def test_get_njoy_from_environ():
    exeold = None
    if "NJOY" in os.environ:
        exeold = os.environ["NJOY"]
        del os.environ["NJOY"]
    os.environ["NJOY"] = "/path/to/my_njoy.exe"
    exe = sandy.get_njoy()
    assert exe == "/path/to/my_njoy.exe"
    del os.environ["NJOY"]
    if exeold:
        os.environ["NJOY"] = exeold

@pytest.mark.njoy
def test_get_njoy_from_environ_error():
    exe = None
    if "NJOY" in os.environ:
        exe = os.environ["NJOY"]
        del os.environ["NJOY"]
    with pytest.raises(Exception):
        sandy.get_njoy()
    if exe:
        os.environ["NJOY"] = exe
    
@pytest.mark.njoy
def test_njoy_process_dryrun():
    """Test default options for njoy.process"""
    endftape = os.path.join(os.path.dirname(__file__), "data", "n-002_He_003.endf")
    input, inputs, outputs = sandy.njoy.process(endftape, dryrun=True)
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
225 1 1 20 32 0 /
293.6 /
1.00E+10 /
0 /
moder
-28 30 /
acer
-21 -28 0 50 70 /
1 0 1 .02 0 /
'sandy runs acer'/
225 293.6 /
1 1 /
/
stop"""
    assert input == text
    assert inputs['tape20'] == endftape
    assert outputs['tape30'] == '2003.pendf'
    assert outputs['tape50'] == '2003.02c'
    assert outputs['tape70'] == '2003.02c.xsd'

@pytest.mark.njoy
def test_njoy_process_no_broadr():
    """Test njoy.process without broadr"""
    endftape = os.path.join(os.path.dirname(__file__), "data", "n-002_He_003.endf")
    input, inputs, outputs = sandy.njoy.process(endftape, dryrun=True, broadr=False)
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
225 1 1 20 32 0 /
293.6 /
1.00E+10 /
0 /
moder
-27 30 /
acer
-21 -27 0 50 70 /
1 0 1 .02 0 /
'sandy runs acer'/
225 293.6 /
1 1 /
/
stop"""
    assert input == text
    assert inputs['tape20'] == endftape
    assert outputs['tape30'] == '2003.pendf'
    assert outputs['tape50'] == '2003.02c'
    assert outputs['tape70'] == '2003.02c.xsd'

@pytest.mark.njoy
def test_njoy_process_no_gaspr():
    """Test njoy.process without gaspr"""
    endftape = os.path.join(os.path.dirname(__file__), "data", "n-002_He_003.endf")
    input, inputs, outputs = sandy.njoy.process(endftape, dryrun=True, broadr=False, gaspr=False)
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
225 1 1 20 32 0 /
293.6 /
1.00E+10 /
0 /
moder
-26 30 /
acer
-21 -26 0 50 70 /
1 0 1 .02 0 /
'sandy runs acer'/
225 293.6 /
1 1 /
/
stop"""
    assert input == text
    assert inputs['tape20'] == endftape
    assert outputs['tape30'] == '2003.pendf'
    assert outputs['tape50'] == '2003.02c'
    assert outputs['tape70'] == '2003.02c.xsd'

@pytest.mark.njoy
def test_njoy_process_no_thermr():
    """Test njoy.process without thermr"""
    endftape = os.path.join(os.path.dirname(__file__), "data", "n-002_He_003.endf")
    input, inputs, outputs = sandy.njoy.process(endftape, dryrun=True, broadr=False, gaspr=False,
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
225 1 1 20 32 0 /
293.6 /
1.00E+10 /
0 /
moder
-25 30 /
acer
-21 -25 0 50 70 /
1 0 1 .02 0 /
'sandy runs acer'/
225 293.6 /
1 1 /
/
stop"""
    assert input == text
    assert inputs['tape20'] == endftape
    assert outputs['tape30'] == '2003.pendf'
    assert outputs['tape50'] == '2003.02c'
    assert outputs['tape70'] == '2003.02c.xsd'

@pytest.mark.njoy
def test_njoy_process_no_acer():
    """Test njoy.process without acer"""
    endftape = os.path.join(os.path.dirname(__file__), "data", "n-002_He_003.endf")
    input, inputs, outputs = sandy.njoy.process(endftape, dryrun=True, broadr=False, gaspr=False,
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
225 1 1 20 32 0 /
293.6 /
1.00E+10 /
0 /
moder
-25 30 /
stop"""
    assert input == text
    assert inputs['tape20'] == endftape
    assert outputs['tape30'] == '2003.pendf'

@pytest.mark.njoy
def test_njoy_process_no_purr():
    """Test njoy.process without acer"""
    endftape = os.path.join(os.path.dirname(__file__), "data", "n-002_He_003.endf")
    input, inputs, outputs = sandy.njoy.process(endftape, dryrun=True, broadr=False, gaspr=False, 
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
-24 30 /
stop"""
    assert input == text
    assert inputs['tape20'] == endftape
    assert outputs['tape30'] == '2003.pendf'

@pytest.mark.njoy
def test_njoy_process_no_heatr():
    """Test njoy.process without heatr"""
    endftape = os.path.join(os.path.dirname(__file__), "data", "n-002_He_003.endf")
    input, inputs, outputs = sandy.njoy.process(endftape, dryrun=True, broadr=False, gaspr=False,
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
-22 30 /
stop"""
    assert input == text
    assert inputs['tape20'] == endftape
    assert outputs['tape30'] == '2003.pendf'

@pytest.mark.njoy
def test_njoy_process_no_keep_pendf():
    """Test njoy.process and do not keep pendf"""
    endftape = os.path.join(os.path.dirname(__file__), "data", "n-002_He_003.endf")
    input, inputs, outputs = sandy.njoy.process(endftape, dryrun=True, broadr=False, gaspr=False,
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
    assert inputs['tape20'] == endftape
    assert not outputs

@pytest.mark.njoy
def test_njoy_process_pendftape():
    """Test njoy.process using argument pendftape (skip reconr)"""
    endftape = os.path.join(os.path.dirname(__file__), "data", "n-002_He_003.endf")
    pendftape = "pendf"
    input, inputs, outputs = sandy.njoy.process(endftape, pendftape=pendftape, dryrun=True, broadr=False, gaspr=False,
                               thermr=False, acer=False, purr=False, heatr=False, keep_pendf=False)
    text = """moder
20 -21 /
moder
99 -22 /
stop"""
    assert input == text
    assert inputs['tape20'] == endftape
    assert inputs['tape99'] == pendftape

@pytest.mark.njoy
def test_njoy_process_temperatures():
    """Test njoy.process for different temperatures"""
    endftape = os.path.join(os.path.dirname(__file__), "data", "n-002_He_003.endf")
    pendftape = "pendf"
    input, inputs, outputs = sandy.njoy.process(endftape, pendftape=pendftape, dryrun=True, gaspr=False,
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
    assert inputs['tape20'] == endftape
    assert inputs['tape99'] == pendftape
    assert not outputs

@pytest.mark.njoy
def test_njoy_process_acer():
    """Test njoy.process for acer at different temperatures"""
    endftape = os.path.join(os.path.dirname(__file__), "data", "n-002_He_003.endf")
    pendftape = "pendf"
    input, inputs, outputs = sandy.njoy.process(endftape, pendftape=pendftape, dryrun=True, broadr=False, gaspr=False,
                               thermr=False, acer=True, purr=False, heatr=False, keep_pendf=False,
                               temperatures=[300, 600.0000, 900.001])
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
    assert inputs['tape20'] == endftape
    assert inputs['tape99'] == pendftape
    assert outputs['tape50'] == '2003.03c'
    assert outputs['tape70'] == '2003.03c.xsd'
    assert outputs['tape51'] == '2003.06c'
    assert outputs['tape71'] == '2003.06c.xsd'
    assert outputs['tape52'] == '2003.09c'
    assert outputs['tape72'] == '2003.09c.xsd'

@pytest.mark.njoy
def test_njoy_process_suffixes():
    """Test njoy.process for acer at different temperatures"""
    endftape = os.path.join(os.path.dirname(__file__), "data", "n-002_He_003.endf")
    pendftape = "pendf"
    input, inputs, outputs = sandy.njoy.process(endftape, pendftape="pendf", dryrun=True, broadr=False, gaspr=False,
                               thermr=False, acer=True, purr=False, heatr=False, keep_pendf=False,
                               temperatures=[300, 600.0000, 900.001], suffixes=["01", "02", "06"])
    text = """moder
20 -21 /
moder
99 -22 /
acer
-21 -22 0 50 70 /
1 0 1 .01 0 /
'sandy runs acer'/
225 300.0 /
1 1 /
/
acer
-21 -22 0 51 71 /
1 0 1 .02 0 /
'sandy runs acer'/
225 600.0 /
1 1 /
/
acer
-21 -22 0 52 72 /
1 0 1 .06 0 /
'sandy runs acer'/
225 900.0 /
1 1 /
/
stop"""
    assert input == text
    assert inputs['tape20'] == endftape
    assert inputs['tape99'] == pendftape
    assert outputs['tape50'] == '2003.01c'
    assert outputs['tape70'] == '2003.01c.xsd'
    assert outputs['tape51'] == '2003.02c'
    assert outputs['tape71'] == '2003.02c.xsd'
    assert outputs['tape52'] == '2003.06c'
    assert outputs['tape72'] == '2003.06c.xsd'

@pytest.mark.njoy
def test_njoy_process_sig0():
    """Test njoy.process for different sig0"""
    endftape = os.path.join(os.path.dirname(__file__), "data", "n-002_He_003.endf")
    pendftape = "pendf"
    input, inputs, outputs = sandy.njoy.process(endftape, pendftape=pendftape, dryrun=True, broadr=False, gaspr=False,
                               thermr=False, acer=False, purr=True, heatr=False, keep_pendf=False,
                               sig0=[1e10, 1E9, 100000000])
    text = """moder
20 -21 /
moder
99 -22 /
purr
-21 -22 -23 /
225 1 3 20 32 0 /
293.6 /
1.00E+10 1.00E+09 1.00E+08 /
0 /
stop"""
    assert input == text
    assert inputs['tape20'] == endftape
    assert inputs['tape99'] == pendftape
    assert not outputs

@pytest.mark.njoy
@pytest.mark.njoy_exe
def test_njoy_process(tmpdir):
    """Test njoy.process for ENDF/B-VIII.0 He-3.
    Check that desired outputs are produced and that xsdir files are correctly updated.
    """
    endftape = os.path.join(os.path.dirname(__file__), "data", "n-002_He_003.endf")
    wdir = str(tmpdir)
    input, inputs, outputs = sandy.njoy.process(endftape, temperatures=[300, 600, 900], suffixes=["03", "06", "15"], tag="_b71", wdir=wdir,
                               thermr=False)
    assert inputs['tape20'] == endftape
    assert outputs['tape30'] == os.path.join(wdir, '2003_b71.pendf')
    assert os.path.isfile(outputs['tape30'])
    assert outputs['tape50'] == os.path.join(wdir, '2003_b71.03c')
    assert os.path.isfile(outputs['tape50'])
    assert outputs['tape70'] == os.path.join(wdir, '2003_b71.03c.xsd')
    assert os.path.isfile(outputs['tape70'])
    assert outputs['tape51'] == os.path.join(wdir, '2003_b71.06c')
    assert os.path.isfile(outputs['tape51'])
    assert outputs['tape71'] == os.path.join(wdir, '2003_b71.06c.xsd')
    assert os.path.isfile(outputs['tape71'])
    assert outputs['tape52'] == os.path.join(wdir, '2003_b71.15c')
    assert os.path.isfile(outputs['tape52'])
    assert outputs['tape72'] == os.path.join(wdir, '2003_b71.15c.xsd')
    assert os.path.isfile(outputs['tape72'])
    for ace in ['2003_b71.03c', '2003_b71.06c', '2003_b71.15c']:
        xsdargs = open(os.path.join(wdir, ace) + ".xsd").read().split()
        assert len(xsdargs) == 10
        assert xsdargs[0] == "2003{}".format(os.path.splitext(ace)[1])
        assert xsdargs[2] == os.path.join(wdir, ace)
        assert xsdargs[3] == "0"

@pytest.mark.njoy
@pytest.mark.njoy_exe
def test_njoy_process_1(tmpdir):
    """Test njoy.process for ENDF/B-VIII.0 He-3.
    Check that no output is produced.
    """
    endftape = os.path.join(os.path.dirname(__file__), "data", "n-002_He_003.endf")
    wdir = str(tmpdir)
    input, inputs, outputs = sandy.njoy.process(endftape, wdir=wdir, thermr=False, acer=False, keep_pendf=False)
    assert not os.listdir(wdir)

@pytest.mark.njoy
@pytest.mark.njoy_exe
def test_njoy_process_2(tmpdir):
    """Test njoy.process for ENDF/B-VIII.0 Co-58m.
    Check that ZA of desired outputs is changed because isotope is metastable.
    """
    endftape = os.path.join(os.path.dirname(__file__), "data", "n-027_Co_058m1.endf")
    wdir = str(tmpdir)
    input, inputs, outputs = sandy.njoy.process(endftape, wdir=wdir, thermr=False, keep_pendf=True, route="1")
    assert inputs['tape20'] == endftape
    assert outputs['tape30'] == os.path.join(wdir, '27458.pendf')
    assert os.path.isfile(outputs['tape30'])
    assert outputs['tape50'] == os.path.join(wdir, '27458.02c')
    assert os.path.isfile(outputs['tape50'])
    assert outputs['tape70'] == os.path.join(wdir, '27458.02c.xsd')
    assert os.path.isfile(outputs['tape70'])
    xsdargs = open(outputs['tape70']).read().split()
    assert len(xsdargs) == 10
    assert xsdargs[0] == "27458.02c"
    assert xsdargs[2] == outputs['tape50']
    assert xsdargs[3] == "1"

@pytest.mark.njoy
@pytest.mark.njoy_exe
def test_njoy_process_addpath(tmpdir):
    """Test add_path keyword"""
    endftape = os.path.join(os.path.dirname(__file__), "data", "n-002_He_003.endf")
    wdir = str(tmpdir)
    input, inputs, outputs = sandy.njoy.process(endftape, wdir=wdir, thermr=False, gaspr=False, heatr=False, purr=False, addpath="")
    text = open(outputs['tape70']).read()
    assert text == '2003.02c 2.989032 2003.02c 0 1 1 7108 0 0 2.530E-08'
    input, inputs, outputs = sandy.njoy.process(endftape, wdir=wdir, thermr=False, gaspr=False, heatr=False, purr=False, addpath="aaa")
    text = open(outputs['tape70']).read()
    assert text == '2003.02c 2.989032 aaa/2003.02c 0 1 1 7108 0 0 2.530E-08'
    input, inputs, outputs = sandy.njoy.process(endftape, wdir=wdir, thermr=False, gaspr=False, heatr=False, purr=False, addpath=None)
    text = open(outputs['tape70']).read()
    assert text == '2003.02c 2.989032 {} 0 1 1 7108 0 0 2.530E-08'.format(outputs['tape50'])
    input, inputs, outputs = sandy.njoy.process(endftape, wdir=wdir, thermr=False, gaspr=False, heatr=False, purr=False)
    text = open(outputs['tape70']).read()
    assert text == '2003.02c 2.989032 {} 0 1 1 7108 0 0 2.530E-08'.format(outputs['tape50'])

@pytest.mark.njoy
def test_moder_1():
    """Test moder with default parameters"""
    text = sandy.njoy._moder_input(-20111111, 21)
    assert text == 'moder\n-20111111 21 /\n'
    with pytest.raises(Exception):
        sandy.njoy._moder_input("aaa", 21)
    with pytest.raises(Exception):
        sandy.njoy._moder_input(-80, 15.5)

@pytest.mark.njoy
def test_reconr_1():
    """Test reconr with default parameters"""
    text = sandy.njoy._reconr_input(-20, -21, 200)
    assert text == "reconr\n-20 -21 /\n'sandy runs njoy'/\n200 0 0 /\n0.001 0. /\n0/\n"

@pytest.mark.njoy
def test_reconr_2():
    """Test reconr parameter err"""
    text = sandy.njoy._reconr_input(-20, -21, 200, err=10.0)
    assert text == "reconr\n-20 -21 /\n'sandy runs njoy'/\n200 0 0 /\n10.0 0. /\n0/\n"

@pytest.mark.njoy
def test_reconr_3():
    """TTest reconr parameter header"""
    text = sandy.njoy._reconr_input(-20, -21, 200, header="aaa")
    assert text == "reconr\n-20 -21 /\n'aaa'/\n200 0 0 /\n0.001 0. /\n0/\n"

@pytest.mark.njoy
def test_broadr_1():
    """Test broadr with default parameters"""
    text = sandy.njoy._broadr_input(-20, -21, -22, 200)
    assert text == 'broadr\n-20 -21 -22 /\n200 1 0 0 0. /\n0.001 /\n293.6 /\n0 /\n'

@pytest.mark.njoy
def test_broadr_2():
    """Test broadr parameter temperatures"""
    text = sandy.njoy._broadr_input(-20, -21, -22, 200, temperatures=[900.51, 1E3])
    assert text == 'broadr\n-20 -21 -22 /\n200 2 0 0 0. /\n0.001 /\n900.5 1000.0 /\n0 /\n'
    with pytest.raises(Exception):
        sandy.njoy._broadr_input(-20, -21, -22, 200, temperatures=["aaa"])

@pytest.mark.njoy
def test_broadr_3():
    """Test broadr parameter err"""
    text = sandy.njoy._broadr_input(-20, -21, -22, 200, err=0.1)
    assert text == 'broadr\n-20 -21 -22 /\n200 1 0 0 0. /\n0.1 /\n293.6 /\n0 /\n'

@pytest.mark.njoy
def test_thermr_1():
    """Test thermr with default parameters"""
    text = sandy.njoy._thermr_input(-20, -21, -22, 200)
    assert text == 'thermr\n-20 -21 -22 /\n0 200 20 1 1 0 0 1 221 0 /\n293.6 /\n0.001 10 /\n'

@pytest.mark.njoy
def test_thermr_2():
    """Test thermr parameter temperatures"""
    text = sandy.njoy._thermr_input(-20, -21, -22, 200, temperatures=[900.51, 1E3])
    assert text == 'thermr\n-20 -21 -22 /\n0 200 20 2 1 0 0 1 221 0 /\n900.5 1000.0 /\n0.001 10 /\n'

@pytest.mark.njoy
def test_thermr_3():
    """Test thermr parameter err"""
    text = sandy.njoy._thermr_input(-20, -21, -22, 200, err=100)
    assert text == 'thermr\n-20 -21 -22 /\n0 200 20 1 1 0 0 1 221 0 /\n293.6 /\n100 10 /\n'

@pytest.mark.njoy
def test_thermr_4():
    """Test thermr parameter angles"""
    text = sandy.njoy._thermr_input(-20, -21, -22, 200, angles=31)
    assert text == 'thermr\n-20 -21 -22 /\n0 200 31 1 1 0 0 1 221 0 /\n293.6 /\n0.001 10 /\n'
    with pytest.raises(Exception):
        sandy.njoy._thermr_input(-20, -21, -22, 200, angles=31.2)
    with pytest.raises(Exception):
        sandy.njoy._thermr_input(-20, -21, -22, 200, angles="aaa")

@pytest.mark.njoy
def test_thermr_5():
    """Test thermr parameter emax"""
    text = sandy.njoy._thermr_input(-20, -21, -22, 200, emax=4)
    assert text == 'thermr\n-20 -21 -22 /\n0 200 20 1 1 0 0 1 221 0 /\n293.6 /\n0.001 4 /\n'

@pytest.mark.njoy
def test_thermr_6():
    """Test thermr parameter iprint"""
    text = sandy.njoy._thermr_input(-20, -21, -22, 200, iprint=True)
    assert text == 'thermr\n-20 -21 -22 /\n0 200 20 1 1 0 0 1 221 1 /\n293.6 /\n0.001 10 /\n'
    with pytest.raises(Exception):
        sandy.njoy._thermr_input(-20, -21, -22, 200, iprint="aa")

@pytest.mark.njoy
def test_purr_1():
    """Test purr with default parameters"""
    text = sandy.njoy._purr_input(-20, -21, -22, 200)
    assert text == 'purr\n-20 -21 -22 /\n200 1 1 20 32 0 /\n293.6 /\n1.00E+10 /\n0 /\n'

@pytest.mark.njoy
def test_purr_2():
    """Test purr parameter temperatures"""
    text = sandy.njoy._purr_input(-20, -21, -22, 200, temperatures=[5, 10])
    assert text == 'purr\n-20 -21 -22 /\n200 2 1 20 32 0 /\n5.0 10.0 /\n1.00E+10 /\n0 /\n'

@pytest.mark.njoy
def test_purr_3():
    """Test purr parameter sig0"""
    text = sandy.njoy._purr_input(-20, -21, -22, 200, sig0=[1e8, 10.123])
    assert text == 'purr\n-20 -21 -22 /\n200 1 2 20 32 0 /\n293.6 /\n1.00E+08 1.01E+01 /\n0 /\n'
    with pytest.raises(Exception):
        sandy.njoy._purr_input(-20, -21, -22, 200, sig0=["aa"])

@pytest.mark.njoy
def test_purr_4():
    """Test purr parameter bins"""
    text = sandy.njoy._purr_input(-20, -21, -22, 200, bins=5)
    assert text == 'purr\n-20 -21 -22 /\n200 1 1 5 32 0 /\n293.6 /\n1.00E+10 /\n0 /\n'
    with pytest.raises(Exception):
        sandy.njoy._purr_input(-20, -21, -22, 200, bins="aaa")
    with pytest.raises(Exception):
        sandy.njoy._purr_input(-20, -21, -22, 200, bins=20.1)

@pytest.mark.njoy
def test_purr_5():
    """Test purr parameter iprint"""
    text = sandy.njoy._purr_input(-20, -21, -22, 200, iprint=True)
    assert text == 'purr\n-20 -21 -22 /\n200 1 1 20 32 1 /\n293.6 /\n1.00E+10 /\n0 /\n'
    with pytest.raises(Exception):
        sandy.njoy._purr_input(-20, -21, -22, 200, iprint="aa")

@pytest.mark.njoy
def test_purr_6():
    """Test purr parameter ladders"""
    text = sandy.njoy._purr_input(-20, -21, -22, 200, ladders=2)
    assert text == 'purr\n-20 -21 -22 /\n200 1 1 20 2 0 /\n293.6 /\n1.00E+10 /\n0 /\n'
    with pytest.raises(Exception):
        sandy.njoy._purr_input(-20, -21, -22, 200, ladders="aaa")
    with pytest.raises(Exception):
        sandy.njoy._purr_input(-20, -21, -22, 200, ladders=20.1)

@pytest.mark.njoy
def test_gaspr_1():
    """Test gaspr with default parameters"""
    text = sandy.njoy._gaspr_input(-20111111, 21, 0)
    assert text == 'gaspr\n-20111111 21 0 /\n'

@pytest.mark.njoy
def test_unresr_1():
    """Test unresr with default parameters"""
    text = sandy.njoy._unresr_input(-20, -21, -22, 200)
    assert text == 'unresr\n-20 -21 -22 /\n200 1 1 0 /\n293.6 /\n1.00E+10 /\n0 /\n'

@pytest.mark.njoy
def test_unresr_2():
    """Test unresr parameter temperatures"""
    text = sandy.njoy._unresr_input(-20, -21, -22, 200, temperatures=[5, 10])
    assert text == 'unresr\n-20 -21 -22 /\n200 2 1 0 /\n5.0 10.0 /\n1.00E+10 /\n0 /\n'
    with pytest.raises(Exception):
        sandy.njoy._unresr_input(-20, -21, -22, 200, temperatures=["aa"])


@pytest.mark.njoy
def test_unresr_3():
    """Test unresr parameter sig0"""
    text = sandy.njoy._unresr_input(-20, -21, -22, 200, sig0=[5, 10])
    assert text == 'unresr\n-20 -21 -22 /\n200 1 2 0 /\n293.6 /\n5.00E+00 1.00E+01 /\n0 /\n'
    with pytest.raises(Exception):
        sandy.njoy._unresr_input(-20, -21, -22, 200, sig0=["aa"])

@pytest.mark.njoy
def test_unresr_4():
    """Test unresr parameter iprint"""
    text = sandy.njoy._unresr_input(-20, -21, -22, 200, iprint=True)
    assert text == 'unresr\n-20 -21 -22 /\n200 1 1 1 /\n293.6 /\n1.00E+10 /\n0 /\n'
    with pytest.raises(Exception):
        sandy.njoy._unresr_input(-20, -21, -22, 200, iprint="aa")

@pytest.mark.njoy
def test_heatr_1():
    """Test heatr with default parameters"""
    text = sandy.njoy._heatr_input(-20, -21, -22, 200, [10, 11])
    assert text == 'heatr\n-20 -21 -22 0 /\n200 2 0 0 0 0 /\n10 11 /\n'
    with pytest.raises(Exception):
        sandy.njoy._unresr_input(-20, -21, -22, 200, ["aa"])

@pytest.mark.njoy
def test_heatr_2():
    """Test heatr parameter local"""
    text = sandy.njoy._heatr_input(-20, -21, -22, 200, [10, 11], local=True)
    assert text == 'heatr\n-20 -21 -22 0 /\n200 2 0 0 1 0 /\n10 11 /\n'
    with pytest.raises(Exception):
        sandy.njoy._heatr_input(-20, -21, -22, 200, [10, 11], local="aa")

@pytest.mark.njoy
def test_heatr_3():
    """Test heatr parameter iprint"""
    text = sandy.njoy._heatr_input(-20, -21, -22, 200, [10, 11], iprint=True)
    assert text == 'heatr\n-20 -21 -22 0 /\n200 2 0 0 0 1 /\n10 11 /\n'
    with pytest.raises(Exception):
        sandy.njoy._heatr_input(-20, -21, -22, 200, [10, 11], iprint="aa")

@pytest.mark.njoy
def test_acer_1():
    """Test acer with default parameters"""
    text = sandy.njoy._acer_input(-20, -21, -60, 80, 200)
    assert text == "acer\n-20 -21 0 -60 80 /\n1 0 1 .00 0 /\n'sandy runs acer'/\n200 293.6 /\n1 1 /\n/\n"

@pytest.mark.njoy
def test_acer_2():
    """Test acer parameter temp"""
    text = sandy.njoy._acer_input(-20, -21, -60, 80, 200, temp=-500)
    assert text == "acer\n-20 -21 0 -60 80 /\n1 0 1 .00 0 /\n'sandy runs acer'/\n200 -500.0 /\n1 1 /\n/\n"
    with pytest.raises(Exception):
        sandy.njoy._acer_input(-20, -21, -60, 80, 200, temp="aa")

@pytest.mark.njoy
def test_acer_3():
    """Test acer parameter iprint"""
    text = sandy.njoy._acer_input(-20, -21, -60, 80, 200, iprint=True)
    assert text == "acer\n-20 -21 0 -60 80 /\n1 1 1 .00 0 /\n'sandy runs acer'/\n200 293.6 /\n1 1 /\n/\n"
    with pytest.raises(Exception):
        sandy.njoy._acer_input(-20, -21, -60, 80, 200, iprint="aa")

@pytest.mark.njoy
def test_acer_4():
    """Test acer parameter itype"""
    text = sandy.njoy._acer_input(-20, -21, -60, 80, 200, itype=2)
    assert text == "acer\n-20 -21 0 -60 80 /\n1 0 2 .00 0 /\n'sandy runs acer'/\n200 293.6 /\n1 1 /\n/\n"
    with pytest.raises(Exception):
        sandy.njoy._acer_input(-20, -21, -60, 80, 200, itype="aa")

@pytest.mark.njoy
def test_acer_5():
    """Test acer parameter suff"""
    text = sandy.njoy._acer_input(-20, -21, -60, 80, 200, suff=5)
    assert text == "acer\n-20 -21 0 -60 80 /\n1 0 1 5 0 /\n'sandy runs acer'/\n200 293.6 /\n1 1 /\n/\n"

@pytest.mark.njoy
def test_acer_6():
    """Test acer parameter header"""
    text = sandy.njoy._acer_input(-20, -21, -60, 80, 200, header="Hi!")
    assert text == "acer\n-20 -21 0 -60 80 /\n1 0 1 .00 0 /\n'Hi!'/\n200 293.6 /\n1 1 /\n/\n"

@pytest.mark.njoy
def test_acer_7():
    """Test acer parameter iprint"""
    text = sandy.njoy._acer_input(-20, -21, -60, 80, 200, photons=False)
    assert text == "acer\n-20 -21 0 -60 80 /\n1 0 1 .00 0 /\n'sandy runs acer'/\n200 293.6 /\n1 0 /\n/\n"
    with pytest.raises(Exception):
        sandy.njoy._acer_input(-20, -21, -60, 80, 200, photons="aa")

@pytest.mark.njoy
def test_process_proton():
    """Test default options for njoy.process_proton"""
    endftape = os.path.join(os.path.dirname(__file__), "data", "O016-p.tendl")
    input, inputs, outputs = sandy.njoy.process_proton(endftape, dryrun=True)
    assert input == "acer\n20 20 0 50 70 /\n1 0 1 .00 0 /\n'sandy runs acer'/\n825 0.0 /\n1 1 /\n/\nstop"
    assert outputs['tape50'] == '8016.00h'
    assert outputs['tape70'] == '8016.00h.xsd'
    assert inputs["tape20"] == endftape

@pytest.mark.njoy
@pytest.mark.njoy_exe
def test_process_proton_2(tmpdir):
    """Test njoy.process for TENDL-2015 O-16.
    Check that desired outputs are produced and that xsdir files are correctly updated.
    """
    endftape = os.path.join(os.path.dirname(__file__), "data", "O016-p.tendl")
    wdir = str(tmpdir)
    input, inputs, outputs = sandy.njoy.process_proton(endftape, wdir=wdir)
    assert input == "acer\n20 20 0 50 70 /\n1 0 1 .00 0 /\n'sandy runs acer'/\n825 0.0 /\n1 1 /\n/\nstop"
    assert outputs['tape50'] == os.path.join(wdir, '8016.00h')
    assert os.path.isfile(outputs['tape50'])
    assert outputs['tape70'] == os.path.join(wdir, '8016.00h.xsd')
    assert os.path.isfile(outputs['tape70'])
    assert inputs["tape20"] == endftape
    xsdargs = open(outputs['tape70']).read().split()
    assert len(xsdargs) == 10
    assert xsdargs[0] == "8016.00h"
    assert xsdargs[2] == outputs['tape50']
    assert xsdargs[3] == "0"

@pytest.mark.njoy
def test_get_suffix():
    """Test function get_suffix"""
    for tmp, ext in sandy.njoy.tmp2ext.items():
        assert sandy.njoy.get_suffix(tmp, 0) == ext
    for tmp, ext in sandy.njoy.tmp2ext.items():
        assert sandy.njoy.get_suffix(tmp, 1) == ext
    for tmp, ext in sandy.njoy.tmp2ext_meta.items():
        assert sandy.njoy.get_suffix(tmp, 1, "aleph") == ext
    for tmp, ext in sandy.njoy.tmp2ext_meta.items():
        assert sandy.njoy.get_suffix(tmp, 0, "aleph") == sandy.njoy.tmp2ext[tmp]
    with pytest.raises(Exception):
        sandy.njoy.get_suffix(150, 0)
    assert sandy.njoy.get_suffix(324, 0) == "03"
    assert sandy.njoy.get_suffix(326, 0) == "35"
    assert sandy.njoy.get_suffix(324, 2, method="aleph") == "31"
    assert sandy.njoy.get_suffix(326, 2, method="aleph") == "32"