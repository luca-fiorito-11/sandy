# -*- coding: utf-8 -*-
"""
Created on Tue Mar 05 15:59:46 2019

@author: Luca Fiorito
"""

__author__ = "Luca Fiorito"

import pytest

import pandas as pd
import numpy as np

import sandy

@pytest.fixture(scope="module")
def basetext1():
    """Random lines taken from the JEFF-3.3 evaluation for Fe-56.
    """
    text = """
JEFF-3.3 Incident Neutron File                                       1 0  0    0
 2.605600+4 5.545440+1          1          0          2          02631 1451    1
 0.000000+0 0.000000+0          0          0          0          62631 1451    2
 1.000000+0 2.000000+8          3          0         10          32631 1451    3
 0.000000+0 0.000000+0          0          0        854        2922631 1451    4
 2.605600+4 5.545440+1          0          0          1          02631 2151    1
 0.000000+0 0.000000+0          0          0          0          02631 0  0    0
 2.605600+4 1.000000+0          0          0          1          02631 2151    2
 0.000000+0 0.000000+0          0          0          0          02631 1  099999
 1.031400+6 2.299110+0 1.031480+6 2.240840+0 1.031630+6 2.400020+02631 3  1  597
 1.031700+6 2.539370+0 1.031780+6 3.004300+0 1.031850+6 3.364950+02631 3  1  598
 1.031900+6 3.424110+0 1.031930+6 3.459720+0 1.032000+6 3.442310+02631 3  1  599
 sfggd;'hfd;nhf;
 -1 0  0    0
 2.000000+7 0.000000+0                                            263134  2  332
 0.000000+0 0.000000+0          6          6          0          1263134  2  333
"""
    return text

@pytest.fixture(scope="module")
def basefile1(basetext1):
    return sandy.formats.endf6._BaseFile.from_text(basetext1)

@pytest.mark.formats
@pytest.mark.endf6
def test_BaseFile_mat(basefile1):
    """Test property mat"""
    assert basefile1.mat == [2631]

@pytest.mark.formats
@pytest.mark.endf6
def test_BaseFile_mf(basefile1):
    """Test property mf"""
    assert basefile1.mf == [1, 2, 3, 34]

@pytest.mark.formats
@pytest.mark.endf6
def test_BaseFile_mt(basefile1):
    """Test property mt"""
    assert basefile1.mt == [1, 2, 151, 451]

@pytest.mark.formats
@pytest.mark.endf6
def test_BaseFile_index_dtype(basefile1):
    """Test that all indices are integer"""
    assert all(isinstance(x, int) for lev in basefile1.index.levels for x in lev)

@pytest.mark.formats
@pytest.mark.endf6
def test_BaseFile_columns(basefile1):
    """Test that columns values are correct"""
    assert basefile1.TEXT.loc[2631,1,451] == ' 2.605600+4 5.545440+1          1          0          2          02631 1451    1\n' + \
                                             ' 0.000000+0 0.000000+0          0          0          0          62631 1451    2\n' + \
                                             ' 1.000000+0 2.000000+8          3          0         10          32631 1451    3\n' + \
                                             ' 0.000000+0 0.000000+0          0          0        854        2922631 1451    4\n'
    assert basefile1.TEXT.loc[2631,2,151] == ' 2.605600+4 5.545440+1          0          0          1          02631 2151    1\n' + \
                                             ' 2.605600+4 1.000000+0          0          0          1          02631 2151    2\n'
    assert basefile1.TEXT.loc[2631,3,1]   == ' 1.031400+6 2.299110+0 1.031480+6 2.240840+0 1.031630+6 2.400020+02631 3  1  597\n' + \
                                             ' 1.031700+6 2.539370+0 1.031780+6 3.004300+0 1.031850+6 3.364950+02631 3  1  598\n' + \
                                             ' 1.031900+6 3.424110+0 1.031930+6 3.459720+0 1.032000+6 3.442310+02631 3  1  599\n'
    assert basefile1.TEXT.loc[2631,34,2]  == ' 2.000000+7 0.000000+0                                            263134  2  332\n' + \
                                             ' 0.000000+0 0.000000+0          6          6          0          1263134  2  333\n'
    assert basefile1.columns == ["TEXT"]

@pytest.mark.formats
@pytest.mark.endf6
def test_BaseFile_empty():
    """Test that BaseFile returns error when empty dataframe is given."""
    with pytest.raises(Exception):
        sandy.formats.endf6.BaseFile(pd.DataFrame())