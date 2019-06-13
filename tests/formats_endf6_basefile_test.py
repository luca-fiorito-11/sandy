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
def basetext2():
    """Random lines taken from the ENDF/B-VII.1 evaluation for He-4.
    MAT number was changed.
    """
    text = """
 0.000000+0 0.000000+0          0          0          0          02631 1  099999
 0.000000+0 0.000000+0          0          0          0          02631 0  0    0
 2.004000+3 3.968219+0          0          0          1          02631 2151    1
 2.004000+3 1.000000+0          0          0          1          02631 2151    2
 1.000000-5 1.000000+5          0          0          0          02631 2151    3
 0.000000+0 2.457900-1          0          0          0          02631 2151    4
 0.000000+0 0.000000+0          0          0          0          02631 2  099999
 0.000000+0 0.000000+0          0          0          0          02631 0  0    0
 2.004000+3 3.968219+0          0          0          0          0 228 3  1    1
 0.000000+0 0.000000+0          0          0          1        152 228 3  1    2
 2.004000+3 3.968219+0          0          0          0          0 228 3  2    1
 0.000000+0 0.000000+0          0          0          1        152 228 3  2    2
        152          2                                             228 3  2    3
 1.000000-5 7.669738-1 2.530000-2 7.669738-1 1.000000+0 7.669738-1 228 3  2    4
"""
    return text

@pytest.fixture(scope="module")
def basefile1(basetext1):
    return sandy.formats.endf6._BaseFile.from_text(basetext1)

    
@pytest.fixture(scope="module")
def basefile2(basetext2):
    return sandy.formats.endf6._BaseFile.from_text(basetext2)

@pytest.mark.formats
@pytest.mark.endf6
def test_BaseFile_mat(basefile1):
    """Test property mat"""
    assert basefile1.mat == [2631]
    assert type(basefile1.mat) == list

@pytest.mark.formats
@pytest.mark.endf6
def test_BaseFile_mf(basefile1):
    """Test property mf"""
    assert basefile1.mf == [1, 2, 3, 34]
    assert type(basefile1.mf) == list

@pytest.mark.formats
@pytest.mark.endf6
def test_BaseFile_mt(basefile1):
    """Test property mt"""
    assert basefile1.mt == [1, 2, 151, 451]
    assert type(basefile1.mt) == list

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
    """Test that BaseFile returns error when empty dataframe is given"""
    with pytest.raises(Exception):
        sandy.formats.endf6.BaseFile(pd.DataFrame())

@pytest.mark.formats
@pytest.mark.endf6
def test_BaseFile_add_sections(basefile1, basefile2):
    """Test BaseFile method add_sections"""
    tape = basefile1.add_sections(basefile2)
    assert tape.TEXT.loc[228,3,1] == basefile2.TEXT.loc[228,3,1]
    assert tape.TEXT.loc[228,3,2] == basefile2.TEXT.loc[228,3,2]
    assert tape.TEXT.loc[2631,1,451] == basefile1.TEXT.loc[2631,1,451]
    assert tape.TEXT.loc[2631,2,151] == basefile2.TEXT.loc[2631,2,151]
    assert tape.TEXT.loc[2631,2,151] != basefile1.TEXT.loc[2631,2,151]
    assert tape.TEXT.loc[2631,3,1] == basefile1.TEXT.loc[2631,3,1]
    assert tape.TEXT.loc[2631,34,2] == basefile1.TEXT.loc[2631,34,2]
    tape1 = basefile2.add_sections(basefile1)
    assert tape1.TEXT.loc[228,3,1] == tape.TEXT.loc[228,3,1]
    assert tape1.TEXT.loc[228,3,2] == tape.TEXT.loc[228,3,2]
    assert tape1.TEXT.loc[2631,1,451] == tape.TEXT.loc[2631,1,451]
    assert tape1.TEXT.loc[2631,2,151] != tape.TEXT.loc[2631,2,151]
    assert tape1.TEXT.loc[2631,2,151] == basefile1.TEXT.loc[2631,2,151]
    assert tape1.TEXT.loc[2631,3,1] == tape.TEXT.loc[2631,3,1]
    assert tape1.TEXT.loc[2631,34,2] == tape.TEXT.loc[2631,34,2]

@pytest.mark.formats
@pytest.mark.endf6
def test_BaseFile_delete_sections(basefile2):
    """Test BaseFile method add_sections"""
    tape = basefile2.delete_sections((None, 3, None))
    assert tape.mat == [2631]
    assert tape.mf == [2]
    assert tape.mt == [151]
    tape1 = basefile2.delete_sections((228, 3, 1))
    assert tape1.mat == [228, 2631]
    assert tape1.mf == [2, 3]
    assert tape1.mt == [2, 151]
    tape2 = basefile2.delete_sections((None, None, 1))
    assert tape2.equals(tape1)
    tape3 = basefile2.delete_sections((228, 3, 1), (2631, None, None))
    assert tape3.mat == [228]
    assert tape3.mf == [3]
    assert tape3.mt == [2]
    tape4 = basefile2.delete_sections((228, 3, 1), (2631, None, None), (None, None, 15))
    assert tape4.equals(tape3)
    with pytest.raises(Exception):
        basefile2.delete_sections((228, None, None), (2631, None, None))

@pytest.mark.formats
@pytest.mark.endf6
def test_BaseFile_get_file_format():
    """Check that method get_file_format works as expected"""
    text = " 2.605600+4 5.545440+1          1          0          2          02631 1451    1\n"
    ftype = sandy.formats.endf6._BaseFile.from_text(text).get_file_format()
    assert ftype == "endf6"
    text = " 2.605600+4 5.545440+1          2          0          2          02631 1451    1\n"
    ftype = sandy.formats.endf6._BaseFile.from_text(text).get_file_format()
    assert ftype == "pendf"
    text = " 2.605600+4 5.545440+1          1          0        -11          02631 1451    1\n"
    ftype = sandy.formats.endf6._BaseFile.from_text(text).get_file_format()
    assert ftype == "errorr"
    text = " 2.605600+4 5.545440+1          1          0        -12          02631 1451    1\n"
    ftype = sandy.formats.endf6._BaseFile.from_text(text).get_file_format()
    assert ftype == "errorr"
    text = " 2.605600+4 5.545440+1          1          0         -1          02631 1451    1\n"
    ftype = sandy.formats.endf6._BaseFile.from_text(text).get_file_format()
    assert ftype == "gendf"