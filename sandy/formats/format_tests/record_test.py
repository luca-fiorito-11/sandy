# -*- coding: utf-8 -*-
"""
Created on Thu Jul 12 11:25:35 2018

@author: Fiorito_L
"""

import pytest

import numpy as np

from ...data import Pu9
from ..records import read_ilist
from ..records import read_dlist
from ..records import read_list
from ..records import read_cont
from ..records import read_tab1
from ..records import write_list
from ..records import write_ilist
from ..records import write_dlist
from ..records import write_cont
from ..records import write_tab1
from ..records import read_control
from ..records import read_text

@pytest.mark.records
@pytest.mark.read
@pytest.mark.ilist
def test_read_ilist():
    string = "       3201          2                                            9437 1452    3"
    L = list(map(int, string[:66].replace("+", "E+").replace("-", "E-").split()))
    assert L == read_ilist(string)[:2]
    assert read_ilist(string)[2:] == [0, 0, 0, 0]

@pytest.mark.records
@pytest.mark.read
@pytest.mark.dlist
def test_read_dlist():
    string = " 1.000000-5 2.868348+0 3.000000-5 2.868348+0 1.000000-3 2.868348+09437 1452    4"
    array = np.array(list(map(float,string[:66].replace("+", "E+").replace("-", "E-").split())))
    assert (array == read_dlist(string)).all()

@pytest.mark.records
@pytest.mark.read
@pytest.mark.write
@pytest.mark.list
def test_read_list():
    (C1, C2, L1, L2, NPL, N2, B), ipos = read_list(Pu9.endf6, 1809)
    assert len(B) == NPL
    text = write_list(C1, C2, L1, L2, N2, B)
    original_text = [x[:66] for x in Pu9.endf6[1809:ipos]]
    assert text == original_text

@pytest.mark.records
@pytest.mark.read
@pytest.mark.cont
def test_read_cont():
    (C1, C2, L1, L2, N1, N2), ipos = read_cont(Pu9.endf6, 1)
    assert ipos == 2
    assert C1 == 94239.0
    assert C2 == 236.9986
    assert L1 == 1
    assert L2 == 1
    assert N1 == 2
    assert N2 == 1

@pytest.mark.records
@pytest.mark.read
@pytest.mark.write
@pytest.mark.tab1
def test_read_tab1():
    (C1, C2, L1, L2, NR, NP, NBT, INT, x, y), ipos = read_tab1(Pu9.endf6, 738)
    mat, mf, mt, ns = read_control(Pu9.endf6[ipos])
    assert mat == 9437
    assert mf == 1
    assert mt == 0
    assert ns == 99999
    text = write_tab1(C1, C2, L1, L2, NBT, INT, x, y)
    original_text = [x[:66] for x in Pu9.endf6[738:ipos]]
    assert text == original_text

@pytest.mark.records
@pytest.mark.read
@pytest.mark.text
def test_read_text():
    string, ipos = read_text(Pu9.endf6, 6)
    assert ipos == 7
    assert string == " JEFF33T3             DIST-DEC17 REV3-DEC17            20171231   "

@pytest.mark.records
@pytest.mark.write
@pytest.mark.cont
def test_write_cont():
    text = write_cont(94239, 236.9986, 1, 1, 2, 1)
    assert len(text) == 1
    assert text[0] == ' 9.423900+4 2.369986+2          1          1          2          1'

@pytest.mark.records
@pytest.mark.write
@pytest.mark.ilist
def test_write_ilist():
    text = write_ilist([3, 4, 5])
    assert len(text) == 1
    assert text[0] == '          3          4          5                                 '

@pytest.mark.records
@pytest.mark.write
@pytest.mark.dlist
def test_write_dlist():
    text = write_dlist([94239, 236.9986, 1, 1, 2, 1, 94239, 236.9986, 1, 1, 2])
    assert len(text) == 2
    assert text[0] == ' 9.423900+4 2.369986+2 1.000000+0 1.000000+0 2.000000+0 1.000000+0'
    assert text[1] == ' 9.423900+4 2.369986+2 1.000000+0 1.000000+0 2.000000+0           '