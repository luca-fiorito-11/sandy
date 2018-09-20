# -*- coding: utf-8 -*-
"""
Created on Thu Jul 12 11:25:35 2018

@author: Fiorito_L
"""

import pytest

import numpy as np
import rwf

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
    assert text[0] == ' 0.00000000 0.00000000          0          0          8          0'
    assert text[-1] == ' 1.63478100 3.55460000                                            '

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
    assert text[0] == ' 0.00000000 0.00000000          0          0          1       3201'
    assert text[-1] == ' 29000000.0 6.49123000 29500000.0 6.53911000 30000000.0 6.62055000'

@pytest.mark.records
@pytest.mark.read
@pytest.mark.text
def test_read_text():
    string, ipos = read_text(Pu9.endf6, 6)
    assert ipos ==  7
    assert string == " JEFF33T3             DIST-DEC17 REV3-DEC17            20171231   "

@pytest.mark.records
@pytest.mark.write
@pytest.mark.cont
def test_write_cont():
    text = write_cont(94239, 236.9986, 1, 1, 2, 1)
    assert len(text) == 1
    assert text[0] == ' 94239.0000 236.998600          1          1          2          1'

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
    assert text[0] == ' 94239.0000 236.998600 1.00000000 1.00000000 2.00000000 1.00000000'
    assert text[1] == ' 94239.0000 236.998600 1.00000000 1.00000000 2.00000000           '

@pytest.mark.records
@pytest.mark.write
@pytest.mark.float
def test_write_float():
    values = {1e11 : ' 1.00000+11',
              1e10 : ' 1.00000+10',
              1e9  : ' 1.000000+9',
              1e8  : ' 100000000.',
              1e7  : ' 10000000.0',
              1e6  : ' 1000000.00',
              1e5  : ' 100000.000',
              1e4  : ' 10000.0000',
              1e3  : ' 1000.00000',
              1e2  : ' 100.000000',
              1e1  : ' 10.0000000',
              1e0  : ' 1.00000000',
              0    : ' 0.000000+0',
              1e-1 : ' 1.000000-1',
              1e-2 : ' 1.000000-2',
              1e-3 : ' 1.000000-3',
              1e-4 : ' 1.000000-4',
              1e-5 : ' 1.000000-5',
              1e-6 : ' 1.000000-6',
              1e-7 : ' 1.000000-7',
              1e-8 : ' 1.000000-8',
              1e-9 : ' 1.000000-9',
              1e-10: ' 1.00000-10',
              1e-11: ' 1.00000-11',
              -1e11 : '-1.00000+11',
              -1e10 : '-1.00000+10',
              -1e9  : '-1.000000+9',
              -1e8  : '-100000000.',
              -1e7  : '-10000000.0',
              -1e6  : '-1000000.00',
              -1e5  : '-100000.000',
              -1e4  : '-10000.0000',
              -1e3  : '-1000.00000',
              -1e2  : '-100.000000',
              -1e1  : '-10.0000000',
              -1e0  : '-1.00000000',
              -1e-1 : '-1.000000-1',
              -1e-2 : '-1.000000-2',
              -1e-3 : '-1.000000-3',
              -1e-4 : '-1.000000-4',
              -1e-5 : '-1.000000-5',
              -1e-6 : '-1.000000-6',
              -1e-7 : '-1.000000-7',
              -1e-8 : '-1.000000-8',
              -1e-9 : '-1.000000-9',
              -1e-10: '-1.00000-10',
              -1e-11: '-1.00000-11',
              }
    for k,v in values.items():
        C1 = np.array(k, dtype=float)
        string = np.array("*"*11)
        rwf.wreal(C1, string)
        string = str(string, 'utf-8')[:11]
        assert v == string