# -*- coding: utf-8 -*-
"""
Created on Thu Jun 14 10:22:58 2018

@author: fiorito_l
"""
import sys, pytest, rwf
import numpy as np
from collections import namedtuple

def read_control(string):
    mat = np.array(0, dtype=int)
    mf = np.array(0, dtype=int)
    mt = np.array(0, dtype=int)
    ns = np.array(0, dtype=int)
    io_status = np.array(0, dtype=int)
    rwf.rcontrol(string[66:], io_status, mat, mf, mt, ns)
    if io_status != 0:
        sys.exit("ERROR : line '{}' is not in CONTROL format".format(string))
    return int(mat), int(mf), int(mt), int(ns)

def read_cont(text, ipos):
    """
    Read ``ENDF-6`` ``CONT`` record in formatted fortran.

    Outputs:
        - :``out``: :
            (tuple) content of ``CONT`` record

    Found error in:
        - n-17-Cl-035.jeff32
        - n-3-Li-007.jeff32
        - n-63-Eu-152.jeff32
        - n-63-Eu-153.jeff32
        - n-64-Gd-155.jeff32
        - n-77-Ir-193.jeff32
        - n-90-Th-229.jeff32
        - n-94-Pu-238.jeff32
        - n-94-Pu-241.jeff32
        - n-94-Pu-242.jeff32
        - n-97-Bk-250.jeff32
        - n-98-Cf-254.jeff32
    """
    CONT = namedtuple('CONT', 'C1 C2 L1 L2 N1 N2')
    string = text[ipos]
    c1 = np.array(0.0, dtype=float)
    c2 = np.array(0.0, dtype=float)
    l1 = np.array(0, dtype=int)
    l2 = np.array(0, dtype=int)
    n1 = np.array(0, dtype=int)
    n2 = np.array(0, dtype=int)
    io_status = np.array(0, dtype=int)
    rwf.rcont(string, io_status, c1, c2, l1, l2, n1, n2)
    if io_status != 0:
        sys.exit("ERROR : line '{}' is not in CONT format".format(string))
    ipos += 1
    return CONT(float(c1), float(c2), int(l1), int(l2), int(n1), int(n2)), ipos

def read_text(text, ipos):
    string = text[ipos]
    try:
        ipos += 1
        return text[ipos-1][:66], ipos
    except:
        sys.exit("ERROR : line '{}' is not in TEXT format".format(string))

def read_ilist(string):
    array = np.zeros(6, dtype=int)
    io_status = np.array(0, dtype=int)
    rwf.rilist(string, io_status, array, 6)
    if io_status != 0:
        sys.exit("ERROR : line '{}' is not in ILIST format".format(string))
    return array.tolist()

def read_dlist(string):
    array = np.zeros(6, dtype=float)
    io_status = np.array(0, dtype=int)
    rwf.rlist(string, io_status, array, 6)
    if io_status != 0:
        sys.exit("ERROR : line '{}' is not in DLIST format".format(string))
    return array.tolist()

def read_tab2(text, ipos):
    """
    Read ENDF-6 **tab2** record.
    This record is used to control the tabulation of a 2-dimensional
    function :math:`y(x,z)`.
    It specifies how many values of z are to be given and how to
    interpolate between the successive values of :math:`z`.
    Tabulated values of :math:`y_l(x)` at each value of :math:`z_l`
    are given in ENDF-6 **tab1** or **list** records following the
    ENDF-6 **tab2** record, with the appropriate value of :math:`z`
    in the field designated as *c2*.

    Outputs:
        - c1 :
            (float) first element of string
        - c2 :
            (float) second element of string
        - l1 :
            (int) third element of string
        - l2 :
            (int) fourth element of string
        - nr :
            (int) number of interpolation ranges
        - nz :
            (int) number of :math:`z` to be given
        - nbt :
            (list of int) upper limits of the interpolation ranges
        - interp :
            (list of int) interpolation functions
    """
    import itertools
    TAB2 = namedtuple('TAB2', 'C1 C2 L1 L2 NR NZ NBT INT')
    C, ipos = read_cont(text, ipos)
    iadd = int(np.ceil(C.N1*2/6))
    tab = map(lambda x:read_ilist(x) ,text[ipos:ipos+iadd])
    tab = list(itertools.chain.from_iterable(tab))[:C.N1*2]
    ipos += iadd
    NBT = tab[::2]
    INT = tab[1::2]
    return TAB2(C.C1, C.C2, C.L1, C.L2, C.N1, C.N2, NBT, INT), ipos

def read_tab1(text, ipos):
    """
    Read ENDF-6 **tab1** record.
    These records are used for one-dimensional tabulated functions such as
    :math:`y(x)`.
    The data needed to specify a one-dimensional tabulated function are
    the interpolation tables *nbt* and *int* for each of the *nr* ranges,
    and the *np* tabulated pairs of :math:`x(n)` and :math:`y(n)`.

    Outputs:
        - c1 :
            (float) first element of string
        - c2 :
            (float) second element of string
        - l1 :
            (int) third element of string
        - l2 :
            (int) fourth element of string
        - nr :
            (int) number of interpolation ranges
        - nz :
            (int) number of :math:`z` to be given
        - nbt :
            (list of int) upper limits of the interpolation ranges
        - interp :
            (list of int) interpolation functions
        - x, y :
            tabulated pairs of :math:`x(n)` and :math:`y(n)`
    """
    import itertools
    TAB1 = namedtuple('TAB1', 'C1 C2 L1 L2 NR NP NBT INT x y')
    TAB2, ipos = read_tab2(text, ipos)
    iadd = int(np.ceil(TAB2.NZ*2/6))
    tab = map(lambda x:read_dlist(x) ,text[ipos:ipos+iadd])
    tab = list(itertools.chain.from_iterable(tab))[:TAB2.NZ*2]
    ipos += iadd
    x = np.array(tab[::2], dtype=float)
    y = np.array(tab[1::2], dtype=float)
    return TAB1(*list(TAB2), x, y), ipos

def read_list(text, ipos):
    import itertools
    LIST = namedtuple('LIST', 'C1 C2 L1 L2 NPL N2 B')
    C, ipos = read_cont(text, ipos)
    iadd = int(np.ceil(C.N1/6))
    tab = map(lambda x:read_dlist(x) ,text[ipos:ipos+iadd])
    tab = list(itertools.chain.from_iterable(tab))[:C.N1]
    ipos += iadd
    return LIST(*list(C), tab), ipos

@pytest.mark.unit_test
@pytest.mark.records
def test_read_ilist():
    string = "       3201          2                                            9437 1452    3"
    L = list(map(int, string[:66].replace("+", "E+").replace("-", "E-").split()))
    assert L == read_ilist(string)[:2]
    assert read_ilist(string)[2:] == [0, 0, 0, 0]

@pytest.mark.unit_test
@pytest.mark.records
def test_read_dlist():
    string = " 1.000000-5 2.868348+0 3.000000-5 2.868348+0 1.000000-3 2.868348+09437 1452    4"
    array = np.array(list(map(float,string[:66].replace("+", "E+").replace("-", "E-").split())))
    assert (array == read_dlist(string)).all()

@pytest.mark.unit_test
@pytest.mark.records
def test_read_cont():
    from sandy.data_test import Pu9
    (C1, C2, L1, L2, N1, N2), ipos = read_cont(Pu9.endf6, 1)
    assert ipos == 2
    assert C1 == 94239.0
    assert C2 == 236.9986
    assert L1 == 1
    assert L2 == 1
    assert N1 == 2
    assert N2 == 1

@pytest.mark.unit_test
@pytest.mark.records
def test_read_tab1():
    from sandy.data_test import Pu9
    (C1, C2, L1, L2, NR, NP, NBT, INT, x, y), ipos = read_tab1(Pu9.endf6, 738)
    mat, mf, mt, ns = read_control(Pu9.endf6[ipos])
    assert mat == 9437
    assert mf == 1
    assert mt == 0
    assert ns == 99999