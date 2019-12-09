# -*- coding: utf-8 -*-
"""
Created on Thu Jun 14 10:22:58 2018

@author: fiorito_l
"""
import sys
from collections import namedtuple
import itertools
import pdb

import numpy as np

import sandy

__author__ = "Luca Fiorito"
__all__ = [
        "read_cont",
        "read_tab1",
        "read_tab2",
        "read_list",
        "write_cont",
        "write_tab1",
        "write_tab2",
        "write_list",
        "write_float",
        "write_integer_list",
        "write_float_list",
        "line_numbers",
        "write_line",
        "write_eol",
        ]


line_pattern = "{:<66}{:4d}{:2d}{:3d}{:5d}"



def read_cont(df, ipos):
    """
    Read ``ENDF-6`` ``CONT`` record in formatted fortran.

    Returns
    -------
    `CONT` : `collections.namedtuple`
        - `C1` : `float`
            first element of string
        - `C2` : `float`
            second element of string
        - `L1` : `int`
            third element of string
        - `L1` : `int`
            fourth element of string
        - `N1` : `int`
            fifth element of string
        - `N2` : `int`
            sixth element of string

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
    series = df.iloc[ipos]
    c1 = float(series.C1)
    c2 = float(series.C2)
    l1 = int(series.L1)
    l2 = int(series.L2)
    n1 = int(series.N1)
    n2 = int(series.N2)
    ipos += 1
    return CONT(c1, c2, l1, l2, n1, n2), ipos



def write_cont(C1, C2, L1, L2, N1, N2):
    """
    Write ENDF-6 **cont** record.

    Returns
    -------
    `list` of `str`
        list of 66-characters-long ENDF-6 formatted string 
    
    Warns
    -----
    This function will produce strings longer than 66 characters if integers 
    `> 99999999999` are given.
    """
    lines = [write_float(C1) + write_float(C2) + "{:11d}{:11d}{:11d}{:11d}".format(L1, L2, N1, N2)]
    return lines



def write_integer_list(lst):
    """
    Write list of integers into ENDF-6 format.

    Returns
    -------
    `list` of `str`
        list of 66-characters-long ENDF-6 formatted string 
    """
    itr = sandy.shared.grouper(map("{:11d}".format, lst), 6, fillvalue=" "*11)
    return ["".join(vals) for vals in itr ]



def write_float_list(lst):
    """
    Write list of floats into ENDF-6 format.

    Returns
    -------
    `list` of `str`
        list of 66-characters-long ENDF-6 formatted string 
    """
    itr = sandy.shared.grouper(map(write_float, lst), 6, fillvalue=" "*11)
    return ["".join(vals) for vals in itr ]



def read_tab2(df, ipos):
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

    Returns
    -------
    `TAB1` : `collections.namedtuple`
        - `C1` : `float`
            first element of string
        - `C2` : `float`
            second element of string
        - `L1` : `int`
            third element of string
        - `L1` : `int`
            fourth element of string
        - `NR` : `int`
            number of interpolation ranges
        - `NZ` : `int`
            number of :math:`z` to be given
        - `NBT` : `list` of `int`
            upper limits of the interpolation ranges
        - `INT` : `list` of `int`
            interpolation functions
    """
    TAB2 = namedtuple('TAB2', 'C1 C2 L1 L2 NR NZ NBT INT')
    C, ipos = read_cont(df, ipos)
    L, ipos = _read_list(df, ipos, C.N1 * 2)
    NBT = [int(x) for x in L[::2]]
    INT = [int(x) for x in L[1::2]]
    return TAB2(C.C1, C.C2, C.L1, C.L2, C.N1, C.N2, NBT, INT), ipos



def write_tab2(C1, C2, L1, L2, N2, NBT, INT):
    """
    Write ENDF-6 **tab2** record.

    Returns
    -------
    `list` of `str`
        list of 66-characters-long ENDF-6 formatted string 
    """
    N1 = len(NBT)
    lines = write_cont(C1, C2, L1, L2, N1, N2)
    lines += write_integer_list(sandy.shared.interwine_lists(NBT, INT))
    return lines



def read_tab1(df, ipos):
    """
    Read ENDF-6 **tab1** record.
    These records are used for one-dimensional tabulated functions such as
    :math:`y(x)`.
    
    The data needed to specify a one-dimensional tabulated function are
    the interpolation tables *nbt* and *int* for each of the *nr* ranges,
    and the *np* tabulated pairs of :math:`x(n)` and :math:`y(n)`.

    Returns
    -------
    `TAB1` : `collections.namedtuple`
        - `C1` : `float`
            first element of string
        - `C2` : `float`
            second element of string
        - `L1` : `int`
            third element of string
        - `L1` : `int`
            fourth element of string
        - `NBT` : `list` of `int`
            upper limits of the interpolation ranges
        - `INT` : `list` of `int`
            interpolation functions
        - `x`,  `y` : `numpy.ndarray`
            tabulated pairs of :math:`x(n)` and :math:`y(n)`
    """
    TAB1 = namedtuple('TAB1', 'C1 C2 L1 L2 NR NP NBT INT x y')
    TAB2, ipos = read_tab2(df, ipos)
    L, ipos = _read_list(df, ipos, TAB2.NZ * 2)
    x = np.array(L[::2], dtype=float)
    y = np.array(L[1::2], dtype=float)
    return TAB1(*list(TAB2), x, y), ipos



def write_tab1(C1, C2, L1, L2, NBT, INT, x, y):
    """
    Write ENDF-6 **tab1** record.

    Returns
    -------
    `list` of `str`
        list of 66-characters-long ENDF-6 formatted string 
    """
    N2 = len(x)
    lines = write_tab2(C1, C2, L1, L2, N2, NBT, INT)
    lines += write_float_list(sandy.shared.interwine_lists(x, y))
    return lines



def _read_list(df, ipos, size):
    iadd = int(np.ceil(size/6))
    vals = df.iloc[ipos:ipos+iadd].values
    tab = list(itertools.chain.from_iterable(vals))[:size]
    ipos += iadd
    return tab, ipos



def read_list(df, ipos):
    LIST = namedtuple('LIST', 'C1 C2 L1 L2 NPL N2 B')
    C, ipos = read_cont(df, ipos)
    L, ipos = _read_list(df, ipos, C.N1)
    return LIST(*list(C), [float(i) for i in L]), ipos



def write_list(C1, C2, L1, L2, N2, B):
    """
    Write ENDF-6 **list** record.

    Outputs:
        - list of string
    """
    NPL = len(B)
    lines = write_cont(C1, C2, L1, L2, NPL, N2)
    lines += write_float_list(B)
    return lines



def write_float(x):
    if abs(x) >= 1e0 and abs(x) < 1e1:
        y = "{:11.8f}".format(x)
    elif abs(x) >= 1E1 and abs(x) < 1E2:
        y = "{:11.7f}".format(x)
    elif abs(x) >= 1E2 and abs(x) < 1E3:
        y = "{:11.6f}".format(x)
    elif abs(x) >= 1E3 and abs(x) < 1E4:
        y = "{:11.5f}".format(x)
    elif abs(x) >= 1E4 and abs(x) < 1E5:
        y = "{:11.4f}".format(x)
    elif abs(x) >= 1E5 and abs(x) < 1E6:
        y = "{:11.3f}".format(x)
    elif abs(x) >= 1E6 and abs(x) < 1E7:
        y = "{:11.2f}".format(x)
    elif abs(x) >= 1E7 and abs(x) < 1E8:
        y = "{:11.1f}".format(x)
    elif abs(x) >= 1E8 and abs(x) < 1E10:
        y = "{:11.0f}".format(x)
    elif x == 0:
        y = "{:11.8f}".format(x)
    elif abs(x) < 1E0 and abs(x) >= 1e-9:
        y = "{:13.6e}".format(x).replace("e-0", "-")
    else:
        y = "{:12.5e}".format(x).replace("e", "")
    return y



def line_numbers(length, istart=1):
    iend = istart + length
    ilines = np.arange(istart, iend)
    return np.where(ilines < 100000, ilines, ilines-99999).tolist()



def write_line(string, mat, mf, mt, iline):
    return line_pattern.format(string, mat, mf, mt, iline)



def write_eol(lines, mat, mf, mt, istart=1):
    """
    Add end-of-line flags MAT, MF, MT and line number to list of strings.

    Returns
    -------
    `str`
        A string that is the sum of the list of strings, i.e. including 
        eol falgs, separated by the newline character `\n`.
    
    Warns
    -----
    This function does not check if the strings are in ENDF-6 format 
    or longer than 66 characters.
    """
    ilines = line_numbers(len(lines), istart=istart)
    return [write_line(string, mat, mf, mt, iline) for string,iline in zip(lines, ilines)]

