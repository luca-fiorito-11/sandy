# -*- coding: utf-8 -*-
"""
Created on Thu Jun 14 10:22:58 2018

@author: fiorito_l
"""
import sys
import numpy as np
from collections import namedtuple
import pdb

import fortranformat as ff
import rwf
from ..settings import SandyError

__author__ = "Luca Fiorito"
__all__ = ["read_control", "read_cont", "write_cont", "read_text", "read_list", 
           "write_list", "read_tab2", "write_tab2", "read_tab1", "write_tab1"]

def wreal(x):
    line = None
    if abs(x) > 1e-1 and abs(x) < 1e1:
        line = ff.FortranRecordWriter('(F11.8)')
        line = line.write([x])
    elif abs(x) >= 1e1 and abs(x) < 1e2:
        line = ff.FortranRecordWriter('(F11.7)')
        line = line.write([x])
    elif abs(x) >= 1e2 and abs(x) < 1e3:
        line = ff.FortranRecordWriter('(F11.6)')
        line = line.write([x])
    elif abs(x) >= 1e3 and abs(x) < 1e4:
        line = ff.FortranRecordWriter('(F11.5)')
        line = line.write([x])
    elif abs(x) >= 1e4 and abs(x) < 1e5:
        line = ff.FortranRecordWriter('(F11.4)')
        line = line.write([x])
    elif abs(x) >= 1e5 and abs(x) < 1e6:
        line = ff.FortranRecordWriter('(F11.3)')
        line = line.write([x])
    elif abs(x) >= 1e6 and abs(x) < 1e7:
        line = ff.FortranRecordWriter('(F11.2)')
        line = line.write([x])
    elif abs(x) >= 1e7 and abs(x) < 1e8:
        line = ff.FortranRecordWriter('(F11.1)')
        line = line.write([x])
    elif abs(x) >= 1e8 and abs(x) < 1e9:
        line = ff.FortranRecordWriter('(F11.0)')
        line = line.write([x])
    else:
        line = ff.FortranRecordWriter('(1PE13.6)')
        line = line.write([x])
        
        if line[11] == '0':
            line = line[1:9] + line[10] + line[12]
            line2 = ff.FortranRecordWriter('A11')
            line = line2.write([line])
        else:
            line = line[0:8] + line[10:13]
            line2 = ff.FortranRecordWriter('A11')
            line = line2.write([line])

    return line

def read_control(string):
    mat = np.array(0, dtype=int)
    mf = np.array(0, dtype=int)
    mt = np.array(0, dtype=int)
    ns = np.array(0, dtype=int)
    io_status = np.array(0, dtype=int)
    line = ff.FortranRecordReader('(I4,I2,I3,I5)')
    line = line.read(string[66:])
    mat = line[0]
    mf = line[1]
    mt = line[2]
    ns =  line[3]
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
    line = ff.FortranRecordReader('(2E11.0,4I11)')
    line = line.read(string)
    c1 = line[0]
    c2 = line[1]
    l1 = line[2]
    l2 = line[3]
    n1 = line[4]
    n2 = line[5]
    ipos += 1
    return CONT(float(c1), float(c2), int(l1), int(l2), int(n1), int(n2)), ipos

def write_cont(C1, C2, L1, L2, N1, N2):
    """
    Write ENDF-6 **cont** record.

    Outputs:
        - list of string
    """
    C1 = np.array(C1, dtype=float)
    C2 = np.array(C2, dtype=float)
    L1 = np.array(L1, dtype=int)
    L2 = np.array(L2, dtype=int)
    N1 = np.array(N1, dtype=int)
    N2 = np.array(N2, dtype=int)
    string = np.array("*"*67)
    c1_str = wreal(C1)
    c2_str = wreal(C2)
    line = ff.FortranRecordWriter('(2A11,4I11)')
    string = line.write([c1_str,c2_str,L1,L2,N1,N2])
    string = string[:66] # byte string coming from f2py must be converted
    return [string]

def read_text(text, ipos):
    string = text[ipos]
    try:
        ipos += 1
        return text[ipos-1][:66], ipos
    except:
        sys.exit("ERROR : line '{}' is not in TEXT format".format(string))

def read_ilist(string):
    array = np.zeros(6, dtype=np.int64)
    line = ff.FortranRecordReader('(6I11)')
    line = line.read(string)
    for idx,value in enumerate(line):
        array[idx] = line[idx]
    return array.tolist()

def write_ilist(b):
    length = len(b)*11
    string = np.array("*"*(length + 1))
    #rwf.wilist(string, b, length)
    total = ""
    for i in range(0,len(b)):
        line = ff.FortranRecordWriter('(I11)')
        line = line.write([b[i]])
        total += line
    string = total
    return [string[0+i:66+i] + ' '*(66-len(string[0+i:66+i])) for i in range(0, len(string), 66)]

def read_dlist(string):
    array = np.zeros(6, dtype=float)
    line = ff.FortranRecordReader('(6E11.0)')
    line = line.read(string)
    for idx,value in enumerate(line):
        array[idx] = line[idx]
    return array.tolist()

def write_dlist(b):
    length = len(b)*11
    string = np.array("*"*(length + 1))
    #rwf.wlist(string, b, length)
    total = ""
    for i in range(0,len(b)):
        #line = ff.FortranRecordWriter('(F11.8)')
        #line.write([b[i]])
        line = wreal(b[i])
        string = line
        total += string
    #string = str(total, 'utf-8')[:length]
    string = total
    return [string[0+i:66+i] + ' '*(66-len(string[0+i:66+i])) for i in range(0, len(string), 66)]

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

def write_tab2(C1, C2, L1, L2, N2, NBT, INT):
    """
    Write ENDF-6 **tab2** record.

    Outputs:
        - list of string
    """
    N1 = len(NBT)
    text = write_cont(C1, C2, L1, L2, N1, N2)
    b = [ z for item in zip(NBT, INT) for z in item ]
    text += write_ilist(b)
    return text

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

def write_tab1(C1, C2, L1, L2, NBT, INT, x, y):
    """
    Write ENDF-6 **tab1** record.

    Outputs:
        - list of string
    """
    N2 = len(x)
    text = write_tab2(C1, C2, L1, L2, N2, NBT, INT)
    b = [ z for item in zip(x, y) for z in item ]
    text += write_dlist(b)
    return text

def read_list(text, ipos):
    import itertools
    LIST = namedtuple('LIST', 'C1 C2 L1 L2 NPL N2 B')
    C, ipos = read_cont(text, ipos)
    iadd = int(np.ceil(C.N1/6))
    tab = map(lambda x:read_dlist(x) ,text[ipos:ipos+iadd])
    tab = list(itertools.chain.from_iterable(tab))[:C.N1]
    ipos += iadd
    return LIST(*list(C), tab), ipos

def write_list(C1, C2, L1, L2, N2, B):
    """
    Write ENDF-6 **list** record.

    Outputs:
        - list of string
    """
    NPL = len(B)
    text = write_cont(C1, C2, L1, L2, NPL, N2)
    text += write_dlist(B)
    return text
