# -*- coding: utf-8 -*-
"""
Created on Wed Feb 14 16:04:33 2018

@author: fiorito_l
"""

import sys
import logging
import numpy as np
import fortranformat as ff
from collections import namedtuple
from sandy.tests import TimeDecorator
import re
list_format_r = ff.FortranRecordReader('(6E11.0)')
ilist_format_r = ff.FortranRecordReader('(6I11)')
ilist_format_w = ff.FortranRecordWriter('(6I11)')

def read_text(text, ipos):
    try:
        out = "".join(["{:>11}".format(x) for x in text.iloc[ipos][:6].fillna('')])
        ipos += 1
        return out, ipos
    except:
        mat=text.MAT.iloc[0]
        mf=text.MF.iloc[0]
        mt=text.MT.iloc[0]
        sys.exit("ERROR: cannot read TEXT at line {} for MAT{}/MF{}/MT{}".format(ipos,mat,mf,mt))

def read_cont(text, ipos):
    CONT = namedtuple('CONT', 'C1 C2 L1 L2 N1 N2')
    try:
        out = text.iloc[ipos][:2].tolist() + text.iloc[ipos][2:6].astype(int).tolist()
        ipos += 1
        return CONT(*out), ipos
    except:
        mat=text.MAT.iloc[0]
        mf=text.MF.iloc[0]
        mt=text.MT.iloc[0]
        sys.exit("ERROR: cannot read CONT at line {} for MAT{}/MF{}/MT{}".format(ipos,mat,mf,mt))

def read_float(item):
    # stripping takes 2 extra secs
#    stripped = item.upper().replace('D', 'E').strip()
    # applying float to the array rather than to each item makes us gain 0.5 secs
    # check for presence of +-" increases time by 2 secs
#    item = item.decode()
#    return float(item[0] + item[1:].replace('+', 'E+').replace('-', 'E-'))
    return item[0] + item[1:].replace('+', 'E+').replace('-', 'E-')

#@TimeDecorator
def read_tab1(text, ipos):
    TAB1 = namedtuple('TAB1', 'C1 C2 L1 L2 NR NP NBT INT x y')
    try:
        C, ipos = read_cont(text, ipos)
        iadd = int(np.ceil(C.N1*2/6))
        tab = text.iloc[ipos:ipos+iadd,:6].values.flatten()[:C.N1*2].astype(int).tolist()
        ipos += iadd
#        i = 0
#        tab = []
#        while i < C.N1*2:
#            out = text.iloc[ipos][:6].fillna(0).astype(int).values
#            tab.extend(out)
#            ipos += 1
#            i += 6
#        tab = tab[:C.N1*2]
        NBT = tab[::2]
        INT = tab[1::2]
        iadd = int(np.ceil(C.N2*2/6))
        tab = text.iloc[ipos:ipos+iadd,:6].values.flatten()[:C.N2*2]
        ipos += iadd
        x = np.array(tab[::2], dtype=float)
        y = np.array(tab[1::2], dtype=float)
        return TAB1(C.C1, C.C2, C.L1, C.L2, C.N1, C.N2, NBT, INT, x, y), ipos
    except:
        sys.exit("ERROR: cannot read TAB1 at '{}'".format(text[ipos]))

def read_tab2(text, ipos):
    TAB2 = namedtuple('TAB2', 'C1 C2 L1 L2 NR NZ NBT INT')
    try:
        C, ipos = read_cont(text, ipos)
        i = 0
        tab = []
        while i < C.N1*2:
            tab.extend(ilist_format_r.read(text[ipos]))
            ipos += 1
            i += 6
        tab = tab[:C.N1*2]
        NBT = tab[::2]
        INT = tab[1::2]
        return TAB2(C.C1, C.C2, C.L1, C.L2, C.N1, C.N2, NBT, INT), ipos
    except:
        sys.exit("ERROR: cannot read TAB2 at '{}'".format(text[ipos]))

def read_list(text, ipos):
    LIST = namedtuple('LIST', 'C1 C2 L1 L2 NPL N2 B')
    try:
        C, ipos = read_cont(text, ipos)
        iadd = int(np.ceil(C.N1/6))
        tab = text.iloc[ipos:ipos+iadd,:6].values.flatten()[:C.N1].tolist()
        ipos += iadd
        return LIST(*list(C), tab), ipos
    except:
        sys.exit("ERROR: cannot read LIST at '{}'".format(text[ipos]))

def write_float(x):
    mantissa, exp = "{:>12.6E}".format(x).split('E')
    exp = int(exp)
    sign = np.sign(exp)
    if exp == -1:
        return "{}{:>1.8f}".format("" if np.sign(float(mantissa)) < 0 else " ", x)
    elif exp == 0:
        return "{}{:>1.8f}".format("" if np.sign(float(mantissa)) < 0 else " ", x)
    elif exp == 1:
        return "{}{:>2.7f}".format("" if np.sign(float(mantissa)) < 0 else " ", x)
    elif exp == 2:
        return "{}{:>3.6f}".format("" if np.sign(float(mantissa)) < 0 else " ", x)
    elif exp == 3:
        return "{}{:>4.5f}".format("" if np.sign(float(mantissa)) < 0 else " ", x)
    elif exp == 4:
        return "{}{:>5.4f}".format("" if np.sign(float(mantissa)) < 0 else " ", x)
    elif exp == 5:
        return "{}{:>6.3f}".format("" if np.sign(float(mantissa)) < 0 else " ", x)
    elif exp == 6:
        return "{}{:>7.2f}".format("" if np.sign(float(mantissa)) < 0 else " ", x)
    elif exp == 7:
        return "{}{:>8.1f}".format("" if np.sign(float(mantissa)) < 0 else " ", x)
    elif exp == 8:
        return "{}{:>9.0f}.".format("" if np.sign(float(mantissa)) < 0 else " ", x)
    elif exp == 9:
        return "{}{:>10.0f}".format("" if np.sign(float(mantissa)) < 0 else " ", x)
    elif exp >=10 or exp <= -10:
        return "{:>8}{}{:>}".format(mantissa[:-1], "-" if sign < 0 else "+", abs(exp))
    else:
        return "{:>9}{}{:>}".format(mantissa, "-" if sign < 0 else "+", abs(exp))

#print(write_float(1e-11))
#print(write_float(1e-10))
#print(write_float(1e-9))
#print(write_float(1e-8))
#print(write_float(1e-7))
#print(write_float(1e-6))
#print(write_float(1e-5))
#print(write_float(1e-4))
#print(write_float(1e-3))
#print(write_float(1e-2))
#print(write_float(1e-1))
#print(write_float(0))
#print(write_float(1e0))
#print(write_float(1e1))
#print(write_float(1e2))
#print(write_float(1e3))
#print(write_float(1e5))
#print(write_float(1e6))
#print(write_float(1e7))
#print(write_float(1e8))
#print(write_float(1e9))
#print(write_float(1e10))
#print(write_float(1e11))

def write_cont(C1, C2, L1, L2, N1, N2):
    from sandy.functions import log10
    exps = np.abs(np.floor(log10(np.abs([C1,C2]))))
    form = "(" + ",".join(list(map(lambda x: 'ES12.6E1' if x<10 else 'ES12.5E2' if x<100 else 'ES12.4E3', exps))) + ",4I11)"
    cont_format_w = ff.FortranRecordWriter(form)
    return [cont_format_w.write((C1, C2, L1, L2, N1, N2)).replace("E","")]

def write_tab1(C1, C2, L1, L2, NBT, INT, x, y):
    from sandy.functions import log10
    tab = [item for pair in zip(NBT, INT) for item in pair]
    tab1 = [item for pair in zip(x, y) for item in pair]
    NR = len(NBT)
    NP = len(y)
    TEXT = write_cont(C1, C2, L1, L2, NR, NP)
    i = 0
    while i < NR*2:
        TEXT.append("".join(map(lambda x:"{:>11}".format(x), tab[i:i+6])))
#        TEXT.append(ilist_format_w.write(tab[i:i+6]))
        i += 6
    i = 0
    while i < NP*2:
        TEXT.append("".join(map(write_float, tab1[i:i+6])))
#        TEXT.append("".join(map(lambda x:"{:>12.5E}".format(x), tab1[i:i+6])).replace("E",""))
#        L = tab1[i:i+6]
#        exps = np.abs(np.floor(log10(np.abs(L))))
#        form = "(" + ",".join(list(map(lambda x: 'ES12.6E1' if x<10 else 'ES12.5E2' if x<100 else 'ES12.4E3', exps))) + ")"
#        list_format_w = ff.FortranRecordWriter(form)
#        TEXT.append(list_format_w.write(L).replace("E",""))
        i += 6
    return TEXT

def write_tab2(C1, C2, L1, L2, NZ, NBT, INT):
    tab = [item for pair in zip(NBT, INT) for item in pair]
    NR = len(NBT)
    TEXT = write_cont(C1, C2, L1, L2, NR, NZ)
    i = 0
    while i < NR*2:
        TEXT.append(ilist_format_w.write(tab[i:i+6]))
        i += 6
    return TEXT

def write_list(C1, C2, L1, L2, N2, B):
    from sandy.functions import log10
    NPL = len(B)
    TEXT = write_cont(C1, C2, L1, L2, NPL, N2)
    i = 0
    while i < NPL:
        TEXT.append("".join(map(lambda x:"{:>12.5E}".format(x), B[i:i+6])).replace("E",""))
#        L = B[i:i+6]
#        exps = np.abs(np.floor(log10(np.abs(L))))
#        form = "(" + ",".join(list(map(lambda x: 'ES12.6E1' if x<10 else 'ES12.5E2' if x<100 else 'ES12.4E3', exps))) + ")"
#        list_format_w = ff.FortranRecordWriter(form)
#        TEXT.append(list_format_w.write(L).replace("E",""))
        i += 6
    return TEXT