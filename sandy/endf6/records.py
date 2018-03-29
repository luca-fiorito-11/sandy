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

text_format_r = ff.FortranRecordReader('(A66)')
cont_format_r = ff.FortranRecordReader('(2E11.0,4I11)')
list_format_r = ff.FortranRecordReader('(6E11.0)')
ilist_format_r = ff.FortranRecordReader('(6I11)')

ilist_format_w = ff.FortranRecordWriter('(6I11)')

def read_text(text, ipos):
#    TEXT = namedtuple('TEXT', 'HL')
    try:
        TEXT = text_format_r.read(text[ipos])[0]
        ipos += 1
        return TEXT, ipos
    except:
        sys.exit("ERROR: cannot read TEXT at '{}'".format(text[ipos]))

def read_cont(text, ipos):
    CONT = namedtuple('CONT', 'C1 C2 L1 L2 N1 N2')
    try:
        C = CONT(*cont_format_r.read(text[ipos]))
        ipos += 1
        return C, ipos
    except:
        sys.exit("ERROR: cannot read CONT at '{}'".format(text[ipos]))

def read_tab1(text, ipos):
    TAB1 = namedtuple('TAB1', 'C1 C2 L1 L2 NR NP NBT INT x y')
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
        i = 0
        tab = []
        while i < C.N2*2:
            tab.extend(list_format_r.read(text[ipos]))
            ipos += 1
            i += 6
        tab = tab[:C.N2*2]
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
        i = 0
        tab = []
        while i < C.N1:
            tab.extend(list_format_r.read(text[ipos]))
            ipos += 1
            i += 6
        tab = tab[:C.N1]
        return LIST(C.C1, C.C2, C.L1, C.L2, C.N1, C.N2, tab), ipos
    except:
        sys.exit("ERROR: cannot read LIST at '{}'".format(text[ipos]))

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
        TEXT.append(ilist_format_w.write(tab[i:i+6]))
        i += 6
    i = 0
    while i < NP*2:
        L = tab1[i:i+6]
        exps = np.abs(np.floor(log10(np.abs(L))))
        form = "(" + ",".join(list(map(lambda x: 'ES12.6E1' if x<10 else 'ES12.5E2' if x<100 else 'ES12.4E3', exps))) + ")"
        list_format_w = ff.FortranRecordWriter(form)
        TEXT.append(list_format_w.write(L).replace("E",""))
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
        L = B[i:i+6]
        exps = np.abs(np.floor(log10(np.abs(L))))
        form = "(" + ",".join(list(map(lambda x: 'ES12.6E1' if x<10 else 'ES12.5E2' if x<100 else 'ES12.4E3', exps))) + ")"
        list_format_w = ff.FortranRecordWriter(form)
        TEXT.append(list_format_w.write(L).replace("E",""))
        i += 6
    return TEXT