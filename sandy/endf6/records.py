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

text_format_w = ff.FortranRecordWriter('(A66,I4,I2,I3,I5)')
cont_format_w = ff.FortranRecordWriter('(2ES11.5E1,4I11,I4,I2,I3,I5)')
list_format_w = ff.FortranRecordWriter('(6ES11.5E1,I4,I2,I3,I5)')
ilist_format_w = ff.FortranRecordWriter('(6I11,I4,I2,I3,I5)')

def read_text(text, ipos):
    TEXT = namedtuple('TEXT', 'HL')
    try:
        text = TEXT(*text_format_r.read(text[ipos]))
        ipos += 1
        return text, ipos
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
        i = 0
        tab = tab[:C.N1*2]
        NBT = tab[::2]
        INT = tab[1::2]
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
        tab = tab[:C.N1*2]
        return LIST(C.C1, C.C2, C.L1, C.L2, C.N1, C.N2, tab), ipos
    except:
        sys.exit("ERROR: cannot read LIST at '{}'".format(text[ipos]))


def write_text(text, ipos):
    TEXT = namedtuple('TEXT', 'HL')
    try:
        text = TEXT(*text_format_r.read(text[ipos]))
        ipos += 1
        return text, ipos
    except:
        sys.exit("ERROR: cannot read TEXT at '{}'".format(text[ipos]))

def write_cont(CONT):
    try:
        return [cont_format_w.write(CONT)]
    except:
        sys.exit("ERROR: cannot write CONT at ")

def write_list(LIST):
    try:
        C, ipos = read_cont(text, ipos)
        i = 0
        tab = []
        while i < C.N1:
            tab.extend(list_format_r.read(text[ipos]))
            ipos += 1
            i += 6
        tab = tab[:C.N1*2]
        return LIST(C.C1, C.C2, C.L1, C.L2, C.N1, C.N2, tab), ipos
    except:
        sys.exit("ERROR: cannot read LIST at '{}'".format(text[ipos]))