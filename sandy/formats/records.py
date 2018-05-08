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
        i = 0
        tab = []
        while i < C.N1*2:
            tab.extend(ilist_format_r.read(text[ipos]))
            ipos += 1
            i += 6
        tab = tab[:C.N1*2]
        NBT = tab[::2]
        INT = tab[1::2]
        iadd = int(np.ceil(C.N2*2/6))
        # 3rd method, 3.14 secs
        string = "".join([ text[ipos+i][:66] for i in range(int(iadd))])[:C.N2*2*11]
        # mapping here makes it 1 sec slower
#        string = "".join(map(lambda x : x[:66], text[ipos:ipos+iadd]))[:C.N2*2*11]
        tab = list(map(read_float, re.findall(".{11}", string)))
        # 6th method, 7.34
#        tab = np.genfromtxt(
#                map(lambda x: x.encode(), text[ipos:ipos+iadd]),
#                delimiter=[11]*6,
#                converters=dict(zip(range(6), [read_float]*6))
#                ).flatten()[:C.N2*2]
        # 5th method, 13.6 secs
#        string = re.sub("(.{11})", "\\1 ", string) # add spaces
#        tab = np.array(re.sub('([0-9])([+-])', '\\1E\\2', string).split(), dtype=float)
#        B = re.sub('([0-9]-)', '\\1E-', A)
        # 4th method, 9 secs
#        tab = [ read_float(x) for x in re.findall(".{11}", string)]
        # 2nd method, 24 secs
#        list_format_r = ff.FortranRecordReader('({}E11.0)'.format(C.N2*2))
#        tab = list_format_r.read(string)
        # 1st method, 38 secs
#        i = 0
#        tab = []
#        while i < C.N2*2:
#            tab.extend(list_format_r.read(text[ipos]))
#            ipos += 1
#            i += 6
#        tab = tab[:C.N2*2]
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
        TEXT.append("".join(map(lambda x:"{:>11}".format(x), tab[i:i+6])))
#        TEXT.append(ilist_format_w.write(tab[i:i+6]))
        i += 6
    i = 0
    while i < NP*2:
        TEXT.append("".join(map(lambda x:"{:>12.5E}".format(x), tab1[i:i+6])).replace("E",""))
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