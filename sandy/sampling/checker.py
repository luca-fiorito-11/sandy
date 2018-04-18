# -*- coding: utf-8 -*-
"""
Created on Thu Apr 12 11:29:54 2018

@author: fiorito_l
"""
import pandas as pd
from sandy.endf6 import files as e6
from sandy import settings
import numpy as np
import sys
import os
from os.path import join, dirname, realpath

def run():
    settings.init_checker()
    for inp in os.listdir(settings.args.smpdir):
        print ("check file ", inp)
        DfXs = e6.extract_xs( e6.endf2df(join(settings.args.smpdir,inp), keep_mf=[1,3]) )
        for mat in DfXs.columns.get_level_values("MAT").unique():
            daughters = [ x for x in range(800,850) if x in DfXs[mat].columns]
            if daughters:
                Sum = DfXs[mat][daughters].sum(axis=1)
                if 107 in DfXs[mat].columns:
                    assert np.isclose(Sum.values, DfXs[mat,107].values).all(), "Redundant MT{} differs from the sum of its components".format(107)
                DfXs[mat,107] = Sum
            daughters = [ x for x in range(750,800) if x in DfXs[mat].columns]
            if daughters:
                Sum = DfXs[mat][daughters].sum(axis=1)
                if 106 in DfXs[mat].columns:
                    assert np.isclose(Sum.values, DfXs[mat,106].values).all(), "Redundant MT{} differs from the sum of its components".format(106)
                DfXs[mat,106] = Sum
            daughters = [ x for x in range(700,750) if x in DfXs[mat].columns]
            if daughters:
                Sum = DfXs[mat][daughters].sum(axis=1)
                if 105 in DfXs[mat].columns:
                    assert np.isclose(Sum.values, DfXs[mat,105].values).all(), "Redundant MT{} differs from the sum of its components".format(105)
                DfXs[mat,105] = Sum
            daughters = [ x for x in range(650,700) if x in DfXs[mat].columns]
            if daughters:
                Sum = DfXs[mat][daughters].sum(axis=1)
                if 104 in DfXs[mat].columns:
                    assert np.isclose(Sum.values, DfXs[mat,104].values).all(), "Redundant MT{} differs from the sum of its components".format(104)
                DfXs[mat,104] = Sum
            daughters = [ x for x in range(600,650) if x in DfXs[mat].columns]
            if daughters:
                Sum = DfXs[mat][daughters].sum(axis=1)
                if 103 in DfXs[mat].columns:
                    assert np.isclose(Sum.values, DfXs[mat,103].values).all(), "Redundant MT{} differs from the sum of its components".format(103)
                DfXs[mat,103] = Sum
            daughters = [ x for x in range(102,118) if x in DfXs[mat].columns]
            if daughters:
                Sum = DfXs[mat][daughters].sum(axis=1)
                if 101 in DfXs[mat].columns:
                    assert np.isclose(Sum.values, DfXs[mat,101].values).all(), "Redundant MT{} differs from the sum of its components".format(101)
                DfXs[mat,101] = Sum
            daughters = [ x for x in (19,20,21,38) if x in DfXs[mat].columns]
            if daughters:
                Sum = DfXs[mat][daughters].sum(axis=1)
                if 18 in DfXs[mat].columns:
                    assert np.isclose(Sum.values, DfXs[mat,18].values).all(), "Redundant MT{} differs from the sum of its components".format(18)
                DfXs[mat,18] = Sum
            daughters = [ x for x in (18,101) if x in DfXs[mat].columns]
            if daughters:
                Sum = DfXs[mat][daughters].sum(axis=1)
                if 27 in DfXs[mat].columns:
                    assert np.isclose(Sum.values, DfXs[mat,27].values).all(), "Redundant MT{} differs from the sum of its components".format(27)
                DfXs[mat,27] = Sum
            daughters = [ x for x in range(50,92) if x in DfXs[mat].columns]
            if daughters:
                Sum = DfXs[mat][daughters].sum(axis=1)
                if 4 in DfXs[mat].columns:
                    assert np.isclose(Sum.values, DfXs[mat,4].values).all(), "Redundant MT{} differs from the sum of its components".format(4)
                DfXs[mat,4] = Sum
            daughters = [ x for x in (4,5,11,16,17,*range(22,38),41,42,44,45) if x in DfXs[mat].columns]
            if daughters:
                Sum = DfXs[mat][daughters].sum(axis=1)
                if 3 in DfXs[mat].columns:
                    assert np.isclose(Sum.values, DfXs[mat,3].values).all(), "Redundant MT{} differs from the sum of its components".format(3)
                DfXs[mat,3] = Sum
            daughters = [ x for x in (2,3) if x in DfXs[mat].columns]
            if daughters:
                Sum = DfXs[mat][daughters].sum(axis=1)
                if 1 in DfXs[mat].columns:
                    assert np.isclose(Sum.values, DfXs[mat,1].values).all(), "Redundant MT{} differs from the sum of its components".format(1)
                DfXs[mat,1] = Sum
            daughters = [ x for x in (455,456) if x in DfXs[mat].columns]
            if daughters:
                Sum = DfXs[mat][daughters].sum(axis=1)
                if 452 in DfXs[mat].columns:
                    assert np.isclose(Sum.values, DfXs[mat,452].values).all(), "Redundant MT{} differs from the sum of its components".format(452)
                DfXs[mat,452] = Sum

if __name__ == '__main__':
    from sandy import __file__ as sd
    from sandy.data_test import __file__ as td
    sd = dirname(realpath(sd))
    td = dirname(realpath(td))
    extra_args = [join(sd, r"sampling\cm242-tmpdir")]
    sys.argv = [sys.argv[0]] + extra_args
    run()