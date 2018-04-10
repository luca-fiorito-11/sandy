# -*- coding: utf-8 -*-
"""
Created on Mon Apr  9 13:46:32 2018

@author: fiorito_l
"""
from os.path import join, dirname, realpath
import sys

def test_H1():
    from sandy.sampling import sampling
    from sandy.data_test import __file__ as td
    from sandy import __file__ as sd
    sd = dirname(realpath(sd))
    td = dirname(realpath(td))
    sys.argv.extend([join(td, r"1-H-1g.jeff33"),
                     "--pendf", join(td, r"1-H-1g.jeff33.pendf"),
                     "--outdir", r"tmpdir",
                     "--njoy", join(sd, r"njoy2012_50.exe"),
                     "--eig", "10",
                     "--samples", "10",
                     "--processes", "1",
                     "-mf", "33",
                     "-e", "1e-5",
                     "-e", "5e-5",
                     "-e", "1e-4",
                     "-e", "5e-4",
                     "-e", "1e-3",
                     "-e", "5e-3",
                     "-e", "1e-2",
                     "-e", "5e-2",
                     "-e", "1e-1",
                     "-e", "5e-1",
                     "-e", "1e0",
                     "-e", "5e0",
                     "-e", "1e1",
                     "-e", "5e1",])
    sampling.run()

def test_U5():
    from sandy.sampling import sampling
    from sandy.data_test import __file__ as td
    from sandy import __file__ as sd
    sd = dirname(realpath(sd))
    td = dirname(realpath(td))
    sys.argv.extend([join(td, r"92-U-235g.jeff33"),
                     "--pendf", join(td, r"92-U-235g.jeff33.pendf"),
                     "--outdir", r"tmpdir",
                     "--njoy", join(sd, r"njoy2012_50.exe"),
                     "--eig", "10",
                     "--samples", "10",
                     "--processes", "1",
                     "-e", "1e-5",
                     "-e", "5e-5",
                     "-e", "1e-4",
                     "-e", "5e-4",
                     "-e", "1e-3",
                     "-e", "5e-3",
                     "-e", "1e-2",
                     "-e", "5e-2",
                     "-e", "1e-1",
                     "-e", "5e-1",
                     "-e", "1e0",
                     "-e", "5e0",
                     "-e", "1e1",
                     "-e", "5e1",])
    sampling.run()


test_U5()