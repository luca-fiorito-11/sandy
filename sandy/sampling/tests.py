# -*- coding: utf-8 -*-
"""
Created on Mon Apr  9 13:46:32 2018

@author: fiorito_l
"""
from os.path import join, dirname, realpath
import os
import sys
import pytest

#def test_H1():
#    from sandy.sampling import sampling
#    from sandy.data_test import __file__ as td
#    from sandy import __file__ as sd
#    sd = dirname(realpath(sd))
#    td = dirname(realpath(td))
#    extra_args = [join(td, r"h1.endf"),
##                 "--pendf", join(td, r"h1.pendf"),
#                 "--outdir", r"h1-tmpdir",
#                 "--njoy", join(sd, r"njoy2012_50.exe"),
#                 "--eig", "10",
#                 "--samples", "10",
#                 "--processes", "1",
##                 "-mf", "33",
#                 "-e", "1e-5",
#                 "-e", "5e-5",
#                 "-e", "1e-4",
#                 "-e", "5e-4",
#                 "-e", "1e-3",
#                 "-e", "5e-3",
#                 "-e", "1e-2",
#                 "-e", "5e-2",
#                 "-e", "1e-1",
#                 "-e", "5e-1",
#                 "-e", "1e0",
#                 "-e", "5e0",
#                 "-e", "1e1",
#                 "-e", "5e1",]
#    sys.argv = [sys.argv[0]] + extra_args
#    sampling.run()
#
#def test_U5():
#    from sandy.sampling import sampling
#    from sandy.data_test import __file__ as td
#    from sandy import __file__ as sd
#    sd = dirname(realpath(sd))
#    td = dirname(realpath(td))
#    extra_args = [join(td, r"u235.endf"),
#                 "--pendf", join(td, r"u235.pendf"),
#                 "--outdir", r"u5-33-tmpdir",
##                 "--njoy", join(sd, r"njoy2012_50.exe"),
#                 "--eig", "10",
#                 "--samples", "3",
#                 "--processes", "1",
#                 "-mf", "33",
#                 "-e", "1e-5",
#                 "-e", "5e-5",
#                 "-e", "1e-4",
#                 "-e", "5e-4",
#                 "-e", "1e-3",
#                 "-e", "5e-3",
#                 "-e", "1e-2",
#                 "-e", "5e-2",
#                 "-e", "1e-1",
#                 "-e", "5e-1",
#                 "-e", "1e0",
#                 "-e", "5e0",
#                 "-e", "1e1",
#                 "-e", "5e1",]
#    sys.argv = [sys.argv[0]] + extra_args
#    sampling.run()
#
#def test_Pu9():
#    from sandy.sampling import sampling
#    from sandy.data_test import __file__ as td
#    from sandy import __file__ as sd
#    sd = dirname(realpath(sd))
#    td = dirname(realpath(td))
#    extra_args = [join(td, r"pu239.endf"),
#                 "--pendf", join(td, r"pu239.pendf"),
#                 "--outdir", r"pu9-tmpdir",
##                 "--njoy", join(sd, r"njoy2012_50.exe"),
#                 "--eig", "10",
#                 "--samples", "100",
#                 "--processes", "1",
#                 "-mf", "31",
#                 "-e", "1e-5",
#                 "-e", "5e-5",
#                 "-e", "1e-4",
#                 "-e", "5e-4",
#                 "-e", "1e-3",
#                 "-e", "5e-3",
#                 "-e", "1e-2",
#                 "-e", "5e-2",
#                 "-e", "1e-1",
#                 "-e", "5e-1",
#                 "-e", "1e0",
#                 "-e", "5e0",
#                 "-e", "1e1",
#                 "-e", "5e1",]
#    sys.argv = [sys.argv[0]] + extra_args
#    sampling.run()
#

def test_Cm242(tmpdir, capsys):
    from sandy.sampling import sampling
    from sandy.data_test import __file__ as td
    td = dirname(realpath(td))
    iargs = [join(td, r"cm242.endf"),
             "--endf6-cov", join(td, r"cm242.endf"),
             "--outdir", str(tmpdir),
             "--processes", str(os.cpu_count()),
             "--eig", "10",
             "--samples", "300",]
    sampling.run(iargs)
    captured = capsys.readouterr()
    with open(join(str(tmpdir), "sandy.stdout"), 'w') as f: f.write(captured.out)
    with open(join(str(tmpdir), "sandy.stderr"), 'w') as f: f.write(captured.err)

def test_Fe56_errorr(tmpdir, capsys):
    from sandy.sampling import sampling
    from sandy.data_test import __file__ as td
    td = dirname(realpath(td))
    iargs = [join(td, r"fe56.pendf"),
             "--errorr-cov", join(td, r"fe56.errorr"),
             "--outdir", str(tmpdir),
             "--processes", str(os.cpu_count()),
             "--eig", "10",
             "--samples", "300",]
    sampling.run(iargs)
    captured = capsys.readouterr()
    with open(join(str(tmpdir), "sandy.stdout"), 'w') as f: f.write(captured.out)
    with open(join(str(tmpdir), "sandy.stderr"), 'w') as f: f.write(captured.err)

#def test():
#    from sandy.sampling import sampling
#    from sandy.data_test import __file__ as td
#    from sandy import __file__ as sd
#    sd = dirname(realpath(sd))
#    td = dirname(realpath(td))
#    iargs = [join(td, r"cm242.endf"),
#             "--endf6-cov", join(td, r"cm242.endf"),
#             "--outdir", join(sd, r"cm242-tmpdir"),
#             "--eig", "10",
#             "--samples", "10",]
#    sampling.run(iargs)
#
#test()
#sys.exit()

def run():
    args = [realpath(__file__),
            "--basetemp=sandy_tests",]
    if len(sys.argv) > 1:
        args += sys.argv[1:]
    pytest.main(args)

if __name__ == "__main__":
    run()