# -*- coding: utf-8 -*-
"""
Created on Mon Apr  9 13:46:32 2018

@author: fiorito_l
"""
import os
import sys
import pytest
import warnings
warnings.filterwarnings("ignore", message="numpy.dtype size changed")
warnings.filterwarnings("ignore", message="numpy.ufunc size changed")

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

#@pytest.mark.xs
#def test_H1(tmpdir):
#    from sandy.sampling import sampling
#    from sandy.data_test import __file__ as td
#    td = os.path.dirname(os.path.realpath(td))
#    iargs = [os.path.join(td, r"h1.pendf"),
#             "--errorr-cov", os.path.join(td, r"h1.errorr"),
#             "--outdir", str(tmpdir),
#             "--processes", str(os.cpu_count()),
#             "--eig", "10",
#             "--samples", "100",
#             "--plotdir", os.path.join(str(tmpdir), r"html_files"),
#             "-p"]
#    sampling.run(iargs)
#
#@pytest.mark.nubar
#def test_Cm242(tmpdir):
#    from sandy.sampling import sampling
#    from sandy.data_test import __file__ as td
#    td = os.path.dirname(os.path.realpath(td))
#    iargs = [os.path.join(td, r"cm242.endf"),
#             "--endf6-cov", os.path.join(td, r"cm242.endf"),
#             "--outdir", str(tmpdir),
#             "--processes", str(os.cpu_count()),
#             "--eig", "10",
#             "--samples", "100",
#             "--plotdir", os.path.join(str(tmpdir), r"html_files"),
#             "-p"]
#    sampling.run(iargs)
#
#@pytest.mark.xs
#def test_Fe56_errorr(tmpdir):
#    from sandy.sampling import sampling
#    from sandy.data_test import __file__ as td
#    td = os.path.dirname(os.path.realpath(td))
#    iargs = [os.path.join(td, r"fe56.pendf"),
#             "--errorr-cov", os.path.join(td, r"fe56.errorr"),
#             "--outdir", str(tmpdir),
#             "--processes", str(os.cpu_count()),
#             "--eig", "10",
#             "--samples", "10",]
#    sampling.run(iargs)
#
#@pytest.mark.xs
#@pytest.mark.slow
#def test_U5_errorr(tmpdir):
#    from sandy.sampling import sampling
#    from sandy.data_test import __file__ as td
#    td = os.path.dirname(os.path.realpath(td))
#    iargs = [os.path.join(td, r"u235.pendf"),
#             "--errorr-cov", os.path.join(td, r"u235.errorr"),
#             "--outdir", str(tmpdir),
#             "--processes", str(os.cpu_count()) if os.cpu_count() < 10 else str(10),
#             "--eig", "10",
#             "--samples", "100",
#             "--plotdir", os.path.join(str(tmpdir), r"html_files"),
#             "-p"]
#    sampling.run(iargs)
#
#def runtests():
#    args = [os.path.realpath(__file__),
#            "--basetemp=sandy_tests",]
#    if len(sys.argv) > 1:
#        args += sys.argv[1:]
#    pytest.main(args)
#
#
#if __name__ == "__main__":
#    runtests()