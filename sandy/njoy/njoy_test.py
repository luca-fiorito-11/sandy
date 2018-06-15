# -*- coding: utf-8 -*-
"""
Created on Wed Jun 13 14:59:53 2018

@author: fiorito_l
"""

#import os, sys, pytest
#import warnings
#
#@pytest.mark.njoy
#def test_libProcessing(tmpdir):
#    from sandy.njoy import parser
#    iargs = [os.path.join(td, r"u235.pendf"),
#             "--errorr-cov", os.path.join(td, r"u235.errorr"),
#             "--outdir", str(tmpdir),
#             "--processes", str(os.cpu_count()) if os.cpu_count() < 10 else str(10),
#             "--eig", "10",
#             "--samples", "100",
#             "--plotdir", os.path.join(str(tmpdir), r"html_files"),
#             "-p"]
#    from sandy import settings
#    inputs = settings.init_njoy()
#    evalLib.from_file(inputs["inputfile"]).run_njoy(inputs)
#
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
#lib = PyNjoy(iwt=-98, capsys=False)
#file = evalFile("1-H-3g.jeff33")
#file.pendfFile = lib.pendf(**file.__dict__)
#lib.errorr(**file.__dict__)
#sys.exit()

#from sandy.data_test import __file__ as td
#ppp=1