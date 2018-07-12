# -*- coding: utf-8 -*-
"""
Created on Wed Jul 11 14:32:43 2018

@author: Fiorito_L
"""
import pytest
import os

from ..njoy import run
from .. import data

@pytest.mark.njoy
@pytest.mark.slow
def test_H1(tmpdir):
    iargs = [
            os.path.join(data.__path__[0], r"h1.endf"),
            "-p",
            "--temps", "300", "600",
#            *"--sig0 1E10 1E4 1E3 1E2 1E1 1E0".split(),
#            "--kerma",
#            "--free-gas",
#            "--ptable",
#            "--gaspr",
#            "--ign", "19",
#            "--iwt", "2",
#            *"--suffixes .03 .06".split(),
            ]
    run(iargs)