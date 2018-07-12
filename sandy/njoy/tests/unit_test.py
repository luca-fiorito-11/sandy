# -*- coding: utf-8 -*-
"""
Created on Wed Jul 11 17:54:04 2018

@author: Fiorito_L
"""

import pytest
import os

from .. import Njoy
from .. import data



@pytest.mark.njoy
@pytest.mark.aaa
def test_pendf(tmpdir):
    nj = Njoy()
    file = os.path.join(data.__path__[0], r"h1.endf")
    nj.get_pendf(tape=file, mat=125, TAG="H1")
    