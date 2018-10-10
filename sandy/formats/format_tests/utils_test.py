# -*- coding: utf-8 -*-
"""
Created on Thu Sep  6 10:08:49 2018

@author: Fiorito_L
"""
import pytest
import os

from .. import Endf6, BaseFile
from ...data import Pu9

__author__ = "Luca Fiorito"

#@pytest.fixture(scope="module")
#def Pu9endf():
#    return Endf6.from_text("\n".join(Pu9.endf6))

#@pytest.mark.formats
#@pytest.mark.basefile
#def test_grep_mt(Pu9endf):
#    text = "\n".join(Pu9.pendf)
#    A = BaseFile.from_text(text, listmt=[45651])
#    pytest.set_trace()    