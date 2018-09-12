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

@pytest.fixture(scope="module")
def Pu9endf():
    return Endf6.from_text("\n".join(Pu9.endf6))

@pytest.mark.formats
@pytest.mark.replace
def test_replace_endf6(Pu9endf):
    file = os.path.join(Pu9.__path__[0], "pu239.pendf")
    Pu9pendf = Endf6.from_file(file)
    df = Pu9endf.add_sections(file, {2: [152], 10: "all"})
    assert (df.loc[9437,10,42] == Pu9pendf.loc[9437,10,42]).all()
    assert not (df.loc[9437,10,42] == Pu9endf.loc[9437,10,42]).all()
    assert (9437,2,152) in df.index
    assert (9437,2,152) not in Pu9endf.index
    assert (df.loc[9437,2,152] == Pu9pendf.loc[9437,2,152]).all()
    assert not (df.loc[9437,3,2] == Pu9pendf.loc[9437,3,2]).all()
    assert (df.loc[9437,3,2] == Pu9endf.loc[9437,3,2]).all()
    
@pytest.mark.formats
@pytest.mark.basefile
def test_grep_mt(Pu9endf):
    text = "\n".join(Pu9.pendf)
    A = BaseFile.from_text(text, listmt=[45651])
    pytest.set_trace()
    df = Pu9endf.add_sections(file, {2: [152], 10: "all"})
    assert (df.loc[9437,10,42] == Pu9pendf.loc[9437,10,42]).all()
    assert not (df.loc[9437,10,42] == Pu9endf.loc[9437,10,42]).all()
    assert (9437,2,152) in df.index
    assert (9437,2,152) not in Pu9endf.index
    assert (df.loc[9437,2,152] == Pu9pendf.loc[9437,2,152]).all()
    assert not (df.loc[9437,3,2] == Pu9pendf.loc[9437,3,2]).all()
    assert (df.loc[9437,3,2] == Pu9endf.loc[9437,3,2]).all()
    