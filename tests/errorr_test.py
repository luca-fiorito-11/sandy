# -*- coding: utf-8 -*-
"""
Created on Thu Jul 12 10:00:32 2018

@author: Fiorito_L
"""

import pytest

from sandy.formats.errorr import Errorr
from sandy.data import H1

@pytest.fixture(scope="module")
def testH1():
    tape = Errorr.from_text("\n".join(H1.errorr))
    assert (tape.index.get_level_values("MAT").unique() == 125).all()
    return tape

@pytest.mark.formats
@pytest.mark.endf6
def test_mat_property(testH1):
    sorted(testH1.mat) == [125]

@pytest.mark.formats
@pytest.mark.endf6
def test_mf_property(testH1):
    sorted(testH1.mf) == [1,2,3]

@pytest.mark.formats
@pytest.mark.endf6
def test_mt_property(testH1):
    sorted(testH1.mt) == [1, 2, 102, 151, 451]

@pytest.mark.formats
@pytest.mark.errorr
@pytest.mark.info
def test_read_info(testH1):
    testH1.read_section(125, 1, 451)

@pytest.mark.formats
@pytest.mark.errorr
@pytest.mark.xs
def test_read_xs(testH1):
    testH1.read_section(125, 3, 102)

@pytest.mark.formats
@pytest.mark.errorr
@pytest.mark.cov
@pytest.mark.xs
def test_read_xs_cov(testH1):
    testH1.read_section(125, 33, 102)

@pytest.mark.formats
@pytest.mark.errorr
@pytest.mark.xs
def test_extract_xs(testH1):
    testH1.get_xs(listmat=[125], listmt=[1,2])

@pytest.mark.formats
@pytest.mark.errorr
@pytest.mark.cov
@pytest.mark.xs
def test_extract_xs_cov(testH1):
    testH1.get_xs_cov()