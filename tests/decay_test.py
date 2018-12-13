# -*- coding: utf-8 -*-
"""
Created on Tue Nov 20 14:02:00 2018

@author: Luca Fiorito
"""

import pytest
import os

import numpy as np

from sandy.data import RDD
import sandy.decay
from sandy.decay import DecayChains, BMatrix, QMatrix
from sandy.formats.endf6 import Endf6

__author__ = "Luca Fiorito"

@pytest.fixture(scope="module")
def testRDD():
    tape = Endf6.from_text("\n".join(RDD.endf6))
    assert tape.index.get_level_values("MAT").unique().size == 3852
    return tape

@pytest.fixture(scope="module")
def DecayChains_small():
    tape = Endf6.from_text("\n".join(RDD.endf6), listmat=[471, 498, 528])
    return DecayChains.from_endf6(tape)

#@pytest.mark.rdd
#@pytest.mark.qmatrix
#def test_qmatrix():
#    Q = sandy.get_jeff_qmatrix()
#    assert (np.isclose(np.diag(Q.values), 1)).all()

@pytest.mark.formats
@pytest.mark.endf6
@pytest.mark.rdd
def test_read_rdd(testRDD):
    H1 = testRDD.read_section(2, 8, 457)
    assert H1["ZA"] == 1001
    assert H1["LIS"] == 0
    assert H1["LISO"] == 0
    assert not H1["DK"]
    U5 = testRDD.read_section(3542, 8, 457)
    assert U5["ZA"] == 92235
    assert U5["LIS"] == 0
    assert U5["LISO"] == 0
    assert len(U5["DK"]) == 2
    assert U5["DK"][0]["RTYP"] == 4.0
    assert U5["DK"][1]["RTYP"] == 6.0

@pytest.mark.rdd
def test_decay_chains(DecayChains_small):
    DC = DecayChains_small.copy()
    assert DC.loc[(DC.daughter == 250560) & (DC.parent == 240560)]["yield"].iloc[0].sum() == 1
    assert DC.loc[(DC.daughter == 240560) & (DC.parent == 240560)]["yield"].iloc[0].sum() == -1
    assert DC.loc[(DC.daughter == 260560) & (DC.parent == 250560)]["yield"].iloc[0].sum() == 1
    assert DC.loc[(DC.daughter == 250560) & (DC.parent == 250560)]["yield"].iloc[0].sum() == -1
    assert DC.loc[(DC.daughter == 260560) & (DC.parent == 260560)]["yield"].iloc[0].sum() == 0

@pytest.mark.rdd
def test_qmatrix(DecayChains_small):
    Q = DecayChains_small.get_qmatrix()
    assert np.isclose(Q.loc[240560,240560], 1)
    assert np.isclose(Q.loc[250560,250560], 1)
    assert np.isclose(Q.loc[260560,260560], 1)
    assert np.isclose(Q.loc[250560,240560], 1)
    assert np.isclose(Q.loc[260560,240560], 1)
    assert np.isclose(Q.loc[260560,250560], 1)
    assert np.isclose(Q.loc[240560,250560], 0)
    assert np.isclose(Q.loc[240560,260560], 0)
    assert np.isclose(Q.loc[250560,260560], 0)

@pytest.mark.rdd
def test_bmatrix(DecayChains_small):
    B = DecayChains_small.get_bmatrix()
    assert (B.values == np.array([[0,0,0],[1,0,0],[0,1,0]])).all()

@pytest.mark.rdd
def test_transition_matrix(DecayChains_small):
    """Test on a toy problem that the transition matrix contains 
    the correct values.
    """
    T = DecayChains_small.get_transition_matrix()
    assert np.isclose(T.loc[240560,240560], -0.0019448574089785222)
    assert np.isclose(T.loc[250560,240560], 0.0019448574089785222)
    assert np.isclose(T.loc[250560,250560], -7.46635393516041e-05)
    assert np.isclose(T.loc[260560,250560], 7.46635393516041e-05)

@pytest.mark.rdd
def test_bmatrix_methods():
    """Test that all accessible methods to bmatrix are working and are equivalent.
    """
    tape = Endf6.from_text("\n".join(RDD.endf6))
    DC = DecayChains.from_endf6(tape)
    B1 = DC.get_bmatrix()
    B2 = BMatrix.from_file(os.path.join(RDD.__path__[0], "RDD.jeff33"))
    B3 = sandy.decay.get_bmatrix()
    assert B1.equals(B2)
    assert B1.equals(B3)
    
@pytest.mark.rdd
def test_qmatrix_methods():
    """Test that all accessible methods to qmatrix are working and are equivalent.
    """
    tape = Endf6.from_text("\n".join(RDD.endf6))
    DC = DecayChains.from_endf6(tape)
    B = DC.get_bmatrix()
    Q1 = B.to_qmatrix()
    Q2 = DC.get_qmatrix()
    Q3 = QMatrix.from_file(os.path.join(RDD.__path__[0], "RDD.jeff33"))
    assert Q1.equals(Q2)
    assert Q1.equals(Q3)