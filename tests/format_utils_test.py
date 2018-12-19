# -*- coding: utf-8 -*-
"""
Created on Tue Dec 18 14:42:46 2018

@author: Luca Fiorito
"""

__author__ = "Luca Fiorito"

import pytest

import pandas as pd
import numpy as np

from sandy.formats.utils import *

@pytest.fixture(scope="module")
def cov4():
    matrix = np.eye(4)*3
    mat = [9228]
    mt = [18]
    E = [1e-5,1,100,2e7]
    MI = pd.MultiIndex.from_product([mat,mt,E], names=["MAT", "MT", "E"])
    return XsCov(matrix, index=MI, columns=MI)

@pytest.mark.formats
@pytest.mark.utils
def test_basecov_eig(cov4):
    eigs = cov4.eig()
    assert len(eigs) == 4
    assert eigs.index.tolist() == [0,1,2,3]
    assert (eigs == 3).all()

@pytest.mark.formats
@pytest.mark.utils
def test_basecov_corr(cov4):
    corr = cov4.corr()
    assert corr.shape == (4,4)
    assert (corr.index == cov4.index).all()
    assert (corr.columns == cov4.columns).all()
    assert corr.index.names == cov4.index.names
    assert corr.columns.names == cov4.columns.names
    assert np.isclose(corr.values,np.eye(4)).all()

