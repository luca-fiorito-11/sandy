# -*- coding: utf-8 -*-
"""
Created on Tue Dec 18 14:42:46 2018

@author: Luca Fiorito
"""

__author__ = "Luca Fiorito"

import pytest

import pandas as pd
import numpy as np

import sandy


@pytest.fixture(scope="module")
def cov4():
    matrix = np.eye(4)*3
    mat = [9228]
    mt = [18]
    E = [1e-5,1,100,2e7]
    MI = pd.MultiIndex.from_product([mat,mt,E], names=["A", "B", "C"])
    return sandy.formats.utils.BaseCov(matrix, index=MI, columns=MI)

@pytest.mark.formats
@pytest.mark.utils
@pytest.mark.cov
def test_BaseCov_eig(cov4):
    """Test eig method of BaseCov."""
    eigs = cov4.eig()
    assert len(eigs) == 4
    assert eigs.index.tolist() == [0,1,2,3]
    assert (eigs == 3).all()

@pytest.mark.formats
@pytest.mark.utils
@pytest.mark.cov
def test_BaseCov_corr(cov4):
    """Test corr method of BaseCov."""
    corr = cov4.corr()
    assert isinstance(corr, cov4.__class__)
    assert corr.shape == (4,4)
    assert (corr.index == cov4.index).all()
    assert (corr.columns == cov4.columns).all()
    assert corr.index.names == cov4.index.names
    assert corr.columns.names == cov4.columns.names
    assert np.allclose(corr.values, np.eye(4))

@pytest.mark.formats
@pytest.mark.utils
@pytest.mark.cov
def test_BaseCov_get_var(cov4):
    """Test get_var method of BaseCov."""
    var = cov4.get_var()
    assert (var.values == 3).all()
    assert isinstance(var, pd.Series)
    assert (var.index == cov4.index).all()
    assert var.name == "VAR"

@pytest.mark.formats
@pytest.mark.utils
@pytest.mark.cov
def test_BaseCov_get_std(cov4):
    """Test get_std method of BaseCov"""
    std = cov4.get_std()
    assert isinstance(std, pd.Series)
    assert (std.values == np.sqrt(3)).all()
    assert (std.index == cov4.index).all()
    assert std.name == "STD"

@pytest.mark.formats
@pytest.mark.utils
@pytest.mark.cov
def test_BaseCov_to_matrix(cov4):
    """Test to_matrix method of BaseCov"""
    cov = cov4.to_matrix()
    assert (cov == cov4.values).all()
    assert isinstance(cov, sandy.formats.utils.Cov)