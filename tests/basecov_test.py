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
    return sandy.BaseCov(matrix, index=MI, columns=MI)



@pytest.fixture(scope="module")
def eigs4(cov4):
    return cov4.eig()

def test_BaseCov_eig_values(eigs4):
    np.testing.assert_equal(eigs4.values, 3)

def test_BaseCov_eig_type(eigs4):
    assert isinstance(eigs4, pd.Series) 

def test_BaseCov_eig_index(eigs4):
    np.testing.assert_array_equal(eigs4.index, [0, 1, 2, 3])

def test_BaseCov_eig_name(eigs4):
    assert eigs4.name == "eigenvalues"



@pytest.fixture(scope="module")
def corr4(cov4):
    return cov4.corr()

def test_BaseCov_corr_type(corr4, cov4):
    assert isinstance(corr4, cov4.__class__)

def test_BaseCov_corr_values(corr4):
    assert corr4.values == pytest.approx(np.eye(4))

def test_BaseCov_corr_index(corr4, cov4):
    assert corr4.index.equals(cov4.index)

def test_BaseCov_corr_index_names(corr4, cov4):
    assert corr.index.names == cov4.index.names

def test_BaseCov_corr_index(corr4, cov4):
    assert corr4.columns.equals(cov4.columns)

def test_BaseCov_corr_index_names(corr4, cov4):
    assert corr4.columns.names == cov4.columns.names



@pytest.fixture(scope="module")
def var4(cov4):
    return cov4.get_var()

def test_BaseCov_get_var_values(var4):
    np.testing.assert_equal(var4.values, 3)

def test_BaseCov_get_var_type(var4):
    assert isinstance(var4, pd.Series)

def test_BaseCov_get_var_index(var4, cov4):
    assert var4.index.equals(cov4.index)

def test_BaseCov_get_var_index_names(var4, cov4):
    assert var4.index.names == cov4.index.names

def test_BaseCov_get_var_name(var4):
    assert var4.name == "VAR"



@pytest.fixture(scope="module")
def std4(cov4):
    return cov4.get_std()

def test_BaseCov_get_std_values(std4):
    assert std4.values == pytest.approx(np.sqrt(3))

def test_BaseCov_get_std_type(std4):
    assert isinstance(std4, pd.Series)

def test_BaseCov_get_std_index(std4, cov4):
    assert std4.index.equals(cov4.index)

def test_BaseCov_get_std_index_names(std4, cov4):
    assert std4.index.names == cov4.index.names

def test_BaseCov_get_std_name(std4):
    assert std4.name == "STD"



def test_BaseCov_to_matrix(cov4):
    cov = cov4.to_matrix()
    np.testing.assert_array_equal(cov, cov4.values)
    assert isinstance(cov, sandy.Cov)



def test_BaseCov_filter_by(cov4):
    C = cov4.filter_by("A", [9228], "C", [1.0, 1e2, 1e3])
    assert C.index.names == cov4.index.names
    assert C.columns.names == C.columns.names
    assert (C.index.values == cov4.index.values).all()
    assert (C.columns.values == cov4.columns.values[1:3]).all()
    assert (C.values == cov4.values[:,1:3]).all()
    assert isinstance(C, cov4.__class__)

def test_BaseCov_filter_by_empty(cov4):
    with pytest.raises(Exception):
        cov4.filter_by("A", [9228], "C", [1e-8, 1e9])

def test_BaseCov_filter_wrong_index(cov4):
    with pytest.raises(Exception):
        cov4.filter_by("A", [9228], "D", [1.0, 1e2, 1e3])