# -*- coding: utf-8 -*-
"""
Created on Tue Nov 20 14:02:00 2018

@author: Fiorito_L
"""

import pytest

import numpy as np

import sandy

__author__ = "Luca Fiorito"

@pytest.mark.ddfy
@pytest.mark.qmatrix
def test_qmatrix():
    Q = sandy.get_jeff_qmatrix()
    assert (np.isclose(np.diag(Q.values), 1)).all()