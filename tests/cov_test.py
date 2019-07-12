# -*- coding: utf-8 -*-
"""
Created on Wed Feb  6 16:37:07 2019

@author: Fiorito_L
"""

__author__ = "Luca Fiorito"

import pytest

import numpy as np
import pandas as pd
import sandy

@pytest.fixture(scope="module")
def testCov():
    cov = np.array([[1.  , 0.  , 1.25, 0.  ],
                    [0.  , 9.  , 1.5 , 0.  ],
                    [1.25, 1.5 , 6.25, 0.  ],
                    [0.  , 0.  , 0.  , 0.  ]])
    return sandy.Cov(cov)

@pytest.mark.formats
@pytest.mark.utils
@pytest.mark.cov
def test_Cov_asym():
    """Check that Cov raises error during initialization when matrix is asymmetric"""
    cov = np.array([[1.  , 0.  , 1.25, 0.  ],
                    [0.  , 9.  , 5.0 , 0.  ],
                    [1.25, 1.5 , 6.25, 0.  ],
                    [0.  , 0.  , 0.  , 0.  ]])
    with pytest.raises(Exception):
        sandy.Cov(cov)

@pytest.mark.formats
@pytest.mark.utils
@pytest.mark.cov
def test_Cov_2D():
    """Check that Cov raises error during initialization when matrix is not 2D"""
    cov = np.array([1.  , 0.  ])
    with pytest.raises(Exception):
        sandy.Cov(cov)

@pytest.mark.formats
@pytest.mark.utils
@pytest.mark.cov
def test_Cov_integer():
    """Check that Cov converts integer matrix into float"""
    cov = np.array([[1  , 0  ],
                    [0  , 9  ]], dtype=int)
    C = sandy.Cov(cov)
    assert C.dtype == float

@pytest.mark.formats
@pytest.mark.utils
@pytest.mark.cov
def test_Cov_not_square():
    """Check that Cov raises error during initialization when matrix is not square"""
    cov = np.array([[1.  , 0.  , 1.25, 0.  ],
                    [0.  , 9.  , 5.0 , 0.  ],
                    [1.25, 1.5 , 6.25, 0.  ]])
    with pytest.raises(Exception):
        sandy.Cov(cov)

@pytest.mark.formats
@pytest.mark.utils
@pytest.mark.cov
def test_Cov_negative_diag():
    """Check that Cov raises error during initialization when matrix has negative variances"""
    cov = np.array([[1.  , 0.  , 1.25, 0.  ],
                    [0.  , -5. , 1.5 , 0.  ],
                    [1.25, 1.5 , 6.25, 0.  ],
                    [0.  , 0.  , 0.  , 0.  ]])
    with pytest.raises(Exception):
        sandy.Cov(cov)

@pytest.mark.formats
@pytest.mark.utils
@pytest.mark.cov
def test_Cov_eig(testCov):
    """Check that Cov method eig works correctly"""
    E, V = testCov.eig()

@pytest.mark.formats
@pytest.mark.utils
@pytest.mark.cov
def test_Cov_corr(testCov):
    """Check that Cov method corr works correctly"""
    C = testCov.corr()
    assert isinstance(C, sandy.Cov)
    assert (C == np.array([[1. , 0. , 0.5, 0. ],
                           [0. , 1. , 0.2, 0. ],
                           [0.5, 0.2, 1. , 0. ],
                           [0. , 0. , 0. , 0. ]])).all()

@pytest.mark.formats
@pytest.mark.utils
@pytest.mark.cov
def test_Cov__reduce_size(testCov):
    """Check that Cov method _reduce_size works correctly"""
    I, C = testCov._reduce_size()
    assert (I == np.array([0, 1, 2])).all()
    assert (C == np.array([[1.  , 0.  , 1.25],
                           [0.  , 9.  , 1.5 ],
                           [1.25, 1.5 , 6.25]])).all()
    assert isinstance(C, sandy.Cov)

@pytest.mark.formats
@pytest.mark.utils
@pytest.mark.cov
def test_Cov__reduce_size_all_zeros():
    """Check that Cov method _reduce_size works correctly"""
    cov = sandy.Cov(np.zeros((3,3)))
    I, C = cov._reduce_size()
    assert not I
    assert not C

@pytest.mark.formats
@pytest.mark.utils
@pytest.mark.cov
def test_Cov__restore_size(testCov):
    """Check that Cov method _restore_size works correctly"""
    I = np.array([0, 1, 2])
    C = np.array([[1.  , 0.  , 1.25],
                  [0.  , 9.  , 1.5 ],
                  [1.25, 1.5 , 6.25]])
    cov = sandy.Cov._restore_size(I, C, 4)
    assert (cov == testCov).all()
    assert isinstance(cov, sandy.Cov)

@pytest.mark.formats
@pytest.mark.utils
@pytest.mark.cov
def test_Cov__restore_size_all_zeros():
    """Check that Cov method _restore_size works correctly"""
    I = np.array([])
    C = sandy.Cov(np.ndarray((0,0)))
    cov = sandy.Cov._restore_size(I, C, 3)
    assert (cov == np.zeros((3,3))).all()
    assert isinstance(cov, sandy.Cov)

@pytest.mark.formats
@pytest.mark.utils
@pytest.mark.cov
def test_Cov_get_L_cholesky():
    """Check that Cov method get_L works correctly with cholesky"""
    C = np.array([[1.  , 0.  ],
                  [0.  , 2   ]])
    cov = sandy.Cov(C)
    L = cov.get_L()
    assert np.allclose(L, np.sqrt(cov))

@pytest.mark.formats
@pytest.mark.utils
@pytest.mark.cov
def test_Cov_get_L_eigendecomp():
    """Check that Cov method get_L works correctly with eidecomp."""
    C = np.array([[1.  , 0.  , 1.25],
                  [0.  , 0.2 , 1.5 ],
                  [1.25, 1.5 , 6.25]])
    cov = sandy.Cov(C)
    L = cov.get_L()
    assert np.allclose(L, np.array([[-1.00842023,  0.        ,  0.        ],
                                    [-0.05742429, -0.62848243,  0.        ],
                                    [-1.22303907, -2.18416873,  0.        ]]))
    CC = sandy.Cov(L.dot(L.T))
    E1 = cov.eig()[0]
    E2 = CC.eig()[0]
    assert np.allclose(E1[:-1], E2[:-1])
    assert np.isclose(E2[-1], 0)

@pytest.mark.formats
@pytest.mark.utils
@pytest.mark.cov
def test_Cov_sampling():
    """Check that Cov method sampling works correctly."""
    C = np.array([[1.  , 0.  , 1.25],
                  [0.  , 0.2 , 1.5 ],
                  [1.25, 1.5 , 6.25]])
    cov = sandy.Cov(C)
    smp = cov.sampling(4, seed=1587)
    assert smp.shape == (3, 4)
    assert np.allclose(smp, np.array([[ 0.99853553, -0.34123088,  0.31211793,  0.79757716],
                                      [-0.09853738, -0.34152333,  0.87597002, -0.92308259],
                                      [ 0.67099231, -1.53322223,  3.36104051, -2.39851248]]))

@pytest.mark.formats
@pytest.mark.utils
@pytest.mark.cov
def test_Cov_sampling_1sample():
    """Check that Cov method sampling works correctly with 1 sample."""
    C = np.array([[1.  , 0.  , 1.25],
                  [0.  , 0.2 , 1.5 ],
                  [1.25, 1.5 , 6.25]])
    cov = sandy.Cov(C)
    smp = cov.sampling(1)
    assert smp.shape == (3, 1)
    assert not isinstance(smp, sandy.Cov)