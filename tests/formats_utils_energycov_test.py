# -*- coding: utf-8 -*-
"""
Created on Thu Feb  7 15:06:06 2019

@author: Fiorito_L
"""

import pytest

import numpy as np
import sandy

@pytest.fixture(scope="module")
def cov1():
    cov1 = np.array([[1.  , 0.  , 1.25, 0.  ],
                     [0.  , 9.  , 1.5 , 0.  ],
                     [1.25, 1.5 , 6.25, 0.  ],
                     [0.  , 0.  , 0.  , 0.  ]])
    return cov1

@pytest.fixture(scope="module")
def text_cov_lb6():
    """JEFF-3.2 covariance for 10-Ne-20g (changed MAT1 from 10020 to 0 and 
    MT1 from 91 to 16)"""
    text = \
    """
 1.002000+4 1.982070+1          0          0          0          1102533 16    1
 0.000000+0 0.000000+0          0         16          0          1102533 16   23
 0.000000+0 0.000000+0          0          6         36          5102533 16   24
 1.000000-5 1.781240+7 2.100000+7 2.500000+7 2.900000+7 1.000000-5102533 16   25
 1.037110+7 1.350000+7 1.700000+7 2.100000+7 2.500000+7 2.900000+7102533 16   26
 2.973560-3 3.358410-3 3.619690-3 3.483760-3 3.282020-3 3.103920-3102533 16   27
 2.473770-3 3.361860-3 3.843270-3 3.931630-3 3.748380-3 3.483920-3102533 16   28
 2.231060-3 3.352350-3 4.109790-3 4.006550-3 3.596800-3 3.193460-3102533 16   29
 1.757740-3 2.928100-3 3.705960-3 3.961660-3 3.761260-3 3.460470-3102533 16   30
    """
    return text

@pytest.fixture(scope="module")
def text_cov_lb5_asym():
    """JEFF-3.1.1 covariance for 9-F-19g (changed MT1 from 16 to 4)"""
    text = \
    """
 9.019000+3 1.883500+1          0          0          0          1 92533  4    1
 0.000000+0 0.000000+0          0          4          0          1 92533  4   19
 0.000000+0 0.000000+0          0          5         43          7 92533  4   20
 1.000000-5 1.098500+7 1.200000+7 1.400000+7 1.600000+7 1.800000+7 92533  4   21
 2.000000+7 0.000000+0 0.000000+0 0.000000+0 0.000000+0 0.000000+0 92533  4   22
 0.000000+0 0.000000+0-2.732300-3-8.431200-4 1.625000-3 2.647600-3 92533  4   23
 4.939500-3 0.000000+0-1.625700-3-3.558100-4 1.459800-3 2.437400-3 92533  4   24
 3.772000-3 0.000000+0-2.314700-4 3.051500-4 1.292200-3 2.152300-3 92533  4   25
 2.538300-3 0.000000+0 3.229600-4 6.356200-4 1.308700-3 1.986300-3 92533  4   26
 2.439300-3 0.000000+0 4.981700-4 8.811700-4 1.559200-3 2.135100-3 92533  4   27
 2.953900-3                                                        92533  4   28
    """
    return text

@pytest.fixture(scope="module")
def ecov_sym(cov1):
    """Test `EnergyCov` initialization"""
    e = [1e-5, 10, 1e6, 2e7]
    ecov = sandy.formats.utils.EnergyCov(cov1, index=e, columns=e)
    assert ecov.index.name == "E"
    assert ecov.columns.name == "E"
    return ecov

@pytest.mark.formats
@pytest.mark.utils
@pytest.mark.cov
def test_energy_cov_change_grid_interp(ecov_sym):
    """Check change_grid method for interpolation"""
    ex = [1e0, 1e2, 2e6, 3e6, 4e6]
    ey = [1e0, 2e1]
    ecov = ecov_sym.change_grid(ex=ex, ey=ey)
    assert type(ecov) == sandy.formats.utils.EnergyCov
    assert ecov.index.values.tolist() == ex
    assert ecov.columns.values.tolist() == ey
    assert (ecov.values == [[ecov_sym.values[0,0], ecov_sym.values[0,1]],
                            [ecov_sym.values[1,0], ecov_sym.values[1,1]],
                            [ecov_sym.values[2,0], ecov_sym.values[2,1]],
                            [ecov_sym.values[2,0], ecov_sym.values[2,1]],
                            [ecov_sym.values[2,0], ecov_sym.values[2,1]]]).all()

@pytest.mark.formats
@pytest.mark.utils
@pytest.mark.cov
def test_energy_cov_change_grid_extrap(ecov_sym):
    """Check change_grid method for extrapolation"""
    ex = ey = [1e-6, 3e7]
    ecov = ecov_sym.change_grid(ex=ex, ey=ey)
    assert (ecov.values == [[0, 0], [0, 0]]).all()

@pytest.mark.formats
@pytest.mark.utils
@pytest.mark.cov
def test_energy_cov_sum_covs_1(ecov_sym):
    """Check sum_covs method"""
    e = [10, 20, 2e6]
    cov = np.array([[ 0.5, 0.2, 0.0],
                    [ 0.2, 1.2, 0.0],
                    [ 0.0, 0.0, 5.0]])
    ecov = sandy.formats.utils.EnergyCov(cov, index=e, columns=e)
    out = sandy.formats.utils.EnergyCov.sum_covs(ecov_sym, ecov)
    outr = sandy.formats.utils.EnergyCov.sum_covs(ecov, ecov_sym)
    assert (out.values == outr.values).all()
    assert (out.index == outr.index).all()
    assert (out.columns == outr.columns).all()
    assert (out.index == np.array([1e-5, 10, 20, 1e6, 2e6, 2e7])).all()
    assert (out.columns == np.array([1e-5, 10, 20, 1e6, 2e6, 2e7])).all()
    assert (out.values == np.array([[ 1.  ,  0.  ,  0.  ,  1.25,  1.25,  0.  ],
                                    [ 0.  ,  9.5 ,  9.2 ,  1.7 ,  1.5 ,  0.  ],
                                    [ 0.  ,  9.2 , 10.2 ,  2.7 ,  1.5 ,  0.  ],
                                    [ 1.25,  1.7 ,  2.7 ,  7.45,  6.25,  0.  ],
                                    [ 1.25,  1.5 ,  1.5 ,  6.25, 11.25,  5.  ],
                                    [ 0.  ,  0.  ,  0.  ,  0.  ,  5.  ,  5.  ]])).all()

@pytest.mark.formats
@pytest.mark.utils
@pytest.mark.cov
def test_energy_cov_sum_covs_2():
    """Check sum_covs method for empty matrices returns empty df"""
    emptydf = sandy.formats.utils.EnergyCov.sum_covs(sandy.formats.utils.EnergyCov(), sandy.formats.utils.EnergyCov())
    assert emptydf.empty

@pytest.mark.formats
@pytest.mark.utils
@pytest.mark.cov
def test_energy_cov_sum_covs_3():
    """Check that sum_covs method for empty args raises error"""
    with pytest.raises(Exception):
        sandy.formats.utils.EnergyCov.sum_covs()

@pytest.mark.formats
@pytest.mark.utils
@pytest.mark.cov
def test_energycov_lb6(text_cov_lb6):
    """Check that MF31/33 LB=6 is processed correctly"""
    tape = sandy.formats.endf6.Endf6.from_text(text_cov_lb6)
    sec = tape.read_section(1025, 33, 16)
    nisec = sec["SUB"][16]["NI"][0]
    ecov = sandy.formats.utils.EnergyCov.from_lb6(nisec["EK"], nisec["EL"], nisec["FKL"])
    assert (ecov.values == np.array([[ 2.973560E-3, 3.358410E-3, 3.619690E-3, 3.483760E-3, 3.282020E-3, 3.103920E-3,          0.],
                                     [ 2.473770E-3, 3.361860E-3, 3.843270E-3, 3.931630E-3, 3.748380E-3, 3.483920E-3,          0.],
                                     [ 2.231060E-3, 3.352350E-3, 4.109790E-3, 4.006550E-3, 3.596800E-3, 3.193460E-3,          0.],
                                     [ 1.757740E-3, 2.928100E-3, 3.705960E-3, 3.961660E-3, 3.761260E-3, 3.460470E-3,          0.],
                                     [          0.,          0.,          0.,          0.,          0.,          0.,          0.]])).all()

@pytest.mark.formats
@pytest.mark.utils
@pytest.mark.cov
def test_energycov_lb5_asym(text_cov_lb5_asym):
    """Check that MF31/33 LB=5 LS=0 is processed correctly"""
    tape = sandy.formats.endf6.Endf6.from_text(text_cov_lb5_asym)
    sec = tape.read_section(925, 33, 4)
    nisec = sec["SUB"][4]["NI"][0]
    ecov = sandy.formats.utils.EnergyCov.from_lb5_asym(nisec["EK"], nisec["FKK"])
    assert (ecov.values == np.array([[0.000000E+0, 0.000000E+0, 0.000000E+0, 0.000000E+0, 0.000000E+0, 0.000000E+0, 0],
                                     [0.000000E+0,-2.732300E-3,-8.431200E-4, 1.625000E-3, 2.647600E-3, 4.939500E-3, 0],
                                     [0.000000E+0,-1.625700E-3,-3.558100E-4, 1.459800E-3, 2.437400E-3, 3.772000E-3, 0],
                                     [0.000000E+0,-2.314700E-4, 3.051500E-4, 1.292200E-3, 2.152300E-3, 2.538300E-3, 0],
                                     [0.000000E+0, 3.229600E-4, 6.356200E-4, 1.308700E-3, 1.986300E-3, 2.439300E-3, 0],
                                     [0.000000E+0, 4.981700E-4, 8.811700E-4, 1.559200E-3, 2.135100E-3, 2.953900E-3, 0],
                                     [0, 0, 0, 0, 0, 0, 0]])).all()