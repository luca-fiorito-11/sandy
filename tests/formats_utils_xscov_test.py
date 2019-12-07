# -*- coding: utf-8 -*-
"""
Created on Fri Jan 25 14:36:53 2019

@author: Fiorito_L
"""
__author__ = "Luca Fiorito"

import pytest

import numpy as np
import pandas as pd
import sandy



@pytest.fixture(scope="module")
def cov_lb1():
    text = \
    """
 9.223500+4 2.330248+2          0          0          0          1922831455    1
 0.000000+0 0.000000+0          0        455          0          1922831455    2
 0.000000+0 0.000000+0          0          1          6          3922831455    3
 1.000000-5 2.025000-3 1.000000+2 2.209000-3 3.000000+7 0.000000+0922831455    4
     """
    tape = sandy.formats.endf6.Endf6.from_text(text)
    return sandy.XsCov.from_endf6(tape)

@pytest.fixture(scope="module")
def cov_lb2():
    text = \
    """
 9.223500+4 2.330248+2          0          0          0          1922831456    1
 3.000000+0 0.000000+0          0        456          0          1922831456    2
 0.000000+0 0.000000+0          0          2          6          3922831456    3
 1.000000-5 1.434000-3 6.000000-1 1.600000-3 3.000000+7 0.000000+0922831456    4
    """
    tape = sandy.formats.endf6.Endf6.from_text(text)
    return sandy.XsCov.from_endf6(tape)

@pytest.fixture(scope="module")
def cov_lb4():
    """JEFF-3.2 covariance for 26-Fe-54g"""
    text = \
    """
 2.605400+4 5.347625+1          0          0          0          1262533103    1
 0.000000+0 0.000000+0          0        103          0          1262533103    2
 0.000000+0 0.000000+0         18          4         72         36262533103   33
 1.000000-5-1.000000+0 5.000000+5-3.100000-1 1.000000+6-4.600000-1262533103   34
 2.000000+6-5.400000-1 3.000000+6-5.500000-1 3.350000+6-5.400000-1262533103   35
 4.000000+6-4.700000-1 5.000000+6-1.000000+0 5.500000+6-1.000000+0262533103   36
 6.000000+6-1.000000+0 7.000000+6-7.300000-1 9.500000+6-1.000000+0262533103   37
 1.000000+7-1.000000+0 1.050000+7-7.700000-1 1.200000+7-6.300000-1262533103   38
 1.400000+7-7.300000-1 1.800000+7-1.000000+0 2.000000+7 0.000000+0262533103   39
 1.000000-5 0.000000+0 5.000000+5 1.621100-1 1.000000+6 7.465900-2262533103   40
 2.000000+6 3.963200-2 3.000000+6 3.935300-2 3.350000+6 3.820200-2262533103   41
 4.000000+6 3.026000-2 5.000000+6 2.705200-2 5.500000+6 2.871100-2262533103   42
 6.000000+6 2.548700-2 7.000000+6 2.332100-2 9.500000+6 2.796500-2262533103   43
 1.000000+7 2.765200-2 1.050000+7 2.916100-2 1.200000+7 5.426700-2262533103   44
 1.400000+7 2.941600-2 1.800000+7 9.267700-2 2.000000+7 0.000000+0262533103   45
    """

@pytest.fixture(scope="module")
def cov_lb5_sym():
    """JEFF-3.1.1 covariance for 40-Zr-90g"""
    text = \
    """
 4.009000+4 8.913240+1          0          0          0          1402533 16    1
 0.000000+0 0.000000+0          0         16          0          1402533 16    2
 0.000000+0 0.000000+0          1          5         21          6402533 16    3
 1.000000-5 1.211560+7 1.400000+7 1.600000+7 1.800000+7 2.000000+7402533 16    4
 0.000000+0 0.000000+0 0.000000+0 0.000000+0 0.000000+0 9.610000-4402533 16    5
 7.102100-4 7.142400-4 6.919200-4 8.410000-4 7.145600-4 7.012200-4402533 16    6
 1.024000-3 8.233600-4 9.610000-4                                 402533 16    7
    """
    tape = sandy.formats.endf6.Endf6.from_text(text)
    return sandy.XsCov.from_endf6(tape)

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
def ecov_sym():
    e = [1e-5, 10, 1e6, 2e7]
    cov = np.array([[1.  , 0.  , 1.25, 0.  ],
                    [0.  , 9.  , 1.5 , 0.  ],
                    [1.25, 1.5 , 6.25, 0.  ],
                    [0.  , 0.  , 0.  , 0.  ]])
    return sandy.formats.utils.EnergyCov(cov, index=e, columns=e)
    
@pytest.mark.formats
@pytest.mark.utils
@pytest.mark.cov
def test_xscov_lb1(cov_lb1):
    """Check that MF31/33 LB=1 is processed correctly"""
    assert (cov_lb1.values == cov_lb1.values.T).all()
    assert cov_lb1.values[0,0] == 0.002025
    assert cov_lb1.values[1,1] == 0.002209
    assert (cov_lb1.values[2] == 0).all()
    assert cov_lb1.values[1,0] == 0

@pytest.mark.formats
@pytest.mark.utils
@pytest.mark.cov
def test_xscov_lb2(cov_lb2):
    """Check that MF31/33 LB=2 is processed correctly"""
    assert (cov_lb2.values == cov_lb2.values.T).all()
    assert cov_lb2.values[0,0] == 2.056356e-06
    assert cov_lb2.values[1,1] == 2.560000e-06
    assert (cov_lb2.values[2] == 0).all()
    assert cov_lb2.values[1,0] == 2.294400e-06

@pytest.mark.formats
@pytest.mark.utils
@pytest.mark.cov
def test_xscov_lb5_sym(cov_lb5_sym):
    """Check that MF31/33 LB=5 LS=1 is processed correctly"""
    assert (cov_lb5_sym.values == cov_lb5_sym.values.T).all()
    assert cov_lb5_sym.values[1,1] == 9.610000E-4
    assert cov_lb5_sym.values[2,1] == 7.102100E-4
    assert cov_lb5_sym.values[3,1] == 7.142400E-4
    assert cov_lb5_sym.values[4,1] == 6.919200E-4
    assert cov_lb5_sym.values[2,2] == 8.410000E-4
    assert cov_lb5_sym.values[3,2] == 7.145600E-4
    assert cov_lb5_sym.values[4,2] == 7.012200E-4
    assert cov_lb5_sym.values[3,3] == 1.024000E-3
    assert cov_lb5_sym.values[4,3] == 8.233600E-4
    assert cov_lb5_sym.values[4,4] == 9.610000E-4
    assert (cov_lb5_sym.values[0] == 0).all()
    assert (cov_lb5_sym.values[5] == 0).all()
   
@pytest.mark.formats
@pytest.mark.utils
@pytest.mark.cov
def test_xscov_wrong_shape(text_cov_lb6):
    """Check that `XsCov` raises error because matrix withwrong different shape"""
    tape = sandy.formats.endf6.Endf6.from_text(text_cov_lb6)
    with pytest.raises(Exception):
        return sandy.XsCov.from_endf6(tape)

@pytest.mark.formats
@pytest.mark.utils
@pytest.mark.cov
def test_xscov_asymmetric(text_cov_lb5_asym):
    """Check that `XsCov` raises error because matrix non symmetric"""
    tape = sandy.formats.endf6.Endf6.from_text(text_cov_lb5_asym)
    with pytest.raises(Exception):
        return sandy.XsCov.from_endf6(tape)

@pytest.mark.formats
@pytest.mark.utils
@pytest.mark.cov
def test_xscov__from_list(cov_lb5_sym, cov_lb1, cov_lb2):
    """Check that `XsCov` can be created from list"""
    V1 = cov_lb5_sym.get_section(4025,  16, 4025,  16)
    V2 = cov_lb1.get_section(9228, 455, 9228, 455)
    C12 = cov_lb2.get_section(9228, 456, 9228, 456)
    L = [[(1,  18), (1,  18), V1],
         [(2, 102), (2, 102), V2],
         [(1,  18), (2, 102), C12]]
    cov = sandy.XsCov._from_list(L)
    assert (cov.loc[(1,18),(1,18)].values == V1.values).all()
    assert (cov.loc[(2,102),(2,102)].values == V2.values).all()
    C12_new = C12.change_grid(ex=V1.index.values, ey=V2.columns.values)
    assert (cov.loc[(1,18),(2,102)].values == C12_new.values).all()
    assert (cov.loc[(2,102),(1,18)].values == C12_new.T.values).all()
    assert cov.index.get_level_values("E").tolist() == \
           cov.index.get_level_values("E").tolist() == \
           V1.index.values.tolist() + V2.index.values.tolist()