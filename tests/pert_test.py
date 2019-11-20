import pytest
from io import StringIO

import numpy as np
import pandas as pd

import sandy

__author__ = "Luca Fiorito"

#####################
# Test initialization
#####################
def test_from_file_1_column():
    vals = '1\n5\n9'
    file = StringIO(vals)
    with pytest.raises(Exception):
        sandy.Pert.from_file(file)

def test_from_file_non_monotonic():
    vals = '1  1\n6  5\n5  2\n9  3'
    file = StringIO(vals)
    with pytest.raises(Exception):
        sandy.Pert.from_file(file)

@pytest.fixture(scope="module")
def pert3():
    vals = '1  1  5\n5  2  1\n9  3  1'
    file = StringIO(vals)
    return sandy.Pert.from_file(file)

def test_from_file_3_columns(pert3):
    # should try and catch the warning
    pass

def test_init_with_series(pert3):
    pert = sandy.Pert(pert3.right)
    assert pert.data.equals(pert3.data)

def test_init_with_dataframe(pert3):
    with pytest.raises(Exception):
        sandy.Pert(pert3.right.to_frame())

def test_init_with_intervalindex(pert3):
    pert = sandy.Pert(pert3.data)
    assert pert.data.equals(pert3.data)

################################
# Test initialization attributes
################################

def test_Pert_type(pert3):
    assert isinstance(pert3, sandy.Pert)

def test_Pert_data_index_type(pert3):
    assert isinstance(pert3.data.index, pd.IntervalIndex)

def test_Pert_data_index_right_values(pert3):
    assert pert3.data.index.right.tolist() == [1, 5, 9]

def test_Pert_data_index_left_values(pert3):
    assert pert3.data.index.left.tolist() == [0, 1, 5]

def test_Pert_data_index_float(pert3):
    assert pert3.data.index.right.values.dtype == float

def test_Pert_data_values(pert3):
    np.testing.assert_array_equal(pert3.data.values, [1,2,3])

def test_Pert_data_values_float(pert3):
    assert pert3.data.values.dtype == float
    

########################
# Test attributes
########################



########################
# Test methods
########################

# def test_Spectrum_selfreshape(spec_const):
    # S = spec_const.reshape(spec_const.right.index)
    # assert np.allclose(S.data.values,spec_const.data.values)

# @pytest.mark.parametrize("eg, flux",
                         # [
                                 # ([30], 500),
                                 # ([6e-12], 0.6),
                                 # ([5e-12], 0.5),
                                 # ([4e-12], 0.4),
                                 # ([1e-11], 1),
                                 # ([18.896380829766173], 499),
                                 # ([1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1e0, 1e1, 20], 500),
                         # ]
                         # )
# def test_Spectrum_reshape(spec_const, eg, flux):
    # S = spec_const.reshape(eg)
    # assert S.right.index.tolist() == eg
    # assert S.flux == pytest.approx(flux)

# @pytest.mark.parametrize("eg, err, errtype",
                         # [
                                 # ([2, 1], True, ValueError),
                                 # ([2, 2], False, None),
                                 # ([-1, 2], True, ValueError),
                         # ]
                         # )
# def test_Spectrum_reshape_error(spec_const, eg, err, errtype):
    # if err:
        # with pytest.raises(errtype):
            # spec_const.reshape(eg)
    # else:
        # spec_const.reshape(eg)