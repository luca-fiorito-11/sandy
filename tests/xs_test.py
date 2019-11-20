
__author__ = "Luca Fiorito"

import pytest
import io

import numpy as np
import pandas as pd

import sandy


text1 = """ 1.001000+3 9.991673-1          0          0          0          0 125 3  1    1
 0.000000+0 0.000000+0          0          0          1          3 125 3  1    2
          3          2                                             125 3  1    3
 1.000000-5 1.00000000 1.000000+0 2.00000000 1.000000+5 3.00000000 125 3  1    4
"""
 
 
@pytest.fixture(scope="module")
def tape1():
    return sandy.Endf6.from_text(text1)

@pytest.fixture(scope="module")
def xs1(tape1):
    return sandy.Xs.from_endf6(tape1)

text4 = """ 1.001000+3 9.991673-1          0          0          0          0 125 3  4    1
 0.000000+0 0.000000+0          0          0          1          2 125 3  4    2
          2          2                                             125 3  4    3
 1.000000+2 1.00000000 1.000000+5 4.00000000                       125 3  4    4
"""

@pytest.fixture(scope="module")
def tape4():
    return sandy.Endf6.from_text("\n".join((text1, text4)))

@pytest.fixture(scope="module")
def xs4(tape4):
    return sandy.Xs.from_endf6(tape4)


def test_from_endf6_non_monotonic_grid():
    text = """ 1.001000+3 9.991673-1          0          0          0          0 125 3  1    1
 0.000000+0 0.000000+0          0          0          1          3 125 3  1    2
          3          2                                             125 3  1    3
 1.000000+3 1.00000000 1.000000+0 2.00000000 1.000000+5 3.00000000 125 3  1    4
"""
    tape = sandy.Endf6.from_text(text)
    with pytest.raises(sandy.Error):
        sandy.Xs.from_endf6(tape)

def test_from_endf6_duplicated_values():
    text = """ 1.001000+3 9.991673-1          0          0          0          0 125 3  1    1
 0.000000+0 0.000000+0          0          0          1          3 125 3  1    2
          3          2                                             125 3  1    3
 1.000000+0 1.00000000 1.000000+0 2.00000000 1.000000+5 3.00000000 125 3  1    4
"""
    tape = sandy.Endf6.from_text(text)
    xs = sandy.Xs.from_endf6(tape)
    assert xs.data.index.size == 2
    np.testing.assert_array_equal(xs.data.index, [1, 1e5])



def test_from_endf6_one_reaction(xs1):
    pass

def test_to_endf6_one_reaction(xs1, tape1):
    tape = xs1.to_endf6(tape1)
    assert tape.iloc[0].TEXT == text1


def test_from_endf6_two_reactions(xs4):
    pass

def test_to_endf6_two_reactions(xs4, tape4):
    tape = xs4.to_endf6(tape4)
    assert tape.loc[(125,3,1)].TEXT == text1
    assert tape.loc[(125,3,4)].TEXT == text4



def test_reshape(xs1):
    xsnew = xs1.reshape([1e-4, 1e0])
    np.testing.assert_array_equal(xsnew.data.index, [1e-5, 1e-4, 1e0, 1e5])
    oldvals = xs1.data.values.T[0]
    vals = np.array([oldvals[0],
            oldvals[0]*(1e-4-1e-5)/(1e0-1e-5)*(oldvals[1]-oldvals[0])+oldvals[0],
            oldvals[1],
            oldvals[2],
    ])
    np.testing.assert_array_almost_equal(xsnew.data.values.T[0], vals)

def test_reshape_extrapolate_left(xs1):
    xsnew = xs1.reshape([1e-6])
    np.testing.assert_array_equal(xsnew.data.index, [1e-6, 1e-5, 1e0, 1e5])
    oldvals = xs1.data.values.T[0]
    vals = np.array([0,
            oldvals[0],
            oldvals[1],
            oldvals[2],
    ])
    np.testing.assert_array_almost_equal(xsnew.data.values.T[0], vals)

def test_reshape_extrapolate_right(xs1):
    xsnew = xs1.reshape([1e8])
    np.testing.assert_array_equal(xsnew.data.index, [1e-5, 1e0, 1e5, 1e8])
    oldvals = xs1.data.values.T[0]
    vals = np.array([oldvals[0],
            oldvals[1],
            oldvals[2],
            0,
    ])
    np.testing.assert_array_almost_equal(xsnew.data.values.T[0], vals)



def test_custom_perturbation(xs1):
    pert = sandy.Pert(pd.Series([1, 1.1, 0.7], index=[1e-5, 1e-1, 1e0]))
    xs1.custom_perturbation(125, 1, pert)