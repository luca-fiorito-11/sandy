"""
Collection of utilities, functions and classes that are requetsed in all code components.
"""

import pdb
import os
import re
from itertools import zip_longest
from scipy.constants import Avogadro

import numpy as np
import pandas as pd
import scipy.interpolate
import numba

__author__ = "Luca Fiorito"
__all__ = [
        "pad_from_beginning",
        "pad_from_beginning_fast",
        "uniform_loggrid",
        ]


MeV_MWs = 1.60217733e-19 # conversion coefficient MeV --> MW s
Amn = 1.00866491578      # molar mass of a neutron in amu, source: P. J. Mohr and B. N. Taylor, "The 1998 CODATA Recommended Values of the Fundamental Physics Constants", Version 3.1


def grouper(iterable, n, fillvalue=None):
    """
    Collect data into fixed-length chunks or blocks
    """
    args = [iter(iterable)] * n
    return zip_longest(*args, fillvalue=fillvalue)

def interwine_lists(a, b):
    return [ z for item in zip(a, b) for z in item ]

def expand_za(za, method="nndc", meta=0):
    z = int(za//1000)
    a = int(za - z*1000)
    if method == "nndc":
        m = 0
        if a >= 300:
            m = 1
            a = a - 300 - m*100
    else:
        m = int(meta)
    return z, a, m

def get_za(z, a, m, method="nndc"):
    za = z*1000 + a + 300 + m*100 if m != 0 and method == "nndc" else z*1000 + a
    return int(za), m

def expand_zam(zam):
    z = int(zam//10000)
    a = int(zam - z*10000)//10
    m = int(zam - z*10000 - a*10)
    return z, a, m

def get_zam(z, a, m):
    zam = z*10000 + a*10 + m
    return int(zam)

def za2zam(za, method="nndc", meta=0):
    return get_zam(*expand_za(za, method=method, meta=meta))

def zam2za(zam, method="nndc"):
    z, a, m = expand_zam(zam)
    return get_za(z, a, m, method=method)

def at2wt(zam, at):
    """
    Given an isotope ZAM id and its atomic density in at/cm/b, convert the 
    latter in g/cm3.
    
    Parameters
    ----------
    zam : `int`
        isotope ZAM id
    at : `float`
        atomic density in at/b/cm
    
    Returns
    -------
    `float`
        weight density in g/cm3
    """
    return at * 1e24 / Avogadro * aleph.common.awr[zam] * aleph.common.Amn

def wt2at(zam, wt):
    """
    Given an isotope ZAM id and its weight density in g/cm3, convert the 
    latter in at/cm/b.
    
    Parameters
    ----------
    zam : `int`
        isotope ZAM id
    wt : `float`
        weight density in g/cm3
    
    Returns
    -------
    `float`
        atomic density in at/b/cm
    """
    z, a, m = expand_zam(zam)
    return wt / 1e24 * Avogadro / aleph.common.awr[zam] / aleph.common.Amn


def uniform_loggrid(xmin, xmax, npoints=100):
    """
    Given lower and upper limits, produce a grid with a number of points
    `npoints` that define equivalent intervals in log scale.

    Parameters
    ----------
    xmin : `float`
        lower bound of the grid structure
    xmax : `float`
        upper bound of the grid structure
    npoints : `int`, optional, default `100`

    Returns
    -------
    `numpy` array
        grid equally spaced in logarithmic scale

    Examples
    --------
    >>> uniform_loggrid(1e-5, 1e7, 25)
    array([1.00000000e-05, 3.16227766e-05, 1.00000000e-04, 3.16227766e-04,
           1.00000000e-03, 3.16227766e-03, 1.00000000e-02, 3.16227766e-02,
           1.00000000e-01, 3.16227766e-01, 1.00000000e+00, 3.16227766e+00,
           1.00000000e+01, 3.16227766e+01, 1.00000000e+02, 3.16227766e+02,
           1.00000000e+03, 3.16227766e+03, 1.00000000e+04, 3.16227766e+04,
           1.00000000e+05, 3.16227766e+05, 1.00000000e+06, 3.16227766e+06,
           1.00000000e+07])
    """
    return 10.0**np.linspace(np.log10(xmin), np.log10(xmax), npoints)


def pad_from_beginning(vals, maxlen=None, value=0., axis=0):
    """
    Convert list of arrays into matrix by backward padding.
    .. note:: this function can be used to put cross sections into one matrix
              by adding zeros before the first value of threshold reactions.

    Parameters
    ----------
    vals : `iterable` of arrays/lists
        values of the matrix
    maxlen : `int`, optional, default is `None`
        length to fill with padding.
        If not given, use the maximum length of the arrays in `vals`
        .. important:: if given, maxlen should be
                       `maxlen <= max([len(v) for v in vals])`

    value : `float`, optional, default is `0.`
        value used for padding
    axis : `int`, optional, either `0` or `1`, default is `0`
        axis along whihc the arrays should be positioned.
        `0` means rows, `1` means columns

    Returns
    -------
    `numpy.array`
        2D `numpy.array` with shape `(len(vals), maxlen)` if `axis=0` and
        `(maxlen, len(vals))` if `axis=1`

    Raises
    ------
    `ValueError`
        if `axis` is neither `0` nor `1`
    `ValueError`
        if `maxlen <= max([len(v) for v in vals])`
    """
    length = len(vals)
    lens = [len(v) for v in vals]                     # only iteration
    maxlen_ = max(lens)
    if maxlen is None:
        pass
    elif maxlen < maxlen_:
        raise ValueError("'maxlen' must be >= '{}'".format(maxlen_))
    else:
        maxlen_ = maxlen
    matrix = np.ones((length, maxlen_), dtype=float)*value
    mask = np.arange(maxlen_)[::-1] < np.array(lens)[:, None]  # key line
    matrix[mask] = np.concatenate(vals)
    if axis == 0:
        return matrix
    elif axis == 1:
        return matrix.T
    else:
        raise ValueError("'axis' can be '0' (rows) or '1' (columns), not '{}'".format(axis))


def pad_from_beginning_fast(vals, maxlen):
    """
    Like `aleph.utils.pad_from_beginning` but faster.
    Keyword arguments `axis` and `values` take the default options.
    .. note:: this function can be used to put cross sections into one matrix
              by adding zeros before the first value of threshold reactions.

    Parameters
    ----------
    vals : `iterable` of arrays/lists
        values of the matrix
    maxlen : `int`
        length to fill with padding.

    Returns
    -------
    `numpy.array`
        2D `numpy.array` with shape `(len(vals), maxlen)`
    """
    length = len(vals)
    matrix = np.zeros((length, maxlen))
    lens = [len(v) for v in vals]                    # only iteration
    mask = np.arange(maxlen)[::-1] < np.array(lens)[:, None]  # key line
    matrix[mask] = np.concatenate(vals)
    return matrix


def reshape_differential(x, y, xnew, kind="slinear"):
    """
    Linearly interpolate array over new energy grid structure.

    Extrapolated values are replaced by zeros.

    Parameters
    ----------
    x : 1d array-like object with at least two entries
        energy grid
    xnew : 1d array-like object with at least two entries
        new energy grid
    y : `numpy.ndarray` with at least two entries and same length as `x`
        array to interpolate
    kind: `str` or `int`, optional
        Specifies the kind of interpolation as a string or as an integer
        specifying the order of the spline interpolator to use. The string has
        to be one of ‘linear’, ‘nearest’, ‘nearest-up’, ‘zero’, ‘slinear’,
        ‘quadratic’, ‘cubic’, ‘previous’, or ‘next’. By default is ‘slinear’.

    Returns
    -------
    `numpy.ndarray` with length `len(xnew)`
        interpolated array
    """
    foo = scipy.interpolate.interp1d(
            x, y,
            axis=0,
            copy=False,
            kind=kind,
            bounds_error=False,
            fill_value=0.,
            assume_sorted=True,
            )
    return foo(xnew)


def reshape_integral(x, y, xnew, left_values="first", right_values=0):
    """
    Interpolate array over new energy grid structure using "bfill" method.
    It is assumed that the values of `y` are  multiplied by the grid bin-width.
    The values of the interpolated array are recalculated proportionally to
    the new grid bin-widths.

    Extrapolated values are replaced by zeros.

    Parameters
    ----------
    x : 1d array-like object with at least two entries
        energy grid
    xnew : 1d array-like object with at least two entries
        new energy grid
    y : `numpy.ndarray` with at least two entries and same length as `x`
        array to interpolate

    Returns
    -------
    `numpy.ndarray` with length `len(xnew)`
        interpolated array
    """
    dx = x.copy()
    dx[1:] = np.ediff1d(x)
    dxnew = xnew.copy()
    dxnew[1:] = np.ediff1d(xnew)
    out = reshape_bfill(
        x,
        y / dx,
        xnew,
        left_values=np.nan,
        right_values=np.nan
        ) * dxnew
    return out


def reshape_bfill(x, y, xnew, left_values="first", right_values=0):
    """
    Interpolate array over new energy grid structure using "bfill" method.

    Right-extrapolated values are replaced by zeros.
    Left-extrapolated values are replaced by `y[0]`.

    Parameters
    ----------
    x : 1d array-like object with at least two entries
        energy grid
    xnew : 1d array-like object with at least two entries
        new energy grid
    y : `numpy.ndarray` with at least two entries and same length as `x`
        array to interpolate

    Returns
    -------
    `numpy.ndarray` with length `len(xnew)`
        interpolated array
    """
    fill_value = (left_values, right_values)
    if left_values == "first":
        fill_value[0] = y[0]
    foo = scipy.interpolate.interp1d(
            x, y,
            axis=0,
            copy=False,
            kind="next",
            bounds_error=False,
            fill_value=fill_value,
            assume_sorted=True,
            )
    return foo(xnew)


def add_delimiter_every_n_characters(string, step, delimiter=" "):
    return delimiter.join(string[i:i+step] for i in range(0, len(string), step))


def add_exp_in_endf6_text(text):
    pattern = re.compile("([0-9\.])([+-])([0-9])")
    return pattern.sub("\g<1>E\g<2>\g<3>", text)


def star(func):
    def inner(*args, **kwargs):
        print("*" * 30)
        func(*args, **kwargs)
        print("*" * 30)
    return inner


def percent(func):
    def inner(*args, **kwargs):
        print("%" * 30)
        func(*args, **kwargs)
        print("%" * 30)
    return inner
