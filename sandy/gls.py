# -*- coding: utf-8 -*-
"""
This module contains all classes and functions specific for processing GLS.
"""
import pandas as pd
import sandy

__author__ = "Aitor Bengoechea"
__all__ = [
        "gls_update",
        ]


def gls_update(y, S, Vx, Vy, x, x_p, threshold=None):
    """
    Perform GlS update for a given variances, vectors and sensitivity.

    Parameters
    ----------
    y: 1D iterable
        Vector in which we are going to apply GLS (MX1)
    Vx : 2D iterable
        Covariance vector (MXM).
    Vy : 2D iterable
        Extra Covariance vector (MXM).
    S : 2D iterable
        Sensitivity square matrix (MXM).
    x : 1D iterable
        x vector.
    x_p : 1D iterable
        x vector perturbed.
    threshold : `int`, optional
            Thereshold to avoid numerical fluctuations. The default is None.

    Returns
    -------
    `pd.Series`
        GLS sensitivity for a given Vy and S.

    Example
    -------
    >>> S = [[1, 2], [3, 4]]
    >>> y = pd.Series([1, 1])
    >>> Vx = sandy.CategoryCov.from_var([1, 1]).data
    >>> Vy = pd.DataFrame([[1, 0], [0, 1]], index=[1, 2], columns=[1, 2])
    >>> x = [1, 1]
    >>> x_p = [2, 2]
    >>> gls_update(y, S, Vx, Vy,  x, x_p)
    0   1.02857e+00
    1   1.25714e+00
    dtype: float64
    """
    y_, index = pd.Series(y).values, pd.Series(y).index
    x_, x_p_ = pd.Series(x).values, pd.Series(x_p).values
    delta = x_p_ - x_
    Vx_ = sandy.CategoryCov(Vx)
    A = Vx_._gls_general_sensitivity(S, Vy, threshold).values
    y_new = y_ + A.dot(delta)
    return pd.Series(y_new, index=index)