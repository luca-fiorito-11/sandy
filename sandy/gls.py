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


def gls_update(y, S, Vx_calc, Vy_extra, y_calc, y_extra, threshold=None):
    """
    Perform GlS update for a given variances, vectors and sensitivity.

    Parameters
    ----------
    y: 1D iterable
        Vector in which we are going to apply GLS (MX1)
    Vx_calc : 2D iterable
        2D covariance matrix of x_prior (MXN).
    Vy_extra : 2D iterable
        2D covariance matrix for y_extra (MXN).
    S : 2D iterable
        2D sensitivity of the model y=f(x) (MXN).
    y_calc : 1D iterable
        1D calculated output using S.dot(x_prior), e.g. calculated CFY (NX1)
    y_extra : 1D iterable
        1D extra info on output (NX1)
    threshold : `int`, optional
            Thereshold to avoid numerical fluctuations. The default is None.

    Returns
    -------
    `pd.Series`
        GLS performation for a vector y given S, Vx_prior, Vy_extra, y_calc
        and y_extra.

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
    y_ = pd.Series(y).values
    index = pd.Series(y).index
    y_calc_ = pd.Series(y_calc).values
    y_extra_ = pd.Series(y_extra).values
    delta = y_extra_ - y_calc_
    Vx_calc_ = sandy.CategoryCov(Vx_calc)
    S_ = pd.DataFrame(S).loc[index, pd.Series(y_calc).index]
    A = Vx_calc_._gls_general_sensitivity(S_, Vy_extra, threshold).values
    y_new = y_ + A.dot(delta)
    return pd.Series(y_new, index=index)
