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


def gls_update(x_prior, S, Vx_prior, Vy_extra, y_calc, y_extra, threshold=None):
    """
    Perform GlS update for a given variances, vectors and sensitivity.

    Parameters
    ----------
    x_prior: 1D iterable
        Vector in which we are going to apply GLS (MX1)
    Vx_prior : 2D iterable
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
    # Model calculus:
    y_calc_ = _y_calc(x_prior, S)
    # Data in a appropiate format
    x_prior_ = pd.Series(x_prior)
    y_extra_ = pd.Series(y_extra)
    y_calc_ = y_calc_.reindex(y_extra_.index)
    delta = y_extra_ - y_calc_
    Vx_prior_ = sandy.CategoryCov(Vx_prior)
    S_ = pd.DataFrame(S).loc[Vx_prior_.data.index, Vx_prior_.data.columns]
    # GLS update
    A = Vx_prior_._gls_general_sensitivity(S_, Vy_extra, threshold).values
    x_post = x_prior_ + A.dot(delta)
    return x_post


def _y_calc(x_prior, S):
    """
    Perform model calculation in GLS.

    Parameters
    ----------
    x_prior: 1D iterable
        Vector in which we are going to apply GLS (MX1).
    S : 2D iterable
        2D sensitivity of the model y=f(x) (MXN).

    Returns
    -------
    y_calc : `pd.Series`
        1D calculated output using S.dot(x_prior), e.g. calculated CFY

    """
    S_ = pd.DataFrame(S)
    x_prior_ = pd.Series(x_prior).reindex(S_.columns).fillna(0)
    y_calc = S_.dot(x_prior_)
    return y_calc
