# -*- coding: utf-8 -*-
"""
This module contains all classes and functions specific for processing GLS.
"""
import pandas as pd
import numpy as np
import sandy

__author__ = "Aitor Bengoechea"
__all__ = [
        "gls_update",
        "_y_calc",
        "chi_individual",
        "chi_diag",
        "chi_square",
        "ishikawa_factor"
        ]

x_prior = [1, 2, 3]
y_extra = pd.Series([2, 3, 4], index=[1, 2, 3])
S = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
Vx_prior = [[0, 0, 0], [0, 3, 0], [0, 0, 8]]
Vy_extra = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
N_e = 1


def gls_update(x_prior, S, Vx_prior, Vy_extra, y_extra, threshold=None):
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
    y_extra : 1D iterable
        1D extra info on output (NX1)
    threshold : `int`, optional
            Thereshold to avoid numerical fluctuations. The default is None.

    Returns
    -------
    `pd.Series`
        GLS apply to a vector x_prior given S, Vx_prior, Vy_extra, y_calc
        and y_extra.

    Example
    -------
    >>> S = [[1, 2], [3, 4]]
    >>> y = pd.Series([1, 1])
    >>> Vx = sandy.CategoryCov.from_var([1, 1]).data
    >>> Vy = pd.DataFrame([[1, 0], [0, 1]], index=[1, 2], columns=[1, 2])
    >>> x = [1, 1]
    >>> x_p = [2, 2]
    >>> gls_update(y, S, Vx, Vy, x_p)
    0   2.00000e-01
    1   4.85714e-01
    dtype: float64

    """
    # Model calculus:
    x_prior_ = pd.Series(x_prior)
    S_ = pd.DataFrame(S).loc[:, x_prior_.index]
    y_calc_ = _y_calc(x_prior, S)
    y_extra_ = pd.Series(y_extra)
    y_calc_ = y_calc_.reindex(y_extra_.index)
    S_ = S_.loc[y_extra_.index, :]
    # Data in a appropriate format
    delta = y_extra_ - y_calc_
    Vx_prior_ = sandy.CategoryCov(Vx_prior)
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

    Example
    -------
    S square matrix:
    >>> _y_calc(x_prior, S)
    0    1
    1    2
    2    3
    dtype: int64

    Different number of row and columns in S:
    >>> S = [[1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 1, 1]]
    >>> _y_calc(x_prior, S)
    0    1
    1    2
    2    3
    3    6
    dtype: int64
    """
    S_ = pd.DataFrame(S)
    x_prior_ = pd.Series(x_prior).reindex(S_.columns).fillna(0)
    y_calc = S_.dot(x_prior_)
    return y_calc


def chi_individual(x_prior, S, Vx_prior, Vy_extra, y_extra):
    """
    Function to calculate individual chi-value measured in sigmas according to
    https://www.oecd-nea.org/jcms/pl_19760/intermediate-report-on-methods-and-approaches-to-provide-feedback-from-nuclear-and-covariance-data-adjustment-for-improvement-of-nuclear-data-files
    (page 9, equation (4.2))

    Parameters
    ----------
    x_prior: 1D iterable
        Vector in which we are going to apply GLS (MX1)
    S : 2D iterable
        2D sensitivity of the model y=f(x) (MXN).
    Vx_prior : 2D iterable
        2D covariance matrix of x_prior (MXN).
    Vy_extra : 2D iterable
        2D covariance matrix for y_extra (MXN).
    y_extra : 1D iterable
        1D extra info on output (NX1).

    Returns
    -------
    `pd.Series`
        individual chi-value measured in sigmas.

    Example
    -------
    >>> chi_individual(x_prior, S, Vx_prior, Vy_extra, y_extra)
    1   1.00000e+00
    2   5.00000e-01
    3   3.33333e-01
    dtype: float64
    """
    Vx_prior_ = sandy.CategoryCov(Vx_prior)
    M = Vx_prior_._gls_M(S, Vy_extra)
    M = np.sqrt(np.diag(M))
    y_calc_ = _y_calc(x_prior, S).values
    y_extra_ = pd.Series(y_extra)
    delta = np.abs(y_extra_.values - y_calc_)
    return pd.Series(delta / M, index=y_extra_.index)


def chi_diag(x_prior, S, Vx_prior, Vy_extra, y_extra):
    """
    Function to calculate diagonal chi-value measured in sigmas
    $\chi_{ind,i}$>>1 according to
    https://www.oecd-nea.org/jcms/pl_19760/intermediate-report-on-methods-and-approaches-to-provide-feedback-from-nuclear-and-covariance-data-adjustment-for-improvement-of-nuclear-data-files
    (page 9, equation (4.3))

    Parameters
    ----------
    x_prior: 1D iterable
        Vector in which we are going to apply GLS (MX1)
    S : 2D iterable
        2D sensitivity of the model y=f(x) (MXN).
    Vx_prior : 2D iterable
        2D covariance matrix of x_prior (MXN).
    Vy_extra : 2D iterable
        2D covariance matrix for y_extra (MXN).
    y_extra : 1D iterable
        1D extra info on output (NX1)

    Returns
    -------
    `pd.Series`
        diagonal chi-value measured in sigmas $\chi_{ind,i}$>>1

    Example
    -------
    >>> chi_diag(x_prior, S, Vx_prior, Vy_extra, y_extra)
    1   1.00000e+00
    2   2.00000e+00
    3   3.00000e+00
    dtype: float64
    """
    Vx_prior_ = sandy.CategoryCov(Vx_prior)
    M = Vx_prior_._gls_M(S, Vy_extra)
    M_inv = sandy.CategoryCov(M).invert().data.values
    M_inv = np.sqrt(np.diag(M_inv))
    y_calc_ = _y_calc(x_prior, S).values
    y_extra_ = pd.Series(y_extra)
    delta = np.abs(y_extra_.values - y_calc_)
    return pd.Series(delta / M_inv, index=y_extra_.index)


def chi_square(x_prior, S, Vx_prior, Vy_extra, y_extra, N_e):
    """
    Function to calculate contribution to chi-square value according to
    https://www.oecd-nea.org/jcms/pl_19760/intermediate-report-on-methods-and-approaches-to-provide-feedback-from-nuclear-and-covariance-data-adjustment-for-improvement-of-nuclear-data-files
    (page 10, equation (4.4))

    Parameters
    ----------
    x_prior: 1D iterable
        Vector in which we are going to apply GLS (MX1)
    S : 2D iterable
        2D sensitivity of the model y=f(x) (MXN).
    Vx_prior : 2D iterable
        2D covariance matrix of x_prior (MXN).
    Vy_extra : 2D iterable
        2D covariance matrix for y_extra (MXN).
    y_extra : 1D iterable
        1D extra info on output (NX1)
    N_e : `int`
        Number of experimental values used in adjustment.

    Returns
    -------
    `pd.Series`
        contribution to chi-square value

    Example
    -------
    >>> chi_square(x_prior, S, Vx_prior, Vy_extra, y_extra, N_e)
    1   1.00000e+00
    2   2.50000e-01
    3   1.11111e-01
    dtype: float64
    """
    Vx_prior_ = sandy.CategoryCov(Vx_prior)
    M = Vx_prior_._gls_M(S, Vy_extra)
    M_inv = sandy.CategoryCov(M).invert().data.values
    y_calc_ = _y_calc(x_prior, S).values
    y_extra_ = pd.Series(y_extra)
    delta = y_extra_.values - y_calc_
    chi_square = delta.T.dot(M_inv) * delta / N_e
    return pd.Series(chi_square, index=y_extra_.index)


def ishikawa_factor(S, Vx_prior, Vy_extra):
    """
    Function to obtain Ishikawa factor according to
    https://www.oecd-nea.org/jcms/pl_19760/intermediate-report-on-methods-and-approaches-to-provide-feedback-from-nuclear-and-covariance-data-adjustment-for-improvement-of-nuclear-data-files
    (page 10, equation (4.5))

    Parameters
    ----------
    S : 2D iterable
        2D sensitivity of the model y=f(x) (MXN).
    Vx_prior : 2D iterable
        2D covariance matrix of x_prior (MXN).
    Vy_extra : 2D iterable
        2D covariance matrix for y_extra (MXN).

    Returns
    -------
    `pd.Series`
        Ishikawa factor.

    Example
    -------
    >>> ishikawa_factor(S, Vx_prior, Vy_extra)
    0   0.00000e+00
    1   3.00000e+00
    2   8.00000e+00
    dtype: float64
    """
    Vx_prior_ = sandy.CategoryCov(Vx_prior)
    Vy_calc = Vx_prior_._gls_Vy_calc(S)
    Vy_values = np.diag(Vy_calc)
    Vy_extra_ = np.diag(pd.DataFrame(Vy_extra).values)
    index = pd.DataFrame(Vy_extra).index
    return pd.Series(Vy_values / Vy_extra_, index=index)


def constrained_ls_update(x_prior, S, Vx_prior, threshold=None):
    """
    Perform Constrained Least-Squares update for a given variances, vectors
    and sensitivity.

    Parameters
    ----------
    x_prior: 1D iterable
        Vector in which we are going to apply GLS (MX1)
    S : 2D iterable
        2D sensitivity of the model y=f(x) (MXN).
    Vx_prior : 2D iterable
        2D covariance matrix of x_prior (MXN).
    threshold : `int`, optional
            Thereshold to avoid numerical fluctuations. The default is None.

    Returns
    -------
    x_post : `pandas.Series`
        Constrained Least-Squares apply to a vector x_prior given S and
        Vx_prior.

    Example
    -------
    >>> S = np.array([[1, 2], [3, 4]])
    >>> Vx_prior = sandy.CategoryCov.from_var([1, 1]).data
    >>> x_prior = np.array([1, 2])
    >>> constrained_ls_update(x_prior, S, Vx_prior)
    0   -4.00000e+00
    1    5.50000e+00
    dtype: float64
    """
    x_prior_ = pd.Series(x_prior)
    S_ = pd.DataFrame(S).loc[:, x_prior_.index]
    # Data in a appropriate format:
    delta = _y_calc(x_prior, S_.T) - x_prior.dot(S_.T)
    delta = delta.reindex(x_prior_.index)
    Vx_prior_ = sandy.CategoryCov(Vx_prior)
    # Constrained Least-Square sensitivity
    A = Vx_prior_._constrained_ls_sensitivity(S, threshold=threshold).values
    x_post = x_prior_ + delta.dot(A)
    return x_post
