# -*- coding: utf-8 -*-
"""
This module contains all classes and functions specific for processing GLS and
Assessment of adjustment.
"""
import pandas as pd
import numpy as np
import scipy.sparse as sps
import sandy

__author__ = "Aitor Bengoechea"
__all__ = [
        "gls_update",
        "_y_calc",
        "chi_individual",
        "chi_diag",
        "chi_square",
        "ishikawa_factor",
        "constrained_gls_update"
        ]

x_prior = [1, 2, 3]
y_extra = pd.Series([2, 3, 4], index=[1, 2, 3])
S = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
Vx_prior = [[0, 0, 0], [0, 3, 0], [0, 0, 8]]
Vy_extra = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
N_e = 1


def gls_update(x_prior, S, Vx_prior, Vy_extra, y_extra, sparse=False,
               rows=None, threshold=None):
    """
    Perform the GlS update of a prior vector, given its prior covariance
    matrix, a lekelyhood matrix and additional info on the model obserbale
    (both values and covariance matrix).
    .. math::
        $$
        x_{post} = x_{prior} + V_{x_{prior}}\cdot S.T \cdot \left(S\cdot V_{x_{prior}}\cdot S.T + V_{y_{extra}}\right)^{-1} \cdot \left(y_{extra} - y_{calc}\right)
        $$

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
    rows : `int`, optional
        Option to use row calculation for matrix calculations. This option
        defines the number of lines to be taken into account in each loop.
        The default is None.
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
    >>> x_prior = pd.Series([1, 1])
    >>> Vx_prior = sandy.CategoryCov.from_var([1, 1]).data
    >>> Vy_extra = pd.DataFrame([[1, 0], [0, 1]])
    >>> y_extra = [2, 2]
    >>> gls_update(x_prior, S, Vx_prior, Vy_extra, y_extra)
    0   2.00000e-01
    1   4.85714e-01
    dtype: float64

    >>> gls_update(x_prior, S, Vx_prior, Vy_extra, y_extra, rows=1)
    0   2.00000e-01
    1   4.85714e-01
    dtype: float64
    """
    # Put data in a appropiate format
    S_ = pd.DataFrame(S)
    x_prior_ = pd.Series(x_prior)
    Vx_prior_ = sandy.CategoryCov(Vx_prior).data
    y_extra_ = pd.Series(y_extra)
    union_index = S_.columns.union(x_prior_.index)
    x_prior_ = x_prior_.reindex(union_index).fillna(0)
    S_ = S_.reindex(columns=union_index).fillna(0)
    Vx_prior_ = Vx_prior_.reindex(index=union_index, columns=union_index).fillna(0)
    Vx_prior_ = sandy.CategoryCov(Vx_prior_)
    # Model calculus:
    y_calc_ = _y_calc(x_prior, S_)
    # Fix model calculus and extra information
    y_calc_ = y_calc_.reindex(y_extra_.index).fillna(0)
    S_ = S_.reindex(index=y_extra_.index).fillna(0)
    delta = y_extra_ - y_calc_
    # GLS update
    A = Vx_prior_._gls_general_sensitivity(S_, Vy_extra, rows=rows,
                                           threshold=threshold).values
    A = sps.csr_matrix(A)
    index = x_prior_.index
    x_prior_ = x_prior_.values
    x_post = x_prior_ + A.dot(delta)
    x_post = pd.Series(x_post, index=index)
    if threshold is not None:
        x_post[abs(x_post) < threshold] = 0
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
    sparse : `bool`, optional
        Option to use sparse matrix for calculations. The default is False.

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
    x_prior_ = pd.Series(x_prior)
    union_index = S_.columns.union(x_prior_.index)
    S_ = S_.reindex(columns=union_index).fillna(0)
    x_prior_ = x_prior_.reindex(union_index).fillna(0)
    index = S_.index
    S_ = sps.csr_matrix(S_.values)
    y_calc = S_.dot(x_prior_.values)
    y_calc = pd.Series(y_calc, index=index)
    return y_calc


def chi_individual(x_prior, S, Vx_prior, Vy_extra, y_extra,
                   rows=None):
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
    rows : `int`, optional
        Option to use row calculation for matrix calculations. This option
        defines the number of lines to be taken into account in each loop.
        The default is None.

    Returns
    -------
    `pd.Series`
        individual chi-value measured in sigmas.

    Results:
    -------
    chi individual >> 1 :
        Inconsistency may exist between |y_extra - y_calc| and covariance
        matrix, S*Vx_prior*S.T, and Vy_extra.

    Example
    -------
    >>> x_prior = [1, 2, 3]
    >>> y_extra = pd.Series([2, 3, 4], index=[1, 2, 3])
    >>> S = pd.DataFrame([[1, 0, 0], [0, 1, 0], [0, 0, 1]], index=[1, 2, 3])
    >>> Vx_prior = [[0, 0, 0], [0, 3, 0], [0, 0, 8]]
    >>> Vy_extra = pd.DataFrame([[1, 0, 0], [0, 1, 0], [0, 0, 1]], index=[1, 2, 3], columns=[1, 2, 3])
    >>> chi_individual(x_prior, S, Vx_prior, Vy_extra, y_extra)
    1   1.00000e+00
    2   5.00000e-01
    3   3.33333e-01
    dtype: float64

    >>> chi_individual(x_prior, S, Vx_prior, Vy_extra, y_extra, rows=1)
    1   1.00000e+00
    2   5.00000e-01
    3   3.33333e-01
    dtype: float64
    """
    # Data in a appropriate format
    Vx_prior_ = sandy.CategoryCov(Vx_prior).data
    y_extra_ = pd.Series(y_extra)
    S_ = pd.DataFrame(S)
    union_index = S_.columns.union(Vx_prior_.index)
    S_ = S_.reindex(columns=union_index).fillna(0)
    Vx_prior_ = sandy.CategoryCov(Vx_prior_.reindex(columns=union_index,
                                                    index=union_index).fillna(0))
    # Model calculations
    y_calc_ = _y_calc(x_prior, S_)
    S_ = S_.reindex(index=y_extra_.index).fillna(0)
    y_calc_ = y_calc_.reindex(index=y_extra_.index).fillna(0).values
    # Chi individual calculations
    G = Vx_prior_._gls_G(S_, Vy_extra, rows=rows)
    G = np.sqrt(np.diag(G))
    delta = np.abs(y_extra_.values - y_calc_)
    return pd.Series(delta / G, index=y_extra_.index)


def chi_diag(x_prior, S, Vx_prior, Vy_extra, y_extra, rows=None):
    """
    Function to calculate diagonal chi-value
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
    rows : `int`, optional
        Option to use row calculation for matrix calculations. This option
        defines the number of lines to be taken into account in each loop.
        The default is None.

    Returns
    -------
    `pd.Series`
        diagonal chi-value 

    Results:
    -------
    chi diagonal >> 1 :
        Inconsistency may exist between |y_extra - y_calc| and covariance
        matrix, S*Vx_prior*S.T, and Vy_extra.

    Example
    -------
    >>> x_prior = [1, 2, 3]
    >>> y_extra = pd.Series([2, 3, 4], index=[1, 2, 3])
    >>> S = pd.DataFrame([[1, 0, 0], [0, 1, 0], [0, 0, 1]], index=[1, 2, 3])
    >>> Vx_prior = [[0, 0, 0], [0, 3, 0], [0, 0, 8]]
    >>> Vy_extra = pd.DataFrame([[1, 0, 0], [0, 1, 0], [0, 0, 1]], index=[1, 2, 3], columns=[1, 2, 3])
    >>> N_e = 1
    >>> chi_diag(x_prior, S, Vx_prior, Vy_extra, y_extra)
    1   1.00000e+00
    2   2.00000e+00
    3   3.00000e+00
    dtype: float64

    >>> chi_diag(x_prior, S, Vx_prior, Vy_extra, y_extra, rows=1)
    1   1.00000e+00
    2   2.00000e+00
    3   3.00000e+00
    dtype: float64
    """
    # Data in a appropriate format
    Vx_prior_ = sandy.CategoryCov(Vx_prior).data
    y_extra_ = pd.Series(y_extra)
    S_ = pd.DataFrame(S)
    union_index = S_.columns.union(Vx_prior_.index)
    S_ = S_.reindex(columns=union_index).fillna(0)
    Vx_prior_ = sandy.CategoryCov(Vx_prior_.reindex(columns=union_index,
                                                    index=union_index).fillna(0))
    # Model calculations
    y_calc_ = _y_calc(x_prior, S_)
    S_ = S_.reindex(index=y_extra_.index).fillna(0)
    y_calc_ = y_calc_.reindex(index=y_extra_.index).fillna(0).values
    G_inv = Vx_prior_._gls_G_inv(S_, Vy_extra=Vy_extra, rows=rows).values
    # Chi diagonal calculations
    G_inv = np.sqrt(np.diag(G_inv))
    delta = np.abs(y_extra_.values - y_calc_)
    return pd.Series(delta / G_inv, index=y_extra_.index)


def chi_square(x_prior, S, Vx_prior, Vy_extra, y_extra, N_e, rows=None):
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
    rows : `int`, optional
        Option to use row calculation for matrix calculations. This option
        defines the number of lines to be taken into account in each loop.
        The default is None.

    Returns
    -------
    `pd.Series`
        contribution to chi-square value

    Results:
    -------
    chi square < 0 :
        The experiment is very effective in the adjustment.

    Example
    -------
    >>> x_prior = [1, 2, 3]
    >>> y_extra = pd.Series([2, 3, 4], index=[1, 2, 3])
    >>> S = pd.DataFrame([[1, 0, 0], [0, 1, 0], [0, 0, 1]], index=[1, 2, 3])
    >>> Vx_prior = [[0, 0, 0], [0, 3, 0], [0, 0, 8]]
    >>> Vy_extra = pd.DataFrame([[1, 0, 0], [0, 1, 0], [0, 0, 1]], index=[1, 2, 3], columns=[1, 2, 3])
    >>> chi_square(x_prior, S, Vx_prior, Vy_extra, y_extra, N_e)
    1   1.00000e+00
    2   2.50000e-01
    3   1.11111e-01
    dtype: float64

    >>> chi_square(x_prior, S, Vx_prior, Vy_extra, y_extra, N_e, rows=1)
    1   1.00000e+00
    2   2.50000e-01
    3   1.11111e-01
    dtype: float64
    """
    # Data in a appropriate format
    Vx_prior_ = sandy.CategoryCov(Vx_prior).data
    y_extra_ = pd.Series(y_extra)
    S_ = pd.DataFrame(S)
    union_index = S_.columns.union(Vx_prior_.index)
    S_ = S_.reindex(columns=union_index).fillna(0)
    Vx_prior_ = sandy.CategoryCov(Vx_prior_.reindex(columns=union_index,
                                                    index=union_index).fillna(0))
    # Model calculations
    y_calc_ = _y_calc(x_prior, S_)
    S_ = S_.reindex(index=y_extra_.index).fillna(0)
    y_calc_ = y_calc_.reindex(index=y_extra_.index).fillna(0).values
    G_inv = Vx_prior_._gls_G_inv(S_, Vy_extra, rows=rows).values
    # Chi square calculations
    delta = y_extra_.values - y_calc_
    chi_square = delta.T.dot(G_inv) * delta / N_e
    return pd.Series(chi_square, index=y_extra_.index)


def ishikawa_factor(S, Vx_prior, Vy_extra, rows=None):
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
    rows : `int`, optional
        Option to use row calculation for matrix calculations. This option
        defines the number of lines to be taken into account in each loop.
        The default is None.

    Returns
    -------
    `pd.Series`
        Ishikawa factor.

    Results:
    -------
    Ishikawa factor << 1 :
        The extra data is not so useful and the data remain unchanged.
    Ishikawa factor >> 1 :
        The extra data very useful and the 'posteriori' covariance will
        be reduced to the same level as the integral parameter covariance.
    Ishikawa factor ~ 1 :
        The experiment is useful and the 'posteriori'  covariance will be
        reduced by approximately half

    Example
    -------
    >>> ishikawa_factor(S, Vx_prior, Vy_extra)
    0   0.00000e+00
    1   3.00000e+00
    2   8.00000e+00
    dtype: float64
    """
    # Data in a appropriate format
    Vx_prior_ = sandy.CategoryCov(Vx_prior).data
    Vy_extra_ = pd.DataFrame(Vy_extra)
    index = Vy_extra_.index
    S_ = pd.DataFrame(S)
    union_index = S_.columns.union(Vx_prior_.index)
    S_ = S_.reindex(columns=union_index).fillna(0)
    Vx_prior_ = sandy.CategoryCov(Vx_prior_.reindex(columns=union_index,
                                                    index=union_index).fillna(0))
    # Model calculations
    Vy_calc_ = Vx_prior_._gls_Vy_calc(S_, rows=rows)
    Vy_calc_ = Vy_calc_.reindex(index=index, columns=index).fillna(0)
    # Ishikawa factor calculations
    Vy_values = np.diag(Vy_calc_)
    Vy_extra_ = np.diag(Vy_extra_.values)
    return pd.Series(Vy_values / Vy_extra_, index=index)


def constrained_gls_update(x_prior, S, Vx_prior, rows=None, threshold=None):
    """
    Perform Constrained Least-Squares update for a given variances, vectors
    and sensitivity:
    .. math::
        $$
        x_{post} = x_{prior} + \left(S.T \cdot x_{prior} - x_{prior} \cdot S.T\right) \cdot \left(S\cdot V_{x_{prior}}\cdot S.T + V_{y_{extra}}\right)^{-1} \cdot S \cdot V_{x_{prior}}
        $$

    Parameters
    ----------
    x_prior: 1D iterable
        Vector in which we are going to apply GLS (MX1)
    S : 2D iterable
        2D sensitivity of the model y=f(x) (MXN).
    Vx_prior : 2D iterable
        2D covariance matrix of x_prior (MXN).
    rows : `int`, optional
        Option to use row calculation for matrix calculations. This option
        defines the number of lines to be taken into account in each loop.
        The default is None.
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
    >>> constrained_gls_update(x_prior, S, Vx_prior)
    0   -4.00000e+00
    1    5.50000e+00
    dtype: float64

    >>> constrained_gls_update(x_prior, S, Vx_prior, rows=1)
    0   -4.00000e+00
    1    5.50000e+00
    dtype: float64
    """
    # Data in a appropriate format
    x_prior_ = pd.Series(x_prior)
    original_index = x_prior_.index
    S_ = pd.DataFrame(S).reindex(columns=x_prior_.index).fillna(0)
    Vx_prior_ = sandy.CategoryCov(Vx_prior).data
    index = S_.index
    union_index = S_.columns.union(Vx_prior_.index)
    S_ = S_.reindex(columns=union_index).fillna(0)
    Vx_prior_ = sandy.CategoryCov(Vx_prior_.reindex(columns=union_index,
                                                    index=union_index).fillna(0))
    x_prior_ = x_prior_.reindex(union_index).fillna(0)
    # Common Model calculation for sparse and no sparse:
    y_calc = _y_calc(x_prior, S_.T).reindex(index).fillna(0)
    A = Vx_prior_._gls_constrained_sensitivity(S_, rows=rows,
                                               threshold=threshold)
    # Constrained gls calculations
    x_prior_sps = sps.coo_matrix(x_prior_)
    S_sps = sps.csc_matrix(S_)
    diff = x_prior_sps.dot(S_sps.T)
    diff = pd.Series(diff.toarray()[0], index=index)   
    delta = y_calc - diff
    delta = sps.coo_matrix(delta)
    A = sps.csc_matrix(A.values)
    x_post = x_prior_sps + delta.dot(A)
    x_post = pd.Series(x_post.toarray()[0], index=original_index)
    if threshold is not None:
        x_post[abs(x_post) < threshold] = 0
    return x_post
