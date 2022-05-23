# -*- coding: utf-8 -*-
"""
This module contains all classes and functions specific for processing GLS and
assessment of adjustment in matrix form.
"""
import pandas as pd
import numpy as np
import scipy.sparse as sps
import sandy

__author__ = "Aitor Bengoechea"
__all__ = [
        "gls_update",
        "gls_cov_update",
        "constrained_gls_update",
        "_y_calc",
        "chi_individual",
        "chi_diag",
        "chi_square",
        "ishikawa_factor",
        ]

x_prior = [1, 2, 3]
y_extra = pd.Series([2, 3, 4], index=[1, 2, 3])
S = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
Vx_prior = [[0, 0, 0], [0, 3, 0], [0, 0, 8]]
Vy_extra = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
N_e = 1


def gls_update(x_prior, S, Vx_prior, Vy_extra, y_extra,
               threshold=None):
    """
    Perform the GlS update of a prior vector, given its prior covariance
    matrix, additional info on the model observable and their covariance
    matrix.
    .. math::
        $$
        x_{post} = x_{prior} + V_{x_{prior}}\cdot S^T \cdot \left(S\cdot V_{x_{prior}}\cdot S^T + V_{y_{extra}}\right)^{-1} \cdot \left(y_{extra} - y_{calc}\right)
        $$

    Parameters
    ----------
    x_prior : 1D iterable
        Vector to be updated (NX1)
    Vx_prior : 2D iterable
        2D covariance matrix of x_prior (NXN).
    Vy_extra : 2D iterable or sigle element 1D iterable
        covariance matrix with the uncertainties of the extra information,
        (MXM) or (1x1).
    S : 2D or 1D iterable
         Sensitivity matrix (MXN) or sensitivity vector(1xN).
    y_extra : 1D iterable
        1D extra info on output (MX1).
    threshold : `int`, optional, default is None
        Thereshold to avoid numerical fluctuations.

    Returns
    -------
    `pd.Series`
        updated vector adjusted with the GLS technique.

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

    >>> S = [1, 2]
    >>> x_prior = pd.Series([1, 1])
    >>> Vx_prior = sandy.CategoryCov.from_var([1, 1]).data
    >>> Vy_extra = [1]
    >>> y_extra = [2]
    >>> gls_update(x_prior, S, Vx_prior, Vy_extra, y_extra)
    0   8.33333e-01
    1   6.66667e-01
    dtype: float64
    """
    # Put data in a appropiate format
    S_ = pd.DataFrame(S).values
    if S_.shape[1] == 1:
        S_ = S_.T
    x_prior_ = pd.Series(x_prior).values
    Vx_prior_ = pd.DataFrame(Vx_prior).values
    index = pd.Series(x_prior).index
    y_extra_ = pd.Series(y_extra).values
    Vy_extra = pd.DataFrame(Vy_extra).values
    # GLS update
    G_inv = _gls_G_inv(Vx_prior_, S_, Vy_extra, threshold=threshold).values
    A = Vx_prior_.dot(S_.T).dot(G_inv)
    y_calc = S_.dot(x_prior_)
    delta = y_extra_ - y_calc
    x_post = x_prior_ + A.dot(delta)
    x_post = pd.Series(x_post, index=index)
    if threshold is not None:
        x_post[abs(x_post) < threshold] = 0
    return x_post


def gls_cov_update(Vx_prior, S, Vy_extra=None, threshold=None):
    """
    Perform GlS update for a given covariance matrix, sensitivity and
    covariance matrix of the extra information, according with
    https://www.tandfonline.com/action/journalInformation?journalCode=tnst20:

    .. math::
        $$
        V_{x_{post}} = V_{x_{prior}} - V_{x_{prior}}\cdot S^T \cdot \left(S\cdot V_{x_{prior}}\cdot S^T + V_{y_{extra}}\right)^{-1} \cdot S \cdot V_{x_{prior}}
        $$

    Parameters
    ----------
    Vx_prior : 2D iterable
        prior covariance matrix to be updated (NxN).
    Vy_extra : 2D iterable or sigle element 1D iterable
        covariance matrix with the uncertainties of the extra information,
        (MXM) or (1x1).
    S : 2D or 1D iterable
        Sensitivity matrix (MXN) or sensitivity vector(1xN).
    threshold : `int`, optional
        Thereshold to avoid numerical fluctuations. The default is None.

    Returns
    -------
    `pd.DataFrame`
       updated covariance matrix adjusted with the GLS technique.

    Notes
    -----
    .. note:: If Vy_extra=None the constraint GLS update technique
    will be performed

    Example
    -------
    >>> S = np.array([[1, 2], [3, 4]])
    >>> cov = sandy.CategoryCov.from_var([1, 1]).data
    >>> Vy = np.diag(pd.Series([1, 1]))
    >>> gls_cov_update(cov, S, Vy)
                 0            1
    0  6.00000e-01 -4.00000e-01
    1 -4.00000e-01  3.14286e-01

    >>> S = np.array([1, 2])
    >>> cov = sandy.CategoryCov.from_var([1, 1]).data
    >>> Vy = [1]
    >>> gls_cov_update(cov, S, Vy)
                 0            1
    0  8.33333e-01 -3.33333e-01
    1 -3.33333e-01  3.33333e-01
    """
    Vx_prior = pd.DataFrame(Vx_prior)
    index, columns = Vx_prior.index, Vx_prior.columns
    s_ = pd.DataFrame(S).values
    G_inv = _gls_G_inv(Vx_prior, s_, Vy_extra=Vy_extra, threshold=threshold).values
    # gls update
    if pd.DataFrame(S).shape[1] == 1:
        A = np.outer(Vx_prior.values.dot(s_) * G_inv, s_)
    else:
        A = Vx_prior.values.dot(s_.T).dot(G_inv).dot(s_)
    diff = A.dot(Vx_prior)
    Vx_post = Vx_prior - diff
    if threshold is not None:
        Vx_post[abs(Vx_post) < threshold] = 0
    return pd.DataFrame(Vx_post, index=index, columns=columns)


def constrained_gls_update(x_prior, y_constraint, S, Vx_prior, threshold=None):
    """
    Perform Constrained Least-Squares update for a vector to be updated,
    its covariance matrix, constraint to be introduced and sensitivity:
    .. math::
        $$
        x_{post}^T = x_{prior}^T + \left(y_{constraint}^T - x_{prior}^T \cdot S^T\right) \cdot \left(S\cdot V_{x_{prior}}\cdot S^T \right)^{-1} \cdot S \cdot V_{x_{prior}}
        $$

    Parameters
    ----------
    x_prior : 1D iterable
        Vector to be updated (NX1)
    Vx_prior : 2D iterable
        2D covariance matrix of x_prior (NXN).
    Vy_extra : 2D iterable or sigle element 1D iterable
        covariance matrix with the uncertainties of the extra information,
        (MXM) or (1x1).
     S : 2D or 1D iterable
         Sensitivity matrix (MXN) or sensitivity vector (1xN).
    y_constraint : 1D iterable
        information that can be used to reﬁne the ﬁt (MX1).
    threshold : `int`, optional, default is None
        Thereshold to avoid numerical fluctuations.

    Returns
    -------
    x_post : `pandas.Series`
        updated vector adjusted with the constrained GLS technique.

    Example
    -------
    >>> S = [[1, 2], [3, 4]]
    >>> x_prior = pd.Series([1, 1])
    >>> Vx_prior = sandy.CategoryCov.from_var([1, 1]).data
    >>> y_constraint = [2, 2]
    >>> constrained_gls_update(x_prior, y_constraint, S, Vx_prior)
    0   -2.00000e+00
    1    2.00000e+00
    dtype: float64

    >>> S = [[1, 2], [3, 4], [1, 0]]
    >>> x_prior = pd.Series([1, 1])
    >>> Vx_prior = sandy.CategoryCov.from_var([1, 1]).data
    >>> y_constraint = [2, 2, 2]
    >>> constrained_gls_update(x_prior, y_constraint, S, Vx_prior)
    0    2.50000e-01
    1   -3.00000e+00
    dtype: float64

    >>> S = [1, 2]
    >>> x_prior = pd.Series([1, 1])
    >>> Vx_prior = sandy.CategoryCov.from_var([1, 1]).data
    >>> y_constraint = [2]
    >>> constrained_gls_update(x_prior, y_constraint, S, Vx_prior)
    0   8.00000e-01
    1   6.00000e-01
    dtype: float64
    """
    # Data in a appropriate format
    x_prior_ = pd.Series(x_prior)
    index = x_prior_.index
    s_ = pd.DataFrame(S).values
    Vx_prior_ = pd.DataFrame(Vx_prior).values
    y_constraint_ = pd.Series(y_constraint).values
    if s_.shape[1] == 1:
        s_ = s_.T
    # Constrained gls calculations
    y_calc = s_.dot(x_prior_.values)
    G_inv = _gls_G_inv(Vx_prior, s_, Vy_extra=None, threshold=threshold).values
    delta = y_constraint_ - y_calc
    x_post = x_prior_.values + delta.dot(G_inv).dot(s_).dot(Vx_prior)
    x_post = pd.Series(x_post, index=index)
    if threshold is not None:
        x_post[abs(x_post) < threshold] = 0
    return x_post


def _gls_G_inv(Vx_prior, s, Vy_extra=None, threshold=None):
    """
    Compute part of the GLS update technique. Output calculated using
    .. math::
        $$
        (S\cdot V_{x_{prior}}\cdot S^T + V_{y_{extra}})^{-1}
        $$

    Parameters
    ----------
    Vx_prior : 2D iterable
        prior covariance matrix to be updated (NxN).
    Vy_extra : 2D iterable or sigle element 1D iterable
        covariance matrix with the uncertainties of the extra information,
        (MXM) or (1x1).
    S : 2D or 1D iterable
        Sensitivity matrix (MXN) or sensitivity vector(1xN).
    threshold : `int`, optional
        Thereshold to avoid numerical fluctuations. The default is None.

    Returns
    -------
    `pandas.DataFrame`
        matrix calculated using the formula inserted in the description

    Notes
    -----
    .. note:: If the matrix computed with $S\cdot V_{x_{prior}}\cdot S^T + V_{y_{extra}}$
    has zero columns or vectors, it will result singular. In this case
    the invertion is performed considering a reduced matrix; then the size of
    the matrix is restored introducing the null values.

    Example
    -------
    >>> S = np.array([[1, 2], [3, 4]])
    >>> cov = sandy.CategoryCov.from_var([1, 1]).data
    >>> Vy = np.diag(pd.Series([1, 1]))
    >>> _gls_G_inv(cov, S, Vy)
                 0            1
    0  7.42857e-01 -3.14286e-01
    1 -3.14286e-01  1.71429e-01

    >>> _gls_G_inv(cov, S)
                 0            1
    0  6.25000e+00 -2.75000e+00
    1 -2.75000e+00  1.25000e+00
    
    >>> S = np.array([1, 2])
    >>> cov = sandy.CategoryCov.from_var([1, 1]).data
    >>> Vy = [1]
    >>> _gls_G_inv(cov, S, Vy)
                0
    0 1.66667e-01
    """
    # GLS_sensitivity:
    cov_ = pd.DataFrame(Vx_prior)
    s_ = pd.DataFrame(s)
    if s_.shape[1]==1:
        s_ = s_.T
    Vy_calc = sandwich(cov_.values, s_.T.values, threshold=threshold)
    if Vy_extra is not None:
        Vy_extra_ = pd.DataFrame(Vy_extra)
        index = Vy_extra_.index
        G = Vy_calc.values + Vy_extra_.values
    else:
        index = Vy_calc.index
        G = Vy_calc
    M_nonzero_idxs, M_reduce = reduce_size(G)
    G_inv = np.linalg.inv(M_reduce)
    G_inv = restore_size(M_nonzero_idxs, G_inv, len(G))
    return pd.DataFrame(G_inv)


def sandwich(cov, s, threshold=None):
    """
    Apply the "sandwich formula" to the covariance matrix passed for a given
    sensitivity. According with http://dx.doi.org/10.1016/j.anucene.2015.10.027,
    the moment propagation equation is implemented as:

       .. math::
           $$
           V_R = S^T\cdot V_P\cdot S
           $$

    Parameters
    ----------
    cov : 2D iterable
        Parameter covariance matrix (NxN)
    s : 1D or 2D iterable
        General sensitivities (Nx1) or (NxM)
    threshold : `int`, optional, default is None
        Thereshold to avoid numerical fluctuations.

    Returns
    -------
    `pandas.DataFrame`
        response covariance matrix obtained with the sandwich formula.

    Examples
    --------
    >>> var = np.array([1, 2, 3])
    >>> s = pd.Series([1, 2, 3])
    >>> cov = sandy.CategoryCov.from_var(var).data
    >>> sandwich(cov, s)
                0
    0 3.60000e+01

    >>> s = np.array([1, 2, 3])
    >>> var = pd.Series([1, 2, 3])
    >>> cov = sandy.CategoryCov.from_var(var).data
    >>> var = sandy.CategoryCov.from_var(s).data
    >>> sandwich(cov, var)
                0           1           2
    0 1.00000e+00 0.00000e+00 0.00000e+00
    1 0.00000e+00 8.00000e+00 0.00000e+00
    2 0.00000e+00 0.00000e+00 2.70000e+01
    """
    s_ = pd.DataFrame(s)
    index = s_.columns
    cov_ = pd.DataFrame(cov)
    Vy_calc = s_.T.values.dot(cov_.values).dot(s_.values)
    sandwich = pd.DataFrame(Vy_calc, index=index, columns=index)
    if threshold is not None:
        sandwich[sandwich < threshold] = 0
    return sandwich


def reduce_size(data):
    """
    Reduces the size of the matrix, erasing the zero values.

    Parameters
    ----------
    data : 'pd.DataFrame'
        Matrix to be reduced.

    Returns
    -------
    nonzero_idxs : `numpy.ndarray`
        The indices of the diagonal that are not null.
    cov_reduced : `pandas.DataFrame`
        The reduced matrix.

    Examples
    --------
    >>> S = pd.DataFrame(np.diag(np.array([1, 2, 3])))
    >>> non_zero_index, reduce_matrix = reduce_size(S)
    >>> assert reduce_matrix.equals(S)
    >>> assert (non_zero_index == range(3)).all()

    >>> S = pd.DataFrame(np.diag(np.array([0, 2, 3])))
    >>> non_zero_index, reduce_matrix = reduce_size(S)
    >>> assert (non_zero_index == np.array([1, 2])).all()
    >>> reduce_matrix
      1 2
    1 2 0
    2 0 3

    >>> S.index = S.columns = ["a", "b", "c"]
    >>> non_zero_index, reduce_matrix = reduce_size(S)
    >>> reduce_matrix
      b c
    b 2 0
    c 0 3
    """
    data_ = pd.DataFrame(data)
    nonzero_idxs = np.flatnonzero(np.diag(data_))
    cov_reduced = data_.iloc[nonzero_idxs, nonzero_idxs]
    return nonzero_idxs, cov_reduced


def restore_size(nonzero_idxs, mat_reduced, dim):
    """
    Restore the size of a matrix.

    Parameters
    ----------
    nonzero_idxs : `numpy.ndarray`
        The indices of the diagonal that are not null.
    mat_reduced : `numpy.ndarray`
        The reduced matrix.
    dim : `int`
        Dimension of the original matrix.

    Returns
    -------
    mat : `pd.DataFrame`
        Matrix of specified dimensions.

    Notes
    -----
    ..notes:: This funtion was developed to be used after using
              `reduce_size`.

    Examples
    --------
    >>> S = pd.DataFrame(np.diag(np.array([0, 2, 3, 0])))
    >>> M_nonzero_idxs, M_reduce = reduce_size(S)
    >>> M_reduce[::] = 1
    >>> restore_size(M_nonzero_idxs, M_reduce.values, len(S))
                0           1           2           3
    0 0.00000e+00 0.00000e+00 0.00000e+00 0.00000e+00
    1 0.00000e+00 1.00000e+00 1.00000e+00 0.00000e+00
    2 0.00000e+00 1.00000e+00 1.00000e+00 0.00000e+00
    3 0.00000e+00 0.00000e+00 0.00000e+00 0.00000e+00

    >>> S = pd.DataFrame(np.diag(np.array([0, 2, 3, 0])), index=[1, 2, 3, 4], columns=[5, 6, 7, 8])
    >>> M_nonzero_idxs, M_reduce = reduce_size(S)
    >>> M_reduce[::] = 1
    >>> restore_size(M_nonzero_idxs, M_reduce.values, len(S))
                0           1           2           3
    0 0.00000e+00 0.00000e+00 0.00000e+00 0.00000e+00
    1 0.00000e+00 1.00000e+00 1.00000e+00 0.00000e+00
    2 0.00000e+00 1.00000e+00 1.00000e+00 0.00000e+00
    3 0.00000e+00 0.00000e+00 0.00000e+00 0.00000e+00
    """
    mat = np.zeros((dim, dim))
    for i, ni in enumerate(nonzero_idxs):
        mat[ni, nonzero_idxs] = mat_reduced[i]
    return pd.DataFrame(mat)


def _y_calc(x_prior, S):
    """
    Perform model calculation in GLS.

    Parameters
    ----------
    x_prior: 1D iterable
        Vector in which we are going to apply GLS (MX1).
    S : 2D iterable
        2D sensitivity of the model y=f(x) (NXM).
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
    idx = S_.index
    S_ = sps.csr_matrix(S_.values)
    y_calc = S_.dot(x_prior_.values)
    y_calc = pd.Series(y_calc, index=idx)
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
    G = sandwich(Vx_prior_.data, S_.T,)
    if Vy_extra is not None:
        Vy_extra_ = pd.DataFrame(Vy_extra)
        G = G + Vy_extra_
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
    G_inv = _gls_G_inv(Vx_prior_.data, S_, Vy_extra=Vy_extra).values
    # Chi diagonal calculations
    G_inv = np.sqrt(np.diag(G_inv))
    delta = np.abs(y_extra_.values - y_calc_)
    return pd.Series(delta / G_inv, index=y_extra_.index)


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
    G_inv = _gls_G_inv(Vx_prior_.data, S_, Vy_extra).values
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
    Vy_calc_ = sandwich(Vx_prior_.data, S_.T)
    Vy_calc_ = Vy_calc_.reindex(index=index, columns=index).fillna(0)
    # Ishikawa factor calculations
    Vy_values = np.diag(Vy_calc_)
    Vy_extra_ = np.diag(Vy_extra_.values)
    return pd.Series(Vy_values / Vy_extra_, index=index)

