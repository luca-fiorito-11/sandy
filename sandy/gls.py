# -*- coding: utf-8 -*-
"""
This module contains all classes and functions specific for processing GLS and
assessment of adjustment. It handles only numpy.ndarray.
"""
import numpy as np
import logging

__author__ = "Aitor Bengoechea"
__all__ = [
        "gls_update",
        "_gls_parameters_update",
        "_gls_cov_update",
        "chi_individual",
        "chi_diag",
        "chi_square",
        "ishikawa_factor",
        ]


def gls_update(x_prior, S, Vx_prior, y_extra, Vy_extra=None):
    r"""
    Perform the GlS update of a prior vector and its related covariance matrix,
    according with
    https://www.tandfonline.com/action/journalInformation?journalCode=tnst2
    .. math::
        $$
        x_{post} = x_{prior} + V_{x_{prior}}\cdot S^T \cdot \left(S\cdot V_{x_{prior}}\cdot S^T + V_{y_{extra}}\right)^{-1} \cdot \left(y_{extra} - y_{calc}\right)\\
        V_{x_{post}} = V_{x_{prior}} - V_{x_{prior}}\cdot S^T \cdot \left(S\cdot V_{x_{prior}}\cdot S^T + V_{y_{extra}}\right)^{-1} \cdot S \cdot V_{x_{prior}}
        $$

    Parameters
    ----------
    x_prior : 1D iterable
        Vector to be updated (N,)
    S : 2D or 1D iterable
         Sensitivity matrix (M, N) or sensitivity vector(N,)
    Vx_prior : 2D iterable
        2D covariance matrix of x_prior (N, N)
    y_extra : 1D iterable
        1D extra info on output (M,)
    Vy_extra : 2D iterable or sigle element 1D iterable, optional, default is `None`
        covariance matrix with the uncertainties of the extra information,
        (M, M) or (1,).

    Returns
    -------
    `numpy.ndarray`
        updated vector adjusted with the GLS technique
    `numpy.ndarray`
        updated covariance matrix adjusted with the GLS technique

    Notes
    -----
    .. note:: If Vy_extra=None the constraint GLS update technique
    will be performed

    Example
    -------
    >>> S = [[1, 2], [3, 4]]
    >>> x_prior = [1, 1]
    >>> Vx_prior = np.diag([1, 1])
    >>> Vy_extra = np.array([[1, 0], [0, 1]])
    >>> y_extra = [2, 2]
    >>> gls_update(x_prior, S, Vx_prior, y_extra, Vy_extra)
    (array([0.2       , 0.48571429]),
     array([[ 0.6       , -0.4       ],
            [-0.4       ,  0.31428571]]))

    >>> y_constraint = y_extra
    >>> gls_update(x_prior, S, Vx_prior, y_constraint)[0]
    array([-2.,  2.])
    """
    x_post = _gls_parameters_update(x_prior, S, Vx_prior, y_extra, Vy_extra)
    Vx_post = _gls_cov_update(Vx_prior, S, Vy_extra)
    return x_post, Vx_post


def _gls_parameters_update(x_prior, S, Vx_prior, y_extra, Vy_extra=None):
    r"""
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
        Vector to be updated (N,)
    S : 2D or 1D iterable
         Sensitivity matrix (M, N) or sensitivity vector(N,)
    Vx_prior : 2D iterable
        2D covariance matrix of x_prior (N, N)
    y_extra : 1D iterable
        1D extra info on output (M,)
    Vy_extra : 2D iterable or sigle element 1D iterable, optional, default is `None`
        covariance matrix with the uncertainties of the extra information,
        (M, M) or (1,).

    Returns
    -------
    `numpy.ndarray`
        updated vector adjusted with the GLS technique.

    Notes
    -----
    .. note:: If Vy_extra=None the constraint GLS update technique
    will be performed

    Example
    -------
    >>> S = [[1, 2], [3, 4]]
    >>> x_prior = [1, 1]
    >>> Vx_prior = np.diag([1, 1])
    >>> Vy_extra = np.array([[1, 0], [0, 1]])
    >>> y_extra = [2, 2]
    >>> _gls_parameters_update(x_prior, S, Vx_prior, y_extra, Vy_extra)
    array([0.2       , 0.48571429])

    >>> _gls_parameters_update(x_prior, S, Vx_prior, y_extra, Vy_extra=None)
    array([-2.,  2.])

    >>> S = [1, 2]
    >>> x_prior = [1, 1]
    >>> Vx_prior = np.diag([1, 1])
    >>> Vy_extra = [1]
    >>> y_extra = [2]
    >>> _gls_parameters_update(x_prior, S, Vx_prior, y_extra, Vy_extra)
    array([0.83333333, 0.66666667])

    >>> _gls_parameters_update(x_prior, S, Vx_prior, y_extra, Vy_extra=None)
    array([0.8, 0.6])
    """
    # Put data in a appropiate format
    S_ = np.array(S)
    x_prior_ = np.array(x_prior)
    Vx_prior_ = np.array(Vx_prior)
    y_extra_ = np.array(y_extra)
    G_inv = _gls_G_inv(Vx_prior_, S_, Vy_extra=Vy_extra)
    ndim = len(S_.shape)
    S_ = np.array([S]) if ndim == 1 else S_
    # GLS update
    A = Vx_prior_.dot(S_.T).dot(G_inv)
    y_calc = S_.dot(x_prior_)
    delta = y_extra_ - y_calc
    x_post = x_prior_ + A.dot(delta)
    return x_post


def _gls_cov_update(Vx_prior, S, Vy_extra=None):
    r"""
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
        prior covariance matrix to be updated (N, N).
    Vy_extra : 2D iterable or sigle element 1D iterable
        covariance matrix with the uncertainties of the extra information,
        (M, M) or (1,).
    S : 2D or 1D iterable
        Sensitivity matrix (M, N) or sensitivity vector(N,).

    Returns
    -------
    `numpy.ndarray`
       updated covariance matrix adjusted with the GLS technique.

    Notes
    -----
    .. note:: If Vy_extra=None the constraint GLS update technique
    will be performed

    Example
    -------
    >>> S = np.array([[1, 2], [3, 4]])
    >>> cov = np.diag([1, 1])
    >>> Vy = np.diag([1, 1])
    >>> _gls_cov_update(cov, S, Vy)
    array([[ 0.6       , -0.4       ],
           [-0.4       ,  0.31428571]])

    >>> S = np.array([1, 2])
    >>> Vy = [1]
    >>> _gls_cov_update(cov, S, Vy)
    array([[ 0.83333333, -0.33333333],
           [-0.33333333,  0.33333333]])

    >>> _gls_cov_update(cov, S)
    array([[ 0.8, -0.4],
           [-0.4,  0.2]])
    """
    # Put data in a appropiate format
    Vx_prior_ = np.array(Vx_prior)
    s_ = np.array(S)
    G_inv = _gls_G_inv(Vx_prior_, s_, Vy_extra=Vy_extra)
    # Gls update
    ndim = len(s_.shape)
    if ndim == 1:
        s_ = s_.reshape(1, -1)
    A = Vx_prior_.dot(s_.T).dot(G_inv).dot(s_)
    diff = A.dot(Vx_prior_)
    Vx_post = Vx_prior_ - diff
    return Vx_post


def _gls_G_inv(Vx_prior, s, Vy_extra=None):
    r"""
    Compute part of the GLS update technique. Output calculated using
    .. math::
        $$
        (S\cdot V_{x_{prior}}\cdot S^T + V_{y_{extra}})^{-1}
        $$

    Parameters
    ----------
    Vx_prior : 2D iterable
        prior covariance matrix to be updated (N, N).
    Vy_extra : 2D iterable or sigle element 1D iterable
        covariance matrix with the uncertainties of the extra information,
        (M, M) or (1,).
    S : 2D or 1D iterable
        Sensitivity matrix (M, N) or sensitivity vector(N,).

    Returns
    -------
    `numpy.ndarray`
        matrix calculated using the formula inserted in the description

    Notes
    -----
    .. note:: If the matrix computed with
    $S\cdot V_{x_{prior}}\cdot S^T + V_{y_{extra}}$ has determinant close to
    zero the (Moore-Penrose) pseudo-inverse will be performed.

    Example
    -------
    >>> S = np.array([[1, 2], [3, 4]])
    >>> cov = np.diag([1, 1])
    >>> Vy = np.diag([1, 1])
    >>> _gls_G_inv(cov, S, Vy)
    array([[ 0.74285714, -0.31428571],
           [-0.31428571,  0.17142857]])

    >>> _gls_G_inv(cov, S)
    array([[ 6.25, -2.75],
           [-2.75,  1.25]])

    >>> S = [1, 2]
    >>> Vy = [1]
    >>> _gls_G_inv(cov, S, Vy)
    array([[0.16666667]])

    >>> _gls_G_inv(cov, S)
    array([[0.2]])
    """
    # GLS sensitivity:
    cov_ = np.array(Vx_prior)
    s_ = np.array(s)
    Vy_calc = sandwich(cov_, s_)
    if Vy_extra is not None:
        Vy_extra_ = np.array(Vy_extra)
        G = Vy_calc + Vy_extra_
    else:
        G = Vy_calc
    if len(G.shape) > 1 and abs(np.linalg.det(G)) < 10e-5:
        msg = "determinant of the matrix (S*Vprior*S^T + Vextra) " +\
                f"is {np.linalg.det(G)}, the inverse becomes unreliable"
        logging.warning(msg)
    # Compute the (Moore-Penrose) pseudo-inverse of a matrix.
    # It will converge to the inverse if the matrix G is not singular
    G_inv = np.linalg.pinv(G)
    return G_inv


def sandwich(cov, s):
    r"""
    Apply the "sandwich formula" to the covariance matrix passed for a given
    sensitivity. According with https://doi.org/10.1155/2013/380284,
    the moment propagation equation is implemented as:

       .. math::
           $$
           V_R = S\cdot V_P\cdot S^T
           $$

    Parameters
    ----------
    cov : 2D iterable
        Parameter covariance matrix (N, N)
    s : 1D or 2D iterable
        General sensitivities (N,) or (M, N)

    Returns
    -------
    `numpy.ndarray`
        2D iterable,
        response covariance matrix obtained with the sandwich formula.

    Examples
    --------
    >>> var = [1, 2, 3]
    >>> s = [1, 2, 3]
    >>> cov = np.diag(var)
    >>> sandwich(cov, s)
    array([[36]])

    >>> s_matrix = np.diag(s)
    >>> sandwich(cov, s_matrix)
    array([[ 1,  0,  0],
           [ 0,  8,  0],
           [ 0,  0, 27]])
    """
    s_ = np.asarray(s)
    cov_ = np.asarray(cov)

    sandwich_results = s_ @ cov_ @ s_.T

    if s_.ndim == 1:
        sandwich_results = sandwich_results.reshape(-1, s_.ndim)

    return sandwich_results


def chi_individual(x_prior, S, Vx_prior, Vy_extra, y_extra):
    r"""
    Function to calculate individual chi-value measured in sigmas according to
    https://www.oecd-nea.org/jcms/pl_19760/intermediate-report-on-methods-and-approaches-to-provide-feedback-from-nuclear-and-covariance-data-adjustment-for-improvement-of-nuclear-data-files
    (page 9, equation (4.2))

    Parameters
    ----------
    x_prior: 1D iterable
        Vector in which we are going to apply GLS (M,)
    S : 2D iterable
        2D sensitivity of the model y=f(x) (M, N).
    Vx_prior : 2D iterable
        2D covariance matrix of x_prior (M, N).
    Vy_extra : 2D iterable
        2D covariance matrix for y_extra (M, N).
    y_extra : 1D iterable
        1D extra info on output (N,).

    Returns
    -------
    `numpy.ndarray`
        individual chi-value measured in sigmas.

    Results:
    -------
    chi individual >> 1 :
        Inconsistency may exist between |y_extra - y_calc| and covariance
        matrix, S*Vx_prior*S.T, and Vy_extra.

    Example
    -------
    >>> x_prior = [1, 2, 3]
    >>> y_extra = [2, 3, 4]
    >>> S = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
    >>> Vx_prior = np.array([[0, 0, 0], [0, 3, 0], [0, 0, 8]])
    >>> Vy_extra = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
    >>> chi_individual(x_prior, S, Vx_prior, Vy_extra, y_extra)
    array([1.        , 0.5       , 0.33333333])
    """
    # Data in a appropriate format
    x_prior_ = np.array(x_prior)
    Vx_prior_ = np.array(Vx_prior)
    y_extra_ = np.array(y_extra)
    S_ = np.array(S)
    # Model calculations
    y_calc_ = S_.dot(x_prior_)
    G = sandwich(Vx_prior_, S_.T)
    if Vy_extra is not None:
        Vy_extra_ = np.array(Vy_extra)
        G = G + Vy_extra_
    G = np.sqrt(np.diag(G))
    delta = np.abs(y_extra_ - y_calc_)
    return delta / G


def chi_diag(x_prior, S, Vx_prior, Vy_extra, y_extra):
    r"""
    Function to calculate diagonal chi-value
    $\chi_{ind,i}$>>1 according to
    https://www.oecd-nea.org/jcms/pl_19760/intermediate-report-on-methods-and-approaches-to-provide-feedback-from-nuclear-and-covariance-data-adjustment-for-improvement-of-nuclear-data-files
    (page 9, equation (4.3))

    Parameters
    ----------
    x_prior: 1D iterable
        Vector in which we are going to apply GLS (M,)
    S : 2D iterable
        2D sensitivity of the model y=f(x) (M, N).
    Vx_prior : 2D iterable
        2D covariance matrix of x_prior (M, N).
    Vy_extra : 2D iterable
        2D covariance matrix for y_extra (M, N).
    y_extra : 1D iterable
        1D extra info on output (N,)

    Returns
    -------
    `numpy.ndarray`
        diagonal chi-value

    Results:
    -------
    chi diagonal >> 1 :
        Inconsistency may exist between |y_extra - y_calc| and covariance
        matrix, S*Vx_prior*S.T, and Vy_extra.

    Example
    -------
    >>> x_prior = [1, 2, 3]
    >>> y_extra = [2, 3, 4]
    >>> S = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
    >>> Vx_prior = np.array([[0, 0, 0], [0, 3, 0], [0, 0, 8]])
    >>> Vy_extra = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
    >>> N_e = 1
    >>> chi_diag(x_prior, S, Vx_prior, Vy_extra, y_extra)
    array([1., 2., 3.])
    """
    # Data in a appropriate format
    Vx_prior_ = np.array(Vx_prior)
    y_extra_ = np.array(y_extra)
    S_ = np.array(S)
    # Model calculations
    y_calc_ = S_.dot(x_prior)
    G_inv = _gls_G_inv(Vx_prior_, S_, Vy_extra=Vy_extra)
    # Chi diagonal calculations
    G_inv = np.sqrt(np.diag(G_inv))
    delta = np.abs(y_extra_ - y_calc_)
    return delta / G_inv


def chi_square(x_prior, S, Vx_prior, Vy_extra, y_extra, N_e):
    r"""
    Function to calculate contribution to chi-square value according to
    https://www.oecd-nea.org/jcms/pl_19760/intermediate-report-on-methods-and-approaches-to-provide-feedback-from-nuclear-and-covariance-data-adjustment-for-improvement-of-nuclear-data-files
    (page 10, equation (4.4))

    Parameters
    ----------
    x_prior: 1D iterable
        Vector in which we are going to apply GLS (M,)
    S : 2D iterable
        2D sensitivity of the model y=f(x) (M, N).
    Vx_prior : 2D iterable
        2D covariance matrix of x_prior (M, N).
    Vy_extra : 2D iterable
        2D covariance matrix for y_extra (M, N).
    y_extra : 1D iterable
        1D extra info on output (N,)
    N_e : `int`
        Number of experimental values used in adjustment.

    Returns
    -------
    `numpy.ndarray`
        contribution to chi-square value

    Results:
    -------
    chi square < 0 :
        The experiment is very effective in the adjustment.

    Example
    -------
    >>> x_prior = [1, 2, 3]
    >>> y_extra = [2, 3, 4]
    >>> S = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
    >>> Vx_prior = [[0, 0, 0], [0, 3, 0], [0, 0, 8]]
    >>> Vy_extra = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    >>> N_e = 1
    >>> chi_square(x_prior, S, Vx_prior, Vy_extra, y_extra, N_e)
    array([1.        , 0.25      , 0.11111111])
    """
    # Data in a appropriate format
    Vx_prior_ = np.array(Vx_prior)
    y_extra_ = np.array(y_extra)
    S_ = np.array(S)
    # Model calculations
    y_calc_ = S_.dot(x_prior)
    G_inv = _gls_G_inv(Vx_prior_, S_, Vy_extra)
    # Chi square calculations
    delta = y_extra_ - y_calc_
    chi_square = delta.T.dot(G_inv) * delta / N_e
    return chi_square


def ishikawa_factor(S, Vx_prior, Vy_extra):
    r"""
    Function to obtain Ishikawa factor according to
    https://www.oecd-nea.org/jcms/pl_19760/intermediate-report-on-methods-and-approaches-to-provide-feedback-from-nuclear-and-covariance-data-adjustment-for-improvement-of-nuclear-data-files
    (page 10, equation (4.5))

    Parameters
    ----------
    S : 2D iterable
        2D sensitivity of the model y=f(x) (M, N).
    Vx_prior : 2D iterable
        2D covariance matrix of x_prior (M, N).
    Vy_extra : 2D iterable
        2D covariance matrix for y_extra (M, N).

    Returns
    -------
    `numpy.ndarray`
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
    >>> S = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    >>> Vx_prior = [[0, 0, 0], [0, 3, 0], [0, 0, 8]]
    >>> Vy_extra = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    >>> ishikawa_factor(S, Vx_prior, Vy_extra)
    array([0., 3., 8.])
    """
    # Data in a appropriate format
    Vx_prior_ = np.array(Vx_prior)
    Vy_extra_ = np.array(Vy_extra)
    S_ = np.array(S)
    # Model calculations
    Vy_calc_ = sandwich(Vx_prior_, S_)
    # Ishikawa factor calculations
    Vy_values = np.diag(Vy_calc_)
    Vy_extra_ = np.diag(Vy_extra_)
    return Vy_values / Vy_extra_
