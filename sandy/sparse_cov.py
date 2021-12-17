# -*- coding: utf-8 -*-
"""
This module contains all classes and functions specific for processing
covariance matrix with scipy.sparse
"""
import pandas as pd
import numpy as np
import scipy
import scipy.sparse as sps
import scipy.sparse.linalg as spsl
import sandy

__author__ = "Aitor Bengoechea"


class sparse_cov:
    def __repr__(self):
        return self.data.__repr__()

    def __init__(self, datos):
        self.data = datos

    @property
    def data(self):
        """
        Covariance matrix as a sparse matrix

        Returns
        -------
        `array(<NxN sparse matrix of type '<class 'numpy.float64'>`
            Covariance matrix as a sparse matrix

        """
        return self._data

    @data.setter
    def data(self, data):
        data_ = np.array(data)
        self._data = sps.csr_matrix(data_)

    @property
    def size(self):
        return self.data.shape[0]

    def eig(self, sort=True):
        """
        Extract eigenvalues and eigenvectors.

        Parameters
        ----------
        sort : `bool`, optional, default is `True`
            flag to return sorted eigenvalues and eigenfunctions

        Returns
        -------
        `np.array`
            real part of eigenvalues sorted in descending order
        `np.array`
            matrix of eigenvectors

        Examples
        --------
        Extract eigenvalues of covariance matrix.
        >>> eigenvalues = sparse_cov([[1, 0.4],[0.4, 1]]).eig()[0]
        >>> comparison = eigenvalues == np.array([1.4, 0.6])
        >>> assert comparison.all() == True

        Extract eigenfunctions of covariance matrix.
        >>> eigenfunctions = sparse_cov([[1, 0.4],[0.4, 1]]).eig()[1].round(2)
        >>> comparison = eigenfunctions == np.array([[ 0.71, -0.71], [ 0.71,  0.71]])
        >>> assert comparison.all() == True
        """
        return get_eig(self.data, sort=sort)

    @property
    def std(self):
        """
        Extract standard deviations.

        Returns
        -------
        `np.array`
            1d array of standard deviations

        Examples
        --------
        >>> comparison = sparse_cov([[1, 0.4],[0.4, 1]]).std.data == np.array([1., 1.])
        >>> assert comparison.all() == True
        """
        cova = self.data
        std = np.sqrt(cova.diagonal())
        return sps.csr_matrix((std))

    @property
    def corr(self):
        """
        Extract correlation matrix.

        Returns
        -------
        `np.array`
            correlation matrix

        Examples
        --------
        >>> correlation = sparse_cov([[4, 2.4],[2.4, 9]]).corr.round(1)
        >>> comparison = correlation == np.array([[1.0, 0.4], [0.4, 1.0]])
        >>> assert comparison.all() == True
        """
        cov = self.data
        coeff = self.std
        corr = cov/coeff.multiply(coeff.T)
        return corr

    def invert(self):
        """
        Method for calculating the inverse matrix.

        Returns
        -------
        `sparse_cov`
            The inverse matrix.

        Examples
        --------
        >>> S = sparse_cov(np.diag(np.array([1, 2, 3])))
        >>> print(np.round(S.invert().data,2))
        (0, 0)	1.0
        (1, 1)	0.5
        (2, 2)	0.33
        """
        cov = sps.coo_matrix(self.data).tocsc()
        lu = spsl.splu(cov)
        eye = np.eye(cov.shape[0])
        cov_inv = lu.solve(eye)
        return self.__class__(cov_inv)

    def decompose(self):
        """
        Extract lower triangular matrix `L` for which `L*L^T == COV`.

        Returns
        -------
        `numpy.ndarray`
            lower triangular matrix

        Example
        -------
        >>> triangular = sparse_cov([[1, 0.4],[0.4, 1]]).decompose()
        >>> comparison = triangular == np.array([[1. , 0. ], [0.4, 1. ]])
        >>> assert comparison.all() == True
        """
        triangular = sps.tril(self.data, format='csc')
        return triangular.toarray()

    @classmethod
    def from_var(cls, var):
        """
        Construct the covariance matrix from the variance vector.

        Parameters
        ----------
        var : 1D iterable
            Variance vector.

        Returns
        -------
        `sparse_cov`
            Object containing the covariance matrix.

        Example
        -------
        >>> S = np.array([[1, 2], [3, 4]])
        >>> sensitivity = sparse_cov.from_var([1, 1])
        >>> print(sensitivity._gls_Vy_calc(S))
         (0, 1)	11.0
         (0, 0)	5.0
         (1, 1)	25.0
         (1, 0)	11.0
        """
        var_ = pd.Series(var).values
        cov = sps.diags(var_)
        return cls(cov.toarray())

    @classmethod
    def from_stdev(cls, std):
        """
        Construct the covariance matrix from the standard deviation vector.

        Parameters
        ----------
        std : `pandas.Series`
            Standard deviations vector.

        Returns
        -------
        `sparse_cov`
            Object containing the covariance matrix.

        Example
        -------
        >>> S = pd.Series(np.array([0, 2, 3]), index=pd.Index([1, 2, 3]))
        >>> print(sparse_cov.from_stdev(S).data)
         (1, 1)	4.0
         (2, 2)	9.0
        """
        std_ = pd.Series(std).values
        var = std_ * std_
        return sparse_cov.from_var(var)

    def _gls_Vy_calc(self, S):
        """
        2D calculated output using
        .. math::
            $$
            S\cdot V_{x_{prior}}\cdot S.T
            $$

        Parameters
        ----------
        S : 2D iterable
            Sensitivity matrix (MXN).

        Returns
        -------
        `pd.DataFrame`
            Covariance matrix `Vy_calc` calculated using
            S.dot(Vx_prior).dot(S.T)

        Example
        -------
        >>> S = np.array([[1, 2], [3, 4]])
        >>> sensitivity = sparse_cov.from_var([1, 1])
        >>> print(sensitivity._gls_Vy_calc(S))
         (0, 1)	11.0
         (0, 0)	5.0
         (1, 1)	25.0
         (1, 0)	11.0
        """
        S_ = pd.DataFrame(S).values
        S_ = sps.csr_matrix(S_)
        Vx_prior = self.data
        Vy_calc = S_.dot(Vx_prior).dot(S_.T)
        return Vy_calc

    def _gls_G(self, S, Vy_extra):
        """
        2D calculated output using
        .. math::
            $$
            S\cdot V_{x_{prior}}\cdot S.T + V_{y_{extra}}
            $$

        Parameters
        ----------
        S : 2D iterable
            Sensitivity matrix (MXN).
        Vy_extra : 2D iterable
            2D covariance matrix for y_extra (MXM).

        Returns
        -------
        `<NxN sparse matrix of type '<class 'numpy.float64'>'`
            Covariance matrix `G` calculated using
            S.dot(Vx_prior).dot(S.T) + Vy_extra

        Example
        -------
        >>> S = np.array([[1, 2], [3, 4]])
        >>> sensitivity = sparse_cov.from_var([1, 1])
        >>> Vy = np.diag(pd.Series([1, 1]))
        >>> print(sensitivity._gls_G(S, Vy))
         (0, 0)	6.0
         (0, 1)	11.0
         (1, 0)	11.0
         (1, 1)	26.0
        """
        Vy_extra_ = pd.DataFrame(Vy_extra).values
        Vy_extra_ = sparse_cov(Vy_extra_).data
        # GLS_sensitivity:
        Vy_calc = self._gls_Vy_calc(S)
        G = Vy_calc + Vy_extra_
        return G

    def _gls_G_inv(self, S, Vy_extra):
        """
        2D calculated output using
        .. math::
            $$
            \left(S\cdot V_{x_{prior}}\cdot S.T + V_{y_{extra}}\right)^{-1}
            $$

        Parameters
        ----------
        S : 2D iterable
            Sensitivity matrix (MXN).
        Vy_extra : 2D iterable
            2D covariance matrix for y_extra (MXM).

        Returns
        -------
        `sparse_cov`
            Covariance matrix `G_inv` calculated using
            (S.dot(Vx_prior).dot(S.T) + Vy_extra)^-1

        Example
        -------
        >>> S = np.array([[1, 2], [3, 4]])
        >>> sensitivity = sparse_cov.from_var([1, 1])
        >>> Vy = np.diag(pd.Series([1, 1]))
        >>> print(np.round(sensitivity._gls_G_inv(S, Vy).data, 2))
          (0, 0)	0.74
          (0, 1)	-0.31
          (1, 0)	-0.31
          (1, 1)	0.17
        """
        G = self._gls_G(S, Vy_extra).toarray()
        return sparse_cov(G).invert()

    def _gls_general_sensitivity(self, S, Vy_extra, threshold=None):
        """
        Method to obtain general sensitivity according to GLS
        .. math::
            $$
            V_{x_{prior}}\cdot S.T \cdot \left(S\cdot V_{x_{prior}}\cdot S.T + V_{y_{extra}}\right)^{-1}
            $$

        Parameters
        ----------
        S : 2D iterable
            Sensitivity matrix (MXN).
        Vy_extra : 2D iterable
            2D covariance matrix for y_extra (MXM).
        threshold : `int`, optional
            threshold to avoid numerical fluctuations. The default is None.

        Returns
        -------
        `GLS`
            GLS sensitivity for a given Vy_extra and S.

        Example
        -------
        >>> S = np.array([[1, 2], [3, 4]])
        >>> sensitivity = sparse_cov.from_var([1, 1])
        >>> Vy = np.diag(pd.Series([1, 1]))
        >>> print(np.round(sensitivity._gls_general_sensitivity(S, Vy), 2))
          (0, 0)	-0.2
          (0, 1)	0.2
          (1, 0)	0.23
          (1, 1)	0.06
        """
        S_ = sps.csr_matrix(pd.DataFrame(S).values)
        Vx_prior = self.data
        # GLS_sensitivity:
        G_inv = self._gls_G_inv(S, Vy_extra).data
        sensitivity = Vx_prior.dot(S_.T).dot(G_inv)
        if threshold is not None:
            sensitivity[sensitivity < threshold] = 0
        return sensitivity


def get_eig(cov, sort=True):
    try:
        E, V = sps.linalg.eig(cov)
    except:  # Because sometimes sps gives a error to send you to this method
        E, V = scipy.linalg.eig(cov.toarray())
    E = pd.Series(E.real)
    V = pd.DataFrame(V.real)
    if sort:
        idx = E.sort_values(ascending=False).index
        E = E.iloc[idx].reset_index(drop=True)
        V = V.iloc[idx].reset_index(drop=True)
    return E.values, V.values
