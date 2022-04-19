import functools

import numpy as np
import scipy
import scipy.linalg
import scipy
import scipy.sparse as sps
import scipy.sparse.linalg as spsl
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
import logging
import tables as tb
import os

import sandy
import pytest
pd.options.display.float_format = '{:.5e}'.format

__author__ = "Luca Fiorito"
__all__ = [
        "CategoryCov",
        "EnergyCov",
        "triu_matrix",
        "corr2cov",
        "random_corr",
        "random_cov",
        "sample_distribution",
        ]

S = np.array([[1, 1, 1],
              [1, 2, 1],
              [1, 3, 1]])
var = np.array([[0, 0, 0],
                [0, 2, 0],
                [0, 0, 3]])
minimal_covtest = pd.DataFrame(
    [[9437, 2, 1e-2, 9437, 2, 1e-2, 0.02],
     [9437, 2, 2e5, 9437, 2, 2e5, 0.09],
     [9437, 2, 1e-2, 9437, 102, 1e-2, 0.04],
     [9437, 2, 2e5, 9437, 102, 2e5, 0.05],
     [9437, 102, 1e-2, 9437, 102, 1e-2, 0.01],
     [9437, 102, 2e5, 9437, 102, 2e5, 0.01]],
    columns=["MAT", "MT", "E", "MAT1", "MT1", 'E1', "VAL"]
    )


def cov33csv(func):
    def inner(*args, **kwargs):
        key = "cov33csv"
        kw = kwargs.copy()
        if key in kw:
            if kw[key]:
                print(f"found argument '{key}', ignore oher arguments")
                out = func(
                    *args,
                    index_col=[0, 1, 2],
                    header=[0, 1, 2],
                    )
                out.index.names = ["MAT", "MT", "E"]
                out.columns.names = ["MAT", "MT", "E"]
                return out
            else:
                del kw[key]
        out = func(*args, **kw)
        return out
    return inner


class _Cov(np.ndarray):
    """Covariance matrix treated as a `numpy.ndarray`.

    Methods
    -------
    corr
        extract correlation matrix
    corr2cov
        produce covariance matrix given correlation matrix and standard
        deviation array
    eig
        get covariance matrix eigenvalues and eigenvectors
    get_L
        decompose and extract lower triangular matrix
    sampling
        draw random samples
    """

    def __new__(cls, arr):
        obj = np.ndarray.__new__(cls, arr.shape, float)
        obj[:] = arr[:]
        if not obj.ndim == 2:
            raise sandy.Error("covariance matrix must have two dimensions")
        if not np.allclose(obj, obj.T):
            raise sandy.Error("covariance matrix must be symmetric")
        if (np.diag(arr) < 0).any():
            raise sandy.Error("covariance matrix must have positive variances")
        return obj

    @staticmethod
    def _up2down(self):
        U = np.triu(self)
        L = np.triu(self, 1).T
        C = U + L
        return C

    def eig(self):
        """
        Extract eigenvalues and eigenvectors.

        Returns
        -------
        `Pandas.Series`
            real part of eigenvalues sorted in descending order
        `np.array`
            matrix of eigenvectors
        """
        E, V = scipy.linalg.eig(self)
        E, V = E.real, V.real
        return E, V

    def corr(self):
        """Extract correlation matrix.

        .. note:: zeros on the covariance matrix diagonal are translated
                  into zeros also on the the correlation matrix diagonal.

        Returns
        -------
        `sandy.formats.utils.Cov`
            correlation matrix
        """
        std = np.sqrt(np.diag(self))
        with np.errstate(divide='ignore', invalid='ignore'):
            coeff = np.true_divide(1, std)
            coeff[~ np.isfinite(coeff)] = 0  # -inf inf NaN
        corr = np.multiply(np.multiply(self.T, coeff).T, coeff)
        return self.__class__(corr)

    def _reduce_size(self):
        """
        Reduces the size of the matrix, erasing the null values.

        Returns
        -------
        nonzero_idxs : numpy.ndarray
            The indices of the diagonal that are not null.
        cov_reduced : sandy.core.cov._Cov
            The reduced matrix.

        """
        nonzero_idxs = np.flatnonzero(np.diag(self))
        cov_reduced = self[nonzero_idxs][:, nonzero_idxs]
        return nonzero_idxs, cov_reduced

    @classmethod
    def _restore_size(cls, nonzero_idxs, cov_reduced, dim):
        """
        Restore the size of the matrix

        Parameters
        ----------
        nonzero_idxs : numpy.ndarray
            The indices of the diagonal that are not null.
        cov_reduced : sandy.core.cov._Cov
            The reduced matrix.
        dim : int
            Dimension of the original matrix.

        Returns
        -------
        cov : sandy.core.cov._Cov
            Matrix of specified dimensions.

        """
        cov = _Cov(np.zeros((dim, dim)))
        for i, ni in enumerate(nonzero_idxs):
            cov[ni, nonzero_idxs] = cov_reduced[i]
        return cov

    def sampling(self, nsmp, seed=None):
        """
        Extract random samples from the covariance matrix, either using
        the cholesky or the eigenvalue decomposition.

        Parameters
        ----------
        nsmp : `int`
            number of samples
        seed : `int`
            seed for the random number generator (default is `None`)

        Returns
        -------
        `np.array`
            2D array of random samples with dimension `(self.shape[0], nsmp)`
        """
        dim = self.shape[0]
        np.random.seed(seed=seed)
        y = np.random.randn(dim, nsmp)
        nonzero_idxs, cov_reduced = self._reduce_size()
        L_reduced = cov_reduced.get_L()
        L = self.__class__._restore_size(nonzero_idxs, L_reduced, dim)
        samples = np.array(L.dot(y))
        return samples

    def get_L(self):
        """
        Extract lower triangular matrix `L` for which `L*L^T == self`.

        Returns
        -------
        `np.array`
            lower triangular matrix
        """
        try:
            L = scipy.linalg.cholesky(
                    self,
                    lower=True,
                    overwrite_a=False,
                    check_finite=False
                    )
        except np.linalg.linalg.LinAlgError:
            E, V = self.eig()
            E[E <= 0] = 0
            Esqrt = np.diag(np.sqrt(E))
            M = V.dot(Esqrt)
            Q, R = scipy.linalg.qr(M.T)
            L = R.T
        return L


class CategoryCov():
    """

    Properties
    ----------
    data
        covariance matrix as a dataframe
    size
        first dimension of the ocvariance matrix
    std
        `pd.Series` of standard deviations

    Methods
    -------
    get_corr
        extract correlation matrix from covariance matrix
    """

    def __repr__(self):
        return self.data.__repr__()

    def __init__(self, *args, **kwargs):
        self.data = pd.DataFrame(*args, **kwargs)

    @property
    def data(self):
        """
        Covariance matrix as a dataframe.

        Attributes
        ----------
        index : `pandas.Index` or `pandas.MultiIndex`
            indices
        columns : `pandas.Index` or `pandas.MultiIndex`
            columns
        values : `numpy.array`
            covariance values as `float`

        Returns
        -------
        `pandas.DataFrame`
            covariance matrix

        Notes
        -----
        ..note :: In the future, another tests will be implemented to check
        that the covariance matrix is symmetric and have positive variances.

        Examples
        --------
        >>> with pytest.raises(TypeError): sandy.CategoryCov(np.array[1])
        >>> with pytest.raises(TypeError): sandy.CategoryCov(np.array([[1, 2], [2, -4]]))
        >>> with pytest.raises(TypeError): sandy.CategoryCov(np.array([[1, 2], [3, 4]]))
        """
        return self._data

    @data.setter
    def data(self, data):
        self._data = pd.DataFrame(data, dtype=float)
        if not len(data.shape) == 2 and data.shape[0] == data.shape[1]:
            raise TypeError("Covariance matrix must have two dimensions")
        if not (np.diag(data) >= 0).all():
            raise TypeError("Covariance matrix must have positive variance")
        # Round to avoid numerical fluctuations
        if not (data.values.round(14) == data.values.T.round(14)).all():
            raise TypeError("Covariance matrix must be symmetric")

    @property
    def size(self):
        return self.data.values.shape[0]

    @property
    def std(self):
        """
        Extract standard deviations.

        Returns
        -------
        `pandas.Series`
            1d array of standard deviations

        Examples
        --------
        >>> sandy.CategoryCov([[1, 0.4],[0.4, 1]]).std
        0   1.00000e+00
        1   1.00000e+00
        Name: std, dtype: float64
        """
        cov = self.to_sparse().diagonal()
        std = np.sqrt(cov)
        return pd.Series(std, index=self.data.index, name="std")

    def eig(self, tolerance=None):
        """
        Extract eigenvalues and eigenvectors.

        Parameters
        ----------
        tolerance : `float`, optional, default is `None`
            replace all eigenvalues smaller than a given tolerance with zeros.
            The replacement condition is implemented as:

            .. math::
                $$
                \frac{e_i}{e_{MAX}} < tolerance
                $$

            Then, a `tolerance=1e-3` will replace all eigenvalues
            1000 times smaller than the largest eigenvalue.
            A `tolerance=0` will replace all negative eigenvalues.

        Returns
        -------
        `Pandas.Series`
            array of eigenvalues
        `pandas.DataFrame`
            matrix of eigenvectors

        Notes
        -----
        .. note:: only the real part of the eigenvalues is preserved
        
        .. note:: the discussion associated to the implementeation
                  of this algorithm is available [here](https://github.com/luca-fiorito-11/sandy/discussions/135)

        Examples
        --------
        Extract eigenvalues of correlation matrix.
        >>> sandy.CategoryCov([[1, 0.4], [0.4, 1]]).eig()[0]
        0   1.40000e+00
        1   6.00000e-01
        Name: eigenvalues, dtype: float64

        Extract eigenvectors of correlation matrix.
        >>> sandy.CategoryCov([[1, 0.4], [0.4, 1]]).eig()[1]
                    0            1
        0 7.07107e-01 -7.07107e-01
        1 7.07107e-01  7.07107e-01

        Extract eigenvalues of covariance matrix.
        >>> sandy.CategoryCov([[0.1, 0.1], [0.1, 1]]).eig()[0]
        0   8.90228e-02
        1   1.01098e+00
        Name: eigenvalues, dtype: float64

        Set up a tolerance.
        >>> sandy.CategoryCov([[0.1, 0.1], [0.1, 1]]).eig(tolerance=0.1)[0]
        0   0.00000e+00
        1   1.01098e+00
        Name: eigenvalues, dtype: float64

        Test with negative eigenvalues.
        >>> sandy.CategoryCov([[1, 2], [2, 1]]).eig()[0]
        0    3.00000e+00
        1   -1.00000e+00
        Name: eigenvalues, dtype: float64

        Replace negative eigenvalues.
        >>> sandy.CategoryCov([[1, 2], [2, 1]]).eig(tolerance=0)[0]
        0   3.00000e+00
        1   0.00000e+00
        Name: eigenvalues, dtype: float64

        Check output size.
        >>> cov = sandy.CategoryCov.random_cov(50, seed=11)
        >>> assert cov.eig()[0].size == cov.data.shape[0] == 50

        >>> sandy.CategoryCov([[1, 0.2, 0.1], [0.2, 2, 0], [0.1, 0, 3]]).eig()[0]
        0   9.56764e-01
        1   2.03815e+00
        2   3.00509e+00
        Name: eigenvalues, dtype: float64

        Real test on H1 file
        >>> endf6 = sandy.get_endf6_file("jeff_33", "xs", 10010)
        >>> ek = sandy.energy_grids.CASMO12
        >>> err = endf6.get_errorr(ek=ek, err=1)
        >>> cov = err.get_cov()
        >>> cov.eig()[0].sort_values(ascending=False).head(7)
        0    3.66411e-01
        1    7.05311e-03
        2    1.55346e-03
        3    1.60175e-04
        4    1.81374e-05
        5    1.81078e-06
        6    1.26691e-07
        Name: eigenvalues, dtype: float64

        >>> assert not (cov.eig()[0] >= 0).all()

        >>> assert (cov.eig(tolerance=0)[0] >= 0).all()
        """
        E, V = scipy.linalg.eig(self.data)
        E = pd.Series(E.real, name="eigenvalues")
        V = pd.DataFrame(V.real)
        if tolerance is not None:
            E[E/E.max() < tolerance] = 0
        return E, V

    def get_corr(self):
        """
        Extract correlation matrix.

        Returns
        -------
        `pandas.DataFrame`
            correlation matrix

        Examples
        --------
        >>> sandy.CategoryCov([[4, 2.4],[2.4, 9]]).get_corr()
                    0           1
        0 1.00000e+00 4.00000e-01
        1 4.00000e-01 1.00000e+00
        """
        cov = self.data.values
        with np.errstate(divide='ignore', invalid='ignore'):
            coeff = np.true_divide(1, self.std.values)
            coeff[~ np.isfinite(coeff)] = 0   # -inf inf NaN
        corr = np.multiply(np.multiply(cov, coeff).T, coeff)
        return pd.DataFrame(
            corr,
            index=self.data.index,
            columns=self.data.columns,
            )

    def invert(self, rows=None):
        """
        Method for calculating the inverse matrix.

        Parameters
        ----------
        tables : `bool`, optional
            Option to use row calculation for matrix calculations. The
            default is False.
        rows : `int`, optional
            Option to use row calculation for matrix calculations. This option
            defines the number of lines to be taken into account in each loop.
            The default is None.

        Returns
        -------
        `CategoryCov`
            The inverse matrix.

        Examples
        --------
        >>> S = sandy.CategoryCov(np.diag(np.array([1, 2, 3])))
        >>> S.invert()
                    0           1           2
        0 1.00000e+00 0.00000e+00 0.00000e+00
        1 0.00000e+00 5.00000e-01 0.00000e+00
        2 0.00000e+00 0.00000e+00 3.33333e-01

        >>> S = sandy.CategoryCov(np.diag(np.array([0, 2, 3])))
        >>> S.invert()
                    0           1           2
        0 0.00000e+00 0.00000e+00 0.00000e+00
        1 0.00000e+00 5.00000e-01 0.00000e+00
        2 0.00000e+00 0.00000e+00 3.33333e-01

        >>> S = sandy.CategoryCov(np.diag(np.array([0, 2, 3])))
        >>> S.invert(rows=1)
                    0           1           2
        0 0.00000e+00 0.00000e+00 0.00000e+00
        1 0.00000e+00 5.00000e-01 0.00000e+00
        2 0.00000e+00 0.00000e+00 3.33333e-01
        """
        index = self.data.index
        columns = self.data.columns
        M_nonzero_idxs, M_reduce = reduce_size(self.data)
        cov = sps.csc_matrix(M_reduce.values)
        rows_ = cov.shape[0] if rows is None else rows
        data = sparse_tables_inv(cov, rows=rows_)
        M_inv = restore_size(M_nonzero_idxs, data, len(self.data))
        M_inv = M_inv.reindex(index=index, columns=columns).fillna(0)
        return self.__class__(M_inv)

    def log2norm_cov(self, mu):
        """
        Transform covariance matrix to the one of the underlying normal
        distribution.

        Parameters
        ----------
        mu : iterable
            The desired mean values of the target lognormal distribution.

        Returns
        -------
        `CategoryCov` of the underlying normal covariance matrix

        Examples
        --------
        >>> cov = CategoryCov(pd.DataFrame([[8, 2, 3], [2, 10, 4], [3, 4, 7]], index=['A', 'B', 'C'], columns=['A', 'B', 'C']))
        >>> cov.log2norm_cov(pd.Series(np.ones(cov.data.shape[0]), index=cov.data.index))
                    A           B           C
        A 2.19722e+00 1.09861e+00 1.38629e+00
        B 1.09861e+00 2.39790e+00 1.60944e+00
        C 1.38629e+00 1.60944e+00 2.07944e+00

        >>> cov = CategoryCov(pd.DataFrame([[8, 2, 3], [2, 10, 4], [3, 4, 7]], index=['A', 'B', 'C'], columns=['A', 'B', 'C']))
        >>> mu = pd.Series([1, 2, .5], index=["A", "B", "C"])
        >>> cov.log2norm_cov(mu)
                    A           B           C
        A 2.19722e+00 6.93147e-01 1.94591e+00
        B 6.93147e-01 1.25276e+00 1.60944e+00
        C 1.94591e+00 1.60944e+00 3.36730e+00

        >>> cov = CategoryCov(pd.DataFrame([[8, 2, 3], [2, 10, 4], [3, 4, 7]], index=['A', 'B', 'C'], columns=['A', 'B', 'C']))
        >>> mu = [1, 2, .5]
        >>> cov.log2norm_cov(mu)
                    A           B           C
        A 2.19722e+00 6.93147e-01 1.94591e+00
        B 6.93147e-01 1.25276e+00 1.60944e+00
        C 1.94591e+00 1.60944e+00 3.36730e+00

        >>> cov = CategoryCov(pd.DataFrame([[8, 2, 3], [2, 10, 4], [3, 4, 7]], index=['A', 'B', 'C'], columns=['A', 'B', 'C']))
        >>> mu = np.array([1, 2, .5])
        >>> cov.log2norm_cov(mu)
                    A           B           C
        A 2.19722e+00 6.93147e-01 1.94591e+00
        B 6.93147e-01 1.25276e+00 1.60944e+00
        C 1.94591e+00 1.60944e+00 3.36730e+00

        Notes
        -----
        ..notes:: Reference for the equation is 10.1016/j.nima.2012.06.036
        .. math::
            $$
            cov(lnx_i, lnx_j) = \ln\left(\frac{cov(x_i,x_j)}{<x_i>\cdot<x_j>}+1\right)
            $$
        """
        mu_ = pd.Series(mu)
        mu_.index = self.data.index
        return self.__class__(np.log(self.sandwich(1 / mu_).data + 1))


    def log2norm_mean(self, mu):
        """
        Transform mean values to the mean values of the undelying normal
        distribution.

        Parameters
        ----------
        mu : iterable
            The target mean values.

        Returns
        -------
        `pd.Series` of the underlyig normal distribution mean values

        Examples
        --------
        >>> cov = CategoryCov(pd.DataFrame([[8, 2, 3], [2, 10, 4], [3, 4, 7]], index=['A', 'B', 'C'], columns=['A', 'B', 'C']))
        >>> mu = pd.Series(np.ones(cov.data.shape[0]), index=cov.data.index)
        >>> cov.log2norm_mean(mu)
        A   -1.09861e+00
        B   -1.19895e+00
        C   -1.03972e+00
        dtype: float64

        >>> cov = CategoryCov(pd.DataFrame([[8, 2, 3], [2, 10, 4], [3, 4, 7]], index=['A', 'B', 'C'], columns=['A', 'B', 'C']))
        >>> cov.log2norm_mean([1, 1, 1])
        A   -1.09861e+00
        B   -1.19895e+00
        C   -1.03972e+00
        dtype: float64

        >>> cov = CategoryCov(pd.DataFrame([[8, 2, 3], [2, 10, 4], [3, 4, 7]], index=['A', 'B', 'C'], columns=['A', 'B', 'C']))
        >>> mu = np.ones(cov.data.shape[0])
        >>> cov.log2norm_mean(mu)
        A   -1.09861e+00
        B   -1.19895e+00
        C   -1.03972e+00
        dtype: float64

        Reindexing example
        
        """
        mu_ = pd.Series(mu)
        mu_.index = self.data.index
        return np.log(mu_**2 / np.sqrt(np.diag(self.data) + mu_**2))

    def sampling(self, nsmp, seed=None, rows=None, pdf='normal',
                 tolerance=None, relative=True):
        """
        Extract perturbation coefficients according to chosen distribution with
        covariance from given covariance matrix. See note for non-normal
        distribution sampling.
        The samples' mean will be 1 or 0 depending on `relative` kwarg.

        Parameters
        ----------
        nsmp : `int`
            number of samples.
        seed : `int`, optional, default is `None`
            seed for the random number generator (by default use `numpy`
            dafault pseudo-random number generator).
        rows : `int`, optional, default is `None`
            option to use row calculation for matrix calculations. This option
            defines the number of lines to be taken into account in each loop.
        pdf : `str`, optional, default is 'normal'
            random numbers distribution.
            Available distributions are:
                * `'normal'`
                * `'uniform'`
                * `'lognormal'`
        tolerance : `float`, optional, default is `None`
            replace all eigenvalues smaller than a given tolerance with zeros.
        relative : `bool`, optional, default is `True`
            flag to switch between relative and absolute covariance matrix
            handling
                * `True`: samples' mean will be 1
                * `False`: samples' mean will be 0

        Returns
        -------
        `sandy.Samples`
            object containing samples

        Notes
        -----
        .. note:: sampling with uniform distribution is performed on
            diagonal covariance matrix, neglecting all correlations.

        .. note:: sampling with lognormal distribution gives a set of samples
            with mean=1 as lognormal distribution can not have mean=0.
            Therefore, `relative` parameter does not apply to it.

        Examples
        --------
        Draw 3 sets of samples using custom seed:
        >>> sandy.CategoryCov([[1, 0.4],[0.4, 1]]).sampling(3, seed=11)
                     0            1
        0 -7.49455e-01 -2.13159e+00
        1  1.28607e+00  1.10684e+00
        2  1.48457e+00  9.00879e-01

        >>> sandy.CategoryCov([[1, 0.4],[0.4, 1]]).sampling(3, seed=11, rows=1)
                     0            1
        0 -7.49455e-01 -2.13159e+00
        1  1.28607e+00  1.10684e+00
        2  1.48457e+00  9.00879e-01

        >>> sample = sandy.CategoryCov([[1, 0.4],[0.4, 1]]).sampling(1000000, seed=11)
        >>> sample.data.cov()
                    0           1
        0 9.98662e-01 3.99417e-01
        1 3.99417e-01 9.98156e-01

        Small negative eigenvalue:
        >>> sandy.CategoryCov([[1, -.2],[-.2, 3]]).sampling(3, seed=11, tolerance=0)
                    0           1
        0 2.74945e+00 5.21505e+00
        1 7.13927e-01 1.07147e+00
        2 5.15435e-01 1.64683e+00

        >>> sandy.CategoryCov([[1, -.2],[-.2, 3]]).sampling(1000000, seed=11, tolerance=0).data.cov()
                     0            1
        0  9.98662e-01 -1.99822e-01
        1 -1.99822e-01  2.99437e+00

        Sampling with different `pdf`:
        >>> sandy.CategoryCov([[1, -.2],[-.2, 3]]).sampling(3, seed=11, pdf='uniform', tolerance=0)
                     0           1
        0 -1.07578e-01 2.34960e+00
        1 -6.64587e-01 5.21222e-01
        2  8.72585e-01 9.12563e-01

        >>> sandy.CategoryCov([[1, .2],[.2, 3]]).sampling(3, seed=11, pdf='lognormal', tolerance=0)
                    0           1
        0 3.03419e+00 1.57919e+01
        1 5.57248e-01 4.74160e-01
        2 4.72366e-01 6.50840e-01

        >>> sandy.CategoryCov([[1, -.2],[-.2, 3]]).sampling(1000000, seed=11, pdf='uniform', tolerance=0).data.cov()
                     0            1
        0  1.00042e+00 -1.58806e-03
        1 -1.58806e-03  3.00327e+00

        >>> sandy.CategoryCov([[1, .2],[.2, 3]]).sampling(1000000, seed=11, pdf='lognormal', tolerance=0).data.cov()
                    0           1
        0 1.00219e+00 1.99199e-01
        1 1.99199e-01 3.02605e+00

        `relative` kwarg usage:
        >>> sandy.CategoryCov([[1, -.2],[-.2, 3]]).sampling(1000000, seed=11, pdf='normal', tolerance=0, relative=True).data.mean(axis=0)
        0   1.00014e+00
        1   9.99350e-01
        dtype: float64

        >>> sandy.CategoryCov([[1, -.2],[-.2, 3]]).sampling(1000000, seed=11, pdf='normal', tolerance=0, relative=False).data.mean(axis=0)
        0    1.41735e-04
        1   -6.49679e-04
        dtype: float64

        >>> sandy.CategoryCov([[1, -.2],[-.2, 3]]).sampling(1000000, seed=11, pdf='uniform', tolerance=0, relative=True).data.mean(axis=0)
        0   9.98106e-01
        1   9.99284e-01
        dtype: float64

        >>> sandy.CategoryCov([[1, -.2],[-.2, 3]]).sampling(1000000, seed=11, pdf='uniform', tolerance=0, relative=False).data.mean(axis=0)
        0   -1.89367e-03
        1   -7.15929e-04
        dtype: float64

        Lognormal distribution sampling indeoendency from `relative` kwarg
        >>> sandy.CategoryCov([[1, .2],[.2, 3]]).sampling(1000000, seed=11, pdf='lognormal', tolerance=0, relative=True).data.mean(axis=0)
        0   9.99902e-01
        1   9.99284e-01
        dtype: float64

        >>> sandy.CategoryCov([[1, .2],[.2, 3]]).sampling(1000000, seed=11, pdf='lognormal', tolerance=0, relative=False).data.mean(axis=0)
        0   9.99902e-01
        1   9.99284e-01
        dtype: float64
        """
        dim = self.data.shape[0]
        pdf_ = pdf if pdf != 'lognormal' else 'normal'
        y = sample_distribution(dim, nsmp, seed=seed, pdf=pdf_) - 1
        y = sps.csc_matrix(y)
        # the covariance matrix to decompose is created depending on the chosen
        # pdf
        if pdf == 'uniform':
            to_decompose = self.__class__(np.diag(np.diag(self.data)))
        elif pdf == 'lognormal':
            ones = np.ones(self.data.shape[0])
            to_decompose = self.log2norm_cov(ones)
        else:
            to_decompose = self
        L = sps.csr_matrix(to_decompose.get_L(rows=rows,
                                              tolerance=tolerance))
        samples = pd.DataFrame(L.dot(y).toarray(), index=self.data.index,
                               columns=list(range(nsmp)))
        if pdf == 'lognormal':
            # mean value of lognormally sampled distributions will be one by
            # defaul
            samples = np.exp(samples.add(self.log2norm_mean(ones), axis=0))
        elif relative:
            samples += 1
        return sandy.Samples(samples.T)

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
        `CategoryCov`
            Object containing the covariance matrix.

        Example
        -------
        >>> S = pd.Series(np.array([0, 2, 3]), index=pd.Index([1, 2, 3]))
        >>> cov = sandy.CategoryCov.from_var(S)
        >>> cov
                    1           2           3
        1 0.00000e+00 0.00000e+00 0.00000e+00
        2 0.00000e+00 2.00000e+00 0.00000e+00
        3 0.00000e+00 0.00000e+00 3.00000e+00

        >>> assert type(cov) is sandy.CategoryCov

        >>> S = sandy.CategoryCov.from_var((1, 2, 3))
        >>> S
                    0           1           2
        0 1.00000e+00 0.00000e+00 0.00000e+00
        1 0.00000e+00 2.00000e+00 0.00000e+00
        2 0.00000e+00 0.00000e+00 3.00000e+00

        >>> assert type(S) is sandy.CategoryCov
        >>> assert type(sandy.CategoryCov.from_var([1, 2, 3])) is sandy.CategoryCov
        """
        var_ = pd.Series(var)
        cov_values = sps.diags(var_.values).toarray()
        cov = pd.DataFrame(cov_values,
                           index=var_.index, columns=var_.index)
        return cls(cov)

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
        `CategoryCov`
            Object containing the covariance matrix.

        Example
        -------
        >>> S = pd.Series(np.array([0, 2, 3]), index=pd.Index([1, 2, 3]))
        >>> cov = sandy.CategoryCov.from_stdev(S)
        >>> cov
                    1           2           3
        1 0.00000e+00 0.00000e+00 0.00000e+00
        2 0.00000e+00 4.00000e+00 0.00000e+00
        3 0.00000e+00 0.00000e+00 9.00000e+00

        >>> assert type(cov) is sandy.CategoryCov

        >>> S = sandy.CategoryCov.from_stdev((1, 2, 3))
        >>> S
                    0           1           2
        0 1.00000e+00 0.00000e+00 0.00000e+00
        1 0.00000e+00 4.00000e+00 0.00000e+00
        2 0.00000e+00 0.00000e+00 9.00000e+00

        >>> assert type(S) is sandy.CategoryCov
        >>> assert type(sandy.CategoryCov.from_stdev([1, 2, 3])) is sandy.CategoryCov
        """
        std_ = pd.Series(std)
        var = std_ * std_
        return cls.from_var(var)

    @classmethod
    def from_stack(cls, data_stack, index, columns, values, rows=10000000,
                   kind='upper'):
        """
        Create a covariance matrix from a stacked dataframe.

        Parameters
        ----------
        data_stack : `pd.Dataframe`
            Stacked dataframe.
        index : 1D iterable, optional
            Index of the final covariance matrix.
        columns : 1D iterable, optional
            Columns of the final covariance matrix.
        values : `str`, optional
            Name of the column where the values are located.
        rows : `int`, optional
            Number of rows to take into account into each loop. The default
            is 10000000.
        kind : `str`, optional
            Select if the stack data represents upper or lower triangular
            matrix. The default is 'upper.

        Returns
        -------
        `sandy.CategoryCov`
            Covarinace matrix.

        Examples
        --------
        If the stack data represents the covariance matrix:
        >>> S = pd.DataFrame(np.array([[1, 1, 1], [1, 2, 1], [1, 1, 1]]))
        >>> S = S.stack().reset_index().rename(columns = {'level_0': 'dim1', 'level_1': 'dim2', 0: 'cov'})
        >>> S = S[S['cov'] != 0]
        >>> sandy.CategoryCov.from_stack(S, index=['dim1'], columns=['dim2'], values='cov', kind='all')
        dim2           0           1           2
        dim1                                    
        0    1.00000e+00 1.00000e+00 1.00000e+00
        1    1.00000e+00 2.00000e+00 1.00000e+00
        2    1.00000e+00 1.00000e+00 1.00000e+00

        If the stack data represents only the upper triangular part of the
        covariance matrix:
        >>> test_1 = sandy.CategoryCov.from_stack(minimal_covtest, index=["MAT", "MT", "E"], columns=["MAT1", "MT1", "E1"], values='VAL').data
        >>> test_1
        	        MAT1	    9437
                    MT1	        2	                    102
                    E1	        1.00000e-02	2.00000e+05	1.00000e-02	2.00000e+05
        MAT	   MT	E				
        9437	2	1.00000e-02	2.00000e-02	0.00000e+00	4.00000e-02	0.00000e+00
                    2.00000e+05	0.00000e+00	9.00000e-02	0.00000e+00	5.00000e-02
              102	1.00000e-02	4.00000e-02	0.00000e+00	1.00000e-02	0.00000e+00
                    2.00000e+05	0.00000e+00	5.00000e-02	0.00000e+00	1.00000e-02

        >>> test_2 = sandy.CategoryCov.from_stack(minimal_covtest, index=["MAT", "MT", "E"], columns=["MAT1", "MT1", "E1"], values='VAL', rows=1).data
        >>> test_2
        	        MAT1	    9437
                    MT1	        2	                    102
                    E1	        1.00000e-02	2.00000e+05	1.00000e-02	2.00000e+05
        MAT	   MT	E				
        9437	2	1.00000e-02	2.00000e-02	0.00000e+00	4.00000e-02	0.00000e+00
                    2.00000e+05	0.00000e+00	9.00000e-02	0.00000e+00	5.00000e-02
              102	1.00000e-02	4.00000e-02	0.00000e+00	1.00000e-02	0.00000e+00
                    2.00000e+05	0.00000e+00	5.00000e-02	0.00000e+00	1.00000e-02

        >>> assert (test_1 == test_2).all().all()

        If the stack data represents only the lower triangular part of the
        covariance matrix:
        >>> test_1 = sandy.CategoryCov.from_stack(minimal_covtest, index=["MAT1", "MT1", "E1"], columns=["MAT", "MT", "E"], values='VAL', kind="lower").data
        >>> test_1
        	        MAT	        9437
                    MT	        2	                    102
                    E	        1.00000e-02	2.00000e+05	1.00000e-02	2.00000e+05
        MAT1  MT1	E1				
        9437	2	1.00000e-02	2.00000e-02	0.00000e+00	4.00000e-02	0.00000e+00
                    2.00000e+05	0.00000e+00	9.00000e-02	0.00000e+00	5.00000e-02
              102	1.00000e-02	4.00000e-02	0.00000e+00	1.00000e-02	0.00000e+00
                    2.00000e+05	0.00000e+00	5.00000e-02	0.00000e+00	1.00000e-02

        >>> test_2 = sandy.CategoryCov.from_stack(minimal_covtest, index=["MAT1", "MT1", "E1"], columns=["MAT", "MT", "E"], values='VAL', kind="lower", rows=1).data
        >>> test_2
        	        MAT 	    9437
                    MT	        2	                    102
                    E	        1.00000e-02	2.00000e+05	1.00000e-02	2.00000e+05
        MAT1  MT1	E1				
        9437	2	1.00000e-02	2.00000e-02	0.00000e+00	4.00000e-02	0.00000e+00
                    2.00000e+05	0.00000e+00	9.00000e-02	0.00000e+00	5.00000e-02
              102	1.00000e-02	4.00000e-02	0.00000e+00	1.00000e-02	0.00000e+00
                    2.00000e+05	0.00000e+00	5.00000e-02	0.00000e+00	1.00000e-02

        >>> assert (test_1 == test_2).all().all()
        """
        cov = segmented_pivot_table(data_stack, rows=rows, index=index,
                                    columns=columns, values=values)
        if kind == 'all':
            return cls(cov)
        else:
            return triu_matrix(cov, kind=kind)

    def _gls_Vy_calc(self, S, rows=None):
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
        rows : `int`, optional
            Option to use row calculation for matrix calculations. This option
            defines the number of lines to be taken into account in each loop.
            The default is None.

        Returns
        -------
        `pd.DataFrame`
            Covariance matrix `Vy_calc` calculated using
            S.dot(Vx_prior).dot(S.T)

        Example
        -------
        >>> S = np.array([[1, 2], [3, 4]])
        >>> cov = sandy.CategoryCov.from_var([1, 1])
        >>> cov._gls_Vy_calc(S)
                      0	          1
        0	5.00000e+00	1.10000e+01
        1	1.10000e+01	2.50000e+01

        >>> cov._gls_Vy_calc(S, rows=1)
                      0	          1
        0	5.00000e+00	1.10000e+01
        1	1.10000e+01	2.50000e+01
        """
        index = pd.DataFrame(S).index
        S_ = pd.DataFrame(S).values
        rows_ = S_.shape[0] if rows is None else rows
        Vy_calc = sparse_tables_dot_multiple([S_, self.data.values,
                                              S_.T], rows=rows_)
        return pd.DataFrame(Vy_calc, index=index, columns=index)

    def _gls_G(self, S, Vy_extra=None, rows=None):
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
        Vy_extra : 2D iterable, optional.
            2D covariance matrix for y_extra (MXM).
        rows : `int`, optional
            Option to use row calculation for matrix calculations. This option
            defines the number of lines to be taken into account in each loop.
            The default is None.

        Returns
        -------
        `pd.DataFrame`
            Covariance matrix `G` calculated using
            S.dot(Vx_prior).dot(S.T) + Vy_extra

        Example
        -------
        >>> S = np.array([[1, 2], [3, 4]])
        >>> cov = sandy.CategoryCov.from_var([1, 1])
        >>> Vy = np.diag(pd.Series([1, 1]))
        >>> cov._gls_G(S, Vy)
                    	0	      1
        0	6.00000e+00	1.10000e+01
        1	1.10000e+01	2.60000e+01

        >>> cov._gls_G(S, Vy, rows=1)
                    	0	      1
        0	6.00000e+00	1.10000e+01
        1	1.10000e+01	2.60000e+01

        >>> cov._gls_G(S)
                      0	          1
        0	5.00000e+00	1.10000e+01
        1	1.10000e+01	2.50000e+01

        >>> cov._gls_G(S, rows=1)
                      0	          1
        0	5.00000e+00	1.10000e+01
        1	1.10000e+01	2.50000e+01
        """
        # GLS_sensitivity:
        Vy_calc = self._gls_Vy_calc(S, rows=rows)
        if Vy_extra is not None:
            # Data in a appropriate format
            Vy_extra_ = sandy.CategoryCov(Vy_extra).data
            index = pd.DataFrame(Vy_extra).index
            Vy_extra_ = Vy_extra_.values
            Vy_calc = Vy_calc.reindex(index=index, columns=index).fillna(0).values
            # Calculations:
            Vy_calc = sps.csr_matrix(Vy_calc)
            Vy_extra_ = sps.csr_matrix(Vy_extra_)
            # G calculation
            G = Vy_calc + Vy_extra_
            G = pd.DataFrame(G.toarray(), index=index, columns=index)
        else:
            G = Vy_calc
        return G

    def _gls_G_inv(self, S, Vy_extra=None, rows=None):
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
        Vy_extra : 2D iterable, optional
            2D covariance matrix for y_extra (MXM).
        rows : `int`, optional
            Option to use row calculation for matrix calculations. This option
            defines the number of lines to be taken into account in each loop.
            The default is None.

        Returns
        -------
        `pd.DataFrame`
            Covariance matrix `G_inv` calculated using
            (S.dot(Vx_prior).dot(S.T) + Vy_extra)^-1

        Example
        -------
        >>> S = np.array([[1, 2], [3, 4]])
        >>> cov = sandy.CategoryCov.from_var([1, 1])
        >>> Vy = np.diag(pd.Series([1, 1]))
        >>> cov._gls_G_inv(S, Vy)
                       0	          1
        0	 7.42857e-01   -3.14286e-01
        1	-3.14286e-01	1.71429e-01

        >>> cov._gls_G_inv(S, Vy, rows=1)
                       0	          1
        0	 7.42857e-01   -3.14286e-01
        1	-3.14286e-01	1.71429e-01

        >>> cov._gls_G_inv(S)
                      0	               1
        0	 6.25000e+00	-2.75000e+00
        1	-2.75000e+00	 1.25000e+00

        >>> cov._gls_G_inv(S, rows=1)
                      0	               1
        0	 6.25000e+00	-2.75000e+00
        1	-2.75000e+00	 1.25000e+00
        """
        if Vy_extra is not None:
            index = pd.DataFrame(Vy_extra).index
            G = self._gls_G(S, Vy_extra=Vy_extra, rows=rows).values
        else:
            index = pd.DataFrame(S).index
            G = self._gls_Vy_calc(S, rows=rows).values
        G_inv = sandy.CategoryCov(G).invert(rows=rows).data.values
        return pd.DataFrame(G_inv, index=index, columns=index)

    def _gls_general_sensitivity(self, S, Vy_extra=None,
                                 rows=None, threshold=None):
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
        rows : `int`, optional
            Option to use row calculation for matrix calculations. This option
            defines the number of lines to be taken into account in each loop.
            The default is None.
        threshold : `int`, optional
            threshold to avoid numerical fluctuations. The default is None.

        Returns
        -------
        `GLS`
            GLS sensitivity for a given Vy_extra and S.

        Example
        -------
        >>> S = np.array([[1, 2], [3, 4]])
        >>> cov = sandy.CategoryCov.from_var([1, 1])
        >>> Vy = np.diag(pd.Series([1, 1]))
        >>> cov._gls_general_sensitivity(S, Vy)
                      0	              1
        0	-2.00000e-01	2.00000e-01
        1	2.28571e-01	    5.71429e-02

        >>> S = pd.DataFrame([[1, 2], [3, 4]], index=[1, 2],columns=[3, 4])
        >>> cov = sandy.CategoryCov.from_var([1, 1])
        >>> Vy = pd.DataFrame([[1, 0], [0, 1]], index=[1, 2], columns=[1, 2])
        >>> cov._gls_general_sensitivity(S, Vy_extra=Vy)
                      1	              2
        3	-2.00000e-01	2.00000e-01
        4	 2.28571e-01	5.71429e-02

        >>> cov._gls_general_sensitivity(S, Vy_extra=Vy, rows=1)
                      1	              2
        3	-2.00000e-01	2.00000e-01
        4	 2.28571e-01	5.71429e-02

        >>> cov._gls_general_sensitivity(S)
            	       1	           2
        3	-2.00000e+00	 1.00000e+00
        4	 1.50000e+00	-5.00000e-01

        >>> cov._gls_general_sensitivity(S, rows=1)
            	       1	           2
        3	-2.00000e+00	 1.00000e+00
        4	 1.50000e+00	-5.00000e-01
        """
        index = pd.DataFrame(S).columns
        columns = pd.DataFrame(S).index
        S_ = pd.DataFrame(S).values
        # GLS_sensitivity:
        G_inv = self._gls_G_inv(S, Vy_extra=Vy_extra, rows=rows).values
        rows_ = S_.shape[0] if rows is None else rows
        sensitivity = sparse_tables_dot_multiple([self.data.values, S_.T,
                                                  G_inv], rows=rows_)
        if threshold is not None:
            sensitivity[abs(sensitivity) < threshold] = 0
        return pd.DataFrame(sensitivity, index=index, columns=columns)

    def _gls_constrained_sensitivity(self, S, rows=None,
                                     threshold=None):
        """
        Method to obtain sensitivity according to constrained Least-Squares:
        .. math::
            $$
            \left(S\cdot V_{x_{prior}}\cdot S.T + V_{y_{extra}}\right)^{-1} \cdot S \cdot V_{x_{prior}}
            $$

        Parameters
        ----------
        S : 2D iterable
            Sensitivity matrix (MXN).
        rows : `int`, optional
            Option to use row calculation for matrix calculations. This option
            defines the number of lines to be taken into account in each loop.
            The default is None.
        threshold : `int`, optional
            threshold to avoid numerical fluctuations. The default is None.

        Returns
        -------
        `pd.DataFrame`
            constrained Least-Squares sensitivity.

        Notes
        -----
        ..note :: This method is equivalent to `_gls_general_sensitivity`
        but for a constrained system

        Example
        -------
        >>> S = np.array([[1, 2], [3, 4]])
        >>> cov = CategoryCov.from_var([1, 1])
        >>> cov._gls_constrained_sensitivity(S)
                      0	              1
        0	-2.00000e+00	1.50000e+00
        1	 1.00000e+00   -5.00000e-01

        >>> cov._gls_constrained_sensitivity(S, rows=1)
                      0	              1
        0	-2.00000e+00	1.50000e+00
        1	 1.00000e+00   -5.00000e-01
        """
        # Data in a appropiate format
        S_ = pd.DataFrame(S)
        index = S_.index
        columns = S_.columns
        G_inv = self._gls_G_inv(S, rows=rows).values
        rows_ = S_.shape[0] if rows is None else rows
        sensitivity = sparse_tables_dot_multiple([G_inv, S_,
                                                 self.data.values],
                                                 rows=rows_)
        if threshold is not None:
            sensitivity[abs(sensitivity) < threshold] = 0
        return pd.DataFrame(sensitivity, index=index, columns=columns)

    def _gls_cov_sensitivity(self, S, Vy_extra=None,
                             rows=None, threshold=None):
        """
        Method to obtain covariance sensitivity according to GLS:
        .. math::
            $$
            V_{x_{prior}}\cdot S^T \cdot \left(S\cdot V_{x_{prior}}\cdot S.T + V_{y_{extra}}\right)^{-1} \cdot S
            $$

        Parameters
        ----------
        S : 2D iterable
            Sensitivity matrix (MXN).
        Vy_extra : 2D iterable
            2D covariance matrix for y_extra (MXM).
        rows : `int`, optional
            Option to use row calculation for matrix calculations. This option
            defines the number of lines to be taken into account in each loop.
            The default is None.
        threshold : `int`, optional
            threshold to avoid numerical fluctuations. The default is None.

        Returns
        -------
        `pd.DataFrame`
            GlS sensitivity for a given Vy and S.

        Example
        -------
        >>> S = np.array([[1, 2], [3, 4]])
        >>> cov = sandy.CategoryCov.from_var([1, 1])
        >>> Vy = np.diag(pd.Series([1, 1]))
        >>> cov._gls_cov_sensitivity(S, Vy)
                      0	          1
        0	4.00000e-01	4.00000e-01
        1	4.00000e-01	6.85714e-01

        >>> S = pd.DataFrame([[1, 2], [3, 4]], index=[1, 2],columns=[3, 4])
        >>> cov = sandy.CategoryCov.from_var([1, 1])
        >>> Vy = pd.DataFrame([[1, 0], [0, 1]], index=[1, 2], columns=[1, 2])
        >>> cov._gls_cov_sensitivity(S, Vy)
                      3	          4
        3	4.00000e-01	4.00000e-01
        4	4.00000e-01	6.85714e-01

        >>> cov._gls_cov_sensitivity(S, Vy, rows=1)
                      3	          4
        3	4.00000e-01	4.00000e-01
        4	4.00000e-01	6.85714e-01
        """
        index = columns = pd.DataFrame(S).columns
        S_ = pd.DataFrame(S).values
        general_sens = self._gls_general_sensitivity(S, Vy_extra=Vy_extra,
                                                     rows=rows,
                                                     threshold=threshold).values
        rows_ = S_.shape[0] if rows is None else rows
        cov_sens = sparse_tables_dot(general_sens, S_, rows=rows_).toarray()
        if threshold is not None:
            cov_sens[abs(cov_sens) < threshold] = 0
        return pd.DataFrame(cov_sens, index=index, columns=columns)

    def gls_update(self, S, Vy_extra=None, rows=None,
                   threshold=None):
        """
        Perform GlS update for a given variance and sensitivity:
        .. math::
            $$
            V_{x_{post}} = V_{x_{prior}} - V_{x_{prior}}\cdot S.T \cdot \left(S\cdot V_{x_{prior}}\cdot S.T + V_{y_{extra}}\right)^{-1} \cdot S \cdot V_{x_{prior}}
            $$

        Parameters
        ----------
        Vy_extra : 2D iterable, optional
            2D covariance matrix for y_extra (MXM).
        S : 2D iterable
            Sensitivity matrix (MXN).
        rows : `int`, optional
            Option to use row calculation for matrix calculations. This option
            defines the number of lines to be taken into account in each loop.
            The default is None.
        threshold : `int`, optional
            Thereshold to avoid numerical fluctuations. The default is None.

        Returns
        -------
        `CategoryCov`
            GLS method apply to a CategoryCov object for a given Vy and S.

        Example
        -------
        >>> S = np.array([[1, 2], [3, 4]])
        >>> cov = sandy.CategoryCov.from_var([1, 1])
        >>> Vy = np.diag(pd.Series([1, 1]))
        >>> cov.gls_update(S, Vy)
                     0            1
        0  6.00000e-01 -4.00000e-01
        1 -4.00000e-01  3.14286e-01

        >>> cov.gls_update(S, Vy, rows=1)
                     0            1
        0  6.00000e-01 -4.00000e-01
        1 -4.00000e-01  3.14286e-01
        """
        index, columns = self.data.index, self.data.columns
        A = self._gls_cov_sensitivity(S, Vy_extra=Vy_extra,
                                      rows=rows, threshold=threshold).values
        rows_ = self.data.shape[0] if rows is None else rows
        Vx_prior = self.to_sparse(method='csc_matrix')
        diff = sparse_tables_dot(A, Vx_prior, rows=rows_)
        # gls update
        Vx_post = Vx_prior - diff
        Vx_post = Vx_post.toarray()
        if threshold is not None:
            Vx_post[abs(Vx_post) < threshold] = 0
        return self.__class__(pd.DataFrame(Vx_post, index=index, columns=columns))

    def constrained_gls_update(self, S, rows=None,
                               threshold=None):
        """
        Perform constrained Least-Squares update for a given sensitivity:
        .. math::
            $$
            V_{x_{post}} = V_{x_{prior}} - V_{x_{prior}}\cdot S.T \cdot \left(S\cdot V_{x_{prior}}\cdot S.T\right)^{-1} \cdot S \cdot V_{x_{prior}}
            $$

        Parameters
        ----------
        S : 2D iterable
            Sensitivity matrix (MXN).
        rows : `int`, optional
            Option to use row calculation for matrix calculations. This option
            defines the number of lines to be taken into account in each loop.
            The default is None.
        threshold : `int`, optional
            Thereshold to avoid numerical fluctuations. The default is None.

        Returns
        -------
        `CategoryCov`
            Constrained Least-squares method apply to a CategoryCov object
            for a given S.

        Notes
        -----
        ..note :: This method is equivalent to `gls_update` but for a
        constrained system

        Example
        -------
        >>> S = np.array([[1, 2], [3, 4]])
        >>> cov = sandy.CategoryCov.from_var([1, 1])
        >>> cov_update = cov.constrained_gls_update(S).data.round(decimals=6)
        >>> assert np.amax(cov_update.values) == 0.0

        >>> cov_update = cov.constrained_gls_update(S, rows=1).data.round(decimals=6)
        >>> assert np.amax(cov_update.values) == 0.0
        """
        return self.gls_update(S, Vy_extra=None, rows=rows, threshold=threshold)

    def sandwich(self, S, rows=None, threshold=None):
        """
        Apply the sandwich formula to the CategoryCov object for a given
        pandas.Series.

        Parameters
        ----------
        S : 1D or 2D iterable
            General sensitivities.
        rows : `int`, optional
            Option to use row calculation for matrix calculations. This option
            defines the number of lines to be taken into account in each loop.
            The default is None.
        threshold : `int`, optional
            Thereshold to avoid numerical fluctuations. The default is None.

        Returns
        -------
        `CategoryCov`
            `CategoryCov` object to which we have applied sandwich
            formula for a given pd.Series

        Warnings
        --------
        The `CategoryCov` object and the sensitivity (S) must have the same
        indices.

        Examples
        --------
        >>> S = np.array([1, 2, 3])
        >>> var = pd.Series([1, 2, 3])
        >>> cov = sandy.CategoryCov.from_var(S)
        >>> cov.sandwich(var)
                    0           1           2
        0 1.00000e+00 0.00000e+00 0.00000e+00
        1 0.00000e+00 8.00000e+00 0.00000e+00
        2 0.00000e+00 0.00000e+00 2.70000e+01

        >>> S = pd.DataFrame(np.diag(np.array([1, 2, 3])), index=pd.Index([1, 2, 3]), columns=pd.Index([4, 5, 6]))
        >>> var = pd.Series([1, 2, 3], index=pd.Index([1, 2, 3]))
        >>> cov = sandy.CategoryCov(S)
        >>> cov.sandwich(var)
        	          1	          2	          3
        1	1.00000e+00	0.00000e+00	0.00000e+00
        2	0.00000e+00	8.00000e+00	0.00000e+00
        3	0.00000e+00	0.00000e+00	2.70000e+01

        >>> S = np.array([1, 2, 3])
        >>> var = pd.Series([1, 2, 3])
        >>> cov = sandy.CategoryCov.from_var(S)
        >>> var = sandy.CategoryCov.from_var(var).data
        >>> cov.sandwich(var)
        	0	1	2
        0	1.00000e+00	0.00000e+00	0.00000e+00
        1	0.00000e+00	8.00000e+00	0.00000e+00
        2	0.00000e+00	0.00000e+00	2.70000e+01
        """
        if pd.DataFrame(S).shape[1] == 1:
            S_ = sandy.CategoryCov.from_var(S).data
            sandwich = self._gls_Vy_calc(S_, rows=rows)
        else:
            S_ = pd.DataFrame(S).T
            sandwich = self._gls_Vy_calc(S_, rows=rows)
        if threshold is not None:
            sandwich[sandwich < threshold] = 0
        return self.__class__(sandwich)

    def plot_corr(self, ax, **kwargs):
        """
        Plot correlation matrix as a color-encoded matrix.

        Parameters
        ----------
        ax : `matplotlib Axes`
            Axes in which to draw the plot, otherwise use the currently-active
            Axes.
        kwargs : `dict`
            keyword arguments for seaborn heatmap plot.

        Returns
        -------
        ax : `matplotlib Axes`
            Axes object with the heatmap.
        """
        add = {"cbar": True, "vmin": -1, "vmax": 1, "cmap": "RdBu"}
        for k, v in kwargs.items():
            add[k] = v
        ax = sns.heatmap(self.get_corr(), ax=ax, **add)
        return ax

    @classmethod
    def corr2cov(cls, corr, std, **kwargs):
        """
        Produce covariance matrix given correlation matrix and standard
        deviation array.
        Same as :obj: `corr2cov` but it returns a :obj: `CategoryCov`.

        Parameters
        ----------
        corr : 2d iterable
            square 2D correlation matrix
        std : 1d iterable
            array of standard deviations
        kwargs: `dict`
            keyword arguments to pass to `data` method.
            They can be used to initialize index and columns of the
            :obj: `CategoryCov`.

        Returns
        -------
        :obj: `CategoryCov`
            covariance matrix

        Examples
        --------
        Initialize index and columns
        >>> idx = ["A", "B", "C"]
        >>> std = np.array([1, 2, 3])
        >>> corr = np.array([[1, 0, 2], [0, 3, 0], [2, 0, 1]])
        >>> sandy.CategoryCov.corr2cov(corr, std, index=idx, columns=idx)
                    A           B           C
        A 1.00000e+00 0.00000e+00 6.00000e+00
        B 0.00000e+00 1.20000e+01 0.00000e+00
        C 6.00000e+00 0.00000e+00 9.00000e+00
        """
        return cls(corr2cov(corr, std), **kwargs)

    @classmethod
    @cov33csv
    def from_csv(cls, file, **kwargs):
        """
        Read covariance matrix from csv file using `pandas.read_csv`.

        Parameters
        ----------
        file: `str`
            csv file containing covariance matrix (with or w/o indices and
            columns)
        kwargs: `dict`
            keyword arguments to pass to `pd.read_csv`

        Returns
        -------
        `CategoryCov`
            object containing covariance matrix

        Examples
        --------
        Read a 2x2 matrix from a string in csv format.
        >>> from io import StringIO
        >>> cov = pd.DataFrame([[1, 0.4],[0.4, 1]])
        >>> string = StringIO(cov.to_csv())
        >>> sandy.CategoryCov.from_csv(string, index_col=0)
                    0           1
        0 1.00000e+00 4.00000e-01
        1 4.00000e-01 1.00000e+00

        Now use `pandas.MultiIndex` as `index` and `columns`.
        This example represents the case of a cross section covariance matrix
        for `MAT=9437`, `MT=18` and two energy points `[1e-5, 1e6]`.
        >>> tuples = [(9437, 18, 1e-5), (9437, 18, 1e6)]
        >>> index = pd.MultiIndex.from_tuples(tuples, names=("MAT", "MT", "E"))
        >>> cov.index = cov.columns = index
        >>> string = StringIO(cov.to_csv())
        >>> pos = [0, 1, 2]
        >>> sandy.CategoryCov.from_csv(string, index_col=pos, header=pos)
        MAT                        9437            
        MT                           18            
        E                         1e-05   1000000.0
        MAT  MT E                                  
        9437 18 1.00000e-05 1.00000e+00 4.00000e-01
                1.00000e+06 4.00000e-01 1.00000e+00
        """
        df = pd.read_csv(file, **kwargs)
        return cls(df)

    @classmethod
    def random_corr(cls, size, correlations=True, seed=None, **kwargs):
        """
        >>> sandy.CategoryCov.random_corr(2, seed=1)
                    0           1
        0 1.00000e+00 4.40649e-01
        1 4.40649e-01 1.00000e+00

        >>> sandy.CategoryCov.random_corr(2, correlations=False, seed=1)
                    0           1
        0 1.00000e+00 0.00000e+00
        1 0.00000e+00 1.00000e+00
        """
        np.random.seed(seed=seed)
        corr = np.eye(size)
        if correlations:
            offdiag = np.random.uniform(-1, 1, size**2).reshape(size, size)
            up = np.triu(offdiag, 1)
        else:
            up = np.zeros([size, size])
        corr += up + up.T
        return cls(corr, **kwargs)

    @classmethod
    def random_cov(cls, size, stdmin=0.0, stdmax=1.0, correlations=True,
                   seed=None, **kwargs):
        """
        Construct a covariance matrix with random values

        Parameters
        ----------
        size : `int`
            Dimension of the original matrix
        stdmin : `float`, default is 0
            minimum value of the uniform standard deviation vector
        stdmax : `float`, default is 1
            maximum value of the uniform standard deviation vector
        correlation : `bool`, default is True
            flag to insert the random correlations in the covariance matrix
        seed : `int`, optional, default is `None`
            seed for the random number generator (by default use `numpy`
            dafault pseudo-random number generator)

        Returns
        -------
        `CategoryCov`
            object containing covariance matrix

        Examples
        --------
        >>> sandy.CategoryCov.random_cov(2, seed=1)
                    0           1
        0 2.15373e-02 5.97134e-03
        1 5.97134e-03 8.52642e-03
        """
        corr = random_corr(size, correlations=correlations, seed=seed)
        std = np.random.uniform(stdmin, stdmax, size)
        return cls.corr2cov(corr, std, **kwargs)

    def to_sparse(self, method='csr_matrix'):
        """
        Method to extract `CategoryCov` values into a sparse matrix

        Parameters
        ----------
        method : `str`, optional
            SciPy 2-D sparse matrix. The default is 'csr_matrix'.

        Methods
        -------
        `csr_matrix`:
            Compressed Sparse Row matrix.
        `bsr_matrix`:
            Block Sparse Row matrix.
        `coo_matrix`:
            A sparse matrix in COOrdinate format.
        `csc_matrix`:
            Compressed Sparse Column matrix.
        `dia_matrix`:
            Sparse matrix with DIAgonal storage.
        `dok_matrix`:
            Dictionary Of Keys based sparse matrix.
        `lil_matrix`:
            Row-based list of lists sparse matrix.

        Returns
        -------
        data_sp : `scipy.sparse.matrix`
            `CategoryCov` instance values stored as a sparse matrix
        """
        data = self.data.values
        if method == 'csr_matrix':
            data_sp = sps.csr_matrix(data)
        elif method == 'bsr_matrix':
            data_sp = sps.bsr_matrix(data)
        elif method == 'coo_matrix':
            data_sp = sps.coo_matrix(data)
        elif method == 'csc_matrix':
            data_sp = sps.csc_matrix(data)
        elif method == 'dia_matrix':
            data_sp = sps.dia_matrix(data)
        elif method == 'dok_matrix':
            data_sp = sps.dok_matrix(data)
        elif method == 'lil_matrix':
            data_sp = sps.lil_matrix(data)
        else:
            raise ValueError('The method does not exist in scipy.sparse')
        return data_sp

    def get_L(self, rows=None, tolerance=None):
        """
        Extract lower triangular matrix `L` for which `L*L^T == self`.

        Parameters
        ----------
        rows : `int`, optional
            Option to use row calculation for matrix calculations. This option
            defines the number of lines to be taken into account in each loop.
            The default is None.
        tolerance : `float`, optional, default is `None`
            replace all eigenvalues smaller than a given tolerance with zeros.

        Returns
        -------
        `pandas.DataFrame`
            Cholesky descomposition low triangular matrix.

        Examples
        --------
        Positive define matrix:
        >>> a = np.array([[4, 12, -16], [12, 37, -43], [-16, -43, 98]])
        >>> sandy.CategoryCov(a).get_L()
                       0	          1	          2
        0	-2.00000e+00	0.00000e+00	0.00000e+00
        1	-6.00000e+00	1.00000e+00	0.00000e+00
        2	 8.00000e+00	5.00000e+00	3.00000e+00

        >>> sandy.CategoryCov(a).get_L(tolerance=0)
                       0	          1	          2
        0	-2.00000e+00	0.00000e+00	0.00000e+00
        1	-6.00000e+00	1.00000e+00	0.00000e+00
        2	 8.00000e+00	5.00000e+00	3.00000e+00

        >>> sandy.CategoryCov(a).get_L(rows=1)
                       0	          1	          2
        0	-2.00000e+00	0.00000e+00	0.00000e+00
        1	-6.00000e+00	1.00000e+00	0.00000e+00
        2	 8.00000e+00	5.00000e+00	3.00000e+00

        Matrix with negative eigenvalues
        >>> sandy.CategoryCov([[1, -2],[-2, 3]]).get_L(rows=1, tolerance=0)
                       0	          1
        0	-1.08204e+00	0.00000e+00
        1	 1.75078e+00	0.00000e+00

        >>> sandy.CategoryCov([[1, -2],[-2, 3]]).get_L(tolerance=0)
                       0	          1
        0	-1.08204e+00	0.00000e+00
        1	 1.75078e+00	0.00000e+00

        Decomposition test:
        >>> L = sandy.CategoryCov(a).get_L()
        >>> L.dot(L.T)
                       0	           1	           2
        0	 4.00000e+00	 1.20000e+01	-1.60000e+01
        1	 1.20000e+01	 3.70000e+01	-4.30000e+01
        2	-1.60000e+01	-4.30000e+01	 9.80000e+01

        Matrix with negative eigenvalues, tolerance of 0:
        >>> L = sandy.CategoryCov([[1, -2],[-2, 3]]).get_L(rows=1, tolerance=0)
        >>> L.dot(L.T)
        	           0	           1
        0	 1.17082e+00	-1.89443e+00
        1	-1.89443e+00	 3.06525e+00
        """
        index = self.data.index
        columns = self.data.columns
        # Reduces the size of the matrix, erasing the zero values
        nonzero_idxs, cov_reduced = reduce_size(self.data)
        # Obtain the eigenvalues and eigenvectors:
        E, V = sandy.CategoryCov(cov_reduced).eig(tolerance=tolerance)
        E = sps.diags(np.sqrt(E)).toarray()
        # Construct the matrix:
        rows_ = cov_reduced.shape[0] if rows is None else rows
        A = sandy.cov.sparse_tables_dot(V, E, rows=rows_).T.toarray()
        # QR decomposition:
        Q, R = scipy.linalg.qr(A)
        L_redu = R.T
        # Original size
        L = restore_size(nonzero_idxs, L_redu, len(self.data)).values
        return pd.DataFrame(L, index=index, columns=columns)


class EnergyCov(CategoryCov):
    """
    Dataframe for a multigroup covariance matrix.

    .. note:: It is assumed that the covariance matrix is defined over
              multi-group energy grids.

              Only 'zero' interpolation is supported.

    Attributes
    ----------
    data : `pandas.DataFrame`
        covariance matrix as a dataframe

    Methods
    -------
    add
    
    change_grid
        
    from_lb1
        
    from_lb2
        
    from_lb5_sym
        
    from_lb5_asym
        
    from_lb6
        
    sum_covs
       
    Raises
    ------
    `sandy.Error`
        if index values are not monotonically increasing
    `sandy.Error`
        if columns values are not monotonically increasing
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    @property
    def data(self):
        """
        Covariance matrix as a dataframe.

        Attributes
        ----------
        index : `pandas.Index`
            indices
        columns : `pandas.MultiIndex`
            indices
        values : `numpy.array`
            covariance values as `float`

        Returns
        -------
        `pandas.DataFrame`
            covariance matrix

        Raises
        ------
        `sandy.Error`
            if `index` or `columns` are not monotonically increasing
        """
        return self._data

    @data.setter
    def data(self, data):
        self._data = pd.DataFrame(data)
        self._data.index = pd.Index(
                self._data.index.values,
                name="E",
                )
        self._data.columns = pd.Index(
                self._data.columns.values,
                name="E",
                )
        if not self._data.index.is_monotonic_increasing:
            raise sandy.Error("index values are not monotonically increasing")
        if not self._data.columns.is_monotonic_increasing:
            raise sandy.Error("columns values are not monotonically "
                              "increasing")

    def change_grid(self, ex, ey, inplace=False):
        """
        Given one energy grid for the x-axis and one energy grid for the
        y-axis, interpolate/extrapolate the covariance matrix over the new
        points using the *forward-filling* method.

        .. important::

            * backward extrapolated values (e.g. below threshold) are replaced
              by 0
            * forward extrapolated values (e.g. above 20 MeV) are replaced by
              the covariance coefficient that refers to the last point in the
              original grid

        Parameters
        ----------
        ex : `iterable`
            energy grid for the x-axis
        ey : `iterable`
            energy grid for the y-axis

        Returns
        -------
        `sandy.EnergyCov`
            Covariance matrix interpolated over the new axes.

        Examples
        --------
        >>> eg = [1e-2, 1e6]
        >>> C = sandy.EnergyCov.random_corr(2, seed=1, index=eg, columns=eg)
        >>> C.change_grid([0, 1, 1e6, 1e7], [0, 1, 1e6, 1e7])
        E            0.00000e+00  1.00000e+00  1.00000e+06  1.00000e+07
        E                                                              
        0.00000e+00  0.00000e+00  0.00000e+00  0.00000e+00  0.00000e+00
        1.00000e+00  0.00000e+00  1.00000e+00  4.40649e-01  4.40649e-01
        1.00000e+06  0.00000e+00  4.40649e-01  1.00000e+00  1.00000e+00
        1.00000e+07  0.00000e+00  4.40649e-01  1.00000e+00  1.00000e+00
        """
        df = self.data.reindex(index=ex, method="ffill") \
                      .reindex(columns=ey, method="ffill") \
                      .fillna(0)
        if not inplace:
            return self.__class__(df)
        self.data = df

    def _plot_matrix(self, ax, xscale='log', yscale='log', cmap='bwr',
                     vmin=-1, vmax=1, emin=1e-5, emax=2e7, **kwargs):
        new_xgrid = np.unique([*self.data.index, *[emin, emax]])
        new_ygrid = np.unique([*self.data.columns, *[emin, emax]])
        data = self.change_grid(ex=new_xgrid, ey=new_ygrid).data
        X, Y = np.meshgrid(data.index.values, data.columns.values)
        qmesh = ax.pcolormesh(
                    X.T,
                    Y.T,
                    data.values,
                    cmap=cmap,
                    vmin=vmin,
                    vmax=vmax,
                    **kwargs,
                    )
        ax.set_xlim([emin, emax])
        ax.set_ylim([emin, emax])
        plt.colorbar(qmesh)
        ax.set_xscale(xscale)
        ax.set_yscale(yscale)
        return ax

    def add(self, cov, inplace=False):
        """
        Add the content of another `EnergyCov` (sum).
        If the energy grids do not match, interpolate.

        Parameters
        ----------
        cov : `sandy.EnergyCov`
            multigroup covariance matrices (axes can be different)
        inplace : `bool`, optional, default is `False`
            flag to operate **inplace**

        Returns
        -------
        `sandy.EnergyCov`
            Multi-group covariance matrix.

        Examples
        --------
        >>> eg = [1e-2, 1e6]
        >>> C = sandy.EnergyCov.random_corr(2, seed=1, index=eg, columns=eg)
        >>> C.add(C)
        E            1.00000e-02  1.00000e+06
        E                                    
        1.00000e-02  2.00000e+00  8.81298e-01
        1.00000e+06  8.81298e-01  2.00000e+00

        >>> eg = [1e-1, 1]
        >>> D = sandy.EnergyCov.random_corr(2, seed=5, index=eg, columns=eg)
        >>> C.add(D)
        E            1.00000e-02  1.00000e-01  1.00000e+00  1.00000e+06
        E                                                              
        1.00000e-02  1.00000e+00  1.00000e+00  1.00000e+00  4.40649e-01
        1.00000e-01  1.00000e+00  2.00000e+00  1.74146e+00  1.18211e+00
        1.00000e+00  1.00000e+00  1.74146e+00  2.00000e+00  1.44065e+00
        1.00000e+06  4.40649e-01  1.18211e+00  1.44065e+00  2.00000e+00

        >>> assert C.add(D).data.equals(D.add(C).data)
        """
        ex = np.unique([*self.data.index, *cov.data.index])
        ey = np.unique([*self.data.columns, *cov.data.columns])
        x = self.change_grid(ex, ey)
        y = cov.change_grid(ex, ey)
        data = x.data.add(y.data)
        if inplace:
            self.data = data
        else:
            return self.__class__(data)

    @classmethod
    def sum_covs(cls, *covs):
        """
        Sum multigroup covariance matrices into a single one.

        Parameters
        ----------
        covs : iterable of `sandy.EnergyCov`
            list of multigroup covariance matrices (axes can be different)

        Returns
        -------
        `sandy.EnergyCov`
            Multi-group covariance matrix.

        Examples
        --------
        Sum two 2x2 correlation matrices with different indices and columns
        >>> eg = [1e-2, 1e6]
        >>> C = sandy.EnergyCov.random_corr(2, seed=1, index=eg, columns=eg)
        >>> eg = [1e-1, 1]
        >>> D = sandy.EnergyCov.random_corr(2, seed=5, index=eg, columns=eg)
        >>> sandy.EnergyCov.sum_covs(C, D)
        E            1.00000e-02  1.00000e-01  1.00000e+00  1.00000e+06
        E                                                              
        1.00000e-02  1.00000e+00  1.00000e+00  1.00000e+00  4.40649e-01
        1.00000e-01  1.00000e+00  2.00000e+00  1.74146e+00  1.18211e+00
        1.00000e+00  1.00000e+00  1.74146e+00  2.00000e+00  1.44065e+00
        1.00000e+06  4.40649e-01  1.18211e+00  1.44065e+00  2.00000e+00
        """
        return functools.reduce(lambda x, y: x.add(y), covs)

    @classmethod
    def from_lb1(cls, evalues, fvalues):
        """Extract square covariance matrix from NI-type sub-subsection data 
        with flag `lb=1`.
        
        Parameters
        ----------
        evalues : iterable
            covariance energy grid (same for both axes)
        fvalues : iterable
            array of F-values (covriance matrix diagonal)
        
        Returns
        -------
        `sandy.EnergyCov`
            Multi-group covariance matrix.
        """
        cov = np.diag(fvalues)
        return cls(cov, index=evalues, columns=evalues)

    @classmethod
    def from_lb2(cls, evalues, fvalues):
        """Extract square covariance matrix from NI-type sub-subsection data 
        with flag `lb=2`.
        
        Parameters
        ----------
        evalues : `iterable`
            covariance energy grid for both axis
        fvalues : `iterable`
            array of F-values
        
        Returns
        -------
        `sandy.formats.utils.EnergyCov`
            Multi-group covariance matrix.
        """
        f = np.array(fvalues)
        cov = f*f.reshape(-1,1)
        return cls(cov, index=evalues, columns=evalues)

    @classmethod
    def from_lb5_sym(cls, evalues, fvalues):
        """Extract square symmetric covariance matrix from NI-type sub-subsection data 
        with flag `lb=5`.
        
        Parameters
        ----------
        evalues : `iterable`
            covariance energy grid for both axis
        fvalues : `iterable`
            array of F-values (flattened upper triangular matrix coefficients)
        
        Returns
        -------
        `sandy.formats.utils.EnergyCov`
            Multi-group covariance matrix.
        """
        ne = len(evalues)
        cov = np.zeros([ne - 1, ne - 1])
        indices = np.triu_indices(ne - 1)
        cov[indices] = np.array(fvalues)
        cov += np.triu(cov, 1).T
        # add zero row and column at the end of the matrix
        cov = np.insert(cov, cov.shape[0], [0]*cov.shape[1], axis=0)
        cov = np.insert(cov, cov.shape[1], [0]*cov.shape[0], axis=1)
        return cls(cov, index=evalues, columns=evalues)

    @classmethod
    def from_lb5_asym(cls, evalues, fvalues):
        """
        Extract square asymmetric covariance matrix from NI-type sub-subsection data 
        with flag `lb=5`.
        
        Parameters
        ----------
        evalues : `iterable`
            covariance energy grid for both axis
        fvalues : `iterable`
            array of F-values (flattened full matrix)
        
        Returns
        -------
        `sandy.formats.utils.EnergyCov`
            Multi-group covariance matrix.
        """
        ne = len(evalues)
        cov = np.array(fvalues).reshape(ne - 1, ne - 1)
        # add zero row and column at the end of the matrix
        cov = np.insert(cov, cov.shape[0], [0]*cov.shape[1], axis=0)
        cov = np.insert(cov, cov.shape[1], [0]*cov.shape[0], axis=1)
        return cls(cov, index=evalues, columns=evalues)

    @classmethod
    def from_lb6(cls, evalues_r, evalues_c, fvalues):
        """Extract covariance matrix from NI-type sub-subsection data 
        with flag `lb6`.
        
        Parameters
        ----------
        evalues_r : `iterable`
            covariance energy grid for row axis
        evalues_c : `iterable`
            covariance energy grid for column axis
        fvalues : `iterable`
            array of F-values (flattened full matrix)
        
        Returns
        -------
        `sandy.formats.utils.EnergyCov`
            Multi-group covariance matrix.
        """
        ner = len(evalues_r)
        nec = len(evalues_c)
        cov = np.array(fvalues).reshape(ner-1, nec-1)
        # add zero row and column at the end of the matrix
        cov = np.insert(cov, cov.shape[0], [0]*cov.shape[1], axis=0)
        cov = np.insert(cov, cov.shape[1], [0]*cov.shape[0], axis=1)
        return cls(cov, index=evalues_r, columns=evalues_c)


class GlobalCov(CategoryCov):

    @classmethod
    def from_list(cls, iterable):
        """
        Extract global cross section/nubar covariance matrix from iterables 
        of `EnergyCovs`.
        
        Parameters
        ----------
        iterable : iterable
            list of tuples/lists/iterables with content `[mat, mt, mat1, mt1, EnergyCov]`
        
        Returns
        -------
        `XsCov` or `pandas.DataFrame`
            global cross section/nubar covariance matrix (empty dataframe if no covariance matrix was found)
        """
        columns = ["KEYS_ROWS", "KEYS_COLS", "COV"]
        # Reindex the cross-reaction matrices
        covs = pd.DataFrame.from_records(iterable).set_axis(columns, axis=1).set_index(columns[:-1]).COV
        for (keys_rows,keys_cols), cov in covs.iteritems():
            if keys_rows == keys_cols: # diagonal terms
                if cov.data.shape[0] != cov.data.shape[1]:
                    raise SandyError("non-symmetric covariance matrix for ({}, {})".format(keys_rows, keys_cols))
                if not np.allclose(cov.data, cov.data.T):
                    raise SandyError("non-symmetric covariance matrix for ({}, {})".format(keys_rows, keys_cols))
            else: # off-diagonal terms
                condition1 = (keys_rows,keys_rows) in covs.index
                condition2 = (keys_cols,keys_cols) in covs.index
                if not (condition1 and condition2):
                    covs[keys_rows,keys_cols] = np.nan
                    logging.warn("skip covariance matrix for ({}, {})".format(keys_rows, keys_cols))
                    continue
                ex = covs[keys_rows,keys_rows].data.index.values
                ey = covs[keys_cols,keys_cols].data.columns.values
                covs[keys_rows,keys_cols] = cov.change_grid(ex, ey)
        covs.dropna(inplace=True)
        if covs.empty:
            logging.warn("covariance matrix is empty")
            return pd.DataFrame()
        # Create index for global matrix
        rows_levels = covs.index.levels[0]
        indexlist = [(*keys,e) for keys in rows_levels for e in covs[(keys,keys)].data.index.values]
        index = pd.MultiIndex.from_tuples(indexlist, names=cls.labels)
        # Create global matrix
        matrix = np.zeros((len(index),len(index)))
        for (keys_rows,keys_cols), cov in covs.iteritems():
            ix = index.get_loc(keys_rows)
            ix1 = index.get_loc(keys_cols)
            matrix[ix.start:ix.stop,ix1.start:ix1.stop] = cov.data
            if keys_rows != keys_cols:
                matrix[ix1.start:ix1.stop,ix.start:ix.stop] = cov.data.T
        return cls(matrix, index=index, columns=index)


def corr2cov(corr, s):
    """
    Produce covariance matrix given correlation matrix and standard
    deviation array.

    Parameters
    ----------
    corr : 2D iterable
        square 2D correlation matrix
    s : 1D iterable
        1D iterable with standard deviations

    Returns
    -------
    `numpy.ndarray`
        square 2D covariance matrix

    Examples
    --------
    Test with integers
    >>> s = np.array([1, 2, 3])
    >>> corr = np.array([[1, 0, 2], [0, 3, 0], [2, 0, 1]])
    >>> corr2cov(corr, s)
    array([[ 1,  0,  6],
           [ 0, 12,  0],
           [ 6,  0,  9]])

    Test with float
    >>> corr2cov(corr, s.astype(float))
    array([[ 1.,  0.,  6.],
           [ 0., 12.,  0.],
           [ 6.,  0.,  9.]])
    """
    s_ = np.diag(s)
    return s_.dot(corr.dot(s_))


def sparse_tables_dot(a, b, rows=1000):
    """
    Function to perform multiplications between matrices stored on local
    disk instead of memory.

    Parameters
    ----------
    a : 2D iterable
        Matrix.
    b : 2D iterable
        Matrix.
    rows : `int`, optional.
        Number of rows to be calculated in each loop. The default is 1000.

    Returns
    -------
    dot_product : "scipy.sparse.csc_matrix"
        The multiplication of 2 matrix.
    """
    a_ = sps.csr_matrix(a)
    b_ = sps.csc_matrix(b)
    l, n = a_.shape[0], b_.shape[1]
    f = tb.open_file('dot.h5', 'w')
    filters = tb.Filters(complevel=5, complib='blosc')
    out_data = f.create_earray(f.root, 'data', tb.Float64Atom(), shape=(0,),
                               filters=filters)
    out_indices = f.create_earray(f.root, 'indices', tb.Int64Atom(), shape=(0,),
                                  filters=filters)
    out_indptr = f.create_earray(f.root, 'indptr', tb.Int64Atom(), shape=(0,),
                                 filters=filters)
    out_indptr.append(np.array([0]))
    max_indptr = 0
    for i in range(0, l, rows):
        res = a_[i:min(i+rows, l), :].dot(b_)
        out_data.append(res.data)
        indices = res.indices
        indptr = res.indptr
        out_indices.append(indices)
        out_indptr.append(max_indptr+indptr[1:])
        max_indptr += indices.shape[0]
        f.flush()
    dot_product = sps.csr_matrix((f.root.data[:], f.root.indices[:],
                                  f.root.indptr[:]), shape=(l, n))
    f.close()
    os.remove('dot.h5')
    return dot_product


def sparse_tables_dot_multiple(matrix_list, rows=1000):
    """
    Function to perform multiplications between matrices stored on local
    disk instead of memory.

    Parameters
    ----------
    matrix_list : 1D iterables
        Iterable with the matrix inside
    rows : `int`, optional.
        Number of rows to be calculated in each loop. The default is 1000.

    Returns
    -------
    matrix : "scipy.sparse.csc_matrix"
        Result of the multiplication of 2 matrix.

    Example
    -------
    >>> S = np.array([[1, 2], [3, 4]])
    >>> cov = sandy.CategoryCov.from_var([1, 1]).data.values
    >>> sparse_tables_dot_multiple([S, cov, S.T], 1)
    array([[ 5., 11.],
           [11., 25.]])
    """
    matrix = matrix_list[0]
    for b in matrix_list[1::]:
        intermediate_matrix = sparse_tables_dot(matrix, b, rows=rows)
        matrix = intermediate_matrix
    return matrix.toarray()


def sparse_tables_inv(a, rows=1000):
    """
    Function to perform matrix inversion stored on local
    disk instead of memory.

    Parameters
    ----------
    a : 2D iterable
        Matrix to be inverted.
    rows : `int`, optional.
        Number of rows to be calculated in each loop. The default is 1000.

    Returns
    -------
    invert_matrix : "numpy.ndarray"
        The inverted matrix.

    Example
    -------
    >>> S = sandy.CategoryCov(np.diag(np.array([1, 2, 3]))).data.values
    >>> sandy.cov.sparse_tables_inv(S, 1).round(2)
    array([[1.  , 0.  , 0.  ],
           [0.  , 0.5 , 0.  ],
           [0.  , 0.  , 0.33]])
    """
    a_ = sps.csc_matrix(a)
    l, n = a_.shape[0], a_.shape[1]
    LU = spsl.splu(a_)
    f = tb.open_file('inv.h5', 'w')
    filters = tb.Filters(complevel=5, complib='blosc')
    out_data = f.create_carray(f.root, 'data', tb.Float64Atom(),
                               shape=(l, n), filters=filters)
    Identity = np.hsplit(sps.diags([1], shape=a_.shape).toarray(), l/rows)
    j = 0
    for i in range(0, l, rows):
        I_split = Identity[j]
        tmpResult = LU.solve(I_split)
        out_data[:, i:min(i+rows, l)] = tmpResult
        f.flush()
        j += 1
    invert_matrix = sps.csr_matrix(out_data).toarray()
    f.close()
    os.remove('inv.h5')
    return invert_matrix


def segmented_pivot_table(data_stack, index, columns, values, rows=10000000):
    """
    Create a pivot table from a stacked dataframe.

    Parameters
    ----------
    data_stack : `pd.Dataframe`
        Stacked dataframe.
    index : 1D iterable, optional
        Index of the final covariance matrix.
    columns : 1D iterable, optional
        Columns of the final covariance matrix.
    values : `str`, optional
        Name of the column where the values are located.
    rows : `int`, optional
        Number of rows to take into account into each loop. The default
        is 10000000.

    Returns
    -------
    pivot_matrix : `pd.DataFrame`
        Covariance matrix created from a stacked data

    Examples
    --------
    >>> S = pd.DataFrame(np.array([[1, 1, 1], [0, 2, 1], [0, 0, 1]]))
    >>> S = S.stack().reset_index().rename(columns = {'level_0': 'dim1', 'level_1': 'dim2', 0: 'cov'})
    >>> sandy.cov.segmented_pivot_table(S, index=['dim1'], columns=['dim2'], values='cov')
    dim2	0	1	2
    dim1
       0	1	1	1
       1	0	2	1
       2	0	0	1

    >>> sandy.cov.segmented_pivot_table(S, index=['dim1'], columns=['dim2'], values='cov', rows=1)
    dim2	0	        1	        2
    dim1
       0    1.00000e+00 1.00000e+00 1.00000e+00
       1    0.00000e+00 2.00000e+00 1.00000e+00
       2    0.00000e+00 0.00000e+00 1.00000e+00
    """
    size = data_stack.shape[0]
    pivot_matrix = []
    for i in range(0, size, rows):
        partial_pivot = data_stack[i: min(i+rows, size)].pivot_table(
            index=index,
            columns=columns,
            values=values,
            fill_value=0,
            aggfunc=np.sum,
            )
        pivot_matrix.append(partial_pivot)
    pivot_matrix = pd.concat(pivot_matrix).fillna(0)
    # Because the default axis to concatenate is the 0, some duplicate
    # index appear with null values. With this groupby, the duplicate axis
    # disappear, keeping the original values.
    pivot_matrix = pivot_matrix.groupby(pivot_matrix.index).sum()
    if len(index) >= 2:
        # Groupby transforms multiindex structure into a tuple. This line
        # reverse the transformation.
        pivot_matrix.index = pd.MultiIndex.from_tuples(
            pivot_matrix.index,
            names=index,
            )
    return pivot_matrix

def triu_matrix(matrix, kind='upper'):
    """
    Given the upper or lower triangular matrix , return the full symmetric
    matrix.

    Parameters
    ----------
    matrix : 2d iterable
        Upper triangular matrix
    kind : `str`, optional
        Select if matrix variable is upper or lower triangular matrix. The
        default is 'upper'

    Returns
    -------
    `pd.Dataframe`
        reconstructed symmetric matrix

    Examples
    --------
    >>> S = pd.DataFrame(np.array([[1, 2, 1], [0, 2, 4], [0, 0, 3]]))
    >>> triu_matrix(S).data
                0           1           2
    0 1.00000e+00 2.00000e+00 1.00000e+00
    1 2.00000e+00 2.00000e+00 4.00000e+00
    2 1.00000e+00 4.00000e+00 3.00000e+00

    Overwrite the lower triangular part of the matrix:
    >>> S = pd.DataFrame(np.array([[1, 2, 1], [-8, 2, 4], [-6, -5, 3]]))
    >>> triu_matrix(S).data
                0           1           2
    0 1.00000e+00 2.00000e+00 1.00000e+00
    1 2.00000e+00 2.00000e+00 4.00000e+00
    2 1.00000e+00 4.00000e+00 3.00000e+00

    Test for lower triangular matrix:
    >>> S = pd.DataFrame(np.array([[3, 0, 0], [5, 2, 0], [1, 2, 1]]))
    >>> triu_matrix(S, kind='lower').data
                0           1           2
    0 3.00000e+00 5.00000e+00 1.00000e+00
    1 5.00000e+00 2.00000e+00 2.00000e+00
    2 1.00000e+00 2.00000e+00 1.00000e+00
    
    Overwrite the upper triangular part of the matrix:
    >>> S = pd.DataFrame(np.array([[3, 5, -9], [5, 2, 8], [1, 2, 1]]))
    >>> triu_matrix(S, kind='lower').data
                0           1           2
    0 3.00000e+00 5.00000e+00 1.00000e+00
    1 5.00000e+00 2.00000e+00 2.00000e+00
    2 1.00000e+00 2.00000e+00 1.00000e+00
    """
    matrix_ = pd.DataFrame(matrix)
    index = matrix_.index
    columns = matrix_.columns
    values = matrix_.values
    if kind == 'upper':    
        index_lower = np.tril_indices(matrix_.shape[0], -1)
        values[index_lower] = values.T[index_lower]
    elif kind == 'lower':
        index_upper = np.triu_indices(matrix_.shape[0], 1)
        values[index_upper] = values.T[index_upper]
    return CategoryCov(pd.DataFrame(values, index=index, columns=columns))


def sample_distribution(dim, nsmp, seed=None, pdf='normal'):
    """
    Extract random independent and identically distributed samples according to
    the chosen distribution with standard deviation=1 and mean=1.

    Parameters
    ----------
    dim : `int`
        Dimension of the matrix from where we obtain the samples.
    nsmp : `int`
        number of samples.
    seed : `int`, optional, default is `None`
        seed for the random number generator (by default use `numpy`
        dafault pseudo-random number generator).
    pdf : `str`, optional
        Random numbers distribution. The default is 'normal'
        Available distributions are:
            * `'normal'`
            * `'uniform'`
            * `'lognormal'`

    Returns
    -------
    y : `np.array`
        Numpy array with the random numbers.
        
    Notes
    -----
    .. note:: the implementation of the lognormal distribution sampling is
            performed with the following equations:
            ..math::
                $$
                \mu_{LN} = \exp(\mu_N + \frac{\sigma_N^2}{2})\\
                \sigma_{LN} = \sqrt{[\exp(2\mu_N + 2\sigma_N^2) - \exp(2\mu_N + \sigma_N^2)]}
                $$
            More details can be found in reference: https://linkinghub.elsevier.com/retrieve/pii/S0168900213008450
                                                    DOI: 10.1016/J.NIMA.2013.06.025
   
    Examples
    --------
    >>> sandy.cov.sample_distribution(2, 3, seed=11)
    array([[ 2.74945474,  0.713927  ,  0.51543487],
           [-1.65331856,  0.99171537,  0.68036864]])

    >>> sandy.cov.sample_distribution(2, 3, seed=11, pdf='uniform')
    array([[-0.10757829, -0.66458659,  0.87258524],
           [ 1.77919399,  0.72357718,  0.94951799]])
    
    >>> sandy.cov.sample_distribution(2, 3, seed=11, pdf='lognormal')
    array([[3.03418551, 0.55724795, 0.4723663 ],
           [0.07764515, 0.70224636, 0.54189439]])

    >>> sandy.cov.sample_distribution(2, 1000000, seed=11).mean().round(5)
    1.00025

    >>> sandy.cov.sample_distribution(2, 1000000, seed=11, pdf='uniform').mean().round(5)
    0.99885

    >>> sandy.cov.sample_distribution(2, 1000000, seed=11, pdf='lognormal').mean().round(5)
    0.99953

    >>> sandy.cov.sample_distribution(2, 1000000, seed=11).std().round(5)
    0.99919

    >>> sandy.cov.sample_distribution(2, 1000000, seed=11, pdf='uniform').std().round(5)
    1.00038

    >>> sandy.cov.sample_distribution(2, 1000000, seed=11, pdf='lognormal').std().round(5)
    0.99756

    >>> np.corrcoef(sandy.cov.sample_distribution(2, 1000000, seed=11)).round(5)
    array([[1.e+00, 5.e-05],
           [5.e-05, 1.e+00]])

    >>> np.corrcoef(sandy.cov.sample_distribution(2, 1000000, seed=11, pdf='uniform')).round(5)
    array([[ 1.0e+00, -9.2e-04],
           [-9.2e-04,  1.0e+00]])

    >>> np.corrcoef(sandy.cov.sample_distribution(2, 1000000, seed=11, pdf='lognormal')).round(5)
    array([[ 1.0e+00, -2.6e-04],
           [-2.6e-04,  1.0e+00]])
    """
    np.random.seed(seed=seed)
    if pdf == 'normal':
        y = np.random.randn(dim, nsmp) + 1
    elif pdf == 'uniform':
        a = np.sqrt(12) / 2
        y = np.random.uniform(-a, a, (dim, nsmp)) + 1
    elif pdf == 'lognormal':
        sl = ml = 1 # target mean and standard deviation of the lognormal distribution
        mn = 2 * np.log(ml) - .5 * np.log(sl**2 + np.exp(2 * np.log(ml))) # required mean of the corresponding normal distibution (note reference in the docstring)
        sn = np.sqrt(2 * (np.log(ml) - mn)) # required standard deviation of the corresponding normal distibution (note reference in the docstring)
        y = np.random.lognormal(mn, sn, (dim, nsmp))
    return y


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


def random_corr(size, correlations=True, seed=None):
    np.random.seed(seed=seed)
    corr = np.eye(size)
    if correlations:
        up = np.triu(np.random.uniform(-1, 1, size**2).reshape(size, size), 1)
    else:
        up = np.zeros([size, size])
    corr += up + up.T
    return corr


def random_cov(size, stdmin=0.0, stdmax=1.0, correlations=True, seed=None):
    corr = random_corr(size, correlations=correlations , seed=seed)
    std = np.random.uniform(stdmin, stdmax, size)
    return corr2cov(corr, std)


def random_ctg_cov(index, stdmin=0.0, stdmax=1.0, correlations=True, seed=None):
    cov = random_cov(len(index), stdmin=stdmin, stdmax=stdmax, correlations=correlations, seed=seed)
    return pd.DataFrame(cov, index=index, columns=index)


def print_matrix(size, triu_matrices):
    coords = list(zip(*np.triu_indices(size)))
    kwargs = {"cbar": False, "vmin": -1, "vmax": 1, "cmap": "RdBu"}
    fig,ax = plt.subplots(size, size, sharex="col", sharey="row")
    for (i,j), m in zip(coords, MM):
        if m is None:
            continue
        ax[i,j] = sns.heatmap(m, ax=ax[i,j], **kwargs)
        if i != j:
            ax[j,i] = sns.heatmap(m.T, ax=ax[j,i], **kwargs)
    return fig, ax
