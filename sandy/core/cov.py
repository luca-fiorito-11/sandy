import functools

import numpy as np
import scipy
import scipy.linalg
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

import sandy
import scipy
from scipy.sparse import csc_matrix
from scipy.sparse.linalg import splu
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
        "get_eig",
        ]

S = np.array([[1, 1, 1],
              [1, 2, 1],
              [1, 3, 1]])
var = np.array([[0, 0, 0],
                [0, 2, 0],
                [0, 0, 3]])


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
        """Extract random samples from the covariance matrix, either using
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

    Attributes
    ----------
    data
        covariance matrix as a dataframe
    std
        series of standard deviations
    corr
        correlation matrix
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
        """
        return self._data

    @data.setter
    def data(self, data):
        self._data = pd.DataFrame(data, dtype=float)
        if not len(data.shape) == 2 and data.shape[0] == data.shape[1]:
            raise TypeError("covariance matrix must have two dimensions")

    @property
    def size(self):
        return self.data.values.shape[0]

    def eig(self, sort=True):
        """
        Extract eigenvalues and eigenvectors.

        Parameters
        ----------
        sort : `bool`, optional, default is `True`
            flag to return sorted eigenvalues and eigenfunctions

        Returns
        -------
        `Pandas.Series`
            real part of eigenvalues sorted in descending order
        `pandas.DataFrame`
            matrix of eigenvectors

        Examples
        --------
        Extract eigenvalues of covariance matrix.
        >>> sandy.CategoryCov([[1, 0.4],[0.4, 1]]).eig()[0]
        0   1.40000e+00
        1   6.00000e-01
        Name: eigenvalues, dtype: float64

        Extract eigenfunctions of covariance matrix.
        >>> sandy.CategoryCov([[1, 0.4],[0.4, 1]]).eig()[1]
                    0            1
        0 7.07107e-01 -7.07107e-01
        1 7.07107e-01  7.07107e-01
        """
        return get_eig(self.data.values, sort=sort)

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
        cov = self.data.values
        std = np.sqrt(np.diag(cov))
        return pd.Series(std, index=self.data.index, name="std")

    @property
    def corr(self):
        """
        Extract correlation matrix.

        Returns
        -------
        `pandas.DataFrame`
            correlation matrix

        Examples
        --------
        >>> sandy.CategoryCov([[4, 2.4],[2.4, 9]]).corr
                    0           1
        0 1.00000e+00 4.00000e-01
        1 4.00000e-01 1.00000e+00
        """
        cov = self.data.values
        with np.errstate(divide='ignore', invalid='ignore'):
            coeff = np.true_divide(1, self.std.values)
            coeff[~ np.isfinite(coeff)] = 0   # -inf inf NaN
        corr = np.multiply(np.multiply(cov, coeff).T, coeff)
        return pd.DataFrame(corr,
                            index=self.data.index,
                            columns=self.data.columns,
                            )

    def _reduce_size(self):
        """
        Reduces the size of the matrix, erasing the null values.

        Returns
        -------
        nonzero_idxs : `numpy.ndarray`
            The indices of the diagonal that are not null.
        cov_reduced : `pandas.DataFrame`
            The reduced matrix.

        Examples
        --------
        >>> S = sandy.CategoryCov(np.diag(np.array([1, 2, 3])))
        >>> non_zero_index, reduce_matrix = S._reduce_size()
        >>> non_zero_index
        array([0, 1, 2], dtype=int64)
        >>> reduce_matrix
                      0	          1	          2
        0	1.00000e+00	0.00000e+00	0.00000e+00
        1	0.00000e+00	2.00000e+00	0.00000e+00
        2	0.00000e+00	0.00000e+00	3.00000e+00

        >>> S = sandy.CategoryCov(np.diag(np.array([0, 2, 3])))
        >>> non_zero_index, reduce_matrix = S._reduce_size()
        >>> non_zero_index
        array([1, 2], dtype=int64)
        >>> reduce_matrix
                      1	          2
        1	2.00000e+00	0.00000e+00
        2	0.00000e+00	3.00000e+00
        """
        nonzero_idxs = np.flatnonzero(np.diag(self.data))
        cov_reduced = self.data.loc[nonzero_idxs, nonzero_idxs]
        return nonzero_idxs, cov_reduced

    @classmethod
    def _restore_size(cls, nonzero_idxs, cov_reduced, dim):
        """
        Restore the size of the matrix

        Parameters
        ----------
        nonzero_idxs : `numpy.ndarray`
            The indices of the diagonal that are not null.
        cov_reduced : `numpy.ndarray`
            The reduced matrix.
        dim : `int`
            Dimension of the original matrix.

        Returns
        -------
        cov : `CategoryCov`
            Matrix of specified dimensions.

        Notes
        -----
        ..notes:: This method was developed to be used after calling
                  `_reduce_size`.

        Examples
        --------
        >>> S = sandy.CategoryCov(np.diag(np.array([0, 2, 3, 0])))
        >>> S
                    0           1           2           3
        0 0.00000e+00 0.00000e+00 0.00000e+00 0.00000e+00
        1 0.00000e+00 2.00000e+00 0.00000e+00 0.00000e+00
        2 0.00000e+00 0.00000e+00 3.00000e+00 0.00000e+00
        3 0.00000e+00 0.00000e+00 0.00000e+00 0.00000e+00

        >>> M_nonzero_idxs, M_reduce = S._reduce_size()
        >>> M_reduce[::] = 1
        >>> M_reduce
                      1	          2
        1	1.00000e+00	1.00000e+00
        2	1.00000e+00	1.00000e+00

        >>> sandy.CategoryCov._restore_size(M_nonzero_idxs, M_reduce.values, len(S.data))
                    0           1           2           3
        0 0.00000e+00 0.00000e+00 0.00000e+00 0.00000e+00
        1 0.00000e+00 1.00000e+00 1.00000e+00 0.00000e+00
        2 0.00000e+00 1.00000e+00 1.00000e+00 0.00000e+00
        3 0.00000e+00 0.00000e+00 0.00000e+00 0.00000e+00
        """
        cov = np.zeros((dim, dim))
        for i, ni in enumerate(nonzero_idxs):
            cov[ni, nonzero_idxs] = cov_reduced[i]
        return cls(cov)

    def invert(self):
        """
        Method for calculating the inverse matrix.

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
        """
        index, columns = self.data.index, self.data.columns
        M_nonzero_idxs, M_reduce = self._reduce_size()
        unit = np.identity(len(M_reduce))
        M_reduce_inv = splu(csc_matrix(M_reduce)).solve(unit)
        data = CategoryCov._restore_size(M_nonzero_idxs, M_reduce_inv,
                                         len(self.data))
        M_inv = pd.DataFrame(data.data,
                             index=index, columns=columns)
        return self.__class__(M_inv)

    def sampling(self, nsmp, seed=None):
        """
        Extract random samples from normali distribution centered in zero
        and with given covariance matrix.

        Parameters
        ----------
        nsmp : `int`
            number of samples
        seed : `int`, optional, default is `None`
            seed for the random number generator (by default use `numpy`
            dafault pseudo-random number generator)

        Returns
        -------
        `sandy.Samples`
            object containing samples

        Examples
        --------
        Draw 3 sets of samples using custom seed.
        >>> sandy.CategoryCov([[1, 0.4],[0.4, 1]]).sampling(3, seed=11)
                     0            1
        0 -1.74945e+00 -3.13159e+00
        1  2.86073e-01  1.06836e-01
        2  4.84565e-01 -9.91209e-02
        """
        np.random.seed(seed=seed)
        dim = self.data.shape[0]
        y = np.random.randn(dim, nsmp)  # normal pdf
        L = self.decompose()
        samples = L.dot(y)
        df = pd.DataFrame(samples,
                          index=self.data.index,
                          columns=list(range(nsmp)),
                          )
        return sandy.Samples(df.T)

    def decompose(self):
        """
        Extract lower triangular matrix `L` for which `L*L^T == COV`.

        Returns
        -------
        `numpy.ndarray`
            lower triangular matrix

        Example
        -------
        >>> sandy.CategoryCov([[1, 0.4],[0.4, 1]]).decompose()
        array([[-1.        ,  0.        ],
               [-0.4       ,  0.91651514]])
        """
        E, V = self.eig(sort=False)
        E[E <= 0] = 0
        M = V.values.dot(np.diag(np.sqrt(E.values)))
        Q, R = scipy.linalg.qr(M.T)
        return R.T

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
        cov = pd.DataFrame(np.diag(var_),
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
        return CategoryCov.from_var(var)

    def _gls_Vy_calc(self, S):
        """
        2D calculated output using S.T.dot(Vx_prior).dot(S)

        Parameters
        ----------
        S : 2D iterable
            Sensitivity square matrix (MXM).

        Returns
        -------
        `pd.DataFrame`
            Vy_calc calculated using S.T.dot(Vx_prior).dot(S)

        Example
        -------
        >>> S = np.array([[1, 2], [3, 4]])
        >>> sensitivity = sandy.CategoryCov.from_var([1, 1])
        >>> sensitivity._gls_Vy_calc(S)
                      0	          1
        0	5.00000e+00	1.10000e+01
        1	1.10000e+01	2.50000e+01
        """
        index = pd.DataFrame(S).columns
        S_ = pd.DataFrame(S).values
        Vx_prior = self.data.values
        Vy_calc = S_.dot(Vx_prior).dot(S_.T)
        return pd.DataFrame(Vy_calc, index=index, columns=index)

    def _gls_M(self, S, Vy_extra):
        """
        2D calculated output using S.dot(Vx_prior).dot(S.T) + Vy_extra

        Parameters
        ----------
        S : 2D iterable
            Sensitivity matrix (MXN).
        Vy_extra : 2D iterable
            2D covariance matrix for y_extra (MXM).

        Returns
        -------
        `pd.DataFrame`
            M calculated using S.T.dot(Vx_prior).dot(S) + Vy_extra

        Example
        -------
        >>> S = np.array([[1, 2], [3, 4]])
        >>> sensitivity = sandy.CategoryCov.from_var([1, 1])
        >>> Vy = np.diag(pd.Series([1, 1]))
        >>> sensitivity._gls_M(S, Vy)
                    	0	      1
        0	6.00000e+00	1.10000e+01
        1	1.10000e+01	2.60000e+01
        """
        index = pd.DataFrame(Vy_extra).index
        columns = pd.DataFrame(Vy_extra).columns
        Vy_extra_ = sandy.CategoryCov(Vy_extra).data.values
        # GLS_sensitivity:
        Vy_calc = self._gls_Vy_calc(S).values
        M = Vy_calc + Vy_extra_
        return pd.DataFrame(M, index=index, columns=columns)

    def _gls_general_sensitivity(self, S, Vy_extra, threshold=None):
        """
        Method to obtain general sensitivity according to GLS

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
        >>> sensitivity = sandy.CategoryCov.from_var([1, 1])
        >>> Vy = np.diag(pd.Series([1, 1]))
        >>> sensitivity._gls_general_sensitivity(S, Vy)
                      0	              1
        0	-2.00000e-01	2.00000e-01
        1	2.28571e-01	    5.71429e-02

        >>> S = pd.DataFrame([[1, 2], [3, 4]], columns=[1, 2],index=[3, 4])
        >>> sensitivity = sandy.CategoryCov.from_var([1, 1])
        >>> Vy = pd.DataFrame([[1, 0], [0, 1]], index=[1, 2], columns=[1, 2])
        >>> sensitivity._gls_general_sensitivity(S, Vy)
                      1	              2
        1	-2.00000e-01	2.00000e-01
        2	2.28571e-01	    5.71429e-02
        """
        index, columns = pd.DataFrame(S).columns, pd.DataFrame(S).columns
        S_ = pd.DataFrame(S).values
        Vx_prior = self.data.values
        # GLS_sensitivity:
        M = self._gls_M(S, Vy_extra).values
        M_inv = sandy.CategoryCov(M).invert()
        sensitivity = Vx_prior.dot(S_.T).dot(M_inv.data.values)
        if threshold is not None:
            sensitivity[sensitivity < threshold] = 0
        return pd.DataFrame(sensitivity, index=index, columns=columns)

    def _gls_cov_sensitivity(self, S, Vy_extra, threshold=None):
        """
        Method to obtain covariance sensitivity according to GLS

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
        `CategoryCov`
            GlS sensitivity for a given Vy and S.

        Example
        -------
        >>> S = np.array([[1, 2], [3, 4]])
        >>> var = sandy.CategoryCov.from_var([1, 1])
        >>> Vy = np.diag(pd.Series([1, 1]))
        >>> var._gls_cov_sensitivity(S, Vy)
                      0	          1
        0	4.00000e-01	4.00000e-01
        1	4.00000e-01	6.85714e-01

        >>> S = pd.DataFrame([[1, 2], [3, 4]], columns =[1, 2],index=[3, 4])
        >>> var = sandy.CategoryCov.from_var([1, 1])
        >>> Vy = pd.DataFrame([[1, 0], [0, 1]], index=[1, 2], columns=[1, 2])
        >>> var._gls_cov_sensitivity(S, Vy)
                      1	          2
        3	4.00000e-01	4.00000e-01
        4	4.00000e-01	6.85714e-01
        """
        index, columns = pd.DataFrame(S).index, pd.DataFrame(S).columns
        S_ = pd.DataFrame(S).values
        general_sens = self._gls_general_sensitivity(S, Vy_extra, threshold=threshold).values
        cov_sens = general_sens.dot(S_)
        if threshold is not None:
            cov_sens[cov_sens < threshold] = 0
        return pd.DataFrame(cov_sens, index=index, columns=columns)

    def gls_update(self,  S, Vy_extra, threshold=None):
        """
        Perform GlS update for a given variance and sensitivity.

        Parameters
        ----------
        Vy_extra : 2D iterable
            2D covariance matrix for y_extra (MXM).
        S : 2D iterable
            Sensitivity matrix (MXN).
        delta : 1D iterable
            Perturbed vector minus non perturbed vector.
        threshold : `int`, optional
            Thereshold to avoid numerical fluctuations. The default is None.

        Returns
        -------
        `CategoryCov`
            GLS method apply to a CategoryCov object for a given Vy and S.

        Example
        -------
        >>> S = np.array([[1, 2], [3, 4]])
        >>> var = sandy.CategoryCov.from_var([1, 1])
        >>> Vy = np.diag(pd.Series([1, 1]))
        >>> var.gls_update(S, Vy)
                     0            1
        0  6.00000e-01 -4.00000e-01
        1 -4.00000e-01  3.14286e-01
        """
        index, columns = self.data.index, self.data.columns
        Vx_prior = self.data.values
        A = self._gls_cov_sensitivity(S, Vy_extra, threshold=threshold).values
        Vx_post = Vx_prior - A.dot(Vx_prior)
        return self.__class__(pd.DataFrame(Vx_post, index=index, columns=columns))

    def sandwich(self, S, threshold=None):
        """
        Apply the sandwich formula to the CategoryCov object for a given
        pandas.Series.

        Parameters
        ----------
        S : 1D or 2D iterable
            General sensitivities.

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
            C = self.data
            sandwich = corr2cov(C, S)
        else:
            S_ = pd.DataFrame(S).T
            sandwich = self._gls_Vy_calc(S_)
        if threshold is not None:
            sandwich[sandwich < threshold] = 0
        return sandwich

    def plot_corr(self, ax, **kwargs):
        add = {"cbar": True, "vmin": -1, "vmax": 1, "cmap": "RdBu"}
        for k, v in kwargs.items():
            add[k] = v
        ax = sns.heatmap(self.corr, ax=ax, **add)
        return ax

    @classmethod
    def corr2cov(cls, corr, std, **kwargs):
        """
        Produce covariance matrix given correlation matrix and standard
        deviation array.

        Parameters
        ----------
        corr : 2d `numpy.ndarray`
            square 2D correlation matrix

        std : 1d `numpy.ndarray`
            array of standard deviations

        Returns
        -------
        `sandy.CategoryCov`
            covariance matrix
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
        >>> sandy.CategoryCov.random_cov(2, seed=1)
                    0           1
        0 2.15373e-02 5.97134e-03
        1 5.97134e-03 8.52642e-03
        """
        corr = random_corr(size, correlations=correlations, seed=seed)
        std = np.random.uniform(stdmin, stdmax, size)
        return cls.corr2cov(corr, std, **kwargs)


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
                dtype=float,
                )
        self._data.columns = pd.Index(
                self._data.columns.values,
                name="E",
                dtype=float,
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
    corr : 2d iterable
        square 2D correlation matrix

    s : 1d iterable
        Diagonal matrix of standard deviations

    Returns
    -------
    `numpy.ndarray`
        covariance matrix

    Examples
    --------
    >>> S = np.array([1, 2, 3])
    >>> var = np.array([[1, 0, 2, 0], [0, 3, 0, 0], [4, 0, 5, 0], [0, 0, 0, 0]])
    >>> corr2cov(var, S)
        0   1   2  3
    0   1   0   6  0
    1   0  12   0  0
    2  12   0  45  0
    3   0   0   0  0

    >>> S = np.array([1, 2, 3])
    >>> var = np.array([[1, 0, 2], [0, 3, 0], [4, 0, 5]])
    >>> corr2cov(var, S)
        0   1   2
    0   1   0   6
    1   0  12   0
    2  12   0  45
    """
    s_ = pd.Series(s)
    corr_ = pd.DataFrame(corr)
    dim = corr_.values.shape[0]
    if len(s_) > dim:
        raise TypeError("The shape of the variables is not correct")
    # Create s diagonal matrix
    if len(s_) != dim:
        # Increase the size of s to the size of corr by filling it with zeros.
        dim_ext = dim - len(s)
        S = np.pad(np.diag(s_), ((0, dim_ext), (0, dim_ext)))
        S = pd.DataFrame(S, index=corr_.index, columns=corr_.columns)
    else:
        S = pd.DataFrame(np.diag(s_), index=s_.index, columns=s_.index)
        corr_ = pd.DataFrame(corr_.values, index=s_.index, columns=s_.index)
    return S.T.dot(corr_).dot(S)


def triu_matrix(arr, size):
    """
    Given the upper triangular values of a **square symmetric** matrix in
    an array, return the full matrix.

    Parameters
    ----------
    arr : 1d `numpy.ndarray`
        array with the upper triangular values of the matrix
    size : `int`
        dimension of the matrix

    Returns
    -------
    2d `numpy.ndarray`
        reconstructed 2d-array with symmetric matrix
    """
    matrix = np.zeros([size, size])
    indices = np.triu_indices(size)
    matrix[indices] = arr
    matrix += np.triu(matrix, 1).T
    return matrix


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
    corr = random_corr(size, correlations=correlations, seed=seed)
    std = np.random.uniform(stdmin, stdmax, size)
    return corr2cov(corr, std)


def random_ctg_cov(index, stdmin=0.0, stdmax=1.0, correlations=True, seed=None):
    cov = random_cov(len(index), stdmin=stdmin, stdmax=stdmax, correlations=correlations, seed=seed)
    return pd.DataFrame(cov, index=index, columns=index)


def get_eig(cov, sort=True):
    E, V = scipy.linalg.eig(cov)
    E = pd.Series(E.real, name="eigenvalues")
    V = pd.DataFrame(V.real)
    if sort:
        idx = E.sort_values(ascending=False).index
        E = E.iloc[idx].reset_index(drop=True)
        V = V.iloc[idx].reset_index(drop=True)
    return E, V


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