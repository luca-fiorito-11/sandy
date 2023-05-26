import functools

import numpy as np
import scipy
import scipy.linalg
import scipy.sparse as sps
import pandas as pd
import logging
import tables as tb
import os
from sandy.gls import sandwich, _gls_cov_update

import sandy
pd.options.display.float_format = '{:.5e}'.format

__author__ = "Luca Fiorito"
__all__ = [
        "CategoryCov",
        "triu_matrix",
        "corr2cov",
        "random_corr",
        "random_cov",
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


class CategoryCov():
    """

    Properties
    ----------
    data
        covariance matrix as a dataframe
    size
        first dimension of the covariance matrix

    Methods
    -------
    corr2cov
        create a covariance matrix given a correlation matrix and a standard
        deviation vector
    from_stack
        create a covariance matrix from a stacked `pd.DataFrame`
    from_stdev
        construct a covariance matrix from a standard deviation vector
    from_var
        construct a covariance matrix from a variance vector
    get_corr
        extract correlation matrix from covariance matrix
    get_eig
        extract eigenvalues and eigenvectors from covariance matrix
    get_L
        extract lower triangular matrix such that $C=L L^T$
    get_std
        extract standard deviations from covariance matrix
    invert
        calculate the inverse of the matrix
    sampling
        extract perturbation coefficients according to chosen distribution
        and covariance matrix
    """

    def __repr__(self):
        return self.data.__repr__()

    def __init__(self, *args, **kwargs):
        self.data = pd.DataFrame(*args, dtype=float, **kwargs)

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
        >>> import pytest
        >>> with pytest.raises(TypeError): sandy.CategoryCov(np.array[1])
        >>> with pytest.raises(TypeError): sandy.CategoryCov(np.array([[1, 2], [2, -4]]))
        >>> with pytest.raises(TypeError): sandy.CategoryCov(np.array([[1, 2], [3, 4]]))
        """
        return self._data

    @data.setter
    def data(self, data):
        self._data = data

        if not len(data.shape) == 2 and data.shape[0] == data.shape[1]:
            raise TypeError("Covariance matrix must have two dimensions")

        if not (np.diag(data) >= 0).all():
            raise TypeError("Covariance matrix must have positive variance")

        sym_limit = 10
        # Round to avoid numerical fluctuations
        if not (data.values.round(sym_limit) == data.values.T.round(sym_limit)).all():
            raise TypeError("Covariance matrix must be symmetric")

    @property
    def size(self):
        return self.data.values.shape[0]

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
        return CategoryCov(corr).corr2cov(std)

    @classmethod
    def from_stdev(cls, std):
        """
        Construct the covariance matrix from the standard deviation vector.

        Parameters
        ----------
        var : 1D iterable
            Standar deviation vector.

        Returns
        -------
        `CategoryCov`
            Object containing the covariance matrix.

        Example
        -------
        Create covariance from stdev in `pd.Series`.
        >>> var = pd.Series(np.array([0, 2, 3]), index=pd.Index(["A", "B", "C"]))
        >>> std = np.sqrt(var)
        >>> cov = sandy.CategoryCov.from_stdev(std)
        >>> cov
                    A           B           C
        A 0.00000e+00 0.00000e+00 0.00000e+00
        B 0.00000e+00 2.00000e+00 0.00000e+00
        C 0.00000e+00 0.00000e+00 3.00000e+00
        """
        std_ = pd.Series(std)
        return cls.from_var(std_ * std_)

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
        Create covariance from variance in `pd.Series`.
        >>> var = pd.Series(np.array([0, 2, 3]), index=pd.Index(["A", "B", "C"]))
        >>> cov = sandy.CategoryCov.from_var(var)
        >>> cov
                    A           B           C
        A 0.00000e+00 0.00000e+00 0.00000e+00
        B 0.00000e+00 2.00000e+00 0.00000e+00
        C 0.00000e+00 0.00000e+00 3.00000e+00

        Create covariance from variance in list.
        >>> sandy.CategoryCov.from_var([1, 2, 3])
                    0           1           2
        0 1.00000e+00 0.00000e+00 0.00000e+00
        1 0.00000e+00 2.00000e+00 0.00000e+00
        2 0.00000e+00 0.00000e+00 3.00000e+00
        """
        var_ = pd.Series(var)
        values = np.diag(var_)
        cov = pd.DataFrame(values, index=var_.index, columns=var_.index)
        return cls(cov)

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

    def get_std(self):
        """
        Extract standard deviations.

        Returns
        -------
        `pandas.Series`
            1d array of standard deviations

        Examples
        --------
        >>> sandy.CategoryCov([[1, 0.4],[0.4, 1]]).get_std()
        0   1.00000e+00
        1   1.00000e+00
        Name: STD, dtype: float64
        """
        cov = self.to_sparse().diagonal()
        std = np.sqrt(cov)
        return pd.Series(std, index=self.data.index, name="STD")

    def get_eig(self, tolerance=None):
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
        >>> sandy.CategoryCov([[1, 0.4], [0.4, 1]]).get_eig()[0]
        0   1.40000e+00
        1   6.00000e-01
        Name: EIG, dtype: float64

        Extract eigenvectors of correlation matrix.
        >>> sandy.CategoryCov([[1, 0.4], [0.4, 1]]).get_eig()[1]
                    0            1
        0 7.07107e-01 -7.07107e-01
        1 7.07107e-01  7.07107e-01

        Extract eigenvalues of covariance matrix.
        >>> sandy.CategoryCov([[0.1, 0.1], [0.1, 1]]).get_eig()[0]
        0   8.90228e-02
        1   1.01098e+00
        Name: EIG, dtype: float64

        Set up a tolerance.
        >>> sandy.CategoryCov([[0.1, 0.1], [0.1, 1]]).get_eig(tolerance=0.1)[0]
        0   0.00000e+00
        1   1.01098e+00
        Name: EIG, dtype: float64

        Test with negative eigenvalues.
        >>> sandy.CategoryCov([[1, 2], [2, 1]]).get_eig()[0]
        0    3.00000e+00
        1   -1.00000e+00
        Name: EIG, dtype: float64

        Replace negative eigenvalues.
        >>> sandy.CategoryCov([[1, 2], [2, 1]]).get_eig(tolerance=0)[0]
        0   3.00000e+00
        1   0.00000e+00
        Name: EIG, dtype: float64

        Check output size.
        >>> cov = sandy.CategoryCov.random_cov(50, seed=11)
        >>> assert cov.get_eig()[0].size == cov.data.shape[0] == 50

        >>> sandy.CategoryCov([[1, 0.2, 0.1], [0.2, 2, 0], [0.1, 0, 3]]).get_eig()[0]
        0   9.56764e-01
        1   2.03815e+00
        2   3.00509e+00
        Name: EIG, dtype: float64

        Real test on H1 file
        >>> endf6 = sandy.get_endf6_file("jeff_33", "xs", 10010)
        >>> ek = sandy.energy_grids.CASMO12
        >>> err = endf6.get_errorr(errorr_kws=dict(ek=ek), err=1)["errorr33"]
        >>> cov = err.get_cov()
        >>> cov.get_eig()[0].sort_values(ascending=False).head(7)
        0    3.66411e-01
        1    7.05311e-03
        2    1.55346e-03
        3    1.60175e-04
        4    1.81374e-05
        5    1.81078e-06
        6    1.26691e-07
        Name: EIG, dtype: float64

        >>> assert not (cov.get_eig()[0] >= 0).all()

        >>> assert (cov.get_eig(tolerance=0)[0] >= 0).all()
        """
        E, V = scipy.linalg.eig(self.data)
        E = pd.Series(E.real, name="EIG")
        V = pd.DataFrame(V.real)
        if tolerance is not None:
            E[E/E.max() < tolerance] = 0
        return E, V

    def get_corr(self):
        """
        Extract correlation matrix.

        Returns
        -------
        df : :obj: `CetgoryCov`
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
            coeff = np.true_divide(1, self.get_std().values)
            coeff[~ np.isfinite(coeff)] = 0   # -inf inf NaN
        corr = np.multiply(np.multiply(cov, coeff).T, coeff)
        df =  pd.DataFrame(
            corr,
            index=self.data.index,
            columns=self.data.columns,
            )
        return self.__class__(df)

    def invert(self):
        """
        Method for calculating the inverse matrix.

        Returns
        -------
        `pandas.DataFrame`
            The inverse matrix.

        Notes
        -----
        Many covariance matrices for nuclear data are ill-defined and might
        have condition numbers that make the matrix inversion process
        impossible.
        To make up for this limitation we produce the (Moore-Penrose)
        pseudo-inverse of a Hermitian matrix, as implemented in `numpy`.
        Default options are used.
        This method does not require any pre-processing of the covariance
        data, e.g. removing zeros from the matrix diagonal or truncating
        eigenvalues.
        
        >>> c = sandy.CategoryCov(np.diag(np.array([1, 2, 3])))
        >>> c.invert()
                    0           1           2
        0 1.00000e+00 0.00000e+00 0.00000e+00
        1 0.00000e+00 5.00000e-01 0.00000e+00
        2 0.00000e+00 0.00000e+00 3.33333e-01

        Test product $A^T A = 1$.
        >>> c = sandy.CategoryCov([[2, 0.5], [0.5, 2]])
        >>> np.testing.assert_array_almost_equal(c.data @ c.invert(), np.eye(2))

        Test inverse of inverse.
        >>> c = sandy.CategoryCov(np.diag(np.array([1, 2, 3])))
        >>> np.testing.assert_array_equal(sandy.CategoryCov(c.invert()).invert(), c.data)
        
        Test on real ND covariance. With previous implementation this test failed.
        >>> c = sandy.get_endf6_file("jeff_33", "xs", 10010).get_errorr(err=1, errorr_kws=dict(mt=102))["errorr33"].get_cov()
        >>> cinv = c.invert()
        >>> a = c.data.values
        >>> b = cinv.values
        >>> np.testing.assert_array_almost_equal(a, a @ b @  a, decimal=4)
        >>> assert (cinv.index == c.data.index).all()
        >>> assert (cinv.columns == c.data.columns).all()
        """
        return pd.DataFrame(
            np.linalg.pinv(self.data, hermitian=True),
            index=self.data.index,
            columns=self.data.columns,
        )  # do not return CategoryCov because variance can be negative

    def sampling(self, nsmp, seed=None, pdf='normal', tolerance=0, relative=True, **kwargs):
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

        .. note:: sampling with relative covariance matrix is performed
            setting all the negative perturbation coefficients equal to 0
            and the ones larger than 2 equal to 2 for normal distribution, or
            adjusting the standard deviations in such a way that negative
            samples are avoided for uniform distribution.

        .. note:: sampling with lognormal distribution gives a set of samples
            with mean=1 since a lognormal distribution cannot have mean=0.
            In this case the `relative` parameter does not apply to it.

        Examples
        --------
        Common parameters.
        >>> seed = 11
        >>> nsmp = 1e5
        
        Create Positive-Definite covariance matrix with small stdev.
        >>> index = columns = ["A", "B"]
        >>> c = pd.DataFrame([[1, 0.4],[0.4, 1]], index=index, columns=index) / 10
        >>> cov = sandy.CategoryCov(c)

        Draw relative samples using different distributions.
        >>> smp_n = cov.sampling(nsmp, seed=seed, pdf='normal')
        >>> smp_ln = cov.sampling(nsmp, seed=seed, pdf='lognormal')
        >>> smp_u = cov.sampling(nsmp, seed=seed, pdf='uniform')

        The sample mean converges to a unit vector.
        >>> np.testing.assert_array_almost_equal(smp_n.get_mean(), [1, 1], decimal=2)
        >>> np.testing.assert_array_almost_equal(smp_ln.get_mean(), [1, 1], decimal=2)
        >>> np.testing.assert_array_almost_equal(smp_u.get_mean(), [1, 1], decimal=2)

        The sample covariance converges to the original one.
        >>> np.testing.assert_array_almost_equal(smp_n.get_cov(), c, decimal=3)
        >>> np.testing.assert_array_almost_equal(smp_ln.get_cov(), c, decimal=3)
        >>> np.testing.assert_array_almost_equal(np.diag(smp_u.get_cov()), np.diag(c), decimal=3)

        Samples are reproducible by setting a seed.
        assert cov.sampling(nsmp, seed=seed, pdf='normal').data.equals(smp_n.data)



        Create Positive-Definite covariance matrix with small stdev (small negative eig).
        >>> c = pd.DataFrame([[1, 1.2],[1.2, 1]], index=index, columns=index) / 10
        >>> cov = sandy.CategoryCov(c)

        >>> smp_0 = cov.sampling(nsmp, seed=seed, pdf='normal', tolerance=0)
        >>> np.testing.assert_array_almost_equal(np.diag(smp_0.get_cov()), np.diag(c), decimal=2)
        >>> smp_inf = cov.sampling(nsmp, seed=seed, pdf='normal', tolerance=np.inf)
        >>> import pytest
        >>> with pytest.raises(Exception):
        ...    raise np.testing.assert_array_almost_equal(np.diag(smp_0.get_cov()), np.diag(c), decimal=1)



        Create Positive-Definite covariance matrix with small stdev (large negative eig).
        >>> c = pd.DataFrame([[1, 4],[4, 1]], index=index, columns=index) / 10
        >>> cov = sandy.CategoryCov(c)
        
        Samples kind of converge only if we set a low tolerance
        >>> smp_0 = cov.sampling(nsmp, seed=seed, pdf='normal', tolerance=0)
        >>> with pytest.raises(Exception):
        ...    raise np.testing.assert_array_almost_equal(np.diag(smp_0.get_cov()), np.diag(c), decimal=1)
        >>> smp_inf = cov.sampling(nsmp, seed=seed, pdf='normal', tolerance=np.inf)
        >>> with pytest.raises(Exception):
        ...    raise np.testing.assert_array_almost_equal(np.diag(smp_0.get_cov()), np.diag(c), decimal=1)



        Create Positive-Definite covariance matrix with large stdev.
        >>> index = columns = ["A", "B"]
        >>> c = pd.DataFrame([[1, 0.4],[0.4, 1]], index=index, columns=index) / 10
        >>> cov = sandy.CategoryCov(c)

        Need to increase the amount of samples
        >>> nsmp = 1e6

        The sample mean still converges to a unit vector.
        >>> np.testing.assert_array_almost_equal(smp_n.get_mean(), [1, 1], decimal=2)
        >>> np.testing.assert_array_almost_equal(smp_ln.get_mean(), [1, 1], decimal=2)
        >>> np.testing.assert_array_almost_equal(smp_u.get_mean(), [1, 1], decimal=2)

        Only the lognormal covariance still converges.
        >>> with pytest.raises(Exception):   
        ...    raise np.testing.assert_array_almost_equal(smp_n.get_cov(), c, decimal=1)
        >>> np.testing.assert_array_almost_equal(smp_ln.get_cov(), c, decimal=2)
        >>> with pytest.raises(Exception):   
        ...    raise np.testing.assert_array_almost_equal(np.diag(smp_u.get_cov()), np.diag(c), decimal=1)
        """
        allowed_pdf = [
            "normal",
            "lognormal",
            "uniform",
            ]
        if pdf not in allowed_pdf:
            raise ValueError("`pdf='lognormal'` not allowed")

        if not relative and pdf=='lognormal':
            raise ValueError("`pdf='lognormal'` and `relative=False` is not a valid combination")
        
        nsmp_ = int(nsmp)

        # -- Draw IID samples with mu=0 and std=1
        np.random.seed(seed=seed)
        if pdf == 'uniform':
            a = np.sqrt(12) / 2
            y = np.random.uniform(-a, a, (self.size, nsmp_))
        else:
            y = np.random.randn(self.size, nsmp_)

        # -- Fix covariance matrix according to distribution
        if pdf == 'uniform':
            # no cross-correlation term is considered
            if relative:
                a = np.sqrt(12) / 2 # upper bound of the distribution y
                std = np.sqrt(np.diag(self.data)) 
                std_modified = np.where(std < 1 / a, std, 1 / a)
                cov = np.diag(std_modified**2)
            else:
                cov = np.diag(np.diag(self.data))
            to_decompose = self.__class__(cov, index=self.data.index, columns=self.data.columns)
        elif pdf == 'lognormal':
            ucov = np.log(self.sandwich(np.eye(self.size)).data + 1).values  # covariance matrix of underlying normal distribution
            to_decompose = self.__class__(ucov, index=self.data.index, columns=self.data.columns)
        else:
            to_decompose = self

        # -- Decompose covariance into lower triangular
        L = to_decompose.get_L(tolerance=tolerance)

        # -- Apply covariance to samples
        sparse = False  # it seems to me that a sparse inner product is only marginally faster
        if sparse:
            y = sps.csc_matrix(y)
            L = scipy.sparse.csr_matrix(L)
            inner = (L @ y).toarray()
        else:
            inner = L @ y

        index = self.data.index
        columns = list(range(nsmp_))
        samples = pd.DataFrame(inner, index=index, columns=columns)

        # -- Fix sample (and sample mean) according to distribution
        if pdf == 'lognormal':
            mu = np.ones(self.size)
            umu = np.log(mu**2 / np.sqrt(np.diag(self.data) + mu**2))   # mean of the underlying normal distribution
            samples = np.exp(samples.add(umu, axis=0))
        elif relative:
            samples += 1
            lower_bound = samples > 0
            upper_bound = samples < 2
            samples = samples.where(lower_bound, 0)
            samples = samples.where(upper_bound, 2)

        return sandy.Samples(samples)

    def gls_cov_update(self, S, Vy_extra=None):
        """
        Perform GlS update for a given covariance matrix, sensitivity and
        covariance matrix of the extra information:
        .. math::
            $$
            V_{x_{post}} = V_{x_{prior}} - V_{x_{prior}}\cdot S^T\cdot \left(S\cdot V_{x_{prior}}\cdot S^T + V_{y_{extra}}\right)^{-1}\cdot S\cdot V_{x_{prior}}
            $$

        Parameters
        ----------
        Vy_extra : 2D iterable or sigle element 1D iterable
            covariance matrix of the extra information,
            (M, M) or (1,).
        S : 2D or 1D iterable
            Sensitivity matrix (M, N) or sensitivity vector (N,).

        Returns
        -------
        `sandy.CategoryCov`
            `CategoryCov` object corresponding to the updated covariance matrix
            adjusted with the GLS technique.

        Notes
        -----
        .. note:: If Vy_extra=None the constraint GLS update technique
        will be performed

        Example
        -------
        >>> S = np.array([[1, 2], [3, 4]])
        >>> cov = sandy.CategoryCov.from_var([1, 1])
        >>> Vy = np.diag(pd.Series([1, 1]))
        >>> cov.gls_cov_update(S, Vy)
                     0            1
        0  6.00000e-01 -4.00000e-01
        1 -4.00000e-01  3.14286e-01

        >>> cov = pd.DataFrame([[1, 0], [0, 1]], index=[2, 3], columns=[2, 3])
        >>> cov = sandy.CategoryCov(cov)
        >>> cov.gls_cov_update(S, Vy)
                    2            3
        2  6.00000e-01 -4.00000e-01
        3 -4.00000e-01  3.14286e-01

        >>> S = pd.DataFrame([[1, 0], [0, 1], [0, 1]], index=[2, 3, 4], columns=[1, 2])
        >>> cov = pd.DataFrame([[1, 0], [0, 1]], index=[2, 3], columns=[2, 3])
        >>> cov = sandy.CategoryCov(cov)
        >>> Vy = np.diag(pd.Series([1, 1, 1]))
        >>> cov.gls_cov_update(S, Vy)
                    2           3
        2 5.00000e-01 0.00000e+00
        3 0.00000e+00 3.33333e-01

        >>> S = np.array([1, 2])
        >>> cov = sandy.CategoryCov.from_var([1, 1])
        >>> Vy = [1]
        >>> cov.gls_cov_update(S, Vy)
                     0            1
        0  8.33333e-01 -3.33333e-01
        1 -3.33333e-01  3.33333e-01
        """
        idx = self.data.index
        S_ = np.array(S)
        if Vy_extra is not None:
            Vy_extra_ = pd.DataFrame(Vy_extra).values
            Vx_post = _gls_cov_update(self.data.values, S_, Vy_extra=Vy_extra_)
        else:
            Vx_post = _gls_cov_update(self.data.values, S_)
        Vx_post = pd.DataFrame(Vx_post, index=idx, columns=idx)
        return self.__class__(Vx_post)

    def sandwich(self, s):
        """
        Apply the "sandwich formula" to the CategoryCov object for a given
        sensitivity. According with http://dx.doi.org/10.1016/j.anucene.2015.10.027,
        the moment propagation equation is implemented as:

           .. math::
               $$
               V_R = S\cdot V_P\cdot S^T
               $$

        Parameters
        ----------
        s : 1D or 2D iterable
            General sensitivities (N,) or (M, N) with N the size of the
            `CategoryCov` object.

        Returns
        -------
        `sandy.CategoryCov`
            `CategoryCov` object corresponding to the response covariance matrix
            obtained with the sandwich formula.

        Examples
        --------
        >>> var = np.array([1, 2, 3])
        >>> s = np.array([[1, 2, 3]])
        >>> assert s.shape == (1, 3)
        >>> cov = sandy.CategoryCov.from_var(var)
        >>> cov.sandwich(s)
                    0
        0 3.60000e+01

        >>> s = np.array([1, 2, 3])
        >>> var = pd.Series([1, 2, 3])
        >>> cov = sandy.CategoryCov.from_var(var)
        >>> sensitivity = np.diag(s)
        >>> cov.sandwich(sensitivity)
                    0           1           2
        0 1.00000e+00 0.00000e+00 0.00000e+00
        1 0.00000e+00 8.00000e+00 0.00000e+00
        2 0.00000e+00 0.00000e+00 2.70000e+01

        >>> s = pd.DataFrame([[1, 0, 1], [0, 1, 1]], index=[2, 3], columns=[2, 3, 4]).T
        >>> cov = pd.DataFrame([[1, 0], [0, 1]], index=[2, 3], columns=[2, 3])
        >>> cov = sandy.CategoryCov(cov)
        >>> cov.sandwich(s)
                    2           3           4
        2 1.00000e+00 0.00000e+00 1.00000e+00
        3 0.00000e+00 1.00000e+00 1.00000e+00
        4 1.00000e+00 1.00000e+00 2.00000e+00
        """
        s_ = pd.DataFrame(s)
        index = s_.index
        sandwich_ = sandwich(self.data.values, s_.values)
        if len(sandwich_.shape) == 0: 
            sandwich_ = [sandwich_]
        sandwich_ = pd.DataFrame(sandwich_, index=index, columns=index)
        return self.__class__(sandwich_)

    def corr2cov(self, std):
        """
        Produce covariance matrix given correlation matrix and standard
        deviation array.
        Same as :obj: `corr2cov` but it works with :obj: `CategoryCov`
        instances.

        Parameters
        ----------
        corr : :obj: `CategoryCov`
            square 2D correlation matrix
        std : 1d iterable
            array of standard deviations

        Returns
        -------
        :obj: `CategoryCov`
            covariance matrix

        Examples
        --------
        Initialize index and columns
        >>> idx = ["A", "B", "C"]
        >>> std = np.array([1, 2, 3])
        >>> corr = sandy.CategoryCov([[1, 0, 2], [0, 3, 0], [2, 0, 1]], index=idx, columns=idx)
        >>> corr.corr2cov(std)
                    A           B           C
        A 1.00000e+00 0.00000e+00 6.00000e+00
        B 0.00000e+00 1.20000e+01 0.00000e+00
        C 6.00000e+00 0.00000e+00 9.00000e+00
        """
        cov = corr2cov(self.data, std)
        index = self.data.index
        columns = self.data.columns
        return self.__class__(cov, index=index, columns=columns)

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
        E, V = sandy.CategoryCov(cov_reduced).get_eig(tolerance=tolerance)
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
