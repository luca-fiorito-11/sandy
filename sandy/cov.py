import numpy as np
import pandas as pd
from numpy import count_nonzero
from numpy.linalg import svd, matrix_rank
from numpy.linalg import norm as matrixnorm
from scipy.sparse import csr_matrix
from scipy.linalg import eig, qr, eigvals
from numpy.random import default_rng
import logging
import os
from scipy.stats import norm
from scipy.stats.qmc import LatinHypercube
from sandy.gls import sandwich


import sandy

pd.options.display.float_format = '{:.5e}'.format

__author__ = "Luca Fiorito"
__all__ = [
        "CategoryCov",
        "triu_matrix",
        "corr2cov",
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


def matrix_summary(A_values):
    """
    Summarizes key properties of a covariance matrix.
    """
    D = np.diag(A_values)
    
    A_nodiag = (A_values - np.diag(D))

    # calculate singular values
    U, s, V = svd(A_values, hermitian=True)
    s_max = np.max(s[s!=0])
    s_min = np.min(s[s!=0])
    svd_error = matrixnorm(A_values - U @ np.diag(s) @ U.T, ord="fro") / matrixnorm(A_values, ord="fro")

    # calculate eigenvalues
    e = eigvals(A_values)
    e_max = np.max(e[e!=0])
    e_min = np.min(e[e!=0])

    # calculate sparsity
    sparsity = 1.0 - ( count_nonzero(A_values) / float(A_values.size) )
    
    # norms
    fro_norm = matrixnorm(A_values, ord="fro")
    diag_norm = matrixnorm(np.diag(D), ord="fro")
    offdiag_norm = matrixnorm(A_nodiag, ord="fro")
    
    # Test symmetry with the Frobenius Norm
    symmetry_error = matrixnorm(A_values - A_values.T, ord='fro') / fro_norm

    std = np.sqrt(np.diag(A_values))
    std_zero = std[std > 0].size / std.size
    std_one = std[std > 1].size / std.size

    summary = {
        "Shape": A_values.shape,
        "Rank": matrix_rank(A_values),
        "Min Variance": np.min(D[D!=0]),
        "Max Variance": np.max(D[D!=0]),
        "Min Covariance": np.min(A_nodiag),
        "Max Covariance": np.max(A_nodiag),
        "Min Singular Values": s_min,
        "Max Singular Values": s_max,
        "Frobenius Norm": fro_norm,
        "Diagonal Norm": diag_norm,
        "Off-Diagonal Norm": offdiag_norm,
        "Condition number": s_max / s_min,
        "Min Eigenvalue": e_min,
        "Max Eigenvalue": e_max,
        "Sparsity": sparsity,
        "STD>0 Fraction": std_zero,
        "STD>1 Fraction": std_one,
        "SVD-approximation Error": svd_error,
        "Symmetry Error": symmetry_error,
    }
    return summary


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
    from_stdev
        construct a covariance matrix from a standard deviation vector
    from_var
        construct a covariance matrix from a variance vector
    get_corr
        extract correlation matrix from covariance matrix
    get_eig
        extract eigenvalues and eigenvectors from covariance matrix
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

        if not np.allclose(data.values, data.values.T):
            raise TypeError("Covariance matrix must be symmetric")

    @property
    def size(self):
        return self.data.values.shape[0]

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

    def summarize(self):
        """
        Summarizes key properties of the covariance matrix
        """
        A_values = self.data.values  # Convert to NumPy for calculations
        summary = matrix_summary(A_values)
        return summary

    def regularize(self, correction):
        """
        Regularizes the covariance matrix by adding a scaled diagonal matrix.
    
        This method adds a regularization term to the diagonal elements of the covariance matrix,
        improving numerical stability for further processing (e.g., inversion or sampling).
        Specifically, it computes a diagonal matrix where each diagonal element is scaled by
        the given `correction` factor and adds it to the covariance matrix.
    
        Parameters
        ----------
        correction: float
            A scalar multiplier for the diagonal elements to be added as a regularization term.
    
        Returns
        -------
        :obj: `~sandy.cov.CategoryCov`
            An instance of the same class with the regularized covariance matrix.
    
        Notes
        -----
        - This method assumes the covariance matrix is square.
        - The matrix is regularized in a symmetric way by design, since only the diagonal is modified.

        Example
        -------
        
        Regularize toy covariance matrix.
        
        >>> import pandas as pd
        >>> import numpy as np
        >>> import sandy
        >>> arrays = [[1, 1], [1, 2]]
        >>> index = pd.MultiIndex.from_arrays(arrays, names=("MT", "Other"))
        >>> cm = sandy.CategoryCov([[2.0, 0.5], [0.5, 3.0]], index=index, columns=index)
        >>> cm_reg = cm.regularize(0.1)
        >>> original = cm.data.values
        >>> regularized = cm_reg.data.values

        Test that the diagonal values are scaled correctly.
        
        >>> np.testing.assert_array_almost_equal(np.diag(regularized), np.diag(original) * 1.1)

        Test that the off-diagonal values didn't change.

        >>> offdiag_mask = ~np.eye(original.shape[0], dtype=bool)
        >>> np.testing.assert_array_equal(regularized[offdiag_mask], original[offdiag_mask])    
        """
        # don't need to pass via numpy. metadata are preserved
        C = self.data.copy()
        D = np.diag(C.values.diagonal() * correction)
        C += D
        return self.__class__(C)

    def correct_lognormal(self):
        """
        Corrects invalid covariance values in the data for lognormal sampling.
    
        In lognormal sampling, covariance matrix elements must satisfy the condition COV + 1 > 0.
        This method identifies values less than -1 and corrects them by setting them to (-1 + ε),
        where ε is the machine epsilon for float64. A warning is logged with information about how
        many invalid values were found, the MT numbers involved, and the smallest offending value.
    
        Returns
        -------
        :obj: `~sandy.cov.CategoryCov`
            An instance of the same class with corrected covariance matrix.
    
        Notes
        -----
        - Only the lower triangle (or symmetric) elements of the matrix are considered for counting.
        - The method assumes a symmetric covariance matrix indexed by a MultiIndex with level "MT".
        - If any invalid values are found (less than -1), a warning is logged indicating:
            - the number of offending values,
            - the smallest offending value,
            - the affected MT numbers (if "MT" is present in the index).    

        Example
        -------
        
        Simple test case.

        >>> import pandas as pd
        >>> import numpy as np
        >>> import sandy
        >>> arrays = [
        ...     [1, 2],
        ...     [1, 1]
        ... ]
        >>> index = pd.MultiIndex.from_arrays(arrays, names=("MT", "Other"))
        >>> cm = sandy.CategoryCov([[1, -1.2], [-1.2, 1]], index=index, columns=index)
        >>> cm_corrected = cm.correct_lognormal().data
        >>> assert np.all(cm_corrected.values >= -1)
        >>> assert cm_corrected.loc[(1, 1), (2, 1)] > -1

        """
        C = self.data.copy()

        # this condition limits covariances to max -100 %
        mask = C.values < -1

        if mask.any():
            size = ( mask.size - mask.diagonal().size ) // 2
            how_many_bad_values = mask.sum() // 2
            smallest_bad_value = C[mask].min().min()

            msg = f"""Condition COV + 1 > 0 for Lognormal sampling is not respected.
    {how_many_bad_values}/{size} covariance coefficients are set to -1+eps.
    The smallest covariance is {smallest_bad_value:.5f}
    """
            if "MT" in C.index.names:
                bad_mts = C.index[np.where(mask)[0]].get_level_values("MT").unique().tolist()
                msg += f"The concerned MT numbers are {bad_mts}."

            logging.warning(msg)

            C[mask] = -1 + np.finfo(np.float64).eps

        return self.__class__(C)

    def transform_lognormal(self):
        """
        Assuming that `self` constains the covariance matrix of a lognormal
        multivariate distribution centered in a unit vector, this method applies
        a transformation to calculate the covariance matrix of the underlying
        normal distribution.
        
        The function is taken from https://doi.org/10.1016/j.nima.2012.06.036

        Returns
        -------
        :obj: `~sandy.cov.CategoryCov`
            Covariance matrix of the underlying normal distribution.

        """
        # don't need to pass via numpy. metadata are preserved
        C = self.data.copy()
        C = np.log(C + 1)
        return self.__class__(C)

    def draw_sample(self, N, lhs=False, verbose=False, seed=None):
        M = self.size

        # -- Prepare index and columns for Samples object
        index = self.data.index
        columns = list(range(N))

        C = self.data.values
        D = C.diagonal()

        # -- Reduce matrix size by removing rows and columns with zero on diag
        nz = np.flatnonzero(D)
        Cr = C[nz][:, nz]

        # -- Decompose covariance (SVD better than QR or cholesky)
        Ur, Sr, Vr = svd(Cr, hermitian=True)  # hermitian is twice faster (U5 from JEFF33, 240 groups)

        Mr = Sr.size

        # -- Get U back to original size
        U = np.zeros((M, Mr))
        U[nz] = Ur

        # -- Draw IID samples with mu=0 and std=1
        seed_ = seed if seed else sandy.get_seed()
        if lhs:
            engine = LatinHypercube(d=Mr, seed=seed_)
            lhd = engine.random(n=N)
            X_ = norm(loc=0, scale=1).ppf(lhd).T
        else:
            rng = default_rng(seed=seed_)
            X_ = rng.standard_normal(size=(Mr, N))

        if verbose:
            summary = pd.Series(self.__class__(C).summarize()).to_string()
            print("======================================================")
            print(f"drawing sample of size N={N} from covariance matrix")
            print(f"seed: {seed_}")
            print("------------------------------------------------------")
            print("Covariance properties:")
            print(summary)
            print("======================================================")

        # -- Apply covariance to samples
        # -- 12 times faster with sparse (U5 from JEFF33, 240 groups)
        X = (csr_matrix(U) @ csr_matrix(np.diag(np.sqrt(Sr))) @ csr_matrix(X_)).todense()

        samples = pd.DataFrame(X, index=index, columns=columns)
        return sandy.Samples(samples)

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
        var = self.data.values.diagonal()
        std = np.sqrt(var)
        return pd.Series(std, index=self.data.index, name="STD")

    def get_eig(self, tolerance=None):
        """
        Compute the eigenvalues and eigenvectors of the dataset.
    
        This method extracts the eigenvalues and eigenvectors of `self.data`, 
        which is assumed to be a square symmetric matrix (e.g., a covariance 
        or correlation matrix). The eigenvalues are sorted in descending order, 
        and small values can optionally be set to zero using a tolerance threshold.
    
        Parameters
        ----------
        tolerance : float, optional (default: None)
            If specified, replaces eigenvalues smaller than a given fraction of 
            the largest eigenvalue with zero. The condition is:
    
            .. math::
                \\frac{e_i}{e_{MAX}} < \text{tolerance}
    
            - A value of `tolerance=1e-3` sets all eigenvalues 1000 times 
              smaller than the largest eigenvalue to zero.
            - A value of `tolerance=0` replaces all negative eigenvalues with zero.
            - If `None`, all eigenvalues are returned as computed.
    
        Returns
        -------
        pd.Series
            A Series of eigenvalues sorted in descending order, named "EIG".
        pd.DataFrame
            A DataFrame containing the corresponding eigenvectors, where each column 
            represents an eigenvector.
    
        Notes
        -----
        - Only the **real part** of the eigenvalues and eigenvectors is preserved.
        - The eigenvalues are not necessarily positive, especially for covariance matrices.
        - The implementation discussion is available [here](https://github.com/luca-fiorito-11/sandy/discussions/135).
    
        Examples
        --------
        Extract eigenvalues of a correlation matrix.

        >>> sandy.CategoryCov([[1, 0.4], [0.4, 1]]).get_eig()[0]
        0   1.40000e+00
        1   6.00000e-01
        Name: EIG, dtype: float64
    
        Extract eigenvectors.

        >>> sandy.CategoryCov([[1, 0.4], [0.4, 1]]).get_eig()[1]
                    0            1
        0  7.07107e-01 -7.07107e-01
        1  7.07107e-01  7.07107e-01
    
        Replace small eigenvalues using a tolerance.

        >>> sandy.CategoryCov([[0.1, 0.1], [0.1, 1]]).get_eig(tolerance=0.1)[0]
        0   0.00000e+00
        1   1.01098e+00
        Name: EIG, dtype: float64
    
        Handle negative eigenvalues.

        >>> sandy.CategoryCov([[1, 2], [2, 1]]).get_eig()[0]
        0    3.00000e+00
        1   -1.00000e+00
        Name: EIG, dtype: float64
    
        Replace negative eigenvalues with zero.

        >>> sandy.CategoryCov([[1, 2], [2, 1]]).get_eig(tolerance=0)[0]
        0   3.00000e+00
        1   0.00000e+00
        Name: EIG, dtype: float64
    
        Example with a covariance matrix.

        >>> sandy.CategoryCov([[1, 0.2, 0.1], [0.2, 2, 0], [0.1, 0, 3]]).get_eig()[0]
        0   9.56764e-01
        1   2.03815e+00
        2   3.00509e+00
        Name: EIG, dtype: float64
    
        Real test on H1 file.

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
    
        Ensure all eigenvalues are non-negative when using `tolerance=0`.

        >>> assert (cov.get_eig(tolerance=0)[0] >= 0).all()
        """
        E, V = eig(self.data)
        E = pd.Series(E.real, name="EIG")
        V = pd.DataFrame(V.real)
    
        if tolerance is not None:
            E[E / E.max() < tolerance] = 0
    
        return E, V

    def get_corr(self):
        """
        Extract correlation matrix.

        Returns
        -------
        df : :obj: `~sandy.cov.CategoryCov`
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

    def sampling(self, nsmp, seed=None, lognormal=True, correction=0.5/100,
                 lhs=False, verbose=False, **kwargs):
        """
        Extract perturbation coefficients from the covariance matrix using either
        a normal or lognormal distribution. Samples are adjusted to ensure physical
        plausibility (e.g., positivity, bounded range).
    
        Parameters
        ----------
        nsmp : int
            Number of samples to draw.
        seed : int, optional
            Seed for the random number generator (default is None).
        lognormal : bool, optional
            If True, use lognormal distribution for sampling. Otherwise, use (truncated) normal.
        correction : float, optional
            Regularization factor added to the diagonal of the covariance matrix to
            ensure positive definiteness (default is 0.5%).
        lhs : bool, optional
            If True, use Latin Hypercube Sampling (default is False).
        verbose : bool, optional
            If True, print progress information during sampling.
    
        Returns
        -------
        :obj:`~sandy.samples.Samples`
            An object containing the sampled perturbation coefficients.
    
        Notes
        -----
        - For normal sampling with relative perturbations, values below 0 or above 2
          are clipped. This truncation can degrade covariance accuracy for large
          uncertainties.
        - For lognormal sampling, values are always positive and the sample mean is
          guaranteed to converge to 1.
    
        Examples
        --------
        Common setup:
    
        >>> seed = 11
        >>> nsmp = 1e5
        >>> index = columns = ["A", "B"]
        >>> c = pd.DataFrame([[1, 0.4],[0.4, 1]], index=index, columns=columns) / 10
        >>> cov = sandy.CategoryCov(c)
    
        Normal sampling:
    
        >>> smp_n = cov.sampling(nsmp, seed=seed, lognormal=False)
        >>> np.testing.assert_array_almost_equal(smp_n.get_mean(), [1, 1], decimal=2)
        >>> np.testing.assert_array_almost_equal(smp_n.get_cov(), c, decimal=2)
    
        Lognormal sampling:
    
        >>> smp_ln = cov.sampling(nsmp, seed=seed, lognormal=True)
        >>> np.testing.assert_array_almost_equal(smp_ln.get_mean(), [1, 1], decimal=2)
        >>> np.testing.assert_array_almost_equal(smp_ln.get_cov(), c, decimal=2)
    
        Samples are reproducible:
    
        >>> assert cov.sampling(nsmp, seed=seed, lognormal=False).data.equals(smp_n.data)
    
        For large variances, normal sampling is truncated and does not reproduce the full covariance:
    
        >>> c = pd.DataFrame([[2, 0],[0, 2]])
        >>> s = sandy.CategoryCov(c).sampling(nsmp, lognormal=False)
        >>> np.testing.assert_array_almost_equal(s.get_mean(), [1, 1], decimal=2)
        >>> assert not np.allclose(s.get_cov(), c, atol=1)  # due to truncation
        >>> assert (s.get_rstd().values < 1).all()
        >>> assert np.linalg.norm(s.get_cov() - c) / np.linalg.norm(c) > 0.5
    
        For lognormal sampling, large variances are not an issue:
    
        >>> s = sandy.CategoryCov(c).sampling(nsmp, lognormal=True)
        >>> np.testing.assert_array_almost_equal(s.get_mean(), [1, 1], decimal=2)
        >>> np.testing.assert_allclose(
        ...     s.get_rstd().values, np.sqrt(np.diag(c)), rtol=0.2
        ... )
        >>> eigvals_original = np.linalg.eigvalsh(c)
        >>> eigvals_sampled = np.linalg.eigvalsh(s.get_cov())
        >>> np.testing.assert_allclose(eigvals_sampled, eigvals_original, rtol=0.2)
        """
        N = int(nsmp)

        if lognormal:
            C = self.correct_lognormal()

            var = C.data.values.diagonal()  # this will be used to adjust the mean
            # mean of the underlying normal distribution
            # https://stats.stackexchange.com/questions/573808/intuition-for-why-mean-of-lognormal-distribution-depends-on-variance-of-normally
            umu = np.log(1 / np.sqrt(var + 1))

            samples = (
                C.transform_lognormal()
                .regularize(correction=correction)
                .draw_sample(N, lhs=lhs, verbose=verbose, seed=seed)
                .apply_function(lambda x: x + umu.reshape(var.size, -1))
                .apply_function(np.exp)
                )

        else:
            samples = (
                self.regularize(correction=correction)
                .draw_sample(N, lhs=lhs, verbose=verbose, seed=seed)
                .apply_function(lambda x: x + 1)
                .truncate_normal()
                )

        return samples

    def sandwich(self, s):
        r"""
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

    def get_L(self, tolerance=None):
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
        Positive define matrix.

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

        >>> sandy.CategoryCov([[1, -2],[-2, 3]]).get_L(tolerance=0)
                       0	          1
        0	-1.08204e+00	0.00000e+00
        1	 1.75078e+00	0.00000e+00

        Decomposition test.

        >>> L = sandy.CategoryCov(a).get_L()
        >>> L.dot(L.T)
                       0	           1	           2
        0	 4.00000e+00	 1.20000e+01	-1.60000e+01
        1	 1.20000e+01	 3.70000e+01	-4.30000e+01
        2	-1.60000e+01	-4.30000e+01	 9.80000e+01

        Matrix with negative eigenvalues, tolerance of 0.

        >>> L = sandy.CategoryCov([[1, -2],[-2, 3]]).get_L(tolerance=0)
        >>> L.dot(L.T)
        	           0	           1
        0	 1.17082e+00	-1.89443e+00
        1	-1.89443e+00	 3.06525e+00
        """
        index = self.data.index
        columns = self.data.columns

        # Obtain the eigenvalues and eigenvectors
        E, V = sandy.CategoryCov(self.data).get_eig(tolerance=tolerance)

        # need sparse because much faster for large matrices (2kx2k from J33 Pu9)
        # with a lot of zero eigs
        # this is exactly equivalent to V.values @ np.diag(np.sqrt(E.values))
        A = (csr_matrix(V.values) @ csr_matrix(np.diag(np.sqrt(E.values)))).todense()
        
        # QR decomposition
        Q, R = qr(A.T)
        L = R.T

        return pd.DataFrame(L, index=index, columns=columns)

    def to_excel(self, file):
        """
        Save the sample dataset to an Excel file.

        This method exports the sample data to an Excel file, writing the dataset to 
        the 'COV' sheet. If the file already exists, the sheet is replaced; otherwise, 
        a new file is created.

        Parameters
        ----------
        file : str
            The path to the Excel file. If the file exists, the 'COV' sheet is replaced; 
            otherwise, a new file is created.

        Returns
        -------
        None
            This function does not return a value but writes the dataset to an Excel file.

        Notes
        -----
        - The dataset is saved in a sheet named 'COV' using the 'openpyxl' engine.
        - If the file exists, the function appends to it, replacing the 'COV' sheet.
        - The index of the DataFrame remains unchanged.
        """
        # Resetting indices messes everything up
        df = self.data
    
        # Determine write mode
        mode = "a" if os.path.exists(file) else "w"
        if_sheet_exists = "replace" if mode == "a" else None
    
        # Write to Excel
        with pd.ExcelWriter(file, mode=mode, engine="openpyxl", if_sheet_exists=if_sheet_exists) as writer:
            df.to_excel(writer, sheet_name="COV")

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
    >>> corr2cov(corr, s).astype(int)
    array([[ 1,  0,  6],
           [ 0, 12,  0],
           [ 6,  0,  9]])

    Test with float
    >>> corr2cov(corr, s.astype(float))
    array([[ 1.,  0.,  6.],
           [ 0., 12.,  0.],
           [ 6.,  0.,  9.]])
    """
    s_ = csr_matrix(np.diag(s))
    # sparse or else it is too slow (8000x8000), and anyways s_ is basically sparse
    return np.array((s_ @ csr_matrix(corr) @ s_).todense())


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
