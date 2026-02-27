"""
Covariance and correlation utilities for nuclear data analysis.

This module provides tools for constructing, transforming, validating,
regularizing, sampling, and manipulating covariance matrices represented
as pandas DataFrames. It is designed around the `CategoryCov` class, which
encapsulates a covariance matrix together with its index/column metadata,
and exposes numerical diagnostics and domain‑specific operations relevant
to nuclear data uncertainty propagation.

Main Features
-------------
- Construction of covariance matrices from variances, standard deviations,
  and correlation matrices.
- Comprehensive covariance diagnostics, including:
    * eigenvalues, singular values, PSD checks
    * negative-eigenvalue mass, condition numbers
    * correlation analysis and sparsity metrics
    * recommended regularization factors
- Regularization and PSD-repair via diagonal loading.
- Conversion between X‑space (lognormal) and Z‑space (normal) covariance
  representations, assuming a unit‑mean lognormal model.
- Covariance‑based sampling:
    * normal sampling with truncation for physical bounds
    * lognormal sampling with proper mean-shifting and exponentiation
    * optional Latin Hypercube Sampling (LHS)
- Sandwich covariance propagation for sensitivity analysis.

All matrix computations use NumPy/SciPy routines with careful numerical
treatment for symmetry, conditioning, and floating‑point stability.

Intended Use
------------
The module is primarily designed for covariance matrices in nuclear data
evaluation workflows, but it is general enough to be used for any structured
covariance matrix with pandas indexing and metadata.

Dependencies
------------
- numpy
- pandas
- scipy (for eigh, LHS, etc.)
- openpyxl (optional, for Excel export)

"""

import numpy as np
import pandas as pd
import logging
import os

from numpy.linalg import norm as matrixnorm
from numpy.linalg import qr



__author__ = "Luca Fiorito"
__all__ = [
        "CategoryCov",
        "corr2cov",
        ]



class CategoryCov():
    """
    Representation of a covariance matrix with rich metadata and numerical tools.

    `CategoryCov` wraps a square covariance matrix stored as a pandas DataFrame,
    preserving index/column labels (including MultiIndex structures). The class
    provides a coherent interface for:

    • validating covariance structure  
    • inspecting index/column composition  
    • performing spectral diagnostics  
    • regularizing and repairing covariance matrices  
    • transforming between normal and lognormal representations  
    • drawing random samples consistent with the covariance  
    • uncertainty propagation through the sandwich formula  

    The goal is to ensure that covariance matrices used in nuclear data
    simulations are numerically sound, physically meaningful, and easy to
    manipulate while retaining metadata such as MT/E group structure.

    Parameters
    ----------
    *args, **kwargs :
        Passed directly to `pandas.DataFrame` to build the underlying matrix.
        Special keyword:
        - covariance_checks : bool, default True
            If True, enforce basic covariance validity:
            * matrix must be square
            * diagonal variances must be non‑negative
            * matrix must be symmetric within numerical tolerance

    Attributes
    ----------
    data : pandas.DataFrame
        The covariance matrix with index and columns preserved.
    size : int
        Dimension of the covariance matrix (number of energy groups/categories).

    Key Methods
    -----------
    summarize() :
        Compute extensive diagnostics on the covariance matrix, including
        eigenvalues, singular values, PSD checks, sparsity, correlation
        statistics, lognormal constraints, and suggested regularization.
    summarize_index() :
        Inspect structure of index/columns, including per‑level summaries
        for MultiIndex objects.

    transform_lognormal() :
        Convert a *relative* X‑space (lognormal) covariance to the underlying
        Z‑space (normal) covariance using Σ_Z = log1p(C), assuming unit‑mean
        lognormal variables.

    correct_lognormal() :
        Ensure the domain condition C_ij > -1 for lognormal sampling is met,
        replacing invalid entries by (-1 + eps) and issuing a warning.

    regularize(correction) :
        Apply diagonal loading (D = correction * diag(A)) to improve
        positive‑semidefiniteness and numerical stability.

    sampling(nsmp, lognormal=True, ...) :
        Generate samples consistent with the covariance matrix using either a
        normal or lognormal distribution. Supports Latin Hypercube Sampling.

    get_std(), get_corr(), get_eig(), get_L() :
        Extract standard deviations, correlation matrix, eigenvalues/eigenvectors,
        and the lower‑triangular Cholesky‑like factor L.

    sandwich(s) :
        Apply the sensitivity‑based covariance propagation formula:
            V_R = S · V_P · Sᵀ

    corr2cov(std) :
        Construct a covariance matrix from a correlation matrix and a vector
        of standard deviations.

    Notes
    -----
    - All computations assume double‑precision floats.
    - Symmetry and PSD properties are treated carefully with numerical
      tolerances appropriate for floating‑point covariance matrices.
    - The class is designed to allow large matrices (up to several thousand
      dimensions), using sparse operations where appropriate.

    Examples
    --------
    Construct a covariance matrix:
    >>> cov = CategoryCov([[1, 0.4], [0.4, 1]], index=["A", "B"], columns=["A", "B"])

    Compute diagnostics:
    >>> summary = cov.summarize()

    Convert from lognormal to normal covariance (unit-mean model):
    >>> cov_Z = cov.transform_lognormal()

    Draw samples:
    >>> smp = cov.sampling(10000, lognormal=True)
    """

    def __repr__(self):
        with pd.option_context("display.float_format", "{:.5e}".format):
            return self.data.__repr__()

    def __init__(self, *args, **kwargs):
        self._covariance_checks = kwargs.pop("covariance_checks", True)     # store the flag so data.setter can access it
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
        
        Test incorrect covariance matrices.

        >>> import sandy, pytest
        >>> with pytest.raises(TypeError): sandy.CategoryCov(np.array([[1, 2], [2, -4]]))
        >>> with pytest.raises(TypeError): sandy.CategoryCov(np.array([[1, 2], [3, 4]]))

        By-pass checks using `covariance_checks=False`.

        >>> c = sandy.CategoryCov(np.array([[1, 2], [2, -4]]), covariance_checks=False)
        >>> c = sandy.CategoryCov(np.array([[1, 2], [3, 4]]), covariance_checks=False)

        """
        return self._data

    @data.setter
    def data(self, data):
        self._data = data
        
        if self._covariance_checks:

            # Shape check
            if self._data.ndim != 2 or self._data.shape[0] != self._data.shape[1]:
                raise TypeError("Covariance matrix must be square (2D, n x n).")

            # Positive variances
            diag = self._data.values.diagonal()
            if not np.all(diag >= 0):
                raise TypeError("Covariance matrix must have non-negative variances on the diagonal.")

            # Symmetry (tolerance)
            if not np.allclose(self._data.values, self._data.values.T, rtol=1e-5, atol=1e-8):
                raise TypeError("Covariance matrix must be symmetric within numerical tolerance.")

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

        >>> import sandy
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

        >>> import sandy
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

    def summarize_index(self):
        # =================================
        # CHECK INDEX/COLUMNS
        # =================================
        idx = self.data.index
        cols = self.data.columns
        
        # Basic structural checks
        index_equals_columns = idx.equals(cols)
        
        # Determine index type
        is_multi = isinstance(idx, pd.MultiIndex)
        n_levels = idx.nlevels if is_multi else 1
        level_names = list(idx.names) if is_multi else [idx.name]
        
        # Per-level summaries
        level_summaries = {}
        if is_multi:
            for i in range(idx.nlevels):
                lvl = idx.get_level_values(i)
                name = idx.names[i]
                unique = lvl.unique().values.tolist()
                level_summary = {
                    "n_unique": len(unique),
                }
                if name != "E":
                    level_summary["values"] = unique
                level_summaries[name] = level_summary

        else:
            name = idx.name
            unique = idx.unique().values.tolist()
            level_summary = {
                "n_unique": len(unique),
            }
            if name != "E":
                level_summary["values"] = unique
            level_summaries[name] = level_summary
        
        summary_index_columns = {
            "Index Equals Columns": index_equals_columns,
            "Index Type": "MultiIndex" if is_multi else "Index",
            "Levels": n_levels,
            "Level Names": level_names,
            "Length": len(idx),
            "Index Is Unique": idx.is_unique,
            "Level Summaries": level_summaries,
        }
        return summary_index_columns


    def summarize(self):
        """
        Summarizes key properties of the covariance matrix

        Returns a dict with diagnostics.
        """
        from scipy.linalg import eigh

        summary_index_columns = self.summarize_index()

        # =================================
        # CHECK VALUES
        # =================================
        A_values = self.data.values  # Convert to NumPy for calculations

        A = np.asarray(A_values)
        n = A.shape[0]
    
        # -----------------------------
        # Basic structural diagnostics
        # -----------------------------
        rank = np.linalg.matrix_rank(A)
        diag_vec = np.diag(A)
        A_diag = np.diag(diag_vec)
        A_nodiag = A - A_diag
        rank_deficiency = n - rank
    
        # -----------------------------
        # Norm diagnostics
        # -----------------------------
        fro_norm = float(matrixnorm(A, ord="fro"))
        diag_fro_norm = float(matrixnorm(A_diag, ord="fro"))
        offdiag_fro_norm = float(matrixnorm(A_nodiag, ord="fro"))
    
        # Symmetry diagnostics
        sym_residual = A - A.T
        symmetry_error = matrixnorm(sym_residual, ord='fro') / fro_norm
    
        # -----------------------------
        # Eigendecomposition (correct method for cov matrices)
        # -----------------------------
        e, U = eigh(A)
        eig_nonzero = e[e != 0]
        eig_negative = e[e < 0]
        e_max = float(np.max(eig_nonzero))
        e_min = float(np.min(eig_nonzero))
    
        # PSD check with tolerance
        tol = 1e-10
        is_psd = bool(e_min >= -tol)
    
        # Negative eigenvalue severity (important!)
        negative_eig_mass = float(-np.sum(eig_negative))  # clipped total negative mass
        neg_mass_ratio_fro = negative_eig_mass / fro_norm
        total_var = np.sum(diag_vec)
        neg_mass_ratio_total_var = negative_eig_mass / total_var
    
        # -----------------------------
        # SVD diagnostic
        # -----------------------------
        U, s, V = np.linalg.svd(A, hermitian=True)
        s_nonzero = s[s!=0]
        s_max = np.max(s_nonzero)
        s_min = np.min(s_nonzero)
    
        # Reconstruction error using eigen-decomposition (since A_sym is symmetric)
        # A ≈ U diag(s) U^T
        A_recon = U @ np.diag(s) @ U.T
        svd_approx_error = matrixnorm(A - A_recon, ord="fro") / fro_norm
    
        # Condition numbers
        tol = 1e-12
        cond = float(s_max / s_min)
        log_condition = float(np.log10(cond)) if cond not in (0, np.inf) else np.inf
        is_nearly_singular = bool(s_min < tol)
    
        # -----------------------------
        # Sparsity, calculate in relative terms (between 0 and 1)
        # -----------------------------
        sparsity = float(1.0 - (np.count_nonzero(A_values) / float(A.size)))
        
        # -----------------------------
        # Variances and off-diagonals
        # -----------------------------
        diag_nonzero = diag_vec[diag_vec != 0]
        min_var = float(np.min(diag_nonzero))
        max_var = float(np.max(diag_nonzero))
    
        # Off-diagonal min/max (consider entire matrix excluding diag)
        off_min = float(np.min(A_nodiag))
        off_max = float(np.max(A_nodiag))
    
        # -----------------------------
        # Std diagnostics
        # -----------------------------
        std = np.sqrt(diag_vec)
        std_pos_frac = float((std > 0).sum() / std.size)
        std_gt1_frac = float((std > 1).sum() / std.size)
    
        # -----------------------------
        # Correlation diagnostics
        # -----------------------------
        Dinv = np.diag(1.0 / np.sqrt(np.clip(diag_vec, 1e-15, None)))
        Corr = Dinv @ A @ Dinv
        I = np.eye(n)    
        Corr_nodiag = Corr - I
        Corr_nodiag_abs = np.abs(Corr_nodiag)
        max_corr = float(np.max(Corr_nodiag))
        min_corr = float(np.min(Corr_nodiag))
        max_abs_corr = float(np.max(Corr_nodiag_abs))
    
        # Percentiles tell you the typical and extreme correlation strength in the matrix
        corr_p90 = float(np.percentile(Corr_nodiag_abs, 90))  # 90% of abs(correlations) below this value
        corr_p95 = float(np.percentile(Corr_nodiag_abs, 95))  # 95% of abs(correlations) below this value
        corr_p99 = float(np.percentile(Corr_nodiag_abs, 99))  # 99% of abs(correlations) below this value


        # -----------------------------
        # Recommended diagonal regularization (lambda)
        # -----------------------------
        # 1. PSD repair lambda
        if e_min < 0:
            lambda_psd = -e_min / min_var
        else:
            lambda_psd = 0.0
        
        # 2. negative eigenvalue mass lambda
        if negative_eig_mass > 0:
            lambda_neg = (negative_eig_mass / n) / min_var
        else:
            lambda_neg = 0.0
        
        # final recommended lambda
        lambda_recommended = max(lambda_psd, lambda_neg, 0.0)

        # -----------------------------
        # Possible Lognormal transformation
        # -----------------------------
        mask = A < -1.0
        # this condition limits covariances to max -100 %
        if mask.any():
            # Count off-diagonal unique pairs only (upper triangle)
            iu = np.triu_indices(n, k=1)
            how_many_bad_values = mask[iu].sum()
            smallest_bad_value = A[mask].min().min()
            largest_bad_value = A[mask].max().max()
        else:
            how_many_bad_values = 0
            smallest_bad_value = None
            largest_bad_value = None

        # -----------------------------
        # Final summary dictionary
        # -----------------------------
        summary = {
            "Shape": A.shape,
            "Rank": rank,
            "Rank Deficiency": rank_deficiency,
    
            "Min Variance": min_var,
            "Max Variance": max_var,
            "Min Covariance (off-diag)": off_min,
            "Max Covariance (off-diag)": off_max,
    
            "Min Singular Value": s_min,
            "Max Singular Value": s_max,
    
            "Min Eigenvalue": e_min,
            "Max Eigenvalue": e_max,
            "Negative Eigenvalue Mass": negative_eig_mass,
            "Negative Eigenvalue Mass (% Fro Norm)": 100 * neg_mass_ratio_fro,
            "Negative Eigenvalue Mass (% Tot Var)": 100 * neg_mass_ratio_total_var,
            "Is PSD (tol=1e-10)": is_psd,
    
            "Condition Number": cond,
            "Log10 Condition Number": log_condition,
            "Is Nearly Singular (tol=1e-12)": is_nearly_singular,
    
            "Frobenius Norm": fro_norm,
            "Diagonal Norm": diag_fro_norm,
            "Off-Diagonal Norm": offdiag_fro_norm,
    
            "Sparsity (% of zeros)": 100 * sparsity,
            "STD>0 (%)": 100 * std_pos_frac,
            "STD>1 (%)": 100 * std_gt1_frac,
    
            "Max Correlation (off-diag)": max_corr,
            "Min Correlation (off-diag)": min_corr,
            "Max |Correlation| (off-diag)": max_abs_corr,
            "Corr p90": corr_p90,
            "Corr p95": corr_p95,
            "Corr p99": corr_p99,
    
            "SVD Approximation Error (Fro, %)": 100 * svd_approx_error,
            "Symmetry Error (Fro, %)": 100 * symmetry_error,

            "Lognormal invalid pairs (Cov_ij < -1)": how_many_bad_values,
            "Smallest Cov not respecting LogN": smallest_bad_value,
            "Largest Cov not respecting LogN": largest_bad_value,

            "Regularization λ (Recommended)": lambda_recommended,
        }
        
        # Dicts print with np.float64(...), this will remove np.float64()
        summary = {k: (v.item() if isinstance(v, np.generic) else v)
               for k, v in summary.items()}

        summary["Covariance Index"] = summary_index_columns

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
        
        Simple test case. First create MultiIndex.

        >>> import sandy
        >>> index_arrays = [[1, 2], [1, 1]]
        >>> index = pd.MultiIndex.from_arrays(index_arrays, names=("MT", "Other"))

        Then, create covariance matrix that violates log1p domain (entries <= -1).

        >>> a = [[1, -1.2], [-1.2, 1]]
        >>> cm = sandy.CategoryCov(a, index=index, columns=index)

        Correct values that will make `transform_lognormal` fail.

        >>> cm_corrected = cm.correct_lognormal().data

        Check that all covariances respect condition.

        >>> arr = cm_corrected.to_numpy()
        >>> assert np.all(1.0 + arr > 0.0)
        >>> assert cm_corrected.loc[(1, 1), (2, 1)] > -1

        """
        C = self.data.copy()

        # this condition limits covariances to max -100 %
        mask = C.values < -1

        if mask.any():
            n = mask.shape[0]
            iu = np.triu_indices(n, k=1)
            how_many_bad_values = int(mask[iu].sum())
            smallest_bad_value = float(C.values[mask].min())
            msg = (
                f"Condition COV + 1 > 0 for Lognormal sampling is not respected.\n"
                f"{how_many_bad_values} off-diagonal covariance coefficients "
                f"are set to -1+eps. Smallest covariance is {smallest_bad_value:.5f}."
            )

            if "MT" in C.index.names:
                rows_bad = np.unique(iu[0][mask[iu]])
                bad_mts = C.index.get_level_values("MT")[rows_bad].unique().tolist()
                msg += f" Concerned MT numbers (rows): {bad_mts}."

            logging.warning(msg)
            
            # use pandas mask to avoid read-only issues with ".values[...] = "
            C = C.mask(mask, -1 + np.finfo(np.float64).eps)

        return self.__class__(C)

    def transform_lognormal(self):
        """
        Assuming that `self` constains the covariance matrix of a lognormal
        multivariate distribution (X-space) centered in a unit vector, this
        method applies a transformation to calculate the covariance matrix 
        of the underlying normal distribution (Z-space).
    
        **Intended use:**
        - The input `self.data` is the covariance of strictly positive variables
          `X = exp(Z - 0.5 * diag(Σ_Z))`, i.e., a *lognormal* model with **unit mean**
          for each component (E[X_i] = 1).
        - The covariance is in **relative units**, i.e. it already corresponds to
          Cov(X_i, X_j) / (μ_i μ_j) with μ_i = 1, so simply Cov(X_i, X_j).
        - The matrix is symmetric, and all entries satisfy **C_ij > -1**, ensuring
          `log1p(C_ij)` is defined.
    
        **Mathematical background**
        ---------------------------
        Under the unit-mean lognormal parameterization
            X_i = exp(Z_i - 0.5 * σ_i^2),  with  Z ~ N(0, Σ_Z),
        the covariance in X-space satisfies
            Cov(X_i, X_j) = exp(Σ_Z,ij) - 1.
        Therefore the inverse mapping is elementwise:
            Σ_Z,ij = log(1 + Cov(X_i, X_j)) = log1p(C_ij).
        The function is taken from https://doi.org/10.1016/j.nima.2012.06.036
    
        Returns
        -------
        :obj:`~sandy.cov.CategoryCov`
            New instance of the same class with the **Normal-space** covariance
            matrix Σ_Z = log1p(C), preserving index/columns metadata.
    
        Notes
        -----
        - This transform is **not** the general lognormal-to-normal covariance mapping
          for arbitrary means. If means are not 1, first normalize the covariance by
          μ_i μ_j and then apply this method:
              Σ_Z,ij = log(1 + Cov(X_i, X_j) / (μ_i μ_j)).
        - If you later regularize, do it in Σ_Z (Normal-space), not C (X-space).

        Examples
        --------
        
        Apply the transformation to a zero matrix.

        >>> import sandy
        >>> a = sandy.CategoryCov([[0, 0], [0, 0]], index=["A", "B"], columns=["C", "D"])
        >>> a_log = a.transform_lognormal()
        
        Check that values are correct.

        >>> expected = [[0, 0], [0, 0]]
        >>> np.testing.assert_array_equal(a_log.data.values, expected)
        
        Check that indices and columns are the same before and after transformation.

        >>> assert a.data.index.equals(a_log.data.index)
        >>> assert a.data.columns.equals(a_log.data.columns)

        Test the transformation for another simple matrix.

        >>> val = np.e - 1
        >>> a = sandy.CategoryCov([[val, val], [val, val]])
        >>> a_log = a.transform_lognormal()

        Check that values are again correct.

        >>> expected = [[1, 1], [1, 1]]
        >>> np.testing.assert_array_equal(a_log.data.values, expected)
        
        """
        # Copy to preserve metadata (index/columns) and avoid mutating `self.data`
        C = self.data.copy()
    
        # Elementwise inverse mapping: Σ_Z = log(1 + C)
        # Use log1p for numerical stability on small values
        C.loc[:, :] = np.log1p(C.values)
    
        return self.__class__(C)

    def draw_sample(self, N, lhs=False, verbose=False, seed=None,):
        """
        Draw `N` multivariate samples from this covariance matrix.
    
        Sampling is performed using an SVD-based factorization:
    
            C = U S Uᵀ
            X = U sqrt(S) Z
    
        where Z are IID standard normal samples (or LHS samples when
        `lhs=True`). Rows/columns corresponding to zero variance are
        automatically excluded and then padded back.
    
        Parameters
        ----------
        N : int
            Number of samples to draw.
        lhs : bool, optional
            Use Latin Hypercube Sampling (LHS) instead of IID Gaussian.
            Default is False (IID Gaussian).
        verbose : bool, optional
            Print diagnostic information. Default is False.
        seed : int or None, optional
            Random seed. If None, a seed is obtained from sandy.get_seed().
    
        Returns
        -------
        :obj:`~sandy.samples.Samples`
            A Samples object containing an (M × N) DataFrame of samples.
    
        Examples
        --------
        
        Draw a random sample.

        >>> import sandy
        >>> vals = [[4, 2.4],[2.4, 9]]
        >>> a = sandy.CategoryCov(vals)
        >>> s = a.draw_sample(3, seed=1)

        Check size.

        >>> assert s.data.shape == (2, 3)
    
        LHS sampling also works.

        >>> s = a.draw_sample(3, lhs=True, seed=1)
        >>> assert s.data.shape == (2, 3)
        
        Test that mean is about zero (gaussian sampling).
        
        >>> s = a.draw_sample(50_000, seed=123)
        >>> mu = s.data.mean(axis=1).round(2).tolist()
        >>> expected = [0.0, 0.0]
        >>> np.testing.assert_allclose(mu, expected, atol=0.05)
        
        Test that the sample estimate covariance matrix approximates the
        input covariance matrix.

        >>> smp_cov = s.get_cov().values.round(1)
        >>> np.testing.assert_array_equal(smp_cov, vals)
        
        Test reproducibility.

        >>> s1 = a.draw_sample(1000, seed=33).data
        >>> s2 = a.draw_sample(1000, seed=33).data
        >>> assert np.allclose(s1, s2)
     
        Test LHS also respects covariance matrix.
        
        >>> s_lhs = a.draw_sample(50_000, lhs=True, seed=123)
        >>> smp_cov = s_lhs.get_cov().values.round(1)
        >>> np.testing.assert_array_equal(smp_cov, vals)

        """

        from scipy.stats import norm
        from scipy.stats.qmc import LatinHypercube

        from sandy import get_seed           # lazy import
        from .samples import Samples

        M = self.size
        C = self.data.to_numpy()

        # -- Prepare index and columns for Samples object
        index = self.data.index
        columns = list(range(N))

        # --- Identify nonzero-variance dimensions
        # -- Reduce matrix size by removing rows and columns with zero on diag
        D = np.diag(C)
        nz = np.flatnonzero(D)
        Cr = C[nz][:, nz]
        
        # here we don't handle a covariance matrix with all zero uncertainties.

        # -- Decompose covariance (SVD better than QR or cholesky)
        Ur, Sr, _ = np.linalg.svd(Cr, hermitian=True)  # hermitian is twice faster (U5 from JEFF33, 240 groups)

        # This is the dimension of non-zero singular values
        Mr = Sr.size

        # -- Get U back to original size
        U = np.zeros((M, Mr))
        U[nz] = Ur

        # --- Generate standard-normal samples, IID samples with mu=0 and std=1
        seed_ = seed if seed is not None else get_seed()

        if lhs:
            engine = LatinHypercube(d=Mr, seed=seed_)
            lhd = engine.random(n=N)
            Z = norm.ppf(lhd).T  # (Mr × N), loc=0, scale=1
        else:
            rng = np.random.default_rng(seed=seed_)
            Z = rng.standard_normal(size=(Mr, N))

        # --- Optional diagnostics
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
        # --- Covariance factor application (fast version)
        sqrtSr = np.sqrt(Sr)              # (Mr,)
        L = U * sqrtSr                    # broadcasting, shape (M × Mr)
        X = L @ Z                         # dense BLAS GEMM, shape (M × N)

        # --- Build Samples object
        samples = pd.DataFrame(X, index=index, columns=columns)
        return Samples(samples)

    def get_std(self):
        """
        Extract standard deviations.

        Returns
        -------
        `pandas.Series`
            1d array of standard deviations

        Examples
        --------

        Extract standard deviation vector from covariance matrix.

        >>> import sandy
        >>> idx = ["A", "B"]
        >>> a = sandy.CategoryCov([[1, 0.4],[0.4, 1]], index=idx, columns=idx)
        >>> std = a.get_std()
        >>> expected = [1, 1]
        >>> np.testing.assert_array_equal(std.values, expected)

        Check that the series name is correct.

        >>> assert std.name == "STD"

        Check that the series indices match those of the covariance matrix.

        >>> assert std.index.equals(a.data.index)

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
        They are reported in ascending order.

        >>> import sandy
        >>> eigs = sandy.CategoryCov([[1, 0.4], [0.4, 1]]).get_eig()[0]
        >>> expected = [0.6, 1.4]
        >>> np.testing.assert_array_equal(eigs, expected)
        >>> assert eigs.name == "EIG"
    
        Extract eigenvectors.

        >>> eigv = sandy.CategoryCov([[1, 0.4], [0.4, 1]]).get_eig()[1]
        >>> val = 0.707106
        >>> np.testing.assert_array_almost_equal(eigv, [[-val, val], [val, val]], decimal=6)

        Replace small eigenvalues using a tolerance.
        
        >>> cov = sandy.CategoryCov([[0.1, 0.1], [0.1, 1]])
        >>> eig, eigv = cov.get_eig()
        >>> eig_t, eigv_t = cov.get_eig(tolerance=0.1)
        >>> expected = [8.90228e-02, 1.01098]
        >>> np.testing.assert_array_almost_equal(eig, expected, decimal=5)
        >>> expected = [0.00000, 1.01098]
        >>> np.testing.assert_array_almost_equal(eig_t, expected, decimal=5)
        >>> np.testing.assert_array_equal(eigv, eigv_t)

        Handle negative eigenvalues.

        >>> cov = sandy.CategoryCov([[1, 2], [2, 1]])
        >>> eig = cov.get_eig()[0]
        >>> expected = [-1, 3]
        >>> np.testing.assert_array_almost_equal(eig, expected, decimal=6)
    
        Replace negative eigenvalues with zero.

        >>> eig_t = cov.get_eig(tolerance=0)[0]
        >>> expected = [0, 3]
        >>> np.testing.assert_array_almost_equal(eig_t, expected, decimal=6)
    
        Example with a covariance matrix (before they were correlation matrices).

        >>> cov = sandy.CategoryCov([[1, 0.2, 0.1], [0.2, 2, 0], [0.1, 0, 3]])
        >>> eig = cov.get_eig()[0]
        >>> expected = [9.56764e-01, 2.03815e+00, 3.00509e+00]
        >>> np.testing.assert_array_almost_equal(eig, expected, decimal=5)
    
        Real test on H1 file.

        >>> endf6 = sandy.get_endf6_file("jeff_33", "xs", 10010, local=True)
        >>> ek = sandy.energy_grids.CASMO12
        >>> err = endf6.get_errorr(errorr_kws=dict(ek=ek), err=1)["errorr33"]
        >>> cov = err.get_cov()
        >>> eig = cov.get_eig()[0].sort_values(ascending=False)

        Check the largest eigenvalues.

        >>> large_eig = eig.head(7)
        >>> expected = [3.66411e-01, 7.05311e-03, 1.55346e-03, 1.60175e-04,
        ...             1.81374e-05, 1.81078e-06, 1.26691e-07]
        >>> np.testing.assert_array_almost_equal(large_eig, expected, decimal=5)

        Check the smallest eigenvalues.

        >>> small_eig = eig.tail(7)
        >>> expected = [-8.27624942e-17, -4.05553200e-13, -1.13485938e-12, -1.79747459e-12,
        ...             -3.37064670e-12, -1.15488784e-11, -3.99319980e-11]
        >>> np.testing.assert_array_almost_equal(small_eig, expected, decimal=7)

        Ensure all eigenvalues are non-negative when using `tolerance=0`.

        >>> eig_t = cov.get_eig(tolerance=0)[0]
        >>> assert (eig_t >= 0).all()
        """
        from scipy.linalg import eigh

        E, V = eigh(self.data)
        E = pd.Series(E, name="EIG")
        V = pd.DataFrame(V)
    
        if tolerance is not None:
            E[E / E.max() < tolerance] = 0
    
        return E, V

    def get_corr(self):
        """
        Return the correlation matrix corresponding to this covariance matrix.
    
        Returns
        -------
        :obj:`~sandy.cov.CategoryCov`
            Correlation matrix with the same index/columns as the covariance.
    
        Notes
        -----
        The correlation is computed as:
    
            Corr[i,j] = Cov[i,j] / (sqrt(Cov[i,i]) * sqrt(Cov[j,j]))
    
        Zero variances produce zero rows/columns in the result.

        Examples
        --------

        Extract correlation matrix.
        
        >>> import sandy
        >>> idx = ["A", "B"]
        >>> a = sandy.CategoryCov([[4, 2.4],[2.4, 9]])
        >>> corr = a.get_corr()

        Test that values are correct.

        >>> expected = [[1, 0.4], [0.4, 1]]
        >>> np.testing.assert_array_almost_equal(corr.data.values, expected, decimal=10)

        Check that indices and columns are the same as in the covariance matrix.

        >>> assert a.data.index.equals(corr.data.index)
        >>> assert a.data.columns.equals(corr.data.columns)
        """
        
        cov = self.data.to_numpy()  # always float
        std = self.get_std().to_numpy()  # always float

        # Compute inverse std safely
        with np.errstate(divide='ignore', invalid='ignore'):
            invstd = 1.0 / std
            invstd[~np.isfinite(invstd)] = 0.0   # handle inf, -inf, NaN

        # Compute correlation by outer product
        corr = cov * np.outer(invstd, invstd)
    
        # Preserve index/columns and return same class
        df = pd.DataFrame(corr, index=self.data.index, columns=self.data.columns)
        return self.__class__(df)

    def get_L(self, tolerance=None):
        """
        Return a lower-triangular matrix L such that L @ L.T approximates the
        covariance matrix. If the matrix is not PSD, eigenvalues below `tolerance`
        are replaced by zero.

        Behavior:
        - If `tolerance` is provided, clipping is delegated to `get_eig(tolerance)`.
        - If `tolerance` is None and the matrix is not PSD, raise ValueError.
    
        Parameters
        ----------
        tolerance : float, optional, default is `None`
            If provided, eigenvalues < tolerance are set to zero by `get_eig`.
            If None, negative eigenvalues cause a ValueError.

        Returns
        -------
        pandas.DataFrame
            Lower-triangular matrix L with shape (n, n) and original index/columns.

        Examples
        --------

        Decompose a positive definite matrix.

        >>> import sandy
        >>> a = np.array([[4, 12, -16], [12, 37, -43], [-16, -43, 98]])
        >>> L = sandy.CategoryCov(a).get_L()
        >>> np.testing.assert_array_almost_equal(L @ L.T, a, decimal=6)

        For a PSD matrix, eigenvalues are non-negative, so tolerance has no effect.

        >>> L_t = sandy.CategoryCov(a).get_L(tolerance=0)
        >>> np.testing.assert_equal(L.values, L_t.values)

        For a non-PSD covariance, negative eigenvalues are clipped to zero.
        
        >>> a = np.array([[1, -2],[-2, 3]])
        >>> L = sandy.CategoryCov(a).get_L(tolerance=0)
        >>> assert not np.allclose(L @ L.T, a)   # it's the PSD projection

        Raise error if tolerance is not given.
        
        >>> import pytest
        >>> with pytest.raises(ValueError) as excinfo:
        ...    sandy.CategoryCov(a).get_L()
        >>> assert "Matrix is not PSD" in str(excinfo.value)
        """
        index = self.data.index
        columns = self.data.columns

        # --- 1) Eigenpairs with possible clipping ---
        eigvals, eigvecs = self.get_eig(tolerance=tolerance)
        # from pandas to numpy
        eigvals_np = eigvals.values.astype(float)
        eigvecs_np = eigvecs.values.astype(float)

        # --- 2) PSD check if no tolerance was provided ---
        if tolerance is None and (eigvals_np < 0).any():
            raise ValueError(
                "Matrix is not PSD. Provide `tolerance` to repair it "
                "or ensure the matrix is PSD."
            )

        # --- 3) Keep only positive eigenvalues (after clipping) ---
        pos_mask = eigvals_np > 0
        # this is the size of the n x r L-matrix, the space without clipped eigenvalues
        r = pos_mask.sum()
        
        # the case where all eigenvalues are clipped is not handled

        # --- 4) Build thin A = V_pos * sqrt(lambda_pos) ---
        V_pos = eigvecs_np[:, pos_mask]          # (n × r)
        sqrt_lam = np.sqrt(eigvals_np[pos_mask]) # (r,)
        A = V_pos * sqrt_lam                     # (n × r)
    
        # --- 5) Thin QR: A = Q R  =>  A A^T = R^T R ---
        Q, R = qr(A.T, mode="reduced")
        L_thin = R.T  # (n × r)
    
        # --- 6) Build FULL n×n L by zero-padding unused columns ---
        L_full = np.zeros((len(index), len(columns)))
        L_full[:, :r] = L_thin  # remaining columns stay zero
    
        return pd.DataFrame(L_full, index=index, columns=columns)

    def sampling(self, nsmp, seed=None, lognormal=True, correction=0.5/100,
                 lhs=False, verbose=False, **kwargs):
        """
        Extract perturbation coefficients from the covariance matrix using either
        a normal or lognormal distribution. Samples are adjusted to ensure physical
        plausibility (e.g., positivity, bounded range).
    
        Parameters
        ----------
        nsmp : int
            Number of samples to draw. Floats (e.g. 1e5) are cast to int.
        seed : int, optional
            Seed for the random number generator (default is None).
        lognormal : bool, optional
            If True, use lognormal sampling. Otherwise, use (truncated) normal.
        correction : float, optional
            Regularization factor passed to `regularize` to improve conditioning
            (default is 0.5%).
        lhs : bool, optional
            If True, use Latin Hypercube Sampling (default is False).
        verbose : bool, optional
            If True, print diagnostic information during sampling.
    
        Returns
        -------
        :obj:`~sandy.samples.Samples`
            An object containing the sampled perturbation coefficients.
    
        Notes
        -----
        - Normal sampling produces relative perturbations around 1 and is then
          truncated to [0, 2] to preserve basic physical bounds; this truncation
          can bias the covariance if uncertainties are large.
        - Lognormal sampling returns strictly positive factors with mean ~ 1 
          by construction (mean-centering in log-space).
    
        Examples
        --------

        Common setup.
    
        >>> import sandy
        >>> seed, nsmp = 11, 100_000
        >>> index = columns = ["A", "B"]
        >>> c = pd.DataFrame([[1, 0.4],[0.4, 1]], index=index, columns=columns) / 10
        >>> cov = sandy.CategoryCov(c)
    
        Normal sampling reproduces the mean and covariance approximately.
    
        >>> smp_n = cov.sampling(nsmp, seed=seed, lognormal=False)
        >>> expected_mean = [1, 1]
        >>> np.testing.assert_array_almost_equal(smp_n.get_mean(), expected_mean, decimal=2)
        >>> np.testing.assert_array_almost_equal(smp_n.get_cov(), c, decimal=2)
    
        Lognormal sampling also reproduces the targets.
    
        >>> smp_ln = cov.sampling(nsmp, seed=seed, lognormal=True)
        >>> np.testing.assert_array_almost_equal(smp_ln.get_mean(), expected_mean, decimal=2)
        >>> np.testing.assert_array_almost_equal(smp_ln.get_cov(), c, decimal=2)
    
        Reproducibility with a fixed seed.
    
        >>> smp_n2 = cov.sampling(nsmp, seed=seed, lognormal=False)
        >>> assert smp_n2.data.equals(smp_n.data)
        >>> smp_ln2 = cov.sampling(nsmp, seed=seed, lognormal=True)
        >>> assert smp_ln2.data.equals(smp_ln.data)
    
        For large variances, normal sampling is truncated and does not reproduce the full covariance.
    
        >>> c = pd.DataFrame([[2, 0],[0, 2]])
        >>> s = sandy.CategoryCov(c).sampling(nsmp, lognormal=False)
        >>> assert not np.allclose(s.get_cov(), c, atol=1)  # due to truncation

        ...but the mean is preserved.
        
        >>> np.testing.assert_array_almost_equal(s.get_mean(), expected_mean, decimal=2)

        ...and the sample standard deviations are smaller than 1.

        >>> assert (s.get_rstd().values < 1).all()

        ...in this particular case, the relative error between sample covariance
        matrix and original one is large.

        >>> rel_err = np.linalg.norm(s.get_cov() - c) / np.linalg.norm(c)
        >>> assert rel_err > 0.5
    
        For lognormal sampling, large variances remain well-behaved.
    
        >>> s = sandy.CategoryCov(c).sampling(nsmp, lognormal=True)

        This is tested by checking any mean shift.

        >>> np.testing.assert_array_almost_equal(s.get_mean(), expected_mean, decimal=2)

        ... and by checking that the standard deviations are preserved.

        >>> expected_std = np.sqrt(np.diag(c))
        >>> np.testing.assert_allclose(s.get_rstd().values, expected_std, rtol=0.2)

        ... and also by checking that the eigenvalues of the covariance matrix are preserved.

        >>> eigvals_original = np.linalg.eigvalsh(c)
        >>> eigvals_sampled = np.linalg.eigvalsh(s.get_cov())
        >>> np.testing.assert_allclose(eigvals_sampled, eigvals_original, rtol=0.2)

        """
        N = int(nsmp)
        if N <= 0:
            raise ValueError(f"'nsmp' must be > 0, got {nsmp}")


        if lognormal:
            # 1) Start from the "relative covariance" C
            C = self.correct_lognormal()

            var = np.diag(C.data.to_numpy())  # this will be used to adjust the mean
            
            
            # 2) Compute log-space mean shift to enforce E[exp(Y)] ≈ 1
            # mean of the underlying normal distribution
            # https://stats.stackexchange.com/questions/573808/intuition-for-why-mean-of-lognormal-distribution-depends-on-variance-of-normally
            umu = -0.5 * np.log1p(var)  # (m,) identical to np.log(1 / np.sqrt(var + 1))
            umu = umu[:, None]          # (m, N) identical to umu.reshape(var.size, -1)

            # 3) Transform covariance to log-space and draw samples
            samples = (
                C.transform_lognormal()
                .regularize(correction=correction)
                .draw_sample(N, lhs=lhs, verbose=verbose, seed=seed)
                .apply_function(lambda x: x + umu)
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

        >>> import sandy
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
        from .gls import sandwich

        s_ = pd.DataFrame(s)
        index = s_.index
        sandwich_ = sandwich(self.data.values, s_.values)
        if len(sandwich_.shape) == 0: 
            sandwich_ = [sandwich_]
        sandwich_ = pd.DataFrame(sandwich_, index=index, columns=index)
        return self.__class__(sandwich_)

    def corr2cov(self, std):
        """
        Convert a correlation matrix into a covariance matrix using
        the vector of standard deviations.
        
        Use function :func:`~sandy.cov.corr2cov`.
    
        Parameters
        ----------
        std : array_like, shape (n,)
            Standard deviations. Index is ignored if `pd.Series`.

        Returns
        -------
        :obj:`~sandy.cov.CategoryCov`
            Covariance matrix with the same index/columns as the correlation.

        Examples
        --------

        Initialize index and columns.

        >>> import sandy
        >>> idx = ["A", "B", "C"]
        >>> std = np.array([1, 2, 3])
        >>> a = [[1, 0, 2], [0, 1, 0], [2, 0, 1]]
        >>> corr = sandy.CategoryCov(a, index=idx, columns=idx)
        >>> cov = corr.corr2cov(std)
        
        Check that indices and columns are the same as for the correlation matrix.

        >>> assert corr.data.index.equals(cov.data.index)
        >>> assert corr.data.columns.equals(cov.data.columns)

        """
        corr_ = self.data.to_numpy()
        s_ = np.asarray(std)

        cov = corr2cov(corr_, s_)

        return self.__class__(cov, index=self.data.index, columns=self.data.columns)

    def to_excel(self, file):
        """
        Save the sample dataset to an Excel file.

        This method exports the sample data to an Excel file, writing the dataset to 
        the 'COV' sheet. If the file already exists, the sheet is replaced; otherwise, 
        a new file is created.
        
        Replaces the sheet if file exists; overwrites invalid files.

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
        
        Examples
        --------
        
        Create excel file.
        
        >>> import sandy, tempfile
        >>> tmp = tempfile.NamedTemporaryFile(suffix=".xlsx", delete=False)
        >>> fn = tmp.name
        >>> tmp.close()
        >>> idx = ["A", "B"]
        >>> a = [[1, .4], [.4, 1]]
        >>> c = sandy.CategoryCov(a, index=idx, columns=idx)
        >>> c.to_excel(fn)

        Check that file is created.

        >>> assert os.path.exists(fn)

        If we write twice to the same file, sheet_name `"COV"` is overwritten.

        >>> b = [[1, .2], [.2, 1]]
        >>> c1 = sandy.CategoryCov(a, index=idx, columns=idx)
        >>> c2 = sandy.CategoryCov(b, index=idx, columns=idx)
        >>> c1.to_excel(fn)
        >>> c2.to_excel(fn)   # Should replace sheet "COV"

        >>> out = pd.read_excel(fn, sheet_name="COV", index_col=0)
        >>> assert out.equals(c2.data)

        Check that indices are preserved.
        
        >>> with tempfile.TemporaryDirectory() as td:
        ...     fn = os.path.join(td, "cov.xlsx")
        ...     c.to_excel(fn)
        ...     out = pd.read_excel(fn, sheet_name="COV", index_col=0)
        ...     assert list(out.index) == ["A", "B"]

        Try with something more similar to a real covariance matrix.
        
        First, prepare a MultiIndex for rows and columns.

        >>> rows = pd.MultiIndex.from_product([["g1","g2"], ["A","B"]], names=["group","label"])
        >>> cols = rows

        Then, build a 4x4 symmetric correlation-like matrix.

        >>> base = np.array([[1.0, 0.2, 0.3, 0.1],
        ...                  [0.2, 1.0, 0.4, 0.2],
        ...                  [0.3, 0.4, 1.0, 0.5],
        ...                  [0.1, 0.2, 0.5, 1.0]])
        >>> corr_df = pd.DataFrame(base, index=rows, columns=cols)
        >>> corr = sandy.CategoryCov(corr_df)
        
        Use a temp file path that already exists (empty) to test overwrite of invalid `.xlsx`.

        >>> tmp = tempfile.NamedTemporaryFile(suffix=".xlsx", delete=False)
        >>> fn = tmp.name
        >>> tmp.close()
        
        Write MultiIndex covariance/correlation to Excel.

        >>> corr.to_excel(fn)
        
        File exists and is a valid Excel now.

        >>> assert os.path.exists(fn)
        
        Read back as MultiIndex: header=[0,1] for columns, index_col=[0,1] for rows.

        >>> out = pd.read_excel(fn, sheet_name="COV", header=[0,1], index_col=[0,1])
        
        Indices and columns are preserved.
        
        >>> assert corr.data.index.equals(out.index)
        >>> assert corr.data.columns.equals(out.columns)
        
        Numeric content is preserved (within float tolerance).

        >>> np.testing.assert_allclose(out.to_numpy(), corr.data.to_numpy())

        """
        # Resetting indices messes everything up
        df = self.data
        
        file_exists = os.path.exists(file)
        
        # If file exists and is NOT a valid Excel file, remove it.
        # This avoids BadZipFile without complex logic.
        if file_exists:
            try:
                # Try opening the file as an Excel file
                pd.ExcelFile(file)   #  will fail on empty/invalid files
            except Exception:
                os.remove(file)
                file_exists = False   # force write mode

        # Determine write mode
        mode = "a" if file_exists else "w"

        # Build writer kwargs
        kwargs_ = dict(mode=mode, engine="openpyxl")

        # Only append mode accepts if_sheet_exists
        if file_exists:
            kwargs_["if_sheet_exists"] = "replace"

        # Write to Excel
        with pd.ExcelWriter(file, **kwargs_) as writer:
            df.to_excel(writer, sheet_name="COV")

def corr2cov(corr, s):
    """
    Convert a correlation matrix into a covariance matrix using
    the vector of standard deviations.

    Covariance is defined as:

        Cov[i, j] = Corr[i, j] * s[i] * s[j]

    This is computed efficiently via outer products.


    Parameters
    ----------
    corr : array_like, shape (n, n)
        Correlation matrix.
    s : array_like, shape (n,)
        Standard deviations.


    Returns
    -------
    `numpy.ndarray`
        Covariance matrix of shape (n, n).

    Examples
    --------

    Basic test.

    >>> s = np.array([1, 2, 3])
    >>> corr = np.array([[1, 0, 2], [0, 1, 0], [2, 0, 1]])
    >>> cov = corr2cov(corr, s).astype(int)
    >>> expected = [[ 1,  0,  6], [ 0, 4,  0], [ 6,  0,  9]]
    >>> np.testing.assert_array_equal(cov, expected)

    Standard deviations must reproduce correctly.
    
    >>> np.testing.assert_array_equal(np.sqrt(np.diag(cov)), s)

    Symmetry check.
    
    >>> np.testing.assert_array_equal(cov, cov.T)
    """
    # also used in sandy.fy for the CEA covariance matrices

    corr_ = np.asarray(corr)
    s_ = np.asarray(s)
    return corr_ * np.outer(s_, s_)
