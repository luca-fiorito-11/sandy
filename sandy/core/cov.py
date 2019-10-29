import pdb
import logging

import numpy as np
import scipy as sp
import pandas as pd

import sandy

__author__ = "Luca Fiorito"
__all__ = [
        "Cov",
        "corr2cov",
        ]

class Cov(np.ndarray):
    """Covariance matrix treated as a `numpy.ndarray`.
    
    Methods
    -------
    corr
        extract correlation matrix
    corr2cov
        produce covariance matrix given correlation matrix and standard deviation 
        array
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
        """Extract eigenvalues and eigenvectors.
        
        Returns
        -------
        `Pandas.Series`
            real part of eigenvalues sorted in descending order
        `np.array`
            matrix of eigenvectors
        """
        E, V = sp.linalg.eig(self)
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
            coeff = np.true_divide( 1, std )
            coeff[ ~ np.isfinite(coeff)] = 0  # -inf inf NaN
        corr = np.multiply(np.multiply(self.T, coeff).T, coeff)
        return self.__class__(corr)

    def _reduce_size(self):
        nonzero_idxs =  np.flatnonzero(np.diag(self))
        cov_reduced = self[nonzero_idxs][:,nonzero_idxs]
        return nonzero_idxs, cov_reduced

    @classmethod
    def _restore_size(cls, nonzero_idxs, cov_reduced, dim):
        cov = Cov(np.zeros((dim, dim)))
        for i,ni in enumerate(nonzero_idxs):
            cov[ni,nonzero_idxs] = cov_reduced[i]
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
        logging.debug("covariance matrix dimension is {} X {}".format(*self.shape))
        dim = self.shape[0]
        np.random.seed(seed=seed)
        y = np.random.randn(dim, nsmp)
        nonzero_idxs, cov_reduced = self._reduce_size()
        L_reduced = cov_reduced.get_L()
        L = self.__class__._restore_size(nonzero_idxs, L_reduced, dim)
        samples = np.array(L.dot(y))
        return samples

    def get_L(self):
        """Extract lower triangular matrix `L` for which `L*L^T == self`.
        
        Returns
        -------
        `np.array`
            lower triangular matrix
        """
        try:
            L = sp.linalg.cholesky(self, lower=True, overwrite_a=False, check_finite=False)
        except np.linalg.linalg.LinAlgError:
            E, V = self.eig()
            E[E<=0] = 0
            Esqrt = np.diag(np.sqrt(E))
            M = V.dot(Esqrt)
            Q, R = sp.linalg.qr(M.T)
            L = R.T
        return L



def corr2cov(corr, s):
    """
    Produce covariance matrix given correlation matrix and standard 
    deviation array.
    
    Parameters
    ----------
    corr : `numpy.array`
        square 2D correlation matrix
    
    std : `numpy.array`
        array of standard deviations

    Returns
    -------
    `Cov`
        covariance matrix
    """
    _corr = Cov(corr)
    _std = std.flatten()
    dim = _corr.shape[0]
    S = np.repeat(_std, dim).reshape(dim, dim)
    cov = S.T * (_corr * S)
    return Cov(cov)



def triu_matrix(arr, size):
    """
    Given the upper triangular values of a **square symmetric** matrix in
    an array, return the full matrix.

    Inputs:
        - arr :
            (1d array) array with the upper triangular values of the matrix
        - size :
            (int) dimension of the matrix

    Outputs:
        - matrix :
            (2d array) reconstructed 2d-array with symmetric matrix
    """
    matrix = np.zeros([size, size])
    indices = np.triu_indices(size)
    matrix[indices] = arr
    matrix += np.triu(matrix, 1).T
    return matrix