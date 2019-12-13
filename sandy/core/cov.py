import pdb
import logging

import numpy as np
import scipy
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

import sandy

__author__ = "Luca Fiorito"
__all__ = [
        "CategoryCov",
        "triu_matrix",
        "corr2cov",
        "random_corr",
        "random_cov",
        "get_eig",
        ]

class _Cov(np.ndarray):
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



class CategoryCov():
    
    def __init__(self, df, xlabel=None, ylabel=None):
        self.data = df
        self.xlabel = xlabel
        self.ylabel = ylabel
    
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
            if 'data' is not a 'pandas.DataFrame'
        """
        return self._data
    
    @data.setter
    def data(self, data):
        if not isinstance(data, pd.DataFrame) and not isinstance(data, np.ndarray):
            raise sandy.Error("'data' is not a 'pandas.DataFrame' or a 'numpy.ndarray'")
        if isinstance(data, np.ndarray):
            df = pd.DataFrame(data, dtype=float)
        else:
            df = data.astype(float)
        self._data = df

    def eig(self, sort=True):
        """
        Extract eigenvalues and eigenvectors.
        
        Returns
        -------
        `Pandas.Series`
            real part of eigenvalues sorted in descending order
        `pandas.DataFrame`
            matrix of eigenvectors
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
        """
        cov = self.data.values
        with np.errstate(divide='ignore', invalid='ignore'):
            coeff = np.true_divide(1, self.std.values)
            coeff[ ~ np.isfinite(coeff)] = 0  # -inf inf NaN
        corr = np.multiply(np.multiply(cov, coeff).T, coeff)
        return pd.DataFrame(corr).reindex_like(self.data)

    def sampling(self, nsmp, seed=None):
        """Extract random samples from the covariance matrix, either using
        the cholesky or the eigenvalue decomposition.

        Parameters
        ----------
        nsmp : `int`
            number of samples
        seed : `int`, optional, default is `None`
            seed for the random number generator (by default use `numpy` 
            dafault pseudo-random number generator)

        Returns
        -------
        `pandas.DataFrame`
            2D array of random samples
        """
        np.random.seed(seed=seed)
        dim = self.data.shape[0]
        y = np.random.randn(dim, nsmp)
        L = self.decompose()
        samples = L.dot(y)
        df = pd.DataFrame(samples, index=self.data.index, columns=list(range(nsmp))).T
        return sandy.Samples(df)

    def decompose(self):
        """
        Extract lower triangular matrix `L` for which `L*L^T == COV`.
        
        Returns
        -------
        `numpy.ndarray`
            lower triangular matrix
        """
        E, V = self.eig(sort=False)
        E[E<=0] = 0
        M = V.values.dot(np.diag(np.sqrt(E.values)))
        Q, R = scipy.linalg.qr(M.T)
        return R.T
    
    def plot_corr(self, ax):
        kwargs = {"cbar" : True, "vmin" : -1, "vmax" : 1, "cmap" : "RdBu"}
        ax = sns.heatmap(self.corr, ax=ax, **kwargs)
        if self.xlabel:
            ax.set_xlabel("{}".format(self.xlabel))
        if self.ylabel:
            ax.set_ylabel("{}".format(self.xlabel))
        return ax
    
    @classmethod
    def corr2cov(cls, corr, std):
        """
        Produce covariance matrix given correlation matrix and standard 
        deviation array.
        
        Parameters
        ----------
        corr : 2d `numpy.ndarray`
            square 2D correlation matrix
        
        s : 1d `numpy.ndarray`
            array of standard deviations
    
        Returns
        -------
        `sandy.CategoryCov`
            covariance matrix
        """
        return cls(corr2cov(corr, std))



def corr2cov(corr, s):
    """
    Produce covariance matrix given correlation matrix and standard 
    deviation array.
    
    Parameters
    ----------
    corr : 2d `numpy.ndarray`
        square 2D correlation matrix
    
    s : 1d `numpy.ndarray`
        array of standard deviations

    Returns
    -------
    `numpy.ndarray`
        covariance matrix
    """
    dim = corr.shape[0]
    S = np.repeat(s, dim).reshape(dim, dim)
    return S.T * (corr * S)



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
    upper = np.triu(np.random.uniform(-1, 1, size**2).reshape(size, size), 1) if correlations else np.zeros([size, size])
    corr += upper + upper.T
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
    kwargs = {"cbar" : False, "vmin" : -1, "vmax" : 1, "cmap" : "RdBu"}
    fig,ax = plt.subplots(size, size, sharex="col", sharey="row")
    for (i,j), m in zip(coords, MM):
        if m is None:
            continue
        ax[i,j] = sns.heatmap(m, ax=ax[i,j], **kwargs)
        if i != j:
            ax[j,i] = sns.heatmap(m.T, ax=ax[j,i], **kwargs)
    return fig, ax