import pdb
import functools

import numpy as np
import pandas as pd

import sandy

__author__ = "Luca Fiorito"
__all__ = [
#        "EnergyCov",
        ]

class EnergyCov(sandy.formats.utils.BaseCov):
    """Dataframe for a multigroup covariance matrix.
        
    .. note:: It is assumed that the covariance matrix is defined over 
              multi-group energy grids.

              Only 'zero' interpolation is supported.

    Attributes
    ----------
    index : `pandas.Index` of `float`
        energy grid for the 1st reaction
        .. important:: index values must be monotonically increasing
    columns : `pandas.Index` of `float`
        energy grid for the 2nd reaction
        .. important:: columns values must be monotonically increasing
    `values` : `numpy array`
        matrix coefficients
    
    Methods
    -------
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
        self.index = pd.Float64Index(self.index, name="E")
        self.columns = pd.Float64Index(self.columns, name="E")
        if not (self.index.values == np.sort(self.index.values)).all():
            raise sandy.Error("index values are not monotonically increasing")
        if not (self.columns.values == np.sort(self.columns.values)).all():
            raise sandy.Error("columns values are not monotonically increasing")
    
    def change_grid(self, ex, ey):
        """Given one energy grid for the x-axis and one energy grid for the 
        y-axis, interpolate/extrapolate the covariance matrix over the new 
        points using the *forward-filling* method.
        
        .. important::
            
            * backward extrapolated values (e.g. below threshold) are replaced by 0
            * forward extrapolated values (e.g. above 20 MeV) are replaced by 
              the covariance coefficient that refers to the last point in the 
              original grid
        
        Parameters
        ----------
        ex : `iterable`
            covariance energy grid for the x-axis (first reaction)
        ey : `iterable`
            covariance energy grid for the y-axis (second reaction)
        
        Returns
        -------
        `sandy.formats.utils.EnergyCov`
            Covariance matrix interpolated over the new axes.
        """
        df = self.reindex(index=ex, method="ffill"). \
                  reindex(columns=ey, method="ffill"). \
                  fillna(0)
        return self.__class__(df)
    
    def _get_mesh(self):
        X, Y = np.meshgrid(self.index.values, self.columns.values)
        return X.T, Y.T
    
    def _plot_matrix(self, xscale='log', yscale='log', cmap='bwr', vmin=-1, vmax=1, **kwargs):
        ax = plt.pcolormesh(*self._get_mesh(), self.values, cmap=cmap, vmin=vmin, vmax=vmax, **kwargs)
        plt.colorbar(ax)
        plt.gca().set_xscale(xscale)
        plt.gca().set_yscale(yscale)

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
        """
        def foo(x, y):
            ex = sorted(set(x.index.tolist() + y.index.tolist()))
            ey = sorted(set(x.columns.tolist() + y.columns.tolist()))
            x_ = x.change_grid(ex, ey)
            y_ = y.change_grid(ex, ey)
            return cls(x_.add(y_))
        df = functools.reduce(lambda x,y: foo(x,y), covs)
        return cls(df)

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
        """
        Extract square symmetric covariance matrix from NI-type sub-subsection
        data with flag `lb=5`.

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
        """Extract square asymmetric covariance matrix from NI-type sub-subsection data 
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