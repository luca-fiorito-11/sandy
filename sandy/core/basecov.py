import pdb
import logging

import numpy as np
import pandas as pd

import sandy

__author__ = "Luca Fiorito"
__all__ = [
        "BaseCov",
        ]

class BaseCov(pd.DataFrame):
    """
    Base covariance class inheriting from `pandas.DataFrame`.
    Must be used as superclass by all other Nuclear Data Covariance Objects.
    
    Attributes
    ----------
    mat : `numpy.array`
        array of unique MAT number found in index
    mt : `numpy.array`
        array of unique MT number found in index
    index : `pandas.MultiIndex`
        index with at least level `"MAT"` and `"MT"`
    columns : `pandas.MultiIndex`
        index with at least level `"MAT"` and `"MT"`
    `values` : `numpy array`
        matrix coefficients

    Methods
    -------
    corr
        get correlation matrix instance with inherited class type
    eig
        get covariance matrix eigenvalues as a `pandas.Series` instance
    from_list
        extract global cross section/nubar covariance matrix from iterables 
    get_var
        get covariance matrix variances as a `pandas.Series` instance
    get_std
        get covariance matrix standard deviations as a `pandas.Series` instance
    to_matrix
        get covariance matrix as a `sandy.formats.utils.Cov` instance
    """

    @property
    def mat(self):
        return self.index.get_level_values("MAT").unique()

    @property
    def mt(self):
        return self.index.get_level_values("MT").unique()
    
    def to_matrix(self):
        """
        Extract dataframe values as a `Cov` instance
        
        Returns
        -------
        `sandy.Cov`
            covariance matrix as a `numpy` array
        """
        return sandy.Cov(self.values)

    def eig(self):
        """
        Extract eigenvalues in descending order.
        
        Returns
        -------
        `pandas.Series`
            sorted list of eigenvalues
        """
#        NE = np.extract(E < 0, E)    # extract negative eigenvalues
#        if len(NE) != 0:
#            neig_max = max(abs(NE))
#            eig_max = max(abs(E))
#            if neig_max/eig_max >= 0.1:
#                logging.warn("found large negative eigenvalues")
        E = self.to_matrix().eig()[0]
        return pd.Series(sorted(E, reverse=True), name='eigenvalues')

    def corr(self):
        """
        Extract correlation matrix.
        
        Returns
        -------
        `sandy.BaseCov` or daughter instance
            correlation matrix
        """
        corr = self.to_matrix().corr()
        return self.__class__(corr, index=self.index, columns=self.columns)

    def check_diagonal(self, verbose=True):
        """
        Check if any of the diagonal elements is negative.
        Return count of negative variances.
        
        Parameters
        ----------
        verbose : `bool`
            If `True` print list of negative variances

        Returns
        -------
        `int`
        """
        var = self.get_var()
        mask = var < 0
        count = mask.sum()
        if verbose and count > 0:
            string = var[mask].to_string()
            logging.warn("found {} negative variances\n{}".format(count, string))
        return count

    def get_var(self):
        """
        Extract diagonal.

        Returns
        -------
        `pandas.Series`
            array of variances
        """
        return pd.Series(np.diag(self.values), index=self.index, name="VAR")

    def get_std(self):
        """
        Extract square root of diagonal.
        
        Returns
        -------
        `pandas.Series`
            array of standard deviations
        """
        return self.get_var().apply(np.sqrt).rename("STD")

    def filter_by(self, index_key, index_values, columns_key, columns_values):
        """
        Filter dataframe based on given index and allowed values.

        .. hint:: use this method to filter the dataframe other than `.loc` as 
                  it returns a `BaseCov` (or daughter) instance.
        
        Parameters
        ----------
        index_key : `str`
            index level name to which to consider, e.g. "MAT", "MT"
        index_values : iterable
            list of values to include in the filtered matrix
        columns_key : `str`
            columns level name to which to consider, e.g. "MAT", "MT"
        columns_values : iterable
            list of values to include in the filtered matrix
        
        Returns
        -------
        `BaseCov` or daughter instance
        """
        index_cond = self.index.get_level_values(index_key).isin(index_values)
        columns_cond = self.index.get_level_values(columns_key).isin(columns_values)
        df = self.loc[index_cond, columns_cond]
        if df.empty:
            raise SandyError("applied filter returned empty matrix")
        return self.__class__(df)

    def _stack_correlations(self):
        corrstack = self.corr().T.reset_index(drop=True).T.reset_index(drop=True).stack()
        index = self.index.to_flat_index()
        pdb.set_trace()
        multiindex = pd.MultiIndex.from_product([index.values, index.values])
        idx0 = pd.DataFrame(multiindex.get_level_values(0).tolist(), columns=["MAT", "MT", "E"])
        idx1 = pd.DataFrame(multiindex.get_level_values(1).tolist(), columns=["MAT", "MT", "E"])
        idx0["MAT1"] = idx1.MAT
        idx0["MT1"] = idx1.MT
        idx0["E1"] = idx1.E
        corrstack.index = idx0.set_index(["MAT","MT","E", "MAT1","MT1", "E1"]).index
        return corrstack

    @classmethod
    def _from_list(cls, iterable):
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
        covs = pd.DataFrame.from_records(iterable, columns=columns).set_index(columns[:-1]).COV
        for (keys_rows,keys_cols), cov in covs.iteritems():
            if keys_rows == keys_cols: # diagonal terms
                if cov.shape[0] != cov.shape[1]:
                    raise SandyError("non-symmetric covariance matrix for ({}, {})".format(keys_rows, keys_cols))
                if not np.allclose(cov, cov.T):
                    raise SandyError("non-symmetric covariance matrix for ({}, {})".format(keys_rows, keys_cols))
            else: # off-diagonal terms
                condition1 = (keys_rows,keys_rows) in covs.index
                condition2 = (keys_cols,keys_cols) in covs.index
                if not (condition1 and condition2):
                    covs[keys_rows,keys_cols] = np.nan
                    logging.warn("skip covariance matrix for ({}, {})".format(keys_rows, keys_cols))
                    continue
                ex = covs[keys_rows,keys_rows].index.values
                ey = covs[keys_cols,keys_cols].columns.values
                covs[keys_rows,keys_cols] = cov.change_grid(ex, ey)
        covs.dropna(inplace=True)
        if covs.empty:
            logging.warn("covariance matrix is empty")
            return pd.DataFrame()
        # Create index for global matrix
        rows_levels = covs.index.levels[0]
        indexlist = [(*keys,e) for keys in rows_levels for e in covs[(keys,keys)].index.values]
        index = pd.MultiIndex.from_tuples(indexlist, names=cls.labels)
        # Create global matrix
        matrix = np.zeros((len(index),len(index)))
        for (keys_rows,keys_cols), cov in covs.iteritems():
            ix = index.get_loc(keys_rows)
            ix1 = index.get_loc(keys_cols)
            matrix[ix.start:ix.stop,ix1.start:ix1.stop] = cov
            if keys_rows != keys_cols:
                matrix[ix1.start:ix1.stop,ix.start:ix.stop] = cov.T
        return cls(matrix, index=index, columns=index)