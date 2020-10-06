import pdb
import logging

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import statsmodels.api as sm

import sandy

__author__ = "Luca Fiorito"
__all__ = [
        "Samples",
        ]


class Samples():

    def __init__(self, df):
        self.data = df

    def __repr__(self):
        return self.data.__repr__()

    @property
    def data(self):
        """
        Dataframe of samples.

        Attributes
        ----------
        index : `pandas.Index` or `pandas.MultiIndex`
            indices
        columns : `pandas.Index`
            samples numbering starting from zero
        values : `numpy.array`
            sample values as `float`

        Returns
        -------
        `pandas.DataFrame`
            tabulated samples

        Raises
        ------
        `sandy.Error`
            if 'data' is not a 'pandas.DataFrame' or `numpy.ndarray`
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

    def std_convergence(self):
        smp = self.data
        rng = range(2, smp.shape[0])
        foo = lambda x: smp.loc[:x].std()
        return pd.DataFrame(map(foo, rng), index=rng)

    def mean_convergence(self):
        smp = self.data
        rng = range(1, smp.shape[0])
        foo = lambda x: smp.loc[:x].mean()
        return pd.DataFrame(map(foo, rng), index=rng)

    def truncate(left=0, left_value=0, right=2, right_value=2, inplace=True):
        pass

    @property
    def condition_number(self):
        """
        Return condition number of samples.

        Notes
        -----
        ..note:: the condition number can help assess multicollinearity.
        """
        # The first step is to normalize the independent variables to have
        # unit length
        X = self.data.T
        norm_x = X.values
        for i, name in enumerate(X):
            norm_x[:, i] = X[name] / np.linalg.norm(X[name])
        norm_xtx = np.dot(norm_x.T, norm_x)
        # Then, we take the square root of the ratio of the biggest to the
        # smallest eigen values
        eigs = np.linalg.eigvals(norm_xtx)
        return np.sqrt(eigs.max() / eigs.min())

    def heatmap(self, vmin=-1, vmax=1, cmap="bwr", **kwargs):
        corr = np.corrcoef(self.data)
        return sns.heatmap(corr, vmin=vmin, vmax=vmax, cmap=cmap, **kwargs)

    def sm_ols(self, Y, normalized=False, intercept=False):
        X = self.data.T
        if normalized:
            X = X.divide(X.mean()).fillna(0)
            Y = Y.divide(Y.mean()).fillna(0)
        if intercept:
            X = sm.add_constant(X)
        model = sm.OLS(Y.values, X.values)
        return model.fit()

    def sensitivities(self, Y, **kwargs):
        res = self.sm_ols(Y, **kwargs)
        return pd.Series(res.params)

    def ics(self, Y, **kwargs):
        return self.sensitivities(Y, **kwargs).sum()

    @classmethod
    def from_csv(cls, file, **kwargs):
        df = pd.read_csv(file, **kwargs)
        return cls(df)
