import logging
import io

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

np.random.seed(1)
minimal_testcase = np.random.randn(4, 3)


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
                    )
                out.data.index.names = ["MAT", "MT", "E"]
                return out
            else:
                del kw[key]
        out = func(*args, **kw)
        return out
    return inner


class Samples():
    """
    Attributes
    ----------
    condition_number

    data


    Methods
    -------
    filter_by

    from_csv

    regression_coefficients

    sm_ols

    """

    def __init__(self, df):
        self.data = pd.DataFrame(df, dtype=float)

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
            samples numbering
        values : `numpy.array`
            sample values as `float`

        Returns
        -------
        `pandas.DataFrame`
            tabulated samples
        """
        return self._data

    @data.setter
    def data(self, data):
        self._data = data
        self._data.index.name = 'SMP'

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
        X = self.data.T.copy()
        norm_x = X.values
        for i, name in enumerate(X):
            norm_x[:, i] = X[name] / np.linalg.norm(X[name])
        norm_xtx = np.dot(norm_x.T, norm_x)
        # Then, we take the square root of the ratio of the biggest to the
        # smallest eigen values
        eigs = np.linalg.eigvals(norm_xtx)
        return np.sqrt(eigs.max() / eigs.min())

    @property
    def mean(self):
        return self.data.mean(axis=1).rename("MEAN")

    @property
    def rstd(self):
        return (self.std / self.mean).rename("RSTD")

    @property
    def std(self):
        return self.data.std(axis=1).rename("STD")

    def filter_by(self, key, value):
        """
        Apply condition to source data and return filtered results.

        Parameters
        ----------
        `key` : `str`
            any label present in the columns of `data`
        `value` : `int` or `float`
            value used as filtering condition

        Returns
        -------
        `sandy.Samples`
            filtered dataframe of samples

        Raises
        ------
        `sandy.Error`
            if applied filter returned empty dataframe

        Notes
        -----
        .. note:: The primary function of this method is to make sure that
                  the filtered dataframe is still returned as a `Samples`
                  object.
        """
        condition = self.data.index.get_level_values(key) == value
        out = self.data.copy()[condition]
        if out.empty:
            raise sandy.Error("applied filter returned empty dataframe")
        return self.__class__(out)

    def _std_convergence(self):
        smp = self.data
        rng = range(2, smp.shape[0])
        foo = lambda x: smp.loc[:x].std()
        return pd.DataFrame(map(foo, rng), index=rng)

    def _mean_convergence(self):
        smp = self.data
        rng = range(1, smp.shape[0])
        foo = lambda x: smp.loc[:x].mean()
        return pd.DataFrame(map(foo, rng), index=rng)

    def _heatmap(self, vmin=-1, vmax=1, cmap="bwr", **kwargs):
        corr = np.corrcoef(self.data)
        return sns.heatmap(corr, vmin=vmin, vmax=vmax, cmap=cmap, **kwargs)

    def sm_ols(self, Y, normalized=False, intercept=False):
        X = self.data.T.copy()
        NX, MX = X.shape
        NY = Y.size
        N = min(NX, NY)
        if NX != NY:
            print(f"X and Y have different size, fit only first {N} samples")
        if normalized:
            X = X.divide(X.mean()).fillna(0)
            Y = Y.divide(Y.mean()).fillna(0)
        if intercept:
            X = sm.add_constant(X)
        model = sm.OLS(Y.iloc[:N].values, X[:N].values)
        out = model.fit()
        return out

    def regression_coefficients(self, Y, **kwargs):
        """
        Calculate regression coefficients from OLS model given an output
        population.

        Parameters
        ----------
        Y : `pandas.Series`
            tabulated output population
        kwargs : keyword arguments, optional
            arguments to pass to method `sm_ols`

        Returns
        -------
        `pandas.DataFrame`
            Dataframe with regression coefficients and standard errors.
        """
        X = self.data
        MX, NX = X.shape
        index = X.index
        res = self.sm_ols(Y, **kwargs)
        params = res.params
        bse = res.bse
        start_at = 0 if params.size == MX else 1
        coeff = pd.DataFrame({
            "coeff": params[start_at:],
            "stderr": bse[start_at:],
            }, index=index)
        return coeff

    @classmethod
    @cov33csv
    def from_csv(cls, file, **kwargs):
        """
        Read samples from csv file,

        Parameters
        ----------
        file : `str`
            csv file.
        **kwargs : `dict`
            keyword options for `pandas.read_csv`.

        Returns
        -------
        `sandy.Samples`
            samples into a sandy object.

        Examples
        --------
        >>> csv = minimal_testcase.to_string()
        >>> sandy.Samples.from_csv(io.StringIO(csv), sep="\s+")
                     0            1            2
        0  1.62435e+00 -6.11756e-01 -5.28172e-01
        1 -1.07297e+00  8.65408e-01 -2.30154e+00
        2  1.74481e+00 -7.61207e-01  3.19039e-01
        3 -2.49370e-01  1.46211e+00 -2.06014e+00

        >>> index = pd.MultiIndex.from_product(
        ...    [[9437], [102], [1e-5, 1e-1, 1e1, 1e6]]
        ... )
        >>> df = minimal_testcase.copy()
        >>> df.index = index
        >>> csv = df.to_csv()
        >>> sandy.Samples.from_csv(io.StringIO(csv), sep="\s+", cov33csv=True)
                                        0            1            2
        MAT  MT  E
        9437 102 1.00000e-05  1.62435e+00 -6.11756e-01 -5.28172e-01
                 1.00000e-01 -1.07297e+00  8.65408e-01 -2.30154e+00
                 1.00000e+01  1.74481e+00 -7.61207e-01  3.19039e-01
                 1.00000e+06 -2.49370e-01  1.46211e+00 -2.06014e+00
        """
        df = pd.read_csv(file, **kwargs)
        return cls(df)

    def to_pert(self, smp=None):
        """
        Samples in the format of `Sandy.Pert`. If the sample number is not
        entered, a `pd.DataFrame` is created with a multiindex, the first level
        being the sample number and the second the energy grid.

        Parameters
        ----------
        smp : `int`, optional
            Number of the sample. The default is None.

        Returns
        -------
        `sandy.Pert` or `pd.DataFrame`
            Samples in `sandy.Pert` format.

        Examples
        --------
        >>> samples = pd.DataFrame([[0.5, 1.5], [0.75, 1.25]] , index=[0, 1])
        >>> samples.columns = pd.MultiIndex.from_arrays([[125, 125], [2, 2], [1e-05, 19970500.0]], names=['MAT', 'MT', 'E'])
        >>> sandy.Samples(samples).to_pert()
        	MAT	        125
            MT	        2
        SMP	          E
          0	1.00000e-05	5.00000e-01
            1.99705e+07	1.50000e+00
          1	1.00000e-05	7.50000e-01
            1.99705e+07	1.25000e+00

        >>> sandy.Samples(samples).to_pert(smp=0)
        MAT                         125
        MT                            2
        ENERGY
        (0.0, 1e-05]        5.00000e-01
        (1e-05, 19970500.0] 1.50000e+00

        >>> samples = pd.DataFrame([[0.5, 1.5], [0.75, 1.25]] , index=[0, 1])
        >>> samples.columns = pd.MultiIndex.from_arrays([[125, 125], [2, 2], [1.00000e-05, 1.00000e-05], [1.00000e-05, 	1.00000e-05], [1e-05, 19970500.0]], names=['MAT', 'MT', 'ELO', 'EHI', 'E'])
        >>> sandy.Samples(samples).to_pert()
        	MAT	        125
            MT	        2
            ELO	        1.00000e-05
            EHI	        1.00000e-05
        SMP	          E
          0	1.00000e-05	5.00000e-01
            1.99705e+07	1.50000e+00
          1	1.00000e-05	7.50000e-01
            1.99705e+07	1.25000e+00

        >>> sandy.Samples(samples).to_pert(smp=0)
        MAT                         125
        MT                            2
        ELO                 1.00000e-05
        EHI                 1.00000e-05
        ENERGY
        (0.0, 1e-05]        5.00000e-01
        (1e-05, 19970500.0] 1.50000e+00

        >>> samples = pd.DataFrame([[0.5, 1.5], [0.75, 1.25]] , index=[0, 1])
        >>> samples.columns = pd.MultiIndex.from_arrays([[125, 125], [2, 2], [1, 1], [1e-05, 19970500.0]], names=['MAT', 'MT', 'L', 'E'])
        >>> sandy.Samples(samples).to_pert()
            MAT	        125
            MT	        2
            L	        1
        SMP	          E
          0	1.00000e-05	5.00000e-01
            1.99705e+07	1.50000e+00
          1	1.00000e-05	7.50000e-01
            1.99705e+07	1.25000e+00

        >>> sandy.Samples(samples).to_pert(smp=0)
        MAT                         125
        MT                            2
        L                             1
        ENERGY
        (0.0, 1e-05]        5.00000e-01
        (1e-05, 19970500.0] 1.50000e+00
        """
        data = self.data.copy()
        levels = list(np.arange(data.columns.nlevels))
        if smp is not None:
            if isinstance(smp, int) and smp in data.index:
                pert = data.loc[smp].to_frame().unstack(level=levels[0:-1])
                pert.columns = pert.columns.droplevel(level=None)
                pert = sandy.Pert(pert)
            else:
                print(f"{smp} is not correct or do not exist")
        else:
            def foo(df):
                df = df.stack()
                df.index = df.index.droplevel(level='SMP')
                return df
            pert = data.groupby('SMP').apply(foo).fillna(0)
        return pert
