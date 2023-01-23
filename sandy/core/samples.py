import numpy as np
import pandas as pd
import scipy

import sandy

__author__ = "Luca Fiorito"
__all__ = [
        "Samples",
        ]


class Samples():
    """
    Container for samples.
    
    Attributes
    ----------
    condition_number
        
    data
        

    Methods
    -------
    filter_by
        
    from_csv
       
    """

    _columnsname = "SMP"

    def __init__(self, df, *args, **kwargs):
        self.data = pd.DataFrame(df, *args, dtype=float, **kwargs)

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
        self._data = data.rename_axis(self.__class__._columnsname, axis=1)

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

    def get_mean(self):
        return self.data.mean(axis=1).rename("MEAN")

    def get_cov(self):
        return self.data.T.cov()

    def get_std(self):
        return self.data.std(axis=1).rename("STD")

    def get_rstd(self):
        return (self.get_std() / self.get_mean()).rename("RSTD")

    def iterate_xs_samples(self):
        levels = sandy.Xs._columnsnames
        df = self.data.unstack(level=levels)
        for n, p in df.groupby(axis=1, level=self._columnsname):
            yield n, p.droplevel(self._columnsname, axis=1)

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
    
    def test_shapiro(self, size=None, pdf="normal"):
        """
        

        Parameters
        ----------
        size : TYPE, optional
            DESCRIPTION. The default is None.
        pdf : TYPE, optional
            DESCRIPTION. The default is "normal".

        Returns
        -------
        TYPE
            DESCRIPTION.

        Examples
        --------
        Generate 5000 xs samples normally and log-normally distributed
        >>> tape = sandy.get_endf6_file("jeff_33", "xs", 10010)
        >>> njoy_kws = dict(err=1, errorr33_kws=dict(mt=102))
        >>> nsmp = 5000
        >>> smp_norm = tape.get_perturbations(nsmp, njoy_kws=njoy_kws, smp_kws=dict(seed33=1, pdf="normal"))[33]
        >>> smp_lognorm = tape.get_perturbations(nsmp, njoy_kws=njoy_kws, smp_kws=dict(seed33=1, pdf="lognormal"))[33]

        The Shapiro-Wilks test proves right for the normal samples and the normal distribution.
        >>> stat_norm = []
        >>> stat_lognorm = []
        >>> for nsmp in [10, 50, 100, 500, 1000, 5000]:
        ...    df = smp_norm.test_shapiro(pdf="normal", size=nsmp)
        ...    idx = df.statistic.idxmin()
        ...    stat_norm.append(df.loc[idx].rename(nsmp))
        ...
        ...    df = smp_norm.test_shapiro(pdf="lognormal", size=nsmp)
        ...    idx = df.statistic.idxmin()
        ...    stat_lognorm.append(df.loc[idx].rename(nsmp))
        >>>
        >>> opts = dict(left_index=True, right_index=True, suffixes=("_norm", "_lognorm"))
        >>> df = pd.DataFrame(stat_norm).merge(pd.DataFrame(stat_lognorm), **opts).rename_axis("# SMP")
        >>> print(df)
               statistic_norm  pvalue_norm  statistic_lognorm  pvalue_lognorm
        # SMP                                                                
        10        9.19610e-01  3.15408e-01        9.40333e-01     6.43714e-41
        50        9.44342e-01  1.84070e-02        9.40333e-01     6.43714e-41
        100       9.78776e-01  1.03195e-01        9.40333e-01     6.43714e-41
        500       9.85038e-01  5.02533e-05        9.40333e-01     6.43714e-41
        1000      9.89657e-01  1.68679e-06        9.40333e-01     6.43714e-41
        5000      9.90618e-01  1.03590e-17        9.40333e-01     6.43714e-41

        The Shapiro-Wilks test proves right for the lognormal samples and the lognormal distribution.
        >>> stat_norm = []
        >>> stat_lognorm = []
        >>> for nsmp in [10, 50, 100, 500, 1000, 5000]:
        ...    df = smp_lognorm.test_shapiro(pdf="normal", size=nsmp)
        ...    idx = df.statistic.idxmin()
        ...    stat_norm.append(df.loc[idx].rename(nsmp))
        ...
        ...    df = smp_lognorm.test_shapiro(pdf="lognormal", size=nsmp)
        ...    idx = df.statistic.idxmin()
        ...    stat_lognorm.append(df.loc[idx].rename(nsmp))
        >>>
        >>> opts = dict(left_index=True, right_index=True, suffixes=("_norm", "_lognorm"))
        >>> df = pd.DataFrame(stat_norm).merge(pd.DataFrame(stat_lognorm), **opts).rename_axis("# SMP")
        >>> print(df)
               statistic_norm  pvalue_norm  statistic_lognorm  pvalue_lognorm
        # SMP                                                                
        10        6.97781e-01  4.38032e-04        9.99319e-01     5.41511e-02
        50        8.19990e-01  2.17409e-06        9.99319e-01     5.41511e-02
        100       7.90146e-01  1.08222e-10        9.99319e-01     5.41511e-02
        500       8.99920e-01  1.35794e-17        9.99319e-01     5.41511e-02
        1000      9.14015e-01  2.26424e-23        9.99319e-01     5.41511e-02
        5000      9.02778e-01  0.00000e+00        9.99319e-01     5.41511e-02
        """
        size_ = size or self.data.shape[1]
        names = ["statistic", "pvalue"]

        data = self.data.loc[:, :size_]
        if pdf.lower() == "lognormal":
            data = np.log(self.data)

        df = pd.DataFrame({idx: scipy.stats.shapiro(row) for idx, row in data.iterrows()}, index=names).T
        return df.rename_axis(data.index.names)
