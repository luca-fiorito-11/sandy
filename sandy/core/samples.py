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
        Generate 5000 xs samples normally, log-normally and uniform distributed
        >>> tape = sandy.get_endf6_file("jeff_33", "xs", 10010)
        >>> njoy_kws = dict(err=1, errorr33_kws=dict(mt=102))
        >>> nsmp = 5000
        >>> seed = 5
        >>>
        >>> smp_norm = tape.get_perturbations(nsmp, njoy_kws=njoy_kws, smp_kws=dict(seed33=seed, pdf="normal"))[33]
        >>> smp_lognorm = tape.get_perturbations(nsmp, njoy_kws=njoy_kws, smp_kws=dict(seed33=seed, pdf="lognormal"))[33]
        >>> smp_uniform = tape.get_perturbations(nsmp, njoy_kws=njoy_kws, smp_kws=dict(seed33=seed, pdf="uniform"))[33]

        In this example we defined the following arbitrary convergence criteria:
            - if the p value is larger than 0.05 we fail to reject the null-hypothesis and we accept the results
            - if the first condition is accepted, we confirm the pdf if the statistics is larger than 0.95
        >>> threshold = 0.95
        >>> pthreshold = 0.05
        >>> def test(smps):
        ...     data = []
        ...     for n in [10, 50, 100, 500, 1000, 5000]:
        ...         for pdf in ("normal", "lognormal"):
        ...             df = smps.test_shapiro(pdf=pdf, size=n)
        ...             idx = df.statistic.idxmin()
        ...             w = df.loc[idx]
        ...             t = "reject" if w.pvalue < pthreshold else (pdf if w.statistic > threshold else "reject")
        ...             data.append({"PDF": pdf, "test":t, "# SMP": n})
        ...     df = pd.DataFrame(data).pivot_table(index="# SMP", columns="PDF", values="test", aggfunc=lambda x: ' '.join(x))
        ...     return df

        The Shapiro-Wilks test proves wrong the normal samples because of the tail truncation.
        >>> print(test(smp_norm))
        PDF   lognormal      normal
        # SMP                      
        10       reject  not-normal
        50       reject      reject
        100      reject      reject
        500      reject      reject
        1000     reject      reject
        5000     reject      reject

        The Shapiro-Wilks test proves right for the lognormal samples and the lognormal distribution.
        >>> print(test(smp_lognorm))
        PDF    lognormal  normal
        # SMP                   
        10     lognormal  reject
        50     lognormal  reject
        100    lognormal  reject
        500    lognormal  reject
        1000   lognormal  reject
        5000   lognormal  reject

        The Shapiro-Wilks gives too low p-values for the uniform samples.
        >>> print(test(smp_uniform))
        PDF   lognormal  normal
        # SMP                  
        10       reject  reject
        50       reject  reject
        100      reject  reject
        500      reject  reject
        1000     reject  reject
        5000     reject  reject
        """
        size_ = size or self.data.shape[1]
        names = ["statistic", "pvalue"]

        data = self.data.loc[:, :size_]
        if pdf.lower() == "lognormal":
            data = np.log(self.data)

        df = pd.DataFrame({idx: scipy.stats.shapiro(row) for idx, row in data.iterrows()}, index=names).T
        return df.rename_axis(data.index.names)
