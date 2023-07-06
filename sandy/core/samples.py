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
    data
        Dataframe of samples.

    Methods
    -------
    get_condition_number
        Return condition number of samples.
    get_corr
        Return correlation matrix of samples.
    get_cov
        Return covariance matrix of samples.
    get_mean
        Return mean vector of samples.
    get_std
        Return standard deviation vector of samples.
    get_rstd
        Return relative standard deviation vector of samples.   
    iterate_xs_samples
        Generator that iterates over each sample (in the form of :func:`sandy.Xs`).
    test_shapiro
        Perform the Shapiro-Wilk test for normality on the samples.
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

    def get_condition_number(self):
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

    def get_corr(self):
        return self.data.T.corr()

    def get_cov(self):
        return self.data.T.cov()

    def get_std(self):
        return self.data.std(axis=1).rename("STD")

    def get_rstd(self):
        return (self.get_std() / self.get_mean()).rename("RSTD")

    def iterate_samples(self, kind):
        """
        Iterate samples one by one and shape them as a :func:`sandy.Xs`
        dataframe, but with mutligroup structure.
        This output should be passed to :func:`sandy.Xs._perturb`.
        The function is called by :func:`sandy.Endf6.apply_perturbations`
        The function is valid also for the `sandy.Edistr` instance but
        The kind must be specified.

        Yields
        ------
        n : `int`
        kind: `str`
            - "xs" for cross setions and nubar
            - "chi" for energy distributions
        s : `pd.DataFrame`
            dataframe of perturbation coefficients with:
                
                - columns: `pd.MultiIndex` depending on the kind
                - index: `pd.IntervalIndex` with multigroup structure

        Notes
        -----
        If samples refer to redundant MT number, the same identical samples
        are passed one level down to the partial MT components.
        For instance:
            - MT=4 samples will be assigned to MT=50-91
            - MT=1 samples will be assigned to MT=2 and MT=3
            - MT=18 samples will be assigned to MT=19-21 and MT=38
        
        ..important:: Assigning samples from redundant MT number to partial
                      components only applies if the partial components do not
                      have their own samples, and it only goes one level deep.

        Examples
        --------
        Get samples fot MT=1
        >>> endf6 = sandy.get_endf6_file('jeff_33', 'xs', 10010)
        >>> smps1 = endf6.get_perturbations(1, njoy_kws=dict(err=1, chi=False, mubar=False, nubar=False, errorr33_kws=dict(mt=1)))[33]

        Copy samples each time to a redundant or partial MT
        >>> smps3 = sandy.Samples(smps1.data.reset_index().assign(MT=3).set_index(["MAT", "MT", "E"]))
        >>> smps18 = sandy.Samples(smps1.data.reset_index().assign(MT=18).set_index(["MAT", "MT", "E"]))
        >>> smps19 = sandy.Samples(smps1.data.reset_index().assign(MT=19).set_index(["MAT", "MT", "E"]))
        >>> smps27 = sandy.Samples(smps1.data.reset_index().assign(MT=27).set_index(["MAT", "MT", "E"]))
        >>> smps4 = sandy.Samples(smps1.data.reset_index().assign(MT=4).set_index(["MAT", "MT", "E"]))
        >>> smps51 = sandy.Samples(smps1.data.reset_index().assign(MT=51).set_index(["MAT", "MT", "E"]))
        >>> smps101 = sandy.Samples(smps1.data.reset_index().assign(MT=101).set_index(["MAT", "MT", "E"]))
        >>> smps452 = sandy.Samples(smps1.data.reset_index().assign(MT=452).set_index(["MAT", "MT", "E"]))

        Check that samples are passed correctly to daughter MTs (only one level deep)
        >>> expected = pd.MultiIndex.from_product([[125], [51]], names=["MAT", "MT"])
        >>> assert next(smps51.iterate_samples("xs"))[1].columns.equals(expected)

        >>> expected = pd.MultiIndex.from_product([[125], [4] + list(sandy.redundant_xs[4])], names=["MAT", "MT"])
        >>> assert next(smps4.iterate_samples("xs"))[1].columns.equals(expected)

        >>> expected = pd.MultiIndex.from_product([[125], [1] + list(sandy.redundant_xs[1])], names=["MAT", "MT"])
        >>> assert next(smps1.iterate_samples("xs"))[1].columns.equals(expected)

        >>> expected = pd.MultiIndex.from_product([[125], [3] + list(sandy.redundant_xs[3])], names=["MAT", "MT"])
        >>> assert next(smps3.iterate_samples("xs"))[1].columns.equals(expected)

        >>> expected = pd.MultiIndex.from_product([[125], [1] + list(sandy.redundant_xs[1])], names=["MAT", "MT"])
        >>> assert next(smps1.iterate_samples("xs"))[1].columns.equals(expected)

        >>> expected = pd.MultiIndex.from_product([[125], [18] + list(sandy.redundant_xs[18])], names=["MAT", "MT"])
        >>> assert next(smps18.iterate_samples("xs"))[1].columns.equals(expected)

        >>> expected = pd.MultiIndex.from_product([[125], [27] + list(sandy.redundant_xs[27])], names=["MAT", "MT"])
        >>> assert next(smps27.iterate_samples("xs"))[1].columns.equals(expected)

        >>> expected = pd.MultiIndex.from_product([[125], [101] + list(sandy.redundant_xs[101])], names=["MAT", "MT"])
        >>> assert next(smps101.iterate_samples("xs"))[1].columns.equals(expected)

        >>> expected = pd.MultiIndex.from_product([[125], [452] + list(sandy.redundant_xs[452])], names=["MAT", "MT"])
        >>> assert next(smps452.iterate_samples("xs"))[1].columns.equals(expected)


        In this example the original covariance contains data for MT=1 and MT=51.
        >>> endf6 = sandy.get_endf6_file('jeff_33', 'xs', 942400)
        >>> smps = endf6.get_perturbations(1, njoy_kws=dict(err=1, chi=False, mubar=False, nubar=False, errorr33_kws=dict(mt=[1, 51])))[33]

        Then, since MT=1 is redundant, samples are passed to its partial components (MT=2 and MT=3).
        >>> expected = pd.MultiIndex.from_product([[9440], [1, 51] + list(sandy.redundant_xs[1])], names=["MAT", "MT"])
        >>> assert next(smps.iterate_samples("xs"))[1].columns.equals(expected)
        
        If case one of the partial components already has samples, i.e., MT=2...
        >>> endf6 = sandy.get_endf6_file('jeff_33', 'xs', 942400)
        >>> smps = endf6.get_perturbations(1, njoy_kws=dict(err=1, chi=False, mubar=False, nubar=False, errorr33_kws=dict(mt=[1, 2, 51])))[33]

        Then the MT=1 samples are not passed to the partial components, which 
        in this case it means that MT=2 is not changed and MT=3 is not created.
        >>> expected = pd.MultiIndex.from_product([[9440], [1, 2, 51]], names=["MAT", "MT"])
        >>> assert next(smps.iterate_samples("xs"))[1].columns.equals(expected)
        """
        if kind == "xs":
            levels = sandy.Xs._columnsnames
        elif kind == "chi":
            levels = sandy.Edistr._columnsnames

        df = self.data.unstack(level=levels)
        
        # iterate over samples
        for n, p in df.groupby(axis=1, level=self._columnsname):
            s = p.droplevel(self._columnsname, axis=1)
            
            if kind == "xs":
                adds = []
                for mat in s.columns.get_level_values("MAT").unique():
                    
                    # sort from MT107 to MT1
                    for k, v in sandy.redundant_xs.items():
                        if not (mat, k) in s.columns:
                            continue
                        daughters = pd.MultiIndex.from_product([[mat], v], names=["MAT", "MT"])
                        # Only give perturbation for redundant xs to daughters if no perturbation
                        # for partial cross section is found
                        if s.columns.intersection(daughters).empty:
                            
                            # This goes only 1 level deep.
                            # Then, MT=1 perturbations will be given to MT=2 and MT=3
                            # without descending any further
                            add = pd.DataFrame(
                                np.tile(s[(mat, k)].values, (daughters.size, 1)).T,
                                index=s.index,
                                columns=daughters,
                                )
                            adds.append(add)
                s = pd.concat([s, *adds], axis=1)
            yield n, s

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
        Perform the Shapiro-Wilk test for normality on the samples.
        The test can be performed also for a lognormal distribution by testing
        for normality the logarithm of the samples.
        
        The Shapiro-Wilk test tests the null hypothesis that the data was
        drawn from a normal distribution.

        Parameters
        ----------
        size : `int`, optional
            number of samples (starting from the first) that need to be
            considered for the test. The default is `None`, i.e., all samples.
        pdf : `str`, optional
            the pdf used to test the samples. Either `"normal"` or
            `"lognormal"`. The default is "normal".

        Returns
        -------
        pd.DataFrame
            Dataframe with Shapriro-Wilk results (statistic and pvalue) for
            each variable considered in the :func:`~Samples` instance.

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
        # >>> print(test(smp_norm))
        PDF   lognormal      normal
        # SMP                      
        10       reject      reject
        50       reject      reject
        100      reject      reject
        500      reject      reject
        1000     reject      reject
        5000     reject      reject

        The Shapiro-Wilks test proves right for the lognormal samples and the lognormal distribution.
        # >>> print(test(smp_lognorm))
        PDF    lognormal  normal
        # SMP                   
        10     lognormal  reject
        50     lognormal  reject
        100    lognormal  reject
        500    lognormal  reject
        1000   lognormal  reject
        5000   lognormal  reject

        The Shapiro-Wilks gives too low p-values for the uniform samples.
        # >>> print(test(smp_uniform))
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

        data = self.data.iloc[:, :size_]
        if pdf.lower() == "lognormal":
            data = np.log(self.data)

        df = pd.DataFrame({idx: scipy.stats.shapiro(row) for idx, row in data.iterrows()}, index=names).T
        return df.rename_axis(data.index.names)
