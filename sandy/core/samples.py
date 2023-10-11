import numpy as np
import pandas as pd

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
        self._data.columns = self._data.columns.astype(int)

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

    def iterate_xs_samples(self):
        """
        Iterate samples one by one and shape them as a :func:`sandy.Xs`
        dataframe, but with mutligroup structure.
        This output should be passed to :func:`sandy.Xs._perturb`.
        The function is called by :func:`sandy.Endf6.apply_perturbations`

        Yields
        ------
        n : `int`
            .
        s : `pd.DataFrame`
            dataframe of perturbation coefficients with:
                
                - columns: `pd.MultiIndex` with levels `"MAT"` and `"MT"`
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
        >>> smps2 = endf6.get_perturbations(1, njoy_kws=dict(err=1, chi=False, mubar=False, nubar=False, errorr33_kws=dict(mt=2)))[33]

        Copy samples each time to a redundant or partial MT
        >>> smps1 = sandy.Samples(smps2.data.reset_index().assign(MT=1).set_index(["MAT", "MT", "E"]))
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
        >>> assert next(smps51.iterate_xs_samples())[1].columns.equals(expected)

        >>> expected = pd.MultiIndex.from_product([[125], [4] + list(sandy.redundant_xs[4])], names=["MAT", "MT"])
        >>> assert next(smps4.iterate_xs_samples())[1].columns.equals(expected)

        >>> expected = pd.MultiIndex.from_product([[125], [1] + list(sandy.redundant_xs[1])], names=["MAT", "MT"])
        >>> assert next(smps1.iterate_xs_samples())[1].columns.equals(expected)

        >>> expected = pd.MultiIndex.from_product([[125], [3] + list(sandy.redundant_xs[3])], names=["MAT", "MT"])
        >>> assert next(smps3.iterate_xs_samples())[1].columns.equals(expected)

        >>> expected = pd.MultiIndex.from_product([[125], [1] + list(sandy.redundant_xs[1])], names=["MAT", "MT"])
        >>> assert next(smps1.iterate_xs_samples())[1].columns.equals(expected)

        >>> expected = pd.MultiIndex.from_product([[125], [18] + list(sandy.redundant_xs[18])], names=["MAT", "MT"])
        >>> assert next(smps18.iterate_xs_samples())[1].columns.equals(expected)

        >>> expected = pd.MultiIndex.from_product([[125], [27] + list(sandy.redundant_xs[27])], names=["MAT", "MT"])
        >>> assert next(smps27.iterate_xs_samples())[1].columns.equals(expected)

        >>> expected = pd.MultiIndex.from_product([[125], [101] + list(sandy.redundant_xs[101])], names=["MAT", "MT"])
        >>> assert next(smps101.iterate_xs_samples())[1].columns.equals(expected)

        >>> expected = pd.MultiIndex.from_product([[125], [452] + list(sandy.redundant_xs[452])], names=["MAT", "MT"])
        >>> assert next(smps452.iterate_xs_samples())[1].columns.equals(expected)


        In this example the original covariance contains data for MT=1 and MT=51.
        >>> endf6 = sandy.get_endf6_file('jeff_33', 'xs', 942400)
        >>> smps = endf6.get_perturbations(1, njoy_kws=dict(err=1, chi=False, mubar=False, nubar=False, errorr33_kws=dict(mt=[1, 51])))[33]

        Then, since MT=1 is redundant, samples are passed to its partial components (MT=2 and MT=3).
        >>> expected = pd.MultiIndex.from_product([[9440], [1, 51] + list(sandy.redundant_xs[1])], names=["MAT", "MT"])
        >>> assert next(smps.iterate_xs_samples())[1].columns.equals(expected)
        
        If case one of the partial components already has samples, i.e., MT=2...
        >>> endf6 = sandy.get_endf6_file('jeff_33', 'xs', 942400)
        >>> smps = endf6.get_perturbations(1, njoy_kws=dict(err=1, chi=False, mubar=False, nubar=False, errorr33_kws=dict(mt=[1, 2, 51])))[33]

        Then the MT=1 samples are not passed to the partial components, which 
        in this case it means that MT=2 is not changed and MT=3 is not created.
        >>> expected = pd.MultiIndex.from_product([[9440], [1, 2, 51]], names=["MAT", "MT"])
        >>> assert next(smps.iterate_xs_samples())[1].columns.equals(expected)
        """
        levels = sandy.Xs._columnsnames
        df = self.data.unstack(level=levels)
        
        # -- Iterate over samples
        for n, p in df.groupby(axis=1, level=self._columnsname):
            s = p.droplevel(self._columnsname, axis=1)
            adds = []
            for mat in s.columns.get_level_values("MAT").unique():
                
                # -- Iterate redundant xs (from MT107 to MT1)
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

    def from_excel(file, beg=None, end=None):
        """
        Read perturbation coefficients (for nubar and xs) from excel file.
        The file format is compatible with what written in
        :obj: `~sandy.cov.CategoryCov.sampling`
        
        Parameters
        ----------
        file : `str`
            excel file name.
        beg : `int`
            first sample to consider. Default is first in dataframe.
        end : `int`
            last samples to conisder. Default is last in dataframe.

        Returns
        -------
        smp : :obj: `Samples`
            Samples dataframe.
        """
        df = pd.read_excel(file, sheet_name="SMP")
        
        # -- Everything before ELEFT is part of the MultiIndex
        loc = df.columns.get_loc("ELEFT")
        idx = df.iloc[:, :loc]

        # -- ELEFT and ERIGHT are combined into an IntervalIndex (multigroup)
        idx.insert(loc, "E", pd.IntervalIndex.from_arrays(df.ELEFT, df.ERIGHT))

        idx = pd.MultiIndex.from_frame(idx)

        df = df.iloc[:, loc+2:].loc[:, beg:end].reset_index(drop=True)
        df.index = idx
        smp = sandy.Samples(df)
        return smp
