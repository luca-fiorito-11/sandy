"""
This module contains all classes and functions dedicated to handling PFNS
data and, more in general, any tabulated energy distribution provided in
MF5 sections.
"""
import logging


import pandas as pd
import numpy as np

import sandy

__author__ = "Luca Fiorito"
__all__ = [
    "Edistr",
    ]

minimal_edistrtest = pd.DataFrame(
    [[0.4, 0],
     [0, 0.2],
     [0, 0.7],
     [0, 0.1],
     [0.6, 0]],
    index=[1e-5, 1e-4, 1, 1e7, 2e7],
    columns=pd.MultiIndex.from_product([[9437],[18],[0],[1e0, 2e0]])
    )


class Edistr():
    """
    Object to store tabulate energy distributions.

    Attributes
    ----------
    data : `pandas.DataFrame`
        dataframe of energy distribution data 

    Methods
    -------
    add_energy_point
        add outgoing energy distributions at additional single incident energy
        by interpolation
    add_energy_points
        add outgoing energy distributions at additional incident energies by
        interpolation
    custom_perturbation
        perturb individual outgoing distribution
    filter_by
        apply condition to source data and return filtered results
    from_endf6
        extract energy distributions from `Endf6` instance
    get_integrals
        calculate the integral of each energy distribution
    get_table
        pivot dataframe of tabulated energy spectra
    normalize
        renormalize each outgoing energy distribution to 1
    to_endf6
        Given the `Edistr` instance, the `Endf6` instance is obtained 
        with the relative information on the energy distribution
    """

    _indexname = ("EOUT")
    _columnsnames = ("MAT", "MT", "K", "EIN")

    def __repr__(self):
        return self.data.__repr__()

    def __init__(self, df, **kwargs):
        self.data = pd.DataFrame(df, **kwargs)

    @property
    def data(self):
        """
        Dataframe of energy distribution data with index as outgoing energy
        and columns as pd.MultiIndex with last level as incident energy

        Returns
        -------
        `pandas.DataFrame`
            tabulated energy distribution

        Notes
        -----
        .. note :: tabulated values are assumed to be interpolated linearly

        Examples
        --------
        >>> Edistr(minimal_edistrtest)
        MAT                9437
        MT                   18
        K                     0
        EIN         1.00000e+00 2.00000e+00
        EOUT
        1.00000e-05 4.00000e-01 0.00000e+00
        1.00000e-04 0.00000e+00 2.00000e-01
        1.00000e+00 0.00000e+00 7.00000e-01
        1.00000e+07 0.00000e+00 1.00000e-01
        2.00000e+07 6.00000e-01 0.00000e+00
        """
        return self._data

    @data.setter
    def data(self, data):
        self._data = data.rename_axis(self.__class__._indexname, axis=0)\
                         .rename_axis(self.__class__._columnsnames, axis=1)
        if not data.index.is_monotonic_increasing:
            raise ValueError("energy grid is not monotonically increasing")

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
        `sandy.Edistr`
            filtered dataframe of energy distributions

        Raises
        ------
        `sandy.Error`
            if applied filter returned empty dataframe

        Notes
        -----
        .. note:: The primary function of this method is to make sure that
                  the filtered dataframe is still returned as a `Edistr`
                  object.

        Examples
        --------
        >>> Edistr(minimal_edistrtest).filter_by("EIN", 2)
            MAT  MT  K         EIN        EOUT       VALUE
        0  9437  18  0 2.00000e+00 1.00000e-04 2.00000e-01
        1  9437  18  0 2.00000e+00 1.00000e+00 7.00000e-01
        2  9437  18  0 2.00000e+00 1.00000e+07 1.00000e-01
        """
        condition = self.data[key] == value
        out = self.data.copy()[condition].reset_index(drop=True)
        if out.empty:
            raise sandy.Error("applied filter returned empty dataframe")
        return self.__class__(out)

    def get_table(self, mat, mt, k):
        """
        Pivot dataframe of tabulated energy spectra.

        Parameters
        ----------
        mat : `int`
            MAT number to filter data
        mt : `int`
            MT number to filter data
        k : `int`
            subsection number to filter data

        Returns
        -------
        `pandas.DataFrame`
            tabulated energy spectra with incident energies as index and
            outgoing energies as columns.

        Examples
        --------
        >>> Edistr(minimal_edistrtest).get_table(9437, 18, 0)
        EOUT         1.00000e-05  1.00000e-04  1.00000e+00  1.00000e+07  2.00000e+07
        EIN                                                                         
        1.00000e+00  4.00000e-01  4.00000e-01  4.00000e-01  5.00000e-01  6.00000e-01
        2.00000e+00  0.00000e+00  2.00000e-01  7.00000e-01  1.00000e-01  0.00000e+00
        """
        data = self.filter_by("MAT", mat) \
                   .filter_by("MT", mt) \
                   .filter_by("K", k) \
                   .data
        return data.pivot_table(
            values="VALUE",
            index="EIN",
            columns="EOUT",
            aggfunc="sum",
            ) \
            .interpolate(method="slinear", axis="columns") \
            .fillna(0) \
            .sort_index(axis="index") \
            .sort_index(axis="columns")

    def add_energy_point(self, mat, mt, k, enew):
        """
        Add outgoing energy distribution at one additional incident energy by
        interpolation.

        Parameters
        ----------
        mat : `int`
            MAT number.
        mt : `int`
            MT number.
        k : `int`
            subsection.
        enew : float
            energy point in eV

        Returns
        -------
        `sandy.Edistr`
            energy distribution with an additional additional incoming energy.

        Examples
        --------
        >>> orig = Edistr(minimal_edistrtest)
        >>> new = orig.add_energy_point(9437, 18, 0, 1.5)
        >>> new
            MAT  MT  K         EIN        EOUT       VALUE
        0  9437  18  0 1.00000e+00 1.00000e-05 4.00000e-01
        1  9437  18  0 1.00000e+00 2.00000e+07 6.00000e-01
        2  9437  18  0 1.50000e+00 1.00000e-05 2.00000e-01
        3  9437  18  0 1.50000e+00 1.00000e-04 3.00000e-01
        4  9437  18  0 1.50000e+00 1.00000e+00 5.50000e-01
        5  9437  18  0 1.50000e+00 1.00000e+07 3.00000e-01
        6  9437  18  0 1.50000e+00 2.00000e+07 3.00000e-01
        7  9437  18  0 2.00000e+00 1.00000e-04 2.00000e-01
        8  9437  18  0 2.00000e+00 1.00000e+00 7.00000e-01
        9  9437  18  0 2.00000e+00 1.00000e+07 1.00000e-01

        >>> new.get_table(9437, 18, 0)
        EOUT         1.00000e-05  1.00000e-04  1.00000e+00  1.00000e+07  2.00000e+07
        EIN                                                                         
        1.00000e+00  4.00000e-01  4.00000e-01  4.00000e-01  5.00000e-01  6.00000e-01
        1.50000e+00  2.00000e-01  3.00000e-01  5.50000e-01  3.00000e-01  3.00000e-01
        2.00000e+00  0.00000e+00  2.00000e-01  7.00000e-01  1.00000e-01  0.00000e+00

        If energy point is already present
        >>> copy = orig.add_energy_point(9437, 18, 0, 1)
        >>> pd.testing.assert_frame_equal(orig.data, copy.data)
        """
        data = self.data
        ein = data["EIN"].unique()
        if enew in ein:
            return self.__class__(data)
        mask = ein > enew
        upper = ein[mask][0]
        mask = ein < enew
        lower = ein[mask][-1]
        new_ein = pd.Index([lower, enew, upper], name="EIN")
        distr_lower = data[data.EIN == lower][["EOUT", "VALUE"]].rename({"VALUE": lower}, axis="columns")
        distr_upper = data[data.EIN == upper][["EOUT", "VALUE"]].rename({"VALUE": upper}, axis="columns")
        tab = distr_lower.merge(distr_upper, how="outer", on="EOUT") \
                         .set_index("EOUT") \
                         .interpolate(method="slinear", axis="index") \
                         .fillna(0) \
                         .reindex(new_ein, axis="columns") \
                         .interpolate(method="slinear", axis="columns") \
                         .fillna(0)
        df = tab[enew].rename("VALUE") \
                      .reset_index() \
                      .assign(MAT=mat, MT=mt, K=k, EIN=enew)
        return self.__class__(pd.concat((data, df)))

    def add_energy_points(self, mat, mt, k, elist):
        """
        Add outgoing energy distributions at additional incident energies by
        interpolation.

        Parameters
        ----------
        mat : `int`
            MAT number.
        mt : `int`
            MT number.
        k : `int`
            subsection.
        elist : iterable of `float`
            energy points in eV

        Returns
        -------
        new : `sandy.Edistr`
            energy distribution with an additional additional incoming
            energies.

        Examples
        --------
        >>> orig = Edistr(minimal_edistrtest)
        >>> new = orig.add_energy_points(9437, 18, 0, [1, 1.5, 1.7])
        >>> new
             MAT  MT  K         EIN        EOUT       VALUE
        0   9437  18  0 1.00000e+00 1.00000e-05 4.00000e-01
        1   9437  18  0 1.00000e+00 2.00000e+07 6.00000e-01
        2   9437  18  0 1.50000e+00 1.00000e-05 2.00000e-01
        3   9437  18  0 1.50000e+00 1.00000e-04 3.00000e-01
        4   9437  18  0 1.50000e+00 1.00000e+00 5.50000e-01
        5   9437  18  0 1.50000e+00 1.00000e+07 3.00000e-01
        6   9437  18  0 1.50000e+00 2.00000e+07 3.00000e-01
        7   9437  18  0 1.70000e+00 1.00000e-05 1.20000e-01
        8   9437  18  0 1.70000e+00 1.00000e-04 2.60000e-01
        9   9437  18  0 1.70000e+00 1.00000e+00 6.10000e-01
        10  9437  18  0 1.70000e+00 1.00000e+07 2.20000e-01
        11  9437  18  0 1.70000e+00 2.00000e+07 1.80000e-01
        12  9437  18  0 2.00000e+00 1.00000e-04 2.00000e-01
        13  9437  18  0 2.00000e+00 1.00000e+00 7.00000e-01
        14  9437  18  0 2.00000e+00 1.00000e+07 1.00000e-01
        """
        new = self.__class__(self.data.copy())
        for e in elist:
            new = new.add_energy_point(mat, mt, k, e)
        return new

    def get_integrals(self):
        """
        Calculate the integral of each energy distribution.

        Returns
        -------
        `pandas.DataFrame`
            dataframe of integrals

        Examples
        --------
        >>> Edistr(minimal_edistrtest).get_integrals()
                                INTEGRAL
        MAT  MT K EIN
        9437 18 0 1.00000e+00 3.00000e+06
                  2.00000e+00 4.50000e+06
        """
        data = self.data
        integrals = data.apply(lambda x: np.trapz(x, x.index))
        df = integrals.to_frame(name="INTEGRAL")
        return df
    def normalize(self):
        """
        Renormalize each outgoing energy distribution to 1.

        Returns
        -------
        `sandy.Edistr`
            renormalized energy distributions

        Examples
        --------
        >>> new = Edistr(minimal_edistrtest).normalize()
        >>> new
        MAT                9437
        MT                   18
        K                     0
        EIN         1.00000e+00 2.00000e+00
        EOUT
        1.00000e-05 1.33333e-07 0.00000e+00
        1.00000e-04 0.00000e+00 4.44444e-08
        1.00000e+00 0.00000e+00 1.55556e-07
        1.00000e+07 0.00000e+00 2.22222e-08
        2.00000e+07 2.00000e-07 0.00000e+00

        >>> new.get_integrals()
                                INTEGRAL
        MAT  MT K EIN
        9437 18 0 1.00000e+00 1.00000e+00
                  2.00000e+00 1.00000e+00
        """
        integrals = self.get_integrals()
        data = self.data
        b = np.array([[val]*data.shape[0]
                     for val in integrals.values.squeeze()]).T
        df = pd.DataFrame(data.values/b, 
                          index=data.index, 
                          columns=data.columns)
        return self.__class__(df)

    def custom_perturbation(self, pert, mat, mt, k, ein_low, ein_high):
        """
        Given a peruration object (fractions), a MAT number, a MT number,
        a subsection number, a lower and an upper incoming energy bound,
        apply the perturbation to the outgoing energy distributions for all
        incident energies comprised within the given boundaries.

        Parameters
        ----------
        pert : `sandy.Pert`
            perturbation object.
        mat : `int`
            MAT number.
        mt : `int`
            MT number.
        k : `int`
            subsection.
        ein_low : TYPE
            lower energy boundary in eV.
        ein_high : TYPE
            upper energy boundary in eV.

        Returns
        -------
        `sandy.Edistr`
            perturbed distributions.

        Notes
        -----
        .. note:: the output distributions are not renormalized.

        .. note:: The energy grid of the perturbation object refers to the
                  outgoing energy distribution.

        .. note:: If the perturbation exceeds 100%, it is truncated.

        Examples
        --------
        >>> orig = Edistr(minimal_edistrtest)
        >>> pert = sandy.Pert([1.3], index=[1e-3])
        >>> orig.custom_perturbation(pert, 9437, 18, 0, 1.5, 2.5)
            MAT  MT  K         EIN        EOUT       VALUE
        0  9437  18  0 1.00000e+00 1.00000e-05 4.00000e-01
        1  9437  18  0 1.00000e+00 2.00000e+07 6.00000e-01
        2  9437  18  0 2.00000e+00 1.00000e-04 2.60000e-01
        3  9437  18  0 2.00000e+00 1.00000e+00 7.00000e-01
        4  9437  18  0 2.00000e+00 1.00000e+07 1.00000e-01
        """
        data = self.data.copy()
        condition = (data.MT == mt) &\
                    (data.MAT == mat) &\
                    (data.K == k) &\
                    (data.EIN < ein_high) &\
                    (data.EIN >= ein_low)
        dfs = []
        dfs.append(data[~condition])
        for ein, df in data[condition].groupby("EIN"):
            series = pert.reshape(df.EOUT).data.loc[df.EOUT]
            # truncate extremes and replace them with boundaries
            px = sandy.Pert(series).truncate()
            df.VALUE *= px.data.values
            dfs.append(df)
        out = pd.concat(dfs)
        return self.__class__(out)

    def _perturb(self, s, normalize=True):
        """
        Apply perturbations to energy distributions.

        Parameters
        ----------
        s : `pandas.DataFrame` or :func:`~sandy.Samples`
            input perturbations or samples.
            If `s` is a `pandas.DataFrame`, its index and columns must have the
            same names and structure as in `self.data`.
            
            .. note:: the energy grid of `s` must be multigroup, i.e., 
                    rendered by a (right-closed) `pd.IntervalIndex`.
            
            If `s` is a :func:`~sandy.Samples` instance

        Returns
        -------
        Edistr : :func:`~Edistr` or `dict` of :func:`~Edistr`
            perturbed chi object if `s` is a `pandas.DataFrame`,
            otherwise dictionary of perturbed chi objects with
            sample numbers as key.

        Examples
        --------
        Get plutonium energy distributions
        >>> e6 = sandy.get_endf6_file("jeff_33", "xs", 942390)
        >>> edistr = sandy.Edistr.from_endf6(e6)

        Apply multiplication coefficient equal to 1 to outgoing energy up to 3e7 eV (upper xs energy limit)
        and incidnet energy of 2e6 eV. 
        >>> index = pd.IntervalIndex.from_breaks([9e-6, 3e7], name="EOUT", closed="right")
        >>> columns = pd.MultiIndex.from_product([[9437], [18], [0], [2e6]], names=["MAT", "MT", "K", "EIN"])
        >>> s = pd.DataFrame(1, index=index, columns=columns)
        >>> ep = edistr._perturb(s, normalize=False)   # Because the normalization changes the values
        >>> assert ep.data.equals(edistr.data)
        
        Apply multiplication coefficients equal to 1.20 to outgoing energy up to 3e7 eV (upper xs energy limit)
        and incidnet energy of 2e6 eV. 
        >>> s = pd.DataFrame(1.20, index=index, columns=columns)
        >>> ep = edistr._perturb(s, normalize=False)    # Because being multily for 1.2 the norm. deletes it
        >>> assert not ep.data.equals(edistr.data)
        >>> assert ep.data.loc[:, ep.data.columns != (9437, 18, 0, 2e6)].equals(edistr.data.loc[:, edistr.data.columns != (9437, 18, 0, 2e6)])
        >>> assert ep.data[(9437, 18, 0, 2e6)].equals(edistr.data[(9437, 18, 0, 2e6)] * 1.20)
        
        """
        x = self.data
        
        # reshape indices (energy)
        idx = s.index.get_indexer(x.index)
        # need to copy, or else it returns a view
        # seed https://pandas.pydata.org/pandas-docs/stable/user_guide/indexing.html#returning-a-view-versus-a-copy
        s_ = s.iloc[idx].copy()
        s_.index = x.index
        s_.loc[idx < 0, :] = 1.  # idx = -1 indicates out of range lines

        # reshape columns (MAT and MT)
        idx = s_.columns.get_indexer(x.columns)
        s_ = s_.iloc[:, idx].copy()
        s_.columns = x.columns
        s_.loc[:, idx < 0] = 1.  # idx = -1 indicates out of range lines
        
        edistr = self.__class__(s_ * x)

        if normalize:
            edistr = edistr.normalize()
        return edistr

    @classmethod
    def from_endf6(cls, endf6):
        """
        Extract energy distributions from `Endf6` instance.

        Parameters
        ----------
        endf6 : `sandy.Endf6`
            object containing the ENDF-6 text

        Returns
        -------
        `sandy.Edistr`
            object with tabulated energy distributions

        Show content of `sandy.Edistr` instance.
        >>> e6 = sandy.get_endf6_file("jeff_33", "xs", 942390)
        >>> sandy.Edistr.from_endf6(e6).data.head().iloc[:, 1:5]  
        MAT                9437
        MT                   18
        K                     0
        EIN         5.00000e+02 1.00000e+03 1.00000e+04 5.00000e+04
        EOUT
        1.00000e-05 1.85342e-12 1.85341e-12 1.85329e-12 1.85274e-12
        2.00000e-05 2.62113e-12 2.62112e-12 2.62095e-12 2.62017e-12
        4.00000e-05 3.70684e-12 3.70683e-12 3.70658e-12 3.70548e-12
        6.00000e-05 4.53994e-12 4.53992e-12 4.53962e-12 4.53827e-12
        8.00000e-05 5.24227e-12 5.24225e-12 5.24190e-12 5.24035e-12
        """
        tape = endf6.filter_by(listmf=[5])
        for mat, mf, mt in tape.data:
            sec = tape.read_section(mat, mf, mt)
            for k, pdistr in sec["PDISTR"].items():
                if pdistr["LF"] != 1:
                    msg = "non-tabulated distribution for " +\
                         f"MAT{mat}/MF{mf}/MT{mt}, subsection {k}"
                    logging.warning(msg)
                    continue
                if list(filter(lambda x: x["INT"] != [2], pdistr["EIN"].values())):
                    msg = "found non-linlin interpolation, skip " +\
                         f"distribution for MAT{mat}/MF{mf}/MT{mt}," +\
                         f" subsection {k}"
                    logging.warning(msg)
                    continue  
                data = []
                for ein, v in sorted(pdistr["EIN"].items()):
                    series  = pd.Series(v["EDISTR"],
                                        index=v["EOUT"],
                                        name=(mat, mt, k, ein)).to_frame()
                    data.append(series)
        if not data:
            raise sandy.Error("no tabulated energy distribution was found")
        df = pd.concat(data, axis=1).interpolate(method='slinear', axis=0) \
                      .fillna(0)
        return cls(df)
    def to_endf6(self, endf6):
        """
        Update energy distributions in `Endf6` instance with those available
        in a `Edistr` instance.

        .. warning:: only chi with `(MAT,MT)` combinations that are originally
                    present in the `Endf6` instance are modififed, the others
                    are discarded.
                    The reason behind this is that to reconstruct a endf6
                    section we need info that is not available in the `Edistr`
                    instance itself.

        Parameters
        ----------
        `endf6` : `sandy.Endf6`
            `Endf6` instance

        Returns
        -------
        `sandy.Endf6`
            `Endf6` instance with updated chi
        """
        data = endf6.data.copy()
        mf = 5
        for (mat, mt, k), edistr in self.data.groupby(axis=1, level=[0,1,2]):
            # Must read original section to extract info not given in `Xs`
            # instance, e.g. QI, QM
            frame = edistr[mat, mt, k]
            if (mat, mf, mt) not in endf6.keys:
                continue
            sec = endf6.read_section(mat, mf, mt)
            sec["PDISTR"][k]["NBT_EIN"] = [len(frame.columns.unique())]
            sec["PDISTR"][k]["E_P"] = frame.columns.unique().values[[0,-1]]
            for ein, fk in frame.iteritems():
                sec["PDISTR"][k]["EIN"][ein]["EOUT"] = fk.index.values
                sec["PDISTR"][k]["EIN"][ein]["EDISTR"] = fk.values
                sec["PDISTR"][k]["EIN"][ein]["NBT"] = [len(fk.index)]
            data[mat, mf, mt] = sandy.mf5.write_mf5(sec)
        return sandy.Endf6(data)
