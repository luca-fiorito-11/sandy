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
    [[9437, 18, 0, 1e0, 1e-5, 0.4],
     [9437, 18, 0, 1e0, 2e7, 0.6],
     [9437, 18, 0, 2e0, 1e-4, 0.2],
     [9437, 18, 0, 2e0, 1e0, 0.7],
     [9437, 18, 0, 2e0, 1e7, 0.1]],
    columns=["MAT", "MT", "K", "EIN", "EOUT", "VALUE"]
    )


class Edistr():
    """
    Object to store tabulate energy distributions.

    Attributes
    ----------
    data : `pandas.DataFrame`
        dataframe of energy distribution data with the following columns

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
        etract energy distributions from `Endf6` instance
    get_integrals
        calculate the integral of each energy distribution
    get_table
        pivot dataframe of tabulated energy spectra
    normalize
        renormalize each outgoing energy distribution to 1
    """

    _labels = ["MAT", "MT", "K", "EIN"]

    def __repr__(self):
        return self.data.__repr__()

    def __init__(self, df, **kwargs):
        self.data = pd.DataFrame(df, **kwargs)

    @property
    def data(self):
        """
        Dataframe of energy distribution data with the following columns:

            - `MAT` : MAT number
            - `MT` : MT number
            - `K` : subsection number
            - `EIN` : incoming energy
            - `EOUT` : outgoing energy
            - `VALUE` : tabulated value of the distribution

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
            MAT  MT  K         EIN        EOUT       VALUE
        0  9437  18  0 1.00000e+00 1.00000e-05 4.00000e-01
        1  9437  18  0 1.00000e+00 2.00000e+07 6.00000e-01
        2  9437  18  0 2.00000e+00 1.00000e-04 2.00000e-01
        3  9437  18  0 2.00000e+00 1.00000e+00 7.00000e-01
        4  9437  18  0 2.00000e+00 1.00000e+07 1.00000e-01
        """
        return self._data

    @data.setter
    def data(self, data):
        self._data = data.sort_values(by=["MAT", "MT", "K", "EIN", "EOUT"]) \
                         .reset_index(drop=True)

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

    def reshape(self, enew):
        """
        Reshape the outgoing energy.

        Parameters
        ----------
        enew : 1D iterable or `float`
            new energy grid.

        Returns
        -------
        `sandy.Edistr`
            Energy distribution instance over new grid.

        Examples
        --------
        >>> orig = Edistr(minimal_edistrtest)
        >>> orig.reshape(1.5)
            MAT  MT  K         EIN        EOUT       VALUE
        0  9437  18  0 1.00000e+00 1.00000e-05 4.00000e-01
        1  9437  18  0 1.00000e+00 1.50000e+00 4.00000e-01
        2  9437  18  0 1.00000e+00 2.00000e+07 6.00000e-01
        3  9437  18  0 2.00000e+00 1.00000e-04 2.00000e-01
        4  9437  18  0 2.00000e+00 1.00000e+00 7.00000e-01
        5  9437  18  0 2.00000e+00 1.50000e+00 7.00000e-01
        6  9437  18  0 2.00000e+00 1.00000e+07 1.00000e-01

        >>> orig.reshape([1.5e-5, 1.5, 15e6])
            MAT  MT  K         EIN        EOUT       VALUE
        0  9437  18  0 1.00000e+00 1.00000e-05 4.00000e-01
        1  9437  18  0 1.00000e+00 1.50000e-05 4.00000e-01
        2  9437  18  0 1.00000e+00 1.50000e+00 4.00000e-01
        3  9437  18  0 1.00000e+00 1.50000e+07 5.50000e-01
        4  9437  18  0 1.00000e+00 2.00000e+07 6.00000e-01
        5  9437  18  0 2.00000e+00 1.00000e-04 2.00000e-01
        6  9437  18  0 2.00000e+00 1.00000e+00 7.00000e-01
        7  9437  18  0 2.00000e+00 1.50000e+00 7.00000e-01
        8  9437  18  0 2.00000e+00 1.00000e+07 1.00000e-01
        """
        enew_ = np.array(enew) if isinstance(enew, np.ndarray) else enew

        def foo(df, enew):
            df_ = df.loc[:, ["EOUT", "VALUE"]].set_index("EOUT")
            enew_ = np.union1d(df_.index.values, enew)
            enew_ = enew_[(enew_ >= df_.index.min()) & (enew_ <= df_.index.max())]
            new_edistr = sandy.shared.reshape_differential(
                                                            df_.index.values,
                                                            df_.values,
                                                            enew_,
                                                          )
            new_edistr = pd.DataFrame(new_edistr, index=enew_,
                                      columns=df_.columns)
            new_edistr.index.name = 'EOUT'
            return new_edistr
        edistr_reshape = self.data.groupby(['MAT', 'MT', 'K', 'EIN'])\
                             .apply(foo, enew_).reset_index()
        return self.__class__(edistr_reshape)

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
            MAT  MT  K         EIN    INTEGRAL
        0  9437  18  0 1.00000e+00 1.00000e+07
        1  9437  18  0 2.00000e+00 4.00000e+06
        """
        keys = ["MAT", "MT", "K", "EIN"]
        data = []
        for (mat, mt, k, ein), chi in self.data.groupby(keys):
            dx = np.diff(chi.EOUT.values)
            values = chi.VALUE.values
            y = (values[1:] + values[:-1]) / 2.
            integral = y.dot(dx)
            data.append({
                "MAT": mat,
                "MT": mt,
                "K": k,
                "EIN": ein,
                "INTEGRAL": integral,
                })
        return pd.DataFrame(data)

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
            MAT  MT  K         EIN        EOUT       VALUE
        0  9437  18  0 1.00000e+00 1.00000e-05 4.00000e-08
        1  9437  18  0 1.00000e+00 2.00000e+07 6.00000e-08
        2  9437  18  0 2.00000e+00 1.00000e-04 5.00000e-08
        3  9437  18  0 2.00000e+00 1.00000e+00 1.75000e-07
        4  9437  18  0 2.00000e+00 1.00000e+07 2.50000e-08

        >>> new.get_integrals()
            MAT  MT  K         EIN    INTEGRAL
        0  9437  18  0 1.00000e+00 1.00000e+00
        1  9437  18  0 2.00000e+00 1.00000e+00
        """
        integrals = self.get_integrals()
        data = self.data.copy()
        keys = ["MAT", "MT", "K", "EIN"]
        out = []
        for (mat, mt, k, ein), chi in data.groupby(keys):
            mask = (integrals.MAT == mat) & \
                   (integrals.MT == mt) & \
                   (integrals.K == k) & \
                   (integrals.EIN == ein)
            chi.loc[:, "VALUE"] /= integrals[mask].INTEGRAL.squeeze()
            out.append(chi)
        df = pd.concat(out)
        return self.__class__(df)

    def custom_perturbation(self, pert, mat=None, mt=None, k=None,
                            ein_low=None, ein_high=None):
        """
        Given a peruration object (fractions), a MAT number, a MT number,
        a subsection number, a lower and an upper incoming energy bound,
        apply the perturbation to the outgoing energy distributions for all
        incident energies comprised within the given boundaries.

        Parameters
        ----------
        pert : `sandy.Pert` or `pd.Series`
            perturbation object.
        mat : `int`, optional
            MAT number. The default is None.
        mt : `int`
            MT number. The default is None.
        k : `int`
            subsection.  The default is None.
        ein_low : `float`, optional
            lower energy boundary in eV.  The default is None.
        ein_high : `float`, optional
            upper energy boundary in eV.  The default is None.

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
        >>> pert = pd.Series([1.3], index=[1e-3])
        >>> orig.custom_perturbation(pert, 9437, 18, 0, 1.5, 2.5)
            MAT  MT  K         EIN        EOUT       VALUE
        0  9437  18  0 1.00000e+00 1.00000e-05 4.00000e-01
        1  9437  18  0 1.00000e+00 2.00000e+07 6.00000e-01
        2  9437  18  0 2.00000e+00 1.00000e-04 2.60000e-01
        3  9437  18  0 2.00000e+00 1.00000e+00 7.00000e-01
        4  9437  18  0 2.00000e+00 1.00000e+07 1.00000e-01
        """
        if isinstance(pert, pd.Series):
            if mat is not None and mt is not None and k is not None and ein_low is not None and ein_high is not None:
                columns = pd.MultiIndex.from_product([[mat], [mt], [k], [ein_low], [ein_high]],
                                                     names=['MAT', 'MT', 'K', 'ELO', 'EHI'])
                df = pd.DataFrame(pert.values, index=pert.index,
                                  columns=columns)
                pert_ = sandy.Pert(df)
            else:
                print("The input do not have enought information")
                return
        else:
            pert_ = sandy.Pert(pert) if not isinstance(pert, sandy.Pert) else pert
        return self._custom_perturbation(pert_)

    def _custom_perturbation(self, pert):
        """
        Apply the perturbation to the outgoing energy distributions for all
        incident energies comprised within the given boundarie in the pert
        variable columns.

        Parameters
        ----------
        pert : `sandy.Pert`
            perturbation object

        Returns
        -------
        `sandy.Edistr`
            perturbed distributions.

        Examples
        --------
        >>> orig = Edistr(minimal_edistrtest)
        >>> pert = pd.Series([1.3], index=[1e-3])
        >>> columns = pd.MultiIndex.from_product([[9437], [18], [0], [1.5], [2.5]], names=['MAT', 'MT', 'K', 'ELO', 'EHI'])
        >>> df = pd.DataFrame(pert.values, index=pert.index, columns=columns)
        >>> pert_ = sandy.Pert(df)
        >>> orig._custom_perturbation(pert_)
            MAT  MT  K         EIN        EOUT       VALUE
        0  9437  18  0 1.00000e+00 1.00000e-05 4.00000e-01
        1  9437  18  0 1.00000e+00 2.00000e+07 6.00000e-01
        2  9437  18  0 2.00000e+00 1.00000e-04 2.60000e-01
        3  9437  18  0 2.00000e+00 1.00000e+00 7.00000e-01
        4  9437  18  0 2.00000e+00 1.00000e+07 1.00000e-01
        """
        energy_grid = self.data.loc[:, 'EOUT'].unique()
        enew = np.union1d(energy_grid, pert.right.index)
        enew = enew[(enew <= energy_grid.max()) & (enew >= energy_grid.min())]
        u_pert = pert.reshape(enew).right

        def foo(df, pert):
            ein = df.loc[:, 'EIN'].unique()[0]
            mat = df.loc[:, 'MAT'].unique()[0]
            mt = df.loc[:, 'MT'].unique()[0]
            col = pert.columns
            mask = (mat == col.get_level_values('MAT')) & \
                   (mt == col.get_level_values('MT')) & \
                   (ein >= col.get_level_values('ELO')) & \
                   (ein <= col.get_level_values('EHI'))
            pert_ = pert.loc[:, mask]
            if not pert_.empty:
                pert_ = pert_.iloc[:, [0]]\
                             .reindex(index=df.loc[:, "EOUT"].values)\
                             .values\
                             .flatten()
                df["VALUE"] = df['VALUE'].values * pert_
            return df
        pert_edistr = self.data.groupby(['MAT', 'MT', 'K', 'EIN'])\
                          .apply(foo, u_pert)
        return self.__class__(pert_edistr)

    def _perturb(self, pert, method=2, normalize=True, **kwargs):
        """Perturb energy distributions given a set of perturbations.
        
        Parameters
        ----------
        pert : pandas.Series
            multigroup perturbations from sandy.EdistrSamples
        method : int
            * 1 : samples outside the range [0, 2*_mean_] are set to _mean_. 
            * 2 : samples outside the range [0, 2*_mean_] are set to 0 or 2*_mean_ respectively if they fall below or above the defined range.
        normalize : bool
            apply normalization
        
        Returns
        -------
        sandy.Edistr
        """
        frame = self.copy()
        for (mat,mt,k),S in self.groupby(["MAT", "MT", "K"]):
            if (mat,mt) not in pert.index:
                continue
            for ein,edistr in S.loc[mat,mt,k].iterrows():
                for (elo,ehi),P in pert.loc[mat,mt].groupby(["ELO","EHI"]):
                    if ein >= elo and ein <= ehi:
                        P = P[elo,ehi]
                        eg = sorted(set(edistr.index) | set(P.index))
                        P = P.reindex(eg).ffill().fillna(0).reindex(edistr.index)
                        if method == 2:
                            P = P.where(P>=-edistr, -edistr)
                            P = P.where(P<=edistr, edistr)
                        elif method == 1:
                            P = np.where(P.abs() <= edistr, P, 0)
                        frame.loc[mat,mt,k,ein] = edistr + P
                        break
        edistr = Edistr(frame)
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
        """
        tape = endf6.filter_by(listmf=[5])
        data = []
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
                for ein, v in sorted(pdistr["EIN"].items()):
                    for eout, val in zip(v["EOUT"], v["EDISTR"]):
                        dct = {
                            "MAT": mat,
                            "MT": mt,
                            "K": k,
                            "EIN": ein,
                            "EOUT": eout,
                            "VALUE": val,
                            }
                        data.append(dct)
        if not data:
            raise sandy.Error("no tabulated energy distribution was found")
        return cls(data)
