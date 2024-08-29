"""
This module contains all classes and functions dedicated to handling PFNS
data and, more in general, any tabulated energy distribution provided in
MF5 sections.
"""
import logging


import pandas as pd
import numpy as np

from .pert import Pert
from .core.endf6 import Endf6
from .sections.mf5 import write_mf5

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
            raise NotImplementedError("applied filter returned empty dataframe")
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

        >>> import sandy
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

        >>> import sandy
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
        
        >>> import sandy
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
        
        >>> import sandy
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
        
        >>> import sandy
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
            px = Pert(series).truncate()
            df.VALUE *= px.data.values
            dfs.append(df)
        out = pd.concat(dfs)
        return self.__class__(out)

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
            raise NotImplementedError("no tabulated energy distribution was found")
        return cls(data)

    def to_endf6(self, endf6):
        """
        Update cross sections in :obj:`~sandy.core.endf6.Endf6` instance with those available in a
        :obj:`~sandy.edistr.Edistr` instance.

        Parameters
        ----------
        `endf6` : :obj:`~sandy.core.endf6.Endf6`
            ENDF6 object.

        Returns
        -------
        endf6new : :obj:`~sandy.core.endf6.Endf6`
            ENDF6 object with updated xs.

        Examples
        --------

        >>> import sandy
        >>> import numpy as np
        >>> endf6 = sandy.get_endf6_file('jeff_33', 'xs', 922350)
        >>> ed = sandy.Edistr.from_endf6(endf6)

        >>> np.testing.assert_equal(
        ... sandy.read_mf5(ed.to_endf6(endf6), 9228, 18),
        ... sandy.read_mf5(endf6, 9228, 18)
        ... )

        >>> s = ed.data.query("MT==18")['VALUE'].copy()
        >>> s.loc[0] = 1234567
        >>> ed.data["VALUE"] = s
        >>> newEndf6 = ed.to_endf6(endf6)

        >>> assert newEndf6.data.keys() == endf6.data.keys()

        >>> assert sandy.read_mf5(newEndf6, 9228, 18)["PDISTR"][0]["EIN"][1e-5]["EDISTR"][0] == 1234567

        >>> np.testing.assert_equal(
        ... sandy.read_mf5(newEndf6, 9228, 18)["PDISTR"][0]["EIN"][1e-5]["EDISTR"][1:],
        ... sandy.read_mf5(endf6, 9228, 18)["PDISTR"][0]["EIN"][1e-5]["EDISTR"][1:]
        ... )

        >>> np.testing.assert_equal(
        ... sandy.read_mf5(newEndf6, 9228, 18)["PDISTR"][0]["EIN"][1e-5]["EOUT"],
        ... sandy.read_mf5(endf6, 9228, 18)["PDISTR"][0]["EIN"][1e-5]["EOUT"]
        ... )
        """
        data = endf6.data.copy()
        mf = 5
        for (mat, mt), df in self.data.groupby(["MAT", "MT"]):
            # Must read original section to extract info not given in `Xs`
            # instance, e.g. QI, QM
            if (mat, mf, mt) not in endf6.data.keys():
                continue
            sec = endf6.read_section(mat, mf, mt)
            for (k, ein), dd in df.groupby(['K', 'EIN']):
                # EOUT rewritten in case points were added to the grid
                sec["PDISTR"][k]["EIN"][ein]["EOUT"] = dd["EOUT"].values
                sec["PDISTR"][k]["EIN"][ein]["EDISTR"] = dd["VALUE"].values
            data[mat, mf, mt] = write_mf5(sec)
        return Endf6(data)
