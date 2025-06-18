# -*- coding: utf-8 -*-
"""
This module contains all classes and functions specific for the cross section
class `Xs` that acts as a container for energy-dependent tabulated cross
section values.
"""
import logging
import functools

import numpy as np
import pandas as pd

from .endf6 import Endf6
from .shared import reshape_differential
from .sections.mf1 import write_mf1
from .sections.mf3 import write_mf3

__author__ = "Luca Fiorito"
__all__ = [
        "Xs",
        "redundant_xs",
        ]

pd.options.display.float_format = '{:.5e}'.format

redundant_xs = {
        107: range(800, 850),
        106: range(750, 800),
        105: range(700, 750),
        104: range(650, 700),
        103: range(600, 650),
        101: range(102, 118),
        18: (19, 20, 21, 38),
        27: (18, 101),
        4: range(50, 92),
        3: (4, 5, 11, 16, 17, *range(22, 38), 41, 42, 44, 45),
        1: (2, 3),
        452: (455, 456)
        }


class Xs():
    """
    Object for energy dependent cross sections.

    Attributes
    ----------
    data : `pd.DataFrame`
        source of energy dependent tabulated cross sections

    Methods
    -------
    custom_perturbation
        Apply a custom perturbation to a given cross section.
    from_endf6
        Extract cross sections/nubar from :obj:`~sandy.endf6.Endf6` instance.
    perturb
        
    reshape
        Interpolate cross sections over new grid structure.
    to_endf6
        Update cross sections in :obj:`~sandy.endf6.Endf6` instance.
    """

    redundant_xs = {
            107: range(800, 850),
            106: range(750, 800),
            105: range(700, 750),
            104: range(650, 700),
            103: range(600, 650),
            101: range(102, 118),
            18: (19, 20, 21, 38),
            27: (18, 101),
            4: range(50, 92),
            3: (4, 5, 11, 16, 17, *range(22, 38), 41, 42, 44, 45),
            1: (2, 3),
            452: (455, 456)
            }

    _indexname = "E"
    _columnsnames = ("MAT", "MT")

    def __repr__(self):
        return self.data.__repr__()

    def __init__(self, *args, **kwargs):
        self.data = pd.DataFrame(*args, dtype=float, **kwargs)

    @property
    def data(self):
        """
        Dataframe of energy-dependent tabulated cross sections.

        Attributes
        ----------
        index : `pd.Index`
            Energy grid in eV.
        columns : `pd.MultiIndex`
            MAT/MT indices.
        values : `np.array`
            Cross sections in barns.

        Returns
        -------
        `pd.DataFrame`
            Tabulated xs.

        Raises
        ------
        `NotImplementedError`
            If energy grid is not monotonically increasing.

        Examples
        --------

        Create custom cross section.

        >>> import sandy
        >>> index = [1e-5, 2e7]
        >>> columns = pd.MultiIndex.from_tuples([(9437, 1)])
        >>> sandy.Xs([1, 2], index=index, columns=columns)
        MAT                9437
        MT                    1
        E                      
        1.00000e-05 1.00000e+00
        2.00000e+07 2.00000e+00
        """
        return self._data

    @data.setter
    def data(self, data):
        self._data = data.rename_axis(self.__class__._indexname, axis=0)\
                         .rename_axis(self.__class__._columnsnames, axis=1)
        if not data.index.is_monotonic_increasing:
            raise NotImplementedError("energy grid is not monotonically increasing")

    def reshape(self, eg):
        """
        Linearly interpolate cross sections over new grid structure.

        Parameters
        ----------
        eg : array-like object
            New energy grid.

        Returns
        -------
        :obj:`~sandy.xs.Xs`
            Cross section instance over new grid.

        Warnings
        --------
        The new cross sections are tabulated over the union between
        the old and the given energy grid.
        """
        df = self.data
        enew = df.index.union(eg).astype("float").values
        xsnew = reshape_differential(
            df.index.values,
            df.values,
            enew,
            )
        df = pd.DataFrame(xsnew, index=enew, columns=df.columns)
        return self.__class__(df)

    def custom_perturbation(self, mat, mt, pert):
        """
        Apply a custom perturbation to a given cross section identified by
        a MAT and MT number.

        Parameters
        ----------
        mat : `int`
            MAT material number of the xs to which perturbations are to be
            applied.
        mt : `int`
            MT reaction number of the xs to which perturbations are to be
            applied.
        pert : :obj:`~sandy.pert.Pert`
            Tabulated perturbations.

        Returns
        -------
        :obj:`~sandy.xs.Xs`
            Cross section instance with given series MAT/MT perturbed.
        """
        if (mat, mt) not in self.data:
            msg = f"could not find MAT{mat}/MT{mt}, " +\
                "perturbation will not be applied"
            logging.warning(msg)
            u_xs = self
        else:
            enew = np.union1d(self.data.index.values, pert.right.index.values)
            u_xs = self.reshape(enew)
            u_pert = pert.reshape(enew)
            u_xs.data[(mat, mt)] = u_xs.data[(mat, mt)] * u_pert.right.values
        return self.__class__(u_xs.data)

    def to_endf6(self, endf6):
        """
        Update cross sections in :obj:`~sandy.endf6.Endf6` instance with those available in a
        :obj:`~sandy.xs.Xs` instance.

        .. warning:: only xs with `(MAT,MT)` combinations that are originally
                     present in the :obj:`~sandy.endf6.Endf6` instance are modififed, the others
                     are discarded.
                     The reason behind this is that to reconstruct a ENDF6
                     section we need info that is not available in the :obj:`~sandy.xs.Xs`
                     instance itself.

        Parameters
        ----------
        `endf6` : :obj:`~sandy.endf6.Endf6`
            ENDF6 object.

        Returns
        -------
        endf6new : :obj:`~sandy.endf6.Endf6`
            ENDF6 object with updated xs.
        """
        endf6new = self._xs_to_endf6(endf6)
        endf6new = self._nubar_to_endf6(endf6new)
        return endf6new

    def _nubar_to_endf6(self, endf6):
        data = endf6.data.copy()
        mf = 1
        for (mat, mt), xs in self.data.items():
            # Must read original section to extract info not given in `Xs`
            if (mat, mf, mt) not in endf6.keys:
                continue
            sec = endf6.read_section(mat, mf, mt)
            if sec["LNU"] != 2:
                raise NotImplementedError("cannot update nubar if not in tabulated format")
            ethresh = sec["E"][0]
            xs = xs.where(xs.index >= ethresh).dropna()
            sec["E"] = xs.index.values
            sec["NU"] = xs.values
            # Assume only 1 interpolation region and it is linear
            sec["NBT"] = [xs.size]
            sec["INT"] = [2]
            data[mat, mf, mt] = write_mf1(sec)
        return Endf6(data)

    def _xs_to_endf6(self, endf6):
        data = endf6.data.copy()
        mf = 3
        for (mat, mt), xs in self.data.items():
            # Must read original section to extract info not given in `Xs`
            # instance, e.g. QI, QM
            if (mat, mf, mt) not in endf6.keys:
                continue
            sec = endf6.read_section(mat, mf, mt)
            # Cut threshold xs
            ethresh = sec["E"][0]
            xs = xs.where(xs.index >= ethresh).dropna()
            sec["E"] = xs.index.values
            sec["XS"] = xs.values
            # Assume all xs have only 1 interpolation region and it is linear
            sec["NBT"] = [xs.size]
            sec["INT"] = [2]
            data[mat, mf, mt] = write_mf3(sec)
        return Endf6(data)

    @classmethod
    def from_endf6(cls, endf6):
        """
        Extract cross sections from :obj:`~sandy.endf6.Endf6` instance.

        .. note:: Xs are linearized on a unique grid.

        .. note:: Missing points are linearly interpolated if inside the energy
                  domain, else zero is assigned.

        .. note:: Duplicate energy points will be removed, only the first one
                  is kept.

        Parameters
        ----------
        `endf6` : :obj:`~sandy.endf6.Endf6`
            ENDF6 object.

        Returns
        -------
        :obj:`~sandy.xs.Xs`
            Xs tabulated data.

        Raises
        ------
        `NotImplementedError`
            If interpolation scheme is not lin-lin.

        Warns
        -----
        `logging.warning`
            if duplicate energy points are found

        Examples
        --------
        Get H1 file and process it to PENDF.

        >>> import sandy
        >>> tape = sandy.get_endf6_file("jeff_33", "xs", 10010)
        >>> pendf = tape.get_pendf(minimal_processing=True)
        
        Show content of `sandy.Xs` instance.

        >>> sandy.Xs.from_endf6(pendf).data.head()
        MAT                 125                        
        MT                  1           2           102
        E                                              
        1.00000e-05 3.71363e+01 2.04363e+01 1.66999e+01
        1.03125e-05 3.68813e+01 2.04363e+01 1.64450e+01
        1.06250e-05 3.66377e+01 2.04363e+01 1.62013e+01
        1.09375e-05 3.64045e+01 2.04363e+01 1.59682e+01
        1.12500e-05 3.61812e+01 2.04363e+01 1.57448e+01
        """
        data = []
        # read cross sections
        tape = endf6.filter_by(listmf=[3])
        keep = "first"
        for mat, mf, mt in tape.data:
            sec = tape.read_section(mat, mf, mt)
            if sec['INT'] != [2]:
                logging.warning(f"skip MAT{mat}/MF{mf}/MT{mt} "
                                "because interpolation schme is not lin-lin")
                continue
            xs = pd.Series(sec["XS"], index=sec["E"], name=(mat, mt)) \
                   .rename_axis("E") \
                   .to_frame()
            mask_duplicates = xs.index.duplicated(keep=keep)
            for energy in xs.index[mask_duplicates]:
                logging.warning("found duplicate energy for "
                                f"MAT{mat}/MF{mf}/MT{mt} "
                                f"at {energy:.5e} MeV, keep only {keep} value")
            xs = xs[~mask_duplicates]
            data.append(xs)
        # read nubar
        tape = endf6.filter_by(listmf=[1], listmt=[452, 455, 456])
        keep = "first"
        for mat, mf, mt in tape.data:
            sec = tape.read_section(mat, mf, mt)
            if sec["LNU"] != 2:
                logging.warning(f"skip MAT{mat}/MF{mf}/MT{mt} "
                                "because not tabulated")
                continue
            if sec['INT'] != [2]:
                logging.warning(f"skip MAT{mat}/MF{mf}/MT{mt} "
                                "because interpolation schme is not lin-lin")
                continue
            xs = pd.Series(sec["NU"], index=sec["E"], name=(mat, mt)) \
                   .rename_axis("E") \
                   .to_frame()
            mask_duplicates = xs.index.duplicated(keep=keep)
            for energy in xs.index[mask_duplicates]:
                logging.warning("found duplicate energy for "
                                f"MAT{mat}/MF{mf}/MT{mt} "
                                f"at {energy:.5e} MeV, keep only {keep} value")
            xs = xs[~mask_duplicates]
            data.append(xs)
        if not data:
            raise NotImplementedError("cross sections were not found")
        # should we sort index?

        def foo(l, r):
            how = "outer"
            return pd.merge(l, r, left_index=True, right_index=True, how=how)

        df = functools.reduce(foo, data) \
                      .interpolate(method='slinear', axis=0) \
                      .fillna(0)
        return cls(df)

    def reconstruct_sums(self, drop=True):
        """
        Reconstruct redundant xs according to ENDF-6 rules in Appendix B.
        Redundant cross sections are available in `dict`
        :obj:`~sandy.xs.redundant_xs`.

        Parameters
        ----------
        drop : `bool`, optional
            keep in output only the MT number originally present.
            The default is True.

        Returns
        -------
        :obj:`~sandy.xs.Xs`
            Cross section instance where reconstruction rules are enforced.

        Examples
        --------
        Get ENDF-6 file for H1, process it in PENDF and extract xs.

        >>> import sandy
        >>> tape = sandy.get_endf6_file("jeff_33", "xs", 10010)
        >>> pendf = tape.get_pendf(minimal_processing=True, err=1)
        >>> xs = sandy.Xs.from_endf6(pendf)

        We introduce a perturbation to the elastic scattering xs.

        >>> xs.data[(125, 2)] *= 2
        >>> assert not xs.data[(125, 1)].equals(xs.data[(125, 2)] + xs.data[(125, 102)])

        Reconstructing xs enforces consistency.

        >>> xs1 = xs.reconstruct_sums(drop=True).data
        >>> assert xs1.columns.equals(xs.data.columns)
        >>> assert xs1[(125, 1)].equals(xs1[(125, 2)] + xs1[(125, 102)])

        We can keep all redundant xs with keyword `drop=True`.

        >>> xs2 = xs.reconstruct_sums(drop=False).data
        >>> assert not xs2.columns.equals(xs.data.columns)
        >>> assert xs2[xs1.columns].equals(xs1)

        >>> assert xs2[(125, 101)].equals(xs2[(125, 102)])
        >>> assert xs2[(125, 27)].equals(xs2[(125, 101)])
        >>> assert xs2[(125, 3)].equals(xs2[(125, 27)])

        This example shows that also the inelastic cross section is correclty reconstructed.

        >>> pendf = sandy.get_endf6_file("jeff_33", "xs", 952410).get_pendf(minimal_processing=True, err=1)
        >>> xs = sandy.Xs.from_endf6(pendf)
        >>> xsr = xs.reconstruct_sums(drop=True)
        >>> assert not xs.data[(9543, 4)].equals(xs.data[9543].loc[:, 50:91].sum(axis=1))
        >>> assert xsr.data[(9543, 4)].equals(xsr.data[9543].loc[:, 50:91].sum(axis=1))
        """
        df = self.data.copy()
        for mat, group_ in df.T.groupby("MAT"):
            group = group_.T
            # starting from the lat redundant cross section, find daughters and sum them
            for parent, daughters in sorted(redundant_xs.items(), reverse=True):
                # it must be df, not group, because df is updated
                mask = df[mat].columns.intersection(daughters)
                if not mask.empty:
                    df[(mat, parent)] = df[mat][mask].sum(axis=1)

            # keep only mts present in the original file
            if drop:
                keep = group[mat].columns
                todrop = df[mat].columns.difference(keep)
                df.drop(
                    pd.MultiIndex.from_product([[mat], todrop]),
                    axis=1,
                    inplace=True,
                    )

        return self.__class__(df)

    def _perturb(self, s):
        """
        Apply perturbations to cross sections.

        Parameters
        ----------
        s : `pd.DataFrame`
            Input perturbations or samples.
            Index and columns of the dataframe must have the same names and
            structure as in `self.data`.
            
            .. note:: The energy grid of `s` must be multigroup, i.e., 
                      rendered by a (right-closed) `pd.IntervalIndex`.

        Returns
        -------
        xs : :obj:`~sandy.xs.Xs` or `dict` of :obj:`~sandy.xs.Xs`
            perturbed cross section object if `s` is a `pd.DataFrame`,
            otherwise dictionary of perturbed cross section objects with
            sample numbers as key.

        Examples
        --------
        Get plutonium cross sections.

        >>> import sandy
        >>> pendf = sandy.get_endf6_file("jeff_33", "xs", 942390).get_pendf(err=1, minimal_processing=True)
        >>> xs = sandy.Xs.from_endf6(pendf)

        Apply multiplication coefficient equal to 1 to elastic and inelastic
        scattering cross sections up to 3e7 eV (upper xs energy limit).

        >>> index = pd.IntervalIndex.from_breaks([1e-5, 3e7], name="E", closed="right")
        >>> columns = pd.MultiIndex.from_product([[9437], [2, 4]], names=["MAT", "MT"])
        >>> s = pd.DataFrame(1, index=index, columns=columns)
        >>> xp = xs._perturb(s)
        >>> assert xp.data.equals(xs.data)
        
        Apply multiplication coefficients equal to 1 and to 2 respectively to
        elastic and inelastic scattering cross sections up to 3e7 eV (upper xs energy limit).

        >>> s = pd.DataFrame([[1, 2]], index=index, columns=columns)
        >>> xp = xs._perturb(s)
        >>> assert not xp.data.equals(xs.data)
        >>> assert xp.data.loc[:, xp.data.columns != (9437, 4)].equals(xp.data.loc[:, xs.data.columns != (9437, 4)])
        >>> assert xp.data[(9437, 4)].equals(xs.data[(9437, 4)] * 2)
        
        Apply multiplication coefficients equal to 1 and to 2 respectively to
        elastic and inelastic scattering cross sections up to 2e7 eV.

        >>> index = pd.IntervalIndex.from_breaks([1e-5, 2e7], name="E", closed="right")
        >>> columns = pd.MultiIndex.from_product([[9437], [2, 4]], names=["MAT", "MT"])
        >>> s = pd.DataFrame([[1, 2]], index=index, columns=columns)
        >>> xp = xs._perturb(s)
        >>> assert not xp.data.equals(xs.data)
        >>> assert xp.data.loc[:, xp.data.columns != (9437, 4)].equals(xp.data.loc[:, xs.data.columns != (9437, 4)])
        >>> assert xp.data.loc[:2e7, (9437, 4)].equals(xs.data.loc[:2e7, (9437, 4)] * 2)
        >>> assert xp.data.loc[2e7:, (9437, 4)].iloc[1:].equals(xs.data.loc[2e7:, (9437, 4)].iloc[1:])
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
        s_ = s_.iloc[:, idx]
        s_.columns = x.columns
        s_.loc[:, idx < 0] = 1.  # idx = -1 indicates out of range lines
        
        xs = self.__class__(s_ * x)
        return xs
    
    @classmethod
    def _from_errorr(cls, errorr):
        """
        Extract cross sections/nubar from ERRORR instance.

        Parameters
        ----------
        errorr : :obj:`~sandy.errorr.Errorr`
            ERRORR object.

        Returns
        -------
        :obj:`~sandy.xs.Xs`
            Dataframe of cross sections in ERRORR file.
        """
        mat = errorr.mat[0]
        eg = errorr.energy_grid
        tape = errorr.filter_by(listmf=[3])
        listxs = []
        for (mat,mf,mt),text in tape.TEXT.iteritems():
            X = tape.read_section(mat, mf, mt)
            xs = pd.Series(
                      X["XS"],
                      index=errorr.energy_grid[:-1],
                      name=(mat,mt)
                      ).rename_axis("E").to_frame()
            listxs.append(xs)
        if not listxs:
            logging.warn("no xs/nubar data was found")
            return pd.DataFrame()
        # Use concat instead of merge because indexes are the same
        frame = pd.concat(listxs, axis=1).reindex(eg, method="ffill")
        return cls(frame)


def xs_perturb_worker(xs, n, s, verbose=False):
    """
    

    Parameters
    ----------
    xs : TYPE
        DESCRIPTION.
    n : TYPE
        DESCRIPTION.
    s : `pd.DataFrame`
        see :func:`~sandy.Xs._perturb`
    verbose : TYPE, optional
        DESCRIPTION. The default is False.

    Returns
    -------
    TYPE
        DESCRIPTION.

    """
    if verbose:
        print(f"Processing xs sample {n}...")
    return xs._perturb(s)
