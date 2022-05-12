# -*- coding: utf-8 -*-
"""
This module contains all classes and functions specific for the cross section
class `Xs` that acts as a container for energy-dependent tabulated cross
section values.
"""
import os
import logging
import functools

import numpy as np
import pandas as pd

import sandy

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
    data : `pandas.DataFrame`
        source of energy dependent tabulated cross sections

    Methods
    -------
    reshape
        Interpolate cross sections over new grid structure
    custom_perturbation
        Apply a custom perturbation to a given cross section
    to_endf6
        Update cross sections in `Endf6` instance
    from_endf6
        Extract cross sections/nubar from `Endf6` instance
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
        index : `pandas.Index`
            energy grid in eV
        columns : `pandas.MultiIndex`
            MAT/MT indices
        values : `numpy.array`
            cross sections in barns

        Returns
        -------
        `pandas.DataFrame`
            tabulated xs

        Raises
        ------
        `sandy.Error`
            if energy grid is not monotonically increasing

        Examples
        --------
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
        self._data.index = self._data.index
        if not data.index.is_monotonic_increasing:
            raise sandy.Error("energy grid is not monotonically increasing")

    def reshape(self, eg):
        """
        Linearly interpolate cross sections over new grid structure.

        Parameters
        ----------
        eg : array-like object
            new energy grid

        Returns
        -------
        `Xs`
            cross section instance over new grid

        Warnings
        --------
        The new cross sections are tabulated over the union between
        the old and the given energy grid
        """
        df = self.data
        enew = df.index.union(eg).astype("float").values
        xsnew = sandy.shared.reshape_differential(
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
            applied
        mt : `int`
            MT reaction number of the xs to which perturbations are to be
            applied
        pert : `sandy.Pert`
            tabulated perturbations

        Returns
        -------
        `Xs`
            cross section instance with given series MAT/MT perturbed
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

    def filter_energies(self, energies):
        mask = self.data.index.isin(energies)
        data = self.data.loc[mask]
        return self.__class__(data)

    def to_endf6(self, endf6):
        """
        Update cross sections in `Endf6` instance with those available in a
        `Xs` instance.

        .. warning:: only xs with `(MAT,MT)` combinations that are originally
                     present in the `Endf6` instance are modififed, the others
                     are discarded.
                     The reason behind this is that to reconstruct a endf6
                     section we need info that is not available in the `Xs`
                     instance itself.

        Parameters
        ----------
        `endf6` : `sandy.Endf6`
            `Endf6` instance

        Returns
        -------
        `sandy.Endf6`
            `Endf6` instance with updated xs
        """
        endf6new = self._xs_to_endf6(endf6)
        endf6new = self._nubar_to_endf6(endf6new)
        return endf6new

    def _nubar_to_endf6(self, endf6):
        data = endf6.data.copy()
        mf = 1
        for (mat, mt), xs in self.data.iteritems():
            # Must read original section to extract info not given in `Xs`
            if (mat, mf, mt) not in endf6.keys:
                continue
            sec = endf6.read_section(mat, mf, mt)
            if sec["LNU"] != 2:
                raise sandy.Error("cannot update nubar if not in tabulated "
                                  "format")
            ethresh = sec["E"][0]
            xs = xs.where(xs.index >= ethresh).dropna()
            sec["E"] = xs.index.values
            sec["NU"] = xs.values
            # Assume only 1 interpolation region and it is linear
            sec["NBT"] = [xs.size]
            sec["INT"] = [2]
            data[mat, mf, mt] = sandy.write_mf1(sec)
        return sandy.Endf6(data)

    def _xs_to_endf6(self, endf6):
        data = endf6.data.copy()
        mf = 3
        for (mat, mt), xs in self.data.iteritems():
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
            data[(mat, mf, mt)] = sandy.write_mf3(sec)
        return sandy.Endf6(data)

    @classmethod
    def from_endf6(cls, endf6):
        """
        Extract cross sections from `Endf6` instance.

        .. note:: xs are linearized on a unique grid.

        .. note:: missing points are linearly interpolated if inside the energy
                  domain, else zero is assigned.

        Parameters
        ----------
        `endf6` : `sandy.Endf6`
            `Endf6` instance

        Returns
        -------
        `sandy.Xs`
            xs tabulated data

        Raises
        ------
        `sandy.Error`
            if interpolation scheme is not lin-lin
        `sandy.Error`
            if requested cross section was not found

        Warns
        -----
        `logging.warning`
            if duplicate energy points are found

        Notes
        -----
        .. note:: Cross sections are linearized on a unique grid.

        .. note:: Missing points are linearly interpolated if inside the energy
                  domain, else zero is assigned.

        .. note:: Duplicate energy points will be removed, only the first one
                  is kept.
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
            raise sandy.Error("cross sections were not found")
        # should we sort index?

        def foo(l, r):
            how = "outer"
            return pd.merge(l, r, left_index=True, right_index=True, how=how)

        df = functools.reduce(foo, data) \
                      .interpolate(method='slinear', axis=0) \
                      .fillna(0)
        return cls(df)

    def _reconstruct_sums(self, drop=True, inplace=False):
        """
        Reconstruct redundant xs.
        """
        df = self.data.copy()
        for mat in self.data.columns.get_level_values("MAT").unique():
            for parent, daughters in sorted(redundant_xs.items(), reverse=True):
                daughters = [x for x in daughters if x in df[mat]]
                if daughters:
                    df[mat,parent] = df[mat][daughters].sum(axis=1)
            # keep only mts present in the original file
            if drop:
                todrop = [x for x in df[mat].columns if x not in self.data[mat].columns]
                cols_to_drop = pd.MultiIndex.from_product([[mat], todrop])
                df.drop(cols_to_drop, axis=1, inplace=True)
        if inplace:
            self.data = df
        else:
            return Xs(df)
#        frame = self.copy()
#        for mat in frame.columns.get_level_values("MAT").unique():
#            for parent, daughters in sorted(Xs.redundant_xs.items(), reverse=True):
#                daughters = [ x for x in daughters if x in frame[mat]]
#                if daughters:
#                    frame[mat,parent] = frame[mat][daughters].sum(axis=1)
#            # keep only mts present in the original file
#            if drop:
#                todrop = [ x for x in frame[mat].columns if x not in self.columns.get_level_values("MT") ]
#                frame.drop(pd.MultiIndex.from_product([[mat], todrop]), axis=1, inplace=True)
#        return Xs(frame)

    def perturb(self, pert):
        """Perturb cross sections/nubar given a set of perturbations.
        
        Parameters
        ----------
        pert : pandas.Series
            multigroup perturbations from sandy.Samples
        method : int
            * 1 : samples outside the range [0, 2*_mean_] are set to _mean_. 
            * 2 : samples outside the range [0, 2*_mean_] are set to 0 or 2*_mean_ respectively if they fall below or above the defined range.
        
        Returns
        -------
        `sandy.formats.utils.Xs`
        """
        frame = self.copy()
        for mat in frame.columns.get_level_values("MAT").unique():
            if mat not in pert.index.get_level_values("MAT"):
                continue
            for mt in frame[mat].columns.get_level_values("MT").unique():
                lmtp = pert.loc[mat].index.get_level_values("MT").unique()
                mtPert = None
                if lmtp.max() == 3 and mt >= 3:
                    mtPert = 3
                elif mt in lmtp:
                    mtPert = mt
                else:
                    for parent, daughters in sorted(self.__class__.redundant_xs.items(), reverse=True):
                        if mt in daughters and not list(filter(lambda x: x in lmtp, daughters)) and parent in lmtp:
                            mtPert = parent
                            break
                if not mtPert:
                    continue
                P = pert.loc[mat,mtPert]
                P = P.reindex(P.index.union(frame[mat,mt].index)).ffill().fillna(1).reindex(frame[mat,mt].index)
                if method == 2:
                    P = P.where(P>0, 0.0)
                    P = P.where(P<2, 2.0)
                elif method == 1:
                    P = P.where((P>0) & (P<2), 1.0)
                xs = frame[mat,mt].multiply(P, axis="index")
                frame[mat,mt] = xs
        return Xs(frame).reconstruct_sums()

    @classmethod
    def _from_errorr(cls, errorr):
        """
        Extract cross sections/nubar from ERRORR instance.

        Parameters
        ----------
        errorr : `sandy.formats.endf6.Errorr`
            ERRORR instance

        Returns
        -------
        `sandy.formats.utils.Xs`
            dataframe of cross sections in ERRORR file
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
        return Xs(frame)

    @classmethod
    def from_file(cls, file, kind="endf6"):
        """
        Read cross sections directly from file.

        Parameters
        ----------
        file : `str`
            file name with relative or absolute path
        kind : `str`, optional, default is `'endf6'`
            type of file

        Returns
        -------
        `sandy.Xs`
            cross sections tabulated data

        Examples
        --------
        >>> file = os.path.join(sandy.data.__path__[0], "h1.pendf")
        >>> sandy.Xs.from_file(file).data.head()
        MAT                 125                        
        MT                  1           2           102
        E                                              
        1.00000e-05 3.71363e+01 2.04363e+01 1.66999e+01
        1.03125e-05 3.68813e+01 2.04363e+01 1.64450e+01
        1.06250e-05 3.66377e+01 2.04363e+01 1.62013e+01
        1.09375e-05 3.64045e+01 2.04363e+01 1.59682e+01
        1.12500e-05 3.61812e+01 2.04363e+01 1.57448e+01
        """
        if kind != "endf6":
            raise ValueError("sandy can only read cross sections from 'endf6' "
                             "files")
        tape = sandy.Endf6.from_file(file)
        return cls.from_endf6(tape)

    def eV2MeV(self):
        """
        Produce dataframe of cross sections with index in MeV instead of eV.

        Returns
        -------
        `pandas.DataFrame`
            dataframe of cross sections with enery index in MeV

        Examples
        --------
        >>> index = [1e-5, 2e7]
        >>> columns = pd.MultiIndex.from_tuples([(9437, 1)])
        >>> sandy.Xs([1, 2], index=index, columns=columns).eV2MeV()
        MAT                9437
        MT                    1
        E                      
        1.00000e-11 1.00000e+00
        2.00000e+01 2.00000e+00
        """
        df = self.data.copy()
        df.index = df.index * 1e-6
        return df
