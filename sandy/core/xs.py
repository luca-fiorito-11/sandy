# -*- coding: utf-8 -*-
"""
Outline
=======
1. Summary_
2. Examples_
3. Routines_

.. _Summary:

Summary
=======
This module contains all classes and functions specific for the cross section 
class `Xs` that acts as a container for energy-dependent tabulated cross 
section values.

.. important:: once created, the `Xs` instance should be modified only 
               using its public API, to gurantee that the information 
               stored in its attributes is preserved.

.. _Examples:

Examples
========

Extract cross sections from ENDF-6 file
---------------------------------------

First, import `sandy` and define the ENDF-6 file name with its absolute path or a path 
relative to the current working directory.

>>> import sandy
>>> file = "path/to/endf6/file.txt"

Second, import the text content of the file into a hierarchical object `Endf6`.

>>> tape = sandy.read_formatted_file(file)

Then, extract all cross sections into a dedicated object.

>>> xs = sandy.Xs.from_endf6(tape)
>>> xs

The last command shows you the content of the object `xs`.

How to apply a custom perturbation to a given cross section
-----------------------------------------------------------

To apply an energy dependent perturbation to a given cross section you must 
first create the perturbation object.
This can be done by reading the perturbation coefficients from a `.csv` file, as 

>>> file = "path/to/perturbation/file.csv"
>>> pert = sandy.Pert.from_file(file, sep=",")

Notice that we specified the column separator as `","`, the default separator 
used in `.csv` files.

To apply the perturbation to the U-238 fission cross section `(mat=9228, mt=18)` 
stored in `xs` you can 

>>> mat, mt = 9237, 18
>>> xspert = xs.custom_perturbation(mat, mt, pert)
>>> xspert

.. _Routines:

How to update a `Endf6` object with a modified cross section
------------------------------------------------------------

Modified cross sections contained in a `Xs` instance `xspert` can be rewritten into 
text ENDF-6 ascii format and used to update a `Endf6` instance.
 
>>> # `tape` is our `Endf6` instance
>>> newtape = xspert.to_endf6(tape) # newtape is the updated `Endf6` instance

Notice that only the cross sections originally present in `tape` are updated.
Cross sections present in `xs` and not in `tape` are neglected


How to create a text file from a `Endf6` object
-----------------------------------------------
>>> output = "path/to/output/file.txt"

>>> string = tape.write_string()
>>> with open(output, 'w') as f:
>>>    f.write()
>>>

Routines
========

"""
import pdb
import logging
import functools

import numpy as np
import pandas as pd

import sandy

__author__ = "Luca Fiorito"
__all__ = [
        "Xs",
        ]



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
    
    redundant_xs = {107 : range(800,850),
                    106 : range(750,800),
                    105 : range(700,750),
                    104 : range(650,700),
                    103 : range(600,650),
                    101 : range(102,118),
                    18 : (19,20,21,38),
                    27 : (18,101),
                    4 : range(50,92),
                    3 : (4,5,11,16,17,*range(22,38),41,42,44,45),
                    1 : (2,3),
                    452 : (455,456)}

    _indexname = "E"
    _columnsnames = ["MAT", "MT"]
    
    def __repr__(self):
        return self.data.head().__repr__()
    
    def __init__(self, df, isotope=None):
        self.data = df
        self.isotope = isotope
    
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
            if `data` is not a `pandas.DataFrame`
        `sandy.Error`
            if energy grid is not monotonically increasing
        """
        return self._data
    
    @data.setter
    def data(self, data):
        if not isinstance(data, pd.DataFrame):
            raise sandy.Error("'data' is not a 'pandas.DataFrame'")
        self._data = data.astype(float)
        self._data.index = self._data.index.astype(float)
        if not data.index.is_monotonic_increasing:
            raise sandy.Error("energy grid is not monotonically increasing")
        self._data.index.name = self.__class__._indexname
        self._data.columns.names = self.__class__._columnsnames
    
    def reshape(self, eg, inplace=False):
        """
        Linearly interpolate cross sections over new grid structure.
        
        Parameters
        ----------
        eg : array-like object
            new energy grid
        inplace : `bool`, optional, default is `False`
            flag to activate inplace replacement
        
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
        xsnew = sandy.shared.reshape_differential(df.index.values, df.values, enew)
        df = pd.DataFrame(xsnew, index=enew, columns=df.columns)
        if inplace:
            self.data = df
        else:
            return Xs(df)

    def custom_perturbation(self, mat, mt, pert, inplace=False):
        """
        Apply a custom perturbation to a given cross section identified by 
        a MAT and MT number.
        
        Parameters
        ----------
        mat : `int`
            MAT material number of the xs to which perturbations are to be applied
        mt : `int`
            MT reaction number of the xs to which perturbations are to be applied
        pert : `sandy.Pert`
            tabulated perturbations
        inplace : `bool`, optional, default is `False`
            flag to activate inplace replacement
        
        Returns
        -------
        `Xs`
            cross section instance with given series MAT/MT perturbed
        """
        if (mat, mt) not in self.data:
            logging.warning("could not find MAT{}/MT{}, perturbation will not be applied".format(mat, mt))
            u_xs = self
        else:
            enew = np.union1d(self.data.index.values, pert.right.index.values)
            u_xs = self.reshape(enew)
            u_pert = pert.reshape(enew)
            u_xs.data[(mat,mt)] = u_xs.data[(mat,mt)]*u_pert.right.values
        if inplace:
            self.data = u_xs.data
        else:
            return Xs(u_xs.data)

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
        data = endf6.data.copy()
        mf = 3
        for (mat,mt),xs in self.data.iteritems():
            # Must read original section to extract info not given in `Xs` instance, e.g. QI, QM
            if (mat,mf,mt) not in endf6.keys:
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
            data[mat, mf, mt] = sandy.write_mf3(sec)
        return sandy.Endf6(data)
        
    @classmethod
    def from_endf6(cls, endf6):
        """
        Extract cross sections from `Endf6` instance.
        
        .. note:: xs are linearized on a unique grid.

        .. note:: missing points are linearly interpolated if inside the energy domain, 
                  else zero is assigned.

        .. note:: 
        
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
        
        .. note:: Missing points are linearly interpolated if inside the energy domain, 
                  else zero is assigned.
        
        .. note:: Duplicate energy points will be removed, only the first one is kept.
        """
        tape = endf6.filter_by(listmf=[3])
        data = []
        keep = "first"
        for mat, mf, mt in tape.data:
            sec = tape.read_section(mat, mf, mt)
            xs = pd.Series(sec["XS"], index=sec["E"], name=(mat, mt)).rename_axis("E").to_frame()
            mask_duplicates = xs.index.duplicated(keep=keep)
            for energy in xs.index[mask_duplicates]:
                logging.warning("found duplicate energy for MAT{}/MF{}/MT{} at {:.5e} MeV, keep only {} value".format(mat, mf, mt, energy, keep))
            xs = xs[~mask_duplicates]
            if sec['INT'] != [2]:
                raise sandy.Error('MAT{}/MF{}/MT{} interpolation scheme is not lin-lin'.format(mat, mf, mt))
            data.append(xs)
        if not data:
            raise sandy.Error("requested cross sections were not found")
        # should we sort index?
        df = functools.reduce(lambda left,right : pd.merge(left, right, left_index=True, right_index=True, how='outer'), data).\
                       interpolate(method='slinear', axis=0). \
                       fillna(0)
        return cls(df)

    def _reconstruct_sums(self, drop=True):
        """
        Reconstruct redundant xs.
        """
        frame = self.copy()
        for mat in frame.columns.get_level_values("MAT").unique():
            for parent, daughters in sorted(Xs.redundant_xs.items(), reverse=True):
                daughters = [ x for x in daughters if x in frame[mat]]
                if daughters:
                    frame[mat,parent] = frame[mat][daughters].sum(axis=1)
            # keep only mts present in the original file
            if drop:
                todrop = [ x for x in frame[mat].columns if x not in self.columns.get_level_values("MT") ]
                frame.drop(pd.MultiIndex.from_product([[mat], todrop]), axis=1, inplace=True)
        return Xs(frame)

    def _perturb(self, pert, method=2, **kwargs):
        """Perturb cross sections/nubar given a set of perturbations.
        
        Parameters
        ----------
        pert : pandas.Series
            multigroup perturbations from sandy.XsSamples
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
        """Extract cross sections/nubar from ERRORR instance.
        
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
   


def _from_file(file, sep=None, **kwargs):
    """
    Initialize `Xs` object reading it from file.
    The given shall contain energy values in the first column and pointwise 
    cross sections in the the reaining columns (each column is a different xs).
    The default column separator is `'\\s+'` or a separator 
    accepted by `numpy.genfromtxt`.
    
    Parameters
    ----------
    file : `str`
        file name (absolute or relative path)
    sep : `str`, optional, default `None`
        column separator. By default it takes the `numpy.genfromtxt` 
        default separator `'\\s+'`
    **kwargs : `numpy.genfromtxt` properties, optional
    
    Returns
    -------
    `Xs`
        Container for pointwise cross sections

    Raises
    ------
    `aleph.Error`
        if the file containing the spectrum has less than two columns
    """
    data = np.genfromtxt(file, dtype=float, delimiter=sep, **kwargs)
    if data.ndim < 2:
        raise sandy.Error("at least 2 columns should be given in the file")
    df = pd.DataFrame(data[:,1], index=data[:,0])
    return Xs(df)