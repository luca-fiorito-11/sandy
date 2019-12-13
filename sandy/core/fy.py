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
        "Fy",
        ]



class Fy():
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

    _columns = ["MAT", "MT", "ZAM", "ZAP", "E", "FY"]
    
    def __repr__(self):
        return self.data.head().__repr__()
    
    def __init__(self, df):
        self.data = df
    
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
        self._data = data[self._columns]
    
    def energy_table(self, mt, mat=None, zam=None):
        df = self.data
        if mat:
            condition = (df.MAT == mat)
        elif zam:
            condition = (df.ZAM == zam)
        else:
            raise sandy.Error("either keyword argument 'mat' or 'zam' must be provided")
        efy = pd.pivot_table(df[condition & (df.MT==mt)], index="E", columns="ZAP", values="FY", fill_value=0)
        return EFy(efy)

    def expand_zap(self):
        zap = pd.DataFrame(map(sandy.shared.expand_zam, self.data.ZAP), columns=["Z", "A", "M"]) 
        zap["ZAP"] = self.data.ZAP.values
        return self.data.merge(zap, left_on="ZAP", right_on="ZAP")

    def filter_by(self, key, value):
        """
        """
        condition = self.data[key] == value
        out = self.data.copy()[condition]
        return self.__class__(out)

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
        tape = endf6.filter_by(listmf=[8], listmt=[454, 459])
        keys = ("MAT", "MT", "ZAM", "ZAP", "E", "FY")
        data = []
        for mat, mf, mt in tape.data:
            sec = tape.read_section(mat, mf, mt)
            zam = sec["ZAM"]
            for e in sec["E"]:
                for zap in sec["E"][e]["ZAP"]:
                    fy = sec["E"][e]["ZAP"][zap]["FY"]
                    values = (mat, mt, zam, zap, e, fy)
                    data.append(dict(zip(keys, values)))
        df = pd.DataFrame(data)
        return cls(df)

class EFy():
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

    _indexname =  "E"
    
    def __repr__(self):
        return self.data.head().__repr__()
    
    def __init__(self, df):
        self.data = df
    
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
        self._data = data
    