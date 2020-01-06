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
This module contains all classes and functions specific for processing fission 
yield data.

.. _Examples:

Examples
========
Examples can be found on this `jupyter notebook
<https://github.com/luca-fiorito-11/sandy/blob/develop/notebooks/fy_notebook.ipynb>`_.


Routines
========

Fy

"""
import pdb
import logging

import numpy as np
import pandas as pd

import sandy

__author__ = "Luca Fiorito"
__all__ = [
        "Fy",
        ]



class Fy():
    """
    Object for fission yield data.
    
    Attributes
    ----------
    data : `pandas.DataFrame`
        source of fission yield data
    
    Methods
    -------
    custom_perturbation
        Apply a custom perturbation to a given cross section
    energy_table
        Interpolate cross sections over new grid structure
    _expand_zam
        Add columns `Z`, `A` and `M` of the fissioning isotope
    _expand_zap
        Add columns `Z`, `A` and `M` of the fission product
    filter_by
        Apply condition to fission yield data
    from_endf6
         Extract fission yields from `Endf6` instance
    to_endf6
        Update fission yield data in `Endf6` instance
    to_hdf5
        Write fission yield data to hdf5 file
    """

    _columns = ["MAT", "MT", "ZAM", "ZAP", "E", "FY"]
    
    def __repr__(self):
        return self.data.head().__repr__()
    
    def __init__(self, df):
        self.data = df
    
    @property
    def data(self):
        """
        Dataframe of fission yield data with the following columns:
            
            - `MAT` : MAT number
            - `MT` : MT number
            - `ZAM` : `Z*1000 + A*10 + M` for the parent (fissioning) nuclide
            - `ZAP` : `Z*1000 + A*10 + M` for the daughter nuclide (fission product)
            - `E` : fissioning energy
            - `FY` : fission yield (fraction)
        
        Returns
        -------
        `pandas.DataFrame`
            tabulated fission yields
        
        Raises
        ------
        `sandy.Error`
            if `data` is not a `pandas.DataFrame`
        """
        return self._data
    
    @data.setter
    def data(self, data):
        if not isinstance(data, pd.DataFrame):
            raise sandy.Error("'data' is not a 'pandas.DataFrame'")
        self._data = data[self._columns]
    
    def energy_table(self, mt, mat=None, zam=None):
        """
        Produce dataframe of tabulated fission yields as a function of energy.
        
        Parameters
        ----------
        `mt` : `int`
            MT number
        `mat` : `int`, optional, default is `None`
            MAT number, if not given select fissioning isotope by ZAM
        `zam` : `int`, optional, default is `None`
            ZAM number, if not given select fissioning isotope by MAT

        Returns
        -------
        `pandas.DataFrame`
            tabulated fission yields
        
        Raises
        ------
        `sandy.Error`
            if neither keyword argument 'mat' nor 'zam' are given
        """
        df = self.data
        if mat:
            condition = (df.MAT == mat)
        elif zam:
            condition = (df.ZAM == zam)
        else:
            raise sandy.Error("either keyword argument 'mat' or 'zam' must be provided")
        efy = pd.pivot_table(df[condition & (df.MT==mt)], index="E", columns="ZAP", values="FY", fill_value=0)
        return efy

    def _expand_zap(self):
        """
        Produce dataframe with three extra columns containing the `Z`, `A` and 
        `M` numbers of the **parent** (fissioning) nuclide.
        
        Returns
        -------
        `pandas.DataFrame`
            dataframe with Z, A and M columns.
        """
        expand_zam = sandy.shared.expand_zam
        zap = pd.DataFrame(map(expand_zam, self.data.ZAP), columns=["Z", "A", "M"], dtype=int) 
        zap["ZAP"] = self.data.ZAP.values
        return self.data.merge(zap, left_on="ZAP", right_on="ZAP")

    def _expand_zam(self):
        """
        Produce dataframe with three extra columns containing the `Z`, `A` and 
        `M` numbers of the **daughter** nuclide (fission product).
        
        Returns
        -------
        `pandas.DataFrame`
            dataframe with Z, A and M columns.
        """
        expand_zam = sandy.shared.expand_zam
        zam = pd.DataFrame(map(expand_zam, self.data.ZAM), columns=["Z", "A", "M"], dtype=int) 
        zam["ZAM"] = self.data.ZAM.values
        return self.data.merge(zam, left_on="ZAM", right_on="ZAM")

    def filter_by(self, key, value):
        """
        Apply condition to source data and return filtered results in a new 
        `sandy.Fy` instance.
        
        Parameters
        ----------
        `key` : `str`
            any label present in the columns of `data`
        `value` : `int` or `float`
            value used as filetring condition
        
        Returns
        -------
        `sandy.Fy`
            filtered dataframe of fission yields
        
        Notes
        -----
        .. note:: The primary function of this method is to make sure that 
                  the filtered dataframe is still returned as a `Fy` object.
        """
        condition = self.data[key] == value
        out = self.data.copy()[condition]
        return self.__class__(out)

    @classmethod
    def from_endf6(cls, endf6):
        """
        Extract fission yields from `Endf6` instance.
        
        Parameters
        ----------
        `endf6` : `sandy.Endf6`
            object containing the ENDF-6 text
        
        Returns
        -------
        `sandy.Fy`
            fission yield object

        Notes
        -----
        .. note:: Both independent and cumulative fission product yields are 
                  loaded, if found.
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
    
    def _to_hdf5(self, file, lib):
        """
        Write fission yield data to hdf5 file.
        
        Parameters
        ----------
        `file` : `str`
            HDF5 filename (relative or absolute)
        `lib` : `str`
            library name
        
        Raises
        ------
        `sandy.Error`
            ...
        
        Warns
        -----
        `logging.warning`
            ...
        
        Notes
        -----
        .. note:: ...
        """
        # to be written
        pass

    def _to_endf6(self, endf6):
        """
        Update fission_yields in `sandy.Endf6` instance with those available 
        in a `Fy` instance.
        
        Parameters
        ----------
        `endf6` : `sandy.Endf6`
            `Endf6` instance
        
        Returns
        -------
        `sandy.Endf6`
            `Endf6` instance with updated fission yields
        
        Warnings
        --------
        .. warning:: to be decided whthre only fy originally present in endf6 
                     file are updated, or all.

        Raises
        ------
        `sandy.Error`
            ...
        
        Warns
        -----
        `logging.warning`
            ...
        
        Notes
        -----
        .. note:: ...
        """
        # to be written
        pass
    
    def _custom_perturbation(self, pert, inplace):
        """
        Apply a custom perturbation to fission yields.
        The perturbations are applied to the correct yields based on the 
        common indices MT, ZAM, ZAP, and E.
        
        Parameters
        ----------
        pert : ...
            ...
        inplace : `bool`, optional, default is `False`
            flag to activate inplace replacement
        
        Returns
        -------
        `Fy`
            fy instance with given yields perturbed
        """
        # to be written
        pass    