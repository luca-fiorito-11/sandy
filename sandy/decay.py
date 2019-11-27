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
This module contains all classes and functions dedicated to the processing and 
analysis of a decay data.



.. _Routines:

Routines
========

endf2hdf
"""

import logging
import pdb
import h5py
import os

import numpy as np
import pandas as pd

import sandy

__author__ = "Luca Fiorito"

__all__ = [
        "DecayData",
        "decay_modes",
        ]


decay_modes = {
        0 : "gamma",
        1 : "beta",
        2 : "e.c.",
        3 : "i.t.",
        4 : "alpha",
        5 : "n",
        6 : "s.f.",
        7 : "p",
        }



class DecayData():
    """
    Container of radioactive nuclide data for several isotopes.
    
    Attributes
    ----------
    data : `dict`
        source of decay data content

    Methods
    -------
    get_bmatrix
        extract B-matrix inro dataframe
    get_decay_chains
        extract decay chains into dataframe
    get_qmatrix
        extract Q-matrix into dataframe
    get_transition_matrix
        extract transition matrix into dataframe
    from_endf6
        extract decay data from `sandy.Endf6` instance
    from_hdf5
        extract decay data from hdf5 file
    to_hdf5
        write decay data to hdf5 file
    """
   
    def __repr__(self):
        return self.data.__repr__()
    
    def __init__(self, dct):
        self.data = dct
    
    @property
    def data(self):
        """
        Dictionary of RDD content.
        
        Returns
        -------
        `dict`
            hierarchical RDD content
        
        Raises
        ------
        `sandy.Error`
            if `data` is not a dictionary 
        """
        return self._data
    
    @data.setter
    def data(self, data):
        if not isinstance(data, dict):
            raise sandy.Error("'data' is not a 'dict'")
        self._data = data
   
    def to_hdf5(self, filename, lib):
        """
        Write decay data to hdf5 file.
        
        Parameters
        ----------
        filename : `str`
            name of the hdf5 file (with absolute or relative path)
        lib : `str`
            name of the library that will be used
        """
        path_in_h5 = '/decay/{}/'.format(lib)
        with h5py.File(filename, 'w') as h5file:
            sandy.tools.recursively_save_dict_contents_to_group(h5file, path_in_h5 , self.data)

    def get_decay_chains(self, **kwargs):
        """
        Extract decay chains into dataframe.
        
        Returns
        -------
        `pandas.DataFrame`
            decay chains dataframe
        
        Raises
        ------
        `AssertionError`
            when key `"decay_mode"` is not present but `"decay_constant"` is larger than 0
        `AssertionError`
            when key `"decay_products"` is not present and decay mode is not spontaneous fission
        """
        columns = ["PARENT", "DAUGHTER", "YIELD", "BR", "LAMBDA"]
        items = []
        for zam,v in sorted(self.data.items()):
            if v["stable"]:
                assert v["decay_constant"] == 0
                continue
            # add also the disappearance of the parent
            items.append((zam, zam, -1., 1., v["decay_constant"]))
            for km, vm in v["decay_modes"].items():
                if "decay_products" not in vm:
                    continue
                for zap,yld in vm["decay_products"].items():
                    # add the production of each daughter
                    items.append((zam, zap, yld, vm["branching_ratio"], v["decay_constant"]))
        return pd.DataFrame(items, columns=columns).sort_values(by=["PARENT", "DAUGHTER"])

    def get_bmatrix(self, **kwargs):
        """
        Extract B-matrix into dataframe.
        
        Returns
        -------
        `pandas.DataFrame`
            B-matrix associated to the given decay chains
        """
        B = self.get_decay_chains(**kwargs) \
                .pivot_table(index="DAUGHTER", columns="PARENT", values="YIELD", aggfunc=np.sum, fill_value=0.0) \
                .astype(float) \
                .fillna(0)
        np.fill_diagonal(B.values, 0)
        return B.reindex(B.columns.values, fill_value=0.0)

    def get_qmatrix(self, keep_neutrons=False, **kwargs):
        """
        Extract Q-matrix dataframe.
        
        Returns
        -------
        `pandas.DataFrame`
            Q-matrix associated to the given decay chains
        """
        B = self.get_bmatrix(**kwargs)
        if not keep_neutrons:
            if 10 in B.index:
                B.drop(index=10, inplace=True)
            if 10 in B.columns:
                B.drop(columns=10, inplace=True)
        C = np.identity(len(B)) - B.values
        Q = np.linalg.pinv(C)
        return pd.DataFrame(Q, index=B.index, columns=B.columns)

    def get_transition_matrix(self):
        """
        Extract transition matrix into dataframe.
        
        Returns
        -------
        `pandas.DataFrame`
            transition matrix associated to the given decay chains
        """
        df = self.get_decay_chains()
        df["YIELD"] *= df["LAMBDA"]*df["BR"]
        T = df.pivot_table(index="DAUGHTER", columns="PARENT", values="YIELD", aggfunc=np.sum). \
               astype(float). \
               fillna(0)
        return T.reindex(T.columns.values, fill_value=0.0)

    @classmethod
    def from_endf6(cls, endf6):
        """
        Extract hierarchical structure of decay data from `sandy.Endf6` instance.
        
        Parameters
        ----------
        tape : `sandy.Endf6`
            instance containing decay data
        
        Returns
        -------
        `dict`
            structured container with RDD.
        
        Raises
        ------
        `sandy.Error`
            if no decay data is found
        """
        tape = endf6.filter_by(listmf=[8], listmt=[457])
        if tape.empty:
            raise sandy.Error("no decay data found in file")
        groups = {}
        for (mat,mf,mt),text in tape.TEXT.iteritems():
            sec = endf6.read_section(mat,mf,mt)
            zam = int(sec["ZA"]*10 + sec["LISO"])
            groups[zam] = {
                    "half_life" : sec["HL"],
                    "decay_constant" : sec["LAMBDA"], 
                    "decay_modes" : {},
                    "stable" : bool(sec["NST"]),
                    "spin" : sec["SPI"],
                    "parity" : sec["PAR"],
                    "decay_energy" : {
                            "beta" : sec["E"][0],
                            "gamma" : sec["E"][1],
                            "alpha" : sec["E"][2],
                            "total" : sum(sec["E"][:3]),
                            },
                    "decay_energy_uncertainties" : {
                            "beta" : sec["DE"][0],
                            "gamma" : sec["DE"][1],
                            "alpha" : sec["DE"][2],
                            },
                    }
            if "DK" not in sec: # Stable isotope
                groups[zam]["stable"] = True
                continue
            for rtyp, dk in sec["DK"].items():
                residual_state = dk["RFS"]
                decay_mode_data = {
                        "decay_products" : get_decay_products(rtyp, zam, residual_state),
                        "branching_ratio" : dk["BR"],
                        }
                groups[zam]["decay_modes"][rtyp] = decay_mode_data
        return cls(groups)

    @classmethod
    def from_hdf5(cls, filename, lib):
        """
        Extract hierarchical structure of decay data from hdf5 file.
        
        Parameters
        ----------
        filename : `str`
            hdf5 filename (absolute or relative)
        lib : `str`
            library ID contained in the hdf5 file
        
        Returns
        -------
        `DecayData`
            decay data object
        """
        with h5py.File(filename, 'r') as h5file:
            data = sandy.tools.recursively_load_dict_contents_from_group(h5file, '/decay/{}/'.format(lib))
        return cls(data)



def expand_decay_type(zam, dectyp):
    """
    Given an individual decay mode, return the decay products of a 
    given isotope.
    
    Parameters
    ----------
    zam : `int`
        ZAM identifier
    dectyp : `int`
        decay mode
    
    Returns
    -------
    `int`
        decay daughter product ZAM identifier
    `float`
        number of emitted neutrons
    `float`
        number of emitted protons
    `float`
        number of emitted alphas
    """
    daughter = zam//10
    neutrons = 0.
    protons = 0.
    alphas = 0.
    if dectyp == 1: # Beta decay
        daughter += 1001 - 1
    elif dectyp == 2: # Electron capture and/or positron emission
        daughter += 1 - 1001
    elif dectyp == 3: # Isomeric transition
        pass
    elif dectyp == 4: # Alpha decay
        daughter -= 2004
        alphas += 1.
    elif dectyp == 5: # Neutron emission
        daughter -= 1
        neutrons += 1.
    elif dectyp == 6: # Spontaneous fission
        pass
    elif dectyp == 7: # Proton emission
        daughter -= 1001
        protons += 1.
    elif dectyp == 0: # Gamma emission (not used in MT457)
        pass
    else: # Unknown decay mode
        raise sandy.SandyError("unknown decay mode {} for ZAM={}...".format(dectyp, zam))
    return daughter*10, neutrons, protons, alphas



def get_decay_products(rtyp, zam, rfs):
    """
    For a given isotope and decay mode, extract a dictionary of decay products.

    Parameters
    ----------
    rtyp : `str`
        string of decay modes where:
            1. Beta decay
            2. Electron capture and/or positron emission
            3. Isomeric transition
            4. Alpha decay
            5. Neutron emission (not delayed neutron decay)
            6. Spontaneous fission
            7. Proton emission
        
        Decay mode combinations are allowed, e.g. "15" means Beta decay 
        followed by neutron emission (delayed neutron decay).
    zam : `int`
        ZAM identifier of the nuclide undergoing decay (parent)
    rfs : `int`
        Isomeric state flag for daughter nuclide, e.g. `rfs=0` is ground state, 
        `rfs=1` is first isomeric state, etc.
    
    Returns
    -------
    `dict`
        dictionary of decay products where the keys are the ZAM identifiers for the products
        and the values are the corresponding yield.       
    """
    daughter = zam + 0
    neutrons = 0.
    protons = 0.
    alphas = 0.
    for dectyp in map(int, rtyp):
        daughter, n, h, a = expand_decay_type(daughter, dectyp)
        neutrons += n
        protons += h
        alphas += a
    daughter = int(daughter + rfs)
    products = {}
    if daughter != zam:
        products[daughter] = 1.0
    if neutrons != 0:
        products[10] = neutrons 
    if protons != 0:
        products[10010] = protons
    if alphas != 0:
        products[20040] = alphas
    return products



def endf2hdf(e6file, h5file, lib):
    """
    Write to disk a HDF5 file that reproduces the content of a RDD file in 
    ENDF6 format.
    
    Parameters
    ----------
    e6file : `str`
        filename (with absolute or relative path) of the ENDF6 file
    h5file : `str`
        filename (with absolute or relative path) of the HDF5 file
    lib : `str`
        library name (it will appear as a hdf5 group)
    """
    endf6 = sandy.read_formatted_file()
    from_endf6(endf6).to_hdf5(h5file, lib)
