# -*- coding: utf-8 -*-
"""
Created on Mon Sep 24 09:29:32 2018

@author: Fiorito_L
"""

import pandas as pd

__author__ = "Luca Fiorito"
__all__ = ["LpcSamples", "EdistrSamples", "FySamples"]

class LpcSamples(pd.DataFrame):
    """samples for Legendre Polynomial coefficients.
    
    Index
    -----
    MAT : `int`
        MAT number
    MT : `int`
        MT number
    L : `int`
        order of Legendre polynomial
    E : `float`
        incoming energy
    
    Columns
    -------
    sample indices
    """    

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.index.names = ["MAT", "MT", "L", "E"]
        ncols = len(self.columns)
        self.columns = range(1, ncols+1)
        self.columns.name = "SMP"



class EdistrSamples(pd.DataFrame):
    """samples for Tabulated energy distributions.
    
    Index
    -----
    MAT : `int`
        MAT number
    MT : `int`
        MT number
    ELO : `float`
        lower bound for incoming energy
    EHI : `float`
        upper bound for incoming energy
    EOUT : `float`
        outgoing neutron energy
    
    Columns
    -------
    sample indices
    """    

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.index.names = ["MAT", "MT", "ELO", "EHI", "EOUT"]
        ncols = len(self.columns)
        self.columns = range(1, ncols+1)
        self.columns.name = "SMP"



class FySamples(pd.DataFrame):
    """Samples for fission yields.
    
    Index
    -----
    MAT : `int`
        MAT number
    MT : `int`
        MT number
    E : `float`
        incoming neutron energy
    ZAM : `int`
        ZZZ * 10000 + A * 10 + M
    
    Columns
    -------
    sample indices
    """    

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.index.names = ["MAT", "MT", "E", "ZAM"]
        ncols = len(self.columns)
        self.columns = range(1, ncols+1)
        self.columns.name = "SMP"        