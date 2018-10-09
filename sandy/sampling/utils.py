# -*- coding: utf-8 -*-
"""
Created on Mon Sep 24 09:29:32 2018

@author: Fiorito_L
"""

import pandas as pd

__author__ = "Luca Fiorito"
__all__ = ["LpcSamples", "EdistrSamples"]

class LpcSamples(pd.DataFrame):
    """samples for Legendre Polynomial coefficients.
    
    index :
        - MAT number
        - MT number
        - order of Legendre polynomial L
        - incominig neutron energy
    
    columns:
        - sample indices
    """    

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.index.names = ["MAT", "MT", "L", "E"]
        ncols = len(self.columns)
        self.columns = range(1, ncols+1)
        self.columns.name = "SMP"



class EdistrSamples(pd.DataFrame):
    """samples for Tabulated energy distributions.
    
    index :
        - MAT number
        - MT number
        - lower bound for incoming energy ELO
        - upper bound for incoming energy EHI
        - outgoing neeutron energy EOUT
    
    columns:
        - sample indices
    """    

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.index.names = ["MAT", "MT", "ELO", "EHI", "EOUT"]
        ncols = len(self.columns)
        self.columns = range(1, ncols+1)
        self.columns.name = "SMP"