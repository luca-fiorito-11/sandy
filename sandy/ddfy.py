# -*- coding: utf-8 -*-
"""
Created on Tue Nov 20 13:44:51 2018

@author: Fiorito_L
"""

from sandy.formats.endf6 import Endf6
from sandy.data import RDD

__author__ = "Luca Fiorito"
__all__ = ["get_jeff_qmatrix"]

def get_jeff_qmatrix():
    """Extract Q-matrix dataframe using JEFF-3.3 deacy data file.
    
    Returns
    -------
    `sandy.QMatrix`
    """
    tape = Endf6.from_text("\n".join(RDD.endf6))
    return tape.get_qmatrix(verbose=True)