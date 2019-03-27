# -*- coding: utf-8 -*-
"""
Created on Tue Nov 20 13:44:51 2018

@author: Luca Fiorito
"""

import logging
import os
import pdb

import numpy as np
import pandas as pd

from sandy.formats.endf6 import Endf6
from sandy.data import RDD

__author__ = "Luca Fiorito"
__all__ = ["get_bmatrix", "get_qmatrix", "DecayChains", "BMatrix", "QMatrix"]



def get_bmatrix(file="jeff"):
    """Extract Q-matrix dataframe from decay data file.

    Parameters
    ----------
    file : `str`
        ENDF-6 file containing decay data (by default use JEFF-3.3)
    
    Returns
    -------
    `BMatrix`
    """
    if file is "jeff":   
        return BMatrix.from_file(os.path.join(RDD.__path__[0], "RDD.jeff33"))



def get_qmatrix(file="jeff"):
    """Extract Q-matrix dataframe from decay data file.

    Parameters
    ----------
    file : `str`
        ENDF-6 file containing decay data (by default use JEFF-3.3)
    
    Returns
    -------
    `QMatrix`
    """
    if file is "jeff":   
        return QMatrix.from_file(os.path.join(RDD.__path__[0], "RDD.jeff33"))



def get_transition_matrix(file="jeff"):
    """Extract transition matrix from decay data file.
    
    Parameters
    ----------
    file : `str`
        ENDF-6 file containing decay data (by default use JEFF-3.3)

    Returns
    -------
    `TMatrix`
    """
    if file is "jeff":
        return TMatrix.from_file(os.path.join(RDD.__path__[0], "RDD.jeff33"))



class DecayChains(pd.DataFrame):
    """Dataframe of decay chains for several isotopes.
    Each row contain a different decay chain.
    

    **Columns**:
        
        - PARENT : (`int`) `ID = ZZZ * 10000 + AAA * 10 + META` of parent nuclide
        - DAUGHTER : (`int`) `ID = ZZZ * 10000 + AAA * 10 + META` of daughter nuclide
        - YIELD : (`float`) branching ratio (between 0 and 1)
        - CONSTANT : (`float`) decay constant
    
    Methods
    -------
    from_endf6
        Extract dataframe of decay chains from endf6 instance
    from_file
        Extract dataframe of decay chains from file
    get_bmatrix
        extract B-matrix inro dataframe
    get_qmatrix
        extract Q-matrix into dataframe
    get_transition_matrix
        extract transition matrix into dataframe
    """

    labels = ["PARENT", "DAUGHTER", "YIELD", "CONSTANT"]
    
    def get_bmatrix(self):
        """Extract B-matrix dataframe.
        
        Returns
        -------
        `sandy.decay.BMatrix`
            B-matrix associated to the given decay chains
        """
        B = self.pivot_table(index="DAUGHTER", columns="PARENT", values="YIELD", aggfunc=np.sum, fill_value=0.0). \
                 astype(float). \
                 fillna(0)
        np.fill_diagonal(B.values, 0)
        return BMatrix(B)

    def get_qmatrix(self):
        """Extract Q-matrix dataframe.
        
        Returns
        -------
        `sandy.decay.QMatrix`
            Q-matrix associated to the given decay chains
        """
        return self.get_bmatrix().to_qmatrix()

    def get_transition_matrix(self):
        """Extract transition matrix into dataframe.
        
        Returns
        -------
        `sandy.decay.TMatrix`
            transition matrix associated to the given decay chains
        """
        df = self.copy()
        df["YIELD"] *= df["CONSTANT"]
        T = df.pivot_table(index="DAUGHTER", columns="PARENT", values="YIELD", aggfunc=np.sum). \
               astype(float). \
               fillna(0)
        return T

    @classmethod
    def from_file(cls, file):
        """Extract dataframe of decay chains from file.
        
        Parameters
        ----------
        file : `str`
            ENDF-6 file containing decay data
        
        Returns
        -------
        `sandy.decay.DecayChains`
            dataframe of decay chains
        """
        tape = Endf6.from_file(file)
        return cls.from_endf6(tape)
        

    @classmethod
    def from_endf6(cls, endf6):
        """Extract dataframe of decay chains from Endf6 instance.
        
        Parameters
        ----------
        tape : `Endf6`
            Endf6 instance containing decay data
        
        Returns
        -------
        `sandy.decay.DecayChains`
            dataframe of decay chains
        """
        tape = endf6.filter_by(listmf=[8], listmt=[457])
        listrdd = []
        for ix,text in tape.TEXT.iteritems():
            X = endf6.read_section(*ix)
            zam = int(X["ZA"]*10 + X["LISO"])
            listrdd += [dict(zip(cls.labels, (zam, zam, 0, 0)))]
            if "DK" not in X: # Stable isotope
                continue
            for rtyp,dk in X["DK"].items():
                parent = zam
                daughter = zam//10
                neutrons = 0; protons = 0; alphas = 0
                for dtype in map(int, rtyp):
                    if dtype == 1: # Beta decay
                        daughter += 1001 - 1
                    elif dtype == 2: # Electron capture and/or positron emission
                        daughter += 1 - 1001
                    elif dtype == 3: # Isomeric transition
                        pass
                    elif dtype == 4: # Alpha decay
                        daughter -= 2004
                        alphas += 1
                    elif dtype == 5: # Neutron emission
                        daughter -= 1
                        neutrons += 1
                    elif dtype == 6: # Spontaneous fission
                        logging.debug("skip spontaneous fission for {}...".format(parent))
                    elif dtype == 7: # Proton emission
                        daughter -= 1001
                        protons += 1
                    elif dtype == 0: # Gamma emission (not used in MT457)
                        pass
                    else: # Unknown decay mode
                        logging.debug("skip unknown decay mode for {}...".format(parent))
                daughter = int(daughter*10 + dk["RFS"])
                if daughter == parent:
                    continue
                # Add production of daughter to the list of decay chains
                d = dict(zip(cls.labels, (parent, daughter, dk["BR"], X["LAMBDA"])))
                listrdd.append(d)
                # Add decay of parent to the list of decay chains
                d = dict(zip(cls.labels, (parent, parent, -dk["BR"], X["LAMBDA"])))
                listrdd.append(d)
                if neutrons > 0: # add neutrons produced by decay
                    d = dict(zip(cls.labels, (parent, 10, neutrons*dk["BR"], X["LAMBDA"])))
                    listrdd.append(d)
                if protons > 0: # add protons produced by decay
                    d = dict(zip(cls.labels, (parent, 10010, protons*dk["BR"], X["LAMBDA"])))
                    listrdd.append(d)
                if alphas > 0: # add alphas produced by decay
                    d = dict(zip(cls.labels, (parent, 20040, alphas*dk["BR"], X["LAMBDA"])))
                    listrdd.append(d)
        if not listrdd:
            logging.warn("no decay path found in file")
            return pd.DataFrame()
        return cls(listrdd)
    


class BMatrix(pd.DataFrame):
    """Dataframe for a B-matrix, that is, the production yields 
    of all decay chains (sum of branching ratios).

    **Index**:
    
        - DAUGHTER : (`int`) `ID = ZZZ * 10000 + AAA * 10 + META` of daughter nuclide
    

    **Columns**:
        
        - PARENT : (`int`) `ID = ZZZ * 10000 + AAA * 10 + META` of parent nuclide
    
    Methods
    -------
    to_qmatrix
        Convert B-matrix into Q-matrix
    from_file
        Extract B-matrix from file
    """
    
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.index.name = "DAUGHTER"
        self.columns.name = "PARENT"

    def to_qmatrix(self, keep_neutrons=False):
        """Convert B-matrix into Q-matrix.
        
        Parameters
        ----------
        keep_neutrons : `bool`
            if `False` (default) remove neutrons from the matrix
        
        Returns
        -------
        `sandy.decay.BMatrix`
            Q-matrix associated to B-matrix instance
        """
        B = self.copy()
        if not keep_neutrons:
            if 10 in B.index:
                B.drop(index=10, inplace=True)
            if 10 in B.columns:
                B.drop(columns=10, inplace=True)
        C = np.identity(len(B)) - B.values
        Q = np.linalg.pinv(C)
        return QMatrix(Q, index=B.index, columns=B.columns)

    @classmethod
    def from_file(cls, file):
        """Extract B-matrix from file.
        
        Parameters
        ----------
        file : `str`
            ENDF-6 file containing decay data
        
        Returns
        -------
        `sandy.decay.BMatrix`
            B-matrix extracted from given file
        """
        B = DecayChains.from_file(file).get_bmatrix()
        return cls(B)



class QMatrix(pd.DataFrame):
    """Dataframe for a Q-matrix.

    **Index**:
    
        - DAUGHTER : (`int`) `ID = ZZZ * 10000 + AAA * 10 + META` of daughter nuclide
    

    **Columns**:
        
        - PARENT : (`int`) `ID = ZZZ * 10000 + AAA * 10 + META` of parent nuclide
    
    Methods
    -------
    from_file
        Extract Q-matrix from file
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.index.name = "DAUGHTER"
        self.columns.name = "PARENT"
    
    @classmethod
    def from_file(cls, file):
        """Extract Q-matrix from file.
        
        Parameters
        ----------
        file : `str`
            ENDF-6 file containing decay data
        
        Returns
        -------
        `sandy.decay.QMatrix`
            Q-matrix extracted from given file
        """
        Q = DecayChains.from_file(file).get_qmatrix()
        return cls(Q)



class TMatrix(pd.DataFrame):
    """Dataframe for a transition matrix.

    **Index**:
    
        - DAUGHTER : (`int`) `ID = ZZZ * 10000 + AAA * 10 + META` of daughter nuclide
    

    **Columns**:
        
        - PARENT : (`int`) `ID = ZZZ * 10000 + AAA * 10 + META` of parent nuclide
    
    Methods
    -------
    from_file
        Extract transition matrix from file
    """
    
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.index.name = "DAUGHTER"
        self.columns.name = "PARENT"

    @classmethod
    def from_file(cls, file):
        """Extract transition matrix from file.
        
        Parameters
        ----------
        file : `str`
            ENDF-6 file containing decay data
        
        Returns
        -------
        `sandy.decay.TMatrix`
            transition matrix extracted from given file
        """
        T = DecayChains.from_file(file).get_transition_matrix()
        return cls(T)