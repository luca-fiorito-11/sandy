"""
This module contains all classes and functions dedicated to handling PFNS data and, more in general, 
any tabulated energy distribution provided in MF8 sections.
"""
import logging
from functools import reduce
import pdb

import pandas as pd
import numpy as np

__author__ = "Luca Fiorito"
__all__ = ["Edistr"]


def from_endf6(endf6):
    """Extract tabulated energy distribution from `Endf6` instance.
    
    Parameters
    ----------
    endf6 : `Endf6`
        `Endf6` instance containing tabulated energy distributions
    
    Returns
    -------
    `Edistr`
        global xs/nubar covariance matrix from ENDF6 file
    """
    tape = endf6.filter_by(listmf=[5])
    data = []
    # Loop MF/MT
    logging.debug("found {} edistr sections".format(len(tape)))
    for (mat,mf,mt), text in tape.TEXT.iteritems():
        X = tape.read_section(mat, mf, mt)
        # Loop subsections
        logging.debug("reading section MAT={}/MF={}/MT={}".format(mat, mf, mt))
        logging.debug("found {} partial distributions".format(len(X["PDISTR"])))
        for k,pdistr in X["PDISTR"].items():
            if pdistr["LF"] != 1:
                logging.warn("non-tabulated distribution for MAT{}/MF{}/MT{}, subsec {}".format(mat, mf, mt, k))
                continue
            if list(filter(lambda x:x["INT"] != [2], pdistr["EIN"].values())):
                logging.warn("found non-linlin interpolation, skip energy distr. for MAT{}/MF{}/MT{}, subsec {}".format(mat, mf, mt, k))
                continue
            logging.debug("\treading subsection {} for MAT{}/MF{}/MT{}".format(k, mat, mf, mt))
            logging.debug("\tfound distributions for {} incident energies".format(len(pdistr["EIN"])))
            for ein, v in sorted(pdistr["EIN"].items()):
                for eout, val in zip(v["EOUT"], v["EDISTR"]):
                    data += [(mat, mt, k, ein, eout, val)]
    if not data:
        logging.warn("no tabulated energy distribution was found")
        return pd.DataFrame()
    df = pd.DataFrame.from_records(data, columns=["MAT", "MT", "K", "EIN", "EOUT", "VALUE"])
    df = pd.pivot_table(df, values="VALUE", index=["MAT", "MT", "K", "EIN"], columns="EOUT").T.sort_index().interpolate(method="slinear").fillna(0).T
    return Edistr(df)



class Edistr(pd.DataFrame):
    """Energy distribution object.
    
    Dataframe components
    --------------------
    index :
        - MAT number
        - MT number
        - number of partial energy distribution (K)
        - incoming neutron energy EIN
    
    columns :
        - outgoing neutron energies EOUT

    Attributes
    ----------
    labels : `list`
        names of the dataframe index
        
    Methods
    -------
    add_points
        add outgoing energy distributions at additional incident energies by interpolation
    integrals
        calculate the integral of each energy distribution
    normalize
        renormalize each outgoing energy distribution to 1
    perturb
        apply perturbation coefficients to energy distribution
    to_stack
        move 'EOUT' from columns to index and transform dataframe into series
    """
    
    labels = ["MAT", "MT", "K", "EIN"]

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.index.names = self.labels
        self.sort_index(inplace=True)
        self.columns.name = "EOUT"

    def to_stack(self):
        """Move 'EOUT' from columns to index and transform dataframe into series.
        
        Returns
        -------
        `pandas.Series`
            stack of energy distribution values
        """
        return self.stack().rename("VALUE")

    def add_points(self, extra_points):
        """Add outgoing energy distributions at additional incident energies by interpolation.
        
        Parameters
        ----------
        extra_points : iterable
            energy points in eV
        
        Returns
        -------
        `sandy.core.pfns.Edistr`
            `Edistr` object with additional incoming energies
        """
        frame = self.copy()
        dflist = []
        for (mat,mt,k),df in frame.groupby(["MAT","MT","K"]):
            grid = sorted((set(df.loc[mat, mt, k].index) | set(extra_points)))
            out = df.loc[mat, mt, k].reindex(grid).interpolate(method='slinear').fillna(0)
            out.index = pd.MultiIndex.from_product([[mat],[mt],[k],grid], names=self.labels)
            dflist.append(out)
        edistr = Edistr(pd.concat(dflist, axis=0))
        return edistr

    def integrals(self):
        """Calculate the integral of each energy distribution.
        
        Returns
        -------
        `pandas.Series`
            series of integrals indexed by MAT, MT, K, EIN
        """
        stack = self.to_stack()
        data = []
        for (mat,mt,k,ein), edistr in stack.groupby(self.labels):
            dx = edistr.index.get_level_values("EOUT")[1:] - edistr.index.get_level_values("EOUT")[:-1]
            y = (edistr.values[1:] + edistr.values[:-1]) / 2
            data.append({"MAT" : mat, "MT" : mt, "K" : k, "EIN" : ein, "INTEGRALS" : y.dot(dx)})
        return pd.DataFrame.from_dict(data).set_index(self.labels).INTEGRALS

    def normalize(self):
        """Renormalize each outgoing energy distribution to 1.
        
        Returns
        -------
        `sandy.core.pfns.Edistr`
            renormalized `Edistr` object
        """
        List = []
        for i,v in self.iterrows():
            dx = v.index.values[1:] - v.index.values[:-1]
            y = (v.values[1:]+v.values[:-1])/2
            List.append(v/y.dot(dx))
        frame = pd.DataFrame(List)
        frame.index = pd.MultiIndex.from_tuples(frame.index)
        return Edistr(frame)

    def perturb(self, pert, method=2, normalize=True, **kwargs):
        """Perturb energy distributions given a set of perturbations.
        
        Parameters
        ----------
        pert : pandas.Series
            multigroup perturbations from sandy.EdistrSamples
        method : int
            * 1 : samples outside the range [0, 2*_mean_] are set to _mean_. 
            * 2 : samples outside the range [0, 2*_mean_] are set to 0 or 2*_mean_ respectively if they fall below or above the defined range.
        normalize : bool
            apply normalization
        
        Returns
        -------
        sandy.Edistr
        """
        frame = self.copy()
        for (mat,mt,k),S in self.groupby(["MAT", "MT", "K"]):
            if (mat,mt) not in pert.index:
                continue
            for ein,edistr in S.loc[mat,mt,k].iterrows():
                for (elo,ehi),P in pert.loc[mat,mt].groupby(["ELO","EHI"]):
                    if ein >= elo and ein <= ehi:
                        P = P[elo,ehi]
                        eg = sorted(set(edistr.index) | set(P.index))
                        P = P.reindex(eg).ffill().fillna(0).reindex(edistr.index)
                        if method == 2:
                            P = P.where(P>=-edistr, -edistr)
                            P = P.where(P<=edistr, edistr)
                        elif method == 1:
                            P = np.where(P.abs() <= edistr, P, 0)
                        frame.loc[mat,mt,k,ein] = edistr + P
                        break
        edistr = Edistr(frame)
        if normalize:
            edistr = edistr.normalize()
        return edistr