# -*- coding: utf-8 -*-
"""
This module contains all classes needed to organize and structure different 
nuclear data types into python objects.

Nuclear Data Objects (NDO)
==========================

The following objects are considered:

    - `Xs` : dataframe of energy dependent cross section dataframe
    - `Lpc` : dataframe of energy dependent Legendre Polynomial Coefficients
    - `Edistr` : dataframe of outgoing energy distributions for multiple incoming energy

Nuclear Data Covariance Objects (NDCO)
======================================

The following covariance objects are considered:
    
    - `BaseCov` : base covariance object to be inherithed by specific covariance objects
    - `XsCov` : dataframe of multigroup cross section covariances
    - `LpcCov` : dataframe of multigroup Legendre Polynomial Coefficients covariances
    - `EdistrCov` : dataframe of outgoing energy distributions covariances

**Assumptions**:

 * All NDCO must inherit from `pandas.DataFrame`
 * All NCDO must reproduce square covariance matrices
 * All NCDO must have the following methods/attributes:

  - `labels` : list of index/columns names
  - `get_samples` : method to draw random samples
  - `from_endf` : classmethod to retrieve data from an `endf6` instance
  - `from_errorr` : classmethod to retrieve data from a `errorr` instance

.. important:: This module must not import module `endf6`.
"""

import logging
import pdb
import os
from functools import reduce

import pandas as pd
import numpy as np
import scipy as sp
import seaborn as sns
import matplotlib.pyplot as plt

from ..functions import gls, div0
from sandy.settings import SandyError, colors

__author__ = "Luca Fiorito"
__all__ = ["BaseFile", "Xs", "Lpc", "Edistr", "EnergyCov", "XsCov", "EdistrCov", "LpcCov", 
           "Cov", "Fy", "FyCov", "Tpd",
           "LpcSamples", "EdistrSamples", "FySamples"]



class Section(dict):
    pass



class BaseFile(pd.DataFrame):
    """This class is to be inherited by  all classes that parse and analyze 
    nuclear data evaluated files in ENDF-6 or derived (ERRORR) formats.
    
    Index
    -----
    MAT : `int`
        MAT number to identify the isotope
    MF : `int`
        MF number to identify the data type
    MT : `int`
        MT number to identify the reaction

    Columns
    -------
    TEXT : `string`
        MAT/MF/MT section reported as a single string
    """

    @classmethod
    def from_file(cls, file, listmat=range(1,10000), listmf=range(1,100), listmt=range(1,1000)):
        """Create instance by reading a file.
        
        Parameters
        ----------
        file : `str`
            file name
        listmat : `iterable`
            list of MAT number (default all)
        listmf : `iterable`
            list of MF number (default all)
        listmt : `iterable`
            list of MT number (default all)
        """
        with open(file) as f: text = f.read()
        out = cls.from_text(text, listmat=listmat, listmf=listmf, listmt=listmt)
        out.TAPE = os.path.abspath(os.path.realpath(os.path.expandvars(file)))
        out.FILENAME = os.path.basename(out.TAPE)
        return out

    @classmethod
    def from_text(cls, text, listmat=None, listmf=None, listmt=None, empty_err=True):
        """Create instance from string.
        
        Parameters
        ----------
        text : `str`
            string containing the evaluated data
        listmat : `iterable`
            list of MAT number (default all)
        listmf : `iterable`
            list of MF number (default all)
        listmt : `iterable`
            list of MT number (default all)
        """
        from io import StringIO
        tape = pd.read_fwf(
                StringIO(text),
                widths = [66, 4, 2, 3],
                names = ["TEXT", "MAT", "MF", "MT"],
                usecols = ["MAT", "MF", "MT"]
                )
        tape["TEXT"] = text.splitlines(True)

        tape = tape.loc[(tape.MAT>0) & (tape.MF>0) & (tape.MT>0)]. \
               groupby(["MAT","MF","MT"]). \
               apply(lambda x: "".join(x.TEXT.values)). \
               rename("TEXT"). \
               to_frame()
       
#        splitters = tape.loc[(tape.MAT==0) & (tape.MF==0) & (tape.MT==0)].index
#        dfs = []; ibeg = 0
#        for iend in splitters:
#            df = tape[ibeg:iend]
#            for (mat,mf,mt),group in df.loc[(tape.MAT>0) & (tape.MF>0) & (tape.MT>0)].groupby(["MAT","MF","MT"]):
#                # Select only desired sections
#                if listmt is not None and mt not in listmt:
#                    continue
#                if listmat is not None and mat not in listmat:
#                    continue
#                if listmf is not None and mf not in listmf:
#                    continue
#                dfs.append({"MAT" : mat, "MF" : mf, "MT" : mt, "TEXT" : "".join(group.TEXT.values)})
#            ibeg = iend
#        if not dfs:
#            raise SandyError("tape is empty")
#        tape = pd.DataFrame.from_dict(dfs).set_index(["MAT","MF","MT"])

        frame = cls(tape).filter_by(listmat=listmat, listmf=listmf, listmt=listmt)
        if frame.empty and empty_err:
            raise SandyError("tape is empty")
        return frame

    def __init__(self, *args, **kwargs):
        kwargs.update({"columns" : ["TEXT"]})
        super().__init__(*args, **kwargs)
        self.index.names = ['MAT', 'MF', 'MT']
        self.sort_index(level=["MAT","MF","MT"], inplace=True)
        if self.index.duplicated().any():
            raise SandyError("found duplicate MAT/MF/MT")
    
    def add_sections(self, file, sect, kind='replace'):
        """Add MF/MT section from one file to an existing dataframe.
        If they already exist, replace them or keep them according to parameter 
        `kind`.
        """
        keep = "first" if kind is "keep" else "last"
        queries = []
        for mf,mtlist in sect.items():
            if mtlist == "all":
                queries.append("(MF=={})".format(mf))
            else:
                for mt in mtlist:
                    queries.append("(MF=={} & MT=={})".format(mf,mt))
        query = " | ".join(queries)
        newdf = BaseFile.from_file(file).query(query)
        if newdf.empty:
            logging.warn("'{}' does not contain requested sections".format(file))
            return self
        outdf = pd.concat([self, newdf])
        outdf = outdf.reset_index()
        outdf = outdf.drop_duplicates(["MAT","MF","MT"], keep=keep)
        outdf = outdf.set_index(["MAT","MF","MT"])
        return self.__class__(outdf)

    def delete_sections(self, sect):
        """Add MF/MT section from one file to an existing dataframe.
        """
        queries = []
        for mf,mtlist in sect.items():
            if mtlist == "all":
                queries.append("(MF!={})".format(mf))
            else:
                for mt in mtlist:
                    queries.append("(MF!={} & MT!={})".format(mf,mt))
        query = " & ".join(queries)
        newdf = self.query(query)
        if newdf.empty:
            raise SandyError("all sections were deleted")
        return self.__class__(newdf)

    def filter_by(self, listmat=None, listmf=None, listmt=None):
        """Filter dataframe based on MAT, MF, MT values.
        """
        _listmat = range(1,10000) if listmat is None else listmat
        _listmf = range(1,10000) if listmf is None else listmf
        _listmt = range(1,10000) if listmt is None else listmt
        cond_mat = self.index.get_level_values("MAT").isin(_listmat)
        cond_mf = self.index.get_level_values("MF").isin(_listmf)
        cond_mt = self.index.get_level_values("MT").isin(_listmt)
        df = self.loc[cond_mat & cond_mf & cond_mt]
        return self.__class__(df)
    
    @property
    def mat(self):
        return self.index.get_level_values("MAT").unique()

    @property
    def mf(self):
        return self.index.get_level_values("MF").unique()

    @property
    def mt(self):
        return self.index.get_level_values("MT").unique()



class Xs(pd.DataFrame):

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

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.index.name = "E"
        self.columns.names = ["MAT", "MT"]

    def reconstruct_sums(self, drop=True):
        """
        Reconstruct redundant xs.
        """
        frame = self.copy()
        for mat in frame.columns.get_level_values("MAT").unique():
            for parent, daughters in sorted(Xs.redundant_xs.items(), reverse=True):
                daughters = [ x for x in daughters if x in frame[mat].columns]
                if daughters:
                    frame[mat,parent] = frame[mat][daughters].sum(axis=1)
            # keep only mts present in the original file
            if drop:
                todrop = [ x for x in frame[mat].columns if x not in self.columns.get_level_values("MT") ]
                frame.drop(pd.MultiIndex.from_product([[mat], todrop]), axis=1, inplace=True)
        return Xs(frame)

    def perturb(self, pert, method=2, **kwargs):
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
        for mat, mt in frame:
            if mat not in pert.index.get_level_values("MAT").unique(): continue
            lmtp = pert.loc[mat].index.get_level_values("MT").unique()
            mtPert = None
            if mt in lmtp:
                mtPert = mt
            else:
                for parent, daughters in sorted(self.__class__.redundant_xs.items(), reverse=True):
                    if mt in daughters and not list(filter(lambda x: x in lmtp, daughters)) and parent in lmtp:
                        mtPert = parent
                        break
            if not mtPert: continue
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

    def _macs(self, E0=0.0253, Elo=1E-5, Ehi=1E1):
        """
        Calculate Maxwellian averaged cross sections.
        """
        from math import sqrt, pi
        from ..integrals.macs import maxw_int, maxw_xs_int
        # add points to the index
        index = set(self.index.values)
        index.update([Elo, Ehi])
        index = np.array(sorted(index))
        index = index[(index >= Elo) & (index <= Ehi)]
        xs = self.reindex(index).interpolate(method='slinear', axis=0).fillna(0)
        data = [[E0,
                 xs.index[i],
                 xs.index[i+1],
                 maxw_int(E0, xs.index[i], xs.index[i+1])
                 ] for i in range(len(xs)-1)]
        dframe = pd.DataFrame(data, columns=["E0", "E1", "E2", "INT"])
        cond = dframe.E1/E0 >= 1e-5
        records = []
        for (mat,mt),x in xs.items():
            data = [[E0,
                     x.index[i],
                     x.iloc[i],
                     x.index[i+1],
                     x.iloc[i+1],
                     maxw_xs_int(E0, x.index[i], x.iloc[i], x.index[i+1], x.iloc[i+1])
                     ] for i in range(len(x)-1)]
            nframe = pd.DataFrame(data, columns=["E0", "E1", "S1", "E2", "S2", "INT"])
            N = nframe[cond].INT.sum(); D = dframe[cond].INT.sum()
            I = N / D * (2/sqrt(pi))
            skipped = "{}/{}".format(sum(cond==False), len(dframe))
            records.append([mat, mt, I, D, Elo, Ehi, E0, skipped])
        return pd.DataFrame(records, columns=["MAT", "MT", "MACS", "FLUX", "Elo", "Ehi", "E0","SKIPPED"])

    @classmethod
    def from_errorr(cls, errorr):
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

class Lpc(pd.DataFrame):
    """Legendre polynomial coefficients for angular distribution of outgoing particles.
    
    Dataframe components
    --------------------
    index :
        - MAT number
        - MT number
        - incoming neutron energy
    
    columns :
        - Legendre coefficients starting from P0
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.index.names = ["MAT", "MT", "E"]
        self.columns = range(self.shape[1])
        self.columns.name = "L"
        self.sort_index(inplace=True)
    
    def to_stack(self):
        """Convert Lpc instance to stack series.

        Returns
        -------
        pandas.Series
        """
        series = self.stack()
        series.name = "VALUE"
        return series

    def to_tab(self, mat, mt, e, mu=np.linspace(-1,1,201)):
        """Return tabulated angular distribution for given MAT, MT and energy point.
        """
        from numpy.polynomial import legendre
        sec = self.loc[mat, mt]
        if (e < min(sec.index)) | (e > max(sec.index)): raise NotImplementedError("Energy is out of range")
        if e not in sec.index:
            eg =sorted(set(sec.index) | {e})
            sec = sec.reindex(eg).interpolate(method="slinear")
        coeff = sec.loc[e].dropna()
        c = (coeff.index.values*2+1)/2 * coeff.values
        adistr = legendre.legval(mu, c)
        return pd.Series(adistr, index=mu, name=(mat,mt,e))

    def add_points(self, extra_points):
        """Add additional entries to Lpc incoming energies.
        """
        points = np.array(sorted(extra_points))
        frame = self.copy()
        List = []
        for (mat,mt),df in frame.groupby(["MAT","MT"]):
            rdf = df.loc[mat,mt]
            mask = (points >= min(rdf.index)) & (points <= max(rdf.index))
            grid = sorted((set(rdf.index) | set(points[mask])))
            rdf = rdf.reindex(grid)
            df_notnan = rdf.dropna(axis="columns", thresh=2).interpolate(method='slinear')
            rdf.update(df_notnan)
            rdf = rdf.reset_index()
            rdf["MAT"] = mat
            rdf["MT"] = mt
            rdf = rdf.set_index(["MAT","MT","E"])
            List.append(rdf)
        return Lpc(pd.concat(List, axis=0))

    def perturb(self, pert, method=2, **kwargs):
        """Perturb Legendre polynomials coefficients given a set of perturbations.
        
        Parameters
        ----------
        pert : pandas.Series
            multigroup perturbations from sandy.LpcSamples
        method : int
            * 1 : samples outside the range [0, 2*_mean_] are set to _mean_. 
            * 2 : samples outside the range [0, 2*_mean_] are set to 0 or 2*_mean_ respectively if they fall below or above the defined range.
        
        Returns
        -------
        sandy.Lpc
        """
        frame = self.copy()
        for (mat,mt),_ in self.groupby(["MAT", "MT"]):
            if (mat,mt) not in pert.index: continue
            lpc = frame.loc[mat,mt]
            prt = pert.loc[mat,mt]
            eprt = prt.index.get_level_values("E").unique().values # get cov energies
            elpc = lpc.index.get_level_values("E").unique().values # get lpc energies
            eg = np.array(sorted(set(eprt) | set(elpc)))
            eg = eg[(eg <= max(elpc)) & (eg >= min(elpc))] # limit to lpc range
            lpc_copy = lpc.reindex(eg)
            df_notnan = lpc_copy.dropna(axis="columns", thresh=2) # cut P columns with less than 2 not-NaN
            df_notnan = df_notnan.interpolate(method='slinear')
            lpc_copy.update(df_notnan)
            for l,_ in prt.groupby("L"):
                P = prt.loc[l].reindex(eg).ffill()
                if method == 2:
                    P = P.where(P>0, 0.0)
                    P = P.where(P<2, 2.0)
                elif method == 1:
                    P = P.where((P>0) & (P<2), 1.0)
                lpc_copy[l] *= P
            lpc_copy = lpc_copy.reset_index()
            lpc_copy["MAT"] = mat
            lpc_copy["MT"] = mt
            lpc_copy = Lpc(lpc_copy.set_index(["MAT","MT","E"]))
            frame.update(lpc_copy)
        return Lpc(frame)
    
    def to_tpd(self):
        """Convert Lpc instance to Tpd instance.
        Keep indexes.
        """
        out = pd.DataFrame([self.to_tab(*ix) for ix in self.index], index=self.index)
        return Tpd(out)
        


class Tpd(pd.DataFrame):
    """Tabulated probability distribution for angular distribution of outgoing particles.
    
    Dataframe components
    --------------------
    index :
        - MAT number
        - MT number
        - incoming neutron energy
    
    columns :
        - cosines
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.index.names = ["MAT", "MT", "E"]
        self.columns.name = "COS"
        self.sort_index(inplace=True)



class Edistr(pd.DataFrame):
    """Energy distribution outgoing particles.
    
    Dataframe components
    --------------------
    index :
        - MAT number
        - MT number
        - number of partial energy distribution (k)
        - incoming neutron energy
    
    columns :
        - outgoing neutron energies
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.index.names = ["MAT", "MT", "K", "EIN"]
        self.sort_index(inplace=True)
        self.columns.name = "EOUT"

    def to_stack(self):
        """Convert Edistr instance to stack series.
        
        Returns
        -------
        pandas.Series
        """
        series = self.stack()
        series.name = "VALUE"
        return series

    def add_points(self, extra_points):
        """Add additional entries to Edistr incoming energies.
        
        Parameters
        ----------
        extra_points : iterable
            energy points in eV
        
        Returns
        -------
        sandy.Edistr
        """
        frame = self.copy()
        List = []
        for (mat,mt,k),df in frame.groupby(["MAT","MT","K"]):
            grid = sorted((set(df.loc[mat, mt, k].index) | set(extra_points)))
            df = df.reset_index().set_index("EIN").reindex(grid).interpolate(method='slinear').fillna(0).reset_index()
            df["MAT"] = np.round(df.MAT.values).astype(int)
            df["MT"] = np.round(df.MT.values).astype(int)
            df["K"] = np.round(df.K.values).astype(int)
            df = df.set_index(["MAT","MT","K","EIN"]).sort_index()
            List.append(df)
        return Edistr(pd.concat(List, axis=0))

    def normalize(self):
        """Normalize each outgoing energy distribution to 1.
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
            if (mat,mt) not in pert.index: continue
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
        if normalize:
            return Edistr(frame).normalize()
        return Edistr(frame)



class Fy(pd.DataFrame):
    """Dataset of independent and/or cumulative fission yields and 
    uncertainties for one or more energies and fissioning isotope.
    
    Index
    -----
    MAT : `int`
        MAT number
    MT : `int`
        MT number
    E : `float`
        incoming neutron energy
    ZAM : `int`
        ZZZ * 10000 + AAA * 10 + META
    
    Columns
    -------
    YI : `float`
        fission yields
    DFY : `float`
        fission yield uncertainties
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.index.names = ["MAT", "MT", "E", "ZAM"]
        self.sort_index(inplace=True)
    
    def filter_by(self, index, value):
        """Delete covariances for indices not equal to given value.
        
        Parameters
        ----------
        index : `str`
            index on which to apply the filter, i.e. MAT, MT, E, ZAM
        value :
            corresponding value
        
        Returns
        -------
        `sandy.Fy`
        """
        mask = self.index.get_level_values(index) == value
        df = self.iloc[mask]
        return self.__class__(df)

    def get_cov(self, mat, mt, energy):
        """Extract absolute covariance matrix.
        
        Returns
        -------
        `sandy.FyCov`
        """
        df = self.filter_by("MAT", mat).filter_by("MT", mt).filter_by("E", energy)
        cov = np.diag(df.DYI**2)
        return FyCov(cov, index=df.index, columns=df.index)
    
    def perturb(self, pert, method=2, **kwargs):
        """Perturb fission yields given a set of perturbations.
        
        Parameters
        ----------
        pert : pandas.Series
            perturbations from sandy.FySamples
        method : int
            * 1 : samples outside the range [0, 2*_mean_] are set to _mean_. 
            * 2 : samples outside the range [0, 2*_mean_] are set to 0 or 2*_mean_ respectively if they fall below or above the defined range.
        
        Returns
        -------
        sandy.Fy
        """
        frame = self.copy()
        for mat, dfmat in frame.groupby("MAT"):
            if mat not in pert.index.get_level_values("MAT").unique():
                continue
            for mt, dfmt in dfmat.groupby("MT"):
                if mt not in pert.loc[mat].index.get_level_values("MT").unique():
                    continue
                for e, dfe in dfmt.groupby("E"):
                    if e not in pert.loc[mat, mt].index.get_level_values("E").unique():
                        continue
                    for zam, dfzam in dfe.groupby("ZAM"):
                        if zam not in pert.loc[mat, mt, e].index.get_level_values("ZAM").unique():
                            continue
                        X = dfzam.YI.values
                        P = pert.loc[mat,mt,e,zam]
                        if method == 2:
                            if P < -X:
                                X = 0
                            elif P > X:
                                X = 2*X
                            else:
                                X += P
                        elif method == 1:
                            if P >= -X and P <= X:
                                X += P
                        frame.loc[mat,mt,e,zam]["YI"] = X
        return Fy(frame)



class BaseCov(pd.DataFrame):
    """Base covariance class inheriting from `pandas.DataFrame`.
    Must be used as superclass by all other Nuclear Data Covariance Objects.
    
    Methods
    -------
    corr
        get correlation matrix instance with inherited class type
    eig
        get covariance matrix eigenvalues as a `pandas.Series` instance
    from_list
        extract global cross section/nubar covariance matrix from iterables 
    get_var
        get covariance matrix variances as a `pandas.Series` instance
    get_std
        get covariance matrix standard deviations as a `pandas.Series` instance
    to_matrix
        get covariance matrix as a `sandy.formats.utils.Cov` instance
    """
    
    def to_matrix(self):
        """Extract dataframe values as a `Cov` instance
        
        Returns
        -------
        `sandy.formats.utils.Cov`
            covariance matrix as a `numpy` array
        """
        return Cov(self.values)

    def eig(self):
        """Extract eigenvalues in descending order.
        
        Returns
        -------
        `pandas.Series`
            sorted list of eigenvalues
        """
#        NE = np.extract(E < 0, E)    # extract negative eigenvalues
#        if len(NE) != 0:
#            neig_max = max(abs(NE))
#            eig_max = max(abs(E))
#            if neig_max/eig_max >= 0.1:
#                logging.warn("found large negative eigenvalues")
        E = self.to_matrix().eig()[0]
        return pd.Series(sorted(E, reverse=True), name='eigenvalues')

    def corr(self):
        """Extract correlation matrix.
        
        Returns
        -------
        `BaseCov` or daughter instance
            correlation matrix
        """
        corr = self.to_matrix().corr()
        return self.__class__(corr, index=self.index, columns=self.columns)

    def check_diagonal(self, verbose=True):
        """Check if any of the diagonal elements is negative.
        Return count of negative variances.
        
        Parameters
        ----------
        verbose : `bool`
            If `True` print list of negative variances

        Returns
        -------
        `int`
        """
        var = self.get_var()
        mask = var < 0
        count = mask.sum()
        if verbose and count > 0:
            string = var[mask].to_string()
            logging.warn("found {} negative variances\n{}".format(count, string))
        return count

    def get_var(self):
        """Extract diagonal.

        Returns
        -------
        `pandas.Series`
            array of variances
        """
        return pd.Series(np.diag(self.values), index=self.index, name="VAR")

    def get_std(self):
        """Extract square root of diagonal.
        
        Returns
        -------
        `pandas.Series`
            array of standard deviations
        """
        return self.get_var().apply(np.sqrt).rename("STD")

    def filter_by(self, index_key, index_values, columns_key, columns_values):
        """Filter dataframe based on given index and allowed values.

        .. hint:: use this method to filter the dataframe other than `.loc` as 
                  it returns a `BaseCov` (or daughter) instance.
        
        Parameters
        ----------
        index : `str`
            index on which to apply the filter, e.g. "MAT", "MT"
        values : `iter`
            list of accepted corresponding value
        
        Returns
        -------
        `BaseCov` or daughter instance
        """
        index_cond = self.index.get_level_values(index_key).isin(index_values)
        columns_cond = self.index.get_level_values(columns_key).isin(columns_values)
        df = self.loc[index_cond, columns_cond]
        if df.empty:
            raise SandyError("applied filter returned empty matrix")
        return self.__class__(df)

    def _stack_correlations(self):
        corrstack = self.corr().T.reset_index(drop=True).T.reset_index(drop=True).stack()
        index = self.index.to_flat_index()
        multiindex = pd.MultiIndex.from_product([index.values, index.values])
        corrstack.index = multiindex
        corrstack.index.names = [self.index.names, self.index.names]
        return corrstack

    @classmethod
    def _from_list(cls, iterable):
        """Extract global cross section/nubar covariance matrix from iterables 
        of `EnergyCovs`.
        
        Parameters
        ----------
        iterable : iterable
            list of tuples/lists/iterables with content `[mat, mt, mat1, mt1, EnergyCov]`
        
        Returns
        -------
        `XsCov`
            global cross section/nubar covariance matrix
        """
        columns = ["KEYS_ROWS", "KEYS_COLS", "COV"]
        # Reindex the cross-reaction matrices
        covs = pd.DataFrame.from_records(iterable, columns=columns).set_index(columns[:-1]).COV
        for (keys_rows,keys_cols), cov in covs.iteritems():
            if keys_rows == keys_cols: # diagonal terms
                if cov.shape[0] != cov.shape[1]:
                    raise SandyError("non-symmetric covariance matrix for ({}, {})".format(keys_rows, keys_cols))
                if not np.allclose(cov, cov.T):
                    raise SandyError("non-symmetric covariance matrix for ({}, {})".format(keys_rows, keys_cols))
            else: # off-diagonal terms
                condition1 = (keys_rows,keys_rows) in covs.index
                condition2 = (keys_cols,keys_cols) in covs.index
                if not (condition1 and condition2):
                    covs[keys_rows,keys_cols] = np.nan
                    logging.warn("skip covariance matrix for ({}, {})".format(keys_rows, keys_cols))
                    continue
                ex = covs[keys_rows,keys_rows].index.values
                ey = covs[keys_cols,keys_cols].columns.values
                covs[keys_rows,keys_cols] = cov.change_grid(ex, ey)
        covs.dropna(inplace=True)
        # Create index for global matrix
        rows_levels = covs.index.levels[0]
        indexlist = [(*keys,e) for keys in rows_levels for e in covs[(keys,keys)].index.values]
        index = pd.MultiIndex.from_tuples(indexlist, names=cls.labels)
        # Create global matrix
        matrix = np.zeros((len(index),len(index)))
        for (keys_rows,keys_cols), cov in covs.iteritems():
            ix = index.get_loc(keys_rows)
            ix1 = index.get_loc(keys_cols)
            matrix[ix.start:ix.stop,ix1.start:ix1.stop] = cov
            if keys_rows != keys_cols:
                matrix[ix1.start:ix1.stop,ix.start:ix.stop] = cov.T
#        pdb.set_trace()
        return cls(matrix, index=index, columns=index)



class XsCov(BaseCov):
    """Dataframe to contain cross section/nubar covariance matrices.
    Covariances can be stored for:
        
        - individual reactions,
        - cross reactions,
        - cross isotopes,
        - cross sections vs nubar

    **Index**:
        
        - MAT : (`int`) MAT number to identify the isotope
        - MT : (`int`) MT number to identify the reaction
        - E : (`float`) energy of the incident particle

    **Columns**:
        
        - MAT : (`int`) MAT number to identify the isotope
        - MT : (`int`) MT number to identify the reaction
        - E : (`float`) energy of the incident particle
    
    **Values**: matrix coefficients
    
    Methods
    -------
    from_endf6
        extract global cross section/nubar covariance matrix from 
        `sandy.formats.endf6.Endf6` instance
    from_errorr
        extract global cross section/nubar covariance matrix from 
        `sandy.formats.errorr.Errorr` instance  
        of `sandy.formats.utils.EnergyCov` instances
    get_samples
        extract perturbations from global cross section/nubar covariance matrix
    get_section
        extract section of global cross section/nubar covariance matrix as a 
        `sandy.formats.utils.EnergyCov` instance
    """

    labels = ["MAT", "MT", "E"]

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.index.names = self.labels
        self.columns.names = self.labels

    def get_samples(self, nsmp, eig=0, seed=None):
        cov = self.to_matrix()
        frame = pd.DataFrame(cov.sampling(nsmp, seed=seed) + 1, index=self.index, columns=range(1,nsmp+1))
        frame.columns.name = 'SMP'
        if eig > 0:
            eigs = cov.eig()[0]
            idxs = np.abs(eigs).argsort()[::-1]
            dim = min(len(eigs), eig)
            eigs_smp = Cov(np.cov(frame.values)).eig()[0]
            idxs_smp = np.abs(eigs_smp).argsort()[::-1]
            print("MF[31,33] eigenvalues:\n{:^10}{:^10}{:^10}".format("EVAL", "SAMPLES","DIFF %"))
            diff = div0(eigs[idxs]-eigs_smp[idxs_smp], eigs[idxs], value=np.NaN)*100.
            E = ["{:^10.2E}{:^10.2E}{:^10.1F}".format(a,b,c) for a,b,c in zip(eigs[idxs][:dim], eigs_smp[idxs_smp][:dim], diff[:dim])]
            print("\n".join(E))
        return frame
    
    def get_section(self, mat, mt, mat1, mt1):
        """Extract section of the global covariance/correlation matrix.
        A section is defined by a unique combination of MAT/MT and MAT1/MT1 numbers.
        
        Parameters
        ----------
        mat : `int`
            MAT number for index
        mt : `int`
            MAT number for index
        mat1 : `int`
            MAT number for columns
        mt1 : `int`
            MT number for columns
        
        Returns
        -------
        `EnergyCov`
            section of the global covariance matrix
        """
        df = self.loc[(mat,mt), (mat1,mt1)]
        return EnergyCov(df)
    
    def _change_energy_grid(self, mat, mt, new_grid):
        df = self.index.to_frame(index=False)
        listdf = []
        for (mat_,mt_),edf in df.groupby(["MAT","MT"]):
            if mat_ == mat and mt_ == mt:
                edf = pd.MultiIndex.from_product([[mat],[mt],new_grid], names=["MAT","MT","E"]).to_frame(index=False)
            listdf.append(edf)
        df = pd.concat(listdf, ignore_index=True)
        index = df.set_index(['MAT', 'MT', "E"]).index
        cov = self.reindex(index=index, method="ffill").reindex(columns=index, method="ffill").fillna(0)
        return self.__class__(cov)

    @classmethod
    def from_endf6(cls, endf6):
        """Extract cross section/nubar covariance from ```Endf6``` instance.
        
        Parameters
        ----------
        endf6 : `Endf6`
            `Endf6` instance containing covariance sections
        
        Returns
        -------
        `XsCov`
            global xs/nubar covariance matrix from ENDF6 file
        """
        tape = endf6.filter_by(listmf=[31,33])
        data = []
        # Loop MF/MT
        logging.debug("found {} covariance sections".format(len(tape)))
        for (mat,mf,mt), text in tape.TEXT.iteritems():
            X = tape.read_section(mat, mf, mt)
            # Loop subsections
            logging.debug("reading section MAT={}/MF={}/MT={}".format(mat, mf, mt))
            logging.debug("found {} subsections".format(len(X["SUB"])))
            for sub in X["SUB"].values():
                mat1 = sub['MAT1'] if sub['MAT1'] != 0 else mat
                mt1 = sub['MT1']
                logging.debug("\treading subsection MAT1={}/MT1={}".format(mat1, mt1))
                logging.debug("\tfound {} NI-type sub-subsection".format(len(sub["NI"])))
                covs = []
                # Loop NI-type covariances
                for i,nisec in sub["NI"].items():
                    logging.debug("\t\treconstruct covariance from NI-type section LB={}".format(nisec["LB"]))
                    if nisec["LB"] == 5:
                        foo = EnergyCov.from_lb5_asym if nisec["LS"] == 0 else EnergyCov.from_lb5_sym
                        cov = foo(nisec["EK"], nisec["FKK"])
                        covs.append(cov)
                    elif nisec["LB"] == 1:
                        cov = EnergyCov.from_lb1(nisec["EK"], nisec["FK"])
                        covs.append(cov)
                    elif nisec["LB"] == 2:
                        cov = EnergyCov.from_lb2(nisec["EK"], nisec["FK"])
                        covs.append(cov)
                    elif nisec["LB"] == 6:
                        cov = EnergyCov.from_lb6(nisec["EK"], nisec["EL"], nisec["FKL"])
                        covs.append(cov)
                    else:
                        logging.warn("skip LB={} covariance for [({}/{}), ({}/{})]".format(nisec["LB"], mat, mt, mat1, mt1))
                        continue
                if len(covs) == 0:
                    logging.debug("\tsubsection MAT1={}/MT1={} did not provide accetable covariances".format(mat1, mt1))
                    continue
                cov = EnergyCov.sum_covs(*covs)
                if cov.all().all():
                    logging.warn("\tempty covariance for [({}/{}), ({}/{})]".format(mat, mt, mat1, mt1))
                    continue
                data.append([(mat,mt), (mat1,mt1), cov])
        if not data:
            logging.warn("no xs covariance was found")
            return pd.DataFrame()
        return cls._from_list(data)

    @classmethod
    def from_errorr(cls, errorr):
        """Extract cross section/nubar covariance from `Errorr` instance.
        
        Parameters
        ----------
        errorr : `Errorr`
            `Errorr` instance containing covariance sections
        
        Returns
        -------
        `XsCov`
            global xs/nubar covariance matrix from ERRORR file
        """
        mat = errorr.mat[0]
        eg = errorr.read_section(mat,1,451)["EG"]
        List = []
        for (mat,mf,mt),text in errorr.TEXT.iteritems():
            if mf not in [31, 33]:
                continue
            X = errorr.read_section(mat,mf,mt)
            for mt1,y in X["RP"].items():
                List.append([mat, X["MT"], mat, mt1, y])
        frame = pd.DataFrame(List, columns=('MAT', 'MT','MAT1', 'MT1', 'COV'))
        mi = [(mat,mt,e) for mat,mt in sorted(set(zip(frame.MAT, frame.MT))) for e in eg]
        index = pd.MultiIndex.from_tuples(mi, names=("MAT", "MT", "E"))
        # initialize union matrix
        matrix = np.zeros((len(index),len(index)))
        for i,row in frame.iterrows():
            ix = index.get_loc((row.MAT,row.MT))
            ix1 = index.get_loc((row.MAT1,row.MT1))
            matrix[ix.start:ix.stop-1,ix1.start:ix1.stop-1] = row.COV
        i_lower = np.tril_indices(len(index), -1)
        matrix[i_lower] = matrix.T[i_lower]  # make the matrix symmetric
        return XsCov(matrix, index=index, columns=index)



class EdistrCov(BaseCov):

    labels = ["MAT", "MT", "ELO", "EHI", "EOUT"]

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.index.names = self.labels
        self.columns.names = self.labels
    
    @property
    def nblocks(self):
        """Number of covariance blocks.
        """
        return self.index.get_level_values("EHI").unique().size

    def get_samples(self, nsmp, **kwargs):
        """Draw samples from probability distribution centered in 0 and with
        absolute covariance in EdistrCov instance.
        
        Parameters
        ----------
        nsmp : `int`
            number of samples
        
        Returns
        -------
        `sandy.EdistrSamples`
        """
        cov = self.to_matrix()
        frame = pd.DataFrame(cov.sampling(nsmp), index=self.index, columns=range(1,nsmp+1))
        frame.columns.name = 'SMP'
        if "eig" in kwargs:
            if kwargs["eig"] > 0:
                eigs = cov.eig()[0]
                idxs = np.abs(eigs).argsort()[::-1]
                dim = min(len(eigs), kwargs["eig"])
                eigs_smp = Cov(np.cov(frame.values)).eig()[0]
                idxs_smp = np.abs(eigs_smp).argsort()[::-1]
                print("MF35 eigenvalues:\n{:^10}{:^10}{:^10}".format("EVAL", "SAMPLES","DIFF %"))
                diff = div0(eigs[idxs]-eigs_smp[idxs_smp], eigs[idxs], value=np.NaN)*100.
                E = ["{:^10.2E}{:^10.2E}{:^10.1F}".format(a,b,c) for a,b,c in zip(eigs[idxs][:dim], eigs_smp[idxs_smp][:dim], diff[:dim])]
                print("\n".join(E))
        return EdistrSamples(frame)
    
    def plot_block_corr(self, mat, mt, block, display=True, **kwargs):
        """Plot block correlation matrix.
        
        Parameters
        ----------
        mat : `int`
            MAT number
            
        mt : `int`
            MT number

        block : `int`
            covarianceblock number (starting from 0)

        display : `bool`
            flag to display figure to screen
        
        kwargs : keyword arguments
            extra arguments to pass to ```pcolor```
        
        Returns
        -------
        `matplotlib.pyplot.Axes`
        """
        cov = self.filter_by("MAT", mat).filter_by("MT", mt)
        ehi = cov.index.get_level_values("EHI").unique()[block]
        cov = cov.filter_by("EHI", ehi)
        index = cov.index.get_level_values("EOUT")
        corr = cov.to_matrix()[1].corr
        fig, ax = plt.subplots()
        ax.set_yscale('log')
        ax.set_xscale('log')
        ax.set(xlabel='energy (eV)', ylabel='energy (eV)')
        c = ax.pcolor(index, index, corr, vmin=-1, vmax=1, cmap="bwr", **kwargs)
        fig.colorbar(c, ax=ax)
        if display:
            plt.tight_layout()
            plt.grid()
            plt.show()
            plt.close()
        return ax

    def plot_block_cov(self, mat, mt, block, display=True, **kwargs):
        """Plot block covariance matrix.
        
        Parameters
        ----------
        mat : `int`
            MAT number
            
        mt : `int`
            MT number

        block : `int`
            covarianceblock number (starting from 0)

        display : `bool`
            flag to display figure to screen
        
        kwargs : keyword arguments
            extra arguments to pass to ```pcolor```
        
        Returns
        -------
        `matplotlib.pyplot.Axes`
        """
        cov = self.filter_by("MAT", mat).filter_by("MT", mt)
        ehi = cov.index.get_level_values("EHI").unique()[block]
        cov = cov.filter_by("EHI", ehi)
        index = cov.index.get_level_values("EOUT")
        fig, ax = plt.subplots()
        ax.set_yscale('log')
        ax.set_xscale('log')
        ax.set(xlabel='energy (eV)', ylabel='energy (eV)')
        c = ax.pcolor(index, index, cov, **kwargs)
        fig.colorbar(c, ax=ax)
        if display:
            plt.tight_layout()
            plt.grid()
            plt.show()
            plt.close()
        return ax


class LpcCov(BaseCov):
    """Dataframe to contain Legenre Polynomials coefficients covariance 
    matrices.
    Covariances can be stored for:
        
        - individual polynomial coefficients,
        - cross polynomial coefficients,
        - cross isotopes,

    **Index**:
        
        - MAT : (`int`) MAT number to identify the isotope
        - MT : (`int`) MT number to identify the reaction
        - L : (`int`) polynomial order
        - E : (`float`) energy of the incident particle

    **Columns**:
        
        - MAT : (`int`) MAT number to identify the isotope
        - MT : (`int`) MT number to identify the reaction
        - L : (`int`) polynomial order
        - E : (`float`) energy of the incident particle
    
    **Values**: matrix coefficients
    
    Methods
    -------
    from_endf6
        Extract global cross section/nubar covariance matrix from 
        `sandy.formats.endf6.Endf6` instance
    from_list
        Extract global cross section/nubar covariance matrix from iterables 
        of `sandy.formats.utils.EnergyCov` instances
    get_samples
        Extract perturbations from global cross section/nubar covariance matrix
    get_section
        Extract section of global cross section/nubar covariance matrix as a 
        `sandy.formats.utils.EnergyCov` instance
    """
    
    labels = ["MAT", "MT", "L", "E"]

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.index.names = self.labels
        self.columns.names = self.labels
    
    @classmethod
    def from_endf6(cls, endf6):
        """Extract global Legendre Polynomials coefficients covariance matrix 
        from `sandy.formats.endf6.Endf6`.
        
        Parameters
        ----------
        endf6 : `Endf6`
            `Endf6` instance containing covariance sections
        
        Returns
        -------
        `XsCov`
        """
        tape = endf6.filter_by(listmf=[34])
        data = []
        # Loop MF/MT
        logging.debug("found {} covariance sections".format(len(tape)))
        for (mat,mf,mt), text in tape.TEXT.iteritems():
            X = tape.read_section(mat, mf, mt)
            # Loop subsections
            logging.debug("reading section MAT={}/MF={}/MT={}".format(mat, mf, mt))
            logging.debug("found {} subsections".format(len(X["REAC"])))
            for (mat1,mt1), rsec in X["REAC"].items():
                if mat1 == 0:
                    mat1 = mat
                logging.debug("\treading subsection MAT1={}/MT1={}".format(mat1, mt1))
                logging.debug("\tfound {} P sub-subsection".format(len(rsec["P"])))
                for (l,l1), psec in rsec["P"].items():
                    logging.debug("\treading sub-subsection for (P{},P{})".format(l,l1))
                    logging.debug("\tfound {} NI-type sub-sub-subsection".format(len(psec["NI"])))
                    covs = []
                    for i,nisec in psec["NI"].items():
                        logging.debug("\t\treconstruct covariance from NI-type section LB={}".format(nisec["LB"]))
                        if nisec["LB"] == 5:
                            foo = EnergyCov.from_lb5_asym if nisec["LS"] == 0 else EnergyCov.from_lb5_sym
                            cov = foo(nisec["EK"], nisec["FKK"])
                            covs.append(cov)
                        elif nisec["LB"] == 1:
                            cov = EnergyCov.from_lb1(nisec["EK"], nisec["FK"])
                            covs.append(cov)
                        elif nisec["LB"] == 2:
                            cov = EnergyCov.from_lb2(nisec["EK"], nisec["FK"])
                            covs.append(cov)
                        elif nisec["LB"] == 6:
                            cov = EnergyCov.from_lb6(nisec["EK"], nisec["EL"], nisec["FKL"])
                            covs.append(cov)
                        else:
                            logging.warn("skip LB={} covariance for [({}/{}), ({}/{})]".format(nisec["LB"], mat, mt, mat1, mt1))
                            continue
                    if len(covs) == 0:
                        logging.debug("\tsubsection MAT1={}/MT1={} did not provide accetable covariances".format(mat1, mt1))
                        continue
                    cov = EnergyCov.sum_covs(*covs)
                    if cov.all().all():
                        logging.warn("\tempty covariance for [({}/{}), ({}/{})]".format(mat, mt, mat1, mt1))
                        continue
                    data.append([(mat, mt, l), (mat1, mt1, l1), cov])
        if not data:
            logging.warn("no lpc covariance was found")
            return pd.DataFrame()
        return cls._from_list(data)
#
#                    if len(covs) == 0:
#                        continue
#                    cov = reduce(lambda x, y: x.add(y, fill_value=0).fillna(0), covs).fillna(0)
#                    eg |= set(cov.index.values)
#                    List.append([mat, mt, l, mat1, mt1, l1, cov])
#        if not List:
#            logging.warn("no MF34 covariance found")
#            return pd.DataFrame()
#        frame = pd.DataFrame(List, columns=('MAT', 'MT', 'L', 'MAT1', 'MT1', 'L1', 'COV'))
#        eg = sorted(eg)
#        frame.COV = frame.COV.apply(lambda x:cov_interp(x, eg))
#        # From here, the method is identical to Errorr.get_cov()
#        # Except that the size of eg is equal to the size of each matrix (we include the value for 2e7)
#        # and that the indexes are different
#        MI = [(mat,mt,l,e) for mat,mt,l in sorted(set(zip(frame.MAT, frame.MT, frame.L))) for e in eg]
#        index = pd.MultiIndex.from_tuples(MI, names=("MAT", "MT", "L", "E"))
#        # initialize union matrix
#        matrix = np.zeros((len(index),len(index)))
#        for i,row in frame.iterrows():
#            ix = index.get_loc((row.MAT,row.MT,row.L))
#            ix1 = index.get_loc((row.MAT1,row.MT1,row.L1))
#            matrix[ix.start:ix.stop,ix1.start:ix1.stop] = row.COV
#        i_lower = np.tril_indices(len(index), -1)
#        matrix[i_lower] = matrix.T[i_lower]  # make the matrix symmetric
#        return LpcCov(matrix, index=index, columns=index)
    
    def plot_std(self, display=True, **kwargs):
        """Plot standard deviations with seaborn.
        
        Parameters
        ----------
        display : `bool`
            flag to display figure to screen
        
        kwargs : keyword arguments
            extra arguments to pass to ```seaborn.lineplot```
        
        Returns
        -------
        `matplotlib.pyplot.Axes`
        """
        std = self.get_std()*100
        df = std.to_frame().reset_index()
        df["L"] = df["L"].astype("category")
        palette = list(colors.keys())[:len(df.L.unique())]
        ax = sns.lineplot(data=df, drawstyle="steps-post", x="E", y="STD", hue="L", palette=palette, style="MT", **kwargs)
        ax.set_xscale("log")
        if (df.STD > 200).any():
            ax.set_yscale("log")
        ax.set(xlabel='energy (eV)', ylabel='stdev (%)')
        if display:
            plt.grid()
            plt.show()
            plt.close()
        return ax

    def filter_p(self, p):
        """Delete covariances for Legendre polynomial coefficients with order higher than `p`.
        
        Parameters
        ----------
        p : `int`
            maximum order of Legendre polynomial coefficients
        
        Returns
        -------
        `sandy.LpcCov`
        """
        mask = self.index.get_level_values("L") <= p
        lpccov = self.iloc[mask, mask]
        return LpcCov(lpccov)

    def get_samples(self, nsmp, **kwargs):
        """Draw samples from probability distribution centered in 1 and with
        relative covariance in LpcCov instance.
        
        Parameters
        ----------
        nsmp : `int`
            number of samples
        
        Returns
        -------
        `sandy.LpcSamples`
        """
        cov = self.to_matrix()
        frame = pd.DataFrame(cov.sampling(nsmp) + 1, index=self.index, columns=range(1,nsmp+1))
        if "eig" in kwargs:
            if kwargs["eig"] > 0:
                eigs = cov.eig()[0]
                idxs = np.abs(eigs).argsort()[::-1]
                dim = min(len(eigs), kwargs["eig"])
                eigs_smp = Cov(np.cov(frame.values)).eig()[0]
                idxs_smp = np.abs(eigs_smp).argsort()[::-1]
                print("MF34 eigenvalues:\n{:^10}{:^10}{:^10}".format("EVAL", "SAMPLES","DIFF %"))
                diff = div0(eigs[idxs]-eigs_smp[idxs_smp], eigs[idxs], value=np.NaN)*100.
                E = ["{:^10.2E}{:^10.2E}{:^10.1F}".format(a,b,c) for a,b,c in zip(eigs[idxs][:dim], eigs_smp[idxs_smp][:dim], diff[:dim])]
                print("\n".join(E))
        return LpcSamples(frame)



class FyCov(BaseCov):
    """Absolute covariance matrix for independent/cumulative fission yields.
    
    Index / Columns
    ---------------
    MAT : `int`
        MAT number
    MT : `int`
        MT number
    E : `float`
        incoming neutron energy
    ZAM : `int`
        ZZZ * 10000 + AAA * 10 + META
    """
    
    labels = ["MAT", "MT", "E", "ZAM"]
    
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.index.names = self.labels
        self.columns.names = self.labels

    def get_samples(self, nsmp, eig=0):
        """Draw samples from probability distribution centered in 0 and with
        absolute covariance in FyCov instance.
        
        Parameters
        ----------
        nsmp : `int`
            number of samples
        eig : `int`
            number of eigenvalues to display
        
        Returns
        -------
        `sandy.FySamples`
        """
        cov = self.to_matrix()
        frame = pd.DataFrame(cov.sampling(nsmp), index=self.index, columns=range(1,nsmp+1))
        if eig > 0:
            eigs = cov.eig()[0]
            idxs = np.abs(eigs).argsort()[::-1]
            dim = min(len(eigs), eig)
            eigs_smp = Cov(np.cov(frame.values)).eig()[0]
            idxs_smp = np.abs(eigs_smp).argsort()[::-1]
            print("MF8 eigenvalues:\n{:^10}{:^10}{:^10}".format("EVAL", "SAMPLES","DIFF %"))
            diff = div0(eigs[idxs]-eigs_smp[idxs_smp], eigs[idxs], value=np.NaN)*100.
            E = ["{:^10.2E}{:^10.2E}{:^10.1F}".format(a,b,c) for a,b,c in zip(eigs[idxs][:dim], eigs_smp[idxs_smp][:dim], diff[:dim])]
            print("\n".join(E))
        return FySamples(frame)



class EnergyCov(BaseCov):
    """Dataframe for a multigroup covariance matrix.
    
    **Index**:
        
        - E : (`float`) energy grid for the 1st reaction

    **Columns**:
        
        - E : (`float`) energy grid for the 2nd reaction
    
    **Values**: matrix coefficients
    
    .. note:: It is assumed that the covariance matrix is defined over 
              multi-group energy grids.

              Only 'zero' interpolation is supported.
    
    Methods
    -------
    change_grid
        
    from_lb1
        
    from_lb2
        
    from_lb5_sym
        
    from_lb5_asym
        
    from_lb6
        
    sum_covs
        
    """
    
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.index = pd.Float64Index(self.index, name="E")
        self.columns = pd.Float64Index(self.columns, name="E")
        if list(self.index) != sorted(self.index):
            raise SandyError("index values are not monotonically increasing")
        if list(self.columns) != sorted(self.columns):
            raise SandyError("columns values are not monotonically increasing")
    
    def change_grid(self, ex, ey):
        """Given one energy grid for the x-axis and one energy grid for the 
        y-axis, interpolate/extrapolate the covariance matrix over the new 
        points using the *forward-filling* method.
        
        .. important::
            
            * backward extrapolated values (e.g. below threshold) are replaced by 0
            * forward extrapolated values (e.g. above 20 MeV) are replaced by 
              the covariance coefficient that refers to the last point in the 
              original grid
        
        Parameters
        ----------
        ex : `iterable`
            covariance energy grid for the x-axis (first reaction)
        ey : `iterable`
            covariance energy grid for the y-axis (second reaction)
        
        Returns
        -------
        `sandy.formats.utils.EnergyCov`
            Covariance matrix interpolated over the new axes.
        """
        df = self.reindex(index=ex, method="ffill"). \
                  reindex(columns=ey, method="ffill"). \
                  fillna(0)
        return self.__class__(df)
    
    def _get_mesh(self):
        X, Y = np.meshgrid(self.index.values, self.columns.values)
        return X.T, Y.T
    
    def _plot_matrix(self, xscale='log', yscale='log', cmap='bwr', vmin=-1, vmax=1, **kwargs):
        ax = plt.pcolormesh(*self._get_mesh(), self.values, cmap=cmap, vmin=vmin, vmax=vmax, **kwargs)
        plt.colorbar(ax)
        plt.gca().set_xscale(xscale)
        plt.gca().set_yscale(yscale)

    @classmethod
    def sum_covs(cls, *covs):
        """Sum mutligroup covariance matrices into a single one.
        
        Parameters
        ----------
        covs : `iterable` of `sandy.formats.utils.EnergyCov`
            list of multigroup covariance matrices (axes can be different)
        
        Returns
        -------
        `sandy.formats.utils.EnergyCov`
            Multi-group covariance matrix.
        """
        def foo(x, y):
            ex = sorted(set(x.index.tolist() + y.index.tolist()))
            ey = sorted(set(x.columns.tolist() + y.columns.tolist()))
            x_ = x.change_grid(ex, ey)
            y_ = y.change_grid(ex, ey)
            return cls(x_.add(y_))
        df = reduce(lambda x,y: foo(x,y), covs)
        return cls(df)

    @classmethod
    def from_lb1(cls, evalues, fvalues):
        """Extract square covariance matrix from NI-type sub-subsection data 
        with flag `lb=1`.
        
        Parameters
        ----------
        evalues : `iterable`
            covariance energy grid for both axis
        fvalues : `iterable`
            array of F-values (covriance matrix diagonal)
        
        Returns
        -------
        `sandy.formats.utils.EnergyCov`
            Multi-group covariance matrix.
        """
        cov = np.diag(fvalues)
        return cls(cov, index=evalues, columns=evalues)

    @classmethod
    def from_lb2(cls, evalues, fvalues):
        """Extract square covariance matrix from NI-type sub-subsection data 
        with flag `lb=2`.
        
        Parameters
        ----------
        evalues : `iterable`
            covariance energy grid for both axis
        fvalues : `iterable`
            array of F-values
        
        Returns
        -------
        `sandy.formats.utils.EnergyCov`
            Multi-group covariance matrix.
        """
        f = np.array(fvalues)
        cov = f*f.reshape(-1,1)
        return cls(cov, index=evalues, columns=evalues)

    @classmethod
    def from_lb5_sym(cls, evalues, fvalues):
        """Extract square symmetric covariance matrix from NI-type sub-subsection data 
        with flag `lb=5`.
        
        Parameters
        ----------
        evalues : `iterable`
            covariance energy grid for both axis
        fvalues : `iterable`
            array of F-values (flattened upper triangular matrix coefficients)
        
        Returns
        -------
        `sandy.formats.utils.EnergyCov`
            Multi-group covariance matrix.
        """
        ne = len(evalues)
        cov = np.zeros([ne - 1, ne - 1])
        indices = np.triu_indices(ne - 1)
        cov[indices] = np.array(fvalues)
        cov += np.triu(cov, 1).T
        # add zero row and column at the end of the matrix
        cov = np.insert(cov, cov.shape[0], [0]*cov.shape[1], axis=0)
        cov = np.insert(cov, cov.shape[1], [0]*cov.shape[0], axis=1)
        return cls(cov, index=evalues, columns=evalues)

    @classmethod
    def from_lb5_asym(cls, evalues, fvalues):
        """Extract square asymmetric covariance matrix from NI-type sub-subsection data 
        with flag `lb=5`.
        
        Parameters
        ----------
        evalues : `iterable`
            covariance energy grid for both axis
        fvalues : `iterable`
            array of F-values (flattened full matrix)
        
        Returns
        -------
        `sandy.formats.utils.EnergyCov`
            Multi-group covariance matrix.
        """
        ne = len(evalues)
        cov = np.array(fvalues).reshape(ne - 1, ne - 1)
        # add zero row and column at the end of the matrix
        cov = np.insert(cov, cov.shape[0], [0]*cov.shape[1], axis=0)
        cov = np.insert(cov, cov.shape[1], [0]*cov.shape[0], axis=1)
        return cls(cov, index=evalues, columns=evalues)

    @classmethod
    def from_lb6(cls, evalues_r, evalues_c, fvalues):
        """Extract covariance matrix from NI-type sub-subsection data 
        with flag `lb6`.
        
        Parameters
        ----------
        evalues_r : `iterable`
            covariance energy grid for row axis
        evalues_c : `iterable`
            covariance energy grid for column axis
        fvalues : `iterable`
            array of F-values (flattened full matrix)
        
        Returns
        -------
        `sandy.formats.utils.EnergyCov`
            Multi-group covariance matrix.
        """
        ner = len(evalues_r)
        nec = len(evalues_c)
        cov = np.array(fvalues).reshape(ner-1, nec-1)
        # add zero row and column at the end of the matrix
        cov = np.insert(cov, cov.shape[0], [0]*cov.shape[1], axis=0)
        cov = np.insert(cov, cov.shape[1], [0]*cov.shape[0], axis=1)
        return cls(cov, index=evalues_r, columns=evalues_c)

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

class FySystem(pd.DataFrame):
    """Dataset of fission yields and uncertainties for a single fissioning 
    system.

    Index
    -----
    ZAM : `int`
        ZZZ * 10000 + AAA * 10 + META
    
    Columns
    -------
    YI : `float`
        fission yields
    DFY : `float`
        fission yield uncertainties
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.index.name = "ZAM"
        self.columns.name = ["YI", "DYI"]
        self.sort_index(inplace=True)
    
    @property
    def acn(self):
        return self.A.values.dot(self.YI.values)

    @property
    def zcn(self):
        return self.Z.values.dot(self.YI.values)

    @property
    def sum_yields(self):
        return self.YI.sum()
    
    def _get_charge_sensitivity(self):
        return self.Z.values

    def _get_mass_sensitivity(self):
        return self.A.values

    def _get_sum_sensitivity(self):
        return np.array([1]*len(self))
    
    def cov_generator(self, mass, charge):
        """Run GLS adjustment to given fys and uncertainties.
        """
        _be = np.array(self.YI.values)
        _cov = np.diag(self.DYI)
        _be, _cov = gls(_be, _cov, self._get_charge_sensitivity(), charge, 1e-3)
        _be, _cov = gls(_be, _cov, self._get_mass_sensitivity(), mass, 1e-3)
        _be, _cov = gls(_be, _cov, self._get_sum_sensitivity(), 2, 1e-3)
        _be, _cov = gls(_be, _cov, self._get_chain_sensitivity(), chain, cov_chain)
        return _be, _cov



class Cov(np.ndarray):
    """Covariance matrix treated as a `numpy.ndarray`.
    
    Methods
    -------
    corr
        extract correlation matrix
    corr2cov
        produce covariance matrix given correlation matrix and standard deviation 
        array
    eig
        get covariance matrix eigenvalues and eigenvectors
    get_L
        decompose and extract lower triangular matrix
    sampling
        draw random samples
    """
    
    def __new__(cls, arr):
        obj = np.ndarray.__new__(cls, arr.shape, float)
        obj[:] = arr[:]
        if not obj.ndim == 2:
            raise SandyError("covariance matrix must have two dimensions")
        if not np.allclose(obj, obj.T):
            raise SandyError("covariance matrix must be symmetric")
        if (np.diag(arr) < 0).any():
            raise SandyError("covariance matrix must have positive variances")
        return obj

    @classmethod
    def corr2cov(cls, corr, std):
        """Extract `Cov` instance given correlation matrix and standard 
        deviation array.
        
        Parameters
        ----------
        corr : `np.array`
            square 2D correlation matrix
        
        std : `np.array`
            array of standard deviations

        Returns
        -------
        `sandy.formats.utils.Cov`
            covariance matrix
        """
        _corr = cls(corr)
        _std = std.flatten()
        dim = _corr.shape[0]
        S = np.repeat(_std, dim).reshape(dim, dim)
        cov = S.T * (_corr * S)
        return cls(cov)

    @staticmethod
    def _up2down(self):
        U = np.triu(self)
        L = np.triu(self, 1).T
        C = U + L
        return C

    def eig(self):
        """Extract eigenvalues and eigenvectors.
        
        Returns
        -------
        `Pandas.Series`
            real part of eigenvalues sorted in descending order
        `np.array`
            matrix of eigenvectors
        """
        E, V = sp.linalg.eig(self)
        E, V = E.real, V.real
        return E, V

    def corr(self):
        """Extract correlation matrix.
        
        .. note:: zeros on the covariance matrix diagonal are translated 
                  into zeros also on the the correlation matrix diagonal.
        
        Returns
        -------
        `sandy.formats.utils.Cov`
            correlation matrix
        """
        std = np.sqrt(np.diag(self))
        with np.errstate(divide='ignore', invalid='ignore'):
            coeff = np.true_divide( 1, std )
            coeff[ ~ np.isfinite(coeff)] = 0  # -inf inf NaN
        corr = np.multiply(np.multiply(self.T, coeff).T, coeff)
        return self.__class__(corr)

    def _reduce_size(self):
        nonzero_idxs =  np.flatnonzero(np.diag(self))
        cov_reduced = self[nonzero_idxs][:,nonzero_idxs]
        return nonzero_idxs, cov_reduced

    @classmethod
    def _restore_size(cls, nonzero_idxs, cov_reduced, dim):
        cov = Cov(np.zeros((dim, dim)))
        for i,ni in enumerate(nonzero_idxs):
            cov[ni,nonzero_idxs] = cov_reduced[i]
        return cov

    def sampling(self, nsmp, seed=None):
        """Extract random samples from the covariance matrix, either using
        the cholesky or the eigenvalue decomposition.

        Parameters
        ----------
        nsmp : `int`
            number of samples
        seed : `int`
            seed for the random number generator (default is `None`)

        Returns
        -------
        `np.array`
            2D array of random samples with dimension `(self.shape[0], nsmp)`
        """
        logging.debug("covariance matrix dimension is {} X {}".format(*self.shape))
        dim = self.shape[0]
        np.random.seed(seed=seed)
        y = np.random.randn(dim, nsmp)
        nonzero_idxs, cov_reduced = self._reduce_size()
        L_reduced = cov_reduced.get_L()
        L = self.__class__._restore_size(nonzero_idxs, L_reduced, dim)
        samples = np.array(L.dot(y))
        return samples

    def get_L(self):
        """Extract lower triangular matrix `L` for which `L*L^T == self`.
        
        Returns
        -------
        `np.array`
            lower triangular matrix
        """
        try:
            L = sp.linalg.cholesky(self, lower=True, overwrite_a=False, check_finite=False)
        except np.linalg.linalg.LinAlgError:
            E, V = self.eig()
            E[E<=0] = 0
            Esqrt = np.diag(np.sqrt(E))
            M = V.dot(Esqrt)
            Q, R = sp.linalg.qr(M.T)
            L = R.T
        return L



def corr2cov(corr, s):
    dim = corr.shape[0]
    S = np.repeat(s, dim).reshape(dim, dim)
    cov = S.T * (corr * S)
    return cov

def triu_matrix(arr, size):
    """
    Given the upper triangular values of a **square symmetric** matrix in
    an array, return the full matrix.

    Inputs:
        - arr :
            (1d array) array with the upper triangular values of the matrix
        - size :
            (int) dimension of the matrix

    Outputs:
        - matrix :
            (2d array) reconstructed 2d-array with symmetric matrix
    """
    matrix = np.zeros([size, size])
    indices = np.triu_indices(size)
    matrix[indices] = arr
    matrix += np.triu(matrix, 1).T
    return matrix
