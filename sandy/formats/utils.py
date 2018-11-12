# -*- coding: utf-8 -*-
"""
Created on Thu Jun 14 09:19:24 2018

@author: fiorito_l
"""
import logging
import pdb
import os

import pandas as pd
import numpy as np
import scipy as sp
import seaborn as sns
import matplotlib.pyplot as plt

from ..sampling.utils import LpcSamples, EdistrSamples
from ..functions import gls, div0
from ..settings import SandyError, colors

__author__ = "Luca Fiorito"
__all__ = ["BaseFile", "Xs", "Lpc", "Edistr", "XsCov", "EdistrCov", "LpcCov", "Cov", "Fy", "Tpd"]

class Section(dict):
    pass

class BaseFile(pd.DataFrame):

    @classmethod
    def from_file(cls, file, listmat=None, listmf=None, listmt=None):
        """
        Read formatted file and call from_text method.
        """
        with open(file) as f: text = f.read()
        out = cls.from_text(text, listmat=listmat, listmf=listmf, listmt=listmt)
        out.TAPE = os.path.abspath(os.path.realpath(os.path.expandvars(file)))
        out.FILENAME = os.path.basename(out.TAPE)
        return out

    @classmethod
    def from_text(cls, text, empty_err=True, listmat=None, listmf=None, listmt=None):
        """
        Read ENDF-6 formatted file and split it into column based on field width:
            TEXT MAT MF MT
              66   4  2  3
        Store list in dataframe with MultiIndex (MAT,MF,MT).
        """
        from io import StringIO
        tape = pd.read_fwf(
                StringIO(text),
                widths = [66, 4, 2, 3],
                names = ["TEXT", "MAT", "MF", "MT"],
                usecols = ["MAT", "MF", "MT"]
                )
        tape["TEXT"] = text.splitlines(True)
        splitters = tape.loc[(tape.MAT==0) & (tape.MF==0) & (tape.MT==0)].index
        dfs = []; ibeg = 0
        for iend in splitters:
            df = tape[ibeg:iend]
            for (mat,mf,mt),group in df.loc[(tape.MAT>0) & (tape.MF>0) & (tape.MT>0)].groupby(["MAT","MF","MT"]):
                # Select only desired sections
                if listmt is not None and mt not in listmt:
                    continue
                if listmat is not None and mat not in listmat:
                    continue
                if listmf is not None and mf not in listmf:
                    continue
                dfs.append({"MAT" : mat, "MF" : mf, "MT" : mt, "TEXT" : "".join(group.TEXT.values)})
            ibeg = iend
        if not dfs:
            raise SandyError("tape is empty")
        tape = pd.DataFrame.from_dict(dfs).set_index(["MAT","MF","MT"])
        if tape.index.duplicated().any():
            raise SandyError("found duplicate MAT/MF/MT")
        frame = cls(tape)
        if frame.empty and empty_err:
            raise SandyError("tape is empty")
        return frame

    def __init__(self, *args, **kwargs):
        kwargs.update({"columns" : ["TEXT"]})
        super().__init__(*args, **kwargs)
        self.index.names = ['MAT', 'MF', 'MT']
        self.sort_index(level=["MAT","MF","MT"], inplace=True)
    
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
        sandy.Xs
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

    def macs(self, E0=0.0253, Elo=1E-5, Ehi=1E1):
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



class BaseCov(pd.DataFrame):
    """Base covariance class inheriting from `pandas.DataFrame`.
    Must be used as superclass by all other covariances.
    """
    
    def to_matrix(self):
        return self.index, Cov(self.values)

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
        """Extract variances in pandas Series.

        Returns
        -------
        `pandas.Series`
        """
        return pd.Series(np.diag(self.values), index=self.index, name="VAR")

    def get_std(self):
        """Extract standard deviations in pandas Series.
        
        Returns
        -------
        `pandas.Series`
        """
        return self.get_var().apply(np.sqrt).rename("STD")

    def filter_by(self, index, value):
        """Delete covariances for indices not equal to given value.
        
        Parameters
        ----------
        index : `str`
            index on which to apply the filter, e.g. "MAT", "MT"
        value : `int`
            corresponding value
        
        Returns
        -------
        `sandy.LpcCov`, `sandy.EdistrCov`, `sandy.XsCov`
        """
        mask = self.index.get_level_values(index) == value
        cov = self.iloc[mask, mask]
        return self.__class__(cov)



class XsCov(BaseCov):
    """
    columns =  (MATi,MTj) ... (MATm,MTn)
    index = E1, E2, ..., El
    """

    def to_matrix(self):
        return self.index, Cov(self.values)

    def get_samples(self, nsmp, **kwargs):
        index, cov = self.to_matrix()
        frame = pd.DataFrame(cov.sampling(nsmp) + 1, index=index, columns=range(1,nsmp+1))
        frame.columns.name = 'SMP'
        if "eig" in kwargs:
            if kwargs["eig"] > 0:
                eigs = cov.eig()[0]
                idxs = np.abs(eigs).argsort()[::-1]
                dim = min(len(eigs), kwargs["eig"])
                eigs_smp = Cov(np.cov(frame.values)).eig()[0]
                idxs_smp = np.abs(eigs_smp).argsort()[::-1]
                print("MF[31,33] eigenvalues:\n{:^10}{:^10}{:^10}".format("EVAL", "SAMPLES","DIFF %"))
                diff = div0(eigs[idxs]-eigs_smp[idxs_smp], eigs[idxs], value=np.NaN)*100.
                E = ["{:^10.2E}{:^10.2E}{:^10.1F}".format(a,b,c) for a,b,c in zip(eigs[idxs][:dim], eigs_smp[idxs_smp][:dim], diff[:dim])]
                print("\n".join(E))
        return frame

    def macs(self, E0=0.0253, Elo=1E-5, Ehi=1E1):
        from ..integrals.macs import maxw_int
        records = []
        for (mat,mt),sec in self.groupby(["MAT","MT"]):
            C = sec[mat,mt].loc[mat,mt]
            E = set(C.index.values)
            E.update([Elo, Ehi])
            E = np.array(sorted(E))
            E = E[(E >= Elo) & (E <= Ehi)]
            C = C.reindex(E).ffill().fillna(0).T.reindex(E).ffill().fillna(0)
            data = [[E0,
                     E[i],
                     E[i+1],
                     maxw_int(E0, E[i], E[i+1])
                     ] for i in range(len(E)-1)]
            dframe = pd.DataFrame(data, columns=["E0", "E1", "E2", "INT"])
            cond = dframe.E1/E0 >= 1e-5
            D = dframe[cond].INT.sum()
            S = dframe[cond].INT / D
            rvar = S.dot(C.values[:-1,:-1][cond][:,cond].dot(S))
            rstd = np.sqrt(rvar)
            skipped = "{}/{}".format(sum(cond==False), len(dframe))
            records.append([mat, mt, rvar, rstd, D, Elo, Ehi, E0, skipped])
        return pd.DataFrame.from_records(records, columns=["MAT", "MT", "VAR", "STD", "FLUX", "Elo", "Ehi", "E0","SKIPPED"])



class EdistrCov(BaseCov):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.index.names = ["MAT", "MT", "ELO", "EHI", "EOUT"]
        self.columns.names = ["MAT", "MT", "ELO", "EHI", "EOUT"]
    
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
        index, cov = self.to_matrix()
        frame = pd.DataFrame(cov.sampling(nsmp), index=index, columns=range(1,nsmp+1))
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

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.index.names = ["MAT", "MT", "L", "E"]
        self.columns.names = ["MAT", "MT", "L", "E"]

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
        index, cov = self.to_matrix()
        frame = pd.DataFrame(cov.sampling(nsmp) + 1, index=index, columns=range(1,nsmp+1))
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



class Fy(pd.DataFrame):
    """Dataset of independent and/or cumulative fission yields and 
    uncertainties for one or more energies and fissioning isotope.
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.index.names = ["MAT", "MT", "E", "ZA", "META"]

    def get_system(self, mat, mt, e):
        frame = self.loc[mat, mt, e]
        frame = frame.reset_index()
        frame["Z"] = frame.ZA//1000
        frame["A"] = frame.ZA - frame.Z*1000
        frame["ZA"] = frame.ZA*10 + frame.META
        frame = frame.set_index("ZA")
        return FySystem(frame)

class FySystem(pd.DataFrame):
    """Dataset of fission yields and uncertainties for a single fissioning 
    system.
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.index.name = "IZAM"
    
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



def up2down(C):
    """
    Given a covariance matrix in input, copy the upper triangular part to the
    lower triangular part.

    Inputs:
        - C :
            (2d-array) input covariance matrix

    Outputs:
        - C1 :
            (2d-array) output covariance matrix
    """
    U = np.triu(C)
    L = np.triu(C, 1).T
    C1 = U + L
    return C1



def corr2cov(corr, s):
    dim = corr.shape[0]
    S = np.repeat(s, dim).reshape(dim, dim)
    cov = S.T * (corr * S)
    cov = up2down(cov)
    return cov


class Cov(np.ndarray):

    def __new__(cls, arr):
        obj = np.ndarray.__new__(cls, arr.shape, arr.dtype)
        obj[:] = arr[:]
        return obj

    @property
    def prefix(self):
        return "COV : "

    @property
    def dim(self):
        """
        Dimension of the covariance.
        """
        length = self.shape[0]
        return length

    @property
    def var(self):
        r"""
        Variance array.
        """
        var = np.diag(np.array(self))
        return var

    @property
    def std(self):
        r"""
        Standard deviation array.
        """
        var = self.var
        if (var < 0).any():
            raise ValueError("Variances must be non-negative")
        std = np.sqrt(var)
        return std

    @property
    def nnegvar(self):
        r"""
        Number of negative variances.
        """
        return np.flatnonzero(self.var < 0).size

    @property
    def nzerovar(self):
        r"""
        Number of zero variances.
        """
        return np.flatnonzero(self.var == 0).size

    def empty_off_diagonals(self):
        r"""
        Remove off-diagonal elements.

        Outputs:
            - :``C``: :
                (``cov.Cov instance``) covariance with empty off-diagonals
        """
        logging.info(self.prefix + "'no_correlations' option is requested, delete off-diagonal terms")
        C = Cov(np.diag(np.diag(self)))
        return C

    def is_symmetric(self):
        r"""
        Check if covariance is symmetric.

        If it is nearly symmetric (rtol=1e-5), then we copy the upper
        triangular part to the lower triangular part and we make it
        symmetric.

        Outputs:
            - :``check``: :
                (boolean) ``True`` if matrix is symmetric, else ``False``
        """
        check = True
        if not (self.T == self).all():
            check = False
            if np.isclose(self.T, self).all():
                check = True
                self[:] = up2down(self)
        return check

    def reduce_size(self):
        """
        Reduce matrix dimensions when zeros are found on the diagonal.

        Outputs:
            * :``nonzero_idxs``: :
                (1d array) positions of the original diagonal matrix where the
                coefficients were not zero
            * :``cov_reduced``: :
                (``cov.Cov`` instance) reduced covariance matrix

        """
        nonzero_idxs =  np.flatnonzero(np.diag(self))
        cov_reduced = self[nonzero_idxs][:,nonzero_idxs]
        return nonzero_idxs, cov_reduced

    def restore_size(self, nonzero_idxs, cov_reduced):
        """
        Restore original matrix dimensions from a reduced matrix and an array
        of positions to convert from reduced to original size.

        Inputs:
            * :``nonzero_idxs``: :
                (1d array) positions of the original diagonal matrix where the
                coefficients were not zero
            * :``cov_reduced``: :
                (``cov.Cov`` instance) reduced covariance matrix

        Outputs:
            * :``cov``: :
                (``cov.Cov`` instance) reduced covariance matrix increased
                to given size according to the indexes given in input
        """
        cov = Cov(np.zeros_like(self))
        for i,ni in enumerate(nonzero_idxs):
            cov[ni,nonzero_idxs] = cov_reduced[i]
        return cov

    def sampling(self, nsmp, pdf='normal'):
        r"""
        Extract random samples from the covariance matrix, either using
        the cholesky or the eigenvalue decomposition.

        Inputs:
            - :``nsmp``: :
                (integer) number of samples

        Outputs:
            - :``samples``: :
                (array) random samples
        """
        logging.debug("covariance matrix dimension is {} X {}".format(*self.shape))
        y = np.random.randn(self.dim, int(nsmp))
        nonzero_idxs, cov_reduced = self.reduce_size()
        nzeros = self.shape[0] - len(nonzero_idxs)
        if nzeros > 0:
            logging.debug("found {} zeros on the diagonal, reduce matrix dimension to {} X {}".format(nzeros, *cov_reduced.shape))
        try:
            L_reduced = cov_reduced.cholesky()
            logging.debug("cholesky decomposition was successful")
        except np.linalg.linalg.LinAlgError as exc:
            logging.debug("cholesky decomposition was not successful, proceed with eigenvalue decomposition")
            L_reduced = cov_reduced.eigendecomp()
        L = self.restore_size(nonzero_idxs, L_reduced)
        samples = np.array(L.dot(y), dtype=float)
        return samples

    @property
    def corr(self):
        r"""
        Correlation matrix.
        """
        from sandy.functions import div0
        if not self.is_symmetric():
            raise ValueError("Covariance matrix must be square and symmetric")
        coeff = div0(1, self.std)
        corr = np.multiply(np.multiply(self.T, coeff).T, coeff)
        return corr

    def cholesky(self):
        r"""
        Perform a Cholesky decomposition of the covariance matrix.

        Outputs:
            - :``L``: :
                (2d array) lower triangular matrix
        """
        from scipy.linalg import cholesky
        L = cholesky(self, lower=True, overwrite_a=False, check_finite=False)
        return L

    def eig(self):
        r"""
        Extract eigenvalues and eigenvectors of the covariance matrix.

        Outputs:
            - :``E``: :
                (1d-array) eigenvalues
            - :``V``: :
                (2d-array) eigenvectors
        """
        E, V = sp.linalg.eig(self)
        E, V = E.real, V.real
        return E, V

    def eigendecomp(self):
        r"""
        Perform an eigenvalue decomposition of the covariance matrix.

        Outputs:
            - :``L``: :
                (2d-array) lower triangular matrix
        """
        from scipy.linalg import qr
        E, V = self.eig()
        NE = np.extract(E < 0, E)    # extract negative eigenvalues
        if len(NE) != 0:
            neig_max = max(abs(NE))
            eig_max = max(abs(E))
            if neig_max/eig_max >= 0.1:
                logging.warn("found large negative eigenvalues")
#            logging.debug(self.prefix + '{} negative eigenvalues were found and replaced with zero'.format(negative_eig.size))
#            pos = sorted(abs(E),reverse=True).index(largest_negative) + 1
#            logging.debug(self.prefix + 'Largest negative eigenvalue ranks {}/{}'.format(pos, E.size))
#            logging.debug(self.prefix + 'eig(-)/eig_max = {}%'.format(largest_negative/max(abs(E))*100.))
        E[E<=0] = 0
        Esqrt = np.diag(np.sqrt(E))
        M = V.dot(Esqrt)
        Q,R = qr(M.T)
        L = R.T
#        logging.debug(self.prefix + "Eigenvalue decomposition was successful")
        return L

    def plot(self):
        r"""
        Plot covariance matrix as a pseudocolor plot of a 2-D array.
        The colorbar is also added to the figure.
        """
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots()
        pcm = ax.matshow(self.corr, vmin=-1, vmax=1, cmap='bwr', aspect='auto')
        # Resize the plot to make space for the colorbar
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, 0.7, box.height])
        # set labels
        ax.set_title('evaluated correlation matrix')
        ax.set_xlabel('energy (eV)')
        ax.set_ylabel('energy (eV)')
        # Plot the colorbar in desired position
        cbaxes = fig.add_axes([0.85, 0.1, 0.03, 0.8])
        plt.colorbar(pcm, cax=cbaxes)
        plt.show()
#        fig.show()

    def dump(self, fname):
        np.savetxt(fname, self, fmt='%.5e')
