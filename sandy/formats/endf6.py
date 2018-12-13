# -*- coding: utf-8 -*-
"""
Created on Mon Jan 16 18:03:13 2017

@author: lfiorito
"""
import pdb
import os
import logging
from collections import Counter
from functools import reduce

import numpy as np
import pandas as pd

from sandy.formats.utils import (
        BaseFile,
        Xs,
        Edistr,
        Lpc,
        Fy,
        XsCov,
        EdistrCov,
        LpcCov,
        triu_matrix,
        corr2cov,
        )
from sandy.settings import SandyError
from sandy.functions import find_nearest
from sandy.utils.utils import TimeDecorator


__author__ = "Luca Fiorito"
__all__ = ["Endf6"]

#def split_endf(text):
#    """
#    Read ENDF-6 formatted file and split it into columns based on field widths:
#        C1 C2 L1 L2 N1 N2 MAT MF MT
#        11 11 11 11 11 11  4   2  3.
#    Store list in dataframe.
#    """
#    from io import StringIO
#    def read_float(x):
#        try:
#            return float(x[0] + x[1:].replace('+', 'E+').replace('-', 'E-'))
#        except:
#            return x
#    widths = [11,11,11,11,11,11,4,2,3]
#    columns = ["C1", "C2", "L1", "L2", "N1", "N2","MAT", "MF", "MT"]
#    converters = dict(zip(columns[:6],[read_float]*6))
#    frame =  pd.read_fwf(StringIO(text), widths=widths, names=columns, converters=converters)
#    return frame.query("MAT>0 & MF>0 & MT>0")
#
#


class Endf6(BaseFile):

    Format = "endf6"

    def read_section(self, mat, mf, mt):
        """ Parse MAT/MF/MT section.
        """
        if mf == 1:
            from .mf1 import read
        elif mf == 3:
            from .mf3 import read
        elif mf == 5:
            from .mf5 import read
        elif mf == 4:
            from .mf4 import read
        elif mf == 8:
            from .mf8 import read
        elif mf == 33 or mf == 31:
            from .mf33 import read
        elif mf == 34:
            from .mf34 import read
        elif mf == 35:
            from .mf35 import read
        else:
            raise SandyError("SANDY cannot parse section MAT{}/MF{}/MT{}".format(mat,mf,mt))
        if (mat,mf,mt) not in self.index:
            raise SandyError("section MAT{}/MF{}/MT{} is not in tape".format(mat,mf,mt))
        return read(self.loc[mat,mf,mt].TEXT)

    def write_string(self, title=" "*66, skip_title=False, skip_fend=False):
        """Collect all rows in `Endf6` and write them into string.
        
        Parameters
        ----------
        title : `str`
            title of the file
        skip_title : `bool`
            do not write the title
        skip_fend : `bool`
            do not write the last FEND line

        Returns
        -------
        `str`
        """
        from .records import write_cont
        tape = self.copy()
        string = ""
        if not skip_title:
            string += "{:<66}{:4}{:2}{:3}{:5}\n".format(title, 1, 0, 0, 0)
        for mat,dfmat in tape.groupby('MAT', sort=True):
            for mf,dfmf in dfmat.groupby('MF', sort=True):
                for mt,dfmt in dfmf.groupby('MT', sort=True):
                    for text in dfmt.TEXT:
                        string += text.encode('ascii', 'replace').decode('ascii')
                    string += "{:<66}{:4}{:2}{:3}{:5}\n".format(*write_cont(*[0]*6), int(mat), int(mf), 0, 99999)
                string += "{:<66}{:4}{:2}{:3}{:5}\n".format(*write_cont(*[0]*6), int(mat), 0, 0, 0)
            string += "{:<66}{:4}{:2}{:3}{:5}\n".format(*write_cont(*[0]*6), 0, 0, 0, 0)
        if not skip_fend:
            string += "{:<66}{:4}{:2}{:3}{:5}".format(*write_cont(*[0]*6), -1, 0, 0, 0)
        return string

    def get_xs(self, listmat=None, listmt=None):
        """ Extract selected cross sections (xs).
        xs are linearized on unique grid.
        Missing points are linearly interpolated (use zero when out of domain).

        Conditions:
            - Interpolation law must be lin-lin
            - No duplicate points on energy grid
        """
        condition = self.index.get_level_values("MF") == 3
        tape = self[condition]
        if listmat is not None:
            conditions = [tape.index.get_level_values("MAT") == x for x in listmat]
            condition = reduce(lambda x,y: np.logical_or(x, y), conditions)
            tape = tape[condition]
        if listmt is not None:
            conditions = [tape.index.get_level_values("MT") == x for x in listmt]
            condition = reduce(lambda x,y: np.logical_or(x, y), conditions)
            tape = tape[condition]
        ListXs = []
        for ix,text in tape.TEXT.iteritems():
            X = self.read_section(*ix)
            xs = pd.Series(X["XS"], index=X["E"], name=(X["MAT"],X["MT"])).rename_axis("E").to_frame()
            duplicates = [x for x, count in Counter(xs.index).items() if count > 1]
            if duplicates:
                raise SandyError('duplicate energy points found for MAT{}/MF{}/MT{}\n'.format(*ix) +
                         '\n'.join(map(str,duplicates)))
            if X['INT'] != [2]:
                raise SandyError('MAT{}/MF{}/MT{} interpolation scheme is not lin-lin'.format(*ix))
            ListXs.append(xs)
        if not ListXs:
            logging.warn("requested cross sections were not found")
            return pd.DataFrame()
        frame = reduce(lambda left,right : pd.merge(left, right, left_index=True, right_index=True, how='outer'), ListXs).sort_index().interpolate(method='slinear', axis=0).fillna(0)
        return Xs(frame)

    def update_xs(self, xsFrame):
        from .mf3 import write
        tape = self.copy()
        mf = 3
        for (mat,mt),xsSeries in xsFrame.iteritems():
            if (mat,mf,mt) not in self.index: continue
            sec = self.read_section(mat,mf,mt)
            # Cut threshold xs
            ethresh = sec["E"][0]
            xsSeries = xsSeries.where(xsSeries.index >= ethresh).dropna()
#            iNotZero = next((i for i,x in enumerate(xsSeries) if x), None)
#            if iNotZero > 0: xsSeries = xsSeries.iloc[iNotZero-1:]
            sec["E"] = xsSeries.index.values
            sec["XS"] = xsSeries.values
            # Assume all xs have only 1 interpolation region and it is linear
            sec["NBT"] = [xsSeries.size]
            sec["INT"] = [2]
            text = write(sec)
            tape.loc[mat,mf,mt].TEXT = text
        return Endf6(tape)

    def get_nubar(self, listmat=None, listmt=None):
        """
        Extract selected nubar.
        nubar are linearized on unique grid.
        Missing points are linearly interpolated (use zero when out of domain).

        Conditions:
            - Interpolation law must be lin-lin
            - No duplicate points on energy grid
        """
        condition = self.index.get_level_values("MF") == 1
        tape = self[condition]
        conditions = [tape.index.get_level_values("MT") == x for x in [452, 455, 456]]
        condition = reduce(lambda x,y: np.logical_or(x, y), conditions)
        tape = tape[condition]
        if listmat is not None:
            conditions = [tape.index.get_level_values("MAT") == x for x in listmat]
            condition = reduce(lambda x,y: np.logical_or(x, y), conditions)
            tape = tape[condition]
        if listmt is not None:
            conditions = [tape.index.get_level_values("MT") == x for x in listmt]
            condition = reduce(lambda x,y: np.logical_or(x, y), conditions)
            tape = tape[condition]
#        query = "MF==1 & (MT==452 | MT==455 | MT==456)"
#        if listmat is not None:
#            query_mats = " | ".join(["MAT=={}".format(x) for x in listmat])
#            query += " & ({})".format(query_mats)
#        if listmt is not None:
#            query_mts = " | ".join(["MT=={}".format(x) for x in listmt])
#            query += " & ({})".format(query_mts)
#        tape = self.query(query)
        ListXs = []
        for ix,text in tape.TEXT.iteritems():
            X = self.read_section(*ix)
            xs = pd.Series(X["NUBAR"], index=X["E"], name=(X["MAT"],X["MT"])).rename_axis("E").to_frame()
            duplicates = [x for x, count in Counter(xs.index).items() if count > 1]
            if duplicates:
                raise SandyError('duplicate energy points found for MAT{}/MF{}/MT{}\n'.format(*ix) +
                         '\n'.join(map(str,duplicates)))
            if X['INT'] != [2]:
                raise SandyError('MAT{}/MF{}/MT{} interpolation scheme is not lin-lin'.format(*ix))
            ListXs.append(xs)
        if not ListXs:
            logging.warn("no fission neutron multiplicity was found")
            return pd.DataFrame()
        frame = reduce(lambda left,right : pd.merge(left, right, left_index=True, right_index=True, how='outer'), ListXs).sort_index().interpolate(method='slinear', axis=0).fillna(0)
        return Xs(frame)

    def update_nubar(self, xsFrame):
        from .mf1 import write
        tape = self.copy()
        mf = 1
        for (mat,mt),S in xsFrame.iteritems():
            if (mat,mf,mt) not in self.index: continue
            sec = self.read_section(mat,mf,mt)
            # Cut threshold xs
            iNotZero = next((i for i,x in enumerate(S) if x), None)
            if iNotZero > 0: S = S.iloc[iNotZero-1:]
            sec["E"] = S.index.values
            sec["NUBAR"] = S.values
            # Assume all xs have only 1 interpolation region and it is linear
            sec["NBT"] = [S.size]
            sec["INT"] = [2]
            text = write(sec)
            tape.loc[mat,mf,mt].TEXT = text
        return Endf6(tape)

    def get_edistr(self, listmat=None, listmt=None, verbose=True):
        """
        Extract chi cov for all MAT,MT,SUB found in tape.
        Return a df with MAT,MT,SUB as index and COV as value
        Each COV is a df with Ein on rows and Eout on columns.
        """
        condition = self.index.get_level_values("MF") == 5
        tape = self[condition]
        if listmat is not None:
            conditions = [tape.index.get_level_values("MAT") == x for x in listmat]
            condition = reduce(lambda x,y: np.logical_or(x, y), conditions)
            tape = tape[condition]
        if listmt is not None:
            conditions = [tape.index.get_level_values("MT") == x for x in listmt]
            condition = reduce(lambda x,y: np.logical_or(x, y), conditions)
            tape = tape[condition]
        edistr_list = []
        for ix,text in tape.TEXT.iteritems():
            X = self.read_section(*ix)
            for k,pdistr in X["PDISTR"].items():
                if pdistr["LF"] != 1:
                    if verbose: logging.warn("non-tabulated distribution for MAT{}/MF{}/MT{}, subsec {}".format(*ix,k))
                    continue
                if list(filter(lambda x:x["INT"] != [2], pdistr["EIN"].values())):
                    if verbose: logging.warn("found non-linlin interpolation, skip energy distr. for MAT{}/MF{}/MT{}, subsec {}".format(*ix,k))
                    continue
                for ein,v in sorted(pdistr["EIN"].items()):
                    columns = pd.MultiIndex.from_tuples([(X["MAT"], X["MT"], k, ein)], names=("MAT", "MT", "K", "EIN"))
                    df = pd.DataFrame(v["EDISTR"], index=v["EOUT"], columns=columns)
                    df.index.name = "EOUT"
                    edistr_list.append(df)
        if not edistr_list:
            logging.warn("no tabulated energy distribution was found")
            return pd.DataFrame()
        frame = reduce(lambda x,y : pd.merge(x, y, left_index=True, right_index=True), edistr_list)
        frame = frame.sort_index().interpolate(method="slinear").fillna(0)
        return Edistr(frame.T)
#        for ix,text in tape.TEXT.iteritems():
#            X = self.read_section(*ix)
#            for k,pdistr in X["PDISTR"].items():
#                if pdistr["LF"] != 1:
#                    if verbose: logging.warn("non-tabulated distribution for MAT{}/MF{}/MT{}, subsec {}".format(*ix,k))
#                    continue
#                if list(filter(lambda x:x["INT"] != [2], pdistr["EIN"].values())):
#                    if verbose: logging.warn("found non-linlin interpolation, skip energy distr. for MAT{}/MF{}/MT{}, subsec {}".format(*ix,k))
#                    continue
#                for ein,v in sorted(pdistr["EIN"].items()):
#                    DictEdistr.update({(X["MAT"], X["MT"], k, ein) : pd.Series(v["EDISTR"], index=v["EOUT"])})
#        if not DictEdistr:
#            logging.warn("no tabulated energy distribution was found")
#            return pd.DataFrame()
#        pdb.set_trace()
#        frame = pd.DataFrame.from_dict(DictEdistr, orient='index').interpolate(method="slinear", axis=1).fillna(0)
#        return Edistr(frame)

    def update_edistr(self, edistrFrame):
        from .mf5 import write
        mf = 5
        tape = self.copy()
        for (mat,mt),S in edistrFrame.groupby(["MAT","MT"]):
            if (mat,mf,mt) not in self.index: continue
            sec = self.read_section(mat,mf,mt)
            for k,S in S.groupby(["K"]):
                if sec["PDISTR"][k]["LF"] != 1: continue
                ein_orig = sorted(sec["PDISTR"][k]["EIN"].keys())
                for ein in S.index.get_level_values("EIN"):
                    edistr = S.loc[mat,mt,k,ein].values
                    eout = S.loc[mat,mt,k,ein].index.values
                    ein_found = find_nearest(ein_orig, ein)[1]
                    mask = np.in1d(eout, sec["PDISTR"][k]["EIN"][ein_found]["EOUT"])
                    edistr = edistr[mask]
                    eout = eout[mask]
                    dict_distr = {"EDISTR" : edistr,
                                  "EOUT" : eout,
                                  "NBT" : [len(eout)],
                                  "INT" : [2]}
                    sec["PDISTR"][k]["EIN"].update({ein : dict_distr})
                sec["PDISTR"][k]["NBT_EIN"] = [len(sec["PDISTR"][k]["EIN"])]
                sec["PDISTR"][k]["INT_EIN"] = [2]
            text = write(sec)
            tape.loc[mat,mf,mt].TEXT = text
        return Endf6(tape)

    def get_edistr_cov(self, listmat=None, listmt=None):
        condition = self.index.get_level_values("MF") == 35
        tape = self[condition]
        if listmat is not None:
            conditions = [tape.index.get_level_values("MAT") == x for x in listmat]
            condition = reduce(lambda x,y: np.logical_or(x, y), conditions)
            tape = tape[condition]
        if listmt is not None:
            conditions = [tape.index.get_level_values("MT") == x for x in listmt]
            condition = reduce(lambda x,y: np.logical_or(x, y), conditions)
            tape = tape[condition]
#        query = "MF==35"
#        if listmat is not None:
#            query_mats = " | ".join(["MAT=={}".format(x) for x in listmat])
#            query += " & ({})".format(query_mats)
#        if listmt is not None:
#            query_mts = " | ".join(["MT=={}".format(x) for x in listmt])
#            query += " & ({})".format(query_mts)
#        tape = self.query(query)
        List = []; eg = set()
        for ix,text in tape.TEXT.iteritems():
            X = self.read_section(*ix)
            mat = X['MAT']; mt = X['MT']
            for sub in X["SUB"].values():
                # Ek grid is one unit longer than covariance.
                Ek = np.array(sub["EK"])
                Fkk = np.array(sub["FKK"])
                NE = sub["NE"]
                cov = triu_matrix(Fkk, NE-1)
                # Normalize covariance matrix dividing by the energy bin.
                dE = 1./(Ek[1:]-Ek[:-1])
                cov = corr2cov(cov, dE)
                # Add zero row and column at the end of the matrix
                cov = np.insert(cov, cov.shape[0], [0]*cov.shape[1], axis=0)
                cov = np.insert(cov, cov.shape[1], [0]*cov.shape[0], axis=1)
                cov = pd.DataFrame(cov, index=Ek, columns=Ek)
                eg |= set(cov.index.values)
                List.append([mat, mt, sub["ELO"], sub["EHI"], cov])
        if not List:
            logging.warn("no energy distribution covariance found")
            return pd.DataFrame()
        frame = pd.DataFrame(List, columns=('MAT', 'MT', 'ELO', 'EHI', 'COV'))
        eg = sorted(eg)
        frame.COV = frame.COV.apply(lambda x:cov_interp(x, eg))
        # From here, the method is identical to Errorr.get_cov()
        # Except that the size of eg is equal to the size of each matrix (we include the value for 2e7)
        # and that the indexes are different
        MI = [(mat,mt,elo,ehi,e) for mat,mt,elo,ehi in sorted(set(zip(frame.MAT, frame.MT, frame.ELO, frame.EHI))) for e in eg]
        index = pd.MultiIndex.from_tuples(MI, names=("MAT", "MT", 'ELO', 'EHI', "EOUT"))
        # initialize union matrix
        matrix = np.zeros((len(index),len(index)))
        for i,row in frame.iterrows():
            ix = index.get_loc((row.MAT,row.MT,row.ELO,row.EHI))
            ix1 = index.get_loc((row.MAT,row.MT,row.ELO,row.EHI))
            matrix[ix.start:ix.stop,ix1.start:ix1.stop] = row.COV
        i_lower = np.tril_indices(len(index), -1)
        matrix[i_lower] = matrix.T[i_lower]  # make the matrix symmetric
        return EdistrCov(matrix, index=index, columns=index)

    def get_lpc(self, listmat=None, listmt=None, verbose=True):
        condition = self.index.get_level_values("MF") == 4
        tape = self[condition]
        if listmat is not None:
            conditions = [tape.index.get_level_values("MAT") == x for x in listmat]
            condition = reduce(lambda x,y: np.logical_or(x, y), conditions)
            tape = tape[condition]
        if listmt is not None:
            conditions = [tape.index.get_level_values("MT") == x for x in listmt]
            condition = reduce(lambda x,y: np.logical_or(x, y), conditions)
            tape = tape[condition]
        DictLpc =  {}
        for ix,text in tape.TEXT.iteritems():
            X = self.read_section(*ix)
            if "LPC" not in X: continue
            if X["LPC"]["INT"] != [2]:
                if verbose:
                    logging.warn("found non-linlin interpolation, skip angular distr. for MAT{}/MF{}/MT{}".format(*ix))
                continue
            for e,v in X["LPC"]["E"].items():
                DictLpc.update({(X["MAT"], X["MT"],e) : pd.Series([1]+v["COEFF"])})
        if not DictLpc:
            logging.warn("no angular distribution in Legendre expansion was found")
            return pd.DataFrame()
        frame = pd.DataFrame.from_dict(DictLpc, orient="index")
        return Lpc(frame)

    def update_lpc(self, lpcFrame):
        from .mf4 import write
        mf = 4
        tape = self.copy()
        for (mat,mt),S in lpcFrame.groupby(["MAT","MT"]):
            if (mat,mf,mt) not in self.index: continue
            sec = self.read_section(mat,mf,mt)
            if "LPC" not in sec: continue
            for e in S.loc[mat,mt].index:
                if e in sec["LPC"]["E"]:
                    T = sec["LPC"]["E"][e]["T"]
                    LT = sec["LPC"]["E"][e]["LT"]
                else:
                    T = LT = 0
                coeff =  S.loc[mat,mt,e].dropna().values[1:]
                if len(coeff) == 0: continue
                dict_distr = {"COEFF" : coeff, "LT" : LT, "T" : T}
                sec["LPC"]["E"].update({e : dict_distr})
            sec["LPC"]["NBT"] = [len(sec["LPC"]["E"])]
            sec["LPC"]["INT"] = [2]
            text = write(sec)
            tape.loc[mat,mf,mt].TEXT = text
        return Endf6(tape)

    def get_lpc_cov(self, listmat=None, listmt=None):
        condition = self.index.get_level_values("MF") == 34
        tape = self[condition]
        if listmat is not None:
            conditions = [tape.index.get_level_values("MAT") == x for x in listmat]
            condition = reduce(lambda x,y: np.logical_or(x, y), conditions)
            tape = tape[condition]
        if listmt is not None:
            conditions = [tape.index.get_level_values("MT") == x for x in listmt]
            condition = reduce(lambda x,y: np.logical_or(x, y), conditions)
            tape = tape[condition]
        List = []; eg = set()
        for ix,text in tape.TEXT.iteritems():
            X = self.read_section(*ix)
            mat = X['MAT']; mt = X['MT']
            for (mat1,mt1),rsec in X["REAC"].items():
                if mat1 == 0: mat1 = mat;
                for (l,l1),psec in rsec["P"].items():
                    covs = []
                    for nisec in psec["NI"].values():
                        if nisec["LB"] == 5:
                            Fkk = np.array(nisec["FKK"])
                            if nisec["LS"] == 0: # to be tested
                                cov = Fkk.reshape(nisec["NE"]-1, nisec["NE"]-1)
                            else:
                                cov = triu_matrix(Fkk, nisec["NE"]-1)
                            # add zero row and column at the end of the matrix
                            cov = np.insert(cov, cov.shape[0], [0]*cov.shape[1], axis=0)
                            cov = np.insert(cov, cov.shape[1], [0]*cov.shape[0], axis=1)
                            e1 = e2 = nisec["EK"]
                        elif nisec["LB"] == 1:
                            cov = np.diag(nisec["FK"])
                            e1 = e2 = nisec["EK"]
                        elif nisec["LB"] == 2:
                            f = np.array(nisec["FK"])
                            cov = f*f.reshape(-1,1)
                            e1 = e2 = nisec["EK"]
                        elif nisec["LB"] == 6:
                            cov = np.array(nisec["FKL"]).reshape(nisec["NER"]-1, nisec["NEC"]-1)
                            # add zero row and column at the end of the matrix
                            cov = np.insert(cov, cov.shape[0], [0]*cov.shape[1], axis=0)
                            cov = np.insert(cov, cov.shape[1], [0]*cov.shape[0], axis=1)
                            e1 = nisec["EK"]
                            e2 = nisec["EL"]
                        else:
                            logging.warn("skipped NI-type covariance with flag LB={} for MAT{}/MF{}/MT{}".format(nisec["LB"], *ix))
                            continue
                        cov = pd.DataFrame(cov, index=e1, columns=e2)
                        covs.append(cov)
                    if len(covs) == 0:
                        continue
                    cov = reduce(lambda x, y: x.add(y, fill_value=0).fillna(0), covs).fillna(0)
                    eg |= set(cov.index.values)
                    List.append([mat, mt, l, mat1, mt1, l1, cov])
        if not List:
            logging.warn("no MF34 covariance found")
            return pd.DataFrame()
        frame = pd.DataFrame(List, columns=('MAT', 'MT', 'L', 'MAT1', 'MT1', 'L1', 'COV'))
        eg = sorted(eg)
        frame.COV = frame.COV.apply(lambda x:cov_interp(x, eg))
        # From here, the method is identical to Errorr.get_cov()
        # Except that the size of eg is equal to the size of each matrix (we include the value for 2e7)
        # and that the indexes are different
        MI = [(mat,mt,l,e) for mat,mt,l in sorted(set(zip(frame.MAT, frame.MT, frame.L))) for e in eg]
        index = pd.MultiIndex.from_tuples(MI, names=("MAT", "MT", "L", "E"))
        # initialize union matrix
        matrix = np.zeros((len(index),len(index)))
        for i,row in frame.iterrows():
            ix = index.get_loc((row.MAT,row.MT,row.L))
            ix1 = index.get_loc((row.MAT1,row.MT1,row.L1))
            matrix[ix.start:ix.stop,ix1.start:ix1.stop] = row.COV
        i_lower = np.tril_indices(len(index), -1)
        matrix[i_lower] = matrix.T[i_lower]  # make the matrix symmetric
        return LpcCov(matrix, index=index, columns=index)

    def update_tpd(self, tpdFrame):
        from .mf4 import write
        mf = 4
        tape = self.copy()
        for (mat,mt),S in tpdFrame.groupby(["MAT","MT"]):
            if (mat,mf,mt) not in self.index: continue
            sec = self.read_section(mat,mf,mt)
            pdb.set_trace()
            if "TAB" in sec: del sec["TAB"]
            if "LPC" in sec: del sec["LPC"]
            sec["TAB"] = {}
            sub = {}
            for e in S.loc[mat,mt].index:
                tab = S.loc[mat,mt,e]
                T = LT = 0 # do not deal with this
                pdb.set_trace()
                distr = {"T" : T, "LT" : LT, "NBT" : [len(tab)], "INT" : [2], "MU" : tab.index, "ADISTR" : tab.values}
                sub.update({e : distr})
            sec["TAB"].update({"E" : sub})
            sec["TAB"].update({"NBT" : [len(sec["TAB"]["E"])]})   
            sec["TAB"].update({"INT" : [2]})   
            text = write(sec)
            tape.loc[mat,mf,mt].TEXT = text
        return Endf6(tape)

    def get_fy(self, listenergy=None, listmat=None, listmt=None):
        """Extract selected fission yields.
        xs are linearized on unique grid.
        """
        condition = self.index.get_level_values("MF") == 8
        tape = self[condition]
        if listmat is not None:
            conditions = [tape.index.get_level_values("MAT") == x for x in listmat]
            condition = reduce(lambda x,y: np.logical_or(x, y), conditions)
            tape = tape[condition]
        if listmt is not None:
            conditions = [tape.index.get_level_values("MT") == x for x in listmt]
            condition = reduce(lambda x,y: np.logical_or(x, y), conditions)
            tape = tape[condition]
        listfy = []
        for ix,text in tape.TEXT.iteritems():
            X = self.read_section(*ix)
            for e,esec in X["E"].items():
                if listenergy is not None:
                    if e not in listenergy:
                        continue
                for fy in esec["FY"].values():
                    zam = fy["ZAFP"]*10 + fy["FPS"]
                    fydict = {"MAT" : ix[0], "MT" : ix[2], "E" : e, "ZAM" : zam, "YI" : fy["YI"], "DYI" : fy["DYI"]}
                    listfy.append(fydict)
        if not listfy:
            logging.warn("requested fission yields were not found")
            return pd.DataFrame()
        frame = pd.DataFrame.from_dict(listfy).set_index(["MAT","MT","E","ZAM"])
        return Fy(frame)

    def update_fy(self, fyFrame):
        """Update fy sections of `Endf6` instance with new data coming from 
        a `Fy` instance.
        
        Parameters
        ----------
        fyFrame : `sandy.Fy`
            tabulated fission yields
        
        Returns
        -------
        `sandy.Endf6`
        """
        from .mf8 import write
        tape = self.copy()
        mf = 8
        for (mat,mt),df in fyFrame.groupby(["MAT","MT"]):
            if (mat,mf,mt) not in self.index:
                continue
            sec = self.read_section(mat,mf,mt)
            for e,esec in sec["E"].items():
                if (mat, mt, e) not in fyFrame.index:
                    continue
                newfy = fyFrame.loc[mat,mt,e]
                for zam, row in newfy.iterrows():
                    if zam in esec["FY"]:
                        sec["E"][e]["FY"][zam]["YI"] = row.YI
            text = write(sec)
            tape.loc[mat,mf,mt].TEXT = text
        return Endf6(tape)

    def update_info(self, descr=None):
        """Update RECORDS item (in DATA column) for MF1/MT451 of each MAT based on the content of the TEXT column.
        """
        from .mf1 import write
        tape = self.copy()
        for mat in sorted(tape.index.get_level_values('MAT').unique()):
            sec = self.read_section(mat,1,451)
            records = pd.DataFrame(sec["RECORDS"], columns=["MF","MT","NC","MOD"]).set_index(["MF","MT"])
            new_records = []
            dfmat=tape.loc[mat]
#            for (mf,mt),text in sorted(tape.loc[mat].query('MT!=451'.format(mat)).TEXT.items()):
            for (mf,mt),text in sorted(dfmat[dfmat.index.get_level_values("MT")!=451].TEXT.items()):
                nc = len(text.splitlines())
                # when copying PENDF sections (MF2/MT152) mod is not present in the dictionary
                try:
                    mod = records.MOD.loc[mf,mt]
                except:
                    mod = 0
                new_records.append((mf,mt,nc,mod))
            if descr is not None:
                sec["TEXT"] = descr
            nc = 4 + len(sec["TEXT"]) + len(new_records) + 1
            mod = records.MOD.loc[1,451]
            new_records = [(1,451,nc,mod)] + new_records
            sec["RECORDS"] = new_records
            text = write(sec)
            tape.loc[mat,1,451].TEXT = text
        return Endf6(tape)

    def delete_cov(self):
        """Delete covariance sections (MF>=30) from Endf6 dataframe.
        """
        tape = self.query("MF<30")
        return Endf6(tape)

    def parse(self):
        mats = self.index.get_level_values("MAT").unique()
        if len(mats) > 1:
            raise NotImplementedError("file contains more than 1 MAT")
        self.mat = self.endf = mats[0]
        if hasattr(self, "tape"):
            self.filename = os.path.basename(self.tape)
        INFO = self.read_section(mats[0], 1 ,451)
        del INFO["TEXT"], INFO["RECORDS"]
        self.__dict__.update(**INFO)
        self.SECTIONS = self.loc[INFO["MAT"]].reset_index()["MF"].unique()
        self.EHRES = 0
        self.THNMAX = - self.EHRES if self.EHRES != 0 else 1.0E6



def cov_interp(df, interp_column, method='zero', axis='both'):
    # interp_column is a list
    frame = df[~df.index.duplicated(keep='first')].T
    frame = frame[~frame.index.duplicated(keep='first')].T
    index = np.unique(list(frame.index) + interp_column)
    columns = np.unique(list(frame.columns) + interp_column)
    if axis in ['rows', 'both']:
        frame = frame.reindex(index).interpolate(method=method).reindex(interp_column)
    if axis in ['columns', 'both']:
        frame = frame.transpose().reindex(columns).interpolate(method=method).reindex(interp_column).transpose()
    return frame.fillna(0)



#def write_decay_data_csv(tape, filename):
#    df = tape.query("MF==8 & MT==457")
#    df['Z'] = df.DATA.apply(lambda x : int(x['ZA']//1000))
#    df['A'] = df.DATA.apply(lambda x : int(x['ZA'] - x['ZA']//1000*1000))
#    df['M'] = df.DATA.apply(lambda x : 'g' if x['LISO'] == 0 else 'm' if x['LISO'] == 1 else 'n')
#    df['HL'] = df.DATA.apply(lambda x : x['HL'])
#    df['DHL'] = df.DATA.apply(lambda x : x['DHL'])
#    df['ELP'] = df.DATA.apply(lambda x : x['E'][0])
#    df['DELP'] = df.DATA.apply(lambda x : x['DE'][0])
#    df['EEM'] = df.DATA.apply(lambda x : x['E'][1])
#    df['DEEM'] = df.DATA.apply(lambda x : x['DE'][2])
#    df['EHP'] = df.DATA.apply(lambda x : x['E'][2])
#    df['DEHP'] = df.DATA.apply(lambda x : x['DE'][2])
#    df[['Z','A','M','HL','DHL','ELP','DELP','EEM','DEEM','EHP','DEHP']].to_csv(filename, index=False)
    
#def get_masses(file):
#    print(file)
#    tape = e6.split2df_byZAM(file)
#    tape["AWR"] = tape["TEXT"].apply(lambda x: float(read_float(x[11:22])))
#    return tape.drop("TEXT", axis=1)
#
#def check_masses(files):
#    return [ get_masses(file) for file in files ]
#
#def check_masses_jeff33(path):
#    """
#    Given a path with ENDF-6 files, compare the masses in each section of each file
#    with those given in the RDD JEFF3.3 file.
#    Produce a csv summary output.
#    """
#    from sandy.data_test import __file__ as td
#    RDD = join(dirname(realpath(td)), r"RDD.jeff33")
#    DfRDD = check_masses([RDD])[0].query("MF==1 & MT==451")#.rename(columns={"AWR": "AWR_RDD"})
#    DfN = pd.concat(check_masses(map(lambda x : join(path, x), listdir(path))))
#    C = DfN.merge(DfRDD.drop(["MF","MT"], axis=1), how="left", on=["ZAM"], suffixes=("_N","_RDD"))
#    C["RATIO"] = (C.AWR_N/C.AWR_RDD).values
#    C["MATCH"] = np.isclose(C["RATIO"], 1, rtol=1e-6)
#    grouped = C[["ZAM","AWR_N"]].set_index("ZAM").groupby("ZAM")
#    df = grouped.agg(['count', 'min', 'max', lambda x:x.value_counts().index[0]])
#    df.columns = df.columns.droplevel(0)
#    df = df.merge(DfRDD[['ZAM','AWR']], how='left', left_index=True, right_on=['ZAM'])
#    df["Z"] = np.floor(df.ZAM/10000)
#    df["A"] = np.floor((df.ZAM - df.Z*10000)/10)
#    df["M"] = np.floor(df.ZAM - df.Z*10000 - df.A*10)
#    df.rename(columns={"<lambda>" : "most freq", "AWR" : "ref"}, inplace=True)
#    df.set_index(['Z','A','M'], inplace=True)
#    df.drop("ZAM", axis=1, inplace=True)
#    df["message"] = (np.isclose(df["min"], df["ref"], rtol=1e-6) & (np.isclose(df["min"], df["ref"], rtol=1e-6)))
#    df["message"] = df.message.replace({True: "OK", False: "WARNING"})
#    df.to_csv("jeff33_masses.csv")
#