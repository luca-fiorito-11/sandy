# -*- coding: utf-8 -*-
"""
Created on Mon Jan 16 18:03:13 2017

@author: lfiorito
"""
import sys, pdb, os, pytest, logging
import numpy as np
import pandas as pd
from warnings import warn


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


class Endf6(pd.DataFrame):

    Format = "endf6"

    @classmethod
    def from_file(cls, file):
        """
        Read ENDF-6 formatted file and call from_text method.
        """
        with open(file) as f: text = f.read()
        out = cls.from_text(text)
        out.TAPE = os.path.abspath(os.path.realpath(os.path.expandvars(file)))
        out.FILENAME = os.path.basename(out.TAPE)
        return out

    @classmethod
    def from_text(cls, text):
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
        splitters = tape.query("MAT==0 & MF==0 & MT==0").index
        dfs = []; ibeg = 0
        for iend in splitters:
            df = tape[ibeg:iend].query("MAT>0 & MF>0 & MT>0").groupby(["MAT","MF","MT"]).sum().reset_index()
            dfs.append(df)
            ibeg = iend
        tape = pd.concat(dfs).set_index(["MAT","MF","MT"])
        dupl = tape.index.duplicated()
        if dupl.any():
            logging.error("found duplicate MAT/MF/MT")
            sys.exit()
        return cls(tape)

    def __init__(self, *args, **kwargs):
        kwargs.update({"columns" : ["TEXT"]})
        super().__init__(*args, **kwargs)
        self.index.names = ['MAT', 'MF', 'MT']

#    def by_ZAM(self):
#        """
#        Change index from MAT,MF,MT to ZAM,MF,MT.
#        Return a pd.DataFrame instance (not Endf6 instance, because most of the methods do not support ZAM).
#        """
#        tape = self.copy().reset_index()
#        text = tape.query("MF==1 & MT==451").TEXT.iloc[0].splitlines()
#        i = 0
#        A, i = read_cont(text, i)
#        B, i = read_cont(text, i)
#        tape = int(A.C1)*10 + B.L1
#        assert False
#        iso["ZAM"] = iso.TEXT.apply(lambda x: int(float(read_float(x[:11]))*10+int(x[103:114]))).values
#        tape =  tape.merge(iso[["MAT","ZAM"]], how="left", on="MAT").drop("MAT", axis=1).set_index(['ZAM','MF','MT']).sort_index()
#        return tape

    def read_section(self, mat, mf, mt):
        """
        Parse MAT/MF/MT section
        """
        if mf == 1:
            from .MF1 import read
        elif mf == 3:
            from .MF3 import read
        elif mf == 5:
            from .MF5 import read
        elif mf == 4:
            from .MF4 import read
        elif mf == 33 or mf == 31:
            from .MF33 import read
        elif mf == 34:
            from .MF34 import read
        elif mf == 35:
            from .MF35 import read
        else:
            sys.exit("ERROR: SANDY cannot parse section MAT{}/MF{}/MT{}".format(mat,mf,mt))
        if (mat,mf,mt) not in self.index:
            raise NotImplementedError("section MAT{}/MF{}/MT{} is not in tape".format(mat,mf,mt))
        return read(self.loc[mat,mf,mt].TEXT)

    def write_string(self, title=" "*66):
        """
        First update MF1/MT451 dictionary.
        Then, Write TEXT column to string.
        """
        from .records import write_cont
        tape = self.copy()
        string = "{:<66}{:4}{:2}{:3}{:5}\n".format(title, 1, 0, 0, 0)
        for mat,dfmat in tape.groupby('MAT', sort=True):
            for mf,dfmf in dfmat.groupby('MF', sort=True):
                for mt,dfmt in dfmf.groupby('MT', sort=True):
                    for text in dfmt.TEXT:
                        string += text.encode('ascii', 'replace').decode('ascii')
                    string += "{:<66}{:4}{:2}{:3}{:5}\n".format(*write_cont(*[0]*6), int(mat), int(mf), 0, 99999)
                string += "{:<66}{:4}{:2}{:3}{:5}\n".format(*write_cont(*[0]*6), int(mat), 0, 0, 0)
            string += "{:<66}{:4}{:2}{:3}{:5}\n".format(*write_cont(*[0]*6), 0, 0, 0, 0)
        string += "{:<66}{:4}{:2}{:3}{:5}".format(*write_cont(*[0]*6), -1, 0, 0, 0)
        return string

    def to_file(self, file, title=" "*66):
        """
        Write TEXT column to file.
        """
        string = self.to_string(title=title)
        with open(file, 'w', encoding="ascii") as f:
            f.write(string)

    def get_xs(self, listmat=None, listmt=None):
        """
        Extract selected cross sections (xs).
        xs are linearized on unique grid.
        Missing points are linearly interpolated (use zero when out of domain).

        Conditions:
            - Interpolation law must be lin-lin
            - No duplicate points on energy grid
        """
        from .utils import Xs
        from collections import Counter
        from functools import reduce
        query = "MF==3"
        if listmat is not None:
            query_mats = " | ".join(["MAT=={}".format(x) for x in listmat])
            query += " & ({})".format(query_mats)
        if listmt is not None:
            query_mts = " | ".join(["MT=={}".format(x) for x in listmt])
            query += " & ({})".format(query_mts)
        tape = self.query(query)
        ListXs = []
        for ix,text in tape.TEXT.iteritems():
            X = self.read_section(*ix)
            xs = pd.Series(X["XS"], index=X["E"], name=(X["MAT"],X["MT"])).rename_axis("E").to_frame()
            duplicates = [x for x, count in Counter(xs.index).items() if count > 1]
            if duplicates:
                sys.exit('ERROR: duplicate energy points found for MAT{}/MF{}/MT{}\n'.format(*ix) +
                         '\n'.join(map(str,duplicates)))
            if X['INT'] != [2]:
                sys.exit('ERROR: MAT{}/MF{}/MT{} interpolation scheme is not lin-lin'.format(*ix))
            ListXs.append(xs)
        if not ListXs:
            warn(UserWarning("no cross section was found"))
            return pd.DataFrame()
        frame = reduce(lambda left,right : pd.merge(left, right, left_index=True, right_index=True, how='outer'), ListXs).sort_index().interpolate(method='slinear', axis=0).fillna(0)
        return Xs(frame)

    def update_xs(self, xsFrame):
        from .MF3 import write
        tape = self.copy()
        mf = 3
        for (mat,mt),S in xsFrame.iteritems():
            if (mat,mf,mt) not in self.index: continue
            sec = self.read_section(mat,mf,mt)
            # Cut threshold xs
            iNotZero = next((i for i,x in enumerate(S) if x), None)
            if iNotZero > 0: S = S.iloc[iNotZero-1:]
            sec["E"] = S.index.values
            sec["XS"] = S.values
            # Assume all xs have only 1 interpolation region and it is linear
            sec["NBT"] = [S.size]
            sec["INT"] = [2]
            text = write(sec)
            tape.loc[mat,mf,mt].TEXT = text
        return Endf6(tape)

    def get_xs_cov(self, listmat=None, listmt=None):
        from .utils import XsCov, triu_matrix
        from functools import reduce
        query = "(MF==33 | MF==31)"
        if listmat is not None:
            query_mats = " | ".join(["MAT=={}".format(x) for x in listmat])
            query += " & ({})".format(query_mats)
        if listmt is not None:
            query_mts = " | ".join(["MT=={}".format(x) for x in listmt])
            query += " & ({})".format(query_mts)
        tape = self.query(query)
        List = []; eg = set()
        for ix,text in tape.TEXT.iteritems():
            X = self.read_section(*ix)
            mat = X['MAT']; mt = X['MT']
            for sub in X["SUB"].values():
                mat1 = sub['MAT1'] if sub['MAT1'] != 0 else mat;
                mt1 = sub['MT1']
                covs = []
                for i,nisec in sub["NI"].items():
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
                        warn("skipped NI-type covariance with flag LB={} for MAT{}/MF{}/MT{}".format(nisec["LB"], *ix), category=Warning)
                        continue
                    cov = pd.DataFrame(cov, index=e1, columns=e2)
                    covs.append(cov)
                if len(covs) == 0:
                    continue
                # covs > 1 for Fe56
                cov = reduce(lambda x, y: x.add(y, fill_value=0).fillna(0), covs).fillna(0)
                eg |= set(cov.index.values)
                List.append([mat, mt, mat1, mt1, cov])
        if not List:
            warn("no MF[31,33] covariance found", category=Warning)
            return pd.DataFrame()
        frame = pd.DataFrame(List, columns=('MAT', 'MT','MAT1', 'MT1', 'COV'))
        eg = sorted(eg)
        frame.COV = frame.COV.apply(lambda x:cov_interp(x, eg))
        # From here, the method is identical to Errorr.get_cov()
        # Except that the size of eg is equal to the size of each matrix (we include the value for 2e7)
        MI = [(mat,mt,e) for mat,mt in sorted(set(zip(frame.MAT, frame.MT))) for e in eg]
        index = pd.MultiIndex.from_tuples(MI, names=("MAT", "MT", "E"))
        # initialize union matrix
        matrix = np.zeros((len(index),len(index)))
        for i,row in frame.iterrows():
            ix = index.get_loc((row.MAT,row.MT))
            ix1 = index.get_loc((row.MAT1,row.MT1))
            matrix[ix.start:ix.stop,ix1.start:ix1.stop] = row.COV
        i_lower = np.tril_indices(len(index), -1)
        matrix[i_lower] = matrix.T[i_lower]  # make the matrix symmetric
        return XsCov(matrix, index=index, columns=index)

    def get_nubar(self, listmat=None, listmt=None):
        """
        Extract selected nubar.
        nubar are linearized on unique grid.
        Missing points are linearly interpolated (use zero when out of domain).

        Conditions:
            - Interpolation law must be lin-lin
            - No duplicate points on energy grid
        """
        from .utils import Xs
        from collections import Counter
        from functools import reduce
        query = "MF==1 & (MT==452 | MT==455 | MT==456)"
        if listmat is not None:
            query_mats = " | ".join(["MAT=={}".format(x) for x in listmat])
            query += " & ({})".format(query_mats)
        if listmt is not None:
            query_mts = " | ".join(["MT=={}".format(x) for x in listmt])
            query += " & ({})".format(query_mts)
        tape = self.query(query)
        ListXs = []
        for ix,text in tape.TEXT.iteritems():
            X = self.read_section(*ix)
            xs = pd.Series(X["NUBAR"], index=X["E"], name=(X["MAT"],X["MT"])).rename_axis("E").to_frame()
            duplicates = [x for x, count in Counter(xs.index).items() if count > 1]
            if duplicates:
                sys.exit('ERROR: duplicate energy points found for MAT{}/MF{}/MT{}\n'.format(*ix) +
                         '\n'.join(map(str,duplicates)))
            if X['INT'] != [2]:
                sys.exit('ERROR: MAT{}/MF{}/MT{} interpolation scheme is not lin-lin'.format(*ix))
            ListXs.append(xs)
        if not ListXs:
            warn(UserWarning("no fission neutron multiplicity was found"))
            return pd.DataFrame()
        frame = reduce(lambda left,right : pd.merge(left, right, left_index=True, right_index=True, how='outer'), ListXs).sort_index().interpolate(method='slinear', axis=0).fillna(0)
        return Xs(frame)

    def update_nubar(self, xsFrame):
        from .MF1 import write
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
        from .utils import Edistr
        query = "MF==5"
        if listmat is not None:
            query_mats = " | ".join(["MAT=={}".format(x) for x in listmat])
            query += " & ({})".format(query_mats)
        if listmt is not None:
            query_mts = " | ".join(["MT=={}".format(x) for x in listmt])
            query += " & ({})".format(query_mts)
        tape = self.query(query)
        DictEdistr =  {}
        for ix,text in tape.TEXT.iteritems():
            X = self.read_section(*ix)
            for k,pdistr in X["PDISTR"].items():
                if pdistr["LF"] != 1:
                    if verbose: warn(UserWarning("WARNING: non-tabulated distribution for MAT{}/MF{}/MT{}, subsec {}".format(*ix,k)))
                    continue
                if list(filter(lambda x:x["INT"] != [2], pdistr["EIN"].values())):
                    if verbose: warn(UserWarning("WARNING: found non-linlin interpolation, skip energy distr. for MAT{}/MF{}/MT{}, subsec {}".format(*ix,k)))
                    continue
                for ein,v in sorted(pdistr["EIN"].items()):
                    DictEdistr.update({(X["MAT"], X["MT"], k, ein) : pd.Series(v["EDISTR"], index=v["EOUT"])})
        if not DictEdistr:
            warn(UserWarning("no tabulated energy distribution was found"))
            return pd.DataFrame()
        frame = pd.DataFrame.from_dict(DictEdistr, orient='index').interpolate(method="slinear", axis=1).fillna(0)
        return Edistr(frame)

    def update_edistr(self, edistrFrame):
        from .MF5 import write
        mf = 5
        tape = self.copy()
        for (mat,mt),S in edistrFrame.groupby(["MAT","MT"]):
            if (mat,mf,mt) not in self.index: continue
            sec = self.read_section(mat,mf,mt)
            for k,S in S.groupby(["K"]):
                if sec["PDISTR"][k]["LF"] != 1: continue
                for ein in S.index.get_level_values("EIN"):
                    dict_distr = {"EDISTR" : S.loc[mat,mt,k,ein].values,
                                  "EOUT" : S.loc[mat,mt,k,ein].index.values,
                                  "NBT" : [S.loc[mat,mt,k,ein].values.size],
                                  "INT" : [2]}
                    sec["PDISTR"][k]["EIN"].update({ein : dict_distr})
                sec["PDISTR"][k]["NBT_EIN"] = [len(sec["PDISTR"][k]["EIN"])]
                sec["PDISTR"][k]["INT_EIN"] = [2]
            text = write(sec)
            tape.loc[mat,mf,mt].TEXT = text
        return Endf6(tape)

    def get_edistr_cov(self, listmat=None, listmt=None):
        from .utils import EdistrCov, triu_matrix, corr2cov
        query = "MF==35"
        if listmat is not None:
            query_mats = " | ".join(["MAT=={}".format(x) for x in listmat])
            query += " & ({})".format(query_mats)
        if listmt is not None:
            query_mts = " | ".join(["MT=={}".format(x) for x in listmt])
            query += " & ({})".format(query_mts)
        tape = self.query(query)
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
            warn(UserWarning("no energy distribution covariance found"))
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
        from .utils import Lpc
        query = "MF==4"
        if listmat is not None:
            query_mats = " | ".join(["MAT=={}".format(x) for x in listmat])
            query += " & ({})".format(query_mats)
        if listmt is not None:
            query_mts = " | ".join(["MT=={}".format(x) for x in listmt])
            query += " & ({})".format(query_mts)
        tape = self.query(query)
        DictLpc =  {}
        for ix,text in tape.TEXT.iteritems():
            X = self.read_section(*ix)
            if "LPC" not in X: continue
            if X["LPC"]["INT"] != [2]:
                if verbose: warn(UserWarning("found non-linlin interpolation, skip angular distr. for MAT{}/MF{}/MT{}".format(*ix)))
                continue
            for e,v in X["LPC"]["E"].items():
                DictLpc.update({(X["MAT"], X["MT"],e) : pd.Series([1]+v["COEFF"])})
        if not DictLpc:
            warn(UserWarning("no angular distribution in Legendre expansion was found"))
            return pd.DataFrame()
        frame = pd.DataFrame.from_dict(DictLpc, orient="index")
        return Lpc(frame)

    def update_lpc(self, lpcFrame):
        from .MF4 import write
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
        from .utils import triu_matrix, LpcCov
        from functools import reduce
        query = "MF==34"
        if listmat is not None:
            query_mats = " | ".join(["MAT=={}".format(x) for x in listmat])
            query += " & ({})".format(query_mats)
        if listmt is not None:
            query_mts = " | ".join(["MT=={}".format(x) for x in listmt])
            query += " & ({})".format(query_mts)
        tape = self.query(query)
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
                            warn("skipped NI-type covariance with flag LB={} for MAT{}/MF{}/MT{}".format(nisec["LB"], *ix), category=Warning)
                            continue
                        cov = pd.DataFrame(cov, index=e1, columns=e2)
                        covs.append(cov)
                    if len(covs) == 0:
                        continue
                    cov = reduce(lambda x, y: x.add(y, fill_value=0).fillna(0), covs).fillna(0)
                    eg |= set(cov.index.values)
                    List.append([mat, mt, l, mat1, mt1, l1, cov])
        if not List:
            warn("no MF34 covariance found", category=Warning)
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

    def update_info(self):
        """
        Update RECORDS item (in DATA column) for MF1/MT451 of each MAT based on the content of the TEXT column.
        """
        from .MF1 import write
        tape = self.copy()
        for mat in sorted(tape.index.get_level_values('MAT').unique()):
            sec = self.read_section(mat,1,451)
            records = pd.DataFrame(sec["RECORDS"], columns=["MF","MT","NC","MOD"]).set_index(["MF","MT"])
            new_records = []
            for (mf,mt),text in sorted(tape.loc[mat].query('MT!=451'.format(mat)).TEXT.items()):
                nc = len(text.splitlines())
                # when copying PENDF sections (MF2/MT152) mod is not present in the dictionary
                try:
                    mod = records.MOD.loc[mf,mt]
                except:
                    mod = 0
                new_records.append((mf,mt,nc,mod))
            nc = 4 + len(sec["TEXT"]) + len(new_records) + 1
            mod = records.MOD.loc[1,451]
            new_records = [(1,451,nc,mod)] + new_records
            sec["RECORDS"] = new_records
            text = write(sec)
            tape.loc[mat,1,451].TEXT = text
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
        self.EHRES = 0
        self.THNMAX = - self.EHRES if self.EHRES != 0 else 1.0E6



def cov_interp(df, interp_column, method='zero', axis='both'):
    # interp_column is a list
    frame = pd.DataFrame(df.copy(), index=df.index.copy(), columns=df.columns.copy())
    index = np.unique(list(frame.index) + interp_column)
    columns = np.unique(list(frame.columns) + interp_column)
    if axis in ['rows', 'both']:
        frame = frame.reindex(index).interpolate(method=method).reindex(interp_column)
    if axis in ['columns', 'both']:
        frame = frame.transpose().reindex(columns).interpolate(method=method).reindex(interp_column).transpose()
    return frame.fillna(0)



#def extract_mu(tape):
##    NLMAX = 0
#    keys = []
#    ListPc = []
#    for chunk in tape.query('MF==4').DATA:
#        if "LPC" in chunk:
#            # check interpolation
#            for e,sub in sorted(chunk["LPC"]["E"].items()):
#                keys.append((chunk["MAT"], chunk["MT"], e))
##                NPL = len(sub["P"])
##                P[i,:NPL] = sub["P"]
##                NLMAX = max(NLMAX, NPL)
#                ListPc.append(sub["P"])
#    index = pd.MultiIndex.from_tuples(keys, names=("MAT", "MT", "E"))
#    DfPc = pd.DataFrame(ListPc, index=index).fillna(0)
#    DfPc.columns = range(1, DfPc.shape[1]+1)
#    DfPc.columns.name = "LegCoeff"
#    # print warning if exceeded NLMAX
#    return DfPc



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


@pytest.fixture(scope="module")
def testPu9():
    from sandy.data_test import Pu9
    tape = Endf6.from_text("\n".join(Pu9.endf6))
    assert (tape.index.get_level_values("MAT").unique() == 9437).all()
    return tape

@pytest.fixture(scope="module")
def testH1():
    from sandy.data_test import H1
    tape = Endf6.from_text("\n".join(H1.pendf))
    assert (tape.index.get_level_values("MAT").unique() == 125).all()
    return tape

@pytest.fixture(scope="module")
def testFe56():
    from sandy.data_test import Fe56
    tape = Endf6.from_text("\n".join(Fe56.endf6))
    assert (tape.index.get_level_values("MAT").unique() == 2631).all()
    return tape

@pytest.mark.formats
@pytest.mark.endf6
@pytest.mark.info
def test_read_info(testPu9):
    S = testPu9.read_section(9437, 1, 451)
    from .MF1 import write
    text = write(S)
    assert testPu9.TEXT.loc[9437,1,451] == text

@pytest.mark.formats
@pytest.mark.endf6
@pytest.mark.info
def test_update_info(testPu9):
    testPu9.loc[9437,3,1].TEXT = "\n".join(testPu9.loc[9437,3,1].TEXT.splitlines()[:10]) + "\n"
    testPu9 = Endf6(testPu9.drop([(9437,3,102)]))
    new = testPu9.update_info()
    recordsold = testPu9.read_section(9437,1,451)["RECORDS"]
    recordsnew = new.read_section(9437,1,451)["RECORDS"]
    assert (3,102,147,1) in recordsold
    assert (3,102,147,1) not in recordsnew
    assert (3,1,188,1) in recordsold
    assert (3,1,10,1) in recordsnew

@pytest.mark.formats
@pytest.mark.endf6
@pytest.mark.xs
def test_read_xs(testPu9):
    S = testPu9.read_section(9437, 3, 102)
    from .MF3 import write
    text = write(S)
    assert testPu9.TEXT.loc[9437,3,102] == text

@pytest.mark.formats
@pytest.mark.endf6
@pytest.mark.xs
def test_extract_xs(testH1):
    testH1.get_xs(listmat=[125], listmt=[1,2,102,4])
    xs = testH1.get_xs(listmat=[9437])
    assert xs.empty

@pytest.mark.formats
@pytest.mark.endf6
@pytest.mark.xs
def test_reconstruct_xs(testH1):
    xs = testH1.get_xs()
    xs[(125,51)] = 1
    xs[(125,2)] = xs[(125,102)] = xs[(125,1)]
    SUM = xs[(125,1)].values*2 + 1
    rec_xs = xs.reconstruct_sums(drop=False)
    assert (rec_xs[(125,4)].values == 1).all()
    assert np.isclose(rec_xs[(125,1)].values,SUM).all()

@pytest.mark.formats
@pytest.mark.endf6
@pytest.mark.xs
def test_update_xs(testH1):
    xs = testH1.get_xs()
    xs[(125,2)] = 1
    new = testH1.update_xs(xs)
    assert (new.read_section(125,3,2)["XS"] == 1).all()
    assert (testH1.read_section(125,3,2)["XS"] != new.read_section(125,3,2)["XS"]).all()

@pytest.mark.formats
@pytest.mark.endf6
@pytest.mark.xs
@pytest.mark.cov
def test_read_xs_cov(testPu9):
    testPu9.read_section(9437, 33, 1)
    testPu9.read_section(9437, 33, 2)
    testPu9.read_section(9437, 33, 18)
    testPu9.read_section(9437, 31, 456)
    testPu9.read_section(9437, 33, 102)

@pytest.mark.formats
@pytest.mark.endf6
@pytest.mark.xs
@pytest.mark.cov
def test_extract_xs_cov(testPu9):
    testPu9.get_xs_cov()

@pytest.mark.formats
@pytest.mark.endf6
@pytest.mark.nubar
def test_read_nubar(testPu9):
    from .MF1 import write
    for mt in (452,455,456):
        S = testPu9.read_section(9437, 1, mt)
        text = write(S)
        assert testPu9.TEXT.loc[9437,1,mt] == text

@pytest.mark.formats
@pytest.mark.endf6
@pytest.mark.nubar
def test_extract_nubar(testPu9):
    testPu9.get_nubar(listmt=[452,456])

@pytest.mark.formats
@pytest.mark.endf6
@pytest.mark.nubar
def test_reconstruct_nubar(testPu9):
    nubar = testPu9.get_nubar()
    nubar[(9437,455)] = 1
    SUM = nubar[(9437,456)].values + 1
    rec_nubar = nubar.reconstruct_sums()
    assert np.isclose(rec_nubar[(9437,452)].values,SUM).all()

@pytest.mark.formats
@pytest.mark.endf6
@pytest.mark.nubar
def test_update_nubar(testPu9):
    nubar = testPu9.get_nubar()
    nubar[(9437,452)] = 1
    new = testPu9.update_nubar(nubar)
    assert (new.read_section(9437,1,452)["NUBAR"] == 1).all()
    assert (testPu9.read_section(9437,1,452)["NUBAR"] != new.read_section(9437,1,452)["NUBAR"]).all()

@pytest.mark.formats
@pytest.mark.endf6
@pytest.mark.lpc
def test_read_lpc(testFe56):
    from .MF4 import write
    S = testFe56.read_section(2631, 4, 2)
    assert S["LTT"] == 3
    assert "LPC" in S
    assert S["LPC"]["NE"] == len(S["LPC"]["E"]) == 1782
    assert "TAB" in S
    assert S["TAB"]["NE"] == len(S["TAB"]["E"]) == 28
    text = write(S)
    assert testFe56.TEXT.loc[2631,4,2] == text
    for mt in range(51,83):
        S = testFe56.read_section(2631, 4, mt)
        text = write(S)
        assert testFe56.TEXT.loc[2631,4,mt] == text

@pytest.mark.formats
@pytest.mark.endf6
@pytest.mark.lpc
def test_extract_lpc(testFe56):
    testFe56.get_lpc()

@pytest.mark.formats
@pytest.mark.endf6
@pytest.mark.lpc
def test_convert_lpc(testFe56):
    Lpc = testFe56.get_lpc()
    C = Lpc.to_tab(2631, 2, 1e-5)
    assert (C == 0.5).all()
    C = Lpc.to_tab(2631, 2, 1e4)
    assert (C >= 0).all()

@pytest.mark.formats
@pytest.mark.endf6
@pytest.mark.lpc
def test_addpoints_lpc(testFe56):
    lpcold = testFe56.get_lpc(listmt=[2])
    lpcnew = lpcold.add_points([1e4, 1e5, 1e6])
    assert 1e4 in lpcnew.loc[2631,2].index
    assert 1e4 not in lpcold.loc[2631,2].index
    assert 1e5 in lpcnew.loc[2631,2].index
    assert 1e5 not in lpcold.loc[2631,2].index
    assert 1e6 in lpcnew.loc[2631,2].index
    assert 1e6 in lpcold.loc[2631,2].index
    lpcnew = lpcold.add_points([])
    assert (lpcnew.index == lpcold.index).all()
    lpcold = testFe56.get_lpc(listmt=[810])
    lpcnew = lpcold.add_points([1e4])
    assert 1e4 not in lpcnew.loc[2631,810].index
    assert 1e4 not in lpcold.loc[2631,810].index

@pytest.mark.formats
@pytest.mark.endf6
@pytest.mark.lpc
def test_update_lpc(testFe56):
    lpc = testFe56.get_lpc()
    new = testFe56.update_lpc(lpc)
    assert new.TEXT[2631,4,2] == testFe56.TEXT[2631,4,2]
    lpc.loc[2631,2,1e5:2e7] = 0
    new = testFe56.update_lpc(lpc)
    assert new.TEXT[2631,4,2] != testFe56.TEXT[2631,4,2]
    new_sec = new.read_section(2631,4,2)
    assert (np.array(new_sec["LPC"]["E"][2e7]["COEFF"]) == 0).all()

@pytest.mark.formats
@pytest.mark.endf6
@pytest.mark.lpc
@pytest.mark.cov
def test_read_lpc_cov(testFe56):
    S = testFe56.read_section(2631, 34, 2)
    assert S["NMT1"] == len(S["REAC"]) == 1
    assert len(S["REAC"][0,2]["P"]) == S["REAC"][0,2]["NL"]*(S["REAC"][0,2]["NL"]+1)//2
    assert len(S["REAC"][0,2]["P"][1,1]["NI"]) == 3

@pytest.mark.formats
@pytest.mark.endf6
@pytest.mark.lpc
@pytest.mark.cov
def test_extract_lpc_cov(testFe56):
    C = testFe56.get_lpc_cov()
    assert C.index.names == C.columns.names == ['MAT', 'MT', 'L', 'E']
    assert (C.values == C.values.T).all()

@pytest.mark.formats
@pytest.mark.endf6
@pytest.mark.chi
def test_read_chi(testPu9):
    S = testPu9.read_section(9437, 5, 18)
    from .MF5 import write
    text = write(S)
    assert testPu9.TEXT.loc[9437,5,18] == text

@pytest.fixture(scope="module")
@pytest.mark.formats
@pytest.mark.endf6
@pytest.mark.chi
def testPu9chi(testPu9):
    return testPu9.get_edistr()

@pytest.mark.formats
@pytest.mark.endf6
@pytest.mark.chi
def test_extract_chi(testPu9chi):
    testPu9chi.add_points([1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1e0, 1e1, 1e2, 1e3, 1e4, 1e5])

@pytest.mark.formats
@pytest.mark.endf6
@pytest.mark.chi
def test_normalize_chi(testPu9chi):
    chi = testPu9chi.normalize()
    for i,v in chi.iterrows():
        dx = v.index.values[1:] - v.index.values[:-1]
        y = (v.values[1:]+v.values[:-1])/2
        assert np.isclose(y.dot(dx), 1, rtol=1e-10)

@pytest.mark.formats
@pytest.mark.endf6
@pytest.mark.chi
def test_update_chi(testPu9, testPu9chi):
    testPu9chi.loc[(9437,18,1,1e-5)] = 1
    new = testPu9.update_edistr(testPu9chi)
    new_sec = new.read_section(9437,5,18)
    assert (new_sec["PDISTR"][1]["EIN"][1e-5]["EDISTR"] == 1).all()

@pytest.mark.formats
@pytest.mark.endf6
@pytest.mark.chi
@pytest.mark.cov
def test_read_chi_cov(testPu9):
    S = testPu9.read_section(9437, 35, 18)
    assert len(S["SUB"]) == S["NK"] == 8
    for k,v in S["SUB"].items():
        assert v["LS"] == 1
        assert v["LB"] == 7

@pytest.mark.formats
@pytest.mark.endf6
@pytest.mark.chi
@pytest.mark.cov
def test_extract_chi_cov(testPu9):
    C = testPu9.get_edistr_cov()
    assert C.index.names == ['MAT', 'MT', 'ELO', 'EHI', 'EOUT']
    assert C.columns.names == ['MAT', 'MT', 'ELO', 'EHI', 'EOUT']
    assert (C.values == C.values.T).all()

@pytest.mark.formats
@pytest.mark.endf6
@pytest.mark.write
def test_write_to_string(testH1):
    string = testH1.write_string()
    newtape = Endf6.from_text(string)
    assert testH1.equals(newtape)