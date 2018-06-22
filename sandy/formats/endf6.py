# -*- coding: utf-8 -*-
"""
Created on Mon Jan 16 18:03:13 2017

@author: lfiorito
"""
import sys, pdb, os, pytest
import numpy as np
import pandas as pd
from copy import deepcopy
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
#def process_endf_section(text, keep_mf=None, keep_mt=None):
#    mf = int(text[70:72])
#    mt = int(text[72:75])
#    if mf == 1 and mt == 451: # read always
#        return read_mf1_mt451(text)
#    if keep_mf:
#        if mf not in keep_mf:
#            return None
#    if keep_mt:
#        if mt not in keep_mt:
#            return None
#    if mf == 1 and mt in (452, 455, 456):
#        return read_mf1_nubar(text)
#    elif mf == 3:
#        return read_mf3_mt(text)
#    elif mf == 5:
#        return read_mf5_mt(text)
#    elif mf == 31 or mf == 33:
#        return read_mf33_mt(text)
#    elif mf == 35:
#        return read_mf35_mt(text)
#    else:
#        return None


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
        tape = tape.query("MAT>0 & MF>0 & MT>0")
        tape = tape.groupby(["MAT","MF","MT"]).sum()
        return cls(tape)

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
        elif mf == 33 or mf == 31:
            from .MF33 import read
        else:
            sys.exit("ERROR: SANDY cannot parse section MAT{}/MF{}/MT{}".format(mat,mf,mt))
        return read(self.loc[mat,mf,mt].TEXT)

#    def process(self, keep_mf=None, keep_mt=None):
#        """
#        Parse TEXT column.
#        """
#        tape = self.copy()
#        tape['DATA'] = tape['TEXT'].apply(process_endf_section, keep_mf=keep_mf, keep_mt=keep_mt)
#        return Endf6(tape)

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
                        string += text.encode('ascii', 'replace').decode('ascii') + "\n"
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

    def get_edistr(self, listmat=None, listmt=None):
        """
        Extract chi cov for all MAT,MT,SUB found in tape.
        Return a df with MAT,MT,SUB as index and COV as value
        Each COV is a df with Ein on rows and Eout on columns.
        """
        query = "MF==5"
        if listmat is not None:
            query_mats = " | ".join(["MAT=={}".format(x) for x in listmat])
            query += " & ({})".format(query_mats)
        if listmt is not None:
            query_mts = " | ".join(["MT=={}".format(x) for x in listmt])
            query += " & ({})".format(query_mts)
        tape = self.query(query)
        DictDf = { "MAT" : [], "MT" : [], "K" : [], "EDISTR" : [] }
        for ix,text in tape.TEXT.iteritems():
            X = self.read_section(*ix)
            for k,pdistr in X["PDISTR"].items():
                if pdistr["LF"] != 1:
                    print("WARNING: non-tabulated distribution for MAT{}/MF{}/MT{}, subsec {}".format(*ix,k))
                    continue
                if list(filter(lambda x:x["INT"] != [2], pdistr["EIN"].values())):
                    print("WARNING: found non-linlin interpolation, skip energy distr. for MAT{}/MF{}/MT{}, subsec {}".format(*ix,k))
                    continue
                D =  { ein : pd.Series(v["EDISTR"], index=v["EOUT"]) for ein,v in sorted(pdistr["EIN"].items()) }
                # merge chi_E(E') distributions on df[E,E'] with unique E' grid
                # Interpolate row-by-row at missing datapoints
                # When interpolation is not possible (edges) fill with 0
                edistr = pd.DataFrame.from_dict(D, orient='index').interpolate(method="slinear", axis=1).fillna(0)
    #            # include default points in E grid and interpolate
    #            ein_new = list(union_grid(chi.index, np.logspace(-5, 7, 13)))
    #            chi = pandas_interpolate(chi, ein_new, method='slinear', axis='rows')
                edistr.index.name = "EIN"
                edistr.columns.name = "EOUT"
                DictDf["MAT"].append(X["MAT"])
                DictDf["MT"].append(X["MT"])
                DictDf["K"].append(k)
                DictDf["EDISTR"].append(edistr)
        DfChi = pd.DataFrame.from_dict(DictDf).set_index(["MAT", "MT", "K"])
        return DfChi

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
                        warn("skipped NI-type covariance with flag LB={} for MAT{}/MF{}/MT{}".fomat(nisec["LB"], *ix), category=Warning)
                        continue
                    cov = pd.DataFrame(cov, index=e1, columns=e2)
                    covs.append(cov)
                if len(covs) == 0:
                    continue
                if len(covs) > 1:
                    pdb.set_trace()
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
        self.THNMAX = - self.EHRES if self.EHRES != 0 else 2.0E6



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








def extract_chi(tape):
    """
    Extract chi cov for all MAT,MT,SUB found in tape.
    Return a df with MAT,MT,SUB as index and COV as value
    Each COV is a df with Ein on rows and Eout on columns.
    """
#    from sandy.functions import union_grid
    DictDf = { "MAT" : [], "MT" : [], "K" : [], "CHI" : [] }
    for chunk in tape.query('MF==5').DATA:
        for k,sub in enumerate(chunk["SUB"]):
            if sub["LF"] != 1:
                continue
            if list(filter(lambda x:x["INT"] != [2], sub["EIN"].values())):
                print("WARNING: found non-linlin interpolation, skip energy distr. for MAT {}, MT {}, subsec {}".format(chunk["MAT"],chunk["MT"],k))
                continue
            chi_dict =  { ein : ssub["PDF"] for ein,ssub in sorted(sub["EIN"].items()) }
            # merge chi_E(E') distributions on df[E,E'] with unique E' grid
            # Interpolate row-by-row at missing datapoints
            # When interpolation is not possible (edges) fill with 0
            chi = pd.DataFrame.from_dict(chi_dict, orient='index')
            chi = chi.interpolate(method="slinear", axis=1).fillna(0)
#            # include default points in E grid and interpolate
#            ein_new = list(union_grid(chi.index, np.logspace(-5, 7, 13)))
#            chi = pandas_interpolate(chi, ein_new, method='slinear', axis='rows')
            chi.index.name = "EIN"
            chi.columns.name = "EOUT"
            DictDf["MAT"].append(chunk["MAT"])
            DictDf["MT"].append(chunk["MT"])
            DictDf["K"].append(k)
            DictDf["CHI"].append(chi)
    DfChi = pd.DataFrame.from_dict(DictDf)
    DfChi.set_index(["MAT", "MT", "K"], inplace=True)
    return DfChi

def extract_mu(tape):
#    NLMAX = 0
    keys = []
    ListPc = []
    for chunk in tape.query('MF==4').DATA:
        if "LPC" in chunk:
            # check interpolation
            for e,sub in sorted(chunk["LPC"]["E"].items()):
                keys.append((chunk["MAT"], chunk["MT"], e))
#                NPL = len(sub["P"])
#                P[i,:NPL] = sub["P"]
#                NLMAX = max(NLMAX, NPL)
                ListPc.append(sub["P"])
    index = pd.MultiIndex.from_tuples(keys, names=("MAT", "MT", "E"))
    DfPc = pd.DataFrame(ListPc, index=index).fillna(0)
    DfPc.columns = range(1, DfPc.shape[1]+1)
    DfPc.columns.name = "LegCoeff"
    # print warning if exceeded NLMAX
    return DfPc

def extract_cov35(tape):
    from sandy.sampling.cov import triu_matrix, corr2cov
    # Get covariances (df) from each MAT, MT and Erange (these are the keys) in a dictionary.
    DictCov = {}
    for chunk in tape.query('MF==35').DATA:
        for sub in chunk["SUB"]:
            # Ek grid is one unit longer than covariance.
            Ek = np.array(sub["Ek"])
            Fkk = np.array(sub["Fkk"])
            NE = sub["NE"]
            cov = triu_matrix(Fkk, NE-1)
            # Normalize covariance matrix dividing by the energy bin.
            dE = 1./(Ek[1:]-Ek[:-1])
            cov = corr2cov(cov, dE)
            # Add zero row and column at the end of the matrix
            cov = np.insert(cov, cov.shape[0], [0]*cov.shape[1], axis=0)
            cov = np.insert(cov, cov.shape[1], [0]*cov.shape[0], axis=1)
            cov = pd.DataFrame(cov, index=Ek, columns=Ek)
            cov.index.name = "Ek"; cov.columns.name = "El"
            DictCov.update({ (chunk["MAT"], chunk["MT"], sub["Elo"], sub["Ehi"]) :
                cov })
    if not DictCov:
        # No MF35 found. Return empty dataframe
        return pd.DataFrame()
    # Collect covs in a list to pass to block_diag (endf6 allows only covs in diag)
    # Recreate indexes with lists of mat, mt, elo, ehi and e.
    from scipy.linalg import block_diag
    covs = []; mat = []; mt = []; einlo = []; einhi = []; eoutlo= []; eouthi = []
    for k,c in sorted(DictCov.items()):
        covs.append(c)
        mat.extend([k[0]]*len(c.index))
        mt.extend([k[1]]*len(c.index))
        einlo.extend([k[2]]*len(c.index))
        einhi.extend([k[3]]*len(c.index))
        eoutlo.extend(c.index)
        eouthi.extend(c.index[1:].insert(-1, c.index[-1]))
    DfCov = block_diag(*covs)
    DfCov = pd.DataFrame(DfCov)
    DfCov['MAT'] = pd.Series(mat)
    DfCov['MT'] = pd.Series(mt)
    DfCov['EINlo'] = pd.Series(einlo)
    DfCov['EINhi'] = pd.Series(einhi)
    DfCov['EOUTlo'] = pd.Series(eoutlo)
    DfCov['EOUThi'] = pd.Series(eouthi)
    DfCov.set_index(['MAT', 'MT', 'EINlo', 'EINhi', 'EOUTlo', 'EOUThi'], inplace=True)
    DfCov.columns = DfCov.index
    return DfCov



def update_chi(tapein, DfChi):
    tape = pd.DataFrame(index=tapein.index.copy(), columns=tapein.columns.copy())
    for k,row in tapein.iterrows():
        tape.loc[k].DATA = deepcopy(row.DATA)
        tape.loc[k].TEXT = deepcopy(row.TEXT)
    for (mat, mt, k),CHI in DfChi.CHI.items():
        if (mat, 5, mt) not in tape.index:
            continue
        if len(tape.DATA.loc[mat,5,mt]["SUB"]) - 1 < k:
            continue
        # Assume that all given chi are linearly interpolated on Ein and Eout
        D = { ein : {"PDF" : PdfSeries, "NBT" : [CHI.shape[1]], "INT" : [2]} for ein,PdfSeries in CHI.iterrows()}
        tape.DATA.loc[mat,5,mt]["SUB"][k].update({'EIN' : D,
                                                  'INT_EIN' : [2],
                                                  'NBT_EIN' : [CHI.shape[0]]})
    return tape



def write_decay_data_csv(tape, filename):
    df = tape.query("MF==8 & MT==457")
    df['Z'] = df.DATA.apply(lambda x : int(x['ZA']//1000))
    df['A'] = df.DATA.apply(lambda x : int(x['ZA'] - x['ZA']//1000*1000))
    df['M'] = df.DATA.apply(lambda x : 'g' if x['LISO'] == 0 else 'm' if x['LISO'] == 1 else 'n')
    df['HL'] = df.DATA.apply(lambda x : x['HL'])
    df['DHL'] = df.DATA.apply(lambda x : x['DHL'])
    df['ELP'] = df.DATA.apply(lambda x : x['E'][0])
    df['DELP'] = df.DATA.apply(lambda x : x['DE'][0])
    df['EEM'] = df.DATA.apply(lambda x : x['E'][1])
    df['DEEM'] = df.DATA.apply(lambda x : x['DE'][2])
    df['EHP'] = df.DATA.apply(lambda x : x['E'][2])
    df['DEHP'] = df.DATA.apply(lambda x : x['DE'][2])
    df[['Z','A','M','HL','DHL','ELP','DELP','EEM','DEEM','EHP','DEHP']].to_csv(filename, index=False)


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
@pytest.mark.info
def test_read_info(testPu9):
    S = testPu9.read_section(9437, 1, 451)
    from .MF1 import write
    text = write(S)
    assert testPu9.TEXT.loc[9437,1,451] == text

@pytest.mark.formats
@pytest.mark.endf6
@pytest.mark.nubar
def test_read_nubar452(testPu9):
    S = testPu9.read_section(9437, 1, 452)
    from .MF1 import write
    text = write(S)
    assert testPu9.TEXT.loc[9437,1,452] == text

@pytest.mark.formats
@pytest.mark.endf6
@pytest.mark.nubar
def test_read_nubar455(testPu9):
    S = testPu9.read_section(9437, 1, 455)
    from .MF1 import write
    text = write(S)
    assert testPu9.TEXT.loc[9437,1,455] == text

@pytest.mark.formats
@pytest.mark.endf6
@pytest.mark.nubar
def test_read_nubar456(testPu9):
    S = testPu9.read_section(9437, 1, 456)
    from .MF1 import write
    text = write(S)
    assert testPu9.TEXT.loc[9437,1,456] == text

@pytest.mark.formats
@pytest.mark.endf6
@pytest.mark.chi
def test_read_chi(testPu9):
    S = testPu9.read_section(9437, 5, 18)
    from .MF5 import write
    text = write(S)
    assert testPu9.TEXT.loc[9437,5,18] == text

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
@pytest.mark.chi
def test_extract_chi(testPu9):
    testPu9.get_edistr()

@pytest.mark.formats
@pytest.mark.endf6
@pytest.mark.xs
@pytest.mark.cov
def test_extract_xs_cov(testPu9):
    testPu9.get_xs_cov()

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
@pytest.mark.write
def test_write_to_string(testH1):
    string = testH1.write_string()
    newtape = Endf6.from_text(string)
    assert testH1.equals(newtape)