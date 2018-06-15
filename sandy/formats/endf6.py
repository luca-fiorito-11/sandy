# -*- coding: utf-8 -*-
"""
Created on Mon Jan 16 18:03:13 2017

@author: lfiorito
"""
import sys, time, pdb, os, re, pytest
import numpy as np
#from . import read_cont, read_tab1, read_tab2, read_list, read_text, read_float
#from . import write_cont, write_tab1, write_list, write_tab2
#from .records2 import read_control, read_cont
import pandas as pd
from copy import deepcopy
from warnings import warn
from sandy.tests import TimeDecorator


def split_endf(text):
    """
    Read ENDF-6 formatted file and split it into columns based on field widths:
        C1 C2 L1 L2 N1 N2 MAT MF MT
        11 11 11 11 11 11  4   2  3.
    Store list in dataframe.
    """
    from io import StringIO
    def read_float(x):
        try:
            return float(x[0] + x[1:].replace('+', 'E+').replace('-', 'E-'))
        except:
            return x
    widths = [11,11,11,11,11,11,4,2,3]
    columns = ["C1", "C2", "L1", "L2", "N1", "N2","MAT", "MF", "MT"]
    converters = dict(zip(columns[:6],[read_float]*6))
    frame =  pd.read_fwf(StringIO(text), widths=widths, names=columns, converters=converters)
    return frame.query("MAT>0 & MF>0 & MT>0")


def process_endf_section(text, keep_mf=None, keep_mt=None):
    mf = int(text[70:72])
    mt = int(text[72:75])
    if mf == 1 and mt == 451: # read always
        return read_mf1_mt451(text)
    if keep_mf:
        if mf not in keep_mf:
            return None
    if keep_mt:
        if mt not in keep_mt:
            return None
    if mf == 1 and mt in (452, 455, 456):
        return read_mf1_nubar(text)
    elif mf == 3:
        return read_mf3_mt(text)
    elif mf == 5:
        return read_mf5_mt(text)
    elif mf == 31 or mf == 33:
        return read_mf33_mt(text)
    elif mf == 35:
        return read_mf35_mt(text)
    else:
        return None


class Endf6(pd.DataFrame):

    @classmethod
    def from_file(cls, file):
        """
        Read ENDF-6 formatted file and call from_text method.
        """
        with open(file) as f:
            text = f.read()
        return cls.from_text(text)

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
        tape = tape.query("MAT>0 & MF>0 & MT>0").groupby(["MAT","MF","MT"]).sum()
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
        else:
            sys.exit("ERROR: SANDY cannot parse section MAT{}/MF{}/MT{}".format(mat,mf,mt))
        return read(self.loc[mat,mf,mt].TEXT)

    def process(self, keep_mf=None, keep_mt=None):
        """
        Parse TEXT column.
        """
        tape = self.copy()
        tape['DATA'] = tape['TEXT'].apply(process_endf_section, keep_mf=keep_mf, keep_mt=keep_mt)
        return Endf6(tape)

#    def to_string(self, title=" "*66):
#        """
#        First update MF1/MT451 dictionary.
#        Then, Write TEXT column to string.
#        """
#        tape = self.update_dict().write_mf1_mt451()
#        string = "{:<66}{:4}{:2}{:3}{:5}\n".format(title, 1, 0, 0, 0)
#        for mat,dfmat in tape.groupby('MAT', sort=True):
#            for mf,dfmf in dfmat.groupby('MF', sort=True):
#                for mt,dfmt in dfmf.groupby('MT', sort=True):
#                    for text in dfmt.TEXT:
#                        string += text.encode('ascii', 'replace').decode('ascii') + "\n"
#                    string += "{:<66}{:4}{:2}{:3}{:5}\n".format(*write_cont(*[0]*6), mat, mf, 0, 99999)
#                string += "{:<66}{:4}{:2}{:3}{:5}\n".format(*write_cont(*[0]*6), mat, 0, 0, 0)
#            string += "{:<66}{:4}{:2}{:3}{:5}\n".format(*write_cont(*[0]*6), 0, 0, 0, 0)
#        string += "{:<66}{:4}{:2}{:3}{:5}".format(*write_cont(*[0]*6), -1, 0, 0, 0)
#        return string

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
#    def get_xs(self):
#        xsList = []
#        for mat in self.index.get_level_values('MAT').unique():
#            query = "(MF==1 & (MT==452 | MT==455 | MT==456))"
#            if self.DATA.loc[mat,1,451]['LRP'] == 2:
#                query += ' | MF==3'
#            for chunk in self.loc[mat].query(query).DATA:
#                if not chunk: continue
#                key = "NUBAR" if chunk["MF"] == 1 else "XS"
#                xs = chunk[key]
#                duplicates = [x for x, count in Counter(xs.index).items() if count > 1]
#                if duplicates:
#                    sys.exit('ERROR: duplicate energy points found for MAT{}/MF{}/MT{}\n'.format(chunk["MAT"],chunk["MF"],chunk["MT"])+
#                             '\n'.join(map(str,duplicates)))
#                if chunk['INT'] != [2]:
#                    sys.exit('ERROR: MAT{}/MF{}/MT{} interpolation scheme is not lin-lin'.format(chunk["MAT"],chunk["MF"],chunk["MT"]))
#                xsList.append(xs.to_frame())
#        xs = reduce(lambda left,right : pd.merge(left, right, left_index=True, right_index=True, how='outer'), xsList).sort_index().interpolate(method='slinear', axis=0).fillna(0)
#        xs.columns.names = ["MAT", "MT"]
#        return Xs(xs)
#        from collections import Counter
#        from functools import reduce
#        xsList = []
#        for mat in self.index.get_level_values('MAT').unique():
#            query = "(MF==1 & (MT==452 | MT==455 | MT==456))"
#            if self.DATA.loc[mat,1,451]['LRP'] == 2:
#                query += ' | MF==3'
#            for chunk in self.loc[mat].query(query).DATA:
#                if not chunk: continue
#                key = "NUBAR" if chunk["MF"] == 1 else "XS"
#                xs = chunk[key]
#                duplicates = [x for x, count in Counter(xs.index).items() if count > 1]
#                if duplicates:
#                    sys.exit('ERROR: duplicate energy points found for MAT{}/MF{}/MT{}\n'.format(chunk["MAT"],chunk["MF"],chunk["MT"])+
#                             '\n'.join(map(str,duplicates)))
#                if chunk['INT'] != [2]:
#                    sys.exit('ERROR: MAT{}/MF{}/MT{} interpolation scheme is not lin-lin'.format(chunk["MAT"],chunk["MF"],chunk["MT"]))
#                xsList.append(xs.to_frame())
#        xs = reduce(lambda left,right : pd.merge(left, right, left_index=True, right_index=True, how='outer'), xsList).sort_index().interpolate(method='slinear', axis=0).fillna(0)
#        xs.columns.names = ["MAT", "MT"]
#        return Xs(xs)

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

    def get_cov(self):
        from sandy.sampling.cov import triu_matrix
        from functools import reduce
        List = []; eg = set()
        for chunk in self.query('MF==33 | MF==31').DATA:
            mat = chunk['MAT']; mt = chunk['MT']
            for sub in chunk["SUB"].values():
                mat1 = sub['MAT1'] if sub['MAT1'] != 0 else mat; mt1 = sub['MT1']
                covs = []
                for nisec in sub["NI"]:
                    if nisec["LB"] == 5:
                        Fkk = np.array(nisec["Fkk"])
                        if nisec["LS"] == 0: # to be tested
                            cov = Fkk.reshape(nisec["NE"]-1, nisec["NE"]-1)
                        else:
                            cov = triu_matrix(Fkk, nisec["NE"]-1)
                        # add zero row and column at the end of the matrix
                        cov = np.insert(cov, cov.shape[0], [0]*cov.shape[1], axis=0)
                        cov = np.insert(cov, cov.shape[1], [0]*cov.shape[0], axis=1)
                        e1 = e2 = nisec["Ek"]
                    elif nisec["LB"] == 1:
                        cov = np.diag(nisec["Fk"])
                        e1 = e2 = nisec["Ek"]
                    elif nisec["LB"] == 2:
                        f = np.array(nisec["Fk"])
                        cov = f*f.reshape(-1,1)
                        e1 = e2 = nisec["Ek"]
                    elif nisec["LB"] == 6:
                        cov = np.array(nisec["Fkl"]).reshape(nisec["NER"]-1, nisec["NEC"]-1)
                        # add zero row and column at the end of the matrix
                        cov = np.insert(cov, cov.shape[0], [0]*cov.shape[1], axis=0)
                        cov = np.insert(cov, cov.shape[1], [0]*cov.shape[0], axis=1)
                        e1 = nisec["Ek"]
                        e2 = nisec["El"]
                    else:
                        warn("skipped NI-type covariance with flag LB={} (MAT{}/MF{}/MT{})".fomat(nisec["LB"], chunk['MAT'], chunk['MF'], chunk['MT']), category=Warning)
                        continue
                    cov = pd.DataFrame(cov, index=e1, columns=e2)
                    covs.append(cov)
                if len(covs) == 0:
                    continue
                if len(covs) > 1:
                    import pdb
                    pdb.set_trace()
                cov = reduce(lambda x, y: x.add(y, fill_value=0).fillna(0), covs).fillna(0)
                eg |= set(cov.index.values)
                List.append([mat, mt, mat1, mt1, cov])
        if not List:
            warn("no MF[31,33] covariances found", category=Warning)
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

    def update_dict(self):
        """
        Update RECORDS item (in DATA column) for MF1/MT451 of each MAT based on the content of the TEXT column.
        """
        tape = pd.DataFrame(index=self.index.copy(), columns=self.columns.copy())
        for k,row in self.iterrows():
            tape.loc[k].DATA = deepcopy(row.DATA)
            tape.loc[k].TEXT = deepcopy(row.TEXT)
        for mat in sorted(tape.index.get_level_values('MAT').unique()):
            chunk = tape.DATA.loc[mat,1,451]
            records = pd.DataFrame(chunk["RECORDS"],
                                   columns=["MF","MT","NC","MOD"]).set_index(["MF","MT"])
            new_records = []
            for (mf,mt),text in sorted(tape.loc[mat].query('MT!=451'.format(mat)).TEXT.items()):
                nc = len(text.splitlines())
                # when copying PENDF sections (MF2/MT152) mod is not present in the dictionary
                try:
                    mod = records.MOD.loc[mf,mt]
                except:
                    mod = 0
                new_records.append((mf,mt,nc,mod))
            nc = 4 + len(chunk["TEXT"]) + len(new_records) + 1
            mod = records.MOD.loc[1,451]
            new_records = [(1,451,nc,mod)] + new_records
            chunk["RECORDS"] = new_records
        return Endf6(tape)

    def write_mf1_mt451(self):
        tape = pd.DataFrame(index=self.index.copy(), columns=self.columns.copy())
        for k,row in self.iterrows():
            tape.loc[k].DATA = deepcopy(row.DATA)
            tape.loc[k].TEXT = deepcopy(row.TEXT)
        for (mat,mf,mt),df in tape.query('MF==1 & MT==451').iterrows():
            TEXT = write_cont(df.DATA["ZA"], df.DATA["AWR"], df.DATA["LRP"], df.DATA["LFI"], df.DATA["NLIB"], df.DATA["NMOD"])
            TEXT += write_cont(df.DATA["ELIS"], df.DATA["STA"], df.DATA["LIS"], df.DATA["LISO"], 0, df.DATA["NFOR"])
            TEXT += write_cont(df.DATA["AWI"], df.DATA["EMAX"], df.DATA["LREL"], 0, df.DATA["NSUB"], df.DATA["NVER"])
            TEXT += write_cont(df.DATA["TEMP"], 0, df.DATA["LDRV"], 0, len(df.DATA["TEXT"]), len(df.DATA["RECORDS"]))
            TEXT += df.DATA["TEXT"]
            TEXT += [ " "*22 + "{:>11}{:>11}{:>11}{:>11}".format(mfnxc,mtnxc,ncnxc,modnxc) for mfnxc,mtnxc,ncnxc,modnxc in df.DATA["RECORDS"]]
            TextOut = []; iline = 1
            for line in TEXT:
                if iline > 99999:
                    iline = 1
                TextOut.append("{:<66}{:4}{:2}{:3}{:5}".format(line, mat, mf, mt, iline))
                iline += 1
            tape.at[(mat,mf,mt),'TEXT'] = "\n".join(TextOut)
        return Endf6(tape)

    def write_mf1_nubar(self):
        tape = pd.DataFrame(index=self.index.copy(), columns=self.columns.copy())
        for k,row in self.iterrows():
            tape.loc[k].DATA = deepcopy(row.DATA)
            tape.loc[k].TEXT = deepcopy(row.TEXT)
        for (mat,mf,mt),df in tape.query('MF==1 & (MT==452 | MT==455 | MT==456)').iterrows():
            TEXT = write_cont(df.DATA["ZA"], df.DATA["AWR"], df.DATA["LDG"], df.DATA["LNU"], 0, 0)
            if df.DATA["MT"] == 455:
                if df.DATA["LDG"] == 0:
                    TEXT += write_list(0, 0, 0, 0, 0, df.DATA["LAMBDAS"])
                elif df.DATA["LDG"] == 1:
                    # Not found in JEFF33 and ENDFB8, hence not implemented
                    pass
            if df.DATA["LNU"] == 1:
                TEXT += write_list(0, 0, 0, 0, 0, df.DATA["C"])
            else:
                TEXT += write_tab1(0, 0, 0, 0, df.DATA["NBT"], df.DATA["INT"],
                               df.DATA["NUBAR"].index, df.DATA["NUBAR"])
            TextOut = []; iline = 1
            for line in TEXT:
                if iline > 99999:
                    iline = 1
                TextOut.append("{:<66}{:4}{:2}{:3}{:5}".format(line, mat, mf, mt, iline))
                iline += 1
            tape.at[(mat,mf,mt),'TEXT'] = "\n".join(TextOut)
        return Endf6(tape)

    def write_mf3_mt(self):
        tape = pd.DataFrame(index=self.index.copy(), columns=self.columns.copy())
        for k,row in self.iterrows():
            tape.loc[k].DATA = deepcopy(row.DATA)
            tape.loc[k].TEXT = deepcopy(row.TEXT)
        for (mat,mf,mt),df in tape.query('MF==3').iterrows():
            if df.DATA is None:
                continue
            TEXT = write_cont(df.DATA["ZA"], df.DATA["AWR"], 0, 0, 0, 0)
            TEXT += write_tab1(df.DATA["QM"], df.DATA["QI"], 0, df.DATA["LR"],
                               df.DATA["NBT"], df.DATA["INT"],
                               df.DATA["XS"].index, df.DATA["XS"])
            TextOut = []; iline = 1
            for line in TEXT:
                if iline > 99999:
                    iline = 1
                TextOut.append("{:<66}{:4}{:2}{:3}{:5}".format(line, mat, mf, mt, iline))
                iline += 1
            tape.at[(mat,mf,mt),'TEXT'] = "\n".join(TextOut)
        return Endf6(tape)

def read_mf2_mt151(text):
    str_list = split_endf(text)
    i = 0
    out = {"MAT" : str_list["MAT"].iloc[0],
           "MF" : str_list["MF"].iloc[0],
           "MT" : str_list["MT"].iloc[0]}
    C, i = read_cont(str_list, i)
    out.update({ "ZA" : C.C1, "AWR" : C.C2, "NIS" : C.N1, "ZAI" : {} })
    for i_iso in range(out["NIS"]): # LOOP ISOTOPES
        C, i = read_cont(str_list, i)
        zai = { "ZAI" : C.C1, "ABN" : C.C2, "LFW" : C.L2, "NER" : C.N1, "ERANGE" : {} }
        for i_erange in range(zai["NER"]): # LOOP ENERGY RANGES
            C, i = read_cont(str_list, i)
            sub = {"EL" : C.C1, "EH" : C.C2, "LRU" : C.L1, "LRF" : C.L2, "NRO" : C.N1, "NAPS" : C.N2}
            if sub["NRO"] != 0: # Tabulated scattering radius
                T, i = read_tab1(str_list, i)
                sub.update({"NBT" : T.NBT, "INT" : T.INT})
                sub["AP"] = pd.Series(T.y, index = T.x).rename_axis("E")
            if sub["LRU"] == 0: # ONLY SCATTERING RADIUS
                C, i = read_cont(str_list, i)
                sub.update({"SPI" : C.C1, "SR" : C.C2, "NLS" : C.N1})
            if sub["LRU"] == 1: # RESOLVED RESONANCES
                if sub["LRF"] in (1,2): # BREIT-WIGNER
                    C, i = read_cont(str_list, i)
                    sub.update({"SPI" : C.C1, "SR" : C.C2, "NLS" : C.N1})
                    L, i = read_list(str_list, i)
                elif sub["LRF"] == 3: # REICH-MOORE
                    C, i = read_cont(str_list, i)
                    sub.update({"SPI" : C.C1, "SR" : C.C2, "LAD" : C.L1, "NLS" : C.N1, "NLSC" : C.N2})
                    L, i = read_list(str_list, i)
                elif sub["LRF"] == 4: # ADLER-ADLER
                    sys.exit("ERROR: SANDY cannot read resonance parameters in Adler-Adler formalism")
                elif sub["LRF"] == 5: # GENERAL R-MATRIX
                    sys.exit("ERROR: General R-matrix formalism no longer available in ENDF-6")
                elif sub["LRF"] == 6: # HYBRID R-FUNCTION
                    sys.exit("ERROR: Hybrid R-function formalism no longer available in ENDF-6")
                elif sub["LRF"] == 7: # LIMITED R-MATRIX
                    C, i = read_cont(str_list, i)
                    sub.update({"IFG" : C.L1, "KRM" : C.L2, "NJS" : C.N1, "KRL" : C.N2})
                    for j in range(sub["NJS"]):
                        L1, i = read_list(str_list, i)
                        L2, i = read_list(str_list, i)
                        L3, i = read_list(str_list, i)
            elif sub["LRU"] == 2: # UNRESOLVED RESONANCES
                if sub["LRF"] == 1 and out["LFW"] == 0: # CASE A
                    C, i = read_cont(str_list, i)
                    sub.update({"SPI" : C.C1, "SR" : C.C2, "LSSF" : C.L1, "NLS" : C.N1})
                    for k in range(sub["NLS"]):
                        L, i = read_list(str_list, i)
                elif sub["LRF"] == 1 and out["LFW"] == 1: # CASE B
                    C, i = read_cont(str_list, i)
                    sub.update({"SPI" : C.C1, "SR" : C.C2, "LSSF" : C.L1, "NE" : C.N1, "NLS" : C.N2})
                    L, i = read_list(str_list, i)
                    for k in range(sub["NLS"]):
                        C, i = read_cont(str_list, i)
                        for l in range(C.N1):
                            L, i = read_list(str_list, i)
                elif sub["LRF"] == 2: # CASE C
                    C, i = read_cont(str_list, i)
                    sub.update({"SPI" : C.C1, "SR" : C.C2, "LSSF" : C.L1, "NLS" : C.N1})
                    for k in range(sub["NLS"]):
                        C, i = read_cont(str_list, i)
                        for l in range(C.N1):
                            L, i = read_list(str_list, i)
            zai["ERANGE"].update({ (sub["EL"],sub["EH"]) : sub })
        out["ZAI"].update({ zai["ZAI"] : zai })
    return out

def read_mf4_mt(text):
    str_list = split_endf(text)
    i = 0
    out = {"MAT" : str_list["MAT"].iloc[0],
           "MF" : str_list["MF"].iloc[0],
           "MT" : str_list["MT"].iloc[0]}
    C, i = read_cont(str_list, i)
    out.update({"ZA" : C.C1, "AWR" : C.C2, "LTT" : C.L2})
    if out["LTT"] in (1,3):
        C, i = read_cont(str_list, i)
        lpc = {"LI" : C.L1, "LCT" : C.L2}
        T2, i = read_tab2(str_list, i)
        lpc.update({"NE" : T2.NZ, "NBT" : T2.NBT, "INT" : T2.INT, "E" : {} })
        for j in range(lpc["NE"]):
            L, i = read_list(str_list, i)
            lpc["E"].update({ L.C2 : {"P" : L.B, "T" : L.C1, "LT" : L.L1}})
        out.update({"LPC" : lpc})
    if out["LTT"] in (2,3):
        C, i = read_cont(str_list, i)
        sub = {"LI" : C.L1, "LCT" : C.L2}
        T2, i = read_tab2(str_list, i)
        sub.update({"NE" : T2.NZ, "NBT" : T2.NBT, "INT" : T2.INT, "E" : {}})
        for i in range(sub["NE"]):
            T1, i = read_tab1(str_list, i)
            distr = pd.Series(T1.y, index = T1.x, name=T1.C2).rename_axis("mu")
            sub["E"].update({ T1.C2 : {"T" : T1.C1, "LT" : T1.L1, "PDF" : distr, "NBT" : T1.NBT, "INT" : T1.INT}})
        out.update({"TAB" : sub})
    if out["LTT"] == 0:
        C, i = read_cont(str_list, i)
        out.update({"ISO" : {"LI" : C.L1, "LCT" : C.L2}})
    return out


def write_mf5_mt(tapein):
    tape = pd.DataFrame(index=tapein.index.copy(), columns=tapein.columns.copy())
    for k,row in tapein.iterrows():
        tape.loc[k].DATA = deepcopy(row.DATA)
        tape.loc[k].TEXT = deepcopy(row.TEXT)
    for (mat,mf,mt),df in tape.loc[(slice(None),5),:].iterrows():
        TEXT = write_cont(df.DATA["ZA"], df.DATA["AWR"], 0, 0, len(df.DATA["SUB"]), 0)
        for sub in df.DATA["SUB"]:
            U = sub['U'] if 'U' in sub else 0
            TEXT += write_tab1(U, 0, 0, sub["LF"],
                               sub["NBT_P"], sub["INT_P"],
                               sub["P"].index, sub["P"])
            if sub["LF"] == 1:
                TEXT += write_tab2(0, 0, 0, 0,
                                   len(sub['EIN']),
                                   sub["NBT_EIN"], sub["INT_EIN"])
                for ein, distr in sorted(sub['EIN'].items()):
                    TEXT += write_tab1(0, ein, 0, 0,
                                       distr["NBT"], distr["INT"],
                                       distr["PDF"].index, distr["PDF"])
            elif sub["LF"] == 5:
                t = sub['Theta']
                TEXT += write_tab1(t.C1, t.C2, t.L1, t.L2,
                                   t.NBT, t.INT,
                                   t.x, t.y)
                t = sub['Tg']
                TEXT += write_tab1(t.C1, t.C2, t.L1, t.L2,
                                   t.NBT, t.INT,
                                   t.x, t.y)
            elif sub["LF"] in (7,9):
                t = sub['Theta']
                TEXT += write_tab1(t.C1, t.C2, t.L1, t.L2,
                                   t.NBT, t.INT,
                                   t.x, t.y)
            elif sub["LF"] == 11:
                t = sub['Ta']
                TEXT += write_tab1(t.C1, t.C2, t.L1, t.L2,
                                   t.NBT, t.INT,
                                   t.x, t.y)
                t = sub['Tb']
                TEXT += write_tab1(t.C1, t.C2, t.L1, t.L2,
                                   t.NBT, t.INT,
                                   t.x, t.y)
            elif sub["LF"] == 12:
                t = sub['TTm']
                TEXT += write_tab1(t.C1, t.C2, t.L1, t.L2,
                                   t.NBT, t.INT,
                                   t.x, t.y)
        TextOut = []; iline = 1
        for line in TEXT:
            if iline > 99999:
                iline = 1
            TextOut.append("{:<66}{:4}{:2}{:3}{:5}".format(line, mat, mf, mt, iline))
            iline += 1
        tape.at[(mat,mf,mt),'TEXT'] = "\n".join(TextOut) + '\n'
    return tape

def read_mf8_mt457(text):
    str_list = split_endf(text)
    i = 0
    out = {"MAT" : str_list["MAT"].iloc[0],
           "MF" : str_list["MF"].iloc[0],
           "MT" : str_list["MT"].iloc[0]}
    C, i = read_cont(str_list, i)
    out.update({"ZA" : C.C1, "AWR" : C.C2, "LIS" : C.L1, "LISO" : C.L2, "NST" :C.N1, "NSP" : C.N2})
    L, i = read_list(str_list, i)
    out.update({"HL" : L.C1, "DHL" : L.C2, "E" : L.B[::2], "DE" : L.B[1::2]})
    L, i = read_list(str_list, i)
    out.update({"SPI" : L.C1, "PAR" : L.C2, "DK" : []})
    if out["NST"] == 0:
        # Update list of decay modes when nuclide is radioactive
        out.update({ "DK" : [ {"RTYP" : RTYP, "RFS" : RFS, "Q" : Q, "DQ" : DQ, "BR" : BR, "DBR" : DBR } for RTYP, RFS, Q, DQ, BR, DBR in  zip(*[iter(L.B)]*6) ] })
    return out

def read_mf8_fy(text):
    str_list = split_endf(text)
    i = 0
    out = {"MAT" : str_list["MAT"].iloc[0],
           "MF" : str_list["MF"].iloc[0],
           "MT" : str_list["MT"].iloc[0]}
    C, i = read_cont(str_list, i)
    out.update({ "ZA" : C.C1, "AWR" : C.C2, "E" : {} })
    for j in range(C.L1):
        L, i = read_list(str_list, i)
        FY = [{ "ZAFP" : ZAFP, "FPS" : FPS, "YI" : YI, "DYI" : DYI } for ZAFP, FPS, YI, DYI in  zip(*[iter(L.B)]*4) ]
        out["E"].update({ L.C1 : { "FY" : FY } })
        if j > 0:
            out["E"][L.C1].update({ "I" : L.L1 })
    return out

def read_mf33_mt(text):
    str_list = split_endf(text)
    i = 0
    out = {"MAT" : str_list["MAT"].iloc[0],
           "MF" : str_list["MF"].iloc[0],
           "MT" : str_list["MT"].iloc[0]}
    C, i = read_cont(str_list, i)
    # Subsections are given as dictionary values.
    # Keys are MAT1*100+MT1
    out.update({"ZA" : C.C1, "AWR" : C.C2, "MTL" : C.L2, "SUB" : {}})
    for j in range(C.N2):
        sub = {}
        C, i = read_cont(str_list, i)
        sub.update({"XMF1" : C.C1, "XLFS1" : C.C2, "MAT1" : C.L1, "MT1" : C.L2})
        NC = C.N1
        NI = C.N2
        NCLIST = []
        for k in range(NC):
            C, i = read_cont(str_list, i)
            subsub = {"LTY" : C.L2}
            if subsub["LTY"] == 0:
                L, i = read_list(str_list, i)
                subsub.update({"E1" : L.C1, "E2" : L.C2,
                               "CI" : L.B[:L.N2], "XMTI" : L.B[L.N2:]})
                NCLIST.append(subsub)
            elif subsub["LTY"] in (1,2,3):
                L, i = read_list(str_list, i)
                subsub.update({"E1" : L.C1, "E2" : L.C2, "MATS" : L.L1, "MTS" : L.L2,
                               "XMFS" : L.B[0], "XLFSS" : L.B[1],
                               "EI" : L.B[2:2+L.N2], "WEI" : L.B[2+L.N2:]})
                NCLIST.append(subsub)
            NCLIST.append(subsub)
        sub.update({"NC" : NCLIST})
        NILIST = []
        for k in range(NI):
            L, i = read_list(str_list, i)
            subsub = {"LB" : L.L2}
            if subsub["LB"] in range(5):
                subsub.update({"LT" : L.L1, "NT" : L.NPL, "NP" : L.N2})
                if subsub["LT"] == 0:
                    subsub.update({"Ek" : L.B[::2], "Fk" : L.B[1::2]})
                else:
                    pdb.set_trace()
                    Nk = subsub["NP"] - subsub["LT"]
                    ARRk = L.B[:Nk]
                    ARRl = L.B[Nk:]
                    subsub.update({"Ek" : ARRk[:Nk/2], "Fk" : ARRk[Nk/2:],
                                   "El" : ARRl[:subsub["LT"]], "Fl" : ARRl[subsub["LT"]:]})
            elif subsub["LB"] == 5:
                subsub.update({"LS" : L.L1, "NT" : L.NPL, "NE" : L.N2,
                               "Ek" : L.B[:L.N2], "Fkk" : L.B[L.N2:]})
            elif subsub["LB"] == 6:
                subsub.update({"NT" : L.NPL, "NER" : L.N2, "NEC" : (L.NPL-1)//L.N2})
                subsub.update({"Ek" : L.B[:subsub["NER"]]})
                subsub.update({"El" : L.B[subsub["NER"]:subsub["NER"]+subsub["NEC"]]})
                subsub.update({"Fkl" : L.B[subsub["NER"]+subsub["NEC"]:]})
            elif subsub["LB"] in (8,9):
                subsub.update({"LT" : L.L1, "NT" : L.NPL, "NP" : L.N2})
                subsub.update({"Ek" : L.B[:subsub["NP"]], "Fk" : L.B[subsub["NP"]:]})
            else:
                pdb.set_trace()
            NILIST.append(subsub)
        sub.update({"NI" : NILIST})
        out["SUB"].update({sub["MAT1"]*1000+sub["MT1"] : sub})
    return out

def read_mf34_mt(str_list):
#    str_list = text.splitlines()
    i = 0
    out = {"MAT" : str_list["MAT"].iloc[0],
           "MF" : str_list["MF"].iloc[0],
           "MT" : str_list["MT"].iloc[0]}
#    out = {"MAT" : int(str_list[i][66:70]),
#           "MF" : int(str_list[i][70:72]),
#           "MT" : int(str_list[i][72:75])}
    C, i = read_cont(str_list, i)
    # Subsections are given as dictionary values.
    # Keys are MAT1*100+MT1
    out.update({"ZA" : C.C1, "AWR" : C.C2, "MTL" : C.L2, "SUB" : {}})
    for j in range(C.N2):
        sub = {}
        C, i = read_cont(str_list, i)
        sub.update({"XMF1" : C.C1, "XLFS1" : C.C2, "MAT1" : C.L1, "MT1" : C.L2})
        NC = C.N1
        NI = C.N2
        NCLIST = []
        for k in range(NC):
            C, i = read_cont(str_list, i)
            subsub = {"LTY" : C.L2}
            if subsub["LTY"] == 0:
                L, i = read_list(str_list, i)
                subsub.update({"E1" : L.C1, "E2" : L.C2,
                               "CI" : L.B[:L.N2], "XMTI" : L.B[L.N2:]})
                NCLIST.append(subsub)
            elif subsub["LTY"] in (1,2,3):
                L, i = read_list(str_list, i)
                subsub.update({"E1" : L.C1, "E2" : L.C2, "MATS" : L.L1, "MTS" : L.L2,
                               "XMFS" : L.B[0], "XLFSS" : L.B[1],
                               "EI" : L.B[2:2+L.N2], "WEI" : L.B[2+L.N2:]})
                NCLIST.append(subsub)
            NCLIST.append(subsub)
        sub.update({"NC" : NCLIST})
        NILIST = []
        for k in range(NI):
            L, i = read_list(str_list, i)
            subsub = {"LB" : L.L2}
            if subsub["LB"] in range(5):
                subsub.update({"LT" : L.L1, "NT" : L.NPL, "NP" : L.N2})
                if subsub["LT"] == 0:
                    subsub.update({"Ek" : L.B[::2], "Fk" : L.B[1::2]})
                else:
                    pdb.set_trace()
                    Nk = subsub["NP"] - subsub["LT"]
                    ARRk = L.B[:Nk]
                    ARRl = L.B[Nk:]
                    subsub.update({"Ek" : ARRk[:Nk/2], "Fk" : ARRk[Nk/2:],
                                   "El" : ARRl[:subsub["LT"]], "Fl" : ARRl[subsub["LT"]:]})
            elif subsub["LB"] == 5:
                subsub.update({"LS" : L.L1, "NT" : L.NPL, "NE" : L.N2,
                               "Ek" : L.B[:L.N2], "Fkk" : L.B[L.N2:]})
            elif subsub["LB"] == 6:
                subsub.update({"NT" : L.NPL, "NER" : L.N2, "NEC" : (L.NPL-1)//L.N2})
                subsub.update({"Ek" : L.B[:subsub["NER"]]})
                subsub.update({"El" : L.B[subsub["NER"]:subsub["NER"]+subsub["NEC"]]})
                subsub.update({"Fkl" : L.B[subsub["NER"]+subsub["NEC"]:]})
            elif subsub["LB"] in (8,9):
                subsub.update({"LT" : L.L1, "NT" : L.NPL, "NP" : L.N2})
                subsub.update({"Ek" : L.B[:subsub["NP"]], "Fk" : L.B[subsub["NP"]:]})
            else:
                pdb.set_trace()
            NILIST.append(subsub)
        sub.update({"NI" : NILIST})
        out["SUB"].update({sub["MAT1"]*1000+sub["MT1"] : sub})
    return out

def read_mf35_mt(str_list):
#    str_list = text.splitlines()
    i = 0
    out = {"MAT" : str_list["MAT"].iloc[0],
           "MF" : str_list["MF"].iloc[0],
           "MT" : str_list["MT"].iloc[0]}
#    out = {"MAT" : int(str_list[i][66:70]),
#           "MF" : int(str_list[i][70:72]),
#           "MT" : int(str_list[i][72:75])}
    C, i = read_cont(str_list, i)
    out.update({ "ZA" : C.C1, "AWR" : C.C2, "NK" : C.N1, "SUB" : []})
    for k in range(out["NK"]):
        L, i = read_list(str_list, i)
        out["SUB"].append({ "Elo" : L.C1, "Ehi" : L.C2, "NE" : L.N2, "Ek" : L.B[:L.N2],
           "Fkk" : L.B[L.N2:] })
    return out



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




class XsCov(pd.DataFrame):
    """
    columns =  (MATi,MTj) ... (MATm,MTn)
    index = E1, E2, ..., El
    """

    pass



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

#def reconstruct_xs(DfXs):
#    for mat in DfXs.columns.get_level_values("MAT").unique():
#        ListMT = deepcopy(DfXs[mat].columns)
#        for mt in DfXs[mat].columns.get_level_values("MT").unique():
#            daughters = [ x for x in range(800,850) if x in DfXs[mat].columns]
#            if daughters:
#                DfXs[mat,107] = DfXs[mat][daughters].sum(axis=1)
#            daughters = [ x for x in range(750,800) if x in DfXs[mat].columns]
#            if daughters:
#                DfXs[mat,106] = DfXs[mat][daughters].sum(axis=1)
#            daughters = [ x for x in range(700,750) if x in DfXs[mat].columns]
#            if daughters:
#                DfXs[mat,105] = DfXs[mat][daughters].sum(axis=1)
#            daughters = [ x for x in range(650,700) if x in DfXs[mat].columns]
#            if daughters:
#                DfXs[mat,104] = DfXs[mat][daughters].sum(axis=1)
#            daughters = [ x for x in range(600,650) if x in DfXs[mat].columns]
#            if daughters:
#                DfXs[mat,103] = DfXs[mat][daughters].sum(axis=1)
#            daughters = [ x for x in range(102,118) if x in DfXs[mat].columns]
#            if daughters:
#                DfXs[mat,101] = DfXs[mat][daughters].sum(axis=1)
#            daughters = [ x for x in (19,20,21,38) if x in DfXs[mat].columns]
#            if daughters:
#                DfXs[mat,18] = DfXs[mat][daughters].sum(axis=1)
#            daughters = [ x for x in (18,101) if x in DfXs[mat].columns]
#            if daughters:
#                DfXs[mat,27] = DfXs[mat][daughters].sum(axis=1)
#            daughters = [ x for x in range(50,92) if x in DfXs[mat].columns]
#            if daughters:
#                DfXs[mat,4] = DfXs[mat][daughters].sum(axis=1)
#            daughters = [ x for x in (4,5,11,16,17,*range(22,38),41,42,44,45) if x in DfXs[mat].columns]
#            if daughters:
#                DfXs[mat,3] = DfXs[mat][daughters].sum(axis=1)
#            daughters = [ x for x in (2,3) if x in DfXs[mat].columns]
#            if daughters:
#                DfXs[mat,1] = DfXs[mat][daughters].sum(axis=1)
#            daughters = [ x for x in (455,456) if x in DfXs[mat].columns]
#            if daughters:
#                DfXs[mat,452] = DfXs[mat][daughters].sum(axis=1)
#        # keep only mts present in the original file
#        todrop = [ x for x in DfXs[mat].columns if x not in ListMT ]
#        if todrop:
#            DfXs.drop(pd.MultiIndex.from_product([[mat], todrop]), axis=1, inplace=True)
#    return DfXs
#
#def update_xs(tapein, DfXs):
#    tape = pd.DataFrame(index=tapein.index.copy(), columns=tapein.columns.copy())
#    for k,row in tapein.iterrows():
#        tape.loc[k].DATA = deepcopy(row.DATA)
#        tape.loc[k].TEXT = deepcopy(row.TEXT)
#    for mat, mt in DfXs:
#        mf = 1 if mt in (452,455,456) else 3
#        name = 'NUBAR' if mt in (452,455,456) else 'XS'
#        if (mat, mf, mt) not in tape.index:
#            continue
#        # Cut threshold xs
#        iNotZero = next((i for i, x in enumerate(DfXs[mat,mt]) if x), None)
#        if iNotZero > 0:
#            SeriesXs = DfXs[mat,mt].iloc[iNotZero-1:]
#        else:
#            SeriesXs = DfXs[mat,mt]
#        # Assume all xs have only 1 interpolation region and it is linear
#        tape.DATA.loc[mat,mf,mt][name] = SeriesXs
#        tape.DATA.loc[mat,mf,mt]["NBT"] = [len(SeriesXs)]
#        tape.DATA.loc[mat,mf,mt]["INT"] = [2]
#    return tape

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


#from sandy.data_test import __path__ as paths
#tape = Endf6.from_file(os.path.join(paths[0], r"u235.endf")).process(keep_mf=[5]).get_chi()

#for i in open("../../list5_jeff").read().splitlines():
#    print("../../../../list5_jeff/"+i)
#    tape = Endf6.from_file("../../../../list5_jeff/"+i).process(keep_mf=[5])
#aaa=1

@pytest.fixture(scope="module")
def testPu9():
    from sandy.data_test import Pu9
    tape = Endf6.from_text("\n".join(Pu9.endf6))
    assert (tape.index.get_level_values("MAT").unique() == 9437).all()
    return tape

@pytest.mark.formats
@pytest.mark.endf6
@pytest.mark.xs
def test_read_xs(testPu9):
    testPu9.read_section(9437, 3, 102)

@pytest.mark.formats
@pytest.mark.endf6
@pytest.mark.info
def test_read_info(testPu9):
    testPu9.read_section(9437, 1, 451)

@pytest.mark.formats
@pytest.mark.endf6
@pytest.mark.nubar
def test_read_nubar(testPu9):
    testPu9.read_section(9437, 1, 452)
    testPu9.read_section(9437, 1, 455)
    testPu9.read_section(9437, 1, 456)

@pytest.mark.formats
@pytest.mark.endf6
@pytest.mark.chi
def test_read_chi(testPu9):
    testPu9.read_section(9437, 5, 18)

@pytest.mark.formats
@pytest.mark.endf6
@pytest.mark.xs
def test_extract_xs():
    from sandy.data_test import Pu9
    tape = Endf6.from_text("\n".join(Pu9.pendf))
    tape.get_xs(listmat=[9437], listmt=[1,2,102,4])
    xs = tape.get_xs(listmat=[125])
    assert xs.empty

@pytest.mark.formats
@pytest.mark.endf6
@pytest.mark.nubar
def test_extract_nubar(testPu9):
    testPu9.get_nubar(listmt=[452,456])

@pytest.mark.formats
@pytest.mark.endf6
@pytest.mark.chi
def test_extract_chi(testPu9):
    testPu9.get_edistr()