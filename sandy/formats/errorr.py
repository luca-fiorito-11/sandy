# -*- coding: utf-8 -*-
"""
Created on Fri May 11 15:08:25 2018

@author: fiorito_l
"""
import pandas as pd
import sys, pytest
import numpy as np

#def process_errorr_section(text, keep_mf=None, keep_mt=None):
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
#    elif mf == 3:
#        return read_mf3_mt(text)
#    elif mf == 33:
#        return read_mf33_mt(text)
#    else:
#        return None


class Errorr(pd.DataFrame):

    Format = "errorr"

    @classmethod
    def from_file(cls, file):
        from .endf6 import Endf6
        return cls(Endf6.from_file(file))

    @classmethod
    def from_text(cls, text):
        from .endf6 import Endf6
        return cls(Endf6.from_text(text))

    def read_section(self, mat, mf, mt):
        """
        Parse MAT/MF/MT section
        """
        if mf == 1:
            from .MF1 import read_errorr as read
        elif mf == 3:
            from .MF3 import read_errorr as read
        elif mf == 33 or mf == 31:
            from .MF33 import read_errorr as read
        else:
            sys.exit("ERROR: SANDY cannot parse section MAT{}/MF{}/MT{}".format(mat,mf,mt))
        return read(self.loc[mat,mf,mt].TEXT)

    def get_xs(self, listmat=None, listmt=None):
        """
        Extract xs from errorr file into Xs instance.
        """
        from .utils import Xs
        query = "MF==3"
        if listmat is not None:
            query_mats = " | ".join(["MAT=={}".format(x) for x in listmat])
            query += " & ({})".format(query_mats)
        if listmt is not None:
            query_mts = " | ".join(["MT=={}".format(x) for x in listmt])
            query += " & ({})".format(query_mts)
        tape = self.query(query)
        mat = self.index.get_level_values("MAT")[0]
        eg = self.read_section(mat,1,451)["EG"]
        ListXs = []
        for ix,text in tape.TEXT.iteritems():
            mat,mf,mt = ix
            X = self.read_section(*ix)
            xs = pd.Series(X["XS"], index=eg[:-1], name=(X["MAT"],X["MT"])).rename_axis("E").to_frame()
            ListXs.append(xs)
        if not ListXs:
            return pd.DataFrame()
        # Use concat instead of merge because indexes are the same
        frame = pd.concat(ListXs, axis=1).reindex(eg, method="ffill")
        return Xs(frame)

    def get_xs_cov(self, listmat=None, listmt=None):
        """
        Extract xs covariances from errorr file into XsCov instance.
        """
        from .utils import XsCov
        query = "(MF==33 | MF==31)"
        if listmat is not None:
            query_mats = " | ".join(["MAT=={}".format(x) for x in listmat])
            query += " & ({})".format(query_mats)
        if listmt is not None:
            query_mts = " | ".join(["MT=={}".format(x) for x in listmt])
            query += " & ({})".format(query_mts)
        tape = self.query(query)
        mat = self.index.get_level_values("MAT")[0]
        eg = self.read_section(mat,1,451)["EG"]
        List = []
        for ix,text in tape.TEXT.iteritems():
            mat,mf,mt = ix
            X = self.read_section(*ix)
            for mt1,y in X["RP"].items():
                List.append([mat, X["MT"], mat, mt1, y])
        frame = pd.DataFrame(List, columns=('MAT', 'MT','MAT1', 'MT1', 'COV'))
        MI = [(mat,mt,e) for mat,mt in sorted(set(zip(frame.MAT, frame.MT))) for e in eg]
        index = pd.MultiIndex.from_tuples(MI, names=("MAT", "MT", "E"))
        # initialize union matrix
        matrix = np.zeros((len(index),len(index)))
        for i,row in frame.iterrows():
            ix = index.get_loc((row.MAT,row.MT))
            ix1 = index.get_loc((row.MAT1,row.MT1))
            matrix[ix.start:ix.stop-1,ix1.start:ix1.stop-1] = row.COV
        i_lower = np.tril_indices(len(index), -1)
        matrix[i_lower] = matrix.T[i_lower]  # make the matrix symmetric
        return XsCov(matrix, index=index, columns=index)

    def get_std(self):
        """
        Extract xs and std from errorr file into dataframe:
            index = energy
            columns = (MAT, MT, DATA)
        """
        xs = self.get_xs()
        cov = self.get_cov()
        stdvals = np.sqrt(np.diag(cov.values))
        xsvals =  xs.values.T.flatten()
        frame = pd.DataFrame.from_dict({"XS" : xsvals, "STD" : stdvals})
        frame.columns.name = "DATA"
        frame.index = cov.index
        frame = frame.unstack(level=["MAT","MT"])
        frame.columns = frame.columns.reorder_levels(["MAT","MT","DATA"])
        return frame


@pytest.fixture(scope="module")
def testH1():
    from sandy.data_test import H1
    tape = Errorr.from_text("\n".join(H1.errorr))
    assert (tape.index.get_level_values("MAT").unique() == 125).all()
    return tape

@pytest.mark.formats
@pytest.mark.errorr
@pytest.mark.info
def test_read_info(testH1):
    testH1.read_section(125, 1, 451)

@pytest.mark.formats
@pytest.mark.errorr
@pytest.mark.xs
def test_read_xs(testH1):
    testH1.read_section(125, 3, 102)

@pytest.mark.formats
@pytest.mark.errorr
@pytest.mark.cov
@pytest.mark.xs
def test_read_xs_cov(testH1):
    testH1.read_section(125, 33, 102)

@pytest.mark.formats
@pytest.mark.errorr
@pytest.mark.xs
def test_extract_xs(testH1):
    testH1.get_xs(listmat=[125], listmt=[1,2])

@pytest.mark.formats
@pytest.mark.errorr
@pytest.mark.cov
@pytest.mark.xs
def test_extract_xs_cov(testH1):
    testH1.get_xs_cov()