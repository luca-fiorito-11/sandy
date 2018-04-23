# -*- coding: utf-8 -*-
"""
Created on Fri Apr 20 16:42:45 2018

@author: lucaf
"""
from sandy.endf6 import files as e6
import pandas as pd
from os.path import join, dirname, realpath
import numpy as np
import scipy as sp
import copy
import sys
from sandy.tests import TimeDecorator
#from e6 import read_float

class DecayDataFile():

    def __init__(self, file):
        self.file = file

    def extract_lambdas(self):
        ListLambdas = []
        for chunk in e6.split(self.file):
            mf = int(chunk[70:72])
            mt = int(chunk[72:75])
            if mf == 8 and mt == 457:
                # I can do it faster without processing the section
                data = e6.read_mf8_mt457(chunk)
                DC = np.log(2.)/data["HL"] if data["HL"] != 0 else 0
                DDC = np.log(2.) * data["DHL"]/data["HL"] * DC if data["HL"] != 0 else 0
                ListLambdas.append([ data["ZA"], data["LISO"], DC, DDC ])
        return pd.DataFrame(ListLambdas, columns=["ZA","LISO","DC", "DDC"])

    @TimeDecorator
    def extract_decaychains(self):
        ListDC = []
        for chunk in e6.split(self.file):
            mf = int(chunk[70:72])
            mt = int(chunk[72:75])
            if mf == 8 and mt == 457:
                data = e6.read_mf8_mt457(chunk)
                za_par = int(data["ZA"])
                liso_par = data["LISO"] # Isomeric state flag for parent nuclide
                if data["NST"] == 1:
                    ListDC.append([za_par, liso_par])
                for dmode in data["DK"]:
                    za_dau = copy.copy(za_par)
                    liso_dau = dmode["RFS"] # Isomeric state flag for daughter nuclide
                    for prtc in map(int,str(dmode["RTYP"]).replace(".","").replace("0", "")):
                        if prtc == 1: # beta decay
                            mode = "Beta-"
                            za_dau += 1000
                        elif prtc == 2: # Electron capture and/or positron emission
                            mode = "Beta+"
                            za_dau -= 1000
                        elif prtc == 3: # Isomeric transition
                            mode = "I.T."
                            pass
                        elif prtc == 4: # Alpha decay
                            mode = "Alpha"
                            za_dau -= 2004
                            ListDC.append([za_par, liso_par, 2004, 0, data["HL"], data["DHL"], dmode["BR"], dmode["DBR"], mode])
                        elif prtc == 5: # Neutron emission (not delayed neutron decay)
                            mode = "n"
                            za_dau -= 1
                            ListDC.append([za_par, liso_par, 1, 0, data["HL"], data["DHL"], dmode["BR"], dmode["DBR"], mode])
                        elif prtc == 6: # Spontaneous fission
                            mode = "S.F."
                            pass
                        elif prtc == 7: # Proton emission
                            mode = "p"
                            za_dau -= 1001
                            ListDC.append([za_par, liso_par, 1001, 0, data["HL"], data["DHL"], dmode["BR"], dmode["DBR"], mode])
                        else:
                            sys.exit("ERROR: Unknown decay mode")
                        ListDC.append([za_par, liso_par, za_dau, liso_dau, data["HL"], data["DHL"], dmode["BR"], dmode["DBR"], mode])
        df = pd.DataFrame(ListDC, columns=("ZA_PARENT", "LISO_PARENT", "ZA_DAUGHTER", "LISO_DAUGHTER", "HL", "DHL", "BR", "DBR", "MODE"))
        df["DC"] = np.log(2.)/df.HL
        df["DDC"] = np.log(2.) * df.DHL/df.HL * df.DC
        return df

    @TimeDecorator
    def extract_bmatrix(self, timelimit=np.inf):
        secyear = 60*60*24*365.25
        tl = secyear*timelimit
#        query = "HL > 0 & HL < {}".format(tl) if with_stable else "HL < {}".format(tl)
        df = self.extract_decaychains()#.query(query)
        df["ZAMP"] = df.ZA_PARENT*10 + df.LISO_PARENT
        df["ZAMD"] = df.ZA_DAUGHTER*10 + df.LISO_DAUGHTER
        IDS = df.ZAMP.unique()
        df.loc[df.HL > tl, ["BR", "DBR"]] = 0
        xt = pd.crosstab(
                index=[df.ZAMP],
                columns=[df.ZAMD],
                values=df.BR,
                aggfunc=np.sum,
                ).reindex(IDS, columns=IDS).fillna(0)
        xt.index = xt.columns = df.groupby(["ZA_PARENT","LISO_PARENT"]).nunique().index
        return xt

    @TimeDecorator
    def extract_qmatrix(self):
        bm = self.extract_bmatrix()
        B = np.identity(len(bm)) - bm.as_matrix()
        Q = np.linalg.pinv(B)
        return pd.DataFrame(Q, index=bm.index, columns=bm.columns)
#        LU = sp.sparse.linalg.spilu(sp.sparse.csc_matrix(B))
#        Q1 = sp.sparse.csc_matrix(LU.solve(np.eye(B.shape[0]))).todense()
#
#        B = sp.sparse.csc_matrix(B)
#        lu = sp.sparse.linalg.splu(B)
#        eye = np.eye(B.shape[0])
#        Q = sp.sparse.csc_matrix(lu.solve(eye))

def test_jeff33_decay_constants():
    """
    Extract decay constants from jeff file and write in csv format.
    """
    from sandy.data_test import __file__ as td
    from sandy import __file__ as sd
    sd = dirname(realpath(sd))
    td = dirname(realpath(td))
    A = DecayDataFile( join(td, r"RDD.jeff33") ).extract_lambdas()
    aaa=1

def test_jeff33_qmatrix():
    """
    Extract qmatrix from jeff file and write in csv format.
    """
    from sandy.data_test import __file__ as td
    from sandy import __file__ as sd
    sd = dirname(realpath(sd))
    td = dirname(realpath(td))
    Q = DecayDataFile( join(td, r"RDD.jeff33") ).extract_qmatrix()
    aaa=1


test_jeff33_qmatrix()