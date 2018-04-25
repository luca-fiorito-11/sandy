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

def B2Q(dfB):
    """
    Convert dataframe B-matrix into dataframe Q-matrix.
    """
    B = np.identity(len(dfB)) - dfB.as_matrix()
    Q = np.linalg.pinv(B)
    return pd.DataFrame(Q, index=dfB.index, columns=dfB.columns)

class RDDFile():

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
    def extract_bmatrix(self, timelimit=np.inf, n=False):
        secyear = 60*60*24*365.25
        tl = secyear*timelimit
#        query = "HL > 0 & HL < {}".format(tl) if with_stable else "HL < {}".format(tl)
        df = self.extract_decaychains()#.query(query)
        df["ZAMP"] = df.ZA_PARENT*10 + df.LISO_PARENT
        df["ZAMD"] = df.ZA_DAUGHTER*10 + df.LISO_DAUGHTER
        if not n:
            df = df[df.ZA_PARENT != 1]
        IDS = df.ZAMP.unique()
        df.loc[df.HL > tl, ["BR", "DBR"]] = 0
        xt = pd.crosstab(
                index=[df.ZAMP],
                columns=[df.ZAMD],
                values=df.BR,
                aggfunc=np.sum,
                ).reindex(IDS, columns=IDS).fillna(0)
        xt.index = xt.columns = df.groupby(["ZA_PARENT","LISO_PARENT"]).nunique().index
        return xt.T

    @TimeDecorator
    def extract_qmatrix(self, timelimit=np.inf, n=False):
        bm = self.extract_bmatrix(timelimit=timelimit)
        return B2Q(bm)
#        B = np.identity(len(bm)) - bm.as_matrix()
#        Q = np.linalg.pinv(B)
#        return pd.DataFrame(Q, index=bm.index, columns=bm.columns)
#        LU = sp.sparse.linalg.spilu(sp.sparse.csc_matrix(B))
#        Q1 = sp.sparse.csc_matrix(LU.solve(np.eye(B.shape[0]))).todense()
#
#        B = sp.sparse.csc_matrix(B)
#        lu = sp.sparse.linalg.splu(B)
#        eye = np.eye(B.shape[0])
#        Q = sp.sparse.csc_matrix(lu.solve(eye))


class FYFile():

    def __init__(self, file):
        self.file = file

    @TimeDecorator
    def extract_fy(self):
        IFY = []; CFY = []
        for chunk in e6.split(self.file):
            mf = int(chunk[70:72])
            mt = int(chunk[72:75])
            if mf == 8 and mt == 454:
                data = e6.read_mf8_fy(chunk)
                for e in data["E"]:
                    IFY.extend([ dict({ "ZAP" : data["ZA"], "E" : e }, **fy) for fy in data["E"][e]["FY"] ])
            elif mf == 8 and mt == 459:
                data = e6.read_mf8_fy(chunk)
                for e in data["E"]:
                    CFY.extend([ dict({ "ZAP" : data["ZA"], "E" : e }, **fy) for fy in data["E"][e]["FY"] ])
        return pd.DataFrame(IFY), pd.DataFrame(CFY)

def check_rdd_fy(RDD, FY):
    from numpy.linalg import lstsq
    from scipy.stats import chisquare
    B = RDDFile( RDD ).extract_bmatrix(timelimit=1000)
    index = B.index.get_level_values("ZA_PARENT")*10 + B.index.get_level_values("LISO_PARENT")
    IFY, CFY = FYFile( FY ).extract_fy()
    # Loop fissioning systems
    for (zap,e), dfi in IFY.groupby(["ZAP", "E"]):
        dfc = CFY.query("ZAP=={} & E=={}".format(zap,e))
        ify = dfi.set_index(dfi.ZAFP*10+dfi.FPS).reindex(index).YI.fillna(0)
        cfy = dfc.set_index(dfc.ZAFP*10+dfc.FPS).reindex(index).YI.fillna(0)
        C = cfy.to_frame().rename(columns={'YI' : 'CY'})
        C['IY'] = ify.values
        C['CY_calc'], res, rank, sing = lstsq((np.identity(len(B)) - B.as_matrix()), C.IY.values)
        C["diff"] = np.abs(C.CY_calc/C.CY-1).values
        C["Z"] = pd.Series(C.index.values/10000, dtype=int).values
        C["A"] = pd.Series((C.index.values-C.Z.values*10000)/10, dtype=int).values
        chi2, p_value = chisquare(C.CY.loc[C.CY>0], f_exp=C.CY_calc.loc[C.CY>0])
        print(zap,e,chi2,p_value)

def test_jeff33_decay_constants():
    """
    Extract decay constants from jeff file and write in csv format.
    """
    from sandy.data_test import __file__ as td
    from sandy import __file__ as sd
    sd = dirname(realpath(sd))
    td = dirname(realpath(td))
    A = RDDFile( join(td, r"RDD.jeff33") ).extract_lambdas()
    aaa=1

def test_jeff33_rdd_fy():
    """
    ISO E CHI^2 PVALUE
    90232.0 400000.0 9.09849076297 1.0
    90232.0 14000000.0 111.981926453 1.0
    92233.0 0.0253 1.76057717685 1.0
    92233.0 400000.0 1.13547668601 1.0
    92233.0 14000000.0 1.98728653841 1.0
    92234.0 400000.0 0.653759104852 1.0
    92235.0 0.0253 0.929158792529 1.0
    92235.0 400000.0 1.58365367609 1.0
    92235.0 14000000.0 4.4516291131 1.0
    92236.0 0.0253 329.909545277 1.0
    92236.0 400000.0 31.4154974493 1.0
    92238.0 400000.0 128.797656121 1.0
    92238.0 14000000.0 6.84650480105 1.0
    93237.0 0.0253 1.8310359522 1.0
    93237.0 400000.0 4.19247484827 1.0
    93238.0 0.0253 10.4958975282 1.0
    93238.0 400000.0 6.82070527688 1.0
    94238.0 0.0253 0.0253452224653 1.0
    94238.0 400000.0 0.196266006335 1.0
    94239.0 0.0253 0.752322322921 1.0
    94239.0 400000.0 2.5143568177 1.0
    94240.0 400000.0 4.61164409854 1.0
    94241.0 0.0253 0.125947777869 1.0
    94241.0 400000.0 18.8692182249 1.0
    94242.0 400000.0 26.8780592723 1.0
    95241.0 0.0253 0.668792606384 1.0
    95241.0 400000.0 1.08959344195 1.0
    95242.0 0.0253 2.79542677863 1.0
    95242.0 400000.0 5.11511820113 1.0
    95243.0 0.0253 3.10548964903 1.0
    95243.0 400000.0 4.94311274336 1.0
    96243.0 0.0253 0.122632005282 1.0
    96243.0 400000.0 0.234333804263 1.0
    96244.0 0.0253 0.11360176429 1.0
    96244.0 400000.0 1.46044874563 1.0
    96245.0 0.0253 0.311911363568 1.0
    96245.0 400000.0 0.509813267668 1.0
    """
    from sandy.data_test import __file__ as td
    check_rdd_fy(join(dirname(realpath(td)), r"RDD.jeff33"),
                 join(dirname(realpath(td)), r"FY.jeff33"))


test_jeff33_rdd_fy()