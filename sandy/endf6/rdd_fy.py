# -*- coding: utf-8 -*-
"""
Created on Fri Apr 20 16:42:45 2018

@author: lucaf
"""
from sandy.endf6 import files as e6
import pandas as pd
from os.path import join, dirname, realpath
import numpy as np
from e6 import read_float

class DecayDataFile():
    
    def __init__(self, file):
        self.file = file

    def extract_lambdas(self):
        ListLambdas = []
        for chunk in e6.split(self.file):
            mf = int(chunk[70:72])
            mt = int(chunk[72:75])
            if mf == 8 and mt == 457:
                # can do it faster without processing the section
                data = e6.process_endf_section(chunk)
                DC = np.log(2.)/data["HL"] if data["HL"] != 0 else 0
                DDC = np.log(2.) * data["DHL"]/data["HL"] * DC if data["HL"] != 0 else 0
                ListLambdas.append([ data["ZA"], data["LISO"], DC, DDC ])
        return pd.DataFrame(ListLambdas, columns=["ZA","LISO","DC", "DDC"])
    
    def extract_qmatrix(self):
        aaa=1
        pass

def test_jeff33_decay_constants():
    from sandy.data_test import __file__ as td
    from sandy import __file__ as sd
    sd = dirname(realpath(sd))
    td = dirname(realpath(td))
    FileRDD = join(td, r"RDD.jeff33")
    RDD = DecayDataFile(FileRDD)
    A = RDD.extract_lambdas()
    aaa=1

def test_jeff33_qmatrix():
    from sandy.data_test import __file__ as td
    from sandy import __file__ as sd
    sd = dirname(realpath(sd))
    td = dirname(realpath(td))
    FileRDD = join(td, r"RDD.jeff33")
    RDD = DecayDataFile(FileRDD)
    A = RDD.extract_lambdas()
    aaa=1


test_jeff33()