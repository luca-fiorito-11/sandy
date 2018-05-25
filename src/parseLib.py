# -*- coding: utf-8 -*-
"""
Created on Fri May 25 16:58:11 2018

@author: fiorito_l
"""
import pandas as pd

class File(dict):

    def __init__(self, file):
        def read_float(x):
            try:
                return float(x[0] + x[1:].replace('+', 'E+').replace('-', 'E-'))
            except:
                return x
        widths = [11,11,11,11,11,11,4,2,3]
        columns = ["C1", "C2", "L1", "L2", "N1", "N2","MAT", "MF", "MT"]
        converters = dict(zip(columns[:6],[read_float]*6))
        tape = pd.read_fwf(file, widths=widths, names=columns, converters=converters).query("MAT>0 & MF>0 & MT>0")
        endf = tape.MAT.iloc[0]
        info = tape.query("MF==1 & MT==451")
        author = info.L2.iloc[4]
        awi = info.C1.iloc[2]
        awr = info.C2.iloc[0]
        dist = info.L1.iloc[5]
        Eval = info.L1.iloc[4]
        if 2 in tape.MF:
            self.update({"ehRes" : tape.query("MF==2 & MT==151").C2.iloc[2]})

        d = {(int(mat),int(mf),int(mt)):"\n".join([lines[i] for i in g.tolist()]) for (mat,mf,mt),g in tape.groupby(["MAT","MF","MT"]).groups.items()}
        tape = pd.DataFrame.from_dict(d, orient='index').rename(columns={0:"TEXT"})
        tape.index = pd.MultiIndex.from_tuples(tape.index, names=["MAT", "MF", "MT"])
        tape["DATA"] = None

File("1-H-3g.jeff33")