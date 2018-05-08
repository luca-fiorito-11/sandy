# -*- coding: utf-8 -*-
"""
Created on Wed Apr 25 11:21:38 2018

@author: fiorito_l
"""

import sandy.endf6.files as e6
from os.path import join, realpath, dirname
from os import listdir
from sandy.endf6.records import read_float
import pandas as pd
import numpy as np


def get_masses(file):
    print(file)
    tape = e6.split2df_byZAM(file)
    tape["AWR"] = tape["TEXT"].apply(lambda x: float(read_float(x[11:22])))
    return tape.drop("TEXT", axis=1)

def check_masses(files):
    return [ get_masses(file) for file in files ]

def check_masses_jeff33(path):
    """
    Given a path with ENDF-6 files, compare the masses in each section of each file
    with those given in the RDD JEFF3.3 file.
    Produce a csv summary output.
    """
    from sandy.data_test import __file__ as td
    RDD = join(dirname(realpath(td)), r"RDD.jeff33")
    DfRDD = check_masses([RDD])[0].query("MF==1 & MT==451")#.rename(columns={"AWR": "AWR_RDD"})
    DfN = pd.concat(check_masses(map(lambda x : join(path, x), listdir(path))))
    C = DfN.merge(DfRDD.drop(["MF","MT"], axis=1), how="left", on=["ZAM"], suffixes=("_N","_RDD"))
    C["RATIO"] = (C.AWR_N/C.AWR_RDD).values
    C["MATCH"] = np.isclose(C["RATIO"], 1, rtol=1e-6)
    grouped = C[["ZAM","AWR_N"]].set_index("ZAM").groupby("ZAM")
    df = grouped.agg(['count', 'min', 'max', lambda x:x.value_counts().index[0]])
    df.columns = df.columns.droplevel(0)
    df = df.merge(DfRDD[['ZAM','AWR']], how='left', left_index=True, right_on=['ZAM'])
    df["Z"] = np.floor(df.ZAM/10000)
    df["A"] = np.floor((df.ZAM - df.Z*10000)/10)
    df["M"] = np.floor(df.ZAM - df.Z*10000 - df.A*10)
    df.rename(columns={"<lambda>" : "most freq", "AWR" : "ref"}, inplace=True)
    df.set_index(['Z','A','M'], inplace=True)
    df.drop("ZAM", axis=1, inplace=True)
    df["message"] = (np.isclose(df["min"], df["ref"], rtol=1e-6) & (np.isclose(df["min"], df["ref"], rtol=1e-6)))
    df["message"] = df.message.replace({True: "OK", False: "WARNING"})
    df.to_csv("jeff33_masses.csv")


#A=get_masses(r"X:\aaa\78-Pt-192g.jeff33")
check_masses_jeff33(r"X:\aaa")
