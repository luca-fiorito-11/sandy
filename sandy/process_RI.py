#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 12 17:28:35 2018

@author: fiorito_l
"""

from sandy.formats import errorr
import pandas as pd
List = []
for file in open("listtt").read().splitlines():
    toadd = errorr.Errorr.from_file(file).process().get_std().iloc[0].unstack("DATA")
    List.append(toadd)
frame = pd.concat(List).query("MT==1|MT==2|MT==18|MT==102")
import pdb
pdb.set_trace()
