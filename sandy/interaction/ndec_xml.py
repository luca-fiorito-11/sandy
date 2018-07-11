# -*- coding: utf-8 -*-
"""
Created on Tue Jul 10 10:09:55 2018

@author: fiorito_l
"""

import pandas as pd
import pdb
from xml.etree import ElementTree as ET

class NDECOutputs(pd.DataFrame):
    
    _COLUMNS = ["id", "format", "file"]

    def __init__(self, *args, **kwargs):
        kwargs.update({"columns" : self._COLUMNS})
        super().__init__(*args, **kwargs)
    
    def to_xml(self):
        outputs = ET.Element("outputs")
        for i,row in self.iterrows():
            attrib = {k : row[k] for k in self._COLUMNS[:-1]}
            ET.SubElement(outputs, "file", attrib=attrib)
        ET.dump(outputs)


class NDECMessages(pd.DataFrame):
    
    _COLUMNS = ["module", "file", "lines", "type", "message"]

    def __init__(self, *args, **kwargs):
        kwargs.update({"columns" : self._COLUMNS})
        super().__init__(*args, **kwargs)
        kwargs.update({"columns" : ["module", "type", "routine", "message"]})
        super().__init__(*args, **kwargs)

    def to_xml(self):
        messages = ET.Element("messages")
        for i,row in self.iterrows():
            attrib = {k : row[k] for k in self._COLUMNS[:-2]}
            ET.SubElement(messages, row["type"], attrib=attrib)
        ET.dump(messages)