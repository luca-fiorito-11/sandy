# -*- coding: utf-8 -*-
"""
Created on Fri Aug 24 10:06:34 2018

@author: Fiorito_L
"""

from .utils import BaseFile
from ..settings import SandyError

__author__ = "Luca Fiorito"
__all__ = ["Groupr"]

class Groupr(BaseFile):

    Format = "groupr"

    def read_section(self, mat, mf, mt):
        """
        Parse MAT/MF/MT section
        """
        if mf == 1:
            from .MF1 import read_groupr as read
        elif mf == 3:
            from .MF3 import read_groupr as read
        else:
            raise SandyError("SANDY cannot parse section MAT{}/MF{}/MT{}".format(mat,mf,mt))
        if (mat,mf,mt) not in self.index:
            raise SandyError("section MAT{}/MF{}/MT{} is not in tape".format(mat,mf,mt))
        return read(self.loc[mat,mf,mt].TEXT)