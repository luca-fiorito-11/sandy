# -*- coding: utf-8 -*-
"""
Created on Mon Jan 16 18:03:13 2017

@author: lfiorito
"""
import sys
import logging
import numpy as np
import fortranformat as ff
from collections import namedtuple

head_format = '(A66)'
cont_format = '(2E11.0,4I11)'
ilist_format = '(6I11)'
list_format = '(6E11.0)'
head_format_r = ff.FortranRecordReader(head_format)
cont_format_r = ff.FortranRecordReader(cont_format)
list_format_r = ff.FortranRecordReader(list_format)
ilist_format_r = ff.FortranRecordReader(ilist_format)

CONT = namedtuple('CONT', 'C1 C2 L1 L2 N1 N2')
LIST = namedtuple('LIST', 'C1 C2 L1 L2 NPL N2 B')
TAB1 = namedtuple('LIST', 'C1 C2 L1 L2 NR NP x y')

def read_cont(text, ipos):
    try:
        cont = CONT(*cont_format_r.read(text[ipos]))
        ipos += 1
        return cont, ipos
    except:
        sys.exit("ERROR: cannot read CONT at '{}'".format(text[ipos]))
            
def read_tab1(text, ipos):
    try:
        cont, ipos = read_cont(text, ipos)
        i = 0
        x = []
        while i < cont.N1:
            x.extend(ilist_format_r.read(text[ipos]))
            ipos += 1
            i += 6
        i = 0
        y = []
        while i < cont.N2:
            y.extend(list_format_r.read(text[ipos]))
            ipos += 1
            i += 6
        return
    except:
        sys.exit("ERROR: cannot read LIST at '{}'".format(text[ipos]))
    
    

def split(file):
    """
    Split ``ENDF-6`` file  into MFMT sections.
    """
    import re
    pattern = ".{74}0.{5}\n?"
    text = open(file).read()
    U = re.split(pattern, text)
    return list(filter(None, U)) # remove empty lines



class Chunk:
    """
    MFMT section
    """

    @property
    def n(self):
        r"""
        Number of lines of the `ENDF-6` file
        """
        return len(self.text)

    @property
    def line(self):
        """
        Return current line.
        """
        return self.text[self.i]

    def __init__(self, text, mat=None, mf=None, mt=None):
        from io import StringIO
        import pandas as pd
        self.mat = int(text[66:70])
        self.mf = int(text[70:72])
        self.mt = int(text[72:75])
        # first argument must be any object with a read() method (such as a StringIO)
        self.text = pd.read_fwf(StringIO(text), widths=6*[11]+[4,2,3,5],
                                names=("c1","c2","l1","l2","n1","n2","mat","mf","mt","ns"),
                                header=None)
        self.i = 0

    def __iter__(self):
        """
        Make the file an iterator.
        """
        return self

    def next(self):
        """
        Yield a new line of the `ENDF-6` file.
        """
        line = self.line
        self.i += 1
        return line

    def move(self, steps):
        """
        Move backward/forward a number of lines in the ``ENDF-6`` text.

        Inputs:
        - :``steps``: :
            (scalar integer) number of lines to jump
        """
        self.i += int(steps)

    def read_cont(self):
        """
        Read ``ENDF-6`` ``CONT`` record in formatted fortran.

        Outputs:
            - :``out``: :
                (tuple) content of ``CONT`` record

        Found error in:
            - n-17-Cl-035.jeff32
            - n-3-Li-007.jeff32
            - n-63-Eu-152.jeff32
            - n-63-Eu-153.jeff32
            - n-64-Gd-155.jeff32
            - n-77-Ir-193.jeff32
            - n-90-Th-229.jeff32
            - n-94-Pu-238.jeff32
            - n-94-Pu-241.jeff32
            - n-94-Pu-242.jeff32
            - n-97-Bk-250.jeff32
            - n-98-Cf-254.jeff32
        """
        try:
            out = (float(self.line.c1), float(self.line.c2), float(self.line.l1),
                   float(self.line.l2), float(self.line.n1), float(self.line.n2))
            self.i += 1
            return out
        except:
            sys.exit("ERROR: line is not CONT.\n{}".format(self.line))


def list2dict(chunks):
    return 1


A=split("H1.txt")
XS=A[2].splitlines()
ii=0
CC, ii = read_cont(XS,ii)
read_tab1(XS,ii)
B=Chunk(A[-1])
import re
line = re.sub('\n', '', A[1])
AAA=re.findall('.{11}', 'line')
import pandas as pd
line = open('012503001').read()
line = re.sub('\n', '', line)
blocks = re.findall('.{11}', line)
df = pd.read_fwf('012503001', widths=6*[11]+[4,2,3,5],
                 names=("C1","C2","L1","L2","N1","N2","MAT","MF","MT","NS"),
                 header=None)
sys.exit()
df.fillna(0, inplace=True)

class File:
    """ General text file """

    @staticmethod
    def isfile(filename):
        """
        Check if file exists.

        Inputs:
            - filename :
                (string) path+name of the file
        """
        from os.path import isfile
        if not isfile(filename):
            logging.error("FILE : File '{}' does not exist".format(filename))
            sys.exit()

    def __init__(self, filename):
        """
        Initialize file.

        Inputs:
            - filename :
                (string) path+name of the file

        ..Important::
            Use encoding="ascii" and errors="surrogateescape" when opening the
            file in order to process files in ASCII compatible encoding such as
            `utf-8` and `latin-1`.
        """
        self.filename = filename
        self.f = open(self.filename, encoding="ascii", errors="surrogateescape")

    @property
    def filename(self):
        """
        Absolute file path + name.
        """
        return self._filename

    @filename.setter
    def filename(self, filename):
        from os.path import expanduser, abspath
        _f = abspath(expanduser(filename))
        File.isfile(_f)
        self._filename = _f

    @property
    def path(self):
        """
        Absolute file path.
        """
        from os.path import split
        return split(self.filename)[0]

    @property
    def name(self):
        """
        File basename.
        """
        from os.path import split
        return split(self.filename)[1]

    def read(self):
        """
        Load the content of input object `stream`.

        Outputs:
            - text :
                (list of strings) text in the file
        """
        try:
            text = self.f.readlines()
            return text
        except:
            logging.error("FILE : Cannot read file '{}'".format(self.filename))
            sys.exit()

    def load_yaml(self, msg="ERROR! Cannot yaml.load input file"):
        """
        Description
        ===========
        Load the content of input object *stream*.

        Outputs
        ======
         - *load*: text in yaml format
        """
        import yaml
        try:
            load = yaml.load(self.f)
            return load
        except yaml.YAMLError as exc:
            sys.exit(msg + " '{}'".format(self.filename))

#import pytest
#
#class TestFile:
#
#    fileyaml = '../test_objects/input.yaml'
#    filenotyaml = '../test_objects/input.not_yaml'
#    filenotexists = 'nonexistingfile'
#
#    def test_is_not_file(self):
#        from files import File
#        with pytest.raises(SystemExit):
#            File.isfile(self.filenotexists)
#
#    def test_is_file(self):
#        from files import File
#        File.isfile(self.fileyaml)
#
#    def test_init_file(self):
#        from files import File
#        F = File(self.fileyaml)
#        assert F.name == 'input.yaml'
#
#    def test_read_file(self):
#        from files import File
#        F = File(self.fileyaml)
#        assert F.read() == open(self.fileyaml).readlines()
#
#    def test_load_yaml_file(self):
#        from files import File
#        F = File(self.fileyaml)
#        load = F.load_yaml()
#        assert 'replace' in load
#        assert 'mat' in load['replace']
#        assert load['replace']['mat'] == 2631
#        assert 'endf1' in load['replace']
#        assert load['replace']['endf1'] == 'file1'
#        assert 'endf2' in load['replace']
#        assert load['replace']['endf2'] == 'file2'
#        assert 'mf' in load['replace']
#        assert load['replace']['mf'] == 3
#        assert 'mt' in load['replace']
#        assert load['replace']['mt'] == 102
#        assert 'emin' in load['replace']
#        assert load['replace']['emin'] == 1e6
#        assert 'emax' in load['replace']
#        assert load['replace']['emax'] == 10e6
#
#    def test_load_notyaml_file(self):
#        from files import File
#        F = File(self.filenotyaml)
#        with pytest.raises(SystemExit):
#            F.load_yaml()