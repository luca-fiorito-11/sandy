# -*- coding: utf-8 -*-
"""
Created on Mon Jan 16 18:03:13 2017

@author: lfiorito
"""
import sys
import logging

def split(file):
    import re
    pattern = ".{74}0.{5}\n?"
    text = open(file).read()
    U = re.split(pattern, text)
    return list(filter(None, U)) # remove empty lines



class Chunk:
    """
    ``ENDF-6`` control file.
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

    def __init__(self, text):
        self.mat = int(text[66:70])
        self.mf = int(text[70:72])
        self.mt = int(text[72:75])
        self.text = list(map(lambda x: x[:66], text.splitlines()))
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


def list2dict(chunks):
    return 1


A=split("H1.txt")
B=Chunk(A[-1])
list2dict(A)

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