#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 21 09:22:55 2017

@author: lfiorito
"""
from os import symlink, mkdir, remove
from os.path import expanduser, abspath, join, basename, dirname, isfile
from shutil import copy, move
import numpy as np
import sys
from sandy.njoy import exe as xnjoy
from sandy.records import Samples
import logging
import pytest

global outobj
outobj = None


class SandySamples(Samples):
    r"""
    Matrix of output samples.
    """
    
    def write_stats(self):
        outobj.write(80*"+"+'\n')
        outobj.write("{:^80}\n".format('SANDY RESULTS'))
        outobj.write(80*"+"+'\n')
        outobj.write('ID    : ' + ' '.join([ '{:^15}'.format(x) for x in self.xx ]) + '\n')
        outobj.write('MEAN  : ' + ' '.join([ '{:^15.6e}'.format(x) for x in self.Mean ]) + '\n')
        outobj.write('STDEV : ' + ' '.join([ '{:^15.6e}'.format(x) for x in self.Std ]) + '\n')
        corr = self.Corr
        for i in range(len(self.xx)):
            if i==0:
                outobj.write('CORR  : ' + ' '.join([ '{:^15.6e}'.format(x) for x in corr[i] ]) + '\n')
            else:
                outobj.write('        ' + ' '.join([ '{:^15.6e}'.format(x) for x in corr[i] ]) + '\n')


class Folder:
    r"""
    Parent class for every working folder.
    """
    
    def __init__(self, path):
        self.CWD = path

    def mkdir(self):
        r"""
        Make directory.
        """
        mkdir(self.CWD)
    
    def move(self, src, name=None):
        r"""
        Move file to sample directory.
        """
        NAME = name if name else basename(src)
        move(src, join(self.CWD, NAME))

    def copy(self, src, name=None):
        r"""
        Copy file to sample directory.
        """
        NAME = name if name else basename(src)
        copy(src, join(self.CWD, NAME))

    def symlink(self, src, name=None):
        r"""
        Soft link file to sample directory.
        """
        NAME = name if name else basename(src)
        symlink(src, join(self.CWD, NAME))
    
    def get_file(self, basename):
        r"""
        Get path + name of file.
        
        Inputs:
            - :``basename``: :
                (string) file's basename
        """
        return join(self.CWD, basename)
    
    def run_njoy(self, *inputs):
        r"""
        Run ``NJOY`` reference calculation for this test using the 
        prewritten ``NJOY`` input files given in argument.
        
        Inputs:
            - :``inputs``: :
                (positional arguments) list of ``NJOY`` input files
        """
        for file in inputs:
            text = open(file).read()
            TEST.run_exe([xnjoy], text.encode(), cwd=self.CWD)

    def run_sandy(self, sandy_input, outdir=None):
        r"""
        Run ``SANDY`` sampling calculation for this test using the 
        prewritten ``SANDY`` input file given in argument.
        
        Inputs:
            - :``sandy_input``: :
                (string) ``SANDY`` input file
        """
        args = ["sandy", sandy_input]
        if outdir:
            args += ["-o", outdir]
        TEST.run_exe(args, cwd=self.CWD)



class Reference(Folder):
    
    def run_reference(self, tape20, *inputs):
        r"""
        Run ``NJOY`` reference calculation for this test using the 
        prewritten ``NJOY`` input files given in argument.
        
        Inputs:
            - :``inputs``: :
                (positional arguments) list of ``NJOY`` input files
        """
        self.mkdir()
        copy(tape20, join(self.CWD, "tape20"))
        self.run_njoy(*inputs)



class Sandy(Folder):
    
    def __init__(self, path, nsmp):
        super().__init__(path)
        self.FOLDERS = [ Isample(ismp, self.CWD) for ismp in range(1, nsmp+1) ]

    @property
    def nsmp(self):
        return len(self.FOLDERS)
    
    def init_samples(self, parms):
        nparms = len(parms)
        smp = SandySamples(parms, np.zeros((nparms, self.nsmp)))
        return smp
        


class Isample(Folder):
    r"""
    Calculation for individual sample.
    Each calculation instance is characterized by a sample number and a 
    working folder.
    """
    
    def __init__(self, ismp, path):
        self.ismp = ismp
        super().__init__(join(path, self.idir))
    
    @property
    def idir(self):
        r"""
        Directory name for the sample.
        """
        return "smp-{}".format(self.ismp)
    


class NJOYOUT(list):
    
    def __init__(self, file):
        with open(file) as f:
            super().__init__(f)
    
    def get_value(self, mf, mt, row, col):
        from sandy.rw_fortran import r_exp_format
        control = "{:>2}{:>3}".format(mf, mt)
        block = list(filter(lambda x: x[70:75]==control, self))
        if len(block) == 0:
            value = 0.
            return value
        line = block[row-1]
        string = line[11*(col-1):11*col]
        formatted = np.array(0.0, dtype='d')
        io_status = np.array(0, dtype='int')
        r_exp_format(string, formatted, io_status)
        if io_status != 0:
            raise NotImplementedError("Cannot read NJOY output value '{}'".format(string))
        value = float(formatted)
        return value



class TEST:
    r"""
    Parent class for every test.
    Most of its properties require attributes that should be defined in the 
    inheriting class.
    """
    
    @property
    def nsmp(self):
        r"""
        Number of samples. Get it from the ``SANDY`` input file.
        """
        from sandy.files import File
        options = File(self.input_sandy).load_yaml()
        return options['samples']
    
    @property
    def CWD(self):
        r"""
        Current working directory.
        """
        return self._path

    @property
    def ref(self):
        r"""
        Directory for reference ``NJOY`` calculation.
        """
        return join(self.tmpdir.strpath, "reference")
    
    @property
    def endf(self):
        r"""
        ``ENDF-6`` path + file.
        """
        return join(self.CWD, self._endf_file)
    
    @property
    def outdir(self):
        r"""
        ``SANDY`` output directory.
        """
        return join(self.tmpdir.strpath, "SANDY_OUT")
    
    @property
    def input_sandy(self):
        r"""
        ``SANDY`` path + input file.
        """
        return join(self.CWD, "input_sandy.inp")
    
    @property
    def njoy_inputs(self):
        r"""
        List of path + name of ``NJOY`` inputs
        """
        return [ join(self.CWD, i) for i in self._njoy_inputs ]

    @staticmethod
    def run_exe(args, stdin=None, cwd=None):
        from subprocess import Popen, PIPE
        stdout = stderr = PIPE
        process = Popen(args, shell=False, cwd=cwd, 
                        stdin=PIPE,
                        stdout=stdout, 
                        stderr=stderr)
        stdoutdata, stderrdata = process.communicate(input=stdin)
        if process.returncode != 0:
            msg = "Process status={}, cannot run '{}'"
            raise NotImplementedError(msg.format(process.returncode, " ".join(args)))
        
    def write_title(self):
        outobj.write(80*'='+'\n')
        outobj.write("{:^80}\n".format(self.__class__.__name__))
        outobj.write(80*'='+'\n')


    

class Test_mf31_pu244(TEST):

    from sandy.tests import mf31_pu244
    _path = dirname(mf31_pu244.__file__)
    _endf_file = "Pu244.dat"
    _njoy_inputs = ["input_moder.inp",
                    "input_reconr.inp",
                    "input_broadr.inp",
                    "input_groupr.inp",
                    "input_errorr.inp",
                    "input_groupr_1g.inp"]
    mf = 3
    mts = [452, 455, 456]
    
    def test_run(self, tmpdir):
        self.tmpdir = tmpdir
        self.TMPDIR = Folder(self.tmpdir.strpath)
        self.REF = Reference(self.ref)
        self.REF.run_reference(self.endf, *self.njoy_inputs[:-1])
        self.REF.output = NJOYOUT(self.REF.get_file("tape27"))
        self.TMPDIR.copy(self.endf)
        self.TMPDIR.OUTDIR = Sandy(self.outdir, self.nsmp)
        self.TMPDIR.run_sandy(self.input_sandy, outdir=self.TMPDIR.OUTDIR.CWD)
        for ISMP in self.TMPDIR.OUTDIR.FOLDERS:
            ifile = ISMP.get_file(self._endf_file + "." + ISMP.idir)
            ISMP.move(ifile, "tape20")
            ISMP.symlink(self.REF.get_file("tape23"))
            ISMP.run_njoy(self.njoy_inputs[0],
                          self.njoy_inputs[3],
                          self.njoy_inputs[4],
                          self.njoy_inputs[5])
            ISMP.output_errorr = NJOYOUT(ISMP.get_file("tape27"))
            ISMP.output_groupr = NJOYOUT(ISMP.get_file("tape26"))
        self.write_title()
        self.process_output_reference()
        self.process_output_sandy()
    
    def process_output_reference(self):
        from sandy.cov import Cov
        means = np.array([ self.REF.output.get_value(3, mt, 2, 1) for mt in self.mts ])
        size = len(means)
        cov = Cov(np.zeros([size, size]))
        indices = np.triu_indices(size)
        data = []
        for mt in self.mts:
            for i in range(1, size+1):
                row = 1 + 3 * i
                data.append(self.REF.output.get_value(33, mt, row, 1))
            size -= 1
        cov[indices] = data
        cov += np.triu(cov, 1).T
        outobj.write(80*"+"+'\n')
        outobj.write("{:^80}\n".format('NJOY RESULTS'))
        outobj.write(80*"+"+'\n')
        outobj.write('ID    : ' + ' '.join([ '{:^15}'.format(x) for x in self.mts ]) + '\n')
        outobj.write('MEAN  : ' + ' '.join([ '{:^15.6e}'.format(x) for x in means ]) + '\n')
        outobj.write('STDEV : ' + ' '.join([ '{:^15.6e}'.format(x) for x in cov.std*means ]) + '\n')
        corr = cov.corr
        for i in range(len(self.mts)):
            if i==0:
                outobj.write('CORR  : ' + ' '.join([ '{:^15.6e}'.format(x) for x in corr[i] ]) + '\n')
            else:
                outobj.write('        ' + ' '.join([ '{:^15.6e}'.format(x) for x in corr[i] ]) + '\n')

    def process_output_sandy(self):
        # CHECK CONVERGENCE
        smp = self.TMPDIR.OUTDIR.init_samples(self.mts)
        for ISMP in self.TMPDIR.OUTDIR.FOLDERS:
            for i,parm in enumerate(self.mts):
                smp[i, ISMP.ismp-1] = ISMP.output_errorr.get_value(3, parm, 2, 1)
            # Sum prompt and delayed to get total
            smp[0] = smp[1] + smp[2]
        smp.write_stats()
        # CHECK NUBAR PHYSICS
        outobj.write(80*"+"+'\n')
        outobj.write("{:^80}\n".format('CHECK CONSERVATION LAWS'))
        outobj.write(80*"+"+'\n')
        FORMAT = "{:^8}" + "".join(["{:^15.5f}" for i in range(len(self.mts)+1)]) + "{:^15}"
        HEAD   = "SAMPLE  "  + "".join(["{:^15}".format(i) for i in self.mts]) + "{:^15}".format("455 + 456")+ "{:^15}\n".format("PASSED")
        outobj.write(HEAD)
        lines = []
        for ISMP in self.TMPDIR.OUTDIR.FOLDERS:
            parms = [ ISMP.output_groupr.get_value(3, parm, 3, 2) for parm in self.mts ]
            SUM = sum(parms[1:])
            passed = 'True' if np.isclose(SUM, parms[0]) else 'False'
            lines.append(FORMAT.format(ISMP.ismp, *[x for x in parms], SUM, passed))
        outobj.write("\n".join(lines) + '\n')



class Test_mf35_cf252(TEST):
    
    from sandy.tests import mf35_cf252
    _path = dirname(mf35_cf252.__file__)
    _endf_file = "DCf252.dat"
    _njoy_inputs = ["input_moder.inp",
                    "input_reconr.inp",
                    "input_broadr.inp",
                    "input_groupr.inp",
                    "input_errorr.inp"]
    mf = 35
    mts = [1.000000E-5, 1.500000E+4, 1.050000E+6]
    
    def test_run(self, tmpdir):
        self.tmpdir = tmpdir
        self.TMPDIR = Folder(self.tmpdir.strpath)
        self.REF = Reference(self.ref)
        self.REF.run_reference(self.endf, *self.njoy_inputs)
        self.REF.output = NJOYOUT(self.REF.get_file("tape27"))
        self.TMPDIR.copy(self.endf)
        self.TMPDIR.OUTDIR = Sandy(self.outdir, self.nsmp)
        self.TMPDIR.run_sandy(self.input_sandy, outdir=self.TMPDIR.OUTDIR.CWD)
        for ISMP in self.TMPDIR.OUTDIR.FOLDERS:
            ifile = ISMP.get_file(self._endf_file + "." + ISMP.idir)
            ISMP.move(ifile, "tape20")
            ISMP.symlink(self.REF.get_file("tape23"))
            ISMP.run_njoy(self.njoy_inputs[0],
                          self.njoy_inputs[3],
                          self.njoy_inputs[4])
            ISMP.output = NJOYOUT(ISMP.get_file("tape27"))
        self.write_title()
        self.process_output_reference()
        self.process_output_sandy()
    
    def process_output_reference(self):
        from sandy.cov import Cov
        means = np.array([ self.REF.output.get_value(5, 18, 2, col) for col in (1,2,3) ])
        size = len(means)
        cov = Cov(np.zeros([size, size]))
        ibeg = 4
        for i in range(size):
            cov[i] = [ self.REF.output.get_value(35, 18, ibeg, col) for col in (1,2,3) ]
            ibeg += 2
        outobj.write(80*"+"+'\n')
        outobj.write("{:^80}\n".format('NJOY RESULTS'))
        outobj.write(80*"+"+'\n')
        outobj.write('ID    : ' + ' '.join([ '{:^15}'.format(x) for x in self.mts ]) + '\n')
        outobj.write('MEAN  : ' + ' '.join([ '{:^15.6e}'.format(x) for x in means ]) + '\n')
        outobj.write('STDEV : ' + ' '.join([ '{:^15.6e}'.format(x) for x in cov.std*means ]) + '\n')
        corr = cov.corr
        for i in range(len(self.mts)):
            if i==0:
                outobj.write('CORR  : ' + ' '.join([ '{:^15.6e}'.format(x) for x in corr[i] ]) + '\n')
            else:
                outobj.write('        ' + ' '.join([ '{:^15.6e}'.format(x) for x in corr[i] ]) + '\n')

    def process_output_sandy(self):
        smp = self.TMPDIR.OUTDIR.init_samples(self.mts)
        for ISMP in self.TMPDIR.OUTDIR.FOLDERS:
            for i in range(len(self.mts)):
                smp[i, ISMP.ismp-1] = ISMP.output.get_value(5, 18, 2, i+1)
        smp.write_stats()



class Test_mf32_th228(TEST):

    from sandy.tests import mf32_th228
    _path = dirname(mf32_th228.__file__)
    _endf_file = "Th228.dat"
    _njoy_inputs = ["input_moder.inp",
                    "input_reconr.inp",
                    "input_errorr.inp"]
    mf = 3
    mts = [1, 2, 102]
    
    def test_run(self, tmpdir):
        """
        Run test following this list of steps:
            * run the reference ``NJOY`` calculation;
            * run ``SANDY`` calculation;
            * for each sandy output:
                - ``SANDY`` sample must be copied to ``tape20``;
                - run ``NJOY`` in appropriate ``SANDY`` output folder;
            * analyze reference output;
            * analyze ``SANDY`` outputs.
        """
        self.tmpdir = tmpdir
        self.TMPDIR = Folder(self.tmpdir.strpath)
        self.REF = Reference(self.ref)
        self.REF.run_reference(self.endf, *self.njoy_inputs)
        self.REF.output = NJOYOUT(self.REF.get_file("tape27"))
        # MUST COPY THE ENDF FILE TO THE TEMPORARY FOLDER
        self.TMPDIR.copy(self.endf)
        self.TMPDIR.OUTDIR = Sandy(self.outdir, self.nsmp)
        self.TMPDIR.run_sandy(self.input_sandy, outdir=self.TMPDIR.OUTDIR.CWD)
        for ISMP in self.TMPDIR.OUTDIR.FOLDERS:
            ifile = ISMP.get_file(self._endf_file + "." + ISMP.idir)
            ISMP.move(ifile, "tape20")
            ISMP.run_njoy(*self.njoy_inputs)
            ISMP.output = NJOYOUT(ISMP.get_file("tape27"))
        self.write_title()
        self.process_output_reference()
        self.process_output_sandy()
    
    def process_output_reference(self):
        from sandy.cov import Cov
        means = np.array([ self.REF.output.get_value(3, mt, 2, 1) for mt in self.mts ])
        size = len(means)
        cov = Cov(np.zeros([size, size]))
        indices = np.triu_indices(size)
        data = []
        for mt in self.mts:
            for i in range(1, size+1):
                row = 1 + 3 * i
                data.append(self.REF.output.get_value(33, mt, row, 1))
            size -= 1
        cov[indices] = data
        cov += np.triu(cov, 1).T
        outobj.write(80*"+"+'\n')
        outobj.write("{:^80}\n".format('NJOY RESULTS'))
        outobj.write(80*"+"+'\n')
        outobj.write('ID    : ' + ' '.join([ '{:^15}'.format(x) for x in self.mts ]) + '\n')
        outobj.write('MEAN  : ' + ' '.join([ '{:^15.6e}'.format(x) for x in means ]) + '\n')
        outobj.write('STDEV : ' + ' '.join([ '{:^15.6e}'.format(x) for x in cov.std*means ]) + '\n')
        corr = cov.corr
        for i in range(len(self.mts)):
            if i==0:
                outobj.write('CORR  : ' + ' '.join([ '{:^15.6e}'.format(x) for x in corr[i] ]) + '\n')
            else:
                outobj.write('        ' + ' '.join([ '{:^15.6e}'.format(x) for x in corr[i] ]) + '\n')

    def process_output_sandy(self):
        r"""
        Loop sample folders and parameters (one-group ``MT1``, ``MT2`` and 
        ``MT102``) to get output samples.
        Then, write statistics.
        """
        smp = self.TMPDIR.OUTDIR.init_samples(self.mts)
        for ISMP in self.TMPDIR.OUTDIR.FOLDERS:
            for i,parm in enumerate(self.mts):
                smp[i, ISMP.ismp-1] = ISMP.output.get_value(3, parm, 2, 1)
        smp.write_stats()
    

def run_tests():
    from sandy.sandy_input import setup_logging
    global outobj
    setup_logging("info") # SET LOGGING TO INFO
    outfile = "sandy_tests.out"
    if isfile(outfile):
        remove(outfile)
    outobj = open(outfile, 'a')
    args = [abspath(__file__)]
    if len(sys.argv) > 1:
        args += sys.argv[1:]
    pytest.main(args)
    outobj.close()
    logging.info("\nThe tests' results can be found in folder 'sandy_tests.out' in the current working directory")
