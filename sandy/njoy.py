# -*- coding: utf-8 -*-
"""
Created on Thu Jun 18 17:36:27 2015

@author: lfiorito
"""
import sys
import numpy as np
import os
import logging
from os.path import dirname, join
import sandy.NJOY


exe = join(dirname(sandy.NJOY.__file__), 'njoy.exe')

class FileNJOY:

    @staticmethod
    def copy_to_tape(src, tape, dst=None):
        """
        Description
        ===========
        Copy `ENDF-6` or `PENDF` inputfile to `NJOY` tape file, e.g. tape20.
        
        Inputs
        ======
         - *src*: `ENDF-6` or `PENDF` inputfile's name
         - *tape*: tape number
         - *dts*: path where the `NJOY` tape file must be copied
        """
        import shutil
        if dst is None:
            dst = ''
        tape = 'tape{}'.format(tape)
        dst = os.path.join(dst, tape)
        shutil.copy(src, dst)
    
    def __init__(self, name):
        r"""
        Initialize `NJOY` input file.
        
        Inputs:
            - name :
                (string) `NJOY` inputfile's name
        """
        from os.path import join
        from sandy.sandy_input import outdir
        self.file = join(outdir, name)
        self.text = ""
    
    def reconr(self, nendf, npendf, mat, reconr_err=0.001, reconr_grid=[], 
               reconr_cards=[], reconr_tempr=0, reconr_errmax=None, 
               reconr_errint=None, **options):
        """
        Define input file for `RECONR` module of `NJOY` and add it to 
        instance text.
        
        Inputs:
            - nendf :
                unit for endf tape
            - npendf :
                unit for pendf tape
            - mat :
                material to be reconstructed.
            - reconr_err :
                fractional reconstruction tolerance used when
                resonance-integral error criterion (see errint)
                is not satisfied (default=0.001).
            - reconr_grid :
                users energy grid points (default None).
            - reconr_cards :
                cards of descriptive comments for mt451 (default None).
            - reconr_tempr :
                reconstruction temperature (deg kelvin) (default=0).
            - reconr_errmax :
                fractional reconstruction tolerance used when 
                resonance-integral error criterion is satisfied 
                (errmax.ge.err, default=10*err).
            - reconr_errint :
                maximum resonance-integral error (in barns)
                per grid point (default=err/20000)
                (note: the max cross section difference for linearization, 
                errlim, and for reconstruction, errmin, are also tied to 
                errint.
                To get maximum accuracy, set errint to a very small number.
                For economical production, use the defaults.)
        """
        self.text += "reconr\n"
        ngrid = len(reconr_grid)
        ncards = len(reconr_cards)
        self.text += "{} {}\n".format(nendf, npendf)              # card 1
        self.text += "'SANDY runs RECONR'/\n"                     # card 2
        for m in np.atleast_1d(mat):
            self.text += "{} {} {}\n".format(m, ncards, ngrid)  # card 3
            string = "{} {}".format(reconr_err, reconr_tempr)
            if reconr_errmax is None:
                string += " /\n"
            else:
                string += " {}".format(reconr_errmax)
                if reconr_errint is None:
                    string += " /\n"
                else:
                    string += " {}\n".format(reconr_errint)
            self.text += string                                   # card 4
            if ncards > 0:
                self.text += "/\n".format(reconr_cards) + "\n"           # card 5
            if ngrid > 0:
                self.text += "\n".join(map("{:.5e}".format, reconr_grid)) + "\n" # card 6
        self.text += "0 /\n"
    
    def stop(self):
        r"""
        Add ``NJOY`` stop command to instance's ``text``.
        """
        self.text += "stop"

    def write(self):
        r"""
        Write ``NJOY`` input file.
        This method is called in debug mode.
        """
        logging.debug("NJOY : write NJOY input to file '{}'".format(self.file))
        with open(self.file, 'w') as f:
            f.write(self.text)

    def run(self, cwd=None):
        """
        Run ``NJOY``.
        
        Pass input file to process's stdin (it must be encoded first).

        ..Important::
            In ``Python 3`` you need to convert string to bytes with a 
            ``encode()`` function
        
        Inputs:
            - ``cwd`` :
                working directory; default is current directory
        """
        from subprocess import Popen, PIPE
        stdin = self.text.encode()
        if logging.getLogger().getEffectiveLevel() <= logging.DEBUG:
            stdout = stderr = None
            self.write()
        else:
            stdout = stderr = PIPE
        process = Popen(exe, shell=True, cwd=cwd,
                        stdin=PIPE, 
                        stdout=stdout, 
                        stderr=stderr)
        stdoutdata, stderrdata = process.communicate(input=stdin)
        if process.returncode != 0:
            msg = "NJOY : Process status={}, cannot run njoy executable"
            logging.error(msg.format(process.returncode))
            sys.exit()
#        
#import pytest
#
#@pytest.fixture()
#def tape20(tmpdir):
#    from njoy import FileNJOY
#    import py
#    filein = py._path.local.LocalPath('../test_objects/file_H1-H2.endf')
#    fileout = tmpdir.join('tape20')
#    FileNJOY.copy_to_tape(str(filein), 20, dst=str(tmpdir))
#    assert filein.readlines() == fileout.readlines()
#    return fileout
#
#@pytest.fixture(scope="module")
#def fnjoy():
#    from njoy import FileNJOY
#    F = FileNJOY('input_njoy')
#    assert F.name == 'input_njoy'
#    return F
#
#
#def test_reconr_command_write_text(fnjoy):
#    grid = [1.5e-1, 1.6e-1]
#    fnjoy.reconr(20, 21, [125, 128], err=0.1, grid=grid)
#    text = fnjoy.text.split('\n')
#    assert text[0] == 'reconr'
#    assert list(map(int, text[1].split())) == [20, 21]
#    assert text[2] == "'SANDY runs RECONR'/"
#    assert list(map(int, text[3].split())) == [125, 0, 2]
#    assert list(map(float, text[4].split()[:-1])) == [0.1, 0]
#    assert float(text[5]) == grid[0]
#    assert float(text[6]) == grid[1]
#    assert list(map(int, text[7].split())) == [128, 0, 2]
#    assert list(map(float, text[8].split()[:-1])) == [0.1, 0]
#    assert float(text[9]) == grid[0]
#    assert float(text[10]) == grid[1]
#    assert text[11] == '0 /'
#
#def test_stop_command(fnjoy):
#    fnjoy.stop()
#    assert fnjoy.text.split('\n')[-1] == 'stop'
#
#def test_write_njoy(fnjoy, tmpdir):
#    fileout = tmpdir.join(fnjoy.name)
#    fnjoy.file = str(fileout)
#    fnjoy.write()
#    assert fnjoy.text == fileout.read()
#
#def test_run_njoy(fnjoy, tape20, tmpdir):
#    fnjoy.run("/home/lfiorito/work/njoy/source_2012/xnjoy-fast", 
#              cwd=str(tmpdir))
#    assert tmpdir.join('output').isfile()
#    assert tmpdir.join('tape21').isfile()
#
#def test_run_njoy_fail(fnjoy, tmpdir):
#    with pytest.raises(SystemExit):
#        fnjoy.run("/home/lfiorito/work/njoy/source_2012/xnjoy-fast", 
#                  cwd=str(tmpdir))