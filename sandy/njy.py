# -*- coding: utf-8 -*-
"""
Created on Wed Mar 28 15:06:23 2018

@author: fiorito_l
"""
import sys
import numpy as np
import os
import logging
from os.path import join

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

    def __init__(self, name='njoy.inp'):
        r"""
        Initialize `NJOY` input file.

        Inputs:
            - name :
                (string) `NJOY` inputfile's name
        """
        self.file = name
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

    def run(self, exe, cwd=None):
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
        process = Popen(exe,
                        shell=True,
                        cwd=cwd,
                        stdin=PIPE,
                        stdout=stdout,
                        stderr=stderr)
        stdoutdata, stderrdata = process.communicate(input=stdin)
        if process.returncode not in [0, 24]:
            msg = "NJOY : Process status={}, cannot run njoy executable"
            logging.error(msg.format(process.returncode))
            sys.exit()