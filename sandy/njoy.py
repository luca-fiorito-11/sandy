# -*- coding: utf-8 -*-
"""
Created on Fri May 25 16:58:11 2018
@author: fiorito_l
"""

__author__ = "Luca Fiorito"

import os
import sys
import shutil
import re
import logging
import pdb
import argparse
import tempfile
import subprocess as sp
import re


import pandas as pd
import numpy as np

from sandy.settings import SandyError
from sandy.utils import which


class EmptyIsConst(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        if len(values) == 0:
            string = option_string.lower()
            if string == "--temps":
                values = [293.6]
            elif string == "--kerma":
                values = [302, 303, 304, 318, 402, 442, 443, 444, 445, 446, 447]
            elif string == "--sig0":
                values = [1E10]
        setattr(namespace, self.dest, values)


def parser(iargs:"list of command line options as string"=None):
    """Parse command line options.
    """
    parser = argparse.ArgumentParser(prog="python -m sandy.njoy",
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('tape',
                        type=lambda x: utils.is_valid_file(parser, x),
                        help="ENDF-6 format file")
    parser.add_argument('-P','--pendftape',
                        type=lambda x: utils.is_valid_file(parser, x),
                        default=argparse.SUPPRESS,
                        metavar="PENDF",
                        help="processed PENDF format file")
    parser.add_argument('-G','--gendftape',
                        type=lambda x: utils.is_valid_file(parser, x),
                        default=argparse.SUPPRESS,
                        metavar="GENDF",
                        help="processed GENDF format file")
    parser.add_argument('--broadr',
                        type=utils.str2bool,
                        default=argparse.SUPPRESS,
                        metavar=True,
                        help=argparse.SUPPRESS)
    parser.add_argument('--purr',
                        type=utils.str2bool,
                        default=argparse.SUPPRESS,
                        metavar=True,
                        help=argparse.SUPPRESS)
    parser.add_argument('--unresr',
                        type=utils.str2bool,
                        default=argparse.SUPPRESS,
                        metavar=False,
                        help=argparse.SUPPRESS)
    parser.add_argument('--thermr',
                        type=utils.str2bool,
                        default=argparse.SUPPRESS,
                        metavar=True,
                        help=argparse.SUPPRESS)
    parser.add_argument('--heatr',
                        type=utils.str2bool,
                        default=argparse.SUPPRESS,
                        metavar=False,
                        help=argparse.SUPPRESS)
    parser.add_argument('--acer',
                        type=utils.str2bool,
                        default=argparse.SUPPRESS,
                        metavar=False,
                        help=argparse.SUPPRESS)
    parser.add_argument('-p','--pendf',
                        action="store_true",
                        help="produce PENDF file")
    parser.add_argument('-a','--ace',
                        action="store_true",
                        help="produce ACE file")
    parser.add_argument('-g','--gendf',
                        action="store_true",
                        help="produce GENDF file")
    parser.add_argument('-e','--errorr',
                        action="store_true",
                        help="produce ERRORR file")
    parser.add_argument('-H','--hendf',
                        action="store_true",
                        help="produce HENDF file")
#    parser.add_argument('--capsys',
#                        action="store_true",
#                        help="capture NJOY stderr and stdout (default=False)")
    parser.add_argument('--exe',
                        default=argparse.SUPPRESS,
                        help="NJOY executable (default=njoy2016).")
    parser.add_argument('--temps',
                        type=float,
                        default=argparse.SUPPRESS,
                        nargs='+',
                        help="list of temperature values")
    parser.add_argument('--sig0',
                        type=float,
                        default=argparse.SUPPRESS,
                        nargs='+',
                        help="list of dilution values")
    parser.add_argument('--kerma',
                        type=int,
                        action=EmptyIsConst,
                        default=argparse.SUPPRESS,
                        nargs='*',
                        help="list of partial kermas (default=[302, 303, 304, 318, 402, 442, 443, 444, 445, 446, 447])")
#    parser.add_argument('--qa',
#                        default=[],
#                        nargs='+',
#                        help="MT number and associated user Q-value (default=None).")
    parser.add_argument('--err',
                        type=float,
                        default=argparse.SUPPRESS,
                        help="fractional tolerance for RECONR and BROADR (default=0.001)")
    parser.add_argument('--free-gas',
                        default=argparse.SUPPRESS,
                        action="store_true",
                        help="compute thermal cross section for free-gas scatterer (THERMR)")
    parser.add_argument('--ptable',
                        default=argparse.SUPPRESS,
                        action="store_true",
                        help="compute probability tables (PURR)")
    parser.add_argument('--gaspr',
                        default=argparse.SUPPRESS,
                        action="store_true",
                        help="compute gas-production cross sections (GASPR)")
    parser.add_argument('--ign',
                        default=argparse.SUPPRESS,
                        help='''neutron group structure option for GROUPR
 - file : read from file
 - 2 : csewg 239-group structure [default]
 - 3 : lanl 30-group structure
 - 4 : anl 27-group structure
 - 5 : rrd 50-group structure
 - 6 : gam-i 68-group structure
 - 7 : gam-ii 100-group structure
 - 8 : laser-thermos 35-group structure
 - 9 : epri-cpm 69-group structure
 - 10 : lanl 187-group structure
 - 11 : lanl 70-group structure
 - 12 : sand-ii 620-group structure
 - 13 : lanl 80-group structure
 - 14 : eurlib 100-group structure
 - 15 : sand-iia 640-group structure
 - 16 : vitamin-e 174-group structure
 - 17 : vitamin-j 175-group structure
 - 19 : ecco 33-group structure
 - 20 : ecco 1968-group structure
 - 21 : tripoli 315-group structure
 - 22 : xmas lwpc 172-group structure
 - 23 : vit-j lwpc 175-group structure
predefined group structures:
 - scale_238''')
    parser.add_argument('--iwt',
                        default=argparse.SUPPRESS,
                        help='''Weight function option for GROUPR
 - file : read from file
 - 2 : constant
 - 3 : 1/e
 - 5 : epri-cell lwr
 - 6 : (thermal) -- (1/e) -- (fission + fusion)
 - 7 : same with t-dep thermal part
 - 8 : thermal--1/e--fast reactor--fission + fusion
 - 9 : claw weight function
 - 10 : claw with t-dependent thermal part
 - 11 : vitamin-e weight function (ornl-5505)
 - 12 : vit-e with t-dep thermal part
predefined functions:
 - jaea_fns_175
''')
    parser.add_argument('--igne',
                        default=argparse.SUPPRESS,
                        help="neutron group structure option for ERRORR (same as --ign)")
    parser.add_argument('--iwte',
                        default=argparse.SUPPRESS,
                        help="weight function option for ERRORR (same as --iwt)")
    parser.add_argument('--suffixes',
                        type=str,
                        default=argparse.SUPPRESS,
                        nargs='+',
                        metavar=".XX",
                        help="suffixes for ACE files, as many as temperature values (default = None).")
    parser.add_argument("-V","--verbose",
                        type=int,
                        default=argparse.SUPPRESS,
                        help="set verbosity level (default = 1)")
    parser.add_argument("-v",
                        '--version',
                        action='version',
                        version='%(prog)s 1.0',
                        help="")
    return parser.parse_args(iargs)

def mat_from_endf(tape):
    """Extract MAT number from file.
    Take it from section MF1/MT451.
    """
    line = open(tape).readlines()[1]
    return int(line[66:70])

def meta_from_endf(tape):
    """Extract metastate from file.
    Take it from section MF1/MT451.
    """
    line = open(tape).readlines()[2]
    return int(line[33:44])

sab = pd.DataFrame.from_records([[48,9237,1,1,241,'uuo2'],
                                  [42,125,0,8,221,'tol'],
                                  [59,1425,0,1,221,'si'],
                                  [37,125,11,2,223,'pol'],
                                  [2,125,0,2,221,'ph'],
                                  [12,128,0,2,221,'pd'],
                                  [75,825,1,1,239,'ouo2'],
                                  [48,825,0,3,221,'osap'],
                                  [51,825,0,1,222,'ohw'],
                                  [46,825,3,1,237,'obeo'],
                                  [3,125,0,2,221,'oh2'],
                                  [13,128,0,2,221,'od2'],
                                  [52,1225,0,1,249,'mg'],
                                  [38,125,0,12,221,'mesi'],
                                  [10,125,1,2,221,'ice'],
                                  [7,125,12,1,225,'hzr'],
                                  [1,125,0,2,222,'lw'],
                                  [8,125,0,2,249,'hca'],
                                  [31,600,1,1,229,'gr'],
                                  [11,128,0,2,221,'dhw'],
                                  [59,2025,0,1,249,'cah'],
                                  [27,425,3,1,233,'bbeo'],
                                  [26,425,2,1,231,'be'],
                                  [60,1325,0,2,221,'asap']],
columns = ['matde','matdp','icoh','natom','mtref','ext'])


def _set_njoyexe(self, njoyexe):
    """Mimic the behavior of the UNIX 'which' command to find the njoy executable.
    """
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)
    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file
    raise NotImplementedError("Not a valid NJOY-2016 executable: '{}'".format(_exe))

    
def _moder_input(nin, nout, **kwargs):
    text = ["moder"]
    text.append("{} {} /".format(nin, nout))
    return "\n".join(text) + "\n"


def _reconr_input(endfin, pendfout, mat,
                  header="sandy runs reconr",
                  err=0.001,
                  **kwargs):
    text = ["reconr"]
    text += ["{:d} {:d} /".format(endfin, pendfout)]
    text += ["'{}'/".format(header)]
    text += ["{:d} 0 0 /".format(mat)]
    text += ["{} 0. /".format(err)]
    text += ["0/"]
    return "\n".join(text) + "\n"

def _broadr_input(endfin, pendfin, pendfout, mat, temperatures=[293.6],
                  err=0.001,
                  **kwargs):
    text = ["broadr"]
    text += ["{:d} {:d} {:d} /".format(endfin, pendfin, pendfout)]
    text += ["{:d} {:d} 0 0 0. /".format(mat, len(temperatures))]
    text += ["{} /".format(err)]
    text += [" ".join(map("{:.1f}".format, temperatures)) + " /"]
    text += ["0 /"]
    return "\n".join(text) + "\n"

def _thermr_input(endfin, pendfin, pendfout, mat, temperatures=[293.6],
                  bins=20, iprint=0,
                  err=0.001, emax=10, **kwargs):
    text = ["thermr"]
    text += ["{:d} {:d} {:d} /".format(endfin, pendfin, pendfout)]
    text += ["0 {:d} {:d} {:d} 1 0 0 1 221 {:d} /".format(mat, bins, len(temperatures), iprint)]
    text += [" ".join(map("{:.1f}".format, temperatures)) + " /"]
    text += ["{} {} /".format(err, emax)]
    return "\n".join(text) + "\n"

def _purr_input(endfin, pendfin, pendfout, mat, temperatures=[293.6], sig0=[1e10], bins=20, ladders=32,
                **kwargs):
    text = ["purr"]
    text += ["{:d} {:d} {:d} /".format(endfin, pendfin, pendfout)]
    text += ["{:d} {:d} {:d} {:d} {:d} /".format(mat, len(temperatures), len(sig0), bins, ladders)]
    text += [" ".join(map("{:.1f}".format, temperatures)) + " /"]
    text += [" ".join(map("{:.2E}".format, sig0)) + " /"]
    text += ["0 /"]
    return "\n".join(text) + "\n"

def _gaspr_input(endfin, pendfin, pendfout,
                 **kwargs):
    text = ["gaspr"]
    text += ["{:d} {:d} {:d} /".format(endfin, pendfin, pendfout)]
    return "\n".join(text) + "\n"

def _unresr_input(endfin, pendfin, pendfout, mat, temperatures=[293.6], sig0=[1e10],
                  **kwargs):
    text = ["unresr"]
    text += ["{:d} {:d} {:d} /".format(endfin, pendfin, pendfout)]
    text += ["{:d} {:d} {:d} 1 /".format(mat, len(temperatures), len(sig0))]
    text += [" ".join(map("{:.1f}".format, temperatures)) + " /"]
    text += [" ".join(map("{:.2E}".format, sig0)) + " /"]
    text += ["0 /"]
    return "\n".join(text) + "\n"

def _heatr_input(endfin, pendfin, pendfout, mat, pks,
                 **kwargs):
    text = ["heatr"]
    text += ["{:d} {:d} {:d} 0 /".format(endfin, pendfin, pendfout)]
    text += ["{:d} {:d} 0 0 0 0 /".format(mat, len(pks))]
    text += [" ".join(map("{:d}".format, pks)) + " /"]
    return "\n".join(text) + "\n"

def _acer_input(endfin, pendfin, aceout, dirout, mat, temp=293.6,
                iopt=1, iprint=0, itype=1, suff=0,
                header="sandy runs acer",
                newfor=1, iopp=1, check=False,
                **kwargs):
    text = ["acer"]
    text += ["{:d} {:d} 0 {:d} {:d} /".format(endfin, pendfin, aceout, dirout)]
    text += ["{:d} {:d} {:d} .{:02d} 0 /".format(iopt, iprint, itype, suff)]
    text += ["'{}'/".format(header)]
    text += ["{:d} {:.1f} /".format(mat, temp)]
    text += ["{:d} {:d} /".format(newfor, iopp)]
    text += ["/"]
    if check:
        text += ["acer"]
        text += ["0 {} 0 0 0 /"]
        text += ["7 %(IPRINT)d 1 -1/"]
        text += ["/"]
    return "\n".join(text) + "\n"

def run_njoy(text, inputs, outputs, exe=None):
        """
        Run njoy executable.
        
        .. Important::

            In `Python 3` you need to convert input string to bytes with a 
            `encode()` function
        
        Parameters
        ----------
        exe : `str` or `None`
            njoy executable: if `None` (default) search in `PATH`
        inputs : `map`
            map of 
        text : `str`
            njoy input file passed to `Popen` as `stdin` (it must be encoded first)
        """
        if not exe:
            for try_exe in ["njoy2016", "njoy", "njoy2012", "xnjoy"]:
                exe = which(try_exe)
                if exe:
                    break
        if not exe:
            raise SandyError("could not find njoy executable")
        # if logging.getLogger().getEffectiveLevel() <= logging.DEBUG:
            # stdout = stderr = None
            # self.write()
        # else:
            # stdout = stderr = PIPE
        stdout = stderr = None
        stdin = text.encode()
        with tempfile.TemporaryDirectory() as tmpdir:
            logging.debug("Create temprary directory '{}'".format(tmpdir))
            for tape,src in inputs.items():
                shutil.copy(src, os.path.join(tmpdir, tape))
            process = sp.Popen(exe,
                               shell=True,
                               cwd=tmpdir,
                               stdin=sp.PIPE, 
                               stdout=stdout, 
                               stderr=stderr)
            stdoutdata, stderrdata = process.communicate(input=stdin)
            if process.returncode != 0:
                raise SandyError("process status={}, cannot run njoy executable".format(process.returncode))
            for tape,dst in inputs.items():
                path = os.path.split(dst)[0]
                if path:
                    os.makedirs(path, exist_ok=True)
                    
                pdb.set_trace()
                shutil.move(os.path.join(tmpdir, tape), dst)

def process(endftape, pendftape=None,
            kermas=[302, 303, 304, 318, 402, 442, 443, 444, 445, 446, 447],
            temperatures=[293.6],
            sig0=[1e10],
            broadr=True,
            thermr=True,
            unresr=False,
            heatr=True,
            gaspr=True,
            purr=True,
            acer=True,
            ext_pendf="pendf", ext_ace="ace", ext_xsd="xsd",
            wdir="", dryrun=False, tag=None, exe=None, keep_pendf=True, route="0",
            **kwargs):
    """Run sequence to process file with njoy.
    
    Parameters
    ----------
    
    Returns
    -------
    `str`
        njoy input text
    """
    inputs = {}
    outputs = {}
    kwargs.update({"temperatures" : temperatures})
    kwargs.update({"sig0" : sig0})
    if "mat" not in kwargs:
        kwargs.update({"mat" : mat_from_endf(endftape)})
    if "suffixes" not in kwargs:
        kwargs.update({"suffixes" : range(len(kwargs["temperatures"]))})
    else:
        if len(kwargs["suffixes"]) != len(kwargs["temperatures"]):
            raise SandyError('suffixes and temperatures must have the same size')
    inputs["tape20"] = endftape
    e = 21
    p = e + 1
    text = _moder_input(20, -e)
    if pendftape:
        inputs["tape99"] = pendftape
        text += _moder_input(99, -p)
    else:
        text += _reconr_input(-e, -p, **kwargs)
    if broadr:
        o = p + 1
        text += _broadr_input(-e, -p, -o, **kwargs)
        p = o
    if thermr:
        o = p + 1 
        text += _thermr_input(0, -p, -o, **kwargs)
        p = o
    if unresr:
        o = p + 1
        text += _unresr_input(-e, -p, -o, **kwargs)
        p = o
    if heatr:
        for i in range(0, len(kermas), 7):
            o = p + 1
            kwargs["pks"] = kermas[i:i+7]
            text += _heatr_input(-e, -p, -o, **kwargs)
            p = o
    if gaspr:
        o = p + 1
        text += _gaspr_input(-e, -p, -o, **kwargs)
        p = o
    if purr:
        o = p + 1
        text += _purr_input(-e, -p, -o, **kwargs)
        p = o
    if keep_pendf:
        o = p + 1
        text += _moder_input(-p, o)
        outputs["tape{}".format(o)] = os.path.join(wdir, "{}.{}".format(tag, ext_pendf))
    if acer:
        for i,(temp,suff) in enumerate(zip(kwargs["temperatures"], kwargs["suffixes"])):
            a = 50 + i
            x = 70 + i
            kwargs["temp"] = temp
            kwargs["suff"] = suff
            text += _acer_input(-e, -p, a, x, **kwargs)
            outputs["tape{}".format(a)] = os.path.join(wdir, "{}.{}".format(tag, ext_ace))
            outputs["tape{}".format(x)] = os.path.join(wdir, "{}.{}".format(tag, ext_xsd))
    text += "stop"
    if dryrun:
        return text
    run_njoy(text, inputs, outputs, exe=exe)
    # If isotope is metatable rewrite ID in xsdir and ace as ID = Z*1000 + 300 + A + META*100.
    # Also change route and filename in xsdir file.
    return text
    if acer:
        meta = meta_from_endf(endftape)
        for i,(temp,suff) in enumerate(zip(kwargs["temperatures"], kwargs["suffixes"])):
            a = 50 + i
            x = 70 + i
            acefile = outputs["tape{}".format(a)]
            xsdfile = outputs["tape{}".format(x)]
            text_xsd = open(xsdfile).read()
            text_xsd = re.sub("route", route, text_xsd)
            text_xsd = re.sub("filename", acefile, text_xsd)
            text_xsd = " ".join(text_xsd.split())
            if meta:
                pattern = '(?P<za>\d{4,6})\.(?P<ext>\d{2}[ct])'
                found = re.search(pattern, text_xsd)
                za = int(found.group("za"))
                ext = found.group("ext")
                za_new = za + meta*100 + 300
                text_xsd = text_xsd.replace("{}.{}".format(za, ext), "{}.{}".format(za_new, ext), 1)
                text_ace = open(acefile).read()
                text_ace = text_ace.replace("{}.{}".format(za, ext), "{}.{}".format(za_new, ext), 1)
                with open(acefile, 'w') as f:
                    f.write(text_ace)
            with open(xsdfile, 'w') as f:
                f.write(text_xsd)
 


 
class Njoy:
    """ 
    """

    def __init__(self,
                 iopt=1,
                 err=0.001,
                 iin=1,
                 icoh=1,
                 iform=1,
                 natom=1,
                 # ACER
                 itype=1,
                 
                 **kwargs):
#        self.sab = os.path.join(os.path.dirname(os.path.abspath(__file__)), "sab.csv")
#        self.wdir = os.path.normpath(wdir)
#        self.evaluationName = os.path.basename(self.wdir)
        self.processes = None
        self.capsys = False
        self.iwt = 4
        self.legendre = 1
        self.legendregg = 6
        self.scatteringLaw = None
        self.eFiss = None
        self.branchingNG = None
        self.branchingN2N = None
        self.gstr = 0
        self.oldlib = None
        self.sgref = 1.0E10
        self.yields = None
        self.iburn = -1
        self.ip1opt = 0
        self.ifprod = 0
        self.jp1 = 0
        self.iverw = 4
        # General options
        self.WDIR = os.getcwd()
        self.folder = os.path.abspath("lib")
        self.IPRINT = 1  # by default set verbosity to max
        self.PLOT = False
        # modules
        self.gaspr = False
        # RECONR/BROADR default options
        self.err = 0.001
        self.THNMAX = 1E6
        # PURR/UNRESR default options
        self.ptable = False
        # THERMR default options
        self.free_gas = False
        self.IIN = 1
        self.ICOH = 1
        self.IFORM = 0
        self.NATOM = 1
        self.MTREF = 221
        # ACER default options
        self.IOPT = 1
        self.ITYPE = 1
        self.NEWFOR = 1
        self.IOPP = 1
        self.HK = ''
        # GROUPR default options
        self.IGN = 2 # neutron group structure for GROUPR
        self.IGG = 0
        self.IWT = 6
        self.LORD = 1
        self.NSTR = None
        self.GSTR = None
        # ERRORR default options
        self.IGNE = 2 # neutron group structure for ERRORR
        self.IWTE = 6
        self.IRELCO = 1
        self.IRESPR = 1
        # HEATR default options
        self.QA = []
        self.LOCAL = 0
        for k,v in kwargs.items(): setattr(self, k, v)
        if not hasattr(self, "exe"): self.exe = "njoy2016"

    @property
    def errmax(self):
        if hasattr(self, "_errmax"):
            return self._errmax
        else:
            return self.err*10.
    
    @errmax.setter
    def errmax(self, errmax):
        self._errmax = errmax
    
    @property
    def exe(self):
        return self._exe
    
    @exe.setter
    def exe(self, exe):
        _exe = utils.which(exe)
        if _exe is None:
            raise NotImplementedError("Not a valid NJOY-2016 executable: '{}'".format(_exe))
        self._exe = _exe

    @property
    def wdir(self):
        if hasattr(self, "_wdir"):
            return self._wdir
        else:
            return os.getcwd()
    
    @wdir.setter
    def wdir(self, wdir):
        self._wdir = wdir

    @property
    def temps(self):
        if hasattr(self, "_temps"):
            return self._temps
        else:
            return [0]

    @temps.setter
    def temps(self, temps):
        """
        Ensures that the temperatures are given correctly.
        """
        if not isinstance(temps, list):
            raise TypeError("{} must be a list".format("temps"))
        if len(temps)> 10:
            raise ValueError("cannot have more than 10 temperatures")
        self._temps = temps

    @property
    def sig0(self):
        if hasattr(self, "_sig0"):
            return self._sig0
        else:
            return [1e10]

    @sig0.setter
    def sig0(self, sig0):
        """
        Ensures that the dilutions are given correctly.
        """
        if not isinstance(sig0, list):
            raise TypeError("{} must be a list".format("sig0"))
        self._sig0 = sig0

    @property
    def suffixes(self):
        if hasattr(self, "_suffixes"):
            suffixes = self._suffixes
        else:
            suffixes = [".00"]
        if len(suffixes) != len(self.temps):
            raise ValueError("'suffixes' must have as many items as 'temps'")
        return suffixes

    @suffixes.setter
    def suffixes(self, suffixes):
        """
        Ensures that the suffixes are as many as the temperature values.
        """
        if not isinstance(suffixes, list):
            raise TypeError("{} must be a list".format("suffixes"))
        self._suffixes = suffixes

    @property
    def kerma(self):
        if hasattr(self, "_kerma"):
            return self._kerma
        else:
            return [302, 303, 304, 318, 402, 442, 443, 444, 445, 446, 447]

    @kerma.setter
    def kerma(self, kerma):
        """
        Ensures that the partial kerma are given correctly.
        """
        if not isinstance(kerma, list):
            raise TypeError("{} must be a list".format("temps"))
        self._temps = kerma

    @staticmethod
    def _moder_input(tapein, tapeout, **kwargs):
        text = ["moder"]
        text.append("{} {} /".format(tapein, tapeout))
        return "\n".join(text) + "\n"

    @staticmethod
    def _reconr_input(endfin, pendfout, mat, err=0.001, errmax=0.01, **kwargs):
        text = ["reconr"]
        text.append("{} {} /".format(endfin, pendfout))
        text.append("' '/")
        text.append("{} 0 0 /".format(mat))
        text.append("{} 0. {} /".format(err, errmax))
        text.append("0/")
        return "\n".join(text) + "\n"

    @staticmethod
    def _broadr_input(endfin, pendfin, pendfout, mat, temps, err=0.001, **kwargs):
        kwargs.update({"ENDFIN" : endfin, "PENDFIN" : pendfin, "PENDFOUT" : pendfout})
        text = ["broadr"]
        text += ["{} {} {} /".format(endfin, pendfin, pendfout)]
        text += ["{} {} 0 0 0. /".format(mat, len(temps))]
        text += ["{} /".format(err)]
        text += [" ".join(["{}".format(tmp) for tmp in temps])]
        text += ["0 /"]
        return "\n".join(text) + "\n"

    @staticmethod
    def _thermr_input(endfin, pendfin, pendfout, mat, temps, iin=1, ico=1,
                      iform=0, natom=1, iprint=1, **kwargs):
        text = ["thermr"]
        text += ["{} {} {} /".format(endfin, pendfin, pendfout)]
        text += ["0 {} 20 {} {} {} {} {} 221 {}".format(mat, len(temps), iin, ico, iform, natom, iprint)]
        text += [" ".join(["{}".format(tmp) for tmp in temps])]
        text += ["0.001 4.0/"]
        return "\n".join(text) + "\n"

    @staticmethod
    def _purr_input(endfin, pendfin, pendfout, mat, temps, sig0, **kwargs):
        text = ["purr"]
        text += ["{} {} {} /".format(endfin, pendfin, pendfout)]
        text += ["{} {} {} 20 32 /".format(mat, len(temps), len(sig0))]
        text += [" ".join(map(str, temps))]
        text += [" ".join(map("{:E}".format, sig0))]
        text += ["0 /"]
        return "\n".join(text) + "\n"

    @staticmethod
    def _gaspr_input(endfin, pendfin, pendfout, **kwargs):
        text = ["purr"]
        text += ["{} {} {} /".format(endfin, pendfin, pendfout)]
        return "\n".join(text) + "\n"

    @staticmethod
    def _unresr_input(endfin, pendfin, pendfout, mat, temps, sig0, **kwargs):
        text = ["unresr"]
        text += ["{} {} {} /".format(endfin, pendfin, pendfout)]
        text += ["{} {} {} 1 /".format(mat, len(temps), len(sig0))]
        text += [" ".join(map(str, temps))]
        text += [" ".join(map("{:E}".format, sig0))]
        text += ["0 /"]
        return "\n".join(text) + "\n"

    @staticmethod
    def _acer_input(nendf, npendf, nace, ndir, mat, temp, suff, iopt=1,
                    iprint=1, itype=1, descr="", newfor=1, iopp=1, **kwargs):
        text = ["acer"]
        text += ["{} {} 0 {} {} /".format(nendf, npendf, nace, ndir)]
        text += ["{} {} {} {} 0 /".format(iopt, iprint, itype, suff)]
        text += [descr + "/"]
        text += ["{} {} /".format(mat, temp)]
        text += ["{} {} /".format(newfor, iopp)]
        text += ["0.001 /"]
#        if kwargs["VERBOSE"] >= 2:
#            text += """acer / Check ACE files
#0 %(NACE) 0 0 0 /
#7 %(IPRINT)d 1 -1/
#/
#"""
        return "\n".join(text) + "\n"

    @staticmethod
    def _heatr_input(endfin, pendfin, pendfout, **kwargs):
        kwargs = dict(dict(vars(self), TEMPS=self.TEMPS), **fileOptions)
        text = ""
        npks = len(kwargs["KERMA"])
        size = 7
        sign = np.sign(pendfout)
        kwargs.update({"ENDFIN" : endfin, "PENDFIN" : pendfin})
        if "QA" not in kwargs: kwargs["QA"] = []
        kwargs.update({"NQA" : len(kwargs["QA"])//2})
        kwargs["TEXTMTA"] = " ".join(map(lambda x: str(int(x)), kwargs["QA"][::2]))
        kwargs["TEXTQA"] = " ".join(map(str, kwargs["QA"][1::2]))
        out = 79
        for i in range(0, npks, size):
            kwargs["PENDFOUT"] = pendfout if i+size >= npks else sign*(out + 1)
            kwargs["NPK"] = len(kwargs["KERMA"][i:i+size])
            kwargs["TEXTPKS"] = " ".join(map(str,kwargs["KERMA"][i:i+size]))
            text += """heatr
%(ENDFIN)d %(PENDFIN)d %(PENDFOUT)d 0 /
%(MAT)d %(NPK)d %(NQA)d 0 %(LOCAL)d %(IPRINT)d /
%(TEXTPKS)s /
""" % kwargs
            if kwargs["NQA"] > 0:
                text += """%(TEXTMTA)s /
%(TEXTQA)s /
""" % kwargs
            kwargs.update({"PENDFIN" : kwargs["PENDFOUT"]})
        return text

    def _groupr_input(self, nendf, npendf, ngin, ngout, **fileOptions):
        from ..formats.records import write_tab1
        kwargs = dict(dict(vars(self), TEMPS=self.TEMPS, SIG0=self.SIG0), **fileOptions)
        kwargs.update({"NENDF" : nendf, "NPENDF" : npendf, "NGIN" : ngin, "NGOUT" : ngout})
        kwargs["NSIG0"] = len(kwargs["SIG0"])
        kwargs["TEXTSIG0"] = " ".join(["%E"%dil for dil in kwargs["SIG0"]])
        kwargs["NTEMPS"] = len(kwargs["TEMPS"])
        kwargs["TEXTTEMPS"] = " ".join(["%E"%tmp for tmp in kwargs["TEMPS"]])
        try:
            ign = int(kwargs["IGN"])
        except ValueError as exc:
            ign = kwargs["IGN"]
        finally:
            kwargs["IGN"] = ign
        if not isinstance(kwargs["IGN"], int):
            if kwargs["IGN"].lower() == "scale_238":
                from ..energy_grids import scale_238 as grid
            else:
                if not os.path.isfile(kwargs["IGN"]):
                    logging.error("file {} not found".format(kwargs["IGN"]))
                    sys.exit()
                grid = np.genfromtxt(kwargs["IGN"])
            ngroups = len(grid) - 1
            kwargs["IGNSTR"] = "{} /\n".format(ngroups) + "\n".join(map(str,grid)) + " /\n"
            kwargs["IGN"] = 1
        try:
            iwt = int(kwargs["IWT"])
        except ValueError as exc:
            iwt = kwargs["IWT"]
        finally:
            kwargs["IWT"] = iwt
        if not isinstance(kwargs["IWT"], int):
            if kwargs["IWT"].lower() == "jaea_fns_175":
                from ..spectra import jaea_fns_175 as spectrum
            else:
                if not os.path.isfile(kwargs["IWT"]):
                    logging.error("file {} not found".format(kwargs["IWT"]))
                    sys.exit()
                spectrum = pd.read_csv(kwargs["IWT"], header=None, names=["E", "F"])
            tab = spectrum.reset_index().values.flatten()
            x = tab[::2]; y = tab[1::2]
            kwargs["IWTSTR"] = "\n".join(write_tab1(0, 0, 0, 0, [len(x)], [2], x, y)) + "/\n"
            kwargs["IWT"] = 1
        text = """groupr
%(NENDF)d %(NPENDF)d %(NGIN)d %(NGOUT)d /
%(MAT)d %(IGN)d %(IGG)d %(IWT)d %(LORD)d %(NTEMPS)d %(NSIG0)d %(IPRINT)d /
'' /
%(TEXTTEMPS)s /
%(TEXTSIG0)s /
""" % kwargs
        if "IGNSTR" in kwargs:
            text += "%(IGNSTR)s" % kwargs
        if "IWTSTR" in kwargs:
            text += "%(IWTSTR)s" % kwargs # card13
        for tmp in kwargs["TEMPS"]:
            text += """3/
6/
0/
"""
        text += "0/\n"
        return text

    def _errorr_input(self, endfin, pendfin, gendfin, errorrout, **fileOptions):
        from ..formats.records import write_tab1
        kwargs = dict(dict(vars(self), TEMPS=self.TEMPS), **fileOptions)
        kwargs.update({"ENDFIN" : endfin, "PENDFIN" : pendfin, "GENDFIN" : gendfin, "ERRORROUT" : errorrout})
        try:
            igne = int(kwargs["IGNE"])
        except ValueError as exc:
            igne = kwargs["IGNE"]
        finally:
            kwargs["IGNE"] = igne
        if not isinstance(kwargs["IGNE"], int):
            if kwargs["IGNE"].lower() == "scale_238":
                from ..energy_grids import scale_238 as grid
            else:
                if not os.path.isfile(kwargs["IGNE"]):
                    logging.error("file {} not found".format(kwargs["IGNE"]))
                    sys.exit()
                grid = np.genfromtxt(kwargs["IGNE"])
            ngroups = len(grid) - 1
            kwargs["IGNESTR"] = "{} /\n".format(ngroups) + "\n".join(map(str,grid)) + " /\n"
            kwargs["IGNE"] = 1
        try:
            iwte = int(kwargs["IWTE"])
        except ValueError as exc:
            iwte = kwargs["IWTE"]
        finally:
            kwargs["IWTE"] = iwte
        if not isinstance(kwargs["IWTE"], int):
            if kwargs["IWTE"].lower() == "jaea_fns_175":
                from ..csv.spectra import jaea_fns_175 as spectrum
            else:
                if not os.path.isfile(kwargs["IWTE"]):
                    logging.error("file {} not found".format(kwargs["IWTE"]))
                    sys.exit()
                spectrum = pd.read_csv(kwargs["IWTE"], header=None, names=["E", "F"])
            tab = spectrum.reset_index().values.flatten()
            x = tab[::2]; y = tab[1::2]
            kwargs["IWTESTR"] = "\n".join(write_tab1(0, 0, 0, 0, [len(x)], [2], x, y)) + "/\n"
            kwargs["IWTE"] = 1
        text = ""
        if 32 in kwargs["SECTIONS"] and 33 not in kwargs["SECTIONS"]:
            mts = [1, 2, 18, 102, 0] if kwargs["LFI"] == 1 else [1, 2, 102, 0]
            kwargs["MTS"] = " /\n".join(map(str,mts)) + " /"
            text += """moder
%(ENDFIN)d 60
errorr
999 / option to insert dummy file 33 data
60 61 /
%(MTS)s
""" % kwargs
            if kwargs["ENDFIN"] < 0:
                text += """moder
61 -62
"""
                kwargs["ENDFIN"] = 61
            else:
                kwargs["ENDFIN"] = -62
        text += """errorr
%(ENDFIN)d %(PENDFIN)d %(GENDFIN)d %(ERRORROUT)d /
%(MAT)d %(IGNE)d %(IWTE)d %(IPRINT)d %(IRELCO)d /
%(IPRINT)d %(TMP)E /
0 33 %(IRESPR)d %(LORD)d /
""" % kwargs
        if "IGNESTR" in kwargs:
            text += "%(IGNESTR)s" % kwargs # card12 + card12b
        if "IWTESTR" in kwargs:
            text += "%(IWTESTR)s" % kwargs # card13
        return text


    def get_pendf(self, tape, mat=None, pendftape=None, tag=None, 
                  **fileOptions):
        """ Run NJOY sequence to produce a PENDF file.
        """
        if mat is None: mat = mat_from_endf(tape)
        fname = os.path.split(tape)[1]
        if tag is None: tag = fname
        mydir = os.path.join(self.wdir, fname)
        os.makedirs(mydir, exist_ok=True)
        DfOutputs = OutputFiles()
        # Build njoy input text
        utils.force_symlink(tape, os.path.join(mydir, "tape20"))
        text = self._moder_input(20, -21)
        if pendftape is None:
            text += self._reconr_input(-21, -22, mat, self.err)
        else:
            utils.force_symlink(pendftape, os.path.join(mydir, "tape30"))
            text += self._moder_input(30, -22)
        e = -21; p = -22
        if hasattr(self, "_temps"):
            text += self._broadr_input(e, p, -23, mat, self.temps, self.err)
            p = -23
        if self.free_gas:
            text += self._thermr_input(0, p, -24, mat, self.temps)
            p = -24
        if hasattr(self, "_sig0"):
            text += self._unresr_input(e, p, -25, mat, self.temps, self.sig0)
            p = -25
        if self.ptable:
            text += self._purr_input(e, p, -26, mat, self.temps, self.sig0)
            p = -26
        if hasattr(self, "_kerma"):
            text += self._heatr_input(e, p, -27, **fileOptions)
            p = -27
        if self.gaspr:
            text += self._gaspr_input(e, p, -28)
            p = -28
        text += self._moder_input(p, 29)
        text += "stop"
        # Run NJOY
        inputfile = os.path.join(mydir, "input_pendf.{}".format(tag))
        DfOutputs = DfOutputs.append({"id" : "input_pendf", "format" : "TEXT", "file" : inputfile}, ignore_index=True)
        with open(inputfile,'w') as f: f.write(text)
        returncode, stdout, stderr = utils.run_process(
                "{} < {}".format(self.exe, inputfile),
                cwd=mydir,
                verbose=True,
                )
        # Process NJOY outputs
        oldfile = os.path.join(mydir, "tape29")
        newfile = os.path.join(mydir, "{}.pendf".format(tag))
        if (os.path.isfile(oldfile) and os.path.getsize(oldfile) > 0):
            shutil.move(oldfile, newfile)
            DfOutputs = DfOutputs.append({"id" : "pendf", "format" : "ENDF", "file" : newfile}, ignore_index=True)
        oldfile = os.path.join(mydir, "output")
        newfile = os.path.join(mydir, "output_pendf.{}".format(tag))
        if (os.path.isfile(oldfile) and os.path.getsize(oldfile) > 0):
            shutil.move(oldfile, newfile)
            DfOutputs = DfOutputs.append({"id" : "output_pendf", "format" : "TEXT", "file" : newfile}, ignore_index=True)
            DfMessages = NjoyOutput.from_file(newfile).get_messages()
        for filename in os.listdir(mydir):
            if filename[:4] == 'tape': os.remove(os.path.join(mydir, filename))
        return DfOutputs, DfMessages

    def get_gendf(self, **fileOptions):
        mydir = os.path.join(self.folder, fileOptions["FILENAME"])
        os.makedirs(mydir, exist_ok=True)

        utils.force_symlink(fileOptions["TAPE"], os.path.join(mydir, "tape20"))
        text = self._moder_input(20, -21, **fileOptions)

        utils.force_symlink(fileOptions["PENDFTAPE"], os.path.join(mydir, "tape29"))
        text += self._moder_input(29, -25, **fileOptions)

        text += self._groupr_input(-21, -25, 0, -26, **fileOptions)
        text += self._moder_input(-26, 30, **fileOptions)
        text += "stop"

        DfOutputs = OutputFiles()
        inputfile = os.path.join(mydir, "input_gendf.{}".format(fileOptions["TAG"]))
        DfOutputs = DfOutputs.append({"id" : "input_gendf", "format" : "TEXT", "file" : inputfile}, ignore_index=True)
        with open(inputfile,'w') as f: f.write(text)
        print(" --- run gendf for {} ---".format(fileOptions["TAG"]))
        returncode, stdout, stderr = run_process("{} < {}".format(self.exe, inputfile), cwd=mydir)

        oldfile = os.path.join(mydir, "tape30")
        newfile = os.path.join(mydir, "{}.gendf".format(fileOptions["TAG"]))
        if os.path.isfile(oldfile):
            if os.path.getsize(oldfile) > 0:
                shutil.move(oldfile, newfile)
                DfOutputs = DfOutputs.append({"id" : "gendf", "format" : "GENDF", "file" : newfile}, ignore_index=True)

        oldfile = os.path.join(mydir, "output")
        newfile = os.path.join(mydir, "output_gendf.{}".format(fileOptions["TAG"]))
        shutil.move(oldfile, newfile)
        DfOutputs = DfOutputs.append({"id" : "output_gendf", "format" : "TEXT", "file" : newfile}, ignore_index=True)
        DfMessages = NjoyOutput.from_file(newfile).get_messages()

        for filename in os.listdir(mydir):
            if filename[:4] == 'tape': os.remove(os.path.join(mydir, filename))
        return DfOutputs, DfMessages


    def get_errorr(self, **fileOptions):
        mydir = os.path.join(self.folder, fileOptions["FILENAME"])
        os.makedirs(mydir, exist_ok=True)

        utils.force_symlink(fileOptions["TAPE"], os.path.join(mydir, "tape20"))
        text = self._moder_input(20, -21, **fileOptions)

#        if "GENDFTAPE" in fileOptions:
#            force_symlink(fileOptions["GENDFTAPE"], os.path.join(mydir, "tape29"))
#            text += self.moder_input(29, 25, **fileOptions)
#            p = 0; g = 25
#        elif "PENDFTAPE" in kwargs:
        utils.force_symlink(fileOptions["PENDFTAPE"], os.path.join(mydir, "tape29"))
        text += self._moder_input(29, -25, **fileOptions)
        p = -25; g = 0

        errorrs = []; er = 30
        for tmp in self.TEMPS:
            text += self._errorr_input(-21, p, g, er, TMP=tmp, **fileOptions)
            errorrs.append(er)
            er += 1
        text += "stop"

        DfOutputs = OutputFiles()
        inputfile = os.path.join(mydir, "input_errorr.{}".format(fileOptions["TAG"]))
        DfOutputs = DfOutputs.append({"id" : "input_errorr", "format" : "TEXT", "file" : inputfile}, ignore_index=True)
        with open(inputfile,'w') as f: f.write(text)
        print(" --- run errorr for {} ---".format(fileOptions["TAG"]))
        returncode, stdout, stderr = run_process("{} < {}".format(self.exe, inputfile), cwd=mydir)

        newfile = os.path.join(mydir, "{}.errorr".format(fileOptions["TAG"]))
        if os.path.isfile(newfile): os.remove(newfile)
        with open(newfile, "a") as f:
            for er in errorrs:
                oldfile = os.path.join(mydir, "tape{}".format(er))
                if not os.path.isfile(oldfile): continue
                with open(oldfile) as g:
                    lines = g.readlines()
                    if er != errorrs[0]:
                        lines = lines[1:]
                    if er != errorrs[-1]:
                        lines = lines[:-1]
                    f.writelines(lines)
        if os.path.getsize(newfile) == 0: os.remove(newfile)
        if os.path.isfile(newfile): DfOutputs = DfOutputs.append({"id" : "errorr", "format" : "ERRORR", "file" : newfile}, ignore_index=True)

        oldfile = os.path.join(mydir, "output")
        newfile = os.path.join(mydir, "output_errorr.{}".format(fileOptions["TAG"]))
        shutil.move(oldfile, newfile)
        DfOutputs = DfOutputs.append({"id" : "output_errorr", "format" : "TEXT", "file" : newfile}, ignore_index=True)
        DfMessages = NjoyOutput.from_file(newfile).get_messages()

        for filename in os.listdir(mydir):
            if filename[:4] == 'tape': os.remove(os.path.join(mydir, filename))
        return DfOutputs, DfMessages


    def get_ace(self, tape, mat=None, pendftape=None, tag=None, 
                **fileOptions):
        """ Run the NJOY sequence to produce a ACE file.
        """
        if mat is None: mat = mat_from_endf(tape)
        fname = os.path.split(tape)[1]
        if tag is None: tag = fname
        mydir = os.path.join(self.wdir, fname)
        os.makedirs(mydir, exist_ok=True)
        if pendftape is None:
            dfout, dferr = self.get_pendf(tape, mat=mat, pendftape=pendftape, tag=tag)
            pendftape = dfout.query("id==\"pendf\"").file.iloc[0]
        else:
            dfout = OutputFiles()
            dferr = NjoyMessages()
        # Build njoy input text
        utils.force_symlink(tape, os.path.join(mydir, "tape20"))
        text = self._moder_input(20, -21)
        utils.force_symlink(pendftape, os.path.join(mydir, "tape29"))
        text += self._moder_input(29, -25)
        aces = []; xsdirs = []
        e = -21; p = -25; a = 40; x = 50
        for tmp,suff in zip(self.temps, self.suffixes):
            text += self._acer_input(e, p, a, x, mat, tmp, suff)
            aces.append(a); xsdirs.append(x)
            a += 1; x += 1
        text += "stop"
        # Run NJOY and process outputs
        inputfile = os.path.join(mydir, "input_acer.{}".format(tag))
        dfout = OutputFiles()
        dfout = dfout.append({"id" : "input_acer", "format" : "ACE", "file" : inputfile}, ignore_index=True)
        with open(inputfile,'w') as f: f.write(text)
        returncode, stdout, stderr = utils.run_process(
                "{} < {}".format(self.exe, inputfile),
                cwd=mydir,
                verbose=True,
                )
        # Process NJOY outputs
        newfile = os.path.join(mydir, "{}.ace".format(tag))
        if os.path.isfile(newfile): os.remove(newfile)
        for a in aces:
            oldfile = os.path.join(mydir, "tape{}".format(a))
            if not os.path.isfile(oldfile): continue
            with open(newfile, "a") as f, open(oldfile) as g: f.write(g.read())
        if (os.path.isfile(newfile) and os.path.getsize(newfile) > 0):
                dfout = dfout.append({"id" : "ace", "format" : "ACE", "file" : newfile}, ignore_index=True)
        newfile = os.path.join(mydir, "{}.xsdir".format(tag))
        if os.path.isfile(newfile): os.remove(newfile)
        with open(newfile, "a") as f:
            lines = []
            for x in xsdirs:
                oldfile = os.path.join(mydir, "tape{}".format(x))
                if not os.path.isfile(oldfile): continue
                with open(oldfile) as g:
                    xargs = g.read().split()
                    xargs[2] = "{}.ace".format(tag)
                    xargs[3] = "0"
                    lines.append(" ".join(xargs))
            f.write("\n".join(lines))
        if (os.path.isfile(newfile) and os.path.getsize(newfile) > 0):
            dfout = dfout.append({"id" : "xsdir", "format" : "TEXT", "file" : newfile}, ignore_index=True)
        oldfile = os.path.join(mydir, "output")
        newfile = os.path.join(mydir, "output_acer.{}".format(tag))
        if (os.path.isfile(oldfile) and os.path.getsize(oldfile) > 0):
            shutil.move(oldfile, newfile)
            dfout = dfout.append({"id" : "output_acer", "format" : "TEXT", "file" : newfile}, ignore_index=True)
            dferr = NjoyOutput.from_file(newfile).get_messages()
        for filename in os.listdir(mydir):
            if filename[:4] == 'tape': os.remove(os.path.join(mydir, filename))
        return dfout, dferr

def run(iargs=None):
    init = parser(iargs)
    nj = Njoy(**vars(init))

    DfOutputs = OutputFiles()
    DfMessages = NjoyMessages()

    options = {"tape" : init.tape}
    if init.pendftape: options.update({"pendftape" : init.pendftape})
    if init.gendftape: options.update({"gendftape" : init.pendftape})
    if init.tag: options.update({"tag" : init.tag})
    if init.mat: options.update({"mat" : init.mat})
    if init.pendf:
        outs, msgs = nj.get_pendf(**options)
        DfOutputs = pd.concat([DfOutputs, outs]).reset_index(drop=True)
        DfMessages = pd.concat([DfMessages, msgs]).reset_index(drop=True)
        options.update({"pendftape" : DfOutputs.query("id==\"pendf\"").file.iloc[0]})

    if init.ace:
        outs, msgs = nj.get_ace(**options)
        DfOutputs = pd.concat([DfOutputs, outs]).reset_index(drop=True)
        DfMessages = pd.concat([DfMessages, msgs]).reset_index(drop=True)

#    if init.gendf:
#        outs, msgs = nj.get_gendf(**tape.__dict__)
#        DfOutputs = pd.concat([DfOutputs, outs]).reset_index(drop=True)
#        DfMessages = pd.concat([DfMessages, msgs]).reset_index(drop=True)
#    try:
#        tape.GENDFTAPE = DfOutputs.query("id==\"gendf\"").file.iloc[0]
#    except:
#        pass
#
#    if init.errorr:
#        outs, msgs = nj.get_errorr(**tape.__dict__)
#        DfOutputs = pd.concat([DfOutputs, outs]).reset_index(drop=True)
#        DfMessages = pd.concat([DfMessages, msgs]).reset_index(drop=True)
