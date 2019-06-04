# -*- coding: utf-8 -*-
"""
Outline
=======
1. Summary_
2. Examples_
3. Routines_

.. _Summary:

Summary
=======
This module contains template inputs for NJOY routines and functions to run them.

Two major functions `process` and `process_protons` are provided to process nuclear data 
files with NJOY into ACE format, respectively for fast neutron-induced and proton-induced 
nuclear data.

Given any nuclear data evaluation file for incident neutrons (fast, not SAB) function `process` 
generates the correspoding ACE filea for a given set of temperatures (one file per temperature).
If no keyword argument is provided, function `process` runs with default options, which include 
NJOY routines RECONR, BROADR, THERMR, HEATR, GASPR, PURR, ACER.
Keyword arguments can be changed to add/remove NJOY routines using `True/False` flags, or to change 
a routine's input parameters.

Major default parmameters:

+------------------+-----------------------------------------------------------+------------------------------+ 
| Parameter        | Value                                                     | Description                  |
+==================+===========================================================+==============================+
| err              | `0.001`                                                   | xs reconstruction tolerance  |
+------------------+-----------------------------------------------------------+------------------------------+ 
| temperatures     | `[293.6]`                                                 | `list` of temperatures (K)   |
+------------------+-----------------------------------------------------------+------------------------------+ 
| bins             | `20`                                                      | # probability bins (PURR)    |
+------------------+-----------------------------------------------------------+------------------------------+ 
| ladders          | `32`                                                      | # resonance ladders (PURR)   |
+------------------+-----------------------------------------------------------+------------------------------+ 
| iprint           | `False`                                                   | output verbosity             |
+------------------+-----------------------------------------------------------+------------------------------+ 
| kermas           | `[302, 303, 304, 318, 402, 442, 443, 444, 445, 446, 447]` | `list` of KERMA factors (MT) |
+------------------+-----------------------------------------------------------+------------------------------+ 

.. _Examples:

Examples
========

.. _Routines:

Routines
========
"""

import os
import shutil
import re
import logging
import pdb
import tempfile
import subprocess as sp

import pandas as pd
import numpy as np

import sandy
from sandy.settings import SandyError
from sandy.formats.endf6 import Endf6

__author__ = "Luca Fiorito"
__all__ = ["process", "process_proton"]

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

def get_suffix(temp):
    """Determine suffix saccording to temperature value using aleph method.
    
    Parameters
    ----------
    temp : `float`
        processing temperature
    
    Returns
    -------
    `int`
        suffix number
    """
    if abs(temp) >= 1000:
        suff = round(abs(temp)/100)
    else:
        suff = round(abs(temp/10/5))*5
        if suff%10 == 0:
            suff = round(suff/10)
    return suff



def _moder_input(nin, nout, **kwargs):
    """Write moder input.
    
    Parameters
    ----------
    nin : `int`
        tape number for input file
    nout : `int`
        tape number for output file
    
    Returns
    -------
    `str`
        moder input text
    """
    text = ["moder"]
    text.append("{:d} {:d} /".format(nin, nout))
    return "\n".join(text) + "\n"


def _reconr_input(endfin, pendfout, mat,
                  header="sandy runs njoy",
                  err=0.001,
                  **kwargs):
    """Write reconr input.
    
    Parameters
    ----------
    endfin : `int`
        tape number for input ENDF-6 file
    pendfout : `int`
        tape number for output PENDF file
    mat : `int`
        MAT number
    header : `str`
        file header (default is "sandy runs njoy")
    err : `float`
        tolerance (default is 0.001)

    Returns
    -------
    `str`
        reconr input text
    """
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
    """Write broadr input.
    
    Parameters
    ----------
    endfin : `int`
        tape number for input ENDF-6 file
    pendfin : `int`
        tape number for input PENDF file
    pendfout : `int`
        tape number for output PENDF file
    mat : `int`
        MAT number
    temperatures : iterable of `float`
        iterable of temperature values in K (default is 293.6 K)
    err : `float`
        tolerance (default is 0.001)

    Returns
    -------
    `str`
        broadr input text
    """
    text = ["broadr"]
    text += ["{:d} {:d} {:d} /".format(endfin, pendfin, pendfout)]
    text += ["{:d} {:d} 0 0 0. /".format(mat, len(temperatures))]
    text += ["{} /".format(err)]
    text += [" ".join(map("{:.1f}".format, temperatures)) + " /"]
    text += ["0 /"]
    return "\n".join(text) + "\n"

def _thermr_input(endfin, pendfin, pendfout, mat,
                  temperatures=[293.6], angles=20, iprint=False,
                  err=0.001, emax=10, **kwargs):
    """Write thermr input for free-gas.
    
    Parameters
    ----------
    endfin : `int`
        tape number for input ENDF-6 file
    pendfin : `int`
        tape number for input PENDF file
    pendfout : `int`
        tape number for output PENDF file
    mat : `int`
        MAT number
    temperatures : iterable of `float`
        iterable of temperature values in K (default is 293.6 K)
    angles : `int`
        number of equi-probable angles (default is 20)
    iprint : `bool`
        print option (default is `False`)
    err : `float`
        tolerance (default is 0.001)
    emax : `float`
        maximum energy for thermal treatment (default is 10 eV)

    Returns
    -------
    `str`
        thermr input text
    """
    text = ["thermr"]
    text += ["{:d} {:d} {:d} /".format(endfin, pendfin, pendfout)]
    text += ["0 {:d} {:d} {:d} 1 0 0 1 221 {:d} /".format(mat, angles, len(temperatures), int(iprint))]
    text += [" ".join(map("{:.1f}".format, temperatures)) + " /"]
    text += ["{} {} /".format(err, emax)]
    return "\n".join(text) + "\n"

def _purr_input(endfin, pendfin, pendfout, mat,
                temperatures=[293.6], sig0=[1e10], bins=20, ladders=32,
                iprint=False,
                **kwargs):
    """Write purr input.
    
    Parameters
    ----------
    endfin : `int`
        tape number for input ENDF-6 file
    pendfin : `int`
        tape number for input PENDF file
    pendfout : `int`
        tape number for output PENDF file
    mat : `int`
        MAT number
    temperatures : iterable of `float`
        iterable of temperature values in K (default is 293.6 K)
    sig0 : iterable of `float`
        iterable of dilution values in barns (default is 1e10 b)
    bins : `int`
        number of probability bins (default is 20)
    ladders : `int`
        number of resonance ladders (default is 32)
    iprint : `bool`
        print option (default is `False`)

    Returns
    -------
    `str`
        purr input text
    """
    text = ["purr"]
    text += ["{:d} {:d} {:d} /".format(endfin, pendfin, pendfout)]
    text += ["{:d} {:d} {:d} {:d} {:d} {:d} /".format(mat, len(temperatures), len(sig0), bins, ladders, int(iprint))]
    text += [" ".join(map("{:.1f}".format, temperatures)) + " /"]
    text += [" ".join(map("{:.2E}".format, sig0)) + " /"]
    text += ["0 /"]
    return "\n".join(text) + "\n"

def _gaspr_input(endfin, pendfin, pendfout,
                 **kwargs):
    """Write gaspr input.
    
    Parameters
    ----------
    endfin : `int`
        tape number for input ENDF-6 file
    pendfin : `int`
        tape number for input PENDF file
    pendfout : `int`
        tape number for output PENDF file

    Returns
    -------
    `str`
        gaspr input text
    """
    text = ["gaspr"]
    text += ["{:d} {:d} {:d} /".format(endfin, pendfin, pendfout)]
    return "\n".join(text) + "\n"

def _unresr_input(endfin, pendfin, pendfout, mat,
                  temperatures=[293.6], sig0=[1e10], iprint=False,
                  **kwargs):
    """Write unresr input.
    
    Parameters
    ----------
    endfin : `int`
        tape number for input ENDF-6 file
    pendfin : `int`
        tape number for input PENDF file
    pendfout : `int`
        tape number for output PENDF file
    mat : `int`
        MAT number
    temperatures : iterable of `float`
        iterable of temperature values in K (default is 293.6 K)
    sig0 : iterable of `float`
        iterable of dilution values in barns (default is 1e10 b)
    iprint : `bool`
        print option (default is `False`)

    Returns
    -------
    `str`
        unresr input text
    """
    text = ["unresr"]
    text += ["{:d} {:d} {:d} /".format(endfin, pendfin, pendfout)]
    text += ["{:d} {:d} {:d} {:d} /".format(mat, len(temperatures), len(sig0), int(iprint))]
    text += [" ".join(map("{:.1f}".format, temperatures)) + " /"]
    text += [" ".join(map("{:.2E}".format, sig0)) + " /"]
    text += ["0 /"]
    return "\n".join(text) + "\n"

def _heatr_input(endfin, pendfin, pendfout, mat, pks,
                 local=False, iprint=False,
                 **kwargs):
    """Write heatr input.
    
    Parameters
    ----------
    endfin : `int`
        tape number for input ENDF-6 file
    pendfin : `int`
        tape number for input PENDF file
    pendfout : `int`
        tape number for output PENDF file
    mat : `int`
        MAT number
    pks : iterable of `int`
        iterable of MT numbers for partial kermas
    local : `bool`
        option to deposit gamma rays locally (default is `False`)
    iprint : `bool`
        print option (default is `False`)

    Returns
    -------
    `str`
        heatr input text
    """
    text = ["heatr"]
    text += ["{:d} {:d} {:d} 0 /".format(endfin, pendfin, pendfout)]
    text += ["{:d} {:d} 0 0 {:d} {:d} /".format(mat, len(pks), int(local), int(iprint))]
    text += [" ".join(map("{:d}".format, pks)) + " /"]
    return "\n".join(text) + "\n"

def _acer_input(endfin, pendfin, aceout, dirout, mat,
                temp=293.6, iprint=False, itype=1, suff=".00",
                header="sandy runs acer",
                photons=True,
                **kwargs):
    """Write acer input for fast data.
    
    Parameters
    ----------
    endfin : `int`
        tape number for input ENDF-6 file
    pendfin : `int`
        tape number for input PENDF file
    aceout : `int`
        tape number for output ACE file
    dirout : `int`
        tape number for output ACE file
    mat : `int`
        MAT number
    temp : `float`
        temperature in K (default is 293.6 K)
    local : `bool`
        option to deposit gamma rays locally (default is `False`)
    iprint : `bool`
        print option (default is `False`)
    itype : `int`
        ace output type: 1, 2, or 3 (default is 1)
    suff : `str`
        id suffix for zaid (default is ".00")
    header : `str`
        descriptive character string of max. 70 characters (default is "sandy runs acer")
    photons : `bool`
        detailed photons (default is `True`)

    Returns
    -------
    `str`
        acer input text
    """
    text = ["acer"]
    text += ["{:d} {:d} 0 {:d} {:d} /".format(endfin, pendfin, aceout, dirout)]
    text += ["1 {:d} {:d} {} 0 /".format(int(iprint), itype, suff)]
    text += ["'{}'/".format(header)]
    text += ["{:d} {:.1f} /".format(mat, temp)]
    text += ["1 {:d} /".format(int(photons))]
    text += ["/"]
    return "\n".join(text) + "\n"

def _errorr_input(endfin, pendfin, errorrout, mat,
                  ign=2, iwt=2, temp=293.6, iprint=False,
                  **kwargs):
    """Write acer input for fast data.
    
    Parameters
    ----------
    endfin : `int`
        tape number for input ENDF-6 file
    pendfin : `int`
        tape number for input PENDF file
    errorrout : `int`
        tape number for output ERRORR file
    mat : `int`
        MAT number
    ign : `int`
        neutron group option (default is 2, csewg 239-group structure)
    iwt : `int`
        weight function option (default is 2, constant)
    temp : `float`
        temperature in K (default is 293.6 K)
    iprint : `bool`
        print option (default is `False`)

    Returns
    -------
    `str`
        acer input text
    """
    text = ["errorr"]
    text += ["{:d} {:d} 0 {:d} 0 /".format(endfin, pendfin, errorrout)]
    text += ["{:d} {:d} {:d} {:d} 1 /".format(mat, ign, iwt, int(iprint))]
    text += ["{:d} {:.1f} /".format(int(iprint), temp)]
    text += ["0 33 /"]
    return "\n".join(text) + "\n"

def _run_njoy(text, inputs, outputs, exe=None):
    """
    Run njoy executable.
    
    .. Important::

        In `Python 3` you need to convert input string to bytes with a 
        `encode()` function
    
    Parameters
    ----------
    exe : `str` or `None`
        njoy executable: if `None` (default) search firts `NJOY` env variable and then 
        in `PATH`
    inputs : `map`
        map of {`tape` : `file`) for input files
    outputs : `map`
        map of {`tape` : `file`) for ouptut files
    text : `str`
        njoy input file passed to `Popen` as `stdin` (it must be encoded first)
    """
    if not exe:
        if "NJOY" in os.environ:
            exe = os.environ["NJOY"]
    if not exe:
        for try_exe in ["njoy2016", "njoy", "njoy2012", "xnjoy"]:
            exe = sandy.tools.which(try_exe)
            if exe:
                break
    if not exe:
        raise SandyError("could not find njoy executable")
    logging.debug("Use NJOY executable '{}'".format(exe))
    stdout = stderr = None
    stdin = text.encode()
    with tempfile.TemporaryDirectory() as tmpdir:
        logging.debug("Create temporary directory '{}'".format(tmpdir))
        for tape,src in inputs.items():
            shutil.copy(src, os.path.join(tmpdir, tape))
        process = sp.Popen(exe,
                           shell=True,
                           cwd=tmpdir,
                           stdin=sp.PIPE, 
                           stdout=stdout, 
                           stderr=stderr)
        stdoutdata, stderrdata = process.communicate(input=stdin)
        logging.debug(stdoutdata)
        logging.debug(stderrdata)
        if process.returncode != 0:
            raise SandyError("process status={}, cannot run njoy executable".format(process.returncode))
        for tape,dst in outputs.items():
            path = os.path.split(dst)[0]
            if path:
                os.makedirs(path, exist_ok=True)
            shutil.move(os.path.join(tmpdir, tape), dst)



def process(endftape, pendftape=None,
            kermas=[302, 303, 304, 318, 402, 442, 443, 444, 445, 446, 447],
            temperatures=[293.6],
            suffixes=None,
            broadr=True,
            thermr=True,
            unresr=False,
            heatr=True,
            gaspr=True,
            purr=True,
            errorr=False,
            acer=True,
            wdir="", dryrun=False, tag="", exe=None, keep_pendf=True,
            route="0", addpath=None,
            **kwargs):
    """Run sequence to process file with njoy.
    
    Parameters
    ----------
    pendftape : `str`
        use this PENDF file and skip module reconr (defult is `None`: run reconr) 
    kermas : iterable of `int`
        MT numbers for partial kermas to pass to heatr.
        Default is:
            
            - `MT=302` : KERMA elastic
            - `MT=303` : KERMA non-elastic
            - `MT=304` : KERMA inelastic
            - `MT=318` : KERMA fission
            - `MT=402` : KERMA radiative capture
            - `MT=442` : total photon KERMA contribution
            - `MT=443` : total kinematic KERMA (kinematic Limit)
            - `MT=444` : total damage energy production cross section
            - `MT=445` : elastic damage energy production cross section
            - `MT=446` : inelastic damage energy production cross section
            - `MT=447` : neutron disappearance damage energy production cross section
        .. note: `MT=301` is the KERMA total (energy balance) and is always calculated
    temperatures : iterable of `float`
        iterable of temperature values in K (default is 293.6 K)
    suffixes : iterable of `int`
        iterable of suffix values for ACE files (default is `None`, use internal routine to determine suffixes)
        .. warning: must have the same number of entries then `temperatures`
        must match the number of entries in `temperatures`
    broadr : `bool`
        option to run module broadr (default is `True`)
    thermr : `bool`
        option to run module thermr (default is `True`)
    unresr : `bool`
        option to run module unresr (default is `False`)
    heatr : `bool`
        option to run module heatr (default is `True`)
    gaspr : `bool`
        option to run module gapr (default is `True`)
    purr : `bool`
        option to run module purr (default is `True`)
    errorr : `bool`
        option to run module errorr (default is `False`)
    acer : `bool`
        option to run module acer (default is `True`)
    wdir : `str`
        working directory (absolute or relative) where all output files are
        saved
        .. note: `wdir` will appear as part of the `filename` in any `xsdir` file if `addpath` is not set
    addpath : `str`
        path to add in xsdir, by default use `wdir`
    dryrun : `bool`
        option to produce the njoy input file without running njoy
    tag : `str`
        tag to append to each output filename beofre the extension (default is `None`)
        .. hint: to process JEFF-3.3 files you could set `tag = "_j33"`
    exe : `str`
        njoy executable (with path)
        .. note: if no executable is given, SANDY looks for a default executable in `PATH` and in env variable `NJOY`
    keep_pendf : `bool`
        save output PENDF file (default is `True`)
    route : `str`
        xsdir "route" parameter (default is "0")
    
    Returns
    -------
    input : `str`
        njoy input text
    inputs : `map`
        map of {`tape` : `file`) for input files
    outputs : `map`
        map of {`tape` : `file`) for ouptut files
    """
    tape = Endf6.from_file(endftape)
    mat = tape.mat[0]
    info = tape.read_section(mat, 1, 451)
    meta = info["LISO"]
    za = int(info["ZA"])
    za_new = za + meta*100 + 300 if meta else za
    inputs = {}
    outputs = {}
    # Only kwargs are passed to NJOY inputs, therefore add temperatures and mat
    kwargs.update({"temperatures" : temperatures, "mat" : mat})
    # Check input args
    if not suffixes:
        suffixes = list(map(get_suffix, temperatures))
    if len(suffixes) != len(temperatures):
        raise SandyError("number of suffixex must match number of temperatures")
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
        o = 30
        text += _moder_input(-p, o)
        outputs["tape{}".format(o)] = os.path.join(wdir, "{}{}.pendf".format(za_new, tag))
    if errorr:
        for i,(temp,suff) in enumerate(zip(temperatures, suffixes)):
            o = 33 + i
            kwargs["temp"] = temp
            kwargs["suff"] = suff = ".{:02d}".format(suff)
            text += _errorr_input(-e, -p, o, **kwargs)
            outputs["tape{}".format(o)] = os.path.join(wdir, "{}{}{}.errorr".format(za_new, tag, suff))
    if acer:
        for i,(temp,suff) in enumerate(zip(temperatures, suffixes)):
            a = 50 + i
            x = 70 + i
            kwargs["temp"] = temp
            kwargs["suff"] = suff = ".{:02d}".format(suff)
            text += _acer_input(-e, -p, a, x, **kwargs)
            outputs["tape{}".format(a)] = os.path.join(wdir, "{}{}{}c".format(za_new, tag, suff))
            outputs["tape{}".format(x)] = os.path.join(wdir, "{}{}{}c.xsd".format(za_new, tag, suff))
    text += "stop"
    if not dryrun:
        _run_njoy(text, inputs, outputs, exe=exe)
        if acer:
            # Change route and filename in xsdir file.
            for i,(temp,suff) in enumerate(zip(temperatures, suffixes)):
                a = 50 + i
                x = 70 + i
                acefile = outputs["tape{}".format(a)]
                if addpath is None:
                    filename = acefile
                else:
                    filename = os.path.basename(acefile)
                    if addpath:
                        filename = os.path.join(addpath, filename)
                xsdfile = outputs["tape{}".format(x)]
                text_xsd = open(xsdfile).read(). \
                                         replace("route", route). \
                                         replace("filename", filename)
                text_xsd = " ".join(text_xsd.split())
                # If isotope is metatable rewrite ZA in xsdir and ace as ZA = Z*1000 + 300 + A + META*100.
                if meta:
                    pattern = '{:d}'.format(za) + '\.(?P<ext>\d{2}[ct])'
                    found = re.search(pattern, text_xsd)
                    ext = found.group("ext")
                    text_xsd = text_xsd.replace("{:d}.{}".format(za, ext), "{:d}.{}".format(za_new, ext), 1)
                    text_ace = open(acefile).read()
                    text_ace = text_ace.replace("{:d}.{}".format(za, ext), "{:d}.{}".format(za_new, ext), 1)
                    with open(acefile, 'w') as f:
                        f.write(text_ace)
                with open(xsdfile, 'w') as f:
                    f.write(text_xsd)
    return text, inputs, outputs

def process_proton(endftape, wdir="", dryrun=False, tag="", exe=None, route="0", **kwargs):
    """Run sequence to process proton file with njoy.
    
    Parameters
    ----------
    wdir : `str`
        working directory (absolute or relative) where all output files are
        saved
        .. note:
            
            `wdir` will appear as part of the `filename` in 
            any `xsdir` file
    dryrun : `bool`
        option to produce the njoy input file without running njoy
    tag : `str`
        tag to append to each output filename beofre the extension (default is `None`)
        .. hint:
            to process JEFF-3.3 files you could set `tag = "_j33"`
    exe : `str`
        njoy executable (with path)
        .. note:
            If no executable is given, SANDY looks for a default executable in `PATH`
    route : `str`
        xsdir "route" parameter (default is "0")
    
    Returns
    -------
    input : `str`
        njoy input text
    inputs : `map`
        map of {`tape` : `file`) for input files
    outputs : `map`
        map of {`tape` : `file`) for ouptut files
    """
    tape = Endf6.from_file(endftape)
    mat = tape.mat[0]
    info = tape.read_section(mat, 1, 451)
    meta = info["LISO"]
    za = int(info["ZA"])
    za_new = za + meta*100 + 300 if meta else za
    inputs = {}
    outputs = {}
    kwargs["mat"] = mat
    inputs["tape20"] = endftape
    kwargs["temp"] = 0
    kwargs["suff"] = suff = ".00"
    text = _acer_input(20, 20, 50, 70, **kwargs)
    outputs["tape50"] = os.path.join(wdir, "{}{}{}h".format(za_new, tag, suff))
    outputs["tape70"] = os.path.join(wdir, "{}{}{}h.xsd".format(za_new, tag, suff))
    text += "stop"
    if not dryrun:
        _run_njoy(text, inputs, outputs, exe=exe)
        # Change route and filename in xsdir file.
        acefile = outputs["tape50"]
        xsdfile = outputs["tape70"]
        text_xsd = open(xsdfile).read(). \
                                 replace("route", route). \
                                 replace("filename", acefile)
        text_xsd = " ".join(text_xsd.split())
        # If isotope is metatable rewrite ZA in xsdir and ace as ZA = Z*1000 + 300 + A + META*100.
        if meta:
            pattern = '{:d}'.format(za) + '\.(?P<ext>\d{2}[ct])'
            found = re.search(pattern, text_xsd)
            ext = found.group("ext")
            text_xsd = text_xsd.replace("{:d}.{}".format(za, ext), "{:d}.{}".format(za_new, ext), 1)
            text_ace = open(acefile).read()
            text_ace = text_ace.replace("{:d}.{}".format(za, ext), "{:d}.{}".format(za_new, ext), 1)
            with open(acefile, 'w') as f:
                f.write(text_ace)
        with open(xsdfile, 'w') as f:
            f.write(text_xsd)
    return text, inputs, outputs