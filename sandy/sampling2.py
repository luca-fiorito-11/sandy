# -*- coding: utf-8 -*-
"""
Created on Fri Apr  8 08:50:02 2022

@author: aitor
"""

import os
import time
import pdb
import shutil
import sys
import argparse
import logging
import tempfile
import multiprocessing as mp
import platform
import pytest
import pdb

import numpy as np
import pandas as pd

import sandy
from sandy.settings import SandyError
from sandy.tools import is_valid_dir, is_valid_file
from sandy import njoy

__author__ = "Aitor Bengoechea"

def parse(iargs=None):
    """Parse command line arguments for sampling option.

    Parameters
    ----------
    iargs : `list` of `str`
        list of strings to parse. The default is taken from `sys.argv`.

    Returns
    -------
    `argparse.Namespace`
        namespace object containing processed given arguments and/or default
        options.
    """
    description = "Produce perturbed files containing sampled parameters that "
    "represent the information\nstored in the evaluated nuclear "
    "data covariances"
    parser = argparse.ArgumentParser(
                        prog="sandy",
                        description=description,
                        formatter_class=argparse.RawTextHelpFormatter,
                        )
    parser.add_argument('file',
                        type=lambda x: is_valid_file(parser, x),
                        help="ENDF-6 or PENDF format file")
    parser.add_argument('--acer',
                        default=False,
                        action="store_true",
                        help="for each perturbed file, produce ACE files\n"
                             "(argument file must be in ENDF-6 format, not PENDF)\n(argument temperature is required)\n(default = False)")
    parser.add_argument('--cov', '-C',
                        type=lambda x: is_valid_file(parser, x),
                        help="file containing covariances")
    parser.add_argument('--cov33csv',
                        type=lambda x: is_valid_file(parser, x),
                        help="file containing xs/nubar covariances in csv "
                             "format")
    parser.add_argument('--debug',
                        default=False,
                        action="store_true",
                        help="turn on debug mode")
    parser.add_argument('--eig',
                        type=int,
                        default=10,
                        metavar="N",
                        help="print the first N eigenvalues of the evaluated covariance matrices\n(default = do not print)")
    parser.add_argument('--energy-sequence', '-E',
                        type=int,
                        metavar="EL",
                        default=49,
                        help=argparse.SUPPRESS)
    parser.add_argument('--errorr',
                        default=False,
                        action="store_true",
                        help="run NJOY module ERRORR to produce covariance "
                             "matrix for xs data (default = False)")
    parser.add_argument('--fission-yields', '-F',
                        default=False,
                        action="store_true",
                        help="input <file> contains fission yields")
    parser.add_argument('--mat',
                        type=int,
                        default=list(range(1, 10000)),
                        action='store',
                        nargs="+",
                        metavar="{1,..,9999}",
                        help="draw samples only from the selected MAT "
                             "sections (default = keep all)")
    parser.add_argument('--max-polynomial', '-P',
                        type=int,
                        help="Maximum order of Legendre polynomial coefficients considered for sampling (default = all)")
    parser.add_argument('--mf',
                        type=int,
                        default=[31, 33, 34, 35],
                        action='store',
                        nargs="+",
                        metavar="{31,33,34,35}",
                        help="draw samples only from the selected MF sections "
                             "(default = keep all)")
    parser.add_argument('--mt',
                        type=int,
                        default=list(range(1, 1000)),
                        action='store',
                        nargs="+",
                        metavar="{1,..,999}",
                        help="draw samples only from the selected MT sections "
                             "(default = keep all)")
    parser.add_argument('--njoy',
                        type=lambda x: is_valid_file(parser, x),
                        default=None,
                        help="NJOY executable "
                             "(default search PATH, and env variable NJOY)")
    parser.add_argument('--outdir', '-D',
                        metavar="DIR",
                        default=os.getcwd(),
                        type=lambda x: is_valid_dir(parser, x, mkdir=True),
                        help="target directory where outputs are stored\n(default = current working directory)\nif it does not exist it will be created")
    parser.add_argument('--outname', '-O',
                        type=str,
                        help="basename for the output files "
                             "(default is the the basename of <file>.)")
    parser.add_argument('--processes', '-N',
                        type=int,
                        default=1,
                        help="number of worker processes (default = 1)")
    parser.add_argument('--samples', '-S',
                        type=int,
                        default=200,
                        help="number of samples (default = 200)")
    parser.add_argument('--seed31',
                        type=int,
                        default=None,
                        metavar="S31",
                        help="seed for random sampling of MF31 covariance "
                             "matrix (default = random)")
    parser.add_argument('--seed33',
                        type=int,
                        default=None,
                        metavar="S33",
                        help="seed for random sampling of MF33 covariance "
                             "matrix (default = random)")
    parser.add_argument('--seed34',
                        type=int,
                        default=None,
                        metavar="S34",
                        help="seed for random sampling of MF34 covariance "
                             "matrix (default = random)")
    parser.add_argument('--seed35',
                        type=int,
                        default=None,
                        metavar="S35",
                        help="seed for random sampling of MF35 covariance "
                             "matrix (default = random)")
    parser.add_argument('--temperatures', '-T',
                        default=[],
                        type=float,
                        action='store',
                        nargs="+",
                        metavar="T",
                        help="for each perturbed file, produce ACE files at "
                             "given temperatures")
    parser.add_argument("--version", "-v",
                        action='version',
                        version='%(prog)s {}'.format(sandy.__version__),
                        help="SANDY's version.")
    init = parser.parse_known_args(args=iargs)[0]
    if init.acer and not init.temperatures:
        parser.error("--acer requires --temperatures")
    if init.acer and sandy.formats.get_file_format(init.file) != "endf6":
        parser.error("--acer requires file in 'endf6' format")
    return init

def get_data(endf6, mf):
    global init
    data = endf6.filter_by(listmf=[mf],
                           listmt=init.mt,
                           listmat=init.mat)
    # Data to the appropiate object:
    if mf == 8:
        module = sandy.Fy.from_endf6(data)
    elif mf == 31 or mf == 33:
        module = sandy.Xs.from_endf6(data)
    return module

def get_fy_samples(endf6):
    """
    Generate fy samples divided by the MAT, MT and energy

    Parameters
    ----------
    endf6 : `sandy.Endf6`
        DESCRIPTION.

    Returns
    -------
    dfperts : `pd.Series`
        `sandy.core.Samples` objects divided by the MAT, MT and energy.

    """
    fy = get_data(endf6, 8)
    if not fy.data.empty:
        dfperts = fy.data\
                    .groupby(['MAT', 'MT', 'E'])["DFY"]\
                    .apply(lambda x: sandy.CategoryCov.from_var(x)
                                          .sampling(init.samples))
#        if init.debug:
#            PertFy.to_csv("perts_mf8.csv")
    return dfperts


def get_nubar_samples(errorr):
    return


def get_xs_samples(errorr):
    global init, ftape
    cov = errorr if init.cov33csv else errorr.get_cov()
    return cov.sampling(init.samples)


def get_lpc_samples(endf6):
    conditions = {"MAT": init.mat, "MT": init.mt, "MF": 34}
    lpc = endf6._filters(conditions)
    return


def get_chi_samples(endf6):
    return


def sample_manager(endf6):
    """
    

    Parameters
    ----------
    endf6 : `sandy.Endf6`
        DESCRIPTION.

    Returns
    -------
    samples : `dict`
        Dictionary containing the `sandy.Samples` objects containing the
        samples.

    """
    global init
    samples = {}
    if 8 in init.mf and 8 in endf6.mf:
        pfy = get_fy_samples(endf6)
        samples[8] = pfy
    if 34 in init.mf and 34 in endf6.mf:
        plpc = get_lpc_samples(endf6)
        samples[34] = plpc
    if 35 in init.mf and 35 in endf6.mf:
        pchi = get_chi_samples(endf6)
        samples[35] = pchi
    if init.cov33csv:
        logging.warning("found argument '--cov33csv', will skip any other"
                        " covariance")
        errorr = sandy.CategoryCov.from_csv(
                    init.cov33csv,
                    index_col=[0, 1, 2],
                    header=[0, 1, 2],
                    )
    elif init.errorr:
        # Obtain errorr file:
        errorr = endf6.get_errorr(njoy=init.njoy)
    # Create nubar and xs samples with errorr file:
    if 31 in init.mf and 31 in endf6.mf:
        pnu = get_nubar_samples(errorr)
        samples[31] = pnu
    if 33 in init.mf and 33 in endf6.mf:
        pxs = get_xs_samples(errorr)
        samples[33] = pxs
    return samples


def custom_perturbation_mf_8(pert_samples, endf6):

    fy = get_data(endf6, 8)
    fypert = fy.perturb(pert_samples)
    newtape = newtape.update_fy(fypert) #Me he quedado aqui
    return


def custom_perturbation_mf_31(pert_samples, endf6):
    global init

    return


def custom_perturbation_mf_33(pert_samples, endf6):
    global init

    return


def custom_perturbation_mf_34(pert_samples, endf6):
    global init

    return


def custom_perturbation_mf_35(pert_samples, endf6):
    global init

    return

def perturbation_manager(samples, formatted_file):
    """
    

    Parameters
    ----------
    samples : `dict`
        Dictionary containing the `sandy.Samples` objects containing the
        sample.

    Returns
    -------
    None.

    """
    global init
    pert_endf6 = formatted_file.copy()
    if 8 in samples:
        mat, mt, e = samples.index.values[0]
        pert_endf6 = custom_perturbation_mf_8(samples, pert_endf6)
    if 31 in samples:
        pert_endf6 = custom_perturbation_mf_31(samples, pert_endf6)
    if 33 in samples:
        pert_endf6 = custom_perturbation_mf_33(samples, pert_endf6)
    if 34 in samples:
        pert_endf6 = custom_perturbation_mf_34(samples, pert_endf6)
    if 35 in samples:
        pert_endf6 = custom_perturbation_mf_35(samples, pert_endf6)
    return pert_endf6


def sampling(iargs=None):
    """
    

    Parameters
    ----------
    iargs : `list` of `str`
        list of strings to parse. The default is taken from `sys.argv`.

    Returns
    -------
    None.

    """
    global init, ftape
    # Command line information and endf6 file information:
    init = parse(iargs)
    ftape = sandy.Endf6.from_file(init.file)
    # Samples:
    samples = sample_manager(ftape)
    # Endf + Pendf:
    pendf = ftape.get_pendf(njoy=init.njoy)
    with tempfile.TemporaryDirectory() as td:
        dst = os.path.join(td, "merged")
        ftape = ftape.merge_pendf(pendf)
        ftape.to_file(dst)
    # Perturbed endf:
    i = 0
    pert_endf6 = perturbation_manager(samples, ftape, i)
    return pert_endf6


def run():
    t0 = time.time()
    try:
        sampling()
    except SandyError as exc:
        logging.error(exc.args[0])
    print("Total running time: {:.2f} sec".format(time.time() - t0))


if __name__ == "__main__":
    run()
