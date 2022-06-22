# -*- coding: utf-8 -*-
"""
Created on Fri Apr  8 08:50:02 2022

@author: aitor
"""

import os
import time
import shutil
import sys
import argparse
import logging
import tempfile
import multiprocessing as mp
import platform
import pytest
from functools import reduce


import numpy as np
import pandas as pd

import sandy
from sandy.settings import SandyError
from sandy.tools import is_valid_dir, is_valid_file
from sandy import njoy

__author__ = "Aitor Bengoechea"

mf_available = [8, 31, 33, 34, 35]


class sandy_object():
    global init

    def __init__(self, ftape, mf, energy_sequence):
        self.Endf6 = ftape
        extra_points = np.logspace(-5, 7, energy_sequence)
        ftape = ftape.filter_by(listmt=init.mt)
        if 8 in mf:
            self.Fy = sandy.Fy.from_endf6(ftape)
        if 31 in mf or 33 in mf:
            self.Xs = sandy.Xs.from_endf6(ftape)
        if 34 in mf:
            self.Lpc = sandy.Lpc.from_endf6(ftape).reshape(extra_points)
        if 35 in mf:
            self.Edistr = sandy.Edistr.from_endf6(ftape).reshape(extra_points)


def parse(iargs=None):
    """
    Parse command line arguments for sampling option.

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
    parser.add_argument('--pdf',
                        type=str,
                        default='normal',
                        metavar="pdf",
                        help="random numbers distribution (default = normal)")
    parser.add_argument('--tolerance',
                        type=float,
                        default=0.0,
                        metavar="tol",
                        help="Replace all eigenvalues smaller than a given"
                             "tolerance with zeros (default = 0.0)")
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


def get_fy_cov(endf6):
    """
    Generate fy samples divided by the MAT, MT and energy

    Parameters
    ----------
    endf6 : `sandy.Endf6`
        Endf6 file containing nfpy.

    Returns
    -------
    dfperts : `pd.Series`
        `sandy.core.Samples` objects divided by the MAT, MT and E.

    """
    global init
    fy = sandy.Fy.from_endf6(endf6.filter_by(listmat=init.mat, listmt=init.mt))
    fy_cov = fy.data.groupby(['MAT', 'MT', 'E'])["DFY"]\
               .apply(lambda x: sandy.CategoryCov.from_stdev(x))
    return fy_cov


def get_cov(endf6):
    """
    Function to obtain the covariance matrices for all mf entered by the user.

    Parameters
    ----------
    endf6 : `sandy.Endf6`
        Endf file.

    Returns
    -------
    cov : `dict`
        Dictionary containing as keys the mf entered by the user and as values
        objects `sandy.CategoryCov`. The dictionary is divided as follows:

                {mf: {`sandy.CategoryCov`}}

    Notes
    -----
    .. note:: FY data is divided in a differente way. The covariance matrix are
    divided according to the mat, mt and energy.

    """
    global init
    cov = {}
    if init.cov:
        cov[init.mf[0]] = sandy.CategoryCov.from_csv(init.cov)
        return cov

    # Covariance of Sandy modules:
    if 8 in init.mf and 8 in endf6.mf:
        cov[8] = get_fy_cov(endf6)

    # Covariance of the endf6 file:
    mf_process = []
    mf_extract = []
    if 34 in init.mf and 34 in endf6.mf:
        mf_extract.append(34)
    if 35 in init.mf and 35 in endf6.mf:
        mf_extract.append(35)
    if init.cov33csv:
        logging.warning("found argument '--cov33csv', will skip any other"
                        "covariance")
        cov[33] = sandy.CategoryCov.from_csv(
                    init.cov33csv,
                    index_col=[0, 1, 2],
                    header=[0, 1, 2],
                    )
    else:
        if 31 in init.mf and 31 in endf6.mf:
            mf_process.append(31)
        if 33 in init.mf and 33 in endf6.mf:
            mf_process.append(33)
        temp = 0 if len(init.temperatures) == 0 else init.temperatures

    # Endf6 covariance:
    if len(mf_process) + len(mf_extract) == 1:
        unique_mf = mf_process + mf_extract
        cov.update({unique_mf[0]: endf6.get_cov(process_mf=mf_process,
                                                mf=mf_extract,
                                                njoy=init.njoy,
                                                temperature=temp)})
    elif len(mf_process) + len(mf_extract) >= 1:
        cov.update(endf6.get_cov(process_mf=mf_process, mf=mf_extract,
                                 njoy=init.njoy, temperature=temp))
    return cov


def sample_manager(endf6):
    """
    Function to obtain the sample covariance matrices for all mf entered
    by the user.

    Parameters
    ----------
    endf6 : `sandy.Endf6`
        Endf6 file.

    Returns
    -------
    Dictionary containing the `sandy.Samples` objects divided by the mf
    (dict keys). The dictionary is divided as follows:

        {mf: {`sandy.Samples`}}

    Notes
    -----
    .. note:: FY data is divided in a differente way. The samples are
    divided according to the mat, mt and energy.
    """
    global init

    # Obtain the covariance:
    cov = get_cov(endf6)

    # Create the samples:
    def sample_creation(individual_cov):
        mf = individual_cov[0]
        mf_cov = individual_cov[1]
        seed = init.seed31 if mf == 31 else init.seed33 if mf == 33 else init.seed34 if mf == 34 else init.seed35 if mf == 35 else None
        if mf == 8:
            samples = mf_cov.apply(lambda x: x.sampling(init.samples,
                                                        seed=seed,
                                                        tolerance=init.tolerance,
                                                        pdf=init.pdf))
        elif mf == 35:
            samples = mf_cov.sampling(init.samples, seed=seed,
                                      tolerance=init.tolerance,
                                      pdf=init.pdf, relative=False)
        else:
            samples = mf_cov.sampling(init.samples, seed=seed,
                                      tolerance=init.tolerance,
                                      pdf=init.pdf)
        return (mf, samples)

    return dict(map(sample_creation, cov.items()))


def perturbation_manager(samples, endf6):
    """
    Function to produce `sandy.Endf6` objects with perturbed data in
    parallelised form if the operating system permits.

    Parameters
    ----------
    samples : `dict`
        Dictionary containing the `sandy.Samples` objects divided by the mf
        (dict keys).
    endf6 : `sandy.Endf6`
        Endf6 file.

    Returns
    -------
    Perturbed files written in the folder selected

    """
    global init, outname

    # Decide if the sampes are going to be created in series or parallel
    if platform.system() == "Windows":
        proc = 1
        logging.info("Running on Windows does not allow parallel "
                     "processing")
    else:
        proc = init.processes

    # Output:
    outname = init.outname if init.outname else os.path.basename(init.file)

    # Sample creation:
    for mat in endf6.to_series().index.get_level_values("MAT").unique():
        tape = endf6.filter_by(listmat=[mat])
        pert_objects = sandy_object(tape, samples.keys(), init.energy_sequence)
        if proc != 1:
            pool = mp.Pool(processes=proc)
        for i in range(1, init.samples + 1):
            if proc == 1:
                pert_by_mf(samples, pert_objects, i, mat)
            else:
                pool.apply_async(pert_by_mf, (samples, pert_objects, i, mat))
    return


def pert_by_mf(samples, pert_objects, i, mat):
    """

    Parameters
    ----------
    samples : `dict`
        Dictionary containing the `sandy.Samples` objects containing the
        sample information.

    Returns
    -------
    Endf6 with perturbed data in a text file

    """
    global init, outname
    t0 = time.time()
    i_ = i-1
    pert_data = [pert_objects.Endf6]

    # FY:
    if 8 in samples:
        sample = samples[8]
        fy = pert_objects.Fy
        for mat, mt, e in sample.index:
            sample_ = sample.loc[(mat, mt, e)].data.iloc[i_]
            fy = fy.custom_perturbation(sample_.values, mat=mat, e=e, mt=mt)
            pert_data.append(fy)

    # XS and Nubar:
    if 31 in samples or 33 in samples:
        if 31 in samples and 33 in samples:
            pert = pd.concat([samples[31].to_pert(smp=i_).data,
                             samples[33].to_pert(smp=i_).data],
                             axis=1)
            pert = sandy.Pert(pert)
        elif 31 in samples:
            pert = samples[31].to_pert(smp=i_)
        elif 33 in samples:
            pert = samples[33].to_pert(smp=i_)

        xs = pert_objects.Xs.custom_perturbation(pert)
        pert_data.append(xs)

    # LPC:
    if 34 in samples:
        pert = samples[34].to_pert(smp=i_)
        lpc = pert_objects.Lpc.custom_perturbation(pert)
        pert_data.append(lpc)

    # Edistr:
    if 35 in samples:
        pert = samples[35].to_pert(smp=i_)
        edistr = pert_objects.Edistr.custom_perturbation(pert)
        pert_data.append(edistr)

    # Perturbed Endf6 file
    pert_endf6 = reduce(lambda x, y: y.to_endf6(x), pert_data)
    print("Created sample {} for MAT {} in {:.2f} sec"
          .format(i, mat, time.time()-t0,))

    # Output files:
    output = os.path.join(init.outdir, '{}-{}'.format(outname, i))
    pert_endf6.to_file(output)
    return


def ace_files():
    """
    Function to decide if ACE files are created in parallel or serial

    """
    global init
    if platform.system() != "Windows":
        pool = mp.Pool(processes=init.processes)
    for i in range(1, init.samples + 1):
        if platform.system() == "Windows":
            _ace_files(i)
        else:
            pool.apply_async(_ace_files, (i))
    return


def _ace_files(i):
    """
    Run NJOY module to produce ACE files

    Parameters
    ----------
    i : `int`
        Number of the sample.
    """
    global init, outname
    smpfile = os.path.join(init.outdir, f'{outname}-{i}')
    kwargs = dict(
        purr=False,
        wdir=init.outdir,
        keep_pendf=False,
        pendftape=smpfile,
        tag=f"_{i}",
        temperatures=init.temperatures,
        err=0.005,
        addpath="",
        )
    fmt = sandy.formats.get_file_format(smpfile)
    if 31 in init.mf or 33 in init.mf:
        kwargs["pendftape"] = smpfile
        inp = init.file
    elif fmt == "endf6":
        inp = smpfile
    input, inputs, outputs = njoy.process(inp, **kwargs)
    return


def sampling(iargs=None):
    """
    Construct multivariate distributions with a unit vector for
    mean and with relative covariances taken from the evaluated files or
    process with NJOY. Perturbation factors are sampled with the same
    multigroup structure of the covariance matrix, and are applied to the
    pointwise data to produce the perturbed files.

    Parameters
    ----------
    iargs : `list` of `str`
        list of strings to parse. The default is taken from `sys.argv`.

    Returns
    -------
    pert_samples : `dict`
        Dictionary containing the `sandy.Endf6` objects with the perturbed data
        for the number of samples entered by the user. The dictionary is
        divided as follows:

                {mat: {Number of the sample: `sandy.Endf6`}}

    ftape: 'sandy.Endf6'
        Original `sandy.Endf6` object withouth perturbed data.
    """
    global init, ftape

    # Command line information and endf6 file information:
    init = parse(iargs)
    ftape = sandy.Endf6.from_file(init.file)
    # Check if the covariance information is available:
    mf_check = ftape.to_series().index.get_level_values("MF")\
                    .intersection(pd.Index(init.mf))
    if len(mf_check) == 0:
        raise sandy.Error("The selected MF was not found")

    # Samples:
    samples = sample_manager(ftape)

    # Endf + Pendf:
    # Check if NJOY has to be run:
    xs_check = ftape.to_series().index.get_level_values("MF")\
                    .intersection(pd.Index([3]))
    if len(xs_check) == 1 and 31 in init.mf or 33 in init.mf:
        pendf = ftape.get_pendf(njoy=init.njoy)
        with tempfile.TemporaryDirectory() as td:
            dst = os.path.join(td, "merged")
            ftape = ftape.merge_pendf(pendf)
            ftape.to_file(dst)

    # Perturbed endf:
    perturbation_manager(samples, ftape)
    
    # Produce ACE files:
    if init.acer:
        ace_files()

    return


def run():
    t0 = time.time()
    try:
        sampling()
    except SandyError as exc:
        logging.error(exc.args[0])
    print("Total running time: {:.2f} sec".format(time.time() - t0))


if __name__ == "__main__":
    run()
