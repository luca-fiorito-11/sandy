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

mf_available = [8, 31, 33, 34, 35]


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
        Endf file

    Returns
    -------
    dfperts : `pd.Series`
        `sandy.core.Samples` objects divided by the MAT, MT and energy.

    """
    global init
    fy = sandy.Fy.from_endf6(endf6.filter_by(listmat=init.mat))
    fy_cov = fy.data.groupby(['MAT', 'MT', 'E'])["DFY"]\
               .apply(lambda x: sandy.CategoryCov.from_var(x)) # o from_stdev
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
    samples : `dict`
        Dictionary containing the `sandy.Samples` objects divided by the mf
        (dict keys). The dictionary is divided as follows:

                {mf: {`sandy.Samples`}}

    Notes
    -----
    .. note:: FY data is divided in a differente way. The samples are
    divided according to the mat, mt and energy.
    """
    global init
    samples = {}
    # Obtain the covariance:
    cov = get_cov(endf6)

    # Create the samples:
    for mf, mf_cov in cov.items():
        seed = init.seed31 if mf == 31 else init.seed33 if mf == 33 else init.seed34 if mf == 34 else init.seed35 if mf == 35 else None
        if mf == 8:
            samples[mf] = mf_cov.apply(lambda x: x.sampling(init.samples,
                                                            seed=seed,
                                                            tolerance=init.tolerance,
                                                            pdf=init.pdf))
        else:
            samples[mf] = mf_cov.sampling(init.samples, seed=seed,
                                          tolerance=init.tolerance,
                                          pdf=init.pdf)
    return samples


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
    pert_samples : `dict`
        Dictionary containing the `sandy.Endf6` objects with the perturbed data
        for the number of samples entered by the user. The dictionary is
        divided as follows:

                {mat: {Number of the sample: `sandy.Endf6`}}

    """
    global init
    pert_samples = {}
    if platform.system() == "Windows":
        proc = 1
        logging.info("Running on Windows does not allow parallel "
                     "processing")
    else:
        proc = init.processes
    for mat in endf6.to_series().index.get_level_values("MAT").unique():
        tape = endf6.filter_by(listmat=[mat])
        seq = range(1, init.samples + 1)
        if proc == 1:
            outs = {i: pert_by_mf(samples, tape, mat, i) for i in seq}
        else:
            pool = mp.Pool(processes=proc)
            outs = {i: pool.apply_async(pert_by_mf, (samples, tape, mat, i))
                    for i in seq}
            outs = {i: out.get() for i, out in outs.items()}
            pool.close()
            pool.join()
        pert_samples.update({mat: outs})
    return pert_samples


def pert_by_mf(samples, endf6, mat, i):
    """

    Parameters
    ----------
    samples : `dict`
        Dictionary containing the `sandy.Samples` objects containing the
        sample information.

    Returns
    -------
    `pert_endf6`: `sandy.Endf6`
        Endf6 with perturbed data.

    """
    global init
    pert_endf6 = sandy.Endf6(endf6.data.copy())
    t0 = time.time()
    i_ = i-1
    if 8 in samples:
        pert_endf6 = custom_perturbation_mf_8(samples[8],
                                              pert_endf6, mat, i_)
    if 31 in samples:
        pert_endf6 = custom_perturbation_mf_31(samples[31].data.iloc[i_],
                                               pert_endf6, mat)
    if 33 in samples:
        pert_endf6 = custom_perturbation_mf_33(samples[33].data.iloc[i_],
                                               pert_endf6, mat, i_)
    if 34 in samples:
        pert_endf6 = custom_perturbation_mf_34(samples[34].data.iloc[i_],
                                               pert_endf6, mat, i_)
    if 35 in samples:
        pert_endf6 = custom_perturbation_mf_35(samples[35].data.iloc[i_],
                                               pert_endf6, mat, i_)
    print("Created sample {} for MAT {} in {:.2f} sec"
          .format(i, mat, time.time()-t0,))
    return pert_endf6


def custom_perturbation_mf_8(sample, endf6, mat, i):
    """
    Perturb the fy data for a selected mat, mt and energy.

    Parameters
    ----------
    sample : `dict`
        Dictionay containing the Samples (`sandy.Samples` objects) divided by
        the mat, mt and enegy.
    endf6 :`sandy.Endf6`
        Endf6 file.
    mat : `int`
        MAT number.
    i : `int`
        Number of the sample.

    Returns
    -------
    `sandt.Endf6`
        Endf6 file with the pertubed fy data.

    """
    fy = sandy.Fy.from_endf6(endf6)
    cond_mat = sample.index.get_level_values("MAT").isin([mat])
    sample_mat = sample.loc[cond_mat]
    for mat, mt, e in sample_mat.index:
        sample_ = sample_mat.loc[(mat, mt, e)].data.iloc[i]
        fy.custom_perturbation(sample_.values, mat=mat, e=e, mt=mt)
    return fy.to_endf6(endf6)


def custom_perturbation_mf_31(sample, endf6, mat):
    """
    Perturb the nubar data for a selected mat.

    Parameters
    ----------
    sample : `dict`
        Sample.
    endf6 :`sandy.Endf6`
        Endf6 file for pointwise nubar process with NJOY (module PENDF).
    mat : `int`
        MAT number.

    Returns
    -------
    `sandt.Endf6`
        Endf6 file with the pertubed nubar data.

    """
    pert = sandy.Pert(pd.Series(sample.values,
                                index=sample.index.get_level_values("E").left))
    xs = sandy.Xs.from_endf6(endf6)
    for mt in pert.index:
        if mt in xs.redundant_xs:
            for mt_redundant in xs.data.columns.get_level_values("MT").intersection(xs.redundant_xs[mt]):
                xs = xs.custom_perturbation(mat,
                                            mt_redundant,
                                            pert[mt])
        else:
            xs = xs.custom_perturbation(mat,
                                        mt,
                                        pert[mt])
    return xs._reconstruct_sums().to_endf6(endf6)


def custom_perturbation_mf_33(sample, endf6, mat, i):
    """
    Perturb the xs data for a selected mat.

    Parameters
    ----------
    sample : `sandy.Samples`
        Sample.
    endf6 :`sandy.Endf6`
        Endf6 file for pointwise nubar process with NJOY (module PENDF).
    mat : `int`
        MAT number.
    i : `int`
        Number of the sample.

    Returns
    -------
    `sandt.Endf6`
        Endf6 file with the pertubed nubar data.
    """
    xs = sandy.Xs.from_endf6(endf6)
    pert = sample.reset_index().query(f"MAT == {mat}").groupby('MT')\
        .apply(lambda x: sandy.Pert(pd.Series(x[i].values,
                                              index=x["E"].values)))
    for mt in pert.index:
        if mt in xs.redundant_xs:
            for mt_redundant in xs.data.columns.get_level_values("MT").intersection(xs.redundant_xs[mt]):
                xs = xs.custom_perturbation(mat,
                                            mt_redundant,
                                            pert[mt])
        else:
            xs = xs.custom_perturbation(mat,
                                        mt,
                                        pert[mt])
    return xs._reconstruct_sums().to_endf6(endf6)


def custom_perturbation_mf_34(sample, endf6, mat, i):
    """
    Perturb the lpc data for a selected mat.

    Parameters
    ----------
    sample : `sandy.Samples`
        Sample.
    endf6 :`sandy.Endf6`
        Endf6 file.
    mat : `int`
        MAT number.
    i : `int`
        Number of the sample.

    Returns
    -------
    `sandt.Endf6`
        Endf6 file with the pertubed lpc data.
    """
    extra_points = np.logspace(-5, 7, init.energy_sequence)
    conditions = {'MT': 2, 'MAT': mat}
    lpc = sandy.Lpc.from_endf6(endf6)._filters(conditions)\
        ._add_points(extra_points)
    pert = sample.reset_index().query(f"MAT == {mat}").groupby('L')\
        .apply(lambda x: sandy.Pert(pd.Series(x[i].values,
                                              index=x["E"])))
    for p in pert.index:
        lpc = lpc.custom_perturbation(mat, 2, p, pert[p])
    return lpc.to_endf6(endf6)


def custom_perturbation_mf_35(sample, endf6, mat, i):
    """
    Perturb the energy distribution data for a selected mat.

    Parameters
    ----------
    sample : `sandy.Samples`
        Sample.
    endf6 :`sandy.Endf6`
        Endf6 file.
    mat : `int`
        MAT number.
    i : `int`
        Number of the sample.

    Returns
    -------
    `sandt.Endf6`
        Endf6 file with the pertubed energy distribution data.
    """
    extra_points = np.logspace(-5, 7, init.energy_sequence)
    edistr = sandy.Edistr.from_endf6(endf6).add_energy_points(9237, 18, 0,
                                                              extra_points)
    pert = sample.reset_index().query(f"MAT == {mat}").groupby(['ELO', 'EHI'])\
        .apply(lambda x: sandy.Pert(pd.Series(x[i].values,
                                              index=x["E"])))
    for elo, ehi in pert.index:
        edistr = edistr.custom_perturbation(pert[elo, ehi],
                                            mat,
                                            18,
                                            0,
                                            elo,
                                            ehi)
    return edistr.to_endf6(endf6)


def _to_file(frame, ismp, mat, outname):
    """
    Write each sample in the outdir and with the outnames introduced by the
    user.

    Parameters
    ----------
    frame : `pd.Series`
        Pandas series with perturbed information divided by mat and the number
        of the sample.
    ismp : `int`
        Number of the sample.
    mat : `int`
        MAT number.
    outname : `str`
        Outname introduced by the user for the sample.
    """
    output = os.path.join(init.outdir, '{}-{}'.format(outname, ismp))
    return frame[ismp, mat].to_file(output)


def to_file(pert_endf6):
    """
    Function to write Endf6 files with perturbed data in ASCII format.

    Parameters
    ----------
    pert_endf6 : `dict`
        Dictionary containing the `sandy.Endf6` objects with the perturbed data
        for the number of samples entered by the user. The dictionary is
        divided as follows:

                {mat: {Number of the sample: `sandy.Endf6`}}

    """
    global init
    frame = pd.DataFrame(pert_endf6)
    frame.index.name = "SMP"
    frame.columns.name = "MAT"
    frame = frame.stack()
    outname = init.outname if init.outname else os.path.split(init.file)[1]
    proc = 1 if platform.system() == "Windows" else init.processes
    for ismp, dfsmp in frame.groupby("SMP"):
        for mat, dfmat in dfsmp.groupby("MAT"):
            if proc == 1:
                _to_file(frame, ismp, mat, outname)
            else:
                mp.pool.apply_async(_to_file, (frame, ismp, mat, outname))
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

    """
    global init, ftape
    # Command line information and endf6 file information:
    init = parse(iargs)
    ftape = sandy.Endf6.from_file(init.file)
    # Check if the covariance information is available:
    mf_check = ftape.to_series().index.get_level_values("MF")\
                    .intersection(pd.Index(init.mf))
    if len(mf_check) == 0:
        print("Covariance information not available in the file")
        return
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
    pert_endf6 = perturbation_manager(samples, ftape)
    to_file(pert_endf6)
    return pert_endf6, samples


def run():
    t0 = time.time()
    try:
        sampling()
    except SandyError as exc:
        logging.error(exc.args[0])
    print("Total running time: {:.2f} sec".format(time.time() - t0))


if __name__ == "__main__":
    run()
