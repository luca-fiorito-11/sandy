# -*- coding: utf-8 -*-
"""
Created on Mon Mar 19 22:51:03 2018

@author: Luca Fiorito
"""

import os
import time
import pdb
import argparse
import logging
import numpy as np
import multiprocessing as mp

import pandas as pd

from sandy.settings import SandyError
from sandy.formats import read_formatted_file, get_file_format
from sandy.formats.endf6 import Endf6
from sandy.formats.utils import FySamples, XsCov
from sandy.core import pfns
from sandy.utils import is_valid_dir, is_valid_file

__author__ = "Luca Fiorito"
__all__ = ["sampling"]


def _sampling_mp(ismp, skip_title=False, skip_fend=False):
    global tape, PertXs, PertNubar, PertEdistr, PertLpc, init
    t0 = time.time()
    mat = tape.index.get_level_values("MAT")[0]
    info = tape.read_section(mat, 1, 451)
    lrp = info["LRP"]
    name = info["TAG"]
    newtape = Endf6(tape.copy())
    extra_points = np.logspace(-5, 7, init.energy_sequence)
    if not PertXs.empty and lrp == 2:
        xs = newtape.get_xs()
        if not xs.empty:
            xs = xs.perturb(PertXs[ismp])
            newtape = newtape.update_xs(xs)
    if not PertNubar.empty:
        nubar = newtape.get_nubar()
        if not nubar.empty:
            nubar = nubar.perturb(PertNubar[ismp])
            newtape = newtape.update_nubar(nubar)
    if not PertEdistr.empty:
        edistr = pfns.from_endf6(newtape)
        if not edistr.empty:
            edistr = edistr.add_points(extra_points).perturb(PertEdistr[ismp])
            newtape = newtape.update_edistr(edistr)
    if not PertLpc.empty:
        lpc = newtape.get_lpc()
        if not lpc.empty:
            lpc = lpc.add_points(extra_points).perturb(PertLpc[ismp])
            newtape = newtape.update_lpc(lpc)
    print("Created sample {} for MAT {} in {:.2f} sec".format(ismp, mat, time.time()-t0,))
    descr = ["perturbed file No.{} created by SANDY".format(ismp)]
    return newtape.delete_cov().update_info(descr=descr).write_string(skip_title=skip_title, skip_fend=skip_fend)



def _sampling_fy_mp(ismp, skip_title=False, skip_fend=False):
    global tape, PertFY, init
    t0 = time.time()
    newtape = Endf6(tape.copy())
    mat = tape.index.get_level_values("MAT")[0]
    fy = newtape.get_fy()
    fynew = fy.perturb(PertFy[ismp])
    newtape = newtape.update_fy(fynew)
    print("Created sample {} for MAT {} in {:.2f} sec".format(ismp, mat, time.time()-t0,))
    descr = ["perturbed file No.{} created by SANDY".format(ismp)]
    return newtape.delete_cov().update_info(descr=descr).write_string(skip_title=skip_title, skip_fend=skip_fend)



def _parse(iargs=None):
    """Parse command line arguments for sampling option.
    """
    parser = argparse.ArgumentParser(prog="python -m sandy.sampling",
                                     description='Run sampling',
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('file',
                        type=lambda x: is_valid_file(parser, x),
                        help="ENDF-6 or PENDF format file")
    parser.add_argument('--cov', '-C',
                       type=lambda x: is_valid_file(parser, x),
                       help="file containing covariances")
    parser.add_argument('--samples', '-S',
                        type=int,
                        default=200,
                        help="number of samples (default = 200)")
    parser.add_argument('--outdir', '-D',
                        metavar="DIR",
                        default=os.getcwd(),
                        type=lambda x: is_valid_dir(parser, x, mkdir=True),
                        help="target directory where outputs are stored\n(default = current working directory)\nif it does not exist it will be created")
    parser.add_argument('--processes','-N',
                        type=int,
                        default=1,
                        help="number of worker processes (default = 1)")
    parser.add_argument('--max-polynomial','-P',
                        type=int,
                        help="Maximum order of Legendre polynomial coefficients considered for sampling (default = all)")
    parser.add_argument('--eig',
                        type=int,
                        default=10,
                        metavar="N",
                        help="print the first N eigenvalues of the evaluated covariance matrices\n(default = do not print)")
    parser.add_argument('--mat',
                        type=int,
                        default=range(0,10000),
                        action='store',
                        nargs="+",
                        metavar="{1,..,9999}",
                        help="draw samples only from the selected MAT sections (default = keep all)")
    parser.add_argument('--mf',
                        type=int,
                        default=[31,33,34,35],
                        action='store',
                        nargs="+",
                        metavar="{31,33,34,35}",
                        help="draw samples only from the selected MF sections (default = keep all)")
    parser.add_argument('--mt',
                        type=int,
                        action='store',
                        nargs="+",
                        metavar="{1,..,999}",
                        help="draw samples only from the selected MT sections (default = keep all)")
    parser.add_argument('--outname','-O',
                        type=str,
                        help="basename for the output files (default is the the basename of <file>.)")
#    parser.add_argument('--verbose',
#                        default=False,
#                        action="store_true",
#                        help="turn on verbosity (default = quiet)")
    parser.add_argument('--debug',
                        default=False,
                        action="store_true",
                        help="turn on debug mode")
    parser.add_argument('--fission-yields','-F',
                        default=False,
                        action="store_true",
                        help="input <file> contains fission yields")
    parser.add_argument('--energy-sequence','-E',
                        type=int,
                        metavar="EL",
                        default=49,
                        help=argparse.SUPPRESS)
#    parser.add_argument("-v",
#                        '--version',
#                        action='version',
#                        version='%(prog)s 1.0',
#                        help="SANDY's version.")
    return parser.parse_known_args(args=iargs)[0]



def sampling(iargs=None):
    """Construct multivariate normal distributions with a unit vector for 
    mean and with relative covariances taken from the evaluated files.
    Perturbation factors are sampled with the same multigroup structure of 
    the covariance matrix, and are applied to the pointwise data to produce 
    the perturbed files.
    """
    global init, PertXs, PertNubar, PertEdistr, PertLpc, PertFy, tape
    init = _parse(iargs)
    # LOAD ENDF6 AND COVARIANCE FILE
    ftape = read_formatted_file(init.file)
    if init.cov:
        covtape = read_formatted_file(init.cov)
        covtype = get_file_format(init.cov)
    else:
        covtape = ftape
        covtype = get_file_format(init.file)
    df = {}
    if init.fission_yields:
        # EXTRACT FY PERTURBATIONS FROM COV FILE
        fy = ftape.get_fy(listmat=init.mat, listmt=init.mt)
        if fy.empty:
            logging.warn("no fission yield section was selected/found")
            return
        index = fy.index.to_frame(index=False)
        dfperts = []
        for mat,dfmat in index.groupby("MAT"):
            for mt,dfmt in dfmat.groupby("MT"):
                for e,dfe in dfmt.groupby("E"):
                    fycov = fy.get_cov(mat, mt, e)
                    pert = fycov.get_samples(init.samples, eig=0)
                    dfperts.append(pert)
        PertFy = FySamples(pd.concat(dfperts))
        if init.debug: PertFy.to_csv("perts_mf8.csv")
        # DELETE LOCAL VARIABLES
        for k in locals().keys():
            del locals()[k]
        # APPLY PERTURBATIONS BY MAT
        for imat,(mat, tape) in enumerate(sorted(ftape.groupby('MAT'))):
            skip_title = False if imat == 0 else True
            skip_fend = False if imat == ftape.index.get_level_values("MAT").unique().size -1 else True
            tape = Endf6(tape)
            kw = dict(skip_title=skip_title, skip_fend=skip_fend)
            if mat not in init.mat:
                out = tape.write_string(**kw)
                outs = {i : out for i in range(1,init.samples+1)}
            else:
                if init.processes == 1:
                    outs = {i : _sampling_fy_mp(i, **kw) for i in range(1,init.samples+1)}
                else:
                    pool = mp.Pool(processes=init.processes)
                    outs = {i : pool.apply_async(_sampling_fy_mp, (i,), kw) for i in range(1,init.samples+1)}
                    outs = {i : out.get() for i,out in outs.items()}
                    pool.close()
                    pool.join()
            df.update({ mat : outs })
    else:
        # EXTRACT XS PERTURBATIONS FROM COV FILE
        PertXs = pd.DataFrame()
        if 33 in init.mf and 33 in covtape.mf:
            listmt = init.mt if init.mt is None else [451].extend(init.mt) # ERRORR needs MF1/MT451 to get the energy grid
            method = XsCov.from_errorr if covtype is "errorr" else XsCov.from_endf6
            xscov = method(covtape.filter_by(listmt=listmt, listmf=[1,33], listmat=init.mat))
            if not xscov.empty:
                count = xscov.check_diagonal()
                if count != 0:
                    logging.warn("MF33 covariances will not be sampled")
                else:
                    PertXs = xscov.get_samples(init.samples, eig=init.eig)
                    if init.debug: PertLpc.to_csv("perts_mf33.csv")
        # EXTRACT NUBAR PERTURBATIONS FROM ENDF6 FILE
        PertNubar = pd.DataFrame()
        if 31 in init.mf and 31 in covtape.mf:
            listmt = init.mt if init.mt is None else [451].extend(init.mt) # ERRORR needs MF1/MT451 to get the energy grid
            method = XsCov.from_errorr if covtype is "errorr" else XsCov.from_endf6
            nubarcov = method(covtape.filter_by(listmt=listmt, listmf=[1,31], listmat=init.mat))
            if not nubarcov.empty:
                count = nubarcov.check_diagonal()
                if count != 0:
                    logging.warn("MF31 covariances will not be sampled")
                else:
                    PertNubar = nubarcov.get_samples(init.samples, eig=init.eig)
                    if init.debug: PertNubar.to_csv("perts_mf31.csv")
        # EXTRACT PERTURBATIONS FROM EDISTR COV FILE
        PertEdistr = pd.DataFrame()
        if 35 in init.mf and 35 in covtape.mf:
            edistrcov = ftape.get_edistr_cov()
            if not edistrcov.empty:
                count = edistrcov.check_diagonal()
                if count != 0:
                    logging.warn("MF35 covariances will not be sampled")
                else:
                    PertEdistr = edistrcov.get_samples(init.samples, eig=init.eig)
                    if init.debug: PertEdistr.to_csv("perts_mf35.csv")
        # EXTRACT PERTURBATIONS FROM LPC COV FILE
        PertLpc = pd.DataFrame()
        if 34 in init.mf and 34 in covtape.mf:
            lpccov = ftape.get_lpc_cov()
            if not lpccov.empty:
                if init.max_polynomial:
                    lpccov = lpccov.filter_p(init.max_polynomial)
                count = lpccov.check_diagonal()
                if count != 0:
                    logging.warn("MF34 covariances will not be sampled")
                else:
                    PertLpc = lpccov.get_samples(init.samples, eig=init.eig)
                    if init.debug: PertLpc.to_csv("perts_mf34.csv")
        if PertLpc.empty and PertEdistr.empty and PertXs.empty and PertNubar.empty:
            logging.warn("no covariance section was selected/found")
            return
        # DELETE LOCAL VARIABLES
        for k in locals().keys():
            del locals()[k]
        # APPLY PERTURBATIONS BY MAT
        for imat,(mat, tape) in enumerate(sorted(ftape.groupby('MAT'))):
            skip_title = False if imat == 0 else True
            skip_fend = False if imat == ftape.index.get_level_values("MAT").unique().size -1 else True
            tape = Endf6(tape)
            kw = dict(skip_title=skip_title, skip_fend=skip_fend)
            info = tape.read_section(mat, 1, 451)
            if info["LRP"] != 2 and 33 in init.mf:
                logging.warn("not a PENDF tape, cross sections will not be perturbed")
            if init.processes == 1:
                outs = {i : _sampling_mp(i, **kw) for i in range(1,init.samples+1)}
            else:
                pool = mp.Pool(processes=init.processes)
                outs = {i : pool.apply_async(_sampling_mp, (i,), kw) for i in range(1,init.samples+1)}
                outs = {i : out.get() for i,out in outs.items()}
                pool.close()
                pool.join()
            df.update({ mat : outs })
    # DUMP TO FILES
    frame = pd.DataFrame(df)
    frame.index.name = "SMP"
    frame.columns.name = "MAT"
    frame = frame.stack()
    outname = init.outname if init.outname else os.path.split(init.file)[1]
    for ismp,dfsmp in frame.groupby("SMP"):
        output = os.path.join(init.outdir, '{}-{}'.format(outname, ismp))
        with open(output, 'w') as f:
            for mat,dfmat in dfsmp.groupby("MAT"):
                f.write(frame[ismp,mat])


def run():
    t0 = time.time()
    try:
        sampling()
    except SandyError as exc:
        logging.error(exc.args[0])
    print("Total running time: {:.2f} sec".format(time.time() - t0))


if __name__ == "__main__":
    run()