# -*- coding: utf-8 -*-
"""
Created on Mon Mar 19 22:51:03 2018

@author: lucaf
"""

import pandas as pd
from .. import settings
from ..formats import Endf6
import numpy as np
import sys, os, time, platform, pdb, pytest, logging
import multiprocessing as mp


def sampling_mp(ismp):
    global tape, PertXs, PertEdistr, PertLpc, init
    t0 = time.time()
    mat = tape.index.get_level_values("MAT")[0]
    info = tape.read_section(mat, 1, 451)
    lrp = info["LRP"]
    name = info["TAG"]
    newtape = Endf6(tape.copy())

    if not PertXs.empty:
        if lrp == 2:
            xs = newtape.get_xs()
            if not xs.empty:
                xs = xs.perturb(PertXs[ismp])
                newtape = newtape.update_xs(xs)
        else:
            nubar = newtape.get_nubar()
            if not nubar.empty:
                nubar = nubar.perturb(PertXs[ismp])
                newtape = newtape.update_nubar(nubar)

    if not PertEdistr.empty:
        edistr = newtape.get_edistr()
        if not edistr.empty:
            edistr = edistr.add_points(init.energy_points).perturb(PertEdistr[ismp])
            newtape = newtape.update_edistr(edistr)

    if not PertLpc.empty:
        lpc = newtape.get_lpc()
        if not lpc.empty:
            lpc = lpc.add_points(init.energy_points).perturb(PertLpc[ismp], verbose=init.verbose)
            newtape = newtape.update_lpc(lpc)

    print("Created sample {} for {} in {:.2f} sec".format(ismp, name, time.time()-t0,))
    return newtape.update_info().write_string()

def sampling(iargs=None):
    from ..formats import Endf6, Errorr
    from . import from_cli
    t0 = time.time()

    global init
    init = from_cli(iargs)

    # LOAD DATA FILE
    ftape = Endf6.from_file(init.file)
    if ftape.empty: sys.exit("ERROR: tape is empty")
    ftape.parse()

    # LOAD COVARIANCE FILE
    if init.errorr_cov:
        covtape = Errorr.from_file(init.errorr_cov)
    elif init.endf6_cov:
        covtape = Endf6.from_file(init.endf6_cov)
    if covtape.empty: sys.exit("ERROR: covtape is empty")

    # EXTRACT PERTURBATIONS FROM XS/NUBAR COV FILE
    global PertXs
    PertXs = pd.DataFrame()
    if 33 in init.mf or 31 in init.mf:
        xscov = covtape.get_xs_cov(listmt=init.mt, listmat=init.mat)
        if not xscov.empty:
            PertXs = xscov.get_samples(init.samples, eig=init.eig)

    # EXTRACT PERTURBATIONS FROM EDISTR COV FILE
    global PertEdistr
    PertEdistr = pd.DataFrame()
    if 35 in init.mf:
        edistrcov = covtape.get_edistr_cov()
        if not edistrcov.empty:
            PertEdistr = edistrcov.get_samples(init.samples, eig=init.eig)

    # EXTRACT PERTURBATIONS FROM LPC COV FILE
    global PertLpc
    PertLpc = pd.DataFrame()
    if 34 in init.mf:
        lpccov = covtape.get_lpc_cov()
        if not lpccov.empty:
            PertLpc = lpccov.get_samples(init.samples, eig=init.eig)

    if PertLpc.empty and PertEdistr.empty and PertXs.empty:
        print("no covariance section was selected/found")
        return

    # APPLY PERTURBATIONS BY MAT
    global tape
    for mat, tape in ftape.groupby('MAT'):
        tape = Endf6(tape)
#        name = tape.read_section(mat, 1, 451)["TAG"]

        if init.processes == 1:
            outs = {i : sampling_mp(i) for i in range(1,init.samples+1)}
        else:
            if platform.system() == "Windows":
                def init_pool(the_tape):
                    global tape
                    tape = the_tape
                pool = mp.Pool(processes=settings.args.processes,
                               initializer=init_pool(tape))
            else:
                pool = mp.Pool(processes=init.processes)
            outs = {i : pool.apply_async(sampling_mp, args=(i,)) for i in range(1,init.samples+1)}
            outs = {i : out.get() for i,out in outs.items()}

        # DUMP TO FILES
        for ismp, string in outs.items():
            outname = init.outname if init.outname else os.path.split(init.file)[1]
            output = os.path.join(init.outdir, '{}-{}'.format(outname, ismp))
            with open(output, 'w') as f: f.write(string)

    # PLOTTING IS OPTIONAL
#    if kwargs["p"]:
#        plotter.run(iargs)
    print("Total running time 'sampling': {:.2f} sec".format(time.time() - t0))