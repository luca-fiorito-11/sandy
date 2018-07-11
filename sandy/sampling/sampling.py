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



@pytest.mark.sampling
@pytest.mark.xs
def test_sample_xs():
    # Test utils.Xs methods "get_samples", "perturb" and "update"
    from ..data_test import H1
    from ..formats import Errorr, Endf6
    from random import randint
    errtape = Errorr.from_text("\n".join(H1.errorr))
    nsmp = 1000
    perts = errtape.get_xs_cov(listmt=[102]).get_samples(nsmp, eig=10)
    pendftape = Endf6.from_text("\n".join(H1.pendf))
    xs = pendftape.get_xs()
    ismp = randint(1, nsmp);
    pert = perts[ismp]
    pxs = xs.perturb(pert)
    newtape = pendftape.update_xs(pxs)
    assert perts.shape[1] == nsmp
    mat = 125; mt = 102
    ratio = pxs/xs.values
    mask1 = np.in1d(ratio[mat,mt].index, pert[mat,mt].index)
    mask2 = np.in1d(pert[mat,mt].index, ratio[mat,mt].index)
    assert np.isclose(ratio[mat,mt].values[mask1], pert[mat,mt].values[mask2]).all()
    assert newtape.loc[125,3,102].TEXT != pendftape.loc[125,3,102].TEXT
    assert newtape.loc[125,3,2].TEXT == pendftape.loc[125,3,2].TEXT

@pytest.mark.sampling
@pytest.mark.errorr
@pytest.mark.xs
def test_H1(tmpdir):
    from ..data_test import H1
    iargs = [os.path.join(H1.__path__[0], r"h1.pendf"),
             "--errorr-cov", os.path.join(H1.__path__[0], r"h1.errorr"),
             "--outdir", str(tmpdir),
             "--processes", str(os.cpu_count()),
             "--eig", "10",
             "--samples", "100",]
#             "--plotdir", os.path.join(str(tmpdir), r"html_files"),
#             "-p"]
    sampling(iargs)

@pytest.mark.sampling
@pytest.mark.endf6
@pytest.mark.nubar
def test_Cm242(tmpdir):
    from ..data_test import Cm242
    iargs = [os.path.join(Cm242.__path__[0], r"cm242.endf"),
             "--endf6-cov", os.path.join(Cm242.__path__[0], r"cm242.endf"),
             "--outdir", str(tmpdir),
             "--processes", str(os.cpu_count()),
             "--eig", "10",
             "--samples", "100",]
#             "--plotdir", os.path.join(str(tmpdir), r"html_files"),
#             "-p"]
    sampling(iargs)

@pytest.mark.sampling
@pytest.mark.errorr
@pytest.mark.xs
def test_Fe56_errorr(tmpdir):
    from ..data_test import Fe56
    iargs = [os.path.join(Fe56.__path__[0], r"fe56.pendf"),
             "--errorr-cov", os.path.join(Fe56.__path__[0], r"fe56.errorr"),
             "--outdir", str(tmpdir),
             "--processes", str(os.cpu_count()),
             "--eig", "10",
             "--samples", "10",]
    sampling(iargs)

@pytest.mark.sampling
@pytest.mark.errorr
@pytest.mark.xs
@pytest.mark.slow
def test_U5_errorr(tmpdir):
    from ..data_test import U5
    iargs = [os.path.join(U5.__path__[0], r"u235.pendf"),
             "--errorr-cov", os.path.join(U5.__path__[0], r"u235.errorr"),
             "--outdir", str(tmpdir),
             "--processes", str(os.cpu_count()) if os.cpu_count() < 10 else str(10),
             "--eig", "10",
             "--samples", "100",]
#             "--plotdir", os.path.join(str(tmpdir), r"html_files"),
#             "-p"]
    sampling(iargs)

@pytest.mark.sampling
@pytest.mark.chi
def test_U5_chi(tmpdir):
    from ..data_test import U5
    iargs = [os.path.join(U5.__path__[0], r"u235.endf"),
             "--endf6-cov", os.path.join(U5.__path__[0], r"u235.endf"),
             "--outdir", str(tmpdir),
             "--processes", str(os.cpu_count()),
             "--eig", "10",
             "--samples", "100",]
#             "--plotdir", os.path.join(str(tmpdir), r"html_files"),
#             "-p"]
    sampling(iargs)

@pytest.mark.sampling
@pytest.mark.lpc
def test_Fe56_lpc(tmpdir):
    from ..data_test import Fe56
    iargs = [os.path.join(Fe56.__path__[0], r"fe56.endf"),
             "--endf6-cov", os.path.join(Fe56.__path__[0], r"fe56.endf"),
             "--outdir", str(tmpdir),
             "--processes", str(os.cpu_count()),
             "--eig", "10",
             "--samples", "100",]
    sampling(iargs)

@pytest.mark.sampling
@pytest.mark.lpc
def test_U238_lpc(tmpdir):
    from ..data_test import U8
    iargs = [os.path.join(U8.__path__[0], r"u238.endf"),
             "--endf6-cov", os.path.join(U8.__path__[0], r"u238.endf"),
             "--outdir", str(tmpdir),
             "-e", "1e-5", "1e-4", "1e-3", "1e-2", "1e-1", "1e0", "1e1", "1e2", "1e3", "1e4", "1e5",
             "--verbose",
             "--processes", "1",
             "--eig", "10",
             "--mf", "34",
             "--samples", "10",]
    sampling(iargs)