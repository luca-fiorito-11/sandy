# -*- coding: utf-8 -*-
"""
Created on Mon Mar 19 22:51:03 2018

@author: lucaf
"""

import pandas as pd
from .. import settings
import numpy as np
import sys, os, time, platform, pdb, pytest
import multiprocessing as mp


def sample_chi(tape, NSMP, **kwargs):
    # perturbations are in absolute values
    DfCov = e6.extract_cov35(tape)
    if DfCov.empty:
        return pd.DataFrame()
    cov = Cov(DfCov.values)
    DfPert = pd.DataFrame(cov.sampling(NSMP),
                          index = DfCov.index,
                          columns = range(1,NSMP+1))
    DfPert.columns.name = 'SMP'
    if "eig" in kwargs:
        if kwargs["eig"] > 0:
            from sandy.functions import div0
            eigs = cov.eig()[0]
            idxs = np.abs(eigs).argsort()[::-1]
            dim = min(len(eigs), kwargs["eig"])
            eigs_smp = Cov(np.cov(DfPert.values)).eig()[0]
            idxs_smp = np.abs(eigs_smp).argsort()[::-1]
            print("MF35 eigenvalues:\n{:^10}{:^10}{:^10}".format("EVAL", "SAMPLES","DIFF %"))
            diff = div0(eigs[idxs]-eigs_smp[idxs_smp], eigs[idxs], value=np.NaN)*100.
            E = ["{:^10.2E}{:^10.2E}{:^10.1F}".format(a,b,c) for a,b,c in zip(eigs[idxs][:dim], eigs_smp[idxs_smp][:dim], diff[:dim])]
            print("\n".join(E))
    return DfPert


def sample_xs(tape, NSMP, **kwargs):
    # perturbations are in relative values
    DfCov = tape.get_cov()
    if DfCov.empty:
        return pd.DataFrame()
    cov = Cov(DfCov.values)
    DfPert = pd.DataFrame( cov.sampling(NSMP) + 1, index=DfCov.index, columns=range(1,NSMP+1))
    DfPert.columns.name = 'SMP'
    if "eig" in kwargs:
        if kwargs["eig"] > 0:
            from sandy.functions import div0
            eigs = cov.eig()[0]
            idxs = np.abs(eigs).argsort()[::-1]
            dim = min(len(eigs), kwargs["eig"])
            eigs_smp = Cov(np.cov(DfPert.values)).eig()[0]
            idxs_smp = np.abs(eigs_smp).argsort()[::-1]
            print("MF[31,33] eigenvalues:\n{:^10}{:^10}{:^10}".format("EVAL", "SAMPLES","DIFF %"))
            diff = div0(eigs[idxs]-eigs_smp[idxs_smp], eigs[idxs], value=np.NaN)*100.
            E = ["{:^10.2E}{:^10.2E}{:^10.1F}".format(a,b,c) for a,b,c in zip(eigs[idxs][:dim], eigs_smp[idxs_smp][:dim], diff[:dim])]
            print("\n".join(E))
    return DfPert


def perturb_chi(tape, PertSeriesChi, **kwargs):
    PertSeriesChi.name = 'SMP'
    Chi = e6.extract_chi(tape)
    matListPert = PertSeriesChi.index.get_level_values("MAT")
    for mat in set(Chi.index.get_level_values("MAT")):
        if mat not in matListPert:
            continue
        mtListPert = PertSeriesChi.loc[mat].index.get_level_values("MT")
        for mt in set(Chi.index.get_level_values("MT")):
            if mt not in mtListPert:
                continue
            Pert = pd.DataFrame(PertSeriesChi).query('MAT=={} & MT=={}'.format(mat, mt)).reset_index()
            try:
                Pert = Pert.pivot(index='EINlo', columns='EOUTlo', values='SMP').ffill(axis='columns')
            except:
                aaa=1
            for k,chi in Chi.loc[mat,mt]['CHI'].iteritems():
                # Add extra energy points
                if "energy_point" in kwargs:
                    chi = chi.reindex(chi.index.union(kwargs["energy_point"])).interpolate(method="slinear").fillna(0)
                P = e6.pandas_interpolate(Pert,
                                          np.unique(chi.index).tolist(),
                                          method='zero',
                                          axis='rows')
                P = e6.pandas_interpolate(P,
                                          np.unique(chi.columns).tolist(),
                                          method='zero',
                                          axis='columns')
                P.index.name = "EIN"
                P.columns.name = "EOUT"
                PertChi = chi.add(P, fill_value=0).applymap(lambda x: x if x >= 0 else 0)
                E = PertChi.columns
                M = PertChi.values
                intergral_array = ((M[:,1:]+M[:,:-1])/2).dot(E[1:]-E[:-1])
                SmpChi = pd.DataFrame( M/intergral_array.reshape(-1,1),
                                       index=PertChi.index,
                                       columns=PertChi.columns)
                Chi.loc[mat,mt,k]['CHI'] = SmpChi
    return e6.write_mf5_mt( e6.update_chi(tape, Chi) )

def samplingold(tape, ismp, PertSeriesNubar=None, PertSeriesRes=None, PertSeriesXs=None, PertSeriesChi=None, **kwargs):
    t0 = time.time()
    if PertSeriesXs is not None:
        ptape = kwargs["pendf"].query("MF==3 | (MF==2 & MT==152)")
    if PertSeriesNubar is not None:
        tape = perturb_xs(tape, PertSeriesNubar, **kwargs)
    if PertSeriesChi is not None:
        tape = perturb_chi(tape, PertSeriesChi, **kwargs)
    if PertSeriesRes is not None:
        tape = perturb_res(tape, PertSeriesRes, **kwargs)
        if PertSeriesXs is not None:
            tmpdir = os.path.join(kwargs["outdir"], "tmp-{}".format(ismp))
            os.makedirs(tmpdir, exist_ok=True)
            output = os.path.join(tmpdir, kwargs["endf6"])
            e6.write_tape(tape, output)
            mat = tape.index.get_level_values("MAT")[0]
            pendf = get_pendf(output, kwargs['njoy'], mat=mat, wd=tmpdir)
            ptape = e6.endf2df(pendf).query("MF==3 | (MF==2 & MT==152)")
    if PertSeriesXs is not None:
        tape = tape.query("MF!=3").append(ptape)
        tape = perturb_xs(tape, PertSeriesXs, **kwargs)
    output = os.path.join(kwargs["outdir"], os.path.basename(kwargs["endf6"]) + '-{}'.format(ismp))
    e6.write_tape(tape, output)
    print("Created file '{}' in {:.2f} sec".format(output, time.time()-t0,))
#    if 33 in kwargs['keep_cov_mf']:
#        if kwargs['pendf'] is not None:
#            ptape = kwargs['pendf'].query("MF==3 | (MF==2 & MT==152)")
#        else:
#            mat = tape.index.get_level_values("MAT")[0]
#            pendf = get_pendf(kwargs['endf6'], kwargs['njoy'], mat=mat, wd=kwargs['outdir'])
#            ptape = e6.endf2df(pendf).query("MF==3 | (MF==2 & MT==152)")
#        tape = tape.query("MF!=3").append(ptape)
#        for mat in tape.index.get_level_values('MAT').unique():
#            tape.DATA.loc[mat,1,451]['LRP'] == 2
#    tape = perturb_xs(tape, PertSeriesXs, **kwargs)

def sampling_sp(ismp, PertSeriesXs, **kwargs):
    global tape
    t0 = time.time()
    tapeout = tape.get_xs().perturb(PertSeriesXs).update_tape(tape).write_mf1_nubar().write_mf3_mt()
    output = os.path.join(kwargs["outdir"], os.path.basename(kwargs["file"]) + '-{}'.format(ismp))
    string = tapeout.to_string()
    print("Created file '{}' in {:.2f} sec".format(output, time.time()-t0,))
    return string, output


def run(iargs=None):
    from sandy.sampling import plotter
    t0 = time.time()

    # SETUP OF SETTINGS
    settings.init_sampling(iargs)
    kwargs = vars(settings.args)

    # LOAD DATA FILE
    global tape
    tape = e6.Endf6.from_file(kwargs["file"]).process()
    if tape.empty:
        sys.exit("ERROR: tape is empty")

    # LOAD COVARIANCE FILE
    if kwargs["errorr_cov"]:
        covtape = Errorr.from_file(settings.args.errorr_cov).process()
    elif kwargs["endf6_cov"]:
        covtape = e6.Endf6.from_file(settings.args.endf6_cov).process()
    if covtape.empty:
        sys.exit("ERROR: covtape is empty")


    # EXTRACT PERTURBATIONS FROM COV FILE
    PertXs = sample_xs(covtape, settings.args.samples, **kwargs)
    PertChi = sample_chi(covtape, settings.args.samples, **kwargs)

    # APPLY PERTURBATIONS
    if settings.args.processes == 1:
        outs = [sampling_sp(i, PertXs[i], **kwargs) for i in range(1,settings.args.samples+1)]
    else:
        if platform.system() == "Windows":
            def init_pool(the_tape):
                global tape
                tape = the_tape
            pool = mp.Pool(processes=settings.args.processes,
                           initializer=init_pool(tape))
        else:
            pool = mp.Pool(processes=settings.args.processes)
        outs = [pool.apply_async(sampling_sp,
                                 args = (i, PertXs[i]),
                                 kwds = {**kwargs}
                                 ) for i in range(1,settings.args.samples+1) ]
        outs = list(map(lambda x:x.get(), outs))

    # DUMP TO FILES
    for string, output in outs:
        with open(output, 'w') as f:
            f.write(string)

    # PLOTTING IS OPTIONAL
    if kwargs["p"]:
        plotter.run(iargs)
    print("Total running time 'sampling': {:.2f} sec".format(time.time() - t0))

def sampling_mp(ismp):
    global tape, PertXs
    t0 = time.time()
    mat = tape.index.get_level_values("MAT")[0]
    info = tape.read_section(mat, 1, 451)
    lrp = info["LRP"]
    name = info["TAG"]
    if lrp == 2:
        try:
            xs = tape.get_xs().perturb(PertXs[ismp])
            tape = tape.update_xs(xs)
        except:
            pass
    else:
        try:
            nubar = tape.get_nubar().perturb(PertXs[ismp])
            tape = tape.update_nubar(nubar)
        except:
            pass
    print("Created sample {} for {} in {:.2f} sec".format(ismp, name, time.time()-t0,))
    return tape.update_info().write_string()

def sampling(iargs=None):
    from ..formats import Endf6, Errorr
    t0 = time.time()
    init = settings.init_sampling(iargs)

    # LOAD DATA FILE
    ftape = Endf6.from_file(init.file)
    if ftape.empty: sys.exit("ERROR: tape is empty")

    # LOAD COVARIANCE FILE
    if init.errorr_cov:
        covtape = Errorr.from_file(init.errorr_cov)
    elif init.endf6_cov:
        covtape = Endf6.from_file(init.endf6_cov)
    if covtape.empty: sys.exit("ERROR: covtape is empty")


    # EXTRACT PERTURBATIONS FROM XS/NUBAR COV FILE
    global PertXs
    try:
        PertXs = covtape.get_xs_cov().get_samples(init.samples, eig=init.eig)
    except:
        PertXs = pd.DataFrame()
    if PertXs.empty: # Add checks for other data (MF35, MF34, ...)
        sys.exit("ERROR: no covariance matrix was found")

    # APPLY PERTURBATIONS BY MAT
    global tape
    for mat, tape in ftape.groupby('MAT'):
        tape = Endf6(tape)
        name = tape.read_section(mat, 1, 451)["TAG"]

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
            output = os.path.join(init.outdir, '{}-{}-{}'.format(name, int(mat), ismp))
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