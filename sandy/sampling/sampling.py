# -*- coding: utf-8 -*-
"""
Created on Mon Mar 19 22:51:03 2018

@author: lucaf
"""

import pandas as pd
import sandy.formats.endf6 as e6
from sandy import settings
from sandy.sampling.cov import Cov
import numpy as np
import sys
import os
import multiprocessing as mp
from copy import deepcopy
import shutil
import time
import matplotlib.pyplot as plt
from sandy.njoy import get_pendf
from sandy.tests import TimeDecorator
from sandy.formats.errorr import Errorr
import re


#To produce correlation matrix
#Index = df_cov_xs.index
#df_cov_xs.update(pd.DataFrame(Cov(df_cov_xs.as_matrix()).corr, index=Index, columns=Index))

def sample_chi(tape, NSMP, **kwargs):
    # perturbations are in absolute values
    DfCov = e6.extract_cov35(tape)
    if DfCov.empty:
        return pd.DataFrame()
    cov = Cov(DfCov.as_matrix())
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
            eigs_smp = Cov(np.cov(DfPert.as_matrix())).eig()[0]
            idxs_smp = np.abs(eigs_smp).argsort()[::-1]
            print("MF35 eigenvalues:\n{:^10}{:^10}{:^10}".format("EVAL", "SAMPLES","DIFF %"))
            diff = div0(eigs[idxs]-eigs_smp[idxs_smp], eigs[idxs], value=np.NaN)*100.
            E = ["{:^10.2E}{:^10.2E}{:^10.1F}".format(a,b,c) for a,b,c in zip(eigs[idxs][:dim], eigs_smp[idxs_smp][:dim], diff[:dim])]
            print("\n".join(E))
    return DfPert


def sample_xs(tape, NSMP, **kwargs):
    # perturbations are in relative values
    from sandy.formats.endf6 import XsCov
    DfCov = XsCov.from_errorr_tape(tape)
    if DfCov.empty:
        return pd.DataFrame()
    cov = Cov(DfCov.as_matrix())
    DfPert = pd.DataFrame( cov.sampling(NSMP) + 1, index=DfCov.index, columns=range(1,NSMP+1))
    DfPert.columns.name = 'SMP'
    if "eig" in kwargs:
        if kwargs["eig"] > 0:
            from sandy.functions import div0
            eigs = cov.eig()[0]
            idxs = np.abs(eigs).argsort()[::-1]
            dim = min(len(eigs), kwargs["eig"])
            eigs_smp = Cov(np.cov(DfPert.as_matrix())).eig()[0]
            idxs_smp = np.abs(eigs_smp).argsort()[::-1]
            print("MF[31,33] eigenvalues:\n{:^10}{:^10}{:^10}".format("EVAL", "SAMPLES","DIFF %"))
            diff = div0(eigs[idxs]-eigs_smp[idxs_smp], eigs[idxs], value=np.NaN)*100.
            E = ["{:^10.2E}{:^10.2E}{:^10.1F}".format(a,b,c) for a,b,c in zip(eigs[idxs][:dim], eigs_smp[idxs_smp][:dim], diff[:dim])]
            print("\n".join(E))
    return DfPert

#def sample_xs(tape, NSMP, **kwargs):
#    # perturbations are in relative values
#    DfCov = e6.extract_cov33(tape)
#    if DfCov.empty:
#        return pd.DataFrame()
#    cov = Cov(DfCov.as_matrix())
#    DfPert = pd.DataFrame( cov.sampling(NSMP) + 1, index=DfCov.index, columns=range(1,NSMP+1))
#    DfPert.columns.name = 'SMP'
#    if "eig" in kwargs:
#        if kwargs["eig"] > 0:
#            from sandy.functions import div0
#            eigs = cov.eig()[0]
#            idxs = np.abs(eigs).argsort()[::-1]
#            dim = min(len(eigs), kwargs["eig"])
#            eigs_smp = Cov(np.cov(DfPert.as_matrix())).eig()[0]
#            idxs_smp = np.abs(eigs_smp).argsort()[::-1]
#            print("MF[31,33] eigenvalues:\n{:^10}{:^10}{:^10}".format("EVAL", "SAMPLES","DIFF %"))
#            diff = div0(eigs[idxs]-eigs_smp[idxs_smp], eigs[idxs], value=np.NaN)*100.
#            E = ["{:^10.2E}{:^10.2E}{:^10.1F}".format(a,b,c) for a,b,c in zip(eigs[idxs][:dim], eigs_smp[idxs_smp][:dim], diff[:dim])]
#            print("\n".join(E))
#    return DfPert

@TimeDecorator
def perturb_xs(tape, PertSeriesXs, **kwargs):
    Xs = e6.Xs.from_tape(tape)
    indexName = Xs.index.name
    # Add extra energy points
    if "energy_point" in kwargs:
        Xs = Xs.reindex(Xs.index.union(kwargs["energy_point"])).interpolate(method="slinear").fillna(0)
    Xs.index.name = indexName
    for mat, mt in Xs:
        if mat not in PertSeriesXs.index.get_level_values("MAT").unique():
            continue
        mtListPert = PertSeriesXs.loc[mat].index.get_level_values("MT").unique()
        if mt in mtListPert:
            mtPert = mt
        elif mt in range(800,850) \
        and not list(filter(lambda x: x in mtListPert, range(800,850))) \
        and 107 in mtListPert:
            mtPert = 107
        elif mt in range(750,800) \
        and not list(filter(lambda x: x in mtListPert, range(750,800))) \
        and 106 in mtListPert:
            mtPert = 106
        elif mt in range(700,750) \
        and not list(filter(lambda x: x in mtListPert, range(700,750))) \
        and 105 in mtListPert:
            mtPert = 105
        elif mt in range(650,700) \
        and not list(filter(lambda x: x in mtListPert, range(650,700))) \
        and 104 in mtListPert:
            mtPert = 104
        elif mt in range(600,650) \
        and not list(filter(lambda x: x in mtListPert, range(600,650))) \
        and 103 in mtListPert:
            mtPert = 103
        elif mt in range(102,118) \
        and not list(filter(lambda x: x in mtListPert, range(102,118))) \
        and 101 in mtListPert:
            mtPert = 101
        elif mt in (19,20,21,38) \
        and not list(filter(lambda x: x in mtListPert, (19,20,21,38))) \
        and 18 in mtListPert:
            mtPert = 18
        elif mt in (18,101) \
        and not list(filter(lambda x: x in mtListPert, (18,101))) \
        and 27 in mtListPert:
            mtPert = 27
        elif mt in range(50,92) \
        and not list(filter(lambda x: x in mtListPert, range(50,92))) \
        and 4 in mtListPert:
            mtPert = 4
        else:
            continue
        P = PertSeriesXs.loc[mat,mtPert]
        P = P.reindex(P.index.union(Xs[mat,mt].index)).ffill().fillna(1).reindex(Xs[mat,mt].index)
        Xs[mat,mt] = Xs[mat,mt].multiply(P, axis="index")
        # Negative values are set to zero
        Xs[mat,mt][Xs[mat,mt] <= 0] = 0
    Xs = e6.reconstruct_xs(Xs)
    return e6.write_mf1_nubar( e6.write_mf3_mt( e6.update_xs(tape, Xs) ) )

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
                M = PertChi.as_matrix()
                intergral_array = ((M[:,1:]+M[:,:-1])/2).dot(E[1:]-E[:-1])
                SmpChi = pd.DataFrame( M/intergral_array.reshape(-1,1),
                                       index=PertChi.index,
                                       columns=PertChi.columns)
                Chi.loc[mat,mt,k]['CHI'] = SmpChi
    return e6.write_mf5_mt( e6.update_chi(tape, Chi) )

def sampling(tape, ismp, PertSeriesNubar=None, PertSeriesRes=None, PertSeriesXs=None, PertSeriesChi=None, **kwargs):
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

def sampling2(ismp, PertSeriesXs, **kwargs):
    global tape
    t0 = time.time()
    tapeout = tape.get_xs().perturb(PertSeriesXs).update_tape(tape)
    tapeout = e6.write_mf1_nubar(tapeout)
    tapeout = e6.write_mf3_mt(tapeout)

    output = os.path.join(kwargs["outdir"], os.path.basename(kwargs["file"]) + '-{}'.format(ismp))
    string = e6.Endf6(tapeout).to_string(output)
    print("Created file '{}' in {:.2f} sec".format(output, time.time()-t0,))
    return string, output


def run(iargs=None):
    t0 = time.time()
    settings.init_sampling(iargs)

    # LOAD DATA FILE
    global tape
    tape = e6.Endf6.from_file(settings.args.file).process()
    if tape.empty:
        sys.exit("ERROR: tape is empty")

    # LOAD COVARIANCE FILE
    if settings.args.errorr_cov:
        covtape = Errorr.from_file(settings.args.errorr_cov).process()#, keep_mf=[3], keep_mt=[102])
    elif settings.args.endf6_cov:
        covtape = e6.Endf6.from_file(settings.args.endf6_cov).process()
    if covtape.empty:
        sys.exit("ERROR: covtape is empty")

    # Further setup of settings
    kwargs = vars(settings.args)

    # EXTRACT PERTURBATIONS FROM COV FILE
    PertXs = sample_xs(covtape, settings.args.samples, **kwargs)
    PertChi = sample_chi(covtape, settings.args.samples, **kwargs)

    # APPLY PERTURBATIONS
    if settings.args.processes == 1:
        outs = [sampling2(i, PertXs[i], **kwargs) for i in range(1,settings.args.samples+1)]
    else:
        pool = mp.Pool(processes=settings.args.processes)
        outs = [pool.apply_async(sampling2,
                                 args = (i, PertXs[i]),
                                 kwds = {**kwargs}
                                 ) for i in range(1,settings.args.samples+1) ]
        outs = list(map(lambda x:x.get(), outs))

    # DUMP TO FILES
    for string, output in outs:
        with open(output, 'w') as f:
            f.write(string)
    print("Total running time: {:.2f} sec".format(time.time() - t0))