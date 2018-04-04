# -*- coding: utf-8 -*-
"""
Created on Mon Mar 19 22:51:03 2018

@author: lucaf
"""

import pandas as pd
from sandy.endf6 import files as e6
from sandy import settings
from sandy.sampling.cov import Cov
import numpy as np
import sys
import os
import multiprocessing as mp
from copy import deepcopy
import shutil
import time

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
    DfCov = e6.extract_cov33(tape)
    if DfCov.empty:
        return pd.DataFrame()
    cov = Cov(DfCov.as_matrix())
    DfPert = pd.DataFrame( cov.sampling(NSMP) + 1,
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
            print("MF[31,33] eigenvalues:\n{:^10}{:^10}{:^10}".format("EVAL", "SAMPLES","DIFF %"))
            diff = div0(eigs[idxs]-eigs_smp[idxs_smp], eigs[idxs], value=np.NaN)*100.
            E = ["{:^10.2E}{:^10.2E}{:^10.1F}".format(a,b,c) for a,b,c in zip(eigs[idxs][:dim], eigs_smp[idxs_smp][:dim], diff[:dim])]
            print("\n".join(E))
    return DfPert

def perturb_xs(tape, PertSeriesXs, **kwargs):
    Xs = e6.extract_xs(tape)
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
            Pert = Pert.pivot(index='EINlo', columns='EOUTlo', values='SMP').ffill(axis='columns')
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

def sampling(tape, output, PertSeriesXs=None, PertSeriesChi=None, ismp=None, **kwargs):
    t0 = time.time()
    itape = deepcopy(tape)
    if PertSeriesChi is not None:
        itape = perturb_chi(itape, PertSeriesChi, **kwargs)
    if PertSeriesXs is not None:
        itape = perturb_xs(itape, PertSeriesXs, **kwargs)
    e6.write_tape(itape, output)
    if ismp is not None:
        print("Created sample {} in file '{}' in {:.2f} sec".format(ismp, output, time.time()-t0,))
    else:
        print("Created file '{}' in {:.2f} sec".format(output, time.time()-t0,))


def run():
    t0 = time.time()
    if len(sys.argv) == 1:
        from sandy.data_test import __file__ as td
        from sandy import __file__ as sd
        from os.path import join
        sd = os.path.dirname(os.path.realpath(sd))
        td = os.path.dirname(os.path.realpath(td))
        sys.argv.extend([join(td, r"94-Pu-239g.jeff33"),
#                         "--pendf", join(td, r"H1.txt.pendf"),
                         "--outdir", r"tmp-dir",
                         "--njoy", join(sd, r"njoy2012_50.exe"),
                         "--eig", "10",
                         "--samples", "100",
                         "--processes", "1",
                         "-e", "1e-5",
                         "-e", "5e-5",
                         "-e", "1e-4",
                         "-e", "5e-4",
                         "-e", "1e-3",
                         "-e", "5e-3",
                         "-e", "1e-2",
                         "-e", "5e-2",
                         "-e", "1e-1",
                         "-e", "5e-1",
                         "-e", "1e0",
                         "-e", "5e0",
                         "-e", "1e1",
                         "-e", "5e1",])
    settings.init_sampling()

    tape = e6.endf2df(settings.args.endf6)#, keep_mf=[3], keep_mt=[102])

    if settings.args.keep_mat:
        query = "|".join([ "MAT=={}".format(x) for x in settings.args.keep_mat])
        tape = tape.query(query)
    if settings.args.keep_cov_mf:
        query = "|".join([ "MF=={}".format(x) for x in settings.args.keep_cov_mf])
        tape = tape.query("MF < 31 | ({})".format(query))
    if settings.args.keep_cov_mt:
        query = "|".join([ "MT=={}".format(x) for x in settings.args.keep_cov_mt])
        tape = tape.query("MF < 31 | ({})".format(query))
    if tape.empty:
        sys.exit("ERROR: tape is empty")
    MATS = list(tape.index.get_level_values("MAT").unique())
    MFS = list(tape.index.get_level_values("MF").unique())

    # Always run. If MF35 is not wanted, then MF35 sections are already removed.
    PertXs = sample_xs(tape, settings.args.samples, **vars(settings.args))
    PertChi = sample_chi(tape, settings.args.samples, **vars(settings.args))


    if 31 in MFS or 35 in MFS:
        outname = os.path.join(settings.args.outdir, os.path.basename(settings.args.endf6) + '-{}')
        if settings.args.processes == 1:
            for ismp in range(1,settings.args.samples+1):
                sampling(tape,
                         outname.format(ismp),
                         PertSeriesXs=PertXs[ismp] if not PertXs.empty else None,
                         PertSeriesChi=PertChi[ismp] if not PertChi.empty else None,
                         ismp=ismp)
        else:
            pool = mp.Pool(processes=settings.args.processes)
            [ pool.apply(sampling,
                         args = (tape, outname.format(ismp)),
                         kwds = {**{"PertSeriesXs"  : PertXs[ismp] if not PertXs.empty else None,
                                    "PertSeriesChi" : PertChi[ismp] if not PertChi.empty else None,
                                    "ismp" : ismp},
                                 **vars(settings.args)}
                         ) for ismp in range(1,settings.args.samples+1) ]

    if 33 in MFS:
        if not settings.args.pendf:
            from sandy.njoy import FileNJOY
            fnjoy = FileNJOY()
            # Write NJOY file
            fnjoy.copy_to_tape(settings.args.endf6, 20, dst=settings.args.outdir)
            fnjoy.reconr(20, 21, mat=tape.index.get_level_values('MAT').unique())
            fnjoy.stop()
            fnjoy.run(settings.args.njoy, cwd=settings.args.outdir)
            settings.args.pendf = os.path.join(settings.args.outdir,
                                               '{}.pendf'.format(os.path.basename(settings.args.endf6)))
            shutil.move(os.path.join(settings.args.outdir, r'tape21'), settings.args.pendf)
            # Remove NJOY junk outputs
            os.unlink(os.path.join(settings.args.outdir, r'tape20'))
            os.unlink(os.path.join(settings.args.outdir, r'output'))
        ptape = e6.endf2df(settings.args.pendf)

        outname = os.path.join(settings.args.outdir, os.path.basename(settings.args.pendf) + '-{}')
        if settings.args.processes == 1:
            for ismp in range(1,settings.args.samples+1):
                sampling(ptape,
                         outname.format(ismp),
                         PertSeriesXs=PertXs[ismp] if not PertXs.empty else None,
                         ismp=ismp,
                         **vars(settings.args))
        else:
            pool = mp.Pool(processes=settings.args.processes)
            [ pool.apply(sampling,
                         args = (ptape, outname.format(ismp)),
                         kwds = {**{"PertSeriesXs" : PertXs[ismp] if not PertXs.empty else None,
                                    },
                                 **vars(settings.args)}
                         ) for ismp in range(1,settings.args.samples+1) ]
    print("Total running time: {:.2f} sec".format(time.time() - t0))

if __name__ == '__main__':
    run()
