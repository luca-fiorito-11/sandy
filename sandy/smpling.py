# -*- coding: utf-8 -*-
"""
Created on Mon Mar 19 22:51:03 2018

@author: lucaf
"""

import pandas as pd
import sandy.endf6.files as e6
from sandy.cov import Cov
import numpy as np
import sandy.settings as settings
import sys
import os
import multiprocessing as mp


#To produce correlation matrix
#Index = df_cov_xs.index
#df_cov_xs.update(pd.DataFrame(Cov(df_cov_xs.as_matrix()).corr, index=Index, columns=Index))

def sample_chi(tape, NSMP):
    # perturbations are in absolute values
    DfCov = e6.extract_cov35(tape)
    DfPert = pd.DataFrame( Cov(DfCov.as_matrix()).sampling(NSMP),
                                 index = DfCov.index,
                                 columns = range(1,NSMP+1))
    DfPert.columns.name = 'SMP'
    return DfPert


def sample_xs(tape, NSMP):
    # perturbations are in relative values
    DfCov = e6.merge_covs( e6.extract_cov33(tape) )
    DfPert = pd.DataFrame( Cov(DfCov.as_matrix()).sampling(NSMP) + 1,
                                 index = DfCov.index,
                                 columns = range(1,NSMP+1))
    DfPert.columns.name = 'SMP'
    return DfPert


def perturb(XsSeries, PertSeries):
    Pert = e6.pandas_interpolate(PertSeries,
                           np.unique(XsSeries.index).tolist(),
                           axis="rows")
    Xs = XsSeries.multiply(Pert, axis="index")
    Xs.name = XsSeries.name # this is mt
    # Negative values are set to zero
    Xs[Xs <= 0] = 0
    return Xs


def sampling(tape, PertSeriesXs, output):
    Xs = e6.extract_xs(tape) # dictionary (keys are MAT) of dataframes
    matListPert = PertSeriesXs.index.get_level_values("MAT")
    for mat in Xs:
        if mat not in matListPert:
            continue
        mtListPert = PertSeriesXs.loc[mat].index.get_level_values("MT")
        mtListXs = Xs[mat].columns
        for mt in sorted(mtListXs, reverse=True):
            if mt in mtListPert:
                mtPert = mt
            elif mt in range(800,850) and 107 in mtListPert:
                mtPert = 107
            elif mt in range(750,800) and 106 in mtListPert:
                mtPert = 106
            elif mt in range(700,750) and 105 in mtListPert:
                mtPert = 105
            elif mt in range(650,700) and 104 in mtListPert:
                mtPert = 104
            elif mt in range(600,650) and 103 in mtListPert:
                mtPert = 103
            elif mt in range(102,118) and 101 in mtListPert:
                mtPert = 101
            elif mt in (19,20,21,38) and 18 in mtListPert:
                mtPert = 18
            elif mt in (18,101) and 27 in mtListPert:
                mtPert = 27
            elif mt in range(50,92) and 4 in mtListPert:
                mtPert = 4
            else:
                continue
            Xs[mat][mt] = perturb(XsSeries = Xs[mat][mt],
                                  PertSeriesXs = PertSeriesXs.loc[mat,mtPert])
#            Xs[mat] = Xs[mat].add(smp, fill_value=0)
        # Redundant XS
        daughters = [ x for x in range(800,850) if x in Xs[mat].columns]
        if daughters:
            Xs[mat][107] = Xs[mat][daughters].sum(axis=1)
        daughters = [ x for x in range(750,800) if x in Xs[mat].columns]
        if daughters:
            Xs[mat][106] = Xs[mat][daughters].sum(axis=1)
        daughters = [ x for x in range(700,750) if x in Xs[mat].columns]
        if daughters:
            Xs[mat][105] = Xs[mat][daughters].sum(axis=1)
        daughters = [ x for x in range(650,700) if x in Xs[mat].columns]
        if daughters:
            Xs[mat][104] = Xs[mat][daughters].sum(axis=1)
        daughters = [ x for x in range(600,650) if x in Xs[mat].columns]
        if daughters:
            Xs[mat][103] = Xs[mat][daughters].sum(axis=1)
        daughters = [ x for x in range(102,118) if x in Xs[mat].columns]
        if daughters:
            Xs[mat][101] = Xs[mat][daughters].sum(axis=1)
        daughters = [ x for x in (19,20,21,38) if x in Xs[mat].columns]
        if daughters:
            Xs[mat][18] = Xs[mat][daughters].sum(axis=1)
        daughters = [ x for x in (18,101) if x in Xs[mat].columns]
        if daughters:
            Xs[mat][27] = Xs[mat][daughters].sum(axis=1)
        daughters = [ x for x in range(50,92) if x in Xs[mat].columns]
        if daughters:
            Xs[mat][4] = Xs[mat][daughters].sum(axis=1)
        daughters = [ x for x in (4,5,11,*range(16,19),*range(22,38),41,42,44,45) if x in Xs[mat].columns]
        if daughters:
            Xs[mat][3] = Xs[mat][daughters].sum(axis=1)
        daughters = [ x for x in (2,3) if x in Xs[mat].columns]
        if daughters:
            Xs[mat][1] = Xs[mat][daughters].sum(axis=1)
        Xs[mat] = Xs[mat][mtListXs] # keep only mt present in the original file
    tape = e6.update_xs(tape, Xs)
    tape = e6.write_mf3_mt(tape)
    e6.write_tape(tape, output)

def sampling2(tape, PertSeriesChi):
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
    tape = e6.update_chi(tape, Chi)
    tape = e6.write_mf3_mt(tape)
    return tape

if __name__ == '__main__':
    __spec__ = None
    #file = "data_test/H1.txt"
    #file = "96-Cm-242g.jeff33"
    #file = "nubar_endfb80.tape"
    #file = "nubar_jeff33.tape"
    #NSMP = 100
    #file = "JEFF33-rdd_all.asc"
    #file = "26-Fe-56g.jeff33"
    if len(sys.argv) == 1:
        sys.argv.extend(["data_test\92-U-235g.jeff33", "--outdir", os.path.join("..","ttt")])
    settings.init()
    tape = e6.endf2df(settings.args.endf6)
    Chi = e6.extract_chi(tape)
    e6.update_chi(tape, Chi)
    PertChi = sample_chi(tape, settings.args.samples)
    PertXs = sample_xs(tape, settings.args.samples)
    sampling2(tape, PertChi[1])
    #df_nu = extract_nu(tape) # dictionary (keys are MAT) of dataframes

    pool = mp.Pool(processes=settings.args.processes)
    # Problem when running python from windows to linux
    outname = os.path.join(settings.args.outdir, os.path.basename(settings.args.endf6) + '-{}')
#    outname = os.path.join(os.path.basename(settings.args.outdir), settings.args.endf6, '{}')
    [pool.apply( sampling, args=(tape, PertXs[ismp], outname.format(ismp))) for ismp in range(1,settings.args.samples+1) ]