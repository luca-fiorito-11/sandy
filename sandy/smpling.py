# -*- coding: utf-8 -*-
"""
Created on Mon Mar 19 22:51:03 2018

@author: lucaf
"""

import pandas as pd
import sandy.endf6.files as e6
from sandy.cov import Cov
import copy
import numpy as np
import sandy.settings as settings
import sys
import os
import multiprocessing as mp



#if __name__ == "__main__":
#    import argparse
#    parser = argparse.ArgumentParser(description='Run SANDY')
#    parser.add_argument('endf6', help="file in ENDF-6 format")
#    parser.add_argument('N-samples', help="number of samples")
#    parser.add_argument('--pendf', help="file in PENDF format")
#    parser.add_argument('-nw', '--no-write', help="disable writing")
#    parser.add_argument("-v", '--version', action='version',
#                          version='%(prog)s 1.0', help="SANDY's version")
#    args = parser.parse_args()
#    file = args.endf6
#    NSMP = args.N_samples
#else:
#    file = "H1.txt"
#    NSMP = 100




#To produce correlation matrix
#Index = df_cov_xs.index
#df_cov_xs.update(pd.DataFrame(Cov(df_cov_xs.as_matrix()).corr, index=Index, columns=Index))



def sample_xs(tape, NSMP):
    DfCov = e6.merge_covs( e6.extract_cov33(tape) )
    # Must add names to indexes
    DfPert = pd.DataFrame( Cov(DfCov.as_matrix()).sampling(NSMP) + 1 ,
                                 index = DfCov.index,
                                 columns = range(1,NSMP+1))
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


def sampling(tape, PertSeries, output):
    Xs = e6.extract_xs(tape) # dictionary (keys are MAT) of dataframes
    matListPert = PertSeries.index.get_level_values("MAT")
    for mat in Xs:
        if mat not in matListPert:
            continue
        mtListPert = PertSeries.loc[mat].index.get_level_values("MT")
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
                                  PertSeries = PertSeries.loc[mat,mtPert])
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
    print (output)
    e6.write_tape(tape, output)

if __name__ == '__main__':
    #file = "data_test/H1.txt"
    #file = "96-Cm-242g.jeff33"
    #file = "nubar_endfb80.tape"
    #file = "nubar_jeff33.tape"
    #NSMP = 100
    #file = "JEFF33-rdd_all.asc"
    #file = "26-Fe-56g.jeff33"
    if len(sys.argv) == 1:
        sys.argv.extend(["data_test/H1.txt", "--outdir", "../ttt"])
    settings.init()
    tape = pd.DataFrame([[int(x[66:70]), int(x[70:72]), int(x[72:75]), x] for x in e6.split(settings.args.endf6)],
            columns=('MAT', 'MF', 'MT','TEXT'))
    tape = tape.set_index(['MAT','MF','MT']).sort_index() # Multi-indexing
    tape['DATA'] = tape['TEXT'].apply(e6.process_endf_section)

    PertXs = sample_xs(tape, settings.args.samples)
    Xs = e6.extract_xs(tape) # dictionary (keys are MAT) of dataframes
    #df_nu = extract_nu(tape) # dictionary (keys are MAT) of dataframes

    pool = mp.Pool(processes=settings.args.processes)
    # Problem when running python from windows to linux
    outname = settings.args.outdir + '/' + os.path.basename(settings.args.endf6) + '-{}'
#    outname = os.path.join(os.path.basename(settings.args.outdir), settings.args.endf6, '{}')
    [pool.apply( sampling, args=(tape, PertXs[ismp],outname.format(ismp))) for ismp in range(1,settings.args.samples+1) ]