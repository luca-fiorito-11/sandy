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
file = "data_test/H1.txt"
#file = "96-Cm-242g.jeff33"
#file = "nubar_endfb80.tape"
#file = "nubar_jeff33.tape"
NSMP = 100
#file = "JEFF33-rdd_all.asc"
#file = "26-Fe-56g.jeff33"
tape = pd.DataFrame([[int(x[66:70]), int(x[70:72]), int(x[72:75]), x] for x in e6.split(file)],
        columns=('MAT', 'MF', 'MT','TEXT'))
tape = tape.set_index(['MAT','MF','MT']).sort_index() # Multi-indexing
tape['DATA'] = tape['TEXT'].apply(e6.process_endf_section)



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

PertXs = sample_xs(tape, NSMP)
#unique_ids = list(set(zip(df_cov_xs.index.get_level_values(0), df_cov_xs.index.get_level_values(1))))
#idx = pd.IndexSlice
#df_samples.loc[idx[:, :, ['C1', 'C3']], idx[:, 'foo']]
XsOrig = e6.extract_xs(tape) # dictionary (keys are MAT) of dataframes
#df_nu = extract_nu(tape) # dictionary (keys are MAT) of dataframes
orig = copy.deepcopy(Xs)
# Add index level with sample number, 0 is the reference
for mat,df in XsOrig.items():
    df.columns = pd.MultiIndex.from_product([[0], df.columns],
                                            names=("ISMP","MT"))

def perturb(XsSeries, PertDf):
    mt = XsSeries.name
    P = e6.pandas_interpolate(PertDf,
                           np.unique(XsSeries.index).tolist(),
                           axis="rows")
    XS = P.multiply(XsSeries, axis="index")
    # Negative values are set to zero
    XS[XS <= 0] = 0
    XS.columns = pd.MultiIndex.from_product([XS.columns, [mt]],
                                            names=("ISMP","MT"))
    return XS

matListPert = PertXs.index.get_level_values("MAT")
for mat in orig:
    if mat not in matListPert:
        continue
    mtListPert = PertXs.loc[mat].index.get_level_values("MT")
    mtListXs = orig[mat].columns
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
        XsSeries = orig[mat][mt]
        PertDf = PertXs.loc[mat,mtPert]
        smp = perturb(XsSeries, PertDf)
        Xs[mat] = Xs[mat].add(smp, fill_value=0)
    # Redundant XS
    for ismp in range(NSMP+1):
        daughters = [ x for x in range(800,850) if x in Xs[mat][ismp].columns]
        if daughters:
            Xs[mat][(ismp,107)] = Xs[mat][ismp][daughters].sum(axis=1)
        daughters = [ x for x in range(750,800) if x in Xs[mat][ismp].columns]
        if daughters:
            Xs[mat][(ismp,106)] = Xs[mat][ismp][daughters].sum(axis=1)
        daughters = [ x for x in range(700,750) if x in Xs[mat][ismp].columns]
        if daughters:
            Xs[mat][(ismp,105)] = Xs[mat][ismp][daughters].sum(axis=1)
        daughters = [ x for x in range(650,700) if x in Xs[mat][ismp].columns]
        if daughters:
            Xs[mat][(ismp,104)] = Xs[mat][ismp][daughters].sum(axis=1)
        daughters = [ x for x in range(600,650) if x in Xs[mat][ismp].columns]
        if daughters:
            Xs[mat][(ismp,103)] = Xs[mat][ismp][daughters].sum(axis=1)
        daughters = [ x for x in range(102,118) if x in Xs[mat][ismp].columns]
        if daughters:
            Xs[mat][(ismp,101)] = Xs[mat][ismp][daughters].sum(axis=1)
        daughters = [ x for x in (19,20,21,38) if x in Xs[mat][ismp].columns]
        if daughters:
            Xs[mat][(ismp,18)] = Xs[mat][ismp][daughters].sum(axis=1)
        daughters = [ x for x in (18,101) if x in Xs[mat][ismp].columns]
        if daughters:
            Xs[mat][(ismp,27)] = Xs[mat][ismp][daughters].sum(axis=1)
        daughters = [ x for x in range(50,92) if x in Xs[mat][ismp].columns]
        if daughters:
            Xs[mat][(ismp,4)] = Xs[mat][ismp][daughters].sum(axis=1)
        daughters = [ x for x in (4,5,11,*range(16,19),*range(22,38),41,42,44,45) if x in Xs[mat][ismp].columns]
        if daughters:
            Xs[mat][(ismp,3)] = Xs[mat][ismp][daughters].sum(axis=1)
        daughters = [ x for x in (2,3) if x in Xs[mat][ismp].columns]
        if daughters:
            Xs[mat][(ismp,1)] = Xs[mat][ismp][daughters].sum(axis=1)
        Xs[mat][ismp] = Xs[mat][ismp][Xs[mat][ismp].columns]

ixs = copy.deepcopy(orig)
for mat, df_samples_bymat in df_samples_xs.groupby(level=0, sort=True):
    for mt, P in df_samples_bymat.groupby(level=1, sort=True):
        if mat not in ixs:
            continue
        if mt not in ixs[mat]:
            continue
        XS = ixs[mat][mt]
        P = pandas_interpolate(P.loc[mat,mt][ismp],
                               list(XS.index.get_level_values("E")),
                               axis="rows")
        XS = XS.multiply(P)
        # Negative values are set to mean
        XS[XS <= 0] = ixs[mat][mt][XS <= 0]
        ixs[mat][mt][XS > 0] = XS[XS > 0]
    # Redundant XS
    columns = ixs[mat].columns
    daughters = [ x for x in range(800,850) if x in ixs[mat]]
    if daughters:
        ixs[mat][107] = ixs[mat][daughters].sum(axis=1)
    daughters = [ x for x in range(750,800) if x in ixs[mat]]
    if daughters:
        ixs[mat][106] = ixs[mat][daughters].sum(axis=1)
    daughters = [ x for x in range(700,750) if x in ixs[mat]]
    if daughters:
        ixs[mat][105] = ixs[mat][daughters].sum(axis=1)
    daughters = [ x for x in range(650,700) if x in ixs[mat]]
    if daughters:
        ixs[mat][104] = ixs[mat][daughters].sum(axis=1)
    daughters = [ x for x in range(600,650) if x in ixs[mat]]
    if daughters:
        ixs[mat][103] = ixs[mat][daughters].sum(axis=1)
    daughters = [ x for x in range(102,118) if x in ixs[mat]]
    if daughters:
        ixs[mat][101] = ixs[mat][daughters].sum(axis=1)
    daughters = [ x for x in (19,20,21,38) if x in ixs[mat]]
    if daughters:
        ixs[mat][18] = ixs[mat][daughters].sum(axis=1)
    daughters = [ x for x in (18,101) if x in ixs[mat]]
    if daughters:
        ixs[mat][27] = ixs[mat][daughters].sum(axis=1)
    daughters = [ x for x in range(50,92) if x in ixs[mat]]
    if daughters:
        ixs[mat][4] = ixs[mat][daughters].sum(axis=1)
    daughters = [ x for x in (4,5,11,*range(16,19),*range(22,38),41,42,44,45) if x in ixs[mat]]
    if daughters:
        ixs[mat][3] = ixs[mat][daughters].sum(axis=1)
    daughters = [ x for x in (2,3) if x in ixs[mat]]
    if daughters:
        ixs[mat][1] = ixs[mat][daughters].sum(axis=1)
    ixs[mat] = ixs[mat][columns]
tape = update_xs(tape, ixs)
tape = write_mf3_mt(tape)
write_tape(tape, "AAA-{}".format(ismp))

for ismp in range(NSMP):
    ixs = copy.deepcopy(df_xs)
    for mat, df_samples_bymat in df_samples_xs.groupby(level=0, sort=True):
        for mt, P in df_samples_bymat.groupby(level=1, sort=True):
            if mat not in ixs:
                continue
            if mt not in ixs[mat]:
                continue
            XS = ixs[mat][mt]
            P = pandas_interpolate(P.loc[mat,mt][ismp],
                                   list(XS.index.get_level_values("E")),
                                   axis="rows")
            XS = XS.multiply(P)
            # Negative values are set to mean
            ixs[mat][mt][XS > 0] = XS[XS > 0]
        # Redundant XS
        columns = ixs[mat].columns
        daughters = [ x for x in range(800,850) if x in ixs[mat]]
        if daughters:
            ixs[mat][107] = ixs[mat][daughters].sum(axis=1)
        daughters = [ x for x in range(750,800) if x in ixs[mat]]
        if daughters:
            ixs[mat][106] = ixs[mat][daughters].sum(axis=1)
        daughters = [ x for x in range(700,750) if x in ixs[mat]]
        if daughters:
            ixs[mat][105] = ixs[mat][daughters].sum(axis=1)
        daughters = [ x for x in range(650,700) if x in ixs[mat]]
        if daughters:
            ixs[mat][104] = ixs[mat][daughters].sum(axis=1)
        daughters = [ x for x in range(600,650) if x in ixs[mat]]
        if daughters:
            ixs[mat][103] = ixs[mat][daughters].sum(axis=1)
        daughters = [ x for x in range(102,118) if x in ixs[mat]]
        if daughters:
            ixs[mat][101] = ixs[mat][daughters].sum(axis=1)
        daughters = [ x for x in (19,20,21,38) if x in ixs[mat]]
        if daughters:
            ixs[mat][18] = ixs[mat][daughters].sum(axis=1)
        daughters = [ x for x in (18,101) if x in ixs[mat]]
        if daughters:
            ixs[mat][27] = ixs[mat][daughters].sum(axis=1)
        daughters = [ x for x in range(50,92) if x in ixs[mat]]
        if daughters:
            ixs[mat][4] = ixs[mat][daughters].sum(axis=1)
        daughters = [ x for x in (4,5,11,*range(16,19),*range(22,38),41,42,44,45) if x in ixs[mat]]
        if daughters:
            ixs[mat][3] = ixs[mat][daughters].sum(axis=1)
        daughters = [ x for x in (2,3) if x in ixs[mat]]
        if daughters:
            ixs[mat][1] = ixs[mat][daughters].sum(axis=1)
        ixs[mat] = ixs[mat][columns]
    tape = update_xs(tape, ixs)
    tape = write_mf3_mt(tape)
    write_tape(tape, "AAA-{}".format(ismp))