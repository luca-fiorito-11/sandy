# -*- coding: utf-8 -*-
"""
Created on Wed Apr  4 21:43:21 2018

@author: lucaf
"""

#E00visitor
#156916
#
#E10visitor
#536933
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
from os.path import join, dirname, realpath
import matplotlib.pyplot as plt

from bokeh.layouts import column
from bokeh.models import CustomJS, ColumnDataSource, Slider, HoverTool

from bokeh.plotting import figure, output_file, show
from bokeh.palettes import Spectral4

def run():
    tools = "box_zoom,pan,save,reset"
    settings.init_plotter()
    mat = settings.args.mat
    mt = settings.args.mt
    mf,key = (1,"NUBAR") if mt in (452,455,456) else (3,"XS")

    if settings.args.original:
        print ("read best estimates from", settings.args.original)
        Data = e6.extract_xs(e6.endf2df(settings.args.original, keep_mf=[mf], keep_mt=[mt]))[mat,mt]
        Data.name = "Data"

    if settings.args.cov:
        print ("read covariances from", settings.args.cov)
        DfCov = e6.extract_cov33(e6.endf2df(settings.args.cov, keep_mf=[mf+30], keep_mt=[mt]))
        DfCov.index = DfCov.index.droplevel("MAT").droplevel("MT")
        DfCov.columns = DfCov.columns.droplevel("MAT").droplevel("MT")
        Runc = pd.Series(np.sqrt(np.diag(DfCov))*100., index=DfCov.index)
        Runc.name = "Runc"

    ListXs = []
    for inp in os.listdir(settings.args.smpdir):
        print ("read file ", inp)
        D = e6.endf2df(join(settings.args.smpdir,inp), keep_mf=[mf], keep_mt=[mt]).DATA.loc[mat,mf,mt][key]
        D.name = inp
        ListXs.append(D)
    Smp = pd.DataFrame(ListXs).T
    Smp['Mean'] = Smp.mean(axis=1)
    Smp['Std'] = Smp.std(axis=1).fillna(0)
    Smp['Rstd'] = (Smp.Std/Smp.Mean*100).replace([-np.inf,np.inf], np.nan).fillna(0)

    Smp = pd.merge_ordered(Smp.reset_index(), Data.to_frame().reset_index(), on="E", how='left').interpolate(method='slinear', axis=0).fillna(0)
    Smp = pd.merge_ordered(Smp, Runc.to_frame().reset_index(), on="E", how='left').interpolate(method='zero', axis=0).fillna(0)
    Smp["Unc"] = (Smp.Runc*Smp.Data/100.).replace([-np.inf,np.inf], 0)
    Smp["Diff_mean"] = ((Smp.Mean/Smp.Data-1)*100.).replace([-np.inf,np.inf], [-100,100])
    Smp["Diff_rstd"] = ((Smp.Rstd/Smp.Runc-1)*100.).replace([-np.inf,np.inf], [-100,100])

    source = ColumnDataSource(Smp)

    x_axis_type = "log"
    y_axis_type = "linear" if mt in (452,455,456) else "log"
    pstd = figure(plot_width=800, plot_height=250, x_axis_type=x_axis_type, y_axis_type="linear", tools=tools)
    pstd.line(x="E", y="Runc", source=source, color=Spectral4[0], alpha=1, legend=r"data")
    pstd.circle(x="E", y="Rstd", source=source, color=Spectral4[1], alpha=1, legend=r"samples")
    pstd.legend.location = "top_left"
    pstd.legend.click_policy="mute"

    x_axis_type = "log"
    y_axis_type = "linear" if mt in (452,455,456) else "log"
    pmean = figure(x_range=pstd.x_range, plot_width=800, plot_height=250, x_axis_type=x_axis_type, y_axis_type=y_axis_type, tools=tools)
    pmean.line(x="E", y="Data", source=source, color=Spectral4[0], alpha=1, legend=r"data")
    pmean.line(x="E", y="Mean", source=source, color=Spectral4[1], alpha=1, legend=r"samples")
    pmean.legend.location = "top_left"
    pmean.legend.click_policy = "mute"


    Ratio = Smp.copy().drop(["Diff_mean", "Diff_rstd", "Rstd", "Runc"], axis=1).set_index("E")
    Ratio = pd.DataFrame(Ratio.values/Ratio.Data.to_frame().values, index=Ratio.index, columns=Ratio.columns)#.replace([-np.inf,np.inf], [np.NaN,100])
    Ratio["UncTop"] = Ratio.Data + Ratio.Unc
    Ratio["UncBottom"] = Ratio.Data - Ratio.Unc
    Ratio["StdTop"] = Ratio.Data + Ratio.Std
    Ratio["StdBottom"] = Ratio.Data - Ratio.Std
    Ratio.drop(["Std", "Unc"], axis=1, inplace=True)

    source = ColumnDataSource(Ratio.iloc[:,range(10)])

    x_axis_type = "log"
    y_axis_type = "linear"
    pratio = figure(x_range=pstd.x_range, plot_width=1600, plot_height=400, x_axis_type=x_axis_type, y_axis_type=y_axis_type, tools=tools)
    pratio.line(x=Ratio.index.values, y=Ratio.Data.values, line_width=2, color=Spectral4[0], alpha=1, legend=r"data")
    pratio.line(x=Ratio.index.values, y=Ratio.Mean.values, line_width=2, color=Spectral4[1], alpha=1, legend=r"mean")
    pratio.line(x=Ratio.index.values, y=Ratio.UncTop.values, line_width=2, line_dash="dotted", color=Spectral4[2], alpha=1, legend=r"stdev")
    pratio.line(x=Ratio.index.values, y=Ratio.UncBottom.values, line_width=2, line_dash="dotted", color=Spectral4[2], alpha=1)
    pratio.line(x=Ratio.index.values, y=Ratio.StdTop.values, line_width=2, line_dash="dotted", color=Spectral4[3], alpha=1, legend=r"unc")
    pratio.line(x=Ratio.index.values, y=Ratio.StdBottom.values, line_width=2, line_dash="dotted", color=Spectral4[3], alpha=1)
    Ratio.drop(["Data", "Mean", "UncTop", "UncBottom", "StdTop", "StdBottom"], axis=1, inplace=True)
    numcols = len(Ratio.columns)
#    mypalette=Spectral4[4:4+numcols]
    names = []
    for name in  Ratio:
        names.append(name)
        l = pratio.line(x=Ratio.index.values,
                        y=Ratio[name].values,
                        line_color='gray',
                        name=name,
                        line_alpha=0.5,
                        hover_line_alpha=1.0,)
        pratio.add_tools(HoverTool(renderers=[l], tooltips=[
                ('Name', name),
                ("E", "@x"),
                ("y", "@y")
            ]))

#    pratio.add_tools(HoverTool(tooltips=[
#        ('Name', '$name'),
#        ("E", "@x"),
#        ("y", "@y")
#        ]))

    pratio.legend.location = "top_left"
    pratio.legend.click_policy = "mute"



    layout = column(pmean, pstd, pratio)
    show(layout)

#    for i in range(1,6):
#        A[i] /= A[0]
#    ListPerts = []
#    rowold = pd.Series()
#    for e,row in A.iterrows():
#        row.drop(0, inplace=True)
#        if rowold.empty:
#            e0 = deepcopy(e)
#        elif not np.allclose(row, rowold):
#            ListPerts.append(( e0, e, *rowold.tolist() ))
#            e0 = deepcopy(e)
#        rowold = deepcopy(row)
#    B = pd.DataFrame.from_records(ListPerts)
#    B.set_index([0, 1], inplace=True)
#    C = B.as_matrix()
#    C = np.insert(C, C.shape[0], C[-1]*C.shape[1], axis=0)
#    C = np.insert(C, C.shape[1], C[:,-1]*C.shape[0], axis=1)
#    C = np.corrcoef(C)
#    x = B.index.get_level_values(0).tolist() + [B.index.get_level_values(1).tolist()[-1]]
#    e6.plot_heatmap(x, x, C,
#                    xscale="log", yscale="log",
#                    vmin=-1, vmax=1,
#                    cmap="bwr",
#                    xlabel=None, ylabel=None, title=None)

if __name__ == '__main__':
    from sandy import __file__ as sd
    from sandy.data_test import __file__ as td
    sd = dirname(realpath(sd))
    td = dirname(realpath(td))
    extra_args = [join(sd, r"sampling\cm242-tmpdir"),
                 "9631",
                 "102",
                 "--original", join(td,r"u235.pendf"),
                 "--cov", join(td,r"cm242.endf"),]
    sys.argv = [sys.argv[0]] + extra_args
    run()