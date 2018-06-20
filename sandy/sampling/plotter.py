# -*- coding: utf-8 -*-
"""
Created on Wed Apr  4 21:43:21 2018

@author: lucaf
"""
import sys
import os

def functi(x):
    from time import sleep
    print("start f(" + x + ")")
    sleep(5)
    print("end   f(" + x + ")")
    return "did " + x

def mycallback(x):
    global blah
    blah = "called back"
    print("My callback " + str(x))



def myerrorcallback(r):
    print("My errorcallback " + str(r))

def process_sp(ismp, inp, text, **kwargs):
    import pandas as pd
    import time
    import sandy.formats.endf6 as e6
    t0 = time.time()
    xs = e6.Endf6.from_text(text).process(keep_mf=[1,3]).get_xs()
    xs.columns = pd.MultiIndex.from_tuples([(mat,mt,ismp) for mat,mt in xs.columns.values], names=["MAT", "MT", "SMP"])
    print("Processed file '{}' in {:.2f} sec".format(inp, time.time()-t0,))
    return xs

def run(iargs=None):
    import pandas as pd
    import sandy.formats.endf6 as e6
    from sandy.formats.errorr import Errorr
    from sandy import settings
    import numpy as np
    import platform
    import multiprocessing as mp
    import time
    from functools import reduce
    from bokeh.layouts import column
    from bokeh.models import CustomJS, ColumnDataSource, Slider, HoverTool, Div
    from bokeh.plotting import figure, output_file, show
    from bokeh.palettes import Spectral4
    from bokeh.io import save

    t0 = time.time()
    settings.init_plotter(iargs)
    # Further setup of settings
    kwargs = vars(settings.args)
    tools = "box_zoom,pan,save,reset"


    # LOAD REFERENCE FILE
    tape = e6.Endf6.from_file(settings.args.file).process(keep_mf=[1,3])
    if tape.empty:
        sys.exit("ERROR: tape is empty")
    xsRef = tape.get_xs()

    # LOAD COVARIANCE FILE
    if settings.args.errorr_cov:
        covtape = Errorr.from_file(settings.args.errorr_cov).process(keep_mf=[31,33])
    elif settings.args.endf6_cov:
        covtape = e6.Endf6.from_file(settings.args.endf6_cov).process(keep_mf=[31,33])
    if covtape.empty:
        sys.exit("ERROR: covtape is empty")
    covRef = covtape.get_cov()

    global DictText
    DictText = {}
    # IMPORT SAMPLES INTO DATAFRAME
    # Keep reading and processing separate to parallelize the processing
    baseoutput = os.path.join(kwargs["outdir"], os.path.basename(kwargs["file"]))
    for ismp in range(1, kwargs["samples"]+1):
        inp = baseoutput + "-" + str(ismp)
        with open(inp) as f:
            DictText.update({ (ismp,inp) :f.read() })

    ListXs = [process_sp(ismp, inp, text, **kwargs) for (ismp,inp),text in DictText.items()]
#    if kwargs["processes"] == 1:
#        ListXs = [process_sp(ismp, inp, **kwargs) for ismp,inp in DictText.keys()]
#    else:
#        if platform.system() == "Windows":
#            def init_pool(the_dicttext):
#                global DictText
#                DictText = the_dicttext
#            pool = mp.Pool(processes=kwargs["processes"],
#                           initializer=init_pool(DictText))
#        else:
#            pool = mp.Pool(processes=kwargs["processes"])
#        queue = Queue()
#        for times in [2,3,0,2]:
#            queue.put(times)
#        def slowstart(q):
#            import os
#            num = q.get()
#            print "slowstart: process id = {0} (sleep({1}))".format(os.getpid(),num)
#            sleep(num)
#        def init_pool(initarg):
#            text
#            global DictText
#            DictText = the_dicttext
#        servers=["s1","s2","s3","s4","s5","s6"]
#        blah = "no callback"
#        with mp.Pool(processes=kwargs["processes"]) as pool:
#            A=[]
#            for server in servers:
#                r = pool.apply_async(functi, (server,),  callback=mycallback, error_callback=myerrorcallback)
#                A.append(r)
#            pool.close()
#            pool.join()
#            print (blah)
#        outs = [pool.apply_async(functi, args=(x)) for x in [1,2,3]]
#                                 args = (ismp, inp),
#                                 ) for ismp,inp in DictText.keys() ]
#        pool.close()
#        ListXs = list(map(lambda x:x.get(), outs))

    xs = reduce(lambda l,r : l.join(r), ListXs).sort_index().interpolate(method='slinear', axis=0).fillna(0)


    width = 1500; height = 300
    for mat in xs.columns.get_level_values('MAT').unique():
        for mt in xs[mat].columns.get_level_values('MT').unique():
            df = xs[mat][mt]
            df.columns = list(map(str, df.columns))
            MEAN = df.mean(axis=1)
            STD = df.std(axis=1).divide(MEAN.values).replace(np.inf, 0).fillna(0) * 100.
            if (STD==0).all(): continue
            df = df.divide(MEAN.values, axis=0).replace(np.inf, 1).fillna(1)
#            mean = pd.Series({j+1:np.mean(df.iloc[:,:3].values) for j in range(1,kwargs["samples"])})
#            std = pd.Series({j+1:np.std(df.iloc[:,:3].values) for j in range(1,kwargs["samples"])})
#            convergence = pd.DataFrame([mean,std], index=["MEAN", "STD"]).T
#            convergence.index.name = "SMP"

            x_axis_type = "log"
            y_axis_type = "linear" if mt in (452,455,456) else "log"
            title = r"Evaluated reaction"
            peval = figure(plot_width=width, plot_height=height, x_axis_type=x_axis_type, y_axis_type=y_axis_type, tools=tools, title=title)
            data = xsRef[mat][mt]
            x = data.index
            y = data.values
            peval.line(x, y, color=Spectral4[1], alpha=.8)
            peval.add_tools(HoverTool(tooltips=[
                    ("E", "@x"),
                    ("XS", "@y")
                    ]))

            x_axis_type = "log"
            y_axis_type = "linear"
            title = r"Standard deviation (%)"
            pstd = figure(x_range=peval.x_range, plot_width=width, plot_height=height, x_axis_type=x_axis_type, y_axis_type=y_axis_type, tools=tools, title=title)
            pstd.circle(x=STD.index, y=STD.values, color=Spectral4[1], alpha=.8, legend=r"samples")
            try:
                data = covRef[mat][mt].loc[mat,mt]
                x = data.index
                y = np.sqrt(np.diag(data.values)) * 100.
                pstd.circle(x, y, color=Spectral4[2], alpha=.8, legend=r"eval")
            except:
                pass
            pstd.legend.location = "top_left"
            pstd.legend.click_policy="mute"
            pstd.add_tools(HoverTool(tooltips=[
                    ("E", "@x"),
                    ("Stdev (%)", "@y")
                    ]))

            x_axis_type = "log"
            y_axis_type = "log"
            title = r"Ratios"
            pratio = figure(x_range=peval.x_range, plot_width=width, plot_height=height, x_axis_type=x_axis_type, y_axis_type=y_axis_type, tools=tools, title=title)
            left = xsRef[mat][mt].rename("DATA").to_frame().reset_index()
            right = MEAN.rename("SMP").to_frame().reset_index()
            data = pd.merge_ordered(left, right, how='right', fill_method='slinear').fillna(0)
            data["RATIO"] = data.SMP.divide(data.DATA.values).replace(np.inf, 0).fillna(0)
            pratio.line(x=data.E, y=data.RATIO, color=Spectral4[3], alpha=.8, legend=r"mean")
            try:
                cov = covRef[mat][mt].loc[mat,mt]
                left = pd.Series(np.sqrt(np.diag(cov.values)) * 100., index=cov.index).rename("DATA").to_frame().reset_index()
                right = STD.rename("SMP").to_frame().reset_index()
                data = pd.merge_ordered(left, right, how='left', fill_method='ffill').fillna(0)
                data["RATIO"] = data.SMP.divide(data.DATA.values).replace(np.inf, 0).fillna(0)
                pratio.circle(x=data.E, y=data.RATIO, color=Spectral4[1], alpha=.8, legend=r"stdev")
            except:
                pass
            pratio.legend.location = "top_left"
            pratio.legend.click_policy="mute"
            pratio.add_tools(HoverTool(tooltips=[
                    ("E", "@x"),
                    ("Ratio", "@y")
                    ]))
            pratio.xaxis.axis_label = 'energy (eV)'

#            source = ColumnDataSource(convergence)
#            pconv = figure(plot_width=1000, plot_height=350, x_axis_type="linear", y_axis_type="linear", tools=tools)
#            pconv.line(x="SMP", y="MEAN", source=source, color=Spectral4[1], alpha=1, legend=r"mean")
#            pconv.line(x="SMP", y="STD", source=source, color=Spectral4[2], alpha=1, legend=r"std")
#            pconv.legend.location = "top_left"

            plots = column(peval, pstd, pratio)
            header = r"""<h1>MAT={} MT={}</h1> #SAMPLES={}

            <hr>

            """.format(mat,mt,kwargs["samples"])
            layout = column(Div(text=header), plots)
            name = "{}-{}-{}.html".format(os.path.basename(baseoutput),mat,mt)
            outfile = os.path.join(kwargs["plotdir"],name)
            output_file(outfile)
            print("save file {}".format(outfile))
            save(layout)
    print("Total running time 'plotter': {:.2f} sec".format(time.time() - t0))



if __name__ == '__main__':
    from sandy import __file__ as sd
    from sandy.data_test import __file__ as td
    sd = os.path.dirname(os.path.realpath(sd))
    td = os.path.dirname(os.path.realpath(td))
    extra_args = [os.path.join(sd, r"sampling\u5-33-tmpdir\pendf"),
                 "9228",
                 "102",
                 "--original", os.path.join(td,r"u235.pendf"),
                 "--cov", os.path.join(td,r"u235.endf"),]
    sys.argv = [sys.argv[0]] + extra_args
    run()