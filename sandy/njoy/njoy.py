# -*- coding: utf-8 -*-
"""
Created on Fri May 25 16:58:11 2018

@author: fiorito_l
"""
import pandas as pd
import numpy as np
import os, sys, time, pdb, shutil, re, pytest
from ..functions import run_process

sab = pd.DataFrame.from_records([[48,9237,1,1,241,'uuo2'],
                                  [42,125,0,8,221,'tol'],
                                  [59,1425,0,1,221,'si'],
                                  [37,125,11,2,223,'pol'],
                                  [2,125,0,2,221,'ph'],
                                  [12,128,0,2,221,'pd'],
                                  [75,825,1,1,239,'ouo2'],
                                  [48,825,0,3,221,'osap'],
                                  [51,825,0,1,222,'ohw'],
                                  [46,825,3,1,237,'obeo'],
                                  [3,125,0,2,221,'oh2'],
                                  [13,128,0,2,221,'od2'],
                                  [52,1225,0,1,249,'mg'],
                                  [38,125,0,12,221,'mesi'],
                                  [10,125,1,2,221,'ice'],
                                  [7,125,12,1,225,'hzr'],
                                  [1,125,0,2,222,'lw'],
                                  [8,125,0,2,249,'hca'],
                                  [31,600,1,1,229,'gr'],
                                  [11,128,0,2,221,'dhw'],
                                  [59,2025,0,1,249,'cah'],
                                  [27,425,3,1,233,'bbeo'],
                                  [26,425,2,1,231,'be'],
                                  [60,1325,0,2,221,'asap']],
columns = ['matde','matdp','icoh','natom','mtref','ext'])


class OutputFiles(pd.DataFrame):

    def __init__(self, *args, **kwargs):
        kwargs.update({"columns" : ["id", "format", "file"]})
        super().__init__(*args, **kwargs)

class NjoyMessages(pd.DataFrame):

    def __init__(self, *args, **kwargs):
        kwargs.update({"columns" : ["module", "type", "routine", "message"]})
        super().__init__(*args, **kwargs)

class NjoyOutput(pd.DataFrame):

    twoline_msg = "-{3}(?P<type>message) from (?P<mod>\w*)-{3}(?P<msg>.*?\n\s{26}.*?)\n"
    oneline_msg = "-{3}(?P<type>message) from (?P<mod>\w*)-{3}(?P<msg>.*?)\n(?!\s{26}.*?\n)"
    twoline_err = "\*{3}(?P<type>error) in (?P<mod>\w*)\*{3}(?P<msg>.*?\n\s{22}.*?)\n"
    oneline_err = "\*{3}(?P<type>error) in (?P<mod>\w*)\*{3}(?P<msg>.*?)\n(?!\s{22}.*?\n)"

    @classmethod
    def from_file(cls, file):
        """
        Read NJOY output file and call from_text method.
        """
        with open(file) as f: text = f.read()
        return cls.from_text(text)

    @classmethod
    def from_text(cls, text):
        """
        Read NJOY output and split it into column based on module.
        """
        splitted = re.split("(moder|acer|groupr|heatr|reconr|broadr|thermr|purr|unresr|errorr)\.\.\.", text)[1:]
        modules = splitted[::2]
        contents = splitted[1::2]
        return cls([*zip(modules,contents)])

    def __init__(self, *args, **kwargs):
        kwargs.update({"columns" : ["module", "content"]})
        super().__init__(*args, **kwargs)

    def get_messages(self):
        found = []
        for i,item in self.iterrows():
            found += [ (item.module, *x) for x in re.findall(self.oneline_msg, item.content)]
            found += [ (item.module, *x) for x in re.findall(self.twoline_msg, item.content)]
            found += [ (item.module, *x) for x in re.findall(self.oneline_err, item.content)]
            found += [ (item.module, *x) for x in re.findall(self.twoline_err, item.content)]
        return NjoyMessages(found)



class Njoy:

    def __init__(self, **kwargs):
#        self.sab = os.path.join(os.path.dirname(os.path.abspath(__file__)), "sab.csv")
#        self.wdir = os.path.normpath(wdir)
#        self.evaluationName = os.path.basename(self.wdir)
        self.processes = None
        self.capsys = False
        self.iwt = 4
        self.legendre = 1
        self.legendregg = 6
        self.scatteringLaw = None
        self.eFiss = None
        self.branchingNG = None
        self.branchingN2N = None
        self.gstr = 0
        self.oldlib = None
        self.sgref = 1.0E10
        self.yields = None
        self.iburn = -1
        self.ip1opt = 0
        self.ifprod = 0
        self.jp1 = 0
        self.iverw = 4
        # General options
        self.EXEC = "njoy2016"
        self.WDIR = os.getcwd()
        self.FOLDER = os.path.abspath("lib")
        self.IPRINT = 1  # by default set verbosity to max
        self.PLOT = False
        self.TEMPS = [293.6]
        # Run options
        self.PENDF = False
        self.GROUPR = False
        self.ACE = False
        self.ERRORR = False
        # modules
        self.RECONR = True
        self.BROADR = True
        self.UNRESR = True
        self.THERMR = False
        self.HEATR = False
        self.GASPR = False
        self.PURR = False
        self.ACER = False
        # RECONR/BROADR default options
        self.ERR = 0.001
        self.ERRMAX = None
        self.THNMAX = 1E6
        # PURR/UNRESR default options
        self.PTABLE = False
        self.SIG0 = [1e10] # default is infinite dilution only
        # THERMR default options
        self.FREE_GAS = False
        self.IIN = 1
        self.ICOH = 1
        self.IFORM = 0
        self.NATOM = 1
        self.MTREF = 221
        # ACER default options
        self.SUFFIXES = [".00"]
        self.IOPT = 1
        self.ITYPE = 1
        self.NEWFOR = 1
        self.IOPP = 1
        self.HK = ''
        # GROUPR default options
        self.IGN = 2 # neutron group structure for GROUPR
        self.IGG = 0
        self.IWT = 6
        self.LORD = 1
        self.NSTR = None
        self.GSTR = None
        # ERRORR default options
        self.IGNE = 2 # neutron group structure for ERRORR
        self.IWTE = 6
        self.IRELCO = 1
        self.IRESPR = 1
        # HEATR default options
        self.KERMA = []
        self.QA = []
        self.LOCAL = 0
        self.__dict__.update({k.upper() : v for k,v in kwargs.items()})
        # Run checks
        if len(self.TEMPS)> 10: raise NotImplementedError("cannot have more than 10 temperatures")

#    def runNjoy(self, inputfile, cwd, timeout=.1):
#        """
#        Run ``NJOY``.
#        Pass input file to process's stdin (it must be encoded first).
#        ..Important::
#            In ``Python 3`` you need to convert string to bytes with a
#            ``encode()`` function
#        """
#        import subprocess as sp
##        pdb.set_trace()
##        if self.CAPSYS:
##            stderr = stdout = sp.PIPE
##        else:
##            stderr = stdout = None
#        process = sp.Popen(self.EXEC,
#                           shell=False,
#                           cwd=cwd,
#                           stdin=open(inputfile),
#                           stdout=sp.PIPE,
#                           stderr=sp.PIPE,)
##        inp = open(inputfile).read().encode()
#        try:
##            stdoutdata, stderrdata = process.communicate(inp, timeout=timeout)
#            stdoutdata, stderrdata = process.communicate(timeout=timeout)
#        except sp.TimeoutExpired as exc:
#            process.kill()
#            stdoutdata, stderrdata = process.communicate()
#            pdb.set_trace()
#            raise NotImplementedError("'{}' took longer than {} seconds to complete and it was killed".format(self.EXEC, timeout))
#        if not self.CAPSYS:
#            print(stdoutdata.decode('utf-8').rstrip())
#        pdb.set_trace()
#        if process.returncode not in [0, 24]:
#            raise NotImplementedError("NJOY exit status {}, cannot run njoy executable".format(process.returncode))

#    def runBash(self, inputfile, cwd):
#        from subprocess import Popen, PIPE
#        if self.CAPSYS:
#            stderr = stdout = PIPE
#        else:
#            stderr = stdout = None
#        command = "bash {}".format(inputfile)
#        process = Popen(command,
#                        shell=True,
#                        cwd=cwd,
#                        stdin=PIPE,
#                        stdout=stdout,
#                        stderr=stderr)
#        stdoutdata, stderrdata = process.communicate()
#        if process.returncode not in [0, 24]:
#            raise NotImplementedError("exit status {} when running \"{}\"".format(process.returncode, command))

    @staticmethod
    def moder_input(tapein, tapeout, **kwargs):
        text = """moder
{} {}
""".format(tapein, tapeout)
        return text

    @staticmethod
    def reconr_input(endfin, pendfout, **kwargs):
        kwargs.update({"ENDFIN" : endfin, "PENDFOUT" : pendfout})
        if not kwargs["ERRMAX"]: kwargs["ERRMAX"] = kwargs["ERR"]*10
        text = """reconr
%(ENDFIN)d %(PENDFOUT)d
'pendf tape from '/
%(MAT)d 1 0/
%(ERR)E  0.  %(ERRMAX)E/
'%(TAG)s with NJOY on %(TIME)s' /
0/
""" % kwargs
        return text

    @staticmethod
    def broadr_input(endfin, pendfin, pendfout, **kwargs):
        kwargs.update({"ENDFIN" : endfin, "PENDFIN" : pendfin, "PENDFOUT" : pendfout})
        kwargs["NTEMPS"] = len(kwargs["TEMPS"])
        kwargs["TEXTTEMPS"] = " ".join(["%E"%tmp for tmp in kwargs["TEMPS"]])
        text = """broadr
%(ENDFIN)d %(PENDFIN)d %(PENDFOUT)d
%(MAT)d %(NTEMPS)d 0 0 0./
%(ERR)E / %(THNMAX)E %(ERR)E/
%(TEXTTEMPS)s/
0/
""" % kwargs
        return text

    @staticmethod
    def thermr_input(endfin, pendfin, pendfout, **kwargs):
        kwargs.update({"ENDFIN" : endfin, "PENDFIN" : pendfin, "PENDFOUT" : pendfout})
        kwargs["NTEMPS"] = len(kwargs["TEMPS"])
        kwargs["TEXTTEMPS"] = " ".join(["%E"%tmp for tmp in kwargs["TEMPS"]])
        text = """thermr
%(ENDFIN)d %(PENDFIN)d %(PENDFOUT)d
0 %(MAT)d 20 %(NTEMPS)d %(IIN)d %(ICOH)d %(IFORM)d %(NATOM)d 221 %(IPRINT)d
%(TEXTTEMPS)s/
0.001 4.0
""" % kwargs
        return text

    @staticmethod
    def purr_input(endfin, pendfin, pendfout, **kwargs):
        kwargs.update({"ENDFIN" : endfin, "PENDFIN" : pendfin, "PENDFOUT" : pendfout})
        kwargs["NSIG0"] = len(kwargs["SIG0"])
        kwargs["TEXTSIG0"] = " ".join(["%E"%dil for dil in kwargs["SIG0"]])
        kwargs["NTEMPS"] = len(kwargs["TEMPS"])
        kwargs["TEXTTEMPS"] = " ".join(["%E"%tmp for tmp in kwargs["TEMPS"]])
        text = """purr
%(ENDFIN)d %(PENDFIN)d %(PENDFOUT)d
%(MAT)d %(NTEMPS)d %(NSIG0)d 20 32/
%(TEXTTEMPS)s/
%(TEXTSIG0)s/
0/
""" % kwargs
        return text

    @staticmethod
    def gaspr_input(endfin, pendfin, pendfout, **kwargs):
        kwargs.update({"ENDFIN" : endfin, "PENDFIN" : pendfin, "PENDFOUT" : pendfout})
        text = """gaspr
%(ENDFIN)d %(PENDFIN)d %(PENDFOUT)d /
""" % kwargs
        return text

    @staticmethod
    def unresr_input(endfin, pendfin, pendfout, **kwargs):
        kwargs.update({"ENDFIN" : endfin, "PENDFIN" : pendfin, "PENDFOUT" : pendfout})
        kwargs["NSIG0"] = len(kwargs["SIG0"])
        kwargs["TEXTSIG0"] = " ".join(["%E"%dil for dil in kwargs["SIG0"]])
        kwargs["NTEMPS"] = len(kwargs["TEMPS"])
        kwargs["TEXTTEMPS"] = " ".join(["%E"%tmp for tmp in kwargs["TEMPS"]])
        text = """unresr
%(ENDFIN)d %(PENDFIN)d %(PENDFOUT)d
%(MAT)d %(NTEMPS)d %(NSIG0)d 1/
%(TEXTTEMPS)s/
%(TEXTSIG0)s/
0/
""" % kwargs
        return text

    @staticmethod
    def acer_input(nendf, npendf, nace, ndir, **kwargs):
        kwargs.update({"NENDF" : nendf, "NPENDF" : npendf, "NACE" : nace, "NDIR" : ndir})
        text = """acer
%(NENDF)d %(NPENDF)d 0 %(NACE)d %(NDIR)d /
%(IOPT)d %(IPRINT)d %(ITYPE)d %(SUFF)s 0/
'%(HK)s'/
%(MAT)d  %(TMP)E /
%(NEWFOR)d %(IOPP)d/
0.001/
""" % kwargs
#        if kwargs["VERBOSE"] >= 2:
#            text += """acer / Check ACE files
#0 %(NACE) 0 0 0 /
#7 %(IPRINT)d 1 -1/
#/
#"""
        return text

    @staticmethod
    def heatr_input(endfin, pendfin, pendfout, **kwargs):
        text = ""
        npks = len(kwargs["KERMA"])
        size = 7
        sign = np.sign(pendfout)
        kwargs.update({"ENDFIN" : endfin, "PENDFIN" : pendfin})
        kwargs.update({"NQA" : len(kwargs["QA"])//2})
        kwargs["TEXTMTA"] = " ".join(map(lambda x: str(int(x)), kwargs["QA"][::2]))
        kwargs["TEXTQA"] = " ".join(map(str, kwargs["QA"][1::2]))
        out = 79
        for i in range(0, npks, size):
            kwargs["PENDFOUT"] = pendfout if i+size >= npks else sign*(out + 1)
            kwargs["NPK"] = len(kwargs["KERMA"][i:i+size])
            kwargs["TEXTPKS"] = " ".join(map(str,kwargs["KERMA"][i:i+size]))
            text += """heatr
%(ENDFIN)d %(PENDFIN)d %(PENDFOUT)d 0 /
%(MAT)d %(NPK)d %(NQA)d 0 %(LOCAL)d %(IPRINT)d /
%(TEXTPKS)s /
""" % kwargs
            if kwargs["NQA"] > 0:
                text += """%(TEXTMTA)s /
%(TEXTQA)s /
""" % kwargs
            kwargs.update({"PENDFIN" : kwargs["PENDFOUT"]})
        return text

    @staticmethod
    def groupr_input(nendf, npendf, ngin, ngout, **kwargs):
        from ..formats.records import write_tab1
        kwargs.update({"NENDF" : nendf, "NPENDF" : npendf, "NGIN" : ngin, "NGOUT" : ngout})
        kwargs["NSIG0"] = len(kwargs["SIG0"])
        kwargs["TEXTSIG0"] = " ".join(["%E"%dil for dil in kwargs["SIG0"]])
        kwargs["NTEMPS"] = len(kwargs["TEMPS"])
        kwargs["TEXTTEMPS"] = " ".join(["%E"%tmp for tmp in kwargs["TEMPS"]])
        if not isinstance(kwargs["IGN"], int):
            if not os.path.isfile(kwargs["IGN"]): raise NotImplementedError("file {} not found".format(kwargs["IGN"]))
            grid = np.genfromtxt(kwargs["IGN"])
            ngroups = len(grid) - 1
            kwargs["IGNSTR"] = "{} /\n".format(ngroups) + "\n".join(map(str,grid)) + " /\n"
            kwargs["IGN"] = 1
        if not isinstance(kwargs["IWT"], int):
            if kwargs["IWT"].lower() == "jaea_fns_175":
                from ..spectra import jaea_fns_175 as spectrum
            else:
                if not os.path.isfile(kwargs["IWT"]): raise NotImplementedError("file {} not found".format(kwargs["IWT"]))
                spectrum = pd.read_csv(kwargs["IWT"], header=None, names=["E", "F"])
            tab = spectrum.reset_index().values.flatten()
            x = tab[::2]; y = tab[1::2]
            kwargs["IWTSTR"] = "\n".join(write_tab1(0, 0, 0, 0, [len(x)], [2], x, y)) + "/\n"
            kwargs["IWT"] = 1
        text = """groupr
%(NENDF)d %(NPENDF)d %(NGIN)d %(NGOUT)d /
%(MAT)d %(IGN)d %(IGG)d %(IWT)d %(LORD)d %(NTEMPS)d %(NSIG0)d %(IPRINT)d /
'' /
%(TEXTTEMPS)s /
%(TEXTSIG0)s /
""" % kwargs
        if "IGNSTR" in kwargs:
            text += "%(IGNSTR)s" % kwargs
        if "IWTSTR" in kwargs:
            text += "%(IWTSTR)s" % kwargs # card13
        for tmp in kwargs["TEMPS"]:
            text += """3/
6/
0/
"""
        text += "0/\n"
        return text

    @staticmethod
    def errorr_input(nendf, npendf, ngout, nout, **kwargs):
        from ..formats.records import write_tab1
        kwargs.update({"NENDF" : nendf, "NPENDF" : npendf, "NGOUT" : ngout, "NOUT" : nout})
        if not isinstance(kwargs["IGNE"], int):
            if not os.path.isfile(kwargs["IGNE"]): raise NotImplementedError("file {} not found".format(kwargs["IGNE"]))
            grid = np.genfromtxt(kwargs["IGNE"])
            ngroups = len(grid) - 1
            kwargs["IGNESTR"] = "{} /\n".format(ngroups) + "\n".join(map(str,grid)) + " /\n"
            kwargs["IGNE"] = 1
        if not isinstance(kwargs["IWTE"], int):
            if kwargs["IWTE"].lower() == "jaea_fns_175":
                from ..csv.spectra import jaea_fns_175 as spectrum
            else:
                if not os.path.isfile(kwargs["IWTE"]): raise NotImplementedError("file {} not found".format(kwargs["IWTE"]))
                spectrum = pd.read_csv(kwargs["IWTE"], header=None, names=["E", "F"])
            tab = spectrum.reset_index().values.flatten()
            x = tab[::2]; y = tab[1::2]
            kwargs["IWTESTR"] = "\n".join(write_tab1(0, 0, 0, 0, [len(x)], [2], x, y)) + "/\n"
            kwargs["IWTE"] = 1
        text = """errorr
%(NENDF)d %(NPENDF)d %(NGOUT)d %(NOUT)d /
%(MAT)d %(IGNE)d %(IWTE)d %(IPRINT)d %(IRELCO)d /
%(IPRINT)d %(TMP)E /
0 33 %(IRESPR)d %(LORD)d /
""" % kwargs
        if "IGNESTR" in kwargs:
            text += "%(IGNESTR)s" % kwargs # card12 + card12b
        if "IWTESTR" in kwargs:
            text += "%(IWTESTR)s" % kwargs # card13
        return text



    def get_pendf(self, **fileOptions):
        from ..functions import force_symlink
        kwargs = dict(self.__dict__, **fileOptions)
        if not os.path.isfile(kwargs["TAPE"]): raise NotImplementedError("evaluation file {} not found".format(kwargs["TAPE"]))
        mydir = os.path.join(kwargs["FOLDER"], kwargs["FILENAME"])
        os.makedirs(mydir, exist_ok=True)
        kwargs["TIME"] = time.ctime(time.time())

        force_symlink(kwargs["TAPE"], os.path.join(mydir, "tape20"))
        text = self.moder_input(20, -21, **kwargs)
        if "PENDFIN" in kwargs:
            if not os.path.isfile(kwargs["PENDFIN"]): raise NotImplementedError("file {} not found".format(kwargs["PENDFIN"]))
            force_symlink(kwargs["PENDFIN"], os.path.join(mydir, "tape30"))
            text += self.moder_input(30, -22, **kwargs)
        else:
            text += self.reconr_input(-21, -22, **kwargs)
        e = -21; p = -22
        if kwargs["BROADR"]:
            text += self.broadr_input(e, p, -23, **kwargs)
            p = -23
        if kwargs["FREE_GAS"]: kwargs["THERMR"] = True
        if kwargs["THERMR"]:
            text += self.thermr_input(0, p, -24, **kwargs)
            p = -24
        if kwargs["UNRESR"]:
            text += self.unresr_input(e, p, -25, **kwargs)
            p = -25
        if kwargs["PTABLE"] > 0: kwargs["PURR"] = True
        if kwargs["PURR"]:
            text += self.purr_input(e, p, -26, **kwargs)
            p = -26
        if len(kwargs["KERMA"]) > 0: kwargs["HEATR"] = True
        if kwargs["HEATR"]:
            text += self.heatr_input(e, p, -27, **kwargs)
            p = -27
        if kwargs["GASPR"]:
            text += self.gaspr_input(e, p, -28, **kwargs)
            p = -28
        text += self.moder_input(p, 29, **kwargs)
        text += "stop"

        DfOutputs = OutputFiles()
        inputfile = os.path.join(mydir, "input_pendf.{}".format(kwargs["TAG"]))
        DfOutputs = DfOutputs.append({"id" : "input_pendf", "format" : "TEXT", "file" : inputfile}, ignore_index=True)
        with open(inputfile,'w') as f: f.write(text)
        print(" --- run pendf for {} ---".format(kwargs["TAG"]))
        returncode, stdout, stderr = run_process("{} < {}".format(kwargs["EXEC"], inputfile), cwd=mydir, verbose=not kwargs["CAPSYS"])

        oldfile = os.path.join(mydir, "tape29")
        newfile = os.path.join(mydir, "{}.pendf".format(kwargs["TAG"]))
        if os.path.isfile(oldfile):
            if os.path.getsize(oldfile) > 0:
                shutil.move(oldfile, newfile)
                DfOutputs = DfOutputs.append({"id" : "pendf", "format" : "ENDF", "file" : newfile}, ignore_index=True)

        oldfile = os.path.join(mydir, "output")
        newfile = os.path.join(mydir, "output_pendf.{}".format(kwargs["TAG"]))
        shutil.move(oldfile, newfile)
        DfOutputs = DfOutputs.append({"id" : "output_pendf", "format" : "TEXT", "file" : newfile}, ignore_index=True)
        DfMessages = NjoyOutput.from_file(newfile).get_messages()

        for filename in os.listdir(mydir):
            if filename[:4] == 'tape': os.remove(os.path.join(mydir, filename))
        return DfOutputs, DfMessages



    def get_gendf(self, **fileOptions):
        from ..functions import force_symlink
        kwargs = dict(self.__dict__, **fileOptions)

        mydir = os.path.join(kwargs["FOLDER"], kwargs["FILENAME"])
        os.makedirs(mydir, exist_ok=True)
        kwargs["TIME"] = time.ctime(time.time())

        force_symlink(kwargs["TAPE"], os.path.join(mydir, "tape20"))
        text = self.moder_input(20, -21, **kwargs)

        force_symlink(kwargs["PENDFTAPE"], os.path.join(mydir, "tape29"))
        text += self.moder_input(29, -25, **kwargs)

        text += self.groupr_input(-21, -25, 0, -26, **kwargs)
        text += self.moder_input(-26, 30, **kwargs)
        text += "stop"

        DfOutputs = OutputFiles()
        inputfile = os.path.join(mydir, "input_gendf.{}".format(kwargs["TAG"]))
        DfOutputs = DfOutputs.append({"id" : "input_gendf", "format" : "TEXT", "file" : inputfile}, ignore_index=True)
        with open(inputfile,'w') as f: f.write(text)
        print(" --- run gendf for {} ---".format(kwargs["TAG"]))
        returncode, stdout, stderr = run_process("{} < {}".format(kwargs["EXEC"], inputfile), cwd=mydir)

        oldfile = os.path.join(mydir, "tape30")
        newfile = os.path.join(mydir, "{}.gendf".format(kwargs["TAG"]))
        if os.path.isfile(oldfile):
            if os.path.getsize(oldfile) > 0:
                shutil.move(oldfile, newfile)
                DfOutputs = DfOutputs.append({"id" : "gendf", "format" : "GENDF", "file" : newfile}, ignore_index=True)

        oldfile = os.path.join(mydir, "output")
        newfile = os.path.join(mydir, "output_gendf.{}".format(kwargs["TAG"]))
        shutil.move(oldfile, newfile)
        DfOutputs = DfOutputs.append({"id" : "output_gendf", "format" : "TEXT", "file" : newfile}, ignore_index=True)
        DfMessages = NjoyOutput.from_file(newfile).get_messages()

        for filename in os.listdir(mydir):
            if filename[:4] == 'tape': os.remove(os.path.join(mydir, filename))
        return DfOutputs, DfMessages



    def get_errorr(self, **fileOptions):
        from ..functions import force_symlink
        kwargs = dict(self.__dict__, **fileOptions)
        mydir = os.path.join(kwargs["FOLDER"], kwargs["FILENAME"])
        os.makedirs(mydir, exist_ok=True)
        kwargs["TIME"] = time.ctime(time.time())

        if not os.path.isfile(kwargs["TAPE"]): raise NotImplementedError("file {} not found".format(kwargs["TAPE"]))
        force_symlink(kwargs["TAPE"], os.path.join(mydir, "tape20"))
        text = self.moder_input(20, -21, **kwargs)

        if "GENDFTAPE" in kwargs:
            if not os.path.isfile(kwargs["GENDFTAPE"]): raise NotImplementedError("file {} not found".format(kwargs["GENDFTAPE"]))
            force_symlink(kwargs["GENDFTAPE"], os.path.join(mydir, "tape29"))
            text += self.moder_input(29, -25, **kwargs)
            p = 0; g = -25
        elif "PENDFTAPE" in kwargs:
            if not os.path.isfile(kwargs["PENDFTAPE"]): raise NotImplementedError("file {} not found".format(kwargs["PENDFTAPE"]))
            force_symlink(kwargs["PENDFTAPE"], os.path.join(mydir, "tape29"))
            text += self.moder_input(29, -25, **kwargs)
            p = -25; g = 0

        errorrs = []; er = 30
        for tmp in kwargs["TEMPS"]:
            text += self.errorr_input(-21, p, g, er, TMP=tmp, **kwargs)
            errorrs.append(er)
            er += 1
        text += "stop"

        DfOutputs = OutputFiles()
        inputfile = os.path.join(mydir, "input_errorr.{}".format(kwargs["TAG"]))
        DfOutputs = DfOutputs.append({"id" : "input_errorr", "format" : "TEXT", "file" : inputfile}, ignore_index=True)
        with open(inputfile,'w') as f: f.write(text)
        print(" --- run errorr for {} ---".format(kwargs["TAG"]))
        returncode, stdout, stderr = run_process("{} < {}".format(kwargs["EXEC"], inputfile), cwd=mydir)

        newfile = os.path.join(mydir, "{}.errorr".format(kwargs["TAG"]))
        if os.path.isfile(newfile): os.remove(newfile)
        with open(newfile, "a") as f:
            for er in errorrs:
                oldfile = os.path.join(mydir, "tape{}".format(er))
                with open(oldfile) as g:
                    if not os.path.isfile(oldfile): continue
                    lines = g.readlines()
                    if er != errorrs[0]:
                        lines = lines[1:]
                    if er != errorrs[-1]:
                        lines = lines[:-1]
                    f.writelines(lines)
        if os.path.isfile(newfile):
            if os.path.getsize(newfile) > 0:
                DfOutputs = DfOutputs.append({"id" : "errorr", "format" : "ERRORR", "file" : newfile}, ignore_index=True)

        oldfile = os.path.join(mydir, "output")
        newfile = os.path.join(mydir, "output_errorr.{}".format(kwargs["TAG"]))
        shutil.move(oldfile, newfile)
        DfOutputs = DfOutputs.append({"id" : "output_errorr", "format" : "TEXT", "file" : newfile}, ignore_index=True)
        DfMessages = NjoyOutput.from_file(newfile).get_messages()

        for filename in os.listdir(mydir):
            if filename[:4] == 'tape': os.remove(os.path.join(mydir, filename))

        return DfOutputs, DfMessages



    def get_ace(self, **fileOptions):
        from ..functions import force_symlink
        kwargs = dict(self.__dict__, **fileOptions)
        mydir = os.path.join(kwargs["FOLDER"], kwargs["FILENAME"])
        os.makedirs(mydir, exist_ok=True)

        if not os.path.isfile(kwargs["TAPE"]): raise NotImplementedError("evaluation file {} not found".format(kwargs["TAPE"]))
        force_symlink(kwargs["TAPE"], os.path.join(mydir, "tape20"))
        text = self.moder_input(20, -21, **kwargs)

        if not os.path.isfile(kwargs["PENDFTAPE"]): raise NotImplementedError("file {} not found".format(kwargs["PENDFTAPE"]))
        force_symlink(kwargs["PENDFTAPE"], os.path.join(mydir, "tape29"))
        text += self.moder_input(29, -25, **kwargs)

        aces = []; xsdirs = []
        e = -21; p = -25; a = 40; x = 50
        for tmp,suff in zip(kwargs["TEMPS"], kwargs["SUFFIXES"]):
            text += self.acer_input(e, p, a, x, TMP=tmp, SUFF=suff, **kwargs)
            aces.append(a); xsdirs.append(x)
            a += 1; x += 1
        text += "stop"

        DfOutputs = OutputFiles()
        inputfile = os.path.join(mydir, "input_acer.{}".format(kwargs["TAG"]))
        DfOutputs = DfOutputs.append({"id" : "input_acer", "format" : "TEXT", "file" : inputfile}, ignore_index=True)
        with open(inputfile,'w') as f: f.write(text)
        print(" --- run acer for {} ---".format(kwargs["TAG"]))
        returncode, stdout, stderr = run_process("{} < {}".format(kwargs["EXEC"], inputfile), cwd=mydir)

        newfile = os.path.join(mydir, "{}.ace".format(kwargs["TAG"]))
        if os.path.isfile(newfile): os.remove(newfile)
        with open(newfile, "a") as f:
            for a in aces:
                oldfile = os.path.join(mydir, "tape{}".format(a))
                if not os.path.isfile(oldfile): continue
                with open(oldfile) as g: f.write(g.read())
        if os.path.isfile(newfile):
            if os.path.getsize(newfile) > 0:
                DfOutputs = DfOutputs.append({"id" : "ace", "format" : "ACE", "file" : newfile}, ignore_index=True)

        newfile = os.path.join(mydir, "{}.xsdir".format(kwargs["TAG"]))
        if os.path.isfile(newfile): os.remove(newfile)
        with open(newfile, "a") as f:
            for x in xsdirs:
                oldfile = os.path.join(mydir, "tape{}".format(x))
                if not os.path.isfile(oldfile): continue
                with open(oldfile) as g:
                    xargs = g.read().split()
                    xargs[0] = ".".join([str(kwargs['ZA'] + 300*kwargs["LISO"]), xargs[0].split('.')[1]])
                    xargs[2] = "{}.ace".format(kwargs["TAG"])
                    xargs[3] = "0"
                    f.write(" ".join(xargs))
        if os.path.isfile(newfile):
            if os.path.getsize(newfile) > 0:
                DfOutputs = DfOutputs.append({"id" : "xsdir", "format" : "TEXT", "file" : newfile}, ignore_index=True)

        oldfile = os.path.join(mydir, "output")
        newfile = os.path.join(mydir, "output_acer.{}".format(kwargs["TAG"]))
        shutil.move(oldfile, newfile)
        DfOutputs = DfOutputs.append({"id" : "output_acer", "format" : "TEXT", "file" : newfile}, ignore_index=True)
        DfMessages = NjoyOutput.from_file(newfile).get_messages()

        for filename in os.listdir(mydir):
            if filename[:4] == 'tape': os.remove(os.path.join(mydir, filename))
        return DfOutputs, DfMessages



class evalLib(pd.DataFrame):

    @classmethod
    def from_file(cls, inputfile):
        """
        Populate dataframe using each row of a given inputfile.
        Example of inputfile:
            ENDF-1  [PENDF-1]
            ENDF-2  [PENDF-2]
            ...
            ENDF-N  [PENDF-2]
        Values in brackets [] are optional.
        """
        lines = open(inputfile).read().splitlines()
        evals = []
        for line in lines:
            groups = line.split() # 0-ENDF [1-PENDF]
            E = evalFile(groups[0])
            print("parsing {}".format(groups[0]))
            if len(groups) > 1:
                E.pendfFile = os.path.abspath(os.path.realpath(groups[1]))
            evals.append(E.to_frame())
        frame = pd.concat(evals, sort=False)
        return cls(frame)

    def run_njoy(self, **njoyOptions):
        import multiprocessing as mp
        nj = Njoy(**njoyOptions)
        pool = mp.Pool(processes=nj.processes)
        outs = [pool.apply_async(nj.run, kwds = {**row.to_dict()}) for i,row in self.iterrows()]
        outs = list(map(lambda x:x.get(), outs))
        if nj.run_acer:
            xsdirFiles = [os.path.join(dp, f) for dp, dn, fn in os.walk(os.path.expanduser(nj.evaluationName)) for f in fn if f.endswith(".xsdir")]
            aceFiles = [os.path.join(dp, f) for dp, dn, fn in os.walk(os.path.expanduser(nj.evaluationName)) for f in fn if f.endswith(".ace")]
            mydir = os.path.join(nj.evaluationName, "ace")
            os.makedirs(mydir, exist_ok=True)
            with open(os.path.join(mydir, os.path.basename(nj.evaluationName)+".xsdir"), 'w') as f:
                for file in sorted(xsdirFiles):
                    f.write(open(file).read() + '\n')
            for file in sorted(aceFiles):
                shutil.move(file, os.path.join(mydir, os.path.basename(file)))


def process_lib():
    from .. import settings
    inputs = settings.init_njoy()
    evalLib.from_file(inputs["inputfile"]).run_njoy(**inputs)

def run(iargs=None):
    from .. import settings
    from ..formats import Endf6
    init = settings.init_njoy(iargs)
    nj = Njoy(**vars(init))

    tape = Endf6.from_file(nj.TAPE)
    tape.parse()

    DfOutputs = OutputFiles()
    DfMessages = NjoyMessages()

    if nj.PENDF:
        outs, msgs = nj.get_pendf(**tape.__dict__)
        DfOutputs = pd.concat([DfOutputs, outs]).reset_index(drop=True)
        DfMessages = pd.concat([DfMessages, msgs]).reset_index(drop=True)
    try:
        nj.PENDFTAPE = DfOutputs.query("id==\"pendf\"").file.iloc[0]
    except:
        pass

    if nj.ACE:
        outs, msgs = nj.get_ace(**tape.__dict__)
        DfOutputs = pd.concat([DfOutputs, outs]).reset_index(drop=True)
        DfMessages = pd.concat([DfMessages, msgs]).reset_index(drop=True)

    if nj.GENDF:
        outs, msgs = nj.get_gendf(**tape.__dict__)
        DfOutputs = pd.concat([DfOutputs, outs]).reset_index(drop=True)
        DfMessages = pd.concat([DfMessages, msgs]).reset_index(drop=True)
    try:
        nj.GENDFTAPE = DfOutputs.query("id==\"gendf\"").file.iloc[0]
    except:
        pass

    if nj.ERRORR:
        outs, msgs = nj.get_errorr(**tape.__dict__)
        DfOutputs = pd.concat([DfOutputs, outs]).reset_index(drop=True)
        DfMessages = pd.concat([DfMessages, msgs]).reset_index(drop=True)



#if __name__ == "__main__":
#    from ..formats import Endf6
#    from ..data_test import H1
#
#    nj = Njoy(QA=[46, -1.6651e6, 47, -1.6651e6, 48, -1.6651e6, 49, -1.6651e6], THERMR=False, TEMPS=[300,600], SUFFIXES=[".03", ".06"])
#    tape = Endf6.from_file(os.path.join(H1.__path__[0], r"h1.endf"))
#    tape.parse()
#    nj.PENDF = nj.get_pendf(**tape.__dict__)
#    nj.ACE, nj.XSDIR = nj.get_ace(**tape.__dict__)
#    evalLib.from_file(inputs["inputfile"]).run_njoy(inputs)

#@pytest.mark.njoy
#def test_njoy(tmpdir):
#    from ..data_test import U5
#    iargs = [os.path.join(U5.__path__[0], r"u235.pendf"),
#             "--errorr-cov", os.path.join(U5.__path__[0], r"u235.errorr"),
#             "--QA", "46", "-1.6651e6", "47", "-1.6651e6", "48", "-1.6651e6", "49", "-1.6651e6",
#             "--PK", "301", "302", "303", "304", "318", "402", "442", "443", "444", "445", "446", "447",
#             "--outdir", str(tmpdir),
#             "--processes", str(os.cpu_count()) if os.cpu_count() < 10 else str(10),
#             "--eig", "10",
#             "--samples", "100",]
##             "--plotdir", os.path.join(str(tmpdir), r"html_files"),
##             "-p"]
#    sampling(iargs)