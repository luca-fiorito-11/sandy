# -*- coding: utf-8 -*-
"""
Created on Fri May 25 16:58:11 2018

@author: fiorito_l
"""
import pandas as pd
import os, sys, time, pdb, shutil, re


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

class evalFile:

    def __init__(self, file):
        def read_float(x):
            try:
                return float(x[0] + x[1:].replace('+', 'E+').replace('-', 'E-'))
            except:
                return x
        widths = [11,11,11,11,11,11,4,2,3]
        columns = ["C1", "C2", "L1", "L2", "N1", "N2","MAT", "MF", "MT"]
        converters = dict(zip(columns[:6],[read_float]*6))
        tape = pd.read_fwf(file, widths=widths, names=columns, converters=converters).query("MAT>0 & MF>0 & MT>0")
        self.filename = os.path.basename(file)
        self.tape = self.endfFile = os.path.abspath(os.path.realpath(file))
        self.endf = self.mat = tape.MAT.iloc[0] # ENDF-6 MAT number
        info = tape.query("MF==1 & MT==451")
        self.za = int(info.C1.iloc[0])
        self.awr = info.C2.iloc[0]
        self.lrp = int(info.L1.iloc[0])
        if info.L2.iloc[0] == 1:
            self.fission = "yes"
        self.nlib = int(info.N1.iloc[0])
        self.mod = int(info.N2.iloc[0])
        self.lis = int(info.L1.iloc[1])
        self.liso = int(info.L2.iloc[1])
        self.metaStable = "true" if self.liso > 0 else "false"
        self.version = int(info.N2.iloc[1]) # Library format
        self.awi = info.C1.iloc[2]
        emax = info.C2.iloc[2]
        self.rel = int(info.L1.iloc[2])
        self.rev = int(info.N2.iloc[2])
        self.nsub = int(info.N1.iloc[2])
        if self.nsub == 10 or self.nsub == 12:
            self.neutron = "yes"
        tag = info.C1.iloc[4]
        try:
            pattern = re.search("(?P<Z>[0-9\s]*?)-(?P<SYM>[\w\s]*?)-(?P<AM>[0-9\s]*)", tag)
            self.tag = self.hmat = pattern.group("SYM").lower().strip() + pattern.group("AM").lower().strip()
        except:
            self.tag = self.hmat = tag.strip()
        if self.liso > 0:
            self.tag += str(self.liso)
        self.lab = info.C2.iloc[4]
        self.eval = info.L1.iloc[4]
        self.author = info.L2.iloc[4]
        self.dist = info.L1.iloc[5]
        self.rdate = info.L2.iloc[5]
        if not tape.query("MF==3 & MT==18").empty:
            self.totalFission = "yes"
#        self.pureAbsorber = "no" # ???
#        self.nis = 1 # ???
#        self.zai = 1 # ???
#        self.gamma = "no"
#        self.AWP0 = "yes"
#        self.dbrcnuclide = info.L1.iloc[4]
#        self.scattering = info.L1.iloc[4]
        if 2 in tape.MF.values:
            self.file2 = "yes"
#            self.resolved = info.L1.iloc[4]
#            self.resonance = info.L1.iloc[4]
#            self.unresolved = info.L1.iloc[4]
            ehRes = tape.query("MF==2 & MT==151").C2.iloc[2]
            if ehRes == emax: ehRes = 0.
            self.ehRes = ehRes
            self.thnmax = - self.ehRes if self.ehRes != 0 else 2.0E6
        if 3 in tape.MF.values:
            self.file3 = "yes"
        if 4 in tape.MF.values:
            self.file4 = "yes"
        if 5 in tape.MF.values:
            self.file5 = "yes"
        if 6 in tape.MF.values:
            self.file6 = "yes"
        if 7 in tape.MF.values:
            self.file7 = "yes"
        if 8 in tape.MF.values:
            self.file8 = "yes"
        if 9 in tape.MF.values or 10 in tape.MF.values:
            self.file9_10 = "yes"
        if 12 in tape.MF.values:
            self.file12 = "yes"
        if list(filter(lambda x:x>30, tape.MF.values)):
            self.covariances = "yes"
        if self.nsub == 12:
            raise NotImplementedError

    def to_frame(self):
        return pd.DataFrame.from_dict(self.__dict__, orient="index").T


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
        self.BROADR = True
        self.THERMR = True
        self.HEATR = False
        self.UNRESR = False
        self.PURR = True
        self.EXEC = "njoy2016"
        self.WDIR = os.getcwd()
        self.FOLDER = os.path.abspath("lib")
        self.IPRINT = 1  # by default set verbosity to max
        self.TEMPS = [293.6]
        # RECONR/BROADR default options
        self.ERR = 0.005
        # PURR/UNRESR default options
        self.SIG0 = [1e10] # default is infinite dilution only
        # THERMR default options
        self.IIN = 1
        self.ICOH = 1
        self.IFORM = 0
        self.NATOM = 1
        self.MTREF = 221
        # ACER default options
        self.SUFFIXES = [".03"]
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
        # SKIP modules
        self.no_thermr = False
        self.no_broadr = False
        self.no_purr = False
        self.no_reconr = False
        self.unresr= False
        self.run_acer = False
        self.run_groupr = False
        self.run_errorr = False
        self.__dict__.update(kwargs)
        # Run checks
        if len(self.TEMPS)> 10: raise NotImplementedError("cannot have more than 10 temperatures")



    def runNjoy(self, inputfile, cwd):
        """
        Run ``NJOY``.
        Pass input file to process's stdin (it must be encoded first).
        ..Important::
            In ``Python 3`` you need to convert string to bytes with a
            ``encode()`` function
        """
        from subprocess import Popen, PIPE
        if self.capsys:
            stderr = stdout = PIPE
        else:
            stderr = stdout = None
        process = Popen(self.EXEC,
                        shell=False,
                        cwd=cwd,
                        stdin=PIPE,
                        stdout=stdout,
                        stderr=stderr)
        inp = open(inputfile).read().encode()
        stdoutdata, stderrdata = process.communicate(inp)
        if process.returncode not in [0, 24]:
            raise NotImplementedError("NJOY exit status {}, cannot run njoy executable".format(process.returncode))

    def runBash(self, inputfile, cwd):
        from subprocess import Popen, PIPE
        if self.capsys:
            stderr = stdout = PIPE
        else:
            stderr = stdout = None
        command = "bash {}".format(inputfile)
        process = Popen(command,
                        shell=True,
                        cwd=cwd,
                        stdin=PIPE,
                        stdout=stdout,
                        stderr=stderr)
        stdoutdata, stderrdata = process.communicate()
        if process.returncode not in [0, 24]:
            raise NotImplementedError("exit status {} when running \"{}\"".format(process.returncode, command))

    @staticmethod
    def moder_input(tapein, tapeout, **kwargs):
        text = """
moder
{} {}
""".format(tapein, tapeout)
        return text

    @staticmethod
    def reconr_input(endfin, pendfout, **kwargs):
        kwargs.update({"ENDFIN" : endfin, "PENDFOUT" : pendfout})
        text = """--
-- *****************************************
-- Reconstruct XS from resonance parameters
-- *****************************************
--
reconr
%(ENDFIN)d %(PENDFOUT)d
'pendf tape from '/
%(MAT)d 1 0/
%(ERR)E  0.  %(ERR)E/
'%(TAG)s with NJOY at %(TIME)s' /
0/
""" % kwargs
        return text

    @staticmethod
    def broadr_input(endfin, pendfin, pendfout, **kwargs):
        kwargs.update({"ENDFIN" : endfin, "PENDFIN" : pendfin, "PENDFOUT" : pendfout})
        text = """--
-- *****************************************
-- Perform Doppler-broadening
-- *****************************************
--
broadr
%(ENDFIN)d %(PENDFIN)d %(PENDFOUT)d
%(MAT)d %(NTEMPS)d 0 0 0./
%(ERR)E %(THNMAX)E %(ERR)E/
%(TEXTTEMPS)s/
0/
""" % kwargs
        return text

    @staticmethod
    def thermr_input(endfin, pendfin, pendfout, **kwargs):
        kwargs.update({"ENDFIN" : endfin, "PENDFIN" : pendfin, "PENDFOUT" : pendfout})
        text = """--
-- *****************************************
-- Add thermal scattering data for free-gas
-- *****************************************
--
thermr
%(ENDFIN)d %(PENDFIN)d %(PENDFOUT)d
0 %(MAT)d 20 %(NTEMPS)d %(IIN)d %(ICOH)d %(IFORM)d %(NATOM)d 221 %(IPRINT)d
%(TEXTTEMPS)s/
0.001 4.0
""" % kwargs
        return text

    @staticmethod
    def purr_input(endfin, pendfin, pendfout, **kwargs):
        kwargs.update({"ENDFIN" : endfin, "PENDFIN" : pendfin, "PENDFOUT" : pendfout})
        text = """
purr
%(ENDFIN)d %(PENDFIN)d %(PENDFOUT)d
%(MAT)d %(NTEMPS)d %(NSIG0)d 20 32/
%(TEXTTEMPS)s/
%(TEXTSIG0)s/
0/
""" % kwargs
        return text

    @staticmethod
    def unresr_input(endfin, pendfin, pendfout, **kwargs):
        kwargs.update({"ENDFIN" : endfin, "PENDFIN" : pendfin, "PENDFOUT" : pendfout})
        text = """
unresr
%(ENDFIN)d %(PENDFIN)d %(PENDFOUT)d
%(MAT)d %(NTEMPS)d %(NSIG0)d 1/
%(TEXTTEMPS)s/
%(TEXTSIG0)s/
0/
""" % kwargs
        return text


    def pendf(self, **fileOptions):
        from ..functions import force_symlink
        kwargs = dict(self.__dict__, **fileOptions)
        print(" --- run pendf for {} ---".format(kwargs["TAG"]))
        if not os.path.isfile(kwargs["TAPE"]): raise NotImplementedError("evaluation file {} not found".format(kwargs["TAPE"]))
        mydir = os.path.join(kwargs["FOLDER"], kwargs["FILENAME"])
        os.makedirs(mydir, exist_ok=True)
        kwargs["NSIG0"] = len(kwargs["SIG0"])
        kwargs["TEXTSIG0"] = " ".join(["%E"%dil for dil in kwargs["SIG0"]])
        kwargs["NTEMPS"] = len(kwargs["TEMPS"])
        kwargs["TEXTTEMPS"] = " ".join(["%E"%tmp for tmp in kwargs["TEMPS"]])
        kwargs["TIME"] = time.ctime(time.time())
        text_bash = "" #"""ln -sf %(endfFile)s tape20\n""" % kwargs
#        text_data = """#!/bin/bash
#set -e
#cat > input_pendf.%(hmat)s << EOF
#moder
#20 -21
#""" % kwargs


        force_symlink(kwargs["TAPE"], os.path.join(mydir, "tape20"))
        text = self.moder_input(20, -21, **kwargs)
        if "PENDF" in kwargs:
#            text_bash += "ln -sf %(pendfFile)s tape30\n" % kwargs
            text += self.moder_input(30, -22, **kwargs)
        text += self.reconr_input(-21, -22, **kwargs)
        e = -21; p = -22
        if kwargs["BROADR"]:
            text += self.broadr_input(e, p, -23, **kwargs)
            p = -23
        if kwargs["THERMR"]:
            text += self.thermr_input(0, p, -24, **kwargs)
            p = -24
        if kwargs["UNRESR"]:
            text += self.unresr_input(e, p, -25, **kwargs)
            p = -25
        if kwargs["PURR"]:
            text += self.purr_input(e, p, -26, **kwargs)
            p = -26
        text += self.moder_input(p, 29, **kwargs)
        text += "stop"

        inputfile = os.path.join(mydir, "input_pendf_{}".format(kwargs["TAG"]))
        with open(inputfile,'w') as f: f.write(text)
        self.runNjoy(inputfile, mydir)
        oldpendf = os.path.join(mydir, "tape29")
        newpendf = os.path.join(mydir, "{}.pendf".format(kwargs["TAG"]))
        if os.path.isfile(oldpendf):
            if os.path.getsize(oldpendf) > 0:
                shutil.move(oldpendf, newpendf)
        shutil.move(os.path.join(mydir, "output"), os.path.join(mydir, "output_pendf.{}".format(kwargs["TAG"])))
        for filename in os.listdir(mydir):
            if filename[:4] == 'tape': os.remove(os.path.join(mydir, filename))
        return newpendf
        text_bash += """
%(NjoyExec)s < input_pendf.%(hmat)s
mv output output_pendf.%(hmat)s
mv tape29 %(filename)s.pendf
rm -f tape*
""" % kwargs
        text_data += text_bash
        inputfile = os.path.join(mydir, "run_pendf_{}.sh".format(kwargs["hmat"]))
        with open(inputfile,'w') as f:
            f.write(text_data)
        self.runBash(inputfile, mydir)
        return os.path.join(mydir, "%(filename)s.pendf" % kwargs)

    def groupr(self, **fileOptions):
        kwargs = dict(self.__dict__, **fileOptions)
        print(" --- run groupr for " + kwargs["hmat"] + " ---")
        if not os.path.isfile(os.path.expandvars(kwargs["endfFile"])): raise NotImplementedError("evaluation file " + kwargs["endfFile"] + " not found")
        if not os.path.isfile(os.path.expandvars(kwargs["pendfFile"])): raise NotImplementedError("evaluation file " + kwargs["pendfFile"] + " not found")
        mydir = os.path.join(kwargs["evaluationName"], kwargs["filename"])
        os.makedirs(mydir, exist_ok=True)
        kwargs["nbDil"] = len(kwargs["sig0"])
        kwargs["textDil"] = " ".join(["%E"%dil for dil in kwargs["sig0"]])
        kwargs["nbTmp"] = len(kwargs["temps"])
        kwargs["textTmp"] = " ".join(["%E"%tmp for tmp in kwargs["temps"]])
        if kwargs["ign"] == -99: # 1 GROUP (1E-5 to 2E7)
            kwargs["ign"] = 1
            kwargs["nstr"] = """1 /
1E-5 2E7 /
"""
        elif kwargs["ign"] == -98: # 1 GROUP (resonance integral)
            kwargs["ign"] = 1
            kwargs["nstr"] = """1 /
5E-1 1e5 /
"""
        kwargs["htime"] = time.ctime(time.time())
        text_bash = """
ln -sf %(endfFile)s tape20
ln -sf %(pendfFile)s tape29

%(NjoyExec)s < input_gendf.%(hmat)s
mv output output_gendf.%(hmat)s
mv tape30 %(filename)s.gendf
rm -f tape*
""" % kwargs
        text_data = """#!/bin/bash
set -e
cat > input_gendf.%(hmat)s << EOF
moder
20 -21
moder
29 -25
groupr
-21 -25 0 -26
%(mat)d %(ign)d %(igg)d %(iwt)d %(lord)d %(nbTmp)d %(nbDil)d %(iprint)d /
'%(hmat)s from %(evaluationName)s at %(htime)s' /
%(textTmp)s /
%(textDil)s /
""" % kwargs
        if kwargs["nstr"]:
            text_data += """%(nstr)s""" % kwargs
        for tmp in kwargs["temps"]:
            text_data += """3/
6/
0/
"""
        text_data += """0/
moder
-26 30
stop
EOF

""" % kwargs
        text_data += text_bash
        inputfile = os.path.join(mydir, "run_gendf_{}.sh".format(kwargs["hmat"]))
        with open(inputfile,'w') as f:
            f.write(text_data)
        self.runBash(inputfile, mydir)
        return os.path.join(mydir, "%(filename)s.gendf" % kwargs)

    def errorr(self, **fileOptions):
        kwargs = dict(self.__dict__, **fileOptions)
        print(" --- run errorr for " + kwargs["hmat"] + " ---")
        if not os.path.isfile(os.path.expandvars(kwargs["endfFile"])): raise NotImplementedError("evaluation file " + kwargs["endfFile"] + " not found")
        mydir = os.path.join(kwargs["evaluationName"], kwargs["filename"])
        os.makedirs(mydir, exist_ok=True)
        if kwargs["igne"] == -99: # 1 GROUP (1E-5 to 2E7)
            kwargs["igne"] = 1
            kwargs["nstr"] = """1 /
1E-5 2E7 /
"""
        elif kwargs["igne"] == -98: # 1 GROUP (resonance integral)
            kwargs["igne"] = 1
            kwargs["nstr"] = """1 /
5E-1 1e5 /
"""
        kwargs["htime"] = time.ctime(time.time())
        for tmp,suff in zip(kwargs["temps"], kwargs["suffixes"]):
            kwargs.update({"tmp":tmp, "suff":suff})
            text_data = """
moder
20 -21
moder
29 -25
"""
            if "gendfFile" in kwargs:
                if not os.path.isfile(os.path.expandvars(kwargs["gendfFile"])): raise NotImplementedError("evaluation file " + kwargs["gendfFile"] + " not found")
                text_data += """
errorr
-21 0 -25 33 /
"""
            elif "pendfFile" in kwargs:
                if not os.path.isfile(os.path.expandvars(kwargs["pendfFile"])): raise NotImplementedError("evaluation file " + kwargs["pendfFile"] + " not found")
                text_data += """
errorr
-21 -25 0 33 /
"""
            else:
                raise NotImplementedError("user must define either pendfFile or gendfFile attribute")
            text_data += """%(mat)d %(igne)d %(iwte)d %(iprint)d %(irelco)d /
%(iprint)d %(tmp)E /
0 33 %(irespr)d %(lord)d /
""" % kwargs
            if kwargs["nstr"]:
                text_data += """%(nstr)s
""" % kwargs
            text_data += """
stop
"""
            text_bash = """#!/bin/bash
set -e
cat > input_errorr%(suff)s.%(hmat)s << EOF
""" % kwargs + text_data + """
EOF

ln -sf %(endfFile)s tape20
""" % kwargs
            if "gendfFile" in kwargs:
                if not os.path.isfile(os.path.expandvars(kwargs["gendfFile"])): raise NotImplementedError("evaluation file " + kwargs["gendfFile"] + " not found")
                text_bash += """
ln -sf %(gendfFile)s tape29
""" % kwargs
            elif "pendfFile" in kwargs:
                if not os.path.isfile(os.path.expandvars(kwargs["pendfFile"])): raise NotImplementedError("evaluation file " + kwargs["pendfFile"] + " not found")
                text_bash += """
ln -sf %(pendfFile)s tape29
""" % kwargs

            text_bash += """
%(NjoyExec)s < input_errorr%(suff)s.%(hmat)s
mv output output_errorr%(suff)s.%(hmat)s
mv tape33 %(filename)s%(suff)s.errorr
rm -f tape*
""" % kwargs
            inputfile = os.path.join(mydir, "run_errorr_{}{}.sh".format(kwargs["hmat"], suff))
            with open(inputfile,'w') as f:
                f.write(text_bash)
            self.runBash(inputfile, mydir)

    def acer(self, **fileOptions):
        kwargs = dict(self.__dict__, **fileOptions)
        print(" --- run acer for " + kwargs["hmat"] + " ---")
        if not os.path.isfile(os.path.expandvars(kwargs["endfFile"])): raise NotImplementedError("evaluation file " + kwargs["endfFile"] + " not found")
        if not os.path.isfile(os.path.expandvars(kwargs["pendfFile"])): raise NotImplementedError("evaluation file " + kwargs["pendfFile"] + " not found")
        mydir = os.path.join(kwargs["evaluationName"], kwargs["filename"])
        kwargs["htime"] = time.ctime(time.time())
        for tmp,suff in zip(kwargs["temps"], kwargs["suffixes"]):
            kwargs.update({"tmp":tmp, "suff":suff})
            text_data = """#!/bin/bash
set -e
cat > input_acer%(suff)s.%(hmat)s << EOF
moder
20 -21
moder
29 -25
acer
-21 -25 0 38 39
1 %(iprint)d 1 %(suff)s 0/
'%(hk)s'/
%(mat)d  %(tmp)E /
%(newfor)d %(iopp)d/
0.001/
-- acer / Check ACE files
-- 0 38 0 40 41
-- 7 1 1 -1/
-- /
stop
EOF

""" % kwargs
            text_data += """
ln -sf %(endfFile)s tape20
ln -sf %(filename)s.pendf tape29
%(NjoyExec)s < input_acer%(suff)s.%(hmat)s
mv output output_acer%(suff)s.%(hmat)s
mv tape38 %(filename)s%(suff)s.ace
mv tape39 %(filename)s%(suff)s.xsdir
rm -f tape*
""" % kwargs
            inputfile = os.path.join(mydir, "run_acer_{}{}.sh".format(kwargs["hmat"], suff))
            with open(inputfile,'w') as f:
                f.write(text_data)
            self.runBash(inputfile, mydir)
            xsdirfile = os.path.join(mydir, kwargs["filename"] + kwargs["suff"] + ".xsdir")
            if os.path.isfile(xsdirfile):
                with open(xsdirfile) as f:
                    xargs = f.read().split()
                xargs[2] = "%(filename)s%(suff)s.ace" % kwargs
                xargs[3] = "0" % kwargs
                xargs[0] = ".".join([str(kwargs['za'] + 300*kwargs["liso"]), xargs[0].split('.')[1]])
                with open(xsdirfile, 'w') as f:
                    f.write(" ".join(xargs))

    def run(self, **fileOptions):
        fileOptions["pendfFile"] = self.pendf(**fileOptions)
        if self.run_acer:
            self.acer(**fileOptions)
        if self.run_groupr:
            self.groupr(**fileOptions)
        if self.run_errorr:
            self.errorr(**fileOptions)



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

if __name__ == "__main__":
#    from sandy.njoy import parser
#    init = settings.init_njoy()
    from ..formats import Endf6
    from ..data_test import H1
    nj = Njoy(BROADR=False, THERMR=False, PURR=False, UNRESR=True)
    tape = Endf6.from_file(os.path.join(H1.__path__[0], r"h1.endf"))
    tape.parse()
    pendf = nj.pendf(**tape.__dict__)
    pdb.set_trace()
#    evalLib.from_file(inputs["inputfile"]).run_njoy(inputs)