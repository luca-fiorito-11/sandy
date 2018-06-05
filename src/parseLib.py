# -*- coding: utf-8 -*-
"""
Created on Fri May 25 16:58:11 2018

@author: fiorito_l
"""
import pandas as pd
import os, re, sys, time, shutil

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



class PyNjoyError(Exception):
    """Exception indicating an error in PyNjoy."""



def TimeDecorator(foo):
    """
    Output the time a function takes
    to execute.
    """
    def wrapper(*args, **kwargs):
        t1 = time.time()
        out = foo(*args, **kwargs)
        t2 = time.time()
        print("Time to run function {}: {} sec".format(foo, t2-t1))
        return out
    return wrapper



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
        self.tape = self.evaluationFile = os.path.abspath(os.path.realpath(file))
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


class PyNjoy:

    def __init__(self, wdir, **kwargs):
#        self.sab = os.path.join(os.path.dirname(os.path.abspath(__file__)), "sab.csv")
#        self.wdir = os.path.normpath(wdir)
#        self.evaluationName = os.path.basename(self.wdir)
        self.wdir = os.getcwd()
        self.evaluationName = os.path.abspath("pynjoy")
        self.iwt = 4
        self.legendre = 1
        self.legendregg = 6
        self.scatteringLaw = None
        self.eFiss = None
        self.branchingNG = None
        self.branchingN2N = None
        self.gstr = 0
        self.oldlib = None
        self.purr = None
        self.sgref = 1.0E10
        self.yields = None
        self.iburn = -1
        self.ip1opt = 0
        self.ifprod = 0
        self.jp1 = 0
        self.iverw = 4
        # General options
        self.iprint = 1  # by default set verbosity to max
        self.temperatures = [293.6]
        # RECONR/BROADR default options
        self.err = 0.005
        # PURR/UNRESR default options
        self.dilutions = [1e10] # default is infinite dilution only
        self.purr = True
        # THERMR default options
        self.iin = 1
        self.icoh = 1
        self.iform = 0
        self.natom = 1
        self.mtref = 221

    def runNjoy(self, inputfile, cwd):
        """
        Run ``NJOY``.
        Pass input file to process's stdin (it must be encoded first).
        ..Important::
            In ``Python 3`` you need to convert string to bytes with a
            ``encode()`` function
        """
        from subprocess import Popen, PIPE
        process = Popen(self.NjoyExec,
                        shell=True,
                        cwd=cwd,
                        stdin=PIPE,
                        stdout=PIPE,
                        stderr=PIPE)
        inp = open(inputfile).read().encode()
        stdoutdata, stderrdata = process.communicate(inp)
        if process.returncode not in [0, 24]:
            raise PyNjoyError("NJOY exit status {}, cannot run njoy executable".format(process.returncode))


    def pendf(self, *args, **kwargs):
        kwargs = dict(self.__dict__, **kwargs)
        print(" --- make pendf for " + kwargs["hmat"] + " ---")
        if not os.path.isfile(os.path.expandvars(kwargs["evaluationFile"])): raise PyNjoyError("evaluation file " + kwargs["evaluationFile"] + " not found")
        mydir = os.path.join(kwargs["evaluationName"], kwargs["filename"])
        if not os.path.isdir(mydir): os.mkdir(mydir)
        if kwargs["dilutions"]:
            kwargs["nbDil"] = len(kwargs["dilutions"])
            kwargs["textDil"] = " ".join(["%E"%dil for dil in kwargs["dilutions"]])
        else:
            kwargs["nbDil"] = 0
            kwargs["textDil"] = ""
        kwargs["nbTmp"] = len(kwargs["temperatures"])
        kwargs["textTmp"] = " ".join(["%E"%tmp for tmp in kwargs["temperatures"]])
        kwargs["htime"] = time.ctime(time.time())

        text_data = """
moder
20 -21
--
-- *********************************************************
-- Reconstruct XS from resonance parameters
-- *********************************************************
--
reconr
-21 -22
'pendf tape from %(evaluationName)s'/
%(mat)d 1 0/
%(err)E  0.  %(err)E/
'%(hmat)s from %(evaluationName)s at %(htime)s' /
0/
--
-- *********************************************************
-- Perform Doppler-broadening
-- *********************************************************
--
broadr
-21 -22 -23
%(mat)d %(nbTmp)d 0 0 0./
%(err)E %(thnmax)E %(err)E/
%(textTmp)s/
0/
""" % kwargs

        if self.scatteringLaw:
            text_data += """
moder
26 -27
--
-- *********************************************************
-- Add thermal scattering data
-- *********************************************************
--
thermr
-27 -23 -25
%(matLaw)d %(mat)d 32 %(nbTmp)d %(typeLaw)d %(elasOpt)d %(iform)d %(nbAtoms)d %(matsab_inc)d 0/
%(textTmp)s/
0.001 4.0
""" % kwargs
        else:
            text_data += """
-- *********************************************************
-- Add thermal scattering data for free-gas
-- *********************************************************
thermr
0 -23 -24
0 %(mat)d 20 %(nbTmp)d %(iin)d %(icoh)d %(iform)d %(natom)d 221 %(iprint)d
%(textTmp)s/
0.001 4.0
""" % kwargs
            if self.dilutions:
                if self.purr:
                    text_data += """
purr
-21 -24 -25
%(mat)d %(nbTmp)d %(nbDil)d 20 32/
%(textTmp)s/
%(textDil)s/
0/
""" % kwargs
                else:
                    text_data += """
unresr
-21 -24 -25
%(mat)d %(nbTmp)d %(nbDil)d 1/
%(textTmp)s/
%(textDil)s/
0/
""" % kwargs
            else:
                text_data += """
moder
-24 -25
""" % kwargs
        text_data += """
moder
-25 29
stop
""" % kwargs
        inputfile = os.path.join(mydir, "njoy.inp")
        with open(inputfile,'w') as f:
            f.write(text_data)
        os.symlink(kwargs["evaluationFile"], os.path.join(mydir, "tape20"))
        self.runNjoy(inputfile, mydir)
#        shutil.move(os.path.join(self.evaluationName, "tape20"), os.path.join(self.evaluationName, "tape20"))
        import pdb
        pdb.set_trace()
        os.system("mv tape29 pendf" + self.hmat)
        os.system("mv file_data file_data_pendf" + self.hmat)
        os.system("mv output out_pendf_" + self.hmat)
        os.system("chmod 644 out_pendf_" + self.hmat)
        for fileName in os.listdir(os.getcwd()):
          if fileName[:4] == 'tape': os.remove(fileName)
        os.chdir(myCwd)



class evalLib(pd.DataFrame):

#    def __init__(self, listFiles, name):
#        self.name = name
#        self.extend([evalFile(file) for file in listFiles])
    @classmethod
    def from_file(cls, file):
        with open(file) as f:
            frame = pd.concat([evalFile(x).to_frame() for x in f.read().splitlines()])
        return cls(frame)

    def to_njoy(self, **kwargs):
        for i,row in self.iterrows():
            if row.nsub == 10:
                text = """from PyNjoy import PyNjoy

lib = PyNjoy()
lib.evaluationName = '%(evaluationName)s'
lib.NjoyExec = '%(NjoyExec)s'

lib.pendf(mat=%(endf)d, hmat='%(tag)s')
lib.acer(mat=%(endf)d, hmat='%(tag)s')
""" % dict(row.to_dict(), **kwargs)
            elif row.nsub == 12:
                sys.exit()
            with open("aaa", 'w') as f: f.write(text)
            aaa=1

lib = PyNjoy("njoy2016")
file = evalFile("1-H-3g.jeff33")
lib.pendf(**file.__dict__)
A = evalLib.from_file("inputs")
A.to_njoy(evaluationName="AAA", NjoyExec="njoy2016")
bbb=1