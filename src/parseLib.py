# -*- coding: utf-8 -*-
"""
Created on Fri May 25 16:58:11 2018

@author: fiorito_l
"""
import pandas as pd
import os, re, sys, time, shutil

class NjoyOptions(dict):

    def __init__(self, **kwargs):
        self.update({"sab" : os.path.join(os.path.dirname(os.path.abspath(__file__)), "sab.csv"),
               "iwt" : 4,
               "legendre" : 1,
               "legendregg" : 6,
               "scatteringLaw" : None,
               "eFiss" : None,
               "branchingNG" : None,
               "branchingN2N" : None,
               "gstr" : 0,
               "oldlib" : None,
               "sgref" : 1.0E10,
               "yields" : None,
               "iburn" : -1,
               "ip1opt" : 0,
               "ifprod" : 0,
               "jp1" : 0,
               "iverw" : 4,
               # General options
               "iprint" : 1,  # by default set verbosity to max
               "temperatures" : [293.6],
               # RECONR/BROADR default options
               "err" : 0.005,
               # PURR/UNRESR default options
               "dilutions" : [1e10], # default is infinite dilution only
               "purr" : True,
               # THERMR default options
               "iin" : 1,
               "icoh" : 1,
               "iform" : 0,
               "natom" : 1,
               "mtref" : 221,
               })
        self.update(kwargs)



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



class evalFile(dict):

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
        self['filename'] = os.path.basename(file)
        self['tape'] = os.path.abspath(os.path.realpath(file))
        self['endf'] = tape.MAT.iloc[0] # ENDF-6 MAT number
        info = tape.query("MF==1 & MT==451")
        self["za"] = int(info.C1.iloc[0])
        self["awr"] = info.C2.iloc[0]
        lrp = int(info.L1.iloc[0])
        if info.L2.iloc[0] == 1:
            self["fission"] = "yes"
        self["nlib"] = int(info.N1.iloc[0])
        self["mod"] = int(info.N2.iloc[0])
        self["lis"] = int(info.L1.iloc[1])
        self["liso"] = int(info.L2.iloc[1])
        self["metaStable"] = "true" if self["liso"] > 0 else "false"
        self["version"] = int(info.N2.iloc[1]) # Library format
        self["awi"] = info.C1.iloc[2]
        emax = info.C2.iloc[2]
        self["rel"] = int(info.L1.iloc[2])
        self["rev"] = int(info.N2.iloc[2])
        self["nsub"] = int(info.N1.iloc[2])
        if self["nsub"] == 10 or self["nsub"] == 12:
            self["neutron"] = "yes"
        tag = info.C1.iloc[4]
        try:
            pattern = re.search("(?P<Z>[0-9\s]*?)-(?P<SYM>[\w\s]*?)-(?P<AM>[0-9\s]*)", tag)
            self["tag"] = pattern.group("SYM").lower().strip() + pattern.group("AM").lower().strip()
        except:
            self["tag"] = tag.strip()
        if self["liso"] > 0:
            self["tag"] += str(self["liso"])
        self["lab"] = info.C2.iloc[4]
        self["eval"] = info.L1.iloc[4]
        self["author"] = info.L2.iloc[4]
        self["dist"] = info.L1.iloc[5]
        self["rdate"] = info.L2.iloc[5]
        if not tape.query("MF==3 & MT==18").empty:
            self["totalFission"] = "yes"
        self["pureAbsorber"] = "no" # ???
        self["nis"] = 1 # ???
        self["zai"] = 1 # ???
        self["gamma"] = "no"
#        self["AWP0"] = "yes"
#        self["dbrcnuclide"] = info.L1.iloc[4]
#        self["scattering"] = info.L1.iloc[4]
        if 2 in tape.MF.values:
            self["file2"] = "yes"
#            self["resolved"] = info.L1.iloc[4]
#            self["resonance"] = info.L1.iloc[4]
#            self["unresolved"] = info.L1.iloc[4]
            ehRes = tape.query("MF==2 & MT==151").C2.iloc[2]
            if ehRes == emax: ehRes = 0.
            self["ehRes"] = ehRes
        if 3 in tape.MF.values:
            self["file3"] = "yes"
        if 4 in tape.MF.values:
            self["file4"] = "yes"
        if 5 in tape.MF.values:
            self["file5"] = "yes"
        if 6 in tape.MF.values:
            self["file6"] = "yes"
        if 7 in tape.MF.values:
            self["file7"] = "yes"
        if 8 in tape.MF.values:
            self["file8"] = "yes"
        if 9 in tape.MF.values or 10 in tape.MF.values:
            self["file9_10"] = "yes"
        if 12 in tape.MF.values:
            self["file12"] = "yes"
        if list(filter(lambda x:x>30, tape.MF.values)):
            self["covariances"] = "yes"

    def to_frame(self):
        return pd.DataFrame.from_dict(self, orient="index").T

    def pendf(self, **kwargs):
        options = dict(self.__dict__, **kwargs)
        print(" --- make pendf for (%tag)s ---" % options)
        myCwd = os.getcwd()
        if not os.path.isfile(os.path.expandvars(self["filename"]): raise PyNjoyError("evaluation file " + self["filename"] + " not found")
        if not os.path.isdir(self.evaluationName): os.mkdir(self.evaluationName)
        textDil=""
        if self.dilutions:
            nbDil = len(self.dilutions)
            textDil = " ".join(["%E"%dil for dil in self.dilutions])
        else:
            nbDil = 0
        nbTmp = len(self.temperatures)
        textTmp = " ".join(["%E"%tmp for tmp in self.temperatures])

        htime = time.ctime(time.time())
        kwargs.update({"textDil": textDil,
                       "nbDil"  : nbDil,
                       "textTmp": textTmp,
                       "nbTmp"  : nbTmp,
                       "htime"  : htime})
#                              "unitLaw": unitLaw, \
#                              "matLaw" : matLaw,  \
#                              "typeLaw": typeLaw, \
#                              "iform"  : iform,   \
#                              "nbAtoms": nbAtoms, \
#                              "elasOpt": elasOpt, \
#                              "matsab_inc" : matsab_inc})
        kwargs = dict(self.__dict__, **kwargs)

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
        %(err)E/
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
        inputfile = os.path.join(self.evaluationName, "file_data")
        with open(inputfile,'w') as f:
            f.write(text_data)
        evalfile = os.path.abspath(self.evaluationFile)
        outdir = os.path.join(self.evaluationName, self.hmat)
        os.symlink(evalfile, os.path.join(outdir, "tape20"))
        self.runNjoy(inputfile, outdir)
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


A = evalLib.from_file("inputs")
A.to_njoy(evaluationName="AAA", NjoyExec="njoy2016")
bbb=1