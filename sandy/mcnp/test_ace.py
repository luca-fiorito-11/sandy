# -*- coding: utf-8 -*-
"""
Created on Thu Jul  5 14:12:12 2018

@author: fiorito_l
"""

import pdb, pytest, os, re
import pandas as pd
from collections import namedtuple

xsdirLine = namedtuple('xsdirLine', 'id awr fname route ftype address tlength rlength nr temp ptab')
xsdirLine.__new__.__defaults__ = (None,) * len(xsdirLine._fields)


class OutputFiles(pd.DataFrame):

    def __init__(self, *args, **kwargs):
        kwargs.update({"columns" : ["id", "format", "file"]})
        super().__init__(*args, **kwargs)

class McnpMessages(pd.DataFrame):

    def __init__(self, *args, **kwargs):
        kwargs.update({"columns" : ["type", "message"]})
        super().__init__(*args, **kwargs)

class McnpOutput(str):

    wrn = "\s{2}(?P<type>warning)\.  (?P<msg>.*?)\n"
    err = "\s{}(?P<type>fatal error)\.  (?P<msg>.*?)\n"

    @classmethod
    def from_file(cls, file):
        """
        Read MCNP output file and call from_text method.
        """
        with open(file) as f: text = f.read()
        return cls(text)

    def get_messages(self):
        found = []
        for x in re.findall(self.wrn, self): found.append(x)
        for x in re.findall(self.err, self): found.append(x)
        return McnpMessages(found)


def run_sphere(line, ERG="14", cwd=None, **kwargs):
    from .templates import sphere, empty_xsdir
    from ..functions import force_symlink, run_process

    DfOutputs = OutputFiles()

    xs = xsdirLine(*line.split())
    XS = " ".join([x for x in xs if x is not None])
    ID = xs.id
    AWR = xs.awr
    ZA = xs.id.split(".")[0]
    inp = sphere.replace("ERG", ERG).replace("XS", XS).replace("ID", ID).replace("AWR", AWR).replace("ZA", ZA)

    dirname = os.path.join(xs.id, "test_sphere{}".format(ERG))
    mydir = os.path.join(os.getcwd(), os.path.join(cwd, dirname) if cwd else dirname)
    os.makedirs(mydir, exist_ok=True)
    inpfile = os.path.join(mydir, "inp")
    with open(inpfile, "w") as f: f.write(inp)
    with open(os.path.join(mydir, "xsdir"), "w") as f: f.write(empty_xsdir)

    if "ace_files" in kwargs:
        for ace in kwargs["ace_files"]: force_symlink(ace, os.path.join(mydir, os.path.basename(ace)))

    outfile = os.path.join(mydir, "sphere{}.out".format(ERG))
    returncode, stdout, stderr = run_process("mcnp6 o={}".format(outfile), cwd=mydir, verbose=True)

    DfOutputs = DfOutputs.append({"id" : "output_sphere{}".format(ERG), "format" : "TEXT", "file" : outfile}, ignore_index=True)
    DfOutputs = DfOutputs.append({"id" : "input_sphere{}".format(ERG), "format" : "TEXT", "file" : inpfile}, ignore_index=True)
    DfMessages = McnpOutput.from_file(outfile).get_messages()

    for filename in os.listdir(mydir):
        if filename[:3] == 'run': os.remove(os.path.join(mydir, filename))
        if filename == 'xsdir': os.remove(os.path.join(mydir, filename))
    return DfOutputs, DfMessages

def run(iargs=None):
    from .. import settings
    init = settings.init_test_ace(iargs)
    for line in open(init.xsdir).read().splitlines():
        run_sphere(line, **vars(init))

    pdb.set_trace()


if __name__ == "__main__":
    run()

@pytest.mark.mcnp
@pytest.mark.ace
def test_ace(tmpdir):
    iargs = ["lib/83-Bi-209g.jeff33/83-bi-209g.xsdir",
             "--ace-files", "lib/83-Bi-209g.jeff33/83-bi-209g.ace",]
    run(iargs)