import os

endf6 = open(os.path.join(__path__[0], "u235.endf")).read().splitlines()
pendf = open(os.path.join(__path__[0], "u235.pendf")).read().splitlines()
errorr = open(os.path.join(__path__[0], "u235.errorr")).read().splitlines()
nfy = open(os.path.join(__path__[0], "u235.nfy")).read().splitlines()