import os

endf6 = open(os.path.join(__path__[0], "pu239.endf")).read().splitlines()
pendf = open(os.path.join(__path__[0], "pu239.pendf")).read().splitlines()