import os

endf6 = open(os.path.join(__path__[0], "u238.endf")).read().splitlines()
pendf = open(os.path.join(__path__[0], "u238.pendf")).read().splitlines()