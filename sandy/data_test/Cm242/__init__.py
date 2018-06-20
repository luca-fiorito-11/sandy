import os

endf6 = open(os.path.join(__path__[0], "cm242.endf")).read().splitlines()
pendf = open(os.path.join(__path__[0], "cm242.pendf")).read().splitlines()