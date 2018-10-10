import os

endf6 = open(os.path.join(__path__[0], "fe56.endf")).read().splitlines()
pendf = open(os.path.join(__path__[0], "fe56.pendf")).read().splitlines()
errorr = open(os.path.join(__path__[0], "fe56.errorr")).read().splitlines()