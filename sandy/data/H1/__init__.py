import os

endf6 = open(os.path.join(__path__[0], "h1.endf")).read().splitlines()
pendf = open(os.path.join(__path__[0], "h1.pendf")).read().splitlines()
errorr = open(os.path.join(__path__[0], "h1.errorr")).read().splitlines()