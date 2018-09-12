from .endf6 import *
from .errorr import *
from .groupr import *
from .utils import *
from .records import *
from ..settings import SandyError

def read_formatted_file(file, listmat=None, listmf=None, listmt=None):
    """Read formatted file and return either Errorr, Gendf or Endf6 instance
    according to flag.
    """
    flag = None
    with open(file) as f:
        for line in f.readlines():
            MAT, MF, MT = read_control(line)[:3]
            if MF == 1 and MT == 451:
                i = 0
                C, i = read_cont([line], i)
                flag = C.N1
                break
    if flag is None:
        raise SandyError("file '{}' not in a known format".format(file))
    if flag == -11 or flag == -12:
        return Errorr.from_file(file, listmat=listmat, listmf=listmf, listmt=listmt)
    elif flag == -1:
        return Gendf.from_file(file, listmat=listmat, listmf=listmf, listmt=listmt)
    else:
        return Endf6.from_file(file, listmat=listmat, listmf=listmf, listmt=listmt)