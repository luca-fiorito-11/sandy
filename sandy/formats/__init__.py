from .endf6 import *
from .errorr import *
from .groupr import *
from .utils import *
from .records import *
from ..settings import SandyError

def get_file_format(file):
    """Given a file in input return the name of the file format.
    
    Parameters
    ----------
    file : `str`
        input file
    
    Returns
    -------
    `str`
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
        ftype = None
    elif flag == -11 or flag == -12:
        ftype = "errorr"
    elif flag == -1:
        ftype = "gendf"
    else:
        ftype = "endf6"
    return ftype



def read_formatted_file(file, listmat=None, listmf=None, listmt=None):
    """Read formatted file and return either Errorr, Gendf or Endf6 instance
    according to file format.
    
    Returns
    -------
    `Errorr` or `Gendf` or `Endf6`
    """
    ftype = get_file_format(file)
    if ftype is "errorr":
        return Errorr.from_file(file, listmat=listmat, listmf=listmf, listmt=listmt)
    elif ftype is "gendf":
        return Gendf.from_file(file, listmat=listmat, listmf=listmf, listmt=listmt)
    elif ftype is "endf6":
        return Endf6.from_file(file, listmat=listmat, listmf=listmf, listmt=listmt)
    else:
        raise SandyError("file '{}' not in a known format".format(file))