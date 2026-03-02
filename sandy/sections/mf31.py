
__author__ = "Luca Fiorito"


def read_mf31(tape, mat, mt):
    from .mf33 import read_mf33
    return read_mf33(tape, mat, mt, mf=31)
