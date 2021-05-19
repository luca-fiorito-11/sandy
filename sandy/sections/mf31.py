import sandy

__author__ = "Luca Fiorito"
__all__ = [
        "read_mf31",
        ]


def read_mf31(tape, mat, mt):
    return sandy.read_mf33(tape, mat, mt, mf=31)
