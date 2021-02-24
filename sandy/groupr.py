import pandas as pd
import numpy as np

import sandy
from sandy.core.endf6 import _FormattedFile

__author__ = "Luca Fiorito"
__all__ = [
        "Groupr",
        ]

pd.options.display.float_format = '{:.5e}'.format


class Groupr(_FormattedFile):
    """
    Container for ERRORR file text grouped by MAT, MF and MT numbers.
    """

    def get_xs(self, mat, mt):
        """

        Parameters
        ----------
        mat : TYPE
            DESCRIPTION.
        mt : TYPE
            DESCRIPTION.

        Returns
        -------
        xs : TYPE
            DESCRIPTION.

        """
        mf1 = read_mf1(self, mat)
        mf3 = read_mf3(self, mat, mt)
        data = np.array([x["DATA"].tolist() for x in mf3["GROUPS"]]).T
        index = pd.IntervalIndex.from_breaks(mf1["EGN"])
        flux = pd.Series(data[0], index=index, name="iwt")
        xs = pd.Series(data[1], index=index, name=(mat, mt))
        return flux, xs


def read_mf1(tape, mat):
    """
    Parse MAT/MF=1/MT=451 section from `sandy.Errorr` object and return
    structured content in nested dcitionaries.

    Parameters
    ----------
    tape : `sandy.Errorr`
        endf6 object containing requested section
    mat : `int`
        MAT number
    mt : `int`
        MT number

    Returns
    -------
    out : `dict`
        Content of the ENDF-6 tape structured as nested `dict`.
    """
    mf = 1
    mt = 451
    df = tape._get_section_df(mat, mf, mt)
    out = {
            "MAT": mat,
            "MF": mf,
            "MT": mt,
            }
    i = 0
    C, i = sandy.read_cont(df, i)
    NZ = C.L2
    NTW = C.N2
    add = {
        "ZA": C.C1,
        "AWR": C.C2,
        "LRP": C.N1,
    }
    out.update(add)
    L, i = sandy.read_list(df, i)
    NGN = L.L1
    NGG = L.L2
    out["TEMPIN"] = L.C1
    out["TITLE"] = L.B[:NTW]
    del L.B[:NTW]
    out["SIGZ"] = L.B[:NZ]
    del L.B[:NZ]
    out["EGN"] = np.array(L.B[:NGN + 1])
    del L.B[:NGN + 1]
    out["EGG"] = np.array(L.B[:NGG + 1])
    return out


def read_mf3(tape, mat, mt):
    """
    Parse MAT/MF=33/MT section from `sandy.Errorr` object and return
    structured content in nested dcitionaries.

    Parameters
    ----------
    tape : `sandy.Errorr`
        endf6 object containing requested section
    mat : `int`
        MAT number
    mt : `int`
        MT number

    Returns
    -------
    out : `dict`
        Content of the ENDF-6 tape structured as nested `dict`.
    """
    mf = 3
    df = tape._get_section_df(mat, mf, mt)
    out = {
            "MAT": mat,
            "MF": mf,
            "MT": mt,
            }
    i = 0
    C, i = sandy.read_cont(df, i)
    NGN = C.N2
    add = {
        "ZA": C.C1,
        "AWR": C.C2,
        "NL": C.L1,
        "LRFLAG": C.N1,
    }
    out.update(add)
    groups = []
    for ig in range(NGN):
        L, i = sandy.read_list(df, i)
        add = {
            "TEMPIN": L.C1,
            "NG2": L.L1,
            "IG2LO": L.L2,
            "IG": L.N2,
            "DATA": np.array(L.B),
            }
        groups.append(add)
    out["GROUPS"] = groups
    return out
