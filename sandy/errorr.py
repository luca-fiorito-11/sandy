import pandas as pd
import numpy as np

import sandy
from sandy.core.endf6 import _FormattedFile

__author__ = "Luca Fiorito"
__all__ = [
        "Errorr",
        ]

pd.options.display.float_format = '{:.5e}'.format


class Errorr(_FormattedFile):
    """
    Container for ERRORR file text grouped by MAT, MF and MT numbers.
    """
    
    def get_energy_grid(self, mat=None):
        mat_ = mat if mat else self.mat[0]
        mf1 = read_mf1(self, mat_)
        return mf1["EG"]

        

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
        index = pd.IntervalIndex.from_breaks(mf1["EG"])
        xs = pd.Series(mf3["XS"], index=index, name=(mat, mt))
        return xs
    
    def get_cov(self):
        """
        Extract cross section/nubar covariance from `Errorr` instance.
        
        Returns
        -------
        `sandy CategoryCov`
            xs/nubar covariance matrix for all cross section/nubar
            MAT/MT in ERRORR file
        """
        eg = self.get_energy_grid()
        data = []
        for mat, mf, mt in self.filter_by(listmf=[31, 33]).data:
            mf33 = sandy.errorr.read_mf33(self, mat, mt)
            for mt1, cov in mf33["COVS"].items():
                # add zero row and column at the end of the matrix
                # (this must be done for ERRORR covariance matrices)
                cov = np.insert(cov, cov.shape[0], [0]*cov.shape[1], axis=0)
                cov = np.insert(cov, cov.shape[1], [0]*cov.shape[0], axis=1)
                idx = pd.MultiIndex.from_product(
                    [[mat], [mt], eg],
                    names=["MAT", "MT", "E"],
                    )
                idx1 = pd.MultiIndex.from_product(
                    [[mat], [mt1], eg],
                    names=["MAT1", "MT1", "E1"],
                    )
                df = pd.DataFrame(cov, index=idx, columns=idx1) \
                       .stack(level=["MAT1", "MT1", "E1"]) \
                       .rename("VAL") \
                       .reset_index()
                data.append(df)
        data = pd.concat(data).pivot_table(
            index=["MAT", "MT", "E"],
            columns=["MAT1", "MT1", "E1"],
            values="VAL",
            )
        return data


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
    add = {
        "ZA": C.C1,
        "AWR": C.C2,
        "LRP": C.N1,
    }
    out.update(add)
    L, i = sandy.read_list(df, i)
    add = {
        "EG": np.array(L.B),
    }
    out.update(add)
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
    L, i = sandy.read_list(df, i)
    add = {
        "XS": np.array(L.B),
    }
    out.update(add)
    return out


def read_mf33(tape, mat, mt):
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
    mf = 33
    df = tape._get_section_df(mat, mf, mt)
    out = {
            "MAT": mat,
            "MF": mf,
            "MT": mt,
            }
    i = 0
    C, i = sandy.read_cont(df, i)
    add = {
        "ZA": C.C1,
        "AWR": C.C2,
    }
    out.update(add)
    reaction_pairs = {}
    for rp in range(C.N2):  # number of reaction pairs
        C, i = sandy.read_cont(df, i)
        MT1 = C.L2
        NG = C.N2
        M = np.zeros((NG, NG))
        while True:
            L, i = sandy.read_list(df, i)
            NGCOL = L.L1
            GROW = L.N2
            GCOL = L.L2
            M[GROW-1, GCOL-1:GCOL+NGCOL-1] = L.B
            if GCOL+NGCOL >= NG and GROW >= NG:
                break
        reaction_pairs[MT1] = M
    out["COVS"] = reaction_pairs
    return out
