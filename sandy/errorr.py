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
        """
        Obtaining the energy grid. 

        Parameters
        ----------
        mat : `int`, optional
            MAT number. The default is None.

        Returns
        -------
        `np.array`
            The energy grid of the `sandy.Errorr` object.

        Examples
        --------
        >>> endf6_2 = sandy.get_endf6_file("jeff_33", "xs", 942410)
        >>> err = endf6_2.get_errorr(ek=sandy.energy_grids.CASMO12, err=1)
        >>> err.get_energy_grid()
        array([1.0000e-05, 3.0000e-02, 5.8000e-02, 1.4000e-01, 2.8000e-01,
               3.5000e-01, 6.2500e-01, 4.0000e+00, 4.8052e+01, 5.5300e+03,
               8.2100e+05, 2.2310e+06, 1.0000e+07])
        """
        mat_ = mat if mat else self.mat[0]
        mf1 = read_mf1(self, mat_)
        return mf1["EG"]

    def get_xs(self, mat, mt):
        """
        Obtain for a given mat and mt the xs values across the energy grid.

        Parameters
        ----------
        mat : `int`
            MAT number
        mt : `int`
            MT number

        Returns
        -------
        xs : `pd.Series`
            For a given mat and mt, the xs values in the energy grid.

        Examples
        --------
        >>> endf6 = sandy.get_endf6_file("jeff_33", "xs", 10010)
        >>> err = endf6.get_errorr(ek=sandy.energy_grids.CASMO12, err=1)
        >>> err.get_xs(125, 102).head()
        (1e-05, 0.03]   6.17622e-01
        (0.03, 0.058]   2.62307e-01
        (0.058, 0.14]   1.77108e-01
        (0.14, 0.28]    1.21068e-01
        (0.28, 0.35]    1.01449e-01
        Name: (125, 102), dtype: float64
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
        data : `sandy CategoryCov`
            xs/nubar covariance matrix for all cross section/nubar
            MAT/MT in ERRORR file.

        Examples
        --------
        >>> endf6 = sandy.get_endf6_file("jeff_33", "xs", 10010)
        >>> err = endf6.get_errorr(ek=[1e-2, 1e3], err=1)
        >>> err.get_cov().data
        		MAT1	    125
                MT1	        1	                    2	                    102
                E1	        1.00000e-02	1.00000e+03	1.00000e-02	1.00000e+03	1.00000e-02	1.00000e+03
        MAT	MT	          E						
        125	 1	1.00000e-02	8.77247e-06	0.00000e+00	8.77393e-06	0.00000e+00	1.10426e-07	0.00000e+00
                1.00000e+03	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00
             2	1.00000e-02	8.77393e-06	0.00000e+00	8.77542e-06	0.00000e+00	0.00000e+00	0.00000e+00
                1.00000e+03	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00
           102	1.00000e-02	1.10426e-07	0.00000e+00	0.00000e+00	0.00000e+00	6.51764e-04	0.00000e+00
                1.00000e+03	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00
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
        data = pd.concat(data)
        return sandy.CategoryCov.from_stack(data)


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
