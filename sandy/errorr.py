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

    def get_energy_grid(self, **kwargs):
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
        >>> err = endf6_2.get_errorr(ek_errorr=sandy.energy_grids.CASMO12, err=1, ek_groupr=sandy.energy_grids.CASMO12)
        >>> err.get_energy_grid()
        array([1.0000e-05, 3.0000e-02, 5.8000e-02, 1.4000e-01, 2.8000e-01,
               3.5000e-01, 6.2500e-01, 4.0000e+00, 4.8052e+01, 5.5300e+03,
               8.2100e+05, 2.2310e+06, 1.0000e+07])

        >>> err.get_energy_grid(mat=9443)
        array([1.0000e-05, 3.0000e-02, 5.8000e-02, 1.4000e-01, 2.8000e-01,
               3.5000e-01, 6.2500e-01, 4.0000e+00, 4.8052e+01, 5.5300e+03,
               8.2100e+05, 2.2310e+06, 1.0000e+07])
        """
        mat_ = kwargs.get('mat', self.mat[0])
        mf1 = read_mf1(self, mat_)
        return mf1["EG"]

    def get_xs(self, **kwargs):
        """
        Obtain the xs values across the energy grid.

        Returns
        -------
        xs : `pd.Series`
            For a given mat and mt, the xs values in the energy grid.

        Examples
        --------
        >>> endf6 = sandy.get_endf6_file("jeff_33", "xs", 10010)
        >>> err = endf6.get_errorr(ek_errorr=sandy.energy_grids.CASMO12, err=1)
        >>> err.get_xs()
        MAT                     125                        
        MT                      1           2           102
        E                                                          
        (1e-05, 0.03]           2.10540e+01 2.04363e+01 6.17622e-01
        (0.03, 0.058]           2.06986e+01 2.04363e+01 2.62307e-01
        (0.058, 0.14]           2.06134e+01 2.04363e+01 1.77108e-01
        (0.14, 0.28]            2.05574e+01 2.04363e+01 1.21068e-01
        (0.28, 0.35]            2.05377e+01 2.04363e+01 1.01449e-01
        (0.35, 0.625]           2.05156e+01 2.04363e+01 7.93598e-02
        (0.625, 4.0]            2.04756e+01 2.04360e+01 3.95521e-02
        (4.0, 48.052]           2.04452e+01 2.04328e+01 1.23376e-02
        (48.052, 5530.0]        2.00727e+01 2.00714e+01 1.31829e-03
        (5530.0, 821000.0]      8.05810e+00 8.05804e+00 6.41679e-05
        (821000.0, 2231000.0]   3.48867e+00 3.48863e+00 3.54246e-05
        (2231000.0, 10000000.0] 1.52409e+00 1.52406e+00 3.44005e-05

        >>> err.get_xs(mt=[1, 2])
        MAT                             125            
        MT                                1           2
        E                                              
        (1e-05, 0.03]           2.10540e+01 2.04363e+01
        (0.03, 0.058]           2.06986e+01 2.04363e+01
        (0.058, 0.14]           2.06134e+01 2.04363e+01
        (0.14, 0.28]            2.05574e+01 2.04363e+01
        (0.28, 0.35]            2.05377e+01 2.04363e+01
        (0.35, 0.625]           2.05156e+01 2.04363e+01
        (0.625, 4.0]            2.04756e+01 2.04360e+01
        (4.0, 48.052]           2.04452e+01 2.04328e+01
        (48.052, 5530.0]        2.00727e+01 2.00714e+01
        (5530.0, 821000.0]      8.05810e+00 8.05804e+00
        (821000.0, 2231000.0]   3.48867e+00 3.48863e+00
        (2231000.0, 10000000.0] 1.52409e+00 1.52406e+00

        >>> err.get_xs(mt=1)
        MAT                             125
        MT                                1
        E                                  
        (1e-05, 0.03]           2.10540e+01
        (0.03, 0.058]           2.06986e+01
        (0.058, 0.14]           2.06134e+01
        (0.14, 0.28]            2.05574e+01
        (0.28, 0.35]            2.05377e+01
        (0.35, 0.625]           2.05156e+01
        (0.625, 4.0]            2.04756e+01
        (4.0, 48.052]           2.04452e+01
        (48.052, 5530.0]        2.00727e+01
        (5530.0, 821000.0]      8.05810e+00
        (821000.0, 2231000.0]   3.48867e+00
        (2231000.0, 10000000.0] 1.52409e+00
        """
        data = []
        listmt_ = kwargs.get('mt', range(1, 10000))
        listmt_ = [listmt_] if isinstance(listmt_, int) else listmt_
        listmat_ = kwargs.get('mat', range(1, 10000))
        listmat_ = [listmat_] if isinstance(listmat_, int) else listmat_
        for mat, mf, mt in self.filter_by(listmf=[3],
                                          listmt=listmt_,
                                          listmat=listmat_).data:
            mf1 = sandy.errorr.read_mf1(self, mat)
            egn = pd.IntervalIndex.from_breaks(mf1["EG"])
            mf3 = sandy.errorr.read_mf3(self, mat, mt)
            columns = pd.MultiIndex.from_tuples([(mat, mt)],
                                                names=["MAT", "MT"])
            index = pd.Index(egn, name="E")
            data.append(pd.DataFrame(mf3["XS"], index=index, columns=columns))
        data = pd.concat(data, axis=1).fillna(0)
        return sandy.Xs(data)

    def get_cov(self, multigroup=True):
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
        >>> err = endf6.get_errorr(ek_errorr=[1e-2, 1e1, 2e7], err=1)
        >>> err.get_cov().data
		MAT1	                   125
        MT1	                       1	                            2	                                102
                E1	               (0.01, 10.0]	(10.0, 20000000.0]	(0.01, 10.0]	(10.0, 20000000.0]	(0.01, 10.0]	(10.0, 20000000.0]
        MAT	 MT	                 E						
        125	  1	      (0.01, 10.0]	8.74838e-06	       4.62556e-05	 8.76101e-06	       4.62566e-05	1.07035e-06	           5.58627e-07
                (10.0, 20000000.0]	4.62556e-05	       2.47644e-04	 4.63317e-05	       2.47650e-04	7.58742e-09	           1.49541e-06
              2	      (0.01, 10.0]	8.76101e-06	       4.63317e-05	 8.77542e-06	       4.63327e-05	0.00000e+00	           0.00000e+00
                (10.0, 20000000.0]	4.62566e-05	       2.47650e-04	 4.63327e-05	       2.47655e-04	0.00000e+00	           0.00000e+00
            102	      (0.01, 10.0]	1.07035e-06	       7.58742e-09	 0.00000e+00	       0.00000e+00	6.51764e-04	           3.40163e-04
                (10.0, 20000000.0]	5.58627e-07	       1.49541e-06	 0.00000e+00	       0.00000e+00	3.40163e-04	           6.70431e-02

        >>> err.get_cov(multigroup=False).data
                MAT1    125
                MT1	    1	                                2	                                102
                E1	    1.00000e-02	1.00000e+01	2.00000e+07	1.00000e-02	1.00000e+01	2.00000e+07	1.00000e-02	1.00000e+01	2.00000e+07
        MAT	 MT	E									
        125	  1	        1.00000e-02	8.74838e-06	4.62556e-05	0.00000e+00	8.76101e-06	4.62566e-05	0.00000e+00	1.07035e-06	5.58627e-07	0.00000e+00
                        1.00000e+01	4.62556e-05	2.47644e-04	0.00000e+00	4.63317e-05	2.47650e-04	0.00000e+00	7.58742e-09	1.49541e-06	0.00000e+00
                        2.00000e+07	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00
              2	        1.00000e-02	8.76101e-06	4.63317e-05	0.00000e+00	8.77542e-06	4.63327e-05	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00
                        1.00000e+01	4.62566e-05	2.47650e-04	0.00000e+00	4.63327e-05	2.47655e-04	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00
                        2.00000e+07	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00
            102	        1.00000e-02	1.07035e-06	7.58742e-09	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	6.51764e-04	3.40163e-04	0.00000e+00
                        1.00000e+01	5.58627e-07	1.49541e-06	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	3.40163e-04	6.70431e-02	0.00000e+00
                        2.00000e+07	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00   
        """
        eg = self.get_energy_grid()
        if multigroup:
            eg = pd.IntervalIndex.from_breaks(eg)
        data = []
        for mat, mf, mt in self.filter_by(listmf=[31, 33]).data:
            mf33 = sandy.errorr.read_mf33(self, mat, mt)
            for mt1, cov in mf33["COVS"].items():
                if not multigroup:
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
        return sandy.CategoryCov.from_stack(data, index=["MAT", "MT", "E"],
                                            columns=["MAT1", "MT1", "E1"],
                                            values='VAL')


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
