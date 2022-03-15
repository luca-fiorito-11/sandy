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
    Container for groupr information grouped by MAT, MF and MT numbers.
    """
    def get_n_energy_grid(self, mat=None):
        """
        Obtaining neutrons energy grid.

        Parameters
        ----------
        mat : `int`, optional
            MAT number. The default is None.

        Returns
        -------
        `np.array`
            The energy grid of the `sandy.Groupr` object.

        Examples
        --------
        >>> endf6 = sandy.get_endf6_file("jeff_33", "xs", 10010)
        >>> groupr = endf6.get_gendf(verborse=True)
        >>> len(groupr.get_n_energy_grid())
        241
        """
        mat_ = mat if mat else self.mat[0]
        mf1 = read_mf1(self, mat_)
        return mf1["EGN"]

    def get_g_energy_grid(self, mat=None):
        """
        Obtain photons energy grid.

        Parameters
        ----------
        mat : `int`, optional
            MAT number. The default is None.

        Returns
        -------
        `np.array`
            The energy grid of the `sandy.Gendf` object.

        Examples
        --------
        >>> endf6 = sandy.get_endf6_file("jeff_33", "xs", 10010)
        >>> groupr = endf6.get_gendf(verborse=True, ep=sandy.energy_grids.CASMO12)
        >>> groupr.get_g_energy_grid()
        array([1.0000e-05, 3.0000e-02, 5.8000e-02, 1.4000e-01, 2.8000e-01,
               3.5000e-01, 6.2500e-01, 4.0000e+00, 4.8052e+01, 5.5300e+03,
               8.2100e+05, 2.2310e+06, 1.0000e+07])
        """
        mat_ = mat if mat else self.mat[0]
        mf1 = read_mf1(self, mat_)
        return mf1["EGG"]

    def get_xs(self):
        """
        Get multigroup cross sections

        Returns
        -------
        xs : `sandy.Xs`
            group to group cross sections

        Examples
        --------
        >>> endf6 = sandy.get_endf6_file("jeff_33", "xs", 10010)
        >>> groupr = endf6.get_gendf(verborse=True, ek=sandy.energy_grids.CASMO12)
        >>> groupr.get_xs()
        MAT                     125                        
        MT                      1           2           102
        E                                                          
        (1e-05, 0.03]           4.74500e+01 4.68507e+01 5.99276e-01
        (0.03, 0.058]           2.66592e+01 2.64039e+01 2.55277e-01
        (0.058, 0.14]           2.33852e+01 2.32133e+01 1.71860e-01
        (0.14, 0.28]            2.18356e+01 2.17186e+01 1.17013e-01
        (0.28, 0.35]            2.13559e+01 2.12616e+01 9.43025e-02
        (0.35, 0.625]           2.10611e+01 2.09845e+01 7.66054e-02
        (0.625, 4.0]            2.06169e+01 2.05790e+01 3.79424e-02
        (4.0, 48.052]           2.04594e+01 2.04475e+01 1.18527e-02
        (48.052, 5530.0]        2.00729e+01 2.00716e+01 1.28270e-03
        (5530.0, 821000.0]      8.05819e+00 8.05812e+00 6.41591e-05
        (821000.0, 2231000.0]   3.48869e+00 3.48866e+00 3.54245e-05
        (2231000.0, 10000000.0] 1.52409e+00 1.52406e+00 3.44005e-05
        """
        data = []
        for mat, mf, mt in self.filter_by(listmf=[3]).data:
            mf1 = sandy.groupr.read_mf1(self, mat)
            egn = pd.IntervalIndex.from_breaks(mf1["EGN"])
            mf3 = sandy.groupr.read_mf3(self, mat, mt)
            xs = np.array([x["DATA"][1].tolist() for x in mf3["GROUPS"]]).T
            columns = pd.MultiIndex.from_tuples([(mat, mt)],
                                                names=["MAT", "MT"])
            index = pd.Index(egn, name="E")
            data.append(pd.DataFrame(xs, index=index, columns=columns))
        data = pd.concat(data, axis=1).fillna(0)
        return sandy.Xs(data)

    def get_flux(self):
        """
        The flux of the multigroup approach

        Returns
        -------
        flux : `pd.DataFrame`
            Dataframe containing the flux.

        Examples
        --------
        >>> endf6 = sandy.get_endf6_file("jeff_33", "xs", 10010)
        >>> groupr = endf6.get_gendf(verborse=True, ek=sandy.energy_grids.CASMO12)
        >>> groupr.get_flux()
        MAT	                    125
        MT	                    1	        2	        102
        E			
        (1e-05, 0.03]	        2.99900e-02	2.99900e-02	2.99900e-02
        (0.03, 0.058]	        2.80000e-02	2.80000e-02	2.80000e-02
        (0.058, 0.14]	        8.20000e-02	8.20000e-02	8.20000e-02
        (0.14, 0.28]	        1.40000e-01	1.40000e-01	1.40000e-01
        (0.28, 0.35]	        7.00000e-02	7.00000e-02	7.00000e-02
        (0.35, 0.625]	        2.75000e-01	2.75000e-01	2.75000e-01
        (0.625, 4.0]	        3.37500e+00	3.37500e+00	3.37500e+00
        (4.0, 48.052]	        4.40520e+01	4.40520e+01	4.40520e+01
        (48.052, 5530.0]	    5.48195e+03	5.48195e+03	5.48195e+03
        (5530.0, 821000.0]	    8.15470e+05	8.15470e+05	8.15470e+05
        (821000.0, 2231000.0]	1.41000e+06	1.41000e+06	1.41000e+06
        (2231000.0, 10000000.0]	7.76900e+06	7.76900e+06	7.76900e+06
        """
        data = []
        for mat, mf, mt in self.filter_by(listmf=[3]).data:
            mf1 = sandy.groupr.read_mf1(self, mat)
            egn = pd.IntervalIndex.from_breaks(mf1["EGN"])
            mf3 = sandy.groupr.read_mf3(self, mat, mt)
            flux = np.array([x["DATA"][0].tolist() for x in mf3["GROUPS"]]).T
            columns = pd.MultiIndex.from_tuples([(mat, mt)], names=["MAT", "MT"])
            index = pd.Index(egn, name= "E")
            data.append(pd.DataFrame(flux, index=index, columns=columns))
        flux = pd.concat(data, axis=1)
        return flux


def read_mf1(tape, mat):
    """
    Parse MAT/MF=1/MT=451 section from `sandy.Gendf` object and return
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

    Examples
    --------
    >>> endf6 = sandy.get_endf6_file("jeff_33", "xs", 10010)
    >>> groupr = endf6.get_gendf(verborse=True, ek=sandy.energy_grids.CASMO12)
    >>> sandy.groupr.read_mf1(groupr, 125)
    {'MAT': 125,
     'MF': 1,
     'MT': 451,
     'ZA': 1001.0,
     'AWR': 0.9991673,
     'LRP': -1,
     'TEMPIN': 293.6,
     'TITLE': [0.0],
     'SIGZ': [10000000000.0],
     'EGN': array([1.0000e-05, 3.0000e-02, 5.8000e-02, 1.4000e-01, 2.8000e-01,
            3.5000e-01, 6.2500e-01, 4.0000e+00, 4.8052e+01, 5.5300e+03,
            8.2100e+05, 2.2310e+06, 1.0000e+07]),
     'EGG': array([0.])}
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

    Examples
    --------
    >>> endf6 = sandy.get_endf6_file("jeff_33", "xs", 10010)
    >>> groupr = endf6.get_gendf(verborse=True, ek=sandy.energy_grids.CASMO12)
    >>> sandy.groupr.read_mf3(groupr, 125, 1)['GROUPS'][0]
    {'TEMPIN': 293.6,
     'NG2': 2,
     'IG2LO': 1,
     'IG': 1,
     'DATA': array([2.999000e-02, 4.745001e+01])}
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
