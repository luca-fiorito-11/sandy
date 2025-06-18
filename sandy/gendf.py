import pandas as pd
import numpy as np

from .endf6 import _FormattedFile
from .xs import Xs
from .records import read_cont, read_list

__author__ = "Luca Fiorito"
__all__ = [
        "Gendf",
        ]

pd.options.display.float_format = '{:.5e}'.format


class Gendf(_FormattedFile):
    """
    Container for gendf information grouped by MAT, MF and MT numbers.
    """

    def get_n_energy_grid(self, **kwargs):
        """
        Obtaining neutrons energy grid.

        Returns
        -------
        `np.array`
            The energy grid of the :obj:`~sandy.gendf.Gendf` object.

        Examples
        --------

        >>> import sandy
        >>> endf6 = sandy.get_endf6_file("jeff_33", "xs", 10010)
        >>> gendf = endf6.get_gendf(verborse=True)
        >>> assert len(gendf.get_n_energy_grid()) == 241

        >>> endf6 = sandy.get_endf6_file("jeff_33", "xs", 10010)
        >>> gendf = endf6.get_gendf(groupr_kws=dict(ek=sandy.energy_grids.CASMO12))
        >>> np.testing.assert_allclose(gendf.get_n_energy_grid(), sandy.energy_grids.CASMO12, atol=1e-14, rtol=1e-14)
        >>> np.testing.assert_allclose(gendf.get_n_energy_grid(mat=125), sandy.energy_grids.CASMO12, atol=1e-14, rtol=1e-14)
        """
        mat_ = kwargs.get('mat', self.mat[0])
        mf1 = read_mf1(self, mat_)
        return mf1["EGN"]

    def get_g_energy_grid(self, **kwargs):
        """
        Obtain photons energy grid.

        Returns
        -------
        `np.array`
            The energy grid of the :obj:`~sandy.gendf.Gendf` object.

        Examples
        --------

        >>> import sandy
        >>> endf6 = sandy.get_endf6_file("jeff_33", "xs", 10010)
        >>> gendf = endf6.get_gendf(groupr_kws=dict(ep=sandy.energy_grids.CASMO12))
        >>> np.testing.assert_allclose(gendf.get_g_energy_grid(), sandy.energy_grids.CASMO12, atol=1e-14, rtol=1e-14)
        >>> np.testing.assert_allclose(gendf.get_g_energy_grid(mat=125), sandy.energy_grids.CASMO12, atol=1e-14, rtol=1e-14)
        """
        mat_ = kwargs.get('mat', self.mat[0])
        mf1 = read_mf1(self, mat_)
        return mf1["EGG"]

    def get_xs(self, **kwargs):
        """
        Get multigroup neutron cross sections

        Returns
        -------
        xs : :obj:`~sandy.xs.Xs`
            multigroup cross sections

        Examples
        --------

        >>> import sandy
        >>> endf6 = sandy.get_endf6_file("jeff_33", "xs", 10010)
        >>> gendf = endf6.get_gendf(minimal_processing=True, err=0.005, temperature=293.6, groupr_kws=dict(ek=sandy.energy_grids.CASMO12))
        >>> gendf.get_xs()
        MAT                             125
        MT                              1           2           102         251
        E
        (1e-05, 0.03]           4.74500e+01 4.68507e+01 5.99276e-01 6.67348e-01
        (0.03, 0.058]           2.66592e+01 2.64039e+01 2.55277e-01 6.67289e-01
        (0.058, 0.14]           2.33852e+01 2.32133e+01 1.71860e-01 6.67323e-01
        (0.14, 0.28]            2.18356e+01 2.17186e+01 1.17013e-01 6.67259e-01
        (0.28, 0.35]            2.13559e+01 2.12616e+01 9.43025e-02 6.67231e-01
        (0.35, 0.625]           2.10611e+01 2.09845e+01 7.66054e-02 6.67220e-01
        (0.625, 4.0]            2.06169e+01 2.05790e+01 3.79424e-02 6.67215e-01
        (4.0, 48.052]           2.04594e+01 2.04475e+01 1.18527e-02 6.67220e-01
        (48.052, 5530.0]        2.00729e+01 2.00716e+01 1.28270e-03 6.67237e-01
        (5530.0, 821000.0]      8.05819e+00 8.05812e+00 6.41591e-05 6.67120e-01
        (821000.0, 2231000.0]   3.48869e+00 3.48866e+00 3.54245e-05 6.66838e-01
        (2231000.0, 10000000.0] 1.52409e+00 1.52406e+00 3.44005e-05 6.65044e-01

        >>> gendf.get_xs(mt=1)
        MAT                             125
        MT                                1
        E                                  
        (1e-05, 0.03]           4.74500e+01
        (0.03, 0.058]           2.66592e+01
        (0.058, 0.14]           2.33852e+01
        (0.14, 0.28]            2.18356e+01
        (0.28, 0.35]            2.13559e+01
        (0.35, 0.625]           2.10611e+01
        (0.625, 4.0]            2.06169e+01
        (4.0, 48.052]           2.04594e+01
        (48.052, 5530.0]        2.00729e+01
        (5530.0, 821000.0]      8.05819e+00
        (821000.0, 2231000.0]   3.48869e+00
        (2231000.0, 10000000.0] 1.52409e+00

        >>> gendf.get_xs(mt=[1, 2])
        MAT                             125            
        MT                                1           2
        E                                              
        (1e-05, 0.03]           4.74500e+01 4.68507e+01
        (0.03, 0.058]           2.66592e+01 2.64039e+01
        (0.058, 0.14]           2.33852e+01 2.32133e+01
        (0.14, 0.28]            2.18356e+01 2.17186e+01
        (0.28, 0.35]            2.13559e+01 2.12616e+01
        (0.35, 0.625]           2.10611e+01 2.09845e+01
        (0.625, 4.0]            2.06169e+01 2.05790e+01
        (4.0, 48.052]           2.04594e+01 2.04475e+01
        (48.052, 5530.0]        2.00729e+01 2.00716e+01
        (5530.0, 821000.0]      8.05819e+00 8.05812e+00
        (821000.0, 2231000.0]   3.48869e+00 3.48866e+00
        (2231000.0, 10000000.0] 1.52409e+00 1.52406e+00

        Use `err=1` or else it takes too long.

        >>> endf6 = sandy.get_endf6_file('jeff_33','xs', 922350)
        >>> gendf = endf6.get_gendf(minimal_processing=True, err=1, groupr_kws=dict(ek=sandy.energy_grids.CASMO12))
        >>> gendf.get_xs(mt=[4, 5])
        MAT                            9228            
        MT                                4           5
        E                                              
        (1e-05, 0.03]           0.00000e+00 0.00000e+00
        (0.03, 0.058]           0.00000e+00 0.00000e+00
        (0.058, 0.14]           0.00000e+00 0.00000e+00
        (0.14, 0.28]            0.00000e+00 0.00000e+00
        (0.28, 0.35]            0.00000e+00 0.00000e+00
        (0.35, 0.625]           0.00000e+00 0.00000e+00
        (0.625, 4.0]            0.00000e+00 0.00000e+00
        (4.0, 48.052]           0.00000e+00 0.00000e+00
        (48.052, 5530.0]        8.60022e-07 3.41812e-08
        (5530.0, 821000.0]      1.22176e+00 7.72174e-04
        (821000.0, 2231000.0]   1.94841e+00 1.17248e-02
        (2231000.0, 10000000.0] 1.40637e+00 7.49994e-03
        """
        data = []
        mat_ = kwargs.get('mat', self.mat[0])
        mf1 = read_mf1(self, mat_)
        egn = pd.IntervalIndex.from_breaks(mf1["EGN"])
        listmt_ = kwargs.get('mt', range(1, 10000))
        listmt_ = [listmt_] if isinstance(listmt_, int) else listmt_
        listmat_ = kwargs.get('mat', range(1, 10000))
        listmat_ = [listmat_] if isinstance(listmat_, int) else listmat_
        for mat, mf, mt in self.filter_by(listmf=[3],
                                          listmt=listmt_,
                                          listmat=listmat_).data:
            mf3 = read_mf3(self, mat, mt)
            lowest_range = mf3["GROUPS"][0]["IG"] - 1
            xs = np.array([x["DATA"][1].tolist() for x in mf3["GROUPS"]])
            xs = np.insert(xs, [0]*lowest_range, 0) if lowest_range != 0 else xs
            columns = pd.MultiIndex.from_tuples([(mat, mt)],
                                                names=["MAT", "MT"])
            index = pd.Index(egn, name="E")
            data.append(pd.DataFrame(xs.T, index=index, columns=columns))
        data = pd.concat(data, axis=1).fillna(0)
        return Xs(data)

    def get_flux(self, **kwargs):
        """
        The flux of the multigroup approach

        Returns
        -------
        flux : `pd.Series`
            Dataframe containing the flux.

        Examples
        --------

        >>> import sandy
        >>> endf6 = sandy.get_endf6_file("jeff_33", "xs", 10010)
        >>> gendf = endf6.get_gendf(minimal_processing=True, err=1, temperature=293.6, groupr_kws=dict(ek=sandy.energy_grids.CASMO12))
        >>> gendf.get_flux()
        (1e-05, 0.03]             2.99900e-02
        (0.03, 0.058]             2.80000e-02
        (0.058, 0.14]             8.20000e-02
        (0.14, 0.28]              1.40000e-01
        (0.28, 0.35]              7.00000e-02
        (0.35, 0.625]             2.75000e-01
        (0.625, 4.0]              3.37500e+00
        (4.0, 48.052]             4.40520e+01
        (48.052, 5530.0]          5.48195e+03
        (5530.0, 821000.0]        8.15470e+05
        (821000.0, 2231000.0]     1.41000e+06
        (2231000.0, 10000000.0]   7.76900e+06
        Name: iwt, dtype: float64

        >>> gendf.get_flux(mat=125, mt=2)
        (1e-05, 0.03]             2.99900e-02
        (0.03, 0.058]             2.80000e-02
        (0.058, 0.14]             8.20000e-02
        (0.14, 0.28]              1.40000e-01
        (0.28, 0.35]              7.00000e-02
        (0.35, 0.625]             2.75000e-01
        (0.625, 4.0]              3.37500e+00
        (4.0, 48.052]             4.40520e+01
        (48.052, 5530.0]          5.48195e+03
        (5530.0, 821000.0]        8.15470e+05
        (821000.0, 2231000.0]     1.41000e+06
        (2231000.0, 10000000.0]   7.76900e+06
        Name: iwt, dtype: float64
        """
        data = []
        mat_ = kwargs.get('mat', self.mat[0])
        mt_ = 1
        mf3 = read_mf3(self, mat_, mt_)
        mf1 = read_mf1(self, mat_)
        lowest_range = mf3["GROUPS"][0]["IG"] - 1
        data = np.array([x["DATA"][0].tolist() for x in mf3["GROUPS"]])
        index = pd.IntervalIndex.from_breaks(mf1["EGN"])
        flux = pd.Series(data.T, index=index, name="iwt")
        return flux


def read_mf1(tape, mat):
    """
    Parse MAT/MF=1/MT=451 section from :obj:`~sandy.gendf.Gendf` object and return
    structured content in nested dcitionaries.

    Parameters
    ----------
    tape : :obj:`~sandy.gendf.Gendf`
        ENDF6 object containing requested section.
    mat : `int`
        MAT number.
    mt : `int`
        MT number.

    Returns
    -------
    out : `dict`
        Content of the ENDF-6 tape structured as nested `dict`.

    Examples
    --------

    >>> import sandy
    >>> endf6 = sandy.get_endf6_file("jeff_33", "xs", 10010)
    >>> gendf = endf6.get_gendf(groupr_kws=dict(ek=sandy.energy_grids.CASMO12))
    >>> mf1 = sandy.gendf.read_mf1(gendf, 125)
    >>> mf1['AWR'] = round(mf1['AWR'], 3)
    >>> mf1
    {'MAT': 125,
     'MF': 1,
     'MT': 451,
     'ZA': 1001.0,
     'AWR': 0.999,
     'LRP': -1,
     'TEMPIN': 0.0,
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
    C, i = read_cont(df, i)
    NZ = C.L2
    NTW = C.N2
    add = {
        "ZA": C.C1,
        "AWR": C.C2,
        "LRP": C.N1,
    }
    out.update(add)
    L, i = read_list(df, i)
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
    Parse MAT/MF=33/MT section from :obj:`~sandy.errorr.Errorr` object and return
    structured content in nested dcitionaries.

    Parameters
    ----------
    tape : :obj:`~sandy.gendf.Gendf`
        ENDF6 object containing requested section.
    mat : `int`
        MAT number.
    mt : `int`
        MT number.

    Returns
    -------
    out : `dict`
        Content of the ENDF-6 tape structured as nested `dict`.

    Examples
    --------

    >>> import sandy
    >>> endf6 = sandy.get_endf6_file("jeff_33", "xs", 10010)
    >>> gendf = endf6.get_gendf(temperature=293.6, err=0.005, minimal_processing=True, groupr_kws=dict(ek=sandy.energy_grids.CASMO12))
    >>> sandy.gendf.read_mf3(gendf, 125, 1)['GROUPS'][0]
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
    C, i = read_cont(df, i)
    NGN = C.N2
    add = {
        "ZA": C.C1,
        "AWR": C.C2,
        "NL": C.L1,
        "LRFLAG": C.N1,
    }
    out.update(add)
    groups = []
    L, i = read_list(df, i)
    add = {
            "TEMPIN": L.C1,  # Material temperature (Kelvin)
            "NG2": L.L1,  # Number of secondary positions
            "IG2LO": L.L2,  # Index to lowest zero group
            "IG": L.N2,  # Group index
            "DATA": np.array(L.B),  # Array containing the flux and the xs
            }
    groups.append(add)
    NGN_file = NGN - L.N2
    for ig in range(NGN_file):
        L, i = read_list(df, i)
        add = {
            "TEMPIN": L.C1,  # Material temperature (Kelvin)
            "NG2": L.L1,  # Number of secondary positions
            "IG2LO": L.L2,  # Index to lowest zero group
            "IG": L.N2,  # Group index
            "DATA": np.array(L.B),  # Array containing the flux and the xs
            }
        groups.append(add)
    out["GROUPS"] = groups
    return out
