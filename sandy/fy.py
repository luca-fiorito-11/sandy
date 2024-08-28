r"""
This module contains all classes and functions specific for processing fission
yield data.

Examples
--------
Get CEA fission yield evaluations and correlation matrices.

>>> import os, sandy
>>> assert os.path.exists(sandy.fy_cea_u235th)
>>> assert os.path.exists(sandy.fy_cea_pu239th)
>>> assert os.path.exists(sandy.fy_cea_u235th_corr)
>>> assert os.path.exists(sandy.fy_cea_pu239th_corr)
"""
import logging

import pandas as pd
import numpy as np
import scipy.sparse as sps
from os.path import join, dirname
import re

from .zam import ELEMENTS
from .core.cov import CategoryCov
from .core.endf6 import Endf6
from .sections.mf8 import write_mf8
from .gls import _gls_parameters_update, ishikawa_factor
from .shared import expand_zam

__author__ = "Luca Fiorito"
__all__ = [
        "Fy",
        "fy_cea_pu239th",
        "fy_cea_pu239th_corr",
        "fy_cea_u235th",
        "fy_cea_u235th_corr",
        ]


fy_cea_pu239th = join(dirname(__file__), 'appendix', 'fission_yields', r"jeff-4t3_cea_pu9_cons_28-09-2023.stn")
fy_cea_pu239th_corr = join(dirname(__file__), 'appendix', 'fission_yields', r"jeff-4t3_cea_pu9th_cons_28-09-2023_ind_corr")
fy_cea_u235th = join(dirname(__file__), 'appendix', 'fission_yields', r"mixt_cea-jeff33_u_235_th_eval_c1.stn")
fy_cea_u235th_corr = join(dirname(__file__), 'appendix', 'fission_yields', r"mixt_cea-jeff33_u_235_th_ind_corr_mat_c1")


minimal_fytest = pd.DataFrame(
    [[9437, 454, 942390, 380900, 0.0253e-5, 0.1 * 2, 0.02],
     [9437, 454, 942390, 551370, 0.0253e-5, 0.9 * 2, 0.09],
     [9437, 454, 942390, 380900, 500e3, 0.4 * 2, 0.04],
     [9437, 454, 942390, 551370, 500e3, 0.5 * 2, 0.05],
     [9437, 454, 942390, 541350, 500e3, 0.1 * 2, 0.01]],
    columns=["MAT", "MT", "ZAM", "ZAP", "E", "FY", "DFY"]
    )
minimal_fytest_2 = pd.DataFrame([
     [9437, 454, 942390, 591480, 500e3, 0.1, 0.04],
     [9437, 454, 942390, 591481, 500e3, 0.1 * 2, 0.05],
     [9437, 454, 942390, 601480, 500e3, 0.1 * 3, 0.01],
     [9437, 459, 942390, 591480, 500e3, 0.4 * 2, 0.04],
     [9437, 459, 942390, 591481, 500e3, 0.5 * 2, 0.05],
     [9437, 459, 942390, 601480, 500e3, 0.1 * 2, 0.01]],
    columns=["MAT", "MT", "ZAM", "ZAP", "E", "FY", "DFY"]
    )

def get_chain_yields():
    r"""
    Import chain yields information from data stored in sandy. The
    information was taken from 'https://www-nds.iaea.org/endf349/la-ur-94-3106.pdf',
    page 18-29.

    Returns
    -------
    df : `pd.Dataframe`
        Information of the url divided into a dataframe.

    Notes
    -----
    ..note :: It is recognized that even smaller independent yields have even
    larger than 100% error. Indeed, independent yields in the range of 1.0e-9 to
    1.0e-12 may be difficult to predict to better than a factor of 100. For
    this reason, all yields less than 1.0e-12 are blanked out.

    Examples
    --------

    >>> import sandy
    >>> chain_yields = sandy.fy.get_chain_yields()
    >>> chain_yields.head()
        A     ZAM        E         CHY        DCHY
    0  66  902270  thermal 8.53000e-10 2.72960e-10
    1  67  902270  thermal 2.16000e-09 6.91200e-10
    2  68  902270  thermal 8.23000e-09 2.63360e-09
    3  69  902270  thermal 2.35000e-08 7.52000e-09
    4  70  902270  thermal 5.19000e-08 1.66080e-08
    """
    errors = {'a': 0.0035, 'b': 0.0050, 'c': 0.0070, 'd': 0.01, 'e': 0.014,
               'f': 0.02, 'g': 0.028, 'h': 0.04,
               'i': 0.06, 'j': 0.08, 'k': 0.11, 'l': 0.16, 'm': 0.23,
               'n': 0.32, 'o': 0.45, 'p': 0.64}
    energy = {'t': "thermal", 'f': "fast", 'h': "high energy",
              's': 'spontaneous fission'}
    files = ['appendix A.txt', 'appendix B.txt', 'appendix C.txt',
             'appendix D.txt', 'appendix E.txt', 'appendix F.txt']
    #
    path = join(dirname(__file__), 'appendix', 'chain_yields')
    df = pd.concat([pd.read_csv(join(path, file), sep=r"\s+", index_col=0) for file in files], axis=1)
    df.columns.name, df.index.name = "ISO", "A"
    df = df.stack().rename("Y").reset_index("ISO")
    # 
    zam_pattern = re.compile("(?P<SYM>[a-z]+)(?P<A>[0-9]+)(?P<M>[m]*)(?P<E>[a-z]+)", flags=re.IGNORECASE)
    df["SYM"] = df.ISO.apply(lambda x: zam_pattern.search(x).group("SYM").title())
    df["A"] = df.ISO.apply(lambda x: zam_pattern.search(x).group("A")).astype(int)
    df["M"] = df.ISO.apply(lambda x: zam_pattern.search(x).group("M")).astype(bool).astype(int)
    df["E"] = df.ISO.apply(lambda x: energy[zam_pattern.search(x).group("E")])
    df["Z"] = df.SYM.apply(lambda x: {v: k for k, v in ELEMENTS.items()}[x])
    df["ZAM"] = df.Z * 10000 + df.A * 10 + df.M
    #
    df["CHY"] = df.Y.apply(lambda x: x[:-1]).astype(float) / 100
    df["DCHY"] = df.Y.apply(lambda x: errors[x[-1]]) * df.CHY
    return df[["ZAM", "E", "CHY", "DCHY"]].sort_values(by=["ZAM", "E"]) \
                                          .reset_index()


class Fy():
    r"""
    Object for fission yield data.

    Attributes
    ----------
    data : `pandas.DataFrame`
        source of fission yield data

    Methods
    -------
    custom_perturbation
        Apply a custom perturbation to a given fission yield
    energy_table
        Interpolate cross sections over new grid structure
    _expand_zam
        Add columns `Z`, `A` and `M` of the fissioning isotope
    _expand_zap
        Add columns `Z`, `A` and `M` of the fission product
    filter_by
        Apply condition to fission yield data
    from_endf6
         Extract fission yields from `Endf6` instance
    to_endf6
        Update fission yield data in `Endf6` instance
    to_hdf5
        Write fission yield data to hdf5 file
    """

    _columns = ["MAT", "MT", "ZAM", "ZAP", "E", "FY", "DFY"]

    def __repr__(self):
        return self.data.__repr__()

    def __init__(self, df, **kwargs):
        self.data = pd.DataFrame(df, **kwargs)

    @property
    def data(self):
        r"""
        Dataframe of fission yield data with the following columns:

            - `MAT` : MAT number
            - `MT` : MT number
            - `ZAM` : `Z*1000 + A*10 + M` for the parent (fissioning) nuclide
            - `ZAP` : `Z*1000 + A*10 + M` for the daughter nuclide
                      (fission product)
            - `E` : fissioning energy
            - `FY` : fission yield (fraction)

        Returns
        -------
        `pandas.DataFrame`
            tabulated fission yields

        Examples
        --------
        >>> Fy(minimal_fytest)
            MAT   MT     ZAM     ZAP           E          FY         DFY
        0  9437  454  942390  380900 2.53000e-07 2.00000e-01 2.00000e-02
        1  9437  454  942390  551370 2.53000e-07 1.80000e+00 9.00000e-02
        2  9437  454  942390  380900 5.00000e+05 8.00000e-01 4.00000e-02
        3  9437  454  942390  551370 5.00000e+05 1.00000e+00 5.00000e-02
        4  9437  454  942390  541350 5.00000e+05 2.00000e-01 1.00000e-02
        """
        return self._data

    @data.setter
    def data(self, data):
        self._data = data[self._columns]

    def energy_table(self, key, by="ZAM", kind="independent"):
        r"""
        Pivot dataframe of tabulated fission yields as a function of energy.
        Columns are determined by keyword argument `'by'`.

        Parameters
        ----------
        key : `int`
            MAT or ZAM number used as columns of the pivot table
        by : `str`, optional, default is `'ZAM'`
            column name used as columns in the pivot table
        kind : `str`, optional, default is `'independent'`
            either `'independent'` or `'cumulative'`

        Returns
        -------
        `pandas.DataFrame`
            tabulated fission yields

        Notes
        -----
        ..note :: if a fission product has a yield at say, thermal energy but
                  not at fast, then a fast yield of zero will be assigned.

        Raises
        ------
        `ValueError`
            if `kind` is neither `'independent'` nor `'cumulative'`

        Examples
        --------
        >>> Fy(minimal_fytest).energy_table(942390)
        ZAP              380900      541350      551370
        E                                              
        2.53000e-07 2.00000e-01 0.00000e+00 1.80000e+00
        5.00000e+05 8.00000e-01 2.00000e-01 1.00000e+00
        """
        df = self.data
        if kind == "independent":
            mt = 454
        elif kind == "cumulative":
            mt = 459
        else:
            msg = "`kind` must be either `'independent'` or `'cumulative'`"
            raise ValueError(msg)
        condition = (df[by] == key) & (df.MT == mt)
        # if i use fill_value in pivot_table it replaces small values
        # (e.g. 2.88210e-12) with zero. Use fillna!!!
        return pd.pivot_table(
                    df[condition],
                    index="E",
                    columns="ZAP",
                    values="FY",
                    ).fillna(0.)

    def _expand_zap(self):
        r"""
        Produce dataframe with three extra columns containing the `Z`, `A` and
        `M` numbers of the **parent** (fissioning) nuclide.

        Returns
        -------
        `pandas.DataFrame`
            dataframe with Z, A and M columns.

        >>> Fy(minimal_fytest)._expand_zap()
             MAT 	 MT 	   ZAM 	   ZAP 	          E 	         FY 	        DFY 	 Z 	 A 	    M
        0 	9437 	454 	942390 	380900 	2.53000e-07 	2.00000e-01 	2.00000e-02 	38 	90 	    0
        1 	9437 	454 	942390 	551370 	2.53000e-07 	1.80000e+00 	9.00000e-02 	55 	137 	0
        2 	9437 	454 	942390 	380900 	5.00000e+05 	8.00000e-01 	4.00000e-02 	38 	90 	    0
        3 	9437 	454 	942390 	551370 	5.00000e+05 	1.00000e+00 	5.00000e-02 	55 	137 	0
        4 	9437 	454 	942390 	541350 	5.00000e+05 	2.00000e-01 	1.00000e-02 	54 	135 	0
        """
        zam = pd.DataFrame(map(expand_zam, self.data.ZAP),
                           columns=["Z", "A", "M"],
                           dtype=int)
        return self.data.assign(Z=zam.Z, A=zam.A, M=zam.M)

    def _expand_zam(self):
        r"""
        Produce dataframe with three extra columns containing the `Z`, `A` and
        `M` numbers of the **daughter** nuclide (fission product).

        Returns
        -------
        `pandas.DataFrame`
            dataframe with Z, A and M columns.

        >>> Fy(minimal_fytest)._expand_zam()
            MAT 	MT 	ZAM 	ZAP 	E 	FY 	DFY 	Z 	A 	M
        0 	9437 	454 	942390 	380900 	2.53000e-07 	2.00000e-01 	2.00000e-02 	94 	239 	0
        1 	9437 	454 	942390 	551370 	2.53000e-07 	1.80000e+00 	9.00000e-02 	94 	239 	0
        2 	9437 	454 	942390 	380900 	5.00000e+05 	8.00000e-01 	4.00000e-02 	94 	239 	0
        3 	9437 	454 	942390 	551370 	5.00000e+05 	1.00000e+00 	5.00000e-02 	94 	239 	0
        4 	9437 	454 	942390 	541350 	5.00000e+05 	2.00000e-01 	1.00000e-02 	94 	239 	0
        """
        zam = pd.DataFrame(map(expand_zam, self.data.ZAM),
                           columns=["Z", "A", "M"],
                           dtype=int)
        return self.data.assign(Z=zam.Z, A=zam.A, M=zam.M)

    def get_mass_yield(self, zam, e):
        r"""
        Obtain mass yields from the following model: ChY = S * IFY

        Parameters
        ----------
        zam : `int`
            ZAM number of the fissioning nuclide.
        e : `float`
            Energy of the fissioning system.

        Returns
        -------
        `pandas.Series`
            mass yield obtained from ChY = S * IFY

        Examples
        --------
        
        >>> import sandy
        >>> tape_nfpy = sandy.get_endf6_file("jeff_33",'nfpy', 922350)
        >>> nfpy = Fy.from_endf6(tape_nfpy)
        >>> out = nfpy.get_mass_yield(922350, 0.0253).loc[148]
        >>> np.testing.assert_almost_equal(out, 0.0169029147)
        """
        # Filter FY data:
        conditions = {'ZAM': zam, "E": e, 'MT': 454}
        fy_data = self._filters(conditions).data.set_index('ZAP').FY
        mass_data = self._filters({'ZAM': zam, "E": e})
        S = mass_data.get_mass_yield_sensitivity()
        mass_yield = S.dot(fy_data)
        return mass_yield.rename('mass yield')

    def get_chain_yield(self, zam, e, decay_data, **kwargs):
        r"""
        Obtain chain yields from the following model: ChY = S * IFY

        Parameters
        ----------
        zam : `int`
            ZAM number of the fissioning nuclide.
        e : `float`
            Energy of the fissioning system.
        decay_data : :obj:`~sandy.decay.DecayData`
            Radioactive nuclide data from where to obtain chain sensitivities.
        kwargs : `dict`
            Keyword arguments for method :obj:`~sandy.decay.DecayData.get_decay_chains`.

        Returns
        -------
        `pd.Series`
            Chain yield obtained from ChY = S * IFY

        Examples
        --------

        >>> import sandy
        >>> zam = [591480, 591481, 601480, 561480, 571480, 571490, 581480]
        >>> decay_minimal = sandy.get_endf6_file("jeff_33", 'decay', zam)
        >>> decay_fytest = sandy.DecayData.from_endf6(decay_minimal)
        >>> tape_nfpy = sandy.get_endf6_file("jeff_33", 'nfpy', 922350)
        >>> nfpy = Fy.from_endf6(tape_nfpy)
        >>> result_value  = float(nfpy.get_chain_yield(922350, 0.0253, decay_fytest).loc[148])  # Convert to native Python float
        >>> assert result_value == 0.01692277272
        """
        # Filter FY data:
        conditions = {'ZAM': zam, "E": e, 'MT': 454}
        fy_data = self._filters(conditions).data.set_index('ZAP').FY
        S = decay_data. get_chain_yield_sensitivity(**kwargs)
        S = S.reindex(columns=fy_data.index).fillna(0)
        chain_yield = S.dot(fy_data)
        return chain_yield.rename('chain yield')

    def get_mass_yield_sensitivity(self):
        r"""
        Obtain the mass yield sensitivity matrix based only on the
        information given in the `Fy` object (no decay data).

        Returns
        -------
        mass_yield_sensitivity : `pandas.DataFrame`
            Mass yield sensitivity matrix. Mass number in index and ZAM of
            fission products in columns.

        Examples
        --------

        >>> import sandy
        >>> zap =pd.Index([551480, 551490, 561480, 561490, 571480, 571490, 581480, 591480, 591481, 601480])
        >>> tape_nfpy = sandy.get_endf6_file("jeff_33",'nfpy','all')
        >>> zam = 922350
        >>> energy = 0.0253
        >>> conditions = {'ZAM': zam, 'E': energy}
        >>> nfpy = Fy.from_endf6(tape_nfpy)._filters(conditions)
        >>> nfpy.get_mass_yield_sensitivity().loc[148, zap]
        551480   1.00000e+00
        551490   0.00000e+00
        561480   1.00000e+00
        561490   0.00000e+00
        571480   1.00000e+00
        571490   0.00000e+00
        581480   1.00000e+00
        591480   1.00000e+00
        591481   1.00000e+00
        601480   1.00000e+00
        Name: 148, dtype: float64
        """
        # Filter FY data:
        fy_data = self.filter_by('MT', 454)._expand_zap()\
                      .set_index('A')[['ZAP']]
        # Create mass yield sensitivity
        groups = fy_data.groupby(fy_data.index)['ZAP'].value_counts().rename("COUNT")
        return groups.reset_index().pivot_table(index='A', columns='ZAP', values="COUNT", aggfunc="sum").fillna(0)

    def custom_perturbation(self, zam, mt, e, pert):
        r"""
        Apply a custom perturbation to a given fission yield.

        Parameters
        ----------
        zam : `int`
            ZAM number of the material to which perturbations are to be
            applied.
        mt : `int`
            MT reaction number of the FY to which perturbations are to be
            applied either 454 or 459.
        e : `float`
            Energy of the fissioning system to which which perturbations
            are to be applied.
        pert : `pd.Series`
            Perturbation coefficients as ratio values. `pert` index should be
            the corrisponding ZAP numbers.

        Returns
        -------
        `Fy`
            Fission yield instance with applied perturbation.

        Examples
        --------

        >>> import sandy
        >>> tape = sandy.get_endf6_file("jeff_33", 'nfpy', 'all')
        >>> nfpy = Fy.from_endf6(tape)
        >>> pert = pd.Series([0.9], index=[551370])
        >>> nfpy_pert = nfpy.custom_perturbation(922350, 459, 0.0253, pert)
        >>> comp = nfpy_pert.data.query('ZAM==922350 & ZAP==551370 & MT==459 & E==0.0253').squeeze().FY
        >>> assert np.setdiff1d(nfpy_pert.data.values, nfpy.data.values) == comp

        >>> pert = pd.Series([0.9, 0.5], index=[541350, 551370])
        >>> nfpy_pert = nfpy.custom_perturbation(922350, 459, 0.0253, pert)
        >>> comp1 = nfpy_pert.data.query('ZAM==922350 & ZAP==551370 & MT==459 & E==0.0253').squeeze().FY
        >>> f1 = nfpy.data.query('ZAM==922350 & ZAP==551370 & MT==459 & E==0.0253').squeeze().FY
        >>> assert f1 * 0.5 == comp1
        >>> comp2 = nfpy_pert.data.query('ZAM==922350 & ZAP==541350 & MT==459 & E==0.0253').squeeze().FY
        >>> f2 = nfpy.data.query('ZAM==922350 & ZAP==541350 & MT==459 & E==0.0253').squeeze().FY
        >>> assert f2 * 0.9 == comp2
        """
        df = self.data.copy()
        idx = df.query(f"ZAM=={zam} & ZAP=={list(pert.index)} & MT=={mt} & E=={e}").index
        fy_pert = df.loc[idx, :].set_index("ZAP").FY * pert
        fy_pert.index = idx
        fy_pert.name = 'FY'
        df.update(fy_pert)
        return self.__class__(df)

    def apply_bmatrix(self, zam, energy, decay_data, keep_fy_index=False):
        r"""
        Perform IFY = (1-B) * CFY equation to calculate IFY in a given zam
        for a given energy and apply into the original data.

        Parameters
        ----------
        zam : `int`
            ZAM number of the material to which calculations are to be
            applied.
        energy : `float`
            Energy to which calculations are to be applied.
        decay_data : :obj:`~sandy.decay.DecayData`
            Radioactive nuclide data for several isotopes.
        keep_fy_index: `bool`, optional, default is `False`
            Option that allows you to output only the CFY results that were
            part of the original :obj:`~sandy.fy.Fy` object. The default is False.

        Returns
        -------
        :obj:`~sandy.fy.Fy`
            Fission yield instance with IFY calculated for a given combination
            of ZAM/e/decay_data.

        Notes
        -----
        .. note:: This method applies a perturbation to certain IFYs,
        since the equation IFY = (1-B) * CFY is not satisfied for all nuclei.

        Examples
        --------

        >>> import sandy
        >>> zam = [591480, 591481, 601480]
        >>> decay_minimal = sandy.get_endf6_file("jeff_33", 'decay', zam)
        >>> decay_fytest = sandy.DecayData.from_endf6(decay_minimal)
        >>> npfy = Fy(minimal_fytest_2)
        >>> npfy_pert = npfy.apply_bmatrix(942390, 5.00000e+05, decay_fytest)
        >>> npfy_pert.data.query("MT==454")
            MAT   MT     ZAM     ZAP           E           FY         DFY
        3  9437  454  942390  591480 5.00000e+05  8.00000e-01 4.00000e-02
        4  9437  454  942390  591481 5.00000e+05  1.00000e+00 5.00000e-02
        5  9437  454  942390  601480 5.00000e+05 -1.60000e+00 6.48074e-02
        6  9437  454  942390  621480 5.00000e+05 -2.00000e-01 1.00000e-02
        
        >>> zam = [591480, 591481, 601480]
        >>> decay_minimal = sandy.get_endf6_file("jeff_33", 'decay', zam)
        >>> decay_fytest = sandy.DecayData.from_endf6(decay_minimal)
        >>> npfy = Fy(minimal_fytest_2)
        >>> npfy_pert = npfy.apply_bmatrix(942390, 5.00000e+05, decay_fytest, keep_fy_index=True)
        >>> npfy_pert.data.query("MT==454")
            MAT   MT     ZAM     ZAP           E           FY         DFY
        3  9437  454  942390  591480 5.00000e+05  8.00000e-01 4.00000e-02
        4  9437  454  942390  591481 5.00000e+05  1.00000e+00 5.00000e-02
        5  9437  454  942390  601480 5.00000e+05 -1.60000e+00 6.48074e-02
        """
        # Obtain the data:
        data = self.data.copy()
        conditions = {'ZAM': zam, 'MT': 459, "E": energy}
        fy_data = self._filters(conditions).data
        mat = fy_data.MAT.iloc[0]
        std = fy_data.set_index('ZAP')['DFY']
        fy_data = fy_data.set_index('ZAP')['FY']
        if keep_fy_index:
            original_index = fy_data.index

        # Get B-matrix
        B = decay_data.get_bmatrix()
        index, columns = B.index, B.columns

        # Creating (1-B) matrix (which is the sensitivity)
        B = sps.csc_matrix(B)
        unit = sps.csc_matrix(sps.identity(B.shape[0]))
        S = pd.DataFrame(
            (unit - B).toarray(),
            index=index,
            columns=columns,
            )  # sensitivity matrix

        # Rest of the data
        mask = (data.ZAM == zam) & (data.MT == 454) & (data.E == energy)
        fy_data = fy_data.reindex(columns).fillna(0)
        data = data.loc[~mask]
        fy_data = fy_data.reindex(columns).fillna(0)
        std = std.reindex(columns).fillna(0)
        cov_data = CategoryCov.from_stdev(std)

        # Apply (1-B) matrix
        ify_calc_values = (S @ fy_data).rename('FY')
        cov_calc_values = np.sqrt(np.diag((cov_data.sandwich(S).data)))
        cov_calc_values = pd.Series(cov_calc_values, index=S.index)
        if keep_fy_index:
            ify_calc_values = ify_calc_values.reindex(original_index).fillna(0)
            cov_calc_values = cov_calc_values.reindex(original_index).fillna(0)
        calc_values = ify_calc_values.reset_index().rename(columns={'DAUGHTER': 'ZAP'})
        calc_values['DFY'] = cov_calc_values.values
        calc_values = calc_values.assign(MAT=mat, ZAM=zam, MT=454, E=energy)
        data = pd.concat([data, calc_values], ignore_index=True)
        return self.__class__(data)

    def apply_qmatrix(self, zam, energy, decay_data, cut_hl=True, keep_fy_index=False):
        r"""
        Perform CFY = Q * IFY equation to calculate CFY in a given zam
        for a given energy and apply into the original data.

        Parameters
        ----------
        zam : `int`
            ZAM number of the material to which calculations are to be
            applied.
        e : `float`
            Energy to which calculations are to be applied.
        decay_data : :obj:`~sandy.decay.DecayData`
            Radioactive nuclide data for several isotopes.
        cut_hl: `bool`, optional, default is `True`
            cut all the decay modes of the nuclides with a half life larger
            than 100 years.
        keep_fy_index : `bool`, optional, default is `False`
            Option that allows you to output only the CFY results that were
            part of the original :obj:`~sandy.fy.Fy` object.

        Returns
        -------
        :obj:`~sandy.fy.Fy`
            Fission yield instance with IFY calculated for a given combination
            of ZAM/e/decay_data.

        Notes
        -----
        .. note:: This method applies a perturbation to certain CFYs,
        since the equation CFY = Q * IFY is not satisfied for all nuclei.

        Examples
        --------

        >>> import sandy
        >>> zam = [591480, 591481, 601480]
        >>> decay_minimal = sandy.get_endf6_file("jeff_33", 'decay', zam)
        >>> decay_fytest = sandy.DecayData.from_endf6(decay_minimal)
        >>> npfy = Fy(minimal_fytest_2)
        >>> npfy_pert = npfy.apply_qmatrix(942390, 5.00000e+05, decay_fytest, cut_hl=False)
        >>> npfy_pert.data.query("MT==459")
            MAT   MT     ZAM     ZAP           E          FY         DFY
        3  9437  459  942390  591480 5.00000e+05 1.00000e-01 4.00000e-02
        4  9437  459  942390  591481 5.00000e+05 2.00000e-01 5.00000e-02
        5  9437  459  942390  601480 5.00000e+05 6.00000e-01 6.48074e-02
        6  9437  459  942390  621480 5.00000e+05 6.00000e-01 6.48074e-02

        >>> npfy_pert = npfy.apply_qmatrix(942390, 5.00000e+05, decay_fytest, cut_hl=True)
        >>> npfy_pert.data.query("MT==459")
            MAT   MT     ZAM     ZAP           E          FY         DFY
        3  9437  459  942390  591480 5.00000e+05 1.00000e-01 4.00000e-02
        4  9437  459  942390  591481 5.00000e+05 2.00000e-01 5.00000e-02
        5  9437  459  942390  601480 5.00000e+05 6.00000e-01 6.48074e-02
        6  9437  459  942390  621480 5.00000e+05 0.00000e+00 0.00000e+00

        >>> zam = [591480, 591481, 601480]
        >>> decay_minimal = sandy.get_endf6_file("jeff_33", 'decay', zam)
        >>> decay_fytest = sandy.DecayData.from_endf6(decay_minimal)
        >>> npfy = Fy(minimal_fytest_2)
        >>> npfy_pert = npfy.apply_qmatrix(942390, 5.00000e+05, decay_fytest, keep_fy_index=True)
        >>> npfy_pert.data.query("MT==459")
            MAT   MT     ZAM     ZAP           E          FY         DFY
        3  9437  459  942390  591480 5.00000e+05 1.00000e-01 4.00000e-02
        4  9437  459  942390  591481 5.00000e+05 2.00000e-01 5.00000e-02
        5  9437  459  942390  601480 5.00000e+05 6.00000e-01 6.48074e-02
        """
        # Obtain the data:
        data = self.data.copy()
        conditions = {'ZAM': zam, 'MT': 454, "E": energy}
        fy_data = self._filters(conditions).data
        mat = fy_data.MAT.iloc[0]
        std = fy_data.set_index('ZAP')['DFY']
        fy_data = fy_data.set_index('ZAP')['FY']
        if keep_fy_index:
            original_index = fy_data.index
        Q = decay_data.get_qmatrix(cut_hl=cut_hl)

        # Put the data in a approppiate format
        mask = (data.ZAM == zam) & (data.MT == 459) & (data.E == energy)
        data = data.loc[~mask]
        fy_data = fy_data.reindex(Q.columns).fillna(0)
        std = std.reindex(Q.columns).fillna(0)
        cov_data = CategoryCov.from_stdev(std)

        # Apply qmatrix
        cfy_calc_values = (Q @ fy_data).rename('FY')
        cov_calc_values = np.sqrt(np.diag(cov_data.sandwich(Q).data))
        cov_calc_values = pd.Series(cov_calc_values, index=Q.index)
        if keep_fy_index:
            cfy_calc_values = cfy_calc_values.reindex(original_index).fillna(0)
            cov_calc_values = cov_calc_values.reindex(original_index).fillna(0)
        calc_values = cfy_calc_values.reset_index().rename(columns={'DAUGHTER': 'ZAP'})
        calc_values['DFY'] = cov_calc_values.values
        calc_values = calc_values.assign(MAT=mat, ZAM=zam, MT=459, E=energy)
        data = pd.concat([data, calc_values], ignore_index=True)
        return self.__class__(data)
    
    def gls_update(self, zam, energy, S, y_extra, Vy_extra=None):
        r"""
        Perform the GLS update of fission yields and their related covariance
        matrix, according with
        https://doi.org/10.1016/j.anucene.2015.10.027.
        .. math::
            $$
            IFY_{post} = IFY_{prior} + V_{IFY_{prior}}\cdot S^T \cdot \left(S\cdot V_{IFY_{prior}}\cdot S^T + V_{y_{extra}}\right)^{-1} \cdot \left(y_{extra} - y_{calc}\right)\\
            V_{IFY_{post}} = V_{IFY_{prior}} - V_{IFY_{prior}}\cdot S^T \cdot \left(S\cdot V_{IFY_{prior}}\cdot S^T + V_{y_{extra}}\right)^{-1} \cdot S \cdot V_{IFY_{prior}}
            $$
            with $x_{post}$ the updated IFY and $x_{prior}$ the prior IFY

        Parameters
        ----------
        zam : `int`
            ZAM number of the material to which calculations are to be
            applied.
        energy : `float`
            Energy to which calculations are to be applied.
        S : `pandas.DataFrame`
            Sensitivity matrix (M, N) or sensitivity vector (N,).
        y_extra : `pandas.Series`
            Extra information on output (M,). The names of the indeces must be
            consistent with the names of the indices of the sensitivity matrix.
        Vy_extra : `pandas.DataFrame` or `pandas.Series`, optional, default is `None`
            covariance matrix with the uncertainties of the extra information,
            (M, M) or (1,). The names of the indeces and columns must be
            consistent with the names of the indices of the sensitivity matrix.

        Returns
        -------
        :obj:`~sandy.fy.Fy`
            `Fy` object updated with GLS for a given zam, energy, design
            sensitivity and new information.
        :obj:`~sandy.core.cov.CategoryCov`
            :obj:`~sandy.core.cov.CategoryCov` object corresponding to the updated covariance matrix
            adjusted with the GLS technique.

        Notes
        -----
        .. note:: If Vy_extra=None the constraint GLS update technique
        will be performed.

        Examples
        --------

        >>> import sandy
        >>> zam = [591480, 591481, 601480]
        >>> decay_minimal = sandy.get_endf6_file("jeff_33", 'decay', zam)
        >>> decay_fytest = sandy.DecayData.from_endf6(decay_minimal)
        >>> mass_numbers = [147, 148, 149]
        >>> S = decay_fytest.get_chain_yield_sensitivity()
        >>> chain_yields = sandy.fy.get_chain_yields()
        >>> y_extra = chain_yields.query(f"A=={mass_numbers} & E=='thermal' & ZAM==922350").set_index("A").CHY
        >>> std_extra_info = chain_yields.query(f"A=={mass_numbers} & E=='thermal'").set_index("A").DCHY
        >>> Vy_extra = sandy.CategoryCov.from_stdev(std_extra_info).data
        >>> tape = sandy.get_endf6_file("jeff_33", "nfpy", 922350)
        >>> nfpy = sandy.Fy.from_endf6(tape)
        >>> nfpy = sandy.Fy(nfpy.data.query(f"ZAP=={zam}"))
        >>> nfpy.gls_update(922350, 0.0253, S, y_extra, Vy_extra)[0]
             MAT   MT     ZAM     ZAP           E          FY         DFY
         0  9228  454  922350  591480 2.53000e-02 8.40795e-04 2.39023e-05
         1  9228  454  922350  591481 2.53000e-02 1.39006e-02 4.22282e-05
         2  9228  454  922350  601480 2.53000e-02 5.12689e-05 5.33871e-06
         3  9228  459  922350  591480 2.53000e-02 1.66100e-02 1.21560e-04
         4  9228  459  922350  591481 2.53000e-02 3.02430e-04 1.01150e-04
         5  9228  459  922350  601480 2.53000e-02 1.69270e-02 1.18300e-04

        >>> nfpy.gls_update(922350, 0.0253, S, y_extra, Vy_extra)[1]
        ZAP          591480       591481       601480
        ZAP
        591480  5.71318e-10 -4.98686e-10 -1.34355e-12
        591481 -4.98686e-10  1.78322e-09 -2.37615e-11
        601480 -1.34355e-12 -2.37615e-11  2.85018e-11

        >>> y_constraint = pd.Series([2], index=[0])
        >>> tape = sandy.get_endf6_file("jeff_33", "nfpy", 922340)
        >>> nfpy = sandy.Fy.from_endf6(tape)
        >>> S = pd.DataFrame(np.ones(len(nfpy.data.query("E==400000 & MT==454").ZAP)), index=nfpy.data.query("E==400000 & MT==454").ZAP.to_list()).T
        >>> nfpy_post = nfpy.gls_update(922340, 400000, S, y_constraint)[0]
        >>> np.testing.assert_almost_equal(sum(nfpy_post.data.query("MT==454").FY), 2)
        """
        fy_data = self.data.query(f"ZAM=={zam} & E=={energy} & MT==454").set_index('ZAP')
        mat = fy_data.MAT.iloc[0]
        x_prior = fy_data.FY
        std = fy_data.DFY
        Vx_prior = CategoryCov.from_stdev(std)
        index = Vx_prior.data.index
        S_ = S.reindex(columns=index).fillna(0)
        y_extra_ = y_extra.reindex(index=S_.index).fillna(0)
        # Perform GLS:
        if Vy_extra is not None:
            S_ = S_.reindex(index=Vy_extra.index).fillna(0)
            y_extra_ = y_extra.reindex(index=S_.index).fillna(0).values
            x_post = _gls_parameters_update(x_prior.values, S_.values,
                                                  Vx_prior.data.values, y_extra_, Vy_extra.values)
        else:
            x_post = _gls_parameters_update(x_prior.values, S_.values,
                                                  Vx_prior.data.values, y_extra_.values)
        Vx_post = Vx_prior.gls_cov_update(S_, Vy_extra)
        # Results in appropriate format:
        data = self.data.query(f"ZAM=={zam} & E=={energy} & MT==459")
        dx_post = np.sqrt(np.diag(Vx_post.data))
        new_info = pd.DataFrame({"ZAP": index, "FY": x_post, "DFY": dx_post})
        new_info = new_info.assign(MAT=mat, ZAM=zam, MT=454, E=energy)
        return self.__class__(pd.concat([new_info, data], ignore_index=True)), Vx_post

    def ishikawa_factor(self, zam, e, Vy_extra,
                        kind='mass yield', decay_data=None):
        r"""
        Ishikawa factor to determine whether the experiment from where we
        obtain model sensitivity is useful to reduce the IFY uncertainty

        Parameters
        ----------
        zam : `int`
            ZAM number of the material to which calculations are to be
            applied.
        e : `float`
            Energy to which calculations are to be applied.
        Vy_extra : 2D iterable
            2D covariance matrix for y_extra (MXM).
        kind : `str`, optional
            Keyword for obtaining sensitivity. The
            default is 'mass yield'.
        decay_data : `DecayData`, optional
            Object to change the model to CFY = Q*IFY or
            Ch_chain = S_chain*IFY, so the sensitivity (S) is Q or S_chain.
            The default is None, so the model is ChY = ChY_mass*IFY
            and the sensitivity is mass yield sensitivity.

        Returns
        -------
        ishikawa : `pd.Series`
            Ishikawa factor.

        Results:
        -------
        Ishikawa factor << 1 :
            The extra data is not so useful and the data remain unchanged.
        Ishikawa factor >> 1 :
            The extra data very useful and the 'posteriori' covariance will
            be reduced to the same level as the integral parameter covariance.
        Ishikawa factor ~ 1 :
            The experiment is useful and the 'posteriori'  covariance will be
            reduced by approximately half

        Examples
        --------

        >>> import sandy
        >>> zam = [591480, 591481, 601480]
        >>> decay_minimal = sandy.get_endf6_file("jeff_33", 'decay', zam)
        >>> decay_fytest = sandy.DecayData.from_endf6(decay_minimal)
        >>> CFY_var_extra = np.diag(pd.Series([1, 1, 1]))
        >>> CFY_var_extra = pd.DataFrame(CFY_var_extra, index=zam, columns=zam)
        >>> npfy = sandy.Fy(minimal_fytest_2)
        >>> npfy.ishikawa_factor(942390, 500e3, CFY_var_extra, kind='cumulative', decay_data=decay_fytest)
        591480   4.00000e-02
        591481   5.00000e-02
        601480   1.00000e-01
        dtype: float64

        >>> A = [147, 148, 149]
        >>> Chain_var_extra = np.diag(pd.Series([1, 1, 1]))
        >>> Chain_var_extra = pd.DataFrame(Chain_var_extra, index=A, columns=A)
        >>> npfy = sandy.Fy(minimal_fytest_2)
        >>> npfy.ishikawa_factor(942390, 500e3, Chain_var_extra)
        147   0.00000e+00
        148   1.00000e-01
        149   0.00000e+00
        dtype: float64
        """
        # Filter FY data:
        conditions = {'ZAM': zam, "E": e}
        fy_data = self._filters(conditions).data.set_index('ZAP')
        Vx_prior = fy_data.query('MT==454').DFY
        Vx_prior = CategoryCov.from_var(Vx_prior).data
        # Find the GLS sensitivity:
        if kind == 'mass yield':
            model_sensitivity_object = self._filters(conditions)
        elif kind == 'cumulative' or 'chain yield':
            model_sensitivity_object = decay_data
        S = _gls_setup(model_sensitivity_object, kind)
        # Perform Ishikawa factor:
        Vy_extra_ = pd.DataFrame(Vy_extra)
        index = Vy_extra_.index
        S = S.reindex(index=index, columns=Vx_prior.index).fillna(0)
        ishikawa = ishikawa_factor(S.values, Vx_prior.values, Vy_extra_.values)
        ishikawa = pd.Series(ishikawa, index=index)
        return ishikawa

    def _filters(self, conditions):
        r"""
        Apply several condition to source data and return filtered results.

        Parameters
        ----------
        conditions : `dict`
            Conditions to apply, where the key is any label present in the
            columns of `data` and the value is filtering condition.

        Returns
        -------
        :obj:`~sandy.fy.Fy`
            filtered dataframe of fission yields

        Examples
        --------

        >>> conditions = {"ZAP":380900, "E":2.53000e-07}
        >>> Fy(minimal_fytest)._filters(conditions)
            MAT   MT     ZAM     ZAP           E          FY         DFY
        0  9437  454  942390  380900 2.53000e-07 2.00000e-01 2.00000e-02
        """
        conditions = dict(conditions)
        out = self
        for keys, values in conditions.items():
            out = out.filter_by(keys, values)
        return out

    def filter_by(self, key, value):
        r"""
        Apply condition to source data and return filtered results.

        Parameters
        ----------
        `key` : `str`
            any label present in the columns of `data`
        `value` : `int` or `float`
            value used as filtering condition

        Returns
        -------
        :obj:`~sandy.fy.Fy`
            filtered dataframe of fission yields

        Notes
        -----
        .. note:: The primary function of this method is to make sure that
                  the filtered dataframe is still returned as a `Fy` object.

        Examples
        --------

        >>> Fy(minimal_fytest).filter_by("ZAP", 380900)
            MAT   MT     ZAM     ZAP           E          FY         DFY
        0  9437  454  942390  380900 2.53000e-07 2.00000e-01 2.00000e-02
        1  9437  454  942390  380900 5.00000e+05 8.00000e-01 4.00000e-02

        >>> Fy(minimal_fytest).filter_by("E", 5e5)
            MAT   MT     ZAM     ZAP           E          FY         DFY
        0  9437  454  942390  380900 5.00000e+05 8.00000e-01 4.00000e-02
        1  9437  454  942390  551370 5.00000e+05 1.00000e+00 5.00000e-02
        2  9437  454  942390  541350 5.00000e+05 2.00000e-01 1.00000e-02
        """
        condition = self.data[key] == value
        out = self.data.copy()[condition].reset_index(drop=True)
        return self.__class__(out)

    @classmethod
    def from_endf6(cls, endf6, verbose=False):
        r"""
        Extract fission yields from `Endf6` instance.

        Parameters
        ----------
        endf6 : :obj:`~sandy.core.endf6.Endf6`
            object containing the ENDF-6 text
        verbose : `bool`, optional, default is `False`
            flag to print information when reading ENDF-6 file

        Returns
        -------
        :obj:`~sandy.fy.Fy`
            fission yield object

        Notes
        -----
        .. note:: Both independent and cumulative fission product yields are
                  loaded, if found.

        Examples
        --------

        >>> import sandy
        >>> tape = sandy.get_endf6_file("jeff_33", "nfpy", 'all')
        >>> fy = sandy.Fy.from_endf6(tape)
        >>> fy.data.query("ZAM==952421 & MT==454 & E==0.0253").head()
                MAT   MT     ZAM    ZAP           E          FY         DFY
        56250  9547  454  952421  10010 2.53000e-02 3.32190e-05 1.17790e-05
        56251  9547  454  952421  10020 2.53000e-02 1.01520e-05 3.52080e-06
        56252  9547  454  952421  10030 2.53000e-02 1.60000e-04 5.00220e-05
        56253  9547  454  952421  20030 2.53000e-02 0.00000e+00 0.00000e+00
        56254  9547  454  952421  20040 2.53000e-02 2.10000e-03 6.64080e-04
        """
        data = []
        dict_zam = {}
        for (mat, mf, mt) in endf6.keys:
            sec = endf6.read_section(mat, mf, mt)
            if mf == 1:
                dict_zam[mat] = int(sec["ZA"] * 10 + sec["LISO"])
            else:
                if verbose:
                    logging.info(f"reading 'MAT={mat}/MT={mt}'...")
                for e in sec["E"]:
                    for zap in sec["E"][e]["ZAP"]:
                        fy = sec["E"][e]["ZAP"][zap]["FY"]
                        dfy = sec["E"][e]["ZAP"][zap]["DFY"]
                        values = (mat, mt, zap, e, fy, dfy)
                        data.append(dict(zip(["MAT", "MT", "ZAP", "E", "FY", "DFY"], values)))
        df_zam = pd.DataFrame([dict_zam]).T.reset_index()
        df_zam.columns = ['MAT', 'ZAM']
        df = pd.DataFrame(data).merge(df_zam, on='MAT')
        return cls(df)

    def to_endf6(self, endf6):
        r"""
        Update fission yields in `Endf6` instance with those available in a
        `Fy` instance.

        .. warning:: only IFY and CFY that are originally
                     present in the `Endf6` instance are modified

        Parameters
        ----------
        `endf6` : :obj:`~sandy.core.endf6.Endf6`
            ENDF6 object.

        Returns
        -------
        :obj:`~sandy.core.endf6.Endf6`
            ENDF6 objects with updated IFY and CFY.

        Examples
        --------

        >>> import sandy
        >>> tape = sandy.get_endf6_file("jeff_33", "nfpy", "all")
        >>> fy = sandy.Fy.from_endf6(tape)
        >>> new_tape = fy.to_endf6(tape)
        >>> new_tape.filter_by(listmat= [9640], listmf=[8], listmt=[454, 459])
        MAT   MF  MT
        9640  8   454     96245.0000 242.960000          2          0  ...
                  459     96245.0000 242.960000          2          0  ...
        dtype: object
        """
        data_endf6 = Endf6(endf6.data.copy())
        mf = 8
        for (mat, mt, e), data_fy in self.data.groupby(['MAT', 'MT', 'E']):
            sec = data_endf6.read_section(mat, mf, mt)
            new_data = data_fy.set_index('ZAP')[['FY', 'DFY']].T.to_dict()
            sec['E'][e]['ZAP'] = new_data
            data_endf6.data[(mat, mf, mt)] = write_mf8(sec)
        return data_endf6



def _gls_setup(model_sensitivity_object, kind):
    """
    A function to obtain the sensitivity for performing GLS.

    Parameters
    ----------
    model_sensitivity_object : `DecayData` or `Fy`
        Object from which the sensitivity (S) of the model is to be derived:
            y_calc = S*x_prior
    kind : `str`
        Keyword for obtaining sensitivity.

    Raises
    ------
    TypeError
        The kind is not implemented in the method.

    Returns
    -------
    S : `pandas.DataFrame`
        The sensitivity for performing GLS.

    """
    if kind == 'cumulative':
        S = model_sensitivity_object.get_qmatrix()
    elif kind == 'chain yield':
        S = model_sensitivity_object.get_chain_yield_sensitivity()
    elif kind == 'mass yield':
        S = model_sensitivity_object.get_mass_yield_sensitivity()
    else:
        raise ValueError('Keyword argument "kind" is not valid')
    return S

