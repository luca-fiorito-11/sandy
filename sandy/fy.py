# -*- coding: utf-8 -*-
"""
This module contains all classes and functions specific for processing fission
yield data.
"""
from tables import NaturalNameWarning
import h5py
import logging
import warnings

import pandas as pd
import numpy as np
import scipy.sparse as sps
import os

import sandy
from sandy.shared import expand_zam

__author__ = "Luca Fiorito"
__all__ = [
        "Fy",
        "fy2hdf",
        ]


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


def import_chain_yields():
    """
    Import chain yields information from information store in sandy. The
    information was taken from 'https://www-nds.iaea.org/endf349/la-ur-94-3106.pdf',
    page 18-29.

    Returns
    -------
    data : `pd.Dataframe`
        Information of the url divided into a dataframe.

    Notes
    -----
    ..note :: Values in ENDF/B-VI have been extended to cover all nuclides in
    the ENDF/B-VI decay files and also four charge units from Zp, but the
    extension does not alter the values in this document for the specific
    nuclides listed in the appendices.
    ..note :: It is recognized that even smaller independent yields have even
    larger than 100% error.Indeed, independent yields in the range of 1.0e-9 to
    1.0e-12 may be difficult to predict to better than a factor of 100. For
    this reason, all yields less than 1.0e-12 are blanked out.

    Examples
    --------
    >>> sandy.fy.import_chain_yields().head()
           ZAM	         Ch	          E	        DCh
     A				
    66	902270	8.53000e-08	2.53000e-02	2.72960e-08
    67	902270	2.16000e-07	2.53000e-02	6.91200e-08
    68	902270	8.23000e-07	2.53000e-02	2.63360e-07
    69	902270	2.35000e-06	2.53000e-02	7.52000e-07
    70	902270	5.19000e-06	2.53000e-02	1.66080e-06
    """
    letters = {'a': 0.0035, 'b': 0.0050, 'c': 0.0070, 'd': 0.01, 'e': 0.014,
               'f': 0.02, 'g': 0.028, 'h': 0.04,
               'i': 0.06, 'j': 0.08, 'k': 0.11, 'l': 0.16, 'm': 0.23,
               'n': 0.32, 'o': 0.45, 'p': 0.64}
    energy = {'t': 2.53e-2, 'f': 2000000.0, 'h': 14000000.0,
              's': 'spontaneous fission'}
    nuclides_z = {'th': 90, 'pa': 91, 'u': 92, 'np': 93, 'pu': 94, 'am': 95,
                  'cm': 96, 'cf': 98, 'es': 99, 'fm': 100}
    names = {'A': ['mass', 'u235t', 'u235f', 'u235h', 'u238f', 'u238h',
                   'pu239t', 'pu239f', 'pu241t', 'u233t', 'th232f'],
             'B': ['mass', 'u233f', 'u233h', 'u236f', 'pu239h', 'pu240f',
                   'pu241f', 'pu242f', 'th232h', 'np237f', 'cf252s'],
             'C': ['mass', 'u234f', 'u237f', 'pu240h', 'u234h', 'u236h',
                   'pu238f', 'am241f', 'am243f', 'np238f', 'cm242f'],
             'D': ['mass', 'th227t', 'th229t', 'pa231f', 'am241t', 'am241h',
                   'am242t', 'cm245t', 'cf249t', 'cf251t', 'es254t'],
             'E': ['mass', 'cf250s', 'cm244s', 'cm248s', 'es253s', 'fm254s',
                   'fm255t', 'fm256s', 'np237h', 'u232t', 'u238s'],
             'F': ['mass', 'cm243t', 'cm246s', 'cm243f', 'cm244f', 'cm246f',
                   'cm248f', 'pu242h', 'np237t', 'pu240t', 'pu242t']}
    files = ['appendix A.txt', 'appendix B.txt', 'appendix C.txt',
             'appendix D.txt', 'appendix E.txt', 'appendix F.txt']
    colspecs = [(0, 3), (3, 13), (13, 23), (23, 33), (33, 43), (43, 53),
                (53, 63), (63, 73), (73, 83), (83, 93), (93, 103)]
    data = pd.DataFrame(columns=['A', 'ZAM', 'Ch', 'E', 'DCh'])
    for file in files:
        df = pd.read_fwf(os.path.join(os.path.dirname(__file__),
                                      'appendix',
                                      'chain yields',
                                      file),
                         colspecs = colspecs,
                         names = names[file[file.find('.')-1]])[1::].rename(columns ={'mass':'A'}).set_index('A')
        df = df.stack().to_frame().reset_index().rename(columns = {'level_1': 'ZAM', 0: 'Ch'})
        df['E'] = [value for index, valor in df['ZAM'].items()
                   for key, value in energy.items() if key in valor[-2:]]
        df['ZAM'] = [valor[::-1].replace(key, '', 1)[::-1]
                     for index, valor in df['ZAM'].items()
                     for key, value in energy.items() if key in valor[-2:]]
        df['ZAM'] = [nuclides_z[valor[:-3]]*10000 + 10*int(valor[-3:])
                     for index, valor in df['ZAM'].items()]
        df['DCh'] = df['Ch'].str[-1].replace(letters, regex=True)
        df['Ch'] = df['Ch'].str[:-1].astype(float)
        df['DCh'] = df['Ch']*df['DCh']
        df['A'] = df['A'].astype(int)
        df = df[df.E != 'spontaneous fission']
        data = data.append(df.reset_index(drop=True))
    data = data.sort_values(['ZAM', 'E'], ascending=[True, False])\
               .reset_index(drop=True).set_index('A')
    return data


class Fy():
    """
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
        """
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
        """
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
        """
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
        expand_zam = sandy.shared.expand_zam
        zam = pd.DataFrame(map(expand_zam, self.data.ZAP),
                           columns=["Z", "A", "M"],
                           dtype=int)
        return self.data.assign(Z=zam.Z, A=zam.A, M=zam.M)

    def _expand_zam(self):
        """
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

    def get_mass_yield(self, zam, e, sparse=False):
        """
        Obtain mass yields from the following model: ChY = S * IFY

        Parameters
        ----------
        zam : `int`
            ZAM number of the fissioning nuclide.
        e : `float`
            Energy of the fissioning system.
        sparse : `bool`, optional
            Option to use sparse matrix for calculations. The default is False.

        Returns
        -------
        `pandas.Series`
            mass yield obtained from ChY = S * IFY

        Examples
        --------
        >>> tape_nfpy = sandy.get_endf6_file("jeff_33",'nfpy','all')
        >>> nfpy = Fy.from_endf6(tape_nfpy)
        >>> nfpy.get_mass_yield(922350, 0.0253).loc[148]
        0.0169029147
        """
        # Filter FY data:
        conditions = {'ZAM': zam, "E": e, 'MT': 454}
        fy_data = self._filters(conditions).data.set_index('ZAP').FY
        chain_data = self._filters({'ZAM': zam, "E": e})
        kind = 'mass yield'
        S = _gls_setup(chain_data, kind)
        chain = sandy._y_calc(fy_data, S, sparse=sparse)
        return chain.rename('mass yield')

    def get_chain_yield(self, zam, e, decay_data, sparse=False, **kwargs):
        """
        Obtain chain yields from the following model: ChY = S * IFY

        Parameters
        ----------
        zam : `int`
            ZAM number of the fissioning nuclide.
        e : `float`
            Energy of the fissioning system.
        decay_data : `sandy.DecayData`
            Radioactive nuclide data from where to obtain chain sensitivities.
        sparse : `bool`, optional
            Option to use sparse matrix for calculations. The default is False.
        kwargs : `dict`
            keyword arguments for method `get_decay_chains`

        Returns
        -------
        `pandas.Series`
            Chain yield obtained from ChY = S * IFY

        Examples
        --------
        >>> zam = [591480, 591481, 601480]
        >>> decay_minimal = sandy.get_endf6_file("jeff_33", 'decay', zam)
        >>> decay_fytest = sandy.DecayData.from_endf6(decay_minimal)
        >>> npfy = sandy.Fy(minimal_fytest_2)
        >>> chain = npfy.get_chain_yield(942390, 500e3, decay_fytest).iloc[0].round(2) 
        >>> assert chain == 0.6
        >>> chain = npfy.get_chain_yield(942390, 500e3, decay_fytest, sparse=True).iloc[0].round(2) 
        >>> assert chain == 0.6
        """
        # Filter FY data:
        conditions = {'ZAM': zam, "E": e, 'MT': 454}
        fy_data = self._filters(conditions).data.set_index('ZAP').FY
        kind = 'chain yield'
        S = _gls_setup(decay_data, kind)
        chain = sandy._y_calc(fy_data, S, sparse=sparse)
        return chain.rename('chain yield')

    def get_mass_yield_sensitivity(self):
        """
        Obtain the mass yield sensitivity matrix based only on the
        information given in the `Fy` object (no decay data).

        Returns
        -------
        mass_yield_sensitivity : `pandas.DataFrame`
            Mass yield sensitivity matrix. Mass number in index and ZAM of
            fission products in columns.

        Examples
        --------
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
        groups = fy_data.groupby(fy_data.index)['ZAP'].value_counts()
        groups = groups.to_frame()\
                       .rename(columns={'ZAP': 'value'})\
                       .reset_index()
        mass_yield_sensitivity = groups.pivot_table(
            index='A',
            columns='ZAP',
            values='value',
            aggfunc="sum"
            ).fillna(0)
        return mass_yield_sensitivity

    def custom_perturbation(self, zam, mt, e, zap, pert):
        """
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
        zap : `int`
            ZAP of the product to which perturbations are to be
            applied.
        pert : `float`
            Perturbation coefficients as ratio values.

        Returns
        -------
        `Fy`
            Fission yield instance with applied perturbation.

        Examples
        --------
        >>> tape = sandy.get_endf6_file("jeff_33", 'nfpy', 'all')
        >>> nfpy = Fy.from_endf6(tape)
        >>> nfpy_pert = nfpy.custom_perturbation(922350, 459, 0.0253, 551370, 0.9)
        >>> comp = nfpy_pert.data.query('ZAM==922350 & ZAP==551370 & MT==459 & E==0.0253').squeeze().FY
        >>> assert np.setdiff1d(nfpy_pert.data.values, nfpy.data.values) == comp
        """
        df = self.data.copy()
        mask = (df.ZAM == zam) & (df.ZAP == zap) & (df.MT == mt) & (df.E == e)
        df.loc[mask, "FY"] = df[mask].squeeze().FY * pert
        return self.__class__(df)

    def apply_bmatrix(self, zam, e, decay_data, sparse=False, rows=None,
                      keep_fy_index=False, threshold=None, fill_negative=None):
        """
        Perform IFY = (1-B)*CFY equation to calculate IFY in a given zam
        for a given energy and apply into the original data.

        Parameters
        ----------
        zam : `int`
            ZAM number of the material to which calculations are to be
            applied.
        e : `float`
            Energy to which calculations are to be applied.
        decay_data : `sandy.DecayData`
            Radioactive nuclide data for several isotopes.
        sparse : `bool`, optional
            Option to use sparse matrix for calculations. The default is False.
        rows : `int`, optional
            Option to use row calculation for matrix calculations. This option
            defines the number of lines to be taken into account in each loop.
            The default is None.
        keep_fy_index : `bool`, optional
            Option that allows you to output only the CFY results that were
            part of the original `sandy.Fy` object. The default is False.
        threshold : `int`, optional
            Optional argument to avoid numerical fluctuations or
            values so small that they do not have to be taken into
            account. The default is None.
         fill_negative: `float`, optional
            Optional argument to avoid negative IFY after applying this method.
            Since the sensitivity of the method is (1-B), if the nucleus in
            question has a small CFY and is a common product in several decays,
            there is a possibility that the calculated IFY will be negative.
            This option allows us to change all negative values to a value
            chosen by the user. The default is None.

        Returns
        -------
        `sandy.Fy`
            Fission yield instance with IFY calculated for a given combination
            of ZAM/e/decay_data.

        Notes
        -----
        .. note:: This method applies a perturbation to certain IFYs,
        since the equation IFY = (1-B)*CFY is not satisfied for all nuclei.
        .. note:: Some cores may have a negative calculated IFY, which is
        impossible. Therefore, there is a `fill_negative` option that allows us
        to change these negative values to a value chosen by the user.

        Examples
        --------
        >>> zam = [591480, 591481, 601480]
        >>> decay_minimal = sandy.get_endf6_file("jeff_33", 'decay', zam)
        >>> decay_fytest = sandy.DecayData.from_endf6(decay_minimal)
        >>> npfy = Fy(minimal_fytest_2)
        >>> npfy_pert = npfy.apply_bmatrix(942390, 5.00000e+05, decay_fytest)
        >>> npfy_pert.data[npfy_pert.data.MT == 454]
        	 MAT	 MT	   ZAM	   ZAP	          E	          FY	        DFY
        3	9437	454	942390	591480	5.00000e+05	 8.00000e-01	4.00000e-02
        4	9437	454	942390	591481	5.00000e+05	 1.00000e+00	5.00000e-02
        5	9437	454	942390	601480	5.00000e+05	-1.60000e+00	1.00000e-01
        6	9437	454	942390	621480	5.00000e+05	-2.00000e-01	1.00000e-02 

        >>> npfy_pert = npfy.apply_bmatrix(942390, 5.00000e+05, decay_fytest, sparse=True)
        >>> npfy_pert.data[npfy_pert.data.MT == 454]
        	 MAT	 MT	   ZAM	   ZAP	          E	          FY	        DFY
        3	9437	454	942390	591480	5.00000e+05	 8.00000e-01	4.00000e-02
        4	9437	454	942390	591481	5.00000e+05	 1.00000e+00	5.00000e-02
        5	9437	454	942390	601480	5.00000e+05	-1.60000e+00	1.00000e-01
        6	9437	454	942390	621480	5.00000e+05	-2.00000e-01	1.00000e-02

        >>> npfy_pert = npfy.apply_bmatrix(942390, 5.00000e+05, decay_fytest, rows=1)
        >>> npfy_pert.data[npfy_pert.data.MT == 454]
        	 MAT	 MT	   ZAM	   ZAP	          E	          FY	        DFY
        3	9437	454	942390	591480	5.00000e+05	 8.00000e-01	4.00000e-02
        4	9437	454	942390	591481	5.00000e+05	 1.00000e+00	5.00000e-02
        5	9437	454	942390	601480	5.00000e+05	-1.60000e+00	1.00000e-01
        6	9437	454	942390	621480	5.00000e+05	-2.00000e-01	1.00000e-02

        >>> npfy_pert = npfy.apply_bmatrix(942390, 5.00000e+05, decay_fytest, fill_negative=0)
        >>> npfy_pert.data[npfy_pert.data.MT == 454]
        	 MAT	 MT	   ZAM	   ZAP	          E	          FY	        DFY
        3	9437	454	942390	591480	5.00000e+05	 8.00000e-01	4.00000e-02
        4	9437	454	942390	591481	5.00000e+05	 1.00000e+00	5.00000e-02
        5	9437	454	942390	601480	5.00000e+05	 0.00000e+00	0.00000e+00
        6	9437	454	942390	621480	5.00000e+05	 0.00000e+00	0.00000e+00

        >>> npfy_pert = npfy.apply_bmatrix(942390, 5.00000e+05, decay_fytest, sparse=True, fill_negative=0)
        >>> npfy_pert.data[npfy_pert.data.MT == 454]
        	 MAT	 MT	   ZAM	   ZAP	          E	          FY	        DFY
        3	9437	454	942390	591480	5.00000e+05	 8.00000e-01	4.00000e-02
        4	9437	454	942390	591481	5.00000e+05	 1.00000e+00	5.00000e-02
        5	9437	454	942390	601480	5.00000e+05	 0.00000e+00	0.00000e+00
        6	9437	454	942390	621480	5.00000e+05	 0.00000e+00	0.00000e+00

        >>> npfy_pert = npfy.apply_bmatrix(942390, 5.00000e+05, decay_fytest, rows=1, fill_negative=0)
        >>> npfy_pert.data[npfy_pert.data.MT == 454]
        	 MAT	 MT	   ZAM	   ZAP	          E	          FY	        DFY
        3	9437	454	942390	591480	5.00000e+05	 8.00000e-01	4.00000e-02
        4	9437	454	942390	591481	5.00000e+05	 1.00000e+00	5.00000e-02
        5	9437	454	942390	601480	5.00000e+05	 0.00000e+00	0.00000e+00
        6	9437	454	942390	621480	5.00000e+05	 0.00000e+00	0.00000e+00
        
        >>> zam = [591480, 591481, 601480]
        >>> decay_minimal = sandy.get_endf6_file("jeff_33", 'decay', zam)
        >>> decay_fytest = sandy.DecayData.from_endf6(decay_minimal)
        >>> npfy = Fy(minimal_fytest_2)
        >>> npfy_pert = npfy.apply_bmatrix(942390, 5.00000e+05, decay_fytest, keep_fy_index=True)
        >>> npfy_pert.data[npfy_pert.data.MT == 454]
             MAT	 MT	   ZAM	   ZAP	          E	         FY	            DFY
        3	9437	454	942390	591480	5.00000e+05	 8.00000e-01	4.00000e-02
        4	9437	454	942390	591481	5.00000e+05	 1.00000e+00	5.00000e-02
        5	9437	454	942390	601480	5.00000e+05	-1.60000e+00	1.00000e-01

        >>> npfy_pert = npfy.apply_bmatrix(942390, 5.00000e+05, decay_fytest, sparse=True, keep_fy_index=True)
        >>> npfy_pert.data[npfy_pert.data.MT == 454]
             MAT	 MT	   ZAM	   ZAP	          E	         FY	            DFY
        3	9437	454	942390	591480	5.00000e+05	 8.00000e-01	4.00000e-02
        4	9437	454	942390	591481	5.00000e+05	 1.00000e+00	5.00000e-02
        5	9437	454	942390	601480	5.00000e+05	-1.60000e+00	1.00000e-01

        >>> npfy_pert = npfy.apply_bmatrix(942390, 5.00000e+05, decay_fytest, rows=1, keep_fy_index=True)
        >>> npfy_pert.data[npfy_pert.data.MT == 454]
             MAT	 MT	   ZAM	   ZAP	          E	         FY	            DFY
        3	9437	454	942390	591480	5.00000e+05	 8.00000e-01	4.00000e-02
        4	9437	454	942390	591481	5.00000e+05	 1.00000e+00	5.00000e-02
        5	9437	454	942390	601480	5.00000e+05	-1.60000e+00	1.00000e-01

        >>> npfy_pert = npfy.apply_bmatrix(942390, 5.00000e+05, decay_fytest, keep_fy_index=True, fill_negative=0)
        >>> npfy_pert.data[npfy_pert.data.MT == 454]
             MAT	 MT	   ZAM	   ZAP	          E	         FY	            DFY
        3	9437	454	942390	591480	5.00000e+05	 8.00000e-01	4.00000e-02
        4	9437	454	942390	591481	5.00000e+05	 1.00000e+00	5.00000e-02
        5	9437	454	942390	601480	5.00000e+05	 0.00000e+00	0.00000e+00

        >>> npfy_pert = npfy.apply_bmatrix(942390, 5.00000e+05, decay_fytest, sparse=True, keep_fy_index=True, fill_negative=0)
        >>> npfy_pert.data[npfy_pert.data.MT == 454]
             MAT	 MT	   ZAM	   ZAP	          E	         FY	            DFY
        3	9437	454	942390	591480	5.00000e+05	 8.00000e-01	4.00000e-02
        4	9437	454	942390	591481	5.00000e+05	 1.00000e+00	5.00000e-02
        5	9437	454	942390	601480	5.00000e+05	 0.00000e+00	0.00000e+00

        >>> npfy_pert = npfy.apply_bmatrix(942390, 5.00000e+05, decay_fytest, rows=1, keep_fy_index=True, fill_negative=0)
        >>> npfy_pert.data[npfy_pert.data.MT == 454]
             MAT	 MT	   ZAM	   ZAP	          E	         FY	            DFY
        3	9437	454	942390	591480	5.00000e+05	 8.00000e-01	4.00000e-02
        4	9437	454	942390	591481	5.00000e+05	 1.00000e+00	5.00000e-02
        5	9437	454	942390	601480	5.00000e+05	 0.00000e+00	0.00000e+00
        """
        # Obtain the data:
        data = self.data.copy()
        conditions = {'ZAM': zam, 'MT': 459, "E": e}
        fy_data = self._filters(conditions).data
        mat = fy_data.MAT.iloc[0]
        cov_data = fy_data.set_index('ZAP')['DFY']
        fy_data = fy_data.set_index('ZAP')['FY']
        if keep_fy_index:
            original_index = fy_data.index
        B = decay_data.get_bmatrix()
        if 10 in B.index:
            B.drop(index=10, inplace=True)
        if 10 in B.columns:
            B.drop(columns=10, inplace=True)
        # Put the data in a appropriate format:
        index, columns = B.index, B.columns
        # Creating (1-B) matrix:
        if sparse or rows is not None:
            B = sps.csc_matrix(B)
            unit = sps.csc_matrix(sps.identity(B.shape[0]))
            C = unit - B
            sensitivity = pd.DataFrame(C.toarray(), index=index, columns=columns)
        else:
            unit = np.identity(len(B))
            C = unit - B.values
            sensitivity = pd.DataFrame(C, index=index, columns=columns)
        # Rest of the data
        mask = (data.ZAM == zam) & (data.MT == 454) & (data.E == e)
        data = data.loc[~mask]
        columns_union = columns.union(fy_data.index)
        sensitivity = sensitivity.reindex(columns=columns_union).fillna(0)
        if len(columns) != len(columns_union):
            sensitivity = sensitivity.reindex(index=columns_union).fillna(0)
        np.fill_diagonal(sensitivity.values, 1)
        fy_data = fy_data.reindex(columns_union).fillna(0)
        cov_data = cov_data.reindex(columns_union).fillna(0)
        cov_data = sandy.CategoryCov.from_var(cov_data)
        # Apply (1-B) matrix
        if rows is not None:
            ify_calc_values = sandy._y_calc(fy_data, sensitivity, sparse=True).rename('FY')
        else:
            ify_calc_values = sandy._y_calc(fy_data, sensitivity, sparse=sparse).rename('FY')
        cov_calc_values = np.diag(cov_data._gls_Vy_calc(sensitivity,
                                                        sparse=sparse,
                                                        rows=rows))
        cov_calc_values = pd.Series(cov_calc_values, index=sensitivity.index)
        if keep_fy_index:
            ify_calc_values = ify_calc_values.reindex(original_index).fillna(0)
            cov_calc_values = cov_calc_values.reindex(original_index).fillna(0)
        if threshold is not None:
            ify_calc_values[abs(ify_calc_values) < threshold] = 0
            cov_calc_values[abs(cov_calc_values) < threshold] = 0
        calc_values = ify_calc_values.reset_index().rename(columns={'DAUGHTER': 'ZAP'})
        calc_values['DFY'] = cov_calc_values.values
        # Change negative values to zero
        if fill_negative is not None:
            calc_values.loc[calc_values.FY < 0, ['FY', 'DFY']] = [0, 0]
        # Calculus in appropiate way:
        calc_values[['MAT', 'ZAM', 'MT', 'E']] = [mat, zam, 454, e]
        data = pd.concat([data, calc_values], ignore_index=True)
        return self.__class__(data)

    def apply_qmatrix(self, zam, energy, decay_data, sparse=False, rows=None,
                      keep_fy_index=False, threshold=None):
        """
        Perform CFY = Q*IFY equation to calculate CFY in a given zam
        for a given energy and apply into the original data.

        Parameters
        ----------
        zam : `int`
            ZAM number of the material to which calculations are to be
            applied.
        e : `float`
            Energy to which calculations are to be applied.
        decay_data : `sandy.DecayData`
            Radioactive nuclide data for several isotopes.
        sparse : `bool`, optional
            Option to use sparse matrix for calculations. The default is False.
        rows : `int`, optional
            Option to use row calculation for matrix calculations. This option
            defines the number of lines to be taken into account in each loop.
            The default is None.
        keep_fy_index: `bool`, optional
            Option that allows you to output only the CFY results that were
            part of the original `sandy.Fy` object. The default is False.
        threshold : `int`, optional
            Optional argument to avoid numerical fluctuations or
            values so small that they do not have to be taken into
            account. The default is None.

        Returns
        -------
        `sandy.Fy`
            Fission yield instance with CFY calculated for a given combination
            of ZAM/e/decay_data.

        Notes
        -----
        .. note:: This method applies a perturbation to certain CFYs,
        since the equation CFY = Q*IFY is not satisfied for all nuclei.

        Examples
        --------
        >>> zam = [591480, 591481, 601480]
        >>> decay_minimal = sandy.get_endf6_file("jeff_33", 'decay', zam)
        >>> decay_fytest = sandy.DecayData.from_endf6(decay_minimal)
        >>> npfy = Fy(minimal_fytest_2)
        >>> npfy_pert = npfy.apply_qmatrix(942390, 5.00000e+05, decay_fytest)
        >>> npfy_pert.data[npfy_pert.data.MT == 459]
             MAT	 MT	   ZAM	   ZAP	          E	         FY 	    DFY
        3	9437	459	942390	591480	5.00000e+05	1.00000e-01	4.00000e-02
        4	9437	459	942390	591481	5.00000e+05	2.00000e-01	5.00000e-02
        5	9437	459	942390	601480	5.00000e+05	6.00000e-01	1.00000e-01
        6	9437	459	942390	621480	5.00000e+05	6.00000e-01	1.00000e-01

        >>> npfy_pert = npfy.apply_qmatrix(942390, 5.00000e+05, decay_fytest, sparse=True)
        >>> npfy_pert.data[npfy_pert.data.MT == 459]
             MAT	 MT	   ZAM	   ZAP	          E	         FY 	    DFY
        3	9437	459	942390	591480	5.00000e+05	1.00000e-01	4.00000e-02
        4	9437	459	942390	591481	5.00000e+05	2.00000e-01	5.00000e-02
        5	9437	459	942390	601480	5.00000e+05	6.00000e-01	1.00000e-01
        6	9437	459	942390	621480	5.00000e+05	6.00000e-01	1.00000e-01

        >>> npfy_pert = npfy.apply_qmatrix(942390, 5.00000e+05, decay_fytest, rows=1)
        >>> npfy_pert.data[npfy_pert.data.MT == 459]
             MAT	 MT	   ZAM	   ZAP	          E	         FY 	    DFY
        3	9437	459	942390	591480	5.00000e+05	1.00000e-01	4.00000e-02
        4	9437	459	942390	591481	5.00000e+05	2.00000e-01	5.00000e-02
        5	9437	459	942390	601480	5.00000e+05	6.00000e-01	1.00000e-01
        6	9437	459	942390	621480	5.00000e+05	6.00000e-01	1.00000e-01

        >>> zam = [591480, 591481, 601480]
        >>> decay_minimal = sandy.get_endf6_file("jeff_33", 'decay', zam)
        >>> decay_fytest = sandy.DecayData.from_endf6(decay_minimal)
        >>> npfy = Fy(minimal_fytest_2)
        >>> npfy_pert = npfy.apply_qmatrix(942390, 5.00000e+05, decay_fytest, keep_fy_index=True)
        >>> npfy_pert.data[npfy_pert.data.MT == 459]
             MAT	 MT	   ZAM	   ZAP	          E	         FY	        DFY
        3	9437	459	942390	591480	5.00000e+05	1.00000e-01	4.00000e-02
        4	9437	459	942390	591481	5.00000e+05	2.00000e-01	5.00000e-02
        5	9437	459	942390	601480	5.00000e+05	6.00000e-01	1.00000e-01

        >>> npfy_pert = npfy.apply_qmatrix(942390, 5.00000e+05, decay_fytest, sparse=True, keep_fy_index=True)
        >>> npfy_pert.data[npfy_pert.data.MT == 459]
             MAT	 MT	   ZAM	   ZAP	          E	         FY	        DFY
        3	9437	459	942390	591480	5.00000e+05	1.00000e-01	4.00000e-02
        4	9437	459	942390	591481	5.00000e+05	2.00000e-01	5.00000e-02
        5	9437	459	942390	601480	5.00000e+05	6.00000e-01	1.00000e-01

        >>> npfy_pert = npfy.apply_qmatrix(942390, 5.00000e+05, decay_fytest, rows=1, keep_fy_index=True)
        >>> npfy_pert.data[npfy_pert.data.MT == 459]
             MAT	 MT	   ZAM	   ZAP	          E	         FY	        DFY
        3	9437	459	942390	591480	5.00000e+05	1.00000e-01	4.00000e-02
        4	9437	459	942390	591481	5.00000e+05	2.00000e-01	5.00000e-02
        5	9437	459	942390	601480	5.00000e+05	6.00000e-01	1.00000e-01        
        """
        # Obtain the data:
        data = self.data.copy()
        conditions = {'ZAM': zam, 'MT': 454, "E": energy}
        fy_data = self._filters(conditions).data
        mat = fy_data.MAT.iloc[0]
        cov_data = fy_data.set_index('ZAP')['DFY']
        fy_data = fy_data.set_index('ZAP')['FY']
        if keep_fy_index:
            original_index = fy_data.index
        Q = decay_data.get_qmatrix()
        # Put the data in a approppiate format:
        union_index = Q.columns.union(fy_data.index)
        mask = (data.ZAM == zam) & (data.MT == 459) & (data.E == energy)
        data = data.loc[~mask]
        Q = Q.reindex(columns=union_index).fillna(0)
        if len(Q.loc[:, (Q == 0).all()].columns) >= 1:
            Q = Q.reindex(index=union_index).fillna(0)
            approximation = Q.loc[:, (Q == 0).all()].columns
            np.fill_diagonal(Q.loc[approximation, approximation].values, 1)
        fy_data = fy_data.reindex(union_index).fillna(0)
        cov_data = cov_data.reindex(union_index).fillna(0)
        if rows is not None:
            cov_data = sandy.CategoryCov.from_var(cov_data, sparse=True)
            # Apply qmatrix
            cfy_calc_values = sandy._y_calc(fy_data, Q, sparse=True).rename('FY')
            cov_calc_values = np.diag(cov_data._gls_Vy_calc(Q, rows=rows))
        else:
            cov_data = sandy.CategoryCov.from_var(cov_data, sparse=sparse)
            # Apply qmatrix
            cfy_calc_values = sandy._y_calc(fy_data, Q, sparse=sparse).rename('FY')
            cov_calc_values = np.diag(cov_data._gls_Vy_calc(Q, sparse=sparse))
        cov_calc_values = pd.Series(cov_calc_values, index=Q.index)
        if keep_fy_index:
            cfy_calc_values = cfy_calc_values.reindex(original_index).fillna(0)
            cov_calc_values = cov_calc_values.reindex(original_index).fillna(0)
        if threshold is not None:
            cfy_calc_values[abs(cfy_calc_values) < threshold] = 0
            cov_calc_values[abs(cov_calc_values) < threshold] = 0
        calc_values = cfy_calc_values.reset_index().rename(columns={'DAUGHTER': 'ZAP'})
        calc_values['DFY'] = cov_calc_values.values
        # Calculus in appropiate way:
        calc_values[['MAT', 'ZAM', 'MT', 'E']] = [mat, zam, 459, energy]
        data = pd.concat([data, calc_values], ignore_index=True)
        return self.__class__(data)

    def cov_update(self, zam, e, Vy_extra=None, decay_data=None,
                   kind='mass yield', sparse=False,
                   rows=None, threshold=None, **kwargs):
        """
        Update the prior IFY covariance matrix using the GLS technique. If the
        user does not enter Vy_extra, the information taken from
        'https://www-nds.iaea.org/endf349/la-ur-94-3106.pdf' (page 18-29) will
        be used in the options kind = 'mass yield' or kind = 'chain yield'.

        Notes
        -----
        ..note :: In option kind = 'cumulative', if Vy_extra is not entered,
        GLS constrained will be used instead of GLS to update the covariance
        matrix.

        Parameters
        ----------
        zam : `int`
            ZAM number of the material to which calculations are to be
            applied.
        e : `float`
            Energy to which calculations are to be applied.
        Vy_extra : 2D iterable, optional
            Extra Covariance matrix (MXM).The default is None.
        kind : `str`, optional
            Keyword for obtaining sensitivity. The default is 'mass yield'.
        decay_data : `DecayData`, optional
            Object to change the model to CFY = Q*IFY or Ch_chain = S_chain*IFY,
            so the sensitivity (S) is Q or S_chain. The default is None,
            so the model is ChY = ChY_mass*IFY and the sensitivity is mass
            yield sensitivity.
        sparse : `bool`, optional
            Option to use sparse matrix for calculations. The default is False.
        rows : `int`, optional
            Option to use row calculation for matrix calculations. This option
            defines the number of lines to be taken into account in each loop.
            The default is None.
        threshold : `int`, optional
            Optional argument to avoid numerical fluctuations or
            values so small that they do not have to be taken into
            account. The default is None.
        kwargs : `dict`
            keyword arguments for method `get_decay_chains`

        Returns
        -------
        `sandy.CategoryCov`
            IFY covariance matrix performed by GLS for a given energy and zam.

        Examples
        --------
        >>> tape_nfpy = sandy.get_endf6_file("endfb_71",'nfpy','all')
        >>> nfpy = Fy.from_endf6(tape_nfpy)
        >>> e = 2.53000e-02
        >>> zam = 922350
        >>> mass = nfpy.cov_update(zam, e, sparse=True).data.round(6)
        >>> mass.loc[401000, 401000], mass.loc[390971, 390971], mass.loc[390970, 390970], mass.loc[400980, 400980], mass.loc[411020, 411020]
        (0.00561, 0.005398, 0.005398, 0.005058, 0.003176)
        """
        if Vy_extra is None and kind != 'cumulative':
            ch_info = sandy.fy.import_chain_yields()
            ch_info = ch_info.loc[(ch_info.E == e) & (ch_info.ZAM == zam)]
            Vy_extra = sandy.CategoryCov.from_var(ch_info.DCh).data/100
        return self.gls_cov_update(zam, e, Vy_extra=Vy_extra, kind=kind,
                                   sparse=sparse, rows=rows,
                                   decay_data=decay_data, threshold=threshold,
                                   **kwargs)

    def gls_cov_update(self, zam, e, Vy_extra=None, kind='mass yield',
                       decay_data=None, sparse=False, rows=None,
                       threshold=None, **kwargs):
        """
        Update the prior IFY covariance matrix using the GLS technique
        described in https://doi.org/10.1016/j.anucene.2015.10.027
        .. math::
            $$
            V_{IFY_{post}} = V_{IFY_{prior}} - V_{IFY_{prior}}\cdot S.T \cdot \left(S\cdot V_{IFY_{prior}}\cdot S.T + V_{y_{extra}}\right)^{-1} \cdot S \cdot V_{IFY_{prior}}
            $$

        Parameters
        ----------
        zam : `int`
            ZAM number of the material to which calculations are to be
            applied.
        e : `float`
            Energy to which calculations are to be applied.
        Vy_extra : 2D iterable, optional
            Extra Covariance matrix (MXM).The default is None.
        kind : `str`, optional
            Keyword for obtaining sensitivity. The default is 'mass yield'.
        decay_data : `DecayData`, optional
            Object to change the model to CFY = Q*IFY or Ch_chain = S_chain*IFY,
            so the sensitivity (S) is Q or S_chain. The default is None,
            so the model is ChY = ChY_mass*IFY and the sensitivity is mass
            yield sensitivity.
        sparse : `bool`, optional
            Option to use sparse matrix for calculations. The default is False.
        rows : `int`, optional
            Option to use row calculation for matrix calculations. This option
            defines the number of lines to be taken into account in each loop.
            The default is None.
        threshold : `int`, optional
            Optional argument to avoid numerical fluctuations or
            values so small that they do not have to be taken into
            account. The default is None.
        kwargs : `dict`
            keyword arguments for method `get_decay_chains`

        Returns
        -------
        `sandy.CategoryCov`
            IFY covariance matrix performed by GLS for a given energy and zam.

        Notes
        -----
        .. note:: only option `kind='cumulative'` is implemented.

        Examples
        --------
        >>> zam = [591480, 591481, 601480]
        >>> decay_minimal = sandy.get_endf6_file("jeff_33", 'decay', zam)
        >>> decay_fytest = sandy.DecayData.from_endf6(decay_minimal)
        >>> CFY_var_extra = np.diag(pd.Series([1, 1, 1, 1]))
        >>> index = pd.Index([591480, 591481, 601480, 621480])
        >>> CFY_var_extra = pd.DataFrame(CFY_var_extra, index=index, columns=index)
        >>> npfy = Fy(minimal_fytest_2)
        >>> npfy.gls_cov_update(942390, 500e3, CFY_var_extra, kind='cumulative', decay_data=decay_fytest).data
        ZAP	         591480	         591481	        601480
        ZAP
        591480	 3.59798e-02   -3.07257e-03	   -6.45241e-04
        591481	-3.07257e-03	4.38149e-02	   -7.98869e-04
        601480	-6.45241e-04   -7.98869e-04	    9.83224e-03

        >>> npfy.gls_cov_update(942390, 500e3, CFY_var_extra, kind='cumulative', sparse=True, decay_data=decay_fytest).data
        ZAP	         591480	         591481	        601480
        ZAP
        591480	 3.59798e-02   -3.07257e-03	   -6.45241e-04
        591481	-3.07257e-03	4.38149e-02	   -7.98869e-04
        601480	-6.45241e-04   -7.98869e-04	    9.83224e-03

        >>> npfy.gls_cov_update(942390, 500e3, CFY_var_extra, kind='cumulative', rows=1, decay_data=decay_fytest).data
        ZAP	         591480	         591481	        601480
        ZAP
        591480	 3.59798e-02   -3.07257e-03	   -6.45241e-04
        591481	-3.07257e-03	4.38149e-02	   -7.98869e-04
        601480	-6.45241e-04   -7.98869e-04	    9.83224e-03

        >>> CFY_var_extra = np.diag(pd.Series([1, 1, 1]))
        >>> CFY_var_extra = pd.DataFrame(CFY_var_extra, index=zam, columns=zam)
        >>> npfy.gls_cov_update(942390, 500e3, CFY_var_extra, kind='cumulative', decay_data=decay_fytest).data
        ZAP	          591480	      591481	      601480
        ZAP
        591480	 3.71119e-02	-1.67096e-03	-3.50901e-04
        591481	-1.67096e-03	 4.55502e-02	-4.34448e-04
        601480	-3.50901e-04	-4.34448e-04	 9.90877e-03

        >>> npfy.gls_cov_update(942390, 500e3, CFY_var_extra, kind='cumulative', decay_data=decay_fytest, sparse=True).data
        ZAP	          591480	      591481	      601480
        ZAP
        591480	 3.71119e-02	-1.67096e-03	-3.50901e-04
        591481	-1.67096e-03	 4.55502e-02	-4.34448e-04
        601480	-3.50901e-04	-4.34448e-04	 9.90877e-03

        >>> npfy.gls_cov_update(942390, 500e3, CFY_var_extra, kind='cumulative', decay_data=decay_fytest, rows=1).data
        ZAP	          591480	      591481	      601480
        ZAP
        591480	 3.71119e-02	-1.67096e-03	-3.50901e-04
        591481	-1.67096e-03	 4.55502e-02	-4.34448e-04
        601480	-3.50901e-04	-4.34448e-04	 9.90877e-03

        >>> zam = [591480, 591481, 601480]
        >>> chain_var_extra = np.diag(pd.Series([1, 1, 1]))
        >>> index = pd.Index([147, 148, 149])
        >>> chain_var_extra = pd.DataFrame(chain_var_extra, index=index, columns=index)
        >>> npfy = sandy.Fy(minimal_fytest_2)
        >>> npfy.gls_cov_update(942390, 500e3, chain_var_extra).data
        ZAP	          591480	      591481	      601480
        ZAP
        591480	 3.85455e-02	-1.81818e-03	-3.63636e-04
        591481	-1.81818e-03	 4.77273e-02	-4.54545e-04
        601480	-3.63636e-04	-4.54545e-04	 9.90909e-03

        >>> npfy.gls_cov_update(942390, 500e3, chain_var_extra, sparse=True).data
        ZAP	          591480	      591481	      601480
        ZAP
        591480	 3.85455e-02	-1.81818e-03	-3.63636e-04
        591481	-1.81818e-03	 4.77273e-02	-4.54545e-04
        601480	-3.63636e-04	-4.54545e-04	 9.90909e-03

        >>> npfy.gls_cov_update(942390, 500e3, chain_var_extra, rows=1).data
        ZAP	          591480	      591481	      601480
        ZAP
        591480	 3.85455e-02	-1.81818e-03	-3.63636e-04
        591481	-1.81818e-03	 4.77273e-02	-4.54545e-04
        601480	-3.63636e-04	-4.54545e-04	 9.90909e-03

        >>> npfy.gls_cov_update(942390, 500e3)
        ZAP	          591480	      591481	      601480
        ZAP
        591480	 2.40000e-02	-2.00000e-02	-4.00000e-03
        591481	-2.00000e-02	 2.50000e-02	-5.00000e-03
        601480	-4.00000e-03	-5.00000e-03	 9.00000e-03        

        >>> npfy.gls_cov_update(942390, 500e3, sparse=True).data
        ZAP	          591480	      591481	      601480
        ZAP
        591480	 2.40000e-02	-2.00000e-02	-4.00000e-03
        591481	-2.00000e-02	 2.50000e-02	-5.00000e-03
        601480	-4.00000e-03	-5.00000e-03	 9.00000e-03

        >>> npfy.gls_cov_update(942390, 500e3, rows=1).data
        ZAP	          591480	      591481	      601480
        ZAP
        591480	 2.40000e-02	-2.00000e-02	-4.00000e-03
        591481	-2.00000e-02	 2.50000e-02	-5.00000e-03
        601480	-4.00000e-03	-5.00000e-03	 9.00000e-03
        """
        # Divide the data type:
        conditions = {'ZAM': zam, "E": e}
        fy_data = self._filters(conditions).data\
                      .set_index('ZAP')[['MT', 'DFY']]
        Vx_prior = fy_data.query('MT==454').DFY
        if rows is not None:
            Vx_prior = sandy.CategoryCov.from_var(Vx_prior, sparse=True).data
        else:
            Vx_prior = sandy.CategoryCov.from_var(Vx_prior, sparse=sparse).data
        # Fix the S with the correct dimension:
        if kind == 'mass yield':
            model_sensitivity_object = self._filters(conditions)
        elif kind == 'cumulative' or 'chain yield':
            model_sensitivity_object = decay_data
        S = _gls_setup(model_sensitivity_object, kind, **kwargs)
        index = Vx_prior.index
        Vx_prior = Vx_prior.reindex(index=S.columns, columns=S.columns).fillna(0)
        Vx_prior = sandy.CategoryCov(Vx_prior)
        # GLS covariance update
        if Vy_extra is not None:
            Vy_extra_ = pd.DataFrame(Vy_extra)
            S = S.reindex(index=Vy_extra_.index).fillna(0)
            Vx_post = Vx_prior.gls_update(S, Vy_extra_, sparse=sparse,
                                          rows=rows,
                                          threshold=threshold).data
        else:
            Vx_post = Vx_prior.constrained_gls_update(S, sparse=sparse,
                                                      rows=rows,
                                                      threshold=threshold).data
        Vx_post = Vx_post.reindex(index=index, columns=index).fillna(0)
        return sandy.CategoryCov(Vx_post)

    def update(self, zam, e, y_extra=None, Vy_extra=None, decay_data=None,
               kind='mass yield', sparse=False, rows=None,
               threshold=None, **kwargs):
        """
        Update the prior IFY and DIFY using the GLS technique for a given zam
        and energy. If the user does not enter Vy_extra and y_extra, the
        information taken from 'https://www-nds.iaea.org/endf349/la-ur-94-3106.pdf'
        (page 18-29) will be used in the options kind = 'mass yield' or
        kind = 'chain yield'.

        Notes
        -----
        ..note :: In option kind = 'cumulative', if Vy_extra and y_extra is not
        entered, GLS constrained will be used instead of GLS to update the
        IFY and DIFY.

        Parameters
        ----------
        zam : `int`
            ZAM number of the material to which calculations are to be
            applied.
        e : `float`
            Energy to which calculations are to be applied.
        y_extra : 1D iterable, optional
            New value of the vector.
        Vy_extra : 2D iterable, optional
            2D covariance matrix for y_extra (MXM).
        kind : `str`, optional
            Keyword for obtaining sensitivity. The default is 'mass yield'.
        decay_data : `DecayData`, optional
            Object to change the model to CFY = Q*IFY or Ch_chain = S_chain*IFY,
            so the sensitivity (S) is Q or S_chain. The default is None,
            so the model is ChY = ChY_mass*IFY and the sensitivity is mass
            yield sensitivity.
        sparse : `bool`, optional
            Option to use sparse matrix for calculations. The default is False.
        rows : `int`, optional
            Option to use row calculation for matrix calculations. This option
            defines the number of lines to be taken into account in each loop.
            The default is None.
        threshold : `int`, optional
            Optional argument to avoid numerical fluctuations or
            values so small that they do not have to be taken into
            account. The default is None.
        kwargs : `dict`
            keyword arguments for method `get_decay_chains`

        Returns
        -------
        `sandy.Fy`
            IFY updated with GLS for a given zam and energy.

        Examples
        --------
        >>> tape_nfpy = sandy.get_endf6_file("endfb_71",'nfpy','all')
        >>> nfpy = sandy.Fy.from_endf6(tape_nfpy)
        >>> e = 2.53000e-02
        >>> zam = 922350
        >>> conditions = {'E': e, 'ZAM': zam, 'MT': 454}
        >>> mass = nfpy.update(zam, e, sparse=True)._filters(conditions).data.set_index('ZAP').round(6)
        >>> mass.FY.loc[521340], mass.FY.loc[401000], mass.FY.loc[541380], mass.FY.loc[380950], mass.FY.loc[380940]
        (0.062211, 0.049734, 0.04811, 0.04531, 0.044811)

        >>> mass.DFY.loc[401000], mass.DFY.loc[390971], mass.DFY.loc[390970], mass.DFY.loc[400980], mass.DFY.loc[411020]
        (0.00561, 0.005398, 0.005398, 0.005058, 0.003176)
        """
        if y_extra is None and Vy_extra is None and kind != 'cumulative':
            ch_info = sandy.fy.import_chain_yields()
            ch_info = ch_info.loc[(ch_info.E == e) & (ch_info.ZAM == zam)]
            y_extra = ch_info.Ch/100
            Vy_extra = sandy.CategoryCov.from_var(ch_info.DCh).data/100
        return self.gls_update(zam, e, y_extra=y_extra,
                               Vy_extra=Vy_extra, kind=kind,
                               sparse=sparse, rows=rows, decay_data=decay_data,
                               threshold=threshold, **kwargs)

    def gls_update(self, zam, e, y_extra=None, Vy_extra=None,
                   kind='mass yield', decay_data=None, sparse=False, rows=None,
                   threshold=None, **kwargs):
        """
        Update IFY for a given zam, energy, decay_data and new information.
        .. math::
            $$
            IFY_{post} = IFY_{prior} + V_{IFY_{prior}}\cdot S.T \cdot \left(S\cdot V_{IFY_{prior}}\cdot S.T + V_{y_{extra}}\right)^{-1} \cdot \left(y_{extra} - y_{calc}\right)
            $$

        Parameters
        ----------
        zam : `int`
            ZAM number of the material to which calculations are to be
            applied.
        e : `float`
            Energy to which calculations are to be applied.
        y_extra : 1D iterable, optional
            New value of the vector.
        Vy_extra : 2D iterable, optional
            2D covariance matrix for y_extra (MXM).
        kind : `str`, optional
            Keyword for obtaining sensitivity. The default is 'mass yield'.
        decay_data : `DecayData`, optional
            Object to change the model to CFY = Q*IFY or Ch_chain = S_chain*IFY,
            so the sensitivity (S) is Q or S_chain. The default is None,
            so the model is ChY = ChY_mass*IFY and the sensitivity is mass
            yield sensitivity.
        sparse : `bool`, optional
            Option to use sparse matrix for calculations. The default is False.
        rows : `int`, optional
            Option to use row calculation for matrix calculations. This option
            defines the number of lines to be taken into account in each loop.
            The default is None.
        threshold : `int`, optional
            Optional argument to avoid numerical fluctuations or
            values so small that they do not have to be taken into
            account. The default is None.
        kwargs : `dict`
            keyword arguments for method `get_decay_chains`

         Returns
        -------
        `sandy.FY`
            IFY updated with GLS for a given zam, energy, decay_data and
            new information.

        Examples
        --------
        >>> zam = [591480, 591481, 601480]
        >>> decay_minimal = sandy.get_endf6_file("jeff_33", 'decay', zam)
        >>> decay_fytest = sandy.DecayData.from_endf6(decay_minimal)
        >>> CFY_extra = pd.Series([0, 0.1, 0.2], index=zam)
        >>> CFY_var_extra = np.diag(pd.Series([1, 1, 1]))
        >>> CFY_var_extra = pd.DataFrame(CFY_var_extra, index=zam, columns=zam)
        >>> npfy = Fy(minimal_fytest_2)
        >>> npfy.gls_update(942390, 500e3, CFY_extra, CFY_var_extra, kind='cumulative', decay_data=decay_fytest)
            MAT   MT     ZAM     ZAP           E          FY         DFY
        0  9437  459  942390  591480 5.00000e+05 8.00000e-01 4.00000e-02
        1  9437  459  942390  591481 5.00000e+05 1.00000e+00 5.00000e-02
        2  9437  459  942390  601480 5.00000e+05 2.00000e-01 1.00000e-02
        3  9437  454  942390  591480 5.00000e+05 8.24199e-02 3.71119e-02
        4  9437  454  942390  591481 5.00000e+05 1.78234e-01 4.55502e-02
        5  9437  454  942390  601480 5.00000e+05 2.96429e-01 9.90877e-03

        >>> npfy.gls_update(942390, 500e3, CFY_extra, CFY_var_extra, kind='cumulative', decay_data=decay_fytest, sparse=True)
            MAT   MT     ZAM     ZAP           E          FY         DFY
        0  9437  459  942390  591480 5.00000e+05 8.00000e-01 4.00000e-02
        1  9437  459  942390  591481 5.00000e+05 1.00000e+00 5.00000e-02
        2  9437  459  942390  601480 5.00000e+05 2.00000e-01 1.00000e-02
        3  9437  454  942390  591480 5.00000e+05 8.24199e-02 3.71119e-02
        4  9437  454  942390  591481 5.00000e+05 1.78234e-01 4.55502e-02
        5  9437  454  942390  601480 5.00000e+05 2.96429e-01 9.90877e-03

        >>> npfy.gls_update(942390, 500e3, CFY_extra, CFY_var_extra, kind='cumulative', decay_data=decay_fytest, rows=1)
            MAT   MT     ZAM     ZAP           E          FY         DFY
        0  9437  459  942390  591480 5.00000e+05 8.00000e-01 4.00000e-02
        1  9437  459  942390  591481 5.00000e+05 1.00000e+00 5.00000e-02
        2  9437  459  942390  601480 5.00000e+05 2.00000e-01 1.00000e-02
        3  9437  454  942390  591480 5.00000e+05 8.24199e-02 3.71119e-02
        4  9437  454  942390  591481 5.00000e+05 1.78234e-01 4.55502e-02
        5  9437  454  942390  601480 5.00000e+05 2.96429e-01 9.90877e-03

        >>> index = pd.Index([147, 148, 149])
        >>> chain_extra = pd.Series([0, 0.1, 0.2], index=index)
        >>> chain_var_extra = np.diag(pd.Series([1, 1, 1]))
        >>> chain_var_extra = pd.DataFrame(chain_var_extra, index=index, columns=index)
        >>> npfy.gls_update(942390, 500e3, chain_extra, chain_var_extra)
            MAT   MT     ZAM     ZAP           E          FY         DFY
        0  9437  459  942390  591480 5.00000e+05 8.00000e-01 4.00000e-02
        1  9437  459  942390  591481 5.00000e+05 1.00000e+00 5.00000e-02
        2  9437  459  942390  601480 5.00000e+05 2.00000e-01 1.00000e-02
        3  9437  454  942390  591480 5.00000e+05 8.18182e-02 3.85455e-02
        4  9437  454  942390  591481 5.00000e+05 1.77273e-01 4.77273e-02
        5  9437  454  942390  601480 5.00000e+05 2.95455e-01 9.90909e-03

        >>> npfy.gls_update(942390, 500e3, chain_extra, chain_var_extra, sparse=True)
            MAT   MT     ZAM     ZAP           E          FY         DFY
        0  9437  459  942390  591480 5.00000e+05 8.00000e-01 4.00000e-02
        1  9437  459  942390  591481 5.00000e+05 1.00000e+00 5.00000e-02
        2  9437  459  942390  601480 5.00000e+05 2.00000e-01 1.00000e-02
        3  9437  454  942390  591480 5.00000e+05 8.18182e-02 3.85455e-02
        4  9437  454  942390  591481 5.00000e+05 1.77273e-01 4.77273e-02
        5  9437  454  942390  601480 5.00000e+05 2.95455e-01 9.90909e-03

        >>> npfy.gls_update(942390, 500e3, chain_extra, chain_var_extra, rows=1)
            MAT   MT     ZAM     ZAP           E          FY         DFY
        0  9437  459  942390  591480 5.00000e+05 8.00000e-01 4.00000e-02
        1  9437  459  942390  591481 5.00000e+05 1.00000e+00 5.00000e-02
        2  9437  459  942390  601480 5.00000e+05 2.00000e-01 1.00000e-02
        3  9437  454  942390  591480 5.00000e+05 8.18182e-02 3.85455e-02
        4  9437  454  942390  591481 5.00000e+05 1.77273e-01 4.77273e-02
        5  9437  454  942390  601480 5.00000e+05 2.95455e-01 9.90909e-03
        """
        data = self.data.copy()
        # Filter FY data:
        conditions = {'ZAM': zam, "E": e}
        fy_data = self._filters(conditions).data
        mat = fy_data.MAT.iloc[0]
        fy_data = fy_data.set_index('ZAP')[['MT', 'FY', 'DFY']]
        # Divide the data:
        mask = (data.ZAM == zam) & (data.MT == 454) & (data.E == e)
        data = data.loc[~mask]
        x_prior = fy_data.query('MT==454').FY
        Vx_prior = fy_data.query('MT==454').DFY
        index = x_prior.index
        if rows is not None:
            Vx_prior = sandy.CategoryCov.from_var(Vx_prior, sparse=True).data
        else:
            Vx_prior = sandy.CategoryCov.from_var(Vx_prior, sparse=sparse).data
        # Find the GLS varibles:
        if kind == 'mass yield':
            model_sensitivity_object = self._filters(conditions)
        elif kind == 'cumulative' or 'chain yield':
            model_sensitivity_object = decay_data
        S = _gls_setup(model_sensitivity_object, kind, **kwargs)
        # Perform GLS:
        if y_extra is None and Vy_extra is None:
            x_post = sandy.constrained_gls_update(x_prior, S, Vx_prior,
                                                  sparse=sparse, rows=rows,
                                                  threshold=threshold)
            Vx_post = self.gls_cov_update(zam, e, kind=kind,
                                          decay_data=decay_data, sparse=sparse,
                                          rows=rows, threshold=threshold).data
        else:
            x_post = sandy.gls_update(x_prior, S, Vx_prior, Vy_extra, y_extra,
                                      sparse=sparse, rows=rows,
                                      threshold=threshold)
            Vx_post = self.gls_cov_update(zam, e, Vy_extra,
                                          kind=kind, decay_data=decay_data,
                                          sparse=sparse, rows=rows,
                                          threshold=threshold).data
        # Results in appropriate format:
        Vx_post = Vx_post.reindex(index=index, columns=index).fillna(0)
        x_post = x_post.reindex(index).fillna(0)
        calc_values = x_post.rename('FY').reset_index()
        calc_values['DFY'] = pd.Series(np.diag(Vx_post)).values
        calc_values[['MAT', 'ZAM', 'MT', 'E']] = [mat, zam, 454, e]
        data = pd.concat([data, calc_values], ignore_index=True)
        return self.__class__(data)

    def ishikawa_factor(self, zam, e, Vy_extra=None,
                        kind='mass yield', decay_data=None, sparse=False,
                        rows=None, threshold=None, **kwargs):
        """
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
        sparse : `bool`, optional
            Option to use sparse matrix for calculations. The default is False.
        rows : `int`, optional
            Option to use row calculation for matrix calculations. This option
            defines the number of lines to be taken into account in each loop.
            The default is None.
        threshold : `int`, optional
            Optional argument to avoid numerical fluctuations or
            values so small that they do not have to be taken into
            account. The default is None.
        kwargs : `dict`
            keyword arguments for method `get_decay_chains`

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
        >>> zam = [591480, 591481, 601480]
        >>> decay_minimal = sandy.get_endf6_file("jeff_33", 'decay', zam)
        >>> decay_fytest = sandy.DecayData.from_endf6(decay_minimal)
        >>> CFY_var_extra = np.diag(pd.Series([1, 1, 1]))
        >>> CFY_var_extra = pd.DataFrame(CFY_var_extra, index=zam, columns=zam)
        >>> npfy = sandy.Fy(minimal_fytest_2)
        >>> npfy.ishikawa_factor(942390, 500e3, Vy_extra=CFY_var_extra, kind='cumulative', decay_data=decay_fytest)
        591480   4.00000e-02
        591481   5.00000e-02
        601480   1.00000e-01
        dtype: float64
        
        >>> npfy.ishikawa_factor(942390, 500e3, Vy_extra=CFY_var_extra, kind='cumulative', decay_data=decay_fytest, sparse=True)
        591480   4.00000e-02
        591481   5.00000e-02
        601480   1.00000e-01
        dtype: float64

        >>> npfy.ishikawa_factor(942390, 500e3, Vy_extra=CFY_var_extra, kind='cumulative', decay_data=decay_fytest, rows=1)
        591480   4.00000e-02
        591481   5.00000e-02
        601480   1.00000e-01
        dtype: float64

        >>> A = [147, 148, 149]
        >>> Chain_var_extra = np.diag(pd.Series([1, 1, 1]))
        >>> Chain_var_extra = pd.DataFrame(Chain_var_extra, index=A, columns=A)
        >>> npfy = sandy.Fy(minimal_fytest_2)
        >>> npfy.ishikawa_factor(942390, 500e3, Vy_extra=Chain_var_extra)
        147   0.00000e+00
        148   1.00000e-01
        149   0.00000e+00
        dtype: float64

        >>> npfy.ishikawa_factor(942390, 500e3,Vy_extra= Chain_var_extra, sparse=True)
        147   0.00000e+00
        148   1.00000e-01
        149   0.00000e+00
        dtype: float64

        >>> npfy.ishikawa_factor(942390, 500e3, Vy_extra=Chain_var_extra, rows=1)
        147   0.00000e+00
        148   1.00000e-01
        149   0.00000e+00
        dtype: float64
        """
        if Vy_extra is None and kind != 'cumulative':
            ch_info = sandy.fy.import_chain_yields()
            ch_info = ch_info.loc[(ch_info.E == e) & (ch_info.ZAM == zam)]
            Vy_extra = sandy.CategoryCov.from_var(ch_info.DCh).data/100
        return self._ishikawa_factor(zam, e, Vy_extra=Vy_extra,
                                     kind=kind, decay_data=decay_data,
                                     sparse=sparse, rows=rows,
                                     threshold=threshold, **kwargs)

    def _ishikawa_factor(self, zam, e, Vy_extra,
                        kind='mass yield', decay_data=None, sparse=False,
                        rows=None, threshold=None, **kwargs):
        """
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
        sparse : `bool`, optional
            Option to use sparse matrix for calculations. The default is False.
        rows : `int`, optional
            Option to use row calculation for matrix calculations. This option
            defines the number of lines to be taken into account in each loop.
            The default is None.
        threshold : `int`, optional
            Optional argument to avoid numerical fluctuations or
            values so small that they do not have to be taken into
            account. The default is None.
        kwargs : `dict`
            keyword arguments for method `get_decay_chains`

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
        >>> zam = [591480, 591481, 601480]
        >>> decay_minimal = sandy.get_endf6_file("jeff_33", 'decay', zam)
        >>> decay_fytest = sandy.DecayData.from_endf6(decay_minimal)
        >>> CFY_var_extra = np.diag(pd.Series([1, 1, 1]))
        >>> CFY_var_extra = pd.DataFrame(CFY_var_extra, index=zam, columns=zam)
        >>> npfy = sandy.Fy(minimal_fytest_2)
        >>> npfy._ishikawa_factor(942390, 500e3, CFY_var_extra, kind='cumulative', decay_data=decay_fytest)
        591480   4.00000e-02
        591481   5.00000e-02
        601480   1.00000e-01
        dtype: float64
        
        >>> npfy._ishikawa_factor(942390, 500e3, CFY_var_extra, kind='cumulative', decay_data=decay_fytest, sparse=True)
        591480   4.00000e-02
        591481   5.00000e-02
        601480   1.00000e-01
        dtype: float64

        >>> npfy._ishikawa_factor(942390, 500e3, CFY_var_extra, kind='cumulative', decay_data=decay_fytest, rows=1)
        591480   4.00000e-02
        591481   5.00000e-02
        601480   1.00000e-01
        dtype: float64

        >>> A = [147, 148, 149]
        >>> Chain_var_extra = np.diag(pd.Series([1, 1, 1]))
        >>> Chain_var_extra = pd.DataFrame(Chain_var_extra, index=A, columns=A)
        >>> npfy = sandy.Fy(minimal_fytest_2)
        >>> npfy._ishikawa_factor(942390, 500e3, Chain_var_extra)
        147   0.00000e+00
        148   1.00000e-01
        149   0.00000e+00
        dtype: float64

        >>> npfy._ishikawa_factor(942390, 500e3, Chain_var_extra, sparse=True)
        147   0.00000e+00
        148   1.00000e-01
        149   0.00000e+00
        dtype: float64

        >>> npfy._ishikawa_factor(942390, 500e3, Chain_var_extra, rows=1)
        147   0.00000e+00
        148   1.00000e-01
        149   0.00000e+00
        dtype: float64
        """
        # Filter FY data:
        conditions = {'ZAM': zam, "E": e}
        fy_data = self._filters(conditions).data.set_index('ZAP')
        Vx_prior = fy_data.query('MT==454').DFY
        if rows is not None:
            Vx_prior = sandy.CategoryCov.from_var(Vx_prior, sparse=True).data
        else:
            Vx_prior = sandy.CategoryCov.from_var(Vx_prior, sparse=sparse).data
        # Find the GLS sensitivity:
        if kind == 'mass yield':
            model_sensitivity_object = self._filters(conditions)
        elif kind == 'cumulative' or 'chain yield':
            model_sensitivity_object = decay_data
        S = _gls_setup(model_sensitivity_object, kind, **kwargs)
        # Perform Ishikawa factor:
        ishikawa = sandy.ishikawa_factor(S, Vx_prior, Vy_extra, sparse=sparse,
                                         rows=rows)
        if threshold is not None:
            ishikawa[abs(ishikawa) < threshold] = 0
        return ishikawa

    def chi_diagonal(self, zam, e, Vy_extra=None, y_extra=None,
                     kind='mass yield', decay_data=None, sparse=False,
                     rows=None, threshold=None, **kwargs):
        """
        Function to calculate diagonal chi-value and know if may exist a
        inconsistency between |y_extra - y_calc_by_model| and the covariance
        matrices, S*Vx_prior*S, Vy_extra and Vx_posterior.

        Parameters
        ----------
        zam : `int`
            ZAM number of the material to which calculations are to be
            applied.
        e : `float`
            Energy to which calculations are to be applied.
        Vy_extra : 2D iterable
            2D covariance matrix for y_extra (MXM).
        y_extra : TYPE
            DESCRIPTION.
        kind : `str`, optional
            Keyword for obtaining sensitivity. The default is 'mass yield'.
        decay_data : `DecayData`, optional
            Object to change the model to CFY = Q*IFY or
            Ch_chain = S_chain*IFY, so the sensitivity (S) is Q or S_chain.
            The default is None, so the model is ChY = ChY_mass*IFY
            and the sensitivity is mass yield sensitivity.
        sparse : `bool`, optional
            Option to use sparse matrix for calculations. The default is False.
        rows : `int`, optional
            Option to use row calculation for matrix calculations. This option
            defines the number of lines to be taken into account in each loop.
            The default is None.
        threshold : `int`, optional
            Optional argument to avoid numerical fluctuations or
            values so small that they do not have to be taken into
            account. The default is None.
        kwargs : `dict`
            keyword arguments for method `get_decay_chains`

        Returns
        -------
        `pd.Series`
            diagonal chi-value

        Results:
        -------
        chi diagonal >> 1 :
            Inconsistency may exist between |y_extra - y_calc| and covariance
            matrix, S*Vx_prior*S.T, and Vy_extra.
        Example
        -------
        >>> zam = [591480, 591481, 601480]
        >>> decay_minimal = sandy.get_endf6_file("jeff_33", 'decay', zam)
        >>> decay_fytest = sandy.DecayData.from_endf6(decay_minimal)
        >>> CFY_var_extra = np.diag(pd.Series([1, 1, 1]))
        >>> CFY_var_extra = pd.DataFrame(CFY_var_extra, index=zam, columns=zam)
        >>> y_extra = pd.Series([1, 1, 1], index=zam)
        >>> npfy = sandy.Fy(minimal_fytest_2)
        >>> npfy._chi_diagonal(942390, 500e3, Vy_extra=CFY_var_extra, y_extra=y_extra, kind='cumulative', decay_data=decay_fytest)
        591480   9.17180e-01
        591481   8.18867e-01
        601480   4.18775e-01
        dtype: float64

        >>> npfy._chi_diagonal(942390, 500e3, Vy_extra=CFY_var_extra, y_extra=y_extra, kind='cumulative', decay_data=decay_fytest, sparse=True)
        591480   9.17180e-01
        591481   8.18867e-01
        601480   4.18775e-01
        dtype: float64

        >>> npfy._chi_diagonal(942390, 500e3, Vy_extra=CFY_var_extra, y_extra=y_extra, kind='cumulative', decay_data=decay_fytest, rows=1)
        591480   9.17180e-01
        591481   8.18867e-01
        601480   4.18775e-01
        dtype: float64

        >>> A = [147, 148, 149]
        >>> Chain_var_extra = np.diag(pd.Series([1, 1, 1]))
        >>> Chain_var_extra = pd.DataFrame(Chain_var_extra, index=A, columns=A)
        >>> npfy = sandy.Fy(minimal_fytest_2)
        >>> y_extra = pd.Series([1, 1, 1], index=A)
        >>> npfy._chi_diagonal(942390, 500e3, Vy_extra=Chain_var_extra, y_extra=y_extra)
        147   1.00000e+00
        148   4.19524e-01
        149   1.00000e+00
        dtype: float64

        >>> npfy._chi_diagonal(942390, 500e3, Vy_extra=Chain_var_extra, y_extra=y_extra, sparse=True)
        147   1.00000e+00
        148   4.19524e-01
        149   1.00000e+00
        dtype: float64

        >>> npfy.chi_diagonal(942390, 500e3, Vy_extra=Chain_var_extra, y_extra=y_extra, rows=1)
        147   1.00000e+00
        148   4.19524e-01
        149   1.00000e+00
        dtype: float64
        """
        if y_extra is None and Vy_extra is None and kind != 'cumulative':
            ch_info = sandy.fy.import_chain_yields()
            ch_info = ch_info.loc[(ch_info.E == e) & (ch_info.ZAM == zam)]
            y_extra = ch_info.Ch/100
            Vy_extra = sandy.CategoryCov.from_var(ch_info.DCh).data/100
        return self._chi_diagonal(zam, e, Vy_extra=Vy_extra, y_extra=y_extra,
                                  kind=kind, decay_data=decay_data,
                                  sparse=sparse, rows=rows,
                                  threshold=threshold, **kwargs)

    def _chi_diagonal(self, zam, e, Vy_extra, y_extra,
                     kind='mass yield', decay_data=None, sparse=False,
                     rows=None, threshold=None, **kwargs):
        """
        Function to calculate diagonal chi-value and know if may exist a
        inconsistency between |y_extra - y_calc_by_model| and the covariance
        matrices, S*Vx_prior*S, Vy_extra and Vx_posterior.

        Parameters
        ----------
        zam : `int`
            ZAM number of the material to which calculations are to be
            applied.
        e : `float`
            Energy to which calculations are to be applied.
        Vy_extra : 2D iterable
            2D covariance matrix for y_extra (MXM).
        y_extra : TYPE
            DESCRIPTION.
        kind : `str`, optional
            Keyword for obtaining sensitivity. The default is 'mass yield'.
        decay_data : `DecayData`, optional
            Object to change the model to CFY = Q*IFY or
            Ch_chain = S_chain*IFY, so the sensitivity (S) is Q or S_chain.
            The default is None, so the model is ChY = ChY_mass*IFY
            and the sensitivity is mass yield sensitivity.
        sparse : `bool`, optional
            Option to use sparse matrix for calculations. The default is False.
        rows : `int`, optional
            Option to use row calculation for matrix calculations. This option
            defines the number of lines to be taken into account in each loop.
            The default is None.
        threshold : `int`, optional
            Optional argument to avoid numerical fluctuations or
            values so small that they do not have to be taken into
            account. The default is None.
        kwargs : `dict`
            keyword arguments for method `get_decay_chains`

        Returns
        -------
        `pd.Series`
            diagonal chi-value

        Results:
        -------
        chi diagonal >> 1 :
            Inconsistency may exist between |y_extra - y_calc| and covariance
            matrix, S*Vx_prior*S.T, and Vy_extra.
        Example
        -------
        >>> zam = [591480, 591481, 601480]
        >>> decay_minimal = sandy.get_endf6_file("jeff_33", 'decay', zam)
        >>> decay_fytest = sandy.DecayData.from_endf6(decay_minimal)
        >>> CFY_var_extra = np.diag(pd.Series([1, 1, 1]))
        >>> CFY_var_extra = pd.DataFrame(CFY_var_extra, index=zam, columns=zam)
        >>> y_extra = pd.Series([1, 1, 1], index=zam)
        >>> npfy = sandy.Fy(minimal_fytest_2)
        >>> npfy._chi_diagonal(942390, 500e3, CFY_var_extra, y_extra, kind='cumulative', decay_data=decay_fytest)
        591480   9.17180e-01
        591481   8.18867e-01
        601480   4.18775e-01
        dtype: float64

        >>> npfy._chi_diagonal(942390, 500e3, CFY_var_extra, y_extra, kind='cumulative', decay_data=decay_fytest, sparse=True)
        591480   9.17180e-01
        591481   8.18867e-01
        601480   4.18775e-01
        dtype: float64

        >>> npfy._chi_diagonal(942390, 500e3, CFY_var_extra, y_extra, kind='cumulative', decay_data=decay_fytest, rows=1)
        591480   9.17180e-01
        591481   8.18867e-01
        601480   4.18775e-01
        dtype: float64

        >>> A = [147, 148, 149]
        >>> Chain_var_extra = np.diag(pd.Series([1, 1, 1]))
        >>> Chain_var_extra = pd.DataFrame(Chain_var_extra, index=A, columns=A)
        >>> npfy = sandy.Fy(minimal_fytest_2)
        >>> y_extra = pd.Series([1, 1, 1], index=A)
        >>> npfy._chi_diagonal(942390, 500e3, Chain_var_extra, y_extra)
        147   1.00000e+00
        148   4.19524e-01
        149   1.00000e+00
        dtype: float64

        >>> npfy._chi_diagonal(942390, 500e3, Chain_var_extra, y_extra, sparse=True)
        147   1.00000e+00
        148   4.19524e-01
        149   1.00000e+00
        dtype: float64

        >>> npfy._chi_diagonal(942390, 500e3, Chain_var_extra, y_extra, rows=1)
        147   1.00000e+00
        148   4.19524e-01
        149   1.00000e+00
        dtype: float64
        """
        # Filter FY data:
        conditions = {'ZAM': zam, "E": e}
        fy_data = self._filters(conditions).data.set_index('ZAP')
        x_prior = fy_data.query('MT==454').FY
        Vx_prior = fy_data.query('MT==454').DFY
        if rows is not None:
            Vx_prior = sandy.CategoryCov.from_var(Vx_prior, sparse=True).data
        else:
            Vx_prior = sandy.CategoryCov.from_var(Vx_prior, sparse=sparse).data
        # Find the GLS sensitivity:
        if kind == 'mass yield':
            model_sensitivity_object = self._filters(conditions)
        elif kind == 'cumulative' or 'chain yield':
            model_sensitivity_object = decay_data
        S = _gls_setup(model_sensitivity_object, kind, **kwargs)
        # Find the diagonal chi-value
        chi_diagonal = sandy.chi_diag(x_prior, S, Vx_prior, Vy_extra, y_extra,
                                      sparse=sparse, rows=rows)
        if threshold is not None:
            chi_diagonal[abs(chi_diagonal) < threshold] = 0
        return chi_diagonal

    def chi_square(self, zam, e, Vy_extra=None, y_extra=None, N_e=1,
                   kind='mass yield', decay_data=None, sparse=False, rows=None,
                   threshold=None, **kwargs):
        """
        Method to characterises the contribution of the experiment to the
        'a posteriori' adjustment.

        Parameters
        ----------
        zam : `int`
            ZAM number of the material to which calculations are to be
            applied.
        e : `float`
            Energy to which calculations are to be applied.
        Vy_extra : 2D iterable
            2D covariance matrix for y_extra (MXM).
        y_extra : 1D iterable
            1D extra info on output (NX1)
        N_e : `int`
            Number of experimental values used in adjustment.The default is 1.
        kind : `str`, optional
            Keyword for obtaining sensitivity. The default is 'mass yield'.
        decay_data : `DecayData`, optional
            Object to change the model to CFY = Q*IFY or
            Ch_chain = S_chain*IFY, so the sensitivity (S) is Q or S_chain.
            The default is None, so the model is ChY = ChY_mass*IFY
            and the sensitivity is mass yield sensitivity.
        sparse : `bool`, optional
            Option to use sparse matrix for calculations. The default is False.
        rows : `int`, optional
            Option to use row calculation for matrix calculations. This option
            defines the number of lines to be taken into account in each loop.
            The default is None.
        threshold : `int`, optional
            Optional argument to avoid numerical fluctuations or
            values so small that they do not have to be taken into
            account. The default is None.
        kwargs : `dict`
            keyword arguments for method `get_decay_chains`

        Returns
        -------
        `pd.Series`
            contribution to chi-square value
    
        Results:
        -------
        chi square < 0 :
            The experiment is very effective in the adjustment.
    
        Example
        -------
        >>> zam = [591480, 591481, 601480]
        >>> decay_minimal = sandy.get_endf6_file("jeff_33", 'decay', zam)
        >>> decay_fytest = sandy.DecayData.from_endf6(decay_minimal)
        >>> CFY_var_extra = np.diag(pd.Series([1, 1, 1]))
        >>> CFY_var_extra = pd.DataFrame(CFY_var_extra, index=zam, columns=zam)
        >>> y_extra = pd.Series([1, 1, 1], index=zam)
        >>> npfy = sandy.Fy(minimal_fytest_2)
        >>> Ne = 1
        >>> npfy.chi_square(942390, 500e3, Vy_extra=CFY_var_extra, y_extra=y_extra, N_e=Ne, kind='cumulative', decay_data=decay_fytest)
        591480   7.68510e-01
        591481   5.98149e-01
        601480   1.19440e-01
        dtype: float64

        >>> npfy.chi_square(942390, 500e3, Vy_extra=CFY_var_extra, y_extra=y_extra, N_e=Ne, kind='cumulative', decay_data=decay_fytest, sparse=True)
        591480   7.68510e-01
        591481   5.98149e-01
        601480   1.19440e-01
        dtype: float64

        >>> npfy.chi_square(942390, 500e3, Vy_extra=CFY_var_extra, y_extra=y_extra, N_e=Ne, kind='cumulative', decay_data=decay_fytest, rows=1)
        591480   7.68510e-01
        591481   5.98149e-01
        601480   1.19440e-01
        dtype: float64

        >>> A = [147, 148, 149]
        >>> Chain_var_extra = np.diag(pd.Series([1, 1, 1]))
        >>> Chain_var_extra = pd.DataFrame(Chain_var_extra, index=A, columns=A)
        >>> npfy = sandy.Fy(minimal_fytest_2)
        >>> y_extra = pd.Series([1, 1, 1], index=A)
        >>> npfy.chi_square(942390, 500e3, Vy_extra=Chain_var_extra, y_extra=y_extra, N_e=Ne)
        147   1.00000e+00
        148   1.45455e-01
        149   1.00000e+00
        dtype: float64

        >>> npfy.chi_square(942390, 500e3, Vy_extra=Chain_var_extra, y_extra=y_extra, N_e=Ne, sparse=True)
        147   1.00000e+00
        148   1.45455e-01
        149   1.00000e+00
        dtype: float64

        >>> npfy.chi_square(942390, 500e3, Vy_extra=Chain_var_extra, y_extra=y_extra, N_e=Ne, rows=1)
        147   1.00000e+00
        148   1.45455e-01
        149   1.00000e+00
        dtype: float64
        """
        if y_extra is None and Vy_extra is None and kind != 'cumulative':
            ch_info = sandy.fy.import_chain_yields()
            ch_info = ch_info.loc[(ch_info.E == e) & (ch_info.ZAM == zam)]
            y_extra = ch_info.Ch/100
            Vy_extra = sandy.CategoryCov.from_var(ch_info.DCh).data/100
        return self._chi_square(zam, e, Vy_extra=Vy_extra, y_extra=y_extra,
                                N_e=N_e, kind=kind, decay_data=decay_data,
                                sparse=sparse, rows=rows, threshold=threshold,
                                **kwargs)

    def _chi_square(self, zam, e, Vy_extra, y_extra, N_e,
                   kind='mass yield', decay_data=None, sparse=False, rows=None,
                   threshold=None, **kwargs):
        """
        Method to characterises the contribution of the experiment to the
        'a posteriori' adjustment.

        Parameters
        ----------
        zam : `int`
            ZAM number of the material to which calculations are to be
            applied.
        e : `float`
            Energy to which calculations are to be applied.
        Vy_extra : 2D iterable
            2D covariance matrix for y_extra (MXM).
        y_extra : 1D iterable
            1D extra info on output (NX1)
        N_e : `int`
            Number of experimental values used in adjustment.
        kind : `str`, optional
            Keyword for obtaining sensitivity. The default is 'mass yield'.
        decay_data : `DecayData`, optional
            Object to change the model to CFY = Q*IFY or
            Ch_chain = S_chain*IFY, so the sensitivity (S) is Q or S_chain.
            The default is None, so the model is ChY = ChY_mass*IFY
            and the sensitivity is mass yield sensitivity.
        sparse : `bool`, optional
            Option to use sparse matrix for calculations. The default is False.
        rows : `int`, optional
            Option to use row calculation for matrix calculations. This option
            defines the number of lines to be taken into account in each loop.
            The default is None.
        threshold : `int`, optional
            Optional argument to avoid numerical fluctuations or
            values so small that they do not have to be taken into
            account. The default is None.
        kwargs : `dict`
            keyword arguments for method `get_decay_chains`

        Returns
        -------
        `pd.Series`
            contribution to chi-square value
    
        Results:
        -------
        chi square < 0 :
            The experiment is very effective in the adjustment.
    
        Example
        -------
        >>> zam = [591480, 591481, 601480]
        >>> decay_minimal = sandy.get_endf6_file("jeff_33", 'decay', zam)
        >>> decay_fytest = sandy.DecayData.from_endf6(decay_minimal)
        >>> CFY_var_extra = np.diag(pd.Series([1, 1, 1]))
        >>> CFY_var_extra = pd.DataFrame(CFY_var_extra, index=zam, columns=zam)
        >>> y_extra = pd.Series([1, 1, 1], index=zam)
        >>> npfy = sandy.Fy(minimal_fytest_2)
        >>> Ne = 1
        >>> npfy._chi_square(942390, 500e3, CFY_var_extra, y_extra, Ne, kind='cumulative', decay_data=decay_fytest)
        591480   7.68510e-01
        591481   5.98149e-01
        601480   1.19440e-01
        dtype: float64

        >>> npfy._chi_square(942390, 500e3, CFY_var_extra, y_extra, Ne, kind='cumulative', decay_data=decay_fytest, sparse=True)
        591480   7.68510e-01
        591481   5.98149e-01
        601480   1.19440e-01
        dtype: float64

        >>> npfy._chi_square(942390, 500e3, CFY_var_extra, y_extra, Ne, kind='cumulative', decay_data=decay_fytest, rows=1)
        591480   7.68510e-01
        591481   5.98149e-01
        601480   1.19440e-01
        dtype: float64

        >>> A = [147, 148, 149]
        >>> Chain_var_extra = np.diag(pd.Series([1, 1, 1]))
        >>> Chain_var_extra = pd.DataFrame(Chain_var_extra, index=A, columns=A)
        >>> npfy = sandy.Fy(minimal_fytest_2)
        >>> y_extra = pd.Series([1, 1, 1], index=A)
        >>> npfy._chi_square(942390, 500e3, Chain_var_extra, y_extra, Ne)
        147   1.00000e+00
        148   1.45455e-01
        149   1.00000e+00
        dtype: float64

        >>> npfy._chi_square(942390, 500e3, Chain_var_extra, y_extra, Ne, sparse=True)
        147   1.00000e+00
        148   1.45455e-01
        149   1.00000e+00
        dtype: float64

        >>> npfy._chi_square(942390, 500e3, Chain_var_extra, y_extra, Ne, rows=1)
        147   1.00000e+00
        148   1.45455e-01
        149   1.00000e+00
        dtype: float64
        """
        # Filter FY data:
        conditions = {'ZAM': zam, "E": e}
        fy_data = self._filters(conditions).data.set_index('ZAP')
        x_prior = fy_data.query('MT==454').FY
        Vx_prior = fy_data.query('MT==454').DFY
        if rows is not None:
            Vx_prior = sandy.CategoryCov.from_var(Vx_prior, sparse=True).data
        else:
            Vx_prior = sandy.CategoryCov.from_var(Vx_prior, sparse=sparse).data
        # Find the GLS sensitivity:
        if kind == 'mass yield':
            model_sensitivity_object = self._filters(conditions)
        elif kind == 'cumulative' or 'chain yield':
            model_sensitivity_object = decay_data
        S = _gls_setup(model_sensitivity_object, kind, **kwargs)
        # Find the diagonal chi-value
        chi_square = sandy.chi_square(x_prior, S, Vx_prior, Vy_extra, y_extra,
                                      N_e, sparse=sparse, rows=rows)
        if threshold is not None:
            chi_square[abs(chi_square) < threshold] = 0
        return chi_square

    def _filters(self, conditions):
        """
        Apply several condition to source data and return filtered results.

        Parameters
        ----------
        conditions : `dict`
            Conditions to apply, where the key is any label present in the
            columns of `data` and the value is filtering condition.

        Returns
        -------
        `sandy.Fy`
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
        """
        Apply condition to source data and return filtered results.

        Parameters
        ----------
        `key` : `str`
            any label present in the columns of `data`
        `value` : `int` or `float`
            value used as filtering condition

        Returns
        -------
        `sandy.Fy`
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
        """
        Extract fission yields from `Endf6` instance.

        Parameters
        ----------
        endf6 : `sandy.Endf6`
            object containing the ENDF-6 text
        verbose : `bool`, optional, default is `False`
            flag to print information when reading ENDF-6 file

        Returns
        -------
        `sandy.Fy`
            fission yield object

        Notes
        -----
        .. note:: Both independent and cumulative fission product yields are
                  loaded, if found.
        """
        tape = endf6.filter_by(listmf=[8], listmt=[454, 459])
        data = []
        for mat, mf, mt in tape.data:
            sec = tape.read_section(mat, mf, mt)
            zam = int(sec["ZA"]*10)
            if verbose:
                logging.info(f"reading 'ZAM={zam}'...")
            for e in sec["E"]:
                for zap in sec["E"][e]["ZAP"]:
                    fy = sec["E"][e]["ZAP"][zap]["FY"]
                    dfy = sec["E"][e]["ZAP"][zap]["DFY"]
                    values = (mat, mt, zam, zap, e, fy, dfy)
                    data.append(dict(zip(cls._columns, values)))
        df = pd.DataFrame(data)
        return cls(df)

    def to_endf6(self, endf6):
        """
        Update cross sections in `Endf6` instance with those available in a
        `Fy` instance.

        .. warning:: only IFY and CFY that are originally
                     present in the `Endf6` instance are modified

        Parameters
        ----------
        `endf6` : `sandy.Endf6`
            `Endf6` instance

        Returns
        -------
        `sandy.Endf6`
            `Endf6` instance with updated IFY and CFY

        Examples
        --------
        >>> tape = sandy.get_endf6_file("jeff_33", "nfpy", "all")
        >>> from_endf = sandy.sections.mf8.read_mf8(tape, 9228, 454)
        >>> text = sandy.sections.mf8._write_fy(from_endf)
        >>> fy = sandy.Fy.from_endf6(tape)
        >>> new_tape = fy.to_endf6(tape)
        >>> new_from_endf = sandy.sections.mf8.read_mf8(new_tape, 9228, 454)
        >>> new_text = sandy.sections.mf8._write_fy(new_from_endf)
        >>> assert new_text == text

        >>> tape = sandy.get_endf6_file("jeff_33", "nfpy", "all")
        >>> from_endf = sandy.sections.mf8.read_mf8(tape, 9228, 459)
        >>> text = sandy.sections.mf8._write_fy(from_endf)
        >>> fy = sandy.Fy.from_endf6(tape)
        >>> new_tape = fy.to_endf6(tape)
        >>> new_from_endf = sandy.sections.mf8.read_mf8(new_tape, 9228, 459)
        >>> new_text = sandy.sections.mf8._write_fy(new_from_endf)
        >>> assert new_text == text
        """
        data_endf6 = endf6.data.copy()
        mf = 8
        for (mat, mt, e), data_fy in self.data.groupby(['MAT', 'MT', 'E']):
            sec = endf6.read_section(mat, mf, mt)
            new_data = data_fy.set_index('ZAP')[['FY', 'DFY']].T.to_dict()
            sec['E'][e]['ZAP'] = new_data
            data_endf6[(mat, mf, mt)] = sandy.write_mf8(sec)
        return sandy.Endf6(data_endf6)

    def to_hdf5(self, file, library):
        """
        Write fission yield data to hdf5 file.

        Parameters
        ----------
        `file` : `str`
            HDF5 file
        `library` : `str`
            library name

        Warnings
        --------
        `logging.warning`
            raise a warning if any hdf5 group key is already in used, still
            the existing group will be replaced

        Notes
        -----
        .. note:: the group key for each set of fission yields contains
                    * library: the lowercase name of the library
                    * fy: key "fy"
                    * kind: "independent" or "cumulative"
                    * ZAM: the ZAM number proceeded by prefix "i"
        .. note:: the energy values in the HDF5 file are in MeV
        """
        warnings.filterwarnings("ignore", category=NaturalNameWarning)
        lib = library
        with h5py.File(file, "a") as f:
            for (zam, mt), df in self.data.groupby(["ZAM", "MT"]):
                kind = "independent" if mt == 454 else "cumulative"
                library = lib.lower()
                key = f"{library}/fy/{kind}/{zam}"
                if key in f:
                    msg = f"hdf5 dataset '{key}' already exists and " +\
                           "will be replaced"
                    logging.warning(msg)
                    del f[key]
                else:
                    logging.info(f"creating hdf5 dataset '{key}'")
                group = f.create_group(key)
                group.attrs["nuclide"] = zam  # redunant
                group.attrs["kind"] = kind    # redunant
                group.attrs["library"] = lib  # redunant
                tab = self.energy_table(zam, by="ZAM", kind=kind)
                tab.index *= 1e-6
                tab.to_hdf(file, key, format="fixed")


def _gls_setup(model_sensitivity_object, kind, **kwargs):
    """
    A function to obtain the sensitivity for performing GLS.

    Parameters
    ----------
    model_sensitivity_object : `DecayData` or `Fy`
        Object from which the sensitivity (S) of the model is to be derived:
            y_calc = S*x_prior
    kind : `str`
        Keyword for obtaining sensitivity.
    index: `pd.Index`, optional
        Optional argument with `sandy.Fy` object index to compare with decay
        radioactive decay data index.

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
        S = model_sensitivity_object.get_chain_yield_sensitivity(**kwargs)
    elif kind == 'mass yield':
        S = model_sensitivity_object.get_mass_yield_sensitivity()
    else:
        raise ValueError('Keyword argument "kind" is not valid')
    return S


def fy2hdf(e6file, h5file, lib):
    """
    Write to disk a HDF5 file that reproduces the content of a FY file in
    ENDF6 format.

    Parameters
    ----------
    e6file : `str`
        ENDF-6 filename
    h5file : `str`
        HDF5 filename
    lib : `str`
        library name (it will appear as a hdf5 group)
    """
    # This function os tested in an ALEPH notebook
    endf6 = sandy.Endf6.from_file(e6file)
    logging.info(f"adding FY to '{lib}' in '{h5file}'")
    Fy.from_endf6(endf6, verbose=True).to_hdf5(h5file, lib)
