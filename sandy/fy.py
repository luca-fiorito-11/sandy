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
        Apply a custom perturbation to a given cross section
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

    def custom_perturbation(self, zam, mt, e, zap, pert):
        """
        Apply a custom perturbation in the fission yields identified by the mt,
        for a given zam in a given energy and for a given product.

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
        >>> nfpy_pert = nfpy.custom_perturbation(922350, 459, 0.0253, 551370, -0.1)
        >>> comp = nfpy_pert.data.query('ZAM==922350 & ZAP==551370 & MT==459 & E==0.0253').squeeze().FY
        >>> assert np.setdiff1d(nfpy_pert.data.values, nfpy.data.values) == comp
        """
        df = self.data.copy()
        mask = (df.ZAM == zam) & (df.ZAP == zap) & (df.MT == mt) & (df.E == e)
        df.loc[mask, "FY"] = df[mask].squeeze().FY * (1 + pert)
        return self.__class__(df)

    def apply_bmatrix(self, zam, e, decay_data):
        """
        Perform IFY = (1-B)CFY equation to calculate IFY in a given zam
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

        Returns
        -------
        `Fy`
            Fission yield instance with IFY calculated for a given combination
            of ZAM/e/decay_data.

        Examples
        --------
        >>> zam = [591480, 591481, 601480]
        >>> decay_minimal = sandy.get_endf6_file("jeff_33", 'decay', zam)
        >>> decay_fytest = sandy.DecayData.from_endf6(decay_minimal)
        >>> npfy = Fy(minimal_fytest_2)
        >>> npfy_pert = npfy.apply_bmatrix(942390, 5.00000e+05, decay_fytest)
        >>> diff = npfy_pert.data[npfy_pert.data.FY != npfy.data.FY].FY
        >>> comp = npfy_pert.data.query('ZAM==942390 & MT==454 & E==500e3').squeeze().FY
        >>> assert comp.values.all() == diff.values.all()

        """
        new_data = self.data.copy()
        fy_data = self.filter_by('ZAM', zam)\
            .filter_by('MT', 459)\
            .filter_by("E", e).data\
            .set_index('ZAP')['FY']
        B = decay_data.get_bmatrix().loc[fy_data.index, fy_data.index]
        # Creating (1-B) matrix:
        unit = np.identity(len(B))
        C = unit - B.values
        sensitivity = pd.DataFrame(C, index=B.index, columns=B.columns)
        # Apply (1-B) matrix
        fy_calc = sensitivity.dot(fy_data)
        mask = (new_data.ZAM == zam) & (new_data.MT == 454) & (new_data.E == e)
        new_data.loc[mask, 'FY'] = fy_calc.values
        return self.__class__(new_data)

    def apply_qmatrix(self, zam, energy, decay_data):
        """
        Perform CFY =Q*IFY equation to calculate CFY in a given zam
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

        Returns
        -------
        `Fy`
            Fission yield instance with IFY calculated for a given combination
            of ZAM/e/decay_data.

        Examples
        --------
        >>> zam = [591480, 591481, 601480]
        >>> decay_minimal = sandy.get_endf6_file("jeff_33", 'decay', zam)
        >>> decay_fytest = sandy.DecayData.from_endf6(decay_minimal)
        >>> npfy = Fy(minimal_fytest_2)
        >>> npfy_pert = npfy.apply_bmatrix(942390, 5.00000e+05, decay_fytest)
        >>> diff = npfy_pert.data[npfy_pert.data.FY != npfy.data.FY].FY
        >>> comp = npfy_pert.data.query('ZAM==942390 & MT==459 & E==500e3').squeeze().FY
        >>> assert comp.values.all() == diff.values.all()

        """
        new_data = self.data.copy()
        fy_data = self.filter_by('ZAM', zam)\
            .filter_by('MT', 454)\
            .filter_by("E", energy).data\
            .set_index('ZAP')['FY']
        Q = decay_data.get_qmatrix().loc[fy_data.index, fy_data.index]
        # Apply qmatrix
        fy_calc = Q.dot(fy_data)
        mask = (new_data.ZAM == zam) & (new_data.MT == 459) & (new_data.E == energy)
        new_data.loc[mask, 'FY'] = fy_calc.values
        return self.__class__(new_data)

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
            zam = sec["ZAM"]
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
        >>> fy = sandy.Fy.from_endf6(tape)
        >>> new_tape = fy.to_endf6(tape)
        >>> assert new_tape.data == tape.data
        """
        data_endf6 = endf6.data
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

#    def _to_endf6(self, endf6):
#        # to be written
#        pass

#    def _custom_perturbation(self, pert):
#        # to be written
#        pass


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
