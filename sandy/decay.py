# -*- coding: utf-8 -*-
"""
This module contains all classes and functions dedicated to the processing and
analysis of a decay data.
"""
import logging
import os  # used in docstrings
import pytest  # used in docstrings
import tempfile  # used in docstrings
import yaml  # used in docstrings
import h5py
import copy
from math import sqrt


import numpy as np
import pandas as pd
from pandas.api.types import CategoricalDtype
import scipy
import scipy.sparse as sps
import scipy.sparse.linalg as spsl
from scipy.sparse import csc_matrix
from scipy.sparse.linalg import splu

import sandy

__author__ = "Luca Fiorito"
__all__ = [
        "DecayData",
        "decay_modes",
        "rdd2hdf",
        ]

pd.options.display.float_format = '{:.5e}'.format


decay_modes = {
        0: "gamma",
        1: "beta",
        2: "e.c.",
        3: "i.t.",
        4: "alpha",
        5: "n",
        6: "s.f.",
        7: "p",
        }


class DecayData():
    """
    Container of radioactive nuclide data for several isotopes.

    Attributes
    ----------
    data : `dict`
        source of decay data content

    Methods
    -------
    from_endf6
        extract decay data from ENDF-6 instance
    from_hdf5
        extract decay data from hdf5 file
    get_bmatrix
        extract B-matrix inro dataframe
    get_decay_chains
        extract decay chains into dataframe
    get_qmatrix
        extract Q-matrix into dataframe
    get_transition_matrix
        extract transition matrix into dataframe
    to_hdf5
        write decay data to hdf5 file
    """

    def __repr__(self):
        return self.data.__repr__()

    def __init__(self, dct):
        self.data = dct

    @property
    def data(self):
        """
        Dictionary of RDD content.

        Returns
        -------
        `dict`
            hierarchical RDD content
        """
        return self._data

    @data.setter
    def data(self, data):
        self._data = data

    def get_nuclides(self):
        return sorted(self.data.keys())

    def get_pn(self):
        pn = {}
        for zam, data in self.data.items():
            if data["stable"]:
                continue
            for decay_mode in data["decay_modes"].values():
                # number_del_neuts = f"{rdtp}".count("5")
                daughters = decay_mode["decay_products"]
                if 10 in daughters:
                    pn[zam] = daughters[10]
        series = pd.Series(pn, name="PN")
        series.index.name = "ZAM"
        return series

    def get_halflives(self):
        thalf = {k: v["half_life"] for k, v in self.data.items()}
        series = pd.Series(thalf, name="T1/2")
        series.index.name = "ZAM"
        return series

    def get_decayconstant(self, with_uncertainty=True):
        """
        Extract decay constant and its uncertainty into dataframe.

        Parameters
        ----------
        with_uncertainty : `bool`, optional, default is 'True'
            makes the method return lamdba and its uncertainty if set equal True

        Returns
        -------
        `pandas.DataFrame` or `pandas.Series`
            dataframe with decay constant and associated uncertainty
            series with decay constant if with_uncertainty = False

        Examples
        --------
        >>> endf6 = sandy.get_endf6_file("jeff_33", "decay", 922350)
        >>> rdd = sandy.DecayData.from_endf6(endf6)
        >>> DC = sandy.DecayData.get_decayconstant(rdd)
        >>> print(DC)
                DLAMBDA      LAMBDA
        922350 2.21715e-20 3.12085e-17

        >>> endf6 = sandy.get_endf6_file("jeff_33", "decay", 942390)
        >>> rdd = sandy.DecayData.from_endf6(endf6)
        >>> DC = sandy.DecayData.get_decayconstant(rdd, False)
        >>> print(DC)
        942390   9.10900e-13
        Name: LAMBDA, dtype: float64
        """
        decay_constant = {zam: {
             "LAMBDA": dic['decay_constant'],
             "DLAMBDA": dic['decay_constant_uncertainty'],
             } for zam, dic in self.data.items()}

        df = pd.DataFrame(decay_constant).T
        if with_uncertainty:
            return df
        else:
            return df.LAMBDA

    def custom_perturbation_decayconstant(self, pert, zam):
        """
        Apply a custom perturbation to the decay constant for a given ZAM.

        Parameters
        ----------
        pert : `float`
            Perturbation coefficient as ratio value.
        zam : `int`
            ZAM number of the material to which perturbation is to be
            applied.

        Returns
        -------
        `sandy.decay.DecayData`
            DacayData instance with given ZAM and decay constant perturbed.

        Examples
        --------
        >>> endf6 = sandy.get_endf6_file("jeff_33", 'decay', [922350, 922380])
        >>> rdd = sandy.DecayData.from_endf6(endf6)
        >>> dc_pert = rdd.custom_perturbation_decayconstant(0.1, 922350)
        >>> dc_pert.data
        {922350: {'half_life': 2.22102e+16,
          'decay_constant': 3.432935762018982e-17,
          'decay_constant_uncertainty': 2.2171470275223715e-20,
          'stable': False,
          'spin': 3.5,
          'parity': -1.0,
          'decay_energy': {'beta': 50671.7,
           'gamma': 163616.0,
           'alpha': 4464600.0,
           'total': 4678887.7},
          'decay_energy_uncertainties': {'beta': 4291.63,
           'gamma': 1708.01,
           'alpha': 163255.0,
           'total': 163320.33067324167},
          'decay_modes': {'40x0': {'decay_products': {902310: 1.0, 20040: 1.0},
            'branching_ratio': 1.0,
            'branching_ratio_uncertainty': 0.0001},
           '60x0': {'decay_products': {},
            'branching_ratio': 7.2e-11,
            'branching_ratio_uncertainty': 2.1e-11}}},
         922380: {'half_life': 1.40996e+17,
          'decay_constant': 4.916076913954618e-18,
          'decay_constant_uncertainty': 3.3008662253228094e-21,
          'stable': False,
          'spin': 0.0,
          'parity': 1.0,
          'decay_energy': {'beta': 10208.3,
           'gamma': 1100.17,
           'alpha': 4258770.0,
           'total': 4270078.47},
          'decay_energy_uncertainties': {'beta': 851.158,
           'gamma': 104.433,
           'alpha': 29639.0,
           'total': 29651.4029548764},
          'decay_modes': {'40x0': {'decay_products': {902340: 1.0, 20040: 1.0},
            'branching_ratio': 0.999999,
            'branching_ratio_uncertainty': 1e-08},
           '60x0': {'decay_products': {},
            'branching_ratio': 5.46e-07,
            'branching_ratio_uncertainty': 1e-08}}}}
        """
        pert_dc = copy.deepcopy(self.data)
        pert_dc[zam]['decay_constant'] = pert_dc[zam]['decay_constant'] * (1 + pert)
        return self.__class__(pert_dc)

    def get_decay_chains(self, skip_parents=False, **kwargs):
        """
        Extract decay chains into dataframe.

        Parameters
        ----------
        skip_parent : `bool`, optional, default is `False`

        Returns
        -------
        `pandas.DataFrame`
            decay chains dataframe

        Examples
        --------
        >>> file = os.path.join(sandy.data.__path__[0], "rdd.endf")
        >>> endf6 = sandy.Endf6.from_file(file)
        >>> rdd = sandy.DecayData.from_endf6(endf6)
        >>> rdd.get_decay_chains()
           PARENT  DAUGHTER        YIELD      LAMBDA
        0   10010     10010  0.00000e+00 0.00000e+00
        1  270600    270600 -1.00000e+00 4.16705e-09
        2  270600    280600  1.00000e+00 4.16705e-09
        3  280600    280600  0.00000e+00 0.00000e+00

        >>> rdd.get_decay_chains(skip_parents=True)
           PARENT  DAUGHTER        YIELD      LAMBDA
        0  270600    280600  1.00000e+00 4.16705e-09
        """
        items = []
        columns = ["PARENT", "DAUGHTER", "YIELD", "LAMBDA"]
        for zam, nucl in sorted(self.data.items()):
            yld = 0. if nucl["stable"] else -1.
            if not skip_parents:   # add also the disappearance of the parent
                add = {
                        "PARENT": zam,
                        "DAUGHTER": zam,
                        "YIELD": yld,
                        "LAMBDA": nucl["decay_constant"]
                        }
                items.append(add)
            if nucl["stable"]:
                continue
            for decay_mode in nucl["decay_modes"].values():
                br = decay_mode["branching_ratio"]
                if "decay_products" not in decay_mode:
                    continue  # S.F.
                for zap, yld in decay_mode["decay_products"].items():
                    # add the production of each daughter
                    add = {
                        "PARENT": zam,
                        "DAUGHTER": zap,
                        "YIELD": yld * br,
                        "LAMBDA": nucl["decay_constant"]
                        }
                    items.append(add)
        df = pd.DataFrame(items) \
               .groupby(["PARENT", "DAUGHTER", "LAMBDA"]).sum().reset_index() \
               .sort_values(by=["PARENT", "DAUGHTER"]) \
               .reset_index(drop=True)[columns]
        return df

    def get_chain_yield_sensitivity(self, **kwargs):
        """
        Extract chain fission yield sensitivity matrix.
        - Columns: nucleus represented by the ZAP (`Z*1000 + A*10 + M`).
        - Index: Mass number(A)
        - values: 1 (in the row (A) of that nucleus if it is stable or in the
        mass number of the products in which it decays) or a fraction
        (if that nucleus has more than one path to decay, the fraction
        represent the probability of decaying along that path. As in the
        previous case, the fraction is located in the mass number of the
        final nucleus).

        Parameters
        ----------
        kwargs : `dict`
            keyword arguments for method `get_decay_chains`

        Returns
        -------
        `pandas.DataFrame`
             associated to the given decay chains

        Examples
        --------
        >>> zam = [10010, 10020, 10030, 10040, 10050, 10060, 922350]
        >>> tape = sandy.get_endf6_file("jeff_33",'decay', zam)
        >>> decay_data = DecayData.from_endf6(tape)
        >>> decay_data.get_chain_yield_sensitivity()
        ZAP	      10010	      10020	      10030	      10040	      10050	      10060	     922350
        A							
        1	1.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00
        2	0.00000e+00	1.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00
        3	0.00000e+00	0.00000e+00	1.00000e+00	1.00000e+00	1.00000e+00	5.00000e-01	0.00000e+00
        4	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	1.00000e+00
        5	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	5.00000e-01	0.00000e+00
        231	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	0.00000e+00	1.00000e+00
        """
        chain = self.get_decay_chains().iloc[:, 0:3]
        chain = chain.loc[(chain.DAUGHTER != 10) & (chain.YIELD >= 0)]\
                     .rename(columns={'PARENT': 'ZAP', 'DAUGHTER': 'A'})
        chain.loc[chain.YIELD == 0, 'YIELD'] = 1
        chain['A'] = chain.A.apply(sandy.zam.expand_zam).apply(lambda x: x[1])
        return chain.pivot_table(index='A', columns='ZAP', values='YIELD').fillna(0)

    def get_bmatrix(self, **kwargs):
        """
        Extract B-matrix into dataframe.

        Parameters
        ----------
        kwargs : `dict`
            keyword arguments for method `get_decay_chains`

        Returns
        -------
        `pandas.DataFrame`
            B-matrix associated to the given decay chains

        Examples
        --------
        >>> file = os.path.join(sandy.data.__path__[0], "rdd.endf")
        >>> endf6 = sandy.Endf6.from_file(file)
        >>> rdd = sandy.DecayData.from_endf6(endf6)
        >>> rdd.get_bmatrix()
        PARENT        10010       270600      280600
        DAUGHTER
        10010    0.00000e+00 0.00000e+00 0.00000e+00
        270600   0.00000e+00 0.00000e+00 0.00000e+00
        280600   0.00000e+00 1.00000e+00 0.00000e+00

        >>> tape = sandy.endf6.get_endf6_file("endfb_71", 'decay', 571480)
        >>> decay_data = sandy.DecayData.from_endf6(tape)
        >>> decay_data.get_bmatrix()
        PARENT 	       10 	       571480 	        581470 	       581480
        DAUGHTER
        10 	    0.00000e+00 	1.50000e-03 	0.00000e+00 	0.00000e+00
        571480 	0.00000e+00 	0.00000e+00 	0.00000e+00 	0.00000e+00
        581470 	0.00000e+00 	1.50000e-03 	0.00000e+00 	0.00000e+00
        581480 	0.00000e+00 	9.98500e-01 	0.00000e+00 	0.00000e+00


        >>> h1 = sandy.endf6.get_endf6_file("endfb_71","decay",551480)
        >>> h2 = sandy.endf6.get_endf6_file("endfb_71","decay",551490)
        >>> h3 = h1.merge(h2)
        >>> rdd = sandy.DecayData.from_endf6(h3)
        >>> rdd.get_bmatrix()
        PARENT 	       10 	         551480 	     551490 	     561460 	     561470 	     561480 	     561490
        DAUGHTER
        10 	    0.00000e+00 	2.18793e-01 	6.88450e-01 	0.00000e+00 	0.00000e+00 	0.00000e+00 	0.00000e+00
        551480 	0.00000e+00 	0.00000e+00 	0.00000e+00 	0.00000e+00 	0.00000e+00 	0.00000e+00 	0.00000e+00
        551490 	0.00000e+00 	0.00000e+00 	0.00000e+00 	0.00000e+00 	0.00000e+00 	0.00000e+00 	0.00000e+00
        561460 	0.00000e+00 	1.72560e-04 	0.00000e+00 	0.00000e+00 	0.00000e+00 	0.00000e+00 	0.00000e+00
        561470 	0.00000e+00 	2.18447e-01 	4.09780e-07 	0.00000e+00 	0.00000e+00 	0.00000e+00 	0.00000e+00
        561480 	0.00000e+00 	7.81380e-01 	6.88450e-01 	0.00000e+00 	0.00000e+00 	0.00000e+00 	0.00000e+00
        561490 	0.00000e+00 	0.00000e+00 	3.11550e-01 	0.00000e+00 	0.00000e+00 	0.00000e+00 	0.00000e+00
        """
        B_data = self.get_decay_chains(**kwargs)
        index_sp = CategoricalDtype(sorted(B_data.DAUGHTER.unique()), ordered=True)
        column_sp = CategoricalDtype(sorted(B_data.PARENT.unique()), ordered=True)
        row = B_data.DAUGHTER.astype(index_sp).cat.codes
        col = B_data.PARENT.astype(column_sp).cat.codes
        sparse_matrix = sps.csr_matrix((B_data.YIELD, (row, col)),
                                       shape=(index_sp.categories.size,
                                       column_sp.categories.size))
        B = pd.DataFrame(sparse_matrix.toarray(),
                         index=index_sp.categories.set_names('DAUGHTER'),
                         columns=column_sp.categories)
        B = B.reindex(columns=index_sp.categories.set_names('PARENT'))\
            .fillna(0)
        np.fill_diagonal(B.values, 0)
        return B

    def get_qmatrix(self, keep_neutrons=False, threshold=None, **kwargs):
        """
        Extract Q-matrix dataframe.

        Optional argument
        -------
        Thereshold: 'int'
            Optional argument to avoid numerical fluctuations or
            values so small that they do not have to be taken into
            account.
        Returns
        -------
        `pandas.DataFrame`
            Q-matrix associated to the given decay chains

        >>> file = os.path.join(sandy.data.__path__[0], "rdd.endf")
        >>> endf6 = sandy.Endf6.from_file(file)
        >>> rdd = sandy.DecayData.from_endf6(endf6)
        >>> out = rdd.get_qmatrix()
        >>> comp = pd.DataFrame([[1, 0, 0],
        ...                      [0, 1, 0],
        ...                      [0, 1, 1]],
        ...                     dtype=float,
        ...                     index=[10010, 270600, 280600],
        ...                     columns=[10010, 270600, 280600])
        >>> comp.index.name = "DAUGHTER"
        >>> comp.columns.name = "PARENT"
        >>> pd.testing.assert_frame_equal(comp, out)


        >>> h1 = sandy.endf6.get_endf6_file("endfb_71","decay",551480)
        >>> h2 = sandy.endf6.get_endf6_file("endfb_71","decay",551490)
        >>> h3 = h1.merge(h2)
        >>> rdd = sandy.DecayData.from_endf6(h3)
        >>> rdd.get_qmatrix()
        PARENT 	     551480 	     551490 	     561460 	     561470 	     561480 	     561490
        DAUGHTER
        551480 	1.00000e+00 	0.00000e+00 	0.00000e+00 	0.00000e+00 	0.00000e+00 	0.00000e+00
        551490 	0.00000e+00 	1.00000e+00 	0.00000e+00 	0.00000e+00 	0.00000e+00 	0.00000e+00
        561460 	1.72560e-04 	0.00000e+00 	1.00000e+00 	0.00000e+00 	0.00000e+00 	0.00000e+00
        561470 	2.18447e-01 	4.09780e-07 	0.00000e+00 	1.00000e+00 	0.00000e+00 	0.00000e+00
        561480 	7.81380e-01 	6.88450e-01 	0.00000e+00 	0.00000e+00 	1.00000e+00 	0.00000e+00
        561490 	0.00000e+00 	3.11550e-01 	0.00000e+00 	0.00000e+00 	0.00000e+00 	1.00000e+00
        """
        B = self.get_bmatrix(**kwargs)
        if not keep_neutrons:
            if 10 in B.index:
                B.drop(index=10, inplace=True)
            if 10 in B.columns:
                B.drop(columns=10, inplace=True)
        index = B.index
        columns = B.columns
        B = sps.csc_matrix(B)
        unit = sps.csc_matrix(sps.identity(B.shape[0]))
        C = spsl.splu(sps.csc_matrix(unit-B))
        qmatrix = pd.DataFrame(C.solve(unit.toarray()), index=index,
                               columns=columns)
        if threshold is not None:
            qmatrix[qmatrix < threshold] = 0
        return qmatrix

    def get_transition_matrix(self):
        """
        Extract transition matrix into dataframe.

        Returns
        -------
        `pandas.DataFrame`
            transition matrix associated to the given decay chains

        Examples
        --------
        >>> file = os.path.join(sandy.data.__path__[0], "rdd.endf")
        >>> endf6 = sandy.Endf6.from_file(file)
        >>> rdd = sandy.DecayData.from_endf6(endf6)
        >>> rdd.get_transition_matrix()
        PARENT        10010        270600      280600
        DAUGHTER
        10010    0.00000e+00  0.00000e+00 0.00000e+00
        270600   0.00000e+00 -4.16705e-09 0.00000e+00
        280600   0.00000e+00  4.16705e-09 0.00000e+00
        """
        df = self.get_decay_chains()
        df["YIELD"] *= df["LAMBDA"]
        T = df.pivot_table(
                index="DAUGHTER",
                columns="PARENT",
                values="YIELD",
                aggfunc=np.sum,
                )\
              .astype(float)\
              .fillna(0)
        return T.reindex(T.columns.values, fill_value=0.0)

    @classmethod
    def from_endf6(cls, endf6, verbose=False):
        """
        Extract hierarchical structure of decay data from `sandy.Endf6`
        instance.

        Parameters
        ----------
        tape : `sandy.Endf6`
            instance containing decay data
        verbose : `bool`, optional, default is `False`
            flag to print information when reading ENDF-6 file

        Returns
        -------
        `dict`
            structured container with RDD.

        Raises
        ------
        `sandy.Error`
            if no decay data is found

        Examples
        --------
        Load test ENDF-6 file with data for H1 and Co60.
        >>> file = os.path.join(sandy.data.__path__[0], "rdd.endf")
        >>> endf6 = sandy.Endf6.from_file(file)
        >>> rdd = sandy.DecayData.from_endf6(endf6)
        >>> print(yaml.dump(rdd))
        !!python/object:sandy.decay.DecayData
        _data:
          10010:
            decay_constant: 0
            decay_constant_uncertainty: 0
            decay_energy:
              alpha: 0.0
              beta: 0.0
              gamma: 0.0
              total: 0.0
            decay_energy_uncertainties:
              alpha: 0.0
              beta: 0.0
              gamma: 0.0
              total: 0.0
            half_life: 0.0
            parity: 1.0
            spin: 0.5
            stable: true
          270600:
            decay_constant: 4.167050502344267e-09
            decay_constant_uncertainty: 6.324352137605637e-13
            decay_energy:
              alpha: 0.0
              beta: 96522.0
              gamma: 2503840.0
              total: 2600362.0
            decay_energy_uncertainties:
              alpha: 0.0
              beta: 202.529
              gamma: 352.186
              total: 406.26712202318316
            decay_modes:
              10x0:
                branching_ratio: 1.0
                branching_ratio_uncertainty: 0.0
                decay_products:
                  280600: 1.0
            half_life: 166340000.0
            parity: 1.0
            spin: 5.0
            stable: false
          280600:
            decay_constant: 0
            decay_constant_uncertainty: 0
            decay_energy:
              alpha: 0.0
              beta: 0.0
              gamma: 0.0
              total: 0.0
            decay_energy_uncertainties:
              alpha: 0.0
              beta: 0.0
              gamma: 0.0
              total: 0.0
            half_life: 0.0
            parity: 1.0
            spin: 0.0
            stable: true
        <BLANKLINE>
        """
        tape = endf6.filter_by(listmf=[8], listmt=[457])
        if tape.is_empty:
            raise sandy.Error("no decay data found in file")
        groups = {}
        for mat, mf, mt in tape.keys:
            sec = endf6.read_section(mat, mf, mt)
            zam = int(sec["ZA"]*10 + sec["LISO"])
            if verbose:
                logging.info(f"reading 'ZAM={zam}'...")
            groups[zam] = {
                    "half_life": sec["HL"],
                    "decay_constant": sec["LAMBDA"],
                    "decay_constant_uncertainty": sec["DLAMBDA"],
                    "stable": bool(sec["NST"]),
                    "spin": sec["SPI"],
                    "parity": sec["PAR"],
                    "decay_energy": {
                            "beta": sec["E"][0],
                            "gamma": sec["E"][1],
                            "alpha": sec["E"][2],
                            "total": sum(sec["E"][:3]),
                            },
                    "decay_energy_uncertainties": {
                            "beta": sec["DE"][0],
                            "gamma": sec["DE"][1],
                            "alpha": sec["DE"][2],
                            "total": sqrt(
                                sum(x**2 for x in sec["DE"][:3])
                                ),
                            },
                    }
            if groups[zam]["stable"]:
                assert groups[zam]["decay_constant"] == 0
                assert "DK" not in sec
                continue
            groups[zam]["decay_modes"] = {}
            for key, dk in sec["DK"].items():
                rtyp = key.split("x")[0]
                residual_state = dk["RFS"]
                decay_mode_data = {
                        "decay_products": get_decay_products(
                                                rtyp,
                                                zam,
                                                residual_state,
                                                ),
                        "branching_ratio": dk["BR"],
                        "branching_ratio_uncertainty": dk["DBR"],
                        }
                groups[zam]["decay_modes"][key] = decay_mode_data
        return cls(groups)

    @classmethod
    def from_hdf5(cls, filename, lib, zam=None):
        """
        Extract hierarchical structure of decay data from hdf5 file.

        Parameters
        ----------
        filename : `str`
            hdf5 filename (absolute or relative)
        lib : `str`
            library ID contained in the hdf5 file
        zam : `int`, optional, default is `None`
            optional selection of individual nuclide (avoid loading all
            library)

        Returns
        -------
        `DecayData`
            decay data object

        Examples
        --------
        Examples are in method `to_hdf5`
        """
        with h5py.File(filename, 'r') as h5file:
            # the last slash is important
            group = f"{lib}/rdd/{zam}/" if zam else f"{lib}/rdd/"
            data = sandy.tools.recursively_load_dict_contents_from_group(
                    h5file,
                    group,
                    )
        if zam:
            data = {zam: data}
        return cls(data)

    def to_hdf5(self, filename, lib, mode="a"):
        """
        Dump decay data to hdf5 file.

        Parameters
        ----------
        filename : `str`
            name of the hdf5 file (with absolute or relative path)
        lib : `str`
            name of the library that will be used

        Notes
        -----
        .. note:: decay data are saved in groups with key `'{lib}/rdd/{zam}'`,
                  where `'{lib}'` and `'{zam}'` are the library and ZAM
                  identifiers.
                  The rest of the contents is the structured following the
                  nested dictionaries.

        Examples
        --------
        Write file into hdf5 format
        >>> file = os.path.join(sandy.data.__path__[0], "rdd.endf")
        >>> endf6 = sandy.Endf6.from_file(file)
        >>> rdd = sandy.DecayData.from_endf6(endf6)
        >>> f = tempfile.TemporaryDirectory()
        >>> path = os.path.join(f.name, "test.h5")
        >>> rdd.to_hdf5(path, "jeff_33")

        Then, make sure `sandy` reads it correctly
        >>> rdd.from_hdf5(path, "jeff_33")
        {10010: {'decay_constant': 0, 'decay_constant_uncertainty': 0, 'decay_energy': {'alpha': 0.0, 'beta': 0.0, 'gamma': 0.0, 'total': 0.0}, 'decay_energy_uncertainties': {'alpha': 0.0, 'beta': 0.0, 'gamma': 0.0, 'total': 0.0}, 'half_life': 0.0, 'parity': 1.0, 'spin': 0.5, 'stable': True}, 270600: {'decay_constant': 4.167050502344267e-09, 'decay_constant_uncertainty': 6.324352137605637e-13, 'decay_energy': {'alpha': 0.0, 'beta': 96522.0, 'gamma': 2503840.0, 'total': 2600362.0}, 'decay_energy_uncertainties': {'alpha': 0.0, 'beta': 202.529, 'gamma': 352.186, 'total': 406.26712202318316}, 'decay_modes': {'10x0': {'branching_ratio': 1.0, 'branching_ratio_uncertainty': 0.0, 'decay_products': {280600: 1.0}}}, 'half_life': 166340000.0, 'parity': 1.0, 'spin': 5.0, 'stable': False}, 280600: {'decay_constant': 0, 'decay_constant_uncertainty': 0, 'decay_energy': {'alpha': 0.0, 'beta': 0.0, 'gamma': 0.0, 'total': 0.0}, 'decay_energy_uncertainties': {'alpha': 0.0, 'beta': 0.0, 'gamma': 0.0, 'total': 0.0}, 'half_life': 0.0, 'parity': 1.0, 'spin': 0.0, 'stable': True}}
        >>> f.cleanup()
        """
        with h5py.File(filename, mode=mode) as h5file:
            for nucl, data in self.data.items():
                group = f"{lib}/rdd/{nucl:d}/"  # the last slash is important
                logging.info(f"dumping RDD for ZAM={nucl} into '{group}'")
                sandy.tools.recursively_save_dict_contents_to_group(
                        h5file,
                        group,
                        data,
                        )


def expand_decay_type(zam, dectyp):
    """
    Given a nuclide and an individual decay mode as in `decay_modes`,
    return:
        - the decay product
        - the number of emitted neutrons
        - the number of emitted protons
        - the number of emitted alphas

    Parameters
    ----------
    zam : `int`
        ZAM identifier
    dectyp : `int`
        decay mode

    Returns
    -------
    `int`
        decay daughter product ZAM identifier
    `float`
        number of emitted neutrons
    `float`
        number of emitted protons
    `float`
        number of emitted alphas

    Notes
    -----
    ..note :: it is assumed that only one nuclide is produced plus neutrons,
              protons and/or alpha particle.
              Other particles such as photons or betas are not considered.
    ..note :: decay modes are taken from the ENDF-6 format manual

    Examples
    --------
    Expand gamma decay (#0)
    >>> d, n, p, a = sandy.decay.expand_decay_type(581480, 0)
    >>> assert d == 581480
    >>> assert n == 0
    >>> assert p == 0
    >>> assert a == 0

    Expand beta decay (#1)
    >>> d, n, p, a = sandy.decay.expand_decay_type(581480, 1)
    >>> assert d == 591480
    >>> assert n == 0
    >>> assert p == 0
    >>> assert a == 0

    Expand electron capture (#2)
    >>> d, n, p, a = sandy.decay.expand_decay_type(581480, 2)
    >>> assert d == 571480
    >>> assert n == 0
    >>> assert p == 0
    >>> assert a == 0

    Expand isomeric transition (#3)
    >>> d, n, p, a = sandy.decay.expand_decay_type(581480, 3)
    >>> assert d == 581480
    >>> assert n == 0
    >>> assert p == 0
    >>> assert a == 0

    Expand alpha decay (#4)
    >>> d, n, p, a = sandy.decay.expand_decay_type(581480, 4)
    >>> assert d == 561440
    >>> assert n == 0
    >>> assert p == 0
    >>> assert a == 1

    Expand neutron decay (#5)
    d, n, p, a = sandy.decay.expand_decay_type(581480, 5)
    assert d == 581470
    assert n == 1
    assert p == 0
    assert a == 0

    Expand spontaneous fission(#6)
    >>> d, n, p, a = sandy.decay.expand_decay_type(581480, 6)
    >>> assert d == 581480
    >>> assert n == 0
    >>> assert p == 0
    >>> assert a == 0

    Expand proton decay (#7):
    >>> d, n, p, a = sandy.decay.expand_decay_type(581480, 7)
    >>> assert d == 571470
    >>> assert n == 0
    >>> assert p == 1
    >>> assert a == 0

    Expand unknown decay:
    >>> with pytest.raises(ValueError):
    ...     sandy.decay.expand_decay_type(581480, 8)
    """
    daughter = zam//10
    neutrons = 0.
    protons = 0.
    alphas = 0.
    if dectyp == 1:  # Beta decay
        daughter += 1001 - 1
    elif dectyp == 2:  # Electron capture and/or positron emission
        daughter += 1 - 1001
    elif dectyp == 3:  # Isomeric transition
        pass
    elif dectyp == 4:  # Alpha decay
        daughter -= 2004
        alphas += 1.
    elif dectyp == 5:  # Neutron emission
        daughter -= 1
        neutrons += 1.
    elif dectyp == 6:  # Spontaneous fission
        pass
    elif dectyp == 7:  # Proton emission
        daughter -= 1001
        protons += 1.
    elif dectyp == 0:  # Gamma emission (not used in MT457)
        pass
    else:  # Unknown decay mode
        raise ValueError(f"unknown decay mode {dectyp} for ZAM={zam}")
    return daughter*10, neutrons, protons, alphas


def get_decay_products(rtyp, zam, meta=0, br=1.):
    """
    For a given parent nuclide and decay mode (individual or composed),
    extract a dictionary of decay products.

    Parameters
    ----------
    rtyp : `str`
        string of decay modes where:
            1. Beta decay
            2. Electron capture and/or positron emission
            3. Isomeric transition
            4. Alpha decay
            5. Neutron emission (not delayed neutron decay)
            6. Spontaneous fission
            7. Proton emission

        Decay mode combinations are allowed, e.g. "15" means Beta decay
        followed by neutron emission (delayed neutron decay).
    zam : `int`
        ZAM identifier of the nuclide undergoing decay (parent)
    meta : `int`, optional, default is `0`
        Isomeric state flag for daughter nuclide, e.g. `meta=0` is ground
        state, `meta=1` is first isomeric state, etc.
    br : `float`
        branching ratio

    Returns
    -------
    `dict`
        dictionary of decay products where the keys are the ZAM identifiers
        for the products and the values are the corresponding yield.

    Examples
    --------
    Extract products of fake decay process including all available decay modes.
    >>> sandy.decay.get_decay_products("01234567", 581480)
    {551420: 1.0, 10: 1.0, 10010: 1.0, 20040: 1.0}

    ...change the metastate of the product
    >>> sandy.decay.get_decay_products("01234567", 581480, meta=1)
    {551421: 1.0, 10: 1.0, 10010: 1.0, 20040: 1.0}

    ...and then use a different braanching ratio
    >>> sandy.decay.get_decay_products("01234567", 581480, br=0.1)
    {551420: 0.1, 10: 0.1, 10010: 0.1, 20040: 0.1}
    """
    daughter = zam + 0
    neutrons = 0.
    protons = 0.
    alphas = 0.
    for dectyp in map(int, rtyp):
        daughter, n, h, a = expand_decay_type(daughter, dectyp)
        neutrons += n
        protons += h
        alphas += a
    daughter = int(daughter + meta)
    products = {}
    if daughter != zam:
        products[daughter] = 1.0 * br
    if neutrons != 0:
        products[10] = neutrons * br
    if protons != 0:
        products[10010] = protons * br
    if alphas != 0:
        products[20040] = alphas * br
    return products


def rdd2hdf(e6file, h5file, lib):
    """
    Write to disk a HDF5 file that reproduces the content of a RDD file in
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
    endf6 = sandy.Endf6.from_file(e6file)
    logging.info(f"adding RDD to '{lib}' in '{h5file}'")
    DecayData.from_endf6(endf6, verbose=True).to_hdf5(h5file, lib)
