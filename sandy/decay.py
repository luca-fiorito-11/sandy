# -*- coding: utf-8 -*-
"""
This module contains all classes and functions dedicated to the processing and
analysis of a decay data.
"""
import logging
import os  # used in docstrings
import tempfile  # used in docstrings
import copy
from math import sqrt


import numpy as np
import pandas as pd
import scipy
from scipy.sparse import csc_matrix
from scipy.sparse.linalg import splu

import sandy

__author__ = "Luca Fiorito"
__all__ = [
        "DecayData",
        "decay_modes",
        "BranchingRatio",
        "HalfLife",
        "DecayEnergy",
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
        """
        Extract probability of neutron emission.

        Returns
        -------
        `pandas.Series`
            panda series with ZAM index and probability of neutrom emission

        Examples
        --------
        >>> endf6 = sandy.get_endf6_file("jeff_33", "decay", 391000)
        >>> rdd = sandy.DecayData.from_endf6(endf6)
        >>> rdd.get_pn()
        ZAM
        391000   1.00000e+00
        Name: PN, dtype: float64
        """
        pn = {}
        for zam, data in self.data.items():
            if data["stable"]:
                continue
            for (rtyp, rfs), decay_mode in data["decay_modes"].items():
                # number_del_neuts = f"{rdtp}".count("5")
                daughters = decay_mode["decay_products"]
                if 10 in daughters:
                    pn[zam] = daughters[10]
        series = pd.Series(pn, name="PN")
        series.index.name = "ZAM"
        return series

    def get_half_life(self, with_uncertainty=True):
        """
        Extract half life and its uncertainty.

        Parameters
        ----------
        with_uncertainty : `bool`, optional, default is 'True'
            makes the method return half lives and uncertainties
            if set equal True, or else return only the half lives

        Returns
        -------
        `sandy.HalfLife`
            object containing half life and associated uncertainty or
            only half life if with_uncertainty=False

        Notes
        -----
        .. note:: if a nuclide is stable, half-life of zero will be assigned,
            according with the value stored in the ENDF6 format.

        Examples
        --------
        >>> endf6 = sandy.get_endf6_file("jeff_33", "decay", [942400, 922350])
        >>> rdd = sandy.DecayData.from_endf6(endf6)
        >>> rdd.get_half_life()
                        HL         DHL
        ZAM                           
        922350 2.22102e+16 1.57788e+13
        942400 2.07108e+11 1.57785e+08

        >>> rdd.get_half_life(with_uncertainty=False)
                        HL
        ZAM               
        922350 2.22102e+16
        942400 2.07108e+11
        
        Stable nuclide:
        >>> endf6 = sandy.get_endf6_file("jeff_33", "decay", 260560)
        >>> rdd = sandy.DecayData.from_endf6(endf6)
        >>> rdd.get_half_life(with_uncertainty=False)
                        HL
        ZAM               
        260560 0.00000e+00
        """
        thalf = {zam: {
             "HL": dic['half_life'],
             "DHL": dic['half_life_uncertainty'],
             } for zam, dic in self.data.items()}
        df = pd.DataFrame(thalf).T
        df.index.name = "ZAM"
        if with_uncertainty:
            return HalfLife(df)
        else:
            return HalfLife(df.HL)


    def get_branching_ratio(self, with_uncertainty=True):
        """
        Extract branching ratios and their uncertainties.

        Parameters
        ----------
        with_uncertainty : `bool`, optional, default is 'True'
            makes the method return branching ratios and uncertainties
            if set equal True, or else return only the branching ratios

        Returns
        -------
        `sandy.BranchingRatio`
            object containing branching ratios and associated uncertainties or
            only branching ratios if with_uncertainty=False

        Examples
        --------
        >>> endf6 = sandy.get_endf6_file("jeff_33", "decay", [942410, 922350])
        >>> rdd = sandy.DecayData.from_endf6(endf6)
        >>> rdd.get_branching_ratio()
                                 BR         DBR
        ZAM    RTYP RFS                        
        922350 4    0   1.00000e+00 1.00000e-04
               6    0   7.20000e-11 2.10000e-11
        942410 4    0   2.44000e-05 0.00000e+00
               1    0   9.99976e-01 0.00000e+00

        >>> rdd.get_branching_ratio(with_uncertainty=False)
                                 BR
        ZAM    RTYP RFS            
        922350 4    0   1.00000e+00
               6    0   7.20000e-11
        942410 4    0   2.44000e-05
               1    0   9.99976e-01

        >>> endf6 = sandy.get_endf6_file("jeff_33", "decay", [942410, 10010, 922350])
        >>> rdd = sandy.DecayData.from_endf6(endf6)
        >>> rdd.get_branching_ratio(with_uncertainty=False)
                                 BR
        ZAM    RTYP RFS            
        922350 4    0   1.00000e+00
               6    0   7.20000e-11
        942410 4    0   2.44000e-05
               1    0   9.99976e-01
               
        Decay at first isomeric state:
        >>> endf6 = sandy.get_endf6_file("jeff_33", "decay", 942390)
        >>> rdd = sandy.DecayData.from_endf6(endf6)
        >>> rdd.get_branching_ratio(with_uncertainty=False)
                                 BR
        ZAM    RTYP RFS            
        942390 4    0   6.00000e-04
                    1   9.99400e-01
               6    0   3.10000e-12

        Stable nuclide:
        >>> endf6 = sandy.get_endf6_file("jeff_33", "decay", 260560)
        >>> rdd = sandy.DecayData.from_endf6(endf6)
        >>> rdd.get_branching_ratio()
        Empty DataFrame
        Columns: [BR, DBR]
        Index: []
        """
        br = []
        zam = []
        rtyp_ = []
        rfs_ = []
        for z, dic in self.data.items():
            if 'decay_modes' in dic.keys():
               for (rtyp, rfs), dk in dic['decay_modes'].items():
                    br.append([
                        dk['branching_ratio'],
                        dk['branching_ratio_uncertainty'],
                        ])
                    rtyp_.append(rtyp)
                    rfs_.append(rfs)
                    zam.append(z)
        tuples = zip(* [zam,
                        rtyp_,
                        rfs_])
        idx = pd.MultiIndex.from_tuples(tuples, names=['ZAM', 'RTYP', 'RFS'])
        df = pd.DataFrame(br, index=idx, columns=['BR', 'DBR'])
        if with_uncertainty:
            return BranchingRatio(df)
        else:
            return BranchingRatio(df.BR)

    def get_decay_energy(self, with_uncertainty=True):
        """
        Extract decay energy and its uncertainty.

        Parameters
        ----------
        with_uncertainty : `bool`, optional, default is 'True'
            makes the method return decay energies and uncertainties
            if set equal True, or else return only the decay energies

        Returns
        -------
        `sandy.DecayEnergy`
            object containing decay energy and associated uncertainty or
            only decay energy if with_uncertainty=False

        Examples
        --------
        >>> endf6 = sandy.get_endf6_file("jeff_33", "decay", [942400, 922350])
        >>> rdd = sandy.DecayData.from_endf6(endf6)
        >>> rdd.get_decay_energy()
                               E          DE
        ZAM    TYPE                         
        922350 alpha 4.46460e+06 1.63255e+05
               beta  5.06717e+04 4.29163e+03
               gamma 1.63616e+05 1.70801e+03
        942400 alpha 5.24303e+06 3.63881e+04
               beta  1.11164e+04 9.02572e+02
               gamma 1.36292e+03 1.33403e+02

        >>> rdd.get_decay_energy(with_uncertainty=False)
                               E
        ZAM    TYPE             
        922350 alpha 4.46460e+06
               beta  5.06717e+04
               gamma 1.63616e+05
        942400 alpha 5.24303e+06
               beta  1.11164e+04
               gamma 1.36292e+03

        Stable nuclide:
        >>> endf6 = sandy.get_endf6_file("jeff_33", "decay", 260560)
        >>> rdd = sandy.DecayData.from_endf6(endf6)
        >>> rdd.get_decay_energy(with_uncertainty=False)
                               E
        ZAM    TYPE             
        260560 alpha 0.00000e+00
               beta  0.00000e+00
               gamma 0.00000e+00
        """
        decay_energy = []
        decay_energy_uncertainty = []
        zam = []
        for z, dic in self.data.items():
            decay_energy.extend([
                dic['decay_energy']['alpha'],
                dic['decay_energy']['beta'],
                dic['decay_energy']['gamma'],
                ])
            decay_energy_uncertainty.extend([
                dic['decay_energy_uncertainties']['alpha'],
                dic['decay_energy_uncertainties']['beta'],
                dic['decay_energy_uncertainties']['gamma'],
                ])
            zam.append(z)
        name = ['alpha', 'beta', 'gamma']
        df = pd.DataFrame(zip(decay_energy, decay_energy_uncertainty),
                          index=pd.MultiIndex.from_product([zam, name], names=['ZAM', 'TYPE']),
                          columns=['E', 'DE'])
        if with_uncertainty:
            return DecayEnergy(df)
        else:
            return DecayEnergy(df.E)

    def get_decay_chains(self, skip_parents=False, cut_hl=False, **kwargs):
        """
        Extract decay chains into dataframe.

        Parameters
        ----------
        skip_parents : `bool`, optional, default is `False`
            flag to skip the parent information
        cut_hl: `bool`, optional, default is `False`
            cut all the dacay modes of the nuclides with an half life larger
            than 100 years

        Returns
        -------
        `pandas.DataFrame`
            decay chains dataframe

        Examples
        --------
        >>> endf6 = sandy.get_endf6_file("jeff_33", 'decay', [10010, 270600, 280600])
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

        Cut the dacay modes of the nuclides with an half life larger
        than 100 years:
        >>> tape = sandy.get_endf6_file("jeff_33", "decay", 601440)
        >>> rdd = sandy.DecayData.from_endf6(tape)
        >>> rdd.get_decay_chains(cut_hl=False)
           PARENT  DAUGHTER        YIELD      LAMBDA
        0  601440     20040  1.00000e+00 9.59169e-24
        1  601440    581400  1.00000e+00 9.59169e-24
        2  601440    601440 -1.00000e+00 9.59169e-24
        
        >>> rdd.get_decay_chains(cut_hl=True)
           PARENT  DAUGHTER        YIELD      LAMBDA
        0  601440     20040  0.00000e+00 9.59169e-24
        1  601440    581400  0.00000e+00 9.59169e-24
        2  601440    601440 -1.00000e+00 9.59169e-24
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
            for (rtyp, rfs), decay_mode in nucl["decay_modes"].items():
                br = decay_mode["branching_ratio"]
                if cut_hl and nucl['half_life'] > 3153600000: # 100 years in seconds
                    br = 0
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
        return chain.pivot_table(index='A', columns='ZAP', values='YIELD',
                                 aggfunc='sum', fill_value=0).astype(float).fillna(0)

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
        >>> endf6 = sandy.get_endf6_file("jeff_33", 'decay', [10010, 270600, 280600])
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

        >>> h1 = sandy.endf6.get_endf6_file("endfb_71", "decay", 551480)
        >>> h2 = sandy.endf6.get_endf6_file("endfb_71", "decay", 551490)
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
        B = self.get_decay_chains(**kwargs) \
                .pivot_table(
                        index="DAUGHTER",
                        columns="PARENT",
                        values="YIELD",
                        aggfunc="sum",
                        fill_value=0.0,
                        )\
                .astype(float)\
                .fillna(0)
        B_reindex = B.reindex(B.index.values, fill_value=0.0, axis=1)
        vals = B_reindex.values
        np.fill_diagonal(vals, 0)
        return pd.DataFrame(vals, index=B_reindex.index, columns=B_reindex.columns)

    def get_qmatrix(self, keep_neutrons=False, threshold=None, **kwargs):
        """
        Extract Q-matrix dataframe.

        Optional argument
        -------
        kwargs : `dict`
            keyword arguments for method `get_decay_chains`
        thereshold: `int`, optional, default is `None`
            argument to avoid numerical fluctuations or
            values so small that they do not have to be taken into
            account
        keep_neutrons : `bool`, optional, default is `False`
            flag to skip the column with neutron data

        Returns
        -------
        `pandas.DataFrame`
            Q-matrix associated to the given decay chains

        Examples
        --------
        >>> endf6 = sandy.get_endf6_file("jeff_33", 'decay', [10010, 270600, 280600])
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

        >>> h1 = sandy.get_endf6_file("endfb_71", "decay", 551480)
        >>> h2 = sandy.get_endf6_file("endfb_71", "decay", 551490)
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

        Cut the dacay modes of the nuclides with an half life larger
        than 100 years:
        >>> tape = sandy.get_endf6_file("jeff_33", "decay", 601440)
        >>> rdd = sandy.DecayData.from_endf6(tape)
        >>> rdd.get_qmatrix(cut_hl=False)
        PARENT        20040       581400      601440
        DAUGHTER                                    
        20040    1.00000e+00 0.00000e+00 1.00000e+00
        581400   0.00000e+00 1.00000e+00 1.00000e+00
        601440   0.00000e+00 0.00000e+00 1.00000e+00

        >>> rdd.get_qmatrix(cut_hl=True)
        PARENT        20040       581400      601440
        DAUGHTER                                    
        20040    1.00000e+00 0.00000e+00 0.00000e+00
        581400   0.00000e+00 1.00000e+00 0.00000e+00
        601440   0.00000e+00 0.00000e+00 1.00000e+00

        Skip column with neutron information:
        >>> tape = sandy.endf6.get_endf6_file("endfb_71", 'decay', 571480)
        >>> rdd = sandy.DecayData.from_endf6(tape)
        >>> rdd.get_qmatrix(keep_neutrons=False)
        PARENT        571480      581470      581480
        DAUGHTER                                    
        571480   1.00000e+00 0.00000e+00 0.00000e+00
        581470   1.50000e-03 1.00000e+00 0.00000e+00
        581480   9.98500e-01 0.00000e+00 1.00000e+00

        >>> rdd.get_qmatrix(keep_neutrons=True)
        PARENT        10          571480      581470      581480
        DAUGHTER                                                
        10       1.00000e+00 1.50000e-03 0.00000e+00 0.00000e+00
        571480   0.00000e+00 1.00000e+00 0.00000e+00 0.00000e+00
        581470   0.00000e+00 1.50000e-03 1.00000e+00 0.00000e+00
        581480   0.00000e+00 9.98500e-01 0.00000e+00 1.00000e+00
        """
        B = self.get_bmatrix(**kwargs)
        if not keep_neutrons:
            if 10 in B.index:
                B.drop(index=10, inplace=True)
            if 10 in B.columns:
                B.drop(columns=10, inplace=True)
        unit = np.identity(len(B))
        C = unit - B.values
        C_inv = splu(csc_matrix(C))
        qmatrix = pd.DataFrame(
            C_inv.solve(unit),
            index=B.index,
            columns=B.columns,
            )
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
        >>> endf6 = sandy.get_endf6_file("jeff_33", 'decay', [10010, 270600, 280600])
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
        >>> endf6 = sandy.get_endf6_file("jeff_33", 'decay', [10010, 270600, 280600])
        >>> rdd = sandy.DecayData.from_endf6(endf6)
        >>> import yaml
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
            decay_energy_uncertainties:
              alpha: 0.0
              beta: 0.0
              gamma: 0.0
            half_life: 0.0
            half_life_uncertainty: 0.0
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
            decay_energy_uncertainties:
              alpha: 0.0
              beta: 202.529
              gamma: 352.186
            decay_modes:
              ? !!python/tuple
              - 1
              - 0
              : branching_ratio: 1.0
                branching_ratio_uncertainty: 0.0
                decay_products:
                  280600: 1.0
            half_life: 166340000.0
            half_life_uncertainty: 25245.5
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
            decay_energy_uncertainties:
              alpha: 0.0
              beta: 0.0
              gamma: 0.0
            half_life: 0.0
            half_life_uncertainty: 0.0
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
                    "half_life_uncertainty": sec["DHL"],
                    "decay_constant": sec["LAMBDA"],
                    "decay_constant_uncertainty": sec["DLAMBDA"],
                    "stable": bool(sec["NST"]),
                    "spin": sec["SPI"],
                    "parity": sec["PAR"],
                    "decay_energy": {
                            "beta": sec["E"][0],
                            "gamma": sec["E"][1],
                            "alpha": sec["E"][2],
                            },
                    "decay_energy_uncertainties": {
                            "beta": sec["DE"][0],
                            "gamma": sec["DE"][1],
                            "alpha": sec["DE"][2],
                            },
                    }
            if groups[zam]["stable"]:
                assert groups[zam]["decay_constant"] == 0
                assert "DK" not in sec
                continue
            groups[zam]["decay_modes"] = {}
            for dk in sec["DK"]:
                rtyp = dk['RTYP']
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
                groups[zam]["decay_modes"][(rtyp, residual_state)] = decay_mode_data
        return cls(groups)

    def to_endf6(self, endf6):
        """
        Update decay data in `Endf6` instance with those available in a
        `DecayData` instance.

        Parameters
        ----------
        `endf6` : `sandy.Endf6`
            `Endf6` instance

        Returns
        -------
        `sandy.Endf6`
            `Endf6` instance with updated decay data

        Examples
        --------
        >>> tape = sandy.get_endf6_file("jeff_33", "decay", 922350)
        >>> rdd = sandy.DecayData.from_endf6(tape)
        >>> new_tape = rdd.to_endf6(tape)
        >>> new_tape
        MAT   MF  MT
        3542  1   451     9.223500+4 2.330250+2         -1          1  ...
                  452     9.223500+4 2.330250+2          0          1  ...
              8   457     92235.0000 233.025000          0          0  ...
        dtype: object
        """
        data = endf6.data.copy()
        tape = endf6.filter_by(listmf=[8], listmt=[457])
        for (mat, mf, mt) in tape.keys:
            sec = tape.read_section(mat, mf, mt)
            zam = int(sec["ZA"] * 10 + sec["LISO"])
            sec["HL"] = self.data[zam]['half_life']
            sec["LAMBDA"] = self.data[zam]['decay_constant']
            sec["DLAMBDA"] = self.data[zam]['decay_constant_uncertainty']
            sec["NST"] = int(self.data[zam]['stable'])
            sec["SPI"] = self.data[zam]['spin']
            sec["PAR"] = self.data[zam]['parity']
            sec["E"][0] = self.data[zam]['decay_energy']['beta']
            sec["E"][1] = self.data[zam]['decay_energy']['gamma']
            sec["E"][2] = self.data[zam]['decay_energy']['alpha']
            sec["DE"][0] = self.data[zam]['decay_energy_uncertainties']['beta']
            sec["DE"][1] = self.data[zam]['decay_energy_uncertainties']['gamma']
            sec["DE"][2] = self.data[zam]['decay_energy_uncertainties']['alpha']
            if 'DK' in sec.keys():
                i = 0
                for (rtyp, rfs), dk in self.data[zam]['decay_modes'].items():
                    sec['DK'][i]['RTYP'] = rtyp
                    sec['DK'][i]['RFS'] = rfs
                    sec['DK'][i]['BR'] = dk['branching_ratio']
                    sec['DK'][i]['DBR'] = dk['branching_ratio_uncertainty']
                    i += 1
            data[mat, mf, mt] = sandy.write_mf8(sec)
        return sandy.Endf6(data)


class _DecayBase():
    """
    Base class to perturb decay data

    Attributes
    ----------
    data
        best estimates and uncertainty or only best estimates as a dataframe

    Methods
    -------
    custom_perturbation
        apply custom perturbation to a given `BranchingRatio`, `DecayEnergy` 
        or `HalfLife` instance.
    """

    def __init__(self, df):
        self.data = pd.DataFrame(df)

    def __repr__(self):
        return self.data.__repr__()

    def custom_perturbation(self, pert):
        """
        Apply a custom perturbation to a given `BranchingRatio`, `DecayEnergy` 
        or `HalfLife` instance.

        Parameters
        ----------
        pert : `pandas.DataFrame`
            dataframe containing perturbation coefficients as ratio values,
            e.g., 1.05 for a perturbation of +5%.
            Depending on the nuclear data to perturb, `pert` index should be:
                * if perturbing branching ratio: "ZAM", "RTYP", "RFS"
                * if perturbing decay energy: "ZAM", "TYPE"
                * if perturbing half life: "ZAM"

        Returns
        -------
        `sandy.BranchingRatio`, `sandy.DecayEnergy` or `sandy.HalfLife`
            branching ratio, decay energy or half life instance with
            given values perturbed

        Examples
        --------
        Perturbation of 5% on the half life of U235:
        >>> endf6 = sandy.get_endf6_file("jeff_33", "decay", 922350)
        >>> rdd = sandy.DecayData.from_endf6(endf6)
        >>> hl = rdd.get_half_life(with_uncertainty=False)
        >>> pert = pd.DataFrame([{"ZAM": 922350, "PERT": 1.05}]).set_index(["ZAM"])
        >>> hl_new = hl.custom_perturbation(pert)
        >>> assert hl_new.data.values == hl.data.values * 1.05

        >>> hl = rdd.get_half_life()
        >>> hl_new = hl.custom_perturbation(pert)
        >>> assert hl_new.data.HL.values == hl.data.HL.values * 1.05

        Perturbation of 5% on the alpha decay energy of U235:
        >>> e = rdd.get_decay_energy(with_uncertainty=False)
        >>> pert = pd.DataFrame([{"ZAM": 922350, "TYPE": "alpha", "PERT": 1.05}]).set_index(["ZAM", "TYPE"])
        >>> e_new = e.custom_perturbation(pert)
        >>> assert e_new.data.E[922350]['alpha'] == e.data.E[922350]['alpha'] * 1.05

        >>> e = rdd.get_decay_energy()
        >>> e_new = e.custom_perturbation(pert)
        >>> assert e_new.data.E[922350]['alpha'] == e.data.E[922350]['alpha'] * 1.05

        Perturbation of 5% on the branching ratio for alpha decay of U235:
        >>> br = rdd.get_branching_ratio(with_uncertainty=False)
        >>> pert = pd.DataFrame([{"ZAM": 922350, "RTYP": 4, "RFS": 0, "PERT": 1.05}]).set_index(["ZAM", "RTYP", "RFS"])
        >>> br_new = br.custom_perturbation(pert)
        >>> assert br_new.data.BR[922350][4][0] == br.data.BR[922350][4][0] * 1.05

        >>> br = rdd.get_branching_ratio()
        >>> br_new = br.custom_perturbation(pert)
        >>> assert br_new.data.BR[922350][4][0] == br.data.BR[922350][4][0] * 1.05
        """
        name = "BR" if isinstance(self, BranchingRatio) else "E" if isinstance(self, DecayEnergy) else "HL"
        df = self.data.merge(pert.reindex(self.data.index).fillna(1), left_index=True, right_index=True)
        df[name] = df.PERT * df[name]
        return self.__class__(df.drop('PERT', axis=1))


class BranchingRatio(_DecayBase):
    """
    Extension of `sandy._DecayBase`. Container of best estimates and
    uncertainties of branching ratios.

    Methods
    -------
    normalize
        apply normalization condition to each row of `BranchingRatio.data`.

    to_decaydata
        update branching ratios in `DecayData` instance with those available in a
        `BranchingRatio` instance.
    """
    
    def normalize(self):
        """
        Normalize branching ratios.

        Returns
        -------
        `sandy.BranchingRatio`
            `BranchingRatio` object with normalized branching ratio values,
            thus respecting the constraint of their sum equal to one.

        Examples
        --------
        >>> endf6 = sandy.get_endf6_file("jeff_33", "decay", [942410, 922350])
        >>> rdd = sandy.DecayData.from_endf6(endf6)
        >>> br = rdd.get_branching_ratio()
        >>> br_norm = br.normalize()
        >>> assert br_norm.data.query("ZAM == 922350").BR.sum() == 1

        >>> br = rdd.get_branching_ratio(with_uncertainty=False)
        >>> br_norm = br.normalize()
        >>> assert br_norm.data.query("ZAM == 922350").sum().values == 1
        
        >>> endf6 = sandy.get_endf6_file("jeff_33", "decay", 942390)
        >>> rdd = sandy.DecayData.from_endf6(endf6)
        >>> br = rdd.get_branching_ratio()
        >>> br_norm = br.normalize()
        >>> assert br_norm.data.query("ZAM == 942390").BR.sum() == 1

        Stable nuclide:
        >>> endf6 = sandy.get_endf6_file("jeff_33", "decay", 260560)
        >>> rdd = sandy.DecayData.from_endf6(endf6)
        >>> br = rdd.get_branching_ratio()
        >>> br.normalize()
        Empty DataFrame
        Columns: [BR, DBR]
        Index: []
        """
        if self.data.empty:
            return self.__class__(self.data)
        foo = lambda x: x / x.sum()  # normalization function
        df = self.data.BR.to_frame().groupby('ZAM', group_keys=False).apply(foo)
        if 'DBR' in self.data.columns:
            df['DBR'] = self.data['DBR']
        return self.__class__(df)

    def to_decaydata(self, rdd):
        """
        Update branching ratios in `DecayData` instance with those available in
        a `BranchingRatio` instance.

        Parameters
        ----------
        `rdd` : `sandy.DecayData`
            `DecayData` instance

        Returns
        -------
        `sandy.DecayData`
            `DecayData` instance with updated branching ratios.

        Examples
        --------
        
        >>> endf6 = sandy.get_endf6_file("jeff_33", "decay", 922350)
        >>> rdd = sandy.DecayData.from_endf6(endf6)
        >>> br = rdd.get_branching_ratio(with_uncertainty=False)
        >>> pert = pd.DataFrame([{"ZAM": 922350, "RTYP": 4, "RFS": 0, "PERT": 1.05}]).set_index(["ZAM", "RTYP", "RFS"])
        >>> br_new = br.custom_perturbation(pert)
        >>> rdd_updated = br_new.to_decaydata(rdd)
        >>> assert rdd_updated.data[922350]['decay_modes'][(4, 0)]['branching_ratio'] == br_new.data.query("ZAM==922350 & RTYP==4 & RFS==0").BR.values
        
        >>> br = rdd.get_branching_ratio()
        >>> br_new = br.custom_perturbation(pert)
        >>> rdd_updated = br_new.to_decaydata(rdd)
        >>> assert rdd_updated.data[922350]['decay_modes'][(4, 0)]['branching_ratio'] == br_new.data.query("ZAM==922350 & RTYP==4 & RFS==0").BR.values
        
        Perturbing only one branching ratio of one nuclide in `DecayData` instance:
        >>> endf6 = sandy.get_endf6_file("jeff_33", "decay", [922350, 942410])
        >>> rdd = sandy.DecayData.from_endf6(endf6)
        >>> br = rdd.get_branching_ratio(with_uncertainty=False)
        >>> br_new = br.custom_perturbation(pert)
        >>> rdd_updated = br_new.to_decaydata(rdd)
        >>> assert rdd_updated.data[922350]['decay_modes'][(4, 0)]['branching_ratio'] == br_new.data.query("ZAM==922350 & RTYP==4 & RFS==0").BR.values
        >>> assert rdd_updated.data[942410]['decay_modes'][(4, 0)]['branching_ratio'] == br_new.data.query("ZAM==942410 & RTYP==4 & RFS==0").BR.values
        
        Perturbing only one branching ratio of each nuclide in `DecayData` instance:
        >>> endf6 = sandy.get_endf6_file("jeff_33", "decay", [922350, 942410])
        >>> rdd = sandy.DecayData.from_endf6(endf6)
        >>> br = rdd.get_branching_ratio(with_uncertainty=False)
        >>> pert = pd.DataFrame([{"ZAM": 922350, "RTYP": 4, "RFS": 0, "PERT": 1.05}, \
                                 {"ZAM": 942410, "RTYP": 4, "RFS": 0, "PERT": 1.02}]).set_index(["ZAM","RTYP", "RFS"])
        >>> br_new = br.custom_perturbation(pert)
        >>> rdd_updated =br_new.to_decaydata(rdd)
        >>> assert rdd_updated.data[922350]['decay_modes'][(4, 0)]['branching_ratio'] == br_new.data.query("ZAM==922350 & RTYP==4 & RFS==0").BR.values
        >>> assert rdd_updated.data[942410]['decay_modes'][(4, 0)]['branching_ratio'] == br_new.data.query("ZAM==942410 & RTYP==4 & RFS==0").BR.values
        
        Perturbing all branching ratios of each nuclide in `DecayData` instance:
        >>> endf6 = sandy.get_endf6_file("jeff_33", "decay", [922350, 942410])
        >>> rdd = sandy.DecayData.from_endf6(endf6)
        >>> br = rdd.get_branching_ratio(with_uncertainty=False)
        >>> pert = pd.DataFrame([{"ZAM": 922350, "RTYP": 4, "RFS": 0, "PERT": 1.05}, \
                                 {"ZAM": 922350, "RTYP": 6, "RFS": 0, "PERT": 0.95}, \
                                 {"ZAM": 942410, "RTYP": 4, "RFS": 0, "PERT": 1.02}, \
                                 {"ZAM": 942410, "RTYP": 1, "RFS": 0, "PERT": 0.99}]).set_index(["ZAM", "RTYP", "RFS"])
        >>> br_new = br.custom_perturbation(pert)
        >>> rdd_updated = br_new.to_decaydata(rdd)
        >>> assert rdd_updated.data[922350]['decay_modes'][(4, 0)]['branching_ratio'] == br_new.data.query("ZAM==922350 & RTYP==4 & RFS==0").BR.values
        >>> assert rdd_updated.data[922350]['decay_modes'][(6, 0)]['branching_ratio'] == br_new.data.query("ZAM==922350 & RTYP==6 & RFS==0").BR.values
        >>> assert rdd_updated.data[942410]['decay_modes'][(4, 0)]['branching_ratio'] == br_new.data.query("ZAM==942410 & RTYP==4 & RFS==0").BR.values
        >>> assert rdd_updated.data[942410]['decay_modes'][(1, 0)]['branching_ratio'] == br_new.data.query("ZAM==942410 & RTYP==1 & RFS==0").BR.values
        """
        rdd_updated = copy.deepcopy(rdd.data)
        for (zam, rtyp, rfs), val in self.data.iterrows():
            rdd_updated[zam]['decay_modes'][(rtyp, rfs)]['branching_ratio'] = val['BR']
        return DecayData(rdd_updated)


class HalfLife(_DecayBase):
    """
    Extension of `sandy._DecayBase`. Container of best estimates and
    uncertainties of half lives.

    Methods
    -------
    to_decaydata
        update half lives in `DecayData` instance with those available in a
        `HalfLife` instance.
    """

    def to_decaydata(self, rdd):
        """
        Update half lives in `DecayData` instance with those available in a
        `HalfLife` instance.
        Decay consants are also recaluclated and updated accordingly.
        Uncertaintoes are ignored.

        Parameters
        ----------
        `rdd` : `sandy.DecayData`
            `DecayData` instance

        Returns
        -------
        `sandy.DecayData`
            `DecayData` instance with updated half lives.

        Examples
        --------
        >>> endf6 = sandy.get_endf6_file("jeff_33", "decay", 922350)
        >>> rdd = sandy.DecayData.from_endf6(endf6)
        >>> hl = rdd.get_half_life(with_uncertainty=False)
        >>> pert = pd.DataFrame([{"ZAM": 922350, "PERT": 1.05}]).set_index(["ZAM"])
        >>> hl_new = hl.custom_perturbation(pert)
        >>> rdd_updated = hl_new.to_decaydata(rdd)
        >>> assert rdd_updated.data[922350]['half_life'] == hl_new.data.values
        >>> assert rdd_updated.data[922350]['decay_constant'] !=rdd.data[922350]['decay_constant']
        
        >>> hl = rdd.get_half_life()
        >>> hl_new = hl.custom_perturbation(pert)
        >>> rdd_updated = hl_new.to_decaydata(rdd)
        >>> assert rdd_updated.data[922350]['half_life'] == hl_new.data.HL.values
        
        Perturbing only half life of one nuclide in `DecayData` instance:
        >>> endf6 = sandy.get_endf6_file("jeff_33", "decay", [922350, 942410])
        >>> rdd = sandy.DecayData.from_endf6(endf6)
        >>> hl = rdd.get_half_life(with_uncertainty=False)
        >>> pert = pd.DataFrame([{"ZAM": 922350, "PERT": 1.05}]).set_index(["ZAM"])
        >>> hl_new = hl.custom_perturbation(pert)
        >>> rdd_updated = hl_new.to_decaydata(rdd)
        >>> assert rdd_updated.data[922350]['half_life'] == hl_new.data.query('ZAM==922350').HL.values
        >>> assert rdd_updated.data[942410]['half_life'] == hl_new.data.query('ZAM==942410').HL.values
        
        Perturbing half life of each nuclide in `DecayData` instance:
        >>> pert = pd.DataFrame([{"ZAM": 922350,"PERT": 1.05},\
                                 {"ZAM": 942410,"PERT": 1.02}]).set_index(["ZAM"])
        >>> hl_new = hl.custom_perturbation(pert)
        >>> rdd_updated = hl_new.to_decaydata(rdd)
        >>> assert rdd_updated.data[922350]['half_life'] == hl_new.data.query('ZAM==922350').HL.values
        >>> assert rdd_updated.data[942410]['half_life'] == hl_new.data.query('ZAM==942410').HL.values
        """
        rdd_updated = copy.deepcopy(rdd.data)
        for zam, val in self.data.iterrows():
            # update half life and recalculate decay constant.
            # Do not update uncertainties.
            rdd_updated[zam]['half_life'] = val['HL']
            if val['HL'] != 0:
                rdd_updated[zam]['decay_constant'] = np.log(2) / val["HL"]
        return DecayData(rdd_updated)


class DecayEnergy(_DecayBase):
    """
    Extension of `sandy._DecayBase`. Container of best estimates and
    uncertainties of decay energies.

    Methods
    -------
    to_decaydata
        update decay energies in `DecayData` instance with those available in a
        `DecayEnergy` instance.
    """

    def to_decaydata(self, rdd):
        """
        Update decay energies in `DecayData` instance with those available in a
        `DecayEnergy` instance.

        Parameters
        ----------
        `rdd` : `sandy.DecayData`
            `DecayData` instance
        Returns
        -------
        `sandy.DecayData`
            `DecayData` instance with updated decay energies.

        Examples
        --------
        >>> endf6 = sandy.get_endf6_file("jeff_33", "decay", 922350)
        >>> rdd = sandy.DecayData.from_endf6(endf6)
        >>> e = rdd.get_decay_energy(with_uncertainty=False)
        >>> pert = pd.DataFrame([{"ZAM": 922350, "TYPE": "alpha", "PERT": 1.05}]).set_index(["ZAM", "TYPE"])
        >>> e_new = e.custom_perturbation(pert)
        >>> rdd_updated = e_new.to_decaydata(rdd)
        >>> assert rdd_updated.data[922350]['decay_energy']['alpha'] == e_new.data.E[922350]['alpha']
        
        >>> e = rdd.get_decay_energy()
        >>> e_new = e.custom_perturbation(pert)
        >>> rdd_updated = e_new.to_decaydata(rdd)
        >>> assert rdd_updated.data[922350]['decay_energy']['alpha'] == e_new.data.E[922350]['alpha']
        
        Perturbing only one decay energy of one nuclide in `DecayData` instance:
        >>> endf6 = sandy.get_endf6_file("jeff_33", "decay", [922350, 942410])
        >>> rdd = sandy.DecayData.from_endf6(endf6)
        >>> e = rdd.get_decay_energy(with_uncertainty=False)
        >>> e_new = e.custom_perturbation(pert)
        >>> rdd_updated =e_new.to_decaydata(rdd)
        >>> assert rdd_updated.data[922350]['decay_energy']['alpha'] == e_new.data.E[922350]['alpha']
        >>> assert rdd_updated.data[942410]['decay_energy']['alpha'] == e_new.data.E[942410]['alpha']
        
        Perturbing one decay energy of each nuclide in `DecayData` instance:
        >>> pert = pd.DataFrame([{"ZAM": 922350, "TYPE": "alpha", "PERT": 1.05}, \
                                 {"ZAM": 942410, "TYPE": "alpha", "PERT": 1.05}]).set_index(["ZAM", "TYPE"])
        >>> e_new = e.custom_perturbation(pert)
        >>> rdd_updated =e_new.to_decaydata(rdd)
        >>> assert rdd_updated.data[922350]['decay_energy']['alpha'] == e_new.data.E[922350]['alpha']
        >>> assert rdd_updated.data[942410]['decay_energy']['alpha'] == e_new.data.E[942410]['alpha']
        
        Perturbing all decay energies of each nuclide in `DecayData` instance:
        >>> pert = pd.DataFrame([{"ZAM": 922350, "TYPE": "alpha", "PERT": 1.05}, \
                                 {"ZAM": 922350, "TYPE": "beta", "PERT": 1.01}, \
                                 {"ZAM": 922350, "TYPE": "gamma", "PERT": 0.97}, \
                                 {"ZAM": 942410, "TYPE": "alpha", "PERT": 1.05}, \
                                 {"ZAM": 942410, "TYPE": "beta", "PERT": 0.98}, \
                                 {"ZAM": 942410, "TYPE": "gamma", "PERT": 1.02}]).set_index(["ZAM", "TYPE"])
        >>> e_new = e.custom_perturbation(pert)
        >>> rdd_updated =e_new.to_decaydata(rdd)
        >>> assert rdd_updated.data[922350]['decay_energy']['alpha'] == e_new.data.E[922350]['alpha']
        >>> assert rdd_updated.data[922350]['decay_energy']['beta'] == e_new.data.E[922350]['beta']
        >>> assert rdd_updated.data[922350]['decay_energy']['gamma'] == e_new.data.E[922350]['gamma']
        >>> assert rdd_updated.data[942410]['decay_energy']['alpha'] == e_new.data.E[942410]['alpha']
        >>> assert rdd_updated.data[942410]['decay_energy']['beta'] == e_new.data.E[942410]['beta']
        >>> assert rdd_updated.data[942410]['decay_energy']['gamma'] == e_new.data.E[942410]['gamma']
        """
        rdd_updated = copy.deepcopy(rdd.data)
        for (zam, typ), val in self.data.iterrows():
            rdd_updated[zam]['decay_energy'][typ] = val['E']
        return DecayData(rdd_updated)


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
    Expand beta decay (#1)
    >>> d, n, p, a = sandy.decay.expand_decay_type(581480, 1)
    >>> assert d == 591480
    >>> assert n == 0
    >>> assert p == 0
    >>> assert a == 0

    Expand electron capture and/or positron emission (#2)
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
    >>> import pytest
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
    rtyp : `int`
        integer of decay modes where:
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

     Notes
     -----
     .. note:: rtyp=0, corresponding to gamma ray decay, is not used in MF=8, MT=457 section,
         according with what reported in the ENDF6 manual.

    Examples
    --------
    Extract products of fake decay process including all available decay modes.
    >>> sandy.decay.get_decay_products(1234567, 581480)
    {551420: 1.0, 10: 1.0, 10010: 1.0, 20040: 1.0}

    ...change the metastate of the product
    >>> sandy.decay.get_decay_products(1234567, 581480, meta=1)
    {551421: 1.0, 10: 1.0, 10010: 1.0, 20040: 1.0}

    ...and then use a different braanching ratio
    >>> sandy.decay.get_decay_products(1234567, 581480, br=0.1)
    {551420: 0.1, 10: 0.1, 10010: 0.1, 20040: 0.1}
    """
    daughter = zam + 0
    neutrons = 0.
    protons = 0.
    alphas = 0.
    for dectyp in map(int, str(rtyp)):
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
