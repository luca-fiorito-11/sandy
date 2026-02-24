import pandas as pd
import numpy as np
import logging

from .endf6 import _FormattedFile
from .cov import CategoryCov
from .xs import Xs
from .records import read_cont, read_list

__author__ = "Luca Fiorito"
__all__ = [
        "Errorr",
        ]

pd.options.display.float_format = '{:.5e}'.format


class Errorr(_FormattedFile):
    """
    Container for ERRORR file text grouped by MAT, MF and MT numbers.
    
    Methods
    -------
    get_cov
        Extract mulitgroup covariance matrix.
    get_energy_grid
        Extract breaks of multi-group energy grid from ERRORR output file.
    get_xs
        Extract multigroup xs values.
    """

    def get_energy_grid(self, **kwargs):
        """
        Extract breaks of multi-group energy grid from ERRORR output file.

        Parameters
        ----------
        mat : `int`, optional
            MAT number. The default is None.

        Returns
        -------
        `np.array`
            The energy grid of the :obj:`~sandy.errorr.Errorr` object.

        Examples
        --------

        This example shows functioning of the `get_cov` method with `MF=31`, `MF=33` and `MF=35`.
        The first case if for `MF=31`.

        >>> import sandy
        >>> e6 = sandy.get_endf6_file("jeff_33", "xs", 10010, local=True)
        >>> ek = sandy.energy_grids.CASMO12
        >>> err = e6.get_errorr(errorr_kws=dict(ek=ek), err=1)['errorr33']
        >>> np.testing.assert_allclose(err.get_energy_grid(), ek, atol=1e-14, rtol=1e-14)
        >>> np.testing.assert_allclose(err.get_energy_grid(mat=125), ek, atol=1e-14, rtol=1e-14)
        """
        mat_ = kwargs.get('mat', self.mat[0])
        mf1 = read_mf1(self, mat_)
        return mf1["EG"]

    def get_xs(self, mts=None, **kwargs):
        """
        Extract multigroup xs/nubar/pfns values.

        Parameters
        ----------
        mts : `list` of `int`, optional
            MT number(s) to extract. Default is `None` (all available MTs).

        Returns
        -------
        xs : `pd.DataFrame`
            MultiIndex DataFrame with XS values for the selected MAT/MT and energy groups.
            Index: energy intervals (E)
            Columns: MultiIndex (MAT, MT)

        Examples
        --------
        This example shows how to use of the `get_xs` method with `MF=31`, `MF=33` and `MF=35`.

        >>> import sandy, pytest
        >>> e6 = sandy.get_endf6_file("jeff_33", "xs", 922350, local=True)
        >>> ek = sandy.energy_grids.CASMO12
        >>> errs = e6.get_errorr(err=1, errorr_kws=dict(ek=ek), groupr_kws=dict(ek=ek))

        The first case is for `MF=33`.

        >>> errs['errorr33'].get_xs().data.iloc[:, :3]
        MAT                            9228                        
        MT                                1           2           4
        E                                                          
        (1e-05, 0.03]           1.33136e+03 1.40944e+01 0.00000e+00
        (0.03, 0.058]           5.26782e+02 1.40005e+01 0.00000e+00
        (0.058, 0.14]           3.41410e+02 1.38400e+01 0.00000e+00
        (0.14, 0.28]            2.50115e+02 1.35755e+01 0.00000e+00
        (0.28, 0.35]            2.22318e+02 1.35530e+01 0.00000e+00
        (0.35, 0.625]           1.31990e+02 1.33804e+01 0.00000e+00
        (0.625, 4.0]            5.56383e+01 1.21305e+01 0.00000e+00
        (4.0, 48.052]           9.39232e+01 1.14789e+01 0.00000e+00
        (48.052, 5530.0]        2.19874e+01 1.18636e+01 8.60022e-07
        (5530.0, 821000.0]      9.53303e+00 6.75610e+00 1.22176e+00
        (821000.0, 2231000.0]   7.06399e+00 3.78697e+00 1.94841e+00
        (2231000.0, 10000000.0] 6.95725e+00 3.82719e+00 1.40637e+00

        Use `mts` to select only specific MT numbers.

        >>> errs['errorr33'].get_xs(mts=[18])
        MAT                            9228
        MT                               18
        E                                  
        (1e-05, 0.03]           1.11449e+03
        (0.03, 0.058]           4.39929e+02
        (0.058, 0.14]           2.80260e+02
        (0.14, 0.28]            1.95476e+02
        (0.28, 0.35]            1.72670e+02
        (0.35, 0.625]           1.02640e+02
        (0.625, 4.0]            3.15194e+01
        (4.0, 48.052]           4.87917e+01
        (48.052, 5530.0]        7.12775e+00
        (5530.0, 821000.0]      1.29401e+00
        (821000.0, 2231000.0]   1.23450e+00
        (2231000.0, 10000000.0] 1.39786e+00

        An error is raised if no requested MT number is found.

        >>> with pytest.raises(ValueError) as exc:
        ...     errs['errorr33'].get_xs(mts=[10])
        >>> assert str(exc.value) == 'No requested MT number [10] was found in ERRORR file'

        The first case is for `MF=31`.

        >>> errs["errorr31"].get_xs()
        MAT                            9228
        MT                              456
        E                                  
        (1e-05, 0.03]           2.40910e+00
        (0.03, 0.058]           2.40910e+00
        (0.058, 0.14]           2.40910e+00
        (0.14, 0.28]            2.40910e+00
        (0.28, 0.35]            2.40910e+00
        (0.35, 0.625]           2.40910e+00
        (0.625, 4.0]            2.40910e+00
        (4.0, 48.052]           2.40910e+00
        (48.052, 5530.0]        2.40930e+00
        (5530.0, 821000.0]      2.44782e+00
        (821000.0, 2231000.0]   2.57604e+00
        (2231000.0, 10000000.0] 3.28604e+00

        The third case is for `MF=35`.

        >>> errs["errorr35"].get_xs()
        MAT                            9228
        MT                               18
        E                                  
        (1e-05, 0.03]           2.14656e-12
        (0.03, 0.058]           3.63873e-12
        (0.058, 0.14]           1.59188e-11
        (0.14, 0.28]            3.95395e-11
        (0.28, 0.35]            2.42436e-11
        (0.35, 0.625]           1.19290e-10
        (0.625, 4.0]            3.10112e-09
        (4.0, 48.052]           1.34526e-07
        (48.052, 5530.0]        1.69897e-04
        (5530.0, 821000.0]      2.43142e-01
        (821000.0, 2231000.0]   4.11862e-01
        (2231000.0, 10000000.0] 3.44826e-01

        """
        # Check if MT numbers (if requested) exist
        if mts is not None:
            requested_mts = set(mts)
            available_mts = set(self.mt)
            if not requested_mts & available_mts:
                raise ValueError(f"No requested MT number {sorted(requested_mts)} was found in ERRORR file")
            
            notfound = requested_mts - available_mts
            if notfound:
                logging.warning(f"The following MT's were not found in ERRORR file: {sorted(notfound)}")

        listmt_ = mts if mts is not None else range(1, 10000)

        listmf_ = [3, 5]   # not 1 because nubar is given in MF3

        filtered_data = self.filter_by(listmf=listmf_, listmt=listmt_)

        data = []
        for mat, mf, mt in filtered_data.data:

            # This reads intro and MG neutron flux
            mf1 = read_mf1(self, mat)
            energy_index  = pd.IntervalIndex.from_breaks(mf1["EG"], name="E")

            # This reads all nubar, xs and pfns
            mf3 = read_mf3(self, mat, mt, mf)
            columns = pd.MultiIndex.from_tuples([(mat, mt)], names=["MAT", "MT"])

            data.append(pd.DataFrame(mf3["XS"], index=energy_index , columns=columns))

        # Concatenate all MT/MAT columns and fill missing energy bins with 0
        data = pd.concat(data, axis=1).fillna(0)
        xs = Xs(data)
        
        return xs

    def get_cov(self, mts=None, covariance_checks=True):
        """
        Extract cross section/nubar covariance from this :obj:`~sandy.errorr.Errorr` instance.

        Parameters
        ----------
        mts : `list`, optional
            List of MT numbers. If `None` (default), keep all MT numbers.
            Use this command if you want to keep only a subsection of the
            MT numbers available in the ERRORR file.
            The output covariance matrix will containt only the MT numbers
            found both in `mts` and in the ERRORR file.
            If some requested MTs are missing, a warning is logged.
            If none are found, a `ValueError` is raised.
        covariance_checks : `bool`, optional, default is `True`
            Perform symmetry and variance checks on the covariance matrix.
            Set to `False` to disable these checks (useful for debugging).
            See also :obj:`~sandy.cov.CategoryCov.data`.

        Returns
        -------
        :obj:`~sandy.cov.CategoryCov`
            Covariance matrix of xs/nubar/pfns for all MAT/MT pairs found in ERRORR file.

        Raises
        ------
        `ValueError`
            If none of the requested MT numbers are found in the ERRORR file.

        Examples
        --------

        Read cross section covariance matrix for simple case (H1).

        >>> import sandy, pytest
        >>> e6 = sandy.get_endf6_file("jeff_33", "xs", 10010, local=True)
        >>> err = e6.get_errorr(errorr_kws=dict(ek=[1e-2, 1e1, 2e7]), err=1, temperature=0.1)['errorr33']
        >>> datamg = err.get_cov().data
        >>> datamg
      		MAT	                   125
          MT	                       1	                            2	                                102
                  E	               (0.01, 10.0]	(10.0, 20000000.0]	(0.01, 10.0]	(10.0, 20000000.0]	(0.01, 10.0]	(10.0, 20000000.0]
          MAT	 MT	                 E						
          125	  1	      (0.01, 10.0]	8.74835e-06 	       4.62555e-05  8.76099e-06 	    4.62566e-05 1.07148e-06	           5.59219e-07
                  (10.0, 20000000.0]	4.62555e-05	       2.47644e-04	 4.63317e-05	       2.47649e-04	7.58743e-09	           1.49541e-06
                2	      (0.01, 10.0]	8.76099e-06	       4.63317e-05	 8.77542e-06	       4.63327e-05	0.00000e+00	           0.00000e+00
                  (10.0, 20000000.0]	4.62566e-05	       2.47649e-04	 4.63327e-05	       2.47655e-04	0.00000e+00	           0.00000e+00
              102	      (0.01, 10.0]	1.07148e-06 	       7.58743e-09	 0.00000e+00	       0.00000e+00	6.51764e-04	           3.40163e-04
                  (10.0, 20000000.0]	5.59219e-07 	       1.49541e-06	 0.00000e+00	       0.00000e+00	3.40163e-04	           6.70430e-02
      
        This example shows the cross correlation among two MT taken from a
        cross section covariance matrix (MF=33) with 3 energy groups at high
        energy.
        There is no correlation in the last two groups. 

        >>> tape = sandy.get_endf6_file("jeff_33", "xs", 641530, local=True)
        >>> out = tape.get_errorr(err=1, errorr33_kws=dict(irespr=0, mt=[1, 51, 52],ek=[1e7,2e7,2.4e7, 2.8e7]))
        >>> cov = out["errorr33"].get_cov().data
        >>> cov.loc[(6428,1)][(6428,51)]
        E                         (10000000.0, 20000000.0]  (20000000.0, 24000000.0]  (24000000.0, 28000000.0]
        E
        (10000000.0, 20000000.0]               8.66651e-03               0.00000e+00               0.00000e+00
        (20000000.0, 24000000.0]               6.81128e-02               0.00000e+00               0.00000e+00
        (24000000.0, 28000000.0]               7.52293e-02               0.00000e+00               0.00000e+00

        This example shows how to use the `get_cov` method with `MF=31`, `MF=33` and `MF=35`.
        The first case is for `MF=31`.
        
        >>> import numpy as np
        >>> e6 = sandy.get_endf6_file("jeff_33", "xs", 922350, local=True)
        >>> err = e6.get_errorr(errorr_kws=dict(ek=[1e-2, 1e1, 2e7]), groupr_kws=dict(ek=[1e-2, 1e1, 2e7]), err=1, xs=False, chi=False, nubar=True, mubar=False)['errorr31']
        >>> datamg = err.get_cov().data
        >>> np.testing.assert_equal(datamg.values, [[3.153674e-05, 1.413344e-05],[1.413344e-05, 1.643044e-05]])

        The second case is for `MF=33`.
        
        >>> e6 = sandy.get_endf6_file("jeff_33", "xs", 922350, local=True)
        >>> err = e6.get_errorr(errorr_kws=dict(ek=[1e-2, 1e1, 2e7]), err=1, xs=True, chi=False, nubar=False, mubar=False)['errorr33']
        >>> datamg = err.get_cov().data
        >>> np.testing.assert_equal(datamg.loc[(9228, 1), (9228, 1)].values, [[2.060002e-04, 6.686222e-08],[6.686222e-08, 7.581125e-05]])

        The third case is for `MF=35`.

        >>> e6 = sandy.get_endf6_file("jeff_33", "xs", 922350, local=True)
        >>> err = e6.get_errorr(errorr_kws=dict(ek=[1e-2, 1e1, 2e7]), groupr_kws=dict(ek=[1e-2, 1e1, 2e7]), err=1, xs=False, chi=True, nubar=False, mubar=False)['errorr35']
        >>> datamg = err.get_cov().data
        >>> np.testing.assert_equal(datamg.values, [[1.750390e-03, 4.450283e-08],[4.450283e-08, 1.622930e-10]])
        
        In some cases, ERRORR produces non-symmetric covariance matrices.
        The generation of a `CategoryCov` can be enforced with `covariance_checks=False`.

        >>> e6 = sandy.get_endf6_file("jeff_33", "xs", 922350, local=True)
        >>> err = e6.get_errorr(errorr_kws=dict(ek=sandy.energy_grids.SCALE238), groupr_kws=dict(ek=sandy.energy_grids.SCALE238), err=1, xs=False, chi=True, nubar=False, mubar=False)['errorr35']
        >>> with pytest.raises(TypeError) as excinfo:
        ...     err.get_cov().data
        >>> cov = err.get_cov(covariance_checks=False).data

        Test selecting only specific MT's.

        >>> err = sandy.get_endf6_file("jeff_33", "xs", 10010, local=True).get_errorr(err=1)["errorr33"]
        >>> cov = err.get_cov()
        >>> np.testing.assert_array_equal(cov.data.index.get_level_values("MT").unique(), [1, 2, 102])
        >>> with pytest.raises(Exception) as exc:
        ...    err.get_cov(mts={4, 452})
        >>> assert str(exc.value) == 'No requested MT number [4, 452] was found in ERRORR file'

        """
        # Check if MT numbers (if requested) exist
        if mts is not None:
            requested_mts = set(mts)
            available_mts = set(self.mt)
            if not requested_mts & available_mts:
                raise ValueError(f"No requested MT number {sorted(requested_mts)} was found in ERRORR file")
            
            notfound = requested_mts - available_mts
            if notfound:
                logging.warning(f"The following MT's were not found in ERRORR file: {sorted(notfound)}")

        # Retrieve multigroup energy grid
        energy_grid = pd.IntervalIndex.from_breaks(self.get_energy_grid())

        # Filter ERRORR data for relevant MF sections
        filtered_data = self.filter_by(listmf=[31, 33, 34, 35]).data

        # initialize global cov matrix with all MAT, MT
        ix = pd.DataFrame(filtered_data.keys(), columns=["MAT", "MF", "MT"])[["MAT", "MT"]]
        ix["IMIN"] = ix.index * energy_grid.size
        ix["IMAX"] = (ix.index + 1) * energy_grid.size

        total_size = ix.shape[0] * energy_grid.size
        c = np.zeros((total_size, total_size))
        
        # Fill covariance matrix
        for mat, mf, mt in filtered_data:
            mf33 = read_mf33(self, mat, mt, 33 if mf == 31 else mf)

            for mt1, subcov in mf33["COVS"].items():
                
                # it seems that when processing MF34 mubar, NJOY keeps 251 for MT
                # but it sets MT1 to 1.
                # Here we manually set MT1 back to 251.
                if mf == 34 and mt1 == 1:
                    mt1 = 251

                i = ix.query("MAT==@mat & MT==@mt").squeeze()
                j = ix.query("MAT==@mat & MT==@mt1").squeeze()

                c[i.IMIN: i.IMAX, j.IMIN: j.IMAX] = subcov
                if mt != mt1:
                    c[j.IMIN: j.IMAX, i.IMIN: i.IMAX] = subcov.T
        
        # Build MultiIndex for rows/columns
        idx = pd.MultiIndex.from_tuples(
            [(mat, mt, e) for i, (mat, mt) in ix[["MAT", "MT"]].iterrows() for e in energy_grid],
            names=["MAT", "MT", "E"],
        )
  
        # Filter by requested MTs if provided
        if mts is not None:
            mask = idx.get_level_values("MT").isin(requested_mts)
            c = c[mask][:, mask]
            idx = idx[mask]
  
        out = CategoryCov(c, index=idx, columns=idx, covariance_checks=covariance_checks)

        return out


def read_mf1(tape, mat):
    """
    Parse MAT/MF=1/MT=451 section from :obj:`~sandy.errorr.Errorr` object and return
    structured content in nested dcitionaries.

    Parameters
    ----------
    tape : :obj:`~sandy.errorr.Errorr`
        endf6 object containing requested section
    mat : `int`
        MAT number

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
    C, i = read_cont(df, i)
    add = {
        "ZA": C.C1,
        "AWR": C.C2,
        "LRP": C.N1,
    }
    out.update(add)
    L, i = read_list(df, i)
    add = {
        "EG": np.array(L.B),
    }
    out.update(add)
    return out


def read_mf3(tape, mat, mt, mf=3):
    """
    Parse MAT/MF=33/MT section from :obj:`~sandy.errorr.Errorr` object and return
    structured content in nested dcitionaries.

    Parameters
    ----------
    tape : :obj:`~sandy.errorr.Errorr`
        ERRORR object containing requested section
    mat : `int`
        MAT number
    mt : `int`
        MT number

    Returns
    -------
    out : `dict`
        Content of the ENDF-6 tape structured as nested `dict`.
    """
    df = tape._get_section_df(mat, mf, mt)
    out = {
            "MAT": mat,
            "MF": mf,
            "MT": mt,
            }
    i = 0
    L, i = read_list(df, i)
    add = {
        "XS": np.array(L.B),
    }
    out.update(add)
    return out


def read_mf33(tape, mat, mt, mf=33):
    """
    Parse MAT/MF=33/MT section from :obj:`~sandy.errorr.Errorr` object and return
    structured content in nested dcitionaries.

    Parameters
    ----------
    tape : :obj:`~sandy.error.Errorr`
        ERRORR object containing requested section
    mat : `int`
        MAT number
    mt : `int`
        MT number
    mf : `int`, optional
        MF number. Default is `33`.

    Notes
    -----
    .. note: ERRORR sections for nubar and xs are all given by NJOY in `MF=33`.
             For PFNS, `MF=35` is used.
             Independently, all sections can be parse by this function.

    Returns
    -------
    out : `dict`
        Content of the ERRORR tape structured as nested `dict`.
    """
    df = tape._get_section_df(mat, mf, mt)
    out = {
            "MAT": mat,
            "MF": mf,
            "MT": mt,
            }
    i = 0
    C, i = read_cont(df, i)
    add = {
        "ZA": C.C1,
        "AWR": C.C2,
    }
    out.update(add)
    reaction_pairs = {}
    for rp in range(C.N2):  # number of reaction pairs
        C, i = read_cont(df, i)
        MT1 = C.L2
        NG = C.N2
        M = np.zeros((NG, NG))
        while True:
            L, i = read_list(df, i)
            NGCOL = L.L1
            GROW = L.N2
            GCOL = L.L2
            M[GROW-1, GCOL-1:GCOL+NGCOL-1] = L.B
            if GROW >= NG:
                break
        reaction_pairs[MT1] = M
    out["COVS"] = reaction_pairs
    return out
