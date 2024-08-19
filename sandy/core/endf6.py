# -*- coding: utf-8 -*-
"""
Created on Wed Dec  4 14:50:33 2019

@author: lfiorito
"""
import io
import os
import shutil
from os.path import dirname, join
from functools import reduce
from tempfile import TemporaryDirectory
import logging
import urllib
from urllib.request import urlopen, Request
from zipfile import ZipFile
import re
import warnings

import multiprocessing as mp
import numpy as np
import pandas as pd
import numpy as np

import sandy
from sandy.libraries import (
    N_FILES_ENDFB_71_IAEA,
    N_FILES_ENDFB_80_IAEA,
    N_FILES_JEFF_32_NEA,
    N_FILES_JEFF_311_IAEA,
    N_FILES_JEFF_33_IAEA,
    N_FILES_JEFF_40T0_NEA,
    N_FILES_JENDL_40U_IAEA,
    N_FILES_TENDL_2023_PSI,
    N_FILES_IRDFF_2_IAEA,
    URL_N_ENDFB_71_IAEA,
    URL_N_JEFF_32_NEA,
    URL_N_JEFF_311_IAEA,
    URL_N_JEFF_33_IAEA,
    URL_N_JEFF_40T0_NEA,
    URL_N_ENDFB_80_IAEA,
    URL_N_JENDL_40U_IAEA,
    URL_N_TENDL_2023_PSI,
    URL_N_IRDFF_2_IAEA,

    NFPY_FILES_ENDFB_71_IAEA,
    NFPY_FILES_ENDFB_80_IAEA,
    NFPY_FILES_JEFF_311_IAEA,
    NFPY_FILES_JEFF_33_IAEA,
    NFPY_FILES_JENDL_40U_IAEA,
    URL_NFPY_ENDFB_71_IAEA,
    URL_NFPY_ENDFB_80_IAEA,
    URL_NFPY_JEFF_311_IAEA,
    URL_NFPY_JEFF_33_IAEA,
    URL_NFPY_JENDL_40U_IAEA,

    DECAY_FILES_ENDFB_71_IAEA,
    DECAY_FILES_ENDFB_80_IAEA,
    DECAY_FILES_JEFF_311_IAEA,
    DECAY_FILES_JEFF_33_IAEA,
    URL_DECAY_ENDFB_71_IAEA,
    URL_DECAY_ENDFB_80_IAEA,
    URL_DECAY_JEFF_311_IAEA,
    URL_DECAY_JEFF_33_IAEA,

    TSL_FILES_ENDFB_71_IAEA,
    TSL_FILES_ENDFB_80_IAEA,
    TSL_FILES_JEFF_33_IAEA,
    TSL_FILES_JENDL_40U_IAEA,
    URL_TSL_JENDL_40U_IAEA,
    URL_TSL_ENDFB_71_IAEA,
    URL_TSL_ENDFB_80_IAEA,
    URL_TSL_JEFF_33_IAEA,

    DXS_FILES_JEFF_33_IAEA,
    DXS_FILES_PROTON_IAEA,
    URL_DXS_JEFF_33_IAEA,
    URL_DXS_PROTON_IAEA
    )


__author__ = "Luca Fiorito"
__all__ = [
        "Endf6",
        "get_endf6_file",
        "get_tsl_index",
        ]

pd.options.display.float_format = '{:.5e}'.format


def get_tsl_index(library):
    """
    Obtain the index information available in the library web page.

    Parameters
    ----------
    library : `str`
        nuclear data library. Available libraries are:
        for 'tsl'
            * `'endfb_71'`
            * `'jeff_33'`
            * `'endfb_80'`
            * `'jendl_40u`
            * `'irdff_ii`

    Raises
    ------
    ValueError
        if library is not among available selection.

    Example
    ------
    >>> sandy.endf6.get_tsl_index("jendl_40u")
     Lib:         JENDL-4.0
     Library:     JENDL-4.0 Japanese evaluated nuclear data library, 2010
     Sub-library: NSUB=12      Thermal Neutron Scattering Data
    --------------------------------------------------------------------------------
       #)  KEY Material     Lab.         Date         Authors
    --------------------------------------------------------------------------------
       1)    1 1-H(H2O)     LANL         EVAL-apr93   MACFARLANE                         20.MeV   tsl_0001_h(h2o).zip 412Kb
       2)    2 1-Para-H     LANL         EVAL-APR93   MacFarlane                         20.MeV   tsl_0002_para-H.zip 91Kb 
       3)    3 1-Ortho-H    LANL         EVAL-APR93   MacFarlane                         20.MeV   tsl_0003_ortho-H.zip 96Kb
       4)    7 1-H(ZrH)     LANL         EVAL-apr93   MACFARLANE                         20.MeV   tsl_0007_h(zrh).zip 448Kb
       5)   11 1-D(D2O)     GA           EVAL-DEC69   KOPPEL,HOUSTON                     20.MeV   tsl_0011_D(D2O).zip 235Kb
       6)   12 1-Para-D     LANL         EVAL-APR93   MacFarlane                         20.MeV   tsl_0012_para-d.zip 92Kb 
       7)   13 1-Ortho-D    LANL         EVAL-APR93   MacFarlane                         20.MeV   tsl_0013_ortho-d.zip 93Kb
       8)   26 4-Be-metal   LANL         EVAL-apr93   MACFARLANE                         20.MeV   tsl_0026_bemetal.zip 419Kb
       9)   27 4-BeO        LANL         EVAL-apr93   MACFARLANE                         20.MeV   tsl_0027_beo.zip 483Kb   
      10)   31 6-Graphite   LANL         EVAL-apr93   MACFARLANE                         20.MeV   tsl_0031_graphite.zip 397Kb
      11)   33 6-l-CH4      LANL         EVAL-APR93   MacFarlane                         20.MeV   tsl_0033_l-ch4.zip 50Kb  
      12)   34 6-s-CH4      LANL         EVAL-APR93   MacFarlane                         20.MeV   tsl_0034_s-ch4.zip 42Kb  
      13)   37 6-H(CH2)     GA           EVAL-DEC69   KOPPEL,HOUSTON,SPREVAK             20.MeV   tsl_0037_H(CH2).zip 72Kb 
      14)   40 6-BENZINE    GA           EVAL-DEC69   KOPPEL,HOUSTON,BORGONOVI           20.MeV   tsl_0040_BENZINE.zip 236Kb
      15)   58 40-Zr(ZrH)   LANL         EVAL-apr93   MACFARLANE                         20.MeV   tsl_0058_zr(zrh).zip 201Kb
    --------------------------------------------------------------------------------
    Total: Materials:15 Size:11Mb Compressed:4Mb
    """
    available_libs = (
            "endfb_71".upper(),
            "endfb_80".upper(),
            "jeff_33".upper(),
            "jendl_40u".upper(),
            "irdff_ii".upper(),
            )
    library_ = library.lower()
    if library_ == "endfb_71":
        index = "https://www-nds.iaea.org/public/download-endf/ENDF-B-VII.1/tsl-index.htm"
    elif library_ == "endfb_80":
        index = "https://www-nds.iaea.org/public/download-endf/ENDF-B-VIII.0/tsl-index.htm"
    elif library_ == "jeff_33":
        index = "https://www-nds.iaea.org/public/download-endf/JEFF-3.3/tsl-index.htm"
    elif library_ == "jendl_40u":
        index = "https://www-nds.iaea.org/public/download-endf/JENDL-4.0u2-20160106/tsl-index.htm"
    elif library_ == "irdff_ii":
        index = "https://www-nds.iaea.org/public/download-endf/IRDFF-II/n-index.htm"
    else:
        raise ValueError(
            f"""library '{library}' is not available.
            Available libraries are: {available_libs}
            """
            )
    user_agent = 'Mozilla/5.0 (Windows; U; Windows NT 5.1; en-US; rv:1.9.0.7) Gecko/2009021910 Firefox/3.0.7'
    headers = {'User-Agent': user_agent, }
    request = Request(index, None, headers)
    response = urlopen(request)
    data = response.read().decode("utf-8")
    # Remove html style:
    data = data[data.find('<pre>')+5:data.find('</pre>')]
    data = re.sub(r'">tsl\S+', '', data)
    data = re.sub(r'<a href="tsl/', '', data)
    print(data.replace('MAT', 'KEY'))
    return


nsubs = {
    4: "decay",
    10: "neutron",
    11: "nfpy",
    10010: "proton",
    }


def get_endf6_file(library, kind, zam, to_file=False):
    """
    Given a library and a nuclide import the corresponding ENDF-6 nuclear
    data file directly from internet.

    Parameters
    ----------
    library : `str`
        nuclear data library. Available libraries are:
        for 'xs':
            * `'endfb_71'`
            * `'endfb_80'`
            * `'irdff_2'`
            * `'jeff_311'`
            * `'jeff_32'`
            * `'jeff_33'`
            * `'jendl_40u'`
            * `'tendl_2023'`
        for 'nfpy':
            * `'endfb_71'`
            * `'jeff_311'`
            * `'jeff_33'`
            * `'endfb_80'`
            * `'jendl_40u'`
        for 'decay':
            * `'endfb_71'`
            * `'jeff_311'`
            * `'jeff_33'`
            * `'endfb_80'`
        for 'tsl': (read the note)
            * `'endfb_71'`
            * `'jeff_33'`
            * `'endfb_80'`
            * `'jendl_40u'`
        for 'dxs':
            * `'jeff_33'`
            * `'proton'`
    kind : `str`
        nuclear data type:
            * 'xs' is a standard neutron-induced nuclear data file
            * 'nfpy' is a Neutron-Induced Fission Product Yields nuclear data
              file
            * 'decay' is a Radioactive Decay Data nuclear data file
            * 'dxs' is displacement cross section data file
            * 'tsl' is a Thermal Neutron Scattering Data file
    zam : `int` or 'all' or iterable
        zam = 'int' (individual nuclides) or iterable (group of nuclides)
            ZAM nuclide identifier $Z \\times 10000 + A \\times 10 + M$ where:
                * $Z$ is the charge number
                * $A$ is the mass number
                * $M$ is the metastate level (0=ground, 1=1st level)
        zam = 'all'
            We obtain the information of all the library. This option is not
            available for 'xs'

    Raises
    ------
    ValueError
        if library is not among available selection.

    ValueError
        if when you select 'xs', you select zam = 'all'

    Notes
    -----
    .. note:: In the `kind='tls'` option, instead of the zam, integers are used.
              If you need help, the `get_tsl_index` function contains all the
              necessary information for the correct choice of these integers.

    Returns
    -------
    `Endf6`
        `Endf6` object with ENDF-6 data for specified library and nuclide.

    Examples
    --------
    Import hydrogen file from JEFF-3.3.
    >>> tape = sandy.get_endf6_file("jeff_33", 'xs', 10010)
    >>> assert type(tape) is sandy.Endf6

    Import hydrogen file from JEFF-3.1.1.
    >>> tape = sandy.get_endf6_file("jeff_311", 'xs', 10010)
    >>> assert type(tape) is sandy.Endf6

    Import hydrogen file from ENDF/B-VII.1.
    >>> tape = sandy.get_endf6_file("endfb_71", 'xs', 10010)
    >>> assert type(tape) is sandy.Endf6

    Import hydrogen file from ENDF/B-VIII.0.
    >>> tape = sandy.get_endf6_file("endfb_80", 'xs', 10010)
    >>> assert type(tape) is sandy.Endf6

    Import hydrogen file from JENDL-4.0u
    >>> tape = sandy.get_endf6_file("jendl_40u", 'xs', 10010)
    >>> assert type(tape) is sandy.Endf6

    Import hydrogen file from TENDL-2023
    >>> tape = sandy.get_endf6_file("tendl_2023", 'xs', 10010)
    >>> assert type(tape) is sandy.Endf6

    Import Neutron-Induced Fission Product Yields for Th-227 from ENDF/B-VII.1.
    >>> tape = sandy.get_endf6_file("endfb_71", 'nfpy', 902270)
    >>> assert type(tape) is sandy.Endf6

    Import Neutron-Induced Fission Product Yields for Th-227 from ENDF/B-VIII.0
    >>> tape = sandy.get_endf6_file("endfb_80", 'nfpy', 902270)
    >>> assert type(tape) is sandy.Endf6

    Import Neutron-Induced Fission Product Yields for Th-227 from JENDL-4.0u
    >>> tape = sandy.get_endf6_file("jendl_40u", 'nfpy', 902270)
    >>> assert type(tape) is sandy.Endf6

    Import Neutron-Induced Fission Product Yields for Th-232 from JEFF-3.1.1
    >>> tape = sandy.get_endf6_file("jeff_311", 'nfpy', 902320)
    >>> assert type(tape) is sandy.Endf6

    Import Neutron-Induced Fission Product Yields for Th-232 from JEFF-3.3
    >>> tape = sandy.get_endf6_file("jeff_33", 'nfpy', 902320)
    >>> assert type(tape) is sandy.Endf6

    Import Radioactive Decay Data for H-1 from JEFF-3.1.1
    >>> tape = sandy.get_endf6_file("jeff_311", 'decay', 10010)
    >>> assert type(tape) is sandy.Endf6

    Import Radioactive Decay Data for H-1 from JEFF-3.3
    >>> tape = sandy.get_endf6_file("jeff_33", 'decay', 10010)
    >>> assert type(tape) is sandy.Endf6

    Import Radioactive Decay Data for H-1 from ENDF/B-VII.1.
    >>> tape = sandy.get_endf6_file("endfb_71", 'decay', 10010)
    >>> assert type(tape) is sandy.Endf6

    Import Radioactive Decay Data for H-1 from ENDF/B-VIII.0.
    >>> tape = sandy.get_endf6_file("endfb_80", 'decay', 10010)
    >>> assert type(tape) is sandy.Endf6

    Import all Neutron-Induced Fission Product Yields from ENDF/B-VII.1.
    >>> tape = sandy.get_endf6_file("endfb_71", 'nfpy', 'all')
    >>> assert type(tape) is sandy.Endf6

    Import a list of Decay Data for JEFF-3.3.
    >>> tape = sandy.get_endf6_file("jeff_33", 'decay', [380900, 551370, 541350])
    >>> assert type(tape) is sandy.Endf6

    Thermal Neutron Scattering Data from ENDF/B-VII.1.
    >>> tape = sandy.get_endf6_file("endfb_71", 'tsl', [1, 2, 3])
    >>> assert type(tape) is sandy.Endf6

    Thermal Neutron Scattering Data from ENDF/B-VIII.0.
    >>> tape = sandy.get_endf6_file("endfb_80", 'tsl', [1, 2, 3])
    >>> assert type(tape) is sandy.Endf6

    Thermal Neutron Scattering Data from JEFF-3.3.
    >>> tape = sandy.get_endf6_file("jeff_33", 'tsl', [1, 2, 3])
    >>> assert type(tape) is sandy.Endf6

    Thermal Neutron Scattering Data from JENDL-4.0u
    >>> tape = sandy.get_endf6_file("jendl_40u", 'tsl', [1, 2, 3])
    >>> assert type(tape) is sandy.Endf6
    
    Import natural Fe for IRDFF-II
    >>> tape = sandy.get_endf6_file("irdff_2", "xs", 260000)
    >>> assert type(tape) is sandy.Endf6
    """
    foo_get = Endf6.from_zipurl
    foo_read = Endf6.read_zipurl
    if kind == 'xs':
        available_libs = (
            "jeff_311".upper(),
            "jeff_32".upper(),
            "jeff_33".upper(),
            "endfb_71".upper(),
            "endfb_80".upper(),
            "jendl_40u".upper(),
            "irdff_2".upper(),
            )
        library_ = library.lower()
        if library_ == "jeff_40t0":  # not allowed anymore since NEA change website
            url = URL_N_JEFF_40T0_NEA
            files = N_FILES_JEFF_40T0_NEA
            foo_read = Endf6.read_url
            foo_get = Endf6.from_url
        elif library_ == "jeff_311":
            url = URL_N_JEFF_311_IAEA
            files = N_FILES_JEFF_311_IAEA
        elif library_ == "jeff_33":
            url = URL_N_JEFF_33_IAEA
            files = N_FILES_JEFF_33_IAEA
        elif library_ == "jeff_32":
            url = URL_N_JEFF_32_NEA
            files = N_FILES_JEFF_32_NEA
            foo_read = Endf6.read_url
            foo_get = Endf6.from_url
        elif library_ == "endfb_71":
            url = URL_N_ENDFB_71_IAEA
            files = N_FILES_ENDFB_71_IAEA
        elif library_ == "endfb_80":
            url = URL_N_ENDFB_80_IAEA
            files = N_FILES_ENDFB_80_IAEA
        elif library_ == "jendl_40u":
            url = URL_N_JENDL_40U_IAEA
            files = N_FILES_JENDL_40U_IAEA
        elif library_ == "tendl_2023":
            url = URL_N_TENDL_2023_PSI
            files = N_FILES_TENDL_2023_PSI
            foo_read = Endf6.read_url
            foo_get = Endf6.from_url
        elif library_ == "irdff_2":
            url = URL_N_IRDFF_2_IAEA
            files = N_FILES_IRDFF_2_IAEA
        else:
            raise ValueError(
                f"""library '{library}' is not available.
                Available libraries are: {available_libs}
                """
                )
    elif kind == 'dxs':
        available_libs = (
            "jeff_33".upper(),
            "proton".upper(),
            )
        library_ = library.lower()
        if library_ == "jeff_33":
            url = URL_DXS_JEFF_33_IAEA
            files = DXS_FILES_JEFF_33_IAEA
            foo_read = Endf6.read_url
            foo_get = Endf6.from_url
        elif library_ == "proton":
            url = URL_DXS_PROTON_IAEA
            files = DXS_FILES_PROTON_IAEA
            foo_read = Endf6.read_url
            foo_get = Endf6.from_url
        else:
            raise ValueError(
                f"""library '{library}' is not available.
                Available libraries are: {available_libs}
                """
                )
    elif kind == 'nfpy':
        available_libs = (
            "endfb_71".upper(),
            "endfb_80".upper(),
            "jendl_40u".upper(),
            "jeff_311".upper(),
            "jeff_33".upper(),
            )
        library_ = library.lower()
        if library_ == "endfb_71":
            url = URL_NFPY_ENDFB_71_IAEA
            files = NFPY_FILES_ENDFB_71_IAEA
        elif library_ == "endfb_80":
            url = URL_NFPY_ENDFB_80_IAEA
            files = NFPY_FILES_ENDFB_80_IAEA
        elif library_ == "jendl_40u":
            url = URL_NFPY_JENDL_40U_IAEA
            files = NFPY_FILES_JENDL_40U_IAEA
        elif library_ == "jeff_311":
            url = URL_NFPY_JEFF_311_IAEA
            files = NFPY_FILES_JEFF_311_IAEA
        elif library_ == "jeff_33":
            url = URL_NFPY_JEFF_33_IAEA
            files = NFPY_FILES_JEFF_33_IAEA
        else:
            raise ValueError(
                f"""library '{library}' is not available.
                Available libraries are: {available_libs}
                """
                    )
    elif kind == 'decay':
        available_libs = (
            "endfb_71".upper(),
            "endfb_80".upper(),
            "jeff_311".upper(),
            "jeff_33".upper(),
            "jendl_40u".upper(),
            )
        library_ = library.lower()
        if library_ == "endfb_71":
            url = URL_DECAY_ENDFB_71_IAEA
            files = DECAY_FILES_ENDFB_71_IAEA
        elif library_ == "endfb_80":
            url = URL_DECAY_ENDFB_80_IAEA
            files = DECAY_FILES_ENDFB_80_IAEA
        elif library_ == "jeff_311":
            url = URL_DECAY_JEFF_311_IAEA
            files = DECAY_FILES_JEFF_311_IAEA
        elif library_ == "jeff_33":
            url = URL_DECAY_JEFF_33_IAEA
            files = DECAY_FILES_JEFF_33_IAEA
        elif library_ == "jendl_40u":
            # it will fail for indivdual files, but it works for the ntire library
            pass
        else:
            raise ValueError(
                f"""library '{library}' is not available.
                Available libraries are: {available_libs}
                """
                    )
    elif kind == 'tsl':
        available_libs = (
            "endfb_71".upper(),
            "endfb_80".upper(),
            "jeff_33".upper(),
            "jendl_40u".upper(),
            )
        library_ = library.lower()
        if library_ == "endfb_71":
            url = URL_TSL_ENDFB_71_IAEA
            files = TSL_FILES_ENDFB_71_IAEA
        elif library_ == "endfb_80":
            url = URL_TSL_ENDFB_80_IAEA
            files = TSL_FILES_ENDFB_80_IAEA
        elif library_ == "jeff_33":
            url = URL_TSL_JEFF_33_IAEA
            files = TSL_FILES_JEFF_33_IAEA
        elif library_ == "jendl_40u":
            url = URL_TSL_JENDL_40U_IAEA
            files = TSL_FILES_JENDL_40U_IAEA
        else:
            raise ValueError(
                f"""library '{library}' is not available.
                Available libraries are: {available_libs}
                """
                    )
    else:
        raise ValueError(f"option 'kind={kind}' is not supported")

    if str(zam).lower() == 'all':
        if kind.lower() == 'xs' or kind.lower() == 'dxs':
            raise ValueError("'all' option is not available for xs")
        elif kind.lower() == "decay":
            # read from local files with all data, otherwise it's too slow
            file = join(dirname(sandy.__file__), 'appendix', 'decay_data', f"rdd.{library_}")
            tape = sandy.Endf6.from_file(file)
        else:  # basically only for fission yield
            text = "".join([foo_read(name, url) for name in files.values()])
            tape = Endf6.from_text(text)
    else:
        if hasattr(zam, "__len__"):
            tapes = map(lambda x: foo_get(files[x], url), zam)
            tape = reduce(lambda x, y: x.add_sections(y.data), tapes)
        else:
            tape = foo_get(files[zam], url)

    if to_file:
        basename = sandy.zam.zam2nuclide(zam, atomic_number=True, sep="-")
        filename = f"{basename}.{library_}"
        logging.info(f"writing nuclear data to file '{filename}'")
        tape.to_file(filename)
    return tape



class _FormattedFile():
    """
    Base class to store ENDF-6 content grouped by `(MAT, MF, MT)`

    Attributes
    ----------
    data

    keys

    kind : `str`
        Kind of ENDF-6 formatted file (`'endf6'`, `'pendf'`, `'gendf'`,
        `'errorr'`) .

    mat : `int`
        MAT number.
    mf : `int`
        MF number.
    mt : `int`
        MT number

    Methods
    -------
    add_sections
        Add text section for given `(MAT, MF, MT)`.
    filter_by
        Filter dataframe based on `(MAT, MF, MT)` lists.
    from_file
        Create dataframe by reading a ENDF-6-formatted file.
    from_text
        Create dataframe from endf6 text in string.
    to_series
        Covert content into `pandas.Series`.
    to_file
        Given a filename write the content of the instance to disk in
        ASCII format.

    Notes
    -----
    This class supports ENDF-6 content from ENDF-6 files, ERRORR files and
    GROUPR files.
    """

    def __repr__(self):
        return self.to_series().__repr__()

    def __init__(self, data, file=None):
        self.data = data
        self.file = file

    @property
    def data(self):
        return self._data

    @data.setter
    def data(self, data):
        if not isinstance(data, dict):
            raise sandy.Error("'data' is not a 'dict'")
        self._data = data

    @property
    def keys(self):
        """
        List of keys `(MAT, MF, MT)` used to identify each tape section.

        Returns
        -------
        `list`
            list of tuples of type `(MAT, MF, MT)` for each section found in
            tape.
        """
        return list(self.data.keys())

    @property
    def _keys(self):
        mat, mf, mt = zip(*self.data.keys())
        return {"MAT": mat, "MF": mf, "MT": mt}

    @property
    def mat(self):
        return sorted(set(self._keys["MAT"]))

    @property
    def mf(self):
        return sorted(set(self._keys["MF"]))

    @property
    def mt(self):
        return sorted(set(self._keys["MT"]))

    def to_series(self, **kwargs):
        series = pd.Series(self.data, **kwargs).sort_index(ascending=True)
        series.index.names = ["MAT", "MF", "MT"]
        return series

    @property
    def is_empty(self):
        return False if self.data else True

    @property
    def kind(self):
        """
        Kind of ENDF-6 formatted file (`'endf6'`, `'pendf'`, `'gendf'`,
        `'errorr'`) .

        Returns
        -------
        `str`
            kind of ENDF-6 formatted file

        Examples
        --------
        >>> assert sandy.get_endf6_file("jeff_33", "decay", 10010).kind == "endf6"
        >>> assert sandy.get_endf6_file("jeff_33", "nfpy", 922350).kind == "endf6"
        >>> assert sandy.get_endf6_file("jeff_33", "xs", 10010).kind == "endf6"
        >>> assert sandy.get_endf6_file("jeff_33", "xs", 10010).get_pendf(err=1).kind == "pendf"
        >>> assert sandy.get_endf6_file("jeff_33", "xs", 10010).get_gendf(err=1).kind == "gendf"
        >>> outs = sandy.get_endf6_file("jeff_33", "xs", 942410).get_errorr(err=1, errorr_kws=dict(mt=18))
        >>> assert outs["errorr31"].kind == "errorr"
        >>> assert outs["errorr33"].kind == "errorr"
        >>> assert outs["errorr34"].kind == "errorr"
        >>> assert outs["errorr35"].kind == "errorr"
        """
        if len(self.mat) > 1:
            msg = "Attribute 'kind' does not work if more than 1 MAT number is"
            "found"
            raise AttributeError(msg)
        mat = self.mat[0]
        text = self.data[(mat, 1, 451)]
        lrp = int(text[22:33])
        nlib = int(text[44:55])
        if nlib == -11 or nlib == -12:
            kind = "errorr"
        elif nlib == -1:
            kind = "gendf"
        else:
            if lrp == 2:
                kind = "pendf"
            elif lrp in [-1, 0, 1]:
                # -1 for decay and nfpy
                # 0 for endf6
                kind = "endf6"
            else:
                kind = "unkwown"
        return kind

    @classmethod
    def from_url(cls, filename, rooturl):
        """
        Given a filename and the url where the file is located,
        extract the ENDF6 data from the file into a sandy.Endf6
        instance.

        Parameters
        ----------
        filename: 'str'
            The filename without path of the zip file to read
        rooturl: 'str'
            The url direction to extract the zip files

        Returns
        -------
        `Endf6`
            `Endf6` object with ENDF-6 data for specified library and nuclide.

        Examples
        --------

        >>> filename = "n-1-H-001.jeff32"
        >>> rooturl = "https://www.oecd-nea.org/dbforms/data/eva/evatapes/jeff_32/"
        >>> sandy.Endf6.from_url(filename, rooturl)
        MAT  MF  MT
        125  1   451     1.001000+3 9.991673-1          0          0  ...
            2   151     1.001000+3 9.991673-1          0          0  ...
            3   1       1.001000+3 9.991673-1          0          0  ...
                2       1.001000+3 9.991673-1          0          0  ...
                102     1.001000+3 9.991673-1          0          0  ...
            4   2       1.001000+3 9.991673-1          0          1  ...
            6   102     1.001000+3 9.991673-1          0          2  ...
            33  1       1.001000+3 9.991673-1          0          0  ...
                2       1.001000+3 9.991673-1          0          0  ...
                102     1.001000+3 9.991673-1          0          0  ...
        dtype: object
        """
        text = cls.read_url(filename, rooturl)
        tape = cls.from_text(text)
        return tape

    @staticmethod
    def read_url(filename, rooturl):
        """
        Given a filename and the url where the file is located,
        extract the ENDF6 data from the file into a string.

        Parameters
        ----------
        filename: 'str'
            The filename without path of the zip file to read
        rooturl: 'str'
            The url direction to extract the zip files

        Returns
        -------
        `str`
            All the endf6 data in a 'str'

        Examples
        --------
        Removed because website stopped working
        #>>> filename = "n-1-H-001.jeff32"
        #>>> rooturl = "https://www.oecd-nea.org/dbforms/data/eva/evatapes/jeff_32/"
        #>>> file = sandy.Endf6.read_url(filename, rooturl)
        #>>> print(file[0:890])
        JEFF-3.2 Release - Neutron File March 2014                             0  0    0
        1.001000+3 9.991673-1          0          0          2          5 125 1451    1
        0.000000+0 0.000000+0          0          0          0          6 125 1451    2
        1.000000+0 2.000000+7          1          0         10         32 125 1451    3
        0.000000+0 0.000000+0          0          0         87         10 125 1451    4
        1-H -  1 LANL       EVAL-OCT05 G.M.Hale                          125 1451    5
                            DIST-DEC06                       20111222    125 1451    6
        ----JEFF32            MATERIAL  125                                125 1451    7
        -----INCIDENT NEUTRON DATA                                         125 1451    8
        ------ENDF-6 FORMAT                                                125 1451    9
        *****************************  JEFF-3.2    *********************** 125 1451   10
        """
        url = f"{rooturl}/{filename}"
        # set a known browser user agent to ensure access
        req = Request(url, headers={'User-Agent': 'Mozilla/5.0'},)
        with urlopen(req) as f:
            text = f.read().decode('utf-8')
        return text

    @classmethod
    def from_zipurl(cls, filename, rooturl):
        """
        Given a filename and the url where the file is located (in
        zipped format), extract the ENDF6 data from the file into a
        sandy.Endf6 instance.

        Parameters
        ----------
        filename: 'str'
            The filename without path of the zip file to read
        rooturl: 'str'
            The url direction to extract the zip files

        Returns
        -------
        `Endf6`
            `Endf6` object with ENDF-6 data for specified library and nuclide.

        Examples
        --------

        >>> filename = "decay_1907_57-La-149.dat"
        >>> rooturl = "https://www-nds.iaea.org/public/download-endf/ENDF-B-VII.1/decay/"
        >>> sandy.Endf6.from_zipurl(filename, rooturl)
        MAT   MF  MT
        1907  1   451     5.714900+4 1.476553+2         -1          0  ...
            8   457     5.714900+4 1.476553+2          0          0  ...
        dtype: object
        """
        text = cls.read_zipurl(filename, rooturl)
        tape = cls.from_text(text)
        return tape

    @staticmethod
    def read_zipurl(filename, rooturl, mode=2):
        """
        Given a filename and the url where the file is located (in
        zipped format), extract the ENDF6 data from the file into
        a string.

        Parameters
        ----------
        filename: 'str'
            The filename without path of the zip file to read
        rooturl: 'str'
            The url direction to extract the zip files

        Returns
        -------
        `str`
            All the endf6 data in a 'str'

        Examples
        --------

        >>> filename = "decay_1907_57-La-149.dat"
        >>> rooturl = "https://www-nds.iaea.org/public/download-endf/ENDF-B-VII.1/decay/"
        >>> file = sandy.Endf6.read_zipurl(filename, rooturl)
        >>> print(file[0:971])
        Retrieved by E4-util: 2012/01/16,13:45:44                            1 0  0    0
        5.714900+4 1.476553+2         -1          0          0          11907 1451    1
        0.000000+0 1.000000+0          0          0          0          61907 1451    2
        0.000000+0 0.000000+0          1          0          4          71907 1451    3
        0.000000+0 0.000000+0          0          0         27          21907 1451    4
        57-La-149  BNL        EVAL-AUG11 Conv. from CGM                  1907 1451    5
        /ENSDF/                                               20111222   1907 1451    6
        ----ENDF/B-VII.1      Material 1907                               1907 1451    7
        -----RADIOACTIVE DECAY DATA                                       1907 1451    8
        ------ENDF-6 FORMAT                                               1907 1451    9
        *********************** Begin Description *********************** 1907 1451   10
        **         ENDF/B-VII.1 RADIOACTIVE DECAY DATA FILE            ** 1907 1451   11
        """
        rootname = os.path.splitext(filename)[0]
        zipurl = f"{rooturl}/{rootname}.zip"

        # set a known browser user agent to ensure access

        if mode == 1:
            req = Request(zipurl, headers={'User-Agent': 'Mozilla/5.0'})
            with urlopen(req) as zipresp:
                with ZipFile(io.BytesIO(zipresp.read())) as zfile:
                    with TemporaryDirectory() as td:
                        zfile.extract(filename, path=td)
                        tmpfile = join(td, filename)
                        with open(tmpfile, "r") as f:
                            text = f.read()

        elif mode == 2:  # as of 7/10/2024 mode 1 is deprecated because of changes on the IAEA websites
            class AppURLopener(urllib.request.FancyURLopener):
                version = "Mozilla/5.0"

            # Instantiate the custom opener
            opener = AppURLopener()

            # Open the URL using the custom opener
            response = opener.open(zipurl)
                
            # Read the response content into a bytes object
            zip_data = response.read()

            # Ensure the response is closed
            response.close()

            # Use the zipfile module to read the zip file from the bytes object
            with ZipFile(io.BytesIO(zip_data)) as zfile:
                with TemporaryDirectory() as td:
                    zfile.extract(filename, path=td)
                    tmpfile = join(td, filename)
                    with open(tmpfile, "r") as f:
                        text = f.read()

        return text

    @classmethod
    def from_file(cls, file):
        """
        Create dataframe by reading a file.

        Parameters
        ----------
        file : `str`
            filename

        Returns
        -------
        `sandy.formats.endf6.BaseFile` or derived instance
            Dataframe containing ENDF6 data grouped by MAT/MF/MT

        Examples
        --------
        Read hydrogen tape from endf-6 formatted file.
        >>> file = "h1.endf"
        >>> sandy.get_endf6_file("jeff_33", "xs", 10010).to_file(file)
        >>> _FormattedFile.from_file(file)
        MAT  MF  MT
        125  1   451     1.001000+3 9.991673-1          0          0  ...
             2   151     1.001000+3 9.991673-1          0          0  ...
             3   1       1.001000+3 9.991673-1          0          0  ...
                 2       1.001000+3 9.991673-1          0          0  ...
                 102     1.001000+3 9.991673-1          0          0  ...
             4   2       1.001000+3 9.991673-1          0          1  ...
             6   102     1.001000+3 9.991673-1          0          2  ...
             33  1       1.001000+3 9.991673-1          0          0  ...
                 2       1.001000+3 9.991673-1          0          0  ...
                 102     1.001000+3 9.991673-1          0          0  ...
        dtype: object

        Read hydrogen tape from text stream.
        >>> stream = io.StringIO(open(file).read())
        >>> _FormattedFile.from_file(stream)
        MAT  MF  MT
        125  1   451     1.001000+3 9.991673-1          0          0  ...
             2   151     1.001000+3 9.991673-1          0          0  ...
             3   1       1.001000+3 9.991673-1          0          0  ...
                 2       1.001000+3 9.991673-1          0          0  ...
                 102     1.001000+3 9.991673-1          0          0  ...
             4   2       1.001000+3 9.991673-1          0          1  ...
             6   102     1.001000+3 9.991673-1          0          2  ...
             33  1       1.001000+3 9.991673-1          0          0  ...
                 2       1.001000+3 9.991673-1          0          0  ...
                 102     1.001000+3 9.991673-1          0          0  ...
        dtype: object
        """
        if isinstance(file, io.StringIO):
            text = file.read()
        else:
            with open(file) as f:
                text = f.read()
        return cls.from_text(text)

    @classmethod
    def from_text(cls, text):
        """
        Create dataframe from endf6 text in string.

        Parameters
        ----------
        text : `str`
            string containing the evaluated data

        Returns
        -------
        `sandy.formats.endf6.BaseFile` or derived instance
            Dataframe containing ENDF6 data grouped by MAT/MF/MT

        Examples
        --------
        Read hydrogen tape from text.
        >>> file = "h1.endf"
        >>> sandy.get_endf6_file("jeff_33", "xs", 10010).to_file(file)
        >>> text = open(file).read()
        >>> _FormattedFile.from_text(text)
        MAT  MF  MT
        125  1   451     1.001000+3 9.991673-1          0          0  ...
             2   151     1.001000+3 9.991673-1          0          0  ...
             3   1       1.001000+3 9.991673-1          0          0  ...
                 2       1.001000+3 9.991673-1          0          0  ...
                 102     1.001000+3 9.991673-1          0          0  ...
             4   2       1.001000+3 9.991673-1          0          1  ...
             6   102     1.001000+3 9.991673-1          0          2  ...
             33  1       1.001000+3 9.991673-1          0          0  ...
                 2       1.001000+3 9.991673-1          0          0  ...
                 102     1.001000+3 9.991673-1          0          0  ...
        dtype: object
        """
        df = pd.read_fwf(
            io.StringIO(text),
            widths=[66, 4, 2, 3],
            names=["TEXT", "MAT", "MF", "MT"],
            dtype={"TEXT": str, "MAT": str, "MF": int, "MT": int},
            na_filter=False,  # speeds up and does not add NaN in empty lines
            # Do not use TEXT because  the parser does not preserve the
            # whitespaces
            usecols=("MAT", "MF", "MT"),
            )
        # use splitlines instead of readlines to remove "\n"
        df["TEXT"] = text.splitlines()
        #
        title = df["TEXT"].iloc[0]
        title_mat = df["MAT"].iloc[0]
        try:
            int(title_mat)
        except ValueError:
            logging.warning(f"wrong MAT number in the file title\n'{title}'")
            df = df.iloc[1:].reset_index(drop=True)
        finally:
            df["MAT"] = df["MAT"].astype(int)
        condition = (df.MT > 0) & (df.MF > 0) & (df.MAT > 0)
        data = df[condition].groupby(["MAT", "MF", "MT"])\
                            .agg({"TEXT": "\n".join})\
                            .TEXT\
                            .to_dict()
        return cls(data)

    def _get_section_df(self, mat, mf, mt):
        """
        
        Examples
        --------
        Check if we can read Pu240 file from JEFF-3.1.1, where a "?" is found in the header.
        >>> tape = sandy.get_endf6_file("jeff_311", "xs", 942400)
        >>> assert "?" in tape.data[(9440, 1, 451)]
        >>> out = tape._get_section_df(9440, 1, 451)
        
        Let's make it fail.
        >>> import pytest
        >>> text = " 94-Pu-240 BRC,CAD    EVAL-JUL04 Bouland Derrien Morillon R?@$¤n  9440 1451    5"
        >>> with pytest.raises(ValueError) as exc_info:
        ...    sandy.Endf6.from_text(text)._get_section_df(9440, 1, 451)
        """
        text = self.data[(mat, mf, mt)]
        delimiters = ["?", "@", "$", "¤"]
        found = False
        for delimiter in delimiters:
            found = delimiter not in text
            if found:
                break
            logging.info(f"Could not parse Endf6 as DataFrame using '{delimiter}' delimiter")
        if not found:
            raise ValueError("Could not find suitable delimiter to parse Endf6 file.")

        def foo(x):
            return sandy.shared.add_delimiter_every_n_characters(
                x[:66],
                11,
                delimiter=delimiter,
            )

        newtext = "\n".join(map(foo, text.splitlines())).replace('"', '*')
        df = pd.read_csv(
            io.StringIO(sandy.shared.add_exp_in_endf6_text(newtext)),
            delimiter=delimiter,
            na_filter=True,
            names=["C1", "C2", "L1", "L2", "N1", "N2"],
        )
        return df

    def add_section(self, mat, mf, mt, text):
        """
        Given MAT, MF and MT add/replace the corresponding section in the
        `Endf6.data`.

        Parameters
        ----------
        mat : `int`
            MAT number
        mf : `int`
            MF number
        mt : `int`
            MT number

        Returns
        -------
        `sandy._FormattedFile` or derived instance
            object with new section

        Examples
        --------
        >>> tape = sandy.Endf6({(9437, 3, 102) : "lorem ipsum"})
        >>> tape.add_section(9999, 1, 1, "dolor sit amet")
        MAT   MF  MT
        9437  3   102       lorem ipsum
        9999  1   1      dolor sit amet
        dtype: object
        """
        d = self.data.copy()
        key = (mat, mf, mt)
        d[key] = text
        return self.__class__(d)

    def add_sections(self, sections):
        d = self.data.copy()
        for (mat, mf, mt), text in sections.items():
            key = (mat, mf, mt)
            d[key] = text
        return self.__class__(d)

    def delete_section(self, mat, mf, mt, raise_error=True):
        """
        Given MAT, MF and MT delete the corresponding section from the
        `Endf6.data`.

        Parameters
        ----------
        mat : `int`
            MAT number
        mf : `int`
            MF number
        mt : `int`
            MT number

        Returns
        -------
        `sandy._FormattedFile` or derived instance
            object without given section

        Examples
        --------
        Delete capture cross section from hydrogen file.
        >>> tape = sandy.get_endf6_file("jeff_33", "xs", 10010)
        >>> new = tape.delete_section(125, 3, 102)
        >>> new
        MAT  MF  MT
        125  1   451     1.001000+3 9.991673-1          0          0  ...
             2   151     1.001000+3 9.991673-1          0          0  ...
             3   1       1.001000+3 9.991673-1          0          0  ...
                 2       1.001000+3 9.991673-1          0          0  ...
             4   2       1.001000+3 9.991673-1          0          1  ...
             6   102     1.001000+3 9.991673-1          0          2  ...
             33  1       1.001000+3 9.991673-1          0          0  ...
                 2       1.001000+3 9.991673-1          0          0  ...
                 102     1.001000+3 9.991673-1          0          0  ...
        dtype: object
        """
        d = self.data.copy()
        key = (mat, mf, mt)
        if key not in d and raise_error is False:
            pass
        else:
            del d[key]
        return self.__class__(d)

    def merge(self, *iterable):
        """
        Given a single `sandy.Endf6` object or an iterable of `sandy.Endf6`
        objects as keyword arguments, add their sections to a copy of the
        `self` instance and return a new `sandy.Endf6` object.
        The new `sandy.Endf6` object contains all MAT/MF/MT sections in `self`
        and in the passed arguments.

        Parameters
        ----------
        iterable : `sandy.Endf6` or iterable of `sandy.Endf6` objects
            The ENDF6 files that will be merged to `self`.

        Returns
        -------
        merged : :func:`_FormattedFile`
            a ENDF6 file containing the MAT/MF/MT sections of `self` and of
            the passed ENDF6 files.

        Notes
        -----
        .. note:: if any section (MAT/MF/MT) already present in the orginal
                  ENDF6 tape also appears in any tape that is being merged,
                  then the original ENDF6 section will be overwritten.

        Examples
        --------
        Merge two files.
        >>> h1 = sandy.get_endf6_file("jeff_33", 'xs', 10010)
        >>> h2 = sandy.get_endf6_file("endfb_71", 'xs', 10020)
        >>> h = h1.merge(h2)
        >>> assert h.to_series()[h1.to_series().index].equals(h1.to_series())
        >>> assert h.to_series()[h2.to_series().index].equals(h2.to_series())

        Merge three files from different libraries.
        >>> h3 = sandy.get_endf6_file("endfb_71", 'xs', 10030)
        >>> h_ = h1.merge(h2, h3).to_series()
        >>> h__ = h.merge(h3).to_series()
        >>> h___ = h1.merge(h2).merge(h3).to_series()
        >>> assert h_.equals(h__) and h_.equals(h___)

        Merge two evaluations for the same nuclide.
        >>> bi_71 = sandy.get_endf6_file("endfb_71", 'xs', 832090)
        >>> bi_33 = sandy.get_endf6_file("jeff_33", 'xs', 832090)
        >>> bi = bi_71.merge(bi_33)
        >>> assert not bi.to_series()[bi_71.to_series().index].equals(bi_71.to_series())
        >>> assert bi.to_series()[bi_33.to_series().index].equals(bi_33.to_series())
        >>> bi = bi_33.merge(bi_71)
        >>> assert bi.to_series()[bi_71.to_series().index].equals(bi_71.to_series())
        >>> assert not bi.to_series()[bi_33.to_series().index].equals(bi_33.to_series())
        """
        tape = reduce(lambda x, y: x.add_sections(y.data), iterable)
        merged = self.add_sections(tape.data)
        return merged

    def filter_by(self,
                  listmat=range(1, 10000),
                  listmf=range(1, 10000),
                  listmt=range(1, 10000)):
        """Filter dataframe based on MAT, MF, MT lists.

        Parameters
        ----------
        listmat : `list` or `None`
            list of requested MAT values (default is `None`: use all MAT)
        listmf : `list` or `None`
            list of requested MF values (default is `None`: use all MF)
        listmt : `list` or `None`
            list of requested MT values (default is `None`: use all MT)

        Returns
        -------
        `sandy._FormattedFile` or derived instance
            Copy of the original instance with filtered MAT, MF and MT sections
        """
        df = self.to_series().to_frame()
        d = df.query("MAT in @listmat and MF in @listmf and MT in @listmt").squeeze(axis=1).to_dict()
        return self.__class__(d)

    def get_value(self, mat, mf, mt, line_number, pos):
        return self._get_section_df(mat, mf, mt)[pos.upper()] \
                   .iloc[line_number - 1]

    def change_value(self, val, mat, mf, mt, line_number, pos, inplace=False,
                     dtype=float):
        items = ("C1", "C2", "L1", "L2", "N1", "N2")
        positions = {y: x for x, y in enumerate(items)}
        step = positions[pos]
        length = 81
        ibeg = length * (line_number - 1) + 11 * step
        iend = length * (line_number - 1) + 11 * (step + 1)
        text = self.data[(mat, mf, mt)]
        new = sandy.write_int(val) if dtype is int else sandy.write_float(val)
        new_text = "".join((text[:ibeg], new, text[iend:]))
        new_tape = self.add_section(mat, mf, mt, new_text, inplace=False)
        if inplace:
            self.data = new_tape.data
        else:
            return new_tape
        print(new_text)

    def write_string(self, title="", tpid=True, fend=True):
        """
        Write `_FormattedFile.data` content to string according to the ENDF-6
        file rules.

        Parameters
        ----------
        title : `str`, optional, default is an empty string
            first line of the file
        tpid : `bool`, optional, defult is `True`
            write TPID line.
            A TPID line is a text line at the beginning of a file,
            ending with `'   1 0  0    0'`.
            
        fend : `bool`, optional, defult is `True`
            write END-OF-FILE line.
            A FEND line is a text line at the end of a file,
            ending with `'  -1 0  0    0'`.

        Returns
        -------
        `str`
            string containing the ENDF-6 information stored in this instance.

        Notes
        -----
        ..note:: no modification is implemented to the actual content of
                 the `Endf6.data` object.

        Examples
        --------
        >>> string = sandy.get_endf6_file("jeff_33", "xs", 10010).write_string()
        >>> print(string[:81 * 4 - 1])
                                                                             1 0  0    0
         1.001000+3 9.991673-1          0          0          2          5 125 1451    1
         0.000000+0 0.000000+0          0          0          0          6 125 1451    2
         1.000000+0 2.000000+7          3          0         10          3 125 1451    3

        if no modification is applied to the `_FormattedFile` content, the
        `write_string` returns an output identical to the file ASCII content.

        Test with `sandy.Errorr` object and title option:
        >>> endf6 = sandy.get_endf6_file("jeff_33", "xs", 10010)
        >>> err = endf6.get_errorr(ek=[1e-2, 1e1, 2e7], err=1)["errorr33"]
        >>> err.to_file("out.err", title="H with ERRORR")
        >>> err_2 = sandy.Errorr.from_file("out.err")
        >>> os.remove("out.err")
        >>> assert err_2.data[(125, 1, 451)] == err.data[(125, 1, 451)]
        >>> assert err_2.data[(125, 3, 102)] == err.data[(125, 3, 102)]
        >>> assert err_2.data[(125, 33, 102)] == err.data[(125, 33, 102)]

        ..note:: differences might appear from the way zeros were handled at
                 the end of ENDF-6 section, or if a different fiel title is
                 given

        How to use keyword `tpid`.
        >>> first = endf6.write_string(tpid=False).splitlines()[0]
        >>> first_tpid = endf6.write_string(tpid=True).splitlines()[0]
        >>> assert first_tpid != first 
        >>> assert " "*66 + "   1 0  0    0" == first_tpid        
        >>> assert endf6.write_string(tpid=False)[0] == endf6.write_string(tpid=True)[0] == ' '
        
        How to use keyword `fend`.
        >>> last = endf6.write_string(fend=False).splitlines()[-1]
        >>> last_fend = endf6.write_string(fend=True).splitlines()[-1]
        >>> assert last_fend != last 
        >>> assert " "*66 + "  -1 0  0    0" == last_fend        
        >>> assert endf6.write_string(fend=False)[-1] == endf6.write_string(fend=True)[-1] == '0'
        
        Check that there is no line concatenation.
        >>> tape = sandy.get_endf6_file("jeff_33", "decay", [10010, 10020])
        >>> assert all(x==80 for x in map(len, tape.write_string().splitlines()))
        """
        string = ""

        # Write title
        if tpid:
            string += sandy.write_line(title, 1, 0, 0, 0)
            string += "\n"

        for mat, dfmat in self.to_series().groupby('MAT', sort=True):
            for mf, dfmf in dfmat.groupby('MF', sort=True):
                for mt, text in dfmf.groupby('MT', sort=True):
                    string += text.squeeze()\
                                  .encode('ascii', 'replace')\
                                  .decode('ascii')
                    string += "\n"
                    string += sandy.write_line("", mat, mf, 0, 99999)
                    string += "\n"
                string += sandy.write_line("", mat, 0, 0, 0)
                string += "\n"
            string += sandy.write_line("", 0, 0, 0, 0)
            string += "\n"

        # Write end-of-file
        if fend:
            string += sandy.write_line("", -1, 0, 0, 0)
        else:
            # remove laast newline
            string = string[:-1]

        return string

    def to_file(self, filename, mode="w", **kwargs):
        """
        Given a filename write the content of a `_FormattedFile` instance to
        disk in ASCII format.

        Parameters
        ----------
        filename : `str`
            The name of the file.
        mode : `str`, optional
            Mode while opening a file. The default is "w".

        Parameters for `write_string`
        -----------------------------
        title : `str`, optional, default is an empty string
            first line of the file

        Returns
        -------
        None.

        """
        text = self.write_string(**kwargs)
        with open(filename, mode) as f:
            f.write(text)


class Endf6(_FormattedFile):
    """
    Container for ENDF-6 file text grouped by MAT, MF and MT numbers.

    Methods
    -------
    get_ace
        Process `Endf6` instance into an ACE file using NJOY.
    get_pendf
        Process `Endf6` instance into a PENDF file using NJOY.
    get_errorr
        Process `Endf6` instance into a Errorr file using NJOY.
    get_id
        Extract ID for a given MAT for a ENDF-6 file.
    get_nsub
        Determine ENDF-6 sub-library type.
    get_records
        Extract tabulated MAT, MF and MT numbers.
    read_section
        Parse MAT/MF/MT section.
    update_intro
        Update MF1/MT451.
    """

    def update_intro(self, **kwargs):
        """
        Method to update MF1/MT451 of each MAT based on the file content
        (concistency is enforced) and user-given keyword arguments.
        
        Parameters
        ----------
        **kwargs : `dict`
            dictionary of elements to be modified in section MF1/MT451 (it
            applies to all MAT numbers).

        Returns
        -------
        :func:`~sandy.Endf6`
            :func:`~sandy.Endf6` with updated MF1/MT451.

        
        Examples
        --------
        Check how many lines of description and how many sections are recorded
        in a file.
        >>> tape = sandy.get_endf6_file("jeff_33", "xs", 10010)
        >>> intro = tape.read_section(125, 1, 451)
        >>> assert len(intro["DESCRIPTION"]) == 87
        >>> assert len(intro["SECTIONS"]) == 10

        By removing sections in the `Endf6` instance, the recorded number of
        sections does not change.
        >>> tape2 = tape.delete_section(125, 33, 1).delete_section(125, 33, 2)
        >>> intro = tape2.read_section(125, 1, 451)
        >>> assert len(intro["DESCRIPTION"]) == 87
        >>> assert len(intro["SECTIONS"]) == 10

        Running `updated intro` updates the recorded number of sections.
        >>> tape2 = tape.delete_section(125, 33, 1).delete_section(125, 33, 2).update_intro()
        >>> intro = tape2.read_section(125, 1, 451)
        >>> assert len(intro["DESCRIPTION"]) == 87
        >>> assert len(intro["SECTIONS"]) == 8

        It can also be used to update the lines of description.
        >>> intro = tape2.update_intro(**dict(DESCRIPTION=[" new description"])).read_section(125, 1, 451)
        >>> print(sandy.write_mf1(intro))
         1001.00000 9.991673-1          0          0          2          5 125 1451    1
         0.00000000 0.00000000          0          0          0          6 125 1451    2
         1.00000000 20000000.0          3          0         10          3 125 1451    3
         0.00000000 0.00000000          0          0          1          8 125 1451    4
         new description                                                   125 1451    5
                                        1        451         13          0 125 1451    6
                                        2        151          4          0 125 1451    7
                                        3          1         35          0 125 1451    8
                                        3          2         35          0 125 1451    9
                                        3        102         35          0 125 1451   10
                                        4          2        196          0 125 1451   11
                                        6        102        201          0 125 1451   12
                                       33        102         21          0 125 1451   13
        """
        tape = self.data.copy()
        for mat, g in self.to_series().groupby("MAT"):
            intro = self.read_section(mat, 1, 451)
            intro.update(**kwargs)
            new_records = [(mf, mt, sec.count('\n') + 1, 0) for (mat, mf, mt), sec in g.items()]
            NWD, NXC = len(intro["DESCRIPTION"]), g.shape[0]
            new_records[0] = (1, 451, NWD+NXC+4, 0)
            intro["SECTIONS"] = new_records
            tape[(mat, 1, 451)] = sandy.write_mf1(intro)
        return self.__class__(tape)

    def get_nsub(self):
        """
        Determine ENDF-6 sub-library type by reading flag "NSUB" of first MAT
        in file.

        Returns
        -------
        `int`
            NSUB value

        Examples
        --------
        assert sandy.get_endf6_file("jeff_33", "xs", 10010).get_nsub() == "neutron"
        assert sandy.get_endf6_file("jeff_33", "xs", 10010).get_nsub().get_pendf(err=1).get_nsub() == "neutron"
        assert sandy.get_endf6_file("jeff_33", "nfpy", 942410).get_nsub() == "nfpy"
        assert sandy.get_endf6_file("jeff_33", "decay", 942410).get_nsub() == "decay"
        assert sandy.get_endf6_file("jeff_33", "dxs", 26000).get_nsub() == "neutron"
        assert sandy.get_endf6_file("proton", "dxs", 26000).get_nsub() == "proton"
        """
        nsub = self.read_section(self.mat[0], 1, 451)["NSUB"]
        return nsubs(nsub)


    def _handle_njoy_inputs(method):
        """
        Decorator to handle keyword arguments for NJOY before running
        the executable.

        Examples
        --------
        Test that `minimal_processing` filters unwanted modules.
        >>> g = sandy.get_endf6_file("jeff_33", "xs", 10010).get_gendf(err=1, minimal_processing=True, temperature=300, dryrun=True)
        >>> assert "broadr" in g and "reconr" in g
        >>> assert "thermr" not in g and "purr" not in g and "heatr" not in g and "unresr" not in g and "gaspr" not in g

        Test `minimal_processing=False`.
        >>> g = sandy.get_endf6_file("jeff_33", "xs", 10010).get_gendf(err=1, temperature=300, dryrun=True)
        >>> assert "broadr" in g and "reconr" in g
        >>> assert "thermr" in g and "purr" in g and "heatr" in g and "gaspr" in g

        Check that for `temperature=0` the calculation stops after RECONR.
        >>> g = sandy.get_endf6_file("jeff_33", "xs", 10010).get_gendf(err=1, dryrun=True)
        >>> assert "reconr" in g
        >>> assert "broadr" not in g and "thermr" not in g and "purr" not in g and "heatr" not in g and "unresr" not in g and "gaspr" not in g
        """
        def inner(
                self,
                temperature=0,
                err=0.001,
                minimal_processing=False,
                verbose=False,
                **kwargs,
                ):
            """
            Parameters
            ----------
            err : TYPE, optional
                 reconstruction tolerance for RECONR, BROADR and THERMR.
                 The default is 0.001.
            minimal_processing: `bool`, optional
                 deactivate modules THERMR, GASPR, HEATR, PURR and UNRESR.
                 The default is False.
            temperature : `float`, optional
                temperature of the cross sections in K. If not given, stop
                the processing after RECONR (before BROADR). The default is 0.
            verbose : `bool`, optional
                flag to print NJOY input file to screen before running the
                executable. The default is False.
            """
            kwds_njoy = kwargs.copy()

            # Handle 'minimal' processing options
            if minimal_processing or float(temperature) == 0:
                kwds_njoy["thermr"] = False
                kwds_njoy["gaspr"] = False
                kwds_njoy["heatr"] = False
                kwds_njoy["purr"] = False
                kwds_njoy["unresr"] = False
            # deactivate modules if temperature is 0
            if temperature == 0:
                kwds_njoy["broadr"] = False
                msg = """Zero or no temperature was requested, NJOY processing will stop after RECONR.
    If you want to process 0K cross sections use `temperature=0.1`.
    """
                logging.info(msg)

            # handle err
            reconr_kws = kwds_njoy.get("reconr_kws", {})
            broadr_kws = kwds_njoy.get("broadr_kws", {})
            thermr_kws = kwds_njoy.get("thermr_kws", {})
            reconr_kws["err"] = broadr_kws["err"] = thermr_kws["err"] = float(err)
            kwds_njoy["reconr_kws"] = reconr_kws
            kwds_njoy["broadr_kws"] = broadr_kws
            kwds_njoy["thermr_kws"] = thermr_kws
            
            kwds_njoy.update(dict(temperatures=[temperature], verbose=verbose))
    
            return method(self, **kwds_njoy)
        return inner

    def _handle_groupr_inputs(method):
        """
        Decorator to handle keyword arguments for NJOY before running
        the executable.

        Examples
        --------
        Test that `minimal_processing` filters unwanted modules.
        >>> g = sandy.get_endf6_file("jeff_33", "xs", 10010).get_gendf(err=1, minimal_processing=True, temperature=300, dryrun=True)
        >>> assert "broadr" in g and "reconr" in g
        >>> assert "thermr" not in g and "purr" not in g and "heatr" not in g and "unresr" not in g and "gaspr" not in g

        Test `minimal_processing=False`.
        >>> g = sandy.get_endf6_file("jeff_33", "xs", 10010).get_gendf(err=1, temperature=300, dryrun=True)
        >>> assert "broadr" in g and "reconr" in g
        >>> assert "thermr" in g and "purr" in g and "heatr" in g and "gaspr" in g

        Check that for `temperature=0` the calculation stops after RECONR.
        >>> g = sandy.get_endf6_file("jeff_33", "xs", 10010).get_gendf(err=1, dryrun=True)
        >>> assert "reconr" in g
        >>> assert "broadr" not in g and "thermr" not in g and "purr" not in g and "heatr" not in g and "unresr" not in g and "gaspr" not in g
        """
        def inner(
                self,
                groupr_kws={},
                **kwargs,
                ):
            """
            Parameters
            ----------
            err : TYPE, optional
                 reconstruction tolerance for RECONR, BROADR and THERMR.
                 The default is 0.001.
            minimal_processing: `bool`, optional
                 deactivate modules THERMR, GASPR, HEATR, PURR and UNRESR.
                 The default is False.
            temperature : `float`, optional
                temperature of the cross sections in K. If not given, stop
                the processing after RECONR (before BROADR). The default is 0.
            verbose : `bool`, optional
                flag to print NJOY input file to screen before running the
                executable. The default is False.
            """
            fission = 18 in self.get_records().query("MF==3").MT.values
            groupr_kws["nubar"] = fission
            fission = 18 in self.get_records().query("MF==5").MT.values
            groupr_kws["chi"] = fission
            groupr_kws["mubar"] = True
            return method(self, groupr_kws=groupr_kws, **kwargs)
        return inner

    def _handle_mf32_alone(method):
        """
        Decorator to handle files with section MF32 without section MF33.

        Examples
        --------
        91-Pa-231 in JEFF-3.3 has MF32 but not MF33.
        >>> tape = sandy.get_endf6_file("jeff_33", "xs", 912310)
        >>> err = tape.get_errorr(chi=False, nubar=True, mubar=False, err=1, xs=True, errorr33_kws=dict(irespr=0))
        >>> assert "errorr33" in err

        91-Pa-231 in JEFF-3.3 has MF32 but not MF33.
        >>> tape = sandy.get_endf6_file("jeff_33", "xs", 912330)
        >>> err = tape.get_errorr(chi=False, nubar=True, mubar=False, err=1, xs=True, errorr33_kws=dict(irespr=0))
        >>> assert "errorr33" in err

        17-Cl-37 in JEFF-3.3 has MF32 but not MF33.
        >>> tape = sandy.get_endf6_file("jeff_33", "xs", 170370)
        >>> err = tape.get_errorr(chi=False, nubar=True, mubar=False, err=1, xs=True, errorr33_kws=dict(irespr=0))
        >>> assert "errorr33" in err

        95-Am-241 in JEFF-3.3 has MF32 but not MF33.
        >>> tape = sandy.get_endf6_file("jeff_33", "xs", 952410)
        >>> err = tape.get_errorr(chi=False, nubar=True, mubar=False, err=1, xs=True, errorr33_kws=dict(irespr=0))
        >>> assert "errorr33" in err
        """
        def inner(
                self,
                **kwargs,
                ):
            """
            Parameters
            ----------
            kwargs : `dict`, optional
                 keyword arguments.
            """
            # input taken from
            # https://www-nds.iaea.org/index-meeting-crp/TM_NDP/docs/OCabellos_2017.pdf
            inst = self
            if 32 in self.mf and 33 not in self.mf:
                inp = sandy.njoy._input_mf32_nomf33 if 18 in self.mt else sandy.njoy._input_mf32_nomf33_no18
                with TemporaryDirectory() as tmpdir:
                    file = os.path.join(tmpdir, "tape20")
                    self.to_file(file)
                    outs = sandy.njoy._run_njoy(inp, file)
                    inst = sandy.Endf6.from_text(outs["errorr33"])
            return method(inst, **kwargs)
        return inner

    def read_section(self, mat, mf, mt, raise_error=True):
        """
        Parse MAT/MF/MT section.

        Parameters
        ----------
        `mat` : int
            MAT number
        `mf` : int
            MF number
        `mt` : int
            MT number

        Returns
        -------
        `dict`
        """
        read_module = f"read_mf{mf}"
        found = hasattr(sandy, read_module)
        if not raise_error and not found:
            return
        foo = eval(f"sandy.{read_module}")
        return foo(self, mat, mt)

    def _update_info(self, descr=None):
        """Update RECORDS item (in DATA column) for MF1/MT451 of each MAT based on the content of the TEXT column.
        """
        from .mf1 import write
        tape = self.copy()
        for mat in sorted(tape.index.get_level_values('MAT').unique()):
            sec = self.read_section(mat,1,451)
            records = pd.DataFrame(sec["RECORDS"], columns=["MF","MT","NC","MOD"]).set_index(["MF","MT"])
            new_records = []
            dfmat=tape.loc[mat]
#            for (mf,mt),text in sorted(tape.loc[mat].query('MT!=451'.format(mat)).TEXT.items()):
            for (mf,mt),text in sorted(dfmat[dfmat.index.get_level_values("MT")!=451].TEXT.items()):
                nc = len(text.splitlines())
                # when copying PENDF sections (MF2/MT152) mod is not present in the dictionary
                try:
                    mod = records.MOD.loc[mf,mt]
                except:
                    mod = 0
                new_records.append((mf,mt,nc,mod))
            if descr is not None:
                sec["TEXT"] = descr
            nc = 4 + len(sec["TEXT"]) + len(new_records) + 1
            mod = records.MOD.loc[1,451]
            new_records = [(1,451,nc,mod)] + new_records
            sec["RECORDS"] = new_records
            text = write(sec)
            tape.loc[mat,1,451].TEXT = text
        return Endf6(tape)

    def get_id(self, method="nndc"):
        """
        Extract ID for a given MAT for a ENDF-6 file.

        Parameters
        ----------
        method : `str`, optional
            Methods adopted to produce the ID. The default is `"nndc"`.
            - if `method='aleph'` the ID is the ZAM identifier
            - else, the ID is the ZA identifier according to the NNDC rules

        Returns
        -------
        ID : `int`
            ID of the ENDF-6 file.

        Notes
        -----
        .. note:: a warning is raised if more than one MAT is found.
                  Only the ID corresponding to the lowest MAT will be returned.
 
        Examples
        --------
        Extract ID for H1 file using NNDC and ALEPH methods
        >>> tape = sandy.get_endf6_file("jeff_33", "xs", 10010)
        >>> assert tape.get_id() == 1001
        >>> assert tape.get_id(method="aleph") == 10010

        Extract ID for Am242m file using NNDC and ALEPH methods
        >>> tape2 = sandy.get_endf6_file("jeff_33", "xs", 952421)
        >>> assert tape2.get_id() == 95642
        >>> assert tape2.get_id(method="ALEPH") == 952421

        >>> assert tape.merge(tape2).get_id() == 1001
        >>> assert tape2.merge(tape).get_id() == 1001
        """
        mat = self.mat[0]
        if len(self.mat) != 1:
            msg = "More than one MAT found, will give ID only for the lowest MAT"
            logging.warning(msg)
        info = self.read_section(mat, 1, 451)
        meta = info["LISO"]
        za = int(info["ZA"])
        zam = za * 10 + meta
        za_new = za + meta * 100 + 300 if meta else za
        ID = zam if method.lower() == "aleph" else za_new
        return ID

    @_handle_njoy_inputs
    def get_ace(self,
                suffix=None,
                pendf=None,
                **kwargs,
                ):
        """
        Process `Endf6` instance into an ACE file using NJOY.

        Parameters
        ----------
        dryrun : `bool`, optional
            Do not run NJOY and return NJOY input. Default is False.
        pendf : :func:`~sandy.Endf6`, optional
            provide manually PENDF object and add it to the processing
            sequence after RECONR and before BROADR. Default is None.
        suffix : `str`, optional
            suffix in the form `".[0-9][0-9]"` to assign to the ACE data.
            If not given, generate automatic suffix according to ALEPH rules.
            Default is None.
        **kwargs : `dict`
            keyword argument to pass to :func:`~sandy.njoy.process_neutron`.

        Returns
        -------
        `dict` of `str`
            output with `'ace'` and `'xsdir'` as keys.

        Examples
        --------
        Check that output is a ace file.
        >>> e6 = sandy.get_endf6_file("jeff_33", "xs", 10010)
        >>> ace = e6.get_ace(temperature=700, err=1, minimal_processing=True)["ace"]
        >>> assert "1001.07c" in ace
        >>> assert "sandy runs acer" in ace
        >>> assert "mat 125" in ace
        
        Check that ace is processed at a different temperature.
        >>> ace = e6.get_ace(temperature=800, err=1, minimal_processing=True)["ace"]
        >>> assert "1001.08c" in ace
        Check xsdir.
        >>> print(outs[xsdir])
        1001.08c    0.999167 filename 0 1   1     3297     0     0 6.894E-08

        Check that ace file follows "nndc" nomenclature for metastable nuclides.
        >>> e6_m = sandy.get_endf6_file("jeff_33", "xs",521291)
        >>> ace_m = e6_m.get_ace(temperature=800, err=1, minimal_processing=True)["ace"]
        >>> assert "52529.08c" in ace_m

        Check that using option `pendf` results in the same output.
        >>> pendf = e6.get_pendf(temperature=0, err=1)
        >>> ace2 = e6.get_ace(temperature=800, err=1, , minimal_processing=True, pendf=pendf)["ace"]
        >>> assert ace == ace2

        Check that the option suffix is used correctly.
        >>> ace = e6.get_ace(temperature=800, suffix=".85", err=1)
        >>> assert "1001.85c" in ace

        Check input pendf file
        >>> import pytest
        >>> with pytest.raises(Exception) as e_info:
        ...    e6.get_ace(pendf=e6)
        """
        if suffix:
            kwargs["suffixes"] = [suffix]

        with TemporaryDirectory() as td:
            endf6file = join(td, "endf6_file")
            self.to_file(endf6file)
            # we don not call to_pendf because we might want to pass a pendf in input
            if pendf:
                if pendf.kind != 'pendf':
                    raise TypeError("kw argument 'pendf' does not contain a PENDF file")
                pendftape = join(td, "pendf_file")
                pendf.to_file(pendftape)
            else:
                pendftape = None
            outputs = sandy.njoy.process_neutron(
                endf6file,
                pendftape=pendftape,
                **kwargs,
                )
        if kwargs.get("dryrun", False):
            return outputs  # this contains the NJOY input
        return {k: outputs[k] for k in ["ace", "xsdir"]}

    @_handle_njoy_inputs
    def get_pendf(self, **kwargs,):
        """
        Process `Endf6` instance into an PENDF file using NJOY.

        Parameters
        ----------
        **kwargs : `dict`
            keyword argument to pass to :func:`~sandy.njoy.process_neutron`.

        Returns
        -------
        pendf : :func:`~sandy.Endf6`
            Pendf object

        Examples
        --------
        Default run.
        >>> endf6 = sandy.get_endf6_file("jeff_33", "xs", 10010)
        >>> out = endf6.get_pendf(verbose=True, temperature=293.6, err=1, minimal_processing=True)
        >>> assert isinstance(out, sandy.Endf6)
        """
        # always deactivate acer
        kwargs["acer"] = False

        with TemporaryDirectory() as td:
            endf6file = join(td, "endf6_file")
            self.to_file(endf6file)
            outputs = sandy.njoy.process_neutron(
                endf6file,
                suffixes=[0],
                **kwargs,
                )
        if kwargs.get("dryrun", False):
            return outputs  # this contains the NJOY input
        pendf = Endf6.from_text(outputs["pendf"])
        return pendf

    @_handle_njoy_inputs
    @_handle_groupr_inputs
    def get_gendf(self, **kwargs,):
        """
        Process `Endf6` instance into a Gendf file using NJOY.

        Parameters
        ----------
        **kwargs : `dict`
            keyword argument to pass to :func:`~sandy.njoy.process_neutron`.

        Returns
        -------
        gendf : :func:`~sandy.Gendf`
            Gendf object

        Examples
        --------
        Default run.
        >>> endf6 = sandy.get_endf6_file("jeff_33", "xs", 10010)
        >>> out = endf6.get_gendf(temperature=293.6, minimal_processing=True)
        >>> assert isinstance(out, sandy.Gendf)

        Test keyword `sigz`
        >>> out = endf6.get_gendf(groupr_kws=dict(sigz=[1e10, 1e2]))
        >>> assert 1e10 in sandy.gendf.read_mf1(out, 125)[sigz]
        >>> assert 1e10 in sandy.gendf.read_mf1(out, 125)[sigz]

        Test keyword `iwt`
        >>> g = endf6.get_gendf(groupr_kws=dict(iwt=3), dryrun=True)
        >>> found = re.search('groupr(.*)moder', g, flags=re.DOTALL).group().splitlines()
        assert "125 2 0 3 0 1 1 0 /" == found[2]

        Test keyword `ign`
        >>> g = endf6.get_gendf(groupr_kws=dict(ign=3), dryrun=True)
        >>> found = re.search('groupr(.*)moder', g, flags=re.DOTALL).group().splitlines()
        assert "125 3 0 2 0 1 1 0 /" == found[2]

        Test keyword `ek`
        >>> g = endf6.get_gendf(groupr_kws=dict(ek=sandy.energy_grids.CASMO12), dryrun=True)
        >>> found = re.search('groupr(.*)moder', g, flags=re.DOTALL).group().splitlines()
        >>> ek = np.array(list(map(float, found[7].replace("/", "").split())))
        >>> assert np.testing.array_allclose(ek, sandy.energy_grids.CASMO12, rtol=1e-14, atol=1e-14)

        Test groupr MFs and MTs for fissile and non-fissile nuclides
        >>> g = endf6.get_gendf(dryrun=True)
        >>> found = re.search('groupr(.*)moder', g, flags=re.DOTALL).group().splitlines()
        assert " ".join(found[6:10]) == '3/ 3 251 / 0/ 0/'
        U-238 test because it contains mubar, xs, chi and nubar
        >>> endf6 = sandy.get_endf6_file('jeff_33','xs', 922380)
        >>> g = endf6.get_gendf(dryrun=True)
        >>> found = re.search('groupr(.*)moder', g, flags=re.DOTALL).group().splitlines()
        >>> assert " ".join(found[6:15]) == '3/ 3 452 / 3 455 / 3 456 / 3 251 / 5/ 5 18 / 0/ 0/'

        Test custom MTs
        >>> endf6 = sandy.get_endf6_file("jeff_33", "xs", 10010)
        >>> g = endf6.get_gendf(dryrun=True, groupr_kws=dict(mt=4))
        >>> found = re.search('groupr(.*)moder', g, flags=re.DOTALL).group().splitlines()
        assert " ".join(found[6:10]) == '3 4 / 3 251 / 0/ 0/'
        >>> g = endf6.get_gendf(dryrun=True, groupr_kws=dict(mt=4))
        >>> found = re.search('groupr(.*)moder', g, flags=re.DOTALL).group().splitlines()
        assert " ".join(found[6:11]) == '3 4 / 3 102 / 3 251 / 0/ 0/'
        """
        groupr_kws = kwargs.get("groupr_kws", {})
        fission = 18 in self.get_records().query("MF==3").MT.values
        groupr_kws["nubar"] = fission
        groupr_kws["chi"] = fission
        groupr_kws["mubar"] = True
        kwargs["groupr_kws"] = groupr_kws

        # always activate groupr
        kwargs["groupr"] = True

        # always deactivate acer
        kwargs["acer"] = False

        with TemporaryDirectory() as td:
            endf6file = os.path.join(td, "endf6_file")
            self.to_file(endf6file)
            outputs = sandy.njoy.process_neutron(
                    endf6file,
                    suffixes=[0],
                    **kwargs,
                    )
        if kwargs.get("dryrun", False):
            return outputs  # this contains the NJOY input
        gendf = sandy.Gendf.from_text(outputs["gendf"])
        return gendf

    @_handle_njoy_inputs
    @_handle_groupr_inputs
    @_handle_mf32_alone
    def get_errorr(self,
                   nubar=True,
                   mubar=True,
                   chi=True,
                   xs=True,
                   **kwargs,
                   ):
        """
        Process `Endf6` instance into a Errorr file using NJOY.

        Parameters
        ----------
        chi : `bool`, optional
            Process the chi covariance (default is `True`)
        mubar : `bool`, optional
            Process the mubar covariance (default is `True`)
        nubar : `bool`, optional
            Process the nubar covariance (default is `True`)
        xs : `bool`, optional
            Process the xs covariance (default is `True`)
        **kwargs : `dict`
            keyword argument to pass to `sandy.njoy.process`.

        spectrum_errorr : iterable, optional
            weight function as a iterable (default is None)

        Returns
        -------
        errorr : `sandy.Errorr` or `None`
            - `Errorr` instance constaining the nuclear data of the ERRORR
              file, if covariance information is found.
            - `None` if no covariance information is found.

        Notes
        -----
        .. note:: method arguments are consistent with those of `get_pendf`

        Examples
        --------
        Default run
        >>> endf6 = sandy.get_endf6_file("jeff_33", "xs", 942410)
        >>> out = endf6.get_errorr(temperature=300, minimal_processing=True, err=1, errorr_kws=dict(ign=3, mt=18))
        Check `ign` and `ek`
        This test check also the type of each output
        >>> assert out["errorr33"].get_xs().data.shape[0] == 30
        >>> assert out["errorr31"].get_xs().data.shape[0] == 30
        >>> assert out["errorr34"].get_xs().data.shape[0] == 30
        >>> assert out["errorr33"].get_xs().data.shape[0] == 30
        Check `ign` and `ek`
        >>> endf6 = sandy.get_endf6_file("jeff_33", "xs", 10010)
        >>> out = endf6.get_errorr(errorr_kws=dict(ek=sandy.energy_grids.CASMO12))
        Check `mt`
        assert out["errorr33"].get_xs().data.squeeze().name == (9443, 2)
        assert out["errorr34"].get_xs().data.squeeze().name == (9443, 251)
        columns = out["errorr31"].get_xs().data.columns
        assert (9443, 452) in columns and (9443, 455) in columns and (9443, 456) in columns

        Check consistency between keywords errorr_kws and errorr33_kws
        >>> ekws = dict(irespr=0, iwt=5, ek=[1e-5, 2e7], mt=(16, 18, 102))
        >>> e6 = sandy.get_endf6_file("jeff_33", "xs", 942410)
        >>> inp1 = e6.get_errorr(temperature=300, dryrun=True, xs=True, chi=False, nubar=False, mubar=False, errorr_kws=ekws)
        >>> inp2 = e6.get_errorr(temperature=300, dryrun=True, xs=True, chi=False, nubar=False, mubar=False, errorr33_kws=ekws)
        >>> inp3 = e6.get_errorr(temperature=300, dryrun=True, xs=True, chi=False, nubar=False, mubar=False)
        >>> assert "groupr" not in inp1 and "groupr" not in inp2 and "groupr" not in inp3
        >>> assert inp1 == inp2 and inp1 != inp3
        Check consistency between keywords errorr_kws and errorr35_kws
        >>> inp1 = e6.get_errorr(temperature=300, dryrun=True, xs=False, chi=True, nubar=False, mubar=False, errorr_kws=ekws)
        >>> inp2 = e6.get_errorr(temperature=300, dryrun=True, xs=False, chi=True, nubar=False, mubar=False, errorr35_kws=ekws)
        >>> inp3 = e6.get_errorr(temperature=300, dryrun=True, xs=False, chi=True, nubar=False, mubar=False)
        >>> assert "groupr" in inp1 and "groupr" in inp2 and "groupr" in inp3
        >>> assert inp1 == inp2 and inp1 != inp3
        Check consistency between keywords errorr_kws and errorr31_kws
        >>> inp1 = e6.get_errorr(temperature=300, dryrun=True, xs=False, chi=False, nubar=True, mubar=False, errorr_kws=ekws)
        >>> inp2 = e6.get_errorr(temperature=300, dryrun=True, xs=False, chi=False, nubar=True, mubar=False, errorr31_kws=ekws)
        >>> inp3 = e6.get_errorr(temperature=300, dryrun=True, xs=False, chi=False, nubar=True, mubar=False)
        >>> assert inp1 == inp2 and inp1 != inp3
        >>> assert "groupr" in inp1 and "groupr" in inp2 and "groupr" in inp3
        Check consistency between keywords errorr_kws and errorr34_kws       
        >>> inp1 = e6.get_errorr(temperature=300, dryrun=True, xs=False, chi=False, nubar=False, mubar=True, errorr_kws=ekws)
        >>> inp2 = e6.get_errorr(temperature=300, dryrun=True, xs=False, chi=False, nubar=False, mubar=True, errorr34_kws=ekws)
        >>> inp3 = e6.get_errorr(temperature=300, dryrun=True, xs=False, chi=False, nubar=False, mubar=True)
        >>> assert inp1 == inp2 and inp1 != inp3
        >>> assert "groupr" in inp1 and "groupr" in inp2 and "groupr" in inp3
        >>> inp1 = e6.get_errorr(temperature=300, dryrun=True, errorr_kws=ekws)
        >>> inp2 = e6.get_errorr(temperature=300, dryrun=True, errorr33_kws=ekws, errorr31_kws=ekws, errorr34_kws=ekws, errorr35_kws=ekws)
        >>> assert inp1 == inp2
        >>> assert "groupr" in inp1 and "groupr" in inp2

        Check default options
        >>> g = sandy.get_endf6_file("jeff_33", "xs", 10010).get_errorr(temperature=300, dryrun=True)
        >>> found = re.search('errorr(.*)', g, flags=re.DOTALL).group().splitlines()
        Check ign(2), iwt (2), iprint (0) and relative (1) options
        >>> assert found[2] == '125 2 2 0 1 /'
        Check temperature (300) option
        >>> assert found[3] == '0 300.0 /'
        Check irespr (1) option
        >>> assert found[4] = '0 33 1/'
        Check options changes
        >>> ekws = dict(ign=3, iwt=5, iprint=True, relative=False)
        >>> g = sandy.get_endf6_file("jeff_33", "xs", 10010).get_errorr(temperature=400, dryrun=True)
        >>> found = re.search('errorr(.*)', g, flags=re.DOTALL).group().splitlines()
        >>> assert found[2] == '125 3 5 1 0 /'
        >>> assert found[3] == '0 400.0 /'
        >>> assert found[4] = '0 33 0/'

        Test spectrum:
        >>> spect = [1.000000e-5, 2.00000000, 3.000000e-2, 2.00000000, 5.800000e-2, 4.00000000, 3, 1]
        >>> out = endf6.get_errorr(spectrum_errorr=spect, ek_errorr=[1.000000e-5, 3.000000e-2, 5.800000e-2, 3], verbose=True, nubar=False, chi=False, mubar=False)
        moder
        20 -21 /
        reconr
        -21 -22 /
        'sandy runs njoy'/
        9237 0 0 /
        0.005 0. /
        0/
        errorr
        -21 -22 0 33 0 /
        9237 1 1 0 1 /
        0 0.0 /
        0 33 1/
        3 /
        1.00000e-05 3.00000e-02 5.80000e-02 3.00000e+00 /
         0.00000000 0.00000000          0          0          1          4
                  4          1                                            
         1.000000-5 2.00000000 3.000000-2 2.00000000 5.800000-2 4.00000000
         3.00000000 1.00000000                                            
        /
        stop
        """
        # Activate specific errorr module according to covariance info and input options
        mf31 = self.get_records().query("MF==31")
        errorr31 = False if mf31.empty else nubar
        mf33 = self.get_records().query("MF==33")
        errorr33 = False if mf33.empty else xs
        mf34 = self.get_records().query("MF==34")
        errorr34 = False if mf34.empty else mubar
        mf35 = self.get_records().query("MF==35")
        errorr35 = False if mf35.empty else chi
        kwargs.update(dict(
            errorr31=errorr31,
            errorr33=errorr33,
            errorr34=errorr34,
            errorr35=errorr35,
            ))

        # Always deactivate acer
        kwargs["acer"] = False

        # keyword arguments in error_kws, if any, overwrite the others
        errorr_kws = kwargs.get("errorr_kws", {})
        errorr31_kws = kwargs.get("errorr31_kws", {})
        errorr31_kws.update(**errorr_kws)
        errorr33_kws = kwargs.get("errorr33_kws", {})
        errorr33_kws.update(**errorr_kws)
        errorr34_kws = kwargs.get("errorr34_kws", {})
        errorr34_kws.update(**errorr_kws)
        errorr35_kws = kwargs.get("errorr35_kws", {})
        errorr35_kws.update(**errorr_kws)
        kwargs.update(dict(
            errorr31_kws=errorr31_kws,
            errorr33_kws=errorr33_kws,
            errorr34_kws=errorr34_kws,
            errorr35_kws=errorr35_kws,
            ))

        with TemporaryDirectory() as td:
            endf6file = join(td, "endf6_file")
            self.to_file(endf6file)
            # update kwargs, or else error because multiple keyword argument
            outputs = sandy.njoy.process_neutron(
                    endf6file,
                    suffixes=[0],
                    **kwargs,
                    )
        if kwargs.get("dryrun", False):
            return outputs  # this contains the NJOY input
        outputs = {k: sandy.Errorr.from_text(v) for k, v in outputs.items() if "errorr" in k}
        return outputs

    def get_records(self):
        """
        Extract MAT, MF and MT combinations avaialbel in the file and 
        report it in tabulated format.

        Returns
        -------
        df : `pd.DataFrame`
            Dataframe with MAT, MF and MT as columns.

        Examples
        --------
        Short test for hydrogen
        >>> sandy.get_endf6_file("jeff_33", "xs", 10010).get_records()
            MAT	MF	MT
        0	125	1	451
        1	125	2	151
        2	125	3	1
        3	125	3	2
        4	125	3	102
        5	125	4	2
        6	125	6	102
        7	125	33	1
        8	125	33	2
        9	125	33	102
        """
        df = self.to_series().rename("TEXT").reset_index().drop("TEXT", axis=1)
        return df

    def get_perturbations(self, nsmp, njoy_kws={}, smp_kws={}, **kwargs,):
        """
        Construct multivariate distributions with a unit vector for 
        mean and with relative covariances taken from the evaluated files
        processed with the NJOY module ERRORR.

        Perturbation factors are sampled with the same multigroup structure of 
        the covariance matrix and are returned by nuclear datatype as a `dict`
        of `pd.Dataframe` instances .

        Parameters
        ----------
        nsmp : TYPE
            DESCRIPTION.
        njoy_kws : TYPE, optional
            DESCRIPTION. The default is {}.
        smp_kws : TYPE, optional
            DESCRIPTION. The default is {}.
        **kwargs : TYPE
            DESCRIPTION.

        Returns
        -------
        smp : TYPE
            DESCRIPTION.

        Examples
        --------
        Generate a couple of samples from the H1 file of JEFF-3.3.
        >>> njoy_kws = dict(err=1, errorr_kws=dict(mt=102))
        >>> tape = sandy.get_endf6_file("jeff_33", "xs", 10010)
        >>> smps = tape.get_perturbations(nsmp=2, njoy_kws=njoy_kws)
        >>> assert len(smps) == 1
        >>> assert isinstance(smps[33], sandy.Samples)
        >>> assert (smps[33].data.index.get_level_values("MT") == 102).all()
        """
        smp = {}

        # -- produce ERRORR files with covariance data
        njoy_kws["mubar"] = False
        njoy_kws["chi"] = False
        outs = self.get_errorr(**njoy_kws)
        filename = "PERT_{}_MF{}.xlsx"
        filenamecov = "COV_{}_MF{}.tape"

        # -- Extract samples from MF31 covariance data
        if "errorr31" in outs:
            outs["errorr31"].to_file(filenamecov.format(self.get_id(), 31))
            xls = filename.format(self.get_id(), 31)
            smp[31] = outs["errorr31"].get_cov().sampling(nsmp, to_excel=xls, seed=smp_kws.get("seed31"), **smp_kws)

        # -- Extract samples from MF33 covariance data
        if "errorr33" in outs:
            outs["errorr33"].to_file(filenamecov.format(self.get_id(), 33))
            xls = filename.format(self.get_id(), 33)
            smp[33] = outs["errorr33"].get_cov().sampling(nsmp, to_excel=xls, seed=smp_kws.get("seed33"), **smp_kws)

        return smp
    
    def apply_perturbations(self, smps, processes=1, pendf=None, njoy_kws={}, **kwargs):
        """
        Apply relative perturbations to the data contained in ENDF6 file.
        At the moment only the procedure for cross sections and nubar is
        implemented.
        Options are included to directly convert perturbed pendf to ace and
        to write data on files.
        
        Parameters
        ----------
        smps : `dict` of :obj: `~sandy.Samples` instances
            relative perturbation coefficients.
        processes : `int`, optional, default is `1`
            number of processes used to complete the task.
            Creation of perturbed pendf files and conversion to ACE 
            format is done in parallel if processes>1.
        pendf : :obj:`sandy.Endf6` instance, optiona, default is `None`
            if given apply pertubrations to provided instance, or else generate
            a pendf file from `self`.
        njoy_kws : `dict`
            keyword argument to pass to "tape.get_pendf()".
        **kwargs : `dict`
            keyword argument to pass to "tape.get_ace()" plus keyword arguments
            to pass to "endf6_perturb_worker".

        Returns
        -------
        A dictionary of endf/pendf file or ace files depending on to_ace.  

        Notes
        -----
        .. note:: ACE file temperature. Two options are implemented:
                  - Generation of pendf file at 0K and broadening to the 
                    required temperature when ACE file is created.
                  - Generation of pendf file at temperature and broadening not 
                    performed when ACE is created. This approach takes into 
                    account implicit effect.

        Examples
        --------
        Example to produce and apply perturbations to Pu-239 xs and nubar.
        >>> tape = sandy.get_endf6_file("jeff_33", "xs", 942390)
        >>> smps = tape.get_perturbations(2, njoy_kws=dict(err=1, chi=False, mubar=False, errorr33_kws=dict(mt=[2, 4, 18]),), smp_kws=dict(seed31=1, seed33=3))
        
        Let's apply both nubar and xs perturbations, then only nubar and then only xs.
        >>> outs_31_33 = tape.apply_perturbations(smps, njoy_kws=dict(err=1), processes=1)
        >>> outs_31 = tape.apply_perturbations({31: smps[31]}, njoy_kws=dict(err=1), processes=1)
        >>> outs_33 = tape.apply_perturbations({33: smps[33]}, njoy_kws=dict(err=1), processes=1)

        Check that files are different for different samples.
        >>> for i in range(2):
        ...    assert(outs_33[i]["endf6"].data == tape.data)
        ...    assert(outs_31[i]["endf6"].data != tape.data)
        ...    assert(outs_31[i]["endf6"].data == outs_31_33[i]["endf6"].data)
        ...    assert(outs_33[i]["pendf"].data != outs_31[i]["pendf"].data)
        ...    assert(outs_33[i]["pendf"].data == outs_31_33[i]["pendf"].data)

        Check that method is consistent only nubar, only xs or both nubar and xs are perturbed.
        >>> assert outs_33[0]["pendf"].data != outs_33[1]["pendf"].data
        >>> assert outs_33[0]["endf6"].data == outs_33[1]["endf6"].data
        >>> assert outs_31[0]["pendf"].data == outs_31[1]["pendf"].data
        >>> assert outs_31[0]["endf6"].data != outs_31[1]["endf6"].data

        Check that redundant nubar is also perturbed.
        >>> nu0 = sandy.Xs.from_endf6(outs_31[0]["endf6"].filter_by(listmt=[452, 455, 456]))
        >>> nu1 = sandy.Xs.from_endf6(outs_31[1]["endf6"].filter_by(listmt=[452, 455, 456]))
        >>> assert not nu0.data[9437, 452].equals(nu1.data[9437, 452])
        >>> assert nu0.data[9437, 455].equals(nu1.data[9437, 455])
        >>> assert not nu0.data[9437, 456].equals(nu1.data[9437, 456])
        
        Check that redundant and partial cross sections are correctly perturbed.
        >>> xs0 = sandy.Xs.from_endf6(outs_33[0]["pendf"].filter_by(listmf=[3]))
        >>> xs1 = sandy.Xs.from_endf6(outs_33[1]["pendf"].filter_by(listmf=[3]))
        >>> assert not xs0.data[9437, 1].equals(xs1.data[9437, 1])
        >>> assert not xs0.data[9437, 2].equals(xs1.data[9437, 2])
        >>> assert not xs0.data[9437, 4].equals(xs1.data[9437, 4])
        >>> assert not xs0.data[9437, 18].equals(xs1.data[9437, 18])
        >>> assert not xs0.data[9437, 51].equals(xs1.data[9437, 51])
        >>> assert xs0.data[9437, 16].equals(xs1.data[9437, 16])
        >>> assert xs0.data[9437, 102].equals(xs1.data[9437, 102])
        >>> assert xs0.data[9437, 103].equals(xs1.data[9437, 103])
        >>> assert xs0.data[9437, 107].equals(xs1.data[9437, 107])

        Check that ENDF6 and PENDF output filenames are correct
        >>> endf6 = sandy.get_endf6_file('jeff_33', 'xs', 10010)
        >>> smps = endf6.get_perturbations(2, njoy_kws=dict(err=0.1))
        >>> outs = endf6.apply_perturbations(smps, to_file=True)
        >>> assert outs[0]["endf6"] == '1001_0.endf6' and os.path.isfile('1001_0.endf6')
        >>> assert outs[0]["pendf"] == '1001_0.pendf' and os.path.isfile('1001_0.endf6')
        >>> assert outs[1]["endf6"] == '1001_1.endf6' and os.path.isfile('1001_1.endf6')
        >>> assert outs[1]["pendf"] == '1001_1.pendf' and os.path.isfile('1001_1.pendf')

        Check that ACE output filenames are correct
        >>> outs = endf6.apply_perturbations(smps, to_file=True, to_ace=True, ace_kws=dict(err=1, temperature=300, purr=False, heatr=False, thermr=False, gaspr=False))
        >>> assert outs[0]["ace"] == '1001_0.03c' and os.path.isfile('1001_0.03c')
        >>> assert outs[0]["xsdir"] == '1001_0.03c.xsd' and os.path.isfile('1001_0.03c.xsd')
        >>> assert outs[1]["ace"] == '1001_1.03c' and os.path.isfile('1001_1.03c')
        >>> assert outs[1]["xsdir"] == '1001_1.03c.xsd' and os.path.isfile('1001_1.03c.xsd')
        
        Check that keyword `pendf` works 
        >>> pendf = endf6.get_pendf(err=1)
        >>> outs1 = endf6.apply_perturbations(smps, njoy_kws=dict(err=1))
        >>> outs2 = endf6.apply_perturbations(smps, pendf=pendf)
        >>> assert outs1[0]["pendf"].write_string() == outs2[0]["pendf"].write_string()
        """

        if 33 not in smps and 31 not in smps:
            logging.info("no perturbation coefficient was found.")
            return

        if pendf:
            pendf_ = pendf
        else:
            pendf_ = self.get_pendf(**njoy_kws)

        data = {}
        if 31 in smps:
            data["pnu"] = smps[31].iterate_xs_samples()
        if 33 in smps:
            data["pxs"] = smps[33].iterate_xs_samples()

        if processes == 1:
            outs = {}

            while True:
                kws = {}

                # -- Iterate perturbation data (xs, nubar)
                for k, v in data.items():
                    item = next(v, False)
                    if not item:
                        break
                    n, s = item
                    kws[k] = s
                if not item:
                    break
                kws.update(**kwargs)
                outs[n] = endf6_perturb_worker(self.data, pendf_.data, n, **kws)

        elif processes > 1:
            pool = mp.Pool(processes=processes)
            outs = {}

            while True:
                kws = {}
                for k, v in data.items():
                    item = next(v, False)
                    if not item:
                        break
                    n, s = item
                    kws[k] = s
                if not item:
                    break
                kws.update(**kwargs)
                outs[n] = pool.apply_async(
                    endf6_perturb_worker,
                    (self.data, pendf_.data, n),
                    kws,
                    )

            outs = {n: out.get() for n, out in outs.items()}
            pool.close()
            pool.join()

        # if we keep ENDF6 and PENDF files in memory, convert them back into
        # sandy Endf6 instances (must do it here because Endf6 object cannot be pickled)
        if not kwargs.get("to_file", False) and not kwargs.get("to_ace", False):
            outs = {k: {k1: sandy.Endf6(v1) for k1, v1 in v.items()} for k, v in outs.items()}
        return outs
    


def endf6_perturb_worker(e6, pendf, n,
                         pxs=None,
                         pnu=None,
                         plpc=None,
                         pedistr=None,
                         verbose=False,
                         to_ace=False,
                         to_file=False,
                         filename="{ZA}_{SMP}",
                         ace_kws={},
                         **kwargs):
    """
    

    Parameters
    ----------
    e6 : TYPE
        DESCRIPTION.
    pendf : TYPE
        DESCRIPTION.
    n : TYPE
        DESCRIPTION.
    pxs : TYPE
        DESCRIPTION.
    verbose : TYPE, optional
        DESCRIPTION. The default is False.
    to_ace : TYPE, optional
        DESCRIPTION. The default is False.
    to_file : TYPE, optional
        DESCRIPTION. The default is False.
    filename : TYPE, optional
        DESCRIPTION. The default is "{ZA}_{SMP:d}".
    ace_kws : TYPE, optional
        DESCRIPTION. The default is {}.
    **kwargs : TYPE
        DESCRIPTION.

    Returns
    -------
    TYPE
        DESCRIPTION.

    """
    # default initialization
    endf6_pert = sandy.Endf6(e6.copy())
    pendf_pert = sandy.Endf6(pendf.copy())


    # filename options, in case we write to file
    mat = endf6_pert.mat[0]
    intro = endf6_pert.read_section(mat, 1, 451)
    za = int(intro["ZA"])
    meta = int(intro["LISO"])
    zam = sandy.zam.za2zam(za, meta=meta, method=False)
    zaid = ace_kws.get("zaid", "nndc")
    if zaid == "nndc":
        za = sandy.zam.zam2za(zam, method=zaid)[0]
    params = dict(
        MAT=mat,
        ZAM=zam,
        META=meta,
        ZA=za,
        SMP=n,
        )
    fn = filename.format(**params)

    # apply nubar perturbation
    if pnu is not None:
        nu = sandy.Xs.from_endf6(endf6_pert.filter_by(listmt=[452, 455, 456]))
        nu_pert = sandy.core.xs.xs_perturb_worker(nu, n, pnu, verbose=verbose)
        endf6_pert = nu_pert.reconstruct_sums(drop=True).to_endf6(endf6_pert).update_intro()

    # apply lpc perturbation
    if plpc is not None:
        pass

    # apply edistr perturbation
    if pedistr is not None:
        pass

    # apply xs perturbation
    if pxs is not None:
        xs = sandy.Xs.from_endf6(pendf_pert)
        xs_pert = sandy.core.xs.xs_perturb_worker(xs, n, pxs, verbose=verbose)
        pendf_pert = xs_pert.reconstruct_sums(drop=True).to_endf6(pendf_pert).update_intro()

    # Run NJOY and convert to ace
    if to_ace:
        temperature = ace_kws.get("temperature", 0)
        suffix = ace_kws.get("suffix", "." + sandy.njoy.get_temperature_suffix(temperature))
        ace = endf6_pert.get_ace(pendf=pendf_pert, **ace_kws)

        if to_file:
            outfiles = {}
            file = f"{fn}{suffix}c"
            with open(file, "w") as f:
                if verbose:
                    print(f"writing to file '{file}'")
                f.write(ace["ace"])
            outfiles["ace"] = file
            file = f"{file}.xsd"
            with open(file, "w") as f:
                if verbose:
                    print(f"writing to file '{file}'")
                f.write(ace["xsdir"])
            outfiles["xsdir"] = file
            return outfiles
        return ace

    else:
        out = {
            "endf6": endf6_pert.data,
            "pendf": pendf_pert.data,
            }

        if to_file:
            outfiles = {}
            file = f"{fn}.endf6"
            if verbose:
                print(f"writing to file '{file}'")
            endf6_pert.to_file(file)
            outfiles["endf6"] = file
            if pendf_pert:
                file = f"{fn}.pendf"
                if verbose:
                    print(f"writing to file '{file}'")
                pendf_pert.to_file(file)
                outfiles["pendf"] = file
            return outfiles

        return out

