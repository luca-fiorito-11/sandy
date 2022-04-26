# -*- coding: utf-8 -*-
"""
Created on Wed Dec  4 14:50:33 2019

@author: lfiorito
"""
import io
import shutil
import os
from functools import reduce
from tempfile import TemporaryDirectory
import logging
from urllib.request import urlopen, Request, urlretrieve
from zipfile import ZipFile
import re
import pytest

import numpy as np
import pandas as pd
import numpy as np

import sandy
from sandy.libraries import (
    N_FILES_ENDFB_71_IAEA,
    N_FILES_JEFF_32_NEA,
    N_FILES_JEFF_33_IAEA,
    N_FILES_JEFF_40T0_NEA,
    N_FILES_ENDFB_80_IAEA,
    N_FILES_JENDL_40U_IAEA,
    URL_N_ENDFB_71_IAEA,
    URL_N_JEFF_32_NEA,
    URL_N_JEFF_33_IAEA,
    URL_N_JEFF_40T0_NEA,
    URL_N_ENDFB_80_IAEA,
    URL_N_JENDL_40U_IAEA,
    URL_NFPY_ENDFB_80_IAEA,
    NFPY_FILES_ENDFB_80_IAEA,
    URL_NFPY_ENDFB_71_IAEA,
    NFPY_FILES_ENDFB_71_IAEA,
    URL_NFPY_JENDL_40U_IAEA,
    NFPY_FILES_JENDL_40U_IAEA,
    URL_NFPY_JEFF_33_IAEA,
    NFPY_FILES_JEFF_33_IAEA,
    URL_DECAY_ENDFB_71_IAEA,
    DECAY_FILES_ENDFB_71_IAEA,
    URL_DECAY_ENDFB_80_IAEA,
    DECAY_FILES_ENDFB_80_IAEA,
    URL_DECAY_JEFF_33_IAEA,
    DECAY_FILES_JEFF_33_IAEA,
    URL_TSL_ENDFB_71_IAEA,
    TSL_FILES_ENDFB_71_IAEA,
    URL_TSL_ENDFB_80_IAEA,
    TSL_FILES_ENDFB_80_IAEA,
    URL_TSL_JEFF_33_IAEA,
    TSL_FILES_JEFF_33_IAEA,
    URL_TSL_JENDL_40U_IAEA,
    TSL_FILES_JENDL_40U_IAEA,
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


def get_endf6_file(library, kind, zam, to_file=False):
    """
    Given a library and a nuclide import the corresponding ENDF-6 nuclear
    data file directly from internet.

    Parameters
    ----------
    library : `str`
        nuclear data library. Available libraries are:
        for `xs`
            * `'endfb_71'`
            * `'jeff_32'`
            * `'jeff_33'`
            * `'jeff_40t0'`
            * `'endfb_80'`
            * `'jendl_40u`
        for 'nfpy'
            * `'endfb_71'`
            * `'jeff_33'`
            * `'endfb_80'`
            * `'jendl_40u`
        for decay:
            * `'endfb_71'`
            * `'jeff_33'`
            * `'endfb_80'`
        for 'tsl' (read the note)
            * `'endfb_71'`
            * `'jeff_33'`
            * `'endfb_80'`
            * `'jendl_40u`
    kind : `str`
        nuclear data type:
            * 'xs' is a standard neutron-induced nuclear data file
            * 'nfpy' is a Neutron-Induced Fission Product Yields nuclear data
              file
            * 'decay' is a Radioactive Decay Data nuclear data file
            * 'tsl' is a Thermal Neutron Scattering Data file
    zam : `int` or 'all' or iterable
        zam = 'int' (individual nuclides) or iterable (group of nuclides)
            ZAM nuclide identifier $Z \\times 10000 + A \\times 10 + M$ where:
                * $Z$ is the charge number
                * $A$ is the mass number
                * $M$ is the metastate level (0=ground, 1=1st level)
        zam = 'all'
            We obtain the information of all the library. Actually,
            the only all libraries available are:
                * for 'nfpy': 'jeff_33'
                * for 'decay': 'jeff_33'
    Raises
    ------
    ValueError
        if library is not among available selection.

    ValueError
        if when you select 'xs', you select zam = 'all'

    Notes
    -----
    .. note:: In the `kind='tls` option, instead of the zam, integers are used.
              If you need help, the `get_tsl_index` function contains all the
              necessary information for the correct choice of these integers.

    Returns
    -------
    `Endf6`
        `Endf6` object with ENDF-6 data for specified library and nuclide.

    Examples
    --------

#    Import hydrogen file from JEFF-4.0T0.
#    >>> tape = sandy.get_endf6_file("jeff_40t0", 'xs', 10010)
#    >>> assert type(tape) is sandy.Endf6

    Import hydrogen file from JEFF-3.3.
    >>> tape = sandy.get_endf6_file("jeff_33", 'xs', 10010)
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

    Import Neutron-Induced Fission Product Yields for Th-227 from ENDF/B-VII.1.
    >>> tape = sandy.get_endf6_file("endfb_71", 'nfpy', 902270)
    >>> assert type(tape) is sandy.Endf6

    Import Neutron-Induced Fission Product Yields for Th-227 from ENDF/B-VIII.0
    >>> tape = sandy.get_endf6_file("endfb_80", 'nfpy', 902270)
    >>> assert type(tape) is sandy.Endf6

    Import Neutron-Induced Fission Product Yields for Th-227 from JENDL-4.0u
    >>> tape = sandy.get_endf6_file("jendl_40u", 'nfpy', 902270)
    >>> assert type(tape) is sandy.Endf6

    Import Neutron-Induced Fission Product Yields for Th-232 from JEFF-3.3
    >>> tape = sandy.get_endf6_file("jeff_33", 'nfpy', 902320)
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

#    Checked, but the test takes too long(~10 min), that's why it is commented.
#    Import all Radioactive Decay Data from ENDF/B-VIII.0.
#    >>> tape = sandy.get_endf6_file("endfb_80", 'decay', 'all')
#    >>> assert type(tape) is sandy.Endf6
    """
    foo_get = Endf6.from_zipurl
    foo_read = Endf6.read_zipurl
    if kind == 'xs':
        available_libs = (
            "jeff_32".upper(),
            "jeff_33".upper(),
            "jeff_40t0".upper(),
            "endfb_71".upper(),
            "endfb_80".upper(),
            "jendl_40u".upper(),
            )
        library_ = library.lower()
        if library_ == "jeff_40t0":
            url = URL_N_JEFF_40T0_NEA
            files = N_FILES_JEFF_40T0_NEA
            foo_read = Endf6.read_url
            foo_get = Endf6.from_url
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
            "jeff_33".upper(),
            )
        library_ = library.lower()
        if library_ == "endfb_71":
            url = URL_DECAY_ENDFB_71_IAEA
            files = DECAY_FILES_ENDFB_71_IAEA
        elif library_ == "endfb_80":
            url = URL_DECAY_ENDFB_80_IAEA
            files = DECAY_FILES_ENDFB_80_IAEA
        elif library_ == "jeff_33":
            url = URL_DECAY_JEFF_33_IAEA
            files = DECAY_FILES_JEFF_33_IAEA
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
    if str(zam).lower() == 'all':
        if kind.lower() == 'xs':
            raise ValueError("'all' option is not available for xs")
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
        """
        """
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
        >>> file = os.path.join(sandy.data.__path__[0], "h1.endf")
        >>> _FormattedFile.from_file(file).kind
        'endf6'

        >>> file = os.path.join(sandy.data.__path__[0], "h1.pendf")
        >>> _FormattedFile.from_file(file).kind
        'pendf'
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
            elif lrp == 1 or lrp == 0:
                kind = "endf6"
            else:
                kind == "unkwown"
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

        >>> filename = "n-1-H-001.jeff32"
        >>> rooturl = "https://www.oecd-nea.org/dbforms/data/eva/evatapes/jeff_32/"
        >>> file = sandy.Endf6.read_url(filename, rooturl)
        >>> print(file[0:890])
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
    def read_zipurl(filename, rooturl):
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
        req = Request(zipurl, headers={'User-Agent': 'Mozilla/5.0'})
        with urlopen(req) as zipresp:
            with ZipFile(io.BytesIO(zipresp.read())) as zfile:
                with TemporaryDirectory() as td:
                    zfile.extract(filename, path=td)
                    tmpfile = os.path.join(td, filename)
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
        >>> file = os.path.join(sandy.data.__path__[0], "h1.endf")
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
        """Create dataframe from endf6 text in string.

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
        >>> file = os.path.join(sandy.data.__path__[0], "h1.endf")
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

    def _get_section_df(self, mat, mf, mt, delimiter="?"):
        """
        """
        text = self.data[(mat, mf, mt)]

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
        >>> file = os.path.join(sandy.data.__path__[0], "h1.endf")
        >>> tape = _FormattedFile.from_file(file)
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

    def delete_sections(self, sections, raise_error=True):
        d = self.data.copy()
        for mat, mf, mt in sections:
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
        merged : TYPE
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
        >>> h1.merge(h2)
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
        128  1   451     1.002000+3 1.996800+0          0          0  ...
             2   151     1.002000+3 1.996800+0          0          0  ...
             3   1       1.002000+3 1.996800+0          0          0  ...
                 2       1.002000+3 1.996800+0          0          0  ...
                 3       1.002000+3 1.996800+0          0          0  ...
                 16      1.002000+3 1.996800+0          0          0  ...
                 102     1.002000+3 1.996800+0          0          0  ...
             4   2       1.002000+3 1.996800+0          0          2  ...
             6   16      1.002000+3 1.996800+0          0          1  ...
             8   102     1.002000+3 1.996800+0          0          0  ...
             9   102     1.002000+3 1.996800+0          0          0  ...
             12  102     1.002000+3 1.996800+0          1          0  ...
             14  102     1.002000+3 1.996800+0          1          0  ...
             33  1       1.002000+3 1.996800+0          0          0  ...
                 2       1.002000+3 1.996800+0          0          0  ...
                 16      1.002000+3 1.996800+0          0          0  ...
                 102     1.002000+3 1.996800+0          0          0  ...
        dtype: object

        Merge three files from different libraries.
        >>> h3 = sandy.get_endf6_file("endfb_71", 'xs', 10030)
        >>> h1.merge(h2, h3)
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
        128  1   451     1.002000+3 1.996800+0          0          0  ...
             2   151     1.002000+3 1.996800+0          0          0  ...
             3   1       1.002000+3 1.996800+0          0          0  ...
                 2       1.002000+3 1.996800+0          0          0  ...
                 3       1.002000+3 1.996800+0          0          0  ...
                 16      1.002000+3 1.996800+0          0          0  ...
                 102     1.002000+3 1.996800+0          0          0  ...
             4   2       1.002000+3 1.996800+0          0          2  ...
             6   16      1.002000+3 1.996800+0          0          1  ...
             8   102     1.002000+3 1.996800+0          0          0  ...
             9   102     1.002000+3 1.996800+0          0          0  ...
             12  102     1.002000+3 1.996800+0          1          0  ...
             14  102     1.002000+3 1.996800+0          1          0  ...
             33  1       1.002000+3 1.996800+0          0          0  ...
                 2       1.002000+3 1.996800+0          0          0  ...
                 16      1.002000+3 1.996800+0          0          0  ...
                 102     1.002000+3 1.996800+0          0          0  ...
        131  1   451     1.003000+3 2.989596+0          0          0  ...
             2   151     1.003000+3 2.989596+0          0          0  ...
             3   1       1.003000+3 2.989596+0          0          0  ...
                 2       1.003000+3 2.989596+0          0          0  ...
                 16      1.003000+3 2.989596+0          0          0  ...
             4   2       1.003000+3 2.989596+0          0          1  ...
                 16      1.003000+3 2.989596+0          0          2  ...
             5   16      1.003000+3 2.989596+0          0          0  ...
        dtype: object

        Merge two evaluations for the same nuclide.
        >>> h2_2 = sandy.get_endf6_file("jeff_32", 'xs', 10020)
        >>> h2.merge(h2_2)
        MAT  MF  MT
        128  1   451     1.002000+3 1.995712+0          0          0  ...
             2   151     1.002000+3 1.995712+0          0          0  ...
             3   1       1.002000+3 1.995712+0          0          0  ...
                 2       1.002000+3 1.995712+0          0          0  ...
                 3       1.002000+3 1.995712+0          0          0  ...
                 16      1.002000+3 1.995712+0          0          0  ...
                 102     1.002000+3 1.995712+0          0          0  ...
             4   2       1.002000+3 1.995712+0          0          1  ...
             6   16      1.002000+3 1.995712+0          0          1  ...
             8   102     1.002000+3 1.996800+0          0          0  ...
             9   102     1.002000+3 1.996800+0          0          0  ...
             12  102     1.002000+3 1.995712+0          1          0  ...
             14  102     1.002000+3 1.995712+0          1          0  ...
             33  1       1.002000+3 1.996800+0          0          0  ...
                 2       1.002000+3 1.996800+0          0          0  ...
                 16      1.002000+3 1.996800+0          0          0  ...
                 102     1.002000+3 1.996800+0          0          0  ...
        dtype: object
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
        `sandy.formats.endf6.BaseFile` or derived instance
            Copy of the original instance with filtered MAT, MF and MT sections
        """
        series = self.to_series()
        cond_mat = series.index.get_level_values("MAT").isin(listmat)
        cond_mf = series.index.get_level_values("MF").isin(listmf)
        cond_mt = series.index.get_level_values("MT").isin(listmt)
        d = series.loc[cond_mat & cond_mf & cond_mt].to_dict()
        return self.__class__(d, file=self.file)

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

    @classmethod
    def _from_old_format(cls, old_endf6):
        """
        Convert old endf6 tape into new one!
        """
        return cls(old_endf6.TEXT.to_dict())

    def write_string(self, title=""):
        """
        Write `_FormattedFile.data` content to string according to the ENDF-6
        file rules.

        Parameters
        ----------
        title : `str`, optional, default is an empty string
            first line of the file

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
        >>> file = os.path.join(sandy.data.__path__[0], "h1.endf")
        >>> string = sandy.Endf6.from_file(file).write_string()
        >>> print(string[:81 * 4 - 1])
                                                                             1 0  0    0
         1.001000+3 9.991673-1          0          0          2          5 125 1451    1
         0.000000+0 0.000000+0          0          0          0          6 125 1451    2
         1.000000+0 2.000000+7          3          0         10          3 125 1451    3

        if no modification is applied to the `_FormattedFile` content, the
        `write_string` returns an output identical to the file ASCII content.
        >>> assert string == open(file).read()

        Test with `sandy.Errorr` object and title option:
        >>> endf6 = sandy.get_endf6_file("jeff_33", "xs", 10010)
        >>> err = endf6.get_errorr(ek=[1e-2, 1e1, 2e7], err=1)
        >>> err.to_file("out.err", title="H with ERRORR")
        >>> err_2 = sandy.Errorr.from_file("out.err")
        >>> os.remove("out.err")
        >>> assert err_2.data[(125, 1, 451)] == err.data[(125, 1, 451)]
        >>> assert err_2.data[(125, 3, 102)] == err.data[(125, 3, 102)]
        >>> assert err_2.data[(125, 33, 102)] == err.data[(125, 33, 102)]

        ..note:: differences might appear from the way zeros were handled at
                 the end of ENDF-6 section, or if a different fiel title is
                 given
        """
        string = sandy.write_line(title, 1, 0, 0, 0)
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
        string += sandy.write_line("", -1, 0, 0, 0)
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
    read_section
        Parse MAT/MF/MT section.
    to_file
        Given a filename write the content of a `Endf6` instance to disk in
        ASCII format.
    to_string
        Write `Endf6.data` content to string according to the ENDF-6 file
        rules.
    write_string
        Write ENDF-6 content to string.
    """

    def _get_nsub(self):
        """
        Determine ENDF-6 sub-library type by reading flag "NSUB" of first MAT
        in file:

            * `NSUB = 10` : Incident-Neutron Data
            * `NSUB = 11` : Neutron-Induced Fission Product Yields

        Returns
        -------
        `int`
            NSUB value
        """
        return self.read_section(self.mat[0], 1, 451)["NSUB"]

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

    def parse(self):
        mats = self.index.get_level_values("MAT").unique()
        if len(mats) > 1:
            raise NotImplementedError("file contains more than 1 MAT")
        self.mat = self.endf = mats[0]
        if hasattr(self, "tape"):
            self.filename = os.path.basename(self.tape)
        INFO = self.read_section(mats[0], 1 ,451)
        del INFO["TEXT"], INFO["RECORDS"]
        self.__dict__.update(**INFO)
        self.SECTIONS = self.loc[INFO["MAT"]].reset_index()["MF"].unique()
        self.EHRES = 0
        self.THNMAX = - self.EHRES if self.EHRES != 0 else 1.0E6

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

    def get_ace(self,
                temperature,
                njoy=None,
                verbose=False,
                pendf=None,
                **kwargs,
                ):
        """
        Process `Endf6` instance into an ACE file using NJOY.

        Parameters
        ----------
        temperature : `float`
            temperature of the cross sections in K.
        njoy : `str`, optional, default is `None`
            NJOY executable, if `None` search in the system path.
        verbose : TYPE, optional, default is `False`
            flag to print NJOY input file to screen before running the
            executable.
        **kwargs : TYPE
            keyword argument to pass to `sandy.njoy.process`.

        Returns
        -------
        outs : `dict`
            output with `'ace'` and `'xsdir'` as keys pointing to the
            filenames of the corresponding ACE and xsdir files generated in
            the run.

        Examples
        --------
        >>> sandy.get_endf6_file("jeff_33", "xs", 10010).get_ace(700)
        {'ace': '1001.07c', 'xsdir': '1001.07c.xsd'}
        """
        outs = {}
        pendftape = None
        with TemporaryDirectory() as td:
            endf6file = os.path.join(td, "endf6_file")
            self.to_file(endf6file)
            if pendf:
                pendftape = os.path.join(td, "pendf_file")
                pendf.to_file(pendftape)
            text, inputs, outputs = sandy.njoy.process(
                endf6file,
                pendftape=pendftape,
                wdir=".",
                keep_pendf=False,
                exe=njoy,
                temperatures=[temperature],
                verbose=verbose,
                **kwargs,
                )
            acefile = outputs["tape50"]
            basename = os.path.split(acefile)[1]
            dest = os.path.join(os.getcwd(), basename)
            outs["ace"] = basename
            shutil.move(acefile, dest)
            xsdfile = outputs["tape70"]
            basename = os.path.split(xsdfile)[1]
            dest = os.path.join(os.getcwd(), basename)
            outs["xsdir"] = basename
            shutil.move(xsdfile, dest)
        return outs

    def get_pendf(self,
                  temperature=0,
                  njoy=None,
                  to_file=False,
                  verbose=False,
                  **kwargs,
                  ):
        """
        Process `Endf6` instance into an PENDF file using NJOY.

        Parameters
        ----------
        temperature : `float`, optional, default is `0`.
            temperature of the cross sections in K.
            If not given, stop the processing after RECONR (before BROADR).
        njoy : `str`, optional, default is `None`
            NJOY executable, if `None` search in the system path.
        to_file : `str`, optional, default is `None`
            if not `None` write processed ERRORR data to file.
            The name of the PENDF file is the keyword argument.
        verbose : `bool`, optional, default is `False`
            flag to print NJOY input file to screen before running the
            executable.
        **kwargs : `dict`
            keyword argument to pass to `sandy.njoy.process`.

        Returns
        -------
        pendf : `sandy.Endf6`
            `Endf6` instance constaining the nuclear data of the PENDF file.

        Examples
        --------
        Process H1 file from ENDF/B-VII.1 into PENDF
        >>> pendf =sandy.get_endf6_file("endfb_71", "xs", 10010).get_pendf(verbose=True)
        moder
        20 -21 /
        reconr
        -21 -22 /
        'sandy runs njoy'/
        125 0 0 /
        0.001 0. /
        0/
        moder
        -22 30 /
        stop

        >>> pendf
        MAT  MF  MT 
        125  1   451     1.001000+3 9.991673-1          2          0  ...
             2   151     1.001000+3 9.991673-1          0          0  ...
             3   1       1.001000+3 9.991673-1          0         99  ...
                 2       1.001000+3 9.991673-1          0          0  ...
                 102     1.001000+3 9.991673-1          0          0  ...
        dtype: object

        >>> pendf.kind
        'pendf'

        Test `to_file`
        >>> endf6 = sandy.get_endf6_file("jeff_33", "xs", 10010)
        >>> out = endf6.get_pendf(to_file="out.pendf")
        >>> assert os.path.isfile('out.pendf')

        """
        if float(temperature) == 0:
            kwargs["broadr"] = False
            kwargs["thermr"] = False
            kwargs["gaspr"] = False
            kwargs["heatr"] = False
            kwargs["purr"] = False
            kwargs["unresr"] = False
            msg = """Zero or no temperature was requested, NJOY processing will stop after RECONR.
If you want to process 0K cross sections use `temperature=0.1`.
"""
            logging.info(msg)
        with TemporaryDirectory() as td:
            endf6file = os.path.join(td, "endf6_file")
            self.to_file(endf6file)
            text, inputs, outputs = sandy.njoy.process(
                endf6file,
                acer=False,
                keep_pendf=True,
                exe=njoy,
                temperatures=[temperature],
                suffixes=[0],
                verbose=verbose,
                **kwargs,
                )
            pendffile = outputs["tape30"]
            pendf = Endf6.from_file(pendffile)
            if to_file:
                shutil.move(pendffile, to_file)
        return pendf

    def merge_pendf(self, pendf):
        """
        Merge endf-6 file content with that of a pendf file.

        Parameters
        ----------
        pendf : `sandy.Endf6`
            `Endf6` object containing pendf tape.

        Returns
        -------
        `sandy.Endf6`
            `Endf6` object with MF3 and MF1MT451 from PENDF

        Notes
        -----
        .. note:: the `Endf6` object in output has attribute `kind='pendf'`.

        .. note:: the `Endf6` object in output contains all sections from the
                  original endf-6 file, but for all MF=3 and MF=1 MT=451.

        Examples
        --------
        >>> file = os.path.join(sandy.data.__path__[0], "h1.endf")
        >>> endf6 = sandy.Endf6.from_file(file)
        >>> file = os.path.join(sandy.data.__path__[0], "h1.pendf")
        >>> pendf = sandy.Endf6.from_file(file)
        >>> merged = endf6.merge_pendf(pendf)
        >>> merged
        MAT  MF  MT
        125  1   451     1.001000+3 9.991673-1          2          0  ...
             2   151     1.001000+3 9.991673-1          0          0  ...
             3   1       1.001000+3 9.991673-1          0         99  ...
                 2       1.001000+3 9.991673-1          0          0  ...
                 102     1.001000+3 9.991673-1          0          0  ...
             4   2       1.001000+3 9.991673-1          0          1  ...
             6   102     1.001000+3 9.991673-1          0          2  ...
             33  1       1.001000+3 9.991673-1          0          0  ...
                 2       1.001000+3 9.991673-1          0          0  ...
                 102     1.001000+3 9.991673-1          0          0  ...
        dtype: object

        The cross section in the merged file come from the pendf.
        >>> assert merged.data[125, 3, 1] == pendf.data[125, 3, 1]
        >>> assert merged.data[125, 3, 1] != endf6.data[125, 3, 1]

        The new file is also a pendf.
        >>> merged.kind
        'pendf'
                """
        if pendf.kind != "pendf":
            raise sandy.Error("given file is not a pendf")
        section_3_endf6 = self.filter_by(listmf=[3]).data
        section_3_pendf = pendf.filter_by(listmf=[3]).data
        section_1451_pendf = pendf.filter_by(listmf=[1], listmt=[451]).data
        return self.delete_sections(section_3_endf6) \
                   .add_sections(section_3_pendf) \
                   .add_sections(section_1451_pendf)

    def get_errorr(self,
                   temperature=0,
                   njoy=None,
                   to_file=None,
                   verbose=False,
                   groupr=False,
                   err=0.005,
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
        groupr : `bool`, optional, default is `False`
            option to generate covariances from a multigroup cross section
            ..note:: this option is activated by default if `nubar=True` or
                     `chi=True`
        mubar : `bool`, optional
            Process the mubar covariance (default is `True`)
        njoy : `str`, optional, default is `None`
            NJOY executable, if `None` search in the system path.
        nubar : `bool`, optional
            Process the nubar covariance (default is `True`)
        temperature : `float`, optional, default is `0`.
            temperature of the cross sections in K.
            If not given, stop the processing after RECONR (before BROADR).
        to_file : `str`, optional, default is `None`
            if not `None` write processed ERRORR data to file.
            The name of the ERRORR file is the keyword argument.
        verbose : `bool`, optional, default is `False`
            flag to print NJOY input file to screen before running the
            executable.
        xs : `bool`, optional
            Process the xs covariance (default is `True`)
        **kwargs : `dict`
            keyword argument to pass to `sandy.njoy.process`.

        Parameters for RECONR
        ---------------------
        err : `float`, optional
            reconstruction tolerance (default is 0.005)

        Parameters for GROUPR
        ---------------------
        ign_groupr : `int`, optional
            neutron group option (default is 2, csewg 239-group structure)
        ek_groupr : iterable, optional
            derived cross section energy bounds
            (default is `[1e-5, 2e7]` if `ign_groupr==1`)
        sigz : iterable of `float`
            sigma zero values. The default is 1.0e10.
        iwt_groupr : `int`, optional
            weight function option (default is 2, constant)
        spectrum_groupr : iterable, optional
            Weight function as a iterable (default is None)

        Parameters for ERRORR
        ---------------------
        ign_errorr : `int`, optional
            neutron group option (default is 2, csewg 239-group structure)
        ek_errorr : iterable, optional
            derived cross section energy bounds
            (default is `[1e-5, 2e7]` if `ign_errorr==1`)
        iwt_errorr : `int`, optional
            weight function option (default is 2, constant)
        spectrum_errorr : iterable, optional
            weight function as a iterable (default is None)
        irespr: `int`, optional
            processing for resonance parameter covariances
            (default is 1, 1% sensitivity method)
        mt: `int` or iterable of `int`, optional
            list of MT reactions to be processed

            .. note:: this list will be used for all covariance types, i.e.,
                      MF31, MF33, MF34, MF35.
                      If this is not the expected behavior, use keyword
                      arguments `nubar`, `xs`, `mubar` and `chi`.

            .. note:: keyword `mt` is currently incompatible with keyword
                      `groupr`.

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
        Test verbose keyword
        >>> endf6 = sandy.get_endf6_file("jeff_33", "xs", 10010)
        >>> out = endf6.get_errorr(ek_errorr=sandy.energy_grids.CASMO12, verbose=True)
        moder
        20 -21 /
        reconr
        -21 -22 /
        'sandy runs njoy'/
        125 0 0 /
        0.005 0. /
        0/
        errorr
        -21 -22 0 33 0 /
        125 1 2 0 1 /
        0 0.0 /
        0 33 1/
        12 /
        1.00000e-05 3.00000e-02 5.80000e-02 1.40000e-01 2.80000e-01 3.50000e-01 6.25000e-01 4.00000e+00 4.80520e+01 5.53000e+03 8.21000e+05 2.23100e+06 1.00000e+07 /
        stop

        Test output type
        >>> assert isinstance(out, sandy.Errorr)

        Test `ign` and `ek`
        >>> assert out.get_xs().data[(125, 1)].size == 12
        
        Test `to_file`
        >>> out = endf6.get_errorr(to_file="out.err")
        >>> assert os.path.isfile('out.err')

        Test groupr and errorr:
        >>> out = endf6.get_errorr(verbose=True, groupr=True)
        moder
        20 -21 /
        reconr
        -21 -22 /
        'sandy runs njoy'/
        125 0 0 /
        0.005 0. /
        0/
        groupr
        -21 -22 0 -23 /
        125 2 0 2 0 1 1 0 /
        'sandy runs groupr' /
        0.0/
        10000000000.0/
        3/
        0/
        0/
        errorr
        -21 0 -23 33 0 /
        125 2 2 0 1 /
        0 0.0 /
        0 33 1/
        stop

        Test groupr and errorr for neutron energy grids:
        >>> out = endf6.get_errorr(ek_errorr=sandy.energy_grids.CASMO12, ek_groupr=sandy.energy_grids.CASMO12, verbose=True, groupr=True)
        moder
        20 -21 /
        reconr
        -21 -22 /
        'sandy runs njoy'/
        125 0 0 /
        0.005 0. /
        0/
        groupr
        -21 -22 0 -23 /
        125 1 0 2 0 1 1 0 /
        'sandy runs groupr' /
        0.0/
        10000000000.0/
        12 /
        1.00000e-05 3.00000e-02 5.80000e-02 1.40000e-01 2.80000e-01 3.50000e-01 6.25000e-01 4.00000e+00 4.80520e+01 5.53000e+03 8.21000e+05 2.23100e+06 1.00000e+07 /
        3/
        0/
        0/
        errorr
        -21 0 -23 33 0 /
        125 1 2 0 1 /
        0 0.0 /
        0 33 1/
        12 /
        1.00000e-05 3.00000e-02 5.80000e-02 1.40000e-01 2.80000e-01 3.50000e-01 6.25000e-01 4.00000e+00 4.80520e+01 5.53000e+03 8.21000e+05 2.23100e+06 1.00000e+07 /
        stop

        Test groupr and errorr for neutron and photons energy grids:
        >>> out = endf6.get_errorr(ek_groupr=sandy.energy_grids.CASMO12, ek_errorr=sandy.energy_grids.CASMO12, ep=sandy.energy_grids.CASMO12, verbose=True, groupr=True)
        moder
        20 -21 /
        reconr
        -21 -22 /
        'sandy runs njoy'/
        125 0 0 /
        0.005 0. /
        0/
        groupr
        -21 -22 0 -23 /
        125 1 1 2 0 1 1 0 /
        'sandy runs groupr' /
        0.0/
        10000000000.0/
        12 /
        1.00000e-05 3.00000e-02 5.80000e-02 1.40000e-01 2.80000e-01 3.50000e-01 6.25000e-01 4.00000e+00 4.80520e+01 5.53000e+03 8.21000e+05 2.23100e+06 1.00000e+07 /
        12 /
        1.00000e-05 3.00000e-02 5.80000e-02 1.40000e-01 2.80000e-01 3.50000e-01 6.25000e-01 4.00000e+00 4.80520e+01 5.53000e+03 8.21000e+05 2.23100e+06 1.00000e+07 /
        3/
        0/
        0/
        errorr
        -21 0 -23 33 0 /
        125 1 2 0 1 /
        0 0.0 /
        0 33 1/
        12 /
        1.00000e-05 3.00000e-02 5.80000e-02 1.40000e-01 2.80000e-01 3.50000e-01 6.25000e-01 4.00000e+00 4.80520e+01 5.53000e+03 8.21000e+05 2.23100e+06 1.00000e+07 /
        stop

        U-238 test because it contains mubar, xs, chi and nubar:
        >>> endf6 = sandy.get_endf6_file('jeff_33','xs', 922380)
        >>> out = endf6.get_errorr(ek_errorr=sandy.energy_grids.CASMO12, ek_groupr=sandy.energy_grids.CASMO12, verbose=True, err=1)
        moder
        20 -21 /
        reconr
        -21 -22 /
        'sandy runs njoy'/
        9237 0 0 /
        1 0. /
        0/
        groupr
        -21 -22 0 -23 /
        9237 1 0 2 0 1 1 0 /
        'sandy runs groupr' /
        0.0/
        10000000000.0/
        12 /
        1.00000e-05 3.00000e-02 5.80000e-02 1.40000e-01 2.80000e-01 3.50000e-01 6.25000e-01 4.00000e+00 4.80520e+01 5.53000e+03 8.21000e+05 2.23100e+06 1.00000e+07 /
        3/
        3 251 'mubar' /
        5/
        5 18 'chi' /
        0/
        0/
        errorr
        -21 0 -23 31 0 /
        9237 1 2 0 1 /
        0 0.0 /
        0 31 1/
        12 /
        1.00000e-05 3.00000e-02 5.80000e-02 1.40000e-01 2.80000e-01 3.50000e-01 6.25000e-01 4.00000e+00 4.80520e+01 5.53000e+03 8.21000e+05 2.23100e+06 1.00000e+07 /
        errorr
        -21 0 -23 33 0 /
        9237 1 2 0 1 /
        0 0.0 /
        0 33 1/
        12 /
        1.00000e-05 3.00000e-02 5.80000e-02 1.40000e-01 2.80000e-01 3.50000e-01 6.25000e-01 4.00000e+00 4.80520e+01 5.53000e+03 8.21000e+05 2.23100e+06 1.00000e+07 /
        errorr
        -21 0 -23 35 0 /
        9237 1 2 0 1 /
        0 0.0 /
        0 35 1/
        12 /
        1.00000e-05 3.00000e-02 5.80000e-02 1.40000e-01 2.80000e-01 3.50000e-01 6.25000e-01 4.00000e+00 4.80520e+01 5.53000e+03 8.21000e+05 2.23100e+06 1.00000e+07 /
        errorr
        -21 0 -23 34 0 /
        9237 1 2 0 1 /
        0 0.0 /
        0 34 1/
        12 /
        1.00000e-05 3.00000e-02 5.80000e-02 1.40000e-01 2.80000e-01 3.50000e-01 6.25000e-01 4.00000e+00 4.80520e+01 5.53000e+03 8.21000e+05 2.23100e+06 1.00000e+07 /
        stop

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


        >>> spect_g = [1.000000e-5, 1.00000000, 3.000000e-2, 2.00000000, 5.800000e-2, 3.00000000, 3, 2]
        >>> out = endf6.get_errorr(spectrum_errorr=spect, spectrum_groupr=spect_g, ek_errorr=[1.000000e-5, 3.000000e-2, 5.800000e-2, 3], ek_groupr=[1.000000e-5, 3.000000e-2, 5.800000e-2, 3], verbose=True, nubar=False, chi=False, mubar=False, groupr=True)
        moder
        20 -21 /
        reconr
        -21 -22 /
        'sandy runs njoy'/
        9237 0 0 /
        0.005 0. /
        0/
        groupr
        -21 -22 0 -23 /
        9237 1 0 1 0 1 1 0 /
        'sandy runs groupr' /
        0.0/
        10000000000.0/
        3 /
        1.00000e-05 3.00000e-02 5.80000e-02 3.00000e+00 /
         0.00000000 0.00000000          0          0          1          4
                  4          1                                            
         1.000000-5 1.00000000 3.000000-2 2.00000000 5.800000-2 3.00000000
         3.00000000 2.00000000                                            
        /
        3/
        0/
        0/
        errorr
        -21 0 -23 33 0 /
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

        Test irespr:
        out = endf6.get_errorr(spectrum_errorr=spect, ek_errorr=[1.000000e-5, 3.000000e-2, 5.800000e-2, 3], verbose=True, nubar=False, chi=False, mubar=False, irespr=0)
        moder
        20 -21 /
        reconr
        -21 -22 /
        'sandy runs njoy'/
        125 0 0 /
        0.005 0. /
        0/
        errorr
        -21 -22 0 33 0 /
        125 1 1 0 1 /
        0 0.0 /
        0 33 0/
        3 /
        1.00000e-05 3.00000e-02 5.80000e-02 3.00000e+00 /
         1.000000-5 2.00000000 3.000000-2 2.00000000 5.800000-2 4.00000000  3.00000000 1.00000000                                            
        /
        stop
    
        Test for MT:
        >>> endf6 = sandy.get_endf6_file("jeff_33", "xs", 10010)
        >>> out = endf6.get_errorr(verbose=True, mt=[1, 2], ek_errorr=sandy.energy_grids.CASMO12)
        moder
        20 -21 /
        reconr
        -21 -22 /
        'sandy runs njoy'/
        125 0 0 /
        0.005 0. /
        0/
        errorr
        -21 -22 0 33 0 /
        125 1 2 0 1 /
        0 0.0 /
        1 33 1/
        2 0 /
        1 2 /
        12 /
        1.00000e-05 3.00000e-02 5.80000e-02 1.40000e-01 2.80000e-01 3.50000e-01 6.25000e-01 4.00000e+00 4.80520e+01 5.53000e+03 8.21000e+05 2.23100e+06 1.00000e+07 /
        stop
        
        Keywords `mt` and `groupr` are compatible (`nubar` and `xs` covariance
        together):
        >>> with pytest.raises(sandy.SandyError):
        ...    sandy.get_endf6_file("jeff_33", "xs", 10010).get_errorr(err=1, mt=1, groupr=True)

        Test content of output `Errorr` file
        >>> out = sandy.get_endf6_file('jeff_33', "xs", 922350).get_errorr(err=1., irespr=0, mubar=False, chi=False)
        >>> keys = [(9228, 1, 451), (9228, 3, 456), (9228, 33, 456), (9228, 3, 1), (9228, 3, 2), (9228, 3, 4), (9228, 3, 16), (9228, 3, 17), (9228, 3, 18), (9228, 3, 37), (9228, 3, 102), (9228, 33, 1), (9228, 33, 2), (9228, 33, 4), (9228, 33, 16), (9228, 33, 17), (9228, 33, 18), (9228, 33, 37), (9228, 33, 102)]
        >>> for key in keys: assert key in out.data
        """
        kwds_njoy = kwargs.copy()
        if float(temperature) == 0:
            kwds_njoy["broadr"] = False
            kwds_njoy["thermr"] = False
            kwds_njoy["gaspr"] = False
            kwds_njoy["heatr"] = False
            kwds_njoy["purr"] = False
            kwds_njoy["unresr"] = False
            kwds_njoy['keep_pendf'] = False
        else:
            kwds_njoy["broadr"] = True
            kwds_njoy["thermr"] = kwds_njoy.get("thermr", False)
            kwds_njoy["gaspr"] = kwds_njoy.get("gaspr", False)
            kwds_njoy["heatr"] = kwds_njoy.get("heatr", False)
            kwds_njoy["purr"] = kwds_njoy.get("purr", False)
            kwds_njoy["unresr"] = kwds_njoy.get("unresr", False)
            kwds_njoy['keep_pendf'] = kwds_njoy.get('keep_pendf', False)

        cov_info = self.covariance_info(nubar=nubar, xs=xs,
                                        mubar=mubar, chi=chi)
        if not np.any(list(cov_info.values())):
            return  # no covariance found or wanted
        kwds_njoy.update(cov_info)

        # Mandatory groupr module activation
        groupr_ = True if (kwds_njoy["nubar"] or kwds_njoy["chi"] or "ek_groupr" in kwds_njoy or "spectrum_groupr" in kwds_njoy) else groupr

        with TemporaryDirectory() as td:
            endf6file = os.path.join(td, "endf6_file")
            self.to_file(endf6file)
            outputs = sandy.njoy.process(
                    endf6file,
                    errorr=True,
                    acer=False,
                    verbose=verbose,
                    temperatures=[temperature],
                    suffixes=[0],
                    err=err,
                    groupr=groupr_,
                    **kwds_njoy,
                    )[2]
            seq = map(sandy.Errorr.from_file, outputs.values())
            errorr = reduce(lambda x, y: x.merge(y), seq)
            if to_file:
                errorr.to_file(to_file)
        return errorr

    def get_gendf(self,
                  temperature=293.6,
                  njoy=None,
                  to_file=None,
                  verbose=False,
                  err=0.005,
                  nubar=False,
                  xs=True,
                  mubar=False,
                  chi=False,
                  **kwargs):
        """
        Process `Endf6` instance into a Gendf file using NJOY.

        Parameters
        ----------
        temperature : `float`, optional, default is `293.6`.
            temperature of the cross sections in K.
            If not given, stop the processing after RECONR (before BROADR).
        njoy : `str`, optional, default is `None`
            NJOY executable, if `None` search in the system path.
        to_file : `str`, optional, default is `None`
            if not `None` write processed GENDF data to file.
            The name of the GENDF file is the keyword argument.
        verbose : `bool`, optional, default is `False`
            flag to print NJOY input file to screen before running the
            executable.
        broadr : `bool`, optional, default is `True`
            option to generate gendf file with Doppler-broadened cross sections
        **kwargs : `dict`
            keyword argument to pass to `sandy.njoy.process`.

        Parameters for RECONR and BROADR
        --------------------------------
        err : `float`, optional
            reconstruction tolerance (default is 0.005)

        Parameters for GROUPR
        ---------------------
        chi : `bool`, optional
            Proccess the chi covariance(default is `False`)
        ign : `int`, optional
            neutron group option (default is 2, csewg 239-group structure)
        iwt_groupr : `int`, optional
            weight function option (default is 2, constant)
        mubar : `bool`, optional
            Proccess multigroup mubar (default is `False`)
        mt: `int` or iterable of `int`, optional
            run groupr only for the selected MT numbers
        nubar : `bool`, optional
            Proccess multigroup nubar (default is `False`)
        nuclide_production : `bool`, optional
            process multigroup activation yields (default is `False`)
        spectrum_groupr : iterable, optional
            Weight function as a iterable (default is None)
        sigz : iterable of `float`
            sigma zero values. The default is 1.0e10.
        xs : `bool`, optional
            Proccess multigroup xs (default is `True`)

        Returns
        -------
        gendf : `sandy.Gendf`
            Gendf object

        Examples
        --------
        Default test
        >>> endf6 = sandy.get_endf6_file("jeff_33", "xs", 10010)
        >>> out = endf6.get_gendf(verbose=True)
        moder
        20 -21 /
        reconr
        -21 -22 /
        'sandy runs njoy'/
        125 0 0 /
        0.005 0. /
        0/
        broadr
        -21 -22 -23 /
        125 1 0 0 0. /
        0.005 /
        293.6 /
        0 /
        groupr
        -21 -23 0 -24 /
        125 2 0 2 0 1 1 0 /
        'sandy runs groupr' /
        293.6/
        10000000000.0/
        3/
        0/
        0/
        moder
        -24 32 /
        stop

        Test keyword `sigz`
        >>> out = endf6.get_gendf(verbose=True, sigz=[1.0e10, 1e2])
        moder
        20 -21 /
        reconr
        -21 -22 /
        'sandy runs njoy'/
        125 0 0 /
        0.005 0. /
        0/
        broadr
        -21 -22 -23 /
        125 1 0 0 0. /
        0.005 /
        293.6 /
        0 /
        groupr
        -21 -23 0 -24 /
        125 2 0 2 0 1 2 0 /
        'sandy runs groupr' /
        293.6/
        10000000000.0 100.0/
        3/
        0/
        0/
        moder
        -24 32 /
        stop

        Test keyword `iwt_groupr`
        >>> out = endf6.get_gendf(verbose=True, iwt_groupr=3)
        moder
        20 -21 /
        reconr
        -21 -22 /
        'sandy runs njoy'/
        125 0 0 /
        0.005 0. /
        0/
        broadr
        -21 -22 -23 /
        125 1 0 0 0. /
        0.005 /
        293.6 /
        0 /
        groupr
        -21 -23 0 -24 /
        125 2 0 3 0 1 1 0 /
        'sandy runs groupr' /
        293.6/
        10000000000.0/
        3/
        0/
        0/
        moder
        -24 32 /
        stop

        Test keyword `ign_groupr`
        >>> out = endf6.get_gendf(verbose=True, ign_groupr=3)
        moder
        20 -21 /
        reconr
        -21 -22 /
        'sandy runs njoy'/
        125 0 0 /
        0.005 0. /
        0/
        broadr
        -21 -22 -23 /
        125 1 0 0 0. /
        0.005 /
        293.6 /
        0 /
        groupr
        -21 -23 0 -24 /
        125 3 0 2 0 1 1 0 /
        'sandy runs groupr' /
        293.6/
        10000000000.0/
        3/
        0/
        0/
        moder
        -24 32 /
        stop

        Test keyword `to_file`
        >>> endf6 = sandy.get_endf6_file("jeff_33", "xs", 10010)
        >>> out = endf6.get_gendf(to_file="out.gendf")
        >>> assert os.path.isfile('out.gendf')

        Test keyword `ek_groupr`
        >>> out = endf6.get_gendf(verbose=True, ek_groupr=sandy.energy_grids.CASMO12)
        moder
        20 -21 /
        reconr
        -21 -22 /
        'sandy runs njoy'/
        125 0 0 /
        0.005 0. /
        0/
        broadr
        -21 -22 -23 /
        125 1 0 0 0. /
        0.005 /
        293.6 /
        0 /
        groupr
        -21 -23 0 -24 /
        125 1 0 2 0 1 1 0 /
        'sandy runs groupr' /
        293.6/
        10000000000.0/
        12 /
        1.00000e-05 3.00000e-02 5.80000e-02 1.40000e-01 2.80000e-01 3.50000e-01 6.25000e-01 4.00000e+00 4.80520e+01 5.53000e+03 8.21000e+05 2.23100e+06 1.00000e+07 /
        3/
        0/
        0/
        moder
        -24 32 /
        stop

        U-238 test because it contains mubar, xs, chi and nubar:
        >>> endf6 = sandy.get_endf6_file('jeff_33','xs', 922380)
        >>> out = endf6.get_gendf(ek_groupr=sandy.energy_grids.CASMO12, verbose=True, err=1, nubar=True, mubar=True, chi=True)
        moder
        20 -21 /
        reconr
        -21 -22 /
        'sandy runs njoy'/
        9237 0 0 /
        1 0. /
        0/
        broadr
        -21 -22 -23 /
        9237 1 0 0 0. /
        1 /
        293.6 /
        0 /
        groupr
        -21 -23 0 -24 /
        9237 1 0 2 0 1 1 0 /
        'sandy runs groupr' /
        293.6/
        10000000000.0/
        12 /
        1.00000e-05 3.00000e-02 5.80000e-02 1.40000e-01 2.80000e-01 3.50000e-01 6.25000e-01 4.00000e+00 4.80520e+01 5.53000e+03 8.21000e+05 2.23100e+06 1.00000e+07 /
        3/
        3 251 'mubar' /
        5/
        5 18 'chi' /
        0/
        0/
        moder
        -24 32 /
        stop

        Test keyword `spectrum_groupr`
        >>> endf6 = sandy.get_endf6_file('jeff_33','xs', 10010)
        >>> spect = [1.000000e-5, 2.00000000, 3.000000e-2, 2.00000000, 5.800000e-2, 4.00000000, 3, 1]
        >>> out = endf6.get_gendf(spectrum_groupr=spect, ek_groupr=[1.000000e-5, 3.000000e-2, 5.800000e-2, 3], verbose=True, nubar=False, chi=False, mubar=False)
        moder
        20 -21 /
        reconr
        -21 -22 /
        'sandy runs njoy'/
        125 0 0 /
        0.005 0. /
        0/
        broadr
        -21 -22 -23 /
        125 1 0 0 0. /
        0.005 /
        293.6 /
        0 /
        groupr
        -21 -23 0 -24 /
        125 1 0 1 0 1 1 0 /
        'sandy runs groupr' /
        293.6/
        10000000000.0/
        3 /
        1.00000e-05 3.00000e-02 5.80000e-02 3.00000e+00 /
         0.00000000 0.00000000          0          0          1          4
                  4          1                                            
         1.000000-5 2.00000000 3.000000-2 2.00000000 5.800000-2 4.00000000
         3.00000000 1.00000000                                            
        /
        3/
        0/
        0/
        moder
        -24 32 /
        stop
        """
        kwds_njoy = kwargs.copy()
        if float(temperature) == 0:
            kwds_njoy["broadr"] = False
            kwds_njoy["thermr"] = False
            kwds_njoy["gaspr"] = False
            kwds_njoy["heatr"] = False
            kwds_njoy["purr"] = False
            kwds_njoy["unresr"] = False
            kwds_njoy['keep_pendf'] = False
        else:
            kwds_njoy["broadr"] = True
            kwds_njoy["thermr"] = kwds_njoy.get("thermr", False)
            kwds_njoy["gaspr"] = kwds_njoy.get("gaspr", False)
            kwds_njoy["heatr"] = kwds_njoy.get("heatr", False)
            kwds_njoy["purr"] = kwds_njoy.get("purr", False)
            kwds_njoy["unresr"] = kwds_njoy.get("unresr", False)
            kwds_njoy['keep_pendf'] = kwds_njoy.get('keep_pendf', False)

        kwds_njoy["acer"] = False
        kwds_njoy["keep_pendf"] = False

        kwds_njoy["nubar"] = nubar
        kwds_njoy["xs"] = xs
        kwds_njoy["chi"] = chi
        kwds_njoy["mubar"] = mubar

        with TemporaryDirectory() as td:
            endf6file = os.path.join(td, "endf6_file")
            self.to_file(endf6file)
            outputs = sandy.njoy.process(
                    endf6file,
                    groupr=True,
                    verbose=verbose,
                    temperatures=[temperature],
                    suffixes=[0],
                    err=err,
                    **kwds_njoy,
                    )[2]  # keep only gendf filename
            gendf_file = outputs["tape32"]
            groupr = sandy.Groupr.from_file(gendf_file)
            if to_file:
                groupr.to_file(to_file)
        return groupr

    def covariance_info(self, nubar=True, xs=True, mubar=True, chi=True):
        """
        Check the covariance information in the formatted file.

        Parameters
        ----------
        nubar : `bool`, optional
            default parameter for MF31 (default is `True`)
            it will overwrite what found in the file
        xs : `bool`, optional
            default parameter for MF33 (default is `True`)
            it will overwrite what found in the file
        mubar : `bool`, optional
            default parameter for MF34 (default is `True`)
            it will overwrite what found in the file
        chi : `bool`, optional
            default parameter for MF35 (default is `True`)
            it will overwrite what found in the file

        Returns
        -------
        cov_info : `dict`
            dictionary reporting if covariances were found.

        Notes
        -----
        .. note:: this method only works with MF31, MF33, MF34 and MF35

        Examples
        --------
        Check file contatining MF31, MF33, MF34 and MF35
        >>> endf6 = sandy.get_endf6_file('jeff_33', 'xs', 922380)
        >>> endf6.covariance_info()
        {'nubar': True, 'xs': True, 'mubar': True, 'chi': True}

        Set all values to `False`
        >>> endf6.covariance_info(xs=False, mubar=False, chi=False, nubar=False)
        {'nubar': False, 'xs': False, 'mubar': False, 'chi': False}

        2nd example without MF34
        >>> endf6 = sandy.get_endf6_file('jeff_33', 'xs', 922350)
        >>> endf6.covariance_info()
        {'nubar': True, 'xs': True, 'mubar': False, 'chi': True}

        If MF34 is not found, setting `mubar=True` won't change anything'
        >>> endf6 = sandy.get_endf6_file('jeff_33', 'xs', 922350)
        >>> endf6.covariance_info(mubar=True)
        {'nubar': True, 'xs': True, 'mubar': False, 'chi': True}

        All infos are `False` if no covariance is found
        >>> endf6 = sandy.get_endf6_file('jeff_33', 'xs', 10030)
        >>> endf6.covariance_info()
        {'nubar': False, 'xs': False, 'mubar': False, 'chi': False}
        """
        supported_mf = [31, 33, 34, 35]
        endf6_cov_mf = self.to_series().index.get_level_values("MF")\
                           .intersection(supported_mf)

        run_nubar = True if 31 in endf6_cov_mf else False
        run_xs = True if 33 in endf6_cov_mf else False
        run_mubar = True if 34 in endf6_cov_mf else False
        run_chi = True if 35 in endf6_cov_mf else False

        cov_info = {
            'nubar': run_nubar if nubar else False,
            'xs': run_xs if xs else False,
            'mubar': run_mubar if mubar else False,
            'chi': run_chi if chi else False,
            }
        return cov_info
