"""
Created on Wed Dec  4 14:50:33 2019

@author: lfiorito
"""
import io
import os
from os.path import dirname, join, splitext
import logging

import numpy as np
import pandas as pd


__author__ = "Luca Fiorito"


nsubs = {
    4: "decay",
    10: "neutron",
    11: "nfpy",
    10010: "proton",
    }


def get_endf6_file(library, kind, zam, to_file=False, local=False):
    """
    Retrieve an ENDF-6 evaluation for a given nuclear data library, data type,
    and nuclide (ZAM), either from the internet or from the local SANDY
    library appendix.

    This function provides a unified interface for downloading or loading:
    - neutron-induced cross sections       (kind="xs")
    - fission product yields               (kind="nfpy")
    - radioactive decay data               (kind="decay")
    - thermal scattering laws              (kind="tsl")
    - displacement cross sections          (kind="dxs")

    The function uses a centralized registry to map each `(kind, library)`
    pair to the appropriate remote URL and corresponding filename mapping.

    Parameters
    ----------
    library : str
        Name of the nuclear data library. Valid options depend on the
        requested `kind` and typically include:
        ``"endfb_71"``, ``"endfb_80"``, ``"endfb_81"``, ``"tendl_2023"``,
        ``"jeff_311"``, ``"jeff_33"``, ``"jeff_40"``, ``"jendl_40u"``,
        ``"jendl_5"``, ``"irdff_2"``, etc.

    kind : str
        Type of nuclear data to retrieve. One of:
            - ``"xs"``   : neutron-induced cross sections
            - ``"nfpy"`` : fission product yields
            - ``"decay"``: radioactive decay data
            - ``"tsl"``  : thermal scattering law
            - ``"dxs"``  : displacement cross sections

    zam : int, iterable of int, or "all"
        ZAM identifier(s) in the format ``Z*10000 + A*10 + M``.
            - integer: a single nuclide
            - iterable: multiple nuclides combined into one ENDF-6 tape
            - ``"all"``: retrieve the entire library (only for ``kind="decay"`` and ``kind="nfpy"``)

    to_file : bool, optional
        If True, write the resulting ENDF-6 file to disk using the
        nuclide name and library name as filename.
        Default is False.

    local : bool, optional
        If True, load files from SANDY's local library cache
        instead of downloading from the corresponding online source.
        Default is False.

    Returns
    -------
    Endf6
        A populated :class:`~sandy.endf6.Endf6` object containing the
        requested ENDF‑6 evaluation.

    Raises
    ------
    ValueError
        If the requested kind is unsupported, or the library is not
        available for that kind.

    ValueError
        If ``zam="all"`` is used with a kind that does not support it.

    Notes
    -----
    - ZAM format: ``Z*10000 + A*10 + M``  
      (M=0 for ground state, 1 for first metastable state, etc.)
    - For ``kind="tsl"``, the `zam` parameter must be replaced by the
      corresponding integer indices.
    - Remote downloads use ZIP files when available and automatically
      fall back to raw URLs.

    Examples
    --------

    Import hydrogen file from JEFF-3.3.

    >>> import sandy
    >>> tape = sandy.get_endf6_file("jeff_33", 'xs', 10010, local=True)
    >>> assert type(tape) is sandy.Endf6

    Import hydrogen file from JEFF-3.1.1.

    >>> tape = sandy.get_endf6_file("jeff_311", 'xs', 10010, local=True)
    >>> assert type(tape) is sandy.Endf6

    Import hydrogen file from ENDF/B-VII.1.

    >>> tape = sandy.get_endf6_file("endfb_71", 'xs', 10010, local=True)
    >>> assert type(tape) is sandy.Endf6

    Import hydrogen file from ENDF/B-VIII.0.

    >>> tape = sandy.get_endf6_file("endfb_80", 'xs', 10010, local=True)
    >>> assert type(tape) is sandy.Endf6

    Import hydrogen file from ENDF/B-VIII.1.

    >>> tape = sandy.get_endf6_file("endfb_81", 'xs', 10010, local=True)
    >>> assert type(tape) is sandy.Endf6

    Import hydrogen file from JENDL-4.0u

    >>> tape = sandy.get_endf6_file("jendl_40u", 'xs', 10010, local=True)
    >>> assert type(tape) is sandy.Endf6

    Import hydrogen file from JEFF-4.0, non local

    >>> tape = sandy.get_endf6_file("jeff_40", 'xs', 10010)
    >>> assert type(tape) is sandy.Endf6

    Import hydrogen file from TENDL-2023, non local

    >>> tape = sandy.get_endf6_file("tendl_2023", 'xs', 10010)
    >>> assert type(tape) is sandy.Endf6

    Import hydrogen file from JENDL-5

    >>> tape = sandy.get_endf6_file("jendl_5", 'xs', 10010, local=True)
    >>> assert type(tape) is sandy.Endf6

    Import Al-27 filr from IRDFF-II

    >>> tape = sandy.get_endf6_file("irdff_2", "xs", 130270, local=True)
    >>> assert type(tape) is sandy.Endf6

    Import a list of Decay Data for JEFF-3.3.

    >>> tape = sandy.get_endf6_file("jeff_33", 'decay', [10010, 270590, 270600], local=True)
    >>> assert type(tape) is sandy.Endf6

    Import all Decay Data for JEFF-3.3.

    >>> tape = sandy.get_endf6_file("jeff_33", 'decay', 'all')
    >>> assert type(tape) is sandy.Endf6

    Import all FIssion Product Yield Data for JEFF-3.3.

    >>> tape = sandy.get_endf6_file("jeff_33", 'nfpy', 'all')
    >>> assert type(tape) is sandy.Endf6
    
    Thermal Neutron Scattering Data from JEFF-3.3.

    >>> tape = sandy.get_endf6_file("jeff_33", 'tsl', [1, 2, 3], local=True)
    >>> assert type(tape) is sandy.Endf6

    """
    from functools import reduce

    from . import __file__ as sandy__file__
    from .zam import zam2nuclide
    from .libraries import (
        N_FILES_ENDFB_71_IAEA,
        N_FILES_ENDFB_80_IAEA,
        N_FILES_ENDFB_81_IAEA,
        N_FILES_JEFF_311_IAEA,
        N_FILES_JEFF_33_IAEA,
        N_FILES_JEFF_40_IAEA,
        N_FILES_TENDL_2023_IAEA,
        N_FILES_JENDL_40U_IAEA,
        N_FILES_JENDL_5_IAEA,
        N_FILES_IRDFF_2_IAEA,
        URL_N_ENDFB_71_IAEA,
        URL_N_JEFF_311_IAEA,
        URL_N_JEFF_33_IAEA,
        URL_N_JEFF_40_IAEA,
        URL_N_ENDFB_80_IAEA,
        URL_N_ENDFB_81_IAEA,
        URL_N_JENDL_40U_IAEA,
        URL_N_JENDL_5_IAEA,
        URL_N_TENDL_2023_IAEA,
        URL_N_IRDFF_2_IAEA,
    
        NFPY_FILES_ENDFB_71_IAEA,
        NFPY_FILES_ENDFB_80_IAEA,
        NFPY_FILES_ENDFB_81_IAEA,
        NFPY_FILES_JEFF_311_IAEA,
        NFPY_FILES_JEFF_33_IAEA,
        NFPY_FILES_JEFF_40_IAEA,
        NFPY_FILES_JENDL_40U_IAEA,
        NFPY_FILES_JENDL_5_IAEA,
        URL_NFPY_ENDFB_71_IAEA,
        URL_NFPY_ENDFB_80_IAEA,
        URL_NFPY_ENDFB_81_IAEA,
        URL_NFPY_JEFF_311_IAEA,
        URL_NFPY_JEFF_33_IAEA,
        URL_NFPY_JEFF_40_IAEA,
        URL_NFPY_JENDL_40U_IAEA,
        URL_NFPY_JENDL_5_IAEA,
    
        DECAY_FILES_ENDFB_71_IAEA,
        DECAY_FILES_ENDFB_80_IAEA,
        DECAY_FILES_ENDFB_81_IAEA,
        DECAY_FILES_JEFF_311_IAEA,
        DECAY_FILES_JEFF_33_IAEA,
        DECAY_FILES_JEFF_40_IAEA,
        DECAY_FILES_JENDL_5_IAEA,
        URL_DECAY_ENDFB_71_IAEA,
        URL_DECAY_ENDFB_80_IAEA,
        URL_DECAY_ENDFB_81_IAEA,
        URL_DECAY_JEFF_311_IAEA,
        URL_DECAY_JEFF_33_IAEA,
        URL_DECAY_JEFF_40_IAEA,
        URL_DECAY_JENDL_5_IAEA,
    
        TSL_FILES_ENDFB_71_IAEA,
        TSL_FILES_ENDFB_80_IAEA,
        TSL_FILES_ENDFB_81_IAEA,
        TSL_FILES_JEFF_33_IAEA,
        TSL_FILES_JEFF_40_IAEA,
        TSL_FILES_JENDL_40U_IAEA,
        TSL_FILES_JENDL_5_IAEA,
        URL_TSL_JENDL_40U_IAEA,
        URL_TSL_JENDL_5_IAEA,
        URL_TSL_ENDFB_71_IAEA,
        URL_TSL_ENDFB_80_IAEA,
        URL_TSL_ENDFB_81_IAEA,
        URL_TSL_JEFF_33_IAEA,
        URL_TSL_JEFF_40_IAEA,
    
        DXS_FILES_JEFF_33_IAEA,
        DXS_FILES_PROTON_IAEA,
        URL_DXS_JEFF_33_IAEA,
        URL_DXS_PROTON_IAEA
        )

    kind_ = kind.lower()
    library_ = library.lower()
    
    foo_get = Endf6.from_zipurl

    maps = {
        "xs": {
            "jeff_311": (URL_N_JEFF_311_IAEA, N_FILES_JEFF_311_IAEA),
            "jeff_33": (URL_N_JEFF_33_IAEA, N_FILES_JEFF_33_IAEA),
            "jeff_40": (URL_N_JEFF_40_IAEA, N_FILES_JEFF_40_IAEA),
            "endfb_71": (URL_N_ENDFB_71_IAEA, N_FILES_ENDFB_71_IAEA),
            "endfb_80": (URL_N_ENDFB_80_IAEA, N_FILES_ENDFB_80_IAEA),
            "endfb_81": (URL_N_ENDFB_81_IAEA, N_FILES_ENDFB_81_IAEA),
            "jendl_40u": (URL_N_JENDL_40U_IAEA, N_FILES_JENDL_40U_IAEA),
            "jendl_5": (URL_N_JENDL_5_IAEA, N_FILES_JENDL_5_IAEA),
            "irdff_2": (URL_N_IRDFF_2_IAEA, N_FILES_IRDFF_2_IAEA),
            "tendl_2023": (URL_N_TENDL_2023_IAEA, N_FILES_TENDL_2023_IAEA),
            },
        "nfpy": {
            "jeff_311": (URL_NFPY_JEFF_311_IAEA, NFPY_FILES_JEFF_311_IAEA),
            "jeff_33": (URL_NFPY_JEFF_33_IAEA, NFPY_FILES_JEFF_33_IAEA),
            "jeff_40": (URL_NFPY_JEFF_40_IAEA, NFPY_FILES_JEFF_40_IAEA),
            "endfb_71": (URL_NFPY_ENDFB_71_IAEA, NFPY_FILES_ENDFB_71_IAEA),
            "endfb_80": (URL_NFPY_ENDFB_80_IAEA, NFPY_FILES_ENDFB_80_IAEA),
            "endfb_81": (URL_NFPY_ENDFB_81_IAEA, NFPY_FILES_ENDFB_81_IAEA),
            "jendl_40u": (URL_NFPY_JENDL_40U_IAEA, NFPY_FILES_JENDL_40U_IAEA),
            "jendl_5": (URL_NFPY_JENDL_5_IAEA, NFPY_FILES_JENDL_5_IAEA),
            },
        "decay": {
            "jeff_311": (URL_DECAY_JEFF_311_IAEA, DECAY_FILES_JEFF_311_IAEA),
            "jeff_33": (URL_DECAY_JEFF_33_IAEA, DECAY_FILES_JEFF_33_IAEA),
            "jeff_40": (URL_DECAY_JEFF_40_IAEA, DECAY_FILES_JEFF_40_IAEA),
            "endfb_71": (URL_DECAY_ENDFB_71_IAEA, DECAY_FILES_ENDFB_71_IAEA),
            "endfb_80": (URL_DECAY_ENDFB_80_IAEA, DECAY_FILES_ENDFB_80_IAEA),
            "endfb_81": (URL_DECAY_ENDFB_81_IAEA, DECAY_FILES_ENDFB_81_IAEA),
            "jendl_5": (URL_DECAY_JENDL_5_IAEA, DECAY_FILES_JENDL_5_IAEA),
            },
        "tsl": {
            "jeff_33": (URL_TSL_JEFF_33_IAEA, TSL_FILES_JEFF_33_IAEA),
            "jeff_40": (URL_TSL_JEFF_40_IAEA, TSL_FILES_JEFF_40_IAEA),
            "endfb_71": (URL_TSL_ENDFB_71_IAEA, TSL_FILES_ENDFB_71_IAEA),
            "endfb_80": (URL_TSL_ENDFB_80_IAEA, TSL_FILES_ENDFB_80_IAEA),
            "endfb_81": (URL_TSL_ENDFB_81_IAEA, TSL_FILES_ENDFB_81_IAEA),
            "jendl_40u": (URL_TSL_JENDL_40U_IAEA, TSL_FILES_JENDL_40U_IAEA),
            "jendl_5": (URL_TSL_JENDL_5_IAEA, TSL_FILES_JENDL_5_IAEA),
            },
        "dxs": {
            "jeff_33": (URL_DXS_JEFF_33_IAEA, DXS_FILES_JEFF_33_IAEA),
            "proton": (URL_DXS_PROTON_IAEA, DXS_FILES_PROTON_IAEA),
            },
        }
    
    if kind_ not in maps:
        ValueError(f"option 'kind={kind_}' is not supported")

    if library_ not in maps[kind_]:
            raise ValueError(
                f"""library '{library}' is not available.
                Available libraries are: {maps[kind_].keys()}
                """
                )
    url, files = maps[kind_][library_]

    local_space = join(dirname(sandy__file__), "appendix", "libraries")
    local_path = join(local_space, library, kind)

    def local_foo_get(file, url):
        file = splitext(file)[0]
        tape = Endf6.from_zipfile(join(local_path, f"{file}.zip"))
        return tape
    
    if local:
        foo_get = local_foo_get


    if str(zam).lower() == 'all':
        if kind_ not in ['decay', "nfpy"]:
            raise ValueError(f"'all' option is not available for kind='{kind_}'")

        # --- fall back on local files with all data, otherwise it's too slow
        foo_get = local_foo_get
        tape = foo_get(f"{kind}_{library_}.dat", url="not used")

    else:

        # --- helper: try remote first, then fall back to local
        def fetch_with_fallback(key):
            try:
                return foo_get(files[key], url)    # remote first
            except Exception as err:
                logging.warning(
                    f"Remote download failed for ZAM={key}: {err}\n"
                    f"→ Falling back to local file in:\n"
                    fr"   {local_path}"
                )
                return local_foo_get(files[key], url=None)
    
        # --- get data
        if hasattr(zam, "__len__"):
            # --- CASE 1: a list of nuclides is passed (all valid ZAM)
            tapes = [fetch_with_fallback(x) for x in zam]
            tape = reduce(lambda x, y: x.add_sections(y.data), tapes)
        else:
            # --- CASE 2: a single nuclide is passed (valid ZAM)
            tape = fetch_with_fallback(zam)

    if to_file:
        basename = zam2nuclide(zam, atomic_number=True, sep="-")
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
            raise TypeError("'data' is not a 'dict'")
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

        >>> import sandy
        >>> assert sandy.get_endf6_file("jeff_33", "decay", 10010, local=True).kind == "endf6"
        >>> assert sandy.get_endf6_file("jeff_33", "nfpy", 922350, local=True).kind == "endf6"
        >>> assert sandy.get_endf6_file("jeff_33", "xs", 10010, local=True).kind == "endf6"
        >>> assert sandy.get_endf6_file("jeff_33", "xs", 10010, local=True).get_pendf(err=1).kind == "pendf"
        >>> assert sandy.get_endf6_file("jeff_33", "xs", 10010, local=True).get_gendf(err=1).kind == "gendf"
        >>> outs = sandy.get_endf6_file("jeff_33", "xs", 942410, local=True).get_errorr(err=1, errorr_kws=dict(mt=18))
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
    def from_zipurl(cls, filename, rooturl):
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
        Test the fallback option (read url, no zip).
        
        >>> import sandy
        >>> filename = "n-1-H-001.jeff32"
        >>> rooturl = "https://www.oecd-nea.org/dbforms/data/eva/evatapes/jeff_32/"
        >>> file = sandy.Endf6.from_zipurl(filename, rooturl)
        >>> print(file.write_string()[0:890])
                                                                             1 0  0    0
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

        Test the zip option (issue with the website. test removed.)

        # >>> filename = "decay_1907_57-La-149.dat"
        # >>> rooturl = "https://www-nds.iaea.org/public/download-endf/ENDF-B-VII.1/decay/"
        # >>> file = sandy.Endf6.from_zipurl(filename, rooturl)
        # >>> print(file.write_string()[0:971])
        # Retrieved by E4-util: 2012/01/16,13:45:44                            1 0  0    0
        # 5.714900+4 1.476553+2         -1          0          0          11907 1451    1
        # 0.000000+0 1.000000+0          0          0          0          61907 1451    2
        # 0.000000+0 0.000000+0          1          0          4          71907 1451    3
        # 0.000000+0 0.000000+0          0          0         27          21907 1451    4
        # 57-La-149  BNL        EVAL-AUG11 Conv. from CGM                  1907 1451    5
        # /ENSDF/                                               20111222   1907 1451    6
        # ----ENDF/B-VII.1      Material 1907                               1907 1451    7
        # -----RADIOACTIVE DECAY DATA                                       1907 1451    8
        # ------ENDF-6 FORMAT                                               1907 1451    9
        # *********************** Begin Description *********************** 1907 1451   10
        # **         ENDF/B-VII.1 RADIOACTIVE DECAY DATA FILE            ** 1907 1451   11

        """
        from urllib.request import urlopen, Request
        from zipfile import ZipFile
        from tempfile import TemporaryDirectory
        import requests

        # ---- Prepare URLs ----
        # it is assumed that filename has an extension like ".dat"
        rootname = splitext(filename)[0]
        zipurl = f"{rooturl}/{rootname}.zip"
        rawurl = f"{rooturl}/{filename}"

        # Standard desktop browser UA avoids 403 on many servers
        headers = {
            "User-Agent": (
                "Mozilla/5.0 (Windows NT 10.0; Win64; x64; rv:91.0) "
                "Gecko/20100101 Firefox/91.0"
            )
        }

        # ============================================================
        # 1) TRY READING THE ZIP FILE
        # ============================================================
        try:
            r = requests.get(zipurl, headers=headers, timeout=60)
            r.raise_for_status()             # <- triggers exception for 404, 403, etc.
            zip_data = r.content
    
            # Open ZIP in memory
            with ZipFile(io.BytesIO(zip_data)) as zfile:
                with TemporaryDirectory() as td:
                    # The zip contains the requested file
                    zfile.extract(filename, path=td)
                    tmpfile = join(td, filename)
                    with open(tmpfile, "r") as f:
                        text = f.read()
                        return cls.from_text(text)

        except Exception:
            # ZIP not found or extraction failed → fallback to RAW
            pass

        # ============================================================
        # 2) FALLBACK: READ DIRECTLY FROM RAW URL
        # ============================================================
        try:
            req = Request(rawurl, headers=headers)
            with urlopen(req) as f:
                text = f.read().decode("utf-8")
                return cls.from_text(text)
    
        except Exception as e:
            raise RuntimeError(
                f"Could not fetch ENDF file from either ZIP or raw URL:\n"
                f"ZIP attempted: {zipurl}\n"
                f"RAW attempted: {rawurl}\n"
                f"Error: {e}"
            )

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
        :obj:`sandy.endf6._FormattedFile` or derived instance
            Dataframe containing ENDF6 data grouped by MAT/MF/MT

        Examples
        --------

        Read hydrogen tape from endf-6 formatted file.

        >>> import sandy
        >>> file = "h1.endf"
        >>> tape = sandy.get_endf6_file("jeff_33", "xs", 10010, local=True)
        >>> tape.to_file(file)
        >>> obj = _FormattedFile.from_file(file)
    
        The returned object must be a formatted ENDF-6 file:
    
        >>> assert isinstance(obj, _FormattedFile)
    
        Check that all known keys are present:
    
        >>> keys = obj.data.keys()
        >>> assert len(keys) == 10
        >>> assert (125, 1, 451) in keys
        >>> assert (125, 2, 151) in keys
        >>> assert (125, 3,   1) in keys
        >>> assert (125, 3,   2) in keys
        >>> assert (125, 3, 102) in keys
        >>> assert (125, 4,   2) in keys
        >>> assert (125, 6, 102) in keys
        >>> assert (125,33,   1) in keys
        >>> assert (125,33,   2) in keys
        >>> assert (125,33, 102) in keys
    
        Reading from a text stream must yield the same result:
    
        >>> import io
        >>> stream = io.StringIO(open(file).read())
        >>> obj2 = _FormattedFile.from_file(stream)
    
        >>> assert obj2.data == obj.data

        """
        if isinstance(file, io.StringIO):
            text = file.read()

        else:
            with open(file) as f:
                text = f.read()
        return cls.from_text(text)

    @classmethod
    def from_zipfile(cls, zip_filename):
        """
        Read and return the contents of the only file inside a ZIP archive.
    
        This function assumes the ZIP file contains exactly one file, typically
        with the same base name as the ZIP itself. The internal file is read
        directly from the ZIP archive without extracting it to disk.
        
        Use case: a ENDF-6 zip downloaded from the IAEA website such as `'n_0125_1-H-1.zip'`, 
        which contains a single file `'n_0125_1-H-1.dat'`
    
        Parameters
        ----------
        zip_filename : str
            Path to the ZIP file on disk.
    
        Returns
        -------
        :obj:`~sandy.endf6.Endf6`
            An instance of decoded text content of the internal file.
    
        Raises
        ------
        FileNotFoundError
            If the ZIP file does not exist.
        zipfile.BadZipFile
            If the file is not a valid ZIP archive.
        IndexError
            If the ZIP file is empty and contains no internal files.
        UnicodeDecodeError
            If the internal file cannot be decoded as UTF‑8 text.
    
        Notes
        -----
        - No temporary directories are created.
        - The internal file is read entirely into memory.
        """
        from zipfile import ZipFile

        # Open the ZIP file from disk
        with ZipFile(zip_filename, "r") as z:
            # Expecting exactly one file inside
            internal_name = z.namelist()[0]
    
            # Read and decode as text
            with z.open(internal_name) as f:
                text = f.read().decode("utf-8")
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

        >>> import sandy
        >>> file = "h1.endf"
        >>> tape = sandy.get_endf6_file("jeff_33", "xs", 10010, local=True)
        >>> tape.to_file(file)
        >>> text = open(file).read()
        >>> obj = _FormattedFile.from_text(text)
        
        The returned object must be a formatted ENDF-6 file (tested in `.frome_file`).
        
        Reading the same text with extra empty lines at top and bottom
        should yield identical parsed data:
        
        >>> text_with_empty = "\\n" * 10 + text + "\\n" * 10
        >>> obj2 = _FormattedFile.from_text(text_with_empty)
        >>> assert obj2.data == obj.data
        """

        # -----------------------------
        # 1. Parse MAT/MF/MT using read_fwf
        # -----------------------------
        df = pd.read_fwf(
            io.StringIO(text),
            widths=[66, 4, 2, 3],
            names=["TEXT", "MAT", "MF", "MT"],
            dtype={"TEXT": str, "MAT": str, "MF": int, "MT": int},
            na_filter=False,  # speeds up and does not add NaN in empty lines
            # Do not use TEXT because  the parser does not preserve the whitespaces
            usecols=("MAT", "MF", "MT"),
            )


        # -----------------------------
        # 2. Rebuild TEXT column manually (preserving whitespace)
        # -----------------------------
        # Use splitlines instead of readlines to remove "\n"
        # The if clause removes empty lines.
        df["TEXT"] = [line for line in text.splitlines() if line.split()]


        # -----------------------------
        # 3. Fix title line if MAT is not integer
        # -----------------------------
        title = df["TEXT"].iloc[0]
        title_mat = df["MAT"].iloc[0]

        try:
            int(title_mat)

        except ValueError:
            logging.warning(f"wrong MAT number in the file title\n'{title}'")
            df = df.iloc[1:].reset_index(drop=True)

        finally:
            df["MAT"] = df["MAT"].astype(int)


        # -----------------------------
        # 4. Compute mask using NumPy (avoids pandas ops → avoids NumExpr)
        # -----------------------------
        mt = df["MT"].to_numpy()
        mf = df["MF"].to_numpy()
        mat = df["MAT"].to_numpy()

        mask = (mt > 0) & (mf > 0) & (mat > 0)        # NumPy ops → no pandas.core.ops

        df2 = df.loc[mask, ["MAT", "MF", "MT", "TEXT"]]
    
        # -----------------------------
        # 5. Manual group-by (avoids pandas.groupby → no NumExpr paths)
        # -----------------------------
        data = {}
        # group rows by (MAT, MF, MT)
        for (mat_v, mf_v, mt_v), group in df2.groupby(["MAT", "MF", "MT"], sort=False):
            data[(mat_v, mf_v, mt_v)] = "\n".join(group["TEXT"].tolist())

        return cls(data)

    def _get_section_df(self, mat, mf, mt):
        """
        
        Examples
        --------

        Check if we can read Pu240 file from JEFF-3.1.1, where a "?" is found in the header.

        >>> import sandy
        >>> tape = sandy.get_endf6_file("jeff_311", "xs", 942400, local=True)
        >>> assert "?" in tape.data[(9440, 1, 451)]
        >>> out = tape._get_section_df(9440, 1, 451)
        
        Let's make it fail.

        >>> import pytest
        >>> text = " 94-Pu-240 BRC,CAD    EVAL-JUL04 Bouland Derrien Morillon R?@$¤n  9440 1451    5"
        >>> with pytest.raises(ValueError) as exc_info:
        ...    sandy.Endf6.from_text(text)._get_section_df(9440, 1, 451)

        """
        from .utils import add_delimiter_every_n_characters, add_exp_in_endf6_text

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
            return add_delimiter_every_n_characters(
                x[:66],
                11,
                delimiter=delimiter,
            )

        newtext = "\n".join(map(foo, text.splitlines())).replace('"', '*')
        df = pd.read_csv(
            io.StringIO(add_exp_in_endf6_text(newtext)),
            delimiter=delimiter,
            na_filter=True,
            names=["C1", "C2", "L1", "L2", "N1", "N2"],
        )
        return df

    def add_section(self, mat, mf, mt, text):
        """
        Add or replace a section identified by (MAT, MF, MT) in the underlying
        ENDF-6 container.
    
        The method returns a **new** instance with the updated content; the
        original object is not modified.
    
        Parameters
        ----------
        mat : int
            MAT number.
        mf : int
            MF number.
        mt : int
            MT number.
        text : str
            ENDF-6 section body to store at the given (MAT, MF, MT).

        Returns
        -------
        :obj:`~sandy.endf6._FormattedFile` or derived instance
            object with new section

        Examples
        --------

        Basic add of a new section and structural checks:

        >>> import sandy
        >>> tape = sandy.Endf6({(9437, 3, 102) : "lorem ipsum"})
        >>> new_tape = tape.add_section(9999, 1, 1, "dolor sit amet")


        The returned object must be an Endf6 instance

        >>> assert isinstance(new_tape, sandy.Endf6)
        
        It must contain exactly the two keys

        >>> keys = new_tape.data.keys()
        >>> assert len(keys) == 2
        
        Values must be preserved and correctly inserted

        >>> assert new_tape.data[(9437, 3, 102)] == "lorem ipsum"
        >>> assert new_tape.data[(9999, 1, 1)] == "dolor sit amet"
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
        Delete the section identified by (MAT, MF, MT) from `Endf6.data`.
    
        The method returns a **new** instance with the section removed; the
        original object is not modified.
    
        Parameters
        ----------
        mat : int
            MAT number.
        mf : int
            MF number.
        mt : int
            MT number.
        raise_error : bool, optional
            If True (default), raise a KeyError when the (MAT, MF, MT) section
            does not exist. If False, return the object unchanged when the key
            is absent.
    
        Returns
        -------
        :obj:`sandy.endf6._FormattedFile`
            A new instance (same concrete class as `self`) without the given section.
    
        Raises
        ------
        KeyError
            If the key does not exist and `raise_error=True`.
    
        Examples
        --------

        Delete capture cross section from hydrogen (JEFF-3.3) and verify
        the key is removed while other sections remain:

        >>> import sandy
        >>> tape = sandy.get_endf6_file("jeff_33", "xs", 10010, local=True)
        >>> new = tape.delete_section(125, 3, 102)
    
        The removed key must not be present

        >>> keys = new.data.keys()
        >>> assert (125, 3, 102) not in keys
    
        Some other known sections must still be present
        >>> assert (125, 1, 451) in keys
        >>> assert (125, 6, 102) in keys
        >>> assert len(keys) == 9
    
        If the section is absent and raise_error=False, no exception is raised:
    
        >>> _ = new.delete_section(125, 99, 999, raise_error=False)
    
        If the section is absent and raise_error=True, a KeyError is raised:
    
        >>> try:
        ...     _ = new.delete_section(125, 99, 999, raise_error=True)
        ...     assert False, "Expected KeyError"
        ... except KeyError:
        ...     pass
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

        >>> import sandy
        >>> h1 = sandy.get_endf6_file("jeff_33", 'xs', 10010, local=True)
        >>> h2 = sandy.get_endf6_file("endfb_71", 'xs', 10020, local=True)
        >>> h = h1.merge(h2)
        >>> assert h.to_series()[h1.to_series().index].equals(h1.to_series())
        >>> assert h.to_series()[h2.to_series().index].equals(h2.to_series())

        Merge three files from different libraries.

        >>> h3 = sandy.get_endf6_file("endfb_71", 'xs', 10030, local=True)
        >>> h_ = h1.merge(h2, h3).to_series()
        >>> h__ = h.merge(h3).to_series()
        >>> h___ = h1.merge(h2).merge(h3).to_series()
        >>> assert h_.equals(h__) and h_.equals(h___)

        Merge two evaluations for the same nuclide.

        >>> bi_71 = sandy.get_endf6_file("endfb_71", 'xs', 832090, local=True)
        >>> bi_33 = sandy.get_endf6_file("jeff_33", 'xs', 832090, local=True)
        >>> bi = bi_71.merge(bi_33)
        >>> assert not bi.to_series()[bi_71.to_series().index].equals(bi_71.to_series())
        >>> assert bi.to_series()[bi_33.to_series().index].equals(bi_33.to_series())
        >>> bi = bi_33.merge(bi_71)
        >>> assert bi.to_series()[bi_71.to_series().index].equals(bi_71.to_series())
        >>> assert not bi.to_series()[bi_33.to_series().index].equals(bi_33.to_series())
        """
        from functools import reduce

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

        >>> import sandy
        >>> string = sandy.get_endf6_file("jeff_33", "xs", 10010, local=True).write_string()
        >>> print(string[:81 * 4 - 1])
                                                                             1 0  0    0
         1.001000+3 9.991673-1          0          0          2          5 125 1451    1
         0.000000+0 0.000000+0          0          0          0          6 125 1451    2
         1.000000+0 2.000000+7          3          0         10          3 125 1451    3

        if no modification is applied to the `_FormattedFile` content, the
        `write_string` returns an output identical to the file ASCII content.

        Test with `sandy.Errorr` object and title option.

        >>> endf6 = sandy.get_endf6_file("jeff_33", "xs", 10010, local=True)
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

        >>> tape = sandy.get_endf6_file("jeff_33", "decay", [10010, 10040], local=True)
        >>> assert all(x==80 for x in map(len, tape.write_string().splitlines()))
        """
        from .records import write_line

        string = ""

        # Write title
        if tpid:
            string += write_line(title, 1, 0, 0, 0)
            string += "\n"

        for mat, dfmat in self.to_series().groupby('MAT', sort=True):
            for mf, dfmf in dfmat.groupby('MF', sort=True):
                for mt, text in dfmf.groupby('MT', sort=True):
                    string += text.squeeze()\
                                  .encode('ascii', 'replace')\
                                  .decode('ascii')
                    string += "\n"
                    string += write_line("", mat, mf, 0, 99999)
                    string += "\n"
                string += write_line("", mat, 0, 0, 0)
                string += "\n"
            string += write_line("", 0, 0, 0, 0)
            string += "\n"

        # Write end-of-file
        if fend:
            string += write_line("", -1, 0, 0, 0)
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
    Source data is found in attribute :obj:`~sandy.endf6.Endf6.data`.

    Methods
    -------
    get_ace
        Process :obj:`~sandy.endf6.Endf6` instance into an ACE file using NJOY.
    get_pendf
        Process :obj:`~sandy.endf6.Endf6` instance into a PENDF file using NJOY.
    get_errorr
        Process :obj:`~sandy.endf6.Endf6` instance into a ERRORR file using NJOY.
    get_id
        Extract ID for a given MAT for a ENDF-6 file.
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
        :obj:`~sandy.endf6.Endf6`
            :obj:`~sandy.endf6.Endf6` instance with updated MF1/MT451.

        
        Examples
        --------

        Check how many lines of description and how many sections are recorded
        in a file.

        >>> import sandy
        >>> tape = sandy.get_endf6_file("jeff_33", "xs", 10010, local=True)
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
        from . import write_mf1

        tape = self.data.copy()
        for mat, g in self.to_series().groupby("MAT"):
            intro = self.read_section(mat, 1, 451)
            intro.update(**kwargs)
            new_records = [(mf, mt, sec.count('\n') + 1, 0) for (mat, mf, mt), sec in g.items()]
            NWD, NXC = len(intro["DESCRIPTION"]), g.shape[0]
            new_records[0] = (1, 451, NWD+NXC+4, 0)
            intro["SECTIONS"] = new_records
            tape[(mat, 1, 451)] = write_mf1(intro)
        return self.__class__(tape)

    def read_section(self, mat, mf, mt, raise_error=True):
        """
        Parse MAT/MF/MT section.

        Parameters
        ----------
        mat : `int`
            MAT number
        mf : `int`
            MF number
        mt : `int`
            MT number
        raise_error : `bool`, optional
            Raise or not error if section is not found. The default is True.

        Returns
        -------
        `dict`
        """
        from importlib import import_module

        modname = f"sandy.sections.mf{mf}"

        try:
            module = import_module(modname)  # e.g. sandy.sections.mf1
        except ModuleNotFoundError:
            if raise_error:
                raise ValueError(f"Unsupported MF={mf}")
            return None

        func_name = f"read_mf{mf}"

        try:
            reader = getattr(module, func_name)
        except AttributeError:
            if raise_error:
                raise ValueError(
                    f"Module '{modname}' does not define '{func_name}'"
                )
            return None

        return reader(self, mat, mt)

    def _update_info(self, descr=None):
        """
        Update RECORDS item (in DATA column) for MF1/MT451 of each MAT based on the content of the TEXT column.
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
            
            - If `method='aleph'` the ID is the ZAM identifier.
            - Else, the ID is the ZA identifier according to the NNDC rules.

        Returns
        -------
        ID : `int`
            ID of the ENDF-6 file.

        Notes
        -----
        .. note:: A warning is raised if more than one MAT is found.
                  Only the ID corresponding to the lowest MAT will be returned.
 
        Examples
        --------

        Extract ID for H1 file using NNDC and ALEPH methods.

        >>> import sandy
        >>> tape = sandy.get_endf6_file("jeff_33", "xs", 10010, local=True)
        >>> assert tape.get_id() == 1001
        >>> assert tape.get_id(method="aleph") == 10010

        Extract ID for Am242m file using NNDC and ALEPH methods.

        >>> tape2 = sandy.get_endf6_file("jeff_33", "xs", 952421, local=True)
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

    def _run_njoy(self, pendf=None, pendftape=None, **njoy_kws):
        """
        Internal helper: run NJOY on this ENDF6 tape.
    
        Handles:
        - writing ENDF6 to tmpdir
        - writing optional PENDF
        - preparing arguments for process_neutron
        - returning outputs
        
        Parameters
        ----------
        pendf : :obj:`sandy.endf6.Endf6`, optional, default is `None`.
            Endf6 object containing the pendf file (treated as a python object)
        pendftape : `str`, optional, default is `None`.
            filename (with path) of the pendf file (in this case it is read from file)

        Notes
        -----
        Keyword argument `pendf` is used to pass aPENDF as `Endf6` object,
        while `pendftape` is used to pass a PENDF as the name of a file written on disk.
        """
        from tempfile import TemporaryDirectory
        from .njoy import process_neutron

        with TemporaryDirectory() as td:
            # 1. Write ENDF6 main tape to temp folder
            endf6file = join(td, "tape20")
            self.to_file(endf6file)

            # 2. Handle optional PENDF input
            pendf_file = None

            if pendf is not None:
                # Used passed PENDF object (used by get_pendf)
                if pendf.kind != "pendf":
                    raise TypeError("'pendf' must contain a PENDF tape")
                pendf_file = join(td, "tape21")
                pendf.to_file(pendf_file)

            elif pendftape is not None:
                # User passed filename + path to an existing PENDF file
                pendf_file = pendftape

            # 3. Run NJOY through sandy
            outputs = process_neutron(
                endf6file,
                pendftape=pendf_file,
                **njoy_kws,
            )

        # 4. Return NJOY output (dict)
        return outputs
    
    def _prepare_groupr_kws(self, **groupr_kws):
        """Helper to prepare groupr options"""
        # -- prepare/augment GROUPR options without mutating the user's dict --
        groupr_kws_ = dict(groupr_kws)
    
        # Decide what GROUPR should process based on available sections
        recs = self.get_records()
        has_fission_xs = 18 in recs.query("MF==3").MT.values   # XS fiss present?
        groupr_kws_.setdefault("nubar", has_fission_xs)        # fission XS implies nubar processing
        has_fission_xs = 18 in recs.query("MF==5").MT.values   # PFNS present?
        groupr_kws_.setdefault("chi",   has_fission_xs)        # PFNS implies chi processing
        groupr_kws_.setdefault("mubar", True)                  # always include mubar

        return groupr_kws_
    
    def _prepare_njoy_kws(
            self,
            temperature=0,
            err=0.001,
            minimal_processing=False,
            verbose=False,
            reconr_kws=None,
            broadr_kws=None,
            thermr_kws=None,
            **njoy_kws,
            ):
        """
        Prepare NJOY keyword arguments for neutron processing.
    
        This helper consolidates module-level keyword arguments, applies implicit
        defaults, and enforces the consistent interpretation of `temperature`,
        `err`, and `minimal_processing`. All user-provided dictionaries are copied
        before modification to avoid mutating caller input.
    
        Parameters
        ----------
        temperature : float, optional
            Target processing temperature in Kelvin. If set to ``0``:
            - BROADR is disabled,
            - all post-BROADR modules are disabled,
            - NJOY processing stops after RECONR.
        err : float, optional
            Reconstruction tolerance for RECONR, BROADR, and THERMR. This value is
            applied only to modules whose keyword dictionaries do *not* already
            define ``"err"`` explicitly.
        minimal_processing : bool, optional
            If True, disable all modules after BROADR (THERMR, GASPR, HEATR,
            PURR, UNRESR).
        verbose : bool, optional
            If True, instruct the NJOY driver to print the generated input deck
            before execution.
        reconr_kws, broadr_kws, thermr_kws : dict or None, optional
            Keyword arguments for RECONR, BROADR, and THERMR respectively.
            If None, an empty dict is assumed. Caller dictionaries are never
            mutated; copies are always created.
        **njoy_kws :
            Additional keyword arguments to forward to the NJOY driver.
    
        Returns
        -------
        dict
            A fully prepared and internally consistent NJOY keyword dictionary,
            containing:
            - module activation flags (possibly modified),
            - `reconr_kws`, `broadr_kws`, `thermr_kws`,
            - global NJOY settings for temperature and verbosity.
    
        Notes
        -----
        - Caller-provided dictionaries are always copied before modification.
        - User-specified ``"err"`` values always override the default.

        Examples
        --------
        Test that `minimal_processing` filters unwanted modules.

        >>> import sandy
        >>> g = sandy.get_endf6_file("jeff_33", "xs", 10010, local=True).get_gendf(err=1, minimal_processing=True, temperature=300, dryrun=True)
        >>> assert "broadr" in g and "reconr" in g
        >>> assert "thermr" not in g and "purr" not in g and "heatr" not in g and "unresr" not in g and "gaspr" not in g

        Test `minimal_processing=False`.

        >>> g = sandy.get_endf6_file("jeff_33", "xs", 10010, local=True).get_gendf(err=1, temperature=300, dryrun=True)
        >>> assert "broadr" in g and "reconr" in g
        >>> assert "thermr" in g and "purr" in g and "heatr" in g and "gaspr" in g

        Check that for `temperature=0` the calculation stops after RECONR.

        >>> g = sandy.get_endf6_file("jeff_33", "xs", 10010, local=True).get_gendf(err=1, dryrun=True)
        >>> assert "reconr" in g
        >>> assert "broadr" not in g and "thermr" not in g and "purr" not in g and "heatr" not in g and "unresr" not in g and "gaspr" not in g
        """
        njoy_kws_ = dict(njoy_kws)  # clean copy

        # minimal processing
        if minimal_processing or float(temperature) == 0:
            njoy_kws_.update(dict(
                thermr=False, gaspr=False, heatr=False, purr=False, unresr=False
            ))
        
        # stop after reconr
        if temperature == 0:
            njoy_kws_["broadr"] = False
            logging.warning(
                "Zero Kelvin requested; NJOY will stop after RECONR. "
                "Use temperature=0.1 for 0K xs processing."
            )
        
        # --- Prepare submodule keyword dicts (no mutation of caller dicts) ---
        reconr_kws_ = (reconr_kws or {}).copy()
        broadr_kws_ = (broadr_kws or {}).copy()
        thermr_kws_ = (thermr_kws or {}).copy()

        # Only apply default 'err' if user did not specify one
        if "err" not in reconr_kws_:
            reconr_kws_["err"] = float(err)
        if "err" not in broadr_kws_:
            broadr_kws_["err"] = float(err)
        if "err" not in thermr_kws_:
            thermr_kws_["err"] = float(err)
        
        njoy_kws_["reconr_kws"] = reconr_kws_
        njoy_kws_["broadr_kws"] = broadr_kws_
        njoy_kws_["thermr_kws"] = thermr_kws_

    
        # --- NJOY global arguments ---
        njoy_kws_["temperatures"] = [temperature]
        njoy_kws_["verbose"] = verbose

        return njoy_kws_

    def get_ace(self, suffix=None, pendf=None, dryrun=False, **njoy_kws):
        """
        Process :obj:`~sandy.endf6.Endf6` instance into an ACE file using NJOY.

        Parameters
        ----------
        dryrun : `bool`, optional
            Do not run NJOY and return only NJOY input. Default is `False`.
        pendf : :obj:`~sandy.endf6.Endf6`, optional
            Provide PENDF object and add it to the processing
            sequence after RECONR and before BROADR. Default is `None`.
        suffix : `str`, optional
            suffix in the form `".[0-9][0-9]"` to assign to the ACE data.
            If not given, generate automatic suffix according to ALEPH rules.
            Default is `None`.
        **kwargs : `dict`
            keyword argument to pass to :obj:`~sandy.njoy.process_neutron`.

        Returns
        -------
        `dict` of `str`
            output with `'ace'` and `'xsdir'` as keys.

        Examples
        --------
        Check that output is a ace file.

        >>> import sandy
        >>> e6 = sandy.get_endf6_file("jeff_33", "xs", 10010, local=True)
        >>> outs = e6.get_ace(temperature=700, err=1, minimal_processing=True)
        >>> ace, xsdir = outs["ace"], outs["xsdir"]
        >>> assert "1001.07c" in ace
        >>> assert "sandy runs acer" in ace
        >>> assert "mat 125" in ace
        
        Check that ace is processed at a different temperature.

        >>> outs = e6.get_ace(temperature=800, err=1, minimal_processing=True)
        >>> ace, xsdir = outs["ace"], outs["xsdir"]
        >>> assert "1001.08c" in ace

        Check xsdir (last entry is a newline).

        >>> assert xsdir[:-1] == '  1001.08c    0.999167 filename 0 1   1     3085     0     0 6.894E-08'

        Check that ace file follows "nndc" nomenclature for metastable nuclides.

        >>> e6_m = sandy.get_endf6_file("jeff_33", "xs", 521291, local=True)
        >>> ace_m = e6_m.get_ace(temperature=800, err=1, minimal_processing=True)["ace"]
        >>> assert "52529.08c" in ace_m

        Check that using option `pendf` results in the same output.

        >>> pendf = e6.get_pendf(temperature=0, err=1)
        >>> ace2 = e6.get_ace(temperature=800, err=1, minimal_processing=True, pendf=pendf)["ace"]
        >>> assert ace == ace2

        Check that the option suffix is used correctly.

        >>> ace = e6.get_ace(temperature=800, suffix=".85", err=1)["ace"]
        >>> assert "1001.85c" in ace

        Check input pendf file

        >>> import pytest
        >>> with pytest.raises(Exception) as e_info:
        ...    e6.get_ace(pendf=e6)
        """
        # avoid mutating the user-supplied kwargs
        njoy_kws_ = self._prepare_njoy_kws(**njoy_kws) | {"dryrun": dryrun, "pendf": pendf}
        njoy_kws_ |= ({"suffixes": [suffix]} if suffix is not None else {})
    
        # --- run via the shared helper ---
        outputs = self._run_njoy(**njoy_kws_)
    
        if dryrun:
            return outputs
    
        return {"ace": outputs["ace"], "xsdir": outputs["xsdir"]}

    def get_pendf(self, dryrun=False, **njoy_kws):
        """
        Process :obj:`~sandy.endf6.Endf6` instance into a PENDF file using NJOY.

        Parameters
        ----------
        dryrun : bool, optional
            If True (default), return the NJOY input deck instead of executing NJOY.
        **njoy_kws : dict
            Additional keyword arguments forwarded to :obj:`~sandy.njoy.process_neutron`.

        Returns
        -------
        Endf6 or dict
        - If `dryrun=True`: returns the NJOY input deck (text blocks).
        - If `dryrun=False`: returns a :obj:`~sandy.endf6.Endf6` PENDF object.

        Examples
        --------

        Default run.

        >>> import sandy
        >>> endf6 = sandy.get_endf6_file("jeff_33", "xs", 10010, local=True)
        >>> out = endf6.get_pendf(verbose=True, temperature=293.6, err=1, minimal_processing=True)
        >>> assert isinstance(out, sandy.Endf6)
        """
        # always deactivate acer, but avoid mutating the user-supplied kwargs
        njoy_kws_ = self._prepare_njoy_kws(**njoy_kws) | {"dryrun": dryrun, "acer": False}
    
        # --- run via the shared helper ---
        outputs = self._run_njoy(**njoy_kws_)
                
        # --- In case of dryrun, 'outputs' contains the text of the NJOY input
        if dryrun:
            return outputs
    
        return Endf6.from_text(outputs["pendf"])

    def get_gendf(self, dryrun=False, groupr_kws=None, **njoy_kws):
        """
        Process the current ENDF‑6 evaluation into a multi‑group GENDF using NJOY.
    
        This method prepares and runs the NJOY processing sequence with GROUPR
        to generate a GENDF representation of the evaluation. It augments
        (without overwriting user input) the GROUPR options to include:
        - ``mubar`` (always enabled)
        - ``nubar`` if fission XS (MF=3/MT=18) are present
        - ``chi``   if PFNS data (MF=5/MT=18) are present
    
        Internally, the method:
        1) writes the ENDF‑6 tape to a temporary working directory,
        2) activates GROUPR (``groupr=True``) and disables ACE (``acer=False``),
        3) forwards all remaining options to :func:`sandy.njoy.process_neutron`,
        4) returns the NJOY deck (when ``dryrun=True``) or parses and returns a
           :class:`sandy.Gendf` object (when ``dryrun=False``).
    
        Parameters
        ----------
        dryrun : bool, optional
            If ``True``, return the generated NJOY input deck (text) 
            instead of executing NJOY. Default is ``False``.
        groupr_kws : dict or None, optional
            Dictionary of keyword arguments passed directly to the NJOY GROUPR
            module (e.g., ``iwt``, ``ign``, ``ek``, ``mt``, ``sigz``).  
            User‑provided values are preserved; defaults are inserted only if not present.
            Default is ``None`
        **njoy_kws : dict
            Additional keyword arguments forwarded to
            :func:`~sandy.njoy.process_neutron`.
            User‑provided values are preserved.
    
        Returns
        -------
        :obj:`~sandy.gendf.Gendf` or str
            - If ``dryrun=True``: text with the NJOY input deck.
            - If ``dryrun=False``: a :class:`~sandy.gendf.Gendf` instance created
              from the produced GENDF text.
    
        Notes
        -----
        - GROUPR is always enabled (``groupr=True``).
        - ACE production is always disabled here (``acer=False``).
        - The working directory is temporary and removed after the run.
        - To control GROUPR specifics (e.g., energy grid, iwt/ign/sigz/ek/mt),
          supply them under ``groupr_kws`` in ``**njoy_kws``.

        Examples
        --------

        Default run.

        >>> import sandy
        >>> endf6 = sandy.get_endf6_file("jeff_33", "xs", 10010, local=True)
        >>> out = endf6.get_gendf(temperature=293.6, minimal_processing=True)
        >>> assert isinstance(out, sandy.Gendf)

        Test keyword `sigz`.

        >>> out = endf6.get_gendf(groupr_kws=dict(sigz=[1e10, 1e2]))
        >>> assert 1e10 in sandy.gendf.read_mf1(out, 125)['SIGZ']
        >>> assert 1e10 in sandy.gendf.read_mf1(out, 125)['SIGZ']

        Test keyword `iwt`.

        >>> import re
        >>> g = endf6.get_gendf(groupr_kws=dict(iwt=3), dryrun=True)
        >>> found = re.search('groupr(.*)moder', g, flags=re.DOTALL).group().splitlines()
        >>> assert "125 2 0 3 0 1 1 0 /" == found[2]

        Test keyword `ign`.

        >>> g = endf6.get_gendf(groupr_kws=dict(ign=3), dryrun=True)
        >>> found = re.search('groupr(.*)moder', g, flags=re.DOTALL).group().splitlines()
        >>> assert "125 3 0 2 0 1 1 0 /" == found[2]

        Test keyword `ek`.

        >>> g = endf6.get_gendf(groupr_kws=dict(ek=sandy.energy_grids.CASMO12), dryrun=True)
        >>> found = re.search('groupr(.*)moder', g, flags=re.DOTALL).group().splitlines()
        >>> ek = np.array(list(map(float, found[7].replace("/", "").split())))
        >>> np.testing.assert_allclose(ek, sandy.energy_grids.CASMO12)

        Test groupr MFs and MTs for fissile and non-fissile nuclides.

        >>> g = endf6.get_gendf(dryrun=True)
        >>> found = re.search('groupr(.*)moder', g, flags=re.DOTALL).group().splitlines()
        >>> assert " ".join(found[6:10]) == '3/ 3 251 / 0/ 0/'

        U-238 test because it contains mubar, xs, chi and nubar.

        >>> endf6 = sandy.get_endf6_file('jeff_33','xs', 922380, local=True)
        >>> g = endf6.get_gendf(dryrun=True)
        >>> found = re.search('groupr(.*)moder', g, flags=re.DOTALL).group().splitlines()
        >>> assert " ".join(found[6:15]) == '3/ 3 452 / 3 455 / 3 456 / 3 251 / 5/ 5 18 / 0/ 0/'

        Test custom MTs.

        >>> endf6 = sandy.get_endf6_file("jeff_33", "xs", 10010, local=True)
        >>> g = endf6.get_gendf(dryrun=True, groupr_kws=dict(mt=4))
        >>> found = re.search('groupr(.*)moder', g, flags=re.DOTALL).group().splitlines()
        >>> assert " ".join(found[6:10]) == '3 4 / 3 251 / 0/ 0/'
        >>> g = endf6.get_gendf(dryrun=True, groupr_kws=dict(mt=[4, 102]))
        >>> found = re.search('groupr(.*)moder', g, flags=re.DOTALL).group().splitlines()
        >>> assert " ".join(found[6:11]) == '3 4 / 3 102 / 3 251 / 0/ 0/'
        """
        from .gendf import Gendf

        # --- start from a clean copy, never mutate the caller's dict ---
        # Always activate GROUPR and never run ACER when producing GENDF
        njoy_kws_ = self._prepare_njoy_kws(**njoy_kws) | {"groupr": True, "acer": False} 
    
        # -- prepare/augment GROUPR options without mutating the user's dict --
        groupr_kws_ = (groupr_kws or {}).copy()
        njoy_kws_["groupr_kws"] = self._prepare_groupr_kws(**groupr_kws_)
    
        # Pass dryrun policy down to NJOY
        njoy_kws_["dryrun"] = dryrun
    
        # --- run via the shared helper ---
        outputs = self._run_njoy(**njoy_kws_)

        # --- In case of dryrun, 'outputs' contains the text of the NJOY input
        if dryrun:
            return outputs
    
        # Parse GENDF text into object
        return Gendf.from_text(outputs["gendf"])

    def get_errorr(
            self,
            nubar=None,  # None means "auto"; bool means user override
            mubar=None,
            chi=None,
            xs=None,
            dryrun=False,
            groupr_kws=None,
            errorr_kws=None,
            errorr31_kws=None,
            errorr33_kws=None,
            errorr34_kws=None,
            errorr35_kws=None,
            **njoy_kws,
            ):
        """
        Generate NJOY ERRORR covariance processing for this ENDF6 evaluation.
    
        This method configures and runs the NJOY ERRORR module (optionally via
        GROUPR when required) to produce covariance matrices for cross sections,
        angular distributions, chi, and/or nubar.
    
        Parameters
        ----------
        nubar : bool or None, optional
            Whether to process MF=31 (nubar covariance).
            - None (default): automatically enable if MF=31 exists.
        mubar : bool or None, optional
            Whether to process MF=34 (mubar covariance).
            - None (default): automatically enable if MF=34 exists.
        chi : bool or None, optional
            Whether to process MF=35 (chi covariance).
            - None (default): automatically enable if MF=35 exists.
        xs : bool or None, optional
            Whether to process MF=33 (cross‑section covariance).
            - None (default): automatically enable if MF=33 exists.
        dryrun : bool, optional
            If True, return the generated NJOY input deck as text instead of
            running NJOY. Useful for testing.
        groupr_kws : dict, optional
            Keyword arguments passed to the GROUPR module (for nubar, chi, mubar).
            User-provided values override automatic defaults.
        errorr_kws : dict, optional
            Base keyword options applied to all ERRORR submodules, unless
            overridden by the specific `errorrXY_kws`.
        errorrXY_kws : dict, optional
            Options passed to individual ERRORR modules:
            - errorr31_kws → MF=31  (nubar)
            - errorr33_kws → MF=33  (cross sections)
            - errorr34_kws → MF=34  (mubar)
            - errorr35_kws → MF=35  (chi)
        **njoy_kws :
            Additional keyword arguments forwarded to the NJOY processing pipeline.
    
        Returns
        -------
        dict or str
            - If `dryrun=True`: return the NJOY input deck as a string.
            - If `dryrun=False`: return a dictionary mapping module names
              ("errorr31", "errorr33", ...) to :obj:`~sandy.errorr.Errorr` objects.
              If no covariance is found/reqeusted, returns empty dict.
    
        Notes
        -----
        - Automatic selection of ERRORR submodules is based on the presence of
          MF=31/33/34/35 in the ENDF file, unless explicitly overridden.
        - GROUPR input flags (nubar/chi/mubar) are prepared via
          `_prepare_groupr_kws`.

        Examples
        --------

        Default run.

        >>> import sandy
        >>> endf6 = sandy.get_endf6_file("jeff_33", "xs", 942410, local=True)
        >>> out = endf6.get_errorr(temperature=300, minimal_processing=True, err=1, errorr_kws=dict(ign=3, mt=18))

        Check `ign` and `ek`. This test checks also the type of each output.

        >>> assert out["errorr33"].get_xs().data.shape[0] == 30
        >>> assert out["errorr31"].get_xs().data.shape[0] == 30
        >>> assert out["errorr34"].get_xs().data.shape[0] == 30
        >>> assert out["errorr33"].get_xs().data.shape[0] == 30

        # Check `ign` and `ek`.

        # >>> endf6 = sandy.get_endf6_file("jeff_33", "xs", 10010, local=True)
        # >>> out = endf6.get_errorr(errorr_kws=dict(ek=sandy.energy_grids.CASMO12))

        Check `mt`.
        >>> assert out["errorr33"].get_xs().data.squeeze().name == (9443, 18)
        >>> assert out["errorr34"].get_xs().data.squeeze().name == (9443, 251)
        >>> columns = out["errorr31"].get_xs().data.columns
        >>> assert (9443, 452) in columns and (9443, 455) in columns and (9443, 456) in columns

        Check consistency between keywords `errorr_kws` and `errorr33_kws`.

        >>> ekws = dict(irespr=0, iwt=5, ek=[1e-5, 2e7], mt=(16, 18, 102))
        >>> e6 = sandy.get_endf6_file("jeff_33", "xs", 942410, local=True)
        >>> inp1 = e6.get_errorr(temperature=300, dryrun=True, xs=True, chi=False, nubar=False, mubar=False, errorr_kws=ekws)
        >>> inp2 = e6.get_errorr(temperature=300, dryrun=True, xs=True, chi=False, nubar=False, mubar=False, errorr33_kws=ekws)
        >>> inp3 = e6.get_errorr(temperature=300, dryrun=True, xs=True, chi=False, nubar=False, mubar=False)
        >>> assert "groupr" not in inp1 and "groupr" not in inp2 and "groupr" not in inp3
        >>> assert inp1 == inp2 and inp1 != inp3

        Check consistency between keywords `errorr_kws` and `errorr35_kws`.

        >>> inp1 = e6.get_errorr(temperature=300, dryrun=True, xs=False, chi=True, nubar=False, mubar=False, errorr_kws=ekws)
        >>> inp2 = e6.get_errorr(temperature=300, dryrun=True, xs=False, chi=True, nubar=False, mubar=False, errorr35_kws=ekws)
        >>> inp3 = e6.get_errorr(temperature=300, dryrun=True, xs=False, chi=True, nubar=False, mubar=False)
        >>> assert "groupr" in inp1 and "groupr" in inp2 and "groupr" in inp3
        >>> assert inp1 == inp2 and inp1 != inp3

        Check consistency between keywords `errorr_kws` and `errorr31_kws`.

        >>> inp1 = e6.get_errorr(temperature=300, dryrun=True, xs=False, chi=False, nubar=True, mubar=False, errorr_kws=ekws)
        >>> inp2 = e6.get_errorr(temperature=300, dryrun=True, xs=False, chi=False, nubar=True, mubar=False, errorr31_kws=ekws)
        >>> inp3 = e6.get_errorr(temperature=300, dryrun=True, xs=False, chi=False, nubar=True, mubar=False)
        >>> assert inp1 == inp2 and inp1 != inp3
        >>> assert "groupr" in inp1 and "groupr" in inp2 and "groupr" in inp3

        Check consistency between keywords `errorr_kws` and `errorr34_kws`.

        >>> inp1 = e6.get_errorr(temperature=300, dryrun=True, xs=False, chi=False, nubar=False, mubar=True, errorr_kws=ekws)
        >>> inp2 = e6.get_errorr(temperature=300, dryrun=True, xs=False, chi=False, nubar=False, mubar=True, errorr34_kws=ekws)
        >>> inp3 = e6.get_errorr(temperature=300, dryrun=True, xs=False, chi=False, nubar=False, mubar=True)
        >>> assert inp1 == inp2 and inp1 != inp3
        >>> assert "groupr" in inp1 and "groupr" in inp2 and "groupr" in inp3
        >>> inp1 = e6.get_errorr(temperature=300, dryrun=True, errorr_kws=ekws)
        >>> inp2 = e6.get_errorr(temperature=300, dryrun=True, errorr33_kws=ekws, errorr31_kws=ekws, errorr34_kws=ekws, errorr35_kws=ekws)
        >>> assert inp1 == inp2
        >>> assert "groupr" in inp1 and "groupr" in inp2

        Check default options.

        >>> import re
        >>> g = sandy.get_endf6_file("jeff_33", "xs", 10010, local=True).get_errorr(temperature=300, dryrun=True)
        >>> found = re.search('errorr(.*)', g, flags=re.DOTALL).group().splitlines()

        Check ign(2), iwt (2), iprint (0) and relative (1) options.

        >>> assert found[2] == '125 2 2 0 1 /'

        Check temperature (300) option.

        >>> assert found[3] == '0 300.0 /'

        Check irespr (1) option.

        >>> assert found[4] == '0 33 1/'

        Check options changes.
        
        >>> ekws = dict(ign=3, iwt=5, iprint=True, relative=False)
        >>> g = sandy.get_endf6_file("jeff_33", "xs", 10010, local=True).get_errorr(temperature=400, errorr_kws=ekws, dryrun=True)
        >>> found = re.search('errorr(.*)', g, flags=re.DOTALL).group().splitlines()
        >>> assert found[2] == '125 3 5 1 0 /'
        >>> assert found[3] == '1 400.0 /'
        >>> assert found[4] == '0 33 1/'

        Test spectrum.

        >>> spect = [1.000000e-5, 2.00000000, 3.000000e-2, 2.00000000, 5.800000e-2, 4.00000000, 3, 1]
        >>> out = endf6.get_errorr(verbose=False, dryrun=True, errorr_kws={"spectrum": spect, "ek": [1.000000e-5, 3.000000e-2, 5.800000e-2, 3]}, nubar=False, chi=False, mubar=False)
        >>> # this test has to be finished

        Example: 91-Pa-231 in JEFF-3.3 has MF32 but not MF33.

        >>> tape = sandy.get_endf6_file("jeff_33", "xs", 912310, local=True)
        >>> err = tape.get_errorr(chi=False, nubar=True, mubar=False, err=1, xs=True, errorr33_kws=dict(irespr=0))
        >>> assert "errorr33" in err

        Example: 91-Pa-231 in JEFF-3.3 has MF32 but not MF33.

        >>> tape = sandy.get_endf6_file("jeff_33", "xs", 912330, local=True)
        >>> err = tape.get_errorr(chi=False, nubar=True, mubar=False, err=1, xs=True, errorr33_kws=dict(irespr=0))
        >>> assert "errorr33" in err

        Example: 17-Cl-37 in JEFF-3.3 has MF32 but not MF33.

        >>> tape = sandy.get_endf6_file("jeff_33", "xs", 170370, local=True)
        >>> err = tape.get_errorr(chi=False, nubar=True, mubar=False, err=1, xs=True, errorr33_kws=dict(irespr=0))
        >>> assert "errorr33" in err

        Example: 95-Am-241 in JEFF-3.3 has MF32 but not MF33.
        
        >>> tape = sandy.get_endf6_file("jeff_33", "xs", 952410, local=True)
        >>> err = tape.get_errorr(chi=False, nubar=True, mubar=False, err=1, xs=True, errorr33_kws=dict(irespr=0))
        >>> assert "errorr33" in err
        """
        from tempfile import TemporaryDirectory

        from .njoy import _input_mf32_nomf33, _input_mf32_nomf33_no18, _run_njoy
        from .errorr import Errorr

        src = self  # change of variables to avoid overwriting self
        # ------------------------------------------------------------------
        # 0. Handle MF32-no-MF33 cases (old @handle_mf32_alone decorator)
        # ------------------------------------------------------------------
        recs = src.get_records()
        if 32 in recs.MF.values and 33 not in recs.MF.values:
            # input taken from
            # https://www-nds.iaea.org/index-meeting-crp/TM_NDP/docs/OCabellos_2017.pdf
            input_fiss = _input_mf32_nomf33
            input_nofiss = _input_mf32_nomf33_no18

            # choose MF32 template depending on whether fission exists
            inp = input_fiss if 18 in recs.MT.values else input_nofiss

            with TemporaryDirectory() as td:
                f20 = os.path.join(td, "tape20")
                src.to_file(f20)
                outs = _run_njoy(inp, f20)

            # Replace self with synthetic ERRORR33-only partial tape
            src = Endf6.from_text(outs["errorr33"])
            recs = src.get_records()


        # --- start from a clean copy, never mutate the caller's dict ---
        njoy_kws_ = src._prepare_njoy_kws(**njoy_kws) | {"dryrun": dryrun, "acer": False}
    
        # -- prepare/augment GROUPR options without mutating the user's dict --
        groupr_kws_ = (groupr_kws or {}).copy()
        njoy_kws_["groupr_kws"] = src._prepare_groupr_kws(**groupr_kws_)

        # -- prepare/augment ERRORR options without mutating the user's dict --
        errorr_kws_ = (errorr_kws or {}).copy()

        # Activate specific errorr module according to covariance info and input options
        has31 = not recs.loc[recs.MF == 31].empty  # nubar cov
        has33 = not recs.loc[recs.MF == 33].empty  # xs cov
        has34 = not recs.loc[recs.MF == 34].empty  # mubar cov
        has35 = not recs.loc[recs.MF == 35].empty  # chi cov

        # Switch off if user provides False
        use31 = has31 if nubar is None else has31 & bool(nubar)
        use33 = has33 if xs    is None else has33 & bool(xs)
        use34 = has34 if mubar is None else has34 & bool(mubar)
        use35 = has35 if chi   is None else has35 & bool(chi)

        # Fan out shared base kwargs without overriding user-provided per-module, rightmost wins
        errorr31_kws_ = errorr_kws_ | (errorr31_kws or {}).copy()
        errorr33_kws_ = errorr_kws_ | (errorr33_kws or {}).copy()
        errorr34_kws_ = errorr_kws_ | (errorr34_kws or {}).copy()
        errorr35_kws_ = errorr_kws_ | (errorr35_kws or {}).copy()

        # Compose process_neutron kwargs (no mutation of caller input)
        njoy_kws_ |= {
            "errorr31": use31,
            "errorr33": use33,
            "errorr34": use34,
            "errorr35": use35,
            "errorr31_kws": errorr31_kws_,
            "errorr33_kws": errorr33_kws_,
            "errorr34_kws": errorr34_kws_,
            "errorr35_kws": errorr35_kws_,
        }
        
        # --- run via the shared helper ---
        outputs = src._run_njoy(**njoy_kws_)

        # --- In case of dryrun, 'outputs' contains the text of the NJOY input
        if dryrun:
            return outputs

        outputs = {k: Errorr.from_text(v) for k, v in outputs.items() if k.startswith("errorr")}
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

        Short test for hydrogen.

        >>> import sandy
        >>> sandy.get_endf6_file("jeff_33", "xs", 10010, local=True).get_records()
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

    def get_perturbations(self, *args, **kwargs,):
        """
        Dispatcher to assign perturbations method: either for rdaioactive
        decay data or cross section.

        Notes
        -----
        .. note :: The perturbation method is selected based on the MT's found
                   in `self`.
        """
        logging.info("########################################################")
        logging.info("                GET PERTURBATIONS                       ")
        logging.info("########################################################")
        # this could have been a decorator...
        if 457 in self.mt:
            out = self.get_perturbations_rdd(*args, **kwargs)

        elif 454 in self.mt:
            out = self.get_perturbations_fy(*args, **kwargs)

        else:
            out = self.get_perturbations_xs(*args, **kwargs)

        return out

    def get_perturbations_xs(self, nsmp, njoy_kws={}, smp_kws={}, **kwargs,):
        """
        Construct multivariate distributions with a unit vector for 
        mean and with relative covariances taken from the evaluated files
        processed with the NJOY module ERRORR.

        Perturbation factors are sampled with the same multigroup structure of 
        the covariance matrix and are returned by nuclear datatype as a `dict`.

        Parameters
        ----------
        nsmp : `int`
            Sample size.
        njoy_kws : `dict`, optional
            Keyword arguments to produce ERRORR file.
            The default is {}.
        smp_kws : `dict`, optional
            Keyword arguments for :obj:`~sandy.cov.CategoryCov.sampling`.
            The default is {}.
        **kwargs : `dict`
            additional keyword arguments.

        Returns
        -------
        smp : `dict` of :obj:`~sandy.samples.Samples`
            Dictionary with sample objects.
            The dictionary keys are `31` and `33`, respectively for cross
            sections and nubar.

        Examples
        --------

        Generate a couple of samples from the H1 file of JEFF-3.3.

        >>> import sandy
        >>> njoy_kws = dict(err=1, errorr_kws=dict(mt=102))
        >>> tape = sandy.get_endf6_file("jeff_33", "xs", 10010, local=True)
        >>> smps = tape.get_perturbations(nsmp=2, njoy_kws=njoy_kws)
        >>> assert len(smps) == 1
        >>> assert isinstance(smps[33], sandy.Samples)
        >>> assert (smps[33].data.index.get_level_values("MT") == 102).all()

        Test get perturbations from MF 35.

        >>> njoy_kws = dict(err=1, errorr_kws=dict(mt=18))
        >>> tape = sandy.get_endf6_file("jeff_33", "xs", 922350, local=True)
        >>> smps = tape.get_perturbations(nsmp=2, njoy_kws=njoy_kws)
        >>> assert len(smps) == 3
        >>> assert isinstance(smps[35], sandy.Samples)
        >>> assert (smps[35].data.index.get_level_values("MT") == 18).all()
        """
        from .samples import summarize_sample

        smp = {}
        
        debug = kwargs.get("verbose", False)

        # -- produce ERRORR files with covariance data
        logging.info(" - Produce ERRORR file with NJOY...")
        njoy_kws["mubar"] = False
        outs = self.get_errorr(**njoy_kws)

        filename = "PERT_{}_MF{}.xlsx"
        filename_err = "ERRORR_{}_MF{}.tape"

        # -- Extract samples from covariance data, iterate over MF31, 33 and 35
        for k, out in outs.items():

            # -- Get MF from keys of get_errorr ouput dictionary
            mf = int(k[-2:])
            logging.info(f" - Processing covariance matrix for MF={mf}...")
            
            # -- Print ERRORR tape to file
            if debug:
                xls = filename_err.format(self.get_id(), mf)
                logging.info(f" - Writing ERRORR file to '{xls}'...")
                out.to_file(xls)

            # -- Extract covariance matrix
            cov = out.get_cov()

            # -- Extract sample
            seed = smp_kws.get(f"seed{mf}")
            smp[mf] = cov.sampling(nsmp, seed=seed, **smp_kws)

            # -- Dump sample to file
            if debug:
                xls = filename.format(self.get_id(), mf)
                logging.info(f" - Writing perturbation file '{xls}'...")
                smp[mf].to_excel(xls)
                # cov.to_excel(xls)
            
                # Write sample and cov stats to Excel
                with pd.ExcelWriter(xls, mode="a", engine="openpyxl", if_sheet_exists="replace") as writer:
                    
                    if nsmp > 1:
                        summary = summarize_sample(smp[mf], cov)
                    else:
                        # don't call the sample summary for nmsp=1 to avoid warnings
                        summary = {}
    
                    df = pd.Series(summary, name="summary")
                    df.to_excel(writer, sheet_name="STATS SMP")

                summary = cov.summarize()
                df = pd.Series(summary, name="summary")
                df.to_excel(writer, sheet_name="STATS COV")

        return smp

    def get_perturbations_rdd(self, nsmp, smp_hl_kws={}, smp_de_kws={}, smp_br_kws={}, fill_zeros=0.05, **kwargs,):
        """
        Construct multivariate distributions with a unit vector for  mean and
        with relative covariances taken from the evaluated radioactive decay
        data files in `self`.

        Perturbation factors are sampled for the decay constants, decay
        energies, and branching ratios and are returned as a `dict`.

        Parameters
        ----------
        nsmp : `int`
            Sample size.
        smp_hl_kws : `dict`, optional
            Keyword arguments for :obj:`~sandy.cov.CategoryCov.sampling` for half-lives.
            The default is {}.
        smp_de_kws : `dict`, optional
            Keyword arguments for :obj:`~sandy.cov.CategoryCov.sampling` for decay energies.
            The default is {}.
        smp_br_kws : `dict`, optional
            Keyword arguments for :obj:`~sandy.cov.CategoryCov.sampling` for branching ratios.
            The default is {}.
        fill_zeros : `float`, optional
            Some decay energy and half-life data carry zero uncertainty in the evaluation.
            This option allows setting a default uncertainty for all such cases,
            e.g., `fill_zeros=0.05` adds a 5% uncertainty to all such cases.
            The default is 0.05.
            Set `fill_zeros=0.0` if you don't want to consider any additional
            uncertainty.
        **kwargs : `dict`
            additional keyword arguments, such as:
                - `rdd`: to pass directly an already processed :obj:`~sandy.decay.DecayData` instance.
                - `verbose`: to activate output verbosity.
        Raises
        ------
        ValueError
            Error if all nuclides are stable. Then there is no variance.

        Returns
        -------
        smp : `dict` of :obj:`~sandy.samples.Samples`
            Dictionary with sample objects.
            The dictionary keys are `'HL'`, `'DE'` and `'BR'`, respectively for
            half-lives, decay energies and bramching ratios.

        
        Notes
        -----
        .. note:: branching ratios are sampled without correlations and must be
                  renormalized.
        .. note:: if branching ratios have zero unceratinty, unit perturbation
                  coefficients are assigned by default.

        Examples
        --------

        Branching ratio coefficients are one if branching ratio uncertainty is
        not given.

        >>> import sandy
        >>> smps = sandy.get_endf6_file("jeff_33", "decay", [10010, 10040, 270600], local=True).get_perturbations(2)
        >>> assert (smps["BR"].data.values == 1).all()

        Raise error if all nuclides are stable.

        >>> import pytest
        >>> with pytest.raises(ValueError) as exc_info:
        ...    sandy.get_endf6_file("jeff_33", "decay", 10010, local=True).get_perturbations(2)
        """
        from .decay import DecayData
        from .cov import CategoryCov
        from .samples import Samples

        debug = kwargs.get("verbose", False)

        # if already available in kwargs, do not extract DecayData again
        rdd = kwargs.get("rdd")
        if not rdd:
            rdd = DecayData.from_endf6(self, verbose=kwargs.get("verbose"))

        #If all nuclides are stable, then there is no variance and CategoryCov fails
        if all([v["stable"] for v in rdd.data.values()]):
            raise ValueError("this method does not work if all nuclides are stable")
            

        # Convert to relative units and replace zero uncertainties with standard value, e.g., 5%
        # --- decay constants ---
        hl = rdd.get_half_life()
        dhl = (hl.data.DHL / hl.data.HL).replace(0, fill_zeros).fillna(0)
        smp_hl = CategoryCov.from_stdev(dhl).sampling(nsmp, **smp_hl_kws)
        
        # --- decay energies ---
        de = rdd.get_decay_energy()
        dde = (de.data.DE / de.data.E).replace(0, fill_zeros).fillna(0)
        smp_de = CategoryCov.from_stdev(dde).sampling(nsmp, **smp_de_kws)
        
        # --- branching ratios ---
        br = rdd.get_branching_ratio()
        # if branching ratios don't have uncertainty, use constant coefficient 1
        if not br.data.DBR.any():
            smp_br = Samples(np.ones([br.data.shape[0], nsmp]), index=br.data.index)
        else:
            dbr = (br.data.DBR / br.data.BR).fillna(0)
            smp_br = CategoryCov.from_stdev(dbr).sampling(nsmp, **smp_br_kws)
        
        if debug:
            xlsx_file = 'PERT_MF8_MT457.xlsx'
            logging.info(f"writing to file '{xlsx_file}'...")
            with pd.ExcelWriter(xlsx_file, engine="openpyxl") as writer:
                smp_hl.data.to_excel(writer, sheet_name='HALF LIFE')
                smp_de.data.to_excel(writer, sheet_name='DECAY ENERGY')
                smp_br.data.to_excel(writer, sheet_name='BRANCHING RATIO')

        smp = {
            "BR": smp_br,
            "DE": smp_de,
            "HL": smp_hl,
            }
        return smp

    def get_perturbations_fy(self, nsmp, smp_kws={}, covariance=None, **kwargs,):
        """
        Construct multivariate distributions with a unit vector for  mean and
        with relative covariances taken from the evaluated fission yield
        data files in `self`.

        Perturbation factors are sampled for the independent fission yields only.        

        Parameters
        ----------
        nsmp : `int`
            Sample size.
        smp_kws : `dict`, optional
            Keyword arguments for :obj:`~sandy.cov.CategoryCov.sampling`.
            The default is {}.
        covariance : `None` or `str`, optional
            Flag to adopt fission yield covariance matrices.
            The only acceptable flag is `covariance='cea'`, which uses
            the covariance evaluation for U-235 and Pu-239 produced by CEA for
            thermal fission yields.
            See :obj:`~sandy.fy.get_cea_fy`.
            The default is `None`.
        **kwargs : `dict`
            Not used.

        Returns
        -------
        smps : `pd.DataFrame`
            Dataframe with perturbation coefficients given per:
                
                - ZAM: fissioning nuclide
                - E: neutron energy
                - ZAP: fission product
                - SMP: sample ID
            
            .. note:: This is different from :obj:`~sandy.endf6.Endf6.get_perturbations_xs`
                      and :obj:`~sandy.endf6.Endf6.get_perturbations_rdd`, which return
                      a :obj:`~sandy.samples.Samples` instance.

        Examples
        --------
        
        Default use case.

        >>> import sandy
        >>> tape = sandy.get_endf6_file("jeff_33", "nfpy", 922350, local=True)
        >>> smps = tape.get_perturbations_fy(2, smp_kws=dict(seed=3))

        Pass already processed fission yield object.

        >>> nfpy = sandy.Fy.from_endf6(tape)
        
        Ensure reproducibility by fixing seed.
        
        >>> smps2 = tape.get_perturbations_fy(2, nfpy=nfpy, smp_kws=dict(seed=3))
        >>> assert smps.equals(smps2)
        
        Test `covariance='cea'` option.
        This is done by checking the sample correlation between nuclides
        `zap=451140` and `461140`, which in the source data is larger than 0.9.

        >>> smps = tape.get_perturbations_fy(50, nfpy=nfpy, covariance=None)
        >>> data = smps.query("ZAP in [451140, 461140] & E==0.0253").pivot_table(index="ZAP", columns="SMP", values="VALS")
        >>> assert np.corrcoef(data)[0, 1] < 0.3
        >>> smps = tape.get_perturbations_fy(50, nfpy=nfpy, covariance='cea')
        >>> data = smps.query("ZAP in [451140, 461140] & E==0.0253").pivot_table(index="ZAP", columns="SMP", values="VALS")
        >>> assert np.corrcoef(data)[0, 1] > 0.9

        """
        import random
        from .cov import CategoryCov              # lazy import to avoid circular import issue
        from .fy import Fy, get_cea_fy           # lazy import to avoid circular import issue
        
        debug = kwargs.get("verbose", False)

        # if already available in kwargs, do not extract fission yields again
        nfpy = kwargs.get("nfpy")
        if not nfpy:
            nfpy = Fy.from_endf6(self, verbose=kwargs.get("verbose"))

        # if seed is given in "smp_kws", it ensures reproducibility
        seed_start = smp_kws.get("seed", random.randrange(2**32 - 1))
        # set the seed that will be used to ensure the same seed generation sequence when calling CategoryCov.sampling
        random.seed(seed_start)

        smps = []
        for (zam, e), fy in nfpy.data.query("MT==454").groupby(["ZAM", "E"]):
            
            if covariance == "cea" and zam in [922350, 942390] and e==0.0253:
                fy, rcov = get_cea_fy(zam)
                
            else:
                rstd = (fy.DFY / fy.FY).fillna(0)  # relative uncertainties
                rcov =  CategoryCov(pd.DataFrame(np.diag(rstd**2), index=fy.ZAP, columns=fy.ZAP))

            # this is a Samples instance, I cannot pass a seed because it would be used for all fissioning systems
            smp = rcov.sampling(nsmp, seed=random.randrange(2**32 - 1))
            # this is not a Samples instance anymore
            smp = (
                smp.data.rename_axis(index="ZAP").
                stack().rename("VALS").reset_index().
                assign(E=e, ZAM=zam)[["ZAM", "E", "ZAP", "SMP", "VALS"]]  # add energy and ZAM and sort keys
                )
            smps.append(smp)

        # stack with all samples for all ZAM, energy and ZAP
        smps = pd.concat(smps, ignore_index=True)

        if debug:
            xlsx_file = 'PERT_MF8_MT454.xlsx'
            logging.info(f"writing to file '{xlsx_file}'...")
            with pd.ExcelWriter(xlsx_file) as writer:
                for zam, smp in smps.groupby("ZAM"):
                    smp.pivot_table(index=["E", "ZAP"], columns="SMP", values="VALS").to_excel(writer, sheet_name=f"{zam}")

        return smps

    def apply_perturbations(self, *args, **kwargs,):
        """
        Dispatcher to assign perturbations method: either for radioactive
        decay data or cross section.

        Notes
        -----
        .. note :: The perturbation method is selected based on the MT's found
                   in `self`.

        Examples
        --------

        The next two examples will mismatch file and samples.
        The output must be `None` if samples and ENDF6 file do not match.

        >>> import sandy
        >>> tape = sandy.get_endf6_file("jeff_33", "xs", 10010, local=True)
        >>> taped = sandy.get_endf6_file("jeff_33", "decay", 10040, local=True)
        
        Mix rdd file with xs samples.

        >>> pendf = tape.get_pendf(err=1)
        >>> smps = tape.get_perturbations(3)
        >>> assert not taped.apply_perturbations(smps, pendf=pendf)

        Mix xs file with rdd samples.

        >>> rdd = sandy.DecayData.from_endf6(taped)
        >>> smps = taped.get_perturbations(2, rdd=rdd)
        >>> assert not tape.apply_perturbations(smps, rdd=rdd)
        """
        logging.info("########################################################")
        logging.info("              APPLY PERTURBATIONS                       ")
        logging.info("########################################################")

        # this could have been a decorator...same as get_perturbations
        if 457 in self.mt:
            out = self.apply_perturbations_rdd(*args, **kwargs)

        elif 454 in self.mt:
            out = self.apply_perturbations_fy(*args, **kwargs)

        else:
            out = self.apply_perturbations_xs(*args, **kwargs)
        return out

    def apply_perturbations_xs(self, smps, processes=1, pendf=None, njoy_kws={}, **kwargs):
        """
        Apply relative perturbations to the cross-section (XS), nubar and
        prompt fission neutron spectrum (chi) data in an ENDF6 file.
    
        This method perturbs reaction cross sections and nubar values based on provided 
        perturbation samples. The process can be performed in parallel for efficiency. 
        If a PENDF file is not provided, it will be generated automatically.
    
        Parameters
        ----------
        smps : dict of :obj:`~sandy.samples.Samples`
            Dictionary containing relative perturbation coefficients for XS and nubar.
            Expected keys:
            - `31`: nubar perturbations
            - `33`: cross-section perturbations
            - `35`: chi perturbations
        processes : int, optional, default=1
            Number of parallel processes. If `processes > 1`, perturbations are applied in parallel.
        pendf : :obj:`~sandy.endf6.Endf6`, optional, default=None
            If provided, perturbations are applied to this PENDF file. 
            Otherwise, a new PENDF file is generated from `self` (more time-consuming).
        njoy_kws : dict, optional
            Dictionary of keyword arguments for `sandy.endf6.Endf6.get_pendf`, 
            used to generate a PENDF file if `pendf` is not provided.
        **kwargs : dict, optional
            Additional options for ACE file generation and arguments passed to :obj:`~sandy.endf6._endf6_perturb_worker`.
    
        Returns
        -------
        dict
            A dictionary (indexed by sample ID) of perturbed ENDF/PENDF files
            or ACE files, depending on `to_ace` and `to_file` options.
            - If `to_file=False` and `to_ace=False`: Returns a dictionary of `sandy.Endf6` objects (both `'endf6'` and `'pendf'`).
            - If `to_file=True`: Saves perturbed files to disk and returns filenames.
            - If `to_ace=True`: Generates ACE files and returns filenames.
    
        Notes
        -----
        - **Temperature Treatment**:
            - By default, perturbations are applied to a 0K PENDF, followed by Doppler broadening.
            - Alternatively, perturbations can be applied directly to a temperature-specific PENDF.
        - **Parallelization**:
            - If `processes=1`, perturbations are applied sequentially.
            - If `processes>1`, a `ProcessPoolExecutor` is used for parallel processing.
        - **Supported Perturbations**:
            - Nubar (`pnu`, MT=31)
            - Cross-sections (`pxs`, MF=3)
            - Chi (`pchi`, MF=5)
    
        Examples
        --------

        Apply perturbations to Pu-239 XS and nubar.
    
        >>> import sandy
        >>> tape = sandy.get_endf6_file("jeff_33", "xs", 942390, local=True)
        >>> smps = tape.get_perturbations(
        ...     2, 
        ...     njoy_kws={"err": 1, "chi": False, "mubar": False, "errorr33_kws": {"mt": [2, 4, 18]}}, 
        ...     smp_kws={"seed31": 1, "seed33": 3}
        ... )
    
        Apply both nubar and XS perturbations.
    
        >>> outs_31_33 = tape.apply_perturbations_xs(smps, njoy_kws={"err": 1}, processes=1)
    
        Apply only nubar perturbations.
    
        >>> outs_31 = tape.apply_perturbations_xs({31: smps[31]}, njoy_kws={"err": 1}, processes=1)
    
        Apply only XS perturbations.
    
        >>> outs_33 = tape.apply_perturbations_xs({33: smps[33]}, njoy_kws={"err": 1}, processes=1)
    
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

        Default use case for chi and xs, perturbed together.

        >>> smps_ = tape.get_perturbations(2, njoy_kws=dict(err=1, nubar=False, mubar=False), smp_kws=dict(seed33=3, seed35=5))
        >>> outs_33_35 = tape.apply_perturbations(smps_, njoy_kws=dict(err=1), processes=1)

        Compare to individual xs and chi perturbations with same seed.

        >>> outs_33_ = tape.apply_perturbations({33: smps_[33]}, njoy_kws=dict(err=1), processes=1)
        >>> outs_35 = tape.apply_perturbations({35: smps_[35]}, njoy_kws=dict(err=1), processes=1)
        
        >>> for i in range(2):
        ...    assert(outs_33_[i]["endf6"].data == tape.data)
        ...    assert(outs_35[i]["endf6"].data != tape.data)
        ...    assert(outs_35[i]["endf6"].data == outs_33_35[i]["endf6"].data)
        ...    assert(outs_33_[i]["pendf"].data != outs_35[i]["pendf"].data)
        ...    assert(outs_33_[i]["pendf"].data == outs_33_35[i]["pendf"].data)

        >>> assert outs_33_[0]["pendf"].data != outs_33_[1]["pendf"].data
        >>> assert outs_33_[0]["endf6"].data == outs_33_[1]["endf6"].data
        >>> assert outs_35[0]["pendf"].data == outs_35[1]["pendf"].data
        >>> assert outs_35[0]["endf6"].data != outs_35[1]["endf6"].data

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

        Check that ENDF6 and PENDF output filenames are correct.

        >>> endf6 = sandy.get_endf6_file('jeff_33', 'xs', 10010, local=True)
        >>> smps = endf6.get_perturbations(2, njoy_kws=dict(err=0.1))
        >>> outs = endf6.apply_perturbations(smps, to_file=True)
        >>> assert outs[0]["endf6"] == '1001_0.endf6' and os.path.isfile('1001_0.endf6')
        >>> assert outs[0]["pendf"] == '1001_0.pendf' and os.path.isfile('1001_0.endf6')
        >>> assert outs[1]["endf6"] == '1001_1.endf6' and os.path.isfile('1001_1.endf6')
        >>> assert outs[1]["pendf"] == '1001_1.pendf' and os.path.isfile('1001_1.pendf')

        Check that ACE output filenames are correct.

        >>> outs = endf6.apply_perturbations(smps, to_file=True, to_ace=True, ace_kws=dict(err=1, temperature=300, purr=False, heatr=False, thermr=False, gaspr=False))
        >>> assert outs[0]["ace"] == '1001_0.03c' and os.path.isfile('1001_0.03c')
        >>> assert outs[0]["xsdir"] == '1001_0.03c.xsd' and os.path.isfile('1001_0.03c.xsd')
        >>> assert outs[1]["ace"] == '1001_1.03c' and os.path.isfile('1001_1.03c')
        >>> assert outs[1]["xsdir"] == '1001_1.03c.xsd' and os.path.isfile('1001_1.03c.xsd')
        
        Check that keyword `pendf` works.

        >>> pendf = endf6.get_pendf(err=1)
        >>> outs1 = endf6.apply_perturbations(smps, njoy_kws=dict(err=1))
        >>> outs2 = endf6.apply_perturbations(smps, pendf=pendf)
        >>> assert outs1[0]["pendf"].write_string() == outs2[0]["pendf"].write_string()

        """
        from concurrent.futures import ProcessPoolExecutor, as_completed

        if 33 not in smps and 31 not in smps and 35 not in smps:
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
        if 35 in smps:
            # At this level, the samples have the same index
            # as a covariance matrix and can be treated as xs
            data["pchi"] = smps[35].iterate_xs_samples()

        # This dict indexed by sample will contain the output of the worker:
        #    - either perturbed endf6 and pendf tape as `Endf6` objects
        #    - or ace files as string
        outs = {}

        if processes == 1:

            logging.info(" - Apply XS perturbations in series...")

            while True:
                kws = {}
                # -- Iterate perturbation data (xs, nubar, chi)
                for k, v in data.items():
                    try:
                        n, s = next(v)  # Get the next (index, sample) from generator
                        kws[k] = s  # Store sample with generator key
                    except StopIteration:
                        break  # Exit the loop immediately if any generator is exhausted
        
                else:  # Only executes if `for` loop completes normally (no `break`)
                    kws.update(kwargs)  # Merge static kwargs

                    # Call `endf6_perturb_worker` directly
                    outs[n] = _endf6_perturb_worker(self.data, pendf_.data, n, **kws)
                    continue  # Continue to the next iteration
        
                break  # If any generator is exhausted, exit the while loop

        elif processes > 1:
            # switched from mp.Pool to ProcessPoolExecutor because compatible with windows

            logging.info(f" - Apply XS perturbations using a pool of {processes} workers...")

            with ProcessPoolExecutor(max_workers=processes) as executor:
                futures = {}
            
                while True:
                    kws = {}
                    # -- Iterate perturbation data (xs, nubar, chi)
                    for k, v in data.items():
                        try:
                            n, s = next(v)  # Get the next (index, sample) from generator
                            kws[k] = s  # Store sample with generator key
                        except StopIteration:
                            break  # Exit the loop immediately if any generator is exhausted
            
                    else:  # Only executes if `for` loop completes normally (no `break`)
                        kws.update(kwargs)  # Merge static kwargs
            
                        # Submit to the process pool, calling `endf6_perturb_worker` directly
                        futures[executor.submit(_endf6_perturb_worker, self.data, pendf_.data, n, **kws)] = n
                        continue  # Continue to the next iteration
            
                    break  # If any generator is exhausted, exit the while loop
            
                for future in as_completed(futures):
                    n = futures[future]
                    try:
                        outs[n] = future.result()
                    except Exception as e:
                        print(f"Error in task {n}: {e}")  # Error handling
        
        # if we keep ENDF6 and PENDF files in memory, convert them back into
        # sandy Endf6 instances (must do it here because Endf6 object cannot be pickled)
        if not kwargs.get("to_file", False) and not kwargs.get("to_ace", False):
            outs = {k: {k1: Endf6(v1) for k1, v1 in v.items()} for k, v in outs.items()}

        return outs

    def apply_perturbations_rdd(self, smps, processes=1, **kwargs):
        """
        Apply relative perturbations to the data contained in
        :obj:`~sandy.endf6.Endf6` instance of radioactive decay data files.

        Parameters
        ----------
        smps : `dict` of :obj:`~sandy.samples.Samples`
            Dictionary with sample objects.
            See output of :obj:`~sandy.endf6.Endf6.get_perturbations_rdd`
        processes : `int`, optional, default is `1`
            Number of processes used to complete the task.
            Creation of ENDF6 files and post-processing is done in parallel if
            `processes>1`.
        **kwargs : `dict`
            Additional keyword arguments, such as:
                - `rdd`: to pass directly an already processed :obj:`~sandy.decay.DecayData` instance.
                - `verbose`: to activate output verbosity.
                - `to_file`: to write output :obj:`~sandy.endf6.Endf6` instances to file.

        Returns
        -------
        outs : `dict` of :obj:`~sandy.endf6.Endf6` or `dict` of `str`
            Depending on whether keyword argument `to_file` is given or not:
                - `to_file=True`: `dict` with filenames, sample ID's are keys
                - `to_file=False`: `dict` with :obj:`~sandy.endf6.Endf6` instances, sample ID's are keys

        Notes
        -----
        .. note :: if `to_file=True`, outputs have names `'decay_data_0'`, `'decay_data_1'`, etc.

        Examples
        --------

        >>> import sandy
        >>> tape = sandy.get_endf6_file("jeff_33", "decay", [10040, 270590, 270600, 571380], local=True)
        >>> rdd = sandy.DecayData.from_endf6(tape)
        >>> smps = tape.get_perturbations(2, rdd=rdd)
        >>> outs = tape.apply_perturbations_rdd(smps, rdd=rdd)
        >>> rdd0 = sandy.DecayData.from_endf6(outs[0])
        
        Check that half-lives are correctly perturbed.

        >>> np.testing.assert_almost_equal(rdd0.data[10040]['half_life'] / rdd.data[10040]['half_life'], smps["HL"].data.loc[10040, 0], decimal=5)
        >>> np.testing.assert_almost_equal(rdd0.data[270600]['half_life'] / rdd.data[270600]['half_life'], smps["HL"].data.loc[270600, 0], decimal=5)
        >>> np.testing.assert_almost_equal(rdd0.data[571380]['half_life'] / rdd.data[571380]['half_life'], smps["HL"].data.loc[571380, 0], decimal=5)
        >>> assert rdd0.data[270590]['half_life'] == rdd.data[270590]['half_life'] == 0   # this is stable
        
        Check that decay constants are also correctly recalculated.

        >>> np.testing.assert_almost_equal(rdd.data[10040]['decay_constant'] / rdd0.data[10040]['decay_constant'], smps["HL"].data.loc[10040, 0], decimal=5)
        >>> np.testing.assert_almost_equal(rdd.data[270600]['decay_constant'] / rdd0.data[270600]['decay_constant'], smps["HL"].data.loc[270600, 0], decimal=5)
        >>> np.testing.assert_almost_equal(rdd.data[571380]['decay_constant'] / rdd0.data[571380]['decay_constant'], smps["HL"].data.loc[571380, 0], decimal=5)
        >>> assert rdd0.data[270590]['decay_constant'] == rdd.data[270590]['decay_constant'] == 0   # this is stable
        
        Check that decay energies are also correctly perturbed.

        >>> np.testing.assert_almost_equal(rdd0.data[270600]['decay_energy']["beta"] / rdd.data[270600]['decay_energy']["beta"], smps["DE"].data.loc[(270600, "beta"), 0], decimal=5)
        >>> np.testing.assert_almost_equal(rdd0.data[270600]['decay_energy']["gamma"] / rdd.data[270600]['decay_energy']["gamma"], smps["DE"].data.loc[(270600, "gamma"), 0], decimal=5)
        >>> assert rdd0.data[270600]['decay_energy']["alpha"] == rdd.data[270600]['decay_energy']["alpha"] == 0
        
        Check that parameters with zero uncertainty are affected by `fill_zeros`.

        >>> np.testing.assert_almost_equal(rdd0.data[10040]['decay_energy']["alpha"] / rdd.data[10040]['decay_energy']["alpha"], smps["DE"].data.loc[(10040, "alpha"), 0], decimal=5)
        >>> smps = tape.get_perturbations(2, rdd=rdd, fill_zeros=0)
        >>> outs = tape.apply_perturbations_rdd(smps, rdd=rdd)
        >>> rdd1 = sandy.DecayData.from_endf6(outs[0])
        >>> assert rdd1.data[10040]['decay_energy']["alpha"] == rdd.data[10040]['decay_energy']["alpha"] != 0

        Check that data is written to file with key `to_file`.

        >>> outs = tape.apply_perturbations_rdd(smps, rdd=rdd, to_file=True)
        >>> assert os.path.exists(outs[0])
        """
        import multiprocessing as mp
        from .decay import DecayData

        # --- PRE-PROCESSING
        if not {"HL", "DE", "BR"}.issubset(smps.keys()):
            logging.info("no (or incomplete) perturbation coefficient was found.")
            return

        # Get nominal decay data. pop it or it will be given twice to _rdd_perturb_worker
        rdd = kwargs.pop("rdd", None)
        if not rdd:
            rdd = DecayData.from_endf6(self, verbose=kwargs.get("verbose"))

        # --- PROCESSING
        if processes == 1:
            outs = {}
            # iterate over columns of samples instance, then samples size is known and key matching is guaranteed
            # assume HL, BR and DE have all the same key matching (sample ids)
            for ismp in smps["HL"].data.columns:
               outs[ismp] = _rdd_perturb_worker(
                   self.data,
                   rdd.data,
                   smps["HL"].data,
                   smps["DE"].data,
                   smps["BR"].data,
                   ismp,
                   **kwargs,
                   )

        elif processes > 1:
            pool = mp.Pool(processes=processes)
            outs = {}
            # iterate over columns of samples instance, then samples size is known and key matching is guaranteed
            # assume HL, BR and DE have all the same key matching (sample ids)
            for ismp in smps["HL"].data.columns:
                outs[ismp] = pool.apply_async(
                    _rdd_perturb_worker,
                    (
                        self.data,
                        rdd.data,
                        smps["HL"].data,
                        smps["DE"].data,
                        smps["BR"].data,
                        ismp,
                        ),
                    kwargs,
                    )
            outs = {n: out.get() for n, out in outs.items()}
            pool.close()
            pool.join()

        # --- POST-PROCESSING
        # if we keep ENDF6 files in memory, convert them back into Endf6 instances
        # (must do it here because Endf6 object cannot be pickled)
        if not kwargs.get("to_file"):
            outs = {k: Endf6(v) for k, v in outs.items()}

        return outs

    def apply_perturbations_fy(self, smps, processes=1, covariance=None, **kwargs):
        """
        Apply relative perturbations to the data contained in
        :obj:`~sandy.endf6.Endf6` instance of fission yield files.

        Parameters
        ----------
        smps : `pd.DataFrame`
            Fission yield sample object.
            See output of :obj:`~sandy.endf6.Endf6.get_perturbations_fy`.
            See also :obj:`~sandy`
        processes : `int`, optional, default is `1`
            Number of processes used to complete the task.
            Creation of ENDF6 files and post-processing is done in parallel if
            `processes>1`.
        covariance : `None` or `str`, optional
            Flag to adopt fission yield covariance matrices.
            The only acceptable flag is `covariance='cea'`.
            This ensures that the nominal values for U-235 and Pu-239 are taken from 
            the CEA evaluations. It must be used if it aws used for the production of `smps`.
            The default is `None`.
        **kwargs : `dict`
            Additional keyword arguments, such as:
                - `nfpy`: to pass directly an already processed :obj:`~sandy.fy.Fy` instance.
                - `verbose`: to activate output verbosity.
                - `to_file`: to write output :obj:`~sandy.endf6.Endf6` instances to file.

        Returns
        -------
        outs : `dict` of :obj:`~sandy.endf6.Endf6` or `dict` of `str`
            Depending on whether keyword argument `to_file` is given or not:
                - `to_file=True`: `dict` with filenames, sample ID's are keys
                - `to_file=False`: `dict` with :obj:`~sandy.endf6.Endf6` instances, sample ID's are keys

        Notes
        -----
        .. note :: if `to_file=True`, outputs have names `'fy_0'`, `'fy_1'`, etc.

        Examples
        --------
        
        Default use case (write data to file).

        >>> import sandy, pytest
        >>> tape = sandy.get_endf6_file("jeff_33", "nfpy", [922350, 922380, 942390], local=True)
        >>> smps = tape.get_perturbations(2, covariance='cea')
        >>> outs = tape.apply_perturbations_fy(smps, covariance='cea', verbose=False, to_file=True)
        
        If the samples were produced with keyword `covariance='cea'`, the same must be 
        used in `apply_perturbations_fy`.

        >>> nfpy = sandy.Fy.from_endf6(tape)
        >>> nfpy_u235 = sandy.Fy.from_endf6(sandy.Endf6.from_file(sandy.fy_cea_u235th))
        >>> nfpy0 = sandy.Fy.from_endf6(sandy.Endf6.from_file(outs[0]))
        
        >>> n = nfpy.data.query("ZAM==922350 and MT==454")
        >>> n0 = nfpy0.data.query("ZAM==922350 and MT==454")
        >>> nu235 = nfpy_u235.data.query("ZAM==922350 and MT==454")

        To match the perturbation values, the ratio must be taken with respect to the 
        CEA nominal values.
        
        >>> sp = n0.set_index("ZAP").FY.divide(nu235.set_index("ZAP").FY).fillna(1)
        >>> p = smps.query("ZAM==922350 and SMP==0").set_index("ZAP").VALS.rename("FY")
        >>> np.testing.assert_array_almost_equal(p, sp, decimal=4)
        
        If `covariance='cea'` was used to produce the samples, at it is not used in 
        `apply_perturbations_fy`, then there is a mismatch between the ZAP numbers of 
        the samples and of the fission yields.

        >>> with pytest.raises(Exception):
        ...    tape.apply_perturbations_fy(smps, verbose=False, to_file=True)
        """
        import multiprocessing as mp
        from .fy import Fy, fy_cea_u235th, fy_cea_pu239th

        # --- PRE-PROCESSING
        # Get nominal fission yield data. pop it or it will be given twice to _fy_perturb_worker
        nfpy = kwargs.pop("nfpy", None)
        if not nfpy:
            nfpy = Fy.from_endf6(self, verbose=kwargs.get("verbose"))
        
        # Change nominal values to CEA values if asked
        # this is needed to ensure that the samples are given for the same nominal values
        if covariance == 'cea':
            tape_u235 = Endf6.from_file(fy_cea_u235th).data if 922350 in nfpy.data.ZAM.values else {}
            tape_pu239 = Endf6.from_file(fy_cea_pu239th).data if 942390 in nfpy.data.ZAM.values else {}
            tape = Endf6({**self.data, **tape_u235, **tape_pu239})
            
            # also the fission yields must be re-extracted
            nfpy = Fy.from_endf6(tape, verbose=kwargs.get("verbose"))

        else:
            tape = self
            

        # --- PROCESSING
        if processes == 1:
            outs = {}
            # iterate over columns of samples instance, then samples size is known and key matching is guaranteed
            # assume HL, BR and DE have all the same key matching (sample ids)
            for ismp, smp in smps.groupby("SMP"):
               outs[ismp] = _fy_perturb_worker(
                   tape.data,
                   nfpy.data,
                   smp,
                   ismp,
                   **kwargs,
                   )

        elif processes > 1:
            pool = mp.Pool(processes=processes)
            outs = {}
            # iterate over columns of samples instance, then samples size is known and key matching is guaranteed
            # assume HL, BR and DE have all the same key matching (sample ids)
            for ismp, smp in smps.groupby("SMP"):
                outs[ismp] = pool.apply_async(
                    _fy_perturb_worker,
                    (
                        tape.data,
                        nfpy.data,
                        smp,
                        ismp,
                        ),
                    kwargs,
                    )
            outs = {n: out.get() for n, out in outs.items()}
            pool.close()
            pool.join()

        # --- POST-PROCESSING
        # if we keep ENDF6 files in memory, convert them back into Endf6 instances
        # (must do it here because Endf6 object cannot be pickled)
        if not kwargs.get("to_file"):
            outs = {k: Endf6(v) for k, v in outs.items()}

        return outs


def _endf6_perturb_worker(
        endf6,
        pendf,
        ismp,
        pxs=None,
        pnu=None,
        plpc=None,
        pchi=None,
        verbose=False,
        to_ace=False,
        to_file=False,
        ace_kws=None,
        **kwargs,
        ):

    """
    Worker to handle ENDF6 neutron data perturbation (xs, nubar, angular and energy distributions).

    Parameters
    ----------
    endf6 : `dict`
        `data` attribute of :obj:`~sandy.endf6.Endf6`.
        It contains the nominal ENDF6 data.
    pendf : `dict`
        `data` attribute of :obj:`~sandy.endf6.Endf6`.
        It contains the nominal PENDF data.
    ismp : `int`
        sample ID.
    pxs : `pd.DataFrame`
        It contains the perturbation coefficients for cross section.
        It corresponds to one single sample (in principle the one with ID `ismp`).
        It should have the same structure as a :obj:`~sandy.xs.Xs` object.
        The default is `None`.
    pnu: `pd.DataFrame`
        It contains the perturbation coefficients for nubar.
        It corresponds to one single sample (in principle the one with ID `ismp`).
        It should have the same structure as a :obj:`~sandy.xs.Xs` object.
        The default is `None`.
    plpc: `pd.DataFrame`
        Not implemented.
    pchi: `pd.DataFrame`
        It contains the perturbation coefficients for chi.
        It corresponds to one single sample (in principle the one with ID `ismp`).
        It should have the same structure as a :obj:`~sandy.xs.Xs` object.
        The default is `None`.
    verbose : `bool`, optional
        Flag to activate verbosity. The default is `False`.
    to_ace : TYPE, optional
        DESCRIPTION. The default is False.
    to_file : `bool`, optional
        Flag to write outputs to file. The default is `False`.
        This key changes the output type.
    ace_kws : `dict` , optional
        Additional keyword arguments for ACE file production. The default is {}.
    **kwargs : `dict`
        Additional keyword arguments (not used).

    Returns
    -------
    `dict`
        
        - if `to_file=False`: a `dict` with keys, values:
            
            - `endf6`: a perturbed :obj:`~sandy.endf6.Endf6` instance of the given ENDF6
            - `pendf`: a perturbed :obj:`~sandy.endf6.Endf6` instance of the given PENDF

        - if `to_file=True`: a `dict` with keys, values:

            - `endf6`: the filename of the perturbed ENDF6
            - `pendf`: the filenmae of the perturbed PENDF

    Examples
    --------

    Test that energy distributions are correctly perturbed.
    Example for Pu239.

    >>> import sandy
    
    Creation of dummy perturbation.

    >>> interval = pd.Interval(left=1e-8, right=10, closed="right")
    >>> idx = pd.MultiIndex.from_tuples([(9437, 18, interval)], names=("MAT", "MT", "E"))
    >>> pert = 1.1
    >>> df = pd.DataFrame([[pert]], index=idx).reset_index()
    >>> smps = sandy.Samples(df.set_index(["MAT", "MT", "E"]))
    >>> smps
    SMP                             0
    MAT  MT E                        
    9437 18 (1e-08, 10.0] 1.10000e+00

    Creation of reference ENDF6 and PENDF.

    >>> ref_endf6 = sandy.get_endf6_file("jeff_33", "xs", 942390, local=True)
    >>> ref_pendf = ref_endf6.get_pendf(err=1)

    Creation of perturbed data modifying the PFNS with the first perturbation sample.
    
    >>> ismp = 0
    >>> perturbed = sandy.endf6._endf6_perturb_worker(ref_endf6.data, ref_pendf.data, ismp, pchi=dict(smps.iterate_xs_samples())[ismp])
    >>> pert_endf6 = sandy.Endf6(perturbed['endf6'])
    >>> pert_pendf = sandy.Endf6(perturbed['pendf'])

    Compare reference and perturbed :obj:`~sandy.edistr.Edistr`.

    >>> ref_edistr = sandy.Edistr.from_endf6(ref_endf6)
    >>> pert_edistr = sandy.Edistr.from_endf6(pert_endf6)

    Test that the perturbation is correct and happened below `ethresh` only.
    This test works for all incident energies.

    >>> np.testing.assert_array_almost_equal(pert_edistr.data.query("EOUT < 10").VALUE, ref_edistr.data.query("EOUT < 10").VALUE * pert)
    >>> np.testing.assert_array_almost_equal(pert_edistr.data.query("EOUT >= 10").VALUE, ref_edistr.data.query("EOUT >= 10").VALUE)


    Test that cross sections are correctly perturbed.
    Example for H1.

    Creation of dummy perturbation for scattering `MT=2`.

    >>> interval = pd.Interval(left=1e-8, right=10, closed="right")
    >>> idx = pd.MultiIndex.from_tuples([(125, 2, interval)], names=("MAT", "MT", "E"))
    >>> pert = 1.2
    >>> df = pd.DataFrame([[pert]], index=idx).reset_index()
    >>> smps = sandy.Samples(df.set_index(["MAT", "MT", "E"]))
    >>> smps
    SMP                            0
    MAT MT E                        
    125 2  (1e-08, 10.0] 1.20000e+00

    Creation of reference ENDF6 and PENDF.

    >>> ref_endf6 = sandy.get_endf6_file("jeff_33", "xs", 10010, local=True)
    >>> ref_pendf = ref_endf6.get_pendf(err=1)

    Creation of perturbed data modifying the PFNS with the first perturbation sample.
    
    >>> ismp = 0
    >>> perturbed = sandy.endf6._endf6_perturb_worker(ref_endf6.data, ref_pendf.data, ismp, pxs=dict(smps.iterate_xs_samples())[ismp])
    >>> pert_endf6 = sandy.Endf6(perturbed['endf6'])
    >>> pert_pendf = sandy.Endf6(perturbed['pendf'])

    Creation of reference and perturbed :obj:`~sandy.xs.Xs`.
    
    >>> ref_xs = sandy.Xs.from_endf6(ref_pendf)
    >>> pert_xs = sandy.Xs.from_endf6(pert_pendf)

    Test that the perturbation is correct.

    >>> np.testing.assert_array_almost_equal(pert_xs.data.query("E<=10")[(125, 2)], ref_xs.data.query("E<=10")[(125, 2)] * 1.2)
    >>> np.testing.assert_array_almost_equal(pert_xs.data.query("E>10")[(125, 2)], ref_xs.data.query("E>10")[(125, 2)])
    >>> np.testing.assert_array_almost_equal(pert_xs.data[(125, 102)], ref_xs.data[(125, 102)])
    >>> assert not np.array_equal(pert_xs.data[(125, 1)], ref_xs.data[(125, 1)])

    """
    from copy import deepcopy
    from .xs import xs_perturb_worker, Xs
    from .edistr import Edistr
    from .zam import za2zam, zam2za

    # --- Initialize data, get them back as Endf6 instances ---
    endf6_pert = Endf6(deepcopy(endf6))
    pendf_pert = Endf6(deepcopy(pendf))

    # apply nubar perturbation
    if pnu is not None:
        nu = Xs.from_endf6(endf6_pert.filter_by(listmt=[452, 455, 456]))
        nu_pert = xs_perturb_worker(nu, ismp, pnu, verbose=verbose)
        endf6_pert = nu_pert.reconstruct_sums(drop=True).to_endf6(endf6_pert).update_intro()

    # apply lpc perturbation
    if plpc is not None:
        pass

    # Apply energy distribution (edistr) perturbation
    if pchi is not None:
        # Applies the same perturbation to all incident particle energies (EIN) and K
        edistr_pert = []
        
        # Group data by EIN and K for processing
        for (ein, k), df in Edistr.from_endf6(endf6_pert).data.groupby(['EIN', 'K']):
            # Prepare dummy energy distribution data as a xs object
            dummy_xs = Xs(
                df.rename({"EOUT": "E"}, axis=1)
                  .set_index(["MAT","MT"])[["E","VALUE"]]
                  .pivot(columns="E").T.droplevel(level=0)
            )

            # Apply perturbation to dummy energy distribution
            dummy_xs_pert = xs_perturb_worker(dummy_xs, ismp, pchi, verbose=verbose)
            
            # Transform xs data into edistr data and append perturbed data
            perturbed_data = (
                dummy_xs_pert.data.stack([1, 0], future_stack=True)  # Use future_stack=True to adopt the new behavior
                .to_frame()
                .reset_index()
                .rename({"E": "EOUT", 0: "VALUE"}, axis=1)
                .assign(K=k, EIN=ein)
                [["MAT", "MT", "K", "EIN", "EOUT", "VALUE"]]
            )
            edistr_pert.append(perturbed_data)

        # Combine and normalize perturbed data, then update ENDF6
        endf6_pert = (
            Edistr(pd.concat(edistr_pert, ignore_index=True))
            .normalize()
            .to_endf6(endf6_pert)
            .update_intro()
            )

    # apply xs perturbation
    if pxs is not None:
        xs = Xs.from_endf6(pendf_pert)
        xs_pert = xs_perturb_worker(xs, ismp, pxs, verbose=verbose)
        pendf_pert = xs_pert.reconstruct_sums(drop=True).to_endf6(pendf_pert).update_intro()


    # --- Return perturbed ENDF6 and PENDF instances as dict
    out_dict = {
        "endf6": endf6_pert.data,
        "pendf": pendf_pert.data,
        }

    if to_ace:
        ace_kws_ = {} if ace_kws is None else ace_kws.copy()

        # --- Add ACE file content and XSDIR file content to the dict
        out_extra = endf6_pert.get_ace(pendf=pendf_pert, **ace_kws_)   # this is a dict {'ace': text, 'xsdir': text}

        # --- I make the dict update explicit
        out_dict.update({
            "ace": out_extra["ace"],
            "xsdir": out_extra["xsdir"],
            })

    
    if to_file:
        # --- The output files basename is hardcoded. It is too complex to maintain a flexible structure
        mat = endf6_pert.mat[0]
        
        intro = endf6_pert.read_section(mat, 1, 451)
        za = int(intro["ZA"])
        meta = int(intro["LISO"])
        zam = za2zam(za, meta=meta, method=False)
        za_nndc = zam2za(zam, method="nndc")[0]
    
        basename = f"{za_nndc}_{ismp}"
        
        # --- Write files to disk and return only the filenames
        out_dict = _write_files_worker(out_dict, basename=basename, verbose=verbose)

    return out_dict


def _write_files_worker(file_dict, basename="output", verbose=False):
    """
    Write ENDF-6, PENDF, and (optionally) ACE and XSDIR contents to disk.

    Parameters
    ----------
    file_dict : Mapping[str, Any]
        Dictionary with optional keys:
          - "endf6": dict-like data to build a `:obj:~sandy.endf6.Endf6` object
          - "pendf": dict-like data to build a `:obj:~sandy.endf6.Endf6` object
          - "ace"  : str content of ACE file (text)
          - "xsdir": str content of XSDIR file (text)

        It is assumed the ENDF-6/PENDF dicts are compatible with `:obj:~sandy.endf6.Endf6`.

    basename : str, default "output"
        Basename (without extension) used for ENDF-6/PENDF files,
        and also as prefix for ACE (*.XXc).

    verbose : bool, default False
        If True, emits progress messages via `logger`.

    Returns
    -------
    Dict[str, str]
        Mapping with the filenames written. Keys may include:
        - "endf6": path to `<basename>.endf6`
        - "pendf": path to `<basename>.pendf`
        - "ace"  : path to `<basename>.<SUFFIX>c` (if "ace" provided)
        - "xsdir": path to `<basename>.<SUFFIX>c` (if "xsdir" provided)

    Raises
    ------
    KeyError
        If required keys "endf6" or "pendf" are missing in `file_dict`.

    Notes
    -----
    - We inferring the ACE/XSDIR suffix from the ACE/XSDIR content.
    - ENDF6 and PENDF are dict and not `:obj:~sandy.endf6.Endf6` to be pickled.

    Examples
    --------
    
    Collect endf6/pendf/ace/xsdir file for test.
    
    >>> import sandy
    >>> endf6 = sandy.get_endf6_file("jeff_33", "xs", 10010, local=True)
    >>> pendf = endf6.get_pendf()
    >>> out_dict = {"endf6": endf6.data, "pendf": pendf.data}
    >>> out_dict |= endf6.get_ace(temperature=0)
    
    Run worker...but first remove outputs.

    >>> from pathlib import Path
    >>> for filename in ['output.00c', 'output.00c.xsd', 'output.endf6', 'output.pendf']:
    ...    p = Path(filename)
    ...    if p.exists():
    ...       p.unlink()

    >>> outfiles = sandy.endf6._write_files_worker(out_dict, verbose=True)
    
    Check that outputs have been created.

    >>> from os.path import isfile
    >>> assert(isfile('output.00c'))
    >>> assert(isfile('output.00c.xsd'))
    >>> assert(isfile('output.endf6'))
    >>> assert(isfile('output.pendf'))

    """
    from .utils import log

    # --- Write files to disk and return only the filenames
    outfiles= {}
    
    def extract_suffix(text):
        """Extract the suffix (example, "03" from 1001.03c) from the ace/xsdir file.
        No fallback if it doesn't work"""
        suffix = text.split()[0].split(".")[1][:2]
        return suffix
    
    if 'ace' in file_dict:
        suffix = extract_suffix(file_dict["ace"])
        file = f"{basename}.{suffix}c"
        log(f" - Writing ACE to file '{file}'", verbose=verbose)

        with open(file, "w") as f:
            f.write(file_dict["ace"])

        outfiles["ace"] = file

    if 'xsdir' in file_dict:
        suffix = extract_suffix(file_dict["xsdir"])
        file = f"{basename}.{suffix}c.xsd"
        log(f" - Writing XSD to file '{file}'", verbose=verbose)

        with open(file, "w") as f:
            f.write(file_dict["xsdir"])

        outfiles["xsdir"] = file
    
    if 'endf6' in file_dict:
        endf6 = Endf6(file_dict["endf6"])
        file = f"{basename}.endf6"
        log(f" - Writing ENDF-6 to file '{file}'", verbose=verbose)
        endf6.to_file(file)
        outfiles["endf6"] = file

    if 'pendf' in file_dict:
        pendf = Endf6(file_dict["pendf"])
        file = f"{basename}.pendf"
        log(f" - Writing PENDF to file '{file}'", verbose=verbose)
        pendf.to_file(file)
        outfiles["pendf"] = file
    
    return outfiles
    

def _rdd_perturb_worker(endf6, rdd, smp_hl, smp_de, smp_br, ismp,
                       verbose=False, to_file=False, **kwargs):
    """
    Worker to handle ENDF6 radioactive decay data perturbation.

    Parameters
    ----------
    endf6 : `dict`
        `data` attribute of :obj:`~sandy.endf6.Endf6`.
        It contains the nominal ENDF6 data.
    rdd : `pd.DataFrame`
        `data` attribute of :obj:`~sandy.decay.DecayData`.
        It contains the nominal decay data.
    smp_hl : `pd.DataFrame`
        `data` attribute of :obj:`~sandy.samples.Samples`.
        It contains the perturbation coefficients for half-lives.
    smp_de : `pd.DataFrame`
        `data` attribute of :obj:`~sandy.samples.Samples`.
        It contains the perturbation coefficients for decay energies.
    smp_br : `pd.DataFrame`
        `data` attribute of :obj:`~sandy.samples.Samples`.
        It contains the perturbation coefficients for branching ratios.
    ismp : `int`
        sample ID.
    verbose : `bool`, optional
        Flag to activate verbosity. The default is False.
    to_file : `bool`, optional
        Flag to write outputs to file. The default is False.
        This key changes the output type.
    **kwargs : `dict`
        Additional keyword arguments (not used).

    Returns
    -------
    `dict`
        Either a dictionary of :obj:`~sandy.endf6.Endf6` instances for each set of
        perturbation coefficients (if `to_file=False`), or a dictionary
        of `str` with the output file name for each set of perturbation
        coefficients.

    Notes
    -----
    .. note: This method is written so that it can be handled by the
             `multiprocess` module (pickling).

    .. note: Branching ratios are renormalized.
    """
    from .decay import DecayData
    from .samples import Samples
    from .utils import log

    endf6_ = Endf6(endf6.copy())
    rdd_ = DecayData(rdd.copy())
    
    smp_hl_ = Samples(smp_hl.copy())
    smp_de_ = Samples(smp_de.copy())
    smp_br_ = Samples(smp_br.copy())
    
    hl_ = rdd_.get_half_life()
    hl_.data["HL"] *= smp_hl_.data[ismp]
    rdd_ = hl_.to_decaydata(rdd_)
    
    de_ = rdd_.get_decay_energy()
    de_.data["E"] *= smp_de_.data[ismp]
    rdd_ = de_.to_decaydata(rdd_)
    
    br_ = rdd_.get_branching_ratio()
    br_.data["BR"] *= smp_br_.data[ismp]
    rdd_ = br_.normalize().to_decaydata(rdd_)
    
    out = rdd_.to_endf6(endf6_)
    
    # Stop here and return dict of Endf6 instance. not Endf6 because it cannot be pickled
    if not to_file:
        return out.data
 
    # continue and return filename where data was written
    file = f"decay_data_{ismp}"
    log(f" - Writing RDD to file '{file}'", verbose=verbose)
    out.to_file(file)

    return file



def _fy_perturb_worker(endf6, fy, smps, ismp,
                       verbose=False, to_file=False, **kwargs):
    """
    Worker to handle ENDF6 fission yield perturbation.

    Parameters
    ----------
    endf6 : `dict`
        `data` attribute of :obj:`~sandy.endf6.Endf6`.
        It contains the nominal ENDF6 data.
    fy : `pd.DataFrame`
        `data` attribute of :obj:`~sandy.fy.Fy`.
        It contains the nominal fission yield data.
        It i sassume they match all the ZAP of the samples.
    smps : `pd.DataFrame`
        It contains the perturbation coefficients for fission yields.
        Columns are `MAT`, `MT`, `E`, `ZAM`, `ZAP`, `SMP`, `VALS`.
        This dataframe is generally produced with `pd.pivot_table`.
    ismp : `int`
        sample ID.
    verbose : `bool`, optional
        Flag to activate verbosity. The default is False.
    to_file : `bool`, optional
        Flag to write outputs to file. The default is False.
        This key changes the output type.
    **kwargs : `dict`
        Additional keyword arguments (not used).
        
    Notes
    -----
    .. note:: It follows the logic of :obj:`~sandy.endf6._endf6_perturb_worker` and
              :obj:`~sandy.endf6._rdd_perturb_worker`.

    Returns
    -------
    `dict`
        Either a dictionary of :obj:`~sandy.endf6.Endf6` instances for each set of
        perturbation coefficients (if `to_file=False`), or a dictionary
        of `str` with the output file name for each set of perturbation
        coefficients.

    Notes
    -----
    .. note: This method is written so that it can be handled by the
             `multiprocess` module (pickling).

    Examples
    --------
    
    Default test: create 1 sample and perturb fission yields for 1 fissioning system.
    
    >>> import sandy
    >>> nsmp = 1   # sample size
    >>> zam, e = 922350, 0.0253
    >>> tape = sandy.get_endf6_file("jeff_33", "nfpy", zam, local=True)
    >>> nfpy = sandy.Fy.from_endf6(tape)
    >>> idx = nfpy.data.query(f"E=={e} & MT==454 & ZAM=={zam}").index
    >>> fy = nfpy.data.loc[idx]
    >>> smps = sandy.CategoryCov(pd.DataFrame(np.diag((fy.DFY/fy.FY)**2), index=fy.ZAP, columns=fy.ZAP).fillna(0)).sampling(nsmp)
    >>> smps = smps.data.rename_axis(index="ZAP").stack().rename("VALS").reset_index().assign(E=e, ZAM=zam)[["ZAM", "E", "ZAP", "SMP", "VALS"]]
    >>> out = sandy.endf6._fy_perturb_worker(tape.data, nfpy.data, smps, nsmp-1, verbose=True, to_file=False)
    >>> out = sandy.Endf6(out)
    
    Silly test: assert the `MT=454` was changed, and `MT=459` was not.

    >>> assert sandy.Fy.from_endf6(out).data.query("MT==459").equals(nfpy.data.query("MT==459"))
    >>> assert not sandy.Fy.from_endf6(out).data.query("MT==454").equals(nfpy.data.query("MT==454"))
    
    Test to check that the output random ENDF6 are perturbed correctly.

    >>> tape = sandy.get_endf6_file("jeff_33", "nfpy", 922350, local=True)
    >>> smps = tape.get_perturbations(2, covariance=None)
    >>> nfpy = sandy.Fy.from_endf6(tape)
    >>> out = sandy.endf6._fy_perturb_worker(tape.data, nfpy.data, smps, 0)
    >>> nfpy0 = sandy.Fy.from_endf6(sandy.Endf6(out))

    Assert that ratio of perturbed to nominal FY's is equal to samples.

    >>> n = nfpy.data.query("ZAM==922350 and MT==454")
    >>> n0 = nfpy0.data.query("ZAM==922350 and MT==454")
    >>> assert not n.equals(n0)
    >>> sp = (n0.set_index(["MAT", "MT", "ZAM", "E", "ZAP"]).FY /  n.set_index(["MAT", "MT", "ZAM", "E", "ZAP"]).FY).fillna(1)
    >>> p = smps.query("ZAM==922350 and SMP==0").VALS
    >>> np.testing.assert_array_almost_equal(p, sp, decimal=4)
    """
    from .fy import Fy  # lazy import to avoid circular import issue
    from .utils import log

    endf6_ = Endf6(endf6.copy())  # this was a dictionary
    fy_ = Fy(fy.copy())    # this was a dataframe

    for (zam, e), smp in smps.groupby(["ZAM", "E"]):
        idx = fy_.data.query(f"ZAM=={zam} & E=={e} & MT==454").index

        # do not assume both FY's and perturbations are sorted, make them match by ZAP
        zap = fy_.data.loc[idx]["ZAP"]
        # update data directly in Fy instance
        fy_.data.loc[idx, "FY"] *= smp.query(f"SMP=={ismp}").set_index("ZAP").loc[zap].VALS.values
        
        # IMPORTANT, this does not update the CFYs, which in random ENDF-6 file are inconsistent with the perturbed IFYs

    out = fy_.to_endf6(endf6_)
    
    # Stop here and return dict of Endf6 instance. not Endf6 because it cannot be pickled
    if not to_file:
        return out.data
 
    # continue and return filename where data was written
    file = f"fy_{ismp}"
    log(f" - Writing NFY to file '{file}'")
    out.to_file(file)

    return file
    
