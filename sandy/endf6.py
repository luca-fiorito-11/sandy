import io
import os
import logging
import time
from typing import Iterable

from .utils import with_optional_warning_suppression

__author__ = "Luca Fiorito"


nsubs = {
    4: "decay",
    10: "neutron",
    11: "nfpy",
    10010: "proton",
}


def get_endf6_file(
        library: str,
        kind: str,
        zam: int | Iterable[int] | str,
        to_file: bool = False,
        local: bool = False,
        verbose: bool =False,
        ) -> "Endf6":


    """
    Retrieve an ENDF‑6 evaluation for a given nuclear data library, data type,
    and nuclide (ZAM). Data can be downloaded from online sources or loaded
    from SANDY's local library appendix.

    This function provides a unified interface for retrieving:

    - neutron-induced cross sections (``kind="xs"``)
    - fission product yields (``kind="nfpy"``)
    - radioactive decay data (``kind="decay"``)
    - thermal scattering laws (``kind="tsl"``)
    - displacement cross sections (``kind="dxs"``)

    A central registry maps each (kind, library) pair to the appropriate
    remote URL and corresponding filename mapping.

    Parameters
    ----------
    library : str
        Name of the nuclear data library. Available libraries depend on
        ``kind`` and typically include:
        ``"endfb_71"``, ``"endfb_80"``, ``"endfb_81"``, ``"tendl_2023"``,
        ``"jeff_311"``, ``"jeff_33"``, ``"jeff_40"``, ``"jendl_40u"``,
        ``"jendl_5"``, ``"irdff_2"``, etc.

    kind : str
        Type of nuclear data to retrieve. Must be one of:
            - ``"xs"``: neutron-induced cross sections
            - ``"nfpy"``: fission product yields
            - ``"decay"``: radioactive decay data
            - ``"tsl"``: thermal scattering laws
            - ``"dxs"``: displacement cross sections

    zam : int, iterable of int, or str
        ZAM identifier(s) in the format ``Z*10000 + A*10 + M``,
        where ``M`` is the metastable state index (0 for ground state).

        - If an integer: load a single nuclide.
        - If an iterable: merge multiple nuclides into one ENDF‑6 tape.
        - If ``"all"``: retrieve the entire decay or NFpy library (only for
          ``kind="decay"`` and ``kind="nfpy"``).

    to_file : bool, optional
        If True, save the resulting ENDF‑6 file to disk.  
        Default is False.

    local : bool, optional
        If True, load files from SANDY's local library cache instead of downloading.
        Default is False.

    verbose : bool, optional
        If True, print progress messages.  
        Default is False.

    Returns
    -------
    Endf6
        A populated :class:`~sandy.endf6.Endf6` object containing the
        requested ENDF‑6 evaluation.

    Raises
    ------
    ValueError
        If the requested ``kind`` is unsupported or if the library is not
        available for that ``kind``.
    ValueError
        If ``zam="all"`` is used with a ``kind`` that does not support it.

    Notes
    -----
    - ZAM format: ``Z*10000 + A*10 + M``  
      Example: Hydrogen-1 → ``10010``.
    - For ``kind="tsl"``, ``zam`` must contain the integer indices used in
      the thermal scattering law library.
    - When possible, remote access uses ZIP archives; otherwise it falls back
      to raw URLs. Local archives for ``"all"`` selections use compressed
      ``.tar.xz`` files.

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

    Import all Decay Data for the supported libraries.
    They are stored locally and use :meth:`~sandy.endf6._FormattedFile.from_xzfile`.

    >>> for lib in ["jeff_311", "jeff_33", "jeff_40", "endfb_71", "endfb_80", "endfb_81", "jendl_5"]:
    ...    tape = sandy.get_endf6_file(lib, 'decay', 'all')
    ...    assert type(tape) is sandy.Endf6

    Import all FIssion Product Yield Data for the supported libraries.
    They are stored locally and use :meth:`~sandy.endf6._FormattedFile.from_xzfile`.

    >>> tape = sandy.get_endf6_file("jeff_33", 'nfpy', 'all')
    >>> for lib in ["jeff_311", "jeff_33", "jeff_40", "endfb_71", "endfb_80", "endfb_81", "jendl_40u"]:
    ...    tape = sandy.get_endf6_file(lib, 'nfpy', 'all')
    ...    assert type(tape) is sandy.Endf6

    Thermal Neutron Scattering Data from JEFF-3.3.

    >>> tape = sandy.get_endf6_file("jeff_33", 'tsl', [1, 2, 3], local=True)
    >>> assert type(tape) is sandy.Endf6

    """
    # ---- IMPORT
    from functools import reduce
    from pathlib import Path
    from os.path import splitext

    from . import __file__ as sandy__file__
    from .zam import zam2nuclide
    from .utils import log
    from ._perturbation_base import log_stage
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
    # ---- SETUP
    method = "get_endf6_file"


    # ---- REGISTRY OF SUPPORTED LIBRARIES 
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



    # ---- VALIDATE INPUT
    kind_ = kind.lower()
    if kind_ not in maps:
        raise ValueError(f"Unsupported kind='{kind_}'. Valid options: {list(maps.keys())}")

    library_ = library.lower()
    if library_ not in maps[kind_]:
        supported = ", ".join(maps[kind_].keys())
        raise ValueError(
            f"Library '{library}' not available for kind='{kind_}'. "
            f"Supported libraries: {supported}"
        )
    url, files = maps[kind_][library_]


    # ---- LOCAL PATHS
    base = Path(sandy__file__).resolve().parent
    local_space = base / "appendix" / "libraries"
    local_space_onefile = base / "appendix" / "onefile_archives"
    local_path = local_space / library_ / kind_


    # ---- HELPERS

    def load_remote_zip(file: str, url: str) -> Endf6:
        return Endf6.from_zipurl(file, url)

    def load_local_zip(file: str, url: str) -> Endf6:
        # mimics load_remote_zip, even if url is not used
        file = splitext(file)[0]
        archive = local_path / f"{file}.zip"
        tape = Endf6.from_zipfile(archive)
        return tape

    def load_local_xz(file: str, url: str) -> Endf6:
        # 'all' decay data and fission yields are stored as xz
        # mimics load_remote_zip, even if url is not used
        file = splitext(file)[0]
        archive = local_space_onefile / f"{file}.tar.xz"
        tape = Endf6.from_xzfile(archive)
        return tape


    # ---- PICK FETCHER
    fetcher = load_local_zip if local else load_remote_zip


    # ---- HANDLE zam == "all"
    if str(zam).lower() == 'all':
        if kind_ not in ('decay', "nfpy"):
            raise ValueError("'zam=\"all\"' is only supported for 'decay' and 'nfpy'.")

        # --- fall back on local files with all data, otherwise it's too slow
        fetcher = load_local_xz  # only stored locally
        filename = f"{kind}_{library_}.dat"
        tape = fetcher(filename, url="not used")

    # ---- STANDARD FETCH WITH FALLBACK
    else:

        # --- helper: try remote first, then fall back to local
        def fetch_with_fallback(key: int) -> Endf6:
            fname = files[key]
            try:
                return fetcher(fname, url)    # remote first
            except Exception as err:
                msg = (
                    f"Remote download failed for ZAM={key}: {err}\n"
                    f"→ Falling back to local file in:\n"
                    fr"   {local_path}"
                )
                warn_logger = logging.getLogger("sandy.warn")
                log(msg, level=logging.WARNING, logger=warn_logger)
                return load_local_zip(files[key], url=None)

        # ---- MULTIPLE ZAMs
        if isinstance(zam, Iterable) and not isinstance(zam, (str, bytes)):
            # --- CASE 1: a list of nuclides is passed (all valid ZAM)
            tapes = [fetch_with_fallback(x) for x in zam]
            combined = reduce(lambda x, y: x.add_sections(y.data), tapes)
            tape = combined

        # ---- SINGLE ZAM
        else:
            tape = fetch_with_fallback(zam)

    # ---- OPTIONAL WRITE TO FILE
    if to_file:
        basename = zam2nuclide(zam, atomic_number=True, sep="-")
        filename = f"{basename}.{library_}"
        msg = f"saving ENDF6 file to '{filename}'"
        log_stage(log, method, None, msg, verbose=verbose)
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
    def mat(self) -> list[int]:
        return sorted({int(mat) for mat in self._keys["MAT"]})


    @property
    def mf(self) -> list[int]:
        return sorted({int(mf) for mf in self._keys["MF"]})

    @property
    def mt(self) -> list[int]:
        return sorted({int(mt) for mt in self._keys["MT"]})

    def to_series(self, **kwargs):
        import pandas as pd
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
        from os.path import splitext, join
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
            # <- triggers exception for 404, 403, etc.
            r.raise_for_status()
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
        >>> obj = sandy.endf6._FormattedFile.from_file(file)

        The returned object must be a formatted ENDF-6 file:

        >>> assert isinstance(obj, sandy.endf6._FormattedFile)

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
        >>> obj2 = sandy.endf6._FormattedFile.from_file(stream)

        >>> assert obj2.data == obj.data

        """
        if isinstance(file, io.StringIO):
            text = file.read()

        else:
            with open(file) as f:
                text = f.read()
        return cls.from_text(text)

    @classmethod
    def from_xzfile(cls, xz_filename):
        """
        Read and return the contents of the only file inside an XZ‑compressed
        tar archive (.tar.xz).
    
        This method assumes the archive contains exactly one file, typically
        with the same base name as the archive itself. The internal file is
        read directly from the compressed tar archive without extracting it
        to disk.
    
        Use case
        --------
        ENDF-6 data distributed as `.tar.xz` archives, such as those containing
        a single `.dat` file.
    
        Parameters
        ----------
        xz_filename : str
            Path to the `.tar.xz` file on disk.
    
        Returns
        -------
        :obj:`~sandy.endf6.Endf6`
            An instance created from the decoded text content of the internal file.
    
        Raises
        ------
        FileNotFoundError
            If the archive file does not exist.
        tarfile.ReadError
            If the file is not a valid tar or tar.xz archive.
        IndexError
            If the archive contains no files.
        UnicodeDecodeError
            If the internal file cannot be decoded as UTF‑8 text.
    
        Notes
        -----
        - No temporary files or directories are created.
        - The internal file is read entirely into memory.
        
        Examples
        --------
        This way of fetching data is used in :func:`~sandy.endf6.get_endf6_file`.
        Here it is testd for the fission yield data of JEFF-3.3.

        >>> import sandy
        >>> from pathlib import Path
        >>> base = Path(sandy.__file__).resolve().parent
        >>> local_space_onefile = base / "appendix" / "onefile_archives"
        >>> file = "nfpy_jeff_33.tar.xz"
        >>> archive = local_space_onefile / file
        >>> tape = Endf6.from_xzfile(archive)        
        """
        from .utils import read_xzfile

        text = read_xzfile(xz_filename, member=0).getvalue().decode("utf-8")

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
        import pandas as pd

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

        # NumPy ops → no pandas.core.ops
        mask = (mt > 0) & (mf > 0) & (mat > 0)

        df2 = df.loc[mask, ["MAT", "MF", "MT", "TEXT"]]

        # -----------------------------
        # 5. Manual group-by (avoids pandas.groupby → no NumExpr paths)
        # -----------------------------
        data = {}
        # group rows by (MAT, MF, MT)
        for (mat_v, mf_v, mt_v), group in df2.groupby(["MAT", "MF", "MT"], sort=False):
            data[(mat_v, mf_v, mt_v)] = "\n".join(group["TEXT"].tolist())

        return cls(data)

    def _get_section_records(self, mat, mf, mt):
        """
        Very fast ENDF fixed-width parser.
        Returns list of string tuples: (C1, C2, L1, L2, N1, N2)
        exactly like a DataFrame row in the slow version.
        """
        from .records import get_records_from_text

        text = self.data[(mat, mf, mt)]

        return get_records_from_text(text)

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
        d = df.query(
            "MAT in @listmat and MF in @listmf and MT in @listmt").squeeze(axis=1).to_dict()
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
            new_records = [(mf, mt, sec.count('\n') + 1, 0)
                           for (mat, mf, mt), sec in g.items()]
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
                raise ValueError(f"Module '{modname}' does not define '{func_name}'")
            return None

        return reader(self, mat, mt)

    def _derive_basename_for_sampling(self, ismp: int) -> str:
        """
        Compute a stable basename like '<ZA(NNDC)>_<ismp>'.
        """
        from .zam import za2zam, zam2za

        mat = self.mat[0]
        intro = self.read_section(mat, 1, 451)

        za = int(intro["ZA"])
        meta = int(intro["LISO"])
        zam = za2zam(za, meta=meta, method=False)
        za_nndc = zam2za(zam, method="nndc")[0]

        return f"{za_nndc}_{ismp}"

    def _update_info(self, descr=None):
        """
        Update RECORDS item (in DATA column) for MF1/MT451 of each MAT based on the content of the TEXT column.
        """
        # ---- IMPORT
        import pandas as pd

        from .mf1 import write

        tape = self.copy()
        for mat in sorted(tape.index.get_level_values('MAT').unique()):
            sec = self.read_section(mat, 1, 451)
            records = pd.DataFrame(sec["RECORDS"], columns=[
                                   "MF", "MT", "NC", "MOD"]).set_index(["MF", "MT"])
            new_records = []
            dfmat = tape.loc[mat]
#            for (mf,mt),text in sorted(tape.loc[mat].query('MT!=451'.format(mat)).TEXT.items()):
            for (mf, mt), text in sorted(dfmat[dfmat.index.get_level_values("MT") != 451].TEXT.items()):
                nc = len(text.splitlines())
                # when copying PENDF sections (MF2/MT152) mod is not present in the dictionary
                try:
                    mod = records.MOD.loc[mf, mt]
                except:
                    mod = 0
                new_records.append((mf, mt, nc, mod))
            if descr is not None:
                sec["TEXT"] = descr
            nc = 4 + len(sec["TEXT"]) + len(new_records) + 1
            mod = records.MOD.loc[1, 451]
            new_records = [(1, 451, nc, mod)] + new_records
            sec["RECORDS"] = new_records
            text = write(sec)
            tape.loc[mat, 1, 451].TEXT = text
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

    def get_mat_zam_mapping(
            self,
            ) -> dict[int, int]:
        """
        Return a dictionary mapping ENDF material numbers (MAT) to their
        corresponding ZAM identifiers.
    
        The ZAM identifier is defined as::
    
            ZAM = ZA * 10 + LISO
    
        where:
            - ``ZA``   = Z*1000 + A (standard ENDF nuclide identifier)
            - ``LISO`` = metastable state index (0 = ground state)
    
        Returns
        -------
        dict[int, int]
            Dictionary where:
            - keys   = MAT numbers in the ENDF6 file
            - values = ZAM identifiers (ZA * 10 + LISO)
    
        Notes
        -----
        ZAM is an ENDF convention combining ZA and metastable state into
        a single integer. Examples:
    
        - Z = 92, A = 235, ground state:
          ZA = 92235 → ZAM = 922350
        - Z = 95, A = 242, metastable 1:
          ZA = 95242 → ZAM = 952421
        
        Examples
        --------
        Standard use.

        >>> import sandy
        >>> tape = sandy.get_endf6_file("jeff_33", "decay", [270590, 270600], local=True)

        This should return a ``dict``.

        >>> assert tape.get_mat_zam_mapping() == {561: 270590, 562: 270600}

        This should return a ``dict``.

        >>> tape = sandy.get_endf6_file("jeff_33", "decay", 270590, local=True)
        >>> assert tape.get_mat_zam_mapping() == {561: 270590}
                
        Test also a metastable nuclide.

        >>> tape = sandy.get_endf6_file("jeff_33", "decay", 591481, local=True)
        >>> assert tape.get_mat_zam_mapping() == {2007: 591481}
        """
        mat_zam_mapping: dict[int, int] = {}
    
        for mat in self.mat:
            info = self.read_section(mat, 1, 451)
            meta = int(info["LISO"])
            za = int(info["ZA"])
            zam = za * 10 + meta
            mat_zam_mapping[mat] = zam

        return mat_zam_mapping

    def get_library(self) -> str:
        mat = self.mat
        if len(mat) != 1:
            raise Exception(
                "'get_library' only works with single MAT. "
                f"Found {len(mat)} of them"
                )

        mat = mat[0]
        if (mat, 1, 451) not in self.data:
            raise Exception(
                f"Section MAT={mat}, MF=1, MT=451 not found"
                )

        descr = self.read_section(mat, 1, 451)["DESCRIPTION"]
        lib = descr[2][4:22].strip()
        return lib

    def get_zam(self) -> int | list[int]:
        """
        Return the ZAM identifiers for all materials in the ENDF6 file.
    
        Each ZAM value is computed as ``ZA * 10 + LISO``, where:
    
        - ``ZA`` is the unique nuclide identifier (Z*1000 + A)
        - ``LISO`` is the metastable state index
    
        Returns
        -------
        int or list[int]
            - If only one material is present, returns a single integer.
            - If multiple materials are present, returns a list of integers.
    
        Notes
        -----
        ZAM is a common ENDF convention for identifying nuclides including
        metastable states. For example:
        - Z=92, A=235, ground state → ZA = 92235 → ZAM = 922350
        - Z=95, A=242, metastable 1 → ZA = 95242 → ZAM = 952421
        
        Examples
        --------
        Standard use.

        >>> import sandy
        >>> tape = sandy.get_endf6_file("jeff_33", "decay", [270590, 270600], local=True)

        It should return a list.

        >>> assert tape.get_zam() == [270590, 270600]

        This should return a scalar.

        >>> tape = sandy.get_endf6_file("jeff_33", "decay", 270590, local=True)
        >>> assert tape.get_zam() == 270590
                
        Test also a metastable nuclide.

        >>> tape = sandy.get_endf6_file("jeff_33", "decay", 591481, local=True)
        >>> assert tape.get_zam() == 591481
        """
        mat_zam_mapping = self.get_mat_zam_mapping()
        zam_list = list(mat_zam_mapping.values())
    
        if len(zam_list) == 1:
            return zam_list[0]

        return zam_list

    def _run_njoy(
            self,
            pendf=None,
            pendftape=None,
            print_njoy_input=False,
            verbose=False,
            **njoy_kws):
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
        - Keyword argument `pendf` is used to pass aPENDF as `Endf6` object,
          while `pendftape` is used to pass a PENDF as the name of a file written on disk.
        - Better to keep the logging minimal here. Only timing the njoy run.
        """
        # ---- IMPORT
        from tempfile import TemporaryDirectory
        from os.path import join

        from .njoy import process_neutron
        from .utils import log
        
        # ---- SETUP
        zam = self.get_zam()
        
        common_msg = f"_run_njoy | ZAM={zam} "

        with TemporaryDirectory() as td:
            # ---- WRITE ENDF6 main tape to temp folder
            endf6file = join(td, "tape20")
            self.to_file(endf6file)

            # ---- HANDLE optional PENDF input
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

            # ---- RUN NJOY through sandy
            t0 = time.perf_counter()

            outputs = process_neutron(
                endf6file,
                pendftape=pendf_file,
                verbose=print_njoy_input,
                **njoy_kws,
            )

            td = time.perf_counter() - t0

        # ---- LOGGING: end
        dt = time.perf_counter() - t0
        msg = (
            f"| NJOY ran in {dt:.3f} s"
        )
        log(common_msg + msg, verbose=verbose)

        # ---- RETURN NJOY output (dict)
        return outputs

    def _prepare_groupr_kws(self, **groupr_kws):
        """Helper to prepare groupr options"""
        # -- prepare/augment GROUPR options without mutating the user's dict --
        groupr_kws_ = dict(groupr_kws)

        # Decide what GROUPR should process based on available sections
        recs = self.get_records()
        has_fission_xs = 18 in recs.query(
            "MF==3").MT.values   # XS fiss present?
        # fission XS implies nubar processing
        groupr_kws_.setdefault("nubar", has_fission_xs)
        has_fission_xs = 18 in recs.query("MF==5").MT.values   # PFNS present?
        # PFNS implies chi processing
        groupr_kws_.setdefault("chi",   has_fission_xs)
        # always include mubar
        groupr_kws_.setdefault("mubar", True)

        return groupr_kws_

    def _prepare_njoy_kws(
            self,
            temperature=0,
            err=0.001,
            minimal_processing=False,
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
        - Better to keep the logging minimal here. Only timing the njoy run.

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
        # ---- IMPORT
        from .utils import log

        njoy_kws_ = dict(njoy_kws)  # clean copy

        # minimal processing
        if minimal_processing or float(temperature) == 0:
            njoy_kws_.update(dict(
                thermr=False, gaspr=False, heatr=False, purr=False, unresr=False
            ))

        # stop after reconr
        if temperature == 0:
            njoy_kws_["broadr"] = False
            msg = (
                "Zero Kelvin requested; NJOY will stop after RECONR. "
                "Use temperature=0.1 for 0K xs processing."
            )
            warn_logger = logging.getLogger("sandy.warn")
            log(msg, level=logging.WARNING, logger=warn_logger)

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
        njoy_kws_ = self._prepare_njoy_kws(
            **njoy_kws) | {"dryrun": dryrun, "pendf": pendf}
        njoy_kws_ |= ({"suffixes": [suffix]} if suffix is not None else {})

        # --- run via the shared helper ---
        outputs = self._run_njoy(**njoy_kws_)

        if dryrun:
            return outputs

        return {"ace": outputs["ace"], "xsdir": outputs["xsdir"]}

    @with_optional_warning_suppression("sandy.warn", default_suppress=False)
    def get_pendf(
            self,
            *,
            dryrun: bool = False,
            print_njoy_input: bool = False,
            suppress_njoy_output: bool = False,
            verbose: bool | int = False,
            **njoy_kws,
            ):
        """
        Generate a PENDF (pointwise ENDF) file from this :class:`~sandy.endf6.Endf6`
        instance using NJOY.
    
        This method wraps :func:`sandy.njoy.process_neutron` via the internal helper
        :meth:`_run_njoy`, and optionally returns the generated NJOY input deck
        (``dryrun=True``).
    
        Parameters
        ----------
        dryrun : bool, optional
            If ``True``, return the NJOY input text deck **without running NJOY**.
            Default is ``False``.
        verbose : bool or int, optional
            Controls the verbosity level passed to NJOY.  
            - ``False``/``0`` → fully silent
            - ``True``/``1`` → status logging  
            Default is ``False``.
        **njoy_kws : dict
            Additional keyword arguments forwarded to
            :func:`sandy.njoy.process_neutron`.  

            Common options include:
    
            - ``temperature`` : float  
            - ``err`` : int  
            - ``minimal_processing`` : bool  
            - ``suffix`` : str  
            - etc.
    
        Returns
        -------
        :class:`sandy.endf6.Endf6`
            If ``dryrun=False``, the parsed :class:`~sandy.endf6.Endf6` PENDF object.
        dict
            If ``dryrun=True``, the NJOY input deck as a ``str``.
    
        Notes
        -----
        This method:
        - Enforces ``acer=False`` (PENDF generation does not require ACER)
        - Use of the shared NJOY runner ``_run_njoy``

        Examples
        --------
        Default run.

        >>> import sandy
        >>> endf6 = sandy.get_endf6_file("jeff_33", "xs", 10010, local=True)
        >>> kws = dict(suppress_njoy_output=True, suppress_warnings=True, temperature=293.6, err=1, minimal_processing=True)
        >>> out = endf6.get_pendf(**kws)
        >>> assert isinstance(out, Endf6)
        """
        # ---- IMPORT
        from subprocess import DEVNULL
        import pprint

        from .utils import log
        from ._perturbation_base import log_stage

        # ---- SETUP
        zam = self.get_zam()
        
        method = "get_pendf"

        # ---- PREPARE KEYWORDS
        # no mutation, _prepare_njoy_kws returns a copy
        njoy_kws_ = self._prepare_njoy_kws(**njoy_kws)
        njoy_kws_["dryrun"] = dryrun
        njoy_kws_["acer"] = False

        # ---- SUPPRESSING NJOY output (optional)
        if suppress_njoy_output:
            msg = "NJOY output to screen is suppressed"
            log_stage(log, method, zam, msg, verbose=verbose)

            njoy_kws_ |= {
                "njoy_output": DEVNULL
                }

        pdict = pprint.pformat(njoy_kws_, indent=2, sort_dicts=True)
        msg = f"augmented NJOY kwargs: {pdict}"
        log_stage(log, method, zam, msg, verbose=verbose)

        # ---- RUN NJOY via the shared helper
        msg = "run NJOY"
        log_stage(log, method, zam, msg, verbose=verbose)

        # --- run via the shared helper ---
        outputs = self._run_njoy(
            print_njoy_input=print_njoy_input,
            verbose=verbose,
            **njoy_kws_,
            )

        # --- In case of dryrun, 'outputs' contains the text of the NJOY input
        if dryrun:
            msg = "dryrun requested — returning NJOY input deck"
            log_stage(log, method, zam, msg, verbose=verbose)
            return outputs

        msg = "parsing NJOY PENDF output into Endf6 structure"
        log_stage(log, method, zam, msg, verbose=verbose)

        return Endf6.from_text(outputs["pendf"])

    @with_optional_warning_suppression("sandy.warn", default_suppress=False)
    def get_gendf(
            self,
            *,
            dryrun: bool = False,
            groupr_kws=None,
            print_njoy_input: bool = False,
            suppress_njoy_output: bool = False,
            verbose: bool | int = False,
            **njoy_kws,
            ):
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

        >>> import sandy, numpy as np
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
        # ---- IMPORT
        from subprocess import DEVNULL
        import pprint

        from .gendf import Gendf
        from .utils import log
        from ._perturbation_base import log_stage

        # ---- SETUP
        zam = self.get_zam()
        
        method = "get_gendf"

        # ---- PREPARE KEYWORDS
        # no mutation, _prepare_njoy_kws returns a copy
        njoy_kws_ = self._prepare_njoy_kws(**njoy_kws)
        njoy_kws_["dryrun"] = dryrun
        njoy_kws_["groupr"] = True
        njoy_kws_["acer"] = False

        # ---- SUPPRESSING NJOY output (optional)
        if suppress_njoy_output:
            msg = "NJOY output to screen is suppressed"
            log_stage(log, method, zam, msg, verbose=verbose)

            njoy_kws_ |= {
                "njoy_output": DEVNULL
                }

        # -- prepare/augment GROUPR options without mutating the user's dict --
        groupr_kws_ = (groupr_kws or {}).copy()
        njoy_kws_["groupr_kws"] = self._prepare_groupr_kws(**groupr_kws_)

        pdict = pprint.pformat(njoy_kws_, indent=2, sort_dicts=True)
        msg = f"augmented NJOY kwargs: {pdict}"
        log_stage(log, method, zam, msg, verbose=verbose)

        # ---- RUN NJOY via the shared helper
        msg = "run NJOY"
        log_stage(log, method, zam, msg, verbose=verbose)

        # Pass dryrun policy down to NJOY
        njoy_kws_["dryrun"] = dryrun

        # --- run via the shared helper ---
        outputs = self._run_njoy(
            print_njoy_input=print_njoy_input,
            verbose=verbose,
            **njoy_kws_,
            )

        # --- In case of dryrun, 'outputs' contains the text of the NJOY input
        if dryrun:
            msg = "dryrun requested — returning NJOY input deck"
            log_stage(log, method, zam, msg, verbose=verbose)
            return outputs

        msg = "parsing NJOY GENDF output into Endf6 structure"
        log_stage(log, method, zam, msg, verbose=verbose)

        # Parse GENDF text into object
        return Gendf.from_text(outputs["gendf"])

    @with_optional_warning_suppression("sandy.warn", default_suppress=False)
    def get_errorr(
            self,
            *,
            nubar: bool | None = None, # None means "auto"; bool means user override
            mubar: bool | None = None,
            chi: bool | None = None,
            xs: bool | None = None,
            dryrun: bool = False,
            groupr_kws: dict | None = None,
            errorr_kws: dict | None = None,
            errorr31_kws: dict | None = None,
            errorr33_kws: dict | None = None,
            errorr34_kws: dict | None = None,
            errorr35_kws: dict | None = None,
            suppress_njoy_output: bool = False,
            suppress_warnings: bool | None = None,
            print_njoy_input: bool = False,
            verbose: bool = False,
            **njoy_kws,
            ):
        """
        Process covariance data for this ENDF‑6 evaluation using NJOY's ERRORR
        module (with optional GROUPR pre‑processing).
    
        Parameters
        ----------
        nubar : bool or None, default None
            Request processing of MF=31 (nubar covariance).
            - ``None`` → enable automatically if MF=31 exists.
            - ``True`` → force enable.
            - ``False`` → force disable.
    
        mubar : bool or None, default None
            Request processing of MF=34 (angular distribution covariance).
            Follows the same auto/override behavior as ``nubar``.
    
        chi : bool or None, default None
            Request processing of MF=35 (χ covariance).
            Follows the same auto/override behavior as ``nubar``.
    
        xs : bool or None, default None
            Request processing of MF=33 (cross‑section covariance).
            Follows the same auto/override behavior as ``nubar``.
    
        dryrun : bool, default False
            If ``True``, return the generated NJOY input deck as a string instead of
            executing NJOY.
    
        groupr_kws : dict, optional
            Keyword arguments forwarded to GROUPR when required
            (e.g., for nubar, mubar, chi).
            User values override internally inferred defaults.
    
        errorr_kws : dict, optional
            Base options applied to all ERRORR calls unless shadowed by a
            module-specific override.
    
        errorr31_kws, errorr33_kws, errorr34_kws, errorr35_kws : dict, optional
            Options passed only to the matching ERRORR submodule:
            - ``errorr31_kws`` → MF=31
            - ``errorr33_kws`` → MF=33
            - ``errorr34_kws`` → MF=34
            - ``errorr35_kws`` → MF=35
    
        suppress_njoy_output : bool, default False
            If ``True``, suppress raw NJOY standard output.
    
        suppress_warnings : bool or None, default None
            Whether to silence warnings raised during ERRORR processing.
            - ``None`` → use package default
            - ``True`` / ``False`` → explicit override
    
        print_njoy_input : bool, default False
            If ``True``, print the generated NJOY deck to screen before running.
    
        verbose : bool, default False
            Enable additional diagnostic information.
    
        **njoy_kws :
            Extra keyword arguments forwarded to the underlying NJOY runner.
    
        Returns
        -------
        dict or str
            If ``dryrun=True``:
                Returns the NJOY input deck as a ``str``.
    
            If ``dryrun=False``:
                Returns a ``dict`` mapping module names to
                :class:`~sandy.errorr.Errorr` objects, e.g.:
    
                ``{"errorr31": Errorr(...), "errorr33": Errorr(...), ...}``
    
                If no covariance channels are available or requested,
                returns an empty ``dict``.
    
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
        >>> out = endf6.get_errorr(temperature=300, minimal_processing=True, err=1, errorr_kws=dict(ign=3, mt=18), suppress_warnings=True, suppress_njoy_output=True)

        This test checks also the type of each output.

        >>> assert out["errorr33"].get_xs().data.shape[0] == 30
        >>> assert out["errorr31"].get_xs().data.shape[0] == 30
        >>> assert out["errorr34"].get_xs().data.shape[0] == 30
        >>> assert out["errorr33"].get_xs().data.shape[0] == 30

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

        Example: 91-Pa-231 in JEFF-3.3 has MF32 but not MF33.

        >>> tape = sandy.get_endf6_file("jeff_33", "xs", 912310, local=True)
        >>> err = tape.get_errorr(chi=False, nubar=True, mubar=False, err=1, xs=True, errorr33_kws=dict(irespr=0), suppress_warnings=True, suppress_njoy_output=True)
        >>> assert "errorr33" in err

        Example: 91-Pa-231 in JEFF-3.3 has MF32 but not MF33.

        >>> tape = sandy.get_endf6_file("jeff_33", "xs", 912330, local=True)
        >>> err = tape.get_errorr(chi=False, nubar=True, mubar=False, err=1, xs=True, errorr33_kws=dict(irespr=0), suppress_warnings=True, suppress_njoy_output=True)
        >>> assert "errorr33" in err

        Example: 17-Cl-37 in JEFF-3.3 has MF32 but not MF33.

        >>> tape = sandy.get_endf6_file("jeff_33", "xs", 170370, local=True)
        >>> err = tape.get_errorr(chi=False, nubar=True, mubar=False, err=1, xs=True, errorr33_kws=dict(irespr=0), suppress_warnings=True, suppress_njoy_output=True)
        >>> assert "errorr33" in err

        Example: 95-Am-241 in JEFF-3.3 has MF32 but not MF33.

        >>> tape = sandy.get_endf6_file("jeff_33", "xs", 952410, local=True)
        >>> err = tape.get_errorr(chi=False, nubar=True, mubar=False, err=1, xs=True, errorr33_kws=dict(irespr=0), suppress_warnings=True, suppress_njoy_output=True)
        >>> assert "errorr33" in err

        Check case when file does not contain covariance data.

        >>> tape = sandy.get_endf6_file("jeff_33", "xs", 10030, local=True)
        >>> outs = tape.get_errorr()
        >>> assert outs == {}

        Check case when file contains covariance data, but they are not requested.

        >>> tape = sandy.get_endf6_file("jeff_33", "xs", 10020, local=True)
        >>> outs = tape.get_errorr(xs=False, suppress_warnings=True)
        >>> assert outs == {}

        """
        # ---- IMPORT
        from tempfile import TemporaryDirectory
        from subprocess import DEVNULL
        import pprint

        from .njoy import _input_mf32_nomf33, _input_mf32_nomf33_no18, _run_njoy
        from .errorr import Errorr
        from .utils import log
        from ._perturbation_base import log_stage

        # ---- SETUP
        src = self  # change of variables to avoid overwriting self
        zam = src.get_zam()
        recs = src.get_records()
        mfs = recs.MF.unique()
        mts = recs.MT.unique()
        
        method = "get_errorr"
        
        msg = f"loaded ENDF records: MF present = {mfs}"
        log_stage(log, method, zam, msg, verbose=verbose)
        
        if set([31, 32, 33, 34, 35]).isdisjoint(mfs):
            msg = "no processable covariance section was found"
            log_stage(log, method, zam, msg, verbose=verbose)
            return {}
            

        # ---- HANDLE MF32-no-MF33 cases
        # this replaces tye (ld @handle_mf32_alone decorator
        if 32 in mfs and 33 not in mfs:
            msg = "detected MF32 without MF33: using synthetic ERRORR33 templates"
            log_stage(log, method, zam, msg, verbose=verbose)

            # input taken from
            # https://www-nds.iaea.org/index-meeting-crp/TM_NDP/docs/OCabellos_2017.pdf
            input_fiss = _input_mf32_nomf33
            input_nofiss = _input_mf32_nomf33_no18

            # choose MF32 template depending on whether fission exists
            inp = input_fiss if 18 in mts else input_nofiss

            addon_msg = "fission" if 18 in mts else "non‑fission"
            msg = f"using MF32 handling template: {addon_msg}"
            log_stage(log, method, zam, msg, verbose=verbose)

            with TemporaryDirectory() as td:
                f20 = os.path.join(td, "tape20")
                src.to_file(f20)
                njoy_output = DEVNULL if suppress_njoy_output else None
                outs = _run_njoy(inp, f20, njoy_output=njoy_output)

            # Replace self with synthetic ERRORR33-only partial tape
            src = Endf6.from_text(outs["errorr33"])
            recs = src.get_records()

        # ---- NORMALIZE dict-like keyword arguments
        msg = "augmenting ERRORR NJOY kwargs"
        log_stage(log, method, zam, msg, verbose=verbose)

        njoy_kws_ = src._prepare_njoy_kws(**njoy_kws)
        njoy_kws_["dryrun"] = dryrun
        njoy_kws_["acer"] = False

        # -- prepare/augment GROUPR options without mutating the user's dict --
        msg = "augmenting GROUPR NJOY kwargs"
        log_stage(log, method, zam, msg, verbose=verbose)

        groupr_kws_ = (groupr_kws or {}).copy()
        njoy_kws_["groupr_kws"] = src._prepare_groupr_kws(**groupr_kws_)

        # -- prepare/augment ERRORR options without mutating the user's dict --
        errorr_kws_ = (errorr_kws or {}).copy()

        # ---- PREPARE NJOY OPTIONS
        # Activate specific errorr module according to covariance info and input options
        has31 = not recs.loc[recs.MF == 31].empty  # nubar cov
        has33 = not recs.loc[recs.MF == 33].empty  # xs cov
        has34 = not recs.loc[recs.MF == 34].empty  # mubar cov
        has35 = not recs.loc[recs.MF == 35].empty  # chi cov

        msg = f"covariance availability: MF31={has31} MF33={has33} MF34={has34} MF35={has35}"
        log_stage(log, method, zam, msg, verbose=verbose)

        # Switch off if user provides False
        use31 = has31 if nubar is None else has31 & bool(nubar)
        use33 = has33 if xs is None else has33 & bool(xs)
        use34 = has34 if mubar is None else has34 & bool(mubar)
        use35 = has35 if chi is None else has35 & bool(chi)

        msg = f"covariance processing: MF31={use31} MF33={use33} MF34={use34} MF35={use35}"
        log_stage(log, method, zam, msg, verbose=verbose)

        if not any([use31, use33, use34, use35]):
            msg = "no processable covariance section was requested"
            log_stage(log, method, zam, msg, verbose=verbose)
            return {}

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

        # ---- SUPPRESSING NJOY output (optional)
        if suppress_njoy_output:
            msg = "NJOY output to screen is suppressed"
            log_stage(log, method, zam, msg, verbose=verbose)

            njoy_kws_ |= {
                "njoy_output": DEVNULL
                }


        pdict = pprint.pformat(njoy_kws_, indent=2, sort_dicts=True)
        msg = f"running NJOY with augmented kwargs: {pdict}"
        log_stage(log, method, zam, msg, verbose=verbose)

        # ---- RUN NJOY via the shared helper
        msg = "run NJOY"
        log_stage(log, method, zam, msg, verbose=verbose)

        outputs = src._run_njoy(
            print_njoy_input=print_njoy_input,
            verbose=verbose,
            **njoy_kws_,
            )   # if dryrun in njoy_kws, this will return the NJOY input without running NJOY

        # --- In case of dryrun, 'outputs' contains the text of the NJOY input
        if dryrun:
            msg = "dryrun requested - returning NJOY input deck"
            log_stage(log, method, zam, msg, verbose=verbose)
            return outputs

        # ---- MAP OUTPUTS into Errorr objects
        msg = "parsing ERRORR outputs"
        log_stage(log, method, zam, msg, verbose=verbose)

        outputs = {k: Errorr.from_text(
            v) for k, v in outputs.items() if k.startswith("errorr")}

        msg = f"produced ERRORR objects: {outputs.keys()}"
        log_stage(log, method, zam, msg, verbose=verbose)

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

    def get_perturbations(
            self,
            *args,
            verbose=False,
            **kwargs
            ):
        """
        Dispatcher to assign perturbations method: either for radioactive
        decay data, fission yields or cross section.

        Parameters
        ----------
        verbose : bool, optional, default is False
            It will log time messages and be bassed to the called method.

        Notes
        -----
        .. note :: The perturbation method is selected based on the MT's found
                   in `self`.
        """
        # ---- IMPORT
        from .utils import log
        from ._perturbation_base import log_stage

        # ---- SETUP
        msg = (
            "########################################################\n"
            "                GET PERTURBATIONS                       \n"
            "########################################################"
            )
        log(msg, verbose=verbose)

        t0 = time.perf_counter()

        # ---- CHOOSE METHOD
        if 457 in self.mt:
            method = "get_perturbations_rdd"
            out = self.get_perturbations_rdd(*args, verbose=verbose, **kwargs)

        elif 454 in self.mt:
            method = "get_perturbations_fy"
            out = self.get_perturbations_fy(*args, verbose=verbose, **kwargs)

        else:
            method = "get_perturbations_xs"
            out = self.get_perturbations_xs(*args, verbose=verbose, **kwargs)

        # ---- LOGGING: end
        dt = time.perf_counter() - t0
        msg = f"finished in {dt:.3f} s"
        log_stage(log, method, None, msg, verbose=verbose)

        return out

    @with_optional_warning_suppression("sandy.warn", default_suppress=True)
    def get_perturbations_xs(
            self,
            nsmp: int,
            njoy_kws: dict | None = None,
            smp_kws: dict | None = None,
            suppress_njoy_output: bool = True,
            suppress_warnings: bool | None = None,
            verbose: bool = False,
            write: bool = True,
            write_errorr: bool = True,
            write_samples: bool = True,
            **kwargs,
            ) -> dict[int, ]:
        """
        Generate multigroup perturbation samples for::
            - cross sections,
            - nubar,
            - secondary neutron angular distributions (not yet implemented),
            - secondary neutron energy distributions 

        using covariance information processed via NJOY ERRORR.
    
        This method extracts covariance matrices from the MF=31, 33, and 35
        sub‑sections produced by :meth:`~sandy.endf6.Endf6.get_errorr`, 
        constructs multivariate distributions with unit mean and relative 
        covariance, and samples perturbation factors with the same multigroup 
        structure as the covariance matrices.
        
        Samples are returned as :obj:`~sandy.samples.Samples` objects,
        grouped by MF into a `dict`.
    
        Optionally, the raw ERRORR tapes (ASCII text) and sample spreadsheets
        (EXCEL) can be written to disk.

        Parameters
        ----------
        nsmp : int
            Number of perturbation samples to generate.
        njoy_kws : dict, optional
            Keyword arguments forwarded to
            :meth:`~sandy.errorr.Endf6.get_errorr` to control the NJOY ERRORR
            processing.
            Keys such as ``errorr31_kws``, ``errorr33_kws`` and
            ``errorr35_kws`` may be used to select MTs or energy grids.
            The dictionary is internally copied and not modified.
        smp_kws : dict, optional
            Additional keyword arguments forwarded to
            :meth:`~sandy.cov.CategoryCov.sampling`. For reproducibility,
            MF‑specific seeds may be passed via keys like ``"seed33"`` or
            ``"seed35"``. The dictionary is internally copied.
        suppress_njoy_output : bool, default=True
            If True, NJOY output is redirected to ``subprocess.DEVNULL``.
        suppress_warnings : bool, default=None
            If True, suppress warnings emitted during ERRORR processing.
            If not given, apply default (True).
        verbose : bool, default=False
            If True, print detailed progress and diagnostic messages.
            This is decoupled from the :meth:`~sandy.errorr.Endf6.get_errorr`
            logging wich is disbaled by default.
            To enable it, argument `verbose=True` shoould aslos be passed to
            `njoy_kws`.
        write : bool, default=True
            Master switch controlling whether output files are written.
            Affects both ERRORR tapes and sample spreadsheets.
        write_errorr : bool, default=True
            Whether to write the raw ``ERRORR_*_MFxx.tape`` files. Ignored
            if ``write=False``.
        write_samples : bool, default=True
            Whether to write the sample spreadsheets
            ``PERT_*_MFxx.xlsx``. Ignored if ``write=False``.
    
        Returns
        -------
        dict of int to :class:`~sandy.samples.Samples`
            A dictionary mapping ENDF MF numbers to
            :class:`~sandy.samples.Samples` objects. Possible keys are 
            (if present):
    
            * ``31`` — Cross section covariances
            * ``33`` — Reaction cross section covariances
            * ``35`` — Fission spectrum & fission‑related covariances
    
            The dictionary is empty if no covariance MF sections are
            available or if all ERRORR channels were disabled in ``njoy_kws``
            (keyword arguments `nubar=False`, `xs=False`, `chi=False`).


        Notes
        -----
        * The method does **not** modify the user‑provided ``njoy_kws`` or
          ``smp_kws`` dictionaries.
        * If the ENDF material contains no covariance sections
          (MF 31–35), the method returns an empty dictionary.
        * When ``write=True``, output files are written in the current
          working directory using the naming convention:
    
          - ``ERRORR_<ZAID>_MF<MF>.tape``
          - ``PERT_<ZAID>_MF<MF>.xlsx``

        * Logging is formatted as  ``{method} | ZAM={ZAM} | {info}``

        Examples
        --------
        Clean up for later testing of `write` keyword.
        
        >>> # Clean up H1 files
        >>> from pathlib import Path
        >>> outdir = Path.cwd()
        >>> p_err_h1 = outdir / "ERRORR_1001_MF33.tape"
        >>> p_pert_h1 = outdir / "PERT_1001_MF33.xlsx"
        >>> if p_err_h1.exists(): p_err_h1.unlink()
        >>> if p_pert_h1.exists(): p_pert_h1.unlink()
        >>> assert (not p_err_h1.exists()) & (not p_pert_h1.exists())

        >>> # Clean up H2 files
        >>> from pathlib import Path
        >>> outdir = Path.cwd()
        >>> p_err_h2 = outdir / "ERRORR_1002_MF33.tape"
        >>> p_pert_h2 = outdir / "PERT_1002_MF33.xlsx"
        >>> if p_err_h2.exists(): p_err_h2.unlink()
        >>> if p_pert_h2.exists(): p_pert_h2.unlink()
        >>> assert (not p_err_h2.exists()) & (not p_pert_h2.exists())

        >>> # Clean up U235 files
        >>> p_err = { mf : outdir / f"ERRORR_92235_MF{mf}.tape" for mf in (31, 33, 35) }
        >>> p_pert = { mf : outdir / f"PERT_92235_MF{mf}.xlsx" for mf in (33, 35, 35) }
        >>> for p in p_err.values():
        ...     if p.exists():
        ...         p.unlink()
        ...     assert not p.exists()

        Test get perturbations for MF 33.
        Generate a couple of samples from the H1 file of JEFF-3.3.

        >>> import sandy, numpy as np
        >>> njoy_kws = dict(err=1, errorr_kws=dict(mt=102))
        >>> sample_size = 2
        >>> tape = sandy.get_endf6_file("jeff_33", "xs", 10010, local=True)
        >>> smps = tape.get_perturbations_xs(sample_size, njoy_kws=njoy_kws)

        Few checks on the output.

        >>> assert len(smps) == 1   # only MF33 is present
        >>> assert isinstance(smps[33], sandy.samples.Samples)
        
        The MT selection worked.

        >>> assert (smps[33].data.index.get_level_values("MT") == 102).all()

        By default, writing is active. Test files are created with correct name.

        >>> assert p_err_h1.exists()
        >>> assert p_pert_h1.exists()

        Test get perturbations from MF 35.
        Generate a couple of samples from the U235 file of JEFF-3.3.

        >>> njoy_kws = dict(err=1, errorr33_kws=dict(mt=18))
        >>> sample_size = 2
        >>> tape = sandy.get_endf6_file("jeff_33", "xs", 922350, local=True)
        >>> smps = tape.get_perturbations_xs(nsmp=sample_size, njoy_kws=njoy_kws, write=False)

        Few checks on the output.

        >>> assert len(smps) == 3   # MF35, MF33 and MF31 all perturbed
        >>> assert isinstance(smps[35], sandy.samples.Samples)
        >>> assert (smps[35].data.index.get_level_values("MT") == 18).all()

        Writing was deactivated. Files should not exist.

        >>> for p in p_err.values():
        ...     assert not p.exists()

        By redirecting njoy outputs to screen `njoy_kws` are locally modified.
        Also `mubar=False` is added.
        Check that they do not mutate outside the method.

        >>> assert njoy_kws == dict(err=1, errorr33_kws=dict(mt=18))
        
        Suppress all MF from `get_errorr`. This is tested on the H2 file of JEFF-3.3.
        
        >>> njoy_kws = dict(err=1, nubar=False, mubar=False, chi=False, xs=False)
        >>> sample_size = 2
        >>> tape = sandy.get_endf6_file("jeff_33", "xs", 10020, local=True)
        >>> smps = tape.get_perturbations_xs(nsmp=sample_size, njoy_kws=njoy_kws, write=True)

        This should nicely return an empty dictionary.

        >>> assert smps == {}

        And no file should be created.

        >>> assert p_err_h1.exists()
        >>> assert p_pert_h1.exists()
        
        The same behavior is expected for a file without covariance data.
        This is tested on the H3 file of JEFF-3.3.

        >>> tape = sandy.get_endf6_file("jeff_33", "xs", 10030, local=True)
        >>> assert 31 not in tape.mf
        >>> assert 32 not in tape.mf
        >>> assert 33 not in tape.mf
        >>> assert 34 not in tape.mf
        >>> assert 35 not in tape.mf

        This should nicely return an empty dictionary.

        >>> sample_size = 2
        >>> njoy_kws = dict(err=1)
        >>> smps = tape.get_perturbations_xs(nsmp=sample_size, njoy_kws=njoy_kws, write=True)
        >>> assert smps == {}

        """
        # ---- IMPORT
        from subprocess import DEVNULL
        from pathlib import Path
        import pprint
        
        from .samples import Samples
        from .utils import log, get_seed
        from ._perturbation_base import log_stage

        # ---- NORMALIZE dict-like keyword arguments
        smp_kws_ = {} if smp_kws is None else smp_kws.copy()
        njoy_kws_ = {} if njoy_kws is None else njoy_kws.copy()
        # do not mutate caller dicts; copy and augment
        # switch off mubar in ERRORR unless user explicitly set it
        njoy_kws_.setdefault("mubar", False)

        # optionally silence NJOY output
        if suppress_njoy_output:
            njoy_kws_["njoy_output"] = DEVNULL

        # ---- SETUP
        zam = self.get_zam()
        method = "get_perturbations_xs"
        
        
        # prepare output directory and basename
        outdir_path = Path.cwd()

        # filename templates (per MF)
        base = str(self.get_id())   # this would be 92235 for U235 and 95642 for Am242m, need string conversion
        fn_errorr = lambda mf: outdir_path / f"ERRORR_{base}_MF{mf}.tape"
        fn_smp    = lambda mf: outdir_path / f"PERT_{base}_MF{mf}.xlsx"

        # ---- PRODUCE ERRORR files with covariance data
        pdict = pprint.pformat(njoy_kws_, indent=2, sort_dicts=True)
        msg = f"run ERRORR via get_errorr({pdict})"
        log_stage(log, method, zam, msg, verbose=verbose)
        
        # do not print NJOY input to screen
        outs = self.get_errorr(**njoy_kws_)

        smp: dict[int, "Samples"] = {}

        # ---- CHECK MF FOUND
        def _parse_mf(k) -> int:
            """"extract XX from key 'errorrXX', like 33 from 'errorr33'"""
            s = str(k)
            return int(s[-2:])
        
        mfs = [ _parse_mf(k) for k in outs.keys() ]
        msg = f"covariance MFs available={mfs}"
        log_stage(log, method, zam, msg, verbose=verbose)

        # -- Extract samples from covariance data, iterate over MF31, 33 and 35
        for mf in mfs:

            key = f"errorr{mf}"
            out = outs[key]

            # ---- WRITE raw ERRORR tape (optional)
            if write and write_errorr:
                err_path = fn_errorr(mf)
                msg = f"MF={mf:d} | writing ERRORR -> '{err_path}'"
                log_stage(log, method, zam, msg, verbose=verbose)
                
                if err_path.exists():
                    msg = f"MF={mf:d} | file exists and will be overwritten"
                    log_stage(log, method, zam, msg, verbose=verbose)

                out.to_file(err_path)

            # ---- EXTRACT covariance matrix
            msg = f"MF={mf:d} | extracting covariance matrix"
            log_stage(log, method, zam, msg, verbose=verbose)

            cov = out.get_cov()   # CategoryCov object

            if cov.data.empty:
                msg = f"MF={mf:d} | covariance matrix is empty"
                log_stage(log, method, zam, msg, verbose=verbose)

            else:
                msg = f"MF={mf:d} | covariance matrix size={cov.data.shape}"
                log_stage(log, method, zam, msg, verbose=verbose)

                mts = cov.data.columns.get_level_values("MT").unique().to_numpy()
                msg = f"MF={mf:d} | found MT numbers={mts}"
                log_stage(log, method, zam, msg, verbose=verbose)

            # ---- EXTRACT sample
            seed_key = f"seed{mf}"
            # generated seed here to be able to log it...don;t let cov.sampling do it
            seed_ = smp_kws_.get(seed_key)
            if seed_ is not None:
                seed = seed_ 
                msg = f"MF={mf:d} | explicit seed provided"
            else:
                seed = get_seed()
                msg = f"MF={mf:d} | explicit seed not provided"
            log_stage(log, method, zam, msg, verbose=verbose)

            # remove seed** from smp_kws_, they are just for the pipeline
            for key in ["seed31", "seed33", "seed34", "seed35"]:
                if key in smp_kws_:
                    smp_kws_.pop(key)

            msg = f"MF={mf:d} | sampling with SMP size={nsmp} via sampling(seed={seed}, {smp_kws_})"
            log_stage(log, method, zam, msg, verbose=verbose)
            
            # draw samples (unit mean, relative covariance)
            smp[mf] = cov.sampling(nsmp, seed=seed, **smp_kws_)  # sandy.Samples


            # ---- WRITE SAMPLE to XLSX with stats (optional)
            if write and write_samples:
                xls_path = fn_smp(mf)
                msg = f"MF={mf:d} | writing samples -> '{xls_path}'"
                log_stage(log, method, zam, msg, verbose=verbose)

                if xls_path.exists():
                    msg = f"MF={mf:d} | file exists and will be overwritten"
                    log_stage(log, method, zam, msg, verbose=verbose)

                # write the long-form sample frame
                smp[mf].to_excel(xls_path)

        msg = f"done | SMP size={nsmp} | MF={mfs} | write={write}"
        log_stage(log, method, zam, msg, verbose=verbose)

        return smp

    @with_optional_warning_suppression("sandy.warn", default_suppress=False)
    def get_perturbations_rdd(
            self,
            nsmp: int,
            *,
            rdd = None,
            fill_zeros_decay_energy: float | None = None,
            fill_zeros_half_life: float | None = None,
            fill_zeros_branching_ratio: float | None = None,
            smp_hl_kws: dict | None = None,
            smp_de_kws: dict | None = None,
            smp_br_kws: dict | None = None,
            verbose: bool = False,
            write: bool = True,
            **kwargs,
            ) -> dict[str, ]:
        """
        Generate perturbation samples for radioactive decay data (half-lives,
        decay energies, and branching ratios) using relative uncertainties derived
        from the evaluated decay data stored in ``self``.
    
        This method constructs multivariate distributions with mean unity and
        covariance matrices derived from ENDF-6 decay uncertainties. The output is
        a dictionary of :class:`~sandy.samples.Samples` objects.
    
        Parameters
        ----------
        nsmp : int
            Number of samples to generate.
    
        rdd : :class:`~sandy.decay.DecayData`, optional
            Pre-loaded decay data. If not provided, it is extracted automatically
            via :meth:`~sandy.decay.DecayData.from_endf6`.
    
        fill_zeros_decay_energy : float or None, optional
            Fill null decay energy uncertainties with a default value.
            Example: ``fill_zeros_decay_energy=0.05`` to add a 5% uncertainty.
    
        fill_zeros_half_life : float or None, optional
            Fill null half life uncertainties with a default value.
    
        fill_zeros_branching_ratio : float or None, optional
            Fill null branching ratio uncertainties with a default value.
    
        smp_hl_kws, smp_de_kws, smp_br_kws : dict, optional
            Keyword arguments passed to :meth:`~sandy.cov.CategoryCov.sampling`
            for half-lives, decay energies, and branching ratios respectively.
            Pass any seed value here using key ``seed``.
    
        verbose : bool, optional
            Enable detailed logging of intermediate steps.
            Default is ``False``.
    
        write : bool, optional
            If ``True``, write sampled perturbations to an Excel file named
            ``PERT_MF8_MF457.xlsx`` in the current working directory.
            Default is ``True``.
    
        Returns
        -------
        dict[str, Samples]
            A dictionary containing perturbation samples with keys:
    
            - ``"HL"`` — half-life perturbations  
            - ``"DE"`` — decay energy perturbations  
            - ``"BR"`` — branching ratio perturbations  
            
            The dictionary is empty if all nuclides are stable.
    
        Notes
        -----
        - Branching ratios are sampled without correlations and must be
          renormalized afterwards.
    
        Examples
        --------
        Sample Co-59 and Co-60 decay data from JEFF-3.3
        
        First read the file.
        
        >>> import sandy, numpy as np
        >>> decay = sandy.get_endf6_file("jeff_33", "decay", [270590, 270600], local=True)

        Draw with a large sample size, to ensure convergence in the checks.
        Also, use seeds for reproducibility.

        >>> kws = {"smp_hl_kws": {"seed": 3}, "smp_de_kws": {"seed": 3}, "smp_br_kws": {"seed": 3}}
        >>> sample_size = 10000
        >>> smps = decay.get_perturbations_rdd(sample_size, write=False, **kws)
        
        Check that the output mapping keys are correct.
        
        >>> assert smps.keys() == {"HL", "DE", "BR"}

        Check mean and std converged for half lives (the sample size should guarantee it).

        >>> smp_mean, smp_std = smps["HL"].get_mean(), smps["HL"].get_std()
        >>> expected = [1, 1]
        >>> np.testing.assert_array_almost_equal(smp_mean.to_numpy(), expected, decimal=4)
        >>> expected = [0, 1.5e-4]
        >>> np.testing.assert_array_almost_equal(smp_std.to_numpy(), expected, decimal=5)

        Check mean and std converged for decay energies (the sample size should guarantee it).

        >>> smp_mean, smp_std = smps["DE"].get_mean(), smps["DE"].get_std()
        >>> expected = [1] * 6
        >>> np.testing.assert_array_almost_equal(smp_mean.to_numpy(), expected, decimal=4)
        >>> expected = [0] * 4 + [0.002098, 0.000141]
        >>> np.testing.assert_array_almost_equal(smp_std.to_numpy(), expected, decimal=5)

        Branching ratios did not contain uncertainty, so they are returned as unperturbed.
        
        >>> assert all(smps["BR"].data.squeeze().to_numpy() == 1)

        The same output can be produced by passing a :class:`~sandy.decay.DecayData` 
        instance.

        >>> rdd = sandy.DecayData.from_endf6(decay)
        >>> smps_rdd = decay.get_perturbations_rdd(sample_size, write=False, rdd=rdd, **kws)
        
        Outputs must be the same. The purpose of the keyword is not to extract 
        once again the decay data, if already done.

        >>> for k in smps:
        ...    assert smps_rdd[k].data.equals(smps[k].data)
        
        Outputs are different if the same seeds are not given.

        >>> smps_2 = decay.get_perturbations_rdd(sample_size, write=False, rdd=rdd)

        >>> assert not smps_2["HL"].data.equals(smps["HL"].data)
        >>> assert not smps_2["DE"].data.equals(smps["DE"].data)

        Branching ratios are not perturbed, so they don't change.

        >>> assert smps_2["BR"].data.equals(smps["BR"].data)

        Let's try the writing option. But first I clean up existing files.

        >>> # Clean up PERT files
        >>> from pathlib import Path
        >>> outdir = Path.cwd()
        >>> p = outdir / "PERT_MF8_MT457.xlsx"
        >>> if p.exists(): p.unlink()
        >>> assert not p.exists()

        Now run again (with a smaller size, not to run too much).
        The PERT file must have been created.

        >>> sample_size = 2
        >>> smps = decay.get_perturbations_rdd(sample_size, write=True, rdd=rdd)
        >>> assert p.exists()

        If all nuclides are stable there is nothing to sample.
        
        >>> decay = sandy.get_endf6_file("jeff_33", "decay", 270590, local=True)
        >>> sample_size = 2
        >>> smps = decay.get_perturbations_rdd(sample_size, suppress_warnings=True)

        It returns an empty ``dict``.

        >>> assert smps == {}


        Here the keys ``fill_zeros_*`` are tested. File for H5 in JEFF-3.3
        gives all decay data without uncertainty.
        Then, we introduce arbitrary uncertainties on all parameters.

        >>> decay = sandy.get_endf6_file("jeff_33", "decay", 10050, local=True)
        >>> sample_size = 1000
        >>> smps = decay.get_perturbations_rdd(sample_size, write=False, fill_zeros_half_life=0.1,
        ...                                    fill_zeros_decay_energy=0.2, fill_zeros_branching_ratio=0.5)

        Samples are produced with variability due to the provided uncertainty.

        >>> import math
        >>> rtol = 0.05
        >>> key = "HL"
        >>> expected_rstd = 0.1
        >>> smp_rst = smps[key].get_rstd().squeeze()  # scalar
        >>> assert math.isclose(expected_rstd, smp_rst, rel_tol=rtol)
        >>> key = "DE"
        >>> expected_rstd = [0.2, 0, 0]
        >>> smp_rst = smps[key].get_rstd().to_numpy()  # array
        >>> assert np.allclose(expected_rstd, smp_rst, rtol=rtol)
        >>> key = "BR"
        >>> expected_rstd = 0.5
        >>> smp_rst = smps[key].get_rstd().squeeze()  # scalar
        >>> assert math.isclose(expected_rstd, smp_rst, rel_tol=rtol)

        Without using the ``fill_zeros_*`` keywords, unit perturbation
        coefficients are returned.

        >>> smps = decay.get_perturbations_rdd(sample_size, write=False)
        >>> assert np.all(smps["HL"].data == 1)
        >>> assert np.all(smps["DE"].data == 1)
        >>> assert np.all(smps["BR"].data == 1)
        """
        # ---- IMPORT
        from pathlib import Path
        import pandas as pd
        import numpy as np

        from .decay import DecayData
        from .cov import CategoryCov
        from .samples import Samples, FILENAME_RDD_PERT
        from .utils import log, get_seed
        from ._perturbation_base import log_stage

        # ---- SETUP
        # there is likely no single zam
        zam = self.get_zam()
        method = "get_perturbations_rdd"

        length = 1 if np.isscalar(zam) else len(zam)
        
        msg = f"found {length} ZAM"
        log_stage(log, method, zam, msg, verbose=verbose)


        # ---- NORMALIZE dict-like arguments
        smp_hl_kws_ = (smp_hl_kws or {}).copy()
        smp_de_kws_ = (smp_de_kws or {}).copy()
        smp_br_kws_ = (smp_br_kws or {}).copy()

        # ---- PREPARE DecayData object
        # if already available in kwargs, do not extract DecayData again
        status = "provided" if rdd is not None else "extracted with DecayData.from_endf6"
        msg = f"DecayData={status}"
        log_stage(log, method, zam, msg, verbose=verbose)

        # pass verbosity
        rdd_ = rdd if rdd is not None else DecayData.from_endf6(self, verbose=verbose)

        if all(v["stable"] for v in rdd_.data.values()):
            msg = "cannot sample perturbations: all nuclides are stable"
            warn_logger = logging.getLogger("sandy.warn")
            log(msg, level=logging.WARNING, logger=warn_logger)
            return {}

        # ---- HELPER fcuntion
        def _sample_category_with_logging(
            *,
            name: str,
            rel_unc: pd.Series,
            fill_zeros: float,
            kwargs: dict,
            ) -> Samples:
            """
            Helper for sampling one decay-category (HL/DE/BR) with logging,
            seed handling and zero-uncertainty replacement.
            """
            # --- ZERO UNCERTAINTY CHECK ---
            zeros_found = (rel_unc == 0).sum()
            msg = f"{name} | null uncertainty in {zeros_found} entries"
            if fill_zeros > 0:
                msg += f" increased to {fill_zeros*100:.1f} %"
                rel_unc = rel_unc.replace(0, fill_zeros)
            log_stage(log, method, zam, msg, verbose=verbose)
        
            # --- Replace NaN generated from divisions ---
            rel_unc = rel_unc.fillna(0)
        
            # --- SEED HANDLING ---
            # generated seed here to be able to log it...don't let cov.sampling do it
            if "seed" in kwargs:
                seed = kwargs.pop("seed")        # user-provided seed
                msg = f"{name} | explicit seed provided"
            else:
                seed = get_seed()                # auto seed
                msg = f"{name} | no explicit seed provided"

            log_stage(log, method, zam, msg, verbose=verbose)

            msg = f"{name} | sampling with SMP size={nsmp} via sampling(seed={seed}, {kwargs})"
            log_stage(log, method, zam, msg, verbose=verbose)

            # --- ACTUAL SAMPLING ---
            return CategoryCov.from_stdev(rel_unc).sampling(nsmp, seed=seed, **kwargs)


        # ---- SAMPLE HALF LIVES
        hl = rdd_.get_half_life()            # sandy.decay.DecayData
        dhl = hl.data.DHL / hl.data.HL       # pd.Series
        fill_hl = fill_zeros_half_life if fill_zeros_half_life is not None else 0
        smp_hl = _sample_category_with_logging(name="HALF LIVES", rel_unc=dhl, fill_zeros=fill_hl, kwargs=smp_hl_kws_)  # sandy.Samples

        # ---- SAMPLE DECAY ENERGY
        de = rdd_.get_decay_energy()            # sandy.decay.DecayEnergy
        dde = de.data.DE / de.data.E            # pd.Series
        fill_de = fill_zeros_decay_energy if fill_zeros_decay_energy is not None else 0
        smp_de = _sample_category_with_logging(name="DECAY ENERGIES", rel_unc=dde, fill_zeros=fill_de, kwargs=smp_de_kws_)

        # ---- SAMPLE BRANCHING RATIO
        br = rdd_.get_branching_ratio()         # sandy.decay.BranchingRatio
        dbr = br.data.DBR / br.data.BR          # pd.Series
        fill_br = fill_zeros_branching_ratio if fill_zeros_branching_ratio is not None else 0
        smp_br = _sample_category_with_logging(name="BRANCHING RATIOS", rel_unc=dbr, fill_zeros=fill_br, kwargs=smp_br_kws_)


        # ---- WRITE SAMPLE to XLSX with stats (optional)
        if write:
            # prepare output directory and basename
            outdir_path = Path.cwd()
            xls_path = outdir_path / FILENAME_RDD_PERT

            msg = f"writing samples -> '{xls_path}'"
            log_stage(log, method, zam, msg, verbose=verbose)

            if xls_path.exists():
                msg = "file exists and will be overwritten"
                log_stage(log, method, zam, msg, verbose=verbose)

            with pd.ExcelWriter(xls_path, engine="openpyxl") as writer:
                sheet_name = 'HALF LIFE'
                msg = f"HALF LIVES written in sheet named '{sheet_name}'"
                log_stage(log, method, zam, msg, verbose=verbose)
                smp_hl.data.to_excel(writer, sheet_name=sheet_name)

                sheet_name = 'DECAY ENERGY'
                msg = f"DECAY ENERGIES written in sheet named '{sheet_name}'"
                log_stage(log, method, zam, msg, verbose=verbose)
                smp_de.data.to_excel(writer, sheet_name=sheet_name)

                sheet_name = 'BRANCHING RATIO'
                msg = f"BRANCHING RATIOS written in sheet named '{sheet_name}'"
                log_stage(log, method, zam, msg, verbose=verbose)
                smp_br.data.to_excel(writer, sheet_name=sheet_name)

        # ---- CREATE OUTPUT MAPPING
        smp = {
            "HL": smp_hl,
            "DE": smp_de,
            "BR": smp_br,
        }

        msg = f"done | SMP size={nsmp} | # ZAM={length} | write={write}"
        log_stage(log, method, zam, msg, verbose=verbose)

        return smp

    @with_optional_warning_suppression("sandy.warn", default_suppress=False)
    def get_perturbations_fy(
            self,
            nsmp: int,
            nfpy = None,
            *,
            covariance: bool = False,
            smp_kws: dict | None = None,
            verbose: bool = False,
            write: bool = True,
            **kwargs,
            ) -> dict[str, ]:
        """
        Generate perturbation samples for independent fission yields (IFYs).
    
        This function builds multivariate distributions for IFY perturbation
        factors, with mean equal to unity, using either:
    
        - diagonal relative variances derived from the evaluated ENDF-6 data, or
        - (when available and explicitly requested) JEFF-4.0 thermal FY
          correlation matrices provided by CEA.
    
        Only MT=454 (independent fission yields) is treated. The returned samples
        represent *relative* perturbation coefficients applied to FY values.
    
        Parameters
        ----------
        nsmp : int
            Number of samples to generate.
    
        nfpy : sandy.Fy, optional
            Precomputed FY object. If not provided, it is extracted from the
            current ENDF-6 tape via ``Fy.from_endf6(self)``.
    
        covariance : bool, optional
            If ``True``, use JEFF-4.0 CEA thermal FY correlation matrices  
            (U‑233, U‑235, Pu‑239, Pu‑241) **when**:
            
            - the library is JEFF-4.0,
            - energy is thermal (0.0253 eV),
            - fissioning nuclide is one of the supported ZAM values.
    
            Otherwise, a diagonal covariance matrix (i.e., uncorrelated
            perturbations with correct variances) is used.  
            Default is ``False``.
    
        smp_kws : dict, optional
            Additional keyword arguments passed to
            :meth:`sandy.cov.CategoryCov.sampling` (e.g. seed specifications).
    
        verbose : bool, optional
            Enable progress logging.
    
        write : bool, optional
            If ``True``, write the generated perturbations to the file
            ``PERT_MF8_MT454.xlsx`` in the current working directory.
    
        Returns
        -------
        smps : dict
            A mapping with one entry:
    
                ``"IFY" → sandy.samples.Samples``
    
            The Samples object contains a dataframe with multi-index
            ``(ZAM, E, ZAP)`` and columns ``SMP`` representing individual samples.
    
            Each entry is a relative perturbation factor (mean ≈ 1).
    
        Notes
        -----
        - Only IFY covariance matrices for JEFF-4.0 thermal evaluations are
          available (U‑233, U‑235, Pu‑239, Pu‑241).
        - For all other cases, perturbations are uncorrelated but preserve FY
          relative standard deviations.
        - Sampling is block-wise per fissioning system (ZAM, E).
        - Seeds may be supplied per (ZAM, E) pair via ``smp_kws={"seed": {...}}``.


        Examples
        --------
        This test suite checks the reproducibility via keywords ``smp_kws={"seed": {}}``
        and ``nfpy``.

        >>> import sandy, numpy as np
        >>> seed_spec = {(922350, 0.0253): 1, (922350, 400e3): 4}
        >>> tape = sandy.get_endf6_file("jeff_33", "nfpy", 922350, local=True)

        After reading the file for one nuclide, a sample is generated.

        >>> sample_size = 2
        >>> smps = tape.get_perturbations_fy(sample_size, smp_kws=dict(seed=seed_spec), write=False)

        The process is repeated also passing the fy data and the same seed specs.

        >>> nfpy = sandy.Fy.from_endf6(tape)
        >>> smps2 = tape.get_perturbations_fy(sample_size, nfpy=nfpy, smp_kws=dict(seed=seed_spec), write=False)

        Since the seed is only given for thermal and fast fission (not high energy), 
        the resulting samples should be the same.

        >>> lower = smps["IFY"].data.query("E<1e7")
        >>> lower2 = smps2["IFY"].data.query("E<1e7")
        >>> assert lower2.equals(lower)

        However, they differ for the high energy fission yields.

        >>> higher = smps["IFY"].data.query("E>1e7")
        >>> higher2 = smps2["IFY"].data.query("E>1e7")
        >>> assert not higher2.equals(higher)



        This test suite checks the ``covariance`` option, which only works for U-235, 
        Pu-239, Pu-241 and U-233 thermal fission of JEFF-4.0
        Test ``covariance`` option.

        This is done by checking the sample correlation between nuclides ``zap=521350``
        and ``zap=531350``, which in the JEFF-4.0 covariance matrix is larger than 0.9
        in absolute value (it is -0.906028).

        The check is done for U-235 for JEFF-4.0. It only works for JEFF-4.0.

        >>> import sandy
        >>> tape = sandy.get_endf6_file("jeff_40", "nfpy", 922350, local=True)
        >>> nfpy = sandy.Fy.from_endf6(tape)

        With the covariance matrix the correlation should be larger than 0.8
        (took some margin for statistical noise).

        >>> sample_size = 50
        >>> smps = tape.get_perturbations_fy(sample_size, nfpy=nfpy, covariance=True, write=False)
        >>> corr = smps["IFY"].data.query("ZAP in [521350, 531350] & E==0.0253").T.corr()
        >>> assert np.abs(corr.iloc[0, 1]) > 0.8

        Without covariance matrix the correlation should be zero, but we accept
        some tolerance because of the small sample size.
        
        >>> sample_size = 50
        >>> smps = tape.get_perturbations_fy(sample_size, nfpy=nfpy, covariance=False, write=False)
        >>> corr = smps["IFY"].data.query("ZAP in [521350, 531350] & E==0.0253").T.corr()
        >>> assert np.abs(corr.iloc[0, 1]) < 0.5



        This test suite checks the writing option. But first I clean up existing files.

        >>> # Clean up PERT files
        >>> from pathlib import Path
        >>> outdir = Path.cwd()
        >>> p = outdir / "PERT_MF8_MT454.xlsx"
        >>> if p.exists(): p.unlink()
        >>> assert not p.exists()

        Now run again (with a smaller size, not to run too much).
        The PERT file must have been created.

        >>> sample_size = 2
        >>> smps = tape.get_perturbations_fy(sample_size, write=True, nfpy=nfpy)
        >>> assert p.exists()



        This test suite checks the sample convergence.

        First, we draw a large number of samples.

        >>> import sandy
        >>> tape = sandy.get_endf6_file("jeff_33", "nfpy", 922350, local=True)
        >>> sample_size = 1000
        >>> nfpy = sandy.Fy.from_endf6(tape)
        >>> smps = tape.get_perturbations_fy(sample_size, nfpy=nfpy, write=False)

        These are the expected results.

        >>> mean = nfpy.data.query("E==0.0253 and MT==454").set_index("ZAP").FY
        >>> std = nfpy.data.query("E==0.0253 and MT==454").set_index("ZAP").DFY
        >>> rstd = (std / mean).fillna(0)

        And these are the sample estimates.

        >>> smp_mean = smps["IFY"].get_mean().reset_index().query("E==0.0253").set_index("ZAP").MEAN
        >>> smp_rstd = smps["IFY"].get_std().reset_index().query("E==0.0253").set_index("ZAP").STD
        >>> smp_std = mean * smp_rstd

        Being perturbations relative, the mean of each one should converge to one.
        We also check that the mean variations across nuclides are minimal by 
        limiting the standard deviation of the statistical estimate.

        >>> assert np.isclose(smp_mean.mean(), 1, rtol=1e-2)
        >>> assert smp_mean.std() < 0.05

        To check the variance convergence we check the relative difference 
        between obtained and expected.

        >>> assert np.sum((smp_rstd - rstd)**2) / np.sum(rstd**2) < 0.05

        Then we also check that the largest variances are captured within 10%.

        >>> top = std.sort_values(ascending=False).head(100)
        >>> assert np.allclose(smp_std.loc[top.index], top, rtol=0.1)


        The convergence is also tested when sampling with covariance data.

        >>> import sandy, numpy as np
        >>> tape = sandy.get_endf6_file("jeff_40", "nfpy", 922350, local=True)
        >>> sample_size = 1000
        >>> smps = tape.get_perturbations_fy(sample_size, covariance=True, write=False)

        These are the expected results from the covariance source.

        >>> corr = sandy.fy.get_jeff40_fy_correlation_matrix(922350)
        >>> fy = sandy.Fy.from_endf6(tape).data.query("E==0.0253 and MT==454")
        >>> mean = fy.set_index("ZAP").FY
        >>> std = fy.set_index("ZAP").DFY
        >>> rstd = (std / mean).fillna(0)

        And these are the sample estimates.

        >>> smp_mean = smps["IFY"].get_mean().reset_index().query("E==0.0253").set_index("ZAP").MEAN
        >>> smp_rstd = smps["IFY"].get_std().reset_index().query("E==0.0253").set_index("ZAP").STD
        >>> smp_std = mean * smp_rstd

        Being perturbations relative, the mean of each one should converge to one.
        We also check that the mean variations across nuclides are minimal by 
        limiting the standard deviation of the statistical estimate.

        >>> assert np.isclose(smp_mean.mean(), 1, rtol=1e-2)
        >>> assert smp_mean.std() < 0.05

        To check the variance convergence we check the relative difference 
        between obtained and expected.

        >>> assert np.sum((smp_rstd - rstd)**2) / np.sum(rstd**2) < 0.05

        Then we also check that the largest variances are captured within 20%.

        >>> top = std.sort_values(ascending=False).head(100)
        >>> assert np.allclose(smp_std.loc[top.index], top, rtol=0.2)
        """
        # ---- IMPORT
        from pathlib import Path
        import pandas as pd
        import numpy as np

        from .cov import CategoryCov, corr2cov
        from .fy import Fy, get_jeff40_fy_correlation_matrix
        from .samples import Samples, FILENAME_FY_PERT
        from .utils import log, get_seed
        from ._perturbation_base import log_stage

        # ---- SETUP
        # there is likely no single zam
        zam = self.get_zam()
        mat_zam_mapping = self.get_mat_zam_mapping()
        zam_mat_mapping = {zam: mat for mat, zam in mat_zam_mapping.items()}
        method = "get_perturbations_fy"

        length = 1 if np.isscalar(zam) else len(zam)


        # ---- NORMALIZE dict-like arguments
        smp_kws_ = (smp_kws or {}).copy()


        # ---- PREPARE Fy object
        # if already available in kwargs, do not extract Fy again
        status = "provided" if nfpy is not None else "extracted with Fy.from_endf6"
        msg = f"Fy={status}"
        log_stage(log, method, zam, msg, verbose=verbose)

        # pass verbosity
        nfpy_ = nfpy if nfpy is not None else Fy.from_endf6(self, verbose=verbose)


        # ---- EXTRACT SEED specification
        seed_spec = smp_kws_.pop("seed", None)


        # ---- HELPER FOR SAMPLING
        def _sample_fy_with_logging(
                zam: int,
                e: float,
                fy,
                lib: str,
                kwargs: dict,
                ) -> Samples:
            """
            Build covariance, sample perturbations, and log.
            Returns a DataFrame with columns [ZAM, E, ZAP, SMP, VALS].
            """
            msg = f"E={e:3E} | processing IFY block"
            log_stage(log, method, zam, msg, verbose=verbose)

            block_id = (zam, e)
            # ---- DETERMINE LOCAL SEED
            if isinstance(seed_spec, dict):
                # Full explicit seed per FY block
                if block_id in seed_spec:
                    local_seed = seed_spec[block_id]
                    msg = f"E={e:3E} | explicit seed provided"
                else:
                    local_seed = get_seed()
                    msg = f"E={e:3E} | no explicit seed provided for this fissioning system"
        
            else:
                # No seeds at all
                local_seed = get_seed()
                msg = f"E={e:3E} | no explicit seed provided"

            log_stage(log, method, zam, msg, verbose=verbose)

            # ---- SELECT COVARIANCE MODEL
            EXPECTED_E = 0.0253
            ALLOWED_ZAM = [922330, 922350, 942390, 942410]
            if covariance and zam in ALLOWED_ZAM  and np.isclose(e, EXPECTED_E) and lib == "JEFF-4.0":

                msg = f"E={e:3E} | using JEFF-4.0 covariance matrix (with correlations) and fission yield data"
                log_stage(log, method, zam, msg, verbose=verbose)

                corr = get_jeff40_fy_correlation_matrix(zam)
            
                # ---- CONVERT correlation → covariance → relative covariance
                abs_cov = corr2cov(corr, fy.DFY)
                rel_cov = np.divide(abs_cov, fy.FY.to_numpy().reshape(-1, 1) @ fy.FY.to_numpy().reshape(1, -1))
                rcov = CategoryCov(rel_cov, index=fy.ZAP, columns=fy.ZAP)

            else:
                if covariance:
                    msg = f"E={e:3E} | covariance is requested but feature is not yet implemented"
                    log_stage(log, method, zam, msg, verbose=verbose)

                msg = f"E={e:3E} | using diagonal matrix (only variance)"
                log_stage(log, method, zam, msg, verbose=verbose)

                rstd = (fy.DFY / fy.FY).fillna(0)
                rcov = CategoryCov(
                    pd.DataFrame(np.diag(rstd**2), index=fy.ZAP, columns=fy.ZAP)
                )
    
            # ---- SAMPLING FOR THIS (ZAM, E)
            msg = f"E={e:3E} | covariance matrix size={rcov.data.shape}"
            log_stage(log, method, zam, msg, verbose=verbose)
            
            msg = f"E={e:3E} | sampling with SMP size={nsmp} via sampling(seed={local_seed}, {kwargs})"
            log_stage(log, method, zam, msg, verbose=verbose)
            smp = rcov.sampling(nsmp, seed=local_seed, **kwargs)

            # ---- FLATTEN into long-form DataFrame (pandas ≥ 2.1)
            smp_block = (
                smp.data.rename_axis(index="ZAP", columns="SMP")
                    .stack(future_stack=True)  # adopt new implementation
                    .rename("VALS")
                    .reset_index()  # -> columns: ["ZAP", "SMP", "VALS"]
                    .assign(E=e, ZAM=zam)[["ZAM", "E", "ZAP", "SMP", "VALS"]]
            )
            
            # ---- ENFORCE SMP integer type and stable ordering
            smp_block["SMP"] = smp_block["SMP"].astype(int)
            smp_block = smp_block.sort_values(["ZAM", "E", "ZAP", "SMP"], kind="mergesort")

            return smp_block


        # ---- LOOP OVER ALL FY BLOCKS
        smp_list = []
        for (zam, e), fy in nfpy_.data.query("MT==454").groupby(["ZAM", "E"]):
            mat = zam_mat_mapping[zam]
            intro_key = mat, 1, 451
            lib = Endf6({intro_key: self.data[intro_key]}).get_library()
            smp_list.append(_sample_fy_with_logging(zam, e, fy, lib, smp_kws_))
    
        smps = Samples(
            pd.concat(smp_list, ignore_index=True)
              .pivot_table(index=["ZAM", "E", "ZAP"], columns="SMP", values="VALS")
            )
 
        # ---- WRITE TO XLSX
        if write:
            xls_path = Path.cwd() / FILENAME_FY_PERT
            msg = f"writing samples -> '{xls_path}'"
            log_stage(log, method, zam, msg, verbose=verbose)

            if xls_path.exists():
                msg = "file exists and will be overwritten"
                log_stage(log, method, zam, msg, verbose=verbose)

            with pd.ExcelWriter(xls_path, engine="openpyxl") as writer:
                smps.data.to_excel(writer, index=True, sheet_name="SMP")
    
        # ---- CREATE OUTPUT MAPPING
        smps = {
            "IFY": smps,
        }

        msg = f"done | SMP size={nsmp} | # ZAM={length} | write={write}"
        log_stage(log, method, zam, msg, verbose=verbose)

        return smps


    def apply_perturbations(
            self,
            smps,
            *args,
            verbose=False,
            **kwargs,
            ):
        """
        Dispatcher to assign perturbations method: either for radioactive
        decay data, fission yields or cross section.
        
        Parameters
        ----------
        verbose : bool, optional, default is False
            It will log time messages and be bassed to the called method.

        Notes
        -----
        .. note :: The perturbation method is selected based on the MT's found
                   in `self`.

        Examples
        --------
        The next two examples will mismatch file and samples.
        The output must be `None` if samples and ENDF6 file do not match.

        >>> import sandy, pytest
        >>> tape_xs = sandy.get_endf6_file("jeff_33", "xs", 10010, local=True)
        >>> tape_d = sandy.get_endf6_file("jeff_33", "decay", 10040, local=True)
        >>> smps_xs = tape_xs.get_perturbations(2)
        >>> smps_d = tape_d.get_perturbations(2)

        Mix rdd file with xs samples. This will result in an error.

        >>> with pytest.raises(Exception):
        ...    not tape_d.apply_perturbations(smps_xs)

        Mix xs file with rdd samples. This will result in an error.

        >>> with pytest.raises(Exception):
        ...    not tape_xs.apply_perturbations(smps_d)

        An error is also raised if ``smps`` is not a mapping.

        >>> with pytest.raises(Exception):
        ...    not tape_xs.apply_perturbations(3)

        ...or if it empty.

        >>> with pytest.raises(Exception):
        ...    not tape_xs.apply_perturbations({})

        """
        # ---- IMPORT
        from collections.abc import Mapping

        from .utils import log
        from ._perturbation_base import log_stage

        # ---- SETUP
        msg = (
            "########################################################\n"
            "              APPLY PERTURBATIONS                       \n"
            "########################################################"
            )
        log(msg, verbose=verbose)

        t0 = time.perf_counter()

        # ---- RUNTIME TYPE & STRUCTURE CHECKS (early and explicit)
        if not isinstance(smps, Mapping):
            raise TypeError(
                f"`smps` must be a mapping (e.g., dict) got {type(smps).__name__}"
            )

        # ---- CHOOSE METHOD
        if 457 in self.mt:
            method = "apply_perturbations_rdd"
            out = self.apply_perturbations_rdd(smps, *args, verbose=verbose, **kwargs)

        elif 454 in self.mt:
            method = "apply_perturbations_fy"
            out = self.apply_perturbations_fy(smps, *args, verbose=verbose, **kwargs)

        else:
            method = "apply_perturbations_xs"
            out = self.apply_perturbations_xs(smps, *args, verbose=verbose, **kwargs)

        # ---- LOGGING: end
        dt = time.perf_counter() - t0
        msg = f"finished in {dt:.3f} s"
        log_stage(log, method, None, msg, verbose=verbose)

        return out

    @with_optional_warning_suppression("sandy.warn", default_suppress=True)
    def apply_perturbations_xs(
            self,
            smps: dict,
            *,
            ace_kws: dict | None = None,
            njoy_kws: dict | None = None,
            pendf=None,
            processes: int | str = 1,
            enable_tqdm: bool | None = None,
            suppress_njoy_output: bool = True,
            suppress_warnings: bool | None = None,
            to_ace: bool = False,
            to_file: bool = False,
            verbose: bool = False,
            **kwargs,
            ):
        """
        Apply relative perturbations to XS (MF=3), nubar (MT=452/455/456), and PFNS chi (MF=5/MT=18)
        for an ENDF6 evaluation, optionally in parallel.

        This method perturbs reaction cross sections and nubar values based on provided 
        perturbation samples. The process can be performed in parallel for efficiency. 
        If a PENDF file is not provided, it will be generated automatically.

        Parameters
        ----------
        smps : dict of :obj:`~sandy.samples.Samples`
            Mapping of MF/MT groups to Samples:
            - 31 → nubar perturbations (pnu)
            - 33 → cross-section perturbations (pxs)
            - 35 → chi perturbations (pchi)
        ace_kws : dict, optional
            Keyword arguments forwarded to `Endf6.get_ace` when `to_ace=True`.
        njoy_kws : dict, optional
            Keyword arguments forwarded to `Endf6.get_pendf` when `pendf` is not provided.
        pendf : :obj:`~sandy.endf6.Endf6`, optional
            If provided, perturbations are applied to this PENDF; otherwise a new PENDF is produced
            from `self` via `get_pendf(**njoy_kws)`.
        processes : int or "auto", optional (default=1)
            Number of worker processes. Use "auto" to pick `os.cpu_count()`.
            • If 1 → run in series (still uses the initializer to set caches).
            • If >1 → run in parallel using a spawn-safe ProcessPoolExecutor.
        enable_tqdm : bool or None, optional
            Control the use of ``tqdm`` progress bars.
            • ``None`` (default): progress bars follow ``verbose``. If
              ``verbose`` is True, ``tqdm`` is enabled; otherwise disabled.
            • ``True``: force ``tqdm`` progress bars on.
            • ``False``: force ``tqdm`` progress bars off.
        suppress_njoy_output : bool, optional, dafault is True
            Suppress NJOY output to screen (stdout or stderr) by redirecting it to DEVNULL.
            This method runs njoy many times (once per sample item). Then the
            default is `True` to avoid logging too much.
            Activate to see what njoy modules are running under the hood, and
            to have a feeling of how fast njoy modules run.
        suppress_warnings : bool, optional, deafult is None
            Standard warning associated to missing temperatures when using
            :obj:`~sandy.endf6.Endf6.get_pendf` are silenced.
            If not given, apply default (True).
        to_ace : bool, optional, default is False
            Flag to request perturbed ACE files.
        to_file : bool, optional, default is False
            Flag to ask workers to write files to disk.
        verbose : bool, optional, default is False
            print status. For printing info on PENDF and ACE creation, use
            `njoy_kws={'verbose': True}` and `ace_kws={'verbose': True}`.
            This key activates logging at INFO level and 

        Returns
        -------
        dict
            Mapping: `sample_id -> { "endf6": ..., "pendf": ..., ["ace": ..., "xsdir": ...] }`
            • If `to_ace=True`: mapping also contains ace/xsdir text as str.
            • If `to_file=True`: values are filenames.
            • If `to_file=False`: `"endf6"`/`"pendf"` are returned as in-memory dicts,
              and are wrapped into :obj:`~sandy.endf6.Endf6`.

        Performance & notes
        -------------------
        - On Windows/macOS (spawn), repeatedly shipping large ENDF6/PENDF dicts to workers is expensive.
          This method uses a per-process **initializer cache** (`init_xs_cache`) and a small **task wrapper**
          (`task_xs`) so that only the **per-sample payloads** (pxs/pnu/pchi for one sample ID) are sent.
        - Being a sort of wrapper, this method requires a high level of logging
          to report on its status. But because of this, we thought it was better to 
          suppress the njoy outputs to screen by default.
        - The function processes the sample IDs that are present in the requested
          perturbation kinds (31/33/35).
          A 1:1 pairing of (pxs, pnu, pchi) for the same `ismp` is requested, or an error is raised.
        - If you pass only one kind (e.g., 33), it will naturally process the sample IDs of that kind only.
        - Logging is formatted as  ``{method} | ZAM={ZAM} | {info}``


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

        First : produce perturbations for XS and nubar on U234 because it runs fast.

        >>> import sandy
        >>> tape = sandy.get_endf6_file("jeff_33", "xs", 922340, local=True)
        >>> sample_size = 2
        >>> njoy_kws={"err": 1, "chi": False, "mubar": False, "errorr33_kws": {"mt": [2, 4, 18]}}
        >>> smp_kws={"seed31": 1, "seed33": 3}
        >>> smps = tape.get_perturbations(sample_size, njoy_kws=njoy_kws, smp_kws=smp_kws, write=False)

        Apply both nubar and XS perturbations.

        >>> kws = dict(njoy_kws={"err": 1}, processes=1)
        >>> outs_31_33 = tape.apply_perturbations_xs(smps, **kws)

        Apply only nubar perturbations.

        >>> outs_31 = tape.apply_perturbations_xs({31: smps[31]}, **kws)

        Apply only XS perturbations.

        >>> outs_33 = tape.apply_perturbations_xs({33: smps[33]}, **kws)

        Check that outputs `'endf6'` and `'pendf'` exist and are the correct type.

        >>> assert isinstance(outs_33[0]["endf6"], sandy.Endf6)
        >>> assert isinstance(outs_33[0]["pendf"], sandy.Endf6)

        Check that files are different for different samples.

        >>> for i in range(2):
        ...    assert(outs_33[i]["endf6"].data == tape.data)   # endf6 did not change if nubar not perturbed
        ...    assert(outs_31[i]["endf6"].data != tape.data)   # endf6 changed if nubar is perturbed
        ...    assert(outs_31[i]["endf6"].data == outs_31_33[i]["endf6"].data)
        ...    assert(outs_33[i]["pendf"].data != outs_31[i]["pendf"].data)
        ...    assert(outs_33[i]["pendf"].data == outs_31_33[i]["pendf"].data)

        Cross-check across samples.

        >>> assert outs_33[0]["pendf"].data != outs_33[1]["pendf"].data
        >>> assert outs_33[0]["endf6"].data == outs_33[1]["endf6"].data
        >>> assert outs_31[0]["pendf"].data == outs_31[1]["pendf"].data
        >>> assert outs_31[0]["endf6"].data != outs_31[1]["endf6"].data

        Check that redundant nubar is also perturbed.

        >>> mat = 9225
        >>> nu0 = sandy.Xs.from_endf6(outs_31[0]["endf6"].filter_by(listmt=[452, 455, 456]))
        >>> nu1 = sandy.Xs.from_endf6(outs_31[1]["endf6"].filter_by(listmt=[452, 455, 456]))
        >>> assert not nu0.data[mat, 456].equals(nu1.data[mat, 456])   # perturbed
        >>> assert not nu0.data[mat, 452].equals(nu1.data[mat, 452])   # reconstructed
        >>> assert nu0.data[mat, 455].equals(nu1.data[mat, 455])       # no covariance

        Check that redundant and partial cross sections are correctly perturbed.

        >>> mat = 9225
        >>> xs0 = sandy.Xs.from_endf6(outs_33[0]["pendf"].filter_by(listmf=[3]))
        >>> xs1 = sandy.Xs.from_endf6(outs_33[1]["pendf"].filter_by(listmf=[3]))
        >>> for mt in [  2,   4,  18]:   # covariance present
        ...    assert not xs0.data[mat, mt].equals(xs1.data[mat, mt])
        >>> for mt in [  51,  52,  53,   # daughter reactions of mt4
        ...              54,  55,  56,  57,  58,  59,  60,  61,  62,  63,  64,  65,  66,  67,
        ...              68,  69,  70,  71,  72,  73,  74,  75,  76,  77,  78,  79,  80,  81,
        ...              82,  83,  84,  85,  86,  87,  88,  89,  90,  91]:
        ...    assert not xs0.data[mat, mt].equals(xs1.data[mat, mt])
        >>> for mt in [  51,  52,  53,   # daughter reactions of mt4
        ...              54,  55,  56,  57,  58,  59,  60,  61,  62,  63,  64,  65,  66,  67,
        ...              68,  69,  70,  71,  72,  73,  74,  75,  76,  77,  78,  79,  80,  81,
        ...              82,  83,  84,  85,  86,  87,  88,  89,  90,  91]:
        ...    assert not xs0.data[mat, mt].equals(xs1.data[mat, mt])
        >>> for mt in [  18,  19,  20,  21,  38]:   # daughter reactions of mt18
        ...    assert not xs0.data[mat, mt].equals(xs1.data[mat, mt])
        >>> for mt in [  16,  17,  37, 102]:   # not perturbed
        ...    assert xs0.data[mat, mt].equals(xs1.data[mat, mt])

        Second : produce perturbations for CHI and XS. Use Pu238 because it contains CHI.

        >>> tape = sandy.get_endf6_file("jeff_33", "xs", 942380, local=True)
        >>> sample_size = 2
        >>> njoy_kws = dict(err=1, nubar=False, mubar=False, errorr33_kws={"mt": [18]})
        >>> smp_kws = dict(seed33=3, seed35=5)  # preserve xs seed
        >>> smps_ = tape.get_perturbations(sample_size, njoy_kws=njoy_kws, smp_kws=smp_kws, write=False)

        Apply both chi and xs perturbations.

        >>> kws = dict(njoy_kws={"err": 1}, processes=1)
        >>> outs_33_35 = tape.apply_perturbations_xs(smps_, **kws)

        Compare to individual xs and chi perturbations with same seed.

        >>> outs_33_ = tape.apply_perturbations_xs({33: smps_[33]}, **kws)
        >>> outs_35 = tape.apply_perturbations_xs({35: smps_[35]}, **kws)

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

        Third : H1 case, check writing to file.

        >>> tape = sandy.get_endf6_file('jeff_33', 'xs', 10010, local=True)
        >>> sample_size = 2
        >>> njoy_kws = dict(err=1)
        >>> smps = tape.get_perturbations(sample_size, njoy_kws=njoy_kws, write=False)

        Check that ENDF6 and PENDF output filenames are correct.

        >>> kws = dict(njoy_kws={"err": 1}, processes=1)
        >>> outs = tape.apply_perturbations_xs(smps, to_file=True)
        >>> assert outs[0]["endf6"] == '1001_0.endf6' and os.path.isfile('1001_0.endf6')
        >>> assert outs[0]["pendf"] == '1001_0.pendf' and os.path.isfile('1001_0.endf6')
        >>> assert outs[1]["endf6"] == '1001_1.endf6' and os.path.isfile('1001_1.endf6')
        >>> assert outs[1]["pendf"] == '1001_1.pendf' and os.path.isfile('1001_1.pendf')

        Check that ACE output filenames are correct.

        >>> kws = dict(to_file=True, to_ace=True, ace_kws=dict(err=1, temperature=300, purr=False, heatr=False, thermr=False, gaspr=False))
        >>> outs = tape.apply_perturbations_xs(smps, **kws)
        >>> assert outs[0]["ace"] == '1001_0.03c' and os.path.isfile('1001_0.03c')
        >>> assert outs[0]["xsdir"] == '1001_0.03c.xsd' and os.path.isfile('1001_0.03c.xsd')
        >>> assert outs[1]["ace"] == '1001_1.03c' and os.path.isfile('1001_1.03c')
        >>> assert outs[1]["xsdir"] == '1001_1.03c.xsd' and os.path.isfile('1001_1.03c.xsd')

        Check that keyword `pendf` works.

        >>> pendf = tape.get_pendf(err=1)
        >>> outs1 = tape.apply_perturbations_xs(smps, njoy_kws=dict(err=1))
        >>> outs2 = tape.apply_perturbations_xs(smps, pendf=pendf)
        >>> assert outs1[0]["pendf"].write_string() == outs2[0]["pendf"].write_string()

        Fourth : Check parallelization vs serial path.

        >>> outs_par = tape.apply_perturbations_xs(smps, njoy_kws=dict(err=1), processes=2, enable_tqdm=False)
        >>> outs_ser = tape.apply_perturbations_xs(smps, njoy_kws=dict(err=1), processes=1)
        >>> for i in range(sample_size):
        ...    assert outs_ser[i]['endf6'].data == outs_par[i]['endf6'].data
        ...    assert outs_ser[i]['pendf'].data == outs_par[i]['pendf'].data

        Check that input kwargs `njoy_kws` and `ace_kws` do not mutate.
        Using `suppress_njoy_output` only mutates a copy of these dictionaries.

        >>> njoy_kws = {"err": 1}
        >>> ace_kws = {"err": 1, "purr": False}
        >>> outs = tape.apply_perturbations_xs(smps, njoy_kws=njoy_kws, ace_kws=ace_kws)
        >>> assert njoy_kws == {"err": 1}
        >>> assert ace_kws == {"err": 1, "purr": False}

        Check that an error is raised if no valid perturbation is present.

        >>> import pytest
        >>> with pytest.raises(Exception):
        ...    tape.apply_perturbations_xs({})

        >>> wrong_smps = {40: "aaa"}
        >>> with pytest.raises(Exception):
        ...    tape.apply_perturbations_xs(wrong_smps)

        """
        # ---- IMPORT
        import os, sys
        from concurrent.futures import ProcessPoolExecutor, as_completed
        from subprocess import DEVNULL

        from tqdm.auto import tqdm
        from tqdm.contrib.logging import logging_redirect_tqdm

        from ._concurrency import spawn_ctx, init_xs_cache, task_xs
        from ._perturbation_base import (
            validate_required_keys,
            validate_sample_ids,
            validate_smps_mapping,
            log_stage,
            )
        from .utils import log

        # ---- LOGGING SETUP
        zam = self.get_zam()
        method = "apply_perturbations_xs"
        # define whether to print warnings or not

        # ---- NORMALIZE dict-like keyword arguments
        njoy_kws_ = {} if njoy_kws is None else njoy_kws.copy()
        ace_kws_ = {} if ace_kws is None else ace_kws.copy()
        if suppress_njoy_output:
            njoy_kws_["njoy_output"] = ace_kws_["njoy_output"] = DEVNULL

        # ---- VALIDATE that at least one perturbation kind is present
        expected_keys = [31, 33, 35]
        validate_smps_mapping(smps)
        present = validate_required_keys(smps, expected_keys, mode="any")
        sample_ids = validate_sample_ids(smps, present)

        sample_size = len(sample_ids)

        msg = f"kinds={present} | SMP size={sample_size}"
        log_stage(log, method, zam, msg, verbose=verbose)

        # ---- PREPARE NOMINAL PENDF (if not provided)
        status = "provided" if pendf is not None else f"generated via get_pendf({njoy_kws_})"
        msg = f"PENDF={status}"
        log_stage(log, method, zam, msg, verbose=verbose)

        pendf_ = pendf if pendf is not None else self.get_pendf(**njoy_kws_)

        # Parallel execution + Windows spawn requires random access to per‑sample DataFrames
        # this cannot be done with a streaming iterator.
        # Need to materialize the iterator using dict()
        data_per_key = {mf: dict(smps[mf].iterate_xs_samples())
                        for mf in present}

        # This dict indexed by sample will contain the output of the worker:
        #    - either perturbed endf6 and pendf tape as `Endf6` objects
        #    - or ace files as string
        outs = {}

        # ---- PROGRESS BAR SETTINGS
        # Decide whether tqdm is enabled
        if enable_tqdm is None:
            # default: tqdm follows verbose
            tqdm_on = bool(verbose)
        else:
            # user override
            tqdm_on = bool(enable_tqdm)
        tqdm_kws = {
            "desc": "XS perturbations",
            "disable": not tqdm_on,
            "file": sys.stderr,
            "dynamic_ncols": True,
        }

        # ---- SERIAL sample production
        if processes in (None, 0, 1):

            msg = "mode=serial"
            log_stage(log, method, zam, msg, verbose=verbose)

            # Set per-process caches in THIS process
            init_xs_cache(self.data, pendf_.data)

            with logging_redirect_tqdm():  # ensure logs go via tqdm.write
                for ismp in tqdm(sample_ids, **tqdm_kws):
                    outs[ismp] = task_xs(
                        ismp,
                        pxs=data_per_key.get(33, {}).get(ismp),
                        pnu=data_per_key.get(31, {}).get(ismp),
                        pchi=data_per_key.get(35, {}).get(ismp),
                        ace_kws=ace_kws_,
                        to_ace=to_ace,
                        to_file=to_file,
                        verbose=verbose,
                    )

        # ---- PARALLEL sample production
        else:
            nprocs = (os.cpu_count() or 1) if processes in (
                "auto", None, 0) else int(processes)

            msg = f"mode=parallel | workers={nprocs:d}"
            log_stage(log, method, zam, msg, verbose=verbose)

            with ProcessPoolExecutor(
                max_workers=nprocs,
                mp_context=spawn_ctx(),            # Windows/macOS spawn-safe
                initializer=init_xs_cache,         # cache nominal dicts once per worker
                initargs=(self.data, pendf_.data),
            ) as ex:
                futures = {
                    ex.submit(
                        task_xs, ismp,
                        pxs=data_per_key.get(33, {}).get(ismp),
                        pnu=data_per_key.get(31, {}).get(ismp),
                        pchi=data_per_key.get(35, {}).get(ismp),
                        ace_kws=ace_kws_,
                        to_ace=to_ace,
                        to_file=to_file,
                        verbose=verbose,
                    ): ismp for ismp in sample_ids
                }

                msg = f"submitting {sample_size} tasks"
                log_stage(log, method, zam, msg, verbose=verbose)

                with logging_redirect_tqdm():  # ensure logs go via tqdm.write
                    for fut in tqdm(as_completed(futures), total=len(futures), **tqdm_kws):
                        ismp = futures[fut]
                        outs[ismp] = fut.result()

        msg = "collected all results"
        log_stage(log, method, zam, msg, verbose=verbose)

        # ---- WRAPPING UP in-memory dicts into Endf6 instances (cannot be done inside workers)
        if not to_file:

            msg = "wrapping worker outputs into Endf6 objects"
            log_stage(log, method, zam, msg, verbose=verbose)

            for ismp, out_dict in outs.items():
                for key, value in out_dict.items():
                    if key in ("endf6", "pendf"):
                        out_dict[key] = Endf6(value)

        # ---- ENSURE DETERMINISTIC ORDERING OF OUTPUT KEYS
        outs = dict(sorted(outs.items()))

        msg = "done"
        log_stage(log, method, zam, msg, verbose=verbose)

        return outs

    @with_optional_warning_suppression("sandy.warn", default_suppress=True)
    def apply_perturbations_rdd(
            self,
            smps: dict,
            *,
            processes: int | str = 1,
            rdd = None,
            enable_tqdm: bool | None = None,
            suppress_warnings: bool | None = None,
            to_file: bool = False,
            verbose: bool = False,
            **kwargs,
            ):
        """
        Apply sampled perturbations to radioactive-decay data (RDD) contained in an
        :class:`~sandy.endf6.Endf6` object and generate perturbed ENDF-6 files.
    
        This method takes the RDD perturbation samples produced by
        :meth:`~sandy.endf6.Endf6.get_perturbations_rdd` and applies them to the
        nominal decay data (MF=8/MT=457). For each sample ID, a perturbed ENDF-6
        tape is created. Depending on ``to_file``, results are either returned as
        in-memory :class:`~sandy.endf6.Endf6` objects or written directly to disk.
    
        Parameters
        ----------
        smps : dict
            Mapping from perturbation type to :class:`sandy.samples.Samples`
            instances. Must contain the keys ``"HL"``, ``"DE"``, and ``"BR"``,
            representing sampled perturbations for half-lives, decay energies, and
            branching ratios, respectively. The Samples objects must share a
            consistent index layout.
    
        processes : int or {"auto"}, optional
            Number of worker processes:
            - ``1`` (default): run in serial mode.
            - ``>1``: parallel execution using ``ProcessPoolExecutor``.
            - ``"auto"``: automatically use all available CPU cores.
    
        rdd : sandy.DecayData, optional
            Precomputed :class:`sandy.decay.DecayData` object. If not provided,
            it is extracted from ``self`` via
            :meth:`sandy.decay.DecayData.from_endf6`.
    
        enable_tqdm : bool or None, optional
            Control the display of ``tqdm`` progress bars.
    
            - ``None`` (default): follow ``verbose``  
              (progress bars shown when ``verbose=True``).
            - ``True``: always show progress bars.
            - ``False``: always disable progress bars.
    
            This option provides fine-grained control over progress-display
            behavior, preventing clutter in non-interactive CLI environments
            (e.g., when running via ``python -m sandy.sampling``), while still
            enabling helpful progress visualization in interactive Python sessions.
    
        suppress_warnings : bool or None, optional
            Whether to suppress warnings emitted during the calculation.
            The default behavior is controlled by the decorator
            :func:`with_optional_warning_suppression`.
    
        to_file : bool, optional
            If ``True``, each perturbed ENDF-6 tape is written to a file named
            ``decay_data_<sampleID>`` in the current working directory.
            If ``False`` (default), perturbed tapes are returned as Endf6 objects.
    
        verbose : bool, optional
            If ``True``, enable detailed progress messages and diagnostics.
    
        **kwargs :
            Additional keyword arguments reserved for future extensions. They are
            currently ignored.
    
        Returns
        -------
        outs : dict
            Dictionary mapping sample IDs to results:
            - If ``to_file=False``: ``{smpID: Endf6}``
            - If ``to_file=True``:  ``{smpID: filepath}``
    
            Output entries are sorted by sample ID.
    
        Notes
        -----
        - Perturbations are multiplicative factors applied to half-lives,
          decay constants (recomputed from half-lives), and decay energies,
          consistent with MF=8/MT=457.
        - Branching ratios are renormalized to unity after perturbation.
        - Parallel and serial execution produce numerically identical results.
        - When ``to_file=True``, output names follow the pattern
          ``decay_data_<sampleID>``.


        Examples
        --------
        Basic usage and testing.

        Produce few perturbations (only 2) to check consistency.

        >>> import sandy, numpy as np
        >>> tape = sandy.get_endf6_file("jeff_33", "decay", [10040, 270590, 270600, 571380], local=True)
        >>> rdd = sandy.DecayData.from_endf6(tape)
        >>> sample_size = 2
        >>> smps = tape.get_perturbations(sample_size, rdd=rdd)
        >>> outs = tape.apply_perturbations_rdd(smps, rdd=rdd)

        Let's extract the decay data from the first perturbed file.

        >>> idx = 0
        >>> rdd0 = sandy.DecayData.from_endf6(outs[idx])



        This first suite of tests checks that the perturbation is propagated correctly 
        down to the new perturbed ENDF-6 file.

        Check that half-lives are correctly perturbed for all unstable nuclides.

        >>> decimal = 5
        >>> for nuclide in (10040, 270600, 571380):
        ...    expected_relpert = smps["HL"].data.loc[nuclide, idx]
        ...    from_file_relpert = rdd0.data[nuclide]['half_life'] / rdd.data[nuclide]['half_life']
        ...    np.testing.assert_almost_equal(from_file_relpert, expected_relpert, decimal=decimal)

        For the stable nuclide (Co59) the halflife remains 0.

        >>> assert rdd0.data[270590]['half_life'] == rdd.data[270590]['half_life'] == 0

        The same happens for the decay constants (they are recalculated).
        The relative perturbation is the same for decay constants and half lives.

        >>> decimal = 5
        >>> for nuclide in (10040, 270600, 571380):
        ...    expected_relpert = smps["HL"].data.loc[nuclide, idx]
        ...    from_file_relpert = rdd.data[nuclide]['decay_constant'] / rdd0.data[nuclide]['decay_constant']

        For the stable nuclide (Co59) the decay constant remains 0.

        >>> assert rdd0.data[270590]['decay_constant'] == rdd.data[270590]['decay_constant'] == 0

        Check that decay energies are also correctly perturbed.
        We only look at Co60.

        >>> nuclide = 270600
        >>> decimal = 5
        >>> for energy in ("beta", "gamma"):
        ...    expected_relpert = smps["DE"].data.loc[(nuclide, energy), 0]
        ...    from_file_relpert = rdd0.data[nuclide]['decay_energy'][energy] / rdd.data[nuclide]['decay_energy'][energy]
        ...    np.testing.assert_almost_equal(expected_relpert, from_file_relpert, decimal=decimal)

        Co60 does not hava alpha decay energy, so it should remain zero.

        >>> assert rdd0.data[nuclide]['decay_energy']["alpha"] == rdd.data[nuclide]['decay_energy']["alpha"] == 0

        Alpha energy perturbation is checked in H4.

        >>> nuclide = 10040
        >>> decimal = 5
        >>> expected_relpert = smps["DE"].data.loc[(nuclide, "alpha"), 0]
        >>> from_file_relpert = rdd0.data[nuclide]['decay_energy']["alpha"] / rdd.data[nuclide]['decay_energy']["alpha"]
        >>> np.testing.assert_almost_equal(expected_relpert, from_file_relpert, decimal=decimal)



        In this test suite we validate the keywords.

        Check that perturbed files are written to disk with key `to_file`.

        >>> # Clean up perturbed files
        >>> from pathlib import Path
        >>> outdir = Path.cwd()
        >>> p = outdir / "decay_data_0"
        >>> if p.exists(): p.unlink()
        >>> assert not p.exists()

        Now run again and write.

        >>> outs = tape.apply_perturbations_rdd(smps, processes=1, to_file=True, rdd=rdd)
        >>> assert p.exists()
        
        Now let's read it up, so we check that the file is written correctly.
        
        >>> tape0_serial = sandy.Endf6.from_file(p)
        


        In this test suite we check that identical results are produced with
        serial and parallel mode, also with and without ``rdd``.
        
        >>> outs_parallel = tape.apply_perturbations_rdd(smps, processes=2, enable_tqdm=False)
        >>> tape0_parallel = outs_parallel[0]
        
        >>> for k in tape0_serial.data:
        ...    assert tape0_serial.data[k] == tape0_parallel.data[k]
        


        In this test suite we ensure convergence of the sample statistics.

        The test is done for a stable nuclide (Co59) and a non-stable nuclide (Co60).
        The sample size should be high enough to guarantee convergence within the selected tolerances.

        >>> decay = sandy.get_endf6_file("jeff_33", "decay", [270590, 270600], local=True)
        >>> sample_size = 100
        >>> smps = decay.get_perturbations_rdd(sample_size, write=False)
        >>> outs = decay.apply_perturbations_rdd(smps, verbose=False)
        
        This is the statistical summary of the sample.

        >>> import pandas as pd, numpy as np
        >>> hl = pd.DataFrame({ismp: sandy.DecayData.from_endf6(outs[ismp]).get_half_life().data["HL"] for ismp in range(sample_size)})
        >>> descr = hl.T.describe().T
        >>> smp_mean, smp_std = descr["mean"], descr["std"]

        The sample reproduces the original data within a given tolerance.
        
        >>> expected = sandy.DecayData.from_endf6(decay).get_half_life().data
        >>> assert np.allclose(smp_mean, expected["HL"], rtol=1e-4)
        >>> assert np.allclose(smp_std, expected["DHL"], rtol=0.2)

        The convergence for branching ratios is tested on Cs134, for which
        branching ratio uncertainties are given.

        >>> import math
        >>> decay = sandy.get_endf6_file("jeff_33", "decay", 551340, local=True)
        >>> sample_size = 100
        >>> smps = decay.get_perturbations_rdd(sample_size, write=False)
        >>> outs = decay.apply_perturbations_rdd(smps, processes=1)
        
        Here the perturbed branching ratios are colleceted from the output files.

        >>> br = pd.DataFrame({ismp: sandy.DecayData.from_endf6(outs[ismp]).get_branching_ratio().data["BR"] for ismp in range(sample_size)})

        The normalization is verified.

        >>> assert np.allclose(br.sum(), 1, atol=1e-6)
        
        Because of the normalization, a large anti-correlation is introduced
        across the branching ratios.

        >>> smp_corr = np.corrcoef(br)[0, 1]
        >>> assert smp_corr < -0.99

        The statistical mean and standard deviation are collected and compared
        to expected values present in the original tape.

        >>> descr = br.T.describe().T
        >>> smp_mean, smp_std = descr["mean"], descr["std"]
        >>> expected = sandy.DecayData.from_endf6(decay).get_branching_ratio().data

        The convergence is verified.

        >>> for k, v in expected.T.items():
        ...     stat = smp_mean.loc[k]
        ...     expected_stat = v["BR"]
        ...     assert math.isclose(stat, expected_stat, rel_tol=1e-5, abs_tol=2e-7)
        ...     stat = smp_std.loc[k]
        ...     expected_stat = v["DBR"]
        ...     # the convergence for small numbers is only tested with absolute tolerance
        ...     assert math.isclose(stat, expected_stat, abs_tol=2e-7)

        If there is only one decay mode, the branching ratio is `1`, and its
        uncertainty is `0`. Even if this is modified on the perturbation stage,
        the variability is removed by the normalization
        
        >>> decay = sandy.get_endf6_file("jeff_33", "decay", 10050, local=True)
        >>> sample_size = 10
        >>> smps = decay.get_perturbations_rdd(sample_size, write=False, fill_zeros_branching_ratio=0.5)
        >>> outs = decay.apply_perturbations_rdd(smps, processes=1)
        >>> for out in outs.values():
        ...     assert np.all(sandy.DecayData.from_endf6(out).get_branching_ratio().data.BR == 1)



        This test suite checks errors.

        Check that an error is raised if no valid perturbation is present.

        >>> import pytest
        >>> with pytest.raises(Exception):
        ...    tape.apply_perturbations_rdd({})

        >>> wrong_smps = {"XS": "aaa"}
        >>> with pytest.raises(Exception):
        ...    tape.apply_perturbations_rdd(wrong_smps)

        """
        # ---- IMPORT
        import sys
        import numpy as np
        from concurrent.futures import ProcessPoolExecutor, as_completed
        from tqdm.auto import tqdm
        from tqdm.contrib.logging import logging_redirect_tqdm

        from .decay import DecayData
        from ._concurrency import spawn_ctx, init_rdd_cache, task_rdd
        from .utils import log
        from ._perturbation_base import (
            validate_required_keys,
            validate_sample_ids,
            validate_smps_mapping,
            log_stage,
            )

        # ---- SETUP
        # there is likely no single zam
        zam = self.get_zam()
        method = "apply_perturbations_rdd"
        length = 1 if np.isscalar(zam) else len(zam)
        msg = f"found {length} ZAM"
        log_stage(log, method, zam, msg, verbose=verbose)

        # ---- VALIDATE that at least one perturbation kind is present
        required_keys = ["HL", "DE", "BR"]
        validate_smps_mapping(smps)
        present = validate_required_keys(smps, required_keys, mode="all")
        sample_ids = validate_sample_ids(smps, present)

        sample_size = len(sample_ids)
        msg = f"SMP size={sample_size}"
        log_stage(log, method, zam, msg, verbose=verbose)
    

        # ---- PREPARE DecayData object (if not provided)
        # if already available in kwargs, do not extract DecayData again
        status = "provided" if rdd is not None else "extracted with DecayData.from_endf6"
        msg = (f"DecayData={status}")
        log_stage(log, method, zam, msg, verbose=verbose)
        # pass verbosity
        rdd_ = rdd if rdd is not None else DecayData.from_endf6(self, verbose=verbose)


        # This dict will contain outputs per sample id (Endf6 dict or filename)
        outs = {}
    
        # ---- PROGRESS BAR SETTINGS
        # Decide whether tqdm is enabled
        if enable_tqdm is None:
            # default: tqdm follows verbose
            tqdm_on = bool(verbose)
        else:
            # user override
            tqdm_on = bool(enable_tqdm)
        tqdm_kws = {
            "desc": "XS perturbations",
            "disable": not tqdm_on,
            "file": sys.stderr,
            "dynamic_ncols": True,
        }

        # ---- SERIAL EXECUTION
        if processes in (None, 0, 1):
            msg = "mode=serial"
            log_stage(log, method, zam, msg, verbose=verbose)
    
            # Initialize per-process cache in THIS process
            init_rdd_cache(self.data, rdd_.data)
    
            with logging_redirect_tqdm():
                for ismp in tqdm(sample_ids, **tqdm_kws):
                    # pass only the single-column frames to minimize payload
                    # pass as dataframes to mimic xs
                    hl_col = smps["HL"].data[ismp].rename("HL").to_frame()
                    de_col = smps["DE"].data[ismp].rename("E").to_frame()  # inconsistency, "DE" in smps, but "E" in DecayData
                    br_col = smps["BR"].data[ismp].rename("BR").to_frame()
    
                    outs[ismp] = task_rdd(
                        ismp,
                        phl=hl_col,
                        pde=de_col,
                        pbr=br_col,
                        to_file=to_file,
                        verbose=verbose,
                    )
    
        # ---- PARALLEL EXECUTION
        else:
            nprocs = (os.cpu_count() or 1) if processes in ("auto", None, 0) else int(processes)
            msg = f"mode=parallel | workers={nprocs:d}"
            log_stage(log, method, zam, msg, verbose=verbose)
    
            with ProcessPoolExecutor(
                max_workers=nprocs,
                mp_context=spawn_ctx(),           # Windows/macOS spawn-safe
                initializer=init_rdd_cache,       # cache nominal dicts once per worker
                initargs=(self.data, rdd_.data),
            ) as ex:
                futures = {}
                for ismp in sample_ids:
                    # pass only the single-column frames to minimize payload
                    # pass as dataframes to mimic xs
                    hl_col = smps["HL"].data[ismp].rename("HL").to_frame()
                    de_col = smps["DE"].data[ismp].rename("E").to_frame()  # inconsistency, "DE" in smps, but "E" in DecayData
                    br_col = smps["BR"].data[ismp].rename("BR").to_frame()
    
                    fut = ex.submit(
                        task_rdd,
                        ismp,
                        phl=hl_col,
                        pde=de_col,
                        pbr=br_col,
                        to_file=to_file,
                        verbose=verbose,
                    )
                    futures[fut] = ismp
    
                msg = f"submitting {sample_size} tasks"
                log_stage(log, method, zam, msg, verbose=verbose)
    
                with logging_redirect_tqdm():
                    for fut in tqdm(as_completed(futures), total=len(futures), **tqdm_kws):
                        ismp = futures[fut]
                        outs[ismp] = fut.result()

        msg = "collected all results"
        log_stage(log, method, zam, msg, verbose=verbose)


        # ---- WRAP in-memory dicts into Endf6 objects (cannot be done inside workers)
        if not to_file:
            msg = "wrapping worker outputs into Endf6 objects"
            log_stage(log, method, zam, msg, verbose=verbose)
            for ismp, value in outs.items():
                # for key, value in outs[ismp].items():   # only 1 key, that is "endf6"
                outs[ismp] = Endf6(value)

        # ---- DETERMINISTIC ORDER
        outs = dict(sorted(outs.items()))
    
        msg = "done"
        log_stage(log, method, zam, msg, verbose=verbose)

        return outs

    @with_optional_warning_suppression("sandy.warn", default_suppress=True)
    def apply_perturbations_fy(
            self,
            smps,
            *,
            processes: int | str = 1,
            nfpy=None,
            enable_tqdm: bool | None = None,
            suppress_warnings: bool | None = None,
            to_file: bool = False,
            verbose: bool = False,
            **kwargs,
            ):
        """
        Apply sampled perturbations to the independent fission yields (IFYs) in an
        :class:`~sandy.endf6.Endf6` object and generate perturbed ENDF-6 files.
    
        This method takes the perturbation samples produced by
        :meth:`~sandy.endf6.Endf6.get_perturbations_fy` and applies them to the
        nominal FY data (MF=8/MT=454). For each sample, a perturbed ENDF-6 tape
        is created. Depending on ``to_file``, the results are either returned as
        new :class:`~sandy.endf6.Endf6` objects or written to disk.
    
        Parameters
        ----------
        smps : dict
            Dictionary produced by :meth:`get_perturbations_fy`. Must contain
            exactly one key ``"IFY"`` mapped to a
            :class:`sandy.samples.Samples` instance. The Samples object must have
            a multi-index ``(ZAM, E, ZAP)`` and columns containing sample IDs.
    
        processes : int or {"auto"}, optional
            Number of worker processes to use.
            - ``1`` (default): run in serial mode
            - ``>1``: parallel execution using ``ProcessPoolExecutor``
            - ``"auto"``: automatically use all available CPU cores
    
        nfpy : sandy.Fy, optional
            Precomputed FY object. If omitted, FY data is extracted from ``self``
            using :class:`sandy.fy.Fy`.
    
        enable_tqdm : bool or None, optional
            Control the display of ``tqdm`` progress bars.
    
            - ``None`` (default): follow the value of ``verbose``  
              (i.e. progress bars are enabled only when ``verbose=True``).
            - ``True``: always show progress bars, regardless of ``verbose``.
            - ``False``: disable progress bars entirely.
    
            This option allows fine‑grained control of progress display, avoiding
            broken or noisy progress bars in non‑interactive environments (e.g.
            when running via ``python -m sandy.sampling``), while still enabling
            useful progress visualization in interactive Python sessions.
    
        suppress_warnings : bool or None, optional
            Control whether warnings emitted inside the method are suppressed.
            If ``None`` (default), suppression behavior follows the decorator
            :func:`with_optional_warning_suppression`.
    
        to_file : bool, optional
            If ``True``, each perturbed ENDF‑6 tape is written to a file named
            ``fy_<sampleID>`` in the current working directory.  
            If ``False`` (default), perturbed tapes are returned as Endf6 objects.
    
        verbose : bool, optional
            Enable verbose diagnostic logging and status messages.
    
        Returns
        -------
        outs : dict
            A mapping from sample ID to result:
    
            - if ``to_file=False``: ``{smpID: Endf6}``
            - if ``to_file=True``:  ``{smpID: filepath}``
    
            Output entries are sorted by sample ID.
    
        Notes
        -----
        - Perturbations are multiplicative factors applied directly to FY values
          in MF=8/MT=454.
        - Parallel and serial modes produce identical numerical results.
        - When ``to_file=True``, filenames follow the template ``fy_<sampleID>``.
        - Sample IDs are taken from the column names of the Samples object.


        Examples
        --------
        Basic usage and testing.

        Produce few perturbations (only 2) to check consistency.

        >>> import sandy, numpy as np
        >>> tape = sandy.get_endf6_file("jeff_33", "nfpy", 922350, local=True)
        >>> fy = sandy.Fy.from_endf6(tape)
        >>> sample_size = 2
        >>> smps = tape.get_perturbations(sample_size, nfpy=fy, write=False)
        >>> outs = tape.apply_perturbations_fy(smps, nfpy=fy)

        Let's extract the fission yields from the first perturbed file.

        >>> idx = 0
        >>> fy0 = sandy.Fy.from_endf6(outs[idx])



        This first suite of tests checks that the perturbation is propagated correctly 
        down to the new perturbed ENDF-6 file.

        Check that IFYs are correctly perturbed for all energies.

        >>> nuclide, mt = 922350, 454
        >>> for e in fy.data.E.unique():
        ...    expected_relpert = smps["IFY"].data.query("E==@e").droplevel(["ZAM", "E"])[idx]
        ...    block0 = fy0.data.query("MT==@mt and E==@e").set_index("ZAP")["FY"]
        ...    block = fy.data.query("MT==@mt and E==@e").set_index("ZAP")["FY"]
        ...    from_file_relpert = block0.div(block).fillna(1)
        ...    np.testing.assert_array_almost_equal(from_file_relpert, expected_relpert, decimal=4)        
        
        
        
        In this test suite we validate the keywords.

        Check that perturbed files are written to disk with key `to_file`.

        >>> # Clean up perturbed files
        >>> from pathlib import Path
        >>> outdir = Path.cwd()
        >>> p = outdir / "fy_0"
        >>> if p.exists(): p.unlink()
        >>> assert not p.exists()

        Now run again and write.

        >>> outs = tape.apply_perturbations_fy(smps, processes=1, to_file=True, nfpy=fy)
        >>> assert p.exists()
        
        Now let's read it up, so we check that the file is written correctly.
        
        >>> tape0_serial = sandy.Endf6.from_file(p)
        


        In this test suite we check that identical results are produced with
        serial and parallel mode, also with and without ``nfpy``.
        
        >>> outs_parallel = tape.apply_perturbations_fy(smps, processes=2, enable_tqdm=False)
        >>> tape0_parallel = outs_parallel[0]
        
        >>> for k in tape0_serial.data:
        ...    assert tape0_serial.data[k] == tape0_parallel.data[k]



        In this test suite we ensure convergence of the sample statistics.

        The sample size should be high enough to guarantee convergence within the selected tolerances.

        >>> import sandy, pandas as pd, numpy as np
        >>> tape = sandy.get_endf6_file("jeff_33", "nfpy", 922350, local=True)
        >>> fy = sandy.Fy.from_endf6(tape)
        >>> sample_size = 100
        >>> smps = tape.get_perturbations(sample_size, nfpy=fy, write=False)
        >>> outs = tape.apply_perturbations_fy(smps, nfpy=fy)

        To process the outputs we build a dataframe with all samples' FY values, indexed by (E, ZAP).

        >>> dict_df = {}
        >>> for ismp in outs:
        ...     fy_ismp = sandy.Fy.from_endf6(outs[ismp])
        ...     dict_df[ismp] = fy_ismp.data.query("MT==454")[["E","ZAP","FY"]].set_index(["E","ZAP"])["FY"]
        >>> fy_stack = sandy.Samples(dict_df)

        We compare the std from the files with that from the perturbation coefficients.
        It should have been passed through without any major problem
        (`decimal=5` is a rather accurate check, because of 1-to-1 equivalence of FYs).

        >>> expected = smps["IFY"].get_rstd().droplevel("ZAM")
        >>> got = fy_stack.get_rstd().fillna(0)
        >>> np.testing.assert_array_almost_equal(got, expected, decimal=5)

        Same for the mean value.

        >>> fy_nominal = fy.data.query("MT==454")[["E","ZAP","FY"]].set_index(["E","ZAP"])["FY"]
        >>> expected = smps["IFY"].get_mean().droplevel("ZAM")
        >>> got = (fy_stack.get_mean() / fy_nominal).fillna(1)
        >>> np.testing.assert_array_almost_equal(got, expected, decimal=5)



        This test suite checks errors.

        Check that an error is raised if no valid perturbation is present.

        >>> import pytest
        >>> with pytest.raises(Exception):
        ...    tape.apply_perturbations_fy({})

        >>> wrong_smps = {"XS": "aaa"}
        >>> with pytest.raises(Exception):
        ...    tape.apply_perturbations_fy(wrong_smps)

        """
        # ---- IMPORTS
        import os, sys
        import numpy as np
        from concurrent.futures import ProcessPoolExecutor, as_completed
        from tqdm.auto import tqdm
        from tqdm.contrib.logging import logging_redirect_tqdm
    
        from .fy import Fy
        from ._concurrency import spawn_ctx, init_fy_cache, task_fy
        from .utils import log
        from ._perturbation_base import (
            validate_required_keys,
            validate_sample_ids,
            validate_smps_mapping,
            log_stage,
            )


        # ---- SETUP
        zam = self.get_zam()
        method = "apply_perturbations_fy"
        length = 1 if np.isscalar(zam) else len(zam)
        msg = f"found {length} ZAM"
        log_stage(log, method, zam, msg, verbose=verbose)


        # ---- VALIDATE that at least one perturbation kind is present
        required_keys = ["IFY"]
        validate_smps_mapping(smps)
        present = validate_required_keys(smps, required_keys, mode="all")
        sample_ids = validate_sample_ids(smps, present)

        # there is only one data type
        smp = smps["IFY"]
        index_names = set(smp.data.index.names)
        required_idx = {"ZAM", "E", "ZAP"}
        if not required_idx == index_names:
            raise Exception(
                f"'smp' must contain index {required_idx}, got {index_names}"
            )

        sample_ids = sorted(smp.data.columns.unique())
        sample_size = len(sample_ids)
        msg = f"SMP size={sample_size}"
        log_stage(log, method, zam, msg, verbose=verbose)
    

        # ---- PREPARE FY OBJECT
        status = "provided" if nfpy is not None else "extracted with Fy.from_endf6"
        msg = f"Fy={status}"
        log_stage(log, method, zam, msg, verbose=verbose)

        nfpy_ = nfpy if nfpy is not None else Fy.from_endf6(self, verbose=verbose)


        # This dict will contain outputs per sample id (Endf6 dict or filename)
        outs = {}
        
        # ---- PROGRESS BAR SETTINGS
        # Decide whether tqdm is enabled
        if enable_tqdm is None:
            # default: tqdm follows verbose
            tqdm_on = bool(verbose)
        else:
            # user override
            tqdm_on = bool(enable_tqdm)
        tqdm_kws = {
            "desc": "XS perturbations",
            "disable": not tqdm_on,
            "file": sys.stderr,
            "dynamic_ncols": True,
        }
    
        # ---- SERIAL EXECUTION
        if processes in (None, 0, 1):
            msg = "mode=serial"
            log_stage(log, method, zam, msg, verbose=verbose)
    
            # Initialize per-process cache in THIS process
            init_fy_cache(self.data, nfpy_.data)
    
            with logging_redirect_tqdm():
                for ismp in tqdm(sample_ids, **tqdm_kws):
                    # pass only the single-column frames to minimize payload
                    # passed as dataframe to mimic xs
                    fy_col = smp.data[ismp].rename("IFY").to_frame()
    
                    outs[ismp] = task_fy(
                        ismp,
                        pfy=fy_col,
                        to_file=to_file,
                        verbose=verbose,
                    )
        # ---- PARALLEL EXECUTION
        else:
            nprocs = (os.cpu_count() or 1) if processes in ("auto", None, 0) else int(processes)
            msg = f"mode=parallel | workers={nprocs:d}"
            log_stage(log, method, zam, msg, verbose=verbose)
    
            with ProcessPoolExecutor(
                max_workers=nprocs,
                mp_context=spawn_ctx(),           # Windows/macOS spawn-safe
                initializer=init_fy_cache,        # cache nominal dicts once per worker
                initargs=(self.data, nfpy_.data),
            ) as ex:
                futures = {}
                for ismp in sample_ids:
                    # pass only the single-column frames to minimize payload
                    # passed as dataframe to mimic xs
                    fy_col = smp.data[ismp].rename("IFY").to_frame()
    
                    fut = ex.submit(
                        task_fy,
                        ismp,
                        pfy=fy_col,
                        to_file=to_file,
                        verbose=verbose,
                    )
                    futures[fut] = ismp
    
                msg = f"submitting {sample_size} tasks"
                log_stage(log, method, zam, msg, verbose=verbose)
    
                with logging_redirect_tqdm():
                    for fut in tqdm(as_completed(futures), total=len(futures), **tqdm_kws):
                        ismp = futures[fut]
                        outs[ismp] = fut.result()

        # ---- COLLECTED RESULTS
        msg = "collected all results"
        log_stage(log, method, zam, msg, verbose=verbose)

        # ---- WRAP INTO Endf6 OBJECTS ----
        if not to_file:
            msg = "wrapping outputs into Endf6 objects"
            log_stage(log, method, zam, msg, verbose=verbose)
    
            for ismp, data in outs.items():
                outs[ismp] = Endf6(data)
    
        # ---- ORDER RESULTS ----
        outs = dict(sorted(outs.items()))
    
        msg = "done"
        log_stage(log, method, zam, msg, verbose=verbose)
    
        return outs
