
import os
from os.path import join
import shutil
import re
import logging
import pdb
from tempfile import TemporaryDirectory
import subprocess as sp

import pandas as pd
import numpy as np

import sandy

__author__ = "Luca Fiorito"
__all__ = [
    "process_neutron",
    "process_proton",
    "get_njoy",
    ]

sab = pd.DataFrame.from_records([[48,9237,1,1,241,'uuo2'],
                                  [42,125,0,8,221,'tol'],
                                  [59,1425,0,1,221,'si'],
                                  [37,125,11,2,223,'pol'],
                                  [2,125,0,2,221,'ph'],
                                  [12,128,0,2,221,'pd'],
                                  [75,825,1,1,239,'ouo2'],
                                  [48,825,0,3,221,'osap'],
                                  [51,825,0,1,222,'ohw'],
                                  [46,825,3,1,237,'obeo'],
                                  [3,125,0,2,221,'oh2'],
                                  [13,128,0,2,221,'od2'],
                                  [52,1225,0,1,249,'mg'],
                                  [38,125,0,12,221,'mesi'],
                                  [10,125,1,2,221,'ice'],
                                  [7,125,12,1,225,'hzr'],
                                  [1,125,0,2,222,'lw'],
                                  [8,125,0,2,249,'hca'],
                                  [31,600,1,1,229,'gr'],
                                  [11,128,0,2,221,'dhw'],
                                  [59,2025,0,1,249,'cah'],
                                  [27,425,3,1,233,'bbeo'],
                                  [26,425,2,1,231,'be'],
                                  [60,1325,0,2,221,'asap']],
            columns = ['matde','matdp','icoh','natom','mtref','ext'])

tmp2ext = {
    293.6: "02",
    300: "03",
    350: "35",
    400: "04",
    450: "45",
    500: "05",
    550: "55",
    600: "06",
    650: "65",
    700: "07",
    750: "75",
    800: "08",
    850: "85",
    900: "09",
    950: "95",
    1000: "10",
    1100: "11",
    1200: "12",
    1300: "13",
    1400: "14",
    1500: "15",
    1800: "18",
    2100: "21",
    2400: "24",
    }

tmp2ext_meta = {
    293.6: "30",
    300: "31",
    350: "32",
    400: "33",
    450: "34",
    500: "36",
    550: "37",
    600: "38",
    650: "39",
    700: "40",
    750: "41",
    800: "42",
    850: "43",
    900: "44",
    950: "46",
    1000: "47",
    1100: "48",
    1200: "49",
    1300: "50",
    1400: "51",
    1500: "52",
    1800: "53",
    2100: "54",
    2400: "55",
    }


NJOY_TOLER = 0.001
NJOY_TEMPERATURES = [293.6]
NJOY_SIG0 = [1e10]
NJOY_THERMR_EMAX = 10

# input taken from
# https://www-nds.iaea.org/index-meeting-crp/TM_NDP/docs/OCabellos_2017.pdf
_input_mf32_nomf33 = """errorr
999/
20 33/
1/
2/
18/
102/
0/
stop"""
# same but no fissile
_input_mf32_nomf33_no18 = """errorr
999/
20 33/
1/
2/
102/
0/
stop"""


def get_njoy():
    """
    Extract njoy executable from system environment variable `NJOY`.

    Returns
    -------
    `string`
        njoy executable
    """
    return os.environ["NJOY"]


def get_temperature_suffix(temperature, meta=False):
    """
    Determine suffix saccording to temperature value.
    The following table is used, in line with the ALEPH manual.
    
    |                  |   EXT |   META |
    |:-----------------|------:|-------:|
    | [275.0, 325.0)   |    03 |     31 |
    | [325.0, 375.0)   |    35 |     32 |
    | [375.0, 425.0)   |    04 |     33 |
    | [425.0, 475.0)   |    45 |     34 |
    | [475.0, 525.0)   |    05 |     36 |
    | [525.0, 575.0)   |    55 |     37 |
    | [575.0, 625.0)   |    06 |     38 |
    | [625.0, 675.0)   |    65 |     39 |
    | [675.0, 725.0)   |    07 |     40 |
    | [725.0, 775.0)   |    75 |     41 |
    | [775.0, 825.0)   |    08 |     42 |
    | [825.0, 875.0)   |    85 |     43 |
    | [875.0, 925.0)   |    09 |     44 |
    | [925.0, 975.0)   |    95 |     46 |
    | [975.0, 1050.0)  |    10 |     47 |
    | [1050.0, 1150.0) |    11 |     48 |
    | [1150.0, 1250.0) |    12 |     49 |
    | [1250.0, 1350.0) |    13 |     50 |
    | [1350.0, 1450.0) |    14 |     51 |
    | [1450.0, 1650.0) |    15 |     52 |
    | [1650.0, 1950.0) |    18 |     53 |
    | [1950.0, 2250.0) |    21 |     54 |

    If temperature is outside the interval range given above, suffix `'00'` is
    returned.
    If temperature is 293.6 K, suffix `'02'` and `'30'` are returned
    respectively for ground and meta states.
    
    Parameters
    ----------
    temperature : `float`
        processing temperature
    meta : `bool`, optional, default is `True`
        `True` if metastable (any), `False` if ground level

    Returns
    -------
    `str`
        suffix

    Examples
    -------- 
    Test temperatures outside range
    >>> assert get_temperature_suffix(3000, True) == get_temperature_suffix(3000, False)
    >>> assert get_temperature_suffix(3000, True) == get_temperature_suffix(0, True)
    >>> assert get_temperature_suffix(3000, True) == get_temperature_suffix(0, False)
    
    Test reference ALEPH temperatures for ground and meta states
    >>> for k, v in sandy.njoy.tmp2ext.items():
    ...    assert v == get_temperature_suffix(k)
    ...    assert v == get_temperature_suffix(k, 0)

    >>> for k, v in sandy.njoy.tmp2ext_meta.items():
    ...    assert v == get_temperature_suffix(k, meta=True)
    ...    assert v == get_temperature_suffix(k, 1)
    ...    assert v == get_temperature_suffix(k, 2)
    """
    closed = "left"
    
    if temperature == 293.6:
        return "30" if meta else "02"

    # Up to 1000 K temperatures are split every 50 degrees.
    splitter = 50
    idx = pd.IntervalIndex.from_breaks(np.arange(300-splitter/2, 1000+splitter/2, splitter).tolist() + [(1100+1000)/2], closed=closed)
    suffix1 = pd.DataFrame({
        "EXT": ["03", "35", "04", "45", "05", "55", "06", "65", "07", "75", "08", "85", "09", "95", "10"],
        "META": ["31", "32", "33", "34", "36", "37", "38", "39", "40", "41", "42", "43", "44", "46", "47"],
    }, index=idx)

    # Between 1000 K and 1500 K temperatures are split every 100 degrees.
    splitter = 100
    idx = pd.IntervalIndex.from_breaks(np.arange(1100-splitter/2, 1500+splitter/2, splitter).tolist() + [(1800+1500)/2], closed=closed)
    suffix2 = pd.DataFrame({
        "EXT": ["11", "12", "13", "14", "15"],
        "META": ["48", "49", "50", "51", "52"],
    }, index=idx)
    suffix = pd.concat([suffix1, suffix2])

    # Between 1500 K and 2400 K temperatures are split every 300 degrees.
    splitter = 300
    idx = pd.IntervalIndex.from_breaks(np.arange(1800-splitter/2, 2400+splitter/2, splitter).tolist() + [2400+splitter/2], closed=closed)
    suffix3 = pd.DataFrame({
        "EXT": ["18", "21", "24"],
        "META": ["53", "54", "55"],
    }, index=idx)

    suffix = pd.concat([suffix1, suffix2, suffix3])
    
    mask = suffix.index.contains(temperature)
    if suffix[mask].empty:
        suff = "00"
        msg = f"extension '{suff}' will be used for temperature '{temperature}'"
        logging.warning(msg)
    else:
        if meta:
            suff = suffix[mask].META.squeeze()
        else:
            suff = suffix[mask].EXT.squeeze()
    return suff


def _moder_input(nin, nout, **kwargs):
    """
    Write moder input.

    Parameters
    ----------
    nin : `int`
        tape number for input file
    nout : `int`
        tape number for output file

    Returns
    -------
    `str`
        moder input text
    """
    text = ["moder"]
    text += [f"{nin:d} {nout:d} /"]
    return "\n".join(text) + "\n"


def _reconr_input(endfin, pendfout, mat,
                  header="sandy runs njoy",
                  err=NJOY_TOLER,
                  **kwargs):
    """
    Write reconr input.

    Parameters
    ----------
    endfin : `int`
        tape number for input ENDF-6 file
    pendfout : `int`
        tape number for output PENDF file
    mat : `int`
        MAT number
    err : `float`
        tolerance (default is 0.001)
    header : `str`
        file header (default is "sandy runs njoy")

    Returns
    -------
    `str`
        reconr input text
    """
    text = ["reconr"]
    text += [f"{endfin:d} {pendfout:d} /"]
    text += [f"'{header}'/"]
    text += [f"{mat:d} 0 0 /"]
    text += [f"{err} 0. /"]
    text += ["0/"]
    return "\n".join(text) + "\n"


def _broadr_input(endfin, pendfin, pendfout, mat,
                  temperatures=NJOY_TEMPERATURES,
                  err=NJOY_TOLER,
                  **kwargs):
    """
    Write broadr input.

    Parameters
    ----------
    endfin : `int`
        tape number for input ENDF-6 file
    pendfin : `int`
        tape number for input PENDF file
    pendfout : `int`
        tape number for output PENDF file
    mat : `int`
        MAT number
    err : `float`
        tolerance (default is 0.001)
    temperatures : iterable of `float`
        iterable of temperature values in K (default is 293.6 K)

    Returns
    -------
    `str`
        broadr input text
    """
    text = ["broadr"]
    text += [f"{endfin:d} {pendfin:d} {pendfout:d} /"]
    ntemps = len(temperatures)
    text += [f"{mat:d} {ntemps:d} 0 0 0. /"]
    text += [f"{err} /"]
    text += [" ".join(map("{:.1f}".format, temperatures)) + " /"]
    text += ["0 /"]
    return "\n".join(text) + "\n"


def _thermr_input(endfin, pendfin, pendfout, mat,
                  temperatures=NJOY_TEMPERATURES,
                  angles=20,
                  iprint=False,
                  err=NJOY_TOLER,
                  emax=NJOY_THERMR_EMAX,
                  **kwargs):
    """
    Write thermr input for free-gas.

    Parameters
    ----------
    endfin : `int`
        tape number for input ENDF-6 file
    pendfin : `int`
        tape number for input PENDF file
    pendfout : `int`
        tape number for output PENDF file
    mat : `int`
        MAT number
    temperatures : iterable of `float`
        iterable of temperature values in K (default is 293.6 K)
    angles : `int`
        number of equi-probable angles (default is 20)
    iprint : `bool`
        print option (default is `False`)
    err : `float`
        tolerance (default is 0.001)
    emax : `float`
        maximum energy for thermal treatment (default is 10 eV)

    Returns
    -------
    `str`
        thermr input text
    """
    text = ["thermr"]
    text += [f"{endfin:d} {pendfin:d} {pendfout:d} /"]
    ntemps = len(temperatures)
    printflag = int(iprint)
    text += [f"0 {mat:d} {angles:d} {ntemps:d} 1 0 0 1 221 {printflag:d} /"]
    text += [" ".join(map("{:.1f}".format, temperatures)) + " /"]
    text += [f"{err} {emax} /"]
    return "\n".join(text) + "\n"


def _purr_input(endfin, pendfin, pendfout, mat,
                temperatures=NJOY_TEMPERATURES,
                sig0=NJOY_SIG0,
                bins=20,
                ladders=32,
                iprint=False,
                **kwargs):
    """
    Write purr input.

    Parameters
    ----------
    endfin : `int`
        tape number for input ENDF-6 file
    pendfin : `int`
        tape number for input PENDF file
    pendfout : `int`
        tape number for output PENDF file
    mat : `int`
        MAT number
    temperatures : iterable of `float`
        iterable of temperature values in K (default is 293.6 K)
    sig0 : iterable of `float`
        iterable of dilution values in barns (default is 1e10 b)
    bins : `int`
        number of probability bins (default is 20)
    ladders : `int`
        number of resonance ladders (default is 32)
    iprint : `bool`
        print option (default is `False`)

    Returns
    -------
    `str`
        purr input text
    """
    text = ["purr"]
    text += [f"{endfin:d} {pendfin:d} {pendfout:d} /"]
    ntemps = len(temperatures)
    nsig0 = len(sig0)
    printflag = int(iprint)
    text += [f"{mat:d} {ntemps:d} {nsig0:d} {bins:d} {ladders:d} {printflag:d} /"]
    text += [" ".join(map("{:.1f}".format, temperatures)) + " /"]
    text += [" ".join(map("{:.2E}".format, sig0)) + " /"]
    text += ["0 /"]
    return "\n".join(text) + "\n"


def _gaspr_input(endfin, pendfin, pendfout,
                 **kwargs):
    """
    Write gaspr input.

    Parameters
    ----------
    endfin : `int`
        tape number for input ENDF-6 file
    pendfin : `int`
        tape number for input PENDF file
    pendfout : `int`
        tape number for output PENDF file

    Returns
    -------
    `str`
        gaspr input text
    """
    text = ["gaspr"]
    text += [f"{endfin:d} {pendfin:d} {pendfout:d} /"]
    return "\n".join(text) + "\n"


def _unresr_input(endfin, pendfin, pendfout, mat,
                  temperatures=NJOY_TEMPERATURES,
                  sig0=NJOY_SIG0,
                  iprint=False,
                  **kwargs):
    """
    Write unresr input.

    Parameters
    ----------
    endfin : `int`
        tape number for input ENDF-6 file
    pendfin : `int`
        tape number for input PENDF file
    pendfout : `int`
        tape number for output PENDF file
    mat : `int`
        MAT number
    temperatures : iterable of `float`
        iterable of temperature values in K (default is 293.6 K)
    sig0 : iterable of `float`
        iterable of dilution values in barns (default is 1e10 b)
    iprint : `bool`
        print option (default is `False`)

    Returns
    -------
    `str`
        unresr input text
    """
    text = ["unresr"]
    text += [f"{endfin:d} {pendfin:d} {pendfout:d} /"]
    ntemps = len(temperatures)
    nsig0 = len(sig0)
    printflag = int(iprint)
    text += [f"{mat:d} {ntemps:d} {nsig0:d} {printflag:d} /"]
    text += [" ".join(map("{:.1f}".format, temperatures)) + " /"]
    text += [" ".join(map("{:.2E}".format, sig0)) + " /"]
    text += ["0 /"]
    return "\n".join(text) + "\n"


def _heatr_input(endfin, pendfin, pendfout, mat, pks,
                 local=False,
                 iprint=False,
                 **kwargs):
    """
    Write heatr input.

    Parameters
    ----------
    endfin : `int`
        tape number for input ENDF-6 file
    pendfin : `int`
        tape number for input PENDF file
    pendfout : `int`
        tape number for output PENDF file
    mat : `int`
        MAT number
    pks : iterable of `int`
        iterable of MT numbers for partial kermas
    local : `bool`
        option to deposit gamma rays locally (default is `False`)
    iprint : `bool`
        print option (default is `False`)

    Returns
    -------
    `str`
        heatr input text
    """
    text = ["heatr"]
    text += [f"{endfin:d} {pendfin:d} {pendfout:d} 0 /"]
    npks = len(pks)
    localflag = int(local)
    printflag = int(iprint)
    text += [f"{mat:d} {npks:d} 0 0 {localflag:d} {printflag:d} /"]
    text += [" ".join(map("{:d}".format, pks)) + " /"]
    return "\n".join(text) + "\n"


def _acer_input(endfin, pendfin, aceout, dirout, mat,
                temperature=NJOY_TEMPERATURES[0],
                iprint=False,
                itype=1,
                suffix=".00",
                header="sandy runs acer",
                photons=True,
                **kwargs):
    """
    Write acer input for fast data.

    Parameters
    ----------
    endfin : `int`
        tape number for input ENDF-6 file
    pendfin : `int`
        tape number for input PENDF file
    aceout : `int`
        tape number for output ACE file
    dirout : `int`
        tape number for output ACE file
    mat : `int`
        MAT number
    header : `str`
        descriptive character string of max. 70 characters
        (default is "sandy runs acer")
    iprint : `bool`
        print option (default is `False`)
    itype : `int`
        ace output type: 1, 2, or 3 (default is 1)
    photons : `bool`
        detailed photons (default is `True`)
    suffix : `str`
        id suffix for zaid (default is ".00")
    temperature : `float`
        temperature in K (default is 293.6 K)

    Returns
    -------
    `str`
        acer input text
    """
    text = ["acer"]
    text += [f"{endfin:d} {pendfin:d} 0 {aceout:d} {dirout:d} /"]
    printflag = int(iprint)
    text += [f"1 {printflag:d} {itype:d} {suffix} 0 /"]
    text += [f"'{header}'/"]
    text += [f"{mat:d} {temperature:.1f} /"]
    photonsflag = int(photons)
    text += [f"1 {photonsflag:d} /"]
    text += ["/"]
    return "\n".join(text) + "\n"


def _errorr_input(endfin, pendfin, gendfin, errorrout, mat,
                  ign=2, ek=None, spectrum=None,
                  iwt=2, relative=True,
                  mt=None, irespr=1,
                  temperature=NJOY_TEMPERATURES[0], mfcov=33,
                  iprint=False,
                  **kwargs):
    """
    Write errorr input.

    Parameters
    ----------
    endfin : `int`
        tape number for input ENDF-6 file
    pendfin : `int`
        tape number for input PENDF file
    gendfin : `int`
        tape number for input GENDF file
    errorrout : `int`
        tape number for output ERRORR file
    mat : `int`
        MAT number
    ek : iterable, optional
        derived cross section energy bounds (default is None)
    ign : `int`, optional
        neutron group option (default is 2, csewg 239-group structure)

        .. note:: this parameter will not be used if keyword argument
                  `ek` is provided.

    iprint : `bool`, optional
        print option (default is `False`)
    irespr: `int`, optional
        processing for resonance parameter covariances (default is 1)
        - 0: area sensitivity method
        - 1: 1% sensitivity method
    iwt : `int`, optional
        weight function option (default is 2, constant)

        .. note:: this parameter will not be used if keyword argument
                  `spect` is provided.

    relative: `bool`
        use relative covariance form (default is `True`)
    mfcov : `int`
        endf covariance file to be processed (default is 33)
    mt : `int` or iterable of `int`, optional
        run errorr only for the selected mt numbers
        (default is `None`, i.e., process all MT)

        .. note:: this parameter will not be used if keyword argument
                  `mfcov!=33`.

    spectrum : iterable, optional
        weight function (default is `None`)
    temperature : `float`, optional
        temperature in K (default is 293.6 K)

    Returns
    -------
    `str`
        errorr input text

    Examples
    --------
    Default test without keyword arguments
    >>> print(sandy.njoy._errorr_input(20, 21, 0, 22, 9237))
    errorr
    20 21 0 22 0 /
    9237 2 2 0 1 /
    0 293.6 /
    0 33 1/

    Test argument `temperature`
    >>> print(sandy.njoy._errorr_input(20, 21, 0, 22, 9440, temperature=600))
    errorr
    20 21 0 22 0 /
    9440 2 2 0 1 /
    0 600.0 /
    0 33 1/

    Test argument `iwt`
    >>> print(sandy.njoy._errorr_input(20, 21, 0, 22, 9237, iwt=6))
    errorr
    20 21 0 22 0 /
    9237 2 6 0 1 /
    0 293.6 /
    0 33 1/

    Test argument `ek`
    >>> print(sandy.njoy._errorr_input(20, 21, 0, 22, 9237, ek=[1e-2, 1e3, 2e5]))
    errorr
    20 21 0 22 0 /
    9237 1 2 0 1 /
    0 293.6 /
    0 33 1/
    2 /
    1.00000e-02 1.00000e+03 2.00000e+05 /

    Test argument `ign`
    >>> print(sandy.njoy._errorr_input(20, 21, 0, 22, 9237, ign=3))
    errorr
    20 21 0 22 0 /
    9237 3 2 0 1 /
    0 293.6 /
    0 33 1/

    Test nubar
    >>> print(sandy.njoy._errorr_input(20, 21, 0, 22, 9237, mfcov=31))
    errorr
    20 21 0 22 0 /
    9237 2 2 0 1 /
    0 293.6 /
    0 31 1/

    Test mubar
    >>> print(sandy.njoy._errorr_input(20, 21, 0, 22, 9237, mfcov=34))
    errorr
    20 21 0 22 0 /
    9237 2 2 0 1 /
    0 293.6 /
    0 34 1/

    Test chi
    >>> print(sandy.njoy._errorr_input(20, 21, 0, 22, 9237, mfcov=35))
    errorr
    20 21 0 22 0 /
    9237 2 2 0 1 /
    0 293.6 /
    0 35 1/

    Test keyword `relative`
    >>> print(sandy.njoy._errorr_input(20, 21, 0, 22, 9237, relative=False))
    errorr
    20 21 0 22 0 /
    9237 2 2 0 0 /
    0 293.6 /
    0 33 1/

    Test keyword `irespr`
    >>> print(sandy.njoy._errorr_input(20, 21, 0, 22, 9237, irespr=0))
    errorr
    20 21 0 22 0 /
    9237 2 2 0 1 /
    0 293.6 /
    0 33 0/

    Test keyword `mt` as `list`
    >>> print(sandy.njoy._errorr_input(20, 21, 0, 22, 9237, mt=[1, 2]))
    errorr
    20 21 0 22 0 /
    9237 2 2 0 1 /
    0 293.6 /
    1 33 1/
    2 0 /
    1 2 /    

    Test keyword `mt` as `int`:
    >>> print(sandy.njoy._errorr_input(20, 21, 0, 22, 9237, mt=2))
    errorr
    20 21 0 22 0 /
    9237 2 2 0 1 /
    0 293.6 /
    1 33 1/
    1 0 /
    2 /    
    """
    irelco = 0 if relative is False else 1
    iread = 1 if (mt is not None and mfcov == 33) else 0
    iwt_ = 1 if spectrum is not None else iwt
    ign_ = 1 if ek is not None else ign
    text = ["errorr"]
    text += [f"{endfin:d} {pendfin:d} {gendfin:d} {errorrout:d} 0 /"]
    printflag = int(iprint)
    text += [f"{mat:d} {ign_:d} {iwt_:d} {printflag:d} {irelco} /"]
    text += [f"{printflag:d} {temperature:.1f} /"]
    text += [f"{iread:d} {mfcov} {irespr:d}/"]
    if iread == 1:  # only specific mts
        mtlist = [mt] if isinstance(mt, int) else mt
        nmt = len(mtlist)
        text += [f"{nmt:d} 0 /"]
        text += [" ".join(map(str, mtlist)) + " /"]
    if ign_ == 1:
        nk = len(ek) - 1
        text += [f"{nk} /"]
        text += [" ".join(map("{:.5e}".format, ek)) + " /"]
    if iwt_ == 1:
        INT = 1               # constant interpolation
        NBT = int(len(spectrum) / 2)  # only 1 interpolation group
        tab1 = "\n".join(sandy.write_tab1(0, 0, 0, 0, [NBT], [INT],
                                          spectrum[::2],
                                          spectrum[1::2]))
        text += [tab1]
        text += ["/"]
    return "\n".join(text) + "\n"


def _groupr_input(endfin, pendfin, gendfout, mat,
                  ign=2, ek=None, igg=0, ep=None,
                  iwt=2, lord=0, sigz=[1e+10],
                  temperature=NJOY_TEMPERATURES[0],
                  spectrum=None, mt=None,
                  iprint=False,
                  mubar=False, chi=False, nubar=False,
                  **kwargs):
    """
    Write GROUPR input

    Parameters
    ----------
    endfin : `int`
        tape number for input ENDF-6 file
    pendfin : `int`
        tape number for input PENDF file
    gendfout : `int`
        tape number for output PENDF file
    mat : `int`
        MAT number
    chi : `bool`, optional
        Process chi (default is `False`)
    ek : iterable, optional
        derived cross section energy bounds (default is None)
    ep : iterable, optional
        derived gamma cross section energy bounds (default is None)
    igg : `int`, optional
        gamma group option (default is 0, no structure)
    ign : `int`, optional
        neutron group option (default is 2, csewg 239-group structure)
    iprint : `bool`, optional
        print option (default is `False`)
    iwt : `int`, optional
        weight function option (default is 2, constant)
        
        .. note:: this parameter will not be used if keyword argument
                  `spectrum` is provided

    lord : `int`, optional
        Legendre order (default is 0)
    mt: `int` or iterable of `int`, optional
        run groupr only for the selected MT numbers
        (default is `None`, i.e., process all MT)
    mubar : `bool`, optional
        Proccess mubar (default is `False`)
    nubar : `bool`, optional
        Proccess nubar (default is `False`)
    nuclide_production : `bool`, optional
        process MF10 (default is `False`)
    sigz : iterable of `float`
        sigma zero values (he default is 1.0e10)
    spectrum : iterable, optional
        weight function (default is `None`)
    temperature : iterable of `float`
        iterable of temperature values in K (default is 293.6 K)

    Returns
    -------
    `str`
        groupr input text

    Examples
    --------
    Default test without keyword arguments
    >>> print(sandy.njoy._groupr_input(20, 21, 22, 9237))
    groupr
    20 21 0 22 /
    9237 2 0 2 0 1 1 0 /
    'sandy runs groupr' /
    293.6/
    10000000000.0/
    3/
    0/
    0/

    Test argument `temperature`
    >>> print(sandy.njoy._groupr_input(20, 21, 22, 9440, temperature=600))
    groupr
    20 21 0 22 /
    9440 2 0 2 0 1 1 0 /
    'sandy runs groupr' /
    600.0/
    10000000000.0/
    3/
    0/
    0/

    Test argument `iwt`
    >>> print(sandy.njoy._groupr_input(20, 21, 22, 9237, iwt=6))
    groupr
    20 21 0 22 /
    9237 2 0 6 0 1 1 0 /
    'sandy runs groupr' /
    293.6/
    10000000000.0/
    3/
    0/
    0/

    Test argument `ign`
    >>> print(sandy.njoy._groupr_input(20, 21, 22, 9237, ign=3))
    groupr
    20 21 0 22 /
    9237 3 0 2 0 1 1 0 /
    'sandy runs groupr' /
    293.6/
    10000000000.0/
    3/
    0/
    0/

    Test argument `igg`
    >>> print(sandy.njoy._groupr_input(20, 21, 22, 9237, igg=3))
    groupr
    20 21 0 22 /
    9237 2 3 2 0 1 1 0 /
    'sandy runs groupr' /
    293.6/
    10000000000.0/
    3/
    0/
    0/
    
    Test argument `ek`
    >>> print(sandy.njoy._groupr_input(20, 21, 0, 22, 9237, ek=[1e-2, 1e3, 2e5]))
    groupr
    20 21 0 0 /
    22 1 0 2 0 1 1 0 /
    'sandy runs groupr' /
    293.6/
    10000000000.0/
    2 /
    1.00000e-02 1.00000e+03 2.00000e+05 /
    3/
    0/
    0/

    Test argument `ep`
    >>> print(sandy.njoy._groupr_input(20, 21, 0, 22, 9237, ep=[1e-2, 1e3, 2e5]))
    groupr
    20 21 0 0 /
    22 9237 1 2 0 1 1 0 /
    'sandy runs groupr' /
    293.6/
    10000000000.0/
    2 /
    1.00000e-02 1.00000e+03 2.00000e+05 /
    3/
    0/
    0/

    Test argument `lord`
    >>> print(sandy.njoy._groupr_input(20, 21, 0, 22, 9237, lord=3))
    groupr
    20 21 0 0 /
    22 9237 0 2 3 1 1 0 /
    'sandy runs groupr' /
    293.6/
    10000000000.0/
    3/
    0/
    0/

    Test mubar:
    >>> print(sandy.njoy._groupr_input(20, 21, 22, 9237, mubar=True))
    groupr
    20 21 0 22 /
    9237 2 0 2 0 1 1 0 /
    'sandy runs groupr' /
    293.6/
    10000000000.0/
    3/
    3 251 /
    0/
    0/

    Test chi:
    >>> print(sandy.njoy._groupr_input(20, 21, 22, 9237, chi=True))
    groupr
    20 21 0 22 /
    9237 2 0 2 0 1 1 0 /
    'sandy runs groupr' /
    293.6/
    10000000000.0/
    3/
    5/
    5 18 /
    0/
    0/

    Test keyword `mt` as `list`
    >>> print(sandy.njoy._groupr_input(20, 21, 0, 22, 9237, mt=[1, 2]))
    groupr
    20 21 0 0 /
    22 9237 0 2 0 1 1 0 /
    'sandy runs groupr' /
    293.6/
    10000000000.0/
    3 1 /
    3 2 /
    0/
    0/       

    Test keyword `mt` as `int`
    >>> print(sandy.njoy._groupr_input(20, 21, 0, 22, 9237, mt=2))
    groupr
    20 21 0 0 /
    22 9237 0 2 0 1 1 0 /
    'sandy runs groupr' /
    293.6/
    10000000000.0/
    3 2 /
    0/
    0/       
    """
    iwt_ = 1 if spectrum is not None else iwt
    ign_ = 1 if ek is not None else ign
    igg_ = 1 if ep is not None else igg
    sigzlist = sigz if hasattr(sigz, "__len__") else [sigz]
    text = ["groupr"]
    text += [f"{endfin:d} {pendfin:d} 0 {gendfout:d} /"]
    nsigz = len(sigzlist)
    printflag = int(iprint)
    text += [f"{mat:d} {ign_:d} {igg_:d} {iwt_:d} {lord:d} 1 {nsigz:d} {printflag:d} /"]
    text += ["'sandy runs groupr' /"]  # run label
    text += [f"{temperature:.1f}/"]
    text += [" ".join(map("{:.1f}".format, sigzlist)) + "/"]
    if ign_ == 1:
        nk = len(ek) - 1
        text += [f"{nk} /"]
        text += [" ".join(map("{:.5e}".format, ek)) + " /"]
    if igg_ == 1:
        pk = len(ep) - 1
        text += [f"{pk} /"]
        text += [" ".join(map("{:.5e}".format, ep)) + " /"]
    if iwt_ == 1:
        INT = 1               # constant interpolation
        NBT = int(len(spectrum) / 2)  # only 1 interpolation group
        tab1 = "\n".join(sandy.write_tab1(0, 0, 0, 0, [NBT], [INT],
                                          spectrum[::2],
                                          spectrum[1::2]))
        text += [tab1]
        text += ["/"]

    if mt is None:
        text += ["3/"]  # by default process all cross sections (MF=3)
    else:
        mtlist = mt if hasattr(mt, "__len__") else [mt]
        for mt_ in mtlist:
            text += [f"3 {mt_:d} /"]
    if nubar:
        text += [f"3 452 /"]
        text += [f"3 455 /"]
        text += [f"3 456 /"]
    if mubar:
        text += [f"3 251 /"]
    if chi:
        text += ["5/"]
        text += ["5 18 /"]
    text += ["0/"]  # terminate list of reactions for this material
    text += ["0/"]  # terminate materials (only 1 allowed)
    return "\n".join(text) + "\n"


def _run_njoy(text, endf, pendf=None, exe=None):
    """
    Run njoy executable for given input.

    Parameters
    ----------
    inputs : `map`
        map of {`tape`: `file`) for input files
    outputs : `map`
        map of {`tape`: `file`) for ouptut files
    text : `str`
        njoy input file passed to `Popen` as `stdin` (it must be encoded first)
    exe : `str`, optional, default is `None`
        njoy executable: if `None` (default) get it from `NJOY` env variable
    """
    if exe is None:
        exe = get_njoy()
    stdout = stderr = None
    stdin = text.encode()
    with TemporaryDirectory() as tmpdir:
        shutil.copy(endf, os.path.join(tmpdir, "tape20"))
        if pendf:
            shutil.copy(pendf, os.path.join(tmpdir, "tape99"))
        process = sp.Popen(exe,
                           shell=True,
                           cwd=tmpdir,
                           stdin=sp.PIPE,
                           stdout=stdout,
                           stderr=stderr)
        stdoutdata, stderrdata = process.communicate(input=stdin)
        logging.debug(stdoutdata)
        logging.debug(stderrdata)
        retrn = process.returncode
        if retrn != 0:
            cwd = os.getcwd()
            dir = join(cwd, "njoy_outputs")
            if os.path.exists(dir):
                shutil.rmtree(dir)
            os.makedirs(dir)
            for filename in os.listdir(tmpdir):
                shutil.copy(join(tmpdir, filename), join(dir, filename))
            with open(join(dir, "input"), "w") as f:
                f.write(stdin.decode("utf-8"))
            msg = f"process status={retrn} when running njoy executable '{exe}'.\nInputs/Outputs were moved to '{dir}'"
            raise ValueError(msg)

        # Move outputs
        tapes = {
            30: "pendf",
            40: "gendf",
            31: "errorr31",
            33: "errorr33",
            34: "errorr34",
            35: "errorr35",
            }
        outputs = {}
        for k, v in tapes.items():
            out = os.path.join(tmpdir, f"tape{k}")
            if os.path.exists(out):
                with open(out, mode="r") as f:
                    outputs[v] = f.read()
        text = ""
        for k in range(50, 70):
            out = os.path.join(tmpdir, f"tape{k}")
            if os.path.exists(out):
                with open(out, mode="r") as f:
                    text += f.read()
        if text:
            outputs["ace"] = text
        text = ""
        for k in range(70, 90):
            out = os.path.join(tmpdir, f"tape{k}")
            if os.path.exists(out):
                with open(out, mode="r") as f:
                    text += f.read()
        if text:
            outputs["xsdir"] = text
        return outputs


def _prepare_njoy_input(
        mat, temperatures, suffixes,
        acer=True,
        broadr=True,
        gaspr=True,
        groupr=False,
        heatr=True,
        purr=True,
        reconr=True,
        thermr=True,
        unresr=False,
        errorr33=False,
        errorr31=False,
        errorr34=False,
        errorr35=False,
        acer_kws={},
        broadr_kws={},
        gaspr_kws={},
        groupr_kws={},
        heatr_kws={},
        purr_kws={},
        reconr_kws={},
        thermr_kws={},
        unresr_kws={},
        errorr31_kws={},
        errorr33_kws={},
        errorr34_kws={},
        errorr35_kws={},
        **kwargs,
        ):
    """
    Parameters
    ----------
    broadr : `bool`, optional, default is `True`
        option to run module broadr
    thermr : `bool`, optional, default is `True`
        option to run module thermr
    unresr : `bool`, optional, default is `False`
        option to run module unresr
    heatr : `bool`, optional, default is `True`
        option to run module heatr
    gaspr : `bool`, optional, default is `True`
        option to run module gapr
    purr : `bool`, optional, default is `True`
        option to run module purr
    groupr : `bool`, optional, default is `False`
        option to run module groupr
    errorr : `bool`, optional, default is `False`
        option to run module errorr
    acer : `bool`, optional, default is `True`
        option to run module acer

    Notes
    -----
    .. note:: the four calls to NJOY module ERRORR (for MF31, MF33, MF34
              and MF35) are treated as if four different NJOY modules were
              to be run.
    """
    
    e = 21
    p = e + 1
    text = _moder_input(20, -e)

    # this part produces a single PENDF file
    if reconr:
        text += _reconr_input(-e, -p,
                              mat=mat, temperatures=temperatures,
                              **reconr_kws)
    else:
        text += _moder_input(99, -p)
    if broadr:
        o = p + 1
        text += _broadr_input(-e, -p, -o,
                              mat=mat, temperatures=temperatures,
                              **broadr_kws)
        p = o
    if thermr:
        o = p + 1
        text += _thermr_input(0, -p, -o,
                              mat=mat, temperatures=temperatures,
                              **thermr_kws)
        p = o
    if unresr:
        o = p + 1
        text += _unresr_input(-e, -p, -o,
                              mat=mat, temperatures=temperatures,
                              **unresr_kws)
        p = o
    if heatr:
        kermas=[302, 303, 304, 318, 402, 442, 443, 444, 445, 446, 447]
        for i in range(0, len(kermas), 7):
            o = p + 1
            pks = kermas[i:i+7]
            text += _heatr_input(-e, -p, -o,
                                 mat=mat, temperatures=temperatures, pks=pks,
                                 **heatr_kws)
            p = o
    if gaspr:
        o = p + 1
        text += _gaspr_input(-e, -p, -o,
                             mat=mat, temperatures=temperatures,
                             **gaspr_kws)
        p = o
    if purr:
        o = p + 1
        text += _purr_input(-e, -p, -o,
                            mat=mat, temperatures=temperatures,
                            **purr_kws)
        p = o
    o = 30
    text += _moder_input(-p, o)
    
    # this part produces a single GENDF file
    if errorr31 or errorr34 or errorr35:
        # groupr is needed by ERRORR31, 34 and 35
        # this is the only place where we make this check
        groupr_ = True
    else:
        groupr_ = groupr
    g = 39
    if groupr_:
        if len(temperatures) > 1:
            logging.info("Multiple temperatures were requested.\nGROUPR will only process the first.")
        temperature = temperatures[0]
        text += _groupr_input(-e, -p, -g, 
                              mat=mat, temperature=temperature,
                              **groupr_kws)
        o = 40
        text += _moder_input(-g, o)

    # this part produces a maximimum of four ERRORR files, one per data type
    if errorr33 or errorr31 or errorr34 or errorr35:
        if len(temperatures) > 1:
            logging.info("Multiple temperatures were requested.\nERRORR will only process the first.")
        temperature = temperatures[0]
    if errorr33:
        # for xs use a GENDF file only if explicitely asked, not just if
        # groupr_=True because chi or nubar were present
        o = errorr33_kws["mfcov"] = 33
        p_ = 0 if groupr else p
        g_ = g if groupr else 0
        text += _errorr_input(-e, -p_, -g_, o,
                              mat=mat, temperature=temperature,
                              **errorr33_kws)
    if errorr31:
        o = errorr31_kws["mfcov"] = 31
        text += _errorr_input(-e, 0, -g, o,
                              mat=mat, temperature=temperature,
                              **errorr31_kws)
    # NJOY's errorr module WILL produce a MF35 covariance tape
    # if the errorr module called before the errorr call to produce a MF34
    if errorr35:
        o = errorr35_kws["mfcov"] = 35
        text += _errorr_input(-e, 0, -g, o,
                              mat=mat, temperature=temperature,
                              **errorr35_kws)
    if errorr34:
        o = errorr34_kws["mfcov"] = 34
        text += _errorr_input(-e, 0, -g, o,
                              mat=mat, temperature=temperature,
                              **errorr34_kws)

    # this part produces multiple ACE files (one per temperature)
    if acer:
        for i, (temperature, suffix) in enumerate(zip(temperatures, suffixes)):
            a = 50 + i
            x = 70 + i
            text += _acer_input(-e, -p, a, x,
                                mat=mat, temperature=temperature, suffix=suffix,
                                **acer_kws)
    text += "stop"
    return text


def process_neutron(
        endftape,
        temperatures,
        pendftape=None,
        suffixes=None,
        zaid="nndc",
        route="0",
        exe=None,
        verbose=True,
        dryrun=False,
        **kwargs,
        ):
    """
    Run sequence to process file with njoy.

    Parameters
    ----------
    pendftape : `str`, optional, default is `None`
        name (with absolute of relative path) of a pendf file.
        If given, skip module reconr and use this PENDF file, else run reconr
    temperatures : iterable of `float`, optional, default is [293.6]
        iterable of temperature values in K
    dryrun : `bool`, optional, default is `False`
        option to produce the njoy input file without running njoy
    exe : `str`, optional, default is `None`
        njoy executable (with path)
        .. note:: if no executable is given, SANDY looks for a default
                  executable in `PATH` and in env variable `NJOY`
    route : `str`, optional, default is `0`
        xsdir "route" parameter
    suffixes : iterable of `int`, optional, default is `None`
        iterable of suffix values for ACE files: if `None` is given,
        use internal routine to determine suffixes
        
        .. warning :: `suffixes` must match the number of entries in
                     `temperatures`

    verbose : `bool`, optional, default is `False`
        flag to print NJOY input to screen before running the executable

    Returns
    -------
    outputs : `map`
        map of {`tape` : `text`) for ouptut files
    """
    tape = sandy.Endf6.from_file(endftape)
    mat = tape.mat[0]
    info = tape.read_section(mat, 1, 451)
    za = int(info["ZA"])
    meta = info["LISO"]

    kwargs["reconr"] = False if pendftape else True

    # Prepare njoy input
    temperatures_ = temperatures if hasattr(temperatures, "__len__") else [temperatures]
    if suffixes:
        suffixes_ = suffixes if hasattr(suffixes, "__len__") else [suffixes]
    else:
        meta_ = 0 if zaid == "nndc" else meta
        suffixes_ = ["." + get_temperature_suffix(t, meta_) for t in temperatures_]
    text = _prepare_njoy_input(mat, temperatures_, suffixes_, **kwargs)
    if verbose:
        print(text)

    # Run njoy
    if dryrun:
        return text
    outputs = _run_njoy(text, endftape, pendftape, exe=exe)
    
    # Minimal output post-processing
    if "xsdir" in outputs:
        outputs["xsdir"] = outputs["xsdir"].replace("route", route)
    if zaid == "nndc":
        za_new = sandy.zam.zam2za(za*10 + meta, method=zaid)[0]
        for s in suffixes_:
            pattern = f"{za}{s}[c]"
            new_pattern = f"{za_new}{s}c"
            if "xsdir" in outputs:
                outputs["xsdir"] = re.sub(pattern, new_pattern, outputs["xsdir"])
            if "ace" in outputs:
                outputs["ace"] = re.sub(pattern, new_pattern, outputs["ace"])
    return outputs


def process_proton(
        endftape,
        wdir="",
        dryrun=False,
        tag="",
        exe=None,
        route="0",
        **kwargs,
        ):
    """Run sequence to process proton file with njoy.

    Parameters
    ----------
    wdir : `str`
        working directory (absolute or relative) where all output files are
        saved
        .. note:

            `wdir` will appear as part of the `filename` in
            any `xsdir` file
    dryrun : `bool`
        option to produce the njoy input file without running njoy
    tag : `str`
        tag to append to each output filename beofre the extension
        (default is `None`)
        .. hint:
            to process JEFF-3.3 files you could set `tag = "_j33"`
    exe : `str`
        njoy executable (with path)
        .. note:
            If no executable is given, SANDY looks for a default executable
            in `PATH`
    route : `str`
        xsdir "route" parameter (default is "0")

    Returns
    -------
    input : `str`
        njoy input text
    inputs : `map`
        map of {`tape` : `file`) for input files
    outputs : `map`
        map of {`tape` : `file`) for ouptut files
    """
    tape = sandy.Endf6.from_file(endftape)
    mat = tape.mat[0]
    info = tape.read_section(mat, 1, 451)
    meta = info["LISO"]
    za = int(info["ZA"])
    za_new = za + meta*100 + 300 if meta else za
    inputs = {}
    outputs = {}
    kwargs["mat"] = mat
    inputs["tape20"] = endftape
    kwargs["temp"] = 0
    kwargs["suff"] = suff = ".00"
    text = _acer_input(20, 20, 50, 70, **kwargs)
    outputs["tape50"] = join(wdir, "{}{}{}h".format(za_new, tag, suff))
    outputs["tape70"] = join(wdir, "{}{}{}h.xsd".format(za_new, tag, suff))
    text += "stop"
    if not dryrun:
        _run_njoy(text, inputs, outputs, exe=exe)
        # Change route and filename in xsdir file.
        acefile = outputs["tape50"]
        xsdfile = outputs["tape70"]
        text_xsd = open(xsdfile).read() \
                                .replace("route", route) \
                                .replace("filename", acefile)
        text_xsd = " ".join(text_xsd.split())
        # If isotope is metatable rewrite ZA in xsdir and ace as
        # ZA = Z*1000 + 300 + A + META*100.
        if meta:
            pattern = f'{za:d}' + r'.(?P<ext>\d{2}[ct])'
            found = re.search(pattern, text_xsd)
            ext = found.group("ext")
            text_xsd = text_xsd.replace("{:d}.{}".format(za, ext), "{:d}.{}".format(za_new, ext), 1)
            text_ace = open(acefile).read()
            text_ace = text_ace.replace("{:d}.{}".format(za, ext), "{:d}.{}".format(za_new, ext), 1)
            with open(acefile, 'w') as f:
                f.write(text_ace)
        with open(xsdfile, 'w') as f:
            f.write(text_xsd)
    return text, inputs, outputs
