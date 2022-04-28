import os

import pandas as pd

import sandy

__author__ = "Luca Fiorito"
__all__ = [
        "ELEMENTS",
        "ATOMIC_NUMBERS",
        "METASTATES",
        "METASTATES_FLIP",
        "NATURAL_ABUNDANCE",
        "abundance_per_element",
        "expand_za",
        "expand_za",
        "za2latex",
        "zam2latex",
        ]


pd.options.display.float_format = '{:.5e}'.format


ELEMENTS = {
    1: 'H',
    2: 'He',
    3: 'Li',
    4: 'Be',
    5: 'B',
    6: 'C',
    7: 'N',
    8: 'O',
    9: 'F',
    10: 'Ne',
    11: 'Na',
    12: 'Mg',
    13: 'Al',
    14: 'Si',
    15: 'P',
    16: 'S',
    17: 'Cl',
    18: 'Ar',
    19: 'K',
    20: 'Ca',
    21: 'Sc',
    22: 'Ti',
    23: 'V',
    24: 'Cr',
    25: 'Mn',
    26: 'Fe',
    27: 'Co',
    28: 'Ni',
    29: 'Cu',
    30: 'Zn',
    31: 'Ga',
    32: 'Ge',
    33: 'As',
    34: 'Se',
    35: 'Br',
    36: 'Kr',
    37: 'Rb',
    38: 'Sr',
    39: 'Y',
    40: 'Zr',
    41: 'Nb',
    42: 'Mo',
    43: 'Tc',
    44: 'Ru',
    45: 'Rh',
    46: 'Pd',
    47: 'Ag',
    48: 'Cd',
    49: 'In',
    50: 'Sn',
    51: 'Sb',
    52: 'Te',
    53: 'I',
    54: 'Xe',
    55: 'Cs',
    56: 'Ba',
    57: 'La',
    58: 'Ce',
    59: 'Pr',
    60: 'Nd',
    61: 'Pm',
    62: 'Sm',
    63: 'Eu',
    64: 'Gd',
    65: 'Tb',
    66: 'Dy',
    67: 'Ho',
    68: 'Er',
    69: 'Tm',
    70: 'Yb',
    71: 'Lu',
    72: 'Hf',
    73: 'Ta',
    74: 'W',
    75: 'Re',
    76: 'Os',
    77: 'Ir',
    78: 'Pt',
    79: 'Au',
    80: 'Hg',
    81: 'Tl',
    82: 'Pb',
    83: 'Bi',
    84: 'Po',
    85: 'At',
    86: 'Rn',
    87: 'Fr',
    88: 'Ra',
    89: 'Ac',
    90: 'Th',
    91: 'Pa',
    92: 'U',
    93: 'Np',
    94: 'Pu',
    95: 'Am',
    96: 'Cm',
    97: 'Bk',
    98: 'Cf',
    99: 'Es',
    100: 'Fm',
    101: 'Md',
    102: 'No',
    103: 'Lr',
    104: 'Rf',
    105: 'Db',
    106: 'Sg',
    107: 'Bh',
    108: 'Hs',
    109: 'Mt',
    110: 'Ds',
    111: 'Rg',
    112: 'Uub',
    113: 'Uut',
    114: 'Uuq',
    115: 'Uup',
    116: 'Uuh',
    117: 'Uus',
    118: 'UUp',
    }

ATOMIC_NUMBERS = {v: k for k, v in ELEMENTS.items()}

METASTATES = {
    0: "g",
    1: "m",
    2: "n",
    3: "o",
    }

METASTATES_FLIP = {v: k for k, v in METASTATES.items()}

NATURAL_ABUNDANCE = {
    10010: 0.99984426,
    10020: 0.00015574,
    20030: 2e-06,
    20040: 0.999998,
    30060: 0.07589,
    30070: 0.92411,
    40090: 1.0,
    50100: 0.1982,
    50110: 0.8018,
    60120: 0.988922,
    60130: 0.011078,
    70140: 0.996337,
    70150: 0.003663,
    80160: 0.9976206,
    80170: 0.000379,
    80180: 0.0020004,
    90190: 1.0,
    100200: 0.9048,
    100210: 0.0027,
    100220: 0.0925,
    110230: 1.0,
    120240: 0.78951,
    120250: 0.1002,
    120260: 0.11029,
    130270: 1.0,
    140280: 0.9222968,
    140290: 0.0468316,
    140300: 0.0308716,
    150310: 1.0,
    160320: 0.9504074,
    160330: 0.0074869,
    160340: 0.0419599,
    160360: 0.0001458,
    170350: 0.757647,
    170370: 0.242353,
    180360: 0.003336,
    180380: 0.000629,
    180400: 0.996035,
    190390: 0.932581,
    190400: 0.000117,
    190410: 0.067302,
    200400: 0.96941,
    200420: 0.00647,
    200430: 0.00135,
    200440: 0.02086,
    200460: 4e-05,
    200480: 0.00187,
    210450: 1.0,
    220460: 0.0825,
    220470: 0.0744,
    220480: 0.7372,
    220490: 0.0541,
    220500: 0.0518,
    230500: 0.0025,
    230510: 0.9975,
    240500: 0.04345,
    240520: 0.83789,
    240530: 0.09501,
    240540: 0.02365,
    250550: 1.0,
    260540: 0.05845,
    260560: 0.91754,
    260570: 0.02119,
    260580: 0.00282,
    270590: 1.0,
    280580: 0.680769,
    280600: 0.262231,
    280610: 0.011399,
    280620: 0.036345,
    280640: 0.009256,
    290630: 0.6915,
    290650: 0.3085,
    300640: 0.4917,
    300660: 0.2773,
    300670: 0.0404,
    300680: 0.1845,
    300700: 0.0061,
    310690: 0.60108,
    310710: 0.39892,
    320700: 0.2052,
    320720: 0.2745,
    320730: 0.0776,
    320740: 0.3652,
    320760: 0.0775,
    330750: 1.0,
    340740: 0.0086,
    340760: 0.0923,
    340770: 0.076,
    340780: 0.2369,
    340800: 0.498,
    340820: 0.0882,
    350790: 0.50686,
    350810: 0.49314,
    360780: 0.00355,
    360800: 0.02286,
    360820: 0.11593,
    360830: 0.115,
    360840: 0.56987,
    360860: 0.17279,
    370850: 0.7217,
    370870: 0.2783,
    380840: 0.0056,
    380860: 0.0986,
    380870: 0.07,
    380880: 0.8258,
    390890: 1.0,
    400900: 0.5145,
    400910: 0.1122,
    400920: 0.1715,
    400940: 0.1738,
    400960: 0.028,
    410930: 1.0,
    420920: 0.14649,
    420940: 0.09187,
    420950: 0.15873,
    420960: 0.16673,
    420970: 0.09582,
    420980: 0.24292,
    421000: 0.09744,
    440960: 0.0554,
    440980: 0.0187,
    440990: 0.1276,
    441000: 0.126,
    441010: 0.1706,
    441020: 0.3155,
    441040: 0.1862,
    451030: 1.0,
    461020: 0.0102,
    461040: 0.1114,
    461050: 0.2233,
    461060: 0.2733,
    461080: 0.2646,
    461100: 0.1172,
    471070: 0.51839,
    471090: 0.48161,
    481060: 0.01245,
    481080: 0.00888,
    481100: 0.1247,
    481110: 0.12795,
    481120: 0.24109,
    481130: 0.12227,
    481140: 0.28754,
    481160: 0.07512,
    491130: 0.04281,
    491150: 0.95719,
    501120: 0.0097,
    501140: 0.0066,
    501150: 0.0034,
    501160: 0.1454,
    501170: 0.0768,
    501180: 0.2422,
    501190: 0.0859,
    501200: 0.3258,
    501220: 0.0463,
    501240: 0.0579,
    511210: 0.5721,
    511230: 0.4279,
    521200: 0.0009,
    521220: 0.0255,
    521230: 0.0089,
    521240: 0.0474,
    521250: 0.0707,
    521260: 0.1884,
    521280: 0.3174,
    521300: 0.3408,
    531270: 1.0,
    541240: 0.00095,
    541260: 0.00089,
    541280: 0.0191,
    541290: 0.26401,
    541300: 0.04071,
    541310: 0.21232,
    541320: 0.26909,
    541340: 0.10436,
    541360: 0.08857,
    551330: 1.0,
    561300: 0.0011,
    561320: 0.001,
    561340: 0.0242,
    561350: 0.0659,
    561360: 0.0785,
    561370: 0.1123,
    561380: 0.717,
    571380: 0.0008881,
    571390: 0.9991119,
    581360: 0.00186,
    581380: 0.00251,
    581400: 0.88449,
    581420: 0.11114,
    591410: 1.0,
    601420: 0.27153,
    601430: 0.12173,
    601440: 0.23798,
    601450: 0.08293,
    601460: 0.17189,
    601480: 0.05756,
    601500: 0.05638,
    621440: 0.0308,
    621470: 0.15,
    621480: 0.1125,
    621490: 0.1382,
    621500: 0.0737,
    621520: 0.2674,
    621540: 0.2274,
    631510: 0.4781,
    631530: 0.5219,
    641520: 0.002,
    641540: 0.0218,
    641550: 0.148,
    641560: 0.2047,
    641570: 0.1565,
    641580: 0.2484,
    641600: 0.2186,
    651590: 1.0,
    661560: 0.00056,
    661580: 0.00095,
    661600: 0.02329,
    661610: 0.18889,
    661620: 0.25475,
    661630: 0.24896,
    661640: 0.2826,
    671650: 1.0,
    681620: 0.00139,
    681640: 0.01601,
    681660: 0.33503,
    681670: 0.22869,
    681680: 0.26978,
    681700: 0.1491,
    691690: 1.0,
    701680: 0.00123,
    701700: 0.02982,
    701710: 0.14086,
    701720: 0.21686,
    701730: 0.16103,
    701740: 0.32025,
    701760: 0.12995,
    711750: 0.97401,
    711760: 0.02599,
    721740: 0.0016,
    721760: 0.0526,
    721770: 0.186,
    721780: 0.2728,
    721790: 0.1362,
    721800: 0.3508,
    731800: 0.0001201,
    731810: 0.9998799,
    741800: 0.0012,
    741820: 0.265,
    741830: 0.1431,
    741840: 0.3064,
    741860: 0.2843,
    751850: 0.374,
    751870: 0.626,
    761840: 0.0002,
    761860: 0.0159,
    761870: 0.0196,
    761880: 0.1324,
    761890: 0.1615,
    761900: 0.2626,
    761920: 0.4078,
    771910: 0.373,
    771930: 0.627,
    781900: 0.00012,
    781920: 0.00782,
    781940: 0.32864,
    781950: 0.33775,
    781960: 0.25211,
    781980: 0.07356,
    791970: 1.0,
    801960: 0.0015,
    801980: 0.1004,
    801990: 0.1694,
    802000: 0.2314,
    802010: 0.1317,
    802020: 0.2974,
    802040: 0.0682,
    812030: 0.29524,
    812050: 0.70476,
    822040: 0.014,
    822060: 0.241,
    822070: 0.221,
    822080: 0.524,
    832090: 1.0,
    902300: 0.0002,
    902320: 0.9998,
    912310: 1.0,
    922340: 5.4e-05,
    922350: 0.007204,
    922380: 0.992742,
 }


def abundance_per_element():
    abundance_per_element = {
        expand_zam(zam)[0]: {} for zam in NATURAL_ABUNDANCE
        }
    for zam, v in NATURAL_ABUNDANCE.items():
        z, a, m = expand_zam(zam)
        abundance_per_element[z][zam] = v
    return abundance_per_element


def expand_za(za, method="nndc", meta=0):
    """
    Expand ZA number into atomic number, mass number and metastate number
    of the correspondent nuclide.

    Parameters
    ----------
    za : `int`
        nuclide ZA indicator
    method : `str`, optional, default is 'nndc'
        method of the representation of the metastable state in the ZA number
    meta : `int`, optional, default is 0
        metastable state of the nuclide

    Returns
    -------
    `int`
        atomic number
    `int`
        mass number
    `int`
        metastable state

    Notes
    -----
    .. note:: if the 'nndc' method is used, the specification of the `meta`
        parameter is not needed but the implementation of this method is
        performed only for nuclides at ground state or first isomeric state

    Examples
    --------
    >>> expand_za(92235)
    (92, 235, 0)

    >>> expand_za(95241)
    (95, 241, 0)

    >>> expand_za(95641, method="nndc")
    (95, 241, 1)

    >>> expand_za(95241, method=False, meta=1)
    (95, 241, 1)

    >>> expand_za(95241, method=False, meta=2)
    (95, 241, 2)
    """
    z = int(za//1000)
    a = int(za - z*1000)
    if method == "nndc":
        m = 0
        if a >= 300:
            m = 1
            a = a - 300 - m*100
    else:
        m = int(meta)
    return z, a, m


def get_za(z, a, m, method="nndc"):
    if m != 0 and method == "nndc":
        za = z*1000 + a + 300 + m*100
    else:
        za = z*1000 + a
    return int(za), m


def expand_zam(zam):
    z = int(zam//10000)
    a = int(zam - z*10000)//10
    m = int(zam - z*10000 - a*10)
    return z, a, m


def get_zam(z, a, m):
    zam = z*10000 + a*10 + m
    return int(zam)


def za2zam(za, method="nndc", meta=0):
    return get_zam(*expand_za(za, method=method, meta=meta))


def zam2za(zam, method="nndc"):
    z, a, m = expand_zam(zam)
    return get_za(z, a, m, method=method)


def za2latex(za):
    z, a, m = expand_za(za)
    string = "$^{" + f"{a}" + "}$"
    sym = ELEMENTS[z]
    string += f"{sym}"
    return string


def zam2latex(zam):
    z, a, m = expand_zam(zam)
    string = "$^{" + f"{a}"
    if m == 1:
        string += "m"
    elif m == 2:
        string += "n"
    string += "}$"
    sym = ELEMENTS[z]
    string += f"{sym}"
    return string


def zam2nuclide(zam, atomic_number=False, sep=""):
    """
    Convert ZAM to string such with symbol and mass, such as `922350` to
    `"U235"` or `952421` to `"Am242m"`.

    Parameters
    ----------
    zam : `int`
        nuclide ZAM indicator
    atomic_number : `bool`, optional, default is `False`
        flag to include the atomic number in the nuclide name
    sep : `str`, optional, default is `''`
        separation character(s) to place between the atomic number
        (if present), the element ID, and the mass number.

    Returns
    -------
    `string`
        nuclide expressed with symbol and mass

    Examples
    --------
    >>> zam2nuclide(922350)
    'U235'

    >>> zam2nuclide(922350, atomic_number=True)
    '92U235'

    >>> zam2nuclide(922350, atomic_number=True, sep="-")
    '92-U-235'

    >>> zam2nuclide(922350, atomic_number=False, sep="-")
    'U-235'

    >>> zam2nuclide(952420)
    'Am242'

    >>> zam2nuclide(952421)
    'Am242m'

    >>> zam2nuclide(952421, atomic_number=True, sep="_")
    '95_Am_242m'

    >>> zam2nuclide(952422)
    'Am242n'
    """
    z, a, m = expand_zam(zam)
    sym = ELEMENTS[z]
    meta = get_meta_letter(m, skip_ground=True)
    out = f"{sym}{sep}{a}{meta}"
    if atomic_number:
        out = f"{z}{sep}{out}"
    return out


def za2nuclide(za, method="nndc", meta=0, atomic_number=False, sep=""):
    """
    Convert ZA to string with symbol and mass, such as 92235 to
    "U235" or 95642 to "Am242m".

    Parameters
    ----------
    za : `int`
        nuclide ZA indicator
    method : `str`, optional, default is 'nndc'
        method of the representation of the metastable state in the ZA number
    meta : `int`, optional, default is 0
        metastable state of the nuclide
    atomic_number : bool, optional, default is False
        flag to include the atomic number in the nuclide name
    sep : `str`, optional, default is ''
        separation character(s) to place between the atomic number
        (if present), the element ID, and the mass number.

    Returns
    -------
    `string`
        nuclide expressed with symbol and mass number

    Notes
    -----
    .. note:: if the 'nndc' method is used, the specification of the `meta`
        parameter is not needed but the implementation of this method is
        performed only for nuclides at ground state or first isomeric state

    Examples
    --------
    >>> za2nuclide(92235, atomic_number=False)
    'U235'

    >>> za2nuclide(92235, atomic_number=True)
    '92U235'

    >>> za2nuclide(92235, atomic_number=True, sep="-")
    '92-U-235'

    >>> za2nuclide(92235, atomic_number=False, sep="-")
    'U-235'

    >>> za2nuclide(95242)
    'Am242'

    >>> za2nuclide(95642, method="nndc")
    'Am242m'

    >>> za2nuclide(95642, method="nndc", atomic_number=True, sep="_")
    '95_Am_242m'

    >>> za2nuclide(95242, method=False, meta=2)
    'Am242n'
    """
    zam = za2zam(za, method=method, meta=meta)
    return zam2nuclide(zam, atomic_number=atomic_number, sep=sep)


def nuclide2zam(nuclide, atomic_number=False, sep=""):
    """
    Convert string with symbol and mass number to ZAM, such as `"U235"` to
    `922350` or `"Am242m"` to `952421`.

    Parameters
    ----------
    nuclide : `str`
        nuclide expressed with symbol and mass.
    atomic_number : `bool`, optional, default is `False`
        flag to pass a string with the atomic number in `nuclide`
    sep : `str`, optional, default is `''`
        separation character(s) placed between the atomic number
        (if present), the element ID, and the mass number.

    Returns
    -------
    zam : `int`
        nuclide ZAM indicator

    Examples
    --------
    >>> nuclide2zam('U235')
    922350

    >>> nuclide2zam('92U235', atomic_number=True)
    922350

    >>> nuclide2zam('92-U-235', atomic_number=True, sep="-")
    922350

    >>> nuclide2zam('Am242')
    952420

    >>> nuclide2zam('Am242m')
    952421

    >>> nuclide2zam('95_Am_242m', atomic_number=True, sep="_")
    952421

    >>> nuclide2zam('Am242n')
    952422
    """
    sym = ""
    num = ""
    if not nuclide[-1].isalpha():
        nuclide_ = nuclide + "g"
    else:
        nuclide_ = nuclide
    if atomic_number and nuclide_[:2].isnumeric():
        nuclide_ = nuclide_[2:]
    elif atomic_number and nuclide_[:1].isnumeric():
        nuclide_ = nuclide_[1:]
    if sep != "":
        nuclide_.replace(sep, "")
    for i in nuclide_[:-1]:
        if i.isalpha():
            sym += i
        if i.isnumeric():
            num += i
    z = ATOMIC_NUMBERS[sym]
    a = int(num)
    m = METASTATES_FLIP[nuclide_[-1]]
    out = get_zam(z, a, m)
    return out


def nuclide2za(nuclide):
    zam = nuclide2zam(nuclide)
    return zam2za(zam)


def get_meta_letter(m, skip_ground=False):
    meta = METASTATES[m]
    if skip_ground and m == 0:
        meta = ""
    return meta

