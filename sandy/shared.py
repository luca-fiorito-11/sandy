"""
Collection of utilities, functions and classes that are requetsed in all code components.
"""

# import re
# from scipy.constants import Avogadro
# import numpy as np
# import scipy.interpolate


__author__ = "Luca Fiorito"


# MeV_MWs = 1.60217733e-19 # conversion coefficient MeV --> MW s
# Amn = 1.00866491578      # molar mass of a neutron in amu, source: P. J. Mohr and B. N. Taylor, "The 1998 CODATA Recommended Values of the Fundamental Physics Constants", Version 3.1



# def expand_za(za, method="nndc", meta=0):
#     z = int(za//1000)
#     a = int(za - z*1000)
#     if method == "nndc":
#         m = 0
#         if a >= 300:
#             m = 1
#             a = a - 300 - m*100
#     else:
#         m = int(meta)
#     return z, a, m

# def get_za(z, a, m, method="nndc"):
#     za = z*1000 + a + 300 + m*100 if m != 0 and method == "nndc" else z*1000 + a
#     return int(za), m

# def expand_zam(zam):
#     z = int(zam//10000)
#     a = int(zam - z*10000)//10
#     m = int(zam - z*10000 - a*10)
#     return z, a, m

# def get_zam(z, a, m):
#     zam = z*10000 + a*10 + m
#     return int(zam)

# def za2zam(za, method="nndc", meta=0):
#     return get_zam(*expand_za(za, method=method, meta=meta))

# def zam2za(zam, method="nndc"):
#     z, a, m = expand_zam(zam)
#     return get_za(z, a, m, method=method)


