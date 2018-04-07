# -*- coding: utf-8 -*-
"""
Created on Wed Mar 21 15:21:37 2018

@author: fiorito_l
"""
import os
import argparse


def is_valid_file(parser, arg, r=True, w=False, x=False):
    if not os.path.isfile(arg):
        parser.error("File {} does not exist".format(arg))
    if r and not os.access(arg, os.R_OK):
        parser.error("File {} is not readable".format(arg))
    if w and not os.access(arg, os.W_OK):
        parser.error("File {} is not writable".format(arg))
    if x and not os.access(arg, os.X_OK):
        parser.error("File {} is not executable".format(arg))
    return arg

def is_valid_dir(parser, arg, mkdir=False):
    if os.path.isdir(arg):
        return arg
    if mkdir:
        os.makedirs(arg, exist_ok=True)
    else:
        parser.error("Directory {} does not exist".format(arg))
    return arg


def init_sampling(ARGS=None):
    global args
    parser = argparse.ArgumentParser(description='Run SANDY')
    parser.add_argument('endf6',
                        type=lambda x: is_valid_file(parser, x),
                        help="ENDF-6 format file.")
    parser.add_argument('--pendf',
                        type=lambda x: is_valid_file(parser, x),
                        help="PENDF format file.")
    parser.add_argument('--samples',
                        type=int,
                        default=100,
                        help="Number of samples.")
    parser.add_argument('--outdir',
                        default=os.getcwd(),
                        type=lambda x: is_valid_dir(parser, x, mkdir=True),
                        help="Target directory were outputs are stored (default = current working directory). If it does not exist it will be created.")
    parser.add_argument('-np','--processes',
                        type=int,
                        default=None,
                        help="Number of worker processes. By default, the number returned by os.cpu_count() is used.")
    parser.add_argument('--eig',
                        type=int,
                        default=0,
                        metavar="N",
                        help="Print the first N eigenvalues of the evaluated covariance matrices (default = 0, do not print).")
    parser.add_argument('-mat','--keep-mat',
                        type=int,
                        action='append',
                        help="Keep only the selected MAT sections (default = keep all). Allowed values range from 1 to 9999. Provide each MAT section as an individual optional argument, e.g. -mat 9228 -mat 125")
    parser.add_argument('-mf','--keep-cov-mf',
                        type=int,
                        action='append',
                        choices=range(31,36),
                        help="Keep only the selected covariance MF sections (default = keep all). Allowed values are [31, 32, 33, 34, 35]. Provide each MF section as an individual optional argument, e.g. -mf 33 -mf 35")
    parser.add_argument('-mt','--keep-cov-mt',
                        type=int,
                        action='append',
                        help="Keep only the selected covariance MT sections (default = keep all). Allowed values range from 1 to 999. Provide each MT section as an individual optional argument, e.g. -mt 18 -mt 102")
    parser.add_argument('-e','--energy-point',
                        type=float,
                        action='append',
                        help="Additional energy points (in eV) to include in the incoming-neutron energy grid (default = None). Provide each energy point as an individual optional argument, e.g. -e 100.0 -e 201.5")
    parser.add_argument('--test',
                        default=None,
                        help="Run SANDY test case overwriting the other input options (default = None).")
    parser.add_argument('--njoy',
                        type=lambda x: is_valid_file(parser, x, x=True),
                        metavar="EXE",
                        help="NJOY exe file.")
#    parser.add_argument('-nw',
#                        '--no-write',
#                        help="Disable writing.")
    parser.add_argument("-v",
                        '--version',
                        action='version',
                        version='%(prog)s 1.0',
                        help="SANDY's version.")
    args = parser.parse_args(args=ARGS)