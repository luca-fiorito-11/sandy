# -*- coding: utf-8 -*-
"""
Created on Wed Mar 21 15:21:37 2018

@author: fiorito_l
"""
import os
import argparse
import pdb

__author__ = "Luca Fiorito"

class SandyError(Exception):
    pass

def str2bool(v):
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

def is_valid_file(parser, arg, r=True, w=False, x=False):
    arg = os.path.abspath(os.path.realpath(os.path.normpath(arg)))
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
    arg = os.path.abspath(os.path.realpath(os.path.normpath(arg)))
    if os.path.isdir(arg):
        return arg
    if mkdir:
        os.makedirs(arg, exist_ok=True)
    else:
        parser.error("Directory {} does not exist".format(arg))
    return arg


def init_plotter(iargs=None):
    global args
    parser = argparse.ArgumentParser(description='Run SANDY xs-plotter')
    parser.add_argument('file',
                        type=lambda x: is_valid_file(parser, x),
                        help="ENDF-6 or PENDF format file.")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--endf6-cov',
                       type=lambda x: is_valid_file(parser, x),
                       help="ENDF-6 file containing covariances.")
    group.add_argument('--errorr-cov',
                       type=lambda x: is_valid_file(parser, x),
                       help="ERRORR file containing covariances.")
    parser.add_argument('--samples',
                        type=int,
                        default=100,
                        help="Number of samples.")
    parser.add_argument('--outdir',
                        default=os.getcwd(),
                        type=lambda x: is_valid_dir(parser, x, mkdir=True),
                        help="Target directory where outputs are stored (default = current working directory). If it does not exist it will be created.")
    parser.add_argument('-np','--processes',
                        type=int,
                        default=1,
                        help="Number of worker processes (default=1).")
    parser.add_argument('--plotdir',
                        default=os.path.join(os.getcwd(),"html_files"),
                        type=lambda x: is_valid_dir(parser, x, mkdir=True),
                        help="Target directory where plots are stored (default = current working directory/html_files). If it does not exist it will be created.")
    args = parser.parse_known_args(args=iargs)[0]
    return vars(args)


def init_macs(iargs=None):
    parser = argparse.ArgumentParser(description=None)
    parser.add_argument('--pendf',
                       type=lambda x: is_valid_file(parser, x),
                       help="PENDF file.")
    parser.add_argument('--errorr',
                       type=lambda x: is_valid_file(parser, x),
                       help="ERRORR file.")
    parser.add_argument('--kT',
                        type=float,
                        default=[0.0253],
                        nargs='+',
                        help="Maxwellian temperature kT in eV (default=[0.0253]).")
    parser.add_argument('-mt','--listmt',
                        type=int,
                        default=[1,2,18,102],
                        nargs='+',
                        help="List of MT sections (default=[1,2,18,102]).")
    args = parser.parse_known_args(args=iargs)[0]
    return args

def init_test_ace(iargs=None):
    parser = argparse.ArgumentParser(description=None)
    parser.add_argument('xsdir',
                       type=lambda x: is_valid_file(parser, x),
                       help="xsdir file")
    parser.add_argument('-a','--ace-files',
                       type=lambda x: is_valid_file(parser, x),
                       metavar="ace",
                       default=argparse.SUPPRESS,
                       action="store",
                       nargs='+',
                       help="List of ACE files")
    args = parser.parse_known_args(args=iargs)[0]
    return args