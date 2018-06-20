# -*- coding: utf-8 -*-
"""
Created on Wed Mar 21 15:21:37 2018

@author: fiorito_l
"""
import os, argparse


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



def init_njoy(iargs=None):
    global args
    parser = argparse.ArgumentParser(description='Run NJOY')
    parser.add_argument('-i','--inputfile',
                        type=lambda x: is_valid_file(parser, x),
                        required=True,
                        help="<Required> List of endf files to be processed (one file per line).")
    parser.add_argument('--processes',
                        type=int,
                        default=1,
                        help="Number of worker processes (default=1).")
    parser.add_argument('--capsys',
                        action="store_true",
                        help="Capture NJOY stderr and stdout (default=False).")
    parser.add_argument('--NjoyExec',
                        default='njoy2016',
                        help="NJOY executable (default=njoy2016).")
    parser.add_argument('--evaluationName',
                        type=lambda x: is_valid_dir(parser, x, mkdir=True),
                        default="lib",
                        metavar="lib",
                        help="Name of the evaluation.")
    parser.add_argument('--iprint',
                        type=int,
                        choices=range(2),
                        default=1,
                        help="NJOY verbosity: 0=min, 1=max (default=1).")
    parser.add_argument('--no-reconr',
                        action="store_true",
                        help="<Developer only> Skip RECONR module in the NJOY sequence.")
    parser.add_argument('--no-broadr',
                        action="store_true",
                        help="Skip BROADR module in the NJOY sequence.")
    parser.add_argument('--no-thermr',
                        action="store_true",
                        help="Skip THERMR module in the NJOY sequence.")
    parser.add_argument('--no-purr',
                        action="store_true",
                        help="Skip PURR module in the NJOY sequence.")
    parser.add_argument('--unresr',
                        action="store_true",
                        help="Replace PURR module with UNRESR in the NJOY sequence.")
    parser.add_argument('--run-groupr',
                        action="store_true",
                        help="Run GROUPR module in the NJOY sequence.")
    parser.add_argument('--run-acer',
                        action="store_true",
                        help="Run ACER module in the NJOY sequence.")
    parser.add_argument('--run-errorr',
                        action="store_true",
                        help="Run ERRORR module in the NJOY sequence.")
    parser.add_argument('--temps',
                        type=float,
                        default = [293.6],
                        nargs='+',
                        help="Temperature values (default=[293.6]).")
    parser.add_argument('--sig0',
                        type=float,
                        default=[1e10],
                        nargs='+',
                        help="Sigma0 values (default=[1e10]).")
    parser.add_argument('--err',
                        type=float,
                        default=0.005,
                        help="Fractional tolerance for RECONR and BROADR (default=0.005).")
    parser.add_argument('--ign',
                        type=int,
                        default=2,
                        help="Neutron group structure option for GROUPR (default=2 : csewg 239-group structure).")
    parser.add_argument('--igne',
                        type=int,
                        default=2,
                        help="Neutron group structure option for ERRORR (default=2 : csewg 239-group structure).")
    parser.add_argument('--iwt',
                        type=int,
                        default=6,
                        help="Weight function option for GROUPR (default=6 : Maxwellian - 1/E - fission - fusion).")
    parser.add_argument('--iwte',
                        type=int,
                        default=6,
                        help="Weight function option for ERRORR (default=6 : Maxwellian - 1/E - fission - fusion).")
    parser.add_argument('--suffixes',
                        type=str,
                        default=[".03"],
                        nargs='+',
                        help="Suffixes for ACE files, as many as temperature values (default=[\".03\"]]).")
    parser.add_argument("-v",
                        '--version',
                        action='version',
                        version='%(prog)s 1.0',
                        help="")
    args = parser.parse_args()
    return vars(args)


def init_sampling(iargs=None):
    global args
    parser = argparse.ArgumentParser(description='Run SANDY')
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
                        help="Target directory were outputs are stored (default = current working directory). If it does not exist it will be created.")
    parser.add_argument('-np','--processes',
                        type=int,
                        default=1,
                        help="Number of worker processes (default=1).")
    parser.add_argument('--eig',
                        type=int,
                        default=0,
                        metavar="N",
                        help="Print the first N eigenvalues of the evaluated covariance matrices (default = 0, do not print).")
    parser.add_argument('--plotdir',
                        default=os.path.join(os.getcwd(),"html_files"),
                        type=lambda x: is_valid_dir(parser, x, mkdir=True),
                        help="Target directory where plots are stored (default = current working directory/html_files). If it does not exist it will be created.")
    parser.add_argument('-p',
                        action='store_true',
                        help="Turn on xs plotting.")
#    parser.add_argument('-mat','--keep-mat',
#                        type=int,
#                        action='append',
#                        help="Keep only the selected MAT sections (default = keep all). Allowed values range from 1 to 9999. Provide each MAT section as an individual optional argument, e.g. -mat 9228 -mat 125")
#    parser.add_argument('-mf','--keep-cov-mf',
#                        type=int,
#                        action='append',
#                        choices=range(31,36),
#                        default=[],
#                        help="Keep only the selected covariance MF sections (default = keep all). Allowed values are [31, 32, 33, 34, 35]. Provide each MF section as an individual optional argument, e.g. -mf 33 -mf 35")
#    parser.add_argument('-mt','--keep-cov-mt',
#                        type=int,
#                        action='append',
#                        metavar="{1,..,999}",
#                        help="Keep only the selected covariance MT sections (default = keep all). Allowed values range from 1 to 999. Provide each MT section as an individual optional argument, e.g. -mt 18 -mt 102")
    parser.add_argument('-e','--energy-point',
                        type=float,
                        action='append',
                        help="Additional energy points (in eV) to include in the incoming-neutron energy grid (default = None). Provide each energy point as an individual optional argument, e.g. -e 100.0 -e 201.5")
    parser.add_argument("-v",
                        '--version',
                        action='version',
                        version='%(prog)s 1.0',
                        help="SANDY's version.")
    args = parser.parse_known_args(args=iargs)[0]
    return args


def init_ri(iargs=None):
    parser = argparse.ArgumentParser(description=None)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--pendf',
                       type=lambda x: is_valid_file(parser, x),
                       help="PENDF file.")
    group.add_argument('--errorr',
                       type=lambda x: is_valid_file(parser, x),
                       help="ERRORR file.")
    args = parser.parse_known_args(args=iargs)[0]
    return args