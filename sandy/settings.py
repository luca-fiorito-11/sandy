# -*- coding: utf-8 -*-
"""
Created on Wed Mar 21 15:21:37 2018

@author: fiorito_l
"""
import os, argparse, pdb, textwrap

class EmptyIsConst(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        if len(values) == 0:
            string = option_string.lower()
            if string == "--temps":
                values = [293.6]
            elif string == "--kerma":
                values = [302, 303, 304, 318, 402, 442, 443, 444, 445, 446, 447]
            elif string == "--sig0":
                values = [1E10]
        setattr(namespace, self.dest, values)

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

def init_njoy(iargs=None):
    parser = argparse.ArgumentParser(description='Run NJOY', formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('tape',
                        type=lambda x: is_valid_file(parser, x),
                        help="ENDF-6 format file.")
    parser.add_argument('--broadr',
                        type=str2bool,
                        default=argparse.SUPPRESS,
                        metavar=True,
                        help=argparse.SUPPRESS)
    parser.add_argument('--purr',
                        type=str2bool,
                        default=argparse.SUPPRESS,
                        metavar=True,
                        help=argparse.SUPPRESS)
    parser.add_argument('--unresr',
                        type=str2bool,
                        default=argparse.SUPPRESS,
                        metavar=False,
                        help=argparse.SUPPRESS)
    parser.add_argument('--thermr',
                        type=str2bool,
                        default=argparse.SUPPRESS,
                        metavar=True,
                        help=argparse.SUPPRESS)
    parser.add_argument('--heatr',
                        type=str2bool,
                        default=argparse.SUPPRESS,
                        metavar=False,
                        help=argparse.SUPPRESS)
    parser.add_argument('--acer',
                        type=str2bool,
                        default=argparse.SUPPRESS,
                        metavar=False,
                        help=argparse.SUPPRESS)
    parser.add_argument('-p','--pendf',
                        action="store_true",
                        help="produce PENDF file")
    parser.add_argument('-a','--ace',
                        action="store_true",
                        help="produce ACE file")
    parser.add_argument('-g','--gendf',
                        action="store_true",
                        help="produce GENDF file")
    parser.add_argument('-e','--errorr',
                        action="store_true",
                        help="produce ERRORR file")
#    parser.add_argument('--capsys',
#                        action="store_true",
#                        help="capture NJOY stderr and stdout (default=False)")
#    parser.add_argument('--NjoyExec',
#                        default='njoy2016',
#                        help="NJOY executable (default=njoy2016).")
#    parser.add_argument('--evaluationName',
#                        type=lambda x: is_valid_dir(parser, x, mkdir=True),
#                        default=argparse.SUPPRESS,
#                        metavar="lib",
#                        help="Name of the evaluation.")
#    parser.add_argument('--iprint',
#                        type=int,
#                        choices=range(2),
#                        default=argparse.SUPPRESS,
#                        help="NJOY verbosity: 0=min, 1=max (default=1).")
    parser.add_argument('--temps',
                        type=float,
                        default = argparse.SUPPRESS,
                        nargs='+',
                        help="list of temperature values")
    parser.add_argument('--sig0',
                        type=float,
                        default = argparse.SUPPRESS,
                        nargs='+',
                        help="list of dilution values")
    parser.add_argument('--kerma',
                        type=int,
                        action=EmptyIsConst,
                        default = argparse.SUPPRESS,
                        nargs='*',
                        help="partial kermas (default=None)")
#    parser.add_argument('--qa',
#                        default=[],
#                        nargs='+',
#                        help="MT number and associated user Q-value (default=None).")
    parser.add_argument('--err',
                        type=float,
                        default=argparse.SUPPRESS,
                        help="fractional tolerance for RECONR and BROADR (default=0.001)")
    parser.add_argument('--free-gas',
                        default=argparse.SUPPRESS,
                        action="store_true",
                        help="compute thermal cross section for free-gas scatterer (THERMR)")
    parser.add_argument('--ptable',
                        default=argparse.SUPPRESS,
                        action="store_true",
                        help="compute probability tables (PURR)")
    parser.add_argument('--gaspr',
                        default=argparse.SUPPRESS,
                        action="store_true",
                        help="compute gas-production cross sections (GASPR)")
    parser.add_argument('--ign',
                        default=argparse.SUPPRESS,
                        help='''Neutron group structure option for GROUPR
 file : read from file
 - 2 : csewg 239-group structure [default]
 - 3 : lanl 30-group structure
 - 4 : anl 27-group structure
 - 5 : rrd 50-group structure
 - 6 : gam-i 68-group structure
 - 7 : gam-ii 100-group structure
 - 8 : laser-thermos 35-group structure
 - 9 : epri-cpm 69-group structure
 - 10 : lanl 187-group structure
 - 11 : lanl 70-group structure
 - 12 : sand-ii 620-group structure
 - 13 : lanl 80-group structure
 - 14 : eurlib 100-group structure
 - 15 : sand-iia 640-group structure
 - 16 : vitamin-e 174-group structure
 - 17 : vitamin-j 175-group structure
 - 19 : ecco 33-group structure
 - 20 : ecco 1968-group structure
 - 21 : tripoli 315-group structure
 - 22 : xmas lwpc 172-group structure
 - 23 : vit-j lwpc 175-group structure''')
    parser.add_argument('--iwt',
                        default=argparse.SUPPRESS,
                        help="Weight function option for GROUPR (default=6 : Maxwellian - 1/E - fission - fusion).")
    parser.add_argument('--igne',
                        default=argparse.SUPPRESS,
                        help="Neutron group structure option for ERRORR (default=2 : csewg 239-group structure).")
    parser.add_argument('--iwte',
                        default=argparse.SUPPRESS,
                        help="Weight function option for ERRORR (default=6 : Maxwellian - 1/E - fission - fusion).")
    parser.add_argument('--plot',
                        action="store_true",
                        help="Activate plotting capabilities.")
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
    return args

def init_sampling(iargs=None):
    global args
    parser = argparse.ArgumentParser(description='Run sampling', formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('file',
                        type=lambda x: is_valid_file(parser, x),
                        help="ENDF-6 or PENDF format file")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--endf6-cov',
                       type=lambda x: is_valid_file(parser, x),
                       help="ENDF-6 file containing covariances")
    group.add_argument('--errorr-cov',
                       type=lambda x: is_valid_file(parser, x),
                       help="ERRORR file containing covariances")
    parser.add_argument('--samples',
                        type=int,
                        default=100,
                        help="number of samples (default = 100)")
    parser.add_argument('--outdir',
                        metavar="DIR",
                        default=os.getcwd(),
                        type=lambda x: is_valid_dir(parser, x, mkdir=True),
                        help="target directory where outputs are stored (default = current working directory)\nif it does not exist it will be created")
    parser.add_argument('-np','--processes',
                        type=int,
                        default=1,
                        help="number of worker processes (default = 1)")
    parser.add_argument('--eig',
                        type=int,
                        default=0,
                        metavar="N",
                        help="print the first N eigenvalues of the evaluated covariance matrices (default = do not print)")
#    parser.add_argument('--plotdir',
#                        default=os.path.join(os.getcwd(),"html_files"),
#                        type=lambda x: is_valid_dir(parser, x, mkdir=True),
#                        help="Target directory where plots are stored (default = current working directory/html_files). If it does not exist it will be created.")
#    parser.add_argument('-p',
#                        action='store_true',
#                        help="Turn on xs plotting.")
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
    parser.add_argument('-e','--energy-points',
                        type=float,
                        metavar="E",
                        default=[],
                        action="store",
                        nargs='*',
                        help="additional energy points (in eV) to include in the incoming-neutron energy grid (default = None)")
    parser.add_argument("-v",
                        '--version',
                        action='version',
                        version='%(prog)s 1.0',
                        help="SANDY's version.")
    args = parser.parse_known_args(args=iargs)[0]
    return args


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