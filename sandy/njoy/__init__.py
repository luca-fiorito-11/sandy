import os, argparse, pdb

from .njoy import Njoy, NjoyOutput
from .. import utils


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


def from_cli(iargs=None):
    parser = argparse.ArgumentParser(prog="python -m sandy.njoy",
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('tape',
                        type=lambda x: utils.is_valid_file(parser, x),
                        help="ENDF-6 format file")
    parser.add_argument('-P','--pendftape',
                        type=lambda x: utils.is_valid_file(parser, x),
                        default=argparse.SUPPRESS,
                        metavar="PENDF",
                        help="processed PENDF format file")
    parser.add_argument('-G','--gendftape',
                        type=lambda x: utils.is_valid_file(parser, x),
                        default=argparse.SUPPRESS,
                        metavar="GENDF",
                        help="processed GENDF format file")
    parser.add_argument('--broadr',
                        type=utils.str2bool,
                        default=argparse.SUPPRESS,
                        metavar=True,
                        help=argparse.SUPPRESS)
    parser.add_argument('--purr',
                        type=utils.str2bool,
                        default=argparse.SUPPRESS,
                        metavar=True,
                        help=argparse.SUPPRESS)
    parser.add_argument('--unresr',
                        type=utils.str2bool,
                        default=argparse.SUPPRESS,
                        metavar=False,
                        help=argparse.SUPPRESS)
    parser.add_argument('--thermr',
                        type=utils.str2bool,
                        default=argparse.SUPPRESS,
                        metavar=True,
                        help=argparse.SUPPRESS)
    parser.add_argument('--heatr',
                        type=utils.str2bool,
                        default=argparse.SUPPRESS,
                        metavar=False,
                        help=argparse.SUPPRESS)
    parser.add_argument('--acer',
                        type=utils.str2bool,
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
    parser.add_argument('-H','--hendf',
                        action="store_true",
                        help="produce HENDF file")
#    parser.add_argument('--capsys',
#                        action="store_true",
#                        help="capture NJOY stderr and stdout (default=False)")
    parser.add_argument('--exe',
                        default=argparse.SUPPRESS,
                        help="NJOY executable (default=njoy2016).")
    parser.add_argument('--temps',
                        type=float,
                        default=argparse.SUPPRESS,
                        nargs='+',
                        help="list of temperature values")
    parser.add_argument('--sig0',
                        type=float,
                        default=argparse.SUPPRESS,
                        nargs='+',
                        help="list of dilution values")
    parser.add_argument('--kerma',
                        type=int,
                        action=EmptyIsConst,
                        default=argparse.SUPPRESS,
                        nargs='*',
                        help="list of partial kermas (default=[302, 303, 304, 318, 402, 442, 443, 444, 445, 446, 447])")
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
                        help='''neutron group structure option for GROUPR
 - file : read from file
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
 - 23 : vit-j lwpc 175-group structure
predefined group structures:
 - scale_238''')
    parser.add_argument('--iwt',
                        default=argparse.SUPPRESS,
                        help='''Weight function option for GROUPR
 - file : read from file
 - 2 : constant
 - 3 : 1/e
 - 5 : epri-cell lwr
 - 6 : (thermal) -- (1/e) -- (fission + fusion)
 - 7 : same with t-dep thermal part
 - 8 : thermal--1/e--fast reactor--fission + fusion
 - 9 : claw weight function
 - 10 : claw with t-dependent thermal part
 - 11 : vitamin-e weight function (ornl-5505)
 - 12 : vit-e with t-dep thermal part
predefined functions:
 - jaea_fns_175
''')
    parser.add_argument('--igne',
                        default=argparse.SUPPRESS,
                        help="neutron group structure option for ERRORR (same as --ign)")
    parser.add_argument('--iwte',
                        default=argparse.SUPPRESS,
                        help="weight function option for ERRORR (same as --iwt)")
    parser.add_argument('--suffixes',
                        type=str,
                        default=argparse.SUPPRESS,
                        nargs='+',
                        metavar=".XX",
                        help="suffixes for ACE files, as many as temperature values (default = None).")
    parser.add_argument("-V","--verbose",
                        type=int,
                        default=argparse.SUPPRESS,
                        help="set verbosity level (default = 1)")
    parser.add_argument("-v",
                        '--version',
                        action='version',
                        version='%(prog)s 1.0',
                        help="")
    args = parser.parse_args(iargs)
    return args