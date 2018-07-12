import argparse
import os


from .sampling import sampling
from ..utils import is_valid_dir, is_valid_file

def from_cli(iargs=None):
    """
    Parse command line arguments and return ```argparse.Namespace``` instance
    holding attributes.
    """
    parser = argparse.ArgumentParser(prog="python -m sandy.sampling",
                                     description='Run sampling',
                                     formatter_class=argparse.RawTextHelpFormatter)
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
                        required=True,
                        help="number of samples")
    parser.add_argument('--outdir',
                        metavar="DIR",
                        default=os.getcwd(),
                        type=lambda x: is_valid_dir(parser, x, mkdir=True),
                        help="target directory where outputs are stored\n(default = current working directory)\nif it does not exist it will be created")
    parser.add_argument('-np','--processes',
                        type=int,
                        default=1,
                        help="number of worker processes (default = 1)")
    parser.add_argument('--eig',
                        type=int,
                        default=0,
                        metavar="N",
                        help="print the first N eigenvalues of the evaluated covariance matrices\n(default = do not print)")
    parser.add_argument('--mat',
                        type=int,
                        action='store',
                        nargs="+",
                        metavar="{1,..,9999}",
                        help="draw samples only from the selected MAT sections (default = keep all)")
    parser.add_argument('--mf',
                        type=int,
                        default=range(41),
                        action='store',
                        nargs="+",
                        metavar="{1,..,40}",
                        help="draw samples only from the selected MF sections (default = keep all)")
    parser.add_argument('--mt',
                        type=int,
                        action='store',
                        nargs="+",
                        metavar="{1,..,999}",
                        help="draw samples only from the selected MT sections (default = keep all)")
    parser.add_argument('--outname',
                        type=str,
                        help="basename for the output files (default is the the basename of <file>.)")
    parser.add_argument('--verbose',
                        default=False,
                        action="store_true",
                        help="turn on verbosity (default = quiet)")
    parser.add_argument('-e','--energy-points',
                        type=float,
                        metavar="E",
                        default=[],
                        action="store",
                        nargs='+',
                        help="additional energy points (in eV) to include in the incoming-neutron energy grid\n(default = None)")
    parser.add_argument("-v",
                        '--version',
                        action='version',
                        version='%(prog)s 1.0',
                        help="SANDY's version.")
    args = parser.parse_known_args(args=iargs)[0]
    return args