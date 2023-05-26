import os
import time
import logging
import argparse
import filecmp

import sandy
from sandy.tools import is_valid_dir, is_valid_file


__author__ = "Luca Fiorito"
__all__ = []


def parse(iargs=None):
    """Parse command line arguments for sampling option.

    Parameters
    ----------
    iargs : `list` of `str` or `None`, default is `None`,
        list of strings to parse.
        The default is taken from `sys.argv`.

    Returns
    -------
    `argparse.Namespace`
        namespace object containing processed given arguments and/or default
        options.
    """
    description = "Produce perturbed files containing sampled parameters that represent the information stored in the evaluated nuclear data covariances."""
    parser = argparse.ArgumentParser(
                        prog="sandy",
                        description=description,
                        formatter_class=argparse.RawTextHelpFormatter,
                        )

    parser.add_argument('file',
                        type=lambda x: is_valid_file(parser, x),
                        help="ENDF-6 file")

    parser.add_argument('--acer',
                        default=False,
                        action="store_true",
                        help="Process each perturbed file into ACE format "
                        "(default = False)\n"
                        "(--temperatures is required)")

    parser.add_argument('--mat',
                        type=int,
                        default=list(range(1, 10000)),
                        action='store',
                        nargs="+",
                        metavar="{1,..,9999}",
                        help="draw samples only from the selected MAT "
                             "sections (default = keep all)")

    parser.add_argument('--mf',
                        type=int,
                        default=[31, 33, 34, 35],
                        action='store',
                        nargs="+",
                        metavar="{31,33,34,35}",
                        help="draw samples only from the selected MF sections "
                             "(default = keep all)")

    parser.add_argument('--mt33',
                        type=int,
                        default=None,
                        action='store',
                        nargs="+",
                        metavar="{1,..,999}",
                        help="draw samples only from the selected MT sections "
                             "(default = keep all)")

    parser.add_argument('--njoy',
                        type=lambda x: is_valid_file(parser, x),
                        default=None,
                        help="NJOY executable "
                             "(default search PATH, and env variable NJOY)")

    parser.add_argument('--outname', '-O',
                        type=str,
                        default="{ZA}_{SMP}",
                        help="name template for the output files\n"
                             "(use formatting options in https://pyformat.info/ ,\n"
                             "available keywords are MAT, ZAM, ZA, META, SMP)")

    parser.add_argument('--processes', '-N',
                        type=int,
                        default=1,
                        help="number of worker processes (default = 1)")

    parser.add_argument('--samples', '-S',
                        type=int,
                        default=200,
                        help="number of samples (default = 200)")

    parser.add_argument('--seed31',
                        type=int,
                        default=None,
                        metavar="S31",
                        help="seed for random sampling of MF31 covariance "
                             "matrix (default = random)")

    parser.add_argument('--seed33',
                        type=int,
                        default=None,
                        metavar="S33",
                        help="seed for random sampling of MF33 covariance "
                             "matrix (default = random)")

    parser.add_argument('--seed34',
                        type=int,
                        default=None,
                        metavar="S34",
                        help="seed for random sampling of MF34 covariance "
                             "matrix (default = random)")

    parser.add_argument('--seed35',
                        type=int,
                        default=None,
                        metavar="S35",
                        help="seed for random sampling of MF35 covariance "
                             "matrix (default = random)")

    parser.add_argument('--temperatures', '-T',
                        default=None,
                        type=float,
                        action='store',
                        nargs="+",
                        metavar="T",
                        help="for each perturbed file, produce ACE files at "
                             "given temperatures")

    parser.add_argument("--version", "-v",
                        action='version',
                        version='%(prog)s {}'.format(sandy.__version__),
                        help="SANDY's version.")

    parser.add_argument('--verbose',
                        default=False,
                        action="store_true",
                        help="print additional details to screen during execution")

    init = parser.parse_known_args(args=iargs)[0]
    if init.acer and not init.temperatures:
        parser.error("--acer requires --temperatures")
    return init


def run(cli="--help"):
    """
    

    Parameters
    ----------
    cli : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    Examples
    --------
    Retrieve ENDF-6 tape and write it to file.
    >>> sandy.get_endf6_file("jeff_33", "xs", 10010).to_file("H1.jeff33")

    Produce perturbed ACE file.
    >>> cli = "H1.jeff33 --acer True --samples 2 --processes 2 --temperatures 900 --seed33 5"
    >>> sandy.sampling.run(cli)

    Check if ACE and XSDIR files have the right content.
    >>> assert "1001.09c" in open("1001_0.09c").read()
    >>> assert "1001.09c" in open("1001_0.09c.xsd").read()
    >>> assert "1001.09c" in open("1001_1.09c").read()
    >>> assert "1001.09c" in open("1001_1.09c.xsd").read()
    >>> assert not filecmp.cmp("1001_0.09c", "1001_1.09c")

    Run the same on a single process.
    >>> cli = "H1.jeff33 --acer True --samples 2 --processes 2 --temperatures 900 --seed33 5 --outname={ZAM}_{SMP}_SP"
    >>> sandy.sampling.run(cli)

    The identical seed ensures consistent results with the previous run.
    >>> assert filecmp.cmp("1001_0.09c", "10010_0_SP.09c")
    >>> assert filecmp.cmp("1001_1.09c", "10010_1_SP.09c")
    >>> assert filecmp.cmp("1001_0.09c.xsd", "10010_0_SP.09c.xsd")
    >>> assert filecmp.cmp("1001_1.09c.xsd", "10010_1_SP.09c.xsd")

    Produce perturbed ENDF6 and PENDF files.
    >>> cli = "H1.jeff33 --samples 2 --processes 2 --outname=H1_{MAT}_{SMP} --mt 102"
    >>> sandy.sampling.run(cli)
    >>> assert os.path.getsize("H1_125_0.pendf") > 0 and os.path.getsize("H1_125_1.pendf") > 0

    >>> assert filecmp.cmp("H1_125_0.endf6", "H1_125_1.endf6")
    >>> assert filecmp.cmp("H1_125_0.endf6", "H1.jeff33")
    """
    # >>> assert not filecmp.cmp("H1_125_0.pendf", "H1_125_1.pendf")
    # >>> cli = "H1.jeff33 --samples 2 --processes 2 --seed33 5 --outname=H1_{SMP}"
    # >>> sandy.sampling.run(cli)
    t0 = time.time()
    
    iargs = parse(cli.split())

    endf6 = sandy.Endf6.from_file(iargs.file)
    njoy_kws = {}
    if iargs.mt33:
        njoy_kws["errorr33_kws"] = dict(mt=iargs.mt33)

    smp_kws = {}
    if iargs.seed31:
        smp_kws["seed31"] = iargs.seed31
    if iargs.seed33:
        smp_kws["seed33"] = iargs.seed33
    if iargs.seed34:
        smp_kws["seed34"] = iargs.seed34
    if iargs.seed34:
        smp_kws["seed35"] = iargs.seed35

    smps = endf6.get_perturbations(iargs.samples, njoy_kws=njoy_kws, smp_kws=smp_kws)
    if iargs.temperatures:
        temperature = iargs.temperatures[0] if hasattr(iargs.temperatures, "__len__") else iargs.temperatures
    else:
        temperature = 0
        
    endf6.apply_perturbations(
        smps,
        processes=iargs.processes,
        to_file=True,
        to_ace=iargs.acer,
        filename=iargs.outname,
        njoy_kws=dict(err=0.005),
        ace_kws=dict(temperature=temperature, err=0.005),
        verbose=iargs.verbose,
    )

    dt = time.time() - t0
    logging.info(f"Total running time: {dt:.2f} sec")


if __name__ == "__main__":
    run()
