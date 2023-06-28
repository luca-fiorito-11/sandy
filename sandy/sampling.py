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
                        help="ENDF-6 file")

    parser.add_argument('--acer',
                        default=False,
                        action="store_true",
                        help="Process each perturbed file into ACE format "
                        "(default = False)\n"
                        "(--temperatures is required)")

    parser.add_argument('--debug',
                        default=False,
                        action="store_true",
                        help="activate debug options (err=1, verbose=True, minimal_processing=True)")

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
                        default=[31, 33,],
                        action='store',
                        nargs="+",
                        metavar="{31,33}",
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
                        default=sandy.get_seed(),
                        metavar="S31",
                        help="seed for random sampling of MF31 covariance "
                             "matrix (default = random)")

    parser.add_argument('--seed33',
                        type=int,
                        default=sandy.get_seed(),
                        metavar="S33",
                        help="seed for random sampling of MF33 covariance "
                             "matrix (default = random)")

    parser.add_argument('--seed34',
                        type=int,
                        default=sandy.get_seed(),
                        metavar="S34",
                        help="seed for random sampling of MF34 covariance "
                             "matrix (default = random)")

    parser.add_argument('--seed35',
                        type=int,
                        default=sandy.get_seed(),
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


    init = parser.parse_known_args(args=iargs)[0]
    if init.acer and not init.temperatures:
        parser.error("--acer requires --temperatures")
    return init


def multi_run(foo):
    """
    Decorator to handle keyword arguments for NJOY before running
    the executable.

    Examples
    --------
    Test that `minimal_processing` filters unwanted modules.
    >>> g = sandy.get_endf6_file("jeff_33", "xs", 10010).get_gendf(err=1, minimal_processing=True, temperature=300, dryrun=True)
    >>> assert "broadr" in g and "reconr" in g
    >>> assert "thermr" not in g and "purr" not in g and "heatr" not in g and "unresr" not in g and "gaspr" not in g

    Test `minimal_processing=False`.
    >>> g = sandy.get_endf6_file("jeff_33", "xs", 10010).get_gendf(err=1, temperature=300, dryrun=True)
    >>> assert "broadr" in g and "reconr" in g
    >>> assert "thermr" in g and "purr" in g and "heatr" in g and "gaspr" in g

    Check that for `temperature=0` the calculation stops after RECONR.
    >>> g = sandy.get_endf6_file("jeff_33", "xs", 10010).get_gendf(err=1, dryrun=True)
    >>> assert "reconr" in g
    >>> assert "broadr" not in g and "thermr" not in g and "purr" not in g and "heatr" not in g and "unresr" not in g and "gaspr" not in g
    """
    def inner(cli=None):
        """
        Parameters
        ----------
        """
        iargs = parse(cli)
        if os.path.isdir(iargs.file):
            path = iargs.file
            for file in os.listdir(path):
                iargs.file = os.path.join(path, file)
                foo(iargs)
        else:
            foo(iargs)
    return inner


@multi_run
def run(iargs):
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
    logging.info(f"processing file: '{iargs.file}'")
    
    err_pendf = 0.01
    err_ace = 0.01
    err_errorr = 0.1
    if iargs.debug:
        err_errorr = err_ace = err_pendf = 1

    endf6 = sandy.Endf6.from_file(iargs.file)

    # ERRORR KEYWORDS
    nubar = bool(31 in iargs.mf) and (31 in endf6.mf)
    xs = bool(33 in iargs.mf) and (33 in endf6.mf)
    mubar = False
    chi = False
    errorr_kws = dict(
        verbose=iargs.debug,
        err=err_errorr,
        xs=xs,
        nubar=nubar,
        chi=chi,
        mubar=mubar,
        groupr_kws=dict(nubar=nubar, chi=chi, mubar=mubar, ign=3),
        errorr_kws=dict(ign=3)
        )
    if iargs.mt33:
        errorr_kws["errorr33_kws"] = dict(mt=iargs.mt33)

    smp_kws = {}
    smp_kws["seed31"] = iargs.seed31
    smp_kws["seed33"] = iargs.seed33
    smp_kws["seed34"] = iargs.seed34
    smp_kws["seed35"] = iargs.seed35

    smps = endf6.get_perturbations(iargs.samples, njoy_kws=errorr_kws, smp_kws=smp_kws)


    if iargs.temperatures:
        temperature = iargs.temperatures[0] if hasattr(iargs.temperatures, "__len__") else iargs.temperatures
    else:
        temperature = 0

    # PENDF KEYWORDS
    pendf_kws = dict(
        verbose=iargs.debug,
        err=err_pendf,
        minimal_processing=iargs.debug,
        )

    # ACE KEYWORDS
    ace_kws = dict(
        verbose=iargs.debug,
        err=err_ace,
        minimal_processing=iargs.debug,
        temperature=temperature,
        purr=False,
        )

        
    endf6.apply_perturbations(
        smps,
        processes=iargs.processes,
        to_file=True,
        to_ace=iargs.acer,
        filename=iargs.outname,
        njoy_kws=pendf_kws,
        ace_kws=ace_kws,
        verbose=iargs.debug,
    )

    dt = time.time() - t0
    logging.info(f"Total running time: {dt:.2f} sec")


if __name__ == "__main__":
    run()
