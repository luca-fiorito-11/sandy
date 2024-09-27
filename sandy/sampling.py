import os
import time
import logging
import argparse
import subprocess as sp

from .endf6 import Endf6
from .tools import is_valid_file
from .utils import get_seed
from .samples import Samples
from . import __version__


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

    parser.add_argument('--fycov',
                        default=False,
                        action="store_true",
                        help="use CEA covariance data for U-235 and Pu-239 thermal fission yields")

    parser.add_argument("--from_perturbations",
                        default=False,
                        nargs=3,
                        help="resume the sampling pipeline reading the "
                             "perturbation coefficients from file\n"
                             "The three entires are:\n"
                             " - directory where perturbation coefficients are stored\n"
                             " - first perturbation coefficient to consider\n"
                             " - last perturbation coefficient to consider")

    parser.add_argument('--loglevel',
                        type=str,
                        default="info",
                        action='store',
                        metavar="{debug, info, warning, error, critical}",
                        help="Set the logger verbosity level.")

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
                        help="draw samples only from the selected MT sections for MF33"
                             "(default = keep all)")

    parser.add_argument('--njoy',
                        type=lambda x: is_valid_file(parser, x),
                        default=None,
                        help="NJOY executable "
                             "(default search PATH, and env variable NJOY)")

    parser.add_argument("--only_perturbations",
                        default=False,
                        action="store_true",
                        help="stop the sampling pipeline after the creation "
                             "of perturbation the coefficients")

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
                        default=get_seed(),
                        metavar="S31",
                        help="seed for random sampling of MF31 covariance "
                             "matrix (default = random)")

    parser.add_argument('--seed33',
                        type=int,
                        default=get_seed(),
                        metavar="S33",
                        help="seed for random sampling of MF33 covariance "
                             "matrix (default = random)")

    parser.add_argument('--seed34',
                        type=int,
                        default=get_seed(),
                        metavar="S34",
                        help="seed for random sampling of MF34 covariance "
                             "matrix (default = random)")

    parser.add_argument('--seed35',
                        type=int,
                        default=get_seed(),
                        metavar="S35",
                        help="seed for random sampling of MF35 covariance "
                             "matrix (default = random)")

    parser.add_argument('--supressnjoy',
                        default=False,
                        action="store_true",
                        help="Supress NJOY ouputs.")

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
                        version='%(prog)s {}'.format(__version__),
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
    
    >>> import sandy, filecmp
    >>> import pandas as pd
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

    Retrieve ENDF-6 tape and write it to file.

    >>> sandy.get_endf6_file("jeff_33", "xs", 10010).to_file("H1.jeff33")

    Produce perturbed ACE file.

    >>> cli = "H1.jeff33 --acer True --samples 2 --processes 2 --temperatures 900 --seed33 5"
    >>> sandy.sampling.run(cli.split())

    Check if ACE and XSDIR files have the right content.

    >>> assert "1001.09c" in open("1001_0.09c").read()
    >>> assert "1001.09c" in open("1001_0.09c.xsd").read()
    >>> assert "1001.09c" in open("1001_1.09c").read()
    >>> assert "1001.09c" in open("1001_1.09c.xsd").read()
    >>> assert not filecmp.cmp("1001_0.09c", "1001_1.09c", shallow=False)

    Run the same on a single process.

    >>> cli = "H1.jeff33 --acer True --samples 2 --processes 2 --temperatures 900 --seed33 5 --outname={ZAM}_{SMP}_SP"
    >>> sandy.sampling.run(cli.split())

    The identical seed ensures consistent results with the previous run.

    >>> assert filecmp.cmp("1001_0.09c", "10010_0_SP.09c")
    >>> assert filecmp.cmp("1001_1.09c", "10010_1_SP.09c")
    >>> assert filecmp.cmp("1001_0.09c.xsd", "10010_0_SP.09c.xsd")
    >>> assert filecmp.cmp("1001_1.09c.xsd", "10010_1_SP.09c.xsd")

    Produce perturbed ENDF6 and PENDF files.

    >>> cli = "H1.jeff33 --samples 2 --processes 2 --outname=H1_{MAT}_{SMP} --mt 102"
    >>> sandy.sampling.run(cli.split())
    >>> assert os.path.getsize("H1_125_0.pendf") > 0 and os.path.getsize("H1_125_1.pendf") > 0

    >>> assert filecmp.cmp("H1_125_0.endf6", "H1_125_1.endf6")
    >>> assert filecmp.cmp("H1_125_0.endf6", "H1.jeff33")
    
    Let's see how the sampling process can be interrupted fater.
    Produce random ENDF-6 and PENDF files for Pu-241 with the standard procedure.

    >>> file = "942410.jeff33"
    >>> sandy.get_endf6_file("jeff_33", "xs", 942410).to_file(file)
    >>> cl = f"{file}" + " --samples 2 -O {SMP}-{ZAM} --seed33 1 --seed31 1 --mt33 2"
    >>> sandy.sampling.run(cl.split())

    Now, let's interrupt the process after that the perturbations are
    created (reproducible with fixed seed).

    >>> smps = sandy.sampling.run((cl + " --only_perturbations").split())

    We can read these perturbation coefficients without the need of regenerating them.

    >>> cl = f"{file} --from_perturbations {os.getcwd()} 1 1 --only_perturbations"
    >>> smps2 = sandy.sampling.run(cl.split())
    >>> assert smps2[33].data.shape[1] == smps2[31].data.shape[1] == 1
    >>> assert smps[33].data.reset_index().MT.unique() == 2
    >>> assert smps[31].data.reset_index().MT.unique().size == 3
    >>> pd.testing.assert_frame_equal(smps2[33].data, smps[33].data[[1]])
    >>> pd.testing.assert_frame_equal(smps2[31].data, smps[31].data[[1]])

    Using the perturbation coefficients from the excel files we generate the
    same random files of the standard pipeline.

    >>> cl = f"{file}" + " -O new_{SMP}-{ZAM} " + f"--from_perturbations {os.getcwd()} 1 1"
    >>> sandy.sampling.run(cl.split())
    >>> assert filecmp.cmp("new_1-942410.endf6", "1-942410.endf6")
    >>> assert filecmp.cmp("new_1-942410.pendf", "1-942410.pendf")

    If no perturbation file exist, the calculation stops.

    >>> file = "741840.jeff33"
    >>> sandy.get_endf6_file("jeff_33", "xs", 741840).to_file(file)
    >>> cl = f"{file} --from_perturbations {os.getcwd()} 1 1 --only_perturbations"
    >>> assert not sandy.sampling.run(cl.split())
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
            return foo(iargs)
    return inner


def running_time(foo):
    """
    Decorator to handle keyword arguments for NJOY before running
    the executable.
    """
    def inner(*args, **kwargs):
        t0 = time.time()
        out = foo(*args, **kwargs)
        dt = time.time() - t0
        logging.info(f"Total running time: {dt:.2f} sec")
        return out
    return inner


@running_time
@multi_run
def run(iargs):
    """
    Run `sandy` sampling sequence.

    Parameters
    ----------
    iargs : `list` of `str`
        arguments of the command line.
        For example, the two following options are identical:
        
            - in python
                .. code-block:: python

                cli = "H1.jeff33 --acer True --samples 1 --processes 1 --temperatures 900 --seed33 5"
                sandy.sampling.run(cli.split())
                
            - from command line
                .. code-block:: sh
    
                    H1.jeff33 --acer True --samples 1 --processes 1 --temperatures 900 --seed33 5

    Returns
    -------
    None.
    
    Examples
    --------
    Default use case for decay data sampling.

    >>> import sandy
    >>> sandy.get_endf6_file("jeff_33", "decay", [10010, 10040, 270600]).to_file("AAA.txt")
    >>> sandy.sampling.run("AAA.txt --samples 3 --processes 1".split())
    >>> assert {'decay_data_0', 'decay_data_1', 'decay_data_2'}.issubset(set(glob.glob("decay_data*")))

    Default use case for fission yield sampling.

    >>> import sandy
    >>> sandy.get_endf6_file("jeff_33", "nfpy", [922350, 922380]).to_file("AAA.txt")
    >>> sandy.sampling.run("AAA.txt --samples 3 --processes 1".split())
    >>> assert {'fy_0', 'fy_1', 'fy_2'}.issubset(set(glob.glob("fy*")))
    """

    loglevels = {
        "debug": logging.DEBUG,
        "info": logging.INFO,
        "warning": logging.WARNING,
        "error": logging.ERROR,
        "critical": logging.CRITICAL,
    }
    logging.getLogger().setLevel(loglevels[iargs.loglevel])
    logging.info(f"processing file: '{iargs.file}'")

    err_pendf = 0.01
    err_ace = 0.01
    err_errorr = 0.1
    if iargs.debug:
        err_errorr = err_ace = err_pendf = 1

    endf6 = Endf6.from_file(iargs.file)
    
    if 457 in endf6.mt:
        # very dirty way to add decay data sampling from command line interface
        smps = endf6.get_perturbations(iargs.samples)
        endf6.apply_perturbations(
            smps,
            processes=iargs.processes,
            to_file=True,
            verbose=iargs.debug,
        )
        return

    if 454 in endf6.mt:
        # very dirty way to add decay data sampling from command line interface
        # Let's use CEA covariance data by default
        covariance = "cea" if iargs.fycov else None
        smps = endf6.get_perturbations(iargs.samples, covariance=covariance)
        endf6.apply_perturbations(
            smps,
            processes=iargs.processes,
            covariance=covariance,
            to_file=True,
            verbose=iargs.debug,
        )
        return       

    njoy_output = sp.DEVNULL if iargs.supressnjoy else None

    # ERRORR KEYWORDS
    nubar = bool(31 in iargs.mf) and (31 in endf6.mf)
    xs = bool(33 in iargs.mf) and (33 in endf6.mf or 32 in endf6.mf)  # this handles together MF32 and MF33
    mubar = False
    chi = False
    errorr_kws = dict(
        verbose=iargs.debug,
        err=err_errorr,
        xs=xs,
        nubar=nubar,
        chi=chi,
        mubar=mubar,
        groupr_kws=dict(nubar=nubar, chi=chi, mubar=mubar, ign=2),
        errorr_kws=dict(ign=2),
        njoy_output=njoy_output
        )
    if iargs.mt33:
        errorr_kws["errorr33_kws"] = dict(mt=iargs.mt33)

    smp_kws = {}
    smp_kws["seed31"] = iargs.seed31
    smp_kws["seed33"] = iargs.seed33
    smp_kws["seed34"] = iargs.seed34
    smp_kws["seed35"] = iargs.seed35


    if iargs.from_perturbations:
        smps = {}
        ID = endf6.get_id()
        file = os.path.join(iargs.from_perturbations[0], "PERT_{}_MF{}.xlsx")
        beg = int(iargs.from_perturbations[1])
        end = int(iargs.from_perturbations[2])

        for mf in [31, 33, 34, 35]:
            xls = file.format(ID, mf)
            if os.path.isfile(xls):
                smps[mf] = Samples.from_excel(xls, beg=beg, end=end)
        if not smps:
            logging.warning(f"No perturbation file was found for {ID}")

    else:
        smps = endf6.get_perturbations(iargs.samples, njoy_kws=errorr_kws, smp_kws=smp_kws)

    if iargs.only_perturbations:
        return smps


    if iargs.temperatures:
        temperature = iargs.temperatures[0] if hasattr(iargs.temperatures, "__len__") else iargs.temperatures
    else:
        temperature = 0

    # PENDF KEYWORDS
    pendf_kws = dict(
        verbose=iargs.debug,
        err=err_pendf,
        minimal_processing=iargs.debug,
        njoy_output=njoy_output,
        )

    # ACE KEYWORDS
    ace_kws = dict(
        verbose=iargs.debug,
        err=err_ace,
        minimal_processing=iargs.debug,
        temperature=temperature,
        purr=False,
        njoy_output=njoy_output,
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

    return


if __name__ == "__main__":
    run()
