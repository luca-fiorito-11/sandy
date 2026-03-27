import os


__author__ = "Luca Fiorito"


def parse(iargs=None):
    """
    Parse command line arguments for sampling option.

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
    import argparse
    from .utils import is_valid_file
    from . import __version__

    description = (
        "Produce perturbed files containing sampled parameters that "
        "represent the information stored in the evaluated nuclear "
        "data covariances."
        )
    parser = argparse.ArgumentParser(
                        prog="sandy",
                        description=description,
                        formatter_class=argparse.RawTextHelpFormatter,
                        )

    parser.add_argument(
        'file',
        help="ENDF-6 file",
        )
    parser.add_argument(
        '--acer',
        default=False,
        action="store_true",
        help=(
            "Process each perturbed file into ACE format "
            "(default = False)\n"
            "(--temperatures is required)"
            )
        )
    parser.add_argument(
        '--fycov',
        default=False,
        action="store_true",
        help=(
            "Use covariance data for U-233, U-235, Pu-239, Pu-241 thermal fission yields"
            "(only for JEFF-4.0 fission yields)"
            )
        )
    parser.add_argument(
        "--from_perturbations",
        default=False,
        nargs=3,
        help=(
            "Resume the sampling pipeline reading the "
            "perturbation coefficients from file\n"
            "The three entries are:\n"
            " - directory where perturbation coefficients are stored\n"
            " - first perturbation coefficient to consider\n"
            " - last perturbation coefficient to consider"
            )
        )
    parser.add_argument(
        "--cov_energy_grid",
        default="csewg239",
        choices=["csewg239", "lanl30", "epri69", "ecco33"],
        help=(
            "Energy grid to process covariance matrix\n"
            "Allowed entries are:\n"
            " - csewg239 (default)\n"
            " - lanl30\n"
            " - epri69\n"
            " - ecco33"
            )
        )
    parser.add_argument(
        '--loglevel',
        type=str,
        default="info",
        action='store',
        metavar="{debug, info, warning, error, critical}",
        help="Set the logger verbosity level."
        )
    parser.add_argument(
        '--mf',
        type=int,
        default=[31, 33, 35],
        action='store',
        nargs="+",
        metavar="{31,33,35}",
        help=(
            "Draw samples only from the selected MF sections "
            "(default = keep all)"
            )
        )
    parser.add_argument(
        '--mt33',
        type=int,
        default=None,
        action='store',
        nargs="+",
        metavar="{1,..,999}",
        help=(
            "Draw samples only from the selected MT sections for MF33"
            "(default = keep all)"
            )
        )
    parser.add_argument(
        '--njoy',
        type=lambda x: is_valid_file(parser, x),
        default=None,
        help=(
            "NJOY executable "
            "(default search PATH, and env variable NJOY)"
            )
        )
    parser.add_argument(
        '--no-verbose',
        default=False,
        action="store_true",
        help="Disable verbose output",
        )

    parser.add_argument(
        "--only_perturbations",
        default=False,
        action="store_true",
        help=(
            "Stop the sampling pipeline after the creation "
            "of perturbation coefficients"
            )
        )
    parser.add_argument(
        '--processes', '-N',
        type=int,
        default=1,
        help="Number of worker processes (default = 1)"
        )
    parser.add_argument(
        '--samples', '-S',
        type=int,
        default=200,
        help="Number of samples (default = 200)"
        )
    parser.add_argument(
        '--seed31',
        type=int,
        default=None,
        metavar="S31",
        help=(
            "Seed for random sampling of MF31 covariance "
            "matrix (default = random)"
            )
        )
    parser.add_argument(
        '--seed33',
        type=int,
        default=None,
        metavar="S33",
        help=(
            "Seed for random sampling of MF33 covariance "
            "matrix (default = random)"
            )
        )
    parser.add_argument(
        '--seed34',
        type=int,
        default=None,
        metavar="S34",
        help=(
            "Seed for random sampling of MF34 covariance "
            "matrix (default = random)"
            )
        )
    parser.add_argument(
        '--seed35',
        type=int,
        default=None,
        metavar="S35",
        help=(
            "Seed for random sampling of MF35 covariance "
            "matrix (default = random)"
            )
        )
    parser.add_argument(
        "--show-njoy",
        default=False,
        action="store_true",
        help="Show NJOY output (suppressed by default)."
    )
    parser.add_argument(
        '--tqdm',
        default=False,
        action="store_true",
        help="Enable tqdm progress bars (off by default; not recommended in non-interactive terminals)",
        )
    parser.add_argument(
        '--temperatures', '-T',
        default=None,
        type=float,
        action='store',
        nargs="+",
        metavar="T",
        help=(
            "For each perturbed file, produce ACE files at "
            "given temperatures"
            )
        )
    parser.add_argument(
        "--version", "-v",
        action='version',
        version='%(prog)s {}'.format(__version__),
        help="SANDY's version."
        )

    init = parser.parse_known_args(args=iargs)[0]

    if init.acer and not init.temperatures:
        parser.error("--acer requires --temperatures")

    return init



def run(cli=None):
    import argparse

    iargs = parse(cli)

    if os.path.isdir(iargs.file):

        results = []
        for name in os.listdir(iargs.file):
            path = os.path.join(iargs.file, name)
            if not os.path.isfile(path):
                continue
            # build a fresh namespace for each file
            iargs_per_file = argparse.Namespace(**vars(iargs).copy())
            iargs_per_file.file = path
            results.append(_process_one_file(iargs_per_file))
        return results

    else:
        return _process_one_file(iargs)



def _process_one_file(
        iargs,
        ):
    """
    Run the end-to-end sampling pipeline for a single ENDF-6 file, including:
    - reading the input file,
    - generating or loading perturbation coefficients,
    - optionally returning only perturbations,
    - and applying perturbations to produce ENDF6/PENDF/ACE outputs.

    The function supports three main data paths based on the content of the input
    ENDF-6 file:
      * Decay data (MF=8/MT=457) → decay sampling pipeline
      * Fission yields (MF=8/MT=454) → fission-yield sampling pipeline
      * Cross sections and related (MF=31/33/35) → ERRORR-driven pipeline with NJOY

    Parameters
    ----------
    iargs : argparse.Namespace
        Parsed command-line arguments with (at least) the
        attributes defined in :func:`~sandy.sampling.parse`
            file : str
                Path to the input ENDF-6 file.
            samples : int
                Number of samples to generate.
            processes : int
                Number of parallel processes to use when applying perturbations.
            no_verbose : bool
                If True, suppress verbose logging.
            loglevel : {"debug","info","warning","error","critical"}
                Logging level (string).
            tqdm : bool
                If True, show progress bars while applying perturbations.
            show_njoy : bool
                If True, print NJOY output to stdout.
            acer : bool
                If True, produce ACE (and XSDIR) files.
            temperatures : float | list[float] | None
                ACE temperatures (first one is used when `acer` is True).
            mf : list[int]
                Requested covariance sections to perturb (e.g., [31, 33, 35]).
            mt33 : int | list[int] | None
                If provided, restrict MF=33 perturbations to specific MT(s).
            cov_energy_grid : {"csewg239", "lanl30", "epri69", "ecco33"}
                Discrete-groups grid to be used by ERRORR/GROUPR (maps to NJOY IGN).
            seed31, seed33, seed34, seed35 : int | None
                Random seeds for independent reproducibility of different MFs.
            fycov : str | None
                Optional covariance selection for FY pipeline (used when MF8/MT454).
            only_perturbations : bool
                If True, return perturbation coefficients only and produce no files.
            from_perturbations : tuple[str, int, int] | None
                (base_dir, beg, end). If provided, load perturbations from
                PERT_{MAT}_MF{31|33|35}.xlsx in `base_dir` for sample IDs in
                [beg..end] rather than generating them.
                - `base_dir` must exist and be a directory
                - `beg` and `end` must be integers with 0 <= beg <= end

    Returns
    -------
    smps : dict[int | str, Samples] or None
        - If `iargs.only_perturbations` is True:
            Returns a dictionary of Samples (perturbation coefficients),
            keyed by MF integers (31, 33, 35) for the cross-section path, or
            by string keys ("IFY") for the FY shortcut, or by string keys
            ("HL", "DE", "BR") for the RDD shortcut, depending on the pipeline 
            branch executed. The dictionary may be empty if no available 
            perturbations were found or produced.
        - Otherwise:
            Returns None after writing the corresponding output files on disk.

    Raises
    ------
    Exception
        If the input path does not exist or is not a file.
        If `--samples` or `--processes` are not positive integers.
        If `--from_perturbations` is provided with invalid directory or range.

    Notes
    -----
    Pipeline selection:
        * Decay shortcut (MF=8/MT=457)
            - Generates decay perturbations, applies them, and writes files
              named like: `decay_data_{i}`.
        * Fission-yield shortcut (MF=8/MT=454)
            - Generates FY perturbations, applies them, and writes files named
              like: `fy_{i}`.
        * ERRORR-based pipeline (MF=31/33/35)
            - If `from_perturbations` is provided, loads perturbations from
              Excel files `PERT_{MAT}_MF{31|33|35}.xlsx`.
            - Otherwise, computes perturbations with NJOY (GROUPR/ERRORR).
            - When `acer` is True, ACE/XSDIR files are generated at the first
              provided temperature (if `temperatures` is a list).

    Examples
    --------
    Basic decay data sampling: check that perturbed ENDF‑6 files are created
    and that they differ from each other.

    >>> import sandy, glob, os, filecmp
    >>> # Create a decay ENDF file
    >>> sandy.get_endf6_file("jeff_33", "decay", [10010, 10040, 270600], local=True).to_file("AAA.txt")
    
    First, remove any old output files.
    
    >>> for f in glob.glob("decay_data_*"):
    ...     try: os.remove(f)
    ...     except FileNotFoundError: pass

    Run sampling.

    >>> cli = "AAA.txt --samples 3 --processes 1 --no-verbose"
    >>> sandy.sampling.run(cli.split())

    Check that expected files were created and are not empty.
    
    >>> outfiles = sorted(glob.glob("decay_data_*"))
    >>> assert len(outfiles) == 3
    >>> assert {'decay_data_0', 'decay_data_1', 'decay_data_2'}.issubset(outfiles)
    >>> assert all(os.path.getsize(f) > 0 for f in outfiles)

    Two samples must differ.
    
    >>> assert not filecmp.cmp(outfiles[0], outfiles[1], shallow=False)



    Basic fission yield sampling: check that perturbed ENDF‑6 files are created
    and that they differ from each other.

    >>> import sandy, glob, os, filecmp
    >>> # Create a fy ENDF file
    >>> sandy.get_endf6_file("jeff_33", "nfpy", [922350, 922380], local=True).to_file("AAA.txt")

    First, remove any old output files.
    
    >>> for f in glob.glob("decay_data_*"):
    ...     try: os.remove(f)
    ...     except FileNotFoundError: pass
 
    Run sampling.

    >>> cli = "AAA.txt --samples 3 --processes 1 --no-verbose"
    >>> sandy.sampling.run(cli.split())

    Check that expected files were created and are not empty.
    
    >>> outfiles = sorted(glob.glob("fy_*"))
    >>> assert len(outfiles) == 3
    >>> assert {'fy_0', 'fy_1', 'fy_2'}.issubset(outfiles)
    >>> assert all(os.path.getsize(f) > 0 for f in outfiles)

    Two samples must differ.
    
    >>> assert not filecmp.cmp(outfiles[0], outfiles[1], shallow=False)

    Check option ``--fycov``, which only works for U-235, Pu-239, U-233 and
    Pu-241 of JEFF-4.0.

    >>> import sandy, numpy as np
    >>> sandy.get_endf6_file("jeff_40", "nfpy", 922350, local=True).to_file("92-U-235.fy_jeff40")
    >>> cli = "92-U-235.fy_jeff40 --samples 100 --only_perturbations --no-verbose --fycov"
    >>> smps_corr = sandy.sampling.run(cli.split())

    When using ``--fycov`` for U-235, the correlation between nuclides ``zap=521350``
    and ``zap=531350`` should approach ``-0.906028`` (depending on the sample size).

    >>> corr = smps_corr["IFY"].data.query("ZAP in [521350, 531350] & E==0.0253").T.corr()
    >>> assert np.abs(corr.iloc[0, 1]) > 0.8

    Without ``--fycov`` the same correlation remains zero.

    >>> cli = "92-U-235.fy_jeff40 --samples 100 --only_perturbations --no-verbose"
    >>> smps_nocorr = sandy.sampling.run(cli.split())
    >>> corr = smps_nocorr["IFY"].data.query("ZAP in [521350, 531350] & E==0.0253").T.corr()
    >>> assert np.abs(corr.iloc[0, 1]) < 0.3



    Default use case for xs sampling.

    >>> import sandy, glob, os, filecmp
    >>> # Create a ENDF file
    >>> sandy.get_endf6_file("jeff_33", "xs", 10010, local=True).to_file("H1.jeff33")
    
    First, remove any old output files.
    
    >>> for f in glob.glob("1001_*"):
    ...     try: os.remove(f)
    ...     except FileNotFoundError: pass

    >>> cli = "H1.jeff33 --acer True --samples 2 --processes 2 --temperatures 900 --seed33 5 --no-verbose"
    >>> sandy.sampling.run(cli.split())

    Check if ACE and XSDIR files have the right content.

    >>> for file in ["1001_0.09c", "1001_0.09c.xsd", "1001_1.09c", "1001_1.09c.xsd"]:
    ...    assert "1001.09c" in open(file).read()

    Two samples must differ.
    
    >>> assert not filecmp.cmp("1001_0.09c", "1001_1.09c", shallow=False)

    Run the same on a single process. But first move files.

    >>> for f in ["1001_0_MP.09c", "1001_0_MP.09c.xsd", "1001_1_MP.09c", "1001_1_MP.09c.xsd"]:
    ...     try: os.remove(f)
    ...     except FileNotFoundError: pass
    >>> for file in ["1001_0.09c", "1001_0.09c.xsd", "1001_0.endf6", "1001_0.pendf"]:
    ...     os.rename(file, f"{file}_MP")
    >>> for file in ["1001_1.09c", "1001_1.09c.xsd", "1001_1.endf6", "1001_1.pendf"]:
    ...     os.rename(file, f"{file}_MP")

    Run sampling.

    >>> cli = "H1.jeff33 --acer True --samples 2 --processes 1 --temperatures 900 --seed33 5 --no-verbose"
    >>> sandy.sampling.run(cli.split())

    The identical seed ensures consistent results with the previous run.

    >>> for file in ["1001_0.09c", "1001_0.09c.xsd", "1001_0.endf6", "1001_0.pendf"]:
    ...    assert filecmp.cmp(file, f"{file}_MP", shallow=False)
    >>> for file in ["1001_1.09c", "1001_1.09c.xsd", "1001_1.endf6", "1001_1.pendf"]:
    ...    assert filecmp.cmp(file, f"{file}_MP", shallow=False)
    
    Produce perturbed ENDF6 and PENDF files. But first remove the files again.

    >>> for f in glob.glob("1001_*"):
    ...     try: os.remove(f)
    ...     except FileNotFoundError: pass

    Run sampling.

    >>> cli = "H1.jeff33 --samples 2 --processes 1 --mt 102 --no-verbose"
    >>> sandy.sampling.run(cli.split())
    
    In this case no ace file was created

    >>> assert not glob.glob("1001*c")
    >>> assert not glob.glob("1001*xsd")
    
    PENDF files were created and differ from each other.

    >>> assert os.path.getsize("1001_0.pendf") > 0 and os.path.getsize("1001_1.pendf") > 0
    >>> assert not filecmp.cmp("1001_0.pendf", "1001_1.pendf", shallow=False)
    
    Since only MT=102 is requested, only the capture and total (reconstructed
    cross section) are changed.

    >>> from sandy.utils import assert_all_diffs_match_any
    >>> patterns = [r".{66} 125 3  1[ 0-9]{5}$",
    ...             r".{66} 125 3102[ 0-9]{5}$"]
    >>> assert_all_diffs_match_any("1001_0.pendf", "1001_1.pendf", patterns=patterns)

    Since the covariance data is given only for MF33, the endf6 outputs are unperturbed.
    They are also identical to the original ENDF-6 file.
    
    >>> assert filecmp.cmp("1001_0.endf6", "1001_1.endf6", shallow=False)
    >>> assert filecmp.cmp("1001_0.endf6", "H1.jeff33", shallow=False)
    
    Let's see how the sampling process can be interrupted.
    Produce random ENDF-6 and PENDF files for Pu-241 with the standard procedure.

    >>> import sandy, glob, os
    >>> sandy.get_endf6_file("jeff_33", "xs", 942410, local=True).to_file("942410.jeff33")
    >>> cli = "942410.jeff33 --samples 2 --seed33 1 --seed31 1 --seed35 1 --mt33 2 --no-verbose"

    Run sampling. But first remove files

    >>> # Remove files
    >>> for f in glob.glob("942410_*"):
    ...     try: os.remove(f)
    ...     except FileNotFoundError: pass
    >>> # Run sampling
    >>> sandy.sampling.run(cli.split())

    Let's rename them. We'll use them later on for comparison.

    >>> outfiles = glob.glob("942410_*")
    >>> for f in outfiles:
    ...     os.rename(f, f"{f}_1step")

    Now, let's interrupt the process after the perturbations are created (reproducible with fixed seed).

    >>> # First delete the perturbation file if it exists.
    >>> try: os.remove("PERT_94241_MF31.xlsx")
    ... except FileNotFoundError: pass

    If run from an interactive terminal, it returns the samples.

    >>> cli = cli + " --only_perturbations"
    >>> smps = sandy.sampling.run(cli.split())

    We can read these perturbation coefficients without the need of regenerating them.

    First sample 0.

    >>> import shlex, numpy as np
    >>> cli = f"942410.jeff33 --from_perturbations '{os.getcwd()}' 0 0 --only_perturbations --no-verbose"
    >>> smps0 = sandy.sampling.run(shlex.split(cli))
    >>> assert all([np.allclose(v.data[[0]], smps0[k].data) for k, v in smps.items()])
    >>> assert all([v.data[[0]].index.equals(smps0[k].data.index) for k, v in smps.items()])

    Then sample 1.

    >>> import shlex, numpy as np
    >>> cli = f"942410.jeff33 --from_perturbations '{os.getcwd()}' 1 1 --only_perturbations --no-verbose"
    >>> smps1 = sandy.sampling.run(shlex.split(cli))
    >>> assert all([np.allclose(v.data[[1]], smps1[k].data) for k, v in smps.items()])
    >>> assert all([v.data[[1]].index.equals(smps1[k].data.index) for k, v in smps.items()])

    Then together.

    >>> import shlex, numpy as np
    >>> cli = f"942410.jeff33 --from_perturbations '{os.getcwd()}' 1 1 --only_perturbations --no-verbose"
    >>> smps1 = sandy.sampling.run(shlex.split(cli))
    >>> assert all([np.allclose(v.data[[1]], smps1[k].data) for k, v in smps.items()])
    >>> assert all([v.data[[1]].index.equals(smps1[k].data.index) for k, v in smps.items()])

    An excel file of perturbations was produced.
    
    >>> assert os.path.exists("PERT_94241_MF31.xlsx")

    Using the perturbation coefficients from the excel files we generate the
    same random files of the standard pipeline.

    >>> cli = f"942410.jeff33 --from_perturbations '{os.getcwd()}' 1 1 --no-verbose"
    >>> sandy.sampling.run(shlex.split(cli))

    This produces only the first sample (id "_0").
    
    >>> for file in outfiles:
    ...    if "_1" in file:
    ...       assert not os.path.isfile(file)
    ...    else:
    ...       assert filecmp.cmp(file, f"{file}_1step", shallow=False)
    

    If no perturbation file exists, the calculation stops.

    >>> import os, shlex, sandy
    >>> file = "741840.jeff33"
    >>> sandy.get_endf6_file("jeff_33", "xs", 741840, local=True).to_file("741840.jeff33")
    >>> cli = f"741840.jeff33 --from_perturbations '{os.getcwd()}' 1 1 --only_perturbations --no-verbose"
    >>> assert not sandy.sampling.run(shlex.split(cli))
    
    Now let's check that MF5 gets perturbed.

    >>> sandy.get_endf6_file("jeff_33", "xs", 942410, local=True).to_file("942410.jeff33")
    >>> cli = "942410.jeff33 --samples 2 --mf 35 --no-verbose"
    >>> sandy.sampling.run(cli.split())

    MF5 / MT18 should be the only modified part in the file.

    >>> from sandy.utils import assert_all_diffs_match_any
    >>> patterns = [r".{66}9443 5 18[ 0-9]{5}$"]
    >>> assert_all_diffs_match_any("94241_0.endf6", "94241_1.endf6", patterns=patterns)
    >>> assert filecmp.cmp("94241_0.pendf", "94241_0.pendf", shallow=False)

    Now let's check that MF1 gets perturbed.

    >>> sandy.get_endf6_file("jeff_33", "xs", 942410, local=True).to_file("942410.jeff33")
    >>> cli = "942410.jeff33 --samples 2 --mf 31 --no-verbose"
    >>> sandy.sampling.run(cli.split())

    MF1 / MT452 / 455 / 456 should be the only modified part in the file.

    >>> from sandy.utils import assert_all_diffs_match_any
    >>> patterns = [r".{66}9443 145[256][ 0-9]{5}$"]
    >>> assert_all_diffs_match_any("94241_0.endf6", "94241_1.endf6", patterns=patterns)
    >>> assert filecmp.cmp("94241_0.pendf", "94241_0.pendf", shallow=False)

    """
    # ---- IMPORT
    import pprint
    import time
    import logging
    
    from .endf6 import Endf6
    from .samples import Samples
    from .utils import log
    from ._perturbation_base import log_stage

    t0 = time.time()

    # ---- HELPER: print config
    def _recap_config(args):
        d = vars(args).copy()  # SAFE COPY
        if d["mt33"] is None:
            d["mt33"] = "all"

        cfg = pprint.pformat(d, indent=2, sort_dicts=True)
        msg = (
            f"\n"
            f"========= SAMPLING CONFIGURATION =========\n"
            f"{cfg}\n"
            f"==========================================\n"
        )
        log(msg)

    # ---- HELPER: compute temperature for ACE
    def _compute_temperature(args):
        if args.temperatures:
            return args.temperatures[0] if hasattr(args.temperatures, "__len__") else float(args.temperatures)
        return 0.0

    # ---- HELPER: build sampling keyword dict
    def _build_smp_kws(args):
        # Seeds are already type=int by argparse; just log if something looks odd
        return dict(
            seed31=args.seed31,
            seed33=args.seed33,
            seed34=args.seed34,
            seed35=args.seed35,
            )

    # ---------------------------------------
    # ---- LOGGING SETUP
    # ---------------------------------------
    method = "sampling"
    loglevels = {
        "debug": logging.DEBUG,
        "info": logging.INFO,
        "warning": logging.WARNING,
        "error": logging.ERROR,
        "critical": logging.CRITICAL,
    }
    logging.getLogger().setLevel(loglevels[iargs.loglevel])

    # Verbosity flag
    verbose = not iargs.no_verbose
    
    # Print config recap (if verbose)
    if verbose:
        _recap_config(iargs)

    # ---------------------------------------
    # ---- BASIC ARG VALIDATION
    # ---------------------------------------
    # Let the decorator handle directory traversal; here we expect a file

    if not os.path.exists(iargs.file):
        msg = f"input path does not exist: '{iargs.file}'"
        raise Exception(msg)

    if not os.path.isfile(iargs.file):
        msg = f"expected a file path, got a directory: '{iargs.file}'",
        raise Exception(msg)

    if not (isinstance(iargs.samples, int) and iargs.samples > 0):
        msg = f"'--samples' must be a positive integer (got {iargs.samples})"
        raise Exception(msg)

    if not (isinstance(iargs.processes, int) and iargs.processes > 0):
        msg = f"'--processes' must be a positive integer (got {iargs.processes})"
        raise Exception(msg)

    msg = f"processing file: '{iargs.file}'"
    log_stage(log, method, None, msg, verbose=verbose)

    # ---- SETUP NJOY tolerances
    # later make these configurable
    err_pendf = 0.01
    err_ace = 0.01
    err_errorr = 0.1


    # ---------------------------------------------
    # ---- READ ENDF6 
    # ---------------------------------------------
    endf6 = Endf6.from_file(iargs.file)
    len_mat = len(endf6.mat)
    len_mf = len(endf6.mf)
    len_mt = len(endf6.mt)

    # Diagnose content
    msg = f"ENDF-6 content: {len_mat} MAT, {len_mf} MF, {len_mt} MT"
    log_stage(log, method, None, msg, verbose=verbose)
    
    # ----  DECAY SHORTCUT
    if 457 in endf6.mt:
        # very dirty way to add decay data sampling from command line interface
        msg = "detected MF=8/MT=457 (decay): using decay-data path"
        log_stage(log, method, None, msg, verbose=verbose)
        
        # logging is already in the method
        smps = endf6.get_perturbations(
            iargs.samples,
            verbose=verbose,
            write=True,
            )

        # Only provide perturbations if requested
        if iargs.only_perturbations:
            msg = "returning perturbations only"
            if smps == {}:
                msg += " (empty dict)"
            
            log_stage(log, method, None, msg, verbose=verbose)

            dt = time.time() - t0
            msg = f"total running time: {dt:.2f} sec"
            log_stage(log, method, None, msg, verbose=verbose)
            return smps

        # ----------------------------------------
        # ---- CHECK IF PERTURBATIONS ARE PRODUCED
        # ----------------------------------------
        if smps == {}:
            msg = "no perturbation was produced"
            log_stage(log, method, None, msg, verbose=verbose)

            dt = time.time() - t0
            msg = f"total running time: {dt:.2f} sec"
            log_stage(log, method, None, msg, verbose=verbose)
            return

        # logging is already in the method
        endf6.apply_perturbations(
            smps,
            processes=iargs.processes,
            to_file=True,
            verbose=verbose,
            enable_tqdm=iargs.tqdm,
        )

        msg = "RDD perturbations applied successfully"
        log_stage(log, method, None, msg, verbose=verbose)

        dt = time.time() - t0
        msg = f"total running time: {dt:.2f} sec"
        log_stage(log, method, None, msg, verbose=verbose)

        return

    # ----  FY SHORTCUT
    if 454 in endf6.mt:
        # very dirty way to add decay data sampling from command line interface
        msg = "detected MF=8/MT=454 (decay): using fission-yield path"
        log_stage(log, method, None, msg, verbose=verbose)

        # logging is already in the method
        smps = endf6.get_perturbations(
            iargs.samples,
            covariance=iargs.fycov,
            verbose=verbose,
            write=True,
            )

        # Only provide perturbations if requested
        if iargs.only_perturbations:
            msg = "returning perturbations only"
            if smps == {}:
                msg += " (empty dict)"
            
            log_stage(log, method, None, msg, verbose=verbose)

            dt = time.time() - t0
            msg = f"total running time: {dt:.2f} sec"
            log_stage(log, method, None, msg, verbose=verbose)
            return smps

        # ----------------------------------------
        # ---- CHECK IF PERTURBATIONS ARE PRODUCED
        # ----------------------------------------
        if smps == {}:
            msg = "no perturbation was produced"
            log_stage(log, method, None, msg, verbose=verbose)

            dt = time.time() - t0
            msg = f"total running time: {dt:.2f} sec"
            log_stage(log, method, None, msg, verbose=verbose)
            return

        # logging is already in the method
        endf6.apply_perturbations(
            smps,
            processes=iargs.processes,
            to_file=True,
            verbose=verbose,
            enable_tqdm=iargs.tqdm,
        )

        msg = "FY perturbations applied successfully"
        log_stage(log, method, None, msg, verbose=verbose)

        dt = time.time() - t0
        msg = f"total running time: {dt:.2f} sec"
        log_stage(log, method, None, msg, verbose=verbose)
        return       

    # -------------------------------------------------
    # ---- ERRORR KEYWORDS & REQUESTED MF
    # -------------------------------------------------
    nubar_req = (31 in iargs.mf)
    xs_req   = (33 in iargs.mf)
    chi_req  = (35 in iargs.mf)
    
    nubar = nubar_req and (31 in endf6.mf)
    xs    = xs_req   and ((33 in endf6.mf) or (32 in endf6.mf))  # MF32 handled with MF33
    mubar = False
    chi   = chi_req and (35 in endf6.mf)

    # Warn for requested-but-missing MFs
    if nubar_req and not nubar:
        msg = "MF=31 was requested but is not available in the file"
        log_stage(log, method, None, msg, verbose=verbose)

    if xs_req and not xs:
        msg = "MF=33 was requested but neither MF=33 nor MF=32 are available in the file"
        log_stage(log, method, None, msg, verbose=verbose)

    if chi_req and not chi:
        msg = "MF=35 was requested but is not available in the file"
        log_stage(log, method, None, msg, verbose=verbose)

    grids = dict(zip(["csewg239", "lanl30", "epri69", "ecco33"], [2, 3, 9, 19]))
    ign = grids[iargs.cov_energy_grid]

    errorr_kws = dict(
        verbose=verbose,
        err=err_errorr,
        xs=xs,
        nubar=nubar,
        chi=chi,
        mubar=mubar,
        groupr_kws=dict(nubar=nubar, chi=chi, mubar=mubar, ign=ign),  # both groupr and errorr take the same IGN
        errorr_kws=dict(ign=ign),
        suppress_njoy_output=not iargs.show_njoy,
        )
    if iargs.mt33:
        errorr_kws["errorr33_kws"] = dict(mt=iargs.mt33)

    smp_kws = _build_smp_kws(iargs)


    # ------------------------------------------------------------
    # ----  LOAD PERTURBATIONS or compute them
    # ------------------------------------------------------------
    if iargs.from_perturbations:
        base_dir, beg, end = iargs.from_perturbations
        
        if not os.path.isdir(base_dir):
            msg = f"--from_perturbations directory does not exist: '{base_dir}'"
            raise Exception(msg)

        if not (str(beg).isdigit() and str(end).isdigit()):
            msg = "--from_perturbations beg/end must be integers"
            raise Exception(msg)

        beg, end = int(beg), int(end)

        if not (beg >= 1 and end >= beg):
            msg = f"--from_perturbations range must satisfy 1 <= beg <= end (got {beg}..{end})",
            raise Exception(msg)
            

        ID = endf6.get_id()
        template = os.path.join(iargs.from_perturbations[0], "PERT_{}_MF{}.xlsx")

        smps = {}
        for mf in [31, 33, 35]:
            xls = template.format(ID, mf)
            if os.path.isfile(xls):
                msg = f"reading perturbations for MF={mf} from file '{xls}' (range {beg}..{end})"
                log_stage(log, method, None, msg, verbose=verbose)

                smps[mf] = Samples.from_excel(xls, beg=beg, end=end)
            
            else:
                msg = f"no perturbation file '{xls}' was found"
                log_stage(log, method, None, msg, verbose=verbose)

        if smps is {}:
            msg = f"No perturbation file was found for MAT={ID}"
            log_stage(log, method, None, msg, verbose=verbose)

    else:

        smps = endf6.get_perturbations(
            iargs.samples,
            njoy_kws=errorr_kws,
            smp_kws=smp_kws,
            verbose=verbose,
            suppress_njoy_output=not iargs.show_njoy,
            write_errorr=False,
            )
    
    # Only provide perturbations if requested
    if iargs.only_perturbations:
        msg = "returning perturbations only"
        if smps == {}:
            msg += " (empty dict)"
        
        log_stage(log, method, None, msg, verbose=verbose)

        dt = time.time() - t0
        msg = f"total running time: {dt:.2f} sec"
        log_stage(log, method, None, msg, verbose=verbose)
        return smps

    # ----------------------------------------
    # ---- CHECK IF PERTURBATIONS ARE PRODUCED
    # ----------------------------------------
    if smps == {}:
        msg = "no perturbation was produced"
        log_stage(log, method, None, msg, verbose=verbose)

        dt = time.time() - t0
        msg = f"total running time: {dt:.2f} sec"
        log_stage(log, method, None, msg, verbose=verbose)
        return

    # --------------------------------------
    # ---- COMPUTE TEMPERATURE
    # --------------------------------------
    temperature = _compute_temperature(iargs)

    # --------------------------------------
    # ---- PENDF / ACE keyword dictionaries
    # --------------------------------------
    pendf_kws = dict(
        verbose=verbose,
        err=err_pendf,
        suppress_njoy_output=not iargs.show_njoy,
        )

    ace_kws = dict(
        verbose=verbose,
        err=err_ace,
        temperature=temperature,
        purr=False,
        suppress_njoy_output=not iargs.show_njoy,
        )


    # --------------------------------------
    # ---- APPLY PERTURBATIONS
    # --------------------------------------
    endf6.apply_perturbations(
        smps,
        processes=iargs.processes,
        to_file=True,
        to_ace=iargs.acer,
        njoy_kws=pendf_kws,
        ace_kws=ace_kws,
        verbose=verbose,
        suppress_njoy_output=not iargs.show_njoy,
        enable_tqdm=iargs.tqdm,
    )

    msg = f"Sampling pipeline completed successfully for '{iargs.file}'"
    log_stage(log, method, None, msg, verbose=verbose)

    dt = time.time() - t0
    msg = f"total running time: {dt:.2f} sec"
    log_stage(log, method, None, msg, verbose=verbose)
    return


if __name__ == "__main__":
    run()
