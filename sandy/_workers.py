# -*- coding: utf-8 -*-
"""
Created on Tue Mar  3 16:10:54 2026

@author: lfiorito
"""
import os
import time

def _endf6_perturb_worker(
        endf6,
        pendf,
        ismp,
        *,
        pxs=None,
        pnu=None,
        plpc=None,
        pchi=None,
        ace_kws: dict|None = None,
        to_ace: bool = False,
        verbose: bool = False,
        to_file: bool = False,
        **kwargs,
        ):

    """
    Worker to handle ENDF6 neutron data perturbation (xs, nubar, angular and energy distributions).

    Parameters
    ----------
    endf6 : `dict`
        `data` attribute of :obj:`~sandy.endf6.Endf6`.
        It contains the nominal ENDF6 data.
    pendf : `dict`
        `data` attribute of :obj:`~sandy.endf6.Endf6`.
        It contains the nominal PENDF data.
    ismp : `int`
        sample ID.
    pxs : `pd.DataFrame`
        It contains the perturbation coefficients for cross section.
        It corresponds to one single sample (in principle the one with ID `ismp`).
        It should have the same structure as a :obj:`~sandy.xs.Xs` object.
        The default is `None`.
    pnu: `pd.DataFrame`
        It contains the perturbation coefficients for nubar.
        It corresponds to one single sample (in principle the one with ID `ismp`).
        It should have the same structure as a :obj:`~sandy.xs.Xs` object.
        The default is `None`.
    plpc: `pd.DataFrame`
        Not implemented.
    pchi: `pd.DataFrame`
        It contains the perturbation coefficients for chi.
        It corresponds to one single sample (in principle the one with ID `ismp`).
        It should have the same structure as a :obj:`~sandy.xs.Xs` object.
        The default is `None`.
    verbose : `bool`, optional
        Flag to activate verbosity. The default is `False`.
    to_ace : TYPE, optional
        DESCRIPTION. The default is False.
    to_file : `bool`, optional
        Flag to write outputs to file. The default is `False`.
        This key changes the output type.
    ace_kws : `dict` , optional
        Additional keyword arguments for ACE file production. The default is {}.
    **kwargs : `dict`
        Additional keyword arguments (not used).

    Returns
    -------
    `dict`
        
        - if `to_file=False`: a `dict` with keys, values:
            
            - `endf6`: the ``dict`` mapping from ``.data`` of a perturbed :obj:`~sandy.endf6.Endf6` instance
            - `pendf`: the ``dict`` mapping from ``.data`` a perturbed :obj:`~sandy.endf6.Endf6` instance

        - if `to_file=True`: a `dict` with keys, values:

            - `endf6`: the filename of the perturbed ENDF6
            - `pendf`: the filenmae of the perturbed PENDF

    Examples
    --------
    Test that energy distributions are correctly perturbed.
    Example for Pu239.

    >>> import sandy, pandas as pd, numpy as np
    
    Creation of dummy perturbation: a perturbation of 10% up to 10 eV (included).

    >>> interval = pd.Interval(left=1e-8, right=10, closed="right")
    >>> idx = pd.MultiIndex.from_tuples([(9437, 18, interval)], names=("MAT", "MT", "E"))
    >>> pert = 1.1
    >>> df = pd.DataFrame([[pert]], index=idx).reset_index()
    >>> smps = sandy.Samples(df.set_index(["MAT", "MT", "E"]))
    >>> smps
    SMP                             0
    MAT  MT E                        
    9437 18 (1e-08, 10.0] 1.10000e+00

    Creation of reference ENDF6 and PENDF.

    >>> ref_endf6 = sandy.get_endf6_file("jeff_33", "xs", 942390, local=True)
    >>> ref_pendf = ref_endf6.get_pendf(err=1)

    Creation of perturbed data modifying the PFNS with the first perturbation
    sample (the only one).
    
    >>> ismp = 0
    >>> pchi = dict(smps.iterate_xs_samples())[ismp]
    >>> perturbed = sandy._workers._endf6_perturb_worker(ref_endf6.data, ref_pendf.data, ismp, pchi=pchi)
    >>> pert_endf6 = sandy.Endf6(perturbed['endf6'])
    >>> pert_pendf = sandy.Endf6(perturbed['pendf'])

    Compare reference and perturbed :class:`~sandy.edistr.Edistr`.

    >>> ref_edistr = sandy.Edistr.from_endf6(ref_endf6)
    >>> pert_edistr = sandy.Edistr.from_endf6(pert_endf6)

    Test that the perturbation is correct and happened below `ethresh` only.
    This test works for all incident energies.

    >>> lower_perturbed = pert_edistr.data.query("EOUT < 10").VALUE
    >>> lower_expected = ref_edistr.data.query("EOUT < 10").VALUE * pert
    >>> np.testing.assert_array_almost_equal(lower_perturbed, lower_expected)

    At energies >10 eV the perturbation did not apply.

    >>> higher_perturbed = pert_edistr.data.query("EOUT >= 10").VALUE
    >>> higher_expected = ref_edistr.data.query("EOUT >= 10").VALUE
    >>> np.testing.assert_array_almost_equal(higher_perturbed, higher_expected)



    Test that cross sections are correctly perturbed.
    Example for H1.

    Creation of a dummy perturbation for scattering `MT=2`: a perturbation of
    20% up to 10 eV (included).

    >>> interval = pd.Interval(left=1e-8, right=10, closed="right")
    >>> idx = pd.MultiIndex.from_tuples([(125, 2, interval)], names=("MAT", "MT", "E"))
    >>> pert = 1.2
    >>> df = pd.DataFrame([[pert]], index=idx).reset_index()
    >>> smps = sandy.Samples(df.set_index(["MAT", "MT", "E"]))
    >>> smps
    SMP                            0
    MAT MT E                        
    125 2  (1e-08, 10.0] 1.20000e+00

    Creation of reference ENDF6 and PENDF.

    >>> ref_endf6 = sandy.get_endf6_file("jeff_33", "xs", 10010, local=True)
    >>> ref_pendf = ref_endf6.get_pendf(err=1)

    Creation of perturbed data modifying the PFNS with the first perturbation
    sample (the only one).
    
    >>> ismp = 0
    >>> pxs = dict(smps.iterate_xs_samples())[ismp]
    >>> perturbed = sandy._workers._endf6_perturb_worker(ref_endf6.data, ref_pendf.data, ismp, pxs=pxs)
    >>> pert_endf6 = sandy.Endf6(perturbed['endf6'])
    >>> pert_pendf = sandy.Endf6(perturbed['pendf'])

    Creation of reference and perturbed :class:`~sandy.xs.Xs`.
    
    >>> ref_xs = sandy.Xs.from_endf6(ref_pendf)
    >>> pert_xs = sandy.Xs.from_endf6(pert_pendf)

    Test that the perturbation is correct.

    >>> lower_perturbed = pert_xs.data.query("E<=10")[(125, 2)]
    >>> lower_expected = ref_xs.data.query("E<=10")[(125, 2)] * 1.2
    >>> np.testing.assert_array_almost_equal(lower_perturbed, lower_expected)

    At energies >10 eV the perturbation did not apply.

    >>> higher_perturbed = pert_xs.data.query("E>10")[(125, 2)]
    >>> higher_expected = ref_xs.data.query("E>10")[(125, 2)]
    >>> np.testing.assert_array_almost_equal(higher_perturbed, higher_expected)

    And the total cross section also changed because reconstructed from the perturbed scattering.

    >>> assert not np.array_equal(pert_xs.data[(125, 1)], ref_xs.data[(125, 1)])

    """
    # ---- IMPORTS (local to keep worker import cost low at module import time)
    from typing import Any
    import pandas as pd

    from .xs import Xs
    from .edistr import Edistr
    from .utils import log
    from .endf6 import Endf6
    from ._perturbation_base import log_stage

    # ---- BASIC INPUT CHECK (lightweight)
    if not isinstance(endf6, dict):
        raise TypeError("'endf6' must be plain dict (sandy.Endf6.data)")
    if not isinstance(pendf, dict):
        raise TypeError("'pendf' must be plain dict (sandy.Endf6.data)")
    for name, inst in zip(
            ["pnu", "pxs", "plpc", "pchi"],
            [pnu, pxs, plpc, pchi],
            ):
        if inst is not None:
            if not isinstance(inst, pd.DataFrame):
                raise TypeError(f"'{name}' must be pd.DataFrame (index=E, columns=[MAT, MT, ...])")

    # ---- NORMALIZE dict-like keyword arguments
    # will be used if to_ace=True
    ace_kws_ = {} if ace_kws is None else ace_kws

    # ---- INITIALIZE working copies
    # Get them back as Endf6 instances (never mutate cached inputs)
    # A full deepcopy is not necessary (the dict values are string and are immutable), and it’s expensive.
    endf6_pert = Endf6(endf6.copy())
    pendf_pert = Endf6(pendf.copy())

    # ---- LOGGING SETUP
    zam = endf6_pert.get_zam()
    method = (
        f"[PID {os.getpid():5d}] XS-worker "
        f"| sample={ismp:4d}"
        )

    # ---- LOGGING: start perturbation
    # we need to get the ZAM before logging to include it in the message
    msg = f"(NU={pnu is not None}, XS={pxs is not None}, PFNS={pchi is not None}) | start"
    log_stage(log, method, zam, msg, verbose=verbose)
    t0 = time.perf_counter()

    # ---- NUBAR PERTURBATION (on ENDF6)
    if pnu is not None:
        nu = Xs.from_endf6(endf6_pert.filter_by(listmt=[452, 455, 456]))
        nu_pert = nu._perturb(pnu)   # internal: expects "single-sample" payload
        endf6_pert = nu_pert.reconstruct_sums(drop=True).to_endf6(endf6_pert).update_intro()


    # ---- ANGULAR DISTR. PERTURBATION (placeholder)
    if plpc is not None:
        # reserved for future implementation
        pass

    # ---- PFNS PERTURBATION (on ENDF6)
    if pchi is not None:
        # Applies the same perturbation to all incident particle energies (EIN) and K blocks
        edistr_pert_block = []
        
        # Group data by EIN and K for processing
        for (ein, k), df in Edistr.from_endf6(endf6_pert).data.groupby(['EIN', 'K']):

            # Turn edistr into an "xs-like" layout to reuse xs._perturb
            dummy_xs = Xs(
                df.rename({"EOUT": "E"}, axis=1)
                  .set_index(["MAT", "MT"])[["E","VALUE"]]
                  .pivot(columns="E").T.droplevel(level=0)
            )

            # Apply perturbation to dummy energy distribution
            dummy_xs_pert = dummy_xs._perturb(pchi)
            
            # Transform xs data back into edistr form
            perturbed_data = (
                dummy_xs_pert.data.stack([1, 0], future_stack=True)  # Use future_stack=True to adopt the new behavior
                .to_frame()
                .reset_index()
                .rename({"E": "EOUT", 0: "VALUE"}, axis=1)
                .assign(K=k, EIN=ein)
                [["MAT", "MT", "K", "EIN", "EOUT", "VALUE"]]
            )
            edistr_pert_block.append(perturbed_data)

        # Combine and normalize perturbed data, then update ENDF6
        endf6_pert = (
            Edistr(pd.concat(edistr_pert_block, ignore_index=True))
            .normalize()
            .to_endf6(endf6_pert)
            .update_intro()
            )

    # ---- XS PERTURBATION (on PENDF)
    if pxs is not None:
        xs = Xs.from_endf6(pendf_pert)
        xs_pert = xs._perturb(pxs)
        pendf_pert = xs_pert.reconstruct_sums(drop=True).to_endf6(pendf_pert).update_intro()

    # ---- ASSEMBLE OUTPUTS
    # Return perturbed ENDF6 and PENDF instances as dict
    # The items in the output dict are {str: dict}
    out_dict: dict[str, Any] = {
        "endf6": endf6_pert.data,  # dict
        "pendf": pendf_pert.data,  # dict
    }

    # ---- OPTIONAL ACE/XSDIR (from ENDF6 + PENDF pair)
    if to_ace:
        
        # ---- LOGGING: processing
        msg = "processing ACE"
        log_stage(log, method, zam, msg, verbose=verbose)

        # Add ACE file content and XSDIR file content to the dict
        out_extra = endf6_pert.get_ace(pendf=pendf_pert, **ace_kws_)   # this is a dict {'ace': text, 'xsdir': text}

        # I make the dict update explicit
        # The added items to the output dict are {str: str}
        out_dict.update({
            "ace": out_extra["ace"],      # str
            "xsdir": out_extra["xsdir"],  # str
            })

    # ---- OPTIONAL WRITE TO FILE
    # return only filename
    if to_file:
        # ---- LOGGING: writing
        msg = "writing to file"
        log_stage(log, method, zam, msg, verbose=verbose)

        # --- The output files basename is hardcoded. It is too complex to maintain a flexible structure
        basename = endf6_pert._derive_basename_for_sampling(ismp)
        
        # --- Write files to disk and return only the filenames
        # The added items to the output dict are {str: str}
        out_dict = _write_files_worker(out_dict, ismp, basename=basename, zam=zam, verbose=verbose)

    # ---- LOGGING: end
    dt = time.perf_counter() - t0
    msg = f"finished in {dt:.3f} s"
    log_stage(log, method, zam, msg, verbose=verbose)

    return out_dict


def _write_files_worker(
        file_dict,
        ismp: int,
        basename: str = "output",
        zam: int|None = None,
        verbose: bool = False,
        ):
    """
    Write ENDF-6, PENDF, and (optionally) ACE and XSDIR contents to disk.

    Parameters
    ----------
    file_dict : Mapping[str, Any]
        Dictionary with optional keys:
          - "endf6": dict-like data to build a `:obj:~sandy.endf6.Endf6` object
          - "pendf": dict-like data to build a `:obj:~sandy.endf6.Endf6` object
          - "ace"  : str content of ACE file (text)
          - "xsdir": str content of XSDIR file (text)

        It is assumed the ENDF-6/PENDF dicts are compatible with `:obj:~sandy.endf6.Endf6`.

    basename : str, default "output"
        Basename (without extension) used for ENDF-6/PENDF files,
        and also as prefix for ACE (*.XXc).

    verbose : bool, default False
        If True, emits progress messages via `logger`.

    Returns
    -------
    Dict[str, str]
        Mapping with the filenames written. Keys may include:
        - "endf6": path to `<basename>.endf6`
        - "pendf": path to `<basename>.pendf`
        - "ace"  : path to `<basename>.<SUFFIX>c` (if "ace" provided)
        - "xsdir": path to `<basename>.<SUFFIX>c` (if "xsdir" provided)

    Raises
    ------
    KeyError
        If required keys "endf6" or "pendf" are missing in `file_dict`.

    Notes
    -----
    - We inferring the ACE/XSDIR suffix from the ACE/XSDIR content.
    - ENDF6 and PENDF are dict and not `:obj:~sandy.endf6.Endf6` to be pickled.

    Examples
    --------
    
    Collect endf6/pendf/ace/xsdir file for test.
    
    >>> import sandy, os, glob
    >>> endf6 = sandy.get_endf6_file("jeff_33", "xs", 10010, local=True)
    >>> pendf = endf6.get_pendf()
    >>> out_dict = {"endf6": endf6.data, "pendf": pendf.data}
    >>> out_dict |= endf6.get_ace(temperature=0)
    
    ``out_dict`` is a dictionary with keys ``endf6``, ``pendf``,
    ``ace`` and ``xsdir``, like the samples produced by
    ``Endf6.get_perturbations_xs``.
    
    Run worker...but first remove outputs.

    >>> # Delete files
    >>> for f in glob.glob("output*"):
    ...     try: os.remove(f)
    ...     except FileNotFoundError: pass
    >>> # Run worker
    >>> ismp = 1
    >>> outfiles = sandy._workers._write_files_worker(out_dict, ismp)

    Check that outputs have been created.

    >>> output_files = glob.glob("output*")
    >>> expected_outputs = {'output.00c', 'output.00c.xsd', 'output.endf6', 'output.pendf'}
    >>> assert expected_outputs.issubset(set(output_files))
    """
    # ---- IMPORT
    from .utils import log
    from .endf6 import Endf6
    from ._perturbation_base import log_stage

    # ---- LOGGING SETUP
    method = (
        f"[PID {os.getpid():5d}] WRITE-worker "
        f"| sample={ismp:4d}"
        )

    # --- Write files to disk and return only the filenames
    outfiles= {}
    
    def extract_suffix(text):
        """Extract the suffix (example, "03" from 1001.03c) from the ace/xsdir file.
        No fallback if it doesn't work"""
        suffix = text.split()[0].split(".")[1][:2]
        return suffix
    
    if 'ace' in file_dict:
        suffix = extract_suffix(file_dict["ace"])
        file = f"{basename}.{suffix}c"
        msg = f"writing ACE to file -> '{file}'"
        log_stage(log, method, zam, msg, verbose=verbose)

        with open(file, "w") as f:
            f.write(file_dict["ace"])

        outfiles["ace"] = file

    if 'xsdir' in file_dict:
        suffix = extract_suffix(file_dict["xsdir"])
        file = f"{basename}.{suffix}c.xsd"
        msg = f"writing XSD to file -> '{file}'"
        log_stage(log, method, zam, msg, verbose=verbose)

        with open(file, "w") as f:
            f.write(file_dict["xsdir"])

        outfiles["xsdir"] = file
    
    if 'endf6' in file_dict:
        endf6 = Endf6(file_dict["endf6"])
        file = f"{basename}.endf6"
        msg = f"writing ENDF-6 to file -> '{file}'"
        log_stage(log, method, zam, msg, verbose=verbose)
        endf6.to_file(file)
        outfiles["endf6"] = file

    if 'pendf' in file_dict:
        pendf = Endf6(file_dict["pendf"])
        file = f"{basename}.pendf"
        msg = f"Writing PENDF to file -> '{file}'"
        log_stage(log, method, zam, msg, verbose=verbose)
        pendf.to_file(file)
        outfiles["pendf"] = file
    
    return outfiles
    

def _rdd_perturb_worker(
        endf6,
        rdd,
        phl,
        pde,
        pbr,
        ismp,
        *,
        verbose: bool = False,
        to_file: bool = False,
        ):
    """
    Worker to handle ENDF6 radioactive decay data perturbation.

    Parameters
    ----------
    endf6 : `dict`
        `data` attribute of :obj:`~sandy.endf6.Endf6`.
        It contains the nominal ENDF6 data for a singl or multi-ZAM decay file.
    rdd : `pd.DataFrame`
        `data` attribute of :obj:`~sandy.decay.DecayData`.
        It contains the nominal decay data.
    phl : `pd.DataFrame`
        `data` attribute of :obj:`~sandy.samples.Samples`.
        It contains the perturbation coefficients for half-lives.
        It corresponds to one single sample (in principle the one with ID `ismp`).
        Compared to :func:`sandy._workers._endf6_perturb_worker` it is a mandatory argument.
    pde : `pd.DataFrame`
        `data` attribute of :obj:`~sandy.samples.Samples`.
        It contains the perturbation coefficients for decay energies.
        It corresponds to one single sample (in principle the one with ID `ismp`).
        Compared to :func:`sandy._workers._endf6_perturb_worker` it is a mandatory argument.
    pbr : `pd.DataFrame`
        `data` attribute of :obj:`~sandy.samples.Samples`.
        It contains the perturbation coefficients for branching ratios.
        It corresponds to one single sample (in principle the one with ID `ismp`).
        Compared to :func:`sandy._workers._endf6_perturb_worker` it is a mandatory argument.
    ismp : `int`
        sample ID.
    verbose : `bool`, optional
        Flag to activate verbosity. The default is False.
    to_file : `bool`, optional
        Flag to write outputs to file. The default is False.
        This key changes the output type.

    Returns
    -------
    `dict`
        
        - if `to_file=False`: a `dict` with keys, values:
            
            - `endf6`: the ``dict`` mapping from ``.data`` of a perturbed :obj:`~sandy.endf6.Endf6` instance

        - if `to_file=True`: a `dict` with keys, values:

            - `endf6`: the filename of the perturbed ENDF6

    Notes
    -----
    - This method is written so that it can be pickled.
    - Branching ratios are renormalized.
    """
    import time
    import pandas as pd

    # ---- IMPORTS (local to keep worker import cost low at module import time)
    from .decay import DecayData
    from .utils import log
    from .endf6 import Endf6
    from ._perturbation_base import log_stage

    # ---- BASIC INPUT CHECK (lightweight)
    if not isinstance(endf6, dict):
        raise TypeError("'endf6' must be plain dict (sandy.Endf6.data)")
    if not isinstance(rdd, dict):
        raise TypeError("'rdd' must be plain dict(sandy.DecayData.data)")
    if not isinstance(phl, pd.DataFrame):
        raise TypeError("'phl' must be a pd.DataFrame (index=ZAM, columns=HL)")
    if not isinstance(pde, pd.DataFrame):
        raise TypeError("'pde' must be a pd.DataFrame (index=[ZAM, TYPE], columns=DE)")
    if not isinstance(pbr, pd.DataFrame):
        raise TypeError("'pbr' must be a pd.DataFrame (index=[ZAM, RTYP, RFS], columns=BR)")

    # ---- SHALLOW COPIES (avoid expensive deep copies when safe)
    # Endf6/DecayData will copy internals appropriately when transformed;
    # we avoid deep copy of large dicts here and leave copying to domain objects.
    endf6_ = Endf6(endf6)           # Endf6 typically copies internally as needed
    rdd_ = DecayData(rdd)           # same for DecayData

    # ---- LOGGING SETUP
    method = (
        f"[PID {os.getpid():5d}] RDD-worker "
        f"| sample={ismp:4d}"
        )
    zam = endf6_.get_zam()

    msg = "start"
    log_stage(log, method, zam, msg, verbose=verbose)
    t0 = time.perf_counter()

    # ---- APPLY HL
    # Multiply half-life
    hl_ = rdd_.get_half_life()
    # Align by index; both are Series with same index (ZAM)
    # Using .mul to ensure proper alignment and preserve index.
    # Guard against missing values (fill_value=1.0 means “no change” where missing)
    key = "HL"
    hl_series = hl_.data[key]
    hl_factor = phl[key]  # becomes a series
    # Some implementations keep HL as float; enforce numeric and safe multiply
    hl_.data[key] = hl_series.mul(hl_factor, fill_value=1.0)
    rdd_ = hl_.to_decaydata(rdd_)

    # ---- APPLY DE
    # Multiply decay energies (column "E" in de_.data; indexed like (ZAM, kind) or similar)
    de_ = rdd_.get_decay_energy()
    key = "E"
    de_series = de_.data[key]
    de_factor = pde[key]  # becomes a series
    de_.data[key] = de_series.mul(de_factor, fill_value=1.0)
    rdd_ = de_.to_decaydata(rdd_)

    # ---- APPLY BR + RENORMALIZE
    # Multiply BR, then renormalize per parent (sum of branches = 1 where applicable).
    br_ = rdd_.get_branching_ratio()
    key = "BR"
    br_series = br_.data[key]
    br_factor = pbr[key]  # becomes a series
    br_.data[key] = br_series.mul(br_factor, fill_value=1.0)

    # Renormalization: rely on br_.normalize()
    rdd_ = br_.normalize().to_decaydata(rdd_)
    
    # ---- WRITE BACK INTO ENDF6
    out = rdd_.to_endf6(endf6_)   # Endf6 object, cannot be pickled

    # ---- RETURN TYPE
    if to_file:
        # ---- LOGGING: writing
        # continue and return filename where data was written
        filename = f"decay_data_{ismp}"

        msg = f"writing to file '{filename}'"
        log_stage(log, method, zam, msg, verbose=verbose)

        out.to_file(filename)

        out = filename   # str

    else:
        out = out.data   # dict

    # ---- LOGGING: end
    dt = time.perf_counter() - t0
    msg = f"finished in {dt:.3f} s"
    log_stage(log, method, zam, msg, verbose=verbose)

    return out



def _fy_perturb_worker(
        endf6,
        fy,
        pfy,
        ismp,
        *,
        verbose=False,
        to_file=False,
        ):
    """
    Worker to handle ENDF6 fission yield perturbation.

    Parameters
    ----------
    endf6 : `dict`
        `data` attribute of :obj:`~sandy.endf6.Endf6`.
        It contains the nominal ENDF6 data.
    fy : `pd.DataFrame`
        `data` attribute of :obj:`~sandy.fy.Fy`.
        It contains the nominal fission yield data.
        It i sassume they match all the ZAP of the samples.
    smps : `pd.DataFrame`
        It contains the perturbation coefficients for fission yields.
        Columns are `MAT`, `MT`, `E`, `ZAM`, `ZAP`, `SMP`, `VALS`.
        This dataframe is generally produced with `pd.pivot_table`.
    ismp : `int`
        sample ID.
    verbose : `bool`, optional
        Flag to activate verbosity. The default is False.
    to_file : `bool`, optional
        Flag to write outputs to file. The default is False.
        This key changes the output type.
    **kwargs : `dict`
        Additional keyword arguments (not used).
        
    Notes
    -----
    .. note:: It follows the logic of :obj:`~sandy.endf6._endf6_perturb_worker` and
              :obj:`~sandy.endf6._rdd_perturb_worker`.

    Returns
    -------
    `dict`
        Either a dictionary of :obj:`~sandy.endf6.Endf6` instances for each set of
        perturbation coefficients (if `to_file=False`), or a dictionary
        of `str` with the output file name for each set of perturbation
        coefficients.

    Notes
    -----
    .. note: This method is written so that it can be handled by the
             `multiprocess` module (pickling).

    Examples
    --------
    
    Default test: create 1 sample and perturb fission yields for 1 fissioning system.
    
    >>> import sandy, pandas as pd, numpy as np

    >>> nsmp = 1   # sample size
    >>> zam, e = 922350, 0.0253
    >>> tape = sandy.get_endf6_file("jeff_33", "nfpy", zam, local=True)
    >>> nfpy = sandy.Fy.from_endf6(tape)
    >>> idx = nfpy.data.query(f"E=={e} & MT==454 & ZAM=={zam}").index
    >>> fy = nfpy.data.loc[idx]
    >>> smps = sandy.CategoryCov(pd.DataFrame(np.diag((fy.DFY/fy.FY)**2), index=fy.ZAP, columns=fy.ZAP).fillna(0)).sampling(nsmp)
    >>> smps = smps.data.rename_axis(index="ZAP").stack().rename("VALS").reset_index().assign(E=e, ZAM=zam)[["ZAM", "E", "ZAP", "SMP", "VALS"]]
    >>> out = sandy._workers._fy_perturb_worker(tape.data, nfpy.data, smps, nsmp-1, verbose=True, to_file=False)
    >>> out = sandy.Endf6(out)
    
    Silly test: assert the `MT=454` was changed, and `MT=459` was not.

    >>> assert sandy.Fy.from_endf6(out).data.query("MT==459").equals(nfpy.data.query("MT==459"))
    >>> assert not sandy.Fy.from_endf6(out).data.query("MT==454").equals(nfpy.data.query("MT==454"))
    
    Test to check that the output random ENDF6 are perturbed correctly.

    >>> tape = sandy.get_endf6_file("jeff_33", "nfpy", 922350, local=True)
    >>> smps = tape.get_perturbations(2, covariance=None)
    >>> nfpy = sandy.Fy.from_endf6(tape)
    >>> out = sandy._workers._fy_perturb_worker(tape.data, nfpy.data, smps, 0)
    >>> nfpy0 = sandy.Fy.from_endf6(sandy.Endf6(out))

    Assert that ratio of perturbed to nominal FY's is equal to samples.

    >>> n = nfpy.data.query("ZAM==922350 and MT==454")
    >>> n0 = nfpy0.data.query("ZAM==922350 and MT==454")
    >>> assert not n.equals(n0)
    >>> sp = (n0.set_index(["MAT", "MT", "ZAM", "E", "ZAP"]).FY /  n.set_index(["MAT", "MT", "ZAM", "E", "ZAP"]).FY).fillna(1)
    >>> p = smps.query("ZAM==922350 and SMP==0").VALS
    >>> np.testing.assert_array_almost_equal(p, sp, decimal=4)
    """
    import time
    import pandas as pd
    import numpy as np

    # ---- IMPORTS (local to keep worker import cost low at module import time)
    from .fy import Fy
    from .utils import log
    from .endf6 import Endf6
    from ._perturbation_base import log_stage

    # ---- BASIC INPUT CHECK (lightweight)
    if not isinstance(endf6, dict):
        raise TypeError("'endf6' must be plain dict (sandy.Endf6.data)")
    if not isinstance(fy, pd.DataFrame):
        raise TypeError("'fy' must be a pd.DataFrame (sandy.Fy.data)")
    if not isinstance(pfy, pd.DataFrame):
        raise TypeError("'pfy' must be a pd.DataFrame (index=[ZAM, E, ZAP], columns=IFY).")

    # ---- SHALLOW COPIES (avoid expensive deep copies when safe)
    endf6_ = Endf6(endf6.copy())  # this was a dictionary
    fy_ = Fy(fy.copy())    # this was a dataframe

    # ---- LOGGING SETUP
    method = (
        f"[PID {os.getpid():5d}] FY-worker "
        f"| sample={ismp:4d}"
        )
    zam = endf6_.get_zam()

    msg = "start"
    log_stage(log, method, zam, msg, verbose=verbose)
    t0 = time.perf_counter()


    # ---- GROUP BY FISSIONING SYSTEM (group by MultiIndex levels)
    grouping_levels = ["ZAM", "E"]
    for (zam, E), block in pfy["IFY"].groupby(grouping_levels):

        # Reduce to a Series indexed only by ZAP
        block_zap = block.droplevel(grouping_levels)

        # Extract matching nominal IFYs (MF8/MT454) in the tape
        idx = fy_.data.query("MT == 454 and ZAM == @zam and E == @E").index

        if len(idx) == 0:
            raise Exception(
                f"No nominal fission yields for ZAM={zam}, E={E}"
            )

        # Ensure ZAP alignment
        # Nominal ZAP order for this (ZAM,E) block
        zap_order = fy_.data.loc[idx, "ZAP"].to_numpy()

        # Reindex perturbation vector to match nominal ZAP ordering
        perts = block_zap.reindex(zap_order).to_numpy()

        # Validate presence of all ZAP perturbations
        found_nan = np.isnan(perts)
        if found_nan.any():
            missing = zap_order[found_nan]
            raise ValueError(
                f"Missing FY perturbations for sample {ismp} at ZAM={zam}, E={E}, "
                f"ZAP={missing.tolist()}"
            )
        
        # IMPORTANT, this does not update the CFYs, which in random ENDF-6 file are inconsistent with the perturbed IFYs
        fy_.data.loc[idx, "FY"] *= perts

    # ---- WRITE BACK INTO ENDF6
    out = fy_.to_endf6(endf6_)   # Endf6 object, cannot be pickled

    # ---- RETURN TYPE
    if to_file:
        # ---- LOGGING: writing
        # continue and return filename where data was written
        filename = f"fy_{ismp}"

        msg = f"writing to file '{filename}'"
        log_stage(log, method, zam, msg, verbose=verbose)

        out.to_file(filename)

        out = filename   # str

    else:
        out = out.data   # dict

    # ---- LOGGING: end
    dt = time.perf_counter() - t0
    msg = f"finished in {dt:.3f} s"
    log_stage(log, method, zam, msg, verbose=verbose)

    return out
