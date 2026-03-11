# perturbation_base.py

"""
Shared utilities for perturbation application in sandy.

This module unifies:
- validation
- sample ID checking
- logging

Used internally by:
- Endf6.apply_perturbations_xs
- Endf6.apply_perturbations_rdd
"""



# -----------------------------------------------------------------------------
# Logging helper
# -----------------------------------------------------------------------------
import numpy as np

def log_stage(log, method, zam, message, verbose=False):
    prefix = f"{method}"

    if zam is None:
        # no ZAM printed at all
        pass

    elif np.isscalar(zam):
        prefix += f" | ZAM={zam}"

    else:
        prefix += " | MULTI-ZAM"

    log(f"{prefix} | {message}", verbose=verbose)


# -----------------------------------------------------------------------------
# Validation helpers
# -----------------------------------------------------------------------------
def validate_smps_mapping(smps):
    from collections.abc import Mapping
    if not isinstance(smps, Mapping):
        raise TypeError(
            f"`smps` must be a mapping (dict-like), got {type(smps).__name__}"
        )
    if len(smps) == 0:
        raise ValueError("`smps` is empty.")


def validate_required_keys(smps, required_keys, *, mode="all"):
    """
    Validate that required keys are present in `smps`, and return the keys that
    are present.

    Parameters
    ----------
    smps : Mapping
        The samples dictionary.
    required_keys : iterable
        Keys that are expected in `smps`.
    mode : {"all", "any"}, default="all"
        - "all": all required keys must be present.
        - "any": at least one required key must be present.

    Returns
    -------
    set
        The subset of required keys that are present in `smps`.

    Raises
    ------
    KeyError
        If validation fails based on the chosen mode.
    """

    smps_keys = set(smps.keys())
    required_keys = set(required_keys)

    present = smps_keys & required_keys

    if mode == "all":
        missing = required_keys - smps_keys
        if missing:
            raise KeyError(
                f"Missing required smps keys: {missing}. "
                f"Expected ALL of: {sorted(required_keys)}"
            )

    elif mode == "any":
        if not present:
            raise KeyError(
                f"None of the required smps keys found. "
                f"Expected ANY of: {sorted(required_keys)}"
            )

    else:
        raise ValueError("`mode` must be either 'all' or 'any'.")

    return present


def validate_sample_ids(smps, keys):
    """
    Ensure all provided smps[k] share:
    - the same set of sample IDs
    - in the same order

    Returns the ordered list of sample IDs.
    """
    # Ensure deterministic order even if keys is a set
    keys = list(keys)

    id_lists = {k: list(smps[k].data.columns) for k in keys}

    # Check all sets equal
    sets = {frozenset(v) for v in id_lists.values()}
    if len(sets) != 1:
        raise ValueError(
            "Sample ID sets mismatch across keys.\n" +
            "\n".join(f"{k}: {id_lists[k]}" for k in keys)
        )

    # Check ordering equal
    first = keys[0]
    order_ok = all(id_lists[k] == id_lists[first] for k in keys[1:])
    if not order_ok:
        raise ValueError(
            "Sample ID order mismatch across keys.\n" +
            "\n".join(f"{k}: {id_lists[k]}" for k in keys)
        )

    return id_lists[first]  # ordered sample IDs
