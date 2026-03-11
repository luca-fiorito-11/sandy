"""
Collection of utilities, functions and classes that are requetsed in all code
components.

This module should not depend on any other sandy modules.
"""
import numpy as np
import logging
from contextlib import contextmanager
import os
import re
from functools import wraps
from contextlib import nullcontext
import inspect


__author__ = "Luca Fiorito"



def add_delimiter_every_n_characters(string, step, delimiter=" "):
    return delimiter.join(string[i:i+step] for i in range(0, len(string), step))


def add_exp_in_endf6_text(text):
    pattern = re.compile(r"([0-9\.])([+-])([0-9])")
    return pattern.sub(r"\g<1>E\g<2>\g<3>", text)


def get_seed():
    """
    Wrapper to `np.random.SeedSequence().entropy`.
    """
    seed = np.random.SeedSequence().entropy
    return seed


def grouper(iterable, n, fillvalue=None):
    """
    Collect data into fixed-length chunks or blocks
    """
    from itertools import zip_longest

    args = [iter(iterable)] * n
    return zip_longest(*args, fillvalue=fillvalue)


def interwine_lists(a, b):
    return [ z for item in zip(a, b) for z in item ]


def is_valid_file(parser, arg, r=True, w=False, x=False):
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
    if os.path.isdir(arg):
        return arg
    if mkdir:
        os.makedirs(arg, exist_ok=True)
    else:
        parser.error("Directory {} does not exist".format(arg))
    return arg


def pad_from_beginning(vals, maxlen=None, value=0., axis=0):
    """
    Convert list of arrays into matrix by backward padding.
    .. note:: this function can be used to put cross sections into one matrix
              by adding zeros before the first value of threshold reactions.

    Parameters
    ----------
    vals : `iterable` of arrays/lists
        values of the matrix
    maxlen : `int`, optional, default is `None`
        length to fill with padding.
        If not given, use the maximum length of the arrays in `vals`
        .. important:: if given, maxlen should be
                       `maxlen <= max([len(v) for v in vals])`

    value : `float`, optional, default is `0.`
        value used for padding
    axis : `int`, optional, either `0` or `1`, default is `0`
        axis along whihc the arrays should be positioned.
        `0` means rows, `1` means columns

    Returns
    -------
    `numpy.array`
        2D `numpy.array` with shape `(len(vals), maxlen)` if `axis=0` and
        `(maxlen, len(vals))` if `axis=1`

    Raises
    ------
    `ValueError`
        if `axis` is neither `0` nor `1`
    `ValueError`
        if `maxlen <= max([len(v) for v in vals])`
    """
    length = len(vals)
    lens = [len(v) for v in vals]                     # only iteration
    maxlen_ = max(lens)
    if maxlen is None:
        pass
    elif maxlen < maxlen_:
        raise ValueError("'maxlen' must be >= '{}'".format(maxlen_))
    else:
        maxlen_ = maxlen
    matrix = np.ones((length, maxlen_), dtype=float)*value
    mask = np.arange(maxlen_)[::-1] < np.array(lens)[:, None]  # key line
    matrix[mask] = np.concatenate(vals)
    if axis == 0:
        return matrix
    elif axis == 1:
        return matrix.T
    else:
        raise ValueError("'axis' can be '0' (rows) or '1' (columns), not '{}'".format(axis))


def pad_from_beginning_fast(vals, maxlen):
    """
    Like `aleph.utils.pad_from_beginning` but faster.
    Keyword arguments `axis` and `values` take the default options.
    .. note:: this function can be used to put cross sections into one matrix
              by adding zeros before the first value of threshold reactions.

    Parameters
    ----------
    vals : `iterable` of arrays/lists
        values of the matrix
    maxlen : `int`
        length to fill with padding.

    Returns
    -------
    `numpy.array`
        2D `numpy.array` with shape `(len(vals), maxlen)`
    """
    length = len(vals)
    matrix = np.zeros((length, maxlen))
    lens = [len(v) for v in vals]                    # only iteration
    mask = np.arange(maxlen)[::-1] < np.array(lens)[:, None]  # key line
    matrix[mask] = np.concatenate(vals)
    return matrix


def reshape_bfill(x, y, xnew, left_values="first", right_values=0):
    """
    Interpolate array over new energy grid structure using "bfill" method.

    Right-extrapolated values are replaced by zeros.
    Left-extrapolated values are replaced by `y[0]`.

    Parameters
    ----------
    x : 1d array-like object with at least two entries
        energy grid
    xnew : 1d array-like object with at least two entries
        new energy grid
    y : `numpy.ndarray` with at least two entries and same length as `x`
        array to interpolate

    Returns
    -------
    `numpy.ndarray` with length `len(xnew)`
        interpolated array
    """
    from scipy.interpolate import interp1d

    fill_value = [left_values, right_values]
    if left_values == "first":
        fill_value[0] = y[0]
    fill_value = tuple(fill_value)
    foo = interp1d(
            x, y,
            axis=0,
            copy=False,
            kind="next",
            bounds_error=False,
            fill_value=fill_value,
            assume_sorted=True,
            )
    return foo(xnew)


def reshape_differential(x, y, xnew):
    """
    Linearly interpolate array over new energy grid structure.

    Extrapolated values are replaced by zeros.

    Parameters
    ----------
    x : 1d array-like object with at least two entries
        energy grid
    xnew : 1d array-like object with at least two entries
        new energy grid
    y : `numpy.ndarray` with at least two entries and same length as `x`
        array to interpolate

    Returns
    -------
    `numpy.ndarray` with length `len(xnew)`
        interpolated array
    """
    from scipy.interpolate import interp1d

    foo = interp1d(
            x, y,
            axis=0,
            copy=False,
            kind="slinear",
            bounds_error=False,
            fill_value=0.,
            assume_sorted=True,
            )
    return foo(xnew)


def reshape_integral(x, y, xnew, left_values="first", right_values=0):
    """
    Interpolate array over new energy grid structure using "bfill" method.
    It is assumed that the values of `y` are  multiplied by the grid bin-width.
    The values of the interpolated array are recalculated proportionally to
    the new grid bin-widths.

    Extrapolated values are replaced by zeros.

    Parameters
    ----------
    x : 1d array-like object with at least two entries
        energy grid
    xnew : 1d array-like object with at least two entries
        new energy grid
    y : `numpy.ndarray` with at least two entries and same length as `x`
        array to interpolate

    Returns
    -------
    `numpy.ndarray` with length `len(xnew)`
        interpolated array
    """
    dx = x.copy()
    dx[1:] = np.ediff1d(x)
    dxnew = xnew.copy()
    dxnew[1:] = np.ediff1d(xnew)
    out = reshape_bfill(
        x,
        y / dx,
        xnew,
        left_values=np.nan,
        right_values=np.nan
        ) * dxnew
    return out


def uniform_loggrid(xmin, xmax, npoints=100):
    """
    Given lower and upper limits, produce a grid with a number of points
    `npoints` that define equivalent intervals in log scale.

    Parameters
    ----------
    xmin : `float`
        lower bound of the grid structure
    xmax : `float`
        upper bound of the grid structure
    npoints : `int`, optional, default `100`

    Returns
    -------
    `numpy` array
        grid equally spaced in logarithmic scale

    Examples
    --------
    >>> uniform_loggrid(1e-5, 1e7, 25)
    array([1.00000000e-05, 3.16227766e-05, 1.00000000e-04, 3.16227766e-04,
           1.00000000e-03, 3.16227766e-03, 1.00000000e-02, 3.16227766e-02,
           1.00000000e-01, 3.16227766e-01, 1.00000000e+00, 3.16227766e+00,
           1.00000000e+01, 3.16227766e+01, 1.00000000e+02, 3.16227766e+02,
           1.00000000e+03, 3.16227766e+03, 1.00000000e+04, 3.16227766e+04,
           1.00000000e+05, 3.16227766e+05, 1.00000000e+06, 3.16227766e+06,
           1.00000000e+07])
    """
    return 10.0**np.linspace(np.log10(xmin), np.log10(xmax), npoints)



def star(func):
    def inner(*args, **kwargs):
        print("*" * 30)
        func(*args, **kwargs)
        print("*" * 30)
    return inner


def percent(func):
    def inner(*args, **kwargs):
        print("%" * 30)
        func(*args, **kwargs)
        print("%" * 30)
    return inner


def which(program):
    """
    Mimic the behavior of the UNIX 'which' command.     
    """
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)
    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file
    return None


def force_symlink(file1, file2):
    """
    Mimic the behavior of the UNIX 'ln -sf' command.    
    """
    try:
        os.symlink(file1, file2)
    except FileExistsError:
        os.remove(file2)
        os.symlink(file1, file2)


def log(msg, *, level=logging.INFO, logger=None, verbose=None):
    """
    Emit `msg` using the provided `logger` (preferred). If no logger provided,
    fall back to a module-scoped logger. Do not print directly.
    """
    lg = logger or logging.getLogger(__name__)
    # If you want to respect a 'verbose' flag, ensure it *only* affects INFO/DEBUG,
    # not WARNING+, so warnings are still emitted unless suppressed by logger filters.
    if verbose is False and level < logging.WARNING:
        return
    lg.log(level, msg)


# @contextmanager
# def suppress_logging_warnings(logger_name="sandy.warn"):
#     logger = logging.getLogger(logger_name)
#     old_level = logger.level
#     try:
#         logger.setLevel(logging.ERROR)   # suppress WARNING and INFO
#         yield
#     finally:
#         logger.setLevel(old_level)


# @contextmanager
# def suppress_logging_warnings(logger_name_prefix: str, *, level=logging.WARNING):
#     """
#     Suppress log records at `level` or higher for any logger whose name
#     starts with `logger_name_prefix` (exact match or a dotted descendant),
#     regardless of handler/propagation configuration.

#     This works by attaching a filter to the *root logger* so that even
#     propagated records are dropped.
#     """
#     root = logging.getLogger()  # root logger

#     class _PrefixMaxLevelFilter(logging.Filter):
#         def filter(self, record: logging.LogRecord) -> bool:
#             name = record.name  # e.g., 'sandy.warn' or 'sandy.warn.sub'
#             targeted = (name == logger_name_prefix or
#                         name.startswith(logger_name_prefix + "."))
#             if targeted and record.levelno >= level:
#                 return False  # drop it
#             return True  # keep it

#     flt = _PrefixMaxLevelFilter()
#     root.addFilter(flt)
#     try:
#         yield
#     finally:
#         root.removeFilter(flt)
        
@contextmanager
def suppress_logging_warnings(
    logger_name_prefix: str,
    *,
    level: int = logging.WARNING,
    attach_to_handlers: bool = True,
):
    """
    Suppress *logging* records at `level` or higher for any logger whose name
    equals `logger_name_prefix` or starts with `logger_name_prefix + "."`.

    It installs the filter on:
      - the root logger (to catch propagated records),
      - the target logger (to catch records handled locally with propagate=False),
      - and optionally on all current handlers of the target logger.

    This avoids missing records that don't reach the root due to local handling.
    """
    root = logging.getLogger()
    target = logging.getLogger(logger_name_prefix)

    class _PrefixMaxLevelFilter(logging.Filter):
        def filter(self, record: logging.LogRecord) -> bool:
            name = record.name
            targeted = (name == logger_name_prefix) or name.startswith(logger_name_prefix + ".")
            if targeted and record.levelno >= level:
                return False
            return True

    flt = _PrefixMaxLevelFilter()

    # Add to root and target logger
    root.addFilter(flt)
    target.addFilter(flt)

    # Optionally attach to all current handlers on the target
    attached_handlers = []
    if attach_to_handlers:
        for h in target.handlers:
            h.addFilter(flt)
            attached_handlers.append(h)

    try:
        yield
    finally:
        # Remove from root and target logger
        root.removeFilter(flt)
        target.removeFilter(flt)
        # Remove from handlers we touched
        for h in attached_handlers:
            try:
                h.removeFilter(flt)
            except Exception:
                pass

def with_optional_warning_suppression(
    default_logger: str,
    *,
    level=logging.WARNING,
    default_suppress: bool | None = None,
):
    """
    Respects an existing `suppress_warnings` parameter on the function:
    - If the caller *explicitly* passes suppress_warnings, honor it.
    - If not, and default_suppress is True, inject suppress_warnings=True.
    - Opens suppression context only when `suppress_warnings` is True.

    Suppresses *logging* warnings only (via `suppress_logging_warnings`).
    """
    def decorator(func):
        sig = inspect.signature(func)

        @wraps(func)
        def wrapper(*args, **kwargs):
            has_param = "suppress_warnings" in sig.parameters
            caller_provided = has_param and ("suppress_warnings" in kwargs)
            if has_param and not caller_provided and default_suppress is True:
                kwargs["suppress_warnings"] = True

            suppress = bool(kwargs.get("suppress_warnings", False))
            suppress_logger = kwargs.get("suppress_logger", None)
            logger_name = suppress_logger or default_logger

            ctx = suppress_logging_warnings(logger_name, level=level) if suppress else nullcontext()
            with ctx:
                return func(*args, **kwargs)
        return wrapper
    return decorator

# def with_optional_warning_suppression(default_logger: str, *, level=logging.WARNING):
#     """
#     Decorator that wraps a function and adds a `suppress_warnings` kwarg.
#     When `suppress_warnings=True`, it suppresses log >= `level` from `default_logger`
#     (or from `suppress_logger` if provided at call time).

#     Usage:
#         @with_optional_warning_suppression("sandy.warn")
#         def get_errorr(..., suppress_warnings=False, suppress_logger=None):
#             ...

#     Parameters
#     ----------
#     default_logger : str
#         Logger name to suppress if `suppress_logger` not given at call time.
#     level : int, optional
#         Minimum level to suppress. Default: logging.WARNING.
#     """
#     def decorator(func):
#         @wraps(func)
#         def wrapper(*args, suppress_warnings=False, suppress_logger=None, **kwargs):
#             logger_name = suppress_logger or default_logger
#             ctx = suppress_logging_warnings(logger_name, level=level) if suppress_warnings else nullcontext()
#             with ctx:
#                 return func(*args, **kwargs)
#         return wrapper
#     return decorator



# def with_optional_warning_suppression(
#     default_logger: str,
#     *,
#     level: int = logging.WARNING,
#     default_suppress: bool | None = None,
#     # If True, we force suppression unless the caller explicitly sets False.
#     # If False/None, we don't force; we only inject when param missing and default_suppress is True.
# ):
#     """
#     Decorator that:
#       * Uses an existing `suppress_warnings` parameter if provided by caller.
#       * Otherwise, injects `suppress_warnings=default_suppress` if not provided and not positional.
#       * Opens a suppression context only when `suppress_warnings` is True.

#     Args:
#         default_logger: Logger name to suppress by default (e.g., "sandy.warn").
#         level: Logging level cutoff (default: WARNING).
#         default_suppress:
#             - True  -> default to suppressing when the caller didn't provide a value.
#             - False/None -> do not change the effective default (use function’s own default).
#     """
#     def decorator(func):
#         sig = inspect.signature(func)

#         @wraps(func)
#         def wrapper(*args, **kwargs):
#             # Identify if function actually has a `suppress_warnings` parameter
#             has_param = "suppress_warnings" in sig.parameters

#             # Bind arguments (partial to avoid enforcing defaults here)
#             ba = sig.bind_partial(*args, **kwargs)

#             # Determine if caller explicitly set it (positional or kw)
#             caller_provided = has_param and ("suppress_warnings" in ba.arguments)

#             # If not provided by caller and decorator wants default True, set it
#             if has_param and not caller_provided and default_suppress is True:
#                 ba.arguments["suppress_warnings"] = True

#             # Extract effective values after our injection (if any)
#             suppress = bool(ba.arguments.get("suppress_warnings", False))

#             # Allow per-call override of logger name, if the function also supports it
#             suppress_logger = ba.arguments.get("suppress_logger", None)
#             logger_name = suppress_logger or default_logger

#             ctx = suppress_logging_warnings(logger_name, level=level) if suppress else nullcontext()

#             # Reconstruct args/kwargs to call original function
#             # Keep original ordering for positional params; push everything else to kwargs
#             # This keeps behavior consistent even if we injected a kwarg.
#             with ctx:
#                 # Respect the original call form: rebuild args from parameters order
#                 # and fill remaining with kwargs
#                 params = list(sig.parameters.values())
#                 new_args = []
#                 new_kwargs = dict()

#                 # Walk original parameters; if they were passed positionally, keep them
#                 arg_pos = 0
#                 for p in params:
#                     if p.kind in (inspect.Parameter.POSITIONAL_ONLY, inspect.Parameter.POSITIONAL_OR_KEYWORD):
#                         if p.name in ba.arguments and arg_pos < len(args):
#                             # Original positional was provided as positional
#                             new_args.append(ba.arguments[p.name])
#                             arg_pos += 1
#                         elif p.name in ba.arguments and p.name in kwargs:
#                             # Was provided as kw originally
#                             new_kwargs[p.name] = ba.arguments[p.name]
#                         elif p.name in ba.arguments and p.name not in kwargs and arg_pos >= len(args):
#                             # We injected it (or caller gave as kw for a parameter
#                             # after all positionals), keep as kw to avoid shifting positions
#                             new_kwargs[p.name] = ba.arguments[p.name]
#                     else:
#                         # VAR_POSITIONAL / KEYWORD_ONLY / VAR_KEYWORD
#                         if p.name in ba.arguments:
#                             new_kwargs[p.name] = ba.arguments[p.name]

#                 # Also include any extra kwargs that bind_partial accepted but the loop missed
#                 for k, v in ba.arguments.items():
#                     if k not in sig.parameters:
#                         new_kwargs[k] = v

#                 return func(*new_args, **new_kwargs)

#         return wrapper
#     return decorator
