"""
Collection of utilities, functions and classes that are requetsed in all code
components.
"""
import numpy as np
import logging
import os
import re



__author__ = "Luca Fiorito"



def add_delimiter_every_n_characters(string, step, delimiter=" "):
    return delimiter.join(string[i:i+step] for i in range(0, len(string), step))


def add_exp_in_endf6_text(text):
    pattern = re.compile("([0-9\.])([+-])([0-9])")
    return pattern.sub("\g<1>E\g<2>\g<3>", text)


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


def log(msg, verbose=False, level=logging.INFO):
    """Log message at given level; 'verbose' only controls INFO-level output."""
    if level == logging.INFO and not verbose:
        return
    logging.log(level, msg)
