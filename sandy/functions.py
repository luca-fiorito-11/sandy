# -*- coding: utf-8 -*-
"""
Created on Thu Jan 12 11:10:49 2017

@author: lfiorito
"""
import numpy as np
import os, pdb

def run_process(cmd, cwd=None, timeout=3600, verbose=True):
    import subprocess as sp
    process = sp.Popen("exec " + cmd,
                       shell=True,
                       cwd=cwd,
                       stdin=None,
                       stdout=sp.PIPE,
                       stderr=sp.PIPE,)
    try:
        stdoutdata, stderrdata = process.communicate(timeout=timeout)
    except sp.TimeoutExpired as exc:
        process.kill()
        stdoutdata, stderrdata = process.communicate()
        stderrdata += (r"'{}' took longer than {} seconds to complete and it was killed".format(cmd, timeout)).encode()
    stdout = stdoutdata.decode('utf-8', errors='ignore').rstrip()
    stderr = stderrdata.decode('utf-8', errors='ignore').rstrip()
    if verbose:
        print(stdout)
    returncode = process.returncode
    if returncode: print("process failed to run\n".format(returncode) + stderr)
    return process.returncode, stdout, stderr

def force_symlink(file1, file2):
    try:
        os.symlink(file1, file2)
    except FileExistsError:
        os.remove(file2)
        os.symlink(file1, file2)

def zero_interp(xx_old, arr, xx_new):
    """
    Reshape array on its first dimension by using zero interpolation.
    If some x-values of the new grid exceed the bounds of the old grid, set
    the corresponing interpolated array to zero in those positions.

    Inputs:
        - xx_old :
            old array of tabulated x-values.
        - xx_new :
            new array of tabulated x-values.

    Outputs:
        - arr_new :
            array reshaped according to xx_new.
    """
    arr = np.array(arr)
    if np.ndim(arr) == 1:
        arr_new = zero_interp_1d(xx_old, arr, xx_new)
    elif np.ndim(arr) == 2:
        arr_new = zero_interp_2d(xx_old, arr, xx_new)
    return arr_new

def zero_interp_1d(xx_old, arr, xx_new):
    """
    Reshape array on its first dimension by using zero interpolation.
    If some x-values of the new grid exceed the bounds of the old grid, set
    the corresponing interpolated array to zero in those positions.

    Inputs:
        - xx_old :
            old array of tabulated x-values.
        - xx_new :
            new array of tabulated x-values.

    Outputs:
        - arr_new :
            array reshaped according to xx_new.
    """
    arr = np.atleast_1d(arr)
    xx_old  = np.atleast_1d(xx_old)
    xx_new  = np.atleast_1d(xx_new)
    arr = arr.tolist()
    arr_new = []
    for x in xx_new:
        if x < xx_old[0] or x > xx_old[-1]:
            arr_new.append(0)
        else:
            mask = xx_old <= x
            xfound = xx_old[mask][-1]
            idx = xx_old.tolist().index(xfound)
            arr_new.append(arr[idx])
    return np.array(arr_new)

def zero_interp_2d(xx_old, arr, xx_new):
    """
    Reshape array on its first dimension by using zero interpolation.
    If some x-values of the new grid exceed the bounds of the old grid, set
    the corresponing interpolated array to zero in those positions.

    Inputs:
        - xx_old :
            old array of tabulated x-values.
        - xx_new :
            new array of tabulated x-values.

    Outputs:
        - arr_new :
            array reshaped according to xx_new.
    """
    arr = np.atleast_2d(arr)
    xx_old  = np.atleast_1d(xx_old)
    xx_new  = np.atleast_1d(xx_new)
    d1, d2 = arr.shape
    arr = arr.tolist()
    arr_new = []
    for x in xx_new:
        if x < xx_old[0] or x > xx_old[-1]:
            arr_new.append([0]*d2)
        else:
            mask = xx_old <= x
            xfound = xx_old[mask][-1]
            idx = xx_old.tolist().index(xfound)
            arr_new.append(arr[idx])
    return np.array(arr_new)


def div0( a, b , value=0):
    """
    Ignore division by zero.

    Inputs:
        - a :
            (array or scalar) numerator
        - b :
            (array or scalar) denominator
        - value :
            (scalar) standard replacer when division by zero is found (default is zero)

    Outputs:
        - c :
            (array or scalar) a/b

    ..Example::

        div0( [-1, 0, 1], 0 ) -> [0, 0, 0]
    """
    with np.errstate(divide='ignore', invalid='ignore'):
        c = np.true_divide( a, b )
        if hasattr(c, "__len__"):
            c[ ~ np.isfinite( c )] = value  # -inf inf NaN
        else:
            if not np.isfinite( c ):
                c = value
    return c

def log10( x, value=0):
    """
    Ignore division by zero.

    Inputs:
        - a :
            numerator
        - b :
            denominator
        - value :
            standard replacer when division by zero is found

    Outputs:
        - c: a/b

    ..Example::

        div0( [-1, 0, 1], 0 ) -> [0, 0, 0]
    """
    with np.errstate(divide='ignore', invalid='ignore'):
        c = np.log10(x)
        c[ ~ np.isfinite( c )] = value  # -inf inf NaN
    return c

def split_by_n( seq, n ):
    """
    Description
    ===========
    A generator to divide a sequence into chunks of n units.

    Inputs
    ======
     - *seq*: sequence to be divided;
     - *n*: size of each unit.
    """
    while seq:
        yield seq[:n]
        seq = seq[n:]


def isnum(item, rais=True):
    """
    Return true if input *item* is a scalar number, otherwise false.
    Raise error if input *rais* is true.
    """
    import numbers
    if isinstance(item, numbers.Number):
        return True
    else:
        if rais:
            raise TypeError("Input 'item' must be a scalar number")
        else:
            return False

def printProgress (iteration, total, prefix = '', suffix = '',
                   decimals = 1, barLength = 100):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        barLength   - Optional  : character length of bar (Int)
    """
    import sys
    formatStr       = "{0:." + str(decimals) + "f}"
    percents        = formatStr.format(100 * (iteration / float(total)))
    filledLength    = int(round(barLength * iteration / float(total)))
    bar             = '█' * filledLength + '-' * (barLength - filledLength)
    sys.stdout.write('\r%s |%s| %s%s %s' % (prefix, bar, percents, '%', suffix)),
    sys.stdout.flush()
    if iteration == total:
        sys.stdout.write('\n')
        sys.stdout.flush()



def union_grid(*xgrids):
    """
    Return union grid from a list of grids.

    Inputs:
        - :``xgrids``: :
            (list) input grids

    Outputs:
        - :``ug``: :
            (``records.Grid`` instance) sorted union grid in 1D-array
    """
    from sandy.records import Grid
    ug = Grid(np.unique(np.concatenate([xx for xx in xgrids])))
    if len(ug) == 0:
        raise NotImplementedError("The union grid is empty")
    return ug



def find_nearest(my_array, my_value, opt='below'):
    array = np.array(my_array)
    value = float(my_value)
    idx = (np.abs(array-value)).argmin()
    if opt == 'above' and idx < len(array):
        idx += 1
    return idx, array[idx]



def contains(item, xmin, xmax):
    r"""
    Check if item is included in a given domain.

    Inputs:
        - item :
            (scalar or iterable) x-value(s) to check
        - xmin :
            (scalar) lower bound of the domain
        - xmax :
            (scalar) upper bound of the domain

    Outputs:
        - mask:
            (boolean array/scalar) True if item is included, otherwise False
    """
    if xmin > xmax:
        raise ValueError("Lower bound 'xmin' must be smaller than upper bound 'xmax'")
    try:
        float(item)
    except (ValueError,TypeError):
        try:
            list(map(float, item))
        except (ValueError,TypeError) as exc:
            raise ValueError("item must be a number or a list of numbers")
    mask = (item >= xmin) & (item <= xmax)
    return mask