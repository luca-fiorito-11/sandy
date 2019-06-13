# -*- coding: utf-8 -*-
"""
"""
import pdb
import os
import time
import ctypes
import h5py
import sys

import numpy as np

__author__ = "Luca Fiorito"
__all__ = []

def save_dict_to_hdf5(dic, filename):
    """
    ....
    """
    with h5py.File(filename, 'w') as h5file:
        recursively_save_dict_contents_to_group(h5file, '/', dic)

def recursively_save_dict_contents_to_group(h5file, path, dic):
    """
    ....
    """
    for key, item in dic.items():
        if isinstance(item, (np.ndarray, np.int64, np.float64, str, bytes, int, float)):
            h5file[path + str(key)] = item
        elif isinstance(item, dict):
            recursively_save_dict_contents_to_group(h5file, path + str(key) + '/', item)
        else:
            raise ValueError('Cannot save %s type'%type(item))

def load_dict_from_hdf5(filename):
    """
    ....
    """
    with h5py.File(filename, 'r') as h5file:
        return recursively_load_dict_contents_from_group(h5file, '/')

def recursively_load_dict_contents_from_group(h5file, path):
    """
    ....
    """
    ans = {}
    for key, item in h5file[path].items():
        try:
            kdict = int(key)
        except ValueError:
            try:
                kdict = float(key)
            except ValueError:
                kdict = key
        if isinstance(item, h5py._hl.dataset.Dataset):
            ans[kdict] = item[()]
        elif isinstance(item, h5py._hl.group.Group):
            ans[kdict] = recursively_load_dict_contents_from_group(h5file, path + key + '/')
    return ans

def str2bool(v):
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

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

def TimeDecorator(foo):
    """
    Output the time a function takes to execute.
    """
    def wrapper(*args, **kwargs):
        t1 = time.time()
        out = foo(*args, **kwargs)
        t2 = time.time()
        print("Time to run function {}: {} sec".format(foo, t2-t1))
        return out
    return wrapper

def mkl_get_max_threads():
    mkl_rt = ctypes.CDLL('libmkl_rt.so')
    return mkl_rt.mkl_get_max_threads()

def mkl_set_num_threads(cores):
    mkl_rt = ctypes.CDLL('libmkl_rt.so')
    return mkl_rt.mkl_set_num_threads(ctypes.byref(ctypes.c_int(cores)))

def query_yes_no(question, default="yes"):
    """
    Ask a yes/no question via `input()` and return their answer.

    Parameters
    ----------
    question : `srt`
        string that is presented to the user.
    default : `str`, optional, default is `"yes"`
        it is the presumed answer if the user just hits <Enter>.
        It must be "yes" (the default), "no" or None (meaning
        an answer is required of the user).

    Returns
    -------
    `bool`
        The "answer" return value is `True` for `"yes"` or `False` for `"no"`.
    
    Raises
    ------
    `ValueError`
        if `default` is not a valid option
    """
    valid = {"yes": True, "y": True, "ye": True,
             "no": False, "n": False}
    if default is None:
        prompt = " [y/n] "
    elif default == "yes":
        prompt = " [Y/n] "
    elif default == "no":
        prompt = " [y/N] "
    else:
        raise ValueError("invalid default answer: '%s'" % default)
    while True:
        sys.stdout.write(question + prompt)
        choice = input().lower()
        if default is not None and choice == '':
            return valid[default]
        elif choice in valid:
            return valid[choice]
        else:
            sys.stdout.write("Please respond with 'yes' or 'no' "
                             "(or 'y' or 'n').\n")