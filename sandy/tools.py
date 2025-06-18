# -*- coding: utf-8 -*-
import os
import time
import ctypes

__author__ = "Luca Fiorito"
__all__ = [
    ]


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
