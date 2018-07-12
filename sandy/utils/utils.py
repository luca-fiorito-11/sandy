# -*- coding: utf-8 -*-
"""
Created on Wed Jul 11 15:36:30 2018

@author: Fiorito_L
"""
import os
import time

def which(program):
    """Mimic the behavior of the UNIX 'which' command.     
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
    """Mimic the behavior of the UNIX 'ln -sf' command.     
    """
    try:
        os.symlink(file1, file2)
    except FileExistsError:
        os.remove(file2)
        os.symlink(file1, file2)

def TimeDecorator(foo):
    """Output the time a function takes to execute.
    """
    def wrapper(*args, **kwargs):
        t1 = time.time()
        out = foo(*args, **kwargs)
        t2 = time.time()
        print("Time to run function {}: {} sec".format(foo, t2-t1))
        return out
    return wrapper