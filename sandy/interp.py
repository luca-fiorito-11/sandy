# -*- coding: utf-8 -*-
"""
Created on Wed Feb  8 11:43:40 2017

@author: lfiorito
"""

import numpy as np

interp_dict = {1 : 'constant', 2 : 'linear', 4 : 'semilogx', 3 : 'semilogy', 
               5 : 'loglog'}

def check_decorator(foo):
    def checker(xx, yy, axis=0):
        xx = np.array(xx, dtype=float)
        yy = np.array(yy, dtype=float)
        if axis > 0 and np.ndim(yy) <= 1:
            raise ValueError("Cannot interpolate item with dim={} on axis={}".format(np.ndim(yy), axis))
        if xx.size != yy.shape[axis]:
            raise ValueError("Different dimensions for interpolation, check 'axis' value")
        return foo(xx, yy, axis)
    return checker

@check_decorator
def constant(xx, yy, axis=0):
    from scipy.interpolate import interp1d
    interp = interp1d(xx, yy, kind='zero', copy= True, axis=axis,
                      bounds_error=False, fill_value=0, assume_sorted=True)
    return interp

@check_decorator
def linear(xx, yy, axis=0):
    from scipy.interpolate import interp1d
    interp =  interp1d(xx, yy, kind='linear', copy= True, axis=axis,
                           bounds_error=False, fill_value=0, 
                           assume_sorted=True)
    return interp

@check_decorator
def semilogx(xx, yy, axis=0):
    from scipy.interpolate import interp1d
    xx = np.atleast_1d(xx)
    yy = np.atleast_1d(yy)
    if (xx <= 0).any():
        raise ValueError("'x' values must be all positive for semilogx interpolation")
    logx = np.log10(xx)
    interp = interp1d(logx, yy, 
                      kind='linear', copy= True, axis=axis,
                      bounds_error=False, fill_value=0, assume_sorted=True)
    def fun(zz):
        zz = np.atleast_1d(zz)
        if (zz <= 0).any():
            raise ValueError("'x' values must be all positive for semilogx interpolation")
        logz = np.log10(zz)
        mask_domain = (xx[0] <= zz) & (zz <= xx[-1])
        vals = interp(logz[mask_domain])
        if axis == 0:
            dims = (zz.size, yy.shape[1]) if len(yy.shape) > 1 else (zz.size,)
            out = np.zeros(dims, float)
            out[mask_domain] = vals
        elif axis == 1:
            dims = (yy.shape[0], zz.size)
            out = np.zeros(dims, float)
            out[:,mask_domain] = vals
        return out
    return fun

@check_decorator
def semilogy(xx, yy, axis=0):
    from scipy.interpolate import interp1d
    xx = np.atleast_1d(xx)
    yy = np.atleast_1d(yy)
    if (yy <= 0).any():
        raise ValueError("'y' values must be all positive for semilogy interpolation")
    logy = np.log10(yy)
    interp = interp1d(xx, logy, 
                      kind='linear', copy= True, axis=axis,
                      bounds_error=False, fill_value=0, assume_sorted=True)
    def fun(zz):
        zz = np.atleast_1d(zz)
        mask_domain = (xx[0] <= zz) & (zz <= xx[-1])
        vals = interp(zz[mask_domain])
        vals = np.power(10.0, vals)
        if axis == 0:
            dims = (zz.size, yy.shape[1]) if len(yy.shape) > 1 else (zz.size,)
            out = np.zeros(dims, float)
            out[mask_domain] = vals
        elif axis == 1:
            dims = (yy.shape[0], zz.size)
            out = np.zeros(dims, float)
            out[:,mask_domain] = vals
        return out
    return fun

@check_decorator
def loglog(xx, yy, axis=0):
    from scipy.interpolate import interp1d
    xx = np.atleast_1d(xx)
    yy = np.atleast_1d(yy)
    if (xx <= 0).any():
        raise ValueError("'x' values must be all positive for semilogx interpolation")
    if (yy <= 0).any():
        raise ValueError("'y' values must be all positive for semilogy interpolation")
    logx = np.log10(xx)
    logy = np.log10(yy)
    interp = interp1d(logx, logy,
                      kind='linear', copy= True, axis=axis,
                      bounds_error=False, fill_value=0, assume_sorted=True)
    def fun(zz):
        zz = np.atleast_1d(zz)
        if (zz <= 0).any():
            raise ValueError("'x' values must be all positive for semilogx interpolation")
        logz = np.log10(zz)
        mask_domain = (xx[0] <= zz) & (zz <= xx[-1])
        vals = interp(logz[mask_domain])
        vals = np.power(10., vals)
        if axis == 0:
            dims = (zz.size, yy.shape[1]) if len(yy.shape) > 1 else (zz.size,)
            out = np.zeros(dims, float)
            out[mask_domain] = vals
        elif axis == 1:
            dims = (yy.shape[0], zz.size)
            out = np.zeros(dims, float)
            out[:,mask_domain] = vals
        return out
    return fun

@check_decorator
def zero(xx, yy, axis=0):
    xx = np.atleast_1d(xx)
    yy = np.atleast_1d(yy)
    def fun(zz):
        zz = np.atleast_1d(zz)
        if axis == 0:
            dims = (zz.size, yy.shape[1]) if len(yy.shape) > 1 else (zz.size,)
        elif axis == 1:
            dims = (yy.shape[0], zz.size)
        out = np.zeros(dims, float)
        return out
    return fun



class InterpolationRegion:

    interps = {1 : 'constant',
               2 : 'linear',
               3 : 'semilogy', 
               4 : 'semilogx',
               5 : 'loglog'}
    
    def __init__(self, xx, scheme):
        self.beg = xx[0]
        self.end = xx[-1]
        self.scheme = scheme
        self.foo = eval(self.interps[self.scheme])

    @property
    def scheme(self):
        r"""
        Interpolation scheme according to the `ENDF-6` rules.
        
        It must be:
            - 1 : constant
            - 2 : linear
            - 3 : semilogy
            - 4 : semilogx
            - 5 : loglog
        """
        return self._scheme
    
    @scheme.setter
    def scheme(self, scheme):
        if scheme not in self.interps.keys():
            raise NotImplementedError("Interpolation scheme '{}' was not recognized".format(scheme))
        self._scheme = scheme
    
    def interp(self, xx):
        foo = "interpolation_function"
        if not hasattr(self, foo):
            raise AttributeError("Attribute '{}' was not set for '{}' object".format(foo, self.__class__))
        yy = self.interpolation_function(xx)
        return yy
    
    def set_interp(self, xx, yy):
        from sandy.records import Grid
        _yy = Grid(yy, dim=None, minsize=None)
        _xx = Grid(xx, size=yy.shape[0])
        mask = self.contains(_xx)
        idxs = np.flatnonzero(mask)
        if len(idxs) < 2:
            raise NotImplementedError("Need at least two taulated points for interpolation")
        self.interpolation_function = self.foo(_xx[mask], _yy[mask])
    
    def contains(self, item):
        r"""
        Check docstring of function `sandy.functions.contains`
        """
        from sandy.functions import contains
        return contains(item, self.beg, self.end)