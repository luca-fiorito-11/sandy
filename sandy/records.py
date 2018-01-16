#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 11 17:25:33 2017

@author: lfiorito
"""
from sandy.endf import ENDF
import numpy as np
import logging



class Grid(np.ndarray):
    r"""
    Tabulated 1d array representing a continuous monotone function such as 
    an energy grid.
    """
    
    def __new__(cls, xx, TYPE=float, dim=1, minsize=2, size=None):
        """
        Inherhit from ``numpy.ndarray``.
        """
        xx = np.atleast_1d(xx).astype(TYPE)
        obj = np.ndarray.__new__(cls, xx.shape, xx.dtype)
        obj[:] = xx
        return obj
    
    def __init__(self, TYPE=float, dim=1, minsize=2, size=None):
        """
        Run check tests to comply with the ``records.Grid`` constraints.
        
        Inputs:
            - :``TYPE``: :
                type of the array
            - :``dim`` :
                (scalar integer) dimension of the array (deafult is 1d)
            - :``minsize``: :
                (scalar integer) minimum number of tabulated points in the 
                array (default is 2)
            - :``size``: :
                (scalar integer) requested array size (by default the check 
                test is not performed)
            
        """
        if dim is not None:
            self.check_dim(dim)
        if minsize is not None:
            self.check_minsize(minsize)
        if size is not None:
            self.check_size(size)

    @property
    def dx(self):
        r"""
        :math:`n-1` intervals of the grid.
        """
        return self[1:] - self[:-1]
    
    def add_extra_points(self, add=[], ieps=True):
        r"""
        Add extra points to grid.
        
        For each point :math:`x` in the existing grid, two extra points are 
        added respectively at :math:`x - \epsilon` and :math:`x + \epsilon`.
        Additional points to add can be given in input.
        
        Inputs:
            - :``add``: :
                (list) additional points
            - :``ieps``: :
                (boolean) flag to add the extra points at :math:`\pm \eps`
        
        Outputs:
            - :``inst``: :
                (``records.Grid`` instance) grid with additional points
        """
        from math import log10, floor
        extra_points = [ x for x in add ]
        if ieps:
            for i,x in enumerate(self):
                if x <= 0:
                    eps = 1e-5
                    if self[i+1] < 1e-5:
                        continue
                else:
                    exp = floor(log10(x))
                    eps = 10**(exp-6)
                if i != 0:
                    extra_points.append(x-eps)
                if i != len(self)-1:
                    extra_points.append(x+eps)
        points = np.unique(self.tolist() + extra_points)
        inst = self.__class__(points)
        return inst
            
    
    def check_dim(self, dim):
        r"""
        Check that the tabulated array dimension is consistent with the given 
        argument, otherwise raise error.
        
        Inputs:
            - :``dim``: :
                (scalar integer) dimension of the array
        """
        if np.ndim(self) != dim:
            raise NotImplementedError("grid array must have dimension equal to {}".format(dim))
    
    def check_minsize(self, minsize):
        r"""
        Check that the size of the tabulated array is larger or equal to a 
        minimum value given as argument, otherwise raise error.
        
        Inputs:
            - :``minsize``: :
                (scalar integer) minimum number of tabulated points in the 
                array
        """
        if self.size < minsize:
            raise NotImplementedError("grid array size must be at least {}".format(minsize))

    def check_size(self, size):
        r"""
        Check that the size of the tabulated array is equal to a value given 
        as argument, otherwise raise error.
        
        Inputs:
            - :``size``: :
                (scalar integer) requested array size
        """
        if self.size != size:
            raise NotImplementedError("grid array size must be {}".format(size))



class List(list):
    
    @classmethod
    def read_endf(cls, endf):
        c1, c2, l1, l2, n1, n2, bb = endf.read_list()
        return cls(bb)
    
    def write_endf(self, c1=0., c2=0., l1=0, l2=0, n2=0, ismp=None):
        text = ENDF.write_list(c1, c2, l1, l2, n2, self)
        return text



class Tab:

    def __init__(self, xx, nbt=None, intscheme=[1], **kwargs):
        self.xx = xx
        self.nbt = nbt
        self.intscheme = intscheme

    @property
    def xx(self):
        r"""
        Array of tabulated x-values.
        """
        return self._xx
    
    @xx.setter
    def xx(self, xx):
        self._xx = Grid(xx)

    @property
    def domain(self):
        r"""
        Domain of the tabulated function along the x-axis.
        """
        return self.xx[0], self.xx[-1]
        
    @property
    def nbt(self):
        r"""
        List of interpolation regions along the first axis.
        """
        return self._nbt
    
    @nbt.setter
    def nbt(self, nbt):
        if not hasattr(nbt, "__len__"):
            if nbt is None:
                nbt = len(self.xx)
            nbt = [nbt]
        if len(nbt) == 0:
            raise NotImplementedError("Number of interpolation regions cannot be zero")
        if nbt[-1] != len(self.xx):
            raise NotImplementedError("Last item of 'nbt' sequence must be equal to to the length of the tabulated function")
        self._nbt = nbt

    @property
    def intscheme(self):
        r"""
        List of interpolation schemes.
        
        Allowed schemes:
            - 1 : constant
            - 2 : linear
            - 3 : semilogy
            - 4 : semilogx 
            - 5 : loglog
        """
        return self._intscheme

    @intscheme.setter
    def intscheme(self, intscheme):
        if not hasattr(intscheme, "__len__"):
            intscheme = [intscheme]
        if len(intscheme) == 0:
            raise NotImplementedError("Number of interpolation schemes cannot be zero")
        if len(intscheme) != len(self.nbt):
            raise NotImplementedError("Number of interpolation schemes must be equal to number of interpolation regions")
        self._intscheme = intscheme
   
    @property
    def regions(self):
        r"""
        List of interpolation regions.
        """
        from sandy.interp import InterpolationRegion as IR
        regions = []
        xbeg = self.xx[0]
        for n,scheme in zip(self.nbt, self.intscheme):
            xend = self.xx[n-1]
            if xbeg == xend:
                raise NotImplementedError("Interpolation over a single point is not allowed")
            mask = (self.xx >= xbeg) & (self.xx <= xend)
            if mask.any():
                region = IR(self.xx[mask], scheme)
                regions.append(region)
            xbeg = xend
        if len(regions) == 0:
            raise NotImplementedError("Found zero interpolation regions")
        return regions

    @property
    def dx(self):
        r"""
        :math:`n-1` intervals :math:`\Delta x`of the x-grid.
        """
        return self.xx.dx
    
    def copy_attributes(self, inst):
        r"""
        Copy to input  instance `inst` all the attributes of `self` instance 
        that are not already present in `inst`.
        
        Inputs:
            - inst :
                (instance) object receiving attributes
        """
        for key in self.__dict__:
            if not hasattr(inst, key):
                attribute = self.__getattribute__(key)
                inst.__setattr__(key, attribute)
    
    def where(self, x):
        r"""
        Get the position of value `x` in attribute `xx`.
        
        Inputs:
            - x :
                (scalar float) searched value

        Outputs:
            - idx :
                (scalar integer) position of `x`in attribute `xx`, raise `error` 
                if not found
        """
        idx = self.xx.tolist().index(x)
        return idx
    
    def contains(self, item):
        r"""
        Check docstring of function `sandy.functions.contains`
        """
        from sandy.functions import contains
        return contains(item, *self.domain)

    def find_nearest(self, x, opt='below'):
        r"""
        """
        from sandy.functions import find_nearest
        ifound, xfound = find_nearest(self.xx, x, opt=opt)
        yfound = self[ifound]
        return ifound, xfound, yfound

    def add_points(self, points, prefix=None):
        """
        Add energy points along the x-axis.
        Do not add points outside of the object domain.
        
        Inputs:
            - :``points``: :
                (iterator or scalar) energy points to be added
            - :``prefix``: :
                (string) prefix for the debug message
        
        Output:
            - :``inst``: :
                (instance) instance with additional points
        """
        from sandy.functions import union_grid
        _points = Grid(points, minsize=1)
        mask = np.in1d(_points, self.xx, invert=True)
        ug = union_grid(_points, self.xx)
        add = _points[mask]
        if len(add) > 0:
            string = prefix if prefix is not None else self.prefix
            logging.info(string + "{} additional energy points to add".format(len(add)))
        for point in add:
            msg = " - E = {:.6e} eV"
            if not self.contains(point):
                msg += " out of domain, do not add"
            logging.debug(msg.format(point))
        mask = self.contains(ug)
        inst = self.reinit(ug[mask])
        return inst

    def get_section(self, xmin=None, xmax=None):
        r"""
        Extract a section of the tabulated array within a given domain.
        
        Inputs:
            - xmin :
                (scalar float) lower bound of the requested domain
            - xmax :
                (scalar float) upper bound of the requested domain
        
        Outputs:
            - inst :
                (`Tab` instance) copy of the `Tab` instance
        """
        _xmin = self.domain[0] if xmin is None else float(xmin)
        _xmax = self.domain[-1] if xmax is None else float(xmax)
        if _xmin >= _xmax:
            raise NotImplementedError("Lower limit 'xmin'={} must be smaller than upper limit 'xmax'={}".format(_xmin, _xmax))
        mask = (_xmin <= self.xx) & (self.xx <= _xmax)
        return mask
    
    def copy(self):
        r"""
        Create a copy of `Tab` instance.
        Copy the array and all its attributes.
        
        Outputs:
            - inst :
                (`Tab` instance) copy of the `Tab` instance
        """
        from copy import deepcopy
        inst = deepcopy(self)
        self.copy_attributes(inst)
        return inst



class Tab1(Tab, np.ndarray):
    r"""
    Initialize 1 or 2d function that can be interpolated on the first 
    dimension.
    """    

    @staticmethod
    def sum_functions(items, subtract=False):
        r"""
        Magic method to sum `Tab1` instances.
        
        The method operates only with instances that use a linear 
        interpolation function over their complete domain - e.g. `PENDF` cross 
        sections.
        
        Inputs:
            - item :
                (`Tab1` instance) instance to sum
            - subtract :
                (boolean) flag to switch from sum to subtraction (default is sum)
        
        Outputs:
            - yy :
                (`Tab1` instance) sum of the instances
        """
        sign = -1 if subtract else 1
        from sandy.functions import union_grid
        if len(items) < 1:
            raise NotImplementedError("Cannot sum less than one or more items")
        if len(items) == 1:
            return items[0]
        for item in items:
            if not isinstance(item, Tab1):
                raise TypeError("Only 'Tab1' instances can be added",)
            if item.intscheme != [2]:
                raise NotImplementedError("Only 'Tab1' with linear-linear interpolation scheme can be added")
        xx = union_grid( *[item.xx for item in items])
        # Loop to sum the 'Tab1' and 'Tab1Smp' items        
        yy = items[0].reinit(xx)
        for item in items[1:]:
            add = item.reinit(xx)
            if yy.ndim > item.ndim:
                yy = yy + sign * add[:,np.newaxis]
            elif yy.ndim < item.ndim:
                yy =  sign * add + yy[:,np.newaxis]
            else:
                yy = yy + sign * add
        return type(yy)(xx, yy, nbt=len(xx), intscheme=[2])
    
    @classmethod
    def SUM(cls, items, subtract=False):
        r"""
        Magic method to sum `Tab1` instances.
        
        The method operates only with instances that use a linear 
        interpolation function over their complete domain - e.g. `PENDF` cross 
        sections.
        
        Inputs:
            - item :
                (`Tab1` instance) instance to sum
            - subtract :
                (boolean) flag to switch from sum to subtraction (default is sum)
        
        Outputs:
            - yy :
                (`Tab1` instance) sum of the instances
        """
        sign = -1 if subtract else 1
        from sandy.functions import union_grid
        if len(items) < 1:
            raise NotImplementedError("Cannot sum less than one or more items")
        if len(items) == 1:
            return items[0]
        for item in items:
            if not isinstance(item, Tab1):
                raise TypeError("Only 'Tab1' instances can be added",)
            if item.intscheme != [2]:
                raise NotImplementedError("Only 'Tab1' with linear-linear interpolation scheme can be added")
        xx = union_grid( *[item.xx for item in items])
        fltr = list(filter(lambda x : hasattr(x, 'nsmp'), items))
        if len(fltr) == 0:
            yy = cls(xx, np.zeros_like(xx), intscheme=[2])
        else:
            nsmp = fltr[0].nsmp
            yy = cls(xx, np.zeros((xx.size, nsmp)), intscheme=[2])
        for item in items:
            # All grid elements are contained in union grid --> no interpolation
            mask = np.in1d(xx, item.xx)
            if yy.ndim == 2 and item.ndim == 1:
                yy[mask] += sign * item[:,np.newaxis]
            else:
                yy[mask] = yy[mask] + sign * item
        return yy

    @classmethod
    def read_endf(cls, endf, c1=None, c2=None, l1=None, l2=None, **kwargs):
        r"""
        Get `Tab1` object from `ENDF-6` file.
        
        Inputs:
            - endf:
                (`endf.ENDF` instance) `ENDF-6` file
        """
        _c1, _c2, _l1, _l2, nr, nx, nbt, intscheme, xx, yy = endf.read_tab1()
        inst = cls(xx, yy, nbt=nbt, intscheme=intscheme)
        if c1 is not None:
            inst.__setattr__(c1, _c1)
        if c2 is not None:
            inst.__setattr__(c2, _c2)
        if l1 is not None:
            inst.__setattr__(l1, _l1)
        if l2 is not None:
            inst.__setattr__(l2, _l2)
        return inst

    def __new__(cls, xx, yy, **kwargs):
        yy = np.atleast_1d(yy)
        if np.ndim(yy) != 1 and np.ndim(yy) != 2:
            raise NotImplementedError("tabulated function must have DIM=1 or DIM=2")
        obj = np.ndarray.__new__(cls, yy.shape, yy.dtype)
        obj[:] = yy
        return obj

    def __init__(self, xx, yy, **kwargs):
        super().__init__(xx, **kwargs)
        if self.xx.size != self.shape[0]:
            raise NotImplementedError("Different size of x- and y-tabulated arrays")

    @property
    def boundaries(self):
        r"""
        Boundaries of the interpolation schemes.
        """
        return Grid(np.unique([ [r.beg, r.end] for r in self.regions ]))

    def integrate(self, xmin, xmax):
        from sandy.functions import find_nearest
        if xmin < self.xx[0]:
            raise NotImplementedError("")
        if xmax > self.xx[-1]:
            raise NotImplementedError("")
        if xmin >= xmax:
            raise NotImplementedError("")
        imin = find_nearest(self.xx, xmin)[0]
        imax = find_nearest(self.xx, xmax, opt="above")[0]
        domain = Grid(self.xx[imin:imax+1])
        domain[0] = xmin
        domain[-1] = xmax
        integr = domain.dx.dot(self[imin:imax])
        return integr

    def interp(self, xx):
        r"""
        Interpolate function on a new grid or point.
        
        Inputs:
            - :``xx``: :
                new grid or point
        
        Outputs:
            - :``yy``: :
                (array or scalar) interpolated function
        """
        is_scalar = np.ndim(xx) < 1
        _xx = Grid(xx, minsize=1)
        yy = np.zeros(_xx.size) if np.ndim(self) == 1 else np.zeros((_xx.size, self.shape[1]))
        for region in self.regions:
            mask = region.contains(_xx)
            if mask.any():
                region.set_interp(self.xx, self)
                yy[mask] = region.interp(_xx[mask])
        # Directly replace the values that are present in self.xx
        # first method
#        for x,y in zip(self.xx, self):
#            if x in _xx:
#                yy[np.in1d(_xx, x)] = y
        # second method
#        found_in_xx = np.in1d(_xx, self.xx)
#        found_in_self = np.in1d(self.xx, _xx)
#        yy[found_in_xx] = self[found_in_self]
        if is_scalar:
            yy = yy[0]
        return yy

    def reinit(self, xx):
        r"""
        Reinitialize function by interpolating on a new grid.
        
        Inputs:
            - xx :
                new grid. It can be an array, a list, a scalar.
        """
        _xx = Grid(xx)
        yy = self.interp(_xx)
        nbt = []; intscheme = []
        for r in self.regions:
            idxs = np.flatnonzero(_xx <= r.end)
            if idxs.size > 0:
                nbt.append(idxs[-1] + 1)
                intscheme.append(r.scheme)
        if len(nbt) == 0:
            raise NotImplementedError("All tabulated points are set outside the interpolation regions")
        nbt[-1] = len(xx)
        inst = self.__class__(_xx, yy, nbt=nbt, intscheme=intscheme)
        self.copy_attributes(inst)
        return inst
    
    def write_endf(self, c1=0, c2=0, l1=0, l2=0, ismp=None):
        r"""
        Write `Tab1` record according to `ENDF-6` format.
        """
        yy = self if ismp is None else self[:,ismp]
        text = ENDF.write_tab1(c1, c2, l1, l2, self.nbt, self.intscheme, 
                               self.xx, yy)
        return text

    def init_smp(self, nsmp, TYPE):
        r"""
        Create a copy of the `self` instance where the tabulated function is 
        repeated for the number of samples requested.
        Also, copy all the attributes contained in `self` to the output instance.
        
        Inputs:
            - nsmp :
                (scalar integer) number of samples
            - TYPE :
                (metaclass) class object that defines the samples
        
        Output:
            - inst :
                (appropriate samples instance) this instance contains a 
                2d-array function
        """
        yy = np.tile(self, (nsmp,1)).T
        inst = TYPE(self.xx, yy, nbt=self.nbt, intscheme=self.intscheme)
        self.copy_attributes(inst)
        return inst
  
    def plot(self, ax, **kwargs):
        from sandy.plot import plot1d
        plot1d(ax, self.xx, self, **kwargs)
    
    def fig(self, xscale="log", yscale="log", xlabel=None, ylabel=None, 
            legend=False, **kwargs):
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots()
        self.plot(ax, **kwargs)
        ax.set_xlim(self.domain)
        ax.set_xscale(xscale)
        ax.set_yscale(yscale)
        if xlabel:
            ax.set_xlabel(xlabel)
        if ylabel:
            ax.set_ylabel(ylabel)
        ax.grid()
        if legend:
            ax.legend(loc='best', fancybox=True, shadow=True)
        fig.show()

    @staticmethod
    def figs(*tabs, xscale='log', yscale='log', 
             xlabel=None, ylabel=None, title=None, 
             labels=[],
             name=None):
        r"""
        """
        from matplotlib.colors import BASE_COLORS
        import sandy.sandy_input
        import matplotlib.pyplot as plt
        from os.path import join
        from sandy.plot import save
        set_legend = True if len(labels) > 0 else False
        if len(labels) < len(tabs):
            labels += [None]*(len(tabs)-len(labels))
        fig, ax = plt.subplots()
        colors = list(BASE_COLORS.keys())
        for i,tab in enumerate(tabs):
            tab.plot(ax, color=colors[i], label=labels[i])
        ax.set_xscale(xscale)
        ax.set_yscale(yscale)
        ax.grid()
        if title:
            ax.set_title(title)
        if xlabel:
            ax.set_xlabel(xlabel)
        if ylabel:
            ax.set_ylabel(ylabel)
        if set_legend:
            ax.legend(loc='best', fancybox=True, shadow=True)
        if name:
            figname = join(sandy.sandy_input.outdir, name)
            save(figname, fig, 'pdf')
        else:
            fig.show()


class Tab2(Tab, list):

    @classmethod
    def read_endf(cls, endf, TYPE=Tab1, **kwargs):
        r"""
        """
        mat, mf, mt, ns = endf.read_control()
        _,_,_,_,nrx,nx,nbt,intscheme = endf.read_tab2()
        xx = []
        zz = []
        for i in range(nx):
            xx.append(endf.read_cont([0,2,3,4,5]))
            endf.back()
            obj = TYPE.read_endf(endf)
            zz.append(obj)
        inst = cls(xx, zz, nbt=nbt, intscheme=intscheme, **kwargs)
        inst.mat = mat
        inst.mf = mf
        inst.mt = mt
        return inst

    def __init__(self, xx, yy, **kwargs):
        super().__init__(xx, **kwargs)
        if not hasattr(yy, "__len__"):
            raise NotImplementedError("tabulated 'Tab2' functions must be an iterable")
        self.extend([y for y in yy])
        if self.xx.size != len(self):
            raise NotImplementedError("Different size of x- and y-tabulated arrays")

    @property
    def mat(self):
        r"""
        ``MAT`` value.
        """
        if not hasattr(self, "_mat"):
            raise NotImplementedError("Cannot use this function because 'mat' was not defined")
        return self._mat

    @mat.setter
    def mat(self, mat):
        self._mat = mat

    @property
    def mf(self):
        r"""
        ``MF`` value.
        """
        if not hasattr(self, "_mf"):
            raise NotImplementedError("Cannot use this function because 'mf' was not defined")
        return self._mf

    @mf.setter
    def mf(self, mf):
        self._mf = mf

    @property
    def mt(self):
        r"""
        ``MT`` value.
        """
        if not hasattr(self, "_mt"):
            raise NotImplementedError("Cannot use this function because 'mt' was not defined")
        return self._mt

    @mt.setter
    def mt(self, mt):
        self._mt = mt

    def items(self):
        r"""
        Generator that yields key :math:`x_i` and item :math:`f(x_i, y)`.
        """
        i = 0
        while i < len(self):
            yield self.xx[i], self[i]
            i = i + 1
    
    def write_endf(self, c1=0, c2=0, l1=0, l2=0, ismp=None):
        r"""
        Write 
        """
        text = ENDF.write_tab2(c1, c2, l1, l2, len(self.nbt), len(self), 
                               self.nbt, self.intscheme)
        for x,fooy in self.items():
            text += fooy.write_endf(c2=x, ismp=ismp)
        return text
    
    def init_smp(self, nsmp, TYPE):
        r"""
        Return a new `Tab2` instance where the tabulated function contains 
        2d-arrays of samples.
        repeated for the number of samples requested.
        Also, copy all the attributes contained in `self` to the output instance.
        
        Inputs:
            - nsmp :
                (scalar integer) number of samples
            - TYPE :
                (metaclass) class object that defines the samples
        
        Output:
            - inst :
                (`Tab2` instance)
        """
        yy = [ item.init_smp(nsmp, TYPE) for item in self ]
        inst = self.__class__(self.xx, yy, nbt=self.nbt, intscheme=self.intscheme)
        self.copy_attributes(inst)
        return inst

    def interp(self, xx):
        r"""
        Interpolate function :math:`f(x,y)` over the first dimension using an 
        input point :math:`x_i`.
        Return the interpolated output function.

        ..Important::
            `records.Tab2` cannot extrapolate

        Inputs:
            - xx :
                (1d-arry) new x-grid

        Outputs:
            - ff : 
                (list) list of objects representing :math:`f(x_i,y)`
        """
        from sandy.functions import union_grid
        is_scalar = np.ndim(xx) < 1
        _xx = Grid(xx, minsize=1)
        if not self.contains(_xx).all():
            raise NotImplementedError("Elements in the new x-grid are out of domain")
        ff = []
        for x in _xx:
            for region in self.regions:
                if not region.contains(x):
                    continue
                ifound, xi_found, fi_found = self.find_nearest(x, 'below')
                if xi_found == x:
                    f = fi_found
                else:
                    mytype = type(fi_found)
                    jfound = ifound + 1
                    xj_found, fj_found = self.xx[jfound], self[jfound]
                    yy = union_grid(fi_found.xx, fj_found.xx)
                    obj = Tab1([xi_found, xj_found], [fi_found.interp(yy), fj_found.interp(yy)], 
                                nbt=2, intscheme=region.scheme)
                    f = mytype(yy, obj.interp(x), nbt=len(yy), intscheme=fi_found.intscheme)
                ff.append(f)
        if len(ff) != len(_xx):
            raise NotImplementedError("The interpolated x and y values have different dimensions")
        if is_scalar:
            ff = ff[0]
        return ff

    def reinit(self, xx):
        r"""
        Reinitialize function by interpolating on a new grid.
        Always add the interpolation range boundaries to the new grid.
        
        ..Important::
            This method copies all the attributes of the `self` instance to 
            the output instance, except for those already present.
        
        Inputs:
            - xx :
                new grid. It can be an array, a list, a scalar, an empty list...
        """
        _xx = Grid(xx)
        yy = self.interp(_xx)
        nbt = []; intscheme = []
        for r in self.regions:
            idxs = np.flatnonzero(_xx <= r.end)
            if idxs.size > 0:
                nbt.append(idxs[-1] + 1)
                intscheme.append(r.scheme)
        if len(nbt) == 0:
            raise NotImplementedError("All tabulated points are set outside the interpolation regions")
        nbt[-1] = len(xx)
        inst = self.__class__(_xx, yy, nbt=nbt, intscheme=intscheme)
        self.copy_attributes(inst)
        return inst
    


class Samples(Tab1):

#    @classmethod
#    def SUM(cls, items, nsmp, subtract=False):
#        r"""
#        Magic method to sum `Samples` (and `Tab1`) instances.
#        
#        The method operates only with instances that use a linear 
#        interpolation function over their complete domain - e.g. `PENDF` cross 
#        sections.
#        
#        Inputs:
#            - item :
#                (list of instances) `Samples` or `Tab1` instances to sum
#            - subtract :
#                (boolean) flag to switch from sum to subtraction (default is sum)
#        
#        Outputs:
#            - yy :
#                (`Samples` instance) sum of the instances
#        """
#        from sandy.functions import union_grid
#        sign = -1 if subtract else 1
#        if len(items) < 1:
#            raise NotImplementedError("Cannot sum less than one or more items")
#        if len(items) == 1:
#            return items[0]
#        for item in items:
#            if not isinstance(item, Tab1):
#                raise TypeError("Only 'Tab1' instances can be added",)
#            if item.intscheme != [2]:
#                raise NotImplementedError("Only 'Tab1' with linear-linear interpolation scheme can be added")
#        xx = union_grid( *[item.xx for item in items])
#        yy = cls(xx, np.zeros((xx.size,nsmp)), intscheme=[2])
#        # Loop to sum the 'Tab1' and 'Tab1Smp' items        
#        for item in items:
#            # All grid elements are contained in union grid --> no interpolation
#            mask = np.in1d(xx, item.xx)
#            yy[mask] += sign * item
#        return yy
    
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        if self.ndim != 2:
            raise NotImplementedError("A 'Sample' object must have (ndim=2), not (ndim={})".format(self.ndim))

    @property
    def nsmp(self):
        r"""
        Number of samples.
        """
        return self.shape[1]
        
    @property
    def Mean(self):
        """
        Return 1d-array of sample mean values.
        """
        from sandy.records import Tab1
        mean = np.array(self).mean(axis=-1)
        return Tab1(self.xx, mean)
    
    @property
    def Std(self):
        """
        Return 1d-array of sample standard deviations.
        """
        from sandy.records import Tab1
        std = np.array(self).std(axis=-1)
        return Tab1(self.xx, std)
    
    @property
    def Rstd(self):
        """
        Return 1d-array of sample relative standard deviations.
        """
        from sandy.records import Tab1
        from sandy.functions import div0
        rstd = div0(self.Std, np.abs(self.Mean))
        return Tab1(self.xx, rstd)
    
    @property
    def Var(self):
        """
        Return 1d-array of sample variances.
        """
        from sandy.records import Tab1
        var = np.array(self).var(axis=-1)
        return Tab1(self.xx, var)

    @property
    def Cov(self):
        """
        Return sample covariance matrix.
        """
        from sandy.mf33 import Cov33
        cov = np.cov(np.array(self))
        return Cov33(self.xx, self.xx, cov)

    @property
    def Rcov(self):
        """
        Return sample relative covariance matrix.
        """
        from sandy.mf33 import Cov33
        cov = np.cov(np.array(self.Rsmp))
        return Cov33(self.xx, self.xx, cov)

    @property
    def Corr(self):
        """
        Return 2d-array of sample correlation matrix.
        """
        from sandy.mf33 import Cov33
        with np.errstate(divide='ignore', invalid='ignore'):
            corr = np.corrcoef(np.array(self))
            corr[ ~ np.isfinite( corr )] = 0
        return Cov33(self.xx, self.xx, corr)
    
    @property
    def Rsmp(self):
        """
        Return 2d-array of sample correlation matrix.
        """
        from sandy.functions import div0
        rsmp = div0(np.array(self), self.Mean.reshape(-1,1))
        return self.__class__(self.xx, rsmp)
    
    def replace_negative(self):
        """
        Replace negative samples with zeros.
        """
        mask = self < 0
        self[mask] = 0.

    def add_mean(self, mean, relative=True):
        r"""
        Add the mean values to samples.
        
        Inputs:
            - mean :
                (1d-array) array of mean values
            - relative :
                (boolean) if True samples are assumed to be in relative units, 
                otherwise False (default True)
        """
        if relative:
            self += 1
            self *= mean.reshape(-1,1)
        else:
            self += mean.reshape(-1,1)
            self[mean==0] = 0 # Set samples to zero when the mean is zero

    def stdmax(self, rthresh, rstd=False):
        """
        Set to zero samples for which the standard deviation is 
        above a threshold.
        
        A large standard deviation, e.g. 50% or more, can have strong 
        non-physical effects on the mean of the samples.
        This effects could be prevented by removing the samples corresponding 
        to large standard deviations.
        
        Inputs:
            - rthresh :
                (scalar) standard deviation threshold above which the 
                samples are rejected
            - rstd :
                (boolean) if `False` calculate absolute standard deviation, 
                otherwise calculate relative standard deviation (default is 
                `False`).
        """
        std = self.Rstd if rstd else self.Std 
        mask = std > rthresh
        emask = self.xx[mask]
        smask = std[mask]
        self[mask] = 0
        return emask, smask

    def write_to_csv(self, file, relative=True):
        """
        Write relative or absolute perturbations to `csv` file.
        
        Inputs: 
            - file :
                (string) csv file name
            - relative :
                (boolean) flag to define if perturbations are in absolute or 
                relative units (default is relative)
        """
        import csv
        prefix = 'Relative' if relative else 'Absolute'
        todump = [[prefix + ' perturbations']]
        line = ['Energy (eV)'] + self.xx.tolist()
        todump.append(line)
        for i,row in enumerate(self.T):
            line = ['smp_{}'.format(i+1)] + row.tolist()
            todump.append(line)
        with open(file, 'w') as fout:
            writer = csv.writer(fout, quoting=csv.QUOTE_ALL)
            writer.writerows(todump)
