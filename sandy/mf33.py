# -*- coding: utf-8 -*-
"""
Created on Tue Jan 31 14:25:11 2017

@author: lfiorito
"""

import sys
import numpy as np
import logging
from sandy.cov import Cov
from sandy.endf import MT as MTstandard
from sandy.endf import MF as MFstandard



class MF(MFstandard):
    """
    ``MF33`` section of an ``ENDF-6`` file.
    All the corresponding ``MT`` subsections are stored are elements of the 
    dictionary.
    """

    def __init__(self, endf):
        self.mat = endf.line.mat
        self.mf = endf.line.mf
        logging.debug(80*"+")
        if self.mf == 31:
            logging.debug("{:^80}".format(" Processing ENDF-6 section MF31 : COVARIANCES OF FISSION NUBAR "))
        elif self.mf == 33:
            logging.debug("{:^80}".format(" Processing ENDF-6 section MF33 : COVARIANCES OF NEUTRON CROSS SECTIONS "))
        logging.debug(80*"+")
        for mt in sorted(endf.INFO.records[self.mf]):
            self[mt] = MT(endf)
            endf.next() # send
        endf.next() # fend
    
    def __setitem__(self, key, item):
        """
        Filter the covariance sections to set in the dictionary according to 
        the following rules:
            * do not set if covariance is empty or belongs to lumped 
              covariances (``MT851`` to ``MT870``)
            * do not set if the reaction is one component of the 
              lumped reaction
            * do not set if ``MT`` is not requested
            * do not set if section is empty
            * do not set if the ``key`` is not in the ``item`` dictionary.
        """
        from sandy.sandy_input import options
        if key in range(851,871):
            logging.warn(item.prefix + "Delete section because SANDY does not process lumped covariances")
            return
        if item.mtl != 0:
            logging.warn(item.prefix + "Delete section because SANDY does not process lumped covariances")
            return
        if "mt" in options:
            if key not in options['mt']:
                logging.info(item.prefix + "Delete section because MT was not requested")
                return
        if len(item) == 0:
            logging.warn(item.prefix + "Delete section because it is empty")
            return
        if key not in item:
            logging.warn(item.prefix + "Delete section because covariance for '({0},{0})' was not found".format(key))
            return
        super().__setitem__(key, item)
        
    
    @property
    def union_grid(self):
        r"""
        Union energy grid for the cross sections/nubar in file `MF33`/`MF31`.
        """
        from sandy.functions import union_grid
        return union_grid(*[self[mt][mt].xx for mt in self ])

    def sorted_items(self):
        r"""
        Generator that yields covariance instances sorted by ```MT`` and 
        ``MTi``.
        
        Outputs:
            - :``key``: :
                (tuple) keyword ``(MT,MTi)``
            - :``item``: :
                (``mf33.Cov33`` instance) covariance object
        """
        mtlist = sorted(self.keys())
        for x1 in mtlist:
            mtilist = sorted(self[x1].keys())
            for x2 in mtilist:
                item = self[x1][x2]
                if x2 not in self:
                    logging.info("MF{0} MT{1} MTi{2} : Skip section because covariances for MT{2} were not found".format(self.mf, x1, x2))
                    continue
                else:
                    key = (x1,x2)
                    yield key, item
    
    def sampling(self):
        """
        Generate a union covariance matrix of ``MF33``.
        Then, extract sample from the union covariance matrix and assign them 
        according to their ``MT`` key.
        """
        from sandy.cov import Ucov
        from sandy.sandy_input import options
        nsmp = options['samples']
        logging.info(self.prefix + "Start sampling procedure")
        data = Ucov({ k : v for k,v in self.sorted_items() })
        if not data:
            logging.error(self.prefix + "No covariance available")
            sys.exit()
        logging.debug(self.prefix + "Draw samples from union covariance matrix")
        for k,s in data.sampling(nsmp).items():
            self[k].smp = s
        if "stdmax" in options:
            self.filter_samples(options['stdmax'])
    
    def filter_samples(self, threshold):
        """
        Delete samples which have standard deviation above a given threshold.
        
        Inputs:
            - :``threshold``: :
                (scalar) threshold of standard deviation
        """
        for mt,sec in sorted(self.items()):
            emask, smask = sec.smp.stdmax(threshold) # This method deletes the samples
            logging.info(sec.prefix + "Stdev > {:.1f}% in {} energy points".format(threshold*100., len(emask)))
            for x,s in zip(emask, smask):
                logging.debug(" - E = {:5e} eV --> stdev = {:.1f}%".format(x, s*100.))
    
    def write_to_csv(self, filename):
        """
        Write covariance matrices to several ``csv`` file, each correspondig to 
        a different ``MT`` reaction.
        
        Inputs:
            - :``filename``: :
                (string) name of the ``csv`` file
        """
        from sandy.sandy_input import outdir
        from os.path import join
        namecov = filename + ".cov_mf{}_mt{}.csv"
        namesmp = filename + ".smp_mf{}_mt{}.csv"
        for mt,sec in sorted(self.items()):
            if hasattr(sec, 'smp'):
                file = join(outdir, namecov.format(self.mf, mt))
                logging.info(sec.prefix + "Writing covariance to file  '{}'".format(file))
                sec[mt].write_to_csv(file)
                file = join(outdir, namesmp.format(self.mf, mt))
                logging.info(sec.prefix + "Writing samples to file     '{}'".format(file))
                sec.smp.write_to_csv(file)
    

        
class MT(dict, MTstandard):
    r"""
    ``MT`` section of file ``MF33``.
    This section contains covariances for cross sections.
    """
    
    def __init__(self, endf):
        """
        Initialize text and control parameters of `MT` section for file `MF33`.
        
        Inputs:
            - endf :
                `ENDF-6` file iterator
        """
        self.mat, self.mf, self.mt, ns = endf.read_control()
        self.text = self.read_text(endf, getback=True)
        self.za, self.awr, self.mtl, self.nl = endf.read_cont([2,4])
        for l in range(self.nl):
            sec = SubSection(endf)
            self[sec.mti] = sec
    
    @property
    def mtl(self):
        r"""
        Non-zero value is used as a flag to indicate that reaction ``MT`` is 
        one component of the evaluator-defined lumped reaction ``MTL``.
        In this case, no covariance information subsections are given for 
        reaction ``MT`` and ``nl=0``.
        """
        return self._mtl
    
    @mtl.setter
    def mtl(self, mtl):
        self._mtl = mtl
    
    @property
    def nl(self):
        r"""
        Number of subsections within a section.
        """
        return self._nl
    
    @nl.setter
    def nl(self, nl):
        logging.debug(self.prefix + "Covariance section contains {} subsections".format(nl))
        self._nl = nl

    def __setitem__(self, key, item):
        r"""
        Apply the following filters before setting an item:
            * the item does not refer to a cross-material section
            * the item does not have an empty list of covariances
            * the item is not a lumped covariance
            * the covariance matrix is not empty
            * ``MT`` is requested
        Optionally, remove cross-reaction correlations if requested.
        """
        from sandy.sandy_input import options
        if item.mati != 0:
            logging.info(item.prefix + "Delete subsection because SANDY does not process cross-material covariances (MAT{},MATi{})  in section MF33".format(self.mat, item.mati))
            return
        if "mt" in options:
            if key not in options['mt'] and key != self.mt:
                logging.debug(item.prefix + "Delete subsection because MT{} was not requested".format(key))
                return
        if 'cross_correlations' in options:
            if not options['cross_correlations'] and key != self.mt:
                logging.info(item.prefix + "Delete cross correlations as requested")
                return
        if not item.covs:
            logging.info(item.prefix + "Delete subsection because no NI-Type covariance was found")
            return
        if key in range(851,871):
            logging.info(item.prefix + "Delete subsection because SANDY does not process lumped covariances")
            return
        if (item.cov == 0).all():
            logging.info(item.prefix + "Delete subsection because covariance matrix is empty")
            return
        super().__setitem__(key, item.cov)



class SubSection:
    r"""
    Subsection of ``MT`` section of ``MF33``.
    It contains the list of ni-type covariances.
    The nc-type covariances are read but not stored.
    """
    
    def __init__(self, endf):
        self.mat, self.mf, self.mt, ns = endf.read_control()
        self.xmfi, self.xlfsi, self.mati, self.mti, self.nc, self.ni = endf.read_cont()
        # READ NC AND NI SECTIONS
        for inc in range(self.nc):
            endf.read_cont()
            endf.read_list()
        self.covs = [ Cov33.read_endf(endf) for ini in range(self.ni) ]
    
    @property
    def covs(self):
        r"""
        List of covariance matrices.
        """
        return self._covs
    
    @covs.setter
    def covs(self, covs):
        self._covs = []
        for cov in covs:
            if cov.lb != 0:
                self._covs.append(cov)
            else:
                logging.debug(self.prefix + "Remove covariance section with 'lb=0'")
    
    @property
    def cov(self):
        """
        Union covariance matrix obtained by merging all the ni-type 
        covariances.
        
        It raises error if ``covs`` is an empty list.
        """
        return Cov33.merge(*self.covs)
    
    @property
    def prefix(self):
        return "MF{0.mf} MT{0.mt} MATi{0.mati} MTi{0.mti} : ".format(self)

    @property
    def xmfi(self):
        r"""
        Floating point equivalent of the ``MF`` for the 2nd energy-dependent 
        cross section of the pair, for which the correlation matrix is given.
        If ``mfi=mf``, ``xmfi=0.0`` or blank.
        """
        return self._xmfi

    @xmfi.setter
    def xmfi(self, xmfi):
        self._xmfi = xmfi

    @property
    def xlfsi(self):
        r"""
        Floating point equivalent for the final excited state of the 2nd 
        energy-dependent cross section.
        If ``mfi=10`` then ``xlfsi=10``, else ``xlfsi=10`` or blank.
        """
        return self._xlfsi

    @xlfsi.setter
    def xlfsi(self, xlfsi):
        self._xlfsi = xlfsi

    @property
    def mati(self):
        r"""
        ``MAT`` value of the covariance second component.
        """
        return self._mati

    @mati.setter
    def mati(self, mati):
        _mati = 0 if mati == self.mat else mati
        self._mati = _mati

    @property
    def mti(self):
        r"""
        ``MT`` value of the covariance second component.
        """
        return self._mti

    @mti.setter
    def mti(self, mti):
        self._mti = mti

    @property
    def nc(self):
        r"""
        Number of nc-type sub-subsections.
        """
        return self._nc

    @nc.setter
    def nc(self, nc):
        if nc != 0:
            logging.debug(self.prefix + "{} nc-type covariances were found".format(nc))
        self._nc = nc

    @property
    def ni(self):
        r"""
        Number of ni-type sub-subsections.
        """
        return self._ni

    @ni.setter
    def ni(self, ni):
        if ni != 0:
            logging.debug(self.prefix + "{} ni-type covariances were found".format(ni))
        self._ni = ni



class Cov33(Cov):
    """
    The same as the inherited object *Cov*, but with some extra methods.
        
    The class methods permit to initialize a *Cov33* instance directly when 
    reading the ENDF-6 file or by merging already existing instances.
    
    Method *reshape* exploits the zero-interpolation law inherent in 
    the covariance matrix.
    """
    
    @classmethod
    def merge(cls, *covs):
        """
        Create a new *Cov33* instance by merging multiple *Cov33* instances.

        The new instance has got attributes:
             - _xx :
                 union grid of the *_xx* attributes of the other instances
             - _yy :
                 union grid of the *_yy* attributes of the other instances
             - _cov :
                 sum of the covariance matrices reshaped according to *_xx* 
                 and *_yy*
        
        Inputs:
            - covs :
                positional arguments containing *Cov33* instances
        """
        from sandy.functions import union_grid
        if len(covs) == 0:
            raise ValueError("At least one or more input arguments must be provided")
        for inst in covs:
            if not isinstance(inst, Cov33):
                raise TypeError("All arguments must be of 'Cov33' type")
        uxx = union_grid(*[inst.xx for inst in covs])
        covs = [ cov.reinit(uxx, uxx) for cov in covs ]
        ucov = np.sum([cov for cov in covs], 0)
        return cls(uxx, uxx, ucov)

    @classmethod
    def read_endf(cls, endf):
        """
        Extract data from ``ENDF-6`` file.
        """
        lb = endf.read_cont([0,1,2,4,5])
        mt = endf.read_control()[2]
        logging.debug(" - lb={}".format(lb))
        endf.back()
        if lb == 0 or lb == 1:
            return cls.lb1(endf)
        elif lb == 2:
            return cls.lb2(endf)
        elif lb == 4:
            return cls.lb4(endf)
        elif lb == 5 or lb == 7:
            return cls.lb5(endf)
        elif lb == 6:
            return cls.lb6(endf)
        elif lb == 8:
            return cls.lb8(endf)
        else:
            logging.error("MF33 MT{} : Covariances with 'lb={}' cannot yet be read by SANDY".format(mt, lb))
            sys.exit()

    @classmethod
    def lb1(cls, endf):
        r"""
        Classmethod to initialize the covariance from the text read in `LB=1` 
        format of `ENDF-6`.
        
        ..ENDF-6 description::
            Fractional components correlated only within each :math:`E_k` 
            interval.
            
        ..math::
            rcov(x,y) = f(E_x) \delta(E_x,E_y)

        Inputs:
            - endf : 
                *endf.ENDF* instance
        
        Outputs:
            - inst :
                *Cov33* instance
        
        Section parameters:
            - ne : 
                number of entries in the array {:math:`E_k`} defining 
                (NE-1) energy intervals.
        """
        from sandy.cov import Cov
        _, _, lt, lb, nt, ne, data = endf.read_list()
        xx = yy = data[::2]
        f = data[1::2]
        cov = Cov(np.diag(f))
        inst = cls(xx, yy, cov, lb)
        return inst

    @classmethod
    def lb2(cls, endf):
        r"""
        Classmethod to initialize the covariance from the text read in `LB=2` 
        format of `ENDF-6`.
        
        ..ENDF-6 description::
            Fractional components correlated over all :math:E k intervals.
            
        ..math::
            rcov(x,y) = f(E_x) f(E_y)

        Inputs:
            - endf : 
                *endf.ENDF* instance
        
        Outputs:
            - inst :
                *Cov33* instance
        
        Section parameters:
            - ne : 
                number of entries in the array {:math:`E_k`} defining 
                (NE-1) energy intervals.
        """
        from sandy.cov import Cov
        _, _, lt, lb, nt, ne, data = endf.read_list()
        xx = yy = data[::2]
        f = np.array(data[1::2])
        cov = Cov(f*f.reshape(-1,1))
        return cls(xx, yy, cov, lb)

    @classmethod
    def lb4(cls, endf):
        r"""
        Classmethod to initialize the covariance from the text read in `LB=4` 
        format of `ENDF-6`.

        Inputs
        
        ENDF-6 description
        ==================
        Fractional components correlated over all :math:`E_l` intervals 
        within each :math:`E_k` interval.
        :math:`rcov(x,y) = f(Ex)*\delta(E_x,E_y)*f'(E_x)*f'(E_y):
        """
        from sandy.functions import union_grid
        _, _, lt, lb, nt, ne, data = endf.read_list()
        xx = data[::2][:-lt]  # first energy grid is until LT points to the end
        yy = data[::2][-lt:]  # second energy grid is the last LT points
        fx = np.array(data[1::2][:-lt])
        fy = np.array(data[1::2][-lt:])
        eg = union_grid(xx, yy)
        cov = fy.reshape(-1,1) * fy
        cov = cls(xx, yy, cov).reinit(eg)
        for i,e1 in enumerate(eg):
            for j,e2 in enumerate(eg):
                if j < i:
                    continue
                f = 0
                for k,e in enumerate(xx[:-1]):
                    if e1 >= e and e1 < xx[k+1] and e2 >= e and e2 < xx[k+1]:
                        f = fx[k]
                cov[i,j] = cov[j,i] = cov[i,j] * f
        return cls(xx, yy, cov, lb)

    @classmethod
    def lb5(cls, endf):
        r"""
        Classmethod to initialize the covariance from the text read in `LB=5` 
        format of `ENDF-6`.
        ..`ENDF-6` description::
            Relative covariance matrix of some cross sections averaged over 
            some energy intervals.
            LS=0 Asymmetric matrix, LS=1 Symmetric matrix.
            
        ..math:
            rcov(x,y) = g(Ex,Ey)

        Inputs:
            - endf : 
                *endf.ENDF* instance
        
        Outputs:
            - inst :
                *Cov33* instance
        
        Section parameters:
            - ne : 
                number of entries in the array {:math:`E_k`} defining 
                (NE-1) energy intervals.
            - ls :
                flag stating if the matrix is asymmetric (`ls=0`) or 
                symmetric (`ls=1`).
        """
        from sandy.cov import triu_matrix
        _, _, ls, lb, nt, ne, data = endf.read_list()
        xx = yy = data[:ne]
        f = np.array(data[ne:])
        if ls == 0: # to be tested
            cov = f.reshape(ne-1, ne-1)
        else:
            cov = triu_matrix(f, ne-1)
        # add zero row and column at the end of the matrix
        cov = np.insert(cov, cov.shape[0], [0]*cov.shape[1], axis=0)
        cov = np.insert(cov, cov.shape[1], [0]*cov.shape[0], axis=1)
        return cls(xx, yy, cov, lb)
    
    @classmethod
    def lb6(cls, endf):
        r"""
        Classmethod to initialize the covariance from the text read in ``LB=6`` 
        format of ``ENDF-6``.
        A relative covariance matrix interrelating the cross sections for 
        two different reaction types or materials with (generally) 
        different energy grids for its rows and columns.
        :math:`rcov(x,y) = g'(Ex,Ey)`
        
        Inputs
        
        Outputs
        
        ENDF-6 description
        ==================
        A relative covariance matrix interrelating the cross sections for 
        two different reaction types or materials with (generally) different 
        energy grids for its rows and columns.
        :math:`rcov(x,y) = g'(Ex,Ey)`
        """
        from sandy.functions import union_grid
        _, _, lt, lb, nt, ner, data = endf.read_list()
        nt = len(data)
        nec = (nt - 1)//ner
        er = data[:ner]
        ec = data[ner:ner+nec]
        f = np.array(data[ner+nec:])
        cov = f.reshape(ner-1, nec-1)
        # add zero row and column at the end of the matrix
        cov = np.insert(cov, cov.shape[0], [0]*cov.shape[1], axis=0)
        cov = np.insert(cov, cov.shape[1], [0]*cov.shape[0], axis=1)
        ug = union_grid(er, ec)
        inst = cls(er, ec, cov, lb)
        inst = inst.reinit(ug)
        return inst
    
    @classmethod
    def lb8(cls, endf):
        r"""
        Classmethod to initialize the covariance from the text read in `LB=8` 
        format of `ENDF-6`.
        ..`ENDF-6` description::
            In general, each :math:`f_k` characterizes an uncorrelated 
            contribution to the absolute variance of the indicated cross 
            section averaged over any energy interval (subgroup) that includes 
            a portion of the energy interval :math:`\Delta E_{k}`.
            The variance contribution :math:`Var(X_{jj})` from an `LB=8` 
            sub-subsection to the processed group variance for the energy 
            group :math:`(E_j , E_{j+1} ) is inversely proportional to its 
            width :math:`\Delta E_j` when :math`(E_j , E_{ j+1} )` lies 
            within :math:`(E_k , E_{k+1} )` and is obtained from the relation:
        ..math::
            Var(X_{jj} ) = f_k \Delta E_k/\Delta E_j
        """
        inst = cls.lb1(endf)
        inst.lb = 8
        return inst
        
    def __new__(cls, xx, yy, cov, *args, **kwargs):
        obj = Cov.__new__(cls, cov)
        return obj
        
    def __init__(self, xx, yy, matrix, lb=100):
        self.xx = xx
        self.yy = yy
        self.lb = lb

    @property
    def lb(self):
        r"""
        Flag whose numerical value determines how the covariance matrix is given.
        """
        return self._lb

    @lb.setter
    def lb(self, lb):
        if lb not in [0, 1, 2, 3, 4, 5, 6, 7, 8, 100]:
            raise NotImplementedError("Invalid 'lb' value for 'Cov33' object")
        self._lb = lb

    @property
    def xx(self):
        r"""
        Array of tabulated x-values.
        """
        return self._xx
    
    @xx.setter
    def xx(self, xx):
        from sandy.records import Grid
        self._xx = Grid(xx, size=self.shape[0])

    @property
    def yy(self):
        r"""
        Array of tabulated y-values.
        """
        return self._yy
    
    @yy.setter
    def yy(self, yy):
        from sandy.records import Grid
        self._yy = Grid(yy, size=self.shape[1])
    
    def corr2cov(self, s):
        dim = self.shape[0]
        S = np.repeat(s, dim).reshape(dim, dim)
        cov = self.__class__(self.xx, self.yy, S.T * (self * S), lb=self.lb)
        cov = cov.up2down()
        return cov

    def up2down(self):
        r"""
        Copy the upper triangular part to the lower triangular part.
        
        Outputs:
            - :``cov``: :
                (2d-array) output covariance matrix
        """
        U = np.triu(self)
        L = np.triu(self, 1).T
        cov = self.__class__(self.xx, self.yy, U + L, lb=self.lb)
        return cov
    
    def reinit(self, xx, yy=None):
        r"""
        Reshape covariance using zero interpolation.
        
        Inputs:
            - xx :
                new array of tabulated energy for 1st covariance dimension.
            - yy : 
                new array of tabulated energy for 2nd covariance dimension.
                If not given it is set equal to *xx*.
        
        Outputs :
            - inst :
                new *Cov33* instance with arguments *xx*, *yy* and *cov*, 
                where *cov* is the reshaped covariance matrix.
        """
        from sandy.functions import zero_interp
        xx = np.atleast_1d(xx)
        if yy is None:
            yy = xx
        yy = np.atleast_1d(yy)
        cov = zero_interp(self.xx, self, xx)
        cov = zero_interp(self.yy, cov.T, yy)
        cov = cov.T
        inst = self.__class__(xx, yy, cov, self.lb)
        return inst
    
    @property
    def Var(self):
        """
        Extract the variance array.

        Outputs:
            - var :
                (1d array) variance array
        """
        from sandy.records import Tab1
        var = np.diag(np.array(self))
        if (var < 0).any():
            raise ValueError("Variances must be non-negative")
        return Tab1(self.xx, var)

    @property
    def Std(self):
        """
        Extract the standard deviation array.

        Outputs:
            - std :
                (1d array) standard deviation array
        """
        from sandy.records import Tab1
        std = np.sqrt(np.array(self.Var))
        return Tab1(self.xx, std)
    
    def write_to_csv(self, file, c1='Energy (eV)', c2='Stdev (%)'):
        """
        Write standard deviation and correlation matrix to `csv` file.
        
        Inputs: 
            - file :
                (string) csv file name
        """
        import csv
        todump = [[c1, c2] + self.xx.tolist()]
        for y,s,row in zip(self.yy, self.std*100., self.corr*1000.):
            todump.append([y] + [s] + row.tolist())
        with open(file, 'w') as fout:
            writer = csv.writer(fout, quoting=csv.QUOTE_ALL)
            writer.writerows(todump)

    def plot(self):
        r"""
        Plot covariance matrix as a pseudocolor plot of a 2-D array.
        The colorbar is also added to the figure.
        """
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots()
        X, Y = np.meshgrid(self.xx, self.yy)
        pcm = ax.pcolormesh(X, Y, self.corr, vmin=-1, vmax=1, cmap='bwr')
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_aspect(1) # Height is 0.5 times the width
        # Resize the plot to make space for the colorbar
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, 0.7, box.height])
        # set labels
        ax.set_title('evaluated correlation matrix')
        ax.set_xlabel('energy (eV)')
        ax.set_ylabel('energy (eV)')
        # Plot the colorbar in desired position 
        cbaxes = fig.add_axes([0.85, 0.1, 0.03, 0.8]) 
        plt.colorbar(pcm, cax=cbaxes)
        fig.show()